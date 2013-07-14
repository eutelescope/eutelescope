// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// built only if GEAR is used
#ifdef USE_GEAR

// ROOT includes:
#include "TVector3.h"

// eutelescope includes ".h"
#include "EUTelHitMaker.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelSparseCluster2Impl.h"
#include "EUTelExceptions.h"
#include "EUTelEtaFunctionImpl.h"
#include "EUTelAlignmentConstant.h"
#include "EUTelReferenceHit.h"


// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IHistogram3D.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerPulseImpl.h>
//#include <TrackerHitImpl2.h>
#include <IMPL/TrackerHitImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>

// system includes <>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>
#include <iostream>
#include <iomanip>
#include <cstdio>

using namespace std;
using namespace marlin;
using namespace gear;
using namespace eutelescope;

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelHitMaker::_hitHistoLocalName          = "HitHistoLocal";
std::string EUTelHitMaker::_hitHistoTelescopeName      = "HitHistoTelescope";
std::string EUTelHitMaker::_densityPlotName            = "DensityPlot";
std::string EUTelHitMaker::_clusterCenterHistoName     = "ClusterCenter";
std::string EUTelHitMaker::_clusterCenterXHistoName    = "ClusterCenterX";
std::string EUTelHitMaker::_clusterCenterYHistoName    = "ClusterCenterY";
std::string EUTelHitMaker::_clusterCenterEtaHistoName  = "ClusterCenterEta";
std::string EUTelHitMaker::_clusterCenterEtaXHistoName = "ClusterCenterEtaX";
std::string EUTelHitMaker::_clusterCenterEtaYHistoName = "ClusterCenterEtaY";
#endif

EUTelHitMaker::EUTelHitMaker () : Processor("EUTelHitMaker"),
_pulseCollectionName(),
_hitCollectionName(),
_etaCollectionNames(),
_preAlignmentCollectionName("preAlignment"),
_preAlignmentCollectionVec(),
_referenceHitCollectionName("referenceHit"),
_referenceHitCollectionVec(),
_etaCorrection(true),
_3DHistoSwitch(true),
_wantLocalCoordinates(false),
_offsetDBFile("offset-db.slcio"),
_cogAlgorithm("NxMPixel"),
_nPixel(9),
_xyCluSize(),
_referenceHitLCIOFile("reference.slcio"),
_iRun(0),
_iEvt(0),
_conversionIdMap(),
_etaMap(),
_etaVersion(0),
_alreadyBookedSensorID(),
_siPlanesParameters(0),
_siPlanesLayerLayout(0),
_aidaHistoMap(),
_histogramSwitch(true),
_orderedSensorIDVec(),
_siOffsetXMap(),
_siOffsetYMap()
{

  // modify processor description
  _description =
    "EUTelHitMaker is responsible to translate cluster centers from the local frame of reference"
    " to the external frame of reference using the GEAR geometry description";

  registerInputCollection(LCIO::TRACKERPULSE,"PulseCollectionName",
                          "Cluster (pulse) collection name",
                          _pulseCollectionName, string ( "" ));

  registerOutputCollection(LCIO::TRACKERHIT,"HitCollectionName",
                           "Hit collection name",
                           _hitCollectionName, string ( "" ));


  registerProcessorParameter("CoGAlgorithm", "Select here how the center of gravity should be calculated.\n"
                             "FULL: using the full cluster\n"
                             "NPixel: using only the first N most significant pixels (set NPixel too)\n"
                             "NxMPixel: using a subframe of the cluster N x M pixels wide (set NxMPixel too).",
                             _cogAlgorithm, string( "NxMPixel" ) );

  registerOptionalParameter("NPixel", "The number of most significant pixels to be used if CoGAlgorithm is \"NPixel\"",
                            _nPixel, static_cast<int>( 9 ) );

  registerOptionalParameter("Enable3DHisto","If true a 3D histo will be filled. It may require large memory",
                            _3DHistoSwitch, static_cast<bool> ( true ) );
  
  registerOptionalParameter("EnableLocalCoordidates","Hit coordinates are calculated in local reference frame of sensor",
                            _wantLocalCoordinates, static_cast<bool> ( false ) );


  vector<int > xyCluSizeExample(2,3);
  registerOptionalParameter("NxMPixel", "The submatrix size to be used for CoGAlgorithm = \"NxMPixel\"",
                            _xyCluSize, xyCluSizeExample, xyCluSizeExample.size());

  vector<string > etaNames;
  etaNames.push_back("xEta");
  etaNames.push_back("yEta");

  registerOptionalParameter("EtaCollectionName",
                            "The name of the collections containing the eta function (x and y respectively)",
                            _etaCollectionNames, etaNames, etaNames.size());

  registerOptionalParameter("EtaSwitch","Enable or disable eta correction",
                            _etaCorrection, static_cast< bool > ( true ) );

  registerOptionalParameter("OffsetDBFile","This is the name of the LCIO file name with the output offset db (add .slcio)",
                            _offsetDBFile, static_cast< string > ( "offset-db.slcio" ) );
 
  registerOptionalParameter("OffsetCollection","This is the name of the preAligment collection determined from hit correlations and stored in the offset db file.",
                            _preAlignmentCollectionName, static_cast< string > ( "preAlignment" ) );
 
  registerOptionalParameter("ReferenceCollection","This is the name of the reference hit collection initialized in this processor. This collection provides the reference vector to correctly determine a plane corresponding to a global hit coordiante.",
                            _referenceHitCollectionName, static_cast< string > ( "referenceHit" ) );
 
  registerOptionalParameter("ReferenceHitFile","This is the file where the reference hit collection is stored",
                            _referenceHitLCIOFile, static_cast< string > ( "reference.slcio" ) );
    
}


void EUTelHitMaker::init() {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

  // transform the algorithm string to small letters
  transform(_cogAlgorithm.begin(), _cogAlgorithm.end(), _cogAlgorithm.begin(), ::tolower);

  // check if Marlin was built with GEAR support or not
#ifndef USE_GEAR

  streamlog_out ( ERROR4 ) <<  "Marlin was not built with GEAR support." << endl
                           <<  "You need to install GEAR and recompile Marlin with -DUSE_GEAR before continue." << endl;

  // I'm thinking if this is the case of throwing an exception or
  // not. This is a really error and not something that can
  // exceptionally happens. Still not sure what to do
  exit(-1);

#else

  // check if the GEAR manager pointer is not null!
  if ( Global::GEAR == 0x0 ) {
    streamlog_out ( ERROR4 ) <<  "The GearMgr is not available, for an unknown reason." << endl;
    exit(-1);
  }

  _siPlanesParameters  = const_cast<SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

  _orderedSensorIDVec.clear();
  _siOffsetXMap.clear();
  _siOffsetYMap.clear();
  for ( int iPlane = 0 ; iPlane < _siPlanesParameters->getSiPlanesNumber() ; ++iPlane ) {
    _orderedSensorIDVec.push_back( _siPlanesLayerLayout->getID( iPlane ) );
//    _siOffsetXMap.push_back( 0.) ;
//    _siOffsetYMap.push_back( 0.) ;
  }

  _histogramSwitch = true;

  DumpReferenceHitDB();
 
#endif

}

void EUTelHitMaker::DumpReferenceHitDB()
{
// create a reference hit collection file (DB)

  LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();
  try {
    lcWriter->open( _referenceHitLCIOFile, LCIO::WRITE_NEW    );
  } catch ( IOException& e ) {
    streamlog_out ( ERROR4 ) << e.what() << endl;
    exit(-1);
  }

  streamlog_out ( MESSAGE5 ) << "Writing to " << _referenceHitLCIOFile << std::endl;

  LCRunHeaderImpl * lcHeader  = new LCRunHeaderImpl;
  lcHeader->setRunNumber( 0 );
  lcWriter->writeRunHeader(lcHeader);
  delete lcHeader;
  LCEventImpl * event = new LCEventImpl;
  event->setRunNumber( 0 );
  event->setEventNumber( 0 );
  LCTime * now = new LCTime;
  event->setTimeStamp( now->timeStamp() );
  delete now;

  LCCollectionVec * referenceHitCollection = new LCCollectionVec( LCIO::LCGENERICOBJECT );
  double refVec[3] ;
  for(size_t ii = 0 ; ii <  _orderedSensorIDVec.size(); ii++)
  {
    EUTelReferenceHit * refhit = new EUTelReferenceHit();
    refhit->setSensorID( _orderedSensorIDVec[ii] );
    refhit->setXOffset( _siPlanesLayerLayout->getSensitivePositionX(ii) );
    refhit->setYOffset( _siPlanesLayerLayout->getSensitivePositionY(ii) );
    refhit->setZOffset( _siPlanesLayerLayout->getSensitivePositionZ(ii) + 0.5*_siPlanesLayerLayout->getSensitiveThickness(ii) );
    
    refVec[0] = 0.;
    refVec[1] = 0.;
    refVec[2] = 1.;
  
    double gRotation[3] = { 0., 0., 0.}; // not rotated
    gRotation[0] = _siPlanesLayerLayout->getLayerRotationXY(ii); // Euler alpha ;
    gRotation[1] = _siPlanesLayerLayout->getLayerRotationZX(ii); // Euler alpha ;
    gRotation[2] = _siPlanesLayerLayout->getLayerRotationZY(ii); // Euler alpha ;
    streamlog_out( DEBUG5 ) << "GEAR rotations: " << gRotation[0] << " " << gRotation[1] << " " <<  gRotation[2] << endl;
    gRotation[0] =  gRotation[0]*3.1415926/180.; // 
    gRotation[1] =  gRotation[1]*3.1415926/180.; //
    gRotation[2] =  gRotation[2]*3.1415926/180.; //

    TVector3 _RotatedVector( refVec[0], refVec[1], refVec[2] );
    TVector3 _Xaxis( 1.0, 0.0, 0.0 );
    TVector3 _Yaxis( 0.0, 1.0, 0.0 );
    TVector3 _Zaxis( 0.0, 0.0, 1.0 );

    if( TMath::Abs( gRotation[2]) > 1e-6 ) 
    {
        _RotatedVector.Rotate(  gRotation[2], _Xaxis ); // in ZY
//        _Zaxis.Rotate(  gRotation[2], _Xaxis  ); // in XY
//        _Yaxis.Rotate(  gRotation[2], _Xaxis  ); // in XY
    }
    if( TMath::Abs( gRotation[1]) > 1e-6 ) 
    {
        _RotatedVector.Rotate(  gRotation[1], _Yaxis ); // in ZX 
//        _Xaxis.Rotate(  gRotation[1], _Yaxis  ); // in XY
//        _Zaxis.Rotate(  gRotation[1], _Yaxis  ); // in XY
    }
    if( TMath::Abs( gRotation[0]) > 1e-6 ) 
    {   
        _RotatedVector.Rotate(  gRotation[0], _Zaxis ); // in XY
//        _Xaxis.Rotate(  gRotation[0], _Zaxis  ); // in XY
//        _Yaxis.Rotate(  gRotation[0], _Zaxis  ); // in XY
    }
 
     refhit->setAlpha( _RotatedVector[0] );
     refhit->setBeta( _RotatedVector[1] );
     refhit->setGamma( _RotatedVector[2] );
     referenceHitCollection->push_back( refhit );
  }
  event->addCollection( referenceHitCollection, _referenceHitCollectionName );

  lcWriter->writeEvent( event );
  delete event;
  lcWriter->close();
}

void EUTelHitMaker::addReferenceHitCollection(LCEvent *event, std::string referenceHitName="referenceHit")
{ 

  LCCollectionVec * referenceHitCollection = new LCCollectionVec( LCIO::LCGENERICOBJECT );
  for(size_t ii = 0 ; ii <  _orderedSensorIDVec.size(); ii++)
  {
    EUTelReferenceHit * refhit = new EUTelReferenceHit();
    refhit->setSensorID( _orderedSensorIDVec[ii] );
    refhit->setXOffset( 0.0 );
    refhit->setYOffset( 0.0 );
    refhit->setZOffset( 0.0 );
 
    refhit->setAlpha( 0.0 );
    refhit->setBeta( 0.0 );
    refhit->setGamma( 0.0 );
    referenceHitCollection->push_back( refhit );
  }
 
  event->addCollection( referenceHitCollection, referenceHitName );

}


void EUTelHitMaker::processRunHeader (LCRunHeader * rdr) {


  auto_ptr<EUTelRunHeaderImpl> header ( new EUTelRunHeaderImpl (rdr) );
  header->addProcessor( type() );


  // this is the right place also to check the geometry ID. This is a
  // unique number identifying each different geometry used at the
  // beam test. The same number should be saved in the run header and
  // in the xml file. If the numbers are different, instead of barely
  // quitting ask the user what to do.

  if ( header->getGeoID() == 0 )
    streamlog_out ( WARNING0 ) <<  "The geometry ID in the run header is set to zero." << endl
                               <<  "This may mean that the GeoID parameter was not set" << endl;


  if ( header->getGeoID() != _siPlanesParameters->getSiPlanesID() ) {
    streamlog_out ( WARNING5 ) <<  "Error during the geometry consistency check: " << endl
                             << "The run header says the GeoID is " << header->getGeoID() << endl
                             << "The GEAR description says is     " << _siPlanesParameters->getSiPlanesID() << endl;
#ifdef EUTEL_INTERACTIVE
    string answer;
    while (true) {
      streamlog_out ( ERROR1 ) << "Type Q to quit now or C to continue using the actual GEAR description anyway [Q/C]" << endl;
      cin >> answer;
      // put the answer in lower case before making the comparison.
      transform( answer.begin(), answer.end(), answer.begin(), ::tolower );
      if ( answer == "q" ) {
        exit(-1);
      } else if ( answer == "c" ) {
        break;
      }
    }
#endif
  }

  // now book histograms plz...
  if ( isFirstEvent() )  book3DHisto();

  // increment the run counter
  ++_iRun;
}


void EUTelHitMaker::processEvent (LCEvent * event) {

   ++_iEvt;


  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event) ;

  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return;
  } else if ( evt->getEventType() == kUNKNOWN ) {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }


    LCCollectionVec * xEtaCollection = 0x0, * yEtaCollection = 0x0;

    if ( _etaCorrection == 1 ) 
    {
      // this means that the user wants to apply the eta correction to
      // the center of gravity, so we need to have the two Eta function
      // collections

      try 
      {
        xEtaCollection = static_cast<LCCollectionVec*> (event->getCollection( _etaCollectionNames[0] )) ;
      }
      catch (DataNotAvailableException& e) 
      {
        streamlog_out ( ERROR1 )  << "The eta collection " << _etaCollectionNames[0] << " is not available" << endl
                                  << "Continuing without eta correction " << endl;
        _etaCorrection = 0;
      }

      try 
      {
        yEtaCollection = static_cast<LCCollectionVec*> (event->getCollection( _etaCollectionNames[1] )) ;
      }
      catch (DataNotAvailableException& e) 
      {
        streamlog_out ( ERROR1 ) << "The eta collection " << _etaCollectionNames[1] << " is not available" << endl
                                 << "Continuing without eta correction " << endl;
        _etaCorrection = 0;
      }
    }


    if ( isFirstEvent() && ( _etaCorrection == 1 )) 
    {

      EUTelEtaFunctionImpl * func = static_cast< EUTelEtaFunctionImpl*> ( xEtaCollection->getElementAt( 0 ) );

      if ( func->getNInt() != 0 ) 
      {
        _etaVersion = 2;
      } else {
        _etaVersion = 1;
      }


      for ( size_t iDet = 0 ; iDet < xEtaCollection->size(); iDet++) 
      {

        EUTelEtaFunctionImpl * xEtaFunc = static_cast<EUTelEtaFunctionImpl*> ( xEtaCollection->getElementAt(iDet) );

        if ( _etaVersion == 2 ) 
        {
          _etaMap[ xEtaFunc->getSensorID() ] = iDet;
        }

      }

    }

    
    LCCollectionVec * pulseCollection   = 0;
    LCCollectionVec * hitCollection     = 0;

    try
    {
      pulseCollection = static_cast<LCCollectionVec*> (event->getCollection( _pulseCollectionName ));
    }
    catch (DataNotAvailableException& e  ) 
    {
      streamlog_out  ( WARNING9 ) <<  "No input collection " << _pulseCollectionName << " found on event " << event->getEventNumber()
                                  << " in run " << event->getRunNumber() << endl;
      return ;
    }


    try
    {
       hitCollection  = static_cast<LCCollectionVec*> (event->getCollection( _hitCollectionName ));
    }
    catch(...)
    {
//     streamlog_out  ( WARNING2 ) <<  "No output collection " << _hitCollectionName << " found on event " << event->getEventNumber()
//                                << " in run " << event->getRunNumber() << " creating new one  " << endl;
       hitCollection = new LCCollectionVec(LCIO::TRACKERHIT);
    }
  
    CellIDDecoder<TrackerPulseImpl>  pulseCellDecoder(pulseCollection);

    int detectorID    = -99; // it's a non sense
    int oldDetectorID = -100;

    int    layerIndex = -99;
    double xZero = 0., yZero = 0., zZero = 0. ;
    double xSize = 0., ySize = 0.;
    double zThickness = 0.;
    double resolution = 0.;
    double xPitch = 0., yPitch = 0.;
    double xPointing[2] = { 1., 0. }, yPointing[2] = { 1., 0. };

    double gRotation[3] = { 0., 0., 0.}; // not rotated

    for ( int iPulse = 0; iPulse < pulseCollection->getNumberOfElements(); iPulse++ ) 
    {
        TrackerPulseImpl     * pulse   = static_cast<TrackerPulseImpl*> ( pulseCollection->getElementAt(iPulse) );
        EUTelVirtualCluster  * cluster;
        ClusterType type = static_cast<ClusterType>(static_cast<int>((pulseCellDecoder(pulse)["type"])));

        
        if ( type == kEUTelDFFClusterImpl ) 
        {

        // digital fixed cluster implementation. Remember it can come from
        // both RAW and ZS data        
            cluster = new EUTelDFFClusterImpl( static_cast<TrackerDataImpl*> (pulse->getTrackerData()) );
        }
        else if ( type == kEUTelFFClusterImpl ) 
        {
           
        // fixed cluster implementation. Remember it can come from
        // both RAW and ZS data
            cluster = new EUTelFFClusterImpl( static_cast<TrackerDataImpl*> (pulse->getTrackerData()) );

        }
        else if ( type == kEUTelSparseClusterImpl ) 
        {

            // ok the cluster is of sparse type, but we also need to know
            // the kind of pixel description used. This information is
            // stored in the corresponding original data collection.

            LCCollectionVec * sparseClusterCollectionVec = dynamic_cast < LCCollectionVec * > (evt->getCollection("original_zsdata"));
            TrackerDataImpl * oneCluster = dynamic_cast<TrackerDataImpl*> (sparseClusterCollectionVec->getElementAt( 0 ));
            CellIDDecoder<TrackerDataImpl > anotherDecoder(sparseClusterCollectionVec);
            SparsePixelType pixelType = static_cast<SparsePixelType> ( static_cast<int> ( anotherDecoder( oneCluster )["sparsePixelType"] ));

            // now we know the pixel type. So we can properly create a new
            // instance of the sparse cluster
        
            if ( pixelType == kEUTelSimpleSparsePixel ) 
            {
                cluster = new EUTelSparseClusterImpl< EUTelSimpleSparsePixel >
                    ( static_cast<TrackerDataImpl *> ( pulse->getTrackerData()  ) );
            } 
            else 
            {
                streamlog_out ( ERROR4 ) << "Unknown pixel type.  Sorry for quitting." << endl;
                throw UnknownDataTypeException("Pixel type unknown");
            }
            
        }
        else if ( type == kEUTelBrickedClusterImpl ) 
        { 
            //!HACK TAKI

            //  bricked cluster implementation.      
            cluster = new EUTelBrickedClusterImpl( static_cast<TrackerDataImpl*> (pulse->getTrackerData()) ); //!HACK TAKI

        } 
        else if ( type == kEUTelAPIXClusterImpl ) 
        {
        
            cluster = new EUTelSparseClusterImpl< EUTelAPIXSparsePixel > ( static_cast<TrackerDataImpl *> ( pulse->getTrackerData()  ) );
            int xsize0,ysize0;
            cluster->getClusterSize(xsize0,ysize0);
        }
        else 
        {
            streamlog_out ( ERROR4 ) <<  "Unknown cluster type. Sorry for quitting" << endl;
            throw UnknownDataTypeException("Cluster type unknown");
        }

      // there could be several clusters belonging to the same
      // detector. So update the geometry information only if this new
      // cluster belongs to a different detector.
      detectorID = cluster->getDetectorID();

      if ( detectorID != oldDetectorID ) 
      {
        oldDetectorID = detectorID;

        // check if this telescope setup has a DUT
        if ( 
                ( _siPlanesParameters->getSiPlanesType() == _siPlanesParameters->TelescopeWithDUT ) 
                &&
                ( _siPlanesLayerLayout->getDUTID() == detectorID ) 
                ) 
        {
          xZero        = _siPlanesLayerLayout->getDUTSensitivePositionX(); // mm
          yZero        = _siPlanesLayerLayout->getDUTSensitivePositionY(); // mm
          zZero        = _siPlanesLayerLayout->getDUTSensitivePositionZ(); // mm
          zThickness   = _siPlanesLayerLayout->getDUTSensitiveThickness(); // mm
          resolution   = _siPlanesLayerLayout->getDUTSensitiveResolution();// mm
          xPitch       = _siPlanesLayerLayout->getDUTSensitivePitchX();    // mm
          yPitch       = _siPlanesLayerLayout->getDUTSensitivePitchY();    // mm
          xSize        = _siPlanesLayerLayout->getDUTSensitiveSizeX();     // mm
          ySize        = _siPlanesLayerLayout->getDUTSensitiveSizeY();     // mm
          xPointing[0] = _siPlanesLayerLayout->getDUTSensitiveRotation1(); // was -1 ;
          xPointing[1] = _siPlanesLayerLayout->getDUTSensitiveRotation2(); // was  0 ;
          yPointing[0] = _siPlanesLayerLayout->getDUTSensitiveRotation3(); // was  0 ;
          yPointing[1] = _siPlanesLayerLayout->getDUTSensitiveRotation4(); // was

          // check if the histos for this sensor ID have been booked
          // already.
          if ( _alreadyBookedSensorID.find( detectorID ) == _alreadyBookedSensorID.end() ) {
            // we need to book now!
            bookHistos( detectorID, true, xEtaCollection, yEtaCollection );
          }

        } else {


          if ( _conversionIdMap.size() != static_cast< unsigned >(_siPlanesParameters->getSiPlanesNumber()) ) {
            // first of all try to see if this detectorID already belong to
            if ( _conversionIdMap.find( detectorID ) == _conversionIdMap.end() ) {
              // this means that this detector ID was not already inserted,
              // so this is the right place to do that
              for ( int iLayer = 0; iLayer < _siPlanesLayerLayout->getNLayers(); iLayer++ ) {
                if ( _siPlanesLayerLayout->getID(iLayer) == detectorID ) {
                  _conversionIdMap.insert( make_pair( detectorID, iLayer ) );
                  break;
                }
              }
            }
          }

          // perfect! The full geometry description is now coming from the
          // GEAR interface. Let's keep the finger xed!

          // check if the histos for this sensor ID have been booked
          // already.
          if ( _alreadyBookedSensorID.find( detectorID ) == _alreadyBookedSensorID.end() ) {
            // we need to book now!
            bookHistos( detectorID, false, xEtaCollection, yEtaCollection );
          }

          layerIndex   = _conversionIdMap[detectorID];
          xZero        = _siPlanesLayerLayout->getSensitivePositionX(layerIndex); // mm
          yZero        = _siPlanesLayerLayout->getSensitivePositionY(layerIndex); // mm
          zZero        = _siPlanesLayerLayout->getSensitivePositionZ(layerIndex); // mm
          zThickness   = _siPlanesLayerLayout->getSensitiveThickness(layerIndex); // mm
          resolution   = _siPlanesLayerLayout->getSensitiveResolution(layerIndex);// mm
          xPitch       = _siPlanesLayerLayout->getSensitivePitchX(layerIndex);    // mm
          yPitch       = _siPlanesLayerLayout->getSensitivePitchY(layerIndex);    // mm
          xSize        = _siPlanesLayerLayout->getSensitiveSizeX(layerIndex);     // mm
          ySize        = _siPlanesLayerLayout->getSensitiveSizeY(layerIndex);     // mm
          xPointing[0] = _siPlanesLayerLayout->getSensitiveRotation1(layerIndex); // was -1 ;
          xPointing[1] = _siPlanesLayerLayout->getSensitiveRotation2(layerIndex); // was  0 ;
          yPointing[0] = _siPlanesLayerLayout->getSensitiveRotation3(layerIndex); // was  0 ;
          yPointing[1] = _siPlanesLayerLayout->getSensitiveRotation4(layerIndex); // was -1 ;

          try
          {
          gRotation[0] = _siPlanesLayerLayout->getLayerRotationXY(layerIndex); // Euler gamma ;
          gRotation[1] = _siPlanesLayerLayout->getLayerRotationZX(layerIndex); // Euler beta  ;
          gRotation[2] = _siPlanesLayerLayout->getLayerRotationZY(layerIndex); // Euler alpha ;

          // input angles are in DEGREEs !!!
          // translate into radians
          gRotation[0] =  gRotation[0]*3.1415926/180.; // 
          gRotation[1] =  gRotation[1]*3.1415926/180.; //
          gRotation[2] =  gRotation[2]*3.1415926/180.; //
          }
          catch(...)
          {
	    streamlog_out ( MESSAGE5 ) << " no sensor rotation is given in the GEAR steering file, assume NONE " << endl;
          }
        }

        if( _preAlignmentCollectionVec != 0 )
        {
           xZero -= _siOffsetXMap[detectorID];
           yZero -= _siOffsetYMap[detectorID];
        }
    
        if (  ( xPointing[0] == xPointing[1] ) && ( xPointing[0] == 0 ) ) {
          streamlog_out ( ERROR4 ) << "Detector " << detectorID << " has a singular rotation matrix. Sorry for quitting" << endl;
        }

        if (  ( yPointing[0] == yPointing[1] ) && ( yPointing[0] == 0 ) ) {
          streamlog_out ( ERROR4 ) << "Detector " << detectorID << " has a singular rotation matrix. Sorry for quitting" << endl;
        }

      }


      //
      // LOCAL coordinate system !!!!!!
      //

      // get the position of the seed pixel. This is in pixel number.
//      int xCluCenter, yCluCenter;
//      cluster->getCenterCoord(xCluCenter, yCluCenter);
      int xCluSeed = 0;
      int yCluSeed = 0;
      cluster->getSeedCoord(xCluSeed, yCluSeed);



      // with the charge center of gravity calculation, we get a shift
      // from the seed pixel center due to the charge distribution. Those
      // two numbers are the correction values in the case the Eta
      // correction is not applied.



      //!HACK TAKI:
      //! In case of a bricked cluster, we have to make sure to get the normal CoG first.
      //! The one without the global seed coordinate correction (caused by pixel rows being skewed).
      //! That one has to be eta-corrected and the global coordinate correction has to be applied on top of that!
      //! So prepare a brickedCluster pointer now:

      EUTelBrickedClusterImpl* p_tmpBrickedCluster = NULL;
      if ( type == kEUTelBrickedClusterImpl )
      {
            p_tmpBrickedCluster = dynamic_cast< EUTelBrickedClusterImpl* >(cluster);
            if (p_tmpBrickedCluster == NULL)
            {
                streamlog_out ( ERROR4 ) << " .COULD NOT CREATE EUTelBrickedClusterImpl* !!!" << endl;
                throw UnknownDataTypeException("COULD NOT CREATE EUTelBrickedClusterImpl* !!!");
            }
      }

      float xShift = 0.;
      float yShift = 0.;

      if ( _cogAlgorithm == "full" )
      {
        if (_etaCorrection == 1 && p_tmpBrickedCluster != NULL)
        {
            p_tmpBrickedCluster->getCenterOfGravityShiftWithOutGlobalSeedCoordinateCorrection(xShift,yShift);
        }
        else cluster->getCenterOfGravityShift( xShift, yShift );
      }
      else if ( _cogAlgorithm == "npixel" )
      {
        if (_etaCorrection == 1 && p_tmpBrickedCluster != NULL)
        {
            p_tmpBrickedCluster->getCenterOfGravityShiftWithOutGlobalSeedCoordinateCorrection(xShift,yShift,_nPixel);
        }
        else
        {
            streamlog_out( DEBUG5 ) << "Center of gravity algorithm: N PIXEL" << endl;
            cluster->getCenterOfGravityShift( xShift, yShift, _nPixel );
        }
      }
      else if ( _cogAlgorithm == "nxmpixel")
      {
        if (_etaCorrection == 1 && p_tmpBrickedCluster != NULL)
        {
            streamlog_out ( WARNING2 ) << " .Bricked Cluster to be asked for CoG Without Global Seed Coord Correction." << endl;
            streamlog_out ( WARNING2 ) << " .This is not implemented for the NxMPixel algo, cause it does not make much sense for a bricked cluster." << endl;
            streamlog_out ( WARNING2 ) << " .Doing FULL instead." << endl;
            p_tmpBrickedCluster->getCenterOfGravityShiftWithOutGlobalSeedCoordinateCorrection(xShift, yShift);
        }
        else
        {
            //will be okay for a brickedClusterImpl! accounted for such a call internally.
            streamlog_out( DEBUG5 ) << "Center of gravity algorithm: NxM PIXEL" << endl;
            cluster->getCenterOfGravityShift( xShift, yShift, _xyCluSize[0], _xyCluSize[1]);
        }
      }
    
      double xCorrection = static_cast<double> (xShift) ;
      double yCorrection = static_cast<double> (yShift) ;
      //!HACK TAKI if (_etaCorrection==1 && p_tmpBrickedCluster != NULL) THEN DO NOT FORGET
      //!     TO APPLY THE GLOBAL SEED COORD CORRECTION ON TOP AFTERWARDS!

      //if etaCorrection is to applied, then it will overwrite the two Corrections in here:
      if ( _etaCorrection == 1 ) {

        EUTelEtaFunctionImpl * xEtaFunc = NULL;
        EUTelEtaFunctionImpl * yEtaFunc = NULL;

        if ( _etaVersion == 1 ) {
          xEtaFunc = static_cast<EUTelEtaFunctionImpl*> ( xEtaCollection->getElementAt(detectorID) );
          yEtaFunc = static_cast<EUTelEtaFunctionImpl*> ( yEtaCollection->getElementAt(detectorID) );
        } else {
          xEtaFunc = static_cast<EUTelEtaFunctionImpl*> ( xEtaCollection->getElementAt( _etaMap[ detectorID ]) );
          yEtaFunc = static_cast<EUTelEtaFunctionImpl*> ( yEtaCollection->getElementAt( _etaMap[ detectorID ]) );
        }

        bool anomalous = false;
        if ( ( xShift >= -0.5 ) && ( xShift <= 0.5 ) ) {
          xCorrection = xEtaFunc->getEtaFromCoG( xShift );
        } else {
          anomalous = true;
        }
        if ( ( yShift >= -0.5 ) && ( yShift <= 0.5 ) ) {
          yCorrection = yEtaFunc->getEtaFromCoG( yShift );
        } else {
          anomalous = true;
        }

        if ( anomalous )
        {
          streamlog_out ( WARNING4 ) << "Found anomalous cluster\n" << ( * cluster ) << endl;
        }

        //! HACK TAKI (see above):
        //! APPLY THE GLOBAL SEED COORD CORRECTION ON TOP NOW
        if (type == kEUTelBrickedClusterImpl)
        {
            int seedX, seedY;
            cluster->getSeedCoord(seedX, seedY);
            if (seedY %2 == 0)
            {
                xCorrection -= 0.5f;
            }
        }

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
        string tempHistoName;
        if ( _histogramSwitch ) {
          tempHistoName =  _clusterCenterEtaHistoName + "_" + to_string( detectorID ) ;
          if ( AIDA::IHistogram2D * histo = dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[ tempHistoName ] )) {
            histo->fill( xCorrection, yCorrection );
          }

          tempHistoName =  _clusterCenterHistoName + "_" + to_string( detectorID ) ;
          if ( AIDA::IHistogram2D * histo = dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[ tempHistoName ] )) {
            histo->fill( xShift, yShift );
          }

          tempHistoName =  _clusterCenterEtaXHistoName + "_" + to_string( detectorID );
          if ( AIDA::IHistogram1D * histo = dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[ tempHistoName ] )) {
            histo->fill( xCorrection );
          }

          tempHistoName = _clusterCenterXHistoName + "_" + to_string( detectorID );
          if ( AIDA::IHistogram1D * histo = dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[ tempHistoName ] )) {
            histo->fill( xShift );
          }

          tempHistoName = _clusterCenterEtaYHistoName + "_" + to_string( detectorID );
          if ( AIDA::IHistogram1D * histo = dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[ tempHistoName ] )) {
            histo->fill( yCorrection );
          }

          tempHistoName = _clusterCenterYHistoName + "_" + to_string( detectorID );
          if ( AIDA::IHistogram1D * histo = dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[ tempHistoName ] )) {
            histo->fill( yShift );
          }
        }
#endif

      }

      // rescale the pixel number in millimeter
//      double xDet = ( static_cast<double> (xCluCenter) + xCorrection + 0.5 ) * xPitch ;
//      double yDet = ( static_cast<double> (yCluCenter) + yCorrection + 0.5 ) * yPitch ;
      double xDet = ( static_cast<double> (xCluSeed) + xCorrection + 0.5 ) * xPitch ;
      double yDet = ( static_cast<double> (yCluSeed) + yCorrection + 0.5 ) * yPitch ;

// check the hack from Havard:
      float xCoG(0.0f), yCoG(0.0f);
      cluster->getCenterOfGravity(xCoG, yCoG);
// SIMON: LIKELY TO BE REMOVED:
/*      if( detectorID == 22 )
      {
        xCoG =-999999.9;
//        yCoG = 999999.9;
        EUTelSparseClusterImpl< EUTelAPIXSparsePixel > *apixCluster = static_cast<EUTelSparseClusterImpl<EUTelAPIXSparsePixel>*> (cluster);
        for (int iPixel=0; iPixel <  apixCluster->size(); iPixel++) 
        {
	   EUTelAPIXSparsePixel apixPixel;
	   apixCluster->getSparsePixelAt(iPixel, &apixPixel);
          if( xCoG < apixPixel.getXCoord() ) xCoG = apixPixel.getXCoord() ;
//          if( yCoG > apixPixel.getYCoord() ) yCoG = apixPixel.getYCoord() ;
        }
      }
  */
      xDet = (xCoG + 0.5) * xPitch;
      yDet = (yCoG + 0.5) * yPitch;

 


      
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
      string tempHistoName;
      if ( _histogramSwitch ) 
      {
        tempHistoName =  _hitHistoLocalName + "_" + to_string( detectorID );
        if ( AIDA::IHistogram2D* histo = dynamic_cast<AIDA::IHistogram2D*>(_aidaHistoMap[ tempHistoName ]) )
        {
            histo->fill(xDet, yDet);
        }
        else 
        {
            streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                      << ".\nDisabling histogramming from now on " << endl;
            _histogramSwitch = false;
        }
      }
#endif


      // 
      // STILL in the LOCAL coordinate system !!!
      // 
       
      double telPos[3];
      // now perform the rotation of the frame of references and put the
      // results already into a 3D array of double to be ready for the
      // setPosition method of TrackerHit
      telPos[0] = xPointing[0] * xDet + xPointing[1] * yDet;
      telPos[1] = yPointing[0] * xDet + yPointing[1] * yDet;

      // now the translation
      // not sure about the sign. At least it is working for the current
      // configuration but we need to double check it
      double sign = 0;
      if      ( xPointing[0] < -0.7 )       sign = -1 ;
      else if ( xPointing[0] > 0.7 )       sign =  1 ;
      else {
          if       ( xPointing[1] < -0.7 )    sign = -1 ;
          else if  ( xPointing[1] > 0.7 )    sign =  1 ;
      }
      //     telPos[0] += xZero - sign * xSize/2;
      telPos[0] +=  ( -1 ) * sign * xSize / 2;   //apply shifts few lines later

      if      ( yPointing[0] < -0.7 )       sign = -1 ;
      else if ( yPointing[0] > 0.7 )       sign =  1 ;
      else {
          if       ( yPointing[1] < -0.7 )    sign = -1 ;
          else if  ( yPointing[1] > 0.7 )    sign =  1 ;
      }
      //      telPos[1] += yZero - sign * ySize/2;
      telPos[1] += ( -1 ) * sign * ySize / 2;   // apply shifts few lines later


      //    telPos[2] = zZero + 0.5 * zThickness;
      telPos[2] = 0.0 ;       // apply shifts few lines later 
      if ( !_wantLocalCoordinates ) {
            // 
            // NOW !!
            // GLOBAL coordinate system !!!
            // 

            // 
            // rotate according to gRotation angles:
            // 

            //      if(  detectorID >= 10 )
            //     test mode: all planes can be tilted  

            _EulerRotation( telPos, gRotation );
            //
            //    finally apply initial shifts:       
            telPos[0] += xZero;
            telPos[1] += yZero;
            telPos[2] += zZero + 0.5 * zThickness;
      }
          
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
      if ( _histogramSwitch ) 
      {
        tempHistoName = _hitHistoTelescopeName + "_" + to_string( detectorID );
        AIDA::IHistogram2D * histo2D = dynamic_cast<AIDA::IHistogram2D*> (_aidaHistoMap[ tempHistoName ] );
        if ( histo2D )
        {
            histo2D->fill( telPos[0], telPos[1] );
        }
        else 
        {
          streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                    << ".\nDisabling histogramming from now on " << endl;
          _histogramSwitch = false;
        }

        if ( _3DHistoSwitch ) 
        {
          AIDA::IHistogram3D * histo3D = dynamic_cast<AIDA::IHistogram3D*> (_aidaHistoMap[ _densityPlotName ] );
          if ( histo3D ) 
          {
              histo3D->fill( telPos[0], telPos[1], telPos[2] );
          }
          else 
          {
            streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                      << ".\nDisabling histogramming from now on " << endl;
            _histogramSwitch = false;
          }
        }
      }
#endif

      // create the new hit
      TrackerHitImpl * hit = new TrackerHitImpl;
//      hit->setDetectorID( detectorID ) ;
      hit->setPosition( &telPos[0] );
      float cov[TRKHITNCOVMATRIX] = {0.,0.,0.,0.,0.,0.};
      double resx = resolution;
      double resy = resolution;
      cov[0] = resx * resx; // cov(x,x)
      cov[2] = resy * resy; // cov(y,y)
      hit->setCovMatrix( cov );
      hit->setType( pulseCellDecoder(pulse)["type"] );

      // prepare a LCObjectVec to store the current cluster
      LCObjectVec clusterVec;
      clusterVec.push_back( pulse->getTrackerData() );

      // add the clusterVec to the hit
      hit->rawHits() = clusterVec;

      // add the new hit to the hit collection
      hitCollection->push_back( hit );
 
      // delete the eutel cluster
      delete cluster;
    }

    try
    { 
      event->getCollection( _hitCollectionName ) ;
    }
    catch(...)
    {
      event->addCollection( hitCollection, _hitCollectionName );
    }

    if ( isFirstEvent() ) _isFirstEvent = false;
}

void EUTelHitMaker::end() 
{
  streamlog_out ( MESSAGE4 )  << "Successfully finished" << endl;
}

void EUTelHitMaker::bookHistos(int sensorID, bool isDUT, LCCollection * xEtaCollection, LCCollection * yEtaCollection) {

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  int layerIndex = 0;
  if ( !isDUT ) {
    layerIndex   = _conversionIdMap[ sensorID ];
  }

  string tempHistoName;
  string basePath = "plane_" + to_string( sensorID ) ;
  AIDAProcessor::tree(this)->mkdir(basePath.c_str());
  basePath = basePath + "/";

  tempHistoName = _hitHistoLocalName + "_" + to_string( sensorID ) ;

  double xMin =  0;
  double xMax =  (isDUT) ? _siPlanesLayerLayout->getDUTSensitiveSizeX() : _siPlanesLayerLayout->getSensitiveSizeX ( layerIndex );

  double yMin =  0;
  double yMax =  (isDUT) ? _siPlanesLayerLayout->getDUTSensitiveSizeY() :_siPlanesLayerLayout->getSensitiveSizeY ( layerIndex );

  int xNBin =  (isDUT) ? _siPlanesLayerLayout->getDUTSensitiveNpixelX() : _siPlanesLayerLayout->getSensitiveNpixelX( layerIndex );
  int yNBin =  (isDUT) ? _siPlanesLayerLayout->getDUTSensitiveNpixelY() : _siPlanesLayerLayout->getSensitiveNpixelY( layerIndex );


  AIDA::IHistogram2D * hitHistoLocal = AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(),
                                                                                                 xNBin, xMin, xMax, yNBin, yMin, yMax );
  if ( hitHistoLocal ) {
    hitHistoLocal->setTitle("Hit map in the detector local frame of reference");
    _aidaHistoMap.insert( make_pair( tempHistoName, hitHistoLocal ) );
  } else {
    streamlog_out ( ERROR1 )  << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                              << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
    _histogramSwitch = false;
  }

  // 2 should be enough because it
  // means that the sensor is wrong
  // by all its size.
  double safetyFactor = 2.0;
  double xPosition = (isDUT) ? _siPlanesLayerLayout->getDUTSensitivePositionX( ) : _siPlanesLayerLayout->getSensitivePositionX( layerIndex );
  double yPosition = (isDUT) ? _siPlanesLayerLayout->getDUTSensitivePositionY( ) : _siPlanesLayerLayout->getSensitivePositionY( layerIndex );
  double xSize     = (isDUT) ? _siPlanesLayerLayout->getDUTSensitiveSizeX ( )    : _siPlanesLayerLayout->getSensitiveSizeX ( layerIndex );
  double ySize     = (isDUT) ? _siPlanesLayerLayout->getDUTSensitiveSizeY ( )    : _siPlanesLayerLayout->getSensitiveSizeY ( layerIndex );
  int xBin         = (isDUT) ? _siPlanesLayerLayout->getDUTSensitiveNpixelX( )   : _siPlanesLayerLayout->getSensitiveNpixelX( layerIndex );
  int yBin         = (isDUT) ? _siPlanesLayerLayout->getDUTSensitiveNpixelY( )   : _siPlanesLayerLayout->getSensitiveNpixelY( layerIndex );

  xMin = safetyFactor * ( xPosition - ( 0.5 * xSize ));
  xMax = safetyFactor * ( xPosition + ( 0.5 * xSize ));

  yMin = safetyFactor * ( yPosition - ( 0.5 * ySize ));
  yMax = safetyFactor * ( yPosition + ( 0.5 * ySize ));

  xNBin = static_cast< int > ( safetyFactor  * xBin );
  yNBin = static_cast< int > ( safetyFactor  * yBin );

  tempHistoName =  _hitHistoTelescopeName + "_" + to_string( sensorID );
  AIDA::IHistogram2D * hitHistoTelescope =
    AIDAProcessor::histogramFactory(this)->createHistogram2D( ( basePath + tempHistoName ).c_str(),
                                                              xNBin, xMin, xMax, yNBin, yMin, yMax );

  if ( hitHistoTelescope ) {
    hitHistoTelescope->setTitle("Hit map in the telescope frame of reference");
    _aidaHistoMap.insert( make_pair ( tempHistoName, hitHistoTelescope ) );
  } else {
    streamlog_out ( ERROR1 )  << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                              << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
    _histogramSwitch = false;
  }


  if ( _etaCorrection && xEtaCollection != NULL && yEtaCollection != NULL ) {

    EUTelEtaFunctionImpl * xEtaFunc = NULL, * yEtaFunc = NULL;

    if ( _etaVersion >= 2 ) {

      xEtaFunc = static_cast< EUTelEtaFunctionImpl * > ( xEtaCollection->getElementAt( _etaMap[ sensorID ] ) );
      yEtaFunc = static_cast< EUTelEtaFunctionImpl * > ( yEtaCollection->getElementAt( _etaMap[ sensorID ] ) );

    } else {

      xEtaFunc = static_cast< EUTelEtaFunctionImpl * > ( xEtaCollection->getElementAt( sensorID  ) );
      yEtaFunc = static_cast< EUTelEtaFunctionImpl * > ( yEtaCollection->getElementAt( sensorID  ) );
    }

      int xNoOfBin = xEtaFunc->getNoOfBin();
      int yNoOfBin = yEtaFunc->getNoOfBin();
      string tempHistoName = _clusterCenterEtaHistoName + "_" + to_string( sensorID ) ;

    AIDA::IHistogram2D * clusterCenterEta =
      AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName ).c_str(),
                                                                1* xNoOfBin, -.5, +.5, 1 * yNoOfBin, -.5, +.5);

    if ( clusterCenterEta ) {
      clusterCenterEta->setTitle("Position of the cluster center (Eta corrected)");
      _aidaHistoMap.insert( make_pair( tempHistoName, clusterCenterEta ) );
    } else {
      streamlog_out ( ERROR1 )  << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
      _histogramSwitch = false;
    }


    tempHistoName = _clusterCenterHistoName + "_" + to_string( sensorID );
    AIDA::IHistogram2D * clusterCenter =
      AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName ).c_str(),
                                                                1* xNoOfBin, -.5, +.5, 1 * yNoOfBin, -.5, +.5);

    if ( clusterCenter ) {
      clusterCenterEta->setTitle("Position of the cluster center");
      _aidaHistoMap.insert( make_pair( tempHistoName, clusterCenter ) );
    } else {
      streamlog_out ( ERROR1 )  << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
      _histogramSwitch = false;
    }


    tempHistoName = _clusterCenterEtaXHistoName + "_" + to_string( sensorID );
    AIDA::IHistogram1D * clusterCenterEtaX =
      AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                1 * xNoOfBin, -0.5, +0.5 );
    if ( clusterCenterEtaX ) {
      clusterCenterEtaX->setTitle("Projection along X of the cluster center (Eta corrected)");
      _aidaHistoMap.insert( make_pair( tempHistoName, clusterCenterEtaX ) );
    } else {
      streamlog_out ( ERROR1 )  << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
      _histogramSwitch = false;
    }

    tempHistoName = _clusterCenterXHistoName + "_" + to_string( sensorID );
    AIDA::IHistogram1D * clusterCenterX =
      AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                1 * xNoOfBin, -0.5, +0.5 );
    if ( clusterCenterX ) {
      clusterCenterX->setTitle("Projection along X of the cluster center");
      _aidaHistoMap.insert( make_pair( tempHistoName, clusterCenterX ) );
    } else {
      streamlog_out ( ERROR1 )  << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
      _histogramSwitch = false;
    }

    tempHistoName = _clusterCenterEtaYHistoName + "_" + to_string( sensorID );

    AIDA::IHistogram1D * clusterCenterEtaY =
      AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                1 * xNoOfBin, -0.5, +0.5 );
    if ( clusterCenterEtaY ) {
      clusterCenterEtaY->setTitle("Projection along Y of the cluster center (Eta corrected)");
      _aidaHistoMap.insert( make_pair( tempHistoName, clusterCenterEtaY ) );
    } else {
      streamlog_out ( ERROR1 )  << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
      _histogramSwitch = false;
    }

    tempHistoName =  _clusterCenterYHistoName + "_" + to_string( sensorID );
    AIDA::IHistogram1D * clusterCenterY =
      AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                1 * xNoOfBin, -0.5, +0.5 );
    if ( clusterCenterY ) {
      clusterCenterY->setTitle("Projection along Y of the cluster center");
      _aidaHistoMap.insert( make_pair( tempHistoName, clusterCenterY ) );
    } else {
      streamlog_out ( ERROR1 )  << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
      _histogramSwitch = false;
    }

    }

  _alreadyBookedSensorID.insert( sensorID );

#endif // AIDA

}

void EUTelHitMaker::book3DHisto() {

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    // we have to found the boundaries of this histograms. Let's take
    // the outer positions in all directions
    double xMin  =      numeric_limits< double >::max();
    double xMax  = -1 * numeric_limits< double >::max();
    int    xNBin = numeric_limits< int >::min();

    double yMin  =      numeric_limits< double >::max();
    double yMax  = -1 * numeric_limits< double >::max();
    int    yNBin = numeric_limits< int >::min();

    for ( int iPlane = 0 ; iPlane < _siPlanesParameters->getSiPlanesNumber(); ++iPlane ) {

      // x axis
      xMin  = min( _siPlanesLayerLayout->getSensitivePositionX( iPlane ) - ( 0.5  * _siPlanesLayerLayout->getSensitiveSizeX( iPlane )), xMin);
      xMax  = max( _siPlanesLayerLayout->getSensitivePositionX( iPlane ) + ( 0.5  * _siPlanesLayerLayout->getSensitiveSizeX( iPlane )), xMax);
      xNBin = max( _siPlanesLayerLayout->getSensitiveNpixelX( iPlane ), xNBin );

      // y axis
      yMin  = min( _siPlanesLayerLayout->getSensitivePositionY( iPlane ) - ( 0.5  * _siPlanesLayerLayout->getSensitiveSizeY( iPlane )), yMin);
      yMax  = max( _siPlanesLayerLayout->getSensitivePositionY( iPlane ) + ( 0.5  * _siPlanesLayerLayout->getSensitiveSizeY( iPlane )), yMax);
      yNBin = max( _siPlanesLayerLayout->getSensitiveNpixelY( iPlane ), yNBin );

    }

    if  ( _siPlanesParameters->getSiPlanesType() == _siPlanesParameters->TelescopeWithDUT ) {
      // x axis
      xMin  = min( _siPlanesLayerLayout->getDUTSensitivePositionX(  ) - ( 0.5  * _siPlanesLayerLayout->getDUTSensitiveSizeX(  )), xMin);
      xMax  = max( _siPlanesLayerLayout->getDUTSensitivePositionX(  ) + ( 0.5  * _siPlanesLayerLayout->getDUTSensitiveSizeX(  )), xMax);
      xNBin = max( _siPlanesLayerLayout->getDUTSensitiveNpixelX(  ), xNBin );

      // y axis
      yMin  = min( _siPlanesLayerLayout->getDUTSensitivePositionY(  ) - ( 0.5  * _siPlanesLayerLayout->getDUTSensitiveSizeY(  )), yMin);
      yMax  = max( _siPlanesLayerLayout->getDUTSensitivePositionY(  ) + ( 0.5  * _siPlanesLayerLayout->getDUTSensitiveSizeY(  )), yMax);
      yNBin = max( _siPlanesLayerLayout->getDUTSensitiveNpixelY(  ), yNBin );

    }

    if ( _3DHistoSwitch ) {
      // since we may still have alignment problem, we have to take a
      // safety factor on the x and y direction especially.
      // here I take something less than 2 because otherwise I will have
      // a 200MB histogram.
      double safetyFactor = 1.2;

      double xDistance = std::abs( xMax - xMin ) ;
      double xCenter   = ( xMax + xMin ) / 2 ;
      xMin  = xCenter - safetyFactor * ( xDistance / 2 );
      xMax  = xCenter + safetyFactor * ( xDistance / 2 );
      xNBin = static_cast< int > ( xNBin * safetyFactor );

      // generate the x axis binning
      vector< double > xAxis;
      double step = xDistance / xNBin;
      for ( int i = 0 ; i < xNBin ; ++i ) {
        xAxis.push_back ( xMin + i * step );
      }

      double yDistance = std::abs( yMax - yMin ) ;
      double yCenter   = ( yMax + yMin ) / 2 ;
      yMin  = yCenter - safetyFactor * ( yDistance / 2 );
      yMax  = yCenter + safetyFactor * ( yDistance / 2 );
      yNBin = static_cast< int > ( yNBin * safetyFactor );

      // generate the y axis binning
      vector< double > yAxis;
      step = yDistance / yNBin;
      for ( int i = 0 ; i < yNBin ; ++i ) {
        yAxis.push_back( yMin + i * step ) ;
      }


      // generate the z axis but not equally spaced!
      double safetyMargin = 10; // this is mm
      vector< double > zAxis;

      vector< double > zPos;
      for ( int i = 0 ; i < _siPlanesParameters->getSiPlanesNumber(); ++i ) {
        zPos.push_back( _siPlanesLayerLayout->getSensitivePositionZ( i ) );
      }

      if ( _siPlanesParameters->getSiPlanesType() == _siPlanesParameters->TelescopeWithDUT ) {
        zPos.push_back(  _siPlanesLayerLayout->getDUTSensitivePositionZ(  ) );
      }

      sort( zPos.begin(), zPos.end() );

      for ( size_t pos = 0; pos < zPos.size(); ++pos ) {
        zAxis.push_back( zPos.at( pos ) - safetyMargin );
        zAxis.push_back( zPos.at( pos ) + safetyMargin );
      }

      AIDA::IHistogram3D * densityPlot = AIDAProcessor::histogramFactory(this)->createHistogram3D( _densityPlotName ,
                                                                                                   "Hit position in the telescope frame of reference",
                                                                                                   xAxis, yAxis, zAxis, "");

      if ( densityPlot ) {
        _aidaHistoMap.insert( make_pair ( _densityPlotName, densityPlot ) ) ;
      } else {
        streamlog_out ( ERROR1 )  << "Problem booking the " << (_densityPlotName) << ".\n"
                                  << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }
    }

#endif // AIDA
}

void EUTelHitMaker::_EulerRotation(double* _telPos, double* _gRotation) {
  
    TVector3 _UnrotatedSensorHit( _telPos[0], _telPos[1], 0. );
    TVector3 _RotatedSensorHit( _telPos[0], _telPos[1], 0. );

    TVector3 _Xaxis( 1.0, 0.0, 0.0 );
    TVector3 _Yaxis( 0.0, 1.0, 0.0 );
    TVector3 _Zaxis( 0.0, 0.0, 1.0 );

    if( TMath::Abs(_gRotation[2]) > 1e-6 ) 
    {
        _RotatedSensorHit.Rotate( _gRotation[2], _Xaxis ); // in ZY
//        _Zaxis.Rotate( _gRotation[2], _Xaxis  ); // in XY
//        _Yaxis.Rotate( _gRotation[2], _Xaxis  ); // in XY
    }
    if( TMath::Abs(_gRotation[1]) > 1e-6 ) 
    {
        _RotatedSensorHit.Rotate( _gRotation[1], _Yaxis ); // in ZX 
//        _Xaxis.Rotate( _gRotation[1], _Yaxis  ); // in XY
//        _Zaxis.Rotate( _gRotation[1], _Yaxis  ); // in XY
    } 
    if( TMath::Abs(_gRotation[0]) > 1e-6 ) 
    {   
        _RotatedSensorHit.Rotate( _gRotation[0], _Zaxis ); // in XY
//        _Xaxis.Rotate( _gRotation[0], _Zaxis  ); // in XY
//        _Yaxis.Rotate( _gRotation[0], _Zaxis  ); // in XY
    }
 
    _telPos[0] = _RotatedSensorHit.X();
    _telPos[1] = _RotatedSensorHit.Y();
    _telPos[2] = _telPos[2] + _RotatedSensorHit.Z();
 
}

#endif // GEAR
