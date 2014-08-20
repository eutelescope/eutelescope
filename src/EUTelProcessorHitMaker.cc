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
#include "EUTelGeometryTelescopeGeoDescription.h"

#include "EUTelProcessorHitMaker.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelExceptions.h"
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
std::string EUTelProcessorHitMaker::_hitHistoLocalName          = "HitHistoLocal";
std::string EUTelProcessorHitMaker::_hitHistoTelescopeName      = "HitHistoTelescope";
std::string EUTelProcessorHitMaker::_densityPlotName            = "DensityPlot";
std::string EUTelProcessorHitMaker::_clusterCenterHistoName     = "ClusterCenter";
std::string EUTelProcessorHitMaker::_clusterCenterXHistoName    = "ClusterCenterX";
std::string EUTelProcessorHitMaker::_clusterCenterYHistoName    = "ClusterCenterY";
#endif

EUTelProcessorHitMaker::EUTelProcessorHitMaker () : Processor("EUTelProcessorHitMaker"),
_pulseCollectionName(),
_hitCollectionName(),
_referenceHitCollectionName("referenceHit"),
_referenceHitCollectionVec(),
_wantLocalCoordinates(false),
_referenceHitLCIOFile("reference.slcio"),
_iRun(0),
_iEvt(0),
_conversionIdMap(),
_alreadyBookedSensorID(),
_aidaHistoMap(),
_histogramSwitch(true),
_orderedSensorIDVec()
{

  // modify processor description
  _description =
    "EUTelProcessorHitMaker is responsible to translate cluster centers from the local frame of reference"
    " to the external frame of reference using the GEAR geometry description";


  registerInputCollection(LCIO::TRACKERPULSE,"PulseCollectionName",
                            "Cluster (pulse) collection name",
                            _pulseCollectionName, string ( "" ));

  registerOutputCollection(LCIO::TRACKERHIT,"HitCollectionName",
                            "Hit collection name",
                            _hitCollectionName, string ( "" ));

  registerOptionalParameter("EnableLocalCoordidates","Hit coordinates are calculated in local reference frame of sensor",
                            _wantLocalCoordinates, static_cast<bool> ( true ) );

  registerOptionalParameter("ReferenceCollection","This is the name of the reference hit collection initialized in this processor. This collection provides the reference vector to correctly determine a plane corresponding to a global hit coordiante.",
                            _referenceHitCollectionName, static_cast< string > ( "referenceHit" ) );
 
  registerOptionalParameter("ReferenceHitFile","This is the file where the reference hit collection is stored",
                            _referenceHitLCIOFile, static_cast< string > ( "reference.slcio" ) );
    
}


void EUTelProcessorHitMaker::init() {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

  // Getting access to geometry description
  std::string name("test.root");
  geo::gGeometry().initializeTGeoDescription(name,false);

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

  _histogramSwitch = true;

  DumpReferenceHitDB();
 
#endif

}


void EUTelProcessorHitMaker::DumpReferenceHitDB()
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

  int Nlayers =  geo::gGeometry().nPlanes();
  for(size_t ii = 0 ; ii <  Nlayers; ii++)
  {
    EUTelReferenceHit * refhit = new EUTelReferenceHit();

    int sensorID = geo::gGeometry().sensorZOrderToID( ii );
    refhit->setSensorID( sensorID );
    refhit->setXOffset( geo::gGeometry().siPlaneXPosition( sensorID ) );
    refhit->setYOffset( geo::gGeometry().siPlaneYPosition(sensorID) );
    refhit->setZOffset( geo::gGeometry().siPlaneZPosition(sensorID) + 0.5*geo::gGeometry().siPlaneZSize(sensorID) );
    
    refVec[0] = 0.;
    refVec[1] = 0.;
    refVec[2] = 1.;
  
    double gRotation[3] = { 0., 0., 0.}; // not rotated
    gRotation[0] = geo::gGeometry().siPlaneZRotation(sensorID); // Euler alpha ;
    gRotation[1] = geo::gGeometry().siPlaneYRotation(sensorID); // Euler alpha ;
    gRotation[2] = geo::gGeometry().siPlaneXRotation(sensorID); // Euler alpha ;
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
    }
    if( TMath::Abs( gRotation[1]) > 1e-6 ) 
    {
        _RotatedVector.Rotate(  gRotation[1], _Yaxis ); // in ZX 
    }
    if( TMath::Abs( gRotation[0]) > 1e-6 ) 
    {   
        _RotatedVector.Rotate(  gRotation[0], _Zaxis ); // in XY
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



void EUTelProcessorHitMaker::addReferenceHitCollection(LCEvent *event, std::string referenceHitName="referenceHit")
{ 

  LCCollectionVec * referenceHitCollection = new LCCollectionVec( LCIO::LCGENERICOBJECT );
 
  int Nlayers =  geo::gGeometry().nPlanes();
  for(size_t ii = 0 ; ii <  Nlayers; ii++)
  {
    EUTelReferenceHit * refhit = new EUTelReferenceHit();
    int sensorID = geo::gGeometry().sensorZOrderToID( ii );
    refhit->setSensorID( sensorID );
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


void EUTelProcessorHitMaker::processRunHeader (LCRunHeader * rdr) {


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



  // increment the run counter
  ++_iRun;
}


void EUTelProcessorHitMaker::processEvent (LCEvent * event) {

    ++_iEvt;

    EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event) ;

    if ( evt->getEventType() == kEORE ) {
      streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
      return;
    } else if ( evt->getEventType() == kUNKNOWN ) {
      streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event." << endl;
    }
    
    LCCollectionVec * pulseCollection   = 0;
    LCCollectionVec * hitCollection     = 0;

    try
    {
      pulseCollection = static_cast<LCCollectionVec*> (event->getCollection( _pulseCollectionName ));
    }
    catch (DataNotAvailableException& e  ) 
    {
      streamlog_out  ( MESSAGE2 ) <<  "No input collection " << _pulseCollectionName << " found on event " << event->getEventNumber()
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
    // prepare an encoder for the hit collection
    CellIDEncoder<TrackerHitImpl> idHitEncoder(EUTELESCOPE::HITENCODING, hitCollection);

    CellIDDecoder<TrackerPulseImpl>  clusterCellDecoder(pulseCollection);
    CellIDDecoder<TrackerDataImpl>   cellDecoder( pulseCollection );

    int oldDetectorID = -100;

    double xZero = 0., yZero = 0., zZero = 0. ;
    double xSize = 0., ySize = 0., zThickness = 0.;
    double resolutionX = 0., resolutionY = 0.;
    double xPitch = 0., yPitch = 0.;
    int xNpixels = 0, yNpixels = 0;

    for ( int iCluster = 0; iCluster < pulseCollection->getNumberOfElements(); iCluster++ ) 
    {
 	TrackerPulseImpl * clusterFrame = dynamic_cast<TrackerPulseImpl*> ( pulseCollection->getElementAt( iCluster ) ); // actual cluster
    
	int sensorID    = clusterCellDecoder(clusterFrame)["sensorID"];
	SparsePixelType clusterType = static_cast<SparsePixelType> ( static_cast<int> (clusterCellDecoder(clusterFrame)["type"] ) );

	
        TrackerDataImpl  * channelList  = dynamic_cast<TrackerDataImpl*> ( clusterFrame->getTrackerData() ); // list of pixels ?
		
        EUTelVirtualCluster * cluster   = new EUTelSparseClusterImpl< EUTelGenericSparsePixel > (static_cast<TrackerDataImpl *> ( channelList ));

      // there could be several clusters belonging to the same
      // detector. So update the geometry information only if this new
      // cluster belongs to a different detector.
      // superseeded with above ;; sensorID = cluster->getDetectorID();

      if ( sensorID != oldDetectorID ) 
      {
          oldDetectorID = sensorID;


          // perfect! The full geometry description is now coming from the
          // GEAR interface. Let's keep the finger xed!

          // check if the histos for this sensor ID have been booked
          // already.
          if ( _alreadyBookedSensorID.find( sensorID ) == _alreadyBookedSensorID.end() ) {
            // we need to book now!
            bookHistos( sensorID );
          }

          xZero        = geo::gGeometry().siPlaneXPosition( sensorID ); // mm
          yZero        = geo::gGeometry().siPlaneYPosition( sensorID ); // mm
          zZero        = geo::gGeometry().siPlaneZPosition( sensorID ); // mm

          resolutionX  = geo::gGeometry().siPlaneXResolution( sensorID );// mm
          resolutionY  = geo::gGeometry().siPlaneYResolution( sensorID );// mm

          xSize        = geo::gGeometry().siPlaneXSize ( sensorID );    // mm
          ySize        = geo::gGeometry().siPlaneYSize ( sensorID );    // mm
          zThickness   = geo::gGeometry().siPlaneZSize ( sensorID );    // mm

          xPitch       = geo::gGeometry().siPlaneXPitch( sensorID );    // mm
          yPitch       = geo::gGeometry().siPlaneYPitch( sensorID );    // mm

          xNpixels     = geo::gGeometry().siPlaneXNpixels( sensorID );    // mm
          yNpixels     = geo::gGeometry().siPlaneYNpixels( sensorID );    // mm
      }


      //
      // LOCAL coordinate system !!!!!!
      //

      // get the position of the seed pixel. This is in pixel number.
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
      if ( clusterType == kEUTelBrickedClusterImpl )
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

      cluster->getCenterOfGravityShift( xShift, yShift );
    
      double xCorrection = static_cast<double> (xShift) ;
      double yCorrection = static_cast<double> (yShift) ;

      // rescale the pixel number in millimeter
      double xDet = ( static_cast<double> (xCluSeed) + xCorrection + 0.5 ) * xPitch ;
      double yDet = ( static_cast<double> (yCluSeed) + yCorrection + 0.5 ) * yPitch ;

// check the hack from Havard:
      float xCoG(0.0f), yCoG(0.0f);
      cluster->getCenterOfGravity(xCoG, yCoG);
      xDet = (xCoG + 0.5) * xPitch;
      yDet = (yCoG + 0.5) * yPitch; 

      streamlog_out(DEBUG1) << "cluster[" << setw(4) << iCluster << "] on sensor[" << setw(3) << sensorID 
                            << "] at [" << setw(8) << setprecision(3) << xCoG << ":" << setw(8) << setprecision(3) << yCoG << "]"
                            << " ->  [" << setw(8) << setprecision(3) << xDet << ":" << setw(8) << setprecision(3) << yDet << "]"
                            << endl;
      
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
      string tempHistoName;
      if ( _histogramSwitch ) 
      {
        tempHistoName =  _hitHistoLocalName + "_" + to_string( sensorID );
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
      telPos[0] = xDet - xSize/2. ;
      telPos[1] = yDet - ySize/2. ; 
      telPos[2] =   0.;

      if ( !_wantLocalCoordinates ) {
            // 
            // NOW !!
            // GLOBAL coordinate system !!!
           
        const double localPos[3] = { telPos[0], telPos[1], telPos[2] };
        geo::gGeometry().local2Master( sensorID, localPos, telPos);

      }
          
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
      if ( _histogramSwitch ) 
      {
        tempHistoName = _hitHistoTelescopeName + "_" + to_string( sensorID );
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

            }
#endif

      // create the new hit
      TrackerHitImpl * hit = new TrackerHitImpl;

      hit->setPosition( &telPos[0] );
      float cov[TRKHITNCOVMATRIX] = {0.,0.,0.,0.,0.,0.};
      double resx = resolutionX;
      double resy = resolutionY;
      cov[0] = resx * resx; // cov(x,x)
      cov[2] = resy * resy; // cov(y,y)
      hit->setCovMatrix( cov );
      hit->setType( clusterType  );

      // prepare a LCObjectVec to store the current cluster
      LCObjectVec clusterVec;
       clusterVec.push_back( channelList );

      // add the clusterVec to the hit
      hit->rawHits() = clusterVec;
      
      // Determine sensorID from the cluster data.
      idHitEncoder["sensorID"] =  sensorID ;

      // set the local/global bit flag property for the hit
      idHitEncoder["properties"] = 0; // init
      if (!_wantLocalCoordinates) idHitEncoder["properties"] = kHitInGlobalCoord;

      // store values
      idHitEncoder.setCellID( hit );

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

void EUTelProcessorHitMaker::end() 
{
  streamlog_out ( MESSAGE4 )  << "Successfully finished" << endl;
}

void EUTelProcessorHitMaker::bookHistos(int sensorID) {

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)


  string tempHistoName;
  string basePath = "plane_" + to_string( sensorID ) ;
  AIDAProcessor::tree(this)->mkdir(basePath.c_str());
  basePath = basePath + "/";

  tempHistoName = _hitHistoLocalName + "_" + to_string( sensorID ) ;

  double xMin =  0;
  double xMax =  geo::gGeometry().siPlaneXSize ( sensorID );   

  double yMin =  0;
  double yMax =  geo::gGeometry().siPlaneYSize ( sensorID ); 

  int xNBin =    geo::gGeometry().siPlaneXNpixels ( sensorID );
  int yNBin =    geo::gGeometry().siPlaneYNpixels ( sensorID );


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
  double safetyFactor = 1.2;
  double xPosition =  geo::gGeometry().siPlaneXPosition( sensorID );
  double yPosition =  geo::gGeometry().siPlaneYPosition( sensorID );
  double xSize     =  geo::gGeometry().siPlaneXSize ( sensorID );
  double ySize     =  geo::gGeometry().siPlaneYSize ( sensorID );
  int xBin         =  geo::gGeometry().siPlaneXNpixels( sensorID );
  int yBin         =  geo::gGeometry().siPlaneYNpixels( sensorID );

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




  _alreadyBookedSensorID.insert( sensorID );

#endif // AIDA

}


#endif // GEAR
