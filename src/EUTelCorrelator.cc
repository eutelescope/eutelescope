// Version $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// since v00-00-09, this processor requires GEAR
// mainly for geometry initialization.

#if defined(USE_GEAR)

// eutelescope includes ".h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelHistogramManager.h"
 
#include "EUTelCorrelator.h"
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

#include <UTIL/LCTime.h>

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <Exceptions.h>


// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"

// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ITree.h>
#include <AIDA/IAxis.h>
#endif

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// lcio includes <.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerRawDataImpl.h>

// system includes <>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cstdio>

using namespace std;
using namespace marlin;
using namespace eutelescope;

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelCorrelator::_clusterXCorrelationHistoName   = "ClusterXCorrelation";
std::string EUTelCorrelator::_clusterYCorrelationHistoName   = "ClusterYCorrelation";
std::string EUTelCorrelator::_hitXCorrelationHistoName       = "HitXCorrelation";
std::string EUTelCorrelator::_hitYCorrelationHistoName       = "HitYCorrelation";

std::string EUTelCorrelator::_clusterXCorrShiftHistoName             = "ClusterXCorrShift";
std::string EUTelCorrelator::_clusterYCorrShiftHistoName             = "ClusterYCorrShift";
std::string EUTelCorrelator::_clusterXCorrShiftProjectionHistoName   = "ClusterXCorrShiftProjection";
std::string EUTelCorrelator::_clusterYCorrShiftProjectionHistoName   = "ClusterYCorrShiftProjection";

std::string EUTelCorrelator::_hitXCorrShiftHistoName             = "HitXCorrShift";
std::string EUTelCorrelator::_hitYCorrShiftHistoName             = "HitYCorrShift";
std::string EUTelCorrelator::_hitXCorrShiftProjectionHistoName   = "HitXCorrShiftProjection";
std::string EUTelCorrelator::_hitYCorrShiftProjectionHistoName   = "HitYCorrShiftProjection";

#endif

EUTelCorrelator::EUTelCorrelator () : Processor("EUTelCorrelator"), 
_histoInfoFileName("histoinfo.xml")
{

  // modify processor description
  _description =
    "EUTelCorrelator fills histograms with correlation plots";

  EVENT::StringVec	      _clusterCollectionVecExample;
  
  registerInputCollections ( LCIO::TRACKERPULSE, "InputClusterCollections",
                            "List of cluster collections",
                            _clusterCollectionVec, _clusterCollectionVecExample);

  registerInputCollection(LCIO::TRACKERHIT,"InputHitCollectionName",
                              "Hit collection name",
                              _inputHitCollectionName, string ( "hit" ) );

  registerProcessorParameter ("ClusterChargeMinimum",
                              "Minimum allowed cluster charge to be taken into account for the correlation plots (default = 2)",
                              _clusterChargeMin, static_cast <int> (2) );

  registerProcessorParameter ("DumpOffset",
                              "Dump the offset X and Y values calculated from the correlation bands (default = true)",
                              _dumpOffset, static_cast <bool> (true) );

  registerProcessorParameter ("Events",
                              "How many events are needed to get reasonable correlation plots (and Offset DB)? (default=1000)",
                              _events, static_cast <int> (1000) );

  registerOptionalParameter ("FixedPlane", "SensorID of fixed plane", _fixedPlaneID, 0);

  registerOptionalParameter("OffsetDBFile","This is the name of the LCIO file name with the output offset db (add .slcio)",
                              _offsetDBFile, static_cast< string > ( "offset-db.slcio" ) );

 
  registerOptionalParameter("ResidualsXMin","Minimal values of the hit residuals in the X direction for a correlation band. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsXMin, std::vector<float > (6, -10.) );

  registerOptionalParameter("ResidualsYMin","Minimal values of the hit residuals in the Y direction for a correlation band. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsYMin, std::vector<float > (6, -10.) );

  registerOptionalParameter("ResidualsXMax","Maximal values of the hit residuals in the X direction for a correlation band. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsXMax, std::vector<float > (6,  10.) );

  registerOptionalParameter("ResidualsYMax","Maximal values of the hit residuals in the Y direction for a correlation band. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsYMax, std::vector<float > (6,  10.) );

  registerOptionalParameter ("MinNumberOfCorrelatedHits",
                             "If there are more then this number of correlated hits (planes->track candidate) (default=5)",
                             _minNumberOfCorrelatedHits, static_cast <int> (5) );

  registerOptionalParameter("HotPixelCollectionName", "This is the name of the hot pixel collection to be saved into the output slcio file",
                             _hotPixelCollectionName, static_cast< string > ( "hotpixel" ));

  registerOptionalParameter("HistogramInfoFilename", "Name of histogram info xml file", _histoInfoFileName, string("histoinfo.xml"));

}


void EUTelCorrelator::init() {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // Getting access to geometry description
  std::string name("test.root");
  geo::gGeometry().initializeTGeoDescription(name,false);


  // clear the sensor ID vector
  _sensorIDVec.clear();

  // clear the sensor ID map
  _sensorIDVecMap.clear();
  _sensorIDtoZOrderMap.clear();


  // clear the sensor ID vector (z-axis order)
  _sensorIDVecZOrder.clear();


  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

  // check if the GEAR manager pointer is not null!
  if ( Global::GEAR == 0x0 ) {
    streamlog_out ( ERROR4 ) <<  "The GearMgr is not available, for an unknown reason." << endl;
    exit(-1);
  }

  _siPlanesParameters  = const_cast<gear::SiPlanesParameters* >  ( &(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*>  ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

  _siPlanesRotations.clear();
  _siPlanesRotations.resize( _siPlanesLayerLayout->getNLayers() );

  _siPlanesPitchX.clear();
  _siPlanesPitchX.resize( _siPlanesLayerLayout->getNLayers() );

  _siPlanesPitchY.clear();
  _siPlanesPitchY.resize( _siPlanesLayerLayout->getNLayers() );

  _siPlanesOffsetX.clear();
  _siPlanesOffsetX.resize( _siPlanesLayerLayout->getNLayers() );

  _siPlanesOffsetY.clear();
  _siPlanesOffsetY.resize( _siPlanesLayerLayout->getNLayers() );


   for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++ ) 
   {
      _siPlanesRotations[iPlane].insert( make_pair( 1, _siPlanesLayerLayout->getSensitiveRotation1(iPlane) ) ); 
      _siPlanesRotations[iPlane].insert( make_pair( 2, _siPlanesLayerLayout->getSensitiveRotation2(iPlane) ) ); 
      _siPlanesRotations[iPlane].insert( make_pair( 3, _siPlanesLayerLayout->getSensitiveRotation3(iPlane) ) ); 
      _siPlanesRotations[iPlane].insert( make_pair( 4, _siPlanesLayerLayout->getSensitiveRotation4(iPlane) ) );

      _siPlanesPitchX[iPlane] = _siPlanesLayerLayout->getSensitivePitchX(iPlane);
      _siPlanesPitchY[iPlane] = _siPlanesLayerLayout->getSensitivePitchY(iPlane);
 
      _siPlanesOffsetX[iPlane] = 0.;
      _siPlanesOffsetY[iPlane] = 0.;
   }
   
 
   _siPlaneZPosition = new double[ _siPlanesLayerLayout->getNLayers() ];
   for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++ ) 
   {
    _siPlaneZPosition[ iPlane ] = _siPlanesLayerLayout->getLayerPositionZ(iPlane);
    int sensorID = _siPlanesLayerLayout->getID( iPlane );

    _sensorIDVec.push_back( sensorID );
    _sensorIDVecMap.insert( make_pair( sensorID, iPlane ) );

    // count number of the sensors to the left of the current one:
    int _sensors_to_the_left = 0;
    for ( int jPlane = 0 ; jPlane < _siPlanesLayerLayout->getNLayers(); jPlane++ ) 
    {
        if( _siPlanesLayerLayout->getLayerPositionZ(jPlane) + 1e-06 < _siPlaneZPosition[ iPlane ] )
        {
            _sensors_to_the_left++;
        }
    }
 
    _sensorIDVecZOrder.push_back( _sensors_to_the_left );
    _sensorIDtoZOrderMap.insert(make_pair( sensorID, _sensors_to_the_left));
 
    _minX[ sensorID ] = 0;
    _minY[ sensorID ] = 0;
    _maxX[ sensorID ] = _siPlanesLayerLayout->getSensitiveNpixelX( iPlane ) - 1;
    _maxY[ sensorID ] = _siPlanesLayerLayout->getSensitiveNpixelY( iPlane ) - 1;

    if(                
            _siPlanesRotations[iPlane][1] ==  0.                     &&
            _siPlanesRotations[iPlane][2] !=  0.                     &&
            _siPlanesRotations[iPlane][3] !=  0.                     &&
            _siPlanesRotations[iPlane][4] ==  0.
            )
    {
        _maxX[ sensorID ] = _siPlanesLayerLayout->getSensitiveNpixelY( iPlane ) - 1;
        _maxY[ sensorID ] = _siPlanesLayerLayout->getSensitiveNpixelX( iPlane ) - 1;        

        _hitMinX[ sensorID ] =  _siPlanesLayerLayout->getSensitivePositionY( iPlane ) - 0.5*_siPlanesLayerLayout->getSensitiveSizeY ( iPlane ) ;
        _hitMaxX[ sensorID ] =  _siPlanesLayerLayout->getSensitivePositionY( iPlane ) + 0.5*_siPlanesLayerLayout->getSensitiveSizeY ( iPlane ) ;
        _hitMinY[ sensorID ] =  _siPlanesLayerLayout->getSensitivePositionX( iPlane ) - 0.5*_siPlanesLayerLayout->getSensitiveSizeX ( iPlane ) ;
        _hitMaxY[ sensorID ] =  _siPlanesLayerLayout->getSensitivePositionX( iPlane ) + 0.5*_siPlanesLayerLayout->getSensitiveSizeX ( iPlane ) ;
    } 
    else
        if(                
            _siPlanesRotations[iPlane][1] !=  0.                     &&
            _siPlanesRotations[iPlane][2] ==  0.                     &&
            _siPlanesRotations[iPlane][3] ==  0.                     &&
            _siPlanesRotations[iPlane][4] !=  0.
            )
    {
        _maxX[ sensorID ] = _siPlanesLayerLayout->getSensitiveNpixelX( iPlane ) - 1;
        _maxY[ sensorID ] = _siPlanesLayerLayout->getSensitiveNpixelY( iPlane ) - 1;        

        _hitMinX[ sensorID ] =  _siPlanesLayerLayout->getSensitivePositionX( iPlane ) - 0.5*_siPlanesLayerLayout->getSensitiveSizeX ( iPlane ) ;
        _hitMaxX[ sensorID ] =  _siPlanesLayerLayout->getSensitivePositionX( iPlane ) + 0.5*_siPlanesLayerLayout->getSensitiveSizeX ( iPlane ) ;
        _hitMinY[ sensorID ] =  _siPlanesLayerLayout->getSensitivePositionY( iPlane ) - 0.5*_siPlanesLayerLayout->getSensitiveSizeY ( iPlane ) ;
        _hitMaxY[ sensorID ] =  _siPlanesLayerLayout->getSensitivePositionY( iPlane ) + 0.5*_siPlanesLayerLayout->getSensitiveSizeY ( iPlane ) ;
    }   
        else
        {
            streamlog_out (WARNING5) << "unknown sensor rotation configuration ?! check your Gear file or ammend the code " << endl;
            _hitMinX[ sensorID ] =  -10000.;
            _hitMaxX[ sensorID ] =   10000.;
            _hitMinY[ sensorID ] =  -10000.;
            _hitMaxY[ sensorID ] =   10000.;            
        }
  
   }

   _outputCorrelatedHitCollectionVec = 0;

  _isInitialize = false;

}

void EUTelCorrelator::processRunHeader (LCRunHeader * rdr) {


  EUTelRunHeaderImpl * runHeader = new EUTelRunHeaderImpl( rdr ) ;

  // this is the right place also to check the geometry ID. This is a
  // unique number identifying each different geometry used at the
  // beam test. The same number should be saved in the run header and
  // in the xml file. If the numbers are different, instead of barely
  // quitting ask the user what to do.

  if ( runHeader->getGeoID() == 0 )
    streamlog_out ( WARNING0 ) <<  "The geometry ID in the run header is set to zero." << endl
                               <<  "This may mean that the GeoID parameter was not set" << endl;


  if ( runHeader->getGeoID() != _siPlanesParameters->getSiPlanesID() ) {
    streamlog_out ( WARNING5 ) <<  "Error during the geometry consistency check: " << endl
                             << "The run header says the GeoID is " << runHeader->getGeoID() << endl
                             << "The GEAR description says is     " << _siPlanesParameters->getSiPlanesID()
                             << endl;

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

#endif // EUTEL_INTERACTIVE
  }

  delete runHeader;

  // increment the run counter
  ++_iRun;
}


void EUTelCorrelator::processEvent (LCEvent * event) {

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

 
  if(_iEvt > _events) return;
      ++_iEvt;


  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event) ;

  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return;
  } else if ( evt->getEventType() == kUNKNOWN ) {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event."
                               << endl;
  }
  // if the Event that we are looking is the first we create files
  // with histograms.
  if ( !_isInitialize ) 
  {

    _hasClusterCollection = false;

    for( size_t i = 0; i < _clusterCollectionVec.size() ; i++ )
    {
      std::string _inputClusterCollectionName = _clusterCollectionVec[i];

       try 
       {
         // let's check if we have cluster collections
         event->getCollection( _inputClusterCollectionName );
 
         _hasClusterCollection = true;
         streamlog_out ( DEBUG5 ) << "found " << i <<  " name " <<   _inputClusterCollectionName.c_str() << endl;
 
       } catch ( lcio::Exception& e ) {
         streamlog_out ( WARNING ) << "NOT found " << i <<  " name " <<   _inputClusterCollectionName.c_str() << endl;
       }
    }

    try 
    {
      // let's check if we have hit collections

      event->getCollection( _inputHitCollectionName ) ;

      _hasHitCollection = true;

    } catch ( lcio::Exception& e ) {

      _hasHitCollection = false;
    }


    // check if we have at least one collection.
    if ( ! _hasClusterCollection && ! _hasHitCollection  &&  !_isInitialize) {

      // this is the case we didn't find any collection in this
      // event, so keep the first event flag to true in order to try
      // again with the next event. 
      _isInitialize = false;
      return;

    } else {
    
      bookHistos();
     

      _isInitialize = true;
    }

  }


  try {

    if ( _hasClusterCollection && !_hasHitCollection) {

      for( size_t i = 0; i < _clusterCollectionVec.size() ; i++ )
      {
         std::string _ExternalInputClusterCollectionName = _clusterCollectionVec[i];

         LCCollectionVec * externalInputClusterCollection   = static_cast<LCCollectionVec*>   (event->getCollection( _ExternalInputClusterCollectionName ));
         CellIDDecoder<TrackerPulseImpl>  pulseCellDecoder( externalInputClusterCollection );

      // we have an external detector where we consider a cluster each
      // time (external cluster) that is correlated with another
      // detector's clusters (internal cluster)

      for ( size_t iExt = 0 ; iExt < externalInputClusterCollection->size() ; ++iExt ) {

        TrackerPulseImpl * externalPulse = static_cast< TrackerPulseImpl * >   ( externalInputClusterCollection->getElementAt( iExt ) );

        EUTelVirtualCluster  * externalCluster;

        ClusterType type = static_cast<ClusterType>  (static_cast<int>((pulseCellDecoder(externalPulse)["type"])));

        // we check that the type of cluster is ok
        
        if ( type == kEUTelDFFClusterImpl ) 
        {
          externalCluster = new EUTelDFFClusterImpl( static_cast<TrackerDataImpl*>     ( externalPulse->getTrackerData()) );
        }
        else if ( type == kEUTelBrickedClusterImpl ) 
        {
          externalCluster = new EUTelBrickedClusterImpl( static_cast<TrackerDataImpl*> ( externalPulse->getTrackerData()) );
        }
        else if ( type == kEUTelFFClusterImpl ) 
        {
          externalCluster = new EUTelFFClusterImpl( static_cast<TrackerDataImpl*>      ( externalPulse->getTrackerData()) );
            
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
          if ( pixelType == kEUTelSimpleSparsePixel ) {
            externalCluster = new EUTelSparseClusterImpl< EUTelSimpleSparsePixel >
              ( static_cast<TrackerDataImpl *> ( externalPulse->getTrackerData()  ) );
          } else {
            streamlog_out ( ERROR4 ) << "Unknown pixel type.  Sorry for quitting." << endl;
            throw UnknownDataTypeException("Pixel type unknown");
          }
 
          if( externalCluster->getTotalCharge() < _clusterChargeMin ) 
          {
              delete externalCluster; 
              continue;
          }

       }
	else  continue;

        int externalSensorID = pulseCellDecoder( externalPulse ) [ "sensorID" ] ;
 
        float externalXCenter;
        float externalYCenter;

        // we catch the coordinates of the external seed

        externalCluster->getCenterOfGravity( externalXCenter, externalYCenter ) ;
        if( externalCluster->getTotalCharge() <= _clusterChargeMin ) 
        {
            delete externalCluster; 
            continue;
        }

        for( size_t i = 0; i < _clusterCollectionVec.size() ; i++ )
        {
          std::string _InternalInputClusterCollectionName = _clusterCollectionVec[i];

         LCCollectionVec * internalInputClusterCollection   = static_cast<LCCollectionVec*>   (event->getCollection( _InternalInputClusterCollectionName ));
         CellIDDecoder<TrackerPulseImpl>  pulseCellDecoder( internalInputClusterCollection );


        for ( size_t iInt = 0;  iInt <  internalInputClusterCollection->size() ; ++iInt ) 
        {

          TrackerPulseImpl * internalPulse = static_cast< TrackerPulseImpl * >  ( internalInputClusterCollection->getElementAt( iInt ) );

          EUTelVirtualCluster  * internalCluster;

          ClusterType type = static_cast<ClusterType>  (static_cast<int>((pulseCellDecoder(internalPulse)["type"])));

          // we check that the type of cluster is ok
          if ( type == kEUTelDFFClusterImpl ) {
            internalCluster = new EUTelDFFClusterImpl( static_cast<TrackerDataImpl*>
                                                       (internalPulse->getTrackerData()) );

          } else if ( type == kEUTelBrickedClusterImpl ) {
            internalCluster = new EUTelBrickedClusterImpl( static_cast<TrackerDataImpl*>
                                                       (internalPulse->getTrackerData()) );

          } else if ( type == kEUTelFFClusterImpl ) {
            internalCluster = new EUTelFFClusterImpl( static_cast<TrackerDataImpl*>
                                                      (internalPulse->getTrackerData()) );

          } else if ( type == kEUTelSparseClusterImpl ) {

            // ok the cluster is of sparse type, but we also need to know
            // the kind of pixel description used. This information is
            // stored in the corresponding original data collection.

            LCCollectionVec * sparseClusterCollectionVec = dynamic_cast < LCCollectionVec * > (evt->getCollection("original_zsdata"));
            TrackerDataImpl * oneCluster = dynamic_cast<TrackerDataImpl*> (sparseClusterCollectionVec->getElementAt( 0 ));
            CellIDDecoder<TrackerDataImpl > anotherDecoder(sparseClusterCollectionVec);
            SparsePixelType pixelType = static_cast<SparsePixelType> ( static_cast<int> ( anotherDecoder( oneCluster )["sparsePixelType"] ));

            // now we know the pixel type. So we can properly create a new
            // instance of the sparse cluster
            if ( pixelType == kEUTelSimpleSparsePixel ) {
              internalCluster = new EUTelSparseClusterImpl< EUTelSimpleSparsePixel >
                ( static_cast<TrackerDataImpl *> ( internalPulse->getTrackerData()  ) );
            } else {
              streamlog_out ( ERROR4 ) << "Unknown pixel type.  Sorry for quitting." << endl;
              throw UnknownDataTypeException("Pixel type unknown");
            }

            if( internalCluster->getTotalCharge() < _clusterChargeMin )
            {
                delete internalCluster;
                continue;
            }
 
          }
	   else  continue;

          if( internalCluster->getTotalCharge() < _clusterChargeMin )
          {
             delete internalCluster;
             continue;
          }
          

          int internalSensorID = pulseCellDecoder( internalPulse ) [ "sensorID" ] ;


          if ( ( internalSensorID != getFixedPlaneID() && externalSensorID == getFixedPlaneID() )
                  ||
                  (_sensorIDtoZOrderMap[internalSensorID] ==  _sensorIDtoZOrderMap[externalSensorID] + 1 )
                  ) 
          {

              float internalXCenter = 0. ;
              float internalYCenter = 0. ;

            // we catch the coordinates of the internal seed

            internalCluster->getCenterOfGravity( internalXCenter, internalYCenter ) ;

            streamlog_out ( DEBUG5 ) << "Filling histo " << externalSensorID << " " << internalSensorID << endl;


            // we input the coordinates in the correlation matrix, one
            // for each type of coordinate: X and Y

            // assume simpliest +1 and -1 only :
            int exPlaneGear = _sensorIDVecMap[externalSensorID];

                   std::vector<double> cluster_offset;
                   std::vector<double> cluCenter;
                   cluCenter.push_back(internalXCenter);
                   cluCenter.push_back(internalYCenter);
                   cluCenter.push_back(externalXCenter);
                   cluCenter.push_back(externalYCenter);
                   cluster_offset = guessSensorOffset(internalSensorID, externalSensorID, cluCenter);
                   
                   _clusterXCorrShiftMatrix[ externalSensorID ][ internalSensorID ]->fill( externalXCenter*_siPlanesPitchX[exPlaneGear]-_siPlanesLayerLayout->getSensitiveSizeX(exPlaneGear)/2., cluster_offset[0]  );
                   _clusterYCorrShiftMatrix[ externalSensorID ][ internalSensorID ]->fill( externalYCenter*_siPlanesPitchY[exPlaneGear]-_siPlanesLayerLayout->getSensitiveSizeY(exPlaneGear)/2., cluster_offset[1]  );

                   if( cluster_offset.size() >3 )
                   { 
                       _clusterXCorrelationMatrix[ externalSensorID ][ internalSensorID ]->fill( externalXCenter, cluster_offset[2]);
                       _clusterYCorrelationMatrix[ externalSensorID ][ internalSensorID ]->fill( externalYCenter, cluster_offset[3]);
                   }
                   else
                   { 
                       _clusterXCorrelationMatrix[ externalSensorID ][ internalSensorID ]->fill( externalXCenter, internalYCenter );
                       _clusterYCorrelationMatrix[ externalSensorID ][ internalSensorID ]->fill( externalYCenter, internalXCenter );
                   }
          } // endif

          delete internalCluster;

        } // internal loop
        } // internal loop of collections 

        delete externalCluster; 
      } // external loop
      } // external loop of collections

    } // endif hasCluster

    if ( _hasHitCollection ) {


      LCCollectionVec* inputHitCollection = static_cast<LCCollectionVec*>( event->getCollection(_inputHitCollectionName) );
      UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder ( EUTELESCOPE::HITENCODING );

      streamlog_out  ( MESSAGE2 ) << "inputHitCollection " << _inputHitCollectionName.c_str() << endl;


      for ( size_t iExt = 0 ; iExt < inputHitCollection->size(); ++iExt ) {
        std::vector<double> trackX;
        std::vector<double> trackY;
        std::vector<int  > iplane;

       trackX.clear();
       trackY.clear();
       iplane.clear();

        // this is the external hit

        TrackerHitImpl* externalHit = static_cast<TrackerHitImpl*>( inputHitCollection->getElementAt(iExt) );
        
	double* externalPosition = const_cast<double*>( externalHit->getPosition() );

        int externalSensorID = hitDecoder( externalHit )["sensorID"]; 

        double etrackPointLocal[]  = { externalPosition[0], externalPosition[1], externalPosition[2] };
        double etrackPointGlobal[] = { externalPosition[0], externalPosition[1], externalPosition[2] };

// do local2Master only if the externalHit is defined in the local frame of the module:  
// skip if already in global 
// steer by a processor flag ?
        geo::gGeometry().local2Master( externalSensorID, etrackPointLocal, etrackPointGlobal );
//

        trackX.push_back( etrackPointGlobal[0]);
        trackY.push_back( etrackPointGlobal[1]);

        iplane.push_back( externalSensorID);

        streamlog_out  ( MESSAGE2 ) << "eplane:"  << externalSensorID << " loc: "  << etrackPointLocal[0]  << " "<< etrackPointLocal[1]  << " "
                                                                      << " glo: "  << etrackPointGlobal[0] << " "<< etrackPointGlobal[1] << " " << endl;

        for ( size_t iInt = 0; iInt < inputHitCollection->size(); ++iInt ) 
        {


          TrackerHitImpl* internalHit = static_cast<TrackerHitImpl*>( inputHitCollection->getElementAt(iInt) );

          double* internalPosition = const_cast<double*>( internalHit->getPosition() );

          int internalSensorID = hitDecoder( internalHit )["sensorID"]; 

          bool ishot = hitContainsHotPixels(internalHit); 

          if( ishot ) continue;
   

          double itrackPointLocal[]  = { internalPosition[0], internalPosition[1], internalPosition[2] };
          double itrackPointGlobal[] = { internalPosition[0], internalPosition[1], internalPosition[2] };
// do local2Master only if the externalHit is defined in the local frame of the module:  
// skip if already in global 
// steer by a processor flag ?
          geo::gGeometry().local2Master( internalSensorID, itrackPointLocal, itrackPointGlobal );
 

          if ( 
                  ( internalSensorID != getFixedPlaneID() && externalSensorID == getFixedPlaneID() )
                   ||
                  _sensorIDtoZOrderMap[internalSensorID] ==  _sensorIDtoZOrderMap[externalSensorID] +1 
                  ) 
            {

            int iz = _sensorIDtoZOrderMap[internalSensorID]; 

            if(
               ((etrackPointGlobal[0]-itrackPointGlobal[0] ) < _residualsXMax[iz]) && (_residualsXMin[iz] < (etrackPointGlobal[0]-itrackPointGlobal[0] ))   
               &&
               ((etrackPointGlobal[1]-itrackPointGlobal[1] ) < _residualsYMax[iz]) && (_residualsYMin[iz] < (etrackPointGlobal[1]-itrackPointGlobal[1] ))
              )
               {

        trackX.push_back(itrackPointGlobal[0]);
        trackY.push_back(itrackPointGlobal[1]);
        iplane.push_back(internalSensorID);

        streamlog_out  ( MESSAGE2 ) << "iplane:"  << internalSensorID << " loc: "  << itrackPointLocal[0]  << " "<< itrackPointLocal[1]  << " "
                                                                      << " glo: "  << itrackPointGlobal[0] << " "<< itrackPointGlobal[1] << " " << endl;

               }
            }

        }

  vector<int>  iplane_unique = iplane;
  vector<int>::iterator p, p_end;
 
  p_end = unique( iplane_unique.begin(), iplane_unique.end());       // remove duplicates

if( static_cast< int >(iplane_unique.size()) > _minNumberOfCorrelatedHits && trackX.size() == trackY.size())
{
      for(size_t i=1;i< trackX.size();i++)
      {
            _hitXCorrelationMatrix[ iplane[0]        ] [ iplane[i]        ] -> fill ( trackX[0]          , trackX[i]           ) ;
            _hitYCorrelationMatrix[ iplane[0]        ] [ iplane[i]        ] -> fill ( trackY[0]          , trackY[i]           ) ;
            // assume all rotations have been done in the hitmaker processor:
            _hitXCorrShiftMatrix[ iplane[0]        ][ iplane[i]        ]->fill( trackX[0]          , trackX[0]          - trackX[i]          );
            _hitYCorrShiftMatrix[ iplane[0]        ][ iplane[i]        ]->fill( trackY[0]          , trackY[0]          - trackY[i]         );
      }
}else{
}
 
      }
    }
  } catch (DataNotAvailableException& e  ) {

    streamlog_out  ( MESSAGE2 ) <<  "No input collection found on event " << event->getEventNumber()
                                << " in run " << event->getRunNumber() << endl;
  }

#endif

}

void EUTelCorrelator::end() {


    if( _hasClusterCollection && !_hasHitCollection)
    {
        streamlog_out( MESSAGE5 ) << "The input CollectionVec contains ClusterCollection, calculating offset values " << endl;
        
        for ( int iin = 0 ; iin < _siPlanesLayerLayout->getNLayers(); iin++ ) 
        {           
            int inPlane = _siPlanesLayerLayout->getID( iin );
            for ( int iex = 0 ; iex < _siPlanesLayerLayout->getNLayers(); iex++ ) 
            {
                int exPlane = _siPlanesLayerLayout->getID( iex );

                if(
                    !( inPlane != getFixedPlaneID() && exPlane == getFixedPlaneID() )
                  )continue;

                if( _clusterXCorrShiftMatrix[ exPlane ][ inPlane ] == 0 ) continue;
                if( _clusterXCorrShiftMatrix[ exPlane ][ inPlane ]->yAxis().bins() <= 0 ) continue;

                float _heighestBinX = 0.;
                for(int ibin = 0; ibin < _clusterXCorrShiftMatrix[ exPlane ][ inPlane ]->yAxis().bins(); ibin++)
                {
                    double xbin =  
                        _clusterXCorrShiftProjection[ inPlane ]->axis().binLowerEdge(ibin)
                        +
                        _clusterXCorrShiftProjection[ inPlane ]->axis().binWidth(ibin)/2.
                        ;
                    double _binValue = _clusterXCorrShiftMatrix[ exPlane ][ inPlane ]->binEntriesY( ibin );

                    _clusterXCorrShiftProjection[ inPlane ]->fill( xbin, _binValue );
                    if( _binValue > _heighestBinX )
                    {
                        _heighestBinX = _binValue;
                    }
                }
                
               
                float _heighestBinY = 0.;
                for(int ibin = 0; ibin < _clusterYCorrShiftMatrix[ exPlane ][ inPlane ]->yAxis().bins(); ibin++)
                {
                    double xbin =  
                        _clusterYCorrShiftProjection[ inPlane ]->axis().binLowerEdge(ibin)
                        +
                        _clusterYCorrShiftProjection[ inPlane ]->axis().binWidth(ibin)/2.
                        ;
                    double _binValue = _clusterYCorrShiftMatrix[ exPlane ][ inPlane ]->binEntriesY( ibin );
                    _clusterYCorrShiftProjection[ inPlane ]->fill( xbin, _binValue );
                    if( _binValue > _heighestBinY )
                    {
                        _heighestBinY = _binValue;
                    }
                }

                // get the highert bin and its neighbours
                // 
                double _correlationBandBinsX     = 0.;
                double _correlationBandCenterX   = 0.;

                for(int ibin = 0; ibin < _clusterXCorrShiftProjection[ inPlane ]->axis().bins(); ibin++)
                {
                    double ybin =  _clusterXCorrShiftProjection[ inPlane ]->binHeight(ibin); 
                   
                    if( ybin < _heighestBinX*0.9 ) continue;
                    double xbin =  
                        _clusterXCorrShiftProjection[ inPlane ]->axis().binLowerEdge(ibin)
                        +
                        _clusterXCorrShiftProjection[ inPlane ]->axis().binWidth(ibin)/2.
                        ;
                    

                    _correlationBandBinsX   += ybin;
                    _correlationBandCenterX += xbin*ybin;
                }

                
                double _correlationBandBinsY     = 0.;
                double _correlationBandCenterY   = 0.;

                for(int ibin = 0; ibin < _clusterYCorrShiftMatrix[ exPlane ][ inPlane ]->yAxis().bins(); ibin++)
                {
                    double ybin =  _clusterYCorrShiftProjection[ inPlane ]->binHeight(ibin); 
                    
                    if( ybin < _heighestBinY*0.9 ) continue;
                    double xbin =  
                        _clusterYCorrShiftProjection[ inPlane ]->axis().binLowerEdge(ibin)
                        +
                        _clusterYCorrShiftProjection[ inPlane ]->axis().binWidth(ibin)/2.
                        ;                    
                  
                   _correlationBandBinsY   += ybin;
                   _correlationBandCenterY += ybin*xbin;
                }

                int inPlaneGear = _sensorIDVecMap[inPlane];

                if( _correlationBandBinsX != 0. ) 
                    _siPlanesOffsetX[ inPlaneGear ] = -1*_correlationBandCenterX/_correlationBandBinsX;
                else
                    _siPlanesOffsetX[ inPlaneGear ] = 0.;    
 
                if( _correlationBandBinsY != 0. ) 
                    _siPlanesOffsetY[ inPlaneGear ] = -1*_correlationBandCenterY/_correlationBandBinsY;
                else
                    _siPlanesOffsetY[ inPlaneGear ] = 0.;    
               
                streamlog_out( MESSAGE5 ) << "Offsets in plane " << std::setw(3) << inPlane << ": dX=" << std::setw(9) << std::setprecision(3) << _siPlanesOffsetX[inPlaneGear] << " um,    dY=" << std::setw(9) << std::setprecision(3) << _siPlanesOffsetY[inPlaneGear] << " um" << endl;

           }
        }
    }
 
    if( _hasHitCollection)
    {
        streamlog_out( MESSAGE5 ) << "The input CollectionVec contains HitCollection, calculating offest values " << endl;
 
        for ( int iin = 0 ; iin < _siPlanesLayerLayout->getNLayers(); iin++ ) 
        {           
            int inPlane = _siPlanesLayerLayout->getID( iin );
            for ( int iex = 0 ; iex < _siPlanesLayerLayout->getNLayers(); iex++ ) 
            {
                int exPlane = _siPlanesLayerLayout->getID( iex );
                if( _hitXCorrShiftMatrix[ exPlane ][ inPlane ] == 0 ) continue;
                if( _hitXCorrShiftMatrix[ exPlane ][ inPlane ]->yAxis().bins() <= 0 ) continue;

                if(
                    !( inPlane != getFixedPlaneID() && exPlane == getFixedPlaneID() )
                  )continue;

                float _heighestBinX = 0.;
                for(int ibin = 0; ibin < _hitXCorrShiftMatrix[ exPlane ][ inPlane ]->yAxis().bins(); ibin++)
                {
                    double xbin =  
                        _hitXCorrShiftProjection[ inPlane ]->axis().binLowerEdge(ibin)
                        +
                        _hitXCorrShiftProjection[ inPlane ]->axis().binWidth(ibin)/2.
                        ;
                    double _binValue = _hitXCorrShiftMatrix[ exPlane ][ inPlane ]->binEntriesY( ibin );
                    _hitXCorrShiftProjection[ inPlane ]->fill( xbin, _binValue );
                    if( _binValue>0)
                    if( _binValue > _heighestBinX )
                    {
                        _heighestBinX = _binValue;
                   }
                }
                
               
                float _heighestBinY = 0.;
                for(int ibin = 0; ibin < _hitYCorrShiftMatrix[ exPlane ][ inPlane ]->yAxis().bins(); ibin++)
                {
                    double xbin =  
                        _hitYCorrShiftProjection[ inPlane ]->axis().binLowerEdge(ibin)
                        +
                        _hitYCorrShiftProjection[ inPlane ]->axis().binWidth(ibin)/2.
                        ;
                    double _binValue = _hitYCorrShiftMatrix[ exPlane ][ inPlane ]->binEntriesY( ibin );
                    _hitYCorrShiftProjection[ inPlane ]->fill( xbin, _binValue );
                    if( _binValue>0)
                    if( _binValue > _heighestBinY )
                    {
                        _heighestBinY = _binValue;
                   }
                }

                
                // get the highert bin and its neighbours
                // 
                double _correlationBandBinsX     = 0.;
                double _correlationBandCenterX   = 0.;

                for(int ibin = 0; ibin < _hitXCorrShiftProjection[ inPlane ]->axis().bins(); ibin++)
                {
                    double ybin =  _hitXCorrShiftProjection[ inPlane ]->binHeight(ibin); 

                    if( ybin < _heighestBinX*0.9 ) continue;
                    double xbin =  
                        _hitXCorrShiftProjection[ inPlane ]->axis().binLowerEdge(ibin)
                        +
                        _hitXCorrShiftProjection[ inPlane ]->axis().binWidth(ibin)/2.
                        ;
                    

                    _correlationBandBinsX   += ybin;
                    _correlationBandCenterX += xbin*ybin;
                }

                
                double _correlationBandBinsY     = 0.;
                double _correlationBandCenterY   = 0.;

                for(int ibin = 0; ibin < _hitYCorrShiftMatrix[ exPlane ][ inPlane ]->yAxis().bins(); ibin++)
                {
                    double ybin =  _hitYCorrShiftProjection[ inPlane ]->binHeight(ibin); 
                    
                    if( ybin < _heighestBinY*0.9  ) continue;
                    double xbin =  
                        _hitYCorrShiftProjection[ inPlane ]->axis().binLowerEdge(ibin)
                        +
                        _hitYCorrShiftProjection[ inPlane ]->axis().binWidth(ibin)/2.
                        ;                    
                  
                   _correlationBandBinsY   += ybin;
                   _correlationBandCenterY += ybin*xbin;
               }

// alter sign:
               _correlationBandCenterX = -_correlationBandCenterX;                
               _correlationBandCenterY = -_correlationBandCenterY;                
               
                streamlog_out( MESSAGE5 ) << "Hit Offset values: " ; 
                streamlog_out ( MESSAGE5 ) << " plane : " << inPlane << " " ;
                streamlog_out ( MESSAGE5 ) << " X offset : "<<  (_correlationBandBinsX == 0. ? 0.: _correlationBandCenterX/_correlationBandBinsX) ; 
                streamlog_out ( MESSAGE5 ) << " Y offset : "<<  (_correlationBandBinsY == 0. ? 0.: _correlationBandCenterY/_correlationBandBinsY) ; 
                streamlog_out( MESSAGE5 ) << endl;


            }
        }
    }

    if( _dumpOffset && !_hasHitCollection )
    {

        // reopen the LCIO file this time in append mode
        LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();

        try 
        {
            lcWriter->open( _offsetDBFile, LCIO::WRITE_NEW );
        }
        catch ( IOException& e ) 
        {
            streamlog_out ( ERROR4 ) << e.what() << endl
                << "Sorry for quitting. " << endl;
            exit(-1);
        }
        
        LCEventImpl *event = new LCEventImpl;
        event->setRunNumber( 0 );
        event->setEventNumber( 0 );
        event->setDetectorName("Offset DB");

        LCTime *now = new LCTime;
        event->setTimeStamp( now->timeStamp() );
        delete now;

        LCCollectionVec * constantsCollection = new LCCollectionVec( LCIO::LCGENERICOBJECT );


        for ( int iin = 0 ; iin < _siPlanesLayerLayout->getNLayers(); iin++ ) 
        {           
            int _sensorID = _siPlanesLayerLayout->getID( iin );
            EUTelAlignmentConstant * constant = new EUTelAlignmentConstant;

            constant->setXOffset( _siPlanesOffsetX[iin]  );
            constant->setYOffset( _siPlanesOffsetY[iin]  ) ;
 
            constant->setSensorID( _sensorID );
            constantsCollection->push_back( constant );
            
        }

        event->addCollection( constantsCollection, "preAlignment" ); 
        lcWriter->writeEvent( event );        
        delete event;
    
        lcWriter->close();        
        delete lcWriter;
    }
  
    streamlog_out ( MESSAGE4 )  << "Successfully finished" << endl;
    delete [] _siPlaneZPosition;
}

void EUTelCorrelator::bookHistos() {

  if ( !_hasClusterCollection && !_hasHitCollection ) return ;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  try {

    streamlog_out ( DEBUG5 ) <<  "Booking histograms" << endl;

        auto_ptr<EUTelHistogramManager> histoMgr( new EUTelHistogramManager( _histoInfoFileName ));
        EUTelHistogramInfo    * histoInfo;
        bool                    isHistoManagerAvailable;

        try {
            isHistoManagerAvailable = histoMgr->init( );
        } catch ( ios::failure& e ) {
            streamlog_out( ERROR5 ) << "I/O problem with " << _histoInfoFileName << "\n"
                    << "Continuing without histogram manager using default settings"    << endl;
            isHistoManagerAvailable = false;
        } catch ( ParseException& e ) {
            streamlog_out( ERROR5 ) << e.what( ) << "\n"
                    << "Continuing without histogram manager using default settings" << endl;
            isHistoManagerAvailable = false;
        }

// declare.initialize:

        int    xBin  =  10 ;    
        double xMin  = -10.;
        double xMax  =  10.;

        int    yBin  =  10 ;    
        double yMin  = -10.;  
        double yMax  =  10.;

        int colNBin   =  64 ;
        double colMin =   0.;
        double colMax =  64.;

        int rowNBin   =  64 ;
        double rowMin =   0.;
        double rowMax =  64.;


    // create all the directories first
    vector< string > dirNames;

    if ( _hasClusterCollection && !_hasHitCollection) {
      dirNames.push_back ("ClusterX");
      dirNames.push_back ("ClusterY");
      dirNames.push_back ("ClusterXShift");
      dirNames.push_back ("ClusterYShift");
    }

    if ( _hasHitCollection ) {
      dirNames.push_back ("HitX");
      dirNames.push_back ("HitY");
      dirNames.push_back ("HitXShift");
      dirNames.push_back ("HitYShift");
    }

    for ( size_t iPos = 0 ; iPos < dirNames.size() ; iPos++ ) {

      AIDAProcessor::tree(this)->mkdir( dirNames[iPos].c_str() ) ;

    }


    string tempHistoName  = "";
    string tempHistoTitle = "";


    for ( size_t r = 0 ; r < _sensorIDVec.size(); ++r ) 
    {

      int row = _sensorIDVec.at( r );
      
      map< unsigned int , AIDA::IHistogram2D * > innerMapXCluster;
      map< unsigned int , AIDA::IHistogram2D * > innerMapYCluster;

      map< unsigned int , AIDA::IHistogram2D * > innerMapXCluShift;
      map< unsigned int , AIDA::IHistogram2D * > innerMapYCluShift;
      map< unsigned int , AIDA::IHistogram1D * > innerMapXCluShiftProjection;
      map< unsigned int , AIDA::IHistogram1D * > innerMapYCluShiftProjection;

      map< unsigned int , AIDA::IHistogram2D * > innerMapXHit;
      map< unsigned int , AIDA::IHistogram2D * > innerMapYHit;

      map< unsigned int , AIDA::IHistogram2D * > innerMapXHitShift;
      map< unsigned int , AIDA::IHistogram2D * > innerMapYHitShift;
      map< unsigned int , AIDA::IHistogram1D * > innerMapXHitShiftProjection;
      map< unsigned int , AIDA::IHistogram1D * > innerMapYHitShiftProjection;



      for ( size_t c = 0 ; c < _sensorIDVec.size(); ++c ) {
 
        int col = _sensorIDVec.at( c );

         if ( 
// optional       (_sensorIDtoZOrderMap[internalSensorID] != 0 && _sensorIDtoZOrderMap[externalSensorID] == 0)
                  ( col != getFixedPlaneID() && row == getFixedPlaneID() )
                  ||
                  (_sensorIDtoZOrderMap[ col ] ==  _sensorIDtoZOrderMap[ row ] + 1 )
                  ) 
          {


          //we create histograms for X and Y Cluster correlation
          if ( _hasClusterCollection && !_hasHitCollection) {

            double safetyFactor = 1.0; // 2 should be enough because it
            // means that the sensor is wrong
            // by all its size.

            /////////////////////////////////////////////////
            // book X
            tempHistoName = "ClusterX/" + _clusterXCorrelationHistoName + "_d" + to_string( row ) + "_d" + to_string( col );
            streamlog_out( DEBUG5 ) << "Booking histo " << tempHistoName << endl;

            histoInfo = histoMgr->getHistogramInfo(_clusterXCorrelationHistoName);
            xBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin :  100  ;    
            xMin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : -1000.;
            xMax  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax :  1000.;
            yBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin :  100  ;    
            yMin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : -1000.;
            yMax  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax :  1000.;

            AIDA::IHistogram2D * histo2D =
              AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(),
                                                                        xBin, xMin, xMax, yBin, yMin, yMax );

            tempHistoTitle =  "ClusterX/" +  _clusterXCorrelationHistoName + "_d" + to_string( row ) + "_d" + to_string( col );
            histo2D->setTitle( tempHistoTitle.c_str() );
            innerMapXCluster[ col  ] =  histo2D ;

            /////////////////////////////////////////////////
            // book Y
            tempHistoName =  "ClusterY/" +  _clusterYCorrelationHistoName + "_d" + to_string( row ) + "_d" + to_string( col );
            streamlog_out( DEBUG5 ) << "Booking histo " << tempHistoName << endl;

            histoInfo = histoMgr->getHistogramInfo(_clusterYCorrelationHistoName);
            xBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 100  ;    
            xMin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : -1000.;
            xMax  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax :  1000.;
            yBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : 100  ;    
            yMin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : -1000.;
            yMax  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax :  1000.;

            histo2D =
              AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(),
                                                                        xBin, xMin, xMax, yBin, yMin, yMax );
            tempHistoTitle =  "ClusterY/" +  _clusterYCorrelationHistoName + "_d" + to_string( row ) + "_d" + to_string( col );
            histo2D->setTitle( tempHistoTitle.c_str()) ;

            innerMapYCluster[ col  ] =  histo2D ;
            
            /////////////////////////////////////////////////
            // book special histos to calculate sensors initial offsets in X and Y
            // book X
            tempHistoName =  "ClusterXShift/" +  _clusterXCorrShiftHistoName + "_d" + to_string( row ) + "_d" + to_string( col );
            streamlog_out( DEBUG5 ) << "Booking histo " << tempHistoName << endl;

            histoInfo = histoMgr->getHistogramInfo(_clusterXCorrShiftHistoName);
            xBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 100   ;    
            xMin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : -1000.;
            xMax  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax :  1000.;
            yBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : 100   ;    
            yMin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : -1000.;
            yMax  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax :  1000.;

            histo2D =
              AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(),
                                                                        xBin, xMin, xMax, yBin, yMin, yMax );
            tempHistoTitle =  "ClusterXShift/" +  _clusterXCorrShiftHistoName + "_d" + to_string( row ) + "_d" + to_string( col );
            histo2D->setTitle( tempHistoTitle.c_str()) ;

            innerMapXCluShift[ col  ] =  histo2D ;

               /////////////////////////////////////////////////
            // book Y
            tempHistoName =  "ClusterYShift/" +  _clusterYCorrShiftHistoName + "_d" + to_string( row ) + "_d" + to_string( col );
            streamlog_out( DEBUG5 ) << "Booking histo " << tempHistoName << endl;

            histoInfo = histoMgr->getHistogramInfo(_clusterYCorrShiftHistoName);
            xBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 100   ;    
            xMin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : -1000.;
            xMax  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax :  1000.;
            yBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : 100   ;    
            yMin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : -1000.;
            yMax  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax :  1000.;

            histo2D =
              AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(),
                                                                        xBin, xMin, xMax, yBin, yMin, yMax );
            tempHistoTitle =  "ClusterYShift/" +  _clusterYCorrShiftHistoName + "_d" + to_string( row ) + "_d" + to_string( col );
            histo2D->setTitle( tempHistoTitle.c_str()) ;

            innerMapYCluShift[ col  ] =  histo2D ;
          }


          // the idea of using ICloud2D instead of H2D is interesting,
          // but because of a problem when the clouds is converted, I
          // prefer to use an histogram with a standard binning.
          //
          // the boundaries of the histos can be read from the GEAR
          // description and for safety multiplied by a safety factor
          // to take into account possible misalignment.

          if ( _hasHitCollection ) {

            double safetyFactor = 1.0; // 2 should be enough because it
            // means that the sensor is wrong
            // by all its size.

          
            tempHistoName  =  "HitX/" +  _hitXCorrelationHistoName + "_d" + to_string( row ) + "_d" + to_string( col );
            streamlog_out( DEBUG5 ) << "Booking histo " << tempHistoName << endl;
            tempHistoTitle =  "HitX/" +  _hitXCorrelationHistoName + "_d" + to_string( row ) + "_d" +  to_string( col );

            histoInfo = histoMgr->getHistogramInfo(_hitXCorrelationHistoName);
            colNBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 100   ;    
            colMin   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : -10.;
            colMax   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax :  10.;
            rowNBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : 100   ;    
            rowMin   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : -10.;
            rowMax   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax :  10.;


            AIDA::IHistogram2D * histo2D =
              AIDAProcessor::histogramFactory( this )->createHistogram2D( tempHistoName.c_str(), 
                      rowNBin, rowMin, rowMax, colNBin, colMin, colMax );
            
            histo2D->setTitle( tempHistoTitle.c_str() );

            innerMapXHit[ col  ] =  histo2D ;


            // now the hit on the Y direction
            colNBin = _maxY[col];
            rowNBin = _maxY[row];
    
            rowMin = safetyFactor * ( _hitMinY[row]);
            rowMax = safetyFactor * ( _hitMaxY[row]);
            colMin = safetyFactor * ( _hitMinY[col]);
            colMax = safetyFactor * ( _hitMaxY[col]);

            tempHistoName =  "HitY/" + _hitYCorrelationHistoName + "_d" + to_string( row ) + "_d" + to_string( col );
            streamlog_out( DEBUG5 ) << "Booking cloud " << tempHistoName << endl;
            tempHistoTitle = "HitY/" + _hitYCorrelationHistoName + "_d" + to_string( row ) + "_d" + to_string( col ) ;
 
            histoInfo = histoMgr->getHistogramInfo(_hitYCorrelationHistoName);
            colNBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 100   ;    
            colMin   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : -10.;
            colMax   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax :  10.;
            rowNBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : 100   ;    
            rowMin   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : -10.;
            rowMax   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax :  10.;

            histo2D =
              AIDAProcessor::histogramFactory( this )->createHistogram2D( tempHistoName.c_str(), 
                       rowNBin, rowMin, rowMax, colNBin, colMin, colMax );

            histo2D->setTitle( tempHistoTitle.c_str() );

            innerMapYHit[ col ] =  histo2D ;

           
            // book special histos to calculate sensors initial offsets in X and Y
            // book X
            tempHistoName =  "HitXShift/" +  _hitXCorrShiftHistoName + "_d" + to_string( row ) + "_d" + to_string( col );

            streamlog_out( DEBUG5 ) << "Booking histo " << tempHistoName << endl;

            histoInfo = histoMgr->getHistogramInfo(_hitXCorrShiftHistoName);
            colNBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 100   ;    
            colMin   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : -10.;
            colMax   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax :  10.;
            rowNBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : 100   ;    
            rowMin   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : -10.;
            rowMax   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax :  10.;

            histo2D = AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(), rowNBin, rowMin, rowMax, colNBin, colMin, colMax );

            tempHistoTitle =  "HitXShift/" +  _hitXCorrShiftHistoName + "_d" + to_string( row ) + "_d" + to_string( col );
            histo2D->setTitle( tempHistoTitle.c_str()) ;
            innerMapXHitShift[ col  ] =  histo2D ;


            // book Y
            tempHistoName =  "HitYShift/" +  _hitYCorrShiftHistoName + "_d" + to_string( row ) + "_d" + to_string( col );

            streamlog_out( DEBUG5 ) << "Booking histo " << tempHistoName << endl;

            histoInfo = histoMgr->getHistogramInfo(_hitYCorrShiftHistoName);
            colNBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 100   ;    
            colMin   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : -10.;
            colMax   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax :  10.;
            rowNBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : 100   ;    
            rowMin   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : -10.;
            rowMax   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax :  10.;

            histo2D = AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(), rowNBin, rowMin, rowMax, colNBin, colMin, colMax );
           
            tempHistoTitle =  "HitYShift/" +  _hitYCorrShiftHistoName + "_d" + to_string( row ) + "_d" + to_string( col );
            histo2D->setTitle( tempHistoTitle.c_str()) ;
            innerMapYHitShift[ col  ] =  histo2D ;
          }
 
        } else {

          if ( _hasClusterCollection && !_hasHitCollection) {
            innerMapXCluster[ col ]  = NULL ;
            innerMapYCluster[ col ]  = NULL ;
 
            innerMapXCluShift[ col ] = NULL ;
            innerMapYCluShift[ col ] = NULL ;            
            innerMapXCluShiftProjection[ col ] = NULL ;
            innerMapYCluShiftProjection[ col ] = NULL ;            
          }

          if ( _hasHitCollection ) {
            innerMapXHit[ col ] = NULL ;
            innerMapYHit[ col ] = NULL ;
 
            innerMapXHitShift[ col ] = NULL ;
            innerMapYHitShift[ col ] = NULL ;            
            innerMapXHitShiftProjection[ col ] = NULL ;
            innerMapYHitShiftProjection[ col ] = NULL ;            
          }

        }
 
      }



      if ( _hasClusterCollection && !_hasHitCollection) 
      {
        _clusterXCorrelationMatrix[ row ] = innerMapXCluster  ;
        _clusterYCorrelationMatrix[ row ] = innerMapYCluster  ;        
 
        _clusterXCorrShiftMatrix[ row ]  = innerMapXCluShift  ;
        _clusterYCorrShiftMatrix[ row ]  = innerMapYCluShift  ;        
 
            // book special histos to calculate sensors initial offsets in X and Y (Projection histograms)
            // book X
            tempHistoName =  "ClusterXShift/" +  _clusterXCorrShiftProjectionHistoName + "_d" + to_string( row ) ;
            streamlog_out( DEBUG5 ) << "Booking histo " << tempHistoName << endl;

            double safetyFactor = 1.0; // 2 should be enough because it
            // means that the sensor is wrong
            // by all its size.
 
            histoInfo = histoMgr->getHistogramInfo(_clusterXCorrShiftProjectionHistoName);
            xBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 100   ;    
            xMin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : -1000.;
            xMax  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax :  1000.;
            yBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : 100   ;    
            yMin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : -1000.;
            yMax  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax :  1000.;

            AIDA::IHistogram1D *
                 histo1D = AIDAProcessor::histogramFactory(this)->createHistogram1D( tempHistoName.c_str(), xBin, xMin, xMax );
            tempHistoTitle =  "ClusterXShift/" +  _clusterXCorrShiftProjectionHistoName + "_d" + to_string( row );
            histo1D->setTitle( tempHistoTitle.c_str()) ;

            _clusterXCorrShiftProjection[ row ]  = histo1D  ;        


            // book Y
            tempHistoName =  "ClusterYShift/" +  _clusterYCorrShiftProjectionHistoName + "_d" + to_string( row ) ;
            streamlog_out( DEBUG5 ) << "Booking histo " << tempHistoName << endl;

            histoInfo = histoMgr->getHistogramInfo(_clusterYCorrShiftProjectionHistoName);
            xBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 100   ;    
            xMin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : -1000.;
            xMax  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax :  1000.;
            yBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : 100   ;    
            yMin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : -1000.;
            yMax  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax :  1000.;


            histo1D = AIDAProcessor::histogramFactory(this)->createHistogram1D( tempHistoName.c_str(), yBin, yMin, yMax );
            tempHistoTitle =  "ClusterYShift/" +  _clusterYCorrShiftProjectionHistoName + "_d" + to_string( row ) ;
            histo1D->setTitle( tempHistoTitle.c_str()) ;

            _clusterYCorrShiftProjection[ row ]  = histo1D  ;        

        
      }

      if ( _hasHitCollection ) 
      {
         _hitXCorrelationMatrix[ row ] = innerMapXHit;
         _hitYCorrelationMatrix[ row ] = innerMapYHit;

         _hitXCorrShiftMatrix[ row ]   = innerMapXHitShift  ;
         _hitYCorrShiftMatrix[ row ]   = innerMapYHitShift  ;        
 
            // book special histos to calculate sensors initial offsets in X and Y (Projection histograms)
            // book X
            tempHistoName =  "HitXShift/" +  _hitXCorrShiftProjectionHistoName + "_d" + to_string( row ) ;

            streamlog_out( DEBUG5 ) << "Booking histo " << tempHistoName << endl;

 
            double safetyFactor = 1.0; // 2 should be enough because it
            // means that the sensor is wrong
            // by all its size.

            histoInfo = histoMgr->getHistogramInfo(_hitXCorrShiftProjectionHistoName);
            xBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin :  100 ;    
            xMin   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : -10.;
            xMax   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax :  10.;

            AIDA::IHistogram1D * histo1D = 0;
            histo1D = AIDAProcessor::histogramFactory(this)->createHistogram1D( tempHistoName.c_str(), xBin, xMin, xMax );
            tempHistoTitle =  "HitXShift/" +  _hitXCorrShiftProjectionHistoName + "_d" + to_string( row );
            histo1D->setTitle( tempHistoTitle.c_str()) ;

            _hitXCorrShiftProjection[ row ]  = histo1D  ;        


            // book Y
            tempHistoName =  "HitYShift/" +  _hitYCorrShiftProjectionHistoName + "_d" + to_string( row ) ;

            streamlog_out( DEBUG5 ) << "Booking histo " << tempHistoName << endl;

            histoInfo = histoMgr->getHistogramInfo(_hitXCorrShiftProjectionHistoName);
            xBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 100  ;    
            xMin   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : -10.;
            xMax   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax :  10.;

            histo1D = AIDAProcessor::histogramFactory(this)->createHistogram1D( tempHistoName.c_str(), xBin, xMin, xMax );
            tempHistoTitle =  "HitYShift/" +  _hitYCorrShiftProjectionHistoName + "_d" + to_string( row ) ;
            histo1D->setTitle( tempHistoTitle.c_str()) ;

            _hitYCorrShiftProjection[ row ]  = histo1D  ;        
     }
      
    }

  } catch (lcio::Exception& e ) {

    streamlog_out ( ERROR1 ) << "No AIDAProcessor initialized. Sorry for quitting..." << endl;
    exit( -1 );

  }
#endif
}


std::vector<double> EUTelCorrelator::guessSensorOffset(int internalSensorID, int externalSensorID, std::vector<double> cluCenter)

{
    double internalXCenter = cluCenter.at(0);
    double internalYCenter = cluCenter.at(1);
    double externalXCenter = cluCenter.at(2);
    double externalYCenter = cluCenter.at(3);

    int inPlaneGear = _sensorIDVecMap[internalSensorID];
    int exPlaneGear = _sensorIDVecMap[externalSensorID];
    
    double xDet_in =  internalXCenter*_siPlanesPitchX[inPlaneGear];
    double yDet_in =  internalYCenter*_siPlanesPitchY[inPlaneGear] ;
    double xDet_ex =  externalXCenter*_siPlanesPitchX[exPlaneGear] ;
    double yDet_ex =  externalYCenter*_siPlanesPitchY[exPlaneGear] ;

    double xCoo_in = internalXCenter;                      
    double yCoo_in = internalYCenter; 
   
 
    // get rotated sensors coordinates (in mm or um)
    double xPos_in =  xDet_in*( _siPlanesRotations[inPlaneGear][1] ) + yDet_in*( _siPlanesRotations[inPlaneGear][2] );
    double yPos_in =  xDet_in*( _siPlanesRotations[inPlaneGear][3] ) + yDet_in*( _siPlanesRotations[inPlaneGear][4] );
    double xPos_ex =  xDet_ex*( _siPlanesRotations[exPlaneGear][1] ) + yDet_ex*( _siPlanesRotations[exPlaneGear][2] );
    double yPos_ex =  xDet_ex*( _siPlanesRotations[exPlaneGear][3] ) + yDet_ex*( _siPlanesRotations[exPlaneGear][4] );

    double xCooPos_in =  xCoo_in*( _siPlanesRotations[inPlaneGear][1] ) + yCoo_in*( _siPlanesRotations[inPlaneGear][2] );
    double yCooPos_in =  xCoo_in*( _siPlanesRotations[inPlaneGear][3] ) + yCoo_in*( _siPlanesRotations[inPlaneGear][4] );

    // get rotated sensor coordinates (only in pixel number: col num)
   
    double sign = 0.;

      sign = 0.;
      if      ( _siPlanesRotations[inPlaneGear][1] < -0.7 )       sign = -1 ;
      else if ( _siPlanesRotations[inPlaneGear][1] > 0.7 )       sign =  1 ;
      else 
      {
        if       ( _siPlanesRotations[inPlaneGear][2] < -0.7 )    sign = -1 ;
        else if  ( _siPlanesRotations[inPlaneGear][2] > 0.7 )    sign =  1 ;
      }
      xPos_in +=  _siPlanesLayerLayout->getSensitivePositionX( inPlaneGear ) - sign*_siPlanesLayerLayout->getSensitiveSizeX ( inPlaneGear )/2. ;

      xCooPos_in = xCooPos_in*sign;

      sign = 0.;
      if      ( _siPlanesRotations[inPlaneGear][3] < -0.7 )       sign = -1 ;
      else if ( _siPlanesRotations[inPlaneGear][3] > 0.7 )       sign =  1 ;
      else 
      {
        if       ( _siPlanesRotations[inPlaneGear][4] < -0.7 )    sign = -1 ;
        else if  ( _siPlanesRotations[inPlaneGear][4] > 0.7 )    sign =  1 ;
      }
      yPos_in +=  _siPlanesLayerLayout->getSensitivePositionY( inPlaneGear ) - sign*_siPlanesLayerLayout->getSensitiveSizeY ( inPlaneGear )/2. ;
      yCooPos_in = yCooPos_in*sign;


      sign = 0.;
      if      ( _siPlanesRotations[exPlaneGear][1] < -0.7 )       sign = -1 ;
      else if ( _siPlanesRotations[exPlaneGear][1] > 0.7 )       sign =  1 ;
      else 
      {
        if       ( _siPlanesRotations[exPlaneGear][2] < -0.7 )    sign = -1 ;
        else if  ( _siPlanesRotations[exPlaneGear][2] > 0.7 )    sign =  1 ;
      }
      xPos_ex +=  _siPlanesLayerLayout->getSensitivePositionX( exPlaneGear ) - sign*_siPlanesLayerLayout->getSensitiveSizeX ( exPlaneGear )/2. ;

      sign = 0.;
      if      ( _siPlanesRotations[exPlaneGear][3] < -0.7 )       sign = -1 ;
      else if ( _siPlanesRotations[exPlaneGear][3] > 0.7 )       sign =  1 ;
      else 
      {
        if       ( _siPlanesRotations[exPlaneGear][4] < -0.7 )    sign = -1 ;
        else if  ( _siPlanesRotations[exPlaneGear][4] > 0.7 )    sign =  1 ;
      }
      yPos_ex +=  _siPlanesLayerLayout->getSensitivePositionY( exPlaneGear ) - sign*_siPlanesLayerLayout->getSensitiveSizeY ( exPlaneGear )/2. ;

                     
      std::vector<double> cluster_offset;

      cluster_offset.push_back( -( xPos_in - xPos_ex ) );
      cluster_offset.push_back( -( yPos_in - yPos_ex ) );
// 
// add also internal sensor X and Y coord
// 
      cluster_offset.push_back(  xCooPos_in  );
      cluster_offset.push_back(  yCooPos_in  );
      
      return cluster_offset;
}


void  EUTelCorrelator::FillHotPixelMap(LCEvent *event)
{
    LCCollectionVec *hotPixelCollectionVec = 0;
    try 
    {
      hotPixelCollectionVec = static_cast< LCCollectionVec* > ( event->getCollection( _hotPixelCollectionName  ) );
      streamlog_out ( DEBUG5 ) << "Hotpixel collection " << _hotPixelCollectionName.c_str() << " found" << endl; 
    }
    catch (...)
    {
      streamlog_out ( MESSAGE5 ) << "Hotpixel collection " << _hotPixelCollectionName.c_str() << " not found" << endl; 
      return;
    }

        CellIDDecoder<TrackerDataImpl> cellDecoder( hotPixelCollectionVec );
	
        for(int i=0; i<  hotPixelCollectionVec->getNumberOfElements(); i++)
        {
           TrackerDataImpl* hotPixelData = dynamic_cast< TrackerDataImpl *> ( hotPixelCollectionVec->getElementAt( i ) );
	   SparsePixelType  type         = static_cast<SparsePixelType> (static_cast<int> (cellDecoder( hotPixelData )["sparsePixelType"]));
	   int sensorID              = static_cast<int > ( cellDecoder( hotPixelData )["sensorID"] );
       }

}
 

bool EUTelCorrelator::hitContainsHotPixels( TrackerHitImpl   * hit) 
{

        try
        {
            LCObjectVec clusterVector = hit->getRawHits();

						 // EUTelVirtualCluster * cluster;
            if ( hit->getType() == kEUTelBrickedClusterImpl ) {
               // fixed cluster implementation. Remember it
               //  can come from
               //  both RAW and ZS data
   
//              cluster = new EUTelBrickedClusterImpl(static_cast<TrackerDataImpl *> ( clusterVector[0] ) );
                
            } else if ( hit->getType() == kEUTelDFFClusterImpl ) {
              
              // fixed cluster implementation. Remember it can come from
              // both RAW and ZS data
//              cluster = new EUTelDFFClusterImpl( static_cast<TrackerDataImpl *> ( clusterVector[0] ) );
            } else if ( hit->getType() == kEUTelFFClusterImpl ) {
              
              // fixed cluster implementation. Remember it can come from
              // both RAW and ZS data
//              cluster = new EUTelFFClusterImpl( static_cast<TrackerDataImpl *> ( clusterVector[0] ) );
            } 
       }
       catch(...)
       { 
          printf("something went wrong in EUTelCorrelator::hitContainsHotPixels \n");
          // if anything went wrong in the above return FALSE, meaning do not skip this hit
          return 0;
       }
 
       // if none of the above worked return FALSE, meaning do not skip this hit
       return 0;

}



#endif // USE_GEAR
