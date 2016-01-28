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


std::string EUTelCorrelator::_hitXCorrShiftHistoName             = "HitXCorrShift";
std::string EUTelCorrelator::_hitYCorrShiftHistoName             = "HitYCorrShift";
std::string EUTelCorrelator::_hitXCorrShiftProjectionHistoName   = "HitXCorrShiftProjection";
std::string EUTelCorrelator::_hitYCorrShiftProjectionHistoName   = "HitYCorrShiftProjection";

#endif

EUTelCorrelator::EUTelCorrelator () : Processor("EUTelCorrelator"), 
_histoInfoFileName("histoinfo.xml"),
_sensorIDVec()
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

  registerProcessorParameter ("Events",
                              "How many events are needed to get reasonable correlation plots (and Offset DB)? (default=1000)",
                              _events, static_cast <int> (1000) );

  registerOptionalParameter ("FixedPlane", "SensorID of fixed plane", _fixedPlaneID, 0);


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

  geo::gGeometry().initializeTGeoDescription(EUTELESCOPE::GEOFILENAME, EUTELESCOPE::DUMPGEOROOT);

  _sensorIDVec.clear();
  _sensorIDVec = geo::gGeometry().sensorIDsVec();
   for(std::vector<int>::iterator it = _sensorIDVec.begin(); it != _sensorIDVec.end(); it++) {
	_sensorIDtoZ.insert( std::make_pair( *it, static_cast<int>(it - _sensorIDVec.begin())) );
  } 

  // clear the sensor ID map
  _sensorIDVecMap.clear();
  _sensorIDtoZOrderMap.clear();


  // clear the sensor ID vector (z-axis order)
  _sensorIDVecZOrder.clear();


  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

 
  for ( size_t iin = 0 ; iin < geo::gGeometry().nPlanes(); iin++ ) 
  {           
    int sensorID = _sensorIDVec.at( iin );
 
    _minX[ sensorID ] = 0;
    _minY[ sensorID ] = 0;
    _maxX[ sensorID ] = geo::gGeometry().siPlaneXNpixels( sensorID ) - 1;
    _maxY[ sensorID ] = geo::gGeometry().siPlaneYNpixels( sensorID ) - 1;

    _maxX[ sensorID ] = geo::gGeometry().siPlaneXNpixels( sensorID ) - 1;
    _maxY[ sensorID ] = geo::gGeometry().siPlaneYNpixels( sensorID ) - 1;        

    _hitMinX[ sensorID ] =  geo::gGeometry().siPlaneXPosition( sensorID ) - 0.5*geo::gGeometry().siPlaneXSize ( sensorID ) ;
    _hitMaxX[ sensorID ] =  geo::gGeometry().siPlaneXPosition( sensorID ) + 0.5*geo::gGeometry().siPlaneXSize ( sensorID ) ;
    _hitMinY[ sensorID ] =  geo::gGeometry().siPlaneYPosition( sensorID ) - 0.5*geo::gGeometry().siPlaneYSize ( sensorID ) ;
    _hitMaxY[ sensorID ] =  geo::gGeometry().siPlaneYPosition( sensorID ) + 0.5*geo::gGeometry().siPlaneYSize ( sensorID ) ;
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


  if ( (unsigned int)runHeader->getGeoID() != geo::gGeometry().getSiPlanesLayoutID() ) {
    streamlog_out ( WARNING5 ) <<  "Error during the geometry consistency check: " << endl
                             << "The run header says the GeoID is " << runHeader->getGeoID() << endl
                             << "The GEAR description says is     " << geo::gGeometry().getSiPlanesLayoutID()
                             << endl;
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
 

/// intialise:
     _hasClusterCollection = false;
     _hasHitCollection     = false;


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

         _hasClusterCollection = false;
         streamlog_out ( WARNING ) << "NOT found " << i <<  " name " <<   _inputClusterCollectionName.c_str() << endl;

         break; 
       }
     }

     try 
     {
       // let's check if we have hit collections
 
       event->getCollection( _inputHitCollectionName ) ;
 
       _hasHitCollection = true;
       streamlog_out ( DEBUG5 ) << "found " <<   " name " <<   _inputHitCollectionName.c_str() << endl;

     } catch ( lcio::Exception& e ) {

       _hasHitCollection = false;
       streamlog_out ( DEBUG5 ) << "NOT found "  <<  " name " <<   _inputHitCollectionName.c_str() << endl;
     }

     // check if we have at least one collection.
     if ( ! _hasClusterCollection && ! _hasHitCollection  ) {

       // this is the case we didn't find any collection in this event
       return;
 
     } 

     // if the Event that we are looking is the first we create files
     // with histograms.
     if ( !_isInitialize ) 
     {
       // book histograms anyway, check that collections exist in the next clause
       bookHistos();
       _isInitialize = true;
     }


//  try {

    if ( _hasClusterCollection && !_hasHitCollection) {

      for( size_t eCol = 0; eCol < _clusterCollectionVec.size() ; eCol++ )
      {
         std::string _ExternalInputClusterCollectionName = _clusterCollectionVec[eCol];

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
           externalCluster = new EUTelSparseClusterImpl< EUTelGenericSparsePixel > ( static_cast<TrackerDataImpl *> ( externalPulse->getTrackerData()  ) );

           if( externalCluster != 0 && externalCluster->getTotalCharge() < _clusterChargeMin ) 
           {
              delete externalCluster; 
              continue;
           }
        }
	else  continue;

        int externalSensorID = pulseCellDecoder( externalPulse ) [ "sensorID" ] ;
 
        streamlog_out ( DEBUG1 ) << "externalSensorID : " << externalSensorID << " externalCluster=" << externalCluster << std::endl;

        float externalXCenter = 0.;
        float externalYCenter = 0.;

        // we catch the coordinates of the external seed

        externalCluster->getCenterOfGravity( externalXCenter, externalYCenter ) ;
        if( externalCluster->getTotalCharge() <= _clusterChargeMin ) 
        {
            delete externalCluster; 
            continue;
        }

        for( size_t iCol = 0; iCol < _clusterCollectionVec.size() ; iCol++ )
        {
          std::string _InternalInputClusterCollectionName = _clusterCollectionVec[iCol];

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

            internalCluster = new EUTelSparseClusterImpl< EUTelGenericSparsePixel > ( static_cast<TrackerDataImpl *> ( internalPulse->getTrackerData()  ) );

            if( internalCluster != 0 && internalCluster->getTotalCharge() < _clusterChargeMin )
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
                  (  _sensorIDtoZ.at(internalSensorID) == _sensorIDtoZ.at(externalSensorID) + 1 )
                  ) 
          {

              float internalXCenter = 0. ;
              float internalYCenter = 0. ;

            // we catch the coordinates of the internal seed

            internalCluster->getCenterOfGravity( internalXCenter, internalYCenter ) ;

            streamlog_out ( DEBUG5 ) << "Filling histo " << externalSensorID << " " << internalSensorID << endl;


            // we input the coordinates in the correlation matrix, one
            // for each type of coordinate: X and Y

            
            streamlog_out( MESSAGE1 )  << " ex " << externalSensorID <<" = [" << externalXCenter << ":" << externalYCenter << "]"
                                       << " in " << internalSensorID <<" = [" << internalXCenter << ":" << internalYCenter << "]" << std::endl;

            _clusterXCorrelationMatrix[ externalSensorID ][ internalSensorID ]->fill( externalXCenter, internalXCenter );
            _clusterYCorrelationMatrix[ externalSensorID ][ internalSensorID ]->fill( externalYCenter, internalYCenter );

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

        if ( hitDecoder( externalHit) ["properties"] != kHitInGlobalCoord ) {
           geo::gGeometry().local2Master( externalSensorID, etrackPointLocal, etrackPointGlobal );
        } else {
           // do nothing, already in global telescope frame 
        }
 
          
 
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

          double itrackPointLocal[]  = { internalPosition[0], internalPosition[1], internalPosition[2] };
          double itrackPointGlobal[] = { internalPosition[0], internalPosition[1], internalPosition[2] };

          if ( hitDecoder( internalHit )["properties"] != kHitInGlobalCoord ) {
             geo::gGeometry().local2Master( internalSensorID, itrackPointLocal, itrackPointGlobal );
          } else {
             // do nothing, already in global telescope frame 
          }

          if ( 
                  ( internalSensorID != getFixedPlaneID() && externalSensorID == getFixedPlaneID() )
                   ||
                  (  _sensorIDtoZ.at(internalSensorID) == _sensorIDtoZ.at(externalSensorID) + 1 )
              ) 
            {

            int iz = _sensorIDtoZ.at( internalSensorID ) ;

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
          int indexPlane = 0;
 
          indexPlane = 0; // should be always the first element, as it's filled in the externalID loop
 
          if( indexPlane >= 0 ) {
            for(int i = 0; i < (int)trackX.size();i++)
            {
              if( i == indexPlane ) continue; // skip as this one is not booked
              _hitXCorrelationMatrix[ iplane[ indexPlane ]        ] [ iplane[i]        ] -> fill ( trackX[ indexPlane ]          , trackX[i]           ) ;
              _hitYCorrelationMatrix[ iplane[ indexPlane ]        ] [ iplane[i]        ] -> fill ( trackY[ indexPlane ]          , trackY[i]           ) ;
              // assume all rotations have been done in the hitmaker processor:
              _hitXCorrShiftMatrix[ iplane[ indexPlane ]        ][ iplane[i]        ]->fill( trackX[ indexPlane ]          , trackX[ indexPlane ]          - trackX[i]          );
              _hitYCorrShiftMatrix[ iplane[ indexPlane ]        ][ iplane[i]        ]->fill( trackY[ indexPlane ]          , trackY[ indexPlane ]          - trackY[i]         );
            }
          }
        }else{
        }
 
      }
    }
//  } catch (DataNotAvailableException& e  ) {
//
//    streamlog_out  ( MESSAGE2 ) <<  "No input collection found on event " << event->getEventNumber()
//                                << " in run " << event->getRunNumber() << endl;
//  }

#endif

}

void EUTelCorrelator::end() {


 
    if( _hasHitCollection)
    {
        streamlog_out( MESSAGE5 ) << "The input CollectionVec contains HitCollection, calculating offset values " << endl;
 
        for ( size_t exx = 0 ; exx < geo::gGeometry().nPlanes(); exx++ ) 
        {           
            int exPlaneID = _sensorIDVec.at( exx );
            if( exPlaneID != getFixedPlaneID() ) continue;
            for ( size_t inn = 0 ; inn < geo::gGeometry().nPlanes(); inn++ ) 
            {
                int inPlaneID = _sensorIDVec.at( inn );
                if( inPlaneID == getFixedPlaneID() ) continue;

                if( _hitXCorrShiftMatrix[ exPlaneID ][ inPlaneID ] == 0 ) continue;
                if( _hitXCorrShiftMatrix[ exPlaneID ][ inPlaneID ]->yAxis().bins() <= 0 ) continue;


                float _heighestBinX = 0.;
                for( int ibin = 0; ibin < _hitXCorrShiftMatrix[ exPlaneID ][ inPlaneID ]->yAxis().bins(); ibin++)
                {
                    double xbin =  
                        _hitXCorrShiftProjection[ inPlaneID ]->axis().binLowerEdge(ibin)
                        +
                        _hitXCorrShiftProjection[ inPlaneID ]->axis().binWidth(ibin)/2.
                        ;
                    double _binValue = _hitXCorrShiftMatrix[ exPlaneID ][ inPlaneID ]->binEntriesY( ibin );
                    _hitXCorrShiftProjection[ inPlaneID ]->fill( xbin, _binValue );
                    if( _binValue>0)
                    if( _binValue > _heighestBinX )
                    {
                        _heighestBinX = _binValue;
                   }
                }
                
               
                float _heighestBinY = 0.;
                for( int ibin = 0; ibin < _hitYCorrShiftMatrix[ exPlaneID ][ inPlaneID ]->yAxis().bins(); ibin++)
                {
                    double xbin =  
                        _hitYCorrShiftProjection[ inPlaneID ]->axis().binLowerEdge(ibin)
                        +
                        _hitYCorrShiftProjection[ inPlaneID ]->axis().binWidth(ibin)/2.
                        ;
                    double _binValue = _hitYCorrShiftMatrix[ exPlaneID ][ inPlaneID ]->binEntriesY( ibin );
                    _hitYCorrShiftProjection[ inPlaneID ]->fill( xbin, _binValue );
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

                for( int ibin = 0; ibin < _hitXCorrShiftProjection[ inPlaneID ]->axis().bins(); ibin++)
                {
                    double ybin =  _hitXCorrShiftProjection[ inPlaneID ]->binHeight(ibin); 

                    if( ybin < _heighestBinX*0.9 ) continue;
                    double xbin =  
                        _hitXCorrShiftProjection[ inPlaneID ]->axis().binLowerEdge(ibin)
                        +
                        _hitXCorrShiftProjection[ inPlaneID ]->axis().binWidth(ibin)/2.
                        ;
                    

                    _correlationBandBinsX   += ybin;
                    _correlationBandCenterX += xbin*ybin;
                }

                
                double _correlationBandBinsY     = 0.;
                double _correlationBandCenterY   = 0.;

                for( int ibin = 0; ibin < _hitYCorrShiftMatrix[ exPlaneID ][ inPlaneID ]->yAxis().bins(); ibin++)
                {
                    double ybin =  _hitYCorrShiftProjection[ inPlaneID ]->binHeight(ibin); 
                    
                    if( ybin < _heighestBinY*0.9  ) continue;
                    double xbin =  
                        _hitYCorrShiftProjection[ inPlaneID ]->axis().binLowerEdge(ibin)
                        +
                        _hitYCorrShiftProjection[ inPlaneID ]->axis().binWidth(ibin)/2.
                        ;                    
                  
                   _correlationBandBinsY   += ybin;
                   _correlationBandCenterY += ybin*xbin;
               }

               
                streamlog_out( MESSAGE5 ) << "Hit Offset values: " ; 
                streamlog_out ( MESSAGE5 ) << " plane : " << inPlaneID << " to plane : " << exPlaneID ;
                streamlog_out ( MESSAGE5 ) << " X offset : "<<  (_correlationBandBinsX == 0. ? 0.: _correlationBandCenterX/_correlationBandBinsX) ; 
                streamlog_out ( MESSAGE5 ) << " Y offset : "<<  (_correlationBandBinsY == 0. ? 0.: _correlationBandCenterY/_correlationBandBinsY) ; 
                streamlog_out( MESSAGE5 ) << endl;


            }
        }
    }

  
    streamlog_out ( MESSAGE4 )  << "Successfully finished" << endl;
}

void EUTelCorrelator::bookHistos() {

  if ( !_hasClusterCollection && !_hasHitCollection ) return ;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  try {

    streamlog_out ( DEBUG5 ) <<  "Booking histograms" << endl;

        std::unique_ptr<EUTelHistogramManager> histoMgr = std::make_unique<EUTelHistogramManager>(_histoInfoFileName);
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
                  ( col != getFixedPlaneID() && row == getFixedPlaneID() )
                  ||
                  (  _sensorIDtoZ.at( col ) == _sensorIDtoZ.at( row ) + 1 )
                  ) 
          {


          //we create histograms for X and Y Cluster correlation
          if ( _hasClusterCollection && !_hasHitCollection) {

            //double safetyFactor = 1.0; // 2 should be enough because it
            // means that the sensor is wrong
            // by all its size.

            /////////////////////////////////////////////////
            // book X
            tempHistoName = "ClusterX/" + _clusterXCorrelationHistoName + "_d" + to_string( row ) + "_d" + to_string( col );
            streamlog_out( DEBUG5 ) << "Booking histo " << tempHistoName << endl;

            histoInfo = histoMgr->getHistogramInfo(_clusterXCorrelationHistoName);
            xBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin :  geo::gGeometry().siPlaneXNpixels(row);
            xMin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin :  0.;
            xMax  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax :  geo::gGeometry().siPlaneXNpixels(row);
            yBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin :  geo::gGeometry().siPlaneXNpixels(col);    
            yMin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin :  0.;
            yMax  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax :  geo::gGeometry().siPlaneXNpixels(col);

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
            xBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : geo::gGeometry().siPlaneYNpixels(row);
            xMin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : 0.;
            xMax  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : geo::gGeometry().siPlaneYNpixels(row);
            yBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : geo::gGeometry().siPlaneYNpixels(col);    
            yMin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : 0.;
            yMax  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax : geo::gGeometry().siPlaneYNpixels(col);

            histo2D =
              AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(),
                                                                        xBin, xMin, xMax, yBin, yMin, yMax );
            tempHistoTitle =  "ClusterY/" +  _clusterYCorrelationHistoName + "_d" + to_string( row ) + "_d" + to_string( col );
            histo2D->setTitle( tempHistoTitle.c_str()) ;

            innerMapYCluster[ col  ] =  histo2D ;
            
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
            colMin   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : -0.5* geo::gGeometry().siPlaneXSize(row);
            colMax   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax :  0.5* geo::gGeometry().siPlaneXSize(row);
            rowNBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : 100   ;    
            rowMin   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : -0.5* geo::gGeometry().siPlaneXSize(col);
            rowMax   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax :  0.5* geo::gGeometry().siPlaneXSize(col);


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
            colMin   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : -0.5* geo::gGeometry().siPlaneYSize(row);
            colMax   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax :  0.5* geo::gGeometry().siPlaneYSize(row);
            rowNBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : 100   ;    
            rowMin   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : -0.5* geo::gGeometry().siPlaneYSize(col);
            rowMax   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax :  0.5* geo::gGeometry().siPlaneYSize(col);

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
            colMin   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : -0.5* geo::gGeometry().siPlaneXSize(row);
            colMax   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax :  0.5* geo::gGeometry().siPlaneXSize(row);
            rowNBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : 100   ;    
            rowMin   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : -0.5* geo::gGeometry().siPlaneXSize(col);
            rowMax   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax :  0.5* geo::gGeometry().siPlaneXSize(col);

            histo2D = AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(), rowNBin, rowMin, rowMax, colNBin, colMin, colMax );

            tempHistoTitle =  "HitXShift/" +  _hitXCorrShiftHistoName + "_d" + to_string( row ) + "_d" + to_string( col );
            histo2D->setTitle( tempHistoTitle.c_str()) ;
            innerMapXHitShift[ col  ] =  histo2D ;


            // book Y
            tempHistoName =  "HitYShift/" +  _hitYCorrShiftHistoName + "_d" + to_string( row ) + "_d" + to_string( col );

            streamlog_out( DEBUG5 ) << "Booking histo " << tempHistoName << endl;

            histoInfo = histoMgr->getHistogramInfo(_hitYCorrShiftHistoName);
            colNBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 100   ;    
            colMin   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : -0.5* geo::gGeometry().siPlaneYSize(row);
            colMax   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax :  0.5* geo::gGeometry().siPlaneYSize(row);
            rowNBin  =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : 100   ;    
            rowMin   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : -0.5* geo::gGeometry().siPlaneYSize(col);
            rowMax   =      ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax :  0.5* geo::gGeometry().siPlaneYSize(col);

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

 
            //double safetyFactor = 1.0; // 2 should be enough because it
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

    
    double xDet_in =  internalXCenter*geo::gGeometry().siPlaneXPitch(internalSensorID);
    double yDet_in =  internalYCenter*geo::gGeometry().siPlaneYPitch(internalSensorID) ;
    double xDet_ex =  externalXCenter*geo::gGeometry().siPlaneXPitch(externalSensorID) ;
    double yDet_ex =  externalYCenter*geo::gGeometry().siPlaneYPitch(externalSensorID) ;

    double xCoo_in = internalXCenter;                      
    double yCoo_in = internalYCenter; 
   
 
    // get rotated sensors coordinates (in mm or um)
    double xPos_in =  xDet_in;
    double yPos_in =  yDet_in;
    double xPos_ex =  xDet_ex;
    double yPos_ex =  yDet_ex;

    double xCooPos_in =  xCoo_in;
    double yCooPos_in =  yCoo_in;

    // get rotated sensor coordinates (only in pixel number: col num)
   
      xPos_in +=  geo::gGeometry().siPlaneXPosition( internalSensorID ) + geo::gGeometry().siPlaneXSize ( internalSensorID )/2. ;

      xCooPos_in = xCooPos_in;

      yPos_in +=  geo::gGeometry().siPlaneYPosition( internalSensorID ) + geo::gGeometry().siPlaneYSize ( internalSensorID )/2. ;

      yCooPos_in = yCooPos_in;


      xPos_ex +=  geo::gGeometry().siPlaneXPosition( externalSensorID ) + geo::gGeometry().siPlaneXSize ( externalSensorID )/2. ;

      yPos_ex +=  geo::gGeometry().siPlaneYPosition( externalSensorID ) + geo::gGeometry().siPlaneYSize ( externalSensorID )/2. ;

                     
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



#endif // USE_GEAR
