// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Silvia Bonfanti, Uni. Insubria  <mailto:silviafisica@gmail.com>
// Author Loretta Negrini, Uni. Insubria  <mailto:loryneg@gmail.com>
// Version $Id: EUTelCorrelator.cc,v 1.16 2009/07/29 09:36:49 gelin Exp $
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
#include "EUTelCorrelator.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelExceptions.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"

// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ITree.h>
#endif

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// lcio includes <.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <UTIL/CellIDDecoder.h>

// system includes <>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace marlin;
using namespace eutelescope;

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelCorrelator::_clusterXCorrelationHistoName   = "ClusterXCorrelationHisto";
std::string EUTelCorrelator::_clusterYCorrelationHistoName   = "ClusterYCorrelationHisto";
std::string EUTelCorrelator::_hitXCorrelationHistoName       = "HitXCorrelatioHisto";
std::string EUTelCorrelator::_hitYCorrelationHistoName       = "HitYCorrelationHisto";
#endif

EUTelCorrelator::EUTelCorrelator () : Processor("EUTelCorrelator") {

  // modify processor description
  _description =
    "EUTelCorrelator fills histograms with correlation plots";

  registerInputCollection(LCIO::TRACKERPULSE,"InputClusterCollectionName",
                          "Cluster (pulse) collection name",
                          _inputClusterCollectionName, string ( "cluster" ));

  registerInputCollection(LCIO::TRACKERHIT,"InputHitCollectionName",
                          "Hit collection name",
                          _inputHitCollectionName, string ( "hit" ));


}


void EUTelCorrelator::init() {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // clear the sensor ID vector
  _sensorIDVec.clear();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

  // check if the GEAR manager pointer is not null!
  if ( Global::GEAR == 0x0 ) {
    streamlog_out ( ERROR4 ) <<  "The GearMgr is not available, for an unknown reason." << endl;
    exit(-1);
  }

  _siPlanesParameters  = const_cast<gear::SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*>
    ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

  _siPlaneZPosition = new double[ _siPlanesLayerLayout->getNLayers() ];
  for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++ ) {
    _siPlaneZPosition[ iPlane ] = _siPlanesLayerLayout->getLayerPositionZ(iPlane);
    int sensorID = _siPlanesLayerLayout->getID( iPlane );

    _sensorIDVec.push_back( sensorID );
    _minX[ sensorID ] = 0;
    _minY[ sensorID ] = 0;
    _maxX[ sensorID ] = _siPlanesLayerLayout->getSensitiveNpixelX( iPlane ) - 1;
    _maxY[ sensorID ] = _siPlanesLayerLayout->getSensitiveNpixelY( iPlane ) - 1;


  }


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
    streamlog_out ( ERROR1 ) <<  "Error during the geometry consistency check: " << endl
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

  if (_iEvt % 10 == 0)
    streamlog_out( MESSAGE4 ) << "Processing event "
                              << setw(6) << setiosflags(ios::right) << event->getEventNumber() << " in run "
                              << setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber()
                              << setfill(' ') << " (Total = " << setw(10) << _iEvt << ")"
                              << resetiosflags(ios::left) << endl;
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
  if ( !_isInitialize ) {

    try {
      // let's check if we have cluster collections
      event->getCollection( _inputClusterCollectionName );

      _hasClusterCollection = true;

    } catch ( lcio::Exception& e ) {

      _hasClusterCollection = false;
    }

    try {
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

    if ( _hasClusterCollection ) {

      LCCollectionVec * inputClusterCollection   = static_cast<LCCollectionVec*>
        (event->getCollection( _inputClusterCollectionName ));

      CellIDDecoder<TrackerPulseImpl>  pulseCellDecoder( inputClusterCollection );

      // we have an external detector where we consider a cluster each
      // time (external cluster) that is correlated with another
      // detector's clusters (internal cluster)

      for ( size_t iExt = 0 ; iExt < inputClusterCollection->size() ; ++iExt ) {

        TrackerPulseImpl * externalPulse = static_cast< TrackerPulseImpl * >
          ( inputClusterCollection->getElementAt( iExt ) );

        EUTelVirtualCluster  * externalCluster;

        ClusterType type = static_cast<ClusterType>
          (static_cast<int>((pulseCellDecoder(externalPulse)["type"])));
        // we check that the type of cluster is ok
        
        if ( type == kEUTelDFFClusterImpl ) {
          externalCluster = new EUTelDFFClusterImpl( static_cast<TrackerDataImpl*>
                                                    ( externalPulse->getTrackerData()) );
        } else if ( type == kEUTelFFClusterImpl ) {
          externalCluster = new EUTelFFClusterImpl( static_cast<TrackerDataImpl*>
                                                    ( externalPulse->getTrackerData()) );
            
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
            externalCluster = new EUTelSparseClusterImpl< EUTelSimpleSparsePixel >
              ( static_cast<TrackerDataImpl *> ( externalPulse->getTrackerData()  ) );
          } else {
            streamlog_out ( ERROR4 ) << "Unknown pixel type.  Sorry for quitting." << endl;
            throw UnknownDataTypeException("Pixel type unknown");
          }

        } else {
          streamlog_out ( ERROR4 ) <<  "Unknown cluster type. Sorry for quitting" << endl;
          throw UnknownDataTypeException("Cluster type unknown");
        }

        int externalSensorID = pulseCellDecoder( externalPulse ) [ "sensorID" ] ;

        float externalXCenter;
        float externalYCenter;

        // we catch the coordinates of the external seed

        externalCluster->getCenterOfGravity( externalXCenter, externalYCenter ) ;

        for ( size_t iInt = 0;  iInt <  inputClusterCollection->size() ; ++iInt ) {

          TrackerPulseImpl * internalPulse = static_cast< TrackerPulseImpl * >
            ( inputClusterCollection->getElementAt( iInt ) );

          EUTelVirtualCluster  * internalCluster;

          ClusterType type = static_cast<ClusterType>
            (static_cast<int>((pulseCellDecoder(internalPulse)["type"])));

          // we check that the type of cluster is ok
          if ( type == kEUTelDFFClusterImpl ) {
            internalCluster = new EUTelDFFClusterImpl( static_cast<TrackerDataImpl*>
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

          } else {
            streamlog_out ( ERROR4 ) <<  "Unknown cluster type. Sorry for quitting" << endl;
            throw UnknownDataTypeException("Cluster type unknown");
          }
          int internalSensorID = pulseCellDecoder( internalPulse ) [ "sensorID" ] ;

          if ( internalSensorID != externalSensorID ) {

            float internalXCenter;
            float internalYCenter;

            // we catch the coordinates of the internal seed

            internalCluster->getCenterOfGravity( internalXCenter, internalYCenter ) ;

            streamlog_out ( DEBUG ) << "Filling histo " << externalSensorID << " " << internalSensorID << endl;


            // we input the coordinates in the correlation matrix, one
            // for each type of coordinate: X and Y

            _clusterXCorrelationMatrix[ externalSensorID ][ internalSensorID ]->
              fill( externalXCenter, internalXCenter );
            _clusterYCorrelationMatrix[ externalSensorID ][ internalSensorID ]->
              fill( externalYCenter, internalYCenter );

          } // endif

          delete internalCluster;

        } // internal loop

        delete externalCluster;
      } // external loop
    } // endif hasCluster



    if ( _hasHitCollection ) {

      LCCollectionVec * inputHitCollection = static_cast< LCCollectionVec *>
        ( event->getCollection( _inputHitCollectionName )) ;

      for ( size_t iExt = 0 ; iExt < inputHitCollection->size(); ++iExt ) {

        // this is the external hit
        TrackerHitImpl * externalHit = static_cast< TrackerHitImpl * > ( inputHitCollection->
                                                                         getElementAt( iExt ) );

        int externalSensorID = guessSensorID( externalHit );

        double * externalPosition;
        externalPosition = (double *) externalHit->getPosition();

        for ( size_t iInt = 0; iInt < inputHitCollection->size(); ++iInt ) {

          TrackerHitImpl  * internalHit = static_cast< TrackerHitImpl * > ( inputHitCollection->
                                                                            getElementAt( iInt ) );

          int internalSensorID = guessSensorID( internalHit );

          if ( internalSensorID != externalSensorID ) {

            double * internalPosition;
            internalPosition = (double *) internalHit->getPosition(  );

            _hitXCorrelationMatrix[ externalSensorID ] [ internalSensorID ] ->
              fill ( externalPosition[0] , internalPosition[0] ) ;

            _hitYCorrelationMatrix[ externalSensorID ] [ internalSensorID ] ->
              fill( externalPosition[1], internalPosition[1] );

          }

        }

      }
    }

  } catch (DataNotAvailableException& e  ) {

    streamlog_out  ( WARNING2 ) <<  "No input collection found on event " << event->getEventNumber()
                                << " in run " << event->getRunNumber() << endl;
  }


#endif

}

void EUTelCorrelator::end() {

  streamlog_out ( MESSAGE4 )  << "Successfully finished" << endl;
  delete [] _siPlaneZPosition;
}

void EUTelCorrelator::bookHistos() {

  if ( !_hasClusterCollection && !_hasHitCollection ) return ;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  try {

    streamlog_out ( MESSAGE4 ) <<  "Booking histograms" << endl;

    // create all the directories first
    vector< string > dirNames;

    if ( _hasClusterCollection ) {
      dirNames.push_back ("ClusterX");
      dirNames.push_back ("ClusterY");
    }

    if ( _hasHitCollection ) {
      dirNames.push_back ("HitX");
      dirNames.push_back ("HitY");
    }

    for ( size_t iPos = 0 ; iPos < dirNames.size() ; iPos++ ) {

      AIDAProcessor::tree(this)->mkdir( dirNames[iPos].c_str() ) ;

    }


    string tempHistoName;
    string tempHistoTitle ;


    for ( size_t r = 0 ; r < _sensorIDVec.size(); ++r ) {

      int row = _sensorIDVec.at( r );
      
      map< unsigned int , AIDA::IHistogram2D * > innerMapXCluster;
      map< unsigned int , AIDA::IHistogram2D * > innerMapYCluster;

      map< unsigned int , AIDA::IHistogram2D * > innerMapXHit;
      map< unsigned int , AIDA::IHistogram2D * > innerMapYHit;


      for ( size_t c = 0 ; c < _sensorIDVec.size(); ++c ) {

        int col = _sensorIDVec.at( c );

        if ( col != row ) {

          //we create histograms for X and Y Cluster correlation
          if ( _hasClusterCollection ) {
            tempHistoName = "ClusterX/" + _clusterXCorrelationHistoName + "_d"
              + to_string( row ) + "_d" + to_string( col );

            streamlog_out( DEBUG ) << "Booking histo " << tempHistoName << endl;

            int     xBin = _maxX[ row ] - _minX[ row ] + 1;
            double  xMin = static_cast<double >(_minX[ row ]) - 0.5;
            double  xMax = static_cast<double >(_maxX[ row ]) + 0.5;
            int     yBin = _maxX[ col ] - _minX[ col ] + 1;
            double  yMin = static_cast<double >(_minX[ col ]) - 0.5;
            double  yMax = static_cast<double >(_maxX[ col ]) + 0.5;

            AIDA::IHistogram2D * histo2D =
              AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(),
                                                                        xBin, xMin, xMax, yBin, yMin, yMax );

            tempHistoTitle = "XClusterCorrelation_d" + to_string( row ) + "_d" + to_string( col );
            histo2D->setTitle( tempHistoTitle.c_str() );
            innerMapXCluster[ col  ] =  histo2D ;
            tempHistoName =  "ClusterY/" +  _clusterYCorrelationHistoName + "_d" + to_string( row ) + "_d" + to_string( col );

            streamlog_out( DEBUG ) << "Booking histo " << tempHistoName << endl;

            xBin = _maxY[ row ] - _minY[ row ] + 1;
            xMin = static_cast<double >(_minY[ row ]) - 0.5;
            xMax = static_cast<double >(_maxY[ row ]) + 0.5;
            yBin = _maxY[ col ] - _minY[ col ] + 1;
            yMin = static_cast<double >(_minY[ col ]) - 0.5;
            yMax = static_cast<double >(_maxY[ col ]) + 0.5;


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

            double safetyFactor = 2.0; // 2 should be enough because it
            // means that the sensor is wrong
            // by all its size.
            double rowMin = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionX( r ) -
                                             ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeX ( r ) ));
            double rowMax = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionX( r ) +
                                             ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeX ( r )));

            double colMin = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionX( c ) -
                                             ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeX ( c )));
            double colMax = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionX( col ) +
                                             ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeX ( c )) );
            
            //lets limit the memory usage
            int colNBin = 512;
            int rowNBin = 512;
            
            if(_siPlanesLayerLayout->getSensitiveNpixelX( c ) < 255)
              colNBin = static_cast< int > ( safetyFactor ) * _siPlanesLayerLayout->getSensitiveNpixelX( c );

            if(_siPlanesLayerLayout->getSensitiveNpixelX( r ) < 255)
              rowNBin = static_cast< int > ( safetyFactor ) * _siPlanesLayerLayout->getSensitiveNpixelX( r );

            tempHistoName =  "HitX/" +  _hitXCorrelationHistoName + "_d" + to_string( row ) + "_d" + to_string( col );
            streamlog_out( DEBUG ) << "Booking histo " << tempHistoName << endl;

            string tempHistoTitle =  "HitX/" +  _hitXCorrelationHistoName + "_d" + to_string( row ) + "_d" +  to_string( col );

            AIDA::IHistogram2D * histo2D =
              AIDAProcessor::histogramFactory( this )->createHistogram2D( tempHistoName.c_str(), colNBin, colMin, colMax,
                                                                          rowNBin, rowMin, rowMax);
            histo2D->setTitle( tempHistoTitle.c_str() );

            innerMapXHit[ col  ] =  histo2D ;

            // now the hit on the Y direction
            rowMin = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionY( r ) -
                                      ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeY ( r ) ));
            rowMax = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionY( r ) +
                                      ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeY ( r )));

            colMin = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionY( c ) -
                                      ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeY ( c )));
            colMax = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionY( c ) +
                                      ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeY ( c )) );
            
            colNBin = 512;
            rowNBin = 512;
            if(_siPlanesLayerLayout->getSensitiveNpixelY( c ) < 255)
              colNBin = static_cast< int > ( safetyFactor ) * _siPlanesLayerLayout->getSensitiveNpixelY( c );

            if(_siPlanesLayerLayout->getSensitiveNpixelY( r ) < 255)
              rowNBin = static_cast< int > ( safetyFactor ) * _siPlanesLayerLayout->getSensitiveNpixelY( r );

            tempHistoName =  "HitY/" + _hitYCorrelationHistoName + "_d" + to_string( row ) + "_d" + to_string( col );
            streamlog_out( DEBUG ) << "Booking cloud " << tempHistoName << endl;
            tempHistoTitle = "YHitCorrelation_d" + to_string( row ) + "_d" + to_string( col ) ;
            histo2D =
              AIDAProcessor::histogramFactory( this )->createHistogram2D( tempHistoName.c_str(), colNBin, colMin, colMax,
                                                                          rowNBin, rowMin, rowMax);
            histo2D->setTitle( tempHistoTitle.c_str() );

            innerMapYHit[  col  ] =  histo2D ;

          }

        } else {

          if ( _hasClusterCollection ) {
            innerMapXCluster[ col ] = NULL ;
            innerMapYCluster[ col ] = NULL ;
            
          }

          if ( _hasHitCollection ) {
            innerMapXHit[ col ] = NULL ;
            innerMapYHit[ col  ] = NULL ;
          }

        }

      }

      if ( _hasClusterCollection ) {
        _clusterXCorrelationMatrix[ row ] = innerMapXCluster  ;
        _clusterYCorrelationMatrix[ row ] = innerMapYCluster  ;

        
      }

      if ( _hasHitCollection ) {
        _hitXCorrelationMatrix[ row ] = innerMapXHit;
        _hitYCorrelationMatrix[ row ] = innerMapYHit;
      }
    }

  } catch (lcio::Exception& e ) {

    streamlog_out ( ERROR1 ) << "No AIDAProcessor initialized. Sorry for quitting..." << endl;
    exit( -1 );

  }
#endif
}


int EUTelCorrelator::guessSensorID( TrackerHitImpl * hit ) {

  int sensorID = -1;
  double minDistance =  numeric_limits< double >::max() ;
  double * hitPosition = const_cast<double * > (hit->getPosition());

  for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); ++iPlane ) {
    double distance = std::abs( hitPosition[2] - _siPlaneZPosition[ iPlane ] );
    if ( distance < minDistance ) {
      minDistance = distance;
      sensorID = _siPlanesLayerLayout->getID( iPlane );
    }
  }
  if ( minDistance > 5 /* mm */ ) {
    // advice the user that the guessing wasn't successful 
    streamlog_out( WARNING3 ) << "A hit was found " << minDistance << " mm far from the nearest plane\n"
      "Please check the consistency of the data with the GEAR file" << endl;
  }

  return sensorID;
}

#endif // USE_GEAR
