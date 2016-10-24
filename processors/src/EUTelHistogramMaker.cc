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

// eutelescope includes ".h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelGenericSparsePixel.h"
#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"
#include "EUTelHistogramManager.h"
#include "EUTelMatrixDecoder.h"
#include "EUTelHistogramMaker.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/AIDAProcessor.h"

// lcio includes <.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <Exceptions.h>

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes <.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>

// system includes <>
#include <string>
#include <sstream>
#include <vector>
#include <memory>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;


// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelHistogramMaker::_clusterSignalHistoName      = "clusterSignal";
std::string EUTelHistogramMaker::_seedSignalHistoName         = "seedSignal";
std::string EUTelHistogramMaker::_hitMapHistoName             = "hitMap";
std::string EUTelHistogramMaker::_seedSNRHistoName            = "seedSNR";
std::string EUTelHistogramMaker::_clusterNoiseHistoName       = "clusterNoise";
std::string EUTelHistogramMaker::_clusterSNRHistoName         = "clusterSNR";
std::string EUTelHistogramMaker::_eventMultiplicityHistoName  = "eventMultiplicity";
std::string EUTelHistogramMaker::_clusterNumberOfHitPixelName  = "numberofhitpixel";
#endif

EUTelHistogramMaker::EUTelHistogramMaker () : Processor("EUTelHistogramMaker") {

  // modify processor description
  _description =
    "EUTelHistogramMaker fills reference and control histograms";

  // first of all we need to register the input collection
  registerInputCollection (LCIO::TRACKERPULSE, "PulseCollectionName",
                           "Input tracker pulse collection",
                           _pulseCollectionName, string("pulse"));

  // since v00-00-09 the noise and status collections are compulsory.
  std::vector<std::string > NoiseCollectionNameVecExample;
  NoiseCollectionNameVecExample.push_back("noise");
  std::vector<std::string > StatusCollectionNameVecExample;
  StatusCollectionNameVecExample.push_back("status");
  //LCIO::TRACKERDATA
  registerInputCollections(LCIO::TRACKERDATA, "NoiseCollectionName", "The name of the noise collections",
                          _noiseCollectionName, NoiseCollectionNameVecExample );

  registerInputCollections(LCIO::TRACKERRAWDATA, "StatusCollectionName","The name of the status collections.",
                          _statusCollectionName, StatusCollectionNameVecExample );

  registerProcessorParameter("HistoInfoFileName", "This is the name of the histogram information file",
                             _histoInfoFileName, string( "histoinfo.xml" ) );

  IntVec clusterNxNExample;
  clusterNxNExample.push_back(3);
  clusterNxNExample.push_back(5);

  registerOptionalParameter("ClusterNxN", "The list of cluster NxN to be filled."
                            "For example 3 means filling the 3x3 histogram spectrum",
                            _clusterSpectraNxNVector, clusterNxNExample);

  IntVec clusterNExample;
  clusterNExample.push_back(4);
  clusterNExample.push_back(9);
  clusterNExample.push_back(14);
  clusterNExample.push_back(19);
  clusterNExample.push_back(25);
  registerOptionalParameter("ClusterN", "The list of cluster N to be filled."
                            "For example 7 means filling the cluster spectra with the 7 most significant pixels",
                            _clusterSpectraNVector, clusterNExample );

  _isFirstEvent = true;

}


void EUTelHistogramMaker::init () {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  _iRun = 0;

  // set to false the geometry ready flag, i.e. force a geometry
  // initialization as soon as possible.
  _isGeometryReady = false;

  // by default fill also the noise related histograms.
  _noiseHistoSwitch = true;


}

void EUTelHistogramMaker::processRunHeader (LCRunHeader * rdr) {

  // to make things easier re-cast the input header to the EUTelRunHeaderImpl
  std::unique_ptr<EUTelRunHeaderImpl> runHeader = std::make_unique<EUTelRunHeaderImpl>(rdr);
  runHeader->addProcessor( type() );

  // increment the run counter
  ++_iRun;
}

void EUTelHistogramMaker::initializeGeometry( LCEvent * event ) {

  // initialize the geometry filling the following information:
  //
  // min, max pixel along X and Y
  // _ancillaryMap
  // sensorIDVec
  //
  // all these information are taken from the noise collection
  streamlog_out( MESSAGE2 ) << "Initializing geometry... "<< endl;

  try {
    for(size_t i = 0; i < _noiseCollectionName.size();++i)
      {
        LCCollectionVec * collection = dynamic_cast< LCCollectionVec * > ( event->getCollection( _noiseCollectionName[i] ) );
        _noOfDetector   += collection->size();
        CellIDDecoder< TrackerDataImpl > decoder( collection );
        for ( size_t iDetector = 0 ; iDetector < collection->size(); ++iDetector ) {
          TrackerDataImpl * data = dynamic_cast< TrackerDataImpl * > ( collection->getElementAt( iDetector ) );

          int sensorID = decoder( data ) [ "sensorID" ];
          // the look up tables
          _ancillaryMap.insert( make_pair( sensorID, iDetector ) );
          _sensorIDVec.push_back( sensorID );

          // the boundaries
          _minX.insert( make_pair( sensorID, decoder( data ) [ "xMin" ] ) );
          _minY.insert( make_pair( sensorID, decoder( data ) [ "yMin" ] ) );
          _maxX.insert( make_pair( sensorID, decoder( data ) [ "xMax" ] ) );
          _maxY.insert( make_pair( sensorID, decoder( data ) [ "yMax" ] ) );

        }
      }
  } catch ( lcio::DataNotAvailableException ) {
    streamlog_out( WARNING2 ) << "Unable to initialize the geometry with the current event. Trying with the next one" << endl;
    _isGeometryReady = false;
    throw SkipEventException( this ) ;
  }

  _isGeometryReady = true;

}


void EUTelHistogramMaker::processEvent (LCEvent * evt) {

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  if ( !_isGeometryReady ) {
    initializeGeometry(evt);
  }


  EUTelEventImpl * eutelEvent = static_cast<EUTelEventImpl*> (evt);
  EventType type              = eutelEvent->getEventType();

  if ( type == kEORE ) {
    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return ;
  } else if ( type == kUNKNOWN ) {
    streamlog_out ( WARNING2 ) <<  "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }

  std::vector<LCCollectionVec *> noiseCollectionVec; //noiseCollectionVec = 0x0, * statusCollectionVec = 0x0;
  std::vector<LCCollectionVec *> statusCollectionVec;
  if ( _noiseHistoSwitch ) {
   
      
    //noiseCollectionVec  = dynamic_cast<LCCollectionVec *> (
    //evt->getCollection( _noiseCollectionName ) );
    for(size_t i = 0; i < _noiseCollectionName.size(); ++i)
      {
        try {
          noiseCollectionVec.push_back(dynamic_cast<LCCollectionVec *> (
                                                                        evt->getCollection( _noiseCollectionName[i] ) ));
        }
        catch (lcio::DataNotAvailableException& e) {
          streamlog_out ( ERROR1 ) <<  e.what() << endl << "Switching off the noise histogram filling and continuing" << endl;
          _noiseHistoSwitch &= false;
        }
        
        try {
          statusCollectionVec.push_back(dynamic_cast<LCCollectionVec *> ( evt->getCollection( _statusCollectionName[i] ) ));
        } catch (lcio::DataNotAvailableException& e) {
          streamlog_out ( ERROR1 ) << e.what() << endl  << "Switching off the noise histogram filling and continuing" << endl;
          _noiseHistoSwitch &= false;
        }
      }
  }

  if ( isFirstEvent() ) {
    bookHistos();
    _isFirstEvent = false;
  }

  try {

    LCCollectionVec * pulseCollectionVec = dynamic_cast<LCCollectionVec*>  (evt->getCollection(_pulseCollectionName));
    CellIDDecoder<TrackerPulseImpl> cellDecoder(pulseCollectionVec);

    // prepare and reset the hit counter
    map<int, int> eventCounterMap;
    for ( int iPulse = 0; iPulse < pulseCollectionVec->getNumberOfElements(); iPulse++ ) {
      TrackerPulseImpl * pulse = dynamic_cast<TrackerPulseImpl*> ( pulseCollectionVec->getElementAt(iPulse) );
      ClusterType        type  = static_cast<ClusterType> ( static_cast<int> ( cellDecoder(pulse)["type"] ));
      SparsePixelType    pixelType = static_cast<SparsePixelType> (0);

      EUTelVirtualCluster * cluster;
      if ( type == kEUTelDFFClusterImpl ) {
        cluster = new EUTelDFFClusterImpl ( static_cast<TrackerDataImpl*> ( pulse->getTrackerData() ) );
      }
      else if ( type == kEUTelBrickedClusterImpl ) {
        cluster = new EUTelBrickedClusterImpl ( static_cast<TrackerDataImpl*> ( pulse->getTrackerData() ) );
      }
      else if ( type == kEUTelFFClusterImpl ) {
        
        cluster = new EUTelFFClusterImpl ( static_cast<TrackerDataImpl*> ( pulse->getTrackerData() ) );
        
      } else if ( type == kEUTelSparseClusterImpl ) {
        LCCollectionVec * sparseClusterCollectionVec = dynamic_cast < LCCollectionVec * > (evt->getCollection("original_zsdata"));
        TrackerDataImpl * oneCluster = dynamic_cast<TrackerDataImpl*> (sparseClusterCollectionVec->getElementAt( 0 ));
        CellIDDecoder<TrackerDataImpl > anotherDecoder(sparseClusterCollectionVec);
        pixelType = static_cast<SparsePixelType> ( static_cast<int> ( anotherDecoder( oneCluster )["sparsePixelType"] ));

        if ( pixelType == kEUTelGenericSparsePixel ) {
          cluster = new EUTelSparseClusterImpl<EUTelGenericSparsePixel >
            ( static_cast<TrackerDataImpl *> ( pulse->getTrackerData() ) );
        } else {
          streamlog_out ( ERROR4 ) << "Unknown pixel type. Sorry for quitting." << endl;
          throw UnknownDataTypeException("Pixel type unknown");
        }

        
      }
        else {

        streamlog_out ( ERROR4) << "Unknown cluster type. Sorry for quitting" << endl;
        throw UnknownDataTypeException("Cluster type unknown");

      }

      int detectorID = cluster->getDetectorID();
      // increment of one unit the event counter for this plane
      eventCounterMap[detectorID]++;

      string tempHistoName = _clusterSignalHistoName + "_d" + to_string( detectorID );
      (dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(cluster->getTotalCharge());

      if(type == kEUTelDFFClusterImpl ) {
        tempHistoName = _clusterNumberOfHitPixelName + "_d" + to_string( detectorID );
        (dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(cluster->getTotalCharge());
      }

      tempHistoName = _seedSignalHistoName + "_d" + to_string( detectorID );
      (dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(cluster->getSeedCharge());

      vector<int >::iterator iter = _clusterSpectraNVector.begin();
      while ( iter != _clusterSpectraNVector.end() ) {
        tempHistoName = _clusterSignalHistoName + to_string( *iter ) + "_d" + to_string( detectorID ) ;
        (dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(cluster->getClusterCharge((*iter)));
        ++iter;
      }

      iter = _clusterSpectraNxNVector.begin();
      while ( iter != _clusterSpectraNxNVector.end() ) {
        tempHistoName = _clusterSignalHistoName + to_string(*iter) + "x" + to_string(*iter) + "_d" + to_string(detectorID);
        (dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(cluster->getClusterCharge((*iter), (*iter)));
        ++iter;
      }


      tempHistoName = _hitMapHistoName + "_d" + to_string(detectorID);
      int xSeed, ySeed;
      cluster->getCenterCoord(xSeed, ySeed);
      (dynamic_cast<AIDA::IHistogram2D*> (_aidaHistoMap[tempHistoName]))->fill(static_cast<double >(xSeed), static_cast<double >(ySeed), 1.);

      if ( _noiseHistoSwitch ) 
      {
        TrackerDataImpl    * noiseMatrix = NULL;
        TrackerRawDataImpl * statusMatrix = NULL;
        bool found = false;
        size_t i = 0;
        while(!found && i < noiseCollectionVec.size())
          {
            // define a noise decoder
            CellIDDecoder<TrackerDataImpl> noiseDecoder(noiseCollectionVec[i]);
            
            // get the noise TrackerDataImpl corresponding to the detector
            // under analysis and the status matrix as well
            
            if( noiseCollectionVec[i]->getNumberOfElements() > _ancillaryMap[detectorID] )
              {
                noiseMatrix  = dynamic_cast<TrackerDataImpl *>    (noiseCollectionVec[i]->getElementAt(_ancillaryMap[detectorID]) );
                statusMatrix = dynamic_cast<TrackerRawDataImpl *> (statusCollectionVec[i]->getElementAt(_ancillaryMap[detectorID]) );
                
                CellIDDecoder< TrackerDataImpl > noiseDecoder2( noiseCollectionVec[i] ) ;
                int sensorID = noiseDecoder2(  noiseMatrix  ) [ "sensorID" ];
                if(sensorID == detectorID)
                  {
                    found = true;
                  }
                else
                  {
                    i++;
                  }
              }
            else
              {
                i++;
              }
          }
        if(!found)
          {
            cout << "no noise and status maps were found for sensor " << detectorID << endl;
            exit(-1);
          }


        CellIDDecoder< TrackerDataImpl > noiseDecoder( noiseCollectionVec[i] ) ;
        EUTelMatrixDecoder noiseMatrixDecoder( noiseDecoder, noiseMatrix);
            
        vector<float > noiseValues;
        if ( type == kEUTelDFFClusterImpl )
          {
            int xClusterSize, yClusterSize;
              cluster->getClusterSize(xClusterSize, yClusterSize);
              
              for ( int yPixel = ySeed - ( yClusterSize / 2 ); yPixel <= ySeed + ( yClusterSize / 2 ); yPixel++ ) {
                for ( int xPixel = xSeed - ( xClusterSize / 2 ); xPixel <= xSeed + ( xClusterSize / 2 ); xPixel++ ) {
                  
                  // always check we are still within the sensor!!!
                  if ( ( xPixel >= noiseMatrixDecoder.getMinX() )  &&  ( xPixel <=  noiseMatrixDecoder.getMaxX()) &&
                       ( yPixel >=  noiseMatrixDecoder.getMinY() )  &&  ( yPixel <=  noiseMatrixDecoder.getMaxY()) ) {
                    unsigned int index = noiseMatrixDecoder.getIndexFromXY(xPixel, yPixel);
                    if(statusMatrix->getADCValues().size() > index )
                    {
                      bool isBad  = ( statusMatrix->getADCValues()[index] == EUTELESCOPE::BADPIXEL );
                      if ( !isBad ) {
                        noiseValues.push_back( 1.0 );
                      } else {
                        noiseValues.push_back( 0. );
                      }
                    } 
                  } else {
                    noiseValues.push_back( 0. );
                  }
                }
              }
          }
                                                 //! not sure if this is 100% correct for bricked clusters here
        else if ( type == kEUTelFFClusterImpl || type == kEUTelBrickedClusterImpl ) {
              int xClusterSize, yClusterSize;
              cluster->getClusterSize(xClusterSize, yClusterSize);
              
              for ( int yPixel = ySeed - ( yClusterSize / 2 ); yPixel <= ySeed + ( yClusterSize / 2 ); yPixel++ ) {
                for ( int xPixel = xSeed - ( xClusterSize / 2 ); xPixel <= xSeed + ( xClusterSize / 2 ); xPixel++ ) {
                  
                  // always check we are still within the sensor!!!
                  if ( ( xPixel >= noiseMatrixDecoder.getMinX() )  &&  ( xPixel <=  noiseMatrixDecoder.getMaxX()) &&
                       ( yPixel >=  noiseMatrixDecoder.getMinY() )  &&  ( yPixel <=  noiseMatrixDecoder.getMaxY()) ) {
                    unsigned int index = noiseMatrixDecoder.getIndexFromXY(xPixel, yPixel);
                    
                    // the corresponding position in the status matrix has to be HITPIXEL
                    // in the EUTelClusteringProcessor, we verify also that
                    // the pixel isHit, but this cannot be done in this
                    // processor, since the status matrix could have been reset
                    //
                    // bool isHit  = ( statusMatrix->getADCValues()[index] ==
                    // EUTELESCOPE::HITPIXEL );
                    //
                    if(statusMatrix->getADCValues().size() > index )
                    {
                       bool isBad  = ( statusMatrix->getADCValues()[index] == EUTELESCOPE::BADPIXEL );
                      if ( !isBad ) {
                        noiseValues.push_back( noiseMatrix->getChargeValues()[index] );
                      } else {
                        noiseValues.push_back( 0. );
                      }
                    }
                  } else {
                    noiseValues.push_back( 0. );
                  }
                }
              }
            } else if ( type == kEUTelSparseClusterImpl ) {
              
              if ( pixelType == kEUTelGenericSparsePixel ) {
                
                auto recasted = dynamic_cast<EUTelSparseClusterImpl<EUTelGenericSparsePixel>*> ( cluster );
		auto& pixelVec = recasted->getPixels();

		for( auto& sparsePixel: pixelVec ) {
                  int index = noiseMatrixDecoder.getIndexFromXY( sparsePixel.getXCoord(), sparsePixel.getYCoord() );
                  noiseValues.push_back( noiseMatrix->getChargeValues()[ index ] );
                }
              }
            }
            
            try {
              cluster->setNoiseValues( noiseValues );
            } catch ( IncompatibleDataSetException& e ) {
              streamlog_out ( ERROR1 )  << e.what() << endl <<  "Continuing without filling the histograms" << endl;
              _noiseHistoSwitch = false;
            }
          }
      
      
      if ( _noiseHistoSwitch ) {
        AIDA::IHistogram1D * histo;

        tempHistoName = _clusterNoiseHistoName + "_d" + to_string( detectorID );
        histo = dynamic_cast<AIDA::IHistogram1D* > ( _aidaHistoMap[tempHistoName] );
        if ( histo ) {
          histo->fill( cluster->getClusterNoise() );
        }

        tempHistoName =  _clusterSNRHistoName + "_d" + to_string( detectorID );
        histo = dynamic_cast<AIDA::IHistogram1D* > ( _aidaHistoMap[tempHistoName] );
        if ( histo ) {
          histo->fill( cluster->getClusterSNR() );
        }

        tempHistoName = _seedSNRHistoName + "_d" + to_string( detectorID );
        histo = dynamic_cast<AIDA::IHistogram1D * > ( _aidaHistoMap[tempHistoName] );
        if ( histo ) {
          histo->fill( cluster->getSeedSNR() );
        }

        tempHistoName = _clusterSNRHistoName + "_d" + to_string( detectorID ) ;
        histo = dynamic_cast<AIDA::IHistogram1D * > ( _aidaHistoMap[tempHistoName] );
        if ( histo ) {
          histo->fill( cluster->getClusterSNR() );
        }

        vector<int >::iterator iter = _clusterSpectraNxNVector.begin();
        while ( iter != _clusterSpectraNxNVector.end() ) {
          tempHistoName = _clusterSNRHistoName + to_string(*iter) + "x" + to_string(*iter) + "_d" + to_string(detectorID);
          histo = dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName] ) ;
          if ( histo ) {
            histo->fill(cluster->getClusterSNR( (*iter), (*iter) ));
          }
          ++iter;
        }

        vector<float > snrs = cluster->getClusterSNR(_clusterSpectraNVector);
        for ( unsigned int i = 0; i < snrs.size() ; i++ ) {
          tempHistoName = _clusterSNRHistoName + to_string( _clusterSpectraNVector[i] ) + "_d" + to_string( detectorID );
          histo = dynamic_cast<AIDA::IHistogram1D * > ( _aidaHistoMap[tempHistoName] ) ;
          if ( histo ) {
            histo->fill( snrs[i] );
          }
        }
      }

      delete cluster;
    }


    // fill the event multiplicity here
    string tempHistoName;
    for ( int iDetector = 0; iDetector < _noOfDetector; iDetector++ ) {
      tempHistoName = _eventMultiplicityHistoName + "_d" + to_string( _sensorIDVec[iDetector] );
      AIDA::IHistogram1D * histo = dynamic_cast<AIDA::IHistogram1D *> ( _aidaHistoMap[tempHistoName] );
      if ( histo ) {
        histo->fill( eventCounterMap[_sensorIDVec[iDetector]] );
      }
    }

  } catch( DataNotAvailableException& e ) {
    streamlog_out ( WARNING2 )  << "No input collection found on event " << evt->getEventNumber()
                                << " in run " << evt->getRunNumber() << endl;
    return ;
  }
#endif

}

void EUTelHistogramMaker::check (LCEvent * /*  evt */) {
  // nothing to check here - could be used to fill check plots in reconstruction processor
}


void EUTelHistogramMaker::end() {

  streamlog_out ( MESSAGE4 ) << "Processor finished successfully." << endl;

}

void EUTelHistogramMaker::bookHistos() {

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  // histograms are grouped in loops and detectors
  streamlog_out ( MESSAGE4 )  << "Booking histograms " << endl;

  std::unique_ptr<EUTelHistogramManager> histoMgr = std::make_unique<EUTelHistogramManager>(_histoInfoFileName);
  EUTelHistogramInfo    * histoInfo;
  bool                    isHistoManagerAvailable;

  try {
    isHistoManagerAvailable = histoMgr->init();
  } catch ( ios::failure& e) {
    streamlog_out ( ERROR1 )  << "I/O problem with " << _histoInfoFileName << endl
                              << "Continuing without histogram manager"   << endl;
    isHistoManagerAvailable = false;
  } catch ( ParseException& e ) {
    streamlog_out ( ERROR1 )  << e.what() << endl  << "Continuing without histogram manager" << endl;
    isHistoManagerAvailable = false;
  }

  string tempHistoName;
  string basePath;

  for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {

    basePath = "detector_" + to_string( _sensorIDVec.at( iDetector ) );
    AIDAProcessor::tree(this)->mkdir(basePath.c_str());
    basePath.append("/");


    tempHistoName =  _clusterSignalHistoName + "_d" + to_string( _sensorIDVec.at( iDetector ) ) ;
    int    clusterNBin  = 1000;
    double clusterMin   = 0.;
    double clusterMax   = 1000.;
    string clusterTitle = "Cluster spectrum with all pixels";
    if ( isHistoManagerAvailable ) {
      histoInfo = histoMgr->getHistogramInfo(_clusterSignalHistoName);
      if ( histoInfo ) {
        streamlog_out ( DEBUG1 )  << (* histoInfo ) << endl;
        clusterNBin = histoInfo->_xBin;
        clusterMin  = histoInfo->_xMin;
        clusterMax  = histoInfo->_xMax;
        if ( histoInfo->_title != "" ) clusterTitle = histoInfo->_title;
      }
    }
    AIDA::IHistogram1D * clusterSignalHisto =
      AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                clusterNBin,clusterMin,clusterMax);
    _aidaHistoMap.insert(make_pair(tempHistoName, clusterSignalHisto));
    clusterSignalHisto->setTitle(clusterTitle.c_str());



    tempHistoName =  _clusterNumberOfHitPixelName + "_d" + to_string( _sensorIDVec.at( iDetector ) ) ;
    
    int    nhitsNBin  = 20;
    double nhitsMin   = 0.;
    double nhitsMax   = 20.;
    string nhitsTitle = "Number of hit pixel inside the digital cluster";
    if ( isHistoManagerAvailable ) {
      histoInfo = histoMgr->getHistogramInfo(_clusterNumberOfHitPixelName);
      if ( histoInfo ) {
        streamlog_out ( DEBUG1 )  << (* histoInfo ) << endl;
        clusterNBin = histoInfo->_xBin;
        clusterMin  = histoInfo->_xMin;
        clusterMax  = histoInfo->_xMax;
        if ( histoInfo->_title != "" ) nhitsTitle = histoInfo->_title;
      }
    }
    AIDA::IHistogram1D * nhitsHisto =
      AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                nhitsNBin,nhitsMin,nhitsMax);
    _aidaHistoMap.insert(make_pair(tempHistoName, nhitsHisto));
    nhitsHisto->setTitle(nhitsTitle.c_str());








    int    clusterSNRNBin  = 300;
    double clusterSNRMin   = 0.;
    double clusterSNRMax   = 200;
    string clusterSNRTitle = "Cluster SNR";
    if ( isHistoManagerAvailable ) {
      histoInfo = histoMgr->getHistogramInfo( _clusterSNRHistoName );
      if ( histoInfo ) {
        streamlog_out ( DEBUG1 )  << (* histoInfo ) << endl;
        clusterSNRNBin = histoInfo->_xBin;
        clusterSNRMin  = histoInfo->_xMin;
        clusterSNRMax  = histoInfo->_xMax;
        if ( histoInfo->_title != "" ) clusterSNRTitle = histoInfo->_title;
      }
    }

    if ( _noiseHistoSwitch ) {
      // cluster SNR
      tempHistoName =  _clusterSNRHistoName + "_d" + to_string( _sensorIDVec.at( iDetector ) );
      AIDA::IHistogram1D * clusterSNRHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                  clusterSNRNBin, clusterSNRMin, clusterSNRMax);
      _aidaHistoMap.insert( make_pair(tempHistoName, clusterSNRHisto) ) ;
      clusterSNRHisto->setTitle(clusterSNRTitle.c_str());
    }


    vector<int >::iterator iter = _clusterSpectraNVector.begin();
    while ( iter != _clusterSpectraNVector.end() ) {

      tempHistoName = _clusterSignalHistoName + to_string( *iter ) + "_d" + to_string( _sensorIDVec.at( iDetector )) ;
      AIDA::IHistogram1D * clusterSignalNHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                  clusterNBin, clusterMin, clusterMax);
      _aidaHistoMap.insert(make_pair(tempHistoName, clusterSignalNHisto) );
      string tempTitle =  "Cluster spectrum with the " + to_string(*iter) + " most significant pixels ";
      clusterSignalNHisto->setTitle(tempTitle.c_str());

      if ( _noiseHistoSwitch ) {
        // this is for the SNR
        tempHistoName = _clusterSNRHistoName + to_string(*iter) + "_d" + to_string(_sensorIDVec.at( iDetector ));
        AIDA::IHistogram1D * clusterSNRNHisto =
          AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                    clusterSNRNBin, clusterSNRMin, clusterSNRMax);

        tempTitle = "Cluster SNR with the " + to_string(*iter) + " most significant pixels";
        _aidaHistoMap.insert( make_pair( tempHistoName, clusterSNRNHisto ) );
        clusterSNRNHisto->setTitle(tempTitle.c_str());
      }

      ++iter;
    }

    iter = _clusterSpectraNxNVector.begin();
    while ( iter != _clusterSpectraNxNVector.end() ) {

      tempHistoName = _clusterSignalHistoName + to_string(*iter) + "x" + to_string(*iter) + "_d" + to_string( _sensorIDVec.at( iDetector ) );
      AIDA::IHistogram1D * clusterSignalNxNHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                  clusterNBin, clusterMin, clusterMax);
      _aidaHistoMap.insert(make_pair(tempHistoName, clusterSignalNxNHisto) );
      string tempTitle =  "Cluster spectrum with " + to_string(*iter) + " by " + to_string(*iter) + " pixels ";
      clusterSignalNxNHisto->setTitle(tempTitle.c_str());

      if ( _noiseHistoSwitch ) {
        // then the SNR
        tempHistoName = _clusterSNRHistoName + to_string(*iter) + "x" + to_string(*iter) + "_d" + to_string(_sensorIDVec.at( iDetector ));
        AIDA::IHistogram1D * clusterSNRNxNHisto =
          AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                    clusterSNRNBin, clusterSNRMin, clusterSNRMax);
        _aidaHistoMap.insert(make_pair(tempHistoName, clusterSNRNxNHisto) );
        tempTitle =  "SNR with " + to_string(*iter) + " by " + to_string(*iter) + " pixels ";
        clusterSNRNxNHisto->setTitle(tempTitle.c_str());
      }
      ++iter;
    }

    tempHistoName = _seedSignalHistoName + "_d" + to_string( _sensorIDVec.at( iDetector ) );

    int    seedNBin  = 500;
    double seedMin   = 0.;
    double seedMax   = 500.;
    string seedTitle = "Seed pixel spectrum";
    if ( isHistoManagerAvailable ) {
      histoInfo = histoMgr->getHistogramInfo( _seedSignalHistoName );
      if ( histoInfo ) {
        streamlog_out ( DEBUG1 )  << (* histoInfo ) << endl;
        seedNBin = histoInfo->_xBin;
        seedMin  = histoInfo->_xMin;
        seedMax  = histoInfo->_xMax;
        if ( histoInfo->_title != "" ) seedTitle = histoInfo->_title;
      }
    }
    AIDA::IHistogram1D * seedSignalHisto =
      AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                seedNBin, seedMin, seedMax);
    _aidaHistoMap.insert(make_pair(tempHistoName, seedSignalHisto));
    seedSignalHisto->setTitle(seedTitle.c_str());


    if ( _noiseHistoSwitch ) {
      // seed SNR
      tempHistoName =  _seedSNRHistoName + "_d" + to_string( _sensorIDVec.at( iDetector ) );

      int    seedSNRNBin  =  300;
      double seedSNRMin   =    0.;
      double seedSNRMax   =  200.;
      string seedSNRTitle = "Seed SNR";
      if ( isHistoManagerAvailable ) {
        histoInfo = histoMgr->getHistogramInfo( _seedSNRHistoName );
        if ( histoInfo ) {
          streamlog_out ( DEBUG1 )  << (* histoInfo ) << endl;
          seedSNRNBin = histoInfo->_xBin;
          seedSNRMin  = histoInfo->_xMin;
          seedSNRMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) seedSNRTitle = histoInfo->_title;
        }
      }
      AIDA::IHistogram1D * seedSNRHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                  seedSNRNBin, seedSNRMin, seedSNRMax);
      _aidaHistoMap.insert( make_pair(tempHistoName, seedSNRHisto ));
      seedSNRHisto->setTitle(seedSNRTitle.c_str());

      tempHistoName =  _clusterNoiseHistoName + "_d" + to_string( _sensorIDVec.at( iDetector ));

      // cluster noise
      int    clusterNoiseNBin  =  300;
      double clusterNoiseMin   =    0.;
      double clusterNoiseMax   =  200.;
      string clusterNoiseTitle = "Cluster noise";
      if ( isHistoManagerAvailable ) {
        histoInfo = histoMgr->getHistogramInfo( _clusterNoiseHistoName );
        if ( histoInfo ) {
          streamlog_out ( DEBUG1 )  << (* histoInfo ) << endl;
          clusterNoiseNBin = histoInfo->_xBin;
          clusterNoiseMin  = histoInfo->_xMin;
          clusterNoiseMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) clusterNoiseTitle = histoInfo->_title;
        }
      }
      AIDA::IHistogram1D * clusterNoiseHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                  clusterNoiseNBin, clusterNoiseMin, clusterNoiseMax);
      _aidaHistoMap.insert( make_pair(tempHistoName, clusterNoiseHisto ));
      clusterNoiseHisto->setTitle(clusterNoiseTitle.c_str());

    }

    tempHistoName =  _hitMapHistoName + "_d" + to_string( _sensorIDVec.at( iDetector ));
    int     xBin = _maxX[_sensorIDVec.at( iDetector )] - _minX[_sensorIDVec.at( iDetector )] + 1;
    double  xMin = static_cast<double >(_minX[_sensorIDVec.at( iDetector )]) - 0.5;
    double  xMax = static_cast<double >(_maxX[_sensorIDVec.at( iDetector )]) + 0.5;
    int     yBin = _maxY[_sensorIDVec.at(iDetector)] - _minY[_sensorIDVec.at( iDetector )] + 1;
    double  yMin = static_cast<double >(_minY[_sensorIDVec.at( iDetector )]) - 0.5;
    double  yMax = static_cast<double >(_maxY[_sensorIDVec.at( iDetector )]) + 0.5;
    AIDA::IHistogram2D * hitMapHisto =
      AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(),
                                                                xBin, xMin, xMax,yBin, yMin, yMax);
    _aidaHistoMap.insert(make_pair(tempHistoName, hitMapHisto));
    hitMapHisto->setTitle("Hit map");

    tempHistoName = _eventMultiplicityHistoName + "_d" + to_string( _sensorIDVec.at( iDetector ) );
    int     eventMultiNBin  = 30;
    double  eventMultiMin   =  0.;
    double  eventMultiMax   = 30.;
    string  eventMultiTitle = "Event multiplicity";
    if ( isHistoManagerAvailable ) {
      histoInfo = histoMgr->getHistogramInfo(  _eventMultiplicityHistoName );
      if ( histoInfo ) {
        streamlog_out ( DEBUG1 )  << (* histoInfo ) << endl;
        eventMultiNBin  = histoInfo->_xBin;
        eventMultiMin   = histoInfo->_xMin;
        eventMultiMax   = histoInfo->_xMax;
        if ( histoInfo->_title != "" ) eventMultiTitle = histoInfo->_title;
      }
    }
    AIDA::IHistogram1D * eventMultiHisto =
      AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                eventMultiNBin, eventMultiMin, eventMultiMax);
    _aidaHistoMap.insert( make_pair(tempHistoName, eventMultiHisto) );
    eventMultiHisto->setTitle( eventMultiTitle.c_str() );
  }

#else
  streamlog_out ( MESSAGE2 )  << "No histogram produced because Marlin doesn't use AIDA" << endl;
#endif

}

