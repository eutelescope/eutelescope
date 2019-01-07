/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelSparseClustering.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"

// eutelescope data specific
#include "EUTelSparseClusterImpl.h"
#include "EUTelTrackerDataInterfacerImpl.h"

// eutelescope geometry
#include "EUTelGenericPixGeoDescr.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

// ROOT includes
#include "TGeoBBox.h"
#include "TGeoShape.h"

// marlin includes ".h"
#include "marlin/AIDAProcessor.h"
#include "marlin/Exceptions.h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ITree.h>
#endif

// system includes <>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace lcio;
using namespace marlin;
using namespace eutelescope;

EUTelSparseClustering::EUTelSparseClustering()
    : Processor("EUTelSparseClustering"), _zsDataCollectionName(""),
      _pulseCollectionName(""), _initialPulseCollectionSize(0), _iRun(0),
      _iEvt(0), _fillHistos(false), _totalClusterMap(), _noOfDetector(0), 
      _excludedPlanes(), _isGeometryReady(false), _sensorIDVec(), _zsInputDataCollectionVec(nullptr),
      _pulseCollectionVec(nullptr), _sparseMinDistanceSquared(2) {

  _description = "EUTelSparseClustering is looking for clusters into "
                 "a calibrated pixel matrix.";

  registerInputCollection(LCIO::TRACKERDATA, "ZSDataCollectionName",
                          "Input of Zero Suppressed data",
                          _zsDataCollectionName, std::string("zsdata"));

  registerOutputCollection(LCIO::TRACKERPULSE, "PulseCollectionName",
                           "Cluster (output) collection name",
                           _pulseCollectionName, std::string("cluster"));

  registerProcessorParameter("HistogramFilling",
                             "Switch on or off the histogram filling",
                             _fillHistos, true);

  registerOptionalParameter(
      "ExcludedPlanes",
      "The list of sensor ids that have to be excluded from the clustering.",
      _excludedPlanes, std::vector<int>());

  registerProcessorParameter(
      "SparseMinDistanceSquared",
      "Minimum distance squared between sparsified pixel ( touching == 2) [integer]",
      _sparseMinDistanceSquared, 2);

  _isFirstEvent = true;
}

void EUTelSparseClustering::init() {

  // usually a good idea to do
  printParameters();

  // init new geometry
  geo::gGeometry().initializeTGeoDescription(EUTELESCOPE::GEOFILENAME,
                                             EUTELESCOPE::DUMPGEOROOT);

  // reset the run and event counters
  _iRun = 0;
  _iEvt = 0;

  // the geometry is not yet initialized, set switch to false
  _isGeometryReady = false;
}

void EUTelSparseClustering::processRunHeader(LCRunHeader *rdr) {
  auto runHeader = std::make_unique<EUTelRunHeaderImpl>(rdr);
  runHeader->addProcessor(type());
  ++_iRun;
}

void EUTelSparseClustering::initializeGeometry(LCEvent *event){

  // set the total number of detector to zero. This number can be different from
  // the one written in the gear description because
  // the input collection can contain only a fraction of all the sensors
  _noOfDetector = 0;
  _sensorIDVec.clear();

  streamlog_out(DEBUG5) << "Initializing geometry" << std::endl;

  try {
    _zsInputDataCollectionVec = dynamic_cast<LCCollectionVec *>(
        event->getCollection(_zsDataCollectionName));
    _noOfDetector += _zsInputDataCollectionVec->getNumberOfElements();
    CellIDDecoder<TrackerDataImpl> cellDecoder(_zsInputDataCollectionVec);

    for(size_t icoll = 0; icoll < _zsInputDataCollectionVec->size(); ++icoll) {
      auto data = dynamic_cast<TrackerDataImpl*>(_zsInputDataCollectionVec->getElementAt(icoll));
      _sensorIDVec.push_back(cellDecoder(data)["sensorID"]);
      _totalClusterMap.insert(std::make_pair(cellDecoder(data)["sensorID"], 0));
    }
  } catch(lcio::DataNotAvailableException&) {
    streamlog_out(DEBUG5) << "Could not find the input collection: "
                          << _zsDataCollectionName.c_str() << " !" << std::endl;
    return;
  }
  _isGeometryReady = true;
}

void EUTelSparseClustering::readCollections(LCEvent *event) {

  try {
    _zsInputDataCollectionVec = dynamic_cast<LCCollectionVec *>(
        event->getCollection(_zsDataCollectionName));
    streamlog_out(DEBUG4) << "_zsInputDataCollectionVec: "
                          << _zsDataCollectionName.c_str() << " found "
                          << std::endl;
  } catch(lcio::DataNotAvailableException &e) {
    streamlog_out(ERROR4) << "_zsInputDataCollectionVec: "
                          << _zsDataCollectionName.c_str()
                          << " not found in event " << event->getEventNumber()
                          << ". This shouldn't happen! Check your input data!"
                          << std::endl;
    throw SkipEventException(this);
  }
}

void EUTelSparseClustering::processEvent(LCEvent *event) {

  //increment event counter
  ++_iEvt;

  if(!_isGeometryReady) {
    initializeGeometry(event);
  }
  readCollections(event);

  if(_fillHistos && isFirstEvent()) {
    bookHistos();
  }

  auto evt = static_cast<EUTelEventImpl*>(event);
  if(evt->getEventType() == kEORE) {
    streamlog_out(DEBUG4) << "EORE found: nothing else to do." << std::endl;
    return;
  } else if(evt->getEventType() == kUNKNOWN) {
    streamlog_out(WARNING2) << "Event number " << evt->getEventNumber()
                            << " is of unknown type. Continue considering it "
                               "as a normal Data Event."
                            << std::endl;
  }

  //prepare a pulse collection to add all clusters found; this can be either a
  //new collection or an already existing in the event
  LCCollectionVec* pulseCollection = nullptr;
  bool pulseCollectionExists = false;
  _initialPulseCollectionSize = 0;
  try {
    pulseCollection = dynamic_cast<LCCollectionVec *>(evt->getCollection(_pulseCollectionName));
    pulseCollectionExists = true;
    _initialPulseCollectionSize = pulseCollection->size();
  } catch (lcio::DataNotAvailableException &e) {
    pulseCollection = new LCCollectionVec(LCIO::TRACKERPULSE);
  }
  if(isFirstEvent()) {
    auto& pulseCollectionParameters = pulseCollection->parameters();
    std::vector<int> sensorIDVec;
    pulseCollectionParameters.getIntVals("sensorIDs", sensorIDVec ); 
    sensorIDVec.insert( sensorIDVec.end(), _sensorIDVec.begin(), _sensorIDVec.end());
    pulseCollectionParameters.setValues("sensorIDs", sensorIDVec );
  }  

  sparseClustering(evt, pulseCollection);

  //if the pulseCollection is not empty, add it to the event
  if(!pulseCollectionExists &&
	((pulseCollection->size() != _initialPulseCollectionSize) || _initialPulseCollectionSize == 0)) {  
    evt->addCollection(pulseCollection, _pulseCollectionName);
  }

  //fill histos if flag set and pulse collection increased
  if((pulseCollection->size() != _initialPulseCollectionSize) && (_fillHistos)) {
    fillHistos(event);
  }

  if(!pulseCollectionExists &&
	(pulseCollection->size() == _initialPulseCollectionSize) && _initialPulseCollectionSize != 0) { 
   delete pulseCollection;
  }
  _isFirstEvent = false;
}

void EUTelSparseClustering::sparseClustering(LCEvent *evt, LCCollectionVec *pulseCollection) {

  //prepare some decoders
  CellIDDecoder<TrackerDataImpl> cellDecoder(_zsInputDataCollectionVec);

  bool isDummyAlreadyExisting = false;
  LCCollectionVec *sparseClusterCollectionVec = nullptr;

  try {
    sparseClusterCollectionVec =
        dynamic_cast<LCCollectionVec *>(evt->getCollection("original_zsdata"));
    isDummyAlreadyExisting = true;
  } catch (lcio::DataNotAvailableException &e) {
    sparseClusterCollectionVec = new LCCollectionVec(LCIO::TRACKERDATA);
    isDummyAlreadyExisting = false;
  }

  CellIDEncoder<TrackerDataImpl> idZSClusterEncoder(
      EUTELESCOPE::ZSCLUSTERDEFAULTENCODING, sparseClusterCollectionVec);

  //prepare an encoder also for the pulse collection
  CellIDEncoder<TrackerPulseImpl> idZSPulseEncoder(
      EUTELESCOPE::PULSEDEFAULTENCODING, pulseCollection);

  //in the zsInputDataCollectionVec we should have one TrackerData for each detector working in ZS mode
  //[START] loop over ZS detectors
  for(size_t iDetector = 0; iDetector < _zsInputDataCollectionVec->size(); iDetector++) {
    // get the TrackerData and guess which kind of sparsified data it contains
    TrackerDataImpl *zsData = dynamic_cast<TrackerDataImpl *>(
        _zsInputDataCollectionVec->getElementAt(iDetector));
    SparsePixelType type = static_cast<SparsePixelType>(
        static_cast<int>(cellDecoder(zsData)["sparsePixelType"]));
    int sensorID = static_cast<int>(cellDecoder(zsData)["sensorID"]);

    //if this is an excluded sensor, go to the next element
    bool foundExcludedSensor = false;
    for(size_t iExcluded = 0; iExcluded < _excludedPlanes.size(); ++iExcluded) {
      if(_excludedPlanes[iExcluded] == sensorID) {
        foundExcludedSensor = true;
      }
    }
    if(foundExcludedSensor)	continue;

    auto sparseData = Utility::getSparseData(zsData, type);
    auto hitPixelVec = sparseData->getPixels();
    std::vector<std::reference_wrapper<EUTelBaseSparsePixel const>> newlyAdded;

    //[START] loop over cluster candidates
    while(!hitPixelVec.empty()) {
      //prepare a TrackerData to store the cluster candidate
      std::unique_ptr<TrackerDataImpl> zsCluster = std::make_unique<TrackerDataImpl>();
      //prepare a reimplementation of sparsified cluster
      auto sparseCluster = Utility::getClusterData(zsCluster.get(), type);

      //take any pixel, e.g. the first one, add it to the cluster as well as 
      //the newly added pixels
      newlyAdded.push_back(hitPixelVec.front());
      sparseCluster->push_back(hitPixelVec.front().get());
      //now remove it from the original collection
      hitPixelVec.erase(hitPixelVec.begin());

      //now process all newly added pixels, initially this is the just previously 
      //added one but in the process of neighbour finding more new pixels are added
      while(!newlyAdded.empty()) {
        bool newlyDone = true;
        //check against all pixels in the hitPixelVec
        for(auto hitVec = hitPixelVec.begin(); hitVec != hitPixelVec.end(); ++hitVec) {
          //get the relevant infos from the newly added pixel
          auto x_add = newlyAdded.front().get().getXCoord();
          auto y_add = newlyAdded.front().get().getYCoord();

          //and the pixel we test against
          auto x_test = (hitVec->get()).getXCoord();
          auto y_test = (hitVec->get()).getYCoord();

          auto dX = x_add - x_test;
          auto dY = y_add - y_test;
          int distance = dX * dX + dY * dY;
          //if they pass the spatial cut, add them
          if(distance <= _sparseMinDistanceSquared) {
            //add them to the cluster as well as to the newly added ones
            newlyAdded.push_back(*hitVec);
            sparseCluster->push_back(*hitVec);
            //and remove it from the original collection
            hitPixelVec.erase(hitVec);
            //for test pixel there might be other neighbours, we still have to check
            newlyDone = false;
            break;
          }
        }//[END] for-loop over all pixels

        //if no neighbours are found, we can delete the pixel from the newly added;
        //tested against _ALL_ non cluster pixels, there are no other pixels
        //which could be neighbours
        if(newlyDone)	newlyAdded.erase(newlyAdded.begin());
        
      }//[END] loop over newly added

      //now process the found cluster
      if(sparseCluster->size() > 0) {
        //set the ID for this zsCluster
        idZSClusterEncoder["sensorID"] = sensorID;
        idZSClusterEncoder["sparsePixelType"] = static_cast<int>(type);
        idZSClusterEncoder["quality"] = 0;
        idZSClusterEncoder.setCellID(zsCluster.get());

        //add it to the cluster collection
        sparseClusterCollectionVec->push_back(zsCluster.get());

        //prepare a pulse for this cluster
        std::unique_ptr<TrackerPulseImpl> zsPulse = std::make_unique<TrackerPulseImpl>();
        idZSPulseEncoder["sensorID"] = sensorID;
        idZSPulseEncoder["type"] = static_cast<int>(kEUTelSparseClusterImpl);
        idZSPulseEncoder.setCellID(zsPulse.get());
        zsPulse->setTrackerData(zsCluster.release());
        pulseCollection->push_back(zsPulse.release());

        //increment the totalClusterMap
        _totalClusterMap[sensorID] += 1;
      } else {
        //cluster candidate is not passing the threshold... forget about them, 
        //the memory should be automatically cleaned by smart ptr's
      }
    }//[END] loop over cluster candidates
  }//[END] loop over ZS detectors

  //if sparseClusterCollectionVec isn't empty, add it to the current event
  if(!isDummyAlreadyExisting) {
    if(sparseClusterCollectionVec->size() != 0) {
      evt->addCollection(sparseClusterCollectionVec, "original_zsdata");
    } else {
      delete sparseClusterCollectionVec;
    }
  }
}

void EUTelSparseClustering::end() {

  streamlog_out(MESSAGE4) << "Successfully finished" << std::endl;

  std::map<int, int>::iterator iter = _totalClusterMap.begin();
  while(iter != _totalClusterMap.end()) {
    streamlog_out(MESSAGE4) << "Found " << iter->second
                            << " clusters on detector " << iter->first
                            << std::endl;
    ++iter;
  }
}

void EUTelSparseClustering::fillHistos(LCEvent *evt) {

  EUTelEventImpl *eutelEvent = static_cast<EUTelEventImpl *>(evt);
  EventType type = eutelEvent->getEventType();

  if(type == kEORE) {
    streamlog_out(DEBUG4) << "EORE found: nothing else to do." << std::endl;
    return;
  } else if(type == kUNKNOWN) {
    // if it is unknown we had already issued a warning to the user at
    // the beginning of the processEvent. If we get to here, it means
    // that the assumption that the event was a data event was
    // correct, so no harm to continue...
  }

  try {
    LCCollectionVec *_pulseCollectionVec = dynamic_cast<LCCollectionVec *>(
        evt->getCollection(_pulseCollectionName));
    CellIDDecoder<TrackerPulseImpl> cellDecoder(_pulseCollectionVec);

    std::map<int, int> eventCounterMap;

    for(int iPulse = _initialPulseCollectionSize;
         iPulse < _pulseCollectionVec->getNumberOfElements(); iPulse++) {
      TrackerPulseImpl *pulse = dynamic_cast<TrackerPulseImpl *>(
          _pulseCollectionVec->getElementAt(iPulse));
      ClusterType type = static_cast<ClusterType>(
          static_cast<int>(cellDecoder(pulse)["type"]));
      int detectorID = static_cast<int>(cellDecoder(pulse)["sensorID"]);
      //FIXME: do we need this check?
      // SparsePixelType pixelType = static_cast<SparsePixelType> (0);

      EUTelSparseClusterImpl<EUTelGenericSparsePixel> *cluster;

      if(type == kEUTelSparseClusterImpl) {
        cluster = new EUTelSparseClusterImpl<EUTelGenericSparsePixel>(
            static_cast<TrackerDataImpl *>(pulse->getTrackerData()));
      } else {
        streamlog_out(ERROR4) << "Unknown cluster type. Sorry for quitting"
                              << std::endl;
        throw UnknownDataTypeException("Cluster type unknown");
      }

      //if this key doesn't exist yet, it will be value initialized, this is
      //desired, for int this is 0!
      eventCounterMap[detectorID]++;

      //if this is an excluded sensor, go to the next element
      bool foundExcludedSensor = false;
      for(size_t iExcluded = 0; iExcluded < _excludedPlanes.size(); ++iExcluded) {
        if(_excludedPlanes[iExcluded] == detectorID) {
          foundExcludedSensor = true;
          break;
        }
      }
      if(foundExcludedSensor)	continue;

      //get the cluster size in X and Y separately and plot it:
      int xPos, yPos, xSize, ySize;
      cluster->getCenterCoord(xPos, yPos);
      cluster->getClusterSize(xSize, ySize);

      //fill the histograms
      (dynamic_cast<AIDA::IHistogram1D *>(_clusterSizeXHistos[detectorID]))
          ->fill(xSize);
      (dynamic_cast<AIDA::IHistogram1D *>(_clusterSizeYHistos[detectorID]))
          ->fill(ySize);
      (dynamic_cast<AIDA::IHistogram2D *>(_hitMapHistos[detectorID]))
          ->fill(static_cast<double>(xPos), static_cast<double>(yPos), 1.);
      (dynamic_cast<AIDA::IHistogram1D *>(_clusterSizeTotalHistos[detectorID]))
          ->fill(static_cast<int>(cluster->size()));
      (dynamic_cast<AIDA::IHistogram1D *>(_clusterSignalHistos[detectorID]))
          ->fill(cluster->getTotalCharge());
         
      delete cluster;
    }

    //fill event multiplicity here
    std::string tempHistoName;
    for(int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
      AIDA::IHistogram1D *histo = dynamic_cast<AIDA::IHistogram1D *>(
          _eventMultiplicityHistos[_sensorIDVec.at(iDetector)]);
      if(histo) {
        histo->fill(eventCounterMap[_sensorIDVec.at(iDetector)]);
      }
    }
  } catch (lcio::DataNotAvailableException &e) {
    return;
  }
}

void EUTelSparseClustering::bookHistos() {

  //histograms are grouped in loops and detectors
  streamlog_out(DEBUG5) << "Booking histograms " << std::endl;
  
  for(size_t iDetector = 0; iDetector < _sensorIDVec.size(); iDetector++) {
	//get detector information
    int sensorID = _sensorIDVec.at(iDetector);
    geo::EUTelGenericPixGeoDescr *geoDescr = (geo::gGeometry().getPixGeoDescr(sensorID));
  	
  	//with sensorID, find out specific values of detector
  	int minX, minY, maxX, maxY;	
    minX = minY = maxX = maxY = 0;
    double sizeX, sizeY, sizeZ;
    sizeX = sizeY = sizeZ = 0;
    geoDescr->getPixelIndexRange(minX, maxX, minY, maxY);
    geoDescr->getSensitiveSize(sizeX, sizeY, sizeZ);

	//create folder for current detector
    std::string basePath = "detector_" + std::to_string(sensorID);
    AIDAProcessor::tree(this)->mkdir(basePath.c_str());
    
    //create 1D histogram: cluster total size
    std::string histName_clusterSizeTotal = basePath+"/totalClusterSize_det"+ std::to_string(sensorID);
    AIDA::IHistogram1D *hist1D_clusterSizeTotal = AIDAProcessor::histogramFactory(this)->
    	createHistogram1D(histName_clusterSizeTotal, 61, -0.5, 60.5);
    hist1D_clusterSizeTotal->setTitle("Total cluster size (in hit pixels); #cluster pixel; count");
    _clusterSizeTotalHistos.insert(std::make_pair(sensorID, hist1D_clusterSizeTotal));

    //create 1D histogram: cluster signal
    std::string histName_clusterSignal = basePath+"/clusterSignal_det"+ std::to_string(sensorID);
    AIDA::IHistogram1D *hist1D_clusterSignal = AIDAProcessor::histogramFactory(this)->
    	createHistogram1D(histName_clusterSignal, 61, -0.5, 60.5);
    hist1D_clusterSignal->setTitle("Total signal per cluster (in detector specific charge unit); charge; count");
    _clusterSignalHistos.insert(std::make_pair(sensorID, hist1D_clusterSignal));
    
   	//create 1D histogram: cluster size along X
   	std::string histName_clusterSizeX = basePath+"/clusterSizeX_det"+ std::to_string(sensorID);
    AIDA::IHistogram1D *hist1D_clusterSizeX = AIDAProcessor::histogramFactory(this)->
    	createHistogram1D(histName_clusterSizeX, 21, -0.5, 20.5);
    hist1D_clusterSizeX->setTitle("Cluster size in x-direction (in pixels); size [# of pixels]; count");
    _clusterSizeXHistos.insert(std::make_pair(sensorID, hist1D_clusterSizeX));
   	
   	//create 1D histogram: cluster size along Y
   	std::string histName_clusterSizeY = basePath+"/clusterSizeY_det"+ std::to_string(sensorID);
    AIDA::IHistogram1D *hist1D_clusterSizeY = AIDAProcessor::histogramFactory(this)->
    	createHistogram1D(histName_clusterSizeY, 21, -0.5, 20.5);
    hist1D_clusterSizeY->setTitle("Cluster size in y-direction (in pixels); size [# of pixels]; count");
    _clusterSizeYHistos.insert(std::make_pair(sensorID, hist1D_clusterSizeY));
   
   	//create 2D map: hit map
   	std::string histName_hitMap = basePath+"/hitMap_det"+ std::to_string(sensorID);
    int xBin = maxX - minX + 1;
    double xMin = static_cast<double>(minX) - 0.5;
    double xMax = static_cast<double>(maxX) + 0.5;
    int yBin = maxY - minY + 1;
    double yMin = static_cast<double>(minY) - 0.5;
    double yMax = static_cast<double>(maxY) + 0.5;
    AIDA::IHistogram2D *hist2D_hitMap = AIDAProcessor::histogramFactory(this)->
    	createHistogram2D(histName_hitMap, xBin, xMin, xMax, yBin, yMin, yMax);
    hist2D_hitMap->setTitle("Pixel index hit map; x index; y index; count");
    _hitMapHistos.insert(std::make_pair(sensorID, hist2D_hitMap));
    
	//create 1D histogram: event multiplicity
   	std::string histName_eventMultiplicity = basePath+"/eventMultiplicity_det"+ std::to_string(sensorID);
    AIDA::IHistogram1D *hist1D_eventMultiplicity = AIDAProcessor::histogramFactory(this)->
    	createHistogram1D(histName_eventMultiplicity, 15, -0.5, 14.5);
    hist1D_eventMultiplicity->setTitle("Event multiplicity; multiplicity; count");
    _eventMultiplicityHistos.insert(std::make_pair(sensorID, hist1D_eventMultiplicity));
  }
  streamlog_out(DEBUG5) << "end of booking histograms " << std::endl;
}
