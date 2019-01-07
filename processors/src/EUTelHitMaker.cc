/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelHitMaker.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"

#include "EUTelBrickedClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelGenericSparseClusterImpl.h"
#include "EUTelGeometricClusterImpl.h"
#include "EUTelSimpleVirtualCluster.h"
#include "EUTelSparseClusterImpl.h"

#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelExceptions.h"

// marlin includes ".h"
#include "marlin/Global.h"
#include "marlin/Processor.h"

// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ITree.h>
#include <marlin/AIDAProcessor.h>
#endif

// lcio includes <.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>

// system includes <>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace std;
using namespace marlin;
using namespace eutelescope;

EUTelHitMaker::EUTelHitMaker()
    : Processor("EUTelHitMaker"), _pulseCollectionName(),
      _hitCollectionName(), _switchLocalCoordinates(false), _histogramSwitch(true),
      _iRun(0), _iEvt(0), _alreadyBookedSensorID() {
 
  _description = "EUTelHitMaker is responsible to translate cluster "
                 "centers from the local frame of reference \n to the external "
                 "frame of reference using the GEAR geometry description";

  registerInputCollection(LCIO::TRACKERPULSE,
			  "PulseCollectionName",
                          "Input cluster collection name",
			  _pulseCollectionName,
                          std::string(""));

  registerOutputCollection(LCIO::TRACKERHIT,
			   "HitCollectionName",
                           "Output hit collection name",
			   _hitCollectionName,
                           std::string(""));

  registerOptionalParameter("EnableLocalCoordidates",
			    "Hit coordinates are calculated in local reference frame of sensor",
			    _switchLocalCoordinates,
			    false);
      
  registerOptionalParameter("PlotHistograms",
			    "Switch for histogram plotting (default: true)",
			    _histogramSwitch,
			    true);
}

void EUTelHitMaker::init() {

  //good to do this
  printParameters();

  //reset run and event counters
  _iRun = 0;
  _iEvt = 0;

  geo::gGeometry().initializeTGeoDescription(EUTELESCOPE::GEOFILENAME,
                                             EUTELESCOPE::DUMPGEOROOT);

  _histogramSwitch = true;
}

void EUTelHitMaker::processRunHeader(LCRunHeader *rdr) {

  std::unique_ptr<EUTelRunHeaderImpl> header =
      std::make_unique<EUTelRunHeaderImpl>(rdr);
  header->addProcessor(type());

  // this is the right place also to check the geometry ID. This is a
  // unique number identifying each different geometry used at the
  // beam test. The same number should be saved in the run header and
  // in the xml file. If the numbers are different, instead of barely
  // quitting, ask the user what to do.

  if(header->getGeoID() == 0)
    streamlog_out(WARNING0)
      << "The geometry ID in the run header is set to zero." << std::endl
      << "This may mean that the GeoID parameter was not set" << std::endl;

  //increment run counter
  ++_iRun;
}

void EUTelHitMaker::processEvent(LCEvent *event) {

  ++_iEvt;

  EUTelEventImpl *evt = static_cast<EUTelEventImpl *>(event);

  if(evt->getEventType() == kEORE) {
    streamlog_out(DEBUG4) << "EORE found: nothing else to do." << std::endl;
    return;
  } else if(evt->getEventType() == kUNKNOWN) {
    streamlog_out(WARNING2) << "Event number " << evt->getEventNumber()
                            << " in run " << evt->getRunNumber()
                            << " is of unknown type. Continue considering it "
                               "as a normal Data Event."
                            << std::endl;
  }

  LCCollectionVec *pulseCollection = nullptr;
  LCCollectionVec *hitCollection = nullptr;

  try {
    pulseCollection = static_cast<LCCollectionVec *>(
        event->getCollection(_pulseCollectionName));
  } catch(DataNotAvailableException &e) {
    streamlog_out(MESSAGE2) << "No input collection " << _pulseCollectionName
                            << " found on event " << event->getEventNumber()
                            << " in run " << event->getRunNumber() << std::endl;
    return;
  }

  try {
    hitCollection = static_cast<LCCollectionVec *>(
        event->getCollection(_hitCollectionName));
  } catch(...) {
    hitCollection = new LCCollectionVec(LCIO::TRACKERHIT);
  }

  //prepare an encoder for the hit collection
  CellIDEncoder<TrackerHitImpl> idHitEncoder(EUTELESCOPE::HITENCODING,
                                             hitCollection);
  CellIDDecoder<TrackerPulseImpl> clusterCellDecoder(pulseCollection);
  CellIDDecoder<TrackerDataImpl> cellDecoder(
      EUTELESCOPE::ZSDATADEFAULTENCODING);

  int oldDetectorID = -100;
  double xSize = 0., ySize = 0.;
  double resolutionX = 0., resolutionY = 0.;
  double xPitch = 0., yPitch = 0.;

  //[START] loop over cluster
  for(int iCluster = 0; iCluster < pulseCollection->getNumberOfElements();
       iCluster++) {
    TrackerPulseImpl *pulse = dynamic_cast<TrackerPulseImpl *>(
        pulseCollection->getElementAt(iCluster));
    TrackerDataImpl *trackerData =
        dynamic_cast<TrackerDataImpl *>(pulse->getTrackerData());

    int sensorID = clusterCellDecoder(pulse)["sensorID"];
    ClusterType clusterType = static_cast<ClusterType>(
        static_cast<int>(clusterCellDecoder(pulse)["type"]));
    SparsePixelType pixelType = static_cast<SparsePixelType>(
        static_cast<int>(cellDecoder(trackerData)["sparsePixelType"]));

    //there could be several clusters belonging to the same
    //detector. So update the geometry information only if this new
    //cluster belongs to a different detector

    if(sensorID != oldDetectorID) {
      oldDetectorID = sensorID;
      //check if the histos for this sensor ID have been booked already.
      if (_alreadyBookedSensorID.find(sensorID) == _alreadyBookedSensorID.end()) {
        bookHistos(sensorID);
      }

      //all values given in mm
      resolutionX = geo::gGeometry().getPlaneXResolution(sensorID);
      resolutionY = geo::gGeometry().getPlaneYResolution(sensorID);
      xSize = geo::gGeometry().getPlaneXSize(sensorID);
      ySize = geo::gGeometry().getPlaneYSize(sensorID);
      xPitch = geo::gGeometry().getPlaneXPitch(sensorID);
      yPitch = geo::gGeometry().getPlaneYPitch(sensorID);
    }

    //LOCAL coordinate system!
    double telPos[3];
    
    //[IF] cluster type
    if(clusterType == kEUTelGenericSparseClusterImpl) {
      float xPos = 0;
      float yPos = 0;

      //for genericSparseCluster: need to know underlying pixel type
      if(pixelType == kEUTelGenericSparsePixel) {
        EUTelGenericSparseClusterImpl<EUTelGenericSparsePixel> cluster(
            trackerData);
        cluster.getCenterOfGravity(xPos, yPos);

        //for non-geometric clusters: getCenterOfGravity will return it in
        //pixel indices space, i.e have to transform into mm via the dimensions
        xPos = (xPos + 0.5) * xPitch - xSize / 2.;
        yPos = (yPos + 0.5) * yPitch - ySize / 2.;
      
      } else if(pixelType == kEUTelGeometricPixel) {
        EUTelGeometricClusterImpl cluster(trackerData);
        cluster.getGeometricCenterOfGravity(xPos, yPos);
      
      } else {
        streamlog_out(ERROR4) << "We do not support pixel type: " << pixelType
                              << " for kEUTelGenericSparseClusterImpl"
                              << std::endl;
        throw UnknownDataTypeException(
            "Pixel type not supported for kEUTelGenericSparseClusterImpl");
      }

      telPos[0] = xPos;
      telPos[1] = yPos;
      telPos[2] = 0;
    }
    //[ELSE] cluster type
	else {
      EUTelSparseClusterImpl<EUTelGenericSparsePixel> *cluster =
          new EUTelSparseClusterImpl<EUTelGenericSparsePixel>(trackerData);

      // get the position of the seed pixel (in pixel number)
      int xCluSeed = 0;
      int yCluSeed = 0;
      cluster->getSeedCoord(xCluSeed, yCluSeed);

      //with the charge center of gravity calculation, we get a shift
      //from the seed pixel center due to the charge distribution. Those
      //two numbers are the correction values in the case the Eta
      //correction is not applied.

      //!HACK TAKI:
      //! In case of a bricked cluster, we have to make sure to get the normal
      //! CoG first.
      //! The one without the global seed coordinate correction (caused by pixel
      //! rows being skewed).
      //! That one has to be eta-corrected and the global coordinate correction
      //! has to be applied on top of that!
      //! So prepare a brickedCluster pointer now:

      EUTelBrickedClusterImpl *p_tmpBrickedCluster = nullptr;
      if(clusterType == kEUTelBrickedClusterImpl) {
        p_tmpBrickedCluster = dynamic_cast<EUTelBrickedClusterImpl *>(cluster);
        if(p_tmpBrickedCluster == nullptr) {
          streamlog_out(ERROR4) << " .COULD NOT CREATE EUTelBrickedClusterImpl* !!!" << std::endl;
          throw UnknownDataTypeException(
              "COULD NOT CREATE EUTelBrickedClusterImpl* !!!");
        }
      }

      float xShift = 0.;
      float yShift = 0.;
      cluster->getCenterOfGravityShift(xShift, yShift);
      double xCorrection = static_cast<double>(xShift);
      double yCorrection = static_cast<double>(yShift);

      //rescale the pixel number in millimeter
      double xDet =
          (static_cast<double>(xCluSeed) + xCorrection + 0.5) * xPitch;
      double yDet =
          (static_cast<double>(yCluSeed) + yCorrection + 0.5) * yPitch;

      //FIXME? check the hack from Havard:
      float xCoG(0.0f), yCoG(0.0f);
      cluster->getCenterOfGravity(xCoG, yCoG);
      xDet = (xCoG + 0.5) * xPitch;
      yDet = (yCoG + 0.5) * yPitch;

      streamlog_out(DEBUG1)
          << "cluster[" << setw(4) << iCluster << "] on sensor[" << setw(3)
          << sensorID << "] at [" << setw(8) << setprecision(3) << xCoG << ":"
          << setw(8) << setprecision(3) << yCoG << "]"
          << " ->  [" << setw(8) << setprecision(3) << xDet << ":" << setw(8)
          << setprecision(3) << yDet << "]" << std::endl;

      //We have calculated the cluster hit position in terms of distance along
      //the X and Y axis.
      //However we still do not have the sensor centre as the origin of the
      //coordinate system.
      //To do this we need to deduct xSize/2 and ySize/2 for the respective
      //cluster X/Y position

      telPos[0] = xDet - xSize / 2.;
      telPos[1] = yDet - ySize / 2.;
      telPos[2] = 0.;

      delete cluster;
      cluster = nullptr;
    }//[END] cluster type

	//plot hits in the EUTelescope local frame; this frame has the
	//coordinate centre at the sensor centre
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
	if(_histogramSwitch) {
      AIDA::IHistogram2D *histo = dynamic_cast<AIDA::IHistogram2D *>(
              _hitLocalHistos[sensorID]);
      if(histo) {	              
        histo->fill(telPos[0], telPos[1]);
      } else {
        streamlog_out(ERROR1)
            << "Not able to retrieve histogram pointer for hitLocal_det" << sensorID
            << ".\nDisabling histogramming from now on " << std::endl;
        _histogramSwitch = false;
      }
    }
#endif

    if(!_switchLocalCoordinates) {
      // NOW !!
      // GLOBAL coordinate system !!!
      const double localPos[3] = {telPos[0], telPos[1], telPos[2]};
      geo::gGeometry().local2Master(sensorID, localPos, telPos);
    }

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    if(_histogramSwitch) {
      AIDA::IHistogram2D *histo = dynamic_cast<AIDA::IHistogram2D *>(
              _hitTelescopeHistos[sensorID]);
      if(histo) {	
        histo->fill(telPos[0], telPos[1]);
      } else {
        streamlog_out(ERROR1)
            << "Not able to retrieve histogram pointer for hitTelescope_det" << sensorID
            << ".\nDisabling histogramming from now on " << std::endl;
        _histogramSwitch = false;
      }
    }
#endif

    //create new hit
    TrackerHitImpl *hit = new TrackerHitImpl;
    hit->setPosition(&telPos[0]);
    float cov[TRKHITNCOVMATRIX] = {0., 0., 0., 0., 0., 0.};
    double resx = resolutionX;
    double resy = resolutionY;
    cov[0] = resx * resx; //cov(x,x)
    cov[2] = resy * resy; //cov(y,y)
    hit->setCovMatrix(cov);
    hit->setType(clusterType);
    hit->setTime(pulse->getTime());

    //prepare a LCObjectVec to store the current cluster
    LCObjectVec clusterVec;
    clusterVec.push_back(trackerData);

    //add the clusterVec to the hit
    hit->rawHits() = clusterVec;

    //determine sensorID from the cluster data
    idHitEncoder["sensorID"] = sensorID;

    //set the local/global bit flag property for the hit
    idHitEncoder["properties"] = 0; // init
    if(!_switchLocalCoordinates) idHitEncoder["properties"] = kHitInGlobalCoord;
    //store values
    idHitEncoder.setCellID(hit);
    
    //add new hit to hit collection
    hitCollection->push_back(hit);
    
  }//[END] loop over cluster

  try {
    event->getCollection(_hitCollectionName);
  } catch(...) {
    event->addCollection(hitCollection, _hitCollectionName);
  }

  if(isFirstEvent()) _isFirstEvent = false;
}

void EUTelHitMaker::end() {
  streamlog_out(MESSAGE4) << "Successfully finished" << endl;
}

void EUTelHitMaker::bookHistos(int sensorID) {

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  std::string basePath = "detector_" + to_string(sensorID);
  AIDAProcessor::tree(this)->mkdir(basePath.c_str());
  basePath = basePath + "/";
  
  //get information
  double xPosition = geo::gGeometry().getPlaneXPosition(sensorID);
  double yPosition = geo::gGeometry().getPlaneYPosition(sensorID);
  double xSize = geo::gGeometry().getPlaneXSize(sensorID);
  double ySize = geo::gGeometry().getPlaneYSize(sensorID);
  int xPixels = geo::gGeometry().getPlaneNumberOfPixelsX(sensorID);
  int yPixels = geo::gGeometry().getPlaneNumberOfPixelsY(sensorID);

  //create 2D histogram: hit map local
  std::string histName_hitLocal = "hitLocal_det" + to_string(sensorID);
  //note: in local frame the origin is at the centre of the sensor; need both, 
  //+/- x/y direction; add/subtract constant to see all hits on histogram
  double constant = 5; 
  double xMin = -xSize/2 - constant;
  double xMax = xSize/2 + constant;
  double yMin = -ySize/2 - constant;
  double yMax = ySize/2 + constant;
  int xNBin = xPixels;
  int yNBin = yPixels;

  AIDA::IHistogram2D *hist2D_hitLocal =
      AIDAProcessor::histogramFactory(this)->createHistogram2D(
       (basePath+histName_hitLocal).c_str(), xNBin, xMin, xMax, yNBin, yMin, yMax);
  if(hist2D_hitLocal) {
    hist2D_hitLocal->setTitle("Hit map (detector local frame); x position [mm]; y position [mm]");
    _hitLocalHistos.insert(std::make_pair(sensorID, hist2D_hitLocal));
  } else {
    streamlog_out(ERROR1) << "Problem booking the "
                          << (basePath + histName_hitLocal) << ".\n"
                          << "Very likely a problem with path name. Switching "
                             "off histogramming and continue w/o"
                          << std::endl;
    _histogramSwitch = false;
  }
  
  //create 2D histogram: hit map telescope
  std::string histName_hitTelescope = "hitTelescope_det" + to_string(sensorID);
  double safetyFactor = 1.2;
  xMin = safetyFactor * (xPosition - (0.5*xSize));
  xMax = safetyFactor * (xPosition + (0.5*xSize));
  yMin = safetyFactor * (yPosition - (0.5*ySize));
  yMax = safetyFactor * (yPosition + (0.5*ySize));
  xNBin = static_cast<int>(safetyFactor * xPixels);
  yNBin = static_cast<int>(safetyFactor * yPixels);
  
  AIDA::IHistogram2D *hist2D_hitTelescope =
      AIDAProcessor::histogramFactory(this)->createHistogram2D(
          (basePath+histName_hitTelescope).c_str(), xNBin, xMin, xMax, yNBin, yMin, yMax);
  if(hist2D_hitTelescope) {
    hist2D_hitTelescope->setTitle("Hit map (telescope frame); x position [mm]; y position [mm]");
    _hitLocalHistos.insert(std::make_pair(sensorID, hist2D_hitTelescope));
  } else {
    streamlog_out(ERROR1) << "Problem booking the "
                          << (basePath + histName_hitTelescope) << ".\n"
                          << "Very likely a problem with path name. Switching "
                             "off histogramming and continue w/o"
                          << std::endl;
    _histogramSwitch = false;
  }
  _alreadyBookedSensorID.insert(sensorID);
#endif //AIDA
}
