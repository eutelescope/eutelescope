/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelPreAligner.h"
#include "EUTelAlignmentConstant.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelVirtualCluster.h"

// marlin includes ".h"
#include "marlin/AIDAProcessor.h"
#include "marlin/Exceptions.h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <IO/LCWriter.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTime.h>

// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ITree.h>
#endif

// system includes <>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>

using namespace std;
using namespace eutelescope;

EUTelPreAligner::EUTelPreAligner() 
	: Processor("EUTelPreAligner") {

  _description = "EUTelPreAligner applies alignment constants to hit collection";

  registerInputCollection(LCIO::TRACKERHIT,
			  "InputHitCollectionName",
                          "The name of the input hit collection",
                          _inputHitCollectionName,
			  std::string("hit"));

  registerProcessorParameter("RequiredEvents",
			     "How many events should be used for an "
			     "approximation to the X,Y shifts "
			     "(pre-alignment)? (default=50000)",
                             _requiredEvents, 
			     50000);
                             
  registerOptionalParameter("FixedPlane",
			    "SensorID of fixed plane",
			    _fixedID,
                            0);

  registerOptionalParameter("AlignmentConstantLCIOFile",
			    "Name of LCIO database file where alignment constants will be stored",
			    _alignmentConstantLCIOFile,
			    std::string("alignment.slcio"));

  registerOptionalParameter("ResidualsXMin",
			    "Minimal values of the hit residuals in the X direction "
			    "for a correlation band (ordered according to z position).",
			    _residualsXMin,
			    std::vector<float>(6, -10.));

  registerOptionalParameter("ResidualsYMin", "Minimal values of the hit residuals in the Y direction "
			    "for a correlation band (ordered according to z position).",
			    _residualsYMin,
			    std::vector<float>(6, -10.));

  registerOptionalParameter("ResidualsXMax",
			    "Maximal values of the hit residuals in the X direction "
			    "for a correlation band (ordered according to z position).",
			    _residualsXMax,
			    std::vector<float>(6, 10.));

  registerOptionalParameter("ResidualsYMax",
			    "Maximal values of the hit residuals in the Y direction "
			    "for a correlation band (ordered according to z position).",
			    _residualsYMax,
			    std::vector<float>(6, 10.));

  registerOptionalParameter("MinNumberOfCorrelatedHits",
                            "If there are more then this number of correlated "
                            "hits (planes->track candidate) (default=5)",
                            _minNumberOfCorrelatedHits,
			    5);

  registerOptionalParameter("PlotHistograms",
			    "Switch for histogram plotting (default: true)",
			    _histogramSwitch,
			    true);

  registerOptionalParameter("DumpGEAR",
			    "Dump alignment into GEAR file instead of prealignment database",
			    _dumpGEAR,
			    false);

  registerOptionalParameter("GEARSuffix",
                            "Suffix for the new GEAR file, set to empty string "
                            "to overwrite old GEAR file (default: _pre)",
                            _GEARFileSuffix,
			    std::string("_pre"));

  registerOptionalParameter("ExcludedPlanes",
                            "The list of sensor IDs that shall be excluded.",
                            _ExcludedPlanes,
			    std::vector<int>());

  registerOptionalParameter("ExcludedPlanesXCoord",
			    "The list of sensor IDs for which the X coordinate shall be excluded.",
			    _ExcludedPlanesXCoord,
			    std::vector<int>());

  registerOptionalParameter("ExcludedPlanesYCoord",
			    "The list of sensor IDs for which the Y coordinate  shall be excluded.",
			    _ExcludedPlanesYCoord,
			    std::vector<int>());
}

void EUTelPreAligner::init() {

  //usally good idea to do
  printParameters();

  //reset run and event counters
  _iRun = 0;
  _iEvt = 0;

  _sensorIDVec = geo::gGeometry().sensorIDsVec();
  _sensorIDtoZOrderMap.clear();
  for(size_t index = 0; index < _sensorIDVec.size(); index++) {
  	int sensorID = _sensorIDVec.at(index);
    _sensorIDtoZOrderMap.insert(std::make_pair(sensorID, static_cast<int>(index)));
  
    if(sensorID != _fixedID) {
      _preAligners.push_back(PreAligner(geo::gGeometry().getPlaneXPitch(sensorID)/10.,
					geo::gGeometry().getPlaneYPitch(sensorID)/10.,
					geo::gGeometry().getPlaneZPosition(sensorID),
					sensorID));
    }	
  }	
}

void EUTelPreAligner::processRunHeader(LCRunHeader *rdr) {

  std::unique_ptr<EUTelRunHeaderImpl> runHeader =
    std::make_unique<EUTelRunHeaderImpl>(rdr);
  runHeader->addProcessor(type());
  ++_iRun;
}

void EUTelPreAligner::processEvent(LCEvent *event) {
  
  if(_histogramSwitch && isFirstEvent()) {
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    bookHistos();
#endif
  }

  ++_iEvt;

  //if number of required events reached, stop
  if(_iEvt > _requiredEvents)
    return;

  EUTelEventImpl *evt = static_cast<EUTelEventImpl *>(event);

  //check event type
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

  try {
    LCCollectionVec *inputCollectionVec = dynamic_cast<LCCollectionVec *>(
        evt->getCollection(_inputHitCollectionName));
    UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder(EUTELESCOPE::HITENCODING);

    std::vector<float> residX;
    std::vector<float> residY;
    std::vector<PreAligner *> prealign;

    //[START] loop over hits in fixed plane
    for(size_t ref = 0; ref < inputCollectionVec->size(); ref++) {

      TrackerHitImpl *refHit =
          dynamic_cast<TrackerHitImpl *>(inputCollectionVec->getElementAt(ref));
      const double *refPos = refHit->getPosition();
      int sensorID = hitDecoder(refHit)["sensorID"];
      
      //identify fixed plane
      if(sensorID != _fixedID)
        continue;

      residX.clear();
      residY.clear();
      prealign.clear();

	  //[START] loop over other hits
      for(size_t iHit = 0; iHit < inputCollectionVec->size(); iHit++) {

        TrackerHitImpl *hit = dynamic_cast<TrackerHitImpl *>(
            inputCollectionVec->getElementAt(iHit));
        const double *pos = hit->getPosition();
        int iHitID = hitDecoder(hit)["sensorID"];

		//if fixed plane, skip
        if(iHitID == _fixedID)
          continue;
          
        bool gotIt(false);

		//[START] loop over preAligners
        for(size_t ii = 0; ii < _preAligners.size(); ii++) {

          PreAligner &pa = _preAligners.at(ii);

          if(pa.getIden() != iHitID)
          	continue;
   
   	      gotIt = true;

          double correlationX = refPos[0] - pos[0];
          double correlationY = refPos[1] - pos[1];
          int idZ = _sensorIDtoZOrderMap[iHitID];

          if((_residualsXMin[idZ] < correlationX) &&
              (correlationX < _residualsXMax[idZ]) &&
              (_residualsYMin[idZ] < correlationY) &&
              (correlationY < _residualsYMax[idZ])) {
            residX.push_back(correlationX);
            residY.push_back(correlationY);
            prealign.push_back(&pa);
          }
          break;
        }//[END] loop over prealigners
        
        if(!gotIt) {
          streamlog_out(ERROR5) << "Mismatched hit at " << pos[2] << endl;
        }
      }//[END] loop over other hits

      if(prealign.size() > static_cast<unsigned int>(_minNumberOfCorrelatedHits) 
      		&& residX.size() == residY.size()) {
      		
      	//[START] loop over prealigners
        for(unsigned int ii = 0; ii < prealign.size(); ii++) {

          prealign.at(ii)->addPoint(residX.at(ii), residY.at(ii));

		#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
          if(_histogramSwitch) {
            (dynamic_cast<AIDA::IHistogram1D *>(
                 _hitXCorr[prealign.at(ii)->getIden()]))
                ->fill(residX.at(ii));
            (dynamic_cast<AIDA::IHistogram1D *>(
                 _hitYCorr[prealign.at(ii)->getIden()]))
                ->fill(residY.at(ii));
          }
		#endif
        }//[END] loop over prealigners
     }
     
    }//[END] loop over hits in fixed plane
  } catch (DataNotAvailableException &e) {
    streamlog_out(WARNING2) << "No input collection " << _inputHitCollectionName
                            << " found on event " << event->getEventNumber()
                            << " in run " << event->getRunNumber() << std::endl;
  }

  if(isFirstEvent())
    _isFirstEvent = false;
}

void EUTelPreAligner::end() {

  LCCollectionVec *constantsCollection =
      new LCCollectionVec(LCIO::LCGENERICOBJECT);

  //[START] loop over sensorID
  for(size_t ii = 0; ii < _sensorIDVec.size(); ii++) {
    bool ifound = false;
    
    //[START] loop over prealigners
    for(size_t jj = 0; jj < _preAligners.size(); jj++) {
      int sensorID = _preAligners.at(jj).getIden();
      if(_sensorIDVec[ii] == sensorID) {
        ifound = true;
        break;
      }
    }//[END] loop over prealigners
    
    if(ifound == false) {
      EUTelAlignmentConstant *constant = new EUTelAlignmentConstant();
      constant->setXOffset(0.0);
      constant->setYOffset(0.0);
      constant->setSensorID(_sensorIDVec[ii]);
      constantsCollection->push_back(constant);
      streamlog_out(MESSAGE5) << (*constant) << endl;
      continue;
    }
  }//[END] loop over sensorID

  //[START] loop over prealigners
  for(size_t ii = 0; ii < _preAligners.size(); ii++) {
    int sensorID = _preAligners.at(ii).getIden();
    std::vector<int>::iterator it =
        find(_ExcludedPlanes.begin(), _ExcludedPlanes.end(), sensorID);
    std::vector<int>::iterator itXCoord = find(
        _ExcludedPlanesXCoord.begin(), _ExcludedPlanesXCoord.end(), sensorID);
    std::vector<int>::iterator itYCoord = find(
        _ExcludedPlanesYCoord.begin(), _ExcludedPlanesYCoord.end(), sensorID);

    EUTelAlignmentConstant *constant = new EUTelAlignmentConstant();
    if(it == _ExcludedPlanes.end()) {
      if(itXCoord == _ExcludedPlanesXCoord.end() &&
          abs(_preAligners.at(ii).getPeakX()) < 1000) {
        constant->setXOffset(-1.0 * _preAligners.at(ii).getPeakX());
      } else {
        constant->setXOffset(0.0);
	  }
      if(itYCoord == _ExcludedPlanesYCoord.end() &&
          abs(_preAligners.at(ii).getPeakY()) < 1000.) {
        constant->setYOffset(-1.0 * _preAligners.at(ii).getPeakY());
      } else {
        constant->setYOffset(0.0);
      }
    } else {
      constant->setXOffset(0.0);
      constant->setYOffset(0.0);
    }
    constant->setSensorID(sensorID);
    constantsCollection->push_back(constant);

    //update the EUTelGeometry description
     double updatedXOff;
     double updatedYOff;
     //Only want to update if needed, otherwise get a pointless warning
     if(it == _ExcludedPlanes.end()) {
       updatedXOff = geo::gGeometry().getPlaneXPosition(sensorID) + _preAligners.at(ii).getPeakX();
       updatedYOff = geo::gGeometry().getPlaneYPosition(sensorID) + _preAligners.at(ii).getPeakY();
     } else {
        updatedXOff = 0;
        updatedYOff = 0;
     }

    double oldZPos = geo::gGeometry().getPlaneZPosition(sensorID);

    geo::gGeometry().alignGlobalPos(sensorID, updatedXOff, updatedYOff, oldZPos);

    streamlog_out(MESSAGE5) << (*constant) << std::endl;
  }//[END] loop over prealigners

  //if we don't dump into the new gear file, write out the old database
  if(!_dumpGEAR) {
    LCWriter *lcWriter = LCFactory::getInstance()->createLCWriter();
    try {
      lcWriter->open(_alignmentConstantLCIOFile, LCIO::WRITE_NEW);
    } catch(IOException &e) {
      streamlog_out(ERROR4) << e.what() << std::endl;
      exit(-1);
    }
    streamlog_out(MESSAGE5) << "Writing to " << _alignmentConstantLCIOFile
                            << std::endl;

    LCRunHeaderImpl *lcHeader = new LCRunHeaderImpl;
    lcHeader->setRunNumber(0);
    lcWriter->writeRunHeader(lcHeader);
    delete lcHeader;
    LCEventImpl *event = new LCEventImpl;
    event->setRunNumber(0);
    event->setEventNumber(0);
    LCTime *now = new LCTime;
    event->setTimeStamp(now->timeStamp());
    delete now;

    streamlog_out(DEBUG5) << " adding collection alignment " << std::endl;
    event->addCollection(constantsCollection, "alignment");
    lcWriter->writeEvent(event);
    delete event;
    lcWriter->close();
  } 
  //write updated GEAR file
  else {
    marlin::StringParameters *MarlinStringParams = marlin::Global::parameters;
    std::string outputFilename =
        (MarlinStringParams->getStringVal("GearXMLFile"))
            .substr(0,
                    (MarlinStringParams->getStringVal("GearXMLFile")).size() -
                        4);
    streamlog_out(MESSAGE5) << "Writing updated GEAR file with filename: "
                            << outputFilename + "_pre.xml" << std::endl;
    geo::gGeometry().writeGEARFile(outputFilename + _GEARFileSuffix + ".xml");

    //if not write out collection to LCIO, need to delete ourselves
    delete constantsCollection;
  }
}

void EUTelPreAligner::bookHistos() {

  streamlog_out(DEBUG5) << "Booking histograms " << std::endl;
   
  //allow any plane to be the fixed reference
    for(size_t i = 0; i < _sensorIDVec.size(); i++) {
      int sensorID = _sensorIDVec.at(i);

	//create folder for current detector
      std::string basePath = "detector_" + to_string(sensorID);
      marlin::AIDAProcessor::tree(this)->mkdir(basePath.c_str());
      basePath.append("/");

	//create 1D histogram: hit correlation X
      std::string histName_hitXCorr = "hitXCorr_fixed_to_" + to_string(sensorID);
      AIDA::IHistogram1D *hist1D_hitXCorr =
          marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(
              (basePath + histName_hitXCorr).c_str(), 100, -10., 10.);
      _hitXCorr.insert(make_pair(sensorID, hist1D_hitXCorr));

	//create 1D histogram: hit correlation Y
      std::string histName_hitYCorr = "hitYCorr_fixed_to_" + to_string(sensorID);
      AIDA::IHistogram1D *hist1D_hitYCorr =
          marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(
              (basePath + histName_hitYCorr).c_str(), 100, -10., 10.);
      _hitYCorr.insert(make_pair(sensorID, hist1D_hitYCorr));
    }
   
  streamlog_out(DEBUG5) << "end of booking histograms " << std::endl;
}
