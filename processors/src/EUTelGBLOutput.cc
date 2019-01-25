/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelGBLOutput.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelTrackerDataInterfacerImpl.h"
#include "EUTelGenericPixGeoDescr.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

// lcio includes <.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <UTIL/CellIDDecoder.h>
#include "IMPL/LCGenericObjectImpl.h"

// system includes <>
#include <algorithm>

using namespace eutelescope;

EUTelGBLOutput::EUTelGBLOutput() : Processor("EUTelGBLOutput") {

  _description = "EUTelGBLOutput generates NTuple containing hits, zero suppressed data "
                 "and track position from GBL processor";
  
  registerInputCollections(LCIO::TRACKERHIT,
			   "InputHitCollections",
			   "Vector with the names of the hit collections (use the right coordinate system)",
			   _inputHitCollections,
			   std::vector<std::string>{});
  
  registerInputCollections(LCIO::TRACKERHIT,
			   "InputZsCollections",
			   "Vector with the names of the zs collections",
			   _inputZsCollections,
			   std::vector<std::string>{});
  
  registerProcessorParameter("onlyEventsWithTracks",
			     "Decide whether to dump hits and zs data only for events with a track",
			     _onlyWithTracks,
			     true);
  
  registerProcessorParameter("tracksLocalSystem",
			     "Decide whether to dump track positions in local coordinate system",
			     _tracksLocalSystem,
			     true);
  
  registerProcessorParameter("dumpHeader",
			     "Decide whether to dump the event header information",
			     _dumpHeader,
			     false);
  
  registerProcessorParameter("OutputPath",
			     "Path/File where root-file should be stored",
			     _path2file,
			     std::string("NTuple.root"));
  
  registerProcessorParameter("OutputPlanes",
			     "IDs for which the information should be dumped (leave empty for all the planes)",
			     _SelectedPlanes,
			     std::vector<int>());
}

void EUTelGBLOutput::init() {
  
  //usually a good idea to do
  printParameters();

  //reset run and event counters
  _nRun = 0;
  _nEvt = 0;

  //prepare TTree  
  _file = new TFile(_path2file.c_str(), "RECREATE");

  _planeID = new std::vector<int>();
  _trackID = new std::vector<int>();
  _xPos = new std::vector<double>();
  _yPos = new std::vector<double>();
  _omega = new std::vector<double>();
  _phi = new std::vector<double>();
  _kinkx = new std::vector<double>();
  _kinky = new std::vector<double>();
  _chi2 = new std::vector<double>();
  _ndof = new std::vector<int>();
 
  _hitXPos = new std::vector<double>();
  _hitYPos = new std::vector<double>();
  _hitZPos = new std::vector<double>();
  _hitSensorId = new std::vector<int>();

  zs_id = new std::vector<int>();
  zs_x = new std::vector<int>();
  zs_y = new std::vector<int>();
  zs_signal = new std::vector<double>();
  zs_time = new std::vector<int>();

  //tree for version number used for TBmon2
  _versionTree = new TTree("version", "version");
  _versionNo = new std::vector<double>();
  _versionTree->Branch("no", &_versionNo);

  //tree for storing track information
  _eutracks = new TTree("Tracks", "Tracks");
  _eutracks->SetAutoSave(1000000000);
  _eutracks->Branch("nTrackParams", &_nTrackParams);
  _eutracks->Branch("eventNumber", &_nEvt);
  _eutracks->Branch("planeID", &_planeID);
  _eutracks->Branch("trackID", &_trackID);
  _eutracks->Branch("triggerID",&_triggerID);
  _eutracks->Branch("timestamp",&_timestamp);
  _eutracks->Branch("xPos", &_xPos);
  _eutracks->Branch("yPos", &_yPos);
  _eutracks->Branch("omega", &_omega);
  _eutracks->Branch("phi", &_phi);
  _eutracks->Branch("kinkx", &_kinkx);
  _eutracks->Branch("kinky", &_kinky);
  _eutracks->Branch("chi2", &_chi2);
  _eutracks->Branch("ndof", &_ndof);
  
  if(_inputHitCollections.size() != 0) {
    //tree for storing hit information
    _euhits = new TTree("Hits", "Hits");
    _euhits->SetAutoSave(1000000000);
    _euhits->Branch("nHits", &_nHits);
    _euhits->Branch("eventNumber", &_nEvt);
    _euhits->Branch("ID", &_hitSensorId);
    _euhits->Branch("xPos", &_hitXPos);
    _euhits->Branch("yPos", &_hitYPos);
    _euhits->Branch("zPos", &_hitZPos);
  }
  
  if(_inputZsCollections.size() != 0) {
  //tree for storing zero suppressed data 
    _zstree = new TTree("ZeroSuppressed", "ZeroSuppressed");
    _zstree->SetAutoSave(1000000000);
    _zstree->Branch("nPixHits", &_nPixHits);
    _zstree->Branch("eventNumber", &_nEvt);
    _zstree->Branch("ID", &zs_id);
    _zstree->Branch("xPos", &zs_x);
    _zstree->Branch("yPos", &zs_y);
    _zstree->Branch("Signal", &zs_signal);
    _zstree->Branch("Time", &zs_time);
  }

  _euhits->AddFriend(_zstree);
  _euhits->AddFriend(_eutracks);

  //initialize geometry
  geo::gGeometry().initializeTGeoDescription(EUTELESCOPE::GEOFILENAME,
                                             EUTELESCOPE::DUMPGEOROOT);

  for (auto dutID : _SelectedPlanes) {
    // Later we need to shift the sensor since in EUTel centre of sensor is 0|0
    // while in TBmon(II) it is in the lower left corner
    geo::EUTelGenericPixGeoDescr *geoDescr =
        geo::gGeometry().getPixGeoDescr(dutID);
    double xSize, ySize;
    geoDescr->getSensitiveSize(xSize, ySize);

    _xSensSize[dutID] = xSize;
    _ySensSize[dutID] = ySize;
	}

}

void EUTelGBLOutput::processRunHeader(LCRunHeader *runHeader) {

  auto eutelHeader = std::make_unique<EUTelRunHeaderImpl>(runHeader);
  eutelHeader->addProcessor(type());
  //increment run counter
  _nRun++;
}

void EUTelGBLOutput::processEvent(LCEvent *event) {
  
  //increment event counter
  _nEvt++;
  
  EUTelEventImpl *euEvent = static_cast<EUTelEventImpl *>(event);
  _evtNr = event->getEventNumber();
  
  //check event type
  if(euEvent->getEventType() == kEORE) {
    streamlog_out(DEBUG5) << "EORE found: nothing else to do." << std::endl;
    return;
  }
 
  //clear all event info containers
  clear();

  //fill track TTree (hardcoded name for track collection!)
  LCCollection* TrackCollection = event->getCollection("TracksCollection");

  _triggerID = event->getParameters().getIntVal("TriggerNumber");
  //FIXME: This is disgusting...
  _timestamp = event->getTimeStamp()%((long64)INT_MAX);
  
  int nTrackParams=0;
  
  //[START] loop over track collection
  for(int itrack = 0; itrack < TrackCollection->getNumberOfElements(); itrack++) {
  
	EVENT::LCGenericObject* trackposition = 
	  static_cast<EVENT::LCGenericObject*>(TrackCollection->getElementAt(itrack));
	int thisID = trackposition->getIntVal(0);
	
	if(_SelectedPlanes.size() == 0 || std::find(std::begin(_SelectedPlanes), std::end(_SelectedPlanes),
						    thisID) != _SelectedPlanes.end()) {
								
	
      _planeID->push_back(thisID);
      _trackID->push_back(trackposition->getIntVal(2));  
      _ndof->push_back(trackposition->getIntVal(1));
      //FIXME: inserting float numbers into a double, since root doesn't want vector of floats
      _chi2->push_back(trackposition->getFloatVal(0));
      
      //[IF] local coordinates
      if(_tracksLocalSystem) {
        double pos[3];
        pos[0] = trackposition->getFloatVal(1);
        pos[1] = trackposition->getFloatVal(2);
        pos[2] = trackposition->getFloatVal(3);
        double pos_loc[3];
        geo::gGeometry().master2Local(thisID, pos, pos_loc);
        _xPos->push_back(pos_loc[0] + _xSensSize.at(thisID) / 2.0); 
        _yPos->push_back(pos_loc[1] + _ySensSize.at(thisID) / 2.0);
      } else {
        _xPos->push_back(trackposition->getFloatVal(1)); 
        _yPos->push_back(trackposition->getFloatVal(2));
      }//[ENDIF]
      
      _omega->push_back(trackposition->getFloatVal(4));
      _phi->push_back(trackposition->getFloatVal(5));
      //FIXME: What happens for the first plane which doesn't have well defined kink angles?
      _kinkx->push_back(trackposition->getFloatVal(6));
      _kinky->push_back(trackposition->getFloatVal(7));
      
      nTrackParams++; 
    }   
  }//[END] loop over track collection

  _nTrackParams = nTrackParams;
  
  //[IF] no tracks needed
  if(!(_onlyWithTracks) || TrackCollection->getNumberOfElements() != 0) {
    
    //[START] loop over hit collections
    for(auto hitName : _inputHitCollections) {
      LCCollection *hitCollection = event->getCollection(hitName);

      int nHit = hitCollection->getNumberOfElements();
      _nHits = nHit;

  	  //[START] loop over hits
      for(int ihit = 0; ihit < hitCollection->getNumberOfElements(); ihit++) {
        TrackerHitImpl *meshit = dynamic_cast<TrackerHitImpl *>(hitCollection->getElementAt(ihit));
        const double *pos = meshit->getPosition();
        UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder(EUTELESCOPE::HITENCODING);
        int thisID = hitDecoder(meshit)["sensorID"];
        
        if(_SelectedPlanes.size() == 0 || std::find(std::begin(_SelectedPlanes), std::end(_SelectedPlanes),
						    thisID) != _SelectedPlanes.end()) {

          double x = pos[0];
          double y = pos[1];
          double z = pos[2];

          _hitSensorId->push_back(thisID);   
          _hitXPos->push_back(x + _xSensSize.at(thisID) / 2.0);
          _hitYPos->push_back(y + _ySensSize.at(thisID) / 2.0);
          _hitZPos->push_back(z);

        } 
      }//[END] loop over hits
    }//[END] loop over hit collections

    //[START] loop over zs collections
    for(auto zsName : _inputZsCollections) {
      LCCollectionVec *zsInputCollectionVec = nullptr;
      try {
        zsInputCollectionVec = dynamic_cast<LCCollectionVec *>(event->getCollection(zsName));
      } catch(DataNotAvailableException &e) {
        streamlog_out(DEBUG2) << "Raw ZS data collection " << zsName
                              << " not found in event " << event->getEventNumber()
                              << "!" << std::endl;
      }

      UTIL::CellIDDecoder<TrackerDataImpl> cellDecoder(zsInputCollectionVec);
      //[START] loop over planes
      for(unsigned int plane = 0; plane < zsInputCollectionVec->size(); plane++) {
        TrackerDataImpl *zsData = dynamic_cast<TrackerDataImpl *>
	  (zsInputCollectionVec->getElementAt(plane));
        SparsePixelType type = static_cast<SparsePixelType>
	  (static_cast<int>(cellDecoder(zsData)["sparsePixelType"]));
        int thisID = cellDecoder(zsData)["sensorID"];
        
        if(type == kEUTelGenericSparsePixel) { 
          if(_SelectedPlanes.size() == 0 || std::find(std::begin(_SelectedPlanes), std::end(_SelectedPlanes),
						      thisID) != _SelectedPlanes.end()) {
            auto sparseData = std::make_unique<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>>(zsData);
            //[START] loop over pixel
            for(auto &thispixel : *sparseData) {
              _nPixHits++;
              zs_id->push_back(thisID);
              zs_x->push_back(thispixel.getXCoord());
              zs_y->push_back(thispixel.getYCoord());
              zs_signal->push_back(static_cast<double>(thispixel.getSignal()));
              zs_time->push_back(static_cast<int>(thispixel.getTime()));
            }//[END] loop over pixel
          }
        } else {
          throw UnknownDataTypeException("Unknown sparsified pixel");
        }
      }//[END] loop over planes  
    }//[END] loop over zs collections
  }//[ENDIF] no tracks needed
  
  //fill event header TTree
  if(_dumpHeader) {
    //FIXME: include here
  }
  
  //fill the TTrees
  //the event number would make it fill this TTree even for events with no tracks
  if(!(_onlyWithTracks) || TrackCollection->getNumberOfElements() != 0) _eutracks->Fill();
  if(_inputHitCollections.size() != 0) _euhits->Fill();
  if(_inputZsCollections.size() != 0) _zstree->Fill();
}

void EUTelGBLOutput::end() {
  //Write version number for TBmon2
  _versionNo->push_back(2.0);
  _versionTree->Fill();
  _file->Write();
}

void EUTelGBLOutput::clear() {
  //clear hittrack
  _planeID->clear();
  _trackID->clear();
  _xPos->clear();
  _yPos->clear();
  _omega->clear();
  _phi->clear();
  _kinkx->clear();
  _kinky->clear();
  _chi2->clear();
  _ndof->clear();
  //clear hits
  _hitSensorId->clear();
  _hitXPos->clear();
  _hitYPos->clear();
  _hitZPos->clear();
  //clear zsdata
  zs_id->clear();
  zs_x->clear();
  zs_y->clear();
  zs_signal->clear();
  zs_time->clear();
  _nPixHits = 0;
}
