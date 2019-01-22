/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelGBL.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"
#include "EUTelPStream.h" // process streams redi::ipstream
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelGenericPixGeoDescr.h"
#include "EUTelTripletGBLUtility.h"

// GBL includes
#include "include/GblTrajectory.h"
#include "include/MilleBinary.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"
#include "marlin/Exceptions.h"
#include "marlin/AIDAProcessor.h"
#include "Mille.h"

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <IO/LCWriter.h>
#include <UTIL/LCTime.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/LCCollectionVec.h>
#include "IMPL/LCGenericObjectImpl.h"

// system includes <>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>

// ROOT includes ".h"
#include <TMath.h>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;

EUTelGBL::EUTelGBL(): Processor("EUTelGBL") {

  _description = "EUTelGBL uses the MILLE program to write data files for MILLEPEDE II.";


  registerInputCollections(LCIO::TRACKERHIT,
			   "hitCollectionName",
			   "Input hit collections name",
			   _hitCollectionName,
			   std::vector<std::string>{"hit"});

  registerProcessorParameter("eBeam",
			     "Beam energy [GeV]",
			     _eBeam,
			     4.0);

  registerProcessorParameter("kappa",
			     "Global factor to Highland formula, 1.0 means HL as is, 1.2 means 20/% additional scattering",
			     _kappa,
			     1.0);

  registerOptionalParameter("excludedPlanes",
			    "Exclude planes from fit (their scattering budget is considered)",
			    _excludedPlanes,
			    std::vector<int>());

  registerOptionalParameter("requiredPlane",
			    "Only tracks with a hit on selected plane are considered",
			    _requiredPlane,
			    -1);

  registerOptionalParameter("upstreamTriplet",
			    "The three sensors used as the upstream triplet",
			    _upstreamTriplet_IDs,
			    std::vector<int>{0,1,2});

  registerOptionalParameter("downstreamTriplet",
			    "The three sensors used as the downstream triplet",
			    _downstreamTriplet_IDs,
			    std::vector<int>{3,4,5});

  registerOptionalParameter("lastUpstreamSensor",
			    "The last plane (z-ordered) which still should be attached to the upstream triplet",
			    _lastUpstreamSensorID,
			    2);

  registerOptionalParameter("resolutionX",
			    "x-resolution of sensors (z-ordered) [mm]",
			    _xResolutionVec,
			    std::vector<float>());

  registerOptionalParameter("resolutionY",
			    "y-resolution of sensors (z-ordered) [mm]",
			    _yResolutionVec,
			    std::vector<float>());

  registerOptionalParameter("fixedPlanes",
			    "Fix sensor planes in the fit according to their sensor ids (it is recommended to fix two telescope planes)",
			    _fixedPlanes,
			    std::vector<int>{0,5});

  registerOptionalParameter("maxTrackCandidatesTotal",
			    "Maximal number of track candidates (Total)",
			    _maxTrackCandidatesTotal,
			    10000000);

  registerOptionalParameter("milleBinaryFilename",
			    "Name of the Millepede binary file",
			    _binaryFilename,
			    std::string{"mille.bin"});

  registerOptionalParameter("alignMode","Number of alignment constants used. Available mode are:"
			    "\n\t\tXYShifts - shifts in X and Y"
			    "\n\t\tXYShiftsRotZ - shifts in X and Y and rotation around the Z axis,"
			    "\n\t\tXYZShiftsRotZ - shifts in X,Y and Z and rotation around the Z axis"
			    "\n\t\tXYZShiftsRotXYZ - all shifts and rotations allowed",
			    _alignModeString,
			    std::string{ "XYShiftsRotZ" });

  registerOptionalParameter("performAlignment",
			    "Set to 0 if you do not want to perform the alignment",
			    _performAlignment,
			    1);

  registerOptionalParameter("suggestAlignmentCuts",
			    "Set to 1 if you want suggestions of which cuts to use on the tracks "
			    "slope/residual to be print out during the alignment",
			    _suggestAlignmentCuts,
			    0);

  registerOptionalParameter("dumpTracks",
			    "Set to 0 if you do not want to dump tracks in an lcio collection "
			    "(necessary to dump in an NTuple)",
			    _dumpTracks,
			    1);

  registerOptionalParameter("fixedXShift",
			    "List of planes which should be fixed in X direction",
			    _FixedXShift,
			    std::vector<int>());

  registerOptionalParameter("fixedYShift",
			    "List of planes which should be fixed in Y direction",
			    _FixedYShift,
			    std::vector<int>());

  registerOptionalParameter("fixedZShift",
			    "List of planes which should be fixed in Z direction",
			    _FixedZShift,
			    std::vector<int>());

  registerOptionalParameter("fixedZRot",
			    "List of planes which should have a fixed Z rotation",
			    _FixedZRot,
			    std::vector<int>());

  registerOptionalParameter("fixedXRot",
			    "List of planes which should have a fixed X rotation",
			    _FixedXRot,
			    std::vector<int>{0,1,2,3,4,5});

  registerOptionalParameter("fixedYRot",
			    "List of planes which should have a fixed Y rotation",
			    _FixedYRot,
			    std::vector<int>{0,1,2,3,4,5});

  registerOptionalParameter("upstreamTripletResidualCut",
			    "Upstream triplet residual cut [mm]",
			    _upstreamTriplet_ResCut,
			    0.30);

  registerOptionalParameter("downstreamTripletResidualCut",
			    "Downstream triplet residual cut [mm]",
			    _downstreamTriplet_ResCut,
			    0.40);

  registerOptionalParameter("upstreamTripletSlopeCut",
			    "Upstream triplet slope cut [mrad]",
			    _upstreamTriplet_SlopeCut,
			    3.);

  registerOptionalParameter("downstreamTripletSlopeCut",
			    "Downstream triplet slope cut [mrad]",
			    _downstreamTriplet_SlopeCut,
			    5.);

  registerOptionalParameter("tripletsMatchingCut",
			    "Upstream-downstream triplet matching cut [mm]",
			    _upDownTripletMatchCut,
			    0.60);

  registerOptionalParameter("SUT_ID",
			    "ID of the plane whose scattering should be investigated (negative for none)",
			    _SUT_ID,
			    -1);

  registerOptionalParameter("DUTCuts",
			    "Cuts in x and y for matching DUT hits [mm]",
			    _dutCuts,
			    std::vector<float>{1.,1.});

  registerOptionalParameter("zMid",
			    "Z Position used for triplet matching (default to center between end of "
			    "first triplet and start of second triplet)",
			    _zMid,
			    -1.0);

  registerOptionalParameter("chi2Cut",
			    "Cut on chi2 over Ndf for the tracks to be passed to Millepede",
			    _chi2Cut,
			    100.);

  registerOptionalParameter("pedeSteerfileName",
			    "Name of the steering file for the pede program",
			    _pedeSteerfileName,
			    std::string{"steer_mille.txt"});
}


void EUTelGBL::init() {
  geo::gGeometry().initializeTGeoDescription(EUTELESCOPE::GEOFILENAME,
                                             EUTELESCOPE::DUMPGEOROOT);
  
  //get vector of sensorIDs ordered along the global z-axis (guaranteed by the framework)
  _sensorIDVec = geo::gGeometry().sensorIDsVec();
  _nPlanes = _sensorIDVec.size();

  bool isStillUpstream = true;
  for(auto& sensorID: _sensorIDVec) {
    _planePosition.emplace_back( geo::gGeometry().getPlaneZPosition(sensorID) );
    auto dz = geo::gGeometry().getPlaneZSize(sensorID);
    auto rad = geo::gGeometry().getPlaneRadiationLength(sensorID);
    _planeRadLength.emplace_back(dz/rad);
    
    //Since we go along z-ordering, we can decide if a sensor still is associated to the up or downstream triplet
    _isSensorUpstream[sensorID] = isStillUpstream;
    if(sensorID == _lastUpstreamSensorID) {
       isStillUpstream = false;
    }

    auto belongsToUp = (std::find(std::begin(_upstreamTriplet_IDs), std::end(_upstreamTriplet_IDs), sensorID) != _upstreamTriplet_IDs.end());
    auto belongsToDown = (std::find(std::begin(_downstreamTriplet_IDs), std::end(_downstreamTriplet_IDs), sensorID) != _downstreamTriplet_IDs.end());

    if(!belongsToUp && !belongsToDown){
      _DUT_IDs.emplace_back(sensorID);
    }
  }
  
  //summary output of configuration (IDs)
  streamlog_out(MESSAGE4)   << "Upstream sensors have following IDs : ";
  for(int id : _upstreamTriplet_IDs) streamlog_out(MESSAGE4) << id << " ";
  streamlog_out(MESSAGE4)   << "\n and downstream sensors are : ";
  for(int id : _downstreamTriplet_IDs) streamlog_out(MESSAGE4) << id << " ";
  streamlog_out(MESSAGE4)   << "\n and DUTs are sensors with following ID : ";
  for(int id : _DUT_IDs) streamlog_out(MESSAGE4) << id << " ";
  streamlog_out(MESSAGE4) << std::endl;

  //compute the total radiation length
  double totalRadLength = 0;
  //add air from the first to the last plane
  const double radLengthAir = 304200.; //X0 for air is 30420cm
  totalRadLength += (_planePosition.back()-_planePosition.front())/radLengthAir;
  //loop over all planes
  for(auto& radLen: _planeRadLength) {
    totalRadLength += radLen; //FIXME: SHOULD WE EXCLUDE THE SUT FROM THIS CALCULATION?
  }

  for(auto& radLen: _planeRadLength) {
      //Paper showed HL predicts too high angle, at least for biased measurement. 
      double tetSi = _kappa*0.0136 * sqrt(radLen) / _eBeam * ( 1 + 0.038*std::log(totalRadLength) );
      _planeWscatSi.emplace_back( 1.0/(tetSi*tetSi), 1.0/(tetSi*tetSi) );
  }

  for(size_t ipl = 0; ipl < _nPlanes-1; ipl++) {
    auto distplane = _planePosition[ipl+1] - _planePosition[ipl]; // in [mm]
    double epsAir =   0.5*distplane  / radLengthAir;
    double tetAir = _kappa*0.0136 * sqrt(epsAir) / _eBeam * ( 1 + 0.038*std::log(totalRadLength) );
    _planeWscatAir.emplace_back( 1.0/(tetAir*tetAir), 1.0/(tetAir*tetAir) );
  }

  //FIXME: Really needed this output?
  streamlog_out(MESSAGE4) << "Assumed beam energy " << _eBeam << " GeV" <<  std::endl;
  streamlog_out(MESSAGE4) << "\n\t_planeRadLength has size: \t" << _planeRadLength.size() 
			  << "\n\t_planeWscatSi has size: \t" << _planeWscatSi.size() 
			  << "\n\t_planeWscatAir has size: \t" << _planeWscatAir.size() << std::endl;

  for(size_t ipl = 0; ipl < _nPlanes; ipl++) {
    auto x_res = _xResolutionVec.at(ipl);
    auto y_res = _yResolutionVec.at(ipl);
    _planeMeasPrec.emplace_back(1.0/(x_res*x_res), 1.0/(y_res*y_res));
  }

  //calculation of center of triplets (only needed if not given by user)
  if( _zMid<0 ) {
    double upstreamTriplet_zEnd = 0;
    double downstreamTriplet_zStart = 0;
    for (size_t ipl = 0; ipl < _nPlanes; ipl++){
      if(_sensorIDVec[ipl] == _upstreamTriplet_IDs.back()) upstreamTriplet_zEnd=_planePosition[ipl];
      else if(_sensorIDVec[ipl] == _downstreamTriplet_IDs.front()) downstreamTriplet_zStart=_planePosition[ipl];
    }
    _zMid = 0.5*(upstreamTriplet_zEnd + downstreamTriplet_zStart);
  }
  
  streamlog_out( MESSAGE4 ) << "Extrapolated tracks from the triplets will be matched in z position = " << _zMid << std::endl;

  //usually a good idea to do
  printParameters ();

  //reset run and event counters
  _iRun = 0;
  _iEvt = 0;
  _printEventCounter = 0;

  //initialise Mille statistics
  _nTotalTracks = 0;

  // booking histograms
  bookHistos(_sensorIDVec);
  gblutil.setParent(this);
  gblutil.bookHistos();
  
  //only for alignment
  if(_performAlignment){ 
    streamlog_out( MESSAGE2 ) << "Initialising Mille..." << std::endl;

    unsigned int reserveSize = 8000;
    milleAlignGBL = std::make_unique<gbl::MilleBinary>( _binaryFilename, reserveSize );

    streamlog_out( MESSAGE2 ) << "Filename for the binary file is: " << _binaryFilename.c_str() << std::endl;

    if(_alignModeString.compare("XYShiftsRotZ") == 0 ) {
      _alignMode = Utility::alignMode::XYShiftsRotZ;
    } else if( _alignModeString.compare("XYShifts") == 0 ) {
      _alignMode = Utility::alignMode::XYShifts;
    } else if( _alignModeString.compare("XYZShiftsRotZ") == 0 ) {
      _alignMode = Utility::alignMode::XYZShiftsRotZ;
    } else if( _alignModeString.compare("XYZShiftsRotXYZ") == 0 ) {
      _alignMode = Utility::alignMode::XYZShiftsRotXYZ;
    } else {
      streamlog_out(ERROR) << "The chosen AlignMode: '" << _alignModeString << "' is invalid. Please correct your steering template and retry!" << std::endl;
      throw InvalidParameterException("AlignMode");
    }
  
    //Creating the steering file here in init, since doing it in the end section creates problem with opening it in EUTelPedeGEAR
    streamlog_out( MESSAGE4 ) << "Generating the steering file for the pede program..." << std::endl;

    std::ofstream steerFile;
    steerFile.open(_pedeSteerfileName.c_str());

    if( steerFile.is_open()) {
      steerFile << "Cfiles" << std::endl;
      steerFile << _binaryFilename << std::endl << std::endl;
      steerFile << "Parameter" << std::endl;

      //[START] loop over all planes
      for(size_t ipl=0; ipl<_nPlanes; ipl++) {

	//check if current plane is excluded, if so skip
        if(std::find(std::begin(_excludedPlanes), std::end(_excludedPlanes), 
		     _sensorIDVec[ipl]) != _excludedPlanes.end()) {
	  continue;
	}
	
	//check if fixed: current plane is not fixed
	if(std::find(std::begin(_fixedPlanes), std::end(_fixedPlanes),
		     _sensorIDVec[ipl]) != _fixedPlanes.end()) {
	  std::map<Utility::alignMode, int> map_alignModes;
	  map_alignModes[Utility::alignMode::XYShifts] = 3;
	  map_alignModes[Utility::alignMode::XYShiftsRotZ] = 4;
	  map_alignModes[Utility::alignMode::XYZShiftsRotZ] = 5;
	  map_alignModes[Utility::alignMode::XYZShiftsRotXYZ] = 7;

	  int stopid = map_alignModes[_alignMode];
	  for(int id = 1 ; id < stopid ; id++) steerFile << (_sensorIDVec[ipl] * 10 + id) << "  0.0 -1.0" << endl;
	} 
	//check if fixed: current plane is fixed
	else {
	  std::vector< std::vector<int> > alignModeArray =  {_FixedXShift, _FixedYShift};
          
	  if(_alignMode == Utility::alignMode::XYZShiftsRotZ) 
	    alignModeArray.push_back(_FixedZRot);
	  if(_alignMode == Utility::alignMode::XYZShiftsRotXYZ || _alignMode == Utility::alignMode::XYZShiftsRotZ) 
	    alignModeArray.push_back(_FixedZShift);
	  if(_alignMode == Utility::alignMode::XYZShiftsRotXYZ) {
	    alignModeArray.push_back(_FixedXRot);
	    alignModeArray.push_back(_FixedYRot);
	  }
	  if(_alignMode == Utility::alignMode::XYShiftsRotZ || _alignMode == Utility::alignMode::XYZShiftsRotXYZ) 
	    alignModeArray.push_back(_FixedZRot);
	  
	  if(alignModeArray.size() == 2 && _alignMode != Utility::alignMode::XYShifts ) continue;
          
	  int counter = 1;
	  for(auto current : alignModeArray) {
	    if(std::find(current.begin(), current.end(), _sensorIDVec[ipl]) == current.end()) {
	      steerFile << (_sensorIDVec[ipl] * 10 + counter) << "  0.0  0.0" << endl;
	    } else {
	      steerFile << (_sensorIDVec[ipl] * 10 + counter) << "  0.0  -1.0" << endl; 
	    }
	    counter++;
	  }
	}
	
      }//[END] loop over all planes
      
      steerFile << std::endl;
      steerFile << "! chiscut 5.0 2.5" << std::endl;
      steerFile << "outlierdownweighting 4" << std::endl;
      steerFile << "dwfractioncut 0.1" << std::endl;
      steerFile << std::endl;
      steerFile << "method inversion 10  0.1" << std::endl;
      // Use 10 OpenMP threads to process the data:
      steerFile << "threads 10 1" << std::endl;
      steerFile << std::endl;
      steerFile << "histprint" << std::endl;
      steerFile << std::endl;
      steerFile << "end" << std::endl;
      steerFile.close();
      
      streamlog_out( MESSAGE2 ) << "File " << _pedeSteerfileName << " written." << std::endl;
    } else {
      streamlog_out( ERROR2 ) << "Could not open steering file." << std::endl;
    }
    // end writing the pede steering file
  }
  streamlog_out( MESSAGE2 ) << "end of init" << std::endl;
}

void EUTelGBL::processRunHeader( LCRunHeader* rdr ) {

  auto header = std::make_unique<EUTelRunHeaderImpl>(rdr);
  header->addProcessor(type());
  //increment the run counter
  ++_iRun;
}

inline Eigen::Matrix<double,5,5> Jac55new( double ds ) {
  /* for GBL:
     Jacobian for straight line track
     track = q/p, x', y', x, y
               0,  1,  2, 3, 4
     */
  Eigen::Matrix<double,5,5> jac = Eigen::Matrix<double,5,5>::Identity();
  //jac.UnitMatrix();
  jac(3,1) = ds; // x = x0 + xp * ds
  jac(4,2) = ds; // y = y0 + yp * ds
  return jac;
}

void EUTelGBL::processEvent( LCEvent * event ) {

  if(_iEvt % 1000 == 0) {
    streamlog_out( MESSAGE2 ) << "Processing event "
      << setw(6) << setiosflags(ios::right)
      << event->getEventNumber() << " in run "
      << setw(6) << setiosflags(ios::right)
      << event->getRunNumber()
      << ", currently having "
      << _nTotalTracks << " tracks "
      << std::endl;
  }
  
  //FIXME?: compiler doesn't like it inside an if clause
  LCCollectionVec* _outputTracks = new LCCollectionVec(LCIO::LCGENERICOBJECT);
  
  if(_nTotalTracks > static_cast<size_t>(_maxTrackCandidatesTotal)) {
    throw StopProcessingException(this);
  }

  EUTelEventImpl* evt = static_cast<EUTelEventImpl*> (event) ;

  if(evt->getEventType() == kEORE) {
    streamlog_out( DEBUG2 ) << "EORE found: nothing else to do." << endl;
    return;
  }

  CellIDDecoder<TrackerHit> hitCellDecoder(EUTELESCOPE::HITENCODING);
  std::vector<EUTelTripletGBLUtility::hit> telescopeHitsVec;
  std::vector<EUTelTripletGBLUtility::hit> dutHitsVec;

  //[START] loop over all input hit collections
  for(auto collectionName : _hitCollectionName) {
    LCCollection* collection = nullptr;
    try {
      collection = event->getCollection(collectionName);
    } 
    catch (DataNotAvailableException& e) {
      throw SkipEventException(this);
    }

    if(_printEventCounter < NO_PRINT_EVENT_COUNTER) {
      streamlog_out(DEBUG2) << "Event " << event->getEventNumber() << " contains " 
			    << collection->getNumberOfElements() << " hits" << endl;
    }
    
    hist1D_nTelescopeHits->fill(collection->getNumberOfElements());

    //[START] loop over all hits in collection
    for(int iHit = 0; iHit < collection->getNumberOfElements(); iHit++) {
      auto hit = static_cast<TrackerHitImpl*>( collection->getElementAt(iHit) );
      auto sensorID = hitCellDecoder(hit)["sensorID"];
      auto hitPosition = hit->getPosition();

      if(std::find(std::begin(_upstreamTriplet_IDs), std::end(_upstreamTriplet_IDs), 
		   sensorID) != _upstreamTriplet_IDs.end() || 
	 std::find(std::begin(_downstreamTriplet_IDs), std::end(_downstreamTriplet_IDs), 
		   sensorID) != _downstreamTriplet_IDs.end()) {
        telescopeHitsVec.emplace_back(hitPosition, sensorID);
      } else {
        dutHitsVec.emplace_back(hitPosition, sensorID);
      }
      if(_printEventCounter < NO_PRINT_EVENT_COUNTER) 
	streamlog_out(DEBUG0) << "Hit on plane " << sensorID << " at " 
			      << hitPosition[0] << "|" << hitPosition[1]  
			      << std::endl;
    } //[END] loop over all hits in collection
  }//[END] loop over all input hit collections

  int numbertracks = 0;
  auto upstreamTripletVec = std::vector<EUTelTripletGBLUtility::triplet>();
  auto downstreamTripletVec = std::vector<EUTelTripletGBLUtility::triplet>();

  gblutil.FindTriplets(telescopeHitsVec, _upstreamTriplet_IDs, _upstreamTriplet_ResCut,
		       _upstreamTriplet_SlopeCut/1000., upstreamTripletVec, false, true);
  gblutil.FindTriplets(telescopeHitsVec, _downstreamTriplet_IDs, _downstreamTriplet_ResCut,
		       _downstreamTriplet_SlopeCut/1000., downstreamTripletVec, false, false);

  if(_printEventCounter < NO_PRINT_EVENT_COUNTER){
    streamlog_out(DEBUG2)  << "UpstreamTriplets:" << std::endl;
    for(auto& triplet: upstreamTripletVec) {
      streamlog_out(DEBUG2) << triplet << std::endl;
    }
    streamlog_out(DEBUG2) << "DownstreamTriplets:" << std::endl;
    for(auto& triplet: downstreamTripletVec) {
      streamlog_out(DEBUG2) << triplet << std::endl;
    }
  }

  hist1D_nUpstreamTriplets->fill(upstreamTripletVec.size());
  hist1D_nDownstreamTriplets->fill(downstreamTripletVec.size());

  auto matchedTripletVec = std::vector<EUTelTripletGBLUtility::track>();
  gblutil.MatchTriplets(upstreamTripletVec, downstreamTripletVec, _zMid,
			_upDownTripletMatchCut, matchedTripletVec);

  if(_printEventCounter < NO_PRINT_EVENT_COUNTER)
    streamlog_out(DEBUG2) << "Matched to:" << std::endl; 

  if(!_DUT_IDs.empty())
    {
      for(auto& track: matchedTripletVec) {    

	/*
	auto& triplet = track.get_upstream();
	auto& driplet = track.get_downstream();
	for(auto dutID: _dut_ids) {
	  if(_is_sensor_upstream[dutID]) {
	    gblutil.AttachDUT(triplet, dutHitsVec, dutID, _dutCuts);
	  } else {
	    gblutil.AttachDUT(driplet, dutHitsVec, dutID, _dutCuts);
	  }
	  }*/

	//TO BE TESTED: to replace upper block
	for(auto dutID: _DUT_IDs) {
	  //either attach DUT to upstream
	  if(_isSensorUpstream[dutID]) {
	    gblutil.AttachDUT(track.get_upstream(), dutHitsVec, dutID, _dutCuts);
	  }
	  //or to downstream
	  else {
	    gblutil.AttachDUT(track.get_downstream(), dutHitsVec, dutID, _dutCuts);
	  }
	}
      }
    }

  //[START] loop over matched tracks
  for(auto& track: matchedTripletVec) 
    {
      //GBL point vector for the trajectory (in [mm])
      //GBL with triplet A as seed
      std::vector<gbl::GblPoint> traj_points;
      //build up trajectory:
      std::vector<unsigned int> labelVec;
      vector<double> sPoint;
      
      //arc length at the first measurement plane is 0
      double s = 0;
      
      Eigen::Matrix2d proL2m = Eigen::Matrix2d::Identity();
      Eigen::Vector2d scat = Eigen::Vector2d::Zero();
      
      auto uptriplet = track.get_upstream();
      auto downtriplet = track.get_downstream();
      //need triplet slope to compute residual
      auto tripletSlope = uptriplet.slope();
      
      std::vector<EUTelTripletGBLUtility::hit> DUThits;
      for(auto it = uptriplet.DUT_begin(); it != uptriplet.DUT_end(); ++it){
	DUThits.emplace_back(it->second);
      }
      for(auto it = downtriplet.DUT_begin(); it != downtriplet.DUT_end(); ++it){
	DUThits.emplace_back(it->second);
      }
      
      if(_printEventCounter < NO_PRINT_EVENT_COUNTER) 
	streamlog_out(DEBUG2) << "Track has " << DUThits.size() << " DUT hits" << std::endl;
      
      //selection of tracks with a hit on a selected/required plane
      if(_requiredPlane!=-1) {
        bool rejectTrack = true;
	for(auto& checkhits: DUThits) {
          if(int(checkhits.plane) == _requiredPlane) {
	    rejectTrack = false;
	    break;
	  }
	}
        if(rejectTrack) continue;
      }
      
      //FIXME: to be used only during alignment. Matrix defined outside if clause to 
      //avoid complaints from compiler. Could be done better
      
      //define alignment derivatives
      Eigen::Matrix2d alDer2;
      Eigen::Matrix<double,2,3> alDer3;
      Eigen::Matrix<double, 2,4> alDer4;
      Eigen::Matrix<double, 3,6> alDer6;
      
      if(_performAlignment){ 
	alDer2(0,0) = 1.0; // dx/dx GBL sign convention
	alDer2(1,0) = 0.0; // dy/dx
	alDer2(0,1) = 0.0; // dx/dy
	alDer2(1,1) = 1.0; // dy/dy
	
	alDer3(0,0) = 1.0; // dx/dx
	alDer3(1,0) = 0.0; // dy/dx
	alDer3(0,1) = 0.0; // dx/dy
	alDer3(1,1) = 1.0; // dy/dy
	
	alDer4(0,0) = 1.0; // dx/dx
	alDer4(1,0) = 0.0; // dy/dx
	alDer4(0,1) = 0.0; // dx/dy
	alDer4(1,1) = 1.0; // dy/dy
	alDer4(0,3) = tripletSlope.x; // dx/dz
	alDer4(1,3) = tripletSlope.y; // dx/dz
	
	alDer6(0,0) = 1.0; // dx/dx
	alDer6(0,1) = 0.0; // dx/dy
	alDer6(0,2) = tripletSlope.x; // dx/dz
	alDer6(0,3) = 0.0; // dx/da
	alDer6(1,0) = 0.0; // dy/dx
	alDer6(1,1) = 1.0; // dy/dy
	alDer6(1,2) = tripletSlope.y; // dy/dz
	alDer6(1,4) = 0.0; // dy/db
	alDer6(2,0) = 0.0; // dz/dx
	alDer6(2,1) = 0.0; // dz/dy
	alDer6(2,2) = 1.0; // dz/dz
	alDer6(2,5) = 0.0; // dz/dg
      }
      
      std::vector<double> rx (_nPlanes, -1.0);
      std::vector<double> ry (_nPlanes, -1.0);
      std::vector<bool> hasHit (_nPlanes, false);
      
      double step = 0.0;
      unsigned int iLabel;
      
      //[START] loop over all planes
      for(size_t ipl=0; ipl<_nPlanes; ++ipl) {

	//add all the planes: up/downstream telescope will have hits, DUTs maybe
	EUTelTripletGBLUtility::hit const *hit = nullptr;
	auto sensorID = _sensorIDVec[ipl];

	if(std::find(_upstreamTriplet_IDs.begin(), _upstreamTriplet_IDs.end(), 
		     sensorID) != _upstreamTriplet_IDs.end()) {
	  hit = &uptriplet.gethit(sensorID);
	} else if(std::find(_downstreamTriplet_IDs.begin(), _downstreamTriplet_IDs.end(), 
			    sensorID) != _downstreamTriplet_IDs.end()) {
	  hit = &downtriplet.gethit(sensorID);
	} else if(uptriplet.has_DUT(sensorID)) {
	  hit = &uptriplet.get_DUT_Hit(sensorID);
	} else if(downtriplet.has_DUT(sensorID)) {
	  hit = &downtriplet.get_DUT_Hit(sensorID);
	}

	//if there is no hit, take plane position from the geo description
	double zz = hit ? hit->z : _planePosition[ipl];// [mm]
	
	//transport matrix in (q/p, x', y', x, y) space
	auto jacPointToPoint = Jac55new( step );
	auto point = gbl::GblPoint( jacPointToPoint );
	s += step;
	
	if(hit) {
	  hasHit[ipl] = true; 
	  //if there is a hit, add a measurement to the point
	  //for excluded plane: want to know if there is a hit, but don't process it here
	  if(std::find(std::begin(_excludedPlanes), std::end(_excludedPlanes), 
		       _sensorIDVec[ipl]) == _excludedPlanes.end()){
	    double xs = uptriplet.getx_at(zz);
	    double ys = uptriplet.gety_at(zz);
	    
	    if(_printEventCounter < NO_PRINT_EVENT_COUNTER)
	      streamlog_out(DEBUG2) << "xs = " << xs << "   ys = " << ys << std::endl;
	    
	    //add residuals as hit to triplet
	    rx[ipl] = (hit->x - xs);
	    ry[ipl] = (hit->y - ys);
	    
	    //fill measurement vector for GBL
	    Eigen::Vector2d meas(rx[ipl], ry[ipl]);
	    point.addMeasurement( proL2m, meas, _planeMeasPrec[ipl] );
	    
	    //for SUT: add local parameter for kink estimation for the planes after the SUT
	    if(_SUT_ID > 0){
	      //FIXME: very sloppy. SUT_zpos should not be retrieved every time
	      double SUT_zpos = geo::gGeometry().getPlaneZPosition(_SUT_ID);
	      double distSUT = _planePosition[ipl] - SUT_zpos; 
	      if(distSUT > 0){
		//FIXME: definition can stay here or should be moved out?
		Eigen::Matrix<double,2,4> addDer = Eigen::Matrix<double,2,4>::Zero();
		double thickness = geo::gGeometry().getPlaneZSize(_SUT_ID);
		addDer(0,0) = (distSUT - thickness/sqrt(12)); //first scatterer in target
		addDer(1,1) = (distSUT - thickness/sqrt(12)); 
		addDer(0,2) = (distSUT + thickness/sqrt(12)); //second scatterer in target
		addDer(1,3) = (distSUT + thickness/sqrt(12)); 
		point.addLocals(addDer);
	      }
	    }

	    //only during alignment
	    if(_performAlignment) {
	      
	      //alignMode: only x and y shifts
	      if( _alignMode == Utility::alignMode::XYShifts ) {
		// global labels for MP:
		std::vector<int> globalLabels(2);
		globalLabels[0] = _sensorIDVec[ipl] * 10 + 1; //x
		globalLabels[1] = _sensorIDVec[ipl] * 10 + 2; //y
		point.addGlobals( globalLabels, alDer2 );
	      } 
	      //alignMode: x,y shifts and rotation z
	      else if( _alignMode == Utility::alignMode::XYShiftsRotZ ) {
		std::vector<int> globalLabels(3);
		globalLabels[0] = _sensorIDVec[ipl] * 10 + 1; //x
		globalLabels[1] = _sensorIDVec[ipl] * 10 + 2; //y
		globalLabels[2] = _sensorIDVec[ipl] * 10 + 3; //rotZ
		alDer3(0,2) = -ys; //dx/dphi
		alDer3(1,2) =  xs; //dy/dphi
		point.addGlobals( globalLabels, alDer3 );
	      } 
	      //alignMode: x,y,z shifts and rotation z
	      else if( _alignMode == Utility::alignMode::XYZShiftsRotZ ) {
		std::vector<int> globalLabels(4);
		globalLabels[0] = _sensorIDVec[ipl] * 10 + 1; //x
		globalLabels[1] = _sensorIDVec[ipl] * 10 + 2; //y
		globalLabels[2] = _sensorIDVec[ipl] * 10 + 3; //rotZ
		globalLabels[3] = _sensorIDVec[ipl] * 10 + 4; //z
		alDer4(0,2) = -ys; //dx/dphi
		alDer4(1,2) =  xs; //dy/dphi
		point.addGlobals( globalLabels, alDer4 );
	      } 
	      //alignMode: x,y,z shifts and rotation x,y,z
	      else if( _alignMode == Utility::alignMode::XYZShiftsRotXYZ ) {
		double deltaz = hit->z - _planePosition[ipl];
		//FIXME: a bit hacky? : deltaz cannot be zero, otherwise this mode doesn't work
		if ( deltaz < 1E-9 ) deltaz = 1E-9;
		std::vector<int> globalLabels(6);
		globalLabels[0] = _sensorIDVec[ipl] * 10 + 1; //x
		globalLabels[1] = _sensorIDVec[ipl] * 10 + 2; //y
		globalLabels[2] = _sensorIDVec[ipl] * 10 + 3; //rotZ
		globalLabels[3] = _sensorIDVec[ipl] * 10 + 4; //z
		globalLabels[4] = _sensorIDVec[ipl] * 10 + 5; //rotX
		globalLabels[5] = _sensorIDVec[ipl] * 10 + 6; //rotY
		alDer6(0,4) = deltaz; //dx/db
		alDer6(0,5) = -ys; //dx/dg
		alDer6(1,3) = -deltaz; //dy/da
		alDer6(1,5) = xs; //dy/dg
		alDer6(2,3) = ys; //dz/da
		alDer6(2,4) = -xs; //dz/db
		point.addGlobals( globalLabels, alDer6 );
	      }
	    }
	  }
	}
	
	//don't add a scatterer for the SUT in order to have an unbiased estimation of the kink
	if(_sensorIDVec[ipl] != _SUT_ID) {
	  point.addScatterer( scat, _planeWscatSi[ipl] );
	}
	sPoint.push_back( s );
	iLabel = sPoint.size();
	labelVec.push_back(iLabel);
	traj_points.push_back(point);

	//fill up with two air scatters in between planes
	if( ipl < _nPlanes-1 ) {
	  double distplane = _planePosition[ipl+1] - _planePosition[ipl];
	  step = 0.21*distplane; //in [mm]
	  auto point_left = gbl::GblPoint( Jac55new( step ) );
	  point_left.addScatterer( scat, _planeWscatAir[ipl] );
	  s += step;
	  traj_points.push_back(point_left);
	  sPoint.push_back( s );
	  step = 0.58*distplane; //in [mm]
	  auto point_right = gbl::GblPoint( Jac55new( step ) );
	  point_right.addScatterer( scat, _planeWscatAir[ipl] );
	  s += step;
	  traj_points.push_back(point_right);
	  sPoint.push_back( s );
	  step = 0.21*distplane; //remaining distance to next plane, in [mm]
	}
	
      }//[END] loop over all planes

      double Chi2;
      int Ndf;
      double lostWeight;
      
      gbl::GblTrajectory traj(traj_points, false); // curvature = false
      traj.fit( Chi2, Ndf, lostWeight);
      
      if(_printEventCounter < NO_PRINT_EVENT_COUNTER){
	streamlog_out(DEBUG4) << "traj with " << traj.getNumPoints() << " points:" << std::endl;
	for( size_t ipl = 0; ipl < sPoint.size(); ++ipl ){
	  streamlog_out(DEBUG4) << "  GBL point " << ipl;
	  streamlog_out(DEBUG4) << "  z " << sPoint[ipl]; 
	  streamlog_out(DEBUG4) << std::endl;
	}
	//FIXME: why is here size hardcoded to 6?
	for( size_t ipl = 0; ipl < 6; ++ipl ){
	  streamlog_out(DEBUG4) << " plane " << ipl << ", lab " << labelVec[ipl];
	  streamlog_out(DEBUG4) << " z " << sPoint[labelVec[ipl]-1];
	  streamlog_out(DEBUG4) << "  dx " << rx[ipl];
	  streamlog_out(DEBUG4) << "  dy " << ry[ipl];
	  streamlog_out(DEBUG4) << "  chi2 " << Chi2;
	  streamlog_out(DEBUG4) << "  ndf " << Ndf;
	  streamlog_out(DEBUG4) << std::endl;
	}
	
	streamlog_out(DEBUG2)  << " Is traj valid? " << traj.isValid() << std::endl;

	//FIXME: we should assign some verbosity to these outputs
	//traj.printPoints();
	//traj.printTrajectory();
	//traj.printData();
	_printEventCounter++;
      }
      
      if(Chi2 / Ndf >_chi2Cut) continue;
      
      hist1D_gblNdfAlign->fill( Ndf );
      hist1D_gblChi2Align->fill( Chi2 / Ndf );
      double probchi = TMath::Prob( Chi2, Ndf );
      hist1D_gblProbAlign->fill( probchi );
      
      //look at fit:
      Eigen::VectorXd localPar;
      Eigen::MatrixXd localCov;
      double prevAngleX = 0;
      double prevAngleY = 0;
      
      //at plane 0:
      unsigned int ndata = 2;
      unsigned int ndim = 2;
      Eigen::VectorXd aResiduals(ndim);
      Eigen::VectorXd aMeasErrors(ndim);
      Eigen::VectorXd aResErrors(ndim);
      Eigen::VectorXd aDownWeights(ndim);
      
      for(size_t ix = 0; ix < _nPlanes; ++ix) {
	int ipos = labelVec[ix];
	traj.getResults( ipos, localPar, localCov );
	
	//FIXME: compiler doesn't like it outside an if clause
	IMPL::LCGenericObjectImpl* thisTrack = new IMPL::LCGenericObjectImpl();
	//track = q/p, x', y', x, y
	//        0,   1,  2,  3, 4
	hist1D_gblAngleX[ix]->fill( localPar[1]*1E3 );
	hist1D_gblAngleY[ix]->fill( localPar[2]*1E3 ); 
	
	//check for an hit on plane
	if(hasHit[ix]) {
	  //check: plane is not excluded
	  if(std::find(std::begin(_excludedPlanes), std::end(_excludedPlanes),
		       _sensorIDVec[ix]) == _excludedPlanes.end()){
	    traj.getMeasResults( ipos, ndata, aResiduals, aMeasErrors,
				 aResErrors, aDownWeights );
	    hist1D_gblResidX[ix]->fill( (aResiduals[0])*1E3 );
	    hist1D_gblResidY[ix]->fill( (aResiduals[1])*1E3 );
	    hist1D_gblPullX[ix]->fill( aResiduals[0]/aResErrors[0] );
	    hist1D_gblPullY[ix]->fill( aResiduals[1]/aResErrors[1] );
	  } 
	  //check: plane is excluded
	  else {
	    EUTelTripletGBLUtility::hit const * hit = nullptr;
	    if(std::find(_upstreamTriplet_IDs.begin(), _upstreamTriplet_IDs.end(),
			 _sensorIDVec[ix]) != _upstreamTriplet_IDs.end()) {
	      hit = &uptriplet.gethit(_sensorIDVec[ix]);
	    } else if(std::find(_downstreamTriplet_IDs.begin(), _downstreamTriplet_IDs.end(),
				_sensorIDVec[ix]) != _downstreamTriplet_IDs.end()) {
	      hit = &downtriplet.gethit(_sensorIDVec[ix]);
	    } else if(uptriplet.has_DUT(_sensorIDVec[ix])) {
	      hit = &uptriplet.get_DUT_Hit(_sensorIDVec[ix]);
	    } else if(downtriplet.has_DUT(_sensorIDVec[ix])) {
	      hit = &downtriplet.get_DUT_Hit(_sensorIDVec[ix]);
	    }
	    hist1D_gblResidX[ix]->fill(hit->x*1E3 - uptriplet.getx_at(_planePosition[ix]) *1E3 - localPar[3]*1E3);
	    hist1D_gblResidY[ix]->fill(hit->y*1E3 - uptriplet.gety_at(_planePosition[ix]) *1E3 - localPar[4]*1E3);
	  }
	}
	
	if(_dumpTracks) { //CHECK ME CAREFULLY
	  thisTrack->setIntVal(0, _sensorIDVec[ix]); //sensor ID is an int
	  thisTrack->setIntVal(1, Ndf); //Ndf is an int
	  thisTrack->setIntVal(2, numbertracks);
	  thisTrack->setFloatVal(0, Chi2); //chi2
	  thisTrack->setFloatVal(1, uptriplet.getx_at(_planePosition[ix]) + localPar[3]); // x track position (global system)
	  thisTrack->setFloatVal(2, uptriplet.gety_at(_planePosition[ix]) + localPar[4]); // y track position (global system)
	  thisTrack->setFloatVal(3, _planePosition[ix]); // z of the plane, FIXME: we should use z hit position if possible
	  thisTrack->setFloatVal(4, localPar[1]*1E3); // FIXME: this is not the incidence angle, just the angle correction to the seed
	  thisTrack->setFloatVal(5, localPar[2]*1E3); // FIXME: this is not the incidence angle, just the angle correction to the seed
	  
	  if(_sensorIDVec[ix]!=_SUT_ID) {
	    thisTrack->setFloatVal(6, (localPar[1] - prevAngleX)*1E3 ); // kink angle in x in mrad
	    thisTrack->setFloatVal(7, (localPar[2] - prevAngleY)*1E3 ); // kink angle in y in mrad
	  } else {
	    thisTrack->setFloatVal(6, (localPar[5]+localPar[7])*1E3 ); 
	    thisTrack->setFloatVal(7, (localPar[6]+localPar[8])*1E3 );  
	  }
	  
	  _outputTracks->push_back(static_cast<EVENT::LCGenericObject*>(thisTrack));
	}
	
	//fill kink angle histograms [mrad]
	if(_sensorIDVec[ix]!=_SUT_ID){
	  hist1D_gblKinkX[ix]->fill( (localPar[1] - prevAngleX)*1E3 );
	  hist1D_gblKinkY[ix]->fill( (localPar[2] - prevAngleY)*1E3 );
	} else {
	  hist1D_gblKinkX[ix]->fill( (localPar[5]+localPar[7])*1E3 );
	  hist1D_gblKinkY[ix]->fill( (localPar[6]+localPar[8])*1E3 );

	  //FIXME: newly added kink maps
	  profile2D_gblSUTKinkXvsXY->fill(uptriplet.getx_at(_planePosition[ix])+localPar[3],
					  uptriplet.gety_at(_planePosition[ix])+localPar[4],
					  fabs(localPar[5]+localPar[7])*1E3 );
	  profile2D_gblSUTKinkYvsXY->fill(uptriplet.getx_at(_planePosition[ix])+localPar[3],
					  uptriplet.gety_at(_planePosition[ix])+localPar[4],
					  fabs(localPar[6]+localPar[8])*1E3 );
	}
	
	prevAngleX = localPar[1];
	prevAngleY = localPar[2];
 }
    
  // do not pass very bad tracks to mille. Only if the alignment is performed
  if(_performAlignment) {
          traj.milleOut( *milleAlignGBL );
  }
  _nTotalTracks ++;
   numbertracks++;
    }//[END] loop over matched tracks
  
  if(_dumpTracks) event->addCollection(_outputTracks,"TracksCollection");
  hist1D_nTracksPerEvent->fill( numbertracks );
  
  //count events
  _iEvt++;
}

void EUTelGBL::end() {

  milleAlignGBL.reset(nullptr);
  //if user wishes alignment cut suggestion
  if(_suggestAlignmentCuts) {
  	gblutil.determineBestCuts();
  }
  streamlog_out( MESSAGE5 ) << "Found " << _nTotalTracks << " tracks in " << _iEvt << " events" << std::endl;
  streamlog_out( MESSAGE5 ) << "Successfully finished" << std::endl;
}

void EUTelGBL::bookHistos(std::vector<int> const & sensorIDVec) {

  try {
    streamlog_out( MESSAGE2 ) << "Booking histograms..." << std::endl;

    hist1D_nTelescopeHits = AIDAProcessor::histogramFactory(this)->
      createHistogram1D("nAllHit", 201, -0.5, 200.5);
    hist1D_nTelescopeHits->setTitle( "Telescope hits/event;telescope hits;events");

    hist1D_nUpstreamTriplets = AIDAProcessor::histogramFactory(this)->
      createHistogram1D("nUpstreamTriplets", 21, -0.5, 20.5);
    hist1D_nUpstreamTriplets->setTitle( "Number of Upstream Triplets per Event; "
					"Number of Upstream Triplets;Events");

    hist1D_nDownstreamTriplets = AIDAProcessor::histogramFactory(this)->
      createHistogram1D("nDownstreamTriplets", 21, -0.5, 20.5);
    hist1D_nDownstreamTriplets->setTitle( "Number of Downstream Triplets per Event; "
					  "Number of Downstream Triplets;Events");

    // GBL:
    hist1D_gblNdfAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D("gblNdf", 16, -0.5, 15.5);
    hist1D_gblNdfAlign->setTitle( "GBL fit NDF;GBL NDF;tracks" );

    hist1D_gblChi2Align = AIDAProcessor::histogramFactory(this)->
      createHistogram1D("gblChi2", 100, 0, 10);
    hist1D_gblChi2Align->setTitle( "GBL fit chi2 / degrees of freedom ;GBL chi2/Ndf ;tracks");

    hist1D_gblProbAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D("gblProb", 100, 0, 1);
    hist1D_gblProbAlign->setTitle( "GBL fit probability;GBL fit probability;tracks");

    AIDAProcessor::tree(this)->mkdir("GBLFit");
    AIDAProcessor::tree(this)->mkdir("GBLFit/Angles");
    AIDAProcessor::tree(this)->mkdir("GBLFit/Residuals");
    AIDAProcessor::tree(this)->mkdir("GBLFit/Pulls");
    AIDAProcessor::tree(this)->mkdir("GBLFit/Kinks");
    
    //[START] loop over sensor IDs
    for(size_t ix = 0; ix<sensorIDVec.size(); ++ix) 
      {
	auto sensorIdString = std::to_string(sensorIDVec[ix]);
	
	std::string histNameAngleX = "GBLFit/Angles/ax"+sensorIdString;
	std::string histNameAngleY = "GBLFit/Angles/ay"+sensorIdString;
	
	hist1D_gblAngleX.push_back(AIDAProcessor::histogramFactory(this)->
				   createHistogram1D( histNameAngleX, 100, -5, 5));
	hist1D_gblAngleX.back()->setTitle( "GBL angle at plane "+sensorIdString
					   +";x angle at plane "+sensorIdString+" [mrad];tracks"); 
	
	hist1D_gblAngleY.push_back(AIDAProcessor::histogramFactory(this)->
				   createHistogram1D( histNameAngleY, 100, -5, 5));
	hist1D_gblAngleY.back()->setTitle( "GBL angle at plane "+sensorIdString
					   +";y angle at plane "+sensorIdString+" [mrad];tracks"); 
	
	std::string histNameResidX = "GBLFit/Residuals/rx"+sensorIdString; 
	std::string histNameResidY = "GBLFit/Residuals/ry"+sensorIdString; 
	
	hist1D_gblResidX.push_back(AIDAProcessor::histogramFactory(this)->
				   createHistogram1D( histNameResidX, 500, -250, 250));
	hist1D_gblResidX.back()->setTitle( "GBL residual at plane "+sensorIdString
					   +";x resid at plane "+sensorIdString+" [#mum];tracks"); 
	
	hist1D_gblResidY.push_back(AIDAProcessor::histogramFactory(this)->
				   createHistogram1D( histNameResidY, 500, -250, 250));
	hist1D_gblResidY.back()->setTitle( "GBL residual at plane "+sensorIdString
					   +";y resid at plane "+sensorIdString+" [#mum];tracks"); 
	
	std::string histNamePullX = "GBLFit/Pulls/px"+sensorIdString; 
	std::string histNamePullY = "GBLFit/Pulls/py"+sensorIdString; 
	
	hist1D_gblPullX.push_back(AIDAProcessor::histogramFactory(this)->
				  createHistogram1D( histNamePullX, 100, -5, 5));
	hist1D_gblPullX.back()->setTitle( "GBL pull at plane "+sensorIdString
					  +";x pull at plane "+sensorIdString+";tracks"); 
	
	hist1D_gblPullY.push_back(AIDAProcessor::histogramFactory(this)->
				  createHistogram1D( histNamePullY, 100, -5, 5));
	hist1D_gblPullY.back()->setTitle( "GBL pull at plane "+sensorIdString
					  +";y pull at plane "+sensorIdString+";tracks"); 
	
	std::string histNameKinkX = "GBLFit/Kinks/kx"+sensorIdString;
	std::string histNameKinkY = "GBLFit/Kinks/ky"+sensorIdString;
	
	hist1D_gblKinkX.push_back(AIDAProcessor::histogramFactory(this)->
				  createHistogram1D( histNameKinkX, 100, -5, 5));
	hist1D_gblKinkX.back()->setTitle( "GBL kink angle at plane "+sensorIdString
					  +";plane "+sensorIdString+" x kink [mrad];tracks");
	
	hist1D_gblKinkY.push_back(AIDAProcessor::histogramFactory(this)->
				  createHistogram1D( histNameKinkY, 100, -5, 5));
	hist1D_gblKinkY.back()->setTitle( "GBL kink angle at plane "+sensorIdString
					  +";plane "+sensorIdString+" y kink [mrad];tracks");

	if(sensorIDVec[ix]==_SUT_ID) {
	  std::string histNameKinkXvsXY = "GBLFit/Kinks/kinkx_vs_xy"+sensorIdString;
	  std::string histNameKinkYvsXY = "GBLFit/Kinks/kinky_vs_xy"+sensorIdString;

	  profile2D_gblSUTKinkXvsXY = AIDAProcessor::histogramFactory(this)->
	    createProfile2D( histNameKinkXvsXY, 120, -12, 12, 60, -6, 6, 0, 100);
	  profile2D_gblSUTKinkXvsXY->setTitle( "GBL kink angle at plane "+sensorIdString
					       +"; x pos [mm]; y pos [mm]; sqrt(<kinkX^{2}>) [mrad]");

	  profile2D_gblSUTKinkYvsXY = AIDAProcessor::histogramFactory(this)->
	    createProfile2D( histNameKinkYvsXY, 120, -12, 12, 60, -6, 6, 0, 100);
	  profile2D_gblSUTKinkYvsXY->setTitle( "GBL kink angle at plane "+sensorIdString
					       +"; x pos [mm]; y pos [mm]; sqrt(<kinkY^{2}>) [mrad]");
	}

      }//[END] loop over sensor IDs
    
    hist1D_nTracksPerEvent = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "ntracksperevent", 21, -0.5, 20.5 );
    hist1D_nTracksPerEvent->setTitle( "Matched Tracks;Track Matches in an Event;Events" );
  
  } catch( lcio::Exception& e ) {
    //FIXME: fill in here something?!
  }
}
