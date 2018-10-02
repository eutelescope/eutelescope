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

// GBL:
#include "include/GblTrajectory.h"
#include "include/MilleBinary.h"

#include "EUTelTripletGBLUtility.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"
#include "marlin/Exceptions.h"
#include "marlin/AIDAProcessor.h"

// marlin util includes
#include "Mille.h"

// aida includes <.h>
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



//------------------------------------------------------------------------------
EUTelAlignGBL::EUTelAlignGBL(): Processor("EUTelGBL") {

  // modify processor description
  _description = "EUTelAlignGBL uses the MILLE program to write data files for MILLEPEDE II.";

  registerInputCollections(LCIO::TRACKERHIT,"hitCollectionName","Input hit collections name",_hitCollectionName, std::vector<std::string>{"corrhits"});
  registerProcessorParameter("eBeam","Beam energy [GeV]",_eBeam, 4.0);
  registerOptionalParameter("excludePlanes","Exclude planes from fit according to their sensor ids",_excludePlanes_sensorIDs ,std::vector<int>());
  registerOptionalParameter("upstreamTriplet","The three sensors used as the upstream triplet", _upstream_triplet_ids, std::vector<int>{0,1,2});
  registerOptionalParameter("downstreamTriplet","The three sensors used as the downstream triplet", _downstream_triplet_ids, std::vector<int>{3,4,5});
  registerOptionalParameter("lastUpstreamSensor","The last plane (z-ordered) which still should be attached to the upstream triplet", _last_upstream_sensor, 2);
  registerOptionalParameter("resolutionX","x-resolution of sensors (z-ordered) [mm]", _x_resolution_vec ,std::vector<float>());
  registerOptionalParameter("resolutionY","y-resolution of sensors (z-ordered) [mm]", _y_resolution_vec ,std::vector<float>());
  registerOptionalParameter("fixedPlanes","Fix sensor planes in the fit according to their sensor ids (it is recommended to fix two telescope planes)",_FixedPlanes_sensorIDs ,std::vector<int>{0,5});
  registerOptionalParameter("maxTrackCandidatesTotal","Maximal number of track candidates (Total)",_maxTrackCandidatesTotal, 10000000);
  registerOptionalParameter("maxTrackCandidates","Maximal number of track candidates",_maxTrackCandidates, 2000);
  registerOptionalParameter("milleBinaryFilename","Name of the Millepede binary file",_binaryFilename, std::string{"mille.bin"});
  registerOptionalParameter("alignMode","Number of alignment constants used. Available mode are:"
                              "\n\t\tXYZShifts - shifts in X and Y"
                              "\n\t\tXYShiftsRotZ - shifts in X and Y and rotation around the Z axis,"
                              "\n\t\tXYZShiftsRotZ - shifts in X,Y and Z and rotation around the Z axis"
                              "\n\t\tXYZShiftsRotXYZ - all shifts and rotations allowed",
                              _alignModeString, std::string{ "XYShiftsRotZ" });
  registerOptionalParameter("performAlignment","Set to 0 if you do not want to perform the alignment",_performAlignment, 1);
  registerOptionalParameter("fixedXShift","List of planes which should be fixed in X direction",_FixedXShift ,std::vector<int>());
  registerOptionalParameter("fixedYShift","List of planes which should be fixed in Y direction",_FixedYShift ,std::vector<int>());
  registerOptionalParameter("fixedZShift","List of planes which should be fixed in Z direction",_FixedZShift ,std::vector<int>());
  registerOptionalParameter("fixedZRot","List of planes which should have a fixed Z rotation",_FixedZRot ,std::vector<int>());
  registerOptionalParameter("fixedYRot","List of planes which should have a fixed Y rotation",_FixedYRot ,std::vector<int>{0,1,2,3,4,5});
  registerOptionalParameter("fixedXRot","List of planes which should have a fixed X rotation",_FixedXRot ,std::vector<int>{0,1,2,3,4,5});
  registerOptionalParameter("upstreamTripletResidualCut", "Upstream triplet residual cut [mm]", _upTriResCut, 0.30);
  registerOptionalParameter("downstreamTripletResidualCut", "Downstream triplet residual cut [mm]", _downTriResCut, 0.40);
  registerOptionalParameter("upstreamTripletSlopeCut", "Upstream triplet slope cut [mrad]", _upSlopeCut, 3.);
  registerOptionalParameter("downstreamTripletSlopeCut", "Downstream triplet slope cut [mrad]", _downSlopeCut, 5.);
  registerOptionalParameter("tripletsMatchingCut", "Upstream-downstream triplet matching cut [mm]", _upDownTripletMatchCut, 0.60);
  registerOptionalParameter("DUTCuts", "Cuts in x and y for matching DUT hits [mm]", _DUTCuts, std::vector<float>{1.,1.});
  registerOptionalParameter("zMid", "Z Position used for triplet matching (default to center between end of first triplet and start of second triplet)", _zMid, -1.0);
  registerOptionalParameter("chi2Cut", "Cut on the chi2 for the tracks to be passed to Millepede", _chi2cut, 0.001);
  registerOptionalParameter("pedeSteerfileName","Name of the steering file for the pede program",_pedeSteerfileName, std::string{"steer_mille.txt"});
  registerProcessorParameter("kappa","Global factor to Highland formula, 1.0 means HL as is, 1.2 means 20/% additional scattering", _kappa, 1.0);
}

//------------------------------------------------------------------------------
void EUTelAlignGBL::init() {
  geo::gGeometry().initializeTGeoDescription(EUTELESCOPE::GEOFILENAME,
                                             EUTELESCOPE::DUMPGEOROOT);

  //This is the vector of sensorIDs ordered alogn the gloabl z-axis, 
  //this is guranteed by the framework 
  _sensorIDVec = geo::gGeometry().sensorIDsVec();
  _nPlanes = _sensorIDVec.size();

  bool isStillUpstream = true;
  for(auto& sensorID: _sensorIDVec) {
    _planePosition.emplace_back( geo::gGeometry().getPlaneZPosition(sensorID) );
    auto z = geo::gGeometry().getPlaneZSize(sensorID);
    auto rad = geo::gGeometry().getPlaneRadiationLength(sensorID);
    _planeRadLength.emplace_back(z/rad);
    
    //Since we go along z-ordering, we can decide if a sensor still is associated to the up or downstream triplet
    _is_sensor_upstream[sensorID] = isStillUpstream;
    if(sensorID == _last_upstream_sensor) {
       isStillUpstream = false;
    }

    auto belongsToUp = (std::find(std::begin(_upstream_triplet_ids), std::end(_upstream_triplet_ids), sensorID) != _upstream_triplet_ids.end());
    auto belongsToDown = (std::find(std::begin(_downstream_triplet_ids), std::end(_downstream_triplet_ids), sensorID) != _downstream_triplet_ids.end());

    if(!belongsToUp && !belongsToDown){
      _dut_ids.emplace_back(sensorID);
    }
  }

  //to compute the total radiation length we will loop over all planes and add radiation
  //length for air from the first to the last plane
  double totalRadLength = (_planePosition.back()-_planePosition.front())/304200; //Radiation length for air is 30420 cm
  for(auto& radLen: _planeRadLength) {
    totalRadLength += radLen;
  }

  for(auto& radLen: _planeRadLength) {
      // Paper showed HL predicts too high angle, at least for biased measurement. 
      double tetSi = _kappa*0.0136 * sqrt(radLen) / _eBeam * ( 1 + 0.038*std::log(totalRadLength) );
      _planeWscatSi.emplace_back( 1.0/(tetSi*tetSi), 1.0/(tetSi*tetSi) );
  }

  for(size_t ipl = 0; ipl < _nPlanes-1; ipl++) {
    auto distplane = _planePosition[ipl+1] - _planePosition[ipl];
    double epsAir =   0.5*distplane  / 304200.;  // in [mm]
    double tetAir = _kappa*0.0136 * sqrt(epsAir) / _eBeam * ( 1 + 0.038*std::log(totalRadLength) );
    _planeWscatAir.emplace_back( 1.0/(tetAir*tetAir), 1.0/(tetAir*tetAir) );
  }

  streamlog_out(MESSAGE4)   << "assumed beam energy " << _eBeam << " GeV" <<  endl;
  streamlog_out(MESSAGE4)   << "Summary:\n\t_planeRadLength has size: \t" << _planeRadLength.size() 
                            << "\n\t_planeWscatSi has size: \t" << _planeWscatSi.size() 
                            << "\n\t_planeWscatAir has size: \t" << _planeWscatAir.size() << '\n';

  for(size_t ipl = 0; ipl < _nPlanes; ipl++) {
    auto x_res = _x_resolution_vec.at(ipl);
    auto y_res = _y_resolution_vec.at(ipl);
    _planeMeasPrec.emplace_back(1.0/x_res/x_res, 1.0/y_res/y_res);
    
  }

  // the user is giving sensor ids for the planes to be fixed.
  // These sensor ids have to be converted to a local index
  // according to the planes positions along the z axis, i.e. we
  // need to fill _FixedPlanes with the z-position which corresponds
  // to the postion in _sensorIDVec
  cout << "FixedPlanes " << _FixedPlanes_sensorIDs.size() << endl;
  for(auto& fixedPlaneID: _FixedPlanes_sensorIDs) {
    auto planeIt = std::find(_sensorIDVec.begin(), _sensorIDVec.end(), fixedPlaneID); 
    _FixedPlanes.emplace_back( planeIt - _sensorIDVec.begin() );
  }

  // same for excluded planes
  cout << "ExcludedPlanes " << _excludePlanes_sensorIDs.size() << endl;
  for(auto& excludedPlaneID: _excludePlanes_sensorIDs) {
    auto planeIt = std::find(_sensorIDVec.begin(), _sensorIDVec.end(), excludedPlaneID);
    _excludePlanes.emplace_back( planeIt - _sensorIDVec.begin() );
  }

  int counter = 0;
  for(auto& sensorID: _sensorIDVec) { 
    if(std::find(_excludePlanes_sensorIDs.begin(), _excludePlanes_sensorIDs.end(), sensorID) 
    == _excludePlanes_sensorIDs.end()) {
      indexconverter.emplace_back(counter++);
    } else {
      indexconverter.emplace_back(-1);
    }
  }

  // Center of triplets calculation
  if(_zMid < 0 ){
    double zEndUpstreamTriplet = 0;
    double zStartDownstreamTriplet = 0;
    for (size_t ipl = 0; ipl < _nPlanes; ipl++){
	  if(_sensorIDVec[ipl] == _upstream_triplet_ids.at(2)) zEndUpstreamTriplet=_planePosition[ipl];
      else if(_sensorIDVec[ipl] == _downstream_triplet_ids.at(0)) zStartDownstreamTriplet=_planePosition[ipl];
    }
    _zMid = 0.5*(zEndUpstreamTriplet + zStartDownstreamTriplet);
  }
  
  streamlog_out( MESSAGE4 ) << "Extrapolated tracks from the triplets will be matched in z position = " << _zMid << endl;

  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;
  _printEventCounter = 0;

  // Initialize number of excluded planes
  _nExcludePlanes = _excludePlanes.size();

  streamlog_out( MESSAGE2 ) << "Number of planes excluded from the alignment fit: " << _nExcludePlanes << endl;

  // Initialise Mille statistics:
  _nMilleTracks = 0;

  // booking histograms
  bookHistos(_sensorIDVec);
  gblutil.setParent(this);
  gblutil.bookHistos();
  
  // to be used only during alignment
  if(_performAlignment){ 
    streamlog_out( MESSAGE2 ) << "Initialising Mille..." << endl;

    unsigned int reserveSize = 8000;
    milleAlignGBL = std::make_unique<gbl::MilleBinary>( _binaryFilename, reserveSize );

    streamlog_out( MESSAGE2 ) << "The filename for the binary file is: " << _binaryFilename.c_str() << endl;

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
    streamlog_out( MESSAGE4 ) << "Generating the steering file for the pede program..." << endl;

    ofstream steerFile;
    steerFile.open(_pedeSteerfileName.c_str());

    if( steerFile.is_open()) {

      // find first and last excluded plane
      auto firstnotexcl = _nPlanes;
      size_t lastnotexcl = 0;

      // loop over all planes:

      for( size_t ipl = 0; ipl < _nPlanes; ipl++) {

    size_t excluded = 0;

    // loop over excluded planes:

    for( size_t jpl = 0; jpl < _nExcludePlanes; jpl++ ) {
      if( ipl == _excludePlanes[jpl] ) excluded = 1;
    }

    if( excluded == 0 && firstnotexcl > ipl ) firstnotexcl = ipl;

    if( excluded == 0 && lastnotexcl < ipl ) lastnotexcl = ipl;

      } // end loop over all planes

      steerFile << "Cfiles" << endl;
      steerFile << _binaryFilename << endl;
      steerFile << endl;

      steerFile << "Parameter" << endl;

      int nfix = 0;

      // loop over all planes:

      for( size_t ipl = 0; ipl < _nPlanes; ipl++) {

    int excluded = 0; // flag for excluded planes

    // check in list of excluded planes:

    for( size_t iex = 0; iex < _nExcludePlanes; iex++) {
      if( ipl == _excludePlanes[iex] )
        excluded = 1;
    }

    cout << "Plane " << _sensorIDVec[ipl] << " exclude = " << excluded << endl;

    if( excluded == 0 ) {

      bool fixed = false;
      for( size_t i = 0;i< _FixedPlanes.size(); i++ ) {
        if( _FixedPlanes[i] == ipl )
          fixed = true;
      }

      // if fixed planes
      // all the steer writing can be cleaner
      // Label corresponds to plane ID*10 + X, where X (last digit) is between 1 and 6
      // 1 = X shift, 2 = Y shift, 3 = Z shift, 4 = ZX rot, 5 = ZY rot, 6 = XY rot
      // For example, the Z shift of plane 55 will be labelled as 55 * 10 + 3 = 553
      
      if( fixed || (_FixedPlanes.size() == 0 && (ipl == firstnotexcl || ipl == lastnotexcl) ) ) {
        nfix++;
        if( _alignMode == Utility::alignMode::XYShifts ) {
          steerFile << (_sensorIDVec[ipl] * 10 + 1) << "  0.0 -1.0" << endl;
          steerFile << (_sensorIDVec[ipl] * 10 + 2) << "  0.0 -1.0" << endl;
        }
        if( _alignMode == Utility::alignMode::XYShiftsRotZ ) {
          steerFile << (_sensorIDVec[ipl] * 10 + 1) << "  0.0 -1.0" << endl; 
          steerFile << (_sensorIDVec[ipl] * 10 + 2) << "  0.0 -1.0" << endl; 
          steerFile << (_sensorIDVec[ipl] * 10 + 3) << "  0.0 -1.0" << endl; 
        }
        if( _alignMode == Utility::alignMode::XYZShiftsRotZ ) {
          steerFile << (_sensorIDVec[ipl] * 10 + 1) << "  0.0 -1.0" << endl;
          steerFile << (_sensorIDVec[ipl] * 10 + 2) << "  0.0 -1.0" << endl;
          steerFile << (_sensorIDVec[ipl] * 10 + 3) << "  0.0 -1.0" << endl;
          steerFile << (_sensorIDVec[ipl] * 10 + 4) << "  0.0 -1.0" << endl;
        }
        if( _alignMode == Utility::alignMode::XYZShiftsRotXYZ ) {
          steerFile << (_sensorIDVec[ipl] * 10 + 1) << "  0.0 -1.0" << endl;
          steerFile << (_sensorIDVec[ipl] * 10 + 2) << "  0.0 -1.0" << endl;
          steerFile << (_sensorIDVec[ipl] * 10 + 3) << "  0.0 -1.0" << endl;
          steerFile << (_sensorIDVec[ipl] * 10 + 4) << "  0.0 -1.0" << endl;
          steerFile << (_sensorIDVec[ipl] * 10 + 5) << "  0.0 -1.0" << endl;
          steerFile << (_sensorIDVec[ipl] * 10 + 6) << "  0.0 -1.0" << endl;
        }
      }

      else {

        if( _alignMode == Utility::alignMode::XYShifts ) {
          if(std::find(_FixedXShift.begin(), _FixedXShift.end(), _sensorIDVec[ipl]) == _FixedXShift.end()) {
            steerFile << (_sensorIDVec[ipl] * 10 + 1) << "  0.0  0.0" << endl;
	      } else {
			steerFile << (_sensorIDVec[ipl] * 10 + 1) << "  0.0  -1.0" << endl; 
		  }
          if(std::find(_FixedYShift.begin(), _FixedYShift.end(), _sensorIDVec[ipl]) == _FixedYShift.end()) {
            steerFile << (_sensorIDVec[ipl] * 10 + 2) << "  0.0  0.0" << endl;
	      } else {
			steerFile << (_sensorIDVec[ipl] * 10 + 2) << "  0.0  -1.0" << endl; 
		  }
        }

        if( _alignMode == Utility::alignMode::XYShiftsRotZ ) {
          if(std::find(_FixedXShift.begin(), _FixedXShift.end(), _sensorIDVec[ipl]) == _FixedXShift.end()) {
            steerFile << (_sensorIDVec[ipl] * 10 + 1) << "  0.0  0.0" << endl;
	      } else {
			steerFile << (_sensorIDVec[ipl] * 10 + 1) << "  0.0  -1.0" << endl; 
		  }
          if(std::find(_FixedYShift.begin(), _FixedYShift.end(), _sensorIDVec[ipl]) == _FixedYShift.end()) {
            steerFile << (_sensorIDVec[ipl] * 10 + 2) << "  0.0  0.0" << endl;
	      } else {
			steerFile << (_sensorIDVec[ipl] * 10 + 2) << "  0.0  -1.0" << endl; 
		  }
          if(std::find(_FixedZRot.begin(), _FixedZRot.end(), _sensorIDVec[ipl]) == _FixedZRot.end()) {
            steerFile << (_sensorIDVec[ipl] * 10 + 3) << "  0.0  0.0" << endl;
	      } else {
			steerFile << (_sensorIDVec[ipl] * 10 + 3) << "  0.0  -1.0" << endl; 
		  }
        }

        if( _alignMode == Utility::alignMode::XYZShiftsRotZ ) {
		  if(std::find(_FixedXShift.begin(), _FixedXShift.end(), _sensorIDVec[ipl]) == _FixedXShift.end()) {
            steerFile << (_sensorIDVec[ipl] * 10 + 1) << "  0.0  0.0" << endl;
	      } else {
			steerFile << (_sensorIDVec[ipl] * 10 + 1) << "  0.0  -1.0" << endl; 
		  }
		  if(std::find(_FixedYShift.begin(), _FixedYShift.end(), _sensorIDVec[ipl]) == _FixedYShift.end()) {
            steerFile << (_sensorIDVec[ipl] * 10 + 2) << "  0.0  0.0" << endl;
	      } else {
			steerFile << (_sensorIDVec[ipl] * 10 + 2) << "  0.0  -1.0" << endl; 
		  }
		  if(std::find(_FixedZShift.begin(), _FixedZShift.end(), _sensorIDVec[ipl]) == _FixedZShift.end()) { //Z rotation in 3 and Z shift in 4??
            steerFile << (_sensorIDVec[ipl] * 10 + 4) << "  0.0  0.0" << endl;
	      } else {
			steerFile << (_sensorIDVec[ipl] * 10 + 4) << "  0.0  -1.0" << endl; 
		  }
          if(std::find(_FixedZRot.begin(), _FixedZRot.end(), _sensorIDVec[ipl]) == _FixedZRot.end()) {
            steerFile << (_sensorIDVec[ipl] * 10 + 3) << "  0.0  0.0" << endl;
	      } else {
			steerFile << (_sensorIDVec[ipl] * 10 + 3) << "  0.0  -1.0" << endl; 
		  }
        }
        if( _alignMode == Utility::alignMode::XYZShiftsRotXYZ ) {
		  if(std::find(_FixedXShift.begin(), _FixedXShift.end(), _sensorIDVec[ipl]) == _FixedXShift.end()) {
            steerFile << (_sensorIDVec[ipl] * 10 + 1) << "  0.0  0.0" << endl;
	      } else {
			steerFile << (_sensorIDVec[ipl] * 10 + 1) << "  0.0  -1.0" << endl; 
		  }
		  if(std::find(_FixedYShift.begin(), _FixedYShift.end(), _sensorIDVec[ipl]) == _FixedYShift.end()) {
            steerFile << (_sensorIDVec[ipl] * 10 + 2) << "  0.0  0.0" << endl;
	      } else {
			steerFile << (_sensorIDVec[ipl] * 10 + 2) << "  0.0  -1.0" << endl; 
		  }
		  if(std::find(_FixedZShift.begin(), _FixedZShift.end(), _sensorIDVec[ipl]) == _FixedZShift.end()) {
            steerFile << (_sensorIDVec[ipl] * 10 + 3) << "  0.0  0.0" << endl;
	      } else {
			steerFile << (_sensorIDVec[ipl] * 10 + 3) << "  0.0  -1.0" << endl; 
		  }
		  if(std::find(_FixedXRot.begin(), _FixedXRot.end(), _sensorIDVec[ipl]) == _FixedXRot.end()) {
            steerFile << (_sensorIDVec[ipl] * 10 + 4) << "  0.0  0.0" << endl; // LABELLING TO BE CHECKED
	      } else {
			steerFile << (_sensorIDVec[ipl] * 10 + 4) << "  0.0  -1.0" << endl; 
		  }
		  if(std::find(_FixedYRot.begin(), _FixedYRot.end(), _sensorIDVec[ipl]) == _FixedYRot.end()) {
            steerFile << (_sensorIDVec[ipl] * 10 + 5) << "  0.0  0.0" << endl; // LABELLING TO BE CHECKED
	      } else {
			steerFile << (_sensorIDVec[ipl] * 10 + 5) << "  0.0  -1.0" << endl; 
		  }
		  if(std::find(_FixedZRot.begin(), _FixedZRot.end(), _sensorIDVec[ipl]) == _FixedZRot.end()) {
            steerFile << (_sensorIDVec[ipl] * 10 + 6) << "  0.0  0.0" << endl; // LABELLING TO BE CHECKED
	      } else {
			steerFile << (_sensorIDVec[ipl] * 10 + 6) << "  0.0  -1.0" << endl; 
		  }
        }

      }// not fixed

    } // end if plane not excluded

      } // end loop over all planes

      steerFile << endl;
      steerFile << "! chiscut 5.0 2.5" << endl;
      steerFile << "outlierdownweighting 4" << endl;
      steerFile << "dwfractioncut 0.1" << endl;
      steerFile << endl;
      steerFile << "method inversion 10  0.1" << endl;
      // Use 10 OpenMP threads to process the data:
      steerFile << "threads 10 1" << endl;
      steerFile << endl;
      steerFile << "histprint" << endl;
      steerFile << endl;
      steerFile << "end" << endl;

      steerFile.close();

      streamlog_out( MESSAGE2 ) << "File " << _pedeSteerfileName << " written." << endl;

    }
    else {
      streamlog_out( ERROR2 ) << "Could not open steering file." << endl;
    }
  
  // end writing the pede steering file
  }
  
  streamlog_out( MESSAGE2 ) << "end of init" << endl;
}

//------------------------------------------------------------------------------
void EUTelAlignGBL::processRunHeader( LCRunHeader* rdr ) {
  auto header = std::make_unique<EUTelRunHeaderImpl>(rdr);
  header->addProcessor( type() ) ;
  // increment the run counter
  ++_iRun;
}

//------------------------------------------------------------------------------
inline Eigen::Matrix<double,5,5> Jac55new( double ds ) {
  /* for GBL:
     Jacobian for straight line track
     track = q/p, x', y', x, y
     0,   1,  2,  3, 4
     */
  Eigen::Matrix<double,5,5> jac = Eigen::Matrix<double,5,5>::Identity();
  //jac.UnitMatrix();
  jac(3,1) = ds; // x = x0 + xp * ds
  jac(4,2) = ds; // y = y0 + yp * ds
  return jac;
}

//------------------------------------------------------------------------------
void EUTelAlignGBL::processEvent( LCEvent * event ) {

  if( _iEvt % 1000 == 0 ) {
    streamlog_out( MESSAGE2 ) << "Processing event "
      << setw(6) << setiosflags(ios::right)
      << event->getEventNumber() << " in run "
      << setw(6) << setiosflags(ios::right)
      << event->getRunNumber()
      << ", currently having "
      << _nMilleTracks << " tracks "
      << endl;
  }

  if( _nMilleTracks > static_cast<size_t>(_maxTrackCandidatesTotal) ) {
    throw StopProcessingException(this);
  }
  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event) ;

  if( evt->getEventType() == kEORE ) {
    streamlog_out( DEBUG2 ) << "EORE found: nothing else to do." << endl;
    return;
  }

  CellIDDecoder<TrackerHit> hitCellDecoder(EUTELESCOPE::HITENCODING);
  std::vector<EUTelTripletGBLUtility::hit> _hitsVec;
  std::vector<EUTelTripletGBLUtility::hit> _DUThitsVec;

  for( size_t i = 0; i < _hitCollectionName.size(); i++ ) {
    LCCollection* collection = nullptr;
    try {
      collection = event->getCollection(_hitCollectionName[i]);
    } catch (DataNotAvailableException& e) {
      throw SkipEventException(this);
    }
    // loop over all hits in collection:

    if(_printEventCounter < NO_PRINT_EVENT_COUNTER) {
      streamlog_out(DEBUG2) << "Event " << event->getEventNumber() << " contains " 
      << collection->getNumberOfElements() << " hits" << endl;
    }
    
    nTelHits->fill(collection->getNumberOfElements());

    for( int iHit = 0; iHit < collection->getNumberOfElements(); iHit++ ) {
      auto hit = static_cast<TrackerHitImpl*>( collection->getElementAt(iHit) );
      auto sensorID = hitCellDecoder(hit)["sensorID"];
      auto hitPosition = hit->getPosition();
      if(sensorID <= 6) {
        _hitsVec.emplace_back(hitPosition, sensorID);
      } else {
        _DUThitsVec.emplace_back(hitPosition, sensorID);
      }
      if(_printEventCounter < NO_PRINT_EVENT_COUNTER) streamlog_out(DEBUG0) << "Hit on plane " << sensorID << " at " << hitPosition[0] << "|" << hitPosition[1]  << '\n';
    } // end loop over all hits in given collection
  }//loop over all input hit collection

  int numbertracks = 0;
  auto tripletVec = std::vector<EUTelTripletGBLUtility::triplet>();
  auto dripletVec = std::vector<EUTelTripletGBLUtility::triplet>();

  gblutil.FindTriplets(_hitsVec, _upstream_triplet_ids, _upTriResCut, _upSlopeCut/1000., tripletVec, false, true);
  gblutil.FindTriplets(_hitsVec, _downstream_triplet_ids, _downTriResCut, _downSlopeCut/1000., dripletVec, false, false);

  if(_printEventCounter < NO_PRINT_EVENT_COUNTER){
    streamlog_out(DEBUG2)  << "Triplets:\n";
    for(auto& triplet: tripletVec) {
      streamlog_out(DEBUG2) << triplet << '\n';
    }
    streamlog_out(DEBUG2) << "Driplets:\n";
    for(auto& driplet: dripletVec) {
      streamlog_out(DEBUG2) << driplet << '\n';
    }
  }

  nUpstreamTriplets->fill( tripletVec.size()  );
  nDownstreamTriplets->fill( dripletVec.size() );

  auto matchedTripletVec = std::vector<EUTelTripletGBLUtility::track>();
  gblutil.MatchTriplets(tripletVec, dripletVec, _zMid, _upDownTripletMatchCut, matchedTripletVec);

  if(_printEventCounter < NO_PRINT_EVENT_COUNTER) streamlog_out(DEBUG2) << "Matched to:\n"; 

  if(!_dut_ids.empty()) {
          for(auto& track: matchedTripletVec) {    
            auto& triplet = track.get_upstream();
            auto& driplet = track.get_downstream();

            for(auto dutID: _dut_ids) {
              if(_is_sensor_upstream[dutID]) {
                gblutil.AttachDUT(triplet, _DUThitsVec, dutID, _DUTCuts);
              } else {
                gblutil.AttachDUT(driplet, _DUThitsVec, dutID, _DUTCuts);
              }
            }
          }
  }
  for(auto& track: matchedTripletVec) {
    // GBL point vector for the trajectory, all in [mm] !!
    // GBL with triplet A as seed
    std::vector<gbl::GblPoint> traj_points;
    // build up trajectory:
    std::vector<unsigned int> ilab;
    vector<double> sPoint;

    // the arc length at the first measurement plane is 0.
    double s = 0;

    Eigen::Matrix2d proL2m = Eigen::Matrix2d::Identity();
    Eigen::Vector2d scat = Eigen::Vector2d::Zero(); //mean is zero

    auto triplet = track.get_upstream();
    auto driplet = track.get_downstream();
    //We need the _T_riplet slope to compute residual
    auto triSlope = triplet.slope();

    std::vector<EUTelTripletGBLUtility::hit> DUThits;
    for(auto it = triplet.DUT_begin(); it != triplet.DUT_end(); ++it){
      DUThits.emplace_back(it->second);
    }
    for(auto it = driplet.DUT_begin(); it != driplet.DUT_end(); ++it){
      DUThits.emplace_back(it->second);
    }

    if(_printEventCounter < NO_PRINT_EVENT_COUNTER) streamlog_out(DEBUG2) << "Track has " << DUThits.size() << " DUT hits\n";

  // to be used only during alignment. Matrix defined outside if clause to avoid complaints from compiler. Could be done better
    Eigen::Matrix2d alDer2; // alignment derivatives
    Eigen::Matrix<double,2,3> alDer3; // alignment derivatives
    Eigen::Matrix<double, 2,4> alDer4; // alignment derivatives
    Eigen::Matrix<double, 3,6> alDer6; // alignment derivatives
    
    if(_performAlignment){ 
      alDer2(0,0) = 1.0; // dx/dx GBL sign convetion
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
      alDer4(0,3) = triSlope.x; // dx/dz
      alDer4(1,3) = triSlope.y; // dx/dz
    
      alDer6 ( 0, 0 ) = 1.0; // dx/dx
	  alDer6 ( 0, 1 ) = 0.0; // dx/dy
	  alDer6 ( 0, 2 ) = triSlope.x; // dx/dz
	  alDer6 ( 0, 3 ) = 0.0; // dx/da
	  alDer6 ( 1, 0 ) = 0.0; // dy/dx
	  alDer6 ( 1, 1 ) = 1.0; // dy/dy
      alDer6 ( 1, 2 ) = triSlope.y; // dy/dz
	  alDer6 ( 1, 4 ) = 0.0; // dy/db
      alDer6 ( 2, 0 ) = 0.0; // dz/dx
	  alDer6 ( 2, 1 ) = 0.0; // dz/dy
	  alDer6 ( 2, 2 ) = 1.0; // dz/dz
	  alDer6 ( 2, 5 ) = 0.0; // dz/dg
    }
	
    std::vector<double> rx (_nPlanes, -1.0);
    std::vector<double> ry (_nPlanes, -1.0);
    std::vector<bool> hasHit (_nPlanes, false);
  
    double step = .0;
    unsigned int iLabel;

    for( size_t ipl = 0; ipl < _nPlanes; ++ipl ) {
      //We have to add all the planes, the up and downstream arm of the telescope will definitely have
      //hits, the DUTs might not though!
      EUTelTripletGBLUtility::hit const * hit = nullptr;
      auto sensorID = _sensorIDVec[ipl];
      if(std::find(_upstream_triplet_ids.begin(), _upstream_triplet_ids.end(), sensorID) != _upstream_triplet_ids.end()) {
        hit = &triplet.gethit(sensorID);
      } else if(std::find(_downstream_triplet_ids.begin(), _downstream_triplet_ids.end(), sensorID) != _downstream_triplet_ids.end()) {
		hit = &driplet.gethit(sensorID);
      }
      else if(triplet.has_DUT(sensorID)) hit = &triplet.get_DUT_Hit(sensorID);
      else if(driplet.has_DUT(sensorID)) hit = &driplet.get_DUT_Hit(sensorID);

      //if there is no hit we take the plane position from the geo description
      double zz = hit ? hit->z : _planePosition[ipl];// [mm]

      //transport matrix in (q/p, x', y', x, y) space
      auto jacPointToPoint = Jac55new( step );
      auto point = gbl::GblPoint( jacPointToPoint );
      s += step;

      //of there is a hit we will add a measurement to the point
      if(hit){
        hasHit[ipl] = true;
        double xs = triplet.getx_at(zz);
        double ys = triplet.gety_at(zz);

        if(_printEventCounter < NO_PRINT_EVENT_COUNTER) streamlog_out(DEBUG2) << "xs = " << xs << "   ys = " << ys << std::endl;

        rx[ipl] = (hit->x - xs); // resid hit-triplet
        ry[ipl] = (hit->y - ys); // resid
    
        Eigen::Vector2d meas;
        meas[0] = rx[ipl]; // fill meas vector for GBL
        meas[1] = ry[ipl];
        point.addMeasurement( proL2m, meas, _planeMeasPrec[ipl] );
        
        // to be used only during alignment
        if(_performAlignment){ 
          if( _alignMode == Utility::alignMode::XYShifts ) { // only x and y shifts
            // global labels for MP:
            std::vector<int> globalLabels(2);
            globalLabels[0] = _sensorIDVec[ipl] * 10 + 1;
            globalLabels[1] = _sensorIDVec[ipl] * 10 + 2;
            point.addGlobals( globalLabels, alDer2 ); // for MillePede alignment
          }
          else if( _alignMode == Utility::alignMode::XYShiftsRotZ ) { // with rot
            std::vector<int> globalLabels(3);
            globalLabels[0] = _sensorIDVec[ipl] * 10 + 1; // x
            globalLabels[1] = _sensorIDVec[ipl] * 10 + 2; // y
            globalLabels[2] = _sensorIDVec[ipl] * 10 + 3; // rot
            alDer3(0,2) = -ys; // dx/dphi
            alDer3(1,2) =  xs; // dy/dphi
            point.addGlobals( globalLabels, alDer3 ); // for MillePede alignment
          }
          else if( _alignMode == Utility::alignMode::XYZShiftsRotZ ) { // with rot and z shift
            std::vector<int> globalLabels(4);
            globalLabels[0] = _sensorIDVec[ipl] * 10 + 1;
            globalLabels[1] = _sensorIDVec[ipl] * 10 + 2;
            globalLabels[2] = _sensorIDVec[ipl] * 10 + 3;
            globalLabels[3] = _sensorIDVec[ipl] * 10 + 4; // z
            alDer4(0,2) = -ys; // dx/dphi
            alDer4(1,2) =  xs; // dy/dphi
            point.addGlobals( globalLabels, alDer4 ); // for MillePede alignment
          }
          else if( _alignMode == Utility::alignMode::XYZShiftsRotXYZ ) { // everything
		    double deltaz = hit->z - _planePosition[ipl];
            // deltaz cannot be zero, otherwise this mode doesn't work
            // 1e-9 in mm is quite small
            if ( deltaz < 1E-9 ) {
		       deltaz = 1E-9;
            }
		    std::vector<int> globalLabels(6);
            globalLabels[0] = _sensorIDVec[ipl] * 10 + 1;
            globalLabels[1] = _sensorIDVec[ipl] * 10 + 2;
            globalLabels[2] = _sensorIDVec[ipl] * 10 + 3;
            globalLabels[3] = _sensorIDVec[ipl] * 10 + 4;
            globalLabels[4] = _sensorIDVec[ipl] * 10 + 5;
            globalLabels[5] = _sensorIDVec[ipl] * 10 + 6;
            alDer6 ( 0, 4 ) = deltaz; // dx/db
            alDer6 ( 0, 5 ) = -ys; // dx/dg
            alDer6 ( 1, 3 ) = -deltaz; // dy/da
            alDer6 ( 1, 5 ) = xs; // dy/dg
            alDer6 ( 2, 3 ) = ys; // dz/da
            alDer6 ( 2, 4 ) = -xs; // dz/db
            point.addGlobals( globalLabels, alDer6 ); // for MillePede alignment
		  }
	    }
      }

      point.addScatterer( scat, _planeWscatSi[ipl] );
      sPoint.push_back( s );
      iLabel = sPoint.size();
      ilab.push_back(iLabel);
      traj_points.push_back(point);

      if( ipl < _nPlanes-1 ) {
        double distplane = _planePosition[ipl+1] - _planePosition[ipl];
        step = 0.21*distplane; // in [mm]
        auto point = gbl::GblPoint( Jac55new( step ) );
        point.addScatterer( scat, _planeWscatAir[ipl] );
        s += step;
        traj_points.push_back(point);
        sPoint.push_back( s );
        step = 0.58*distplane; // in [mm]
        auto point1 = gbl::GblPoint( Jac55new( step ) );
        point1.addScatterer( scat, _planeWscatAir[ipl] );
        s += step;
        traj_points.push_back(point1);
        sPoint.push_back( s );
        step = 0.21*distplane; // remaining distance to next plane, in [mm]
      }
    } // loop over planes

    // monitor what we put into GBL:
    double xA = triplet.getx_at(_zMid);
    double yA = triplet.gety_at(_zMid);
    double xB = driplet.getx_at(_zMid);
    double yB = driplet.gety_at(_zMid);

    selxHistGBLAlign->fill( xA ); // triplet at mid
    selyHistGBLAlign->fill( yA );
    selaxHistGBLAlign->fill( triSlope.x*1E3 );//track slope
    selayHistGBLAlign->fill( triSlope.y*1E3 );
    seldxHistGBLAlign->fill( (xB-xA)*1E3 ); // triplet-driplet match
    seldyHistGBLAlign->fill( (yB-yA)*1E3 );
    selkxHistGBLAlign->fill( track.kink_x()*1E3 ); // triplet-driplet kink
    selkyHistGBLAlign->fill( track.kink_y()*1E3 );

    //We need to fill the interpolation to plane 1 and extrapolation to 3, 4, 5 (6, 7, 8, ...)
    seldxSensorHistGBLAlign.at(0)->fill( rx[1]*1E3 );
    seldySensorHistGBLAlign.at(0)->fill( ry[1]*1E3 );
    for(size_t ipl = 3; ipl < _nPlanes; ++ipl) {//should be done for all the planes
      if(hasHit[ipl]){
        seldxSensorHistGBLAlign.at(ipl-2)->fill( rx[ipl]*1E3 );
        seldySensorHistGBLAlign.at(ipl-2)->fill( ry[ipl]*1E3 );
      }
    }

    double Chi2;
    int Ndf;
    double lostWeight;

    gbl::GblTrajectory traj(traj_points, false); // curvature = false
    traj.fit( Chi2, Ndf, lostWeight );
    //traj.getLabels(ilab); // instead pushback sPoint.size() when adding plane

    if(_printEventCounter < NO_PRINT_EVENT_COUNTER){
      streamlog_out(DEBUG4) << "traj with " << traj.getNumPoints() << " points:" << endl;
      for( size_t ipl = 0; ipl < sPoint.size(); ++ipl ){
        streamlog_out(DEBUG4) << "  GBL point " << ipl;
        streamlog_out(DEBUG4) << "  z " << sPoint[ipl]; 
        streamlog_out(DEBUG4) << endl;
      }
      for( size_t ipl = 0; ipl < 6; ++ipl ){
        streamlog_out(DEBUG4) << " plane " << ipl << ", lab " << ilab[ipl];
        streamlog_out(DEBUG4) << " z " << sPoint[ilab[ipl]-1];
        streamlog_out(DEBUG4) << "  dx " << rx[ipl];
        streamlog_out(DEBUG4) << "  dy " << ry[ipl];
        streamlog_out(DEBUG4) << "  chi2 " << Chi2;
        streamlog_out(DEBUG4) << "  ndf " << Ndf;
        streamlog_out(DEBUG4) << endl;
      }

      streamlog_out(DEBUG2)  << " Is traj valid? " << traj.isValid() << std::endl;
      //traj.printPoints();
      //traj.printTrajectory();
      //traj.printData();
      _printEventCounter++;
    }

    //cout << " chi2 " << Chi2 << ", ndf " << Ndf << endl;

    gblndfHistGBLAlign->fill( Ndf );
    gblchi2HistGBLAlign->fill( Chi2 / Ndf );
    double probchi = TMath::Prob( Chi2, Ndf );
    gblprbHistGBLAlign->fill( probchi );

    // look at fit:
    Eigen::VectorXd localPar;
    Eigen::MatrixXd localCov;

    std::array<double,8> ax;
    std::array<double,8> ay;
    size_t k = 0;

    // at plane 0:
    unsigned int ndata = 2;
    unsigned int ndim = 2;
    Eigen::VectorXd aResiduals(ndim);
    Eigen::VectorXd aMeasErrors(ndim);
    Eigen::VectorXd aResErrors(ndim);
    Eigen::VectorXd aDownWeights(ndim);

    for(size_t ix = 0; ix < _nPlanes; ++ix) {
      int ipos = ilab[ix];
      traj.getResults( ipos, localPar, localCov );

      //track = q/p, x', y', x, y
      //        0,   1,  2,  3, 4
      gblAxHist[ix]->fill( localPar[1]*1E3 );
      gblAyHist[ix]->fill( localPar[2]*1E3 );
      gblDxHist[ix]->fill( localPar[3]*1E3 );
      gblDyHist[ix]->fill( localPar[4]*1E3 );
      
      if(hasHit[ix]){
        traj.getMeasResults( ipos, ndata, aResiduals, aMeasErrors, aResErrors, aDownWeights );
        gblRxHist[ix]->fill( (aResiduals[0])*1E3 );
        gblRyHist[ix]->fill( (aResiduals[1])*1E3 );
        gblPxHist[ix]->fill( aResiduals[0]/aResErrors[0] );
        gblPyHist[ix]->fill( aResiduals[1]/aResErrors[1] );
      }

      ax[k] = localPar[1];
      ay[k] = localPar[2];
      ++k;
    }

    for(size_t ix = 0; ix < _nPlanes-1; ++ix) {
      gblKinkXHist[ix]->fill( (ax[ix+1] - ax[ix])*1E3 ); // kink at planes [mrad]
      gblKinkYHist[ix]->fill( (ay[ix+1] - ay[ix])*1E3 ); // kink at planes [mrad]
    }

    // do not pass very bad tracks to mille. Only if the alignment is performed
    if(_performAlignment){
      if(probchi > _chi2cut) {
          traj.milleOut( *milleAlignGBL );
            _nMilleTracks ++;
      }
    }
    numbertracks++;
    
  }
  ntracksperevent->fill( numbertracks );

  // count events:
  ++_iEvt;
  if( isFirstEvent() ) _isFirstEvent = false;
}

//------------------------------------------------------------------------------
void EUTelAlignGBL::end() {
}//end

//------------------------------------------------------------------------------
void EUTelAlignGBL::bookHistos(std::vector<int> const & sensorIDVec) {

  try {
    streamlog_out( MESSAGE2 ) << "Booking histograms..." << endl;

    nTelHits = AIDAProcessor::histogramFactory(this)->createHistogram1D( "nallhit", 201, -0.5, 200.5 );
    nTelHits->setTitle( "Telescope hits/event;telescope hits;events" );

    nUpstreamTriplets = AIDAProcessor::histogramFactory(this)->createHistogram1D( "nupstreamtriplets", 21, -0.5, 20.5 );
    nUpstreamTriplets->setTitle( "Number of Upstream Triplets per Event; Number of Upstream Triplets;Events" );

    nDownstreamTriplets = AIDAProcessor::histogramFactory(this)->createHistogram1D( "ndownstreamtriplets", 21, -0.5, 20.5 );
    nDownstreamTriplets->setTitle( "Number of Downstream Triplets per Event; Number of Downstream Triplets;Events" );

    // GBL:
    AIDAProcessor::tree(this)->mkdir("AllTracks");
    selxHistGBLAlign = AIDAProcessor::histogramFactory(this)->createHistogram1D( "AllTracks/selx", 150, -15, 15 );
    selxHistGBLAlign->setTitle( "extrapolated triplet x at triplet matching point, sel GBL;x [mm];tracks" );

    selyHistGBLAlign = AIDAProcessor::histogramFactory(this)->createHistogram1D( "AllTracks/sely", 100, -10, 10 );
    selyHistGBLAlign->setTitle( "extrapolated triplet y at triplet matching point, sel GBL;y [mm];tracks" );

    selaxHistGBLAlign = AIDAProcessor::histogramFactory(this)->createHistogram1D( "AllTracks/selax", 100, -25, 25 );
    selaxHistGBLAlign->setTitle( "track angle x, sel GBL;x angle [mrad];tracks" );

    selayHistGBLAlign = AIDAProcessor::histogramFactory(this)->createHistogram1D( "AllTracks/selay", 100, -25, 25 );
    selayHistGBLAlign->setTitle( "track angle y, sel GBL;y angle [mrad];tracks" );

    seldxHistGBLAlign = AIDAProcessor::histogramFactory(this)->createHistogram1D( "AllTracks/seldx", 100, -5000, 5000 );
    seldxHistGBLAlign->setTitle( "driplet-triplet x residual at matching point, sel GBL;#Deltax [#mum];tracks" );

    seldyHistGBLAlign = AIDAProcessor::histogramFactory(this)->createHistogram1D( "AllTracks/seldy", 100, -5000, 5000 );
    seldyHistGBLAlign->setTitle( "driplet-triplet y residual at matching point, sel GBL;#Deltay [#mum];tracks" );

    selkxHistGBLAlign = AIDAProcessor::histogramFactory(this)->createHistogram1D( "AllTracks/selkx", 100, -25, 25 );
    selkxHistGBLAlign->setTitle( "triplet-driplet kink x, sel GBL;kink x [mrad];tracks" );

    selkyHistGBLAlign = AIDAProcessor::histogramFactory(this)->createHistogram1D( "AllTracks/selky", 100, -25, 25 );
    selkyHistGBLAlign->setTitle( "triplet-driplet kink y, sel GBL;kink y [mrad];tracks" );

    AIDAProcessor::tree(this)->mkdir("AllTracks/Residuals");
    //Preparing the Triplet interpolation to plane 1 and extrapolation to all planes but 0 & 2 (which are)
    //the ones the triplet was constructed from
    for(size_t ix = 1; ix < sensorIDVec.size(); ++ix) {
      if(ix == 2) continue;
      auto sensorIdString = std::to_string(sensorIDVec[ix]);
      std::string histNameSelX = "AllTracks/Residuals/seldx"+sensorIdString;
      std::string histNameSelY = "AllTracks/Residuals/seldy"+sensorIdString;

      seldxSensorHistGBLAlign.push_back(AIDAProcessor::histogramFactory(this)->createHistogram1D( histNameSelX, 100, -1000, 1000 ));
      seldxSensorHistGBLAlign.back()->setTitle( "triplet resid x at "+sensorIdString+", sel GBL;#Deltax [#mum];tracks" );

      seldySensorHistGBLAlign.push_back(AIDAProcessor::histogramFactory(this)->createHistogram1D( histNameSelY, 100, -1000, 1000 ));
      seldySensorHistGBLAlign.back()->setTitle( "triplet resid y at "+sensorIdString+", sel GBL;#Deltay [#mum];tracks" );
    }

    gblndfHistGBLAlign = AIDAProcessor::histogramFactory(this)->createHistogram1D( "gblndf", 16, -0.5, 15.5 );
    gblndfHistGBLAlign->setTitle( "GBL fit NDF;GBL NDF;tracks" );

    gblchi2HistGBLAlign = AIDAProcessor::histogramFactory(this)->createHistogram1D( "gblchi2", 100, 0, 100 );
    gblchi2HistGBLAlign->setTitle( "GBL fit chi2 / degrees of freedom ;GBL chi2 / Ndf ;tracks" );

    gblprbHistGBLAlign = AIDAProcessor::histogramFactory(this)->createHistogram1D( "gblprb", 100, 0, 1 );
    gblprbHistGBLAlign->setTitle( "GBL fit probability;GBL fit probability;tracks" );

    AIDAProcessor::tree(this)->mkdir("GBLFit");
    AIDAProcessor::tree(this)->mkdir("GBLFit/Angles");
    AIDAProcessor::tree(this)->mkdir("GBLFit/Residuals");
    AIDAProcessor::tree(this)->mkdir("GBLFit/Pulls");
    AIDAProcessor::tree(this)->mkdir("GBLFit/Shifts");
    for(size_t ix = 0; ix < sensorIDVec.size(); ++ix) {
      auto sensorIdString = std::to_string(sensorIDVec[ix]);

      std::string histNameAngleX = "GBLFit/Angles/ax"+sensorIdString;
      std::string histNameAngleY = "GBLFit/Angles/ay"+sensorIdString;

      gblAxHist.push_back(AIDAProcessor::histogramFactory(this)->createHistogram1D( histNameAngleX, 100, -5, 5));
      gblAxHist.back()->setTitle( "GBL angle at plane "+sensorIdString+";x angle at plane "+sensorIdString+" [mrad];tracks" ); 

      gblAyHist.push_back(AIDAProcessor::histogramFactory(this)->createHistogram1D( histNameAngleY, 100, -5, 5));
      gblAyHist.back()->setTitle( "GBL angle at plane "+sensorIdString+";y angle at plane "+sensorIdString+" [mrad];tracks" ); 

      std::string histNameResidX = "GBLFit/Residuals/rx"+sensorIdString; 
      std::string histNameResidY = "GBLFit/Residuals/ry"+sensorIdString; 

      gblRxHist.push_back(AIDAProcessor::histogramFactory(this)->createHistogram1D( histNameResidX, 500, -250, 250));
      gblRxHist.back()->setTitle( "GBL residual at plane "+sensorIdString+";x resid at plane "+sensorIdString+" [#mum];tracks" ); 

      gblRyHist.push_back(AIDAProcessor::histogramFactory(this)->createHistogram1D( histNameResidY, 500, -250, 250));
      gblRyHist.back()->setTitle( "GBL residual at plane "+sensorIdString+";y resid at plane "+sensorIdString+" [#mum];tracks" ); 

      std::string histNamePullX = "GBLFit/Pulls/px"+sensorIdString; 
      std::string histNamePullY = "GBLFit/Pulls/py"+sensorIdString; 

      gblPxHist.push_back(AIDAProcessor::histogramFactory(this)->createHistogram1D( histNamePullX, 100, -5, 5));
      gblPxHist.back()->setTitle( "GBL pull at plane "+sensorIdString+";x pull at plane "+sensorIdString+";tracks" ); 

      gblPyHist.push_back(AIDAProcessor::histogramFactory(this)->createHistogram1D( histNamePullY, 100, -5, 5));
      gblPyHist.back()->setTitle( "GBL pull at plane "+sensorIdString+";y pull at plane "+sensorIdString+";tracks" ); 

      std::string histNameShiftX = "GBLFit/Shifts/dx"+sensorIdString;
      std::string histNameShiftY = "GBLFit/Shifts/dy"+sensorIdString;

      gblDxHist.push_back(AIDAProcessor::histogramFactory(this)->createHistogram1D( histNameShiftX, 500, -250, 250));
      gblDxHist.back()->setTitle( "GBL shift at plane "+sensorIdString+";x shift at plane "+sensorIdString+" [#mum];tracks" );

      gblDyHist.push_back(AIDAProcessor::histogramFactory(this)->createHistogram1D( histNameShiftY, 500, -250, 250));
      gblDyHist.back()->setTitle( "GBL shift at plane "+sensorIdString+";y shift at plane "+sensorIdString+" [#mum];tracks" );
    }
    
    AIDAProcessor::tree(this)->mkdir("GBLFit/Kinks");
    for(size_t ix = 1; ix < sensorIDVec.size(); ++ix) {
      auto sensorIdString = std::to_string(sensorIDVec[ix]);

      std::string histNameKinkX = "GBLFit/Kinks/kx"+sensorIdString;
      std::string histNameKinkY = "GBLFit/Kinks/ky"+sensorIdString;

      gblKinkXHist.push_back(AIDAProcessor::histogramFactory(this)->createHistogram1D( histNameKinkX, 100, -5, 5));
      gblKinkXHist.back()->setTitle( "GBL kink angle at plane "+sensorIdString+";plane "+sensorIdString+" x kink [mrad];tracks" );

      gblKinkYHist.push_back(AIDAProcessor::histogramFactory(this)->createHistogram1D( histNameKinkY, 100, -5, 5));
      gblKinkYHist.back()->setTitle( "GBL kink angle at plane "+sensorIdString+";plane "+sensorIdString+" y kink [mrad];tracks" );
    }

    ntracksperevent = AIDAProcessor::histogramFactory(this)->createHistogram1D( "ntracksperevent", 21, -0.5, 20.5 );
    ntracksperevent->setTitle( "Matched Tracks;Track Matches in an Event;Events" );

  }//try
  catch( lcio::Exception& e ) {

  }
}

