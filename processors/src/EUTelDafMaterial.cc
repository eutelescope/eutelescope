// Author Havard Gjersdal, UiO(haavagj@fys.uio.no)
/*!
 * This is a track fitting processor for the Eutelescope package. 
 *
 * It preforms track finding and fitting on a supplied hit collection.
 *
 * The track finder works by propagating all hits to plane 0, currently assuming straight
 * line fits, then running a cluster finder. Hit clusters above some set value are considered
 * track candidates.
 *
 * This track candidate is then fitted using a implementation of a Deterministic Annealing
 * Filter (DAF), that in short is a Kalman Filter running iteratively over a set of weighted
 * measurements, reweighing the measurements after each fit based on the residuals and a
 * supplied chi2 cut off.
 *
 * This package uses the Eigen library for linear algebra. This package is very quick when
 * compiled properly, but very slow when compiled for debugging. Make sure to compile
 * properly before running productions.
 *
 * Running 'cmake -i' inside the build folder, and then when it asks
 * Variable Name: CMAKE_CXX_FLAGS_RELEASE
 * Description: Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files).                          
 *
 * enter:
 * New Value (Enter to keep current value): -O3 -msse2 -ftree-vectorize -DNDEBUG
 *
 * When it asks
 * Variable Name: CMAKE_BUILD_TYPE
 * enter:
 * New Value (Enter to keep current value): Release
 *
 * If youc cpu supports it, you could try -msse4 or -msse3 aswell.
 */

// built only if GEAR and MARLINUTIL are used
#if defined(USE_GEAR)
// eutelescope includes ".h"
#include "EUTelDafMaterial.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelExceptions.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"
#include "marlin/Exceptions.h"
#include "marlin/AIDAProcessor.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IProfile1D.h>
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

// system includes <>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;


EUTelDafMaterial::EUTelDafMaterial () : EUTelDafBase("EUTelDafMaterial"){
    //Child spesific params and description
  dafParams();
}

void EUTelDafMaterial::dafParams(){
  _description = "This processor preforms track reconstruction. The tracks are used to estimate material thicknesses and measurement resolutions.";
  registerOptionalParameter("MinCol", "Minimum allowed local X of hit. All clusters containing hits outside will be discarded from alignment.", _colMin, std::vector<int>());
  registerOptionalParameter("MaxCol", "Maximum allowed local X of hit. All clusters containing hits outside will be discarded from alignment.", _colMax, std::vector<int>());
  registerOptionalParameter("MinRow", "Minimum allowed local Y of hit. All clusters containing hits outside will be discarded from alignment.", _rowMin, std::vector<int>());
  registerOptionalParameter("MaxRow", "Maximum allowed local Y of hit. All clusters containing hits outside will be discarded from alignment.", _rowMax, std::vector<int>());

  registerOptionalParameter("ResidualXMin", "Min X for residual cuts. All sensors shoule be included, same order as parameter DutPlanes",
			    _resXMin, std::vector<float>());
  registerOptionalParameter("ResidualXMax", "Max X for residual cuts. All sensors shoule be included, same order as parameter DutPlanes",
			    _resXMax, std::vector<float>());
  registerOptionalParameter("ResidualYMin", "Min Y for residual cuts. All sensors shoule be included, same order as parameter DutPlanes",
			    _resYMin, std::vector<float>());
  registerOptionalParameter("ResidualYMax", "Max Y for residual cuts. All sensors shoule be included, same order as parameter DutPlanes",
			    _resYMax, std::vector<float>());

  registerOptionalParameter("RadLengthIndex", "Plane Index for thickness estimator",
			    _radLengthIndex, std::vector<int>());
  registerOptionalParameter("ResXIndex", "Plane Index for sigma X estimator",
			    _resXIndex, std::vector<int>());
  registerOptionalParameter("ResYIndex", "Plane Index for sigma Y estimator",
			    _resYIndex, std::vector<int>());


  registerOptionalParameter("ShiftXIndex", "Plane Index for shift X estimator",
			    _shiftXIndex, std::vector<int>());
  registerOptionalParameter("ShiftYIndex", "Plane Index for shift Y estimator",
			    _shiftYIndex, std::vector<int>());
  registerOptionalParameter("ScaleXIndex", "Plane Index for scale X estimator",
			    _scaleXIndex, std::vector<int>());
  registerOptionalParameter("ScaleYIndex", "Plane Index for scale Y estimator",
			    _scaleYIndex, std::vector<int>());
  registerOptionalParameter("ZRotIndex", "Plane Index for Z Rot estimator",
			    _zRotIndex, std::vector<int>());
  registerOptionalParameter("ZPosIndex", "Plane Index for Z Pos estimator",
			    _zPosIndex, std::vector<int>());
}

void EUTelDafMaterial::dafInit() {
  for(size_t ii = 0; ii < _dutPlanes.size(); ii++){
    int iden = _dutPlanes.at(ii);
    int xMin = _resXMin.size() > ii ? _resXMin.at(ii) : -9999999;
    int xMax = _resXMax.size() > ii ? _resXMax.at(ii) : 9999999;
    int yMin = _resYMin.size() > ii ? _resYMin.at(ii) : -9999999;
    int yMax = _resYMax.size() > ii ? _resYMax.at(ii) : 9999999;
    streamlog_out ( MESSAGE5 ) << "xMin " << xMin << " " << xMax << endl;
    streamlog_out ( MESSAGE5 ) << "yMin " << yMin << " " << yMax << endl;
    _resX[iden] = make_pair(xMin, xMax);
    _resY[iden] = make_pair(yMin, yMax);
  }

  for( size_t ii = 0; ii< _system.planes.size(); ii++){
    //Include DUTs
    if( find(_dutPlanes.begin(), _dutPlanes.end(), _system.planes.at(ii).getSensorID()) != _dutPlanes.end()){ 
      //_system.planes.at(ii).include(); 
    }
  }
  _matest.init(_eBeam, _system.planes.size());
  _angleX = 0;
  _angleY = 0;
}


int EUTelDafMaterial::checkDutResids(daffitter::TrackCandidate<float, 4>& track){
  int nHits(0);
  for( size_t ii = 0; ii < _system.planes.size() ; ii++){
    daffitter::FitPlane<float>& plane = _system.planes.at(ii);
    int iden = plane.getSensorID();
    if( find(_dutPlanes.begin(), _dutPlanes.end(), iden) == _dutPlanes.end()){ continue; }
    std::pair<float, float> resX = _resX[iden];
    std::pair<float, float> resY = _resY[iden];
    track.indexes.at(ii) = -1;

    for(size_t w = 0; w < plane.meas.size(); w++){
      bool includeMeas = true;
      daffitter::Measurement<float>& meas = plane.meas.at(w);
      daffitter::TrackEstimate<float, 4>& estim = track.estimates.at(ii);
      track.weights.at(ii)(w) = 0.0;
      //Resids 
      if( (estim.getX() - meas.getX()) < resX.first) { includeMeas = false; }
      if( (estim.getX() - meas.getX()) > resX.second){ includeMeas = false; }
      if( (estim.getY() - meas.getY()) < resY.first) { includeMeas = false; }
      if( (estim.getY() - meas.getY()) > resY.second){ includeMeas = false; }
      if( includeMeas){
	  //meas.goodRegion() ){
	nHits++;
	track.weights.at(ii)(w) = 1.0;
	track.indexes.at(ii) = w;
	break;
      }
    }
  }
  return(nHits);
}


void EUTelDafMaterial::dafEvent (LCEvent * /*event*/) {
  //Check that there is only one DUT cluster per plane to supress showers
  bool shower(false);
  for(size_t pl = 0; pl < _system.planes.size(); pl++){
    if( find(_dutPlanes.begin(), _dutPlanes.end(), _system.planes.at(pl).getSensorID()) == _dutPlanes.end()){ continue;}
    if( _system.planes.at(pl).getSensorID() == 12) { continue; } 
    if( _system.planes.at(pl).meas.size() > 1){ shower = true; }
  }
  if( shower) { return; }

  //Check found tracks, count how many passes per event
  int matAccept(0);
  for(size_t ii = 0; ii < _system.getNtracks(); ii++ ){
    //run track fitte
    _nCandidates++;
    _system.fitPlanesInfoDaf(_system.tracks.at(ii));
    _system.weightToIndex(_system.tracks.at(ii));
    _system.fitInfoFWBiased(_system.tracks.at(ii));
    _system.getChi2BiasedInfo(_system.tracks.at(ii));
    //Check resids, intime, angles
    if(not checkTrack( _system.tracks.at(ii))) { continue;};
    //int inTimeHits = checkInTime( _system.tracks.at(ii));
    for(size_t pl = 0 ; pl < _system.planes.size(); pl++){
      if( find(_dutPlanes.begin(), _dutPlanes.end(), _system.planes.at(pl).getSensorID()) != _dutPlanes.end()){ 
     	_system.planes.at(pl).include(); 
      }
    }
    _system.weightToIndex(_system.tracks.at(ii));
    int inTimeHits = checkDutResids( _system.tracks.at(ii));
    if( inTimeHits >= _nDutHits) {
      std::vector<Measurement<float> > track;
      for(size_t plane = 0; plane < _system.tracks.at(ii).indexes.size(); plane++ ){
	int hitIndex = _system.tracks.at(ii).indexes.at(plane);
	if( hitIndex >= 0 ){
	  Measurement<float> m = _system.planes.at(plane).meas.at(hitIndex);
	  //if(m.goodRegion()){ track.push_back( Measurement<float>(m.getX(), m.getY(), _system.planes.at(plane).getMeasZ(), true, plane)); }
	  if(m.goodRegion()){ track.push_back( m ); }
	  //track.push_back( m );
	}
	if(plane ==  _system.planes.size()/2){
	  _angleX += _system.tracks.at(ii).estimates.at(plane).getXdz();
	  _angleY += _system.tracks.at(ii).estimates.at(plane).getYdz();
	}
      }
      //Fill plots
      if(_histogramSwitch){ 
	fillDetailPlots( _system.tracks.at(ii) ); 
	fillPlots( _system.tracks.at(ii) ); 
      }
      _matest.addTrack(track);
      _nTracks++;
      
      if(matAccept > 0){
	streamlog_out ( MESSAGE5 ) << "Got more than one track with DUT match!"
				   << endl;
	matAccept++;
      }
    }
    for(size_t pl = 0 ; pl < _system.planes.size(); pl++){
      if( find(_dutPlanes.begin(), _dutPlanes.end(), _system.planes.at(pl).getSensorID()) != _dutPlanes.end()){ 
	_system.planes.at(pl).exclude(); 
      }
    }
  }
}

void EUTelDafMaterial::dafEnd() {
  //Add planes to material estimator fitter
  for( size_t ii = 0; ii< _system.planes.size(); ii++){
    _system.planes.at(ii).include();
    FitPlane<float>& plane = _system.planes.at(ii);
    _matest.system.addPlane(_system.planes.at(ii).getSensorID(),_system.planes.at(ii).getZpos(), 
			    plane.getSigmaX(), plane.getSigmaY(), plane.getScatterThetaSqr(), plane.isExcluded() );
  }

  //These tracks have already been found, use very wide cuts.
  _matest.system.setClusterRadius( 10000.0f);
  _matest.system.setNominalXdz(0.0f);
  _matest.system.setNominalYdz(0.0f);
  _matest.system.setChi2OverNdofCut( 100.0f);
  _matest.system.setDAFChi2Cut( 100.0f);
  _matest.system.init();
      
  //_matest.plot((char*) "/home/haavagj/noalign.root");

  //Initialize EstMat from configuration vectors. These are the alignment parameters
  // for(size_t ii = 0; ii < _xShift.size(); ii++){ _matest.xShift.at(ii) = _xShift.at(ii); }
  // for(size_t ii = 0; ii < _yShift.size(); ii++){ _matest.yShift.at(ii) = _yShift.at(ii); }
  // for(size_t ii = 0; ii < _xScale.size(); ii++){ _matest.xScale.at(ii) = _xScale.at(ii); }
  // for(size_t ii = 0; ii < _yScale.size(); ii++){ _matest.yScale.at(ii) = _yScale.at(ii); }
  // for(size_t ii = 0; ii < _zRot.size(); ii++){ _matest.zRot.at(ii) = _zRot.at(ii); }
  // for(size_t ii = 0; ii < _system.planes.size(); ii++){ _matest.zPos.at(ii) = _system.planes.at(ii).getZpos(); }

  // for(size_t ii = 0; ii < _zPos.size(); ii++){ 
  //   cout << "pl " << ii << " ZPos = " << _zPos.at(ii) << endl;
  //   _matest.zPos.at(ii) = _zPos.at(ii); 
  //   _matest.movePlaneZ(ii, _zPos.at(ii));
  // }

  //Initialize the EstMat free parameters from configuration
  _matest.resXIndex = _resXIndex;
  _matest.resYIndex = _resYIndex;
  _matest.radLengthsIndex = _radLengthIndex;
  _matest.xShiftIndex = _shiftXIndex;
  _matest.yShiftIndex = _shiftYIndex;
  _matest.xScaleIndex = _scaleXIndex;
  _matest.yScaleIndex = _scaleYIndex;
  _matest.zRotIndex = _zRotIndex;
  _matest.zPosIndex = _zPosIndex;

  //Initial state of resolutoins ans radiation lengths
  for(size_t ii = 0; ii < _radLength.size(); ii++){
    _matest.radLengths.at(ii) = _radLength.at(ii); 
    _matest.resX.at(ii) = _system.planes.at(ii).getSigmaX();
    _matest.resY.at(ii) = _system.planes.at(ii).getSigmaY();
  }
  
  streamlog_out ( MESSAGE5 ) << "Starting estimation assuming beam energy of "
			     << _eBeam << std::endl;
  //_matest.plot((char*) "/home/haavagj/preestmat.root");
  
  //Start minimization
  // Minimizer* minimize = new FwBw(_matest); //FWBW
  // Minimizer* minimize = new SDR(true,false,false,_matest); //SDR1
  // Minimizer* minimize = new SDR(false,true,false,_matest); //SDR2
  //Minimizer* minimize = new SDR(true,true,false,_matest); //SDR3, 
  //Minimizer* minimize = new FakeChi2(_matest);
  //Minimizer* minimize = new FakeAbsDev(_matest);
  //_matest.simplexSearch(minimize, 3000, 30);
  
  FwBw* minimize = new FwBw(_matest);
  _matest.quasiNewtonHomeMade(minimize, 400);
  
  //Use this for alignment only. 
  // Minimizer* minimize = new Chi2(_matest); //Alignment
  //_matest.quasiNewtonHomeMade(minimize, 1000);
  
  //Plot the resulting pulls, residuals, chi2s
  //_matest.plot((char*) "/home/haavagj/estmat.root");
}
#endif // USE_GEAR
