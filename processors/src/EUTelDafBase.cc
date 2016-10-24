// Version: $Id$

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
#include "EUTelDafBase.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelExceptions.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelReferenceHit.h"

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
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <IO/LCWriter.h>
#include <UTIL/LCTime.h>
#include <UTIL/CellIDDecoder.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/SimTrackerHitImpl.h>

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

#include <TVector3.h>
#include <Eigen/Geometry> 

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;


EUTelDafBase::EUTelDafBase(std::string name) : marlin::Processor(name) {
  //Universal DAF params

  // input collection
  std::vector<std::string > HitCollectionNameVecExample;
  HitCollectionNameVecExample.push_back("hit");
  registerProcessorParameter( "HitCollectionName", "Names of input hit collections", _hitCollectionName,HitCollectionNameVecExample);
  EVENT::StringVec _mcCollectionExample;
  registerProcessorParameter ("mccollections", "List of hit collections. First one is INPUT collection, every subsequent corresponds to applying alignment collection",
			      _mcCollectionStr, _mcCollectionExample );
  registerProcessorParameter ("clusterfinder", "Name of the clusterfinder which should be used, available are: simpleCluster and combinatorialKF",
				 _clusterFinderName, std::string("simpleCluster"));
  registerOptionalParameter("AlignmentCollectionNames", "Names of alignment collections, should be in same order as application",
			    _alignColNames, std::vector<std::string>());

  // Reference collection
  registerOptionalParameter("ReferenceCollection","reference hit collection name ", _referenceHitCollectionName, static_cast <string> ("referenceHit") );
  registerOptionalParameter("UseReferenceCollection","Do you want the reference hit collection to be used for coordinate transformations?",
			    _useReferenceHitCollection, static_cast< bool   > ( true ));
  
  //Tracker system options
  registerOptionalParameter("MakePlots", "Should plots be made and filled?", _histogramSwitch, static_cast<bool>(false));
  registerProcessorParameter("TelescopePlanes","List of sensor IDs for the telescope planes. These planes are used for the track finder, and track fitter.", _telPlanes ,std::vector<int>());
  registerOptionalParameter("DutPlanes",
			    "List of sensor IDs for the DUT planes. Used to make the decision on whether ro accept the track or not. These planes are not used in track finder, and not in the track fitter unless option 'useDutsInFit' is set."
			    , _dutPlanes ,std::vector<int>());
  registerOptionalParameter("Ebeam", "Beam energy [GeV], used to calculate amount of scatter", _eBeam,  static_cast < float > (120.0));
  registerOptionalParameter("TelResolutionX", "Sigma of telescope resolution in the global X plane,", _telResX,  static_cast < float > (5.3));
  registerOptionalParameter("TelResolutionY", "Sigma of telescope resolution in the global Y plane,", _telResY,  static_cast < float > (5.3));
  registerOptionalParameter("DutResolutionX", "Sigma of telescope resolution in the global X plane,", _dutResX,  static_cast < float > (115.4));
  registerOptionalParameter("DutResolutionY", "Sigma of telescope resolution in the global Y plane,", _dutResY,  static_cast < float > (14.4));

  //Material and resolution
  registerOptionalParameter("RadiationLengths","Radiation lengths of planes, ordered by z-pos..", _radLength, std::vector<float>());
  registerOptionalParameter("ResolutionX","Sigma resolution of planes, ordered by z-pos.", _sigmaX, std::vector<float>());
  registerOptionalParameter("ResolutionY","Sigma resolution of planes, ordered by z-pos.", _sigmaY, std::vector<float>());

  //alignment corrections
  // registerOptionalParameter("XShift","X translation of planes, ordered by z-pos..", _xShift, std::vector<float>());
  // registerOptionalParameter("YShift","Y translation of planes, ordered by z-pos..", _yShift, std::vector<float>());
  // registerOptionalParameter("XScale","X scale of planes, ordered by z-pos..", _xScale, std::vector<float>());
  // registerOptionalParameter("YScale","Y scale of planes, ordered by z-pos..", _yScale, std::vector<float>());
  // registerOptionalParameter("ZRot","Z rotation of planes, ordered by z-pos..", _zRot, std::vector<float>());
  // registerOptionalParameter("ZPos","Z position of planes, ordered by gear z-pos..", _zPos, std::vector<float>());

  //Track finder options
  //registerOptionalParameter("FinderRadius","Track finding: The maximum allowed distance between to hits in the xy plane for inclusion in track candidate", _clusterRadius, static_cast<float>(300.0));
  registerOptionalParameter("FinderRadius","Track finding: The maximum allowed normalized distance between to hits in the xy plane for inclusion in track candidate.", _normalizedRadius, static_cast<float>(300.0));
  registerOptionalParameter("Chi2Cutoff","DAF fitter: The cutoff value for a measurement to be included in the fit.", _chi2cutoff, static_cast<float>(300.0f));
  registerOptionalParameter("RequireNTelPlanes","How many telescope planes do we require to be included in the fit?",_nSkipMax ,static_cast <float> (0.0f));
  registerOptionalParameter("NominalDxdz", "dx/dz assumed by track finder", _nXdz, static_cast<float>(0.0f));
  registerOptionalParameter("NominalDydz", "dy/dz assumed by track finder", _nYdz, static_cast<float>(0.0f));
  registerOptionalParameter("MaxXdxDeviance", "maximum devianve for dx/dz in CKF track finder", _nXdzMaxDeviance, static_cast<float>(0.01f));
  registerOptionalParameter("MaxYdxDeviance", "maximum devianve for dy/dz in CKF track finder", _nYdzMaxDeviance, static_cast<float>(0.01f));

  
  //Track quality parameters
  registerOptionalParameter("MaxChi2OverNdof", "Maximum allowed global chi2/ndof", _maxChi2, static_cast<float> ( 9999.0));
  registerOptionalParameter("NDutHits", "How many DUT hits do we need in order to accept track?", _nDutHits, static_cast <int>(0));
}

bool EUTelDafBase::defineSystemFromData(){
  //Find three measurements per plane, use these three to define the plane as a point and a normal vector
  bool gotIt = true;
  for(size_t plane = 0; plane < _system.planes.size(); plane++){
    daffitter::FitPlane<float>& pl = _system.planes.at(plane);
    bool gotPlane = false;

    //If plane is not DUT and not Tel, it does not have hits, and can not be initialized
    if( find(_telPlanes.begin(), _telPlanes.end(), pl.getSensorID()) == _telPlanes.end() and
	find(_dutPlanes.begin(), _dutPlanes.end(), pl.getSensorID()) == _dutPlanes.end()){
      continue;
    }
    
    for(size_t meas = 0; meas < pl.meas.size(); meas++){
      if( _nRef.at(plane) > 2){ gotPlane = true; continue; }
      if( _nRef.at(plane) == 0 ){
	pl.setRef0( Eigen::Matrix<float, 3, 1>(pl.meas.at(meas).getX(), pl.meas.at(meas).getY(), pl.meas.at(meas).getZ()));
	_nRef.at(plane)++;
	gotPlane = false;
	continue;
      }

      if( fabs(pl.meas.at(meas).getX() - pl.getRef0()(0) ) < 500) { continue; }
      if( fabs(pl.meas.at(meas).getY() - pl.getRef0()(1) ) < 500) { continue; }

      if( _nRef.at(plane) == 1 ){
	pl.setRef1( Eigen::Matrix<float, 3, 1>(pl.meas.at(meas).getX(), pl.meas.at(meas).getY(), pl.meas.at(meas).getZ()));
	_nRef.at(plane)++;
	gotPlane = false;
	continue;
      }
      if( fabs(pl.meas.at(meas).getX() - pl.getRef1()(0) ) < 500) {  continue; }
      if( fabs(pl.meas.at(meas).getY() - pl.getRef1()(1) ) < 500) {  continue; }
      if( _nRef.at(plane) == 2 ){
	pl.setRef2( Eigen::Matrix<float, 3, 1>(pl.meas.at(meas).getX(), pl.meas.at(meas).getY(), pl.meas.at(meas).getZ()));
	_nRef.at(plane)++;
	getPlaneNorm(pl);
	streamlog_out ( MESSAGE5 ) << "Initialized plane " << pl.getSensorID() << endl;
	gotPlane = true;
	continue;
      }
    }
    if(not gotPlane) { gotIt = false;}
  }
  return(gotIt);
}

void EUTelDafBase::gearRotate(size_t index, size_t gearIndex){
  daffitter::FitPlane<float>& pl = _system.planes.at(index);

  double gRotation[3] = { 0., 0., 0.};
  gRotation[0] = _siPlanesLayerLayout->getLayerRotationXY(gearIndex); // Euler alpha
  gRotation[1] = _siPlanesLayerLayout->getLayerRotationZX(gearIndex); // Euler beta
  gRotation[2] = _siPlanesLayerLayout->getLayerRotationZY(gearIndex); // Euler gamma
  // transform into radians
  gRotation[0] =  gRotation[0]*3.1415926/180.; 
  gRotation[1] =  gRotation[1]*3.1415926/180.; 
  gRotation[2] =  gRotation[2]*3.1415926/180.; 
  
  //Reference points define plane
  //Transform ref, cloning hitmaker logic
  TVector3 ref0 = TVector3(0.0, 0.0, 0.0);
  TVector3 ref1 = TVector3(10.0, 0.0, 0.0);
  TVector3 ref2 = TVector3(0.0, 10.0, 0.0);

  double zZero        = _siPlanesLayerLayout->getSensitivePositionZ(gearIndex); // mm
  //double zThickness   = _siPlanesLayerLayout->getSensitiveThickness(gearIndex); // mm

  double nomZ = zZero;// + 0.5 * zThickness;

  if( TMath::Abs( gRotation[0] ) > 1e-6 ){
    ref0.RotateZ( gRotation[0] );
    ref1.RotateZ( gRotation[0] );
    ref2.RotateZ( gRotation[0] );
  }
  if( TMath::Abs( gRotation[1] )> 1e-6 ) {
    ref0.RotateY( gRotation[1] );
    ref1.RotateY( gRotation[1] );
    ref2.RotateY( gRotation[1] );
  }
  if( TMath::Abs( gRotation[2] ) > 1e-6 ){
    ref0.RotateX( gRotation[2] );
    ref1.RotateX( gRotation[2] );
    ref2.RotateX( gRotation[2] );
  }

  pl.setRef0( Eigen::Vector3f( ref0.X() * 1000.0f, ref0.Y() * 1000.0f, (ref0.Z() + nomZ) * 1000.0f ));
  pl.setRef1( Eigen::Vector3f( ref1.X() * 1000.0f, ref1.Y() * 1000.0f, (ref1.Z() + nomZ) * 1000.0f ));
  pl.setRef2( Eigen::Vector3f( ref2.X() * 1000.0f, ref2.Y() * 1000.0f, (ref2.Z() + nomZ) * 1000.0f ));

  //Tracks are propagated to glob xy plane => Errors are in glob xy plane. scales like cosine
  //Errors not corrected for xy rotation
  pl.scaleErrors( std::fabs(std::cos(gRotation[1])), std::fabs(std::cos(gRotation[0])));
  getPlaneNorm(pl);
} 

Eigen::Vector3f EUTelDafBase::applyAlignment(EUTelAlignmentConstant* alignment, Eigen::Vector3f point){
  Eigen::Vector3f outpoint;
  double alpha = alignment->getAlpha();
  double beta  = alignment->getBeta();  
  double gamma = alignment->getGamma();
  double z_sensor = point(2);

  // sync with updated sign convention of the rotation angles:
  //
  outpoint(0) =                point(0) + (-1)*gamma * point(1) +      beta  * (point(2) - z_sensor) ;
  outpoint(1) =        gamma * point(0) +              point(1) + (-1)*alpha * (point(2) - z_sensor) ;
  outpoint(2) = (-1) * beta  * point(0) +      alpha * point(1) +              (point(2) - z_sensor) ;

  // second the shift
  outpoint(0) -= 1000.0f * alignment->getXOffset();
  outpoint(1) -= 1000.0f * alignment->getYOffset();
  outpoint(2) -= 1000.0f * alignment->getZOffset();
           
  outpoint(2) += z_sensor ;

  return(outpoint);
}

void EUTelDafBase::alignRotate(std::string collectionName, LCEvent* event) {
  LCCollectionVec * alignmentCollectionVec;
  try {
    alignmentCollectionVec     = dynamic_cast < LCCollectionVec * > (event->getCollection(collectionName));
  } catch (DataNotAvailableException& e) {
    throw runtime_error("Unable to open alignment collection " + collectionName);
  }
  for( size_t plane = 0; plane < _system.planes.size() ; plane++){
    daffitter::FitPlane<float>& pl = _system.planes.at(plane);
    int iden = pl.getSensorID();
    for ( size_t ii = 0; ii < alignmentCollectionVec->size(); ++ii ) {
      try{
	EUTelAlignmentConstant * alignment = static_cast< EUTelAlignmentConstant * >
	  ( alignmentCollectionVec->getElementAt( ii ) );
	if( alignment->getSensorID() != iden) { continue; }
	pl.setRef0( applyAlignment(alignment, pl.getRef0()) );
	pl.setRef1( applyAlignment(alignment, pl.getRef1()) );
	pl.setRef2( applyAlignment(alignment, pl.getRef2()) );
	//Errors not corrected for xy rotation
	getPlaneNorm(pl);
	pl.scaleErrors(alignment->getAlpha() + 1.0f, alignment->getBeta() + 1.0f);
	pl.print();
      }
      catch(...){
        streamlog_out(WARNING) << "Could not find sensor in " <<  collectionName.c_str() << " at " << ii << endl;
      }
    }
  }
}
void EUTelDafBase::getPlaneNorm(daffitter::FitPlane<float>& pl){
  Eigen::Vector3f l1 = pl.getRef1() - pl.getRef0();
  Eigen::Vector3f l2 = pl.getRef2() - pl.getRef0();
  //Calculate plane normal vector from ref points
  pl.setPlaneNorm( l2.cross(l1));
}

void EUTelDafBase::init() {
  printParameters ();

  _iRun = 0; _iEvt = 0; _nTracks = 0; _nCandidates =0;
  n_passedNdof =0; n_passedChi2OverNdof = 0; n_passedIsnan = 0;
  n_failedNdof =0; n_failedChi2OverNdof = 0; n_failedIsnan = 0;
  _initializedSystem = false;

  //Set-uo the cluster finder
  if(_clusterFinderName == "simpleCluster") {
	_trackFinderType = simpleCluster;
  } else if(_clusterFinderName == "combinatorialKF") {
	_trackFinderType = combinatorialKF;
  } else {
	throw std::runtime_error("DAF-Fitter: Choosen cluster finder: "+_clusterFinderName+"does not exist");
  }

  //Geometry description
  _siPlanesParameters  = const_cast<gear::SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

  //Use map to sort planes by z
  _zSort.clear();
  for(int plane = 0; plane < _siPlanesLayerLayout->getNLayers(); plane++){
    _zSort[ _siPlanesLayerLayout->getLayerPositionZ(plane) * 1000.0 ] = plane;
    _indexIDMap[ _siPlanesLayerLayout->getID( plane )] = plane;
  }

  //Add Planes to tracker system,
  map<float, int>::iterator zit = _zSort.begin();
  size_t index(0), nActive(0);
  for( ; zit != _zSort.end(); index++, zit++){
    _nRef.push_back(3);
    int sensorID = _siPlanesLayerLayout->getID( (*zit).second );
    //Read sensitive as 0, in case the two are different
    float zPos  = _siPlanesLayerLayout->getSensitivePositionZ( (*zit).second )* 1000.0;
    //+ 0.5 * 1000.0 *  _siPlanesLayerLayout->getSensitiveThickness( (*zit).second) ; // Do not move plane to center of plane, use front.
    //Figure out what kind of plane we are dealing with
    float errX(0.0f), errY(0.0f);
    bool excluded = true;

    // Get scatter using x / x0
    float radLength = _siPlanesLayerLayout->getLayerThickness( (*zit).second ) /  _siPlanesLayerLayout->getLayerRadLength( (*zit).second );
    radLength += _siPlanesLayerLayout->getSensitiveThickness( (*zit).second ) /  _siPlanesLayerLayout->getSensitiveRadLength( (*zit).second );

    streamlog_out ( MESSAGE5 ) << " zPos:      " << zPos << " " << radLength ;
    streamlog_out ( MESSAGE5 ) << " sen thick: " << _siPlanesLayerLayout->getSensitiveThickness( (*zit).second ) ;
    streamlog_out ( MESSAGE5 ) << " sens rad:  " << _siPlanesLayerLayout->getSensitiveRadLength( (*zit).second ) << endl;
    if( _radLength.size() > index){
      radLength = _radLength.at(index);
    } else {
      _radLength.push_back(radLength);
    }

    float scatter = getScatterThetaVar( radLength );
    
    streamlog_out ( MESSAGE5 ) << " radlength: "<< radLength << endl;
    streamlog_out ( MESSAGE5 ) << " scatter: "<< scatter << endl;
    
    //Is current plane a telescope plane?
    if( find(_telPlanes.begin(), _telPlanes.end(), sensorID) != _telPlanes.end()){
      _nRef.at(index) = 0;
      errX = _telResX; errY = _telResY;
      excluded = false;
    }
    //Is current plane a DUT plane?
    if( find(_dutPlanes.begin(), _dutPlanes.end(), sensorID) != _dutPlanes.end()){
      _nRef.at(index) = 0;
      errX = _dutResX; errY = _dutResY;
      scatter = getScatterThetaVar( radLength );
    }
    //If plane is neither Tel nor Dut, all we need is the zpos and scatter amount.

    //Add plane to tracker system
    if(not excluded){ nActive++;}
    _system.addPlane(sensorID, zPos , errX, errY, scatter, excluded);
    gearRotate(index, (*zit).second);
  }
 
  //Prepare track finder
  switch( _trackFinderType ) {
	case combinatorialKF: _system.setCKFChi2Cut(_normalizedRadius*_normalizedRadius); break;
	case simpleCluster: _system.setClusterRadius(_normalizedRadius); break;
  }
  _system.setNominalXdz(_nXdz); //What is the tangent angle of the beam? (Probably zero)
  _system.setNominalYdz(_nYdz);
  _system.setXdzMaxDeviance(_nXdzMaxDeviance); //How far og the nominal angle can the first two measurements be?
  _system.setYdzMaxDeviance(_nYdzMaxDeviance);

  //Prepare and preallocate memory for track fitter
  _system.setChi2OverNdofCut(_maxChi2);
  _system.setDAFChi2Cut(_chi2cutoff);
  _system.init();

  //Fuzzy assignment by DAF might make a plane only partially included, This means ndof is
  //not a integer. Everything above (ndof - 0.5) is assumed to include at least ndof degrees of
  //freedom. 
  _ndofMin = -4 + _nSkipMax * 2 - 0.5;
  streamlog_out ( MESSAGE5 ) << "NDOF min is " << _ndofMin << endl;
  if( _ndofMin < 0.5) {
    streamlog_out ( ERROR5 ) << "Too few active planes(" << nActive << ") when " << _nSkipMax << " planes can be skipped." 
			     << "Please check your configuration." << endl;
    exit(1);
  }
  dafInit();
    
  if(_histogramSwitch) {
    bookHistos(); 
    bookDetailedHistos();
  }

  //Define region for edge masking
  for(size_t ii = 0; ii < _dutPlanes.size(); ii++){
    int iden = _dutPlanes.at(ii);
    int xMin = _colMin.size() > ii ? _colMin.at(ii) : -9999999;
    int xMax = _colMax.size() > ii ? _colMax.at(ii) : 9999999;
    int yMin = _rowMin.size() > ii ? _rowMin.at(ii) : -9999999;
    int yMax = _rowMax.size() > ii ? _rowMax.at(ii) : 9999999;
    _colMinMax[iden] = make_pair(xMin, xMax);
    _rowMinMax[iden] = make_pair(yMin, yMax);
  }
}

void EUTelDafBase::processRunHeader (LCRunHeader * rdr) {
  std::unique_ptr<EUTelRunHeaderImpl> header = std::make_unique<EUTelRunHeaderImpl>(rdr);
  header->addProcessor(type());
  ++_iRun;
}


float EUTelDafBase::getScatterThetaVar( float radLength ){
  //From pdg live
  float scatterTheta = 0.0136f/_eBeam * sqrt( radLength ) *  (1.0f + 0.038f * std::log(radLength) );
  return(scatterTheta * scatterTheta);
}

size_t EUTelDafBase::getPlaneIndex(float zPos){
  //Get plane index from z-position of hit
  map<float,int>::iterator it = _zSort.begin();
  size_t index(0);
  bool foundIt(false);
  for(;it != _zSort.end(); index++, it++){
    if( fabs((*it).first - zPos) < 30000.0){ 
      foundIt = true; break;
    }
  }
  if(not foundIt){ 
    streamlog_out ( ERROR5 ) << "Found hit at z=" << zPos << " , not able to assign to any plane!" << endl; 
    return(-1);
  }
  return(index);
}

void EUTelDafBase::readHitCollection(LCEvent* event){
  //Dump LCIO hit collection to tracker system
  //Extract hits from collection, add to tracker system
  streamlog_out ( DEBUG5 ) << " readHitCollection: " << _hitCollectionName.size() << " collections to read " << endl;
  for(size_t i =0;i < _hitCollectionName.size();i++) {
    streamlog_out ( DEBUG5 ) << " hit collection name: " << _hitCollectionName[i] << " found for event " << event->getEventNumber();
    try {
      _hitCollection = event->getCollection(_hitCollectionName[i]);
    } catch (DataNotAvailableException& e){
      streamlog_out ( WARNING2 ) << "No input collection " << _hitCollectionName[i] << " found for event " << event->getEventNumber()
				 << " in run " << event->getRunNumber() << endl;
      throw SkipEventException( this );
    }
    //Add all hits in collection to corresponding plane
    streamlog_out ( DEBUG5 ) << " hit collection size : " << _hitCollection->getNumberOfElements() << endl;
    
    for ( int iHit = 0; iHit < _hitCollection->getNumberOfElements(); iHit++ ) {
      TrackerHitImpl* hit = static_cast<TrackerHitImpl*> ( _hitCollection->getElementAt(iHit) );
      double pos[3]  = {0.,0.,0.};
      bool region    = true;
      int planeIndex = -1;
      
      if( _mcCollectionStr.size() > 0 ){
	_mcCollection = dynamic_cast < LCCollectionVec * > (event->getCollection(  _mcCollectionStr[i] ));
	SimTrackerHitImpl* simhit = 0;
	if(_mcCollection != 0 ) simhit = static_cast<SimTrackerHitImpl*> ( _mcCollection->getElementAt(iHit) );
	if(simhit != 0 ){
	  UTIL::CellIDDecoder<SimTrackerHitImpl> simHitDecoder (_mcCollection);
	  const double * simpos = simhit->getPosition();
          pos[0]=simpos[0];
          pos[1]=simpos[1];
          pos[2]=simpos[2];
	  int planeID = simHitDecoder(simhit)["sensorID"];
	  planeIndex = _indexIDMap[planeID];
	}
	streamlog_out ( DEBUG5 ) << " SIM: simhit="<< ( simhit != 0 ) <<" add point [" << planeIndex << "] "<< 
	  static_cast< float >(pos[0]) * 1000.0f << " " << static_cast< float >(pos[1]) * 1000.0f << " " <<  static_cast< float >(pos[2]) * 1000.0f << endl;
      } else if (hit != 0 ){
	const double * hitpos = hit->getPosition();
	pos[0]=hitpos[0];
	pos[1]=hitpos[1];
	pos[2]=hitpos[2];
	UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder ( EUTELESCOPE::HITENCODING );
	int planeID  = hitDecoder(hit)["sensorID"];
	planeIndex = _indexIDMap[planeID];
	streamlog_out ( DEBUG5 ) << " REAL: add point [" << planeIndex << "] "<< 
	  static_cast< float >(pos[0]) * 1000.0f << " " << static_cast< float >(pos[1]) * 1000.0f << " " <<  static_cast< float >(pos[2]) * 1000.0f << endl;
      }

      if(planeIndex >=0 ){ 
	streamlog_out ( DEBUG5 ) << " add point [" << planeIndex << "] "<< 
	  static_cast< float >(pos[0]) * 1000.0f << " " << static_cast< float >(pos[1]) * 1000.0f << " " <<  static_cast< float >(pos[2]) * 1000.0f << endl;
	_system.addMeasurement( planeIndex,
				static_cast< float >( pos[0] ) * 1000.0f,
				static_cast< float >( pos[1] ) * 1000.0f,
				static_cast< float >( pos[2] ) * 1000.0f,  region, iHit );
      }
    }
  }
}

int EUTelDafBase::checkInTime(daffitter::TrackCandidate<float, 4>& track){
  size_t nMatches(0);
  for( size_t ii = 0; ii < _system.planes.size() ; ii++){
    daffitter::FitPlane<float>& plane = _system.planes.at(ii);
    int sensorID = plane.getSensorID();
    //Check if any DUT plane is "In time" with the 
    if( find(_dutPlanes.begin(), _dutPlanes.end(), sensorID) == _dutPlanes.end()){
      continue;
    }
    //In timeness can be checked by seeing if the plane has assigned weight
    for(size_t w = 0; w < plane.meas.size(); w++){
      if( track.weights.at(ii)(w) < 0.5f ) {  continue; }
      nMatches++;
      break;
    }
  }
  return(nMatches);
}

void EUTelDafBase::processEvent(LCEvent * event){
  //Called once per event, read data, fit, save
  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);
  if(event->getEventNumber() % 1000 == 0){
    streamlog_out ( MESSAGE ) << "Accepted " << _nTracks <<" tracks at event " << event->getEventNumber() << endl;
    if(not _initializedSystem){
      streamlog_out ( MESSAGE ) << "System not initialized " << event->getEventNumber() << endl;
    }
  }
  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG2 ) << "EORE found: nothing else to do." << endl;
    return;
  }
  if( isFirstEvent() ){
    //Resolution needs to be updated w.r.t. alignment and GEAR rotations.
    for(size_t ii = 0; ii < _alignColNames.size(); ii++){
      alignRotate(_alignColNames.at(ii), event);
    }
    
    //Use user defined resolutions if supplied.
    if(_sigmaX.size() != _sigmaY.size()){
      streamlog_out ( MESSAGE )<< "Differing lengths of resolution X and Y, only filling shortest vector. Check config." << endl;
    }
    size_t nResolutions = _sigmaX.size();
    if(_sigmaY.size() < nResolutions){ nResolutions = _sigmaY.size(); }
    
    if(nResolutions > _system.planes.size()){
      streamlog_out ( MESSAGE ) << "More resolutions than planes, check config." << endl;
      nResolutions = _system.planes.size();
    }
    for(size_t ii = 0; ii < nResolutions; ii++){
      _system.planes.at(ii).setSigmas( _sigmaX.at(ii), _sigmaY.at(ii));
    }
  }
  
  if ( _useReferenceHitCollection ){
    try {
      _referenceHitVec = dynamic_cast < LCCollectionVec * > (event->getCollection( _referenceHitCollectionName));
    } catch (...){
      streamlog_out ( ERROR5 ) <<  "Reference Hit Collection " << _referenceHitCollectionName.c_str()
			       << " could not be retrieved for event " << event->getEventNumber()
			       << "! Please check your steering files! " << endl;
    }
  }

  //get MC collections if exists
  if( _mcCollectionStr.size() > 0 ) {
    for(unsigned int i=0; i < _mcCollectionStr.size(); i++)  {
      _mcCollection = dynamic_cast < LCCollectionVec * > (event->getCollection(  _mcCollectionStr[i] ));
      streamlog_out( DEBUG5 ) << "Collection " << i << " " << _mcCollectionStr[i].c_str() << " at " << _mcCollection << endl;
    }
  } else {
    _mcCollection = NULL;
  }

  //Prepare tracker system for new data, new tracks
  _system.clear();
  
  //Dump hit collection to collection sorted by plane
  readHitCollection(event);
  
  if( not _initializedSystem ){
    //If the system is not initialized, try to finish initialization from event data.
    _initializedSystem = defineSystemFromData();
    //If initialization is not done, we need data from more events.
    if(not _initializedSystem) { return; }
    streamlog_out(MESSAGE1) << "Initialized system at event " << event->getEventNumber() << std::endl;
    for(size_t ii = 0; ii < _system.planes.size(); ii++){
      _system.planes.at(ii).print();
    }
  }
  
  //Run track finder
  switch( _trackFinderType ) {
	case combinatorialKF: _system.combinatorialKF(); break;
	case simpleCluster: _system.clusterTracker(); break;
  }
 
  //Child specific actions
  dafEvent(event);

  streamlog_out(MESSAGE1) << " dafEvent is OVER " <<std::endl;

  if(event->getEventNumber() % 1000 == 0){
    streamlog_out ( MESSAGE5 ) << "Accepted " << _nTracks <<" tracks at event " << event->getEventNumber() << endl;
  }
}

bool EUTelDafBase::checkTrack(daffitter::TrackCandidate<float,4>& track){
  //Check the track quality
  if( track.ndof < _ndofMin) {n_failedNdof++; return(false); }
  n_passedNdof++;
  if( (track.chi2 / track.ndof) > _maxChi2 ) {n_failedChi2OverNdof++; return(false); }
  n_passedChi2OverNdof++;
  if( isnan(track.ndof)) {n_failedIsnan++; return(false); }
  n_passedIsnan++;

  return(true);
}

void EUTelDafBase::fillPlots(daffitter::TrackCandidate<float,4>& track){
  _aidaHistoMap["chi2"]->fill( track.chi2);
  _aidaHistoMap["logchi2"]->fill( std::log10(track.chi2));
  _aidaHistoMap["ndof"]->fill( track.ndof);
  _aidaHistoMap["chi2overndof"]->fill( track.chi2 / track.ndof);
  //Fill plots per plane
  for( size_t ii = 0; ii < _system.planes.size() ; ii++){
    daffitter::FitPlane<float>& plane = _system.planes.at(ii);
    char iden[4];
    sprintf(iden, "%d", plane.getSensorID());
    string bname = static_cast< string >("pl") + iden + "_";
    //Plot resids, angles for all hits with > 50% includion in track.
    //This should be one measurement per track

    daffitter::TrackEstimate<float,4>& estim = track.estimates.at(ii);
    for(size_t w = 0; w < plane.meas.size(); w++){
      if( track.weights.at(ii)(w) < 0.5f ) {  continue; }
      daffitter::Measurement<float>& meas = plane.meas.at(w);
      //Resids 
      _aidaHistoMap[bname + "residualX"]->fill( (estim.getX() - meas.getX())*1e-3 );
      _aidaHistoMap[bname + "residualY"]->fill( (estim.getY() - meas.getY())*1e-3 );

      //Resids 
      _aidaHistoMapProf1D[bname + "residualdXvsX"]->fill(estim.getX(), estim.getX() - meas.getX() );
      _aidaHistoMapProf1D[bname + "residualdYvsX"]->fill(estim.getX(), estim.getY() - meas.getY() );
      _aidaHistoMapProf1D[bname + "residualdXvsY"]->fill(estim.getY(), estim.getX() - meas.getX() );
      _aidaHistoMapProf1D[bname + "residualdYvsY"]->fill(estim.getY(), estim.getY() - meas.getY() );
      _aidaHistoMapProf1D[bname + "residualdZvsX"]->fill(estim.getX(), plane.getMeasZ() - meas.getZ()  );
      _aidaHistoMapProf1D[bname + "residualdZvsY"]->fill(estim.getY(), plane.getMeasZ() - meas.getZ()  );
      _aidaHistoMap2D[bname + "residualmeasZvsmeasX"]->fill(  meas.getZ()/1000., meas.getX()  );
      _aidaHistoMap2D[bname + "residualmeasZvsmeasY"]->fill(  meas.getZ()/1000., meas.getY()  );
      _aidaHistoMap2D[bname + "residualfitZvsmeasX"]->fill( plane.getMeasZ()/1000., meas.getX() );
      _aidaHistoMap2D[bname + "residualfitZvsmeasY"]->fill( plane.getMeasZ()/1000., meas.getY() );
 
      _aidaHistoMap2D[ "AllResidmeasZvsmeasX"]->fill(  meas.getZ()/1000., meas.getX()  );
      _aidaHistoMap2D[ "AllResidmeasZvsmeasY"]->fill(  meas.getZ()/1000., meas.getY()  );
      _aidaHistoMap2D[ "AllResidfitZvsmeasX"]->fill( plane.getMeasZ()/1000., meas.getX() );
      _aidaHistoMap2D[ "AllResidfitZvsmeasY"]->fill( plane.getMeasZ()/1000., meas.getY() );
      //Angles
      _aidaHistoMap[bname + "dxdz"]->fill( estim.getXdz() );
      _aidaHistoMap[bname + "dydz"]->fill( estim.getYdz() );
      if( ii != 4) { continue; }
      _aidaZvHitX->fill(estim.getX(), meas.getZ() - plane.getZpos());
      _aidaZvFitX->fill(estim.getX(), (plane.getMeasZ() - plane.getZpos()) - (meas.getZ() - plane.getZpos()));
      _aidaZvHitY->fill(estim.getY(), meas.getZ() - plane.getZpos());
      _aidaZvFitY->fill(estim.getY(), (plane.getMeasZ() - plane.getZpos()) - (meas.getZ() - plane.getZpos()));
    }
  }
}

void EUTelDafBase::fillDetailPlots(daffitter::TrackCandidate<float,4>& track){
  for( size_t ii = 0; ii < _system.planes.size() ; ii++){
    daffitter::FitPlane<float>& plane = _system.planes.at(ii);

    daffitter::TrackEstimate<float,4>& estim = track.estimates.at(ii);

    char iden[4];
    sprintf(iden, "%d", plane.getSensorID());
    string bname = static_cast< string >("pl") + iden + "_";

    //Plot resids, angles for all hits with > 50% includion in track.
    //This should be one measurement per track
    for(size_t w = 0; w < plane.meas.size(); w++){
      daffitter::Measurement<float>& meas = plane.meas.at(w);
      if( track.weights.at(ii)(w) < 0.5f ) {  continue; }
      //Resids 
      float resX = ( estim.getX() - meas.getX() );
      resX *= resX;
      resX /= plane.getSigmaX() *  plane.getSigmaX() + estim.cov(0,0);
      float resY = ( estim.getY() - meas.getY() );
      resY *= resY;
      resY /= plane.getSigmaY() *  plane.getSigmaY() + estim.cov(1,1);
      _aidaHistoMap[bname + "hitChi2"]->fill( resX + resY );
      
      _aidaHistoMap[bname + "sigmaX"]->fill( sqrt(estim.cov(0,0)) );
      _aidaHistoMap[bname + "sigmaY"]->fill( sqrt(estim.cov(1,1)) );
      
      float pullX =  ( estim.getX() - meas.getX() ) / sqrt(plane.getSigmaX() * plane.getSigmaX() + estim.cov(0,0));
      float pullY =  ( estim.getY() - meas.getY() ) / sqrt(plane.getSigmaY() * plane.getSigmaY() + estim.cov(1,1));
      _aidaHistoMap[bname + "pullX"]->fill( pullX );
      _aidaHistoMap[bname + "pullY"]->fill( pullY );
    }
  }
}

void EUTelDafBase::bookHistos(){

  int maxNdof = -4 + _system.planes.size() * 2 + 1;
  _aidaHistoMap["chi2"] = AIDAProcessor::histogramFactory(this)->createHistogram1D("chi2", 100, 0, maxNdof * _maxChi2);
  _aidaHistoMap["logchi2"] = AIDAProcessor::histogramFactory(this)->createHistogram1D("logchi2", 100, 0, std::log10(maxNdof * _maxChi2));
  if( _aidaHistoMap["chi2"] == NULL){
    streamlog_out ( ERROR2 ) << "Problem with histo booking. Check paths!" << std::endl;
    _histogramSwitch = false;
    return;
  }
  _aidaHistoMap["ndof"] = AIDAProcessor::histogramFactory(this)->createHistogram1D("ndof", maxNdof * 10, 0, maxNdof);
  _aidaHistoMap["chi2overndof"] = AIDAProcessor::histogramFactory(this)->createHistogram1D("Chi2OverNdof", maxNdof * 10, 0, _maxChi2);

  _aidaZvFitX = AIDAProcessor::histogramFactory(this)->createHistogram2D("ZvHitX", 20,  -5000.0,  5000.0, 20,   -100.0,   100.0);
  _aidaZvHitX = AIDAProcessor::histogramFactory(this)->createHistogram2D("ZvFitX", 20, -10000.0, 10000.0, 20, -10000.0, 10000.0);
  _aidaZvFitY = AIDAProcessor::histogramFactory(this)->createHistogram2D("ZvHitY", 20, -10000.0, 10000.0, 20,   -100.0,   100.0);
  _aidaZvHitY = AIDAProcessor::histogramFactory(this)->createHistogram2D("ZvFitY", 20, -10000.0, 10000.0, 20, -10000.0, 10000.0);

  _aidaHistoMap2D["AllResidmeasZvsmeasX"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( "AllResidmeasZvsmeasX",14 ,-80., 60., 20 ,-10000., 10000.);
  _aidaHistoMap2D["AllResidmeasZvsmeasY"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( "AllResidmeasZvsmeasY",14 ,-80., 60., 20 ,-10000., 10000.);
  _aidaHistoMap2D["AllResidfitZvsmeasX"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( "AllResidfitZvsmeasX",14 ,-80., 60., 20 ,-10000., 10000.);
  _aidaHistoMap2D["AllResidfitZvsmeasY"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( "AllResidfitZvsmeasY",14 ,-80., 60., 20 ,-10000., 10000.);


  for( size_t ii = 0; ii < _system.planes.size() ; ii++)
    {
      daffitter::FitPlane<float>& plane = _system.planes.at(ii);
      char iden[4];
      sprintf(iden, "%d", plane.getSensorID());
      string bname = static_cast< string >("pl") + iden + "_";
      //Resids

      double residminX = -0.3;
      double residmaxX =  0.3;
 
      _aidaHistoMap[bname + "mcresidualX"] =  AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "mcresidualX", 200,  residminX, residmaxX );
      _aidaHistoMap[bname + "mcresidualY"] =  AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "mcresidualY", 200,  residminX, residmaxX );

      _aidaHistoMap[bname + "residualX"] =  AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "residualX",600, residminX, residmaxX);
      _aidaHistoMap[bname + "residualY"] =  AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "residualY",600, residminX, residmaxX);
      //Resids 2D // profiles
      _aidaHistoMapProf1D[bname+"residualdXvsX"]= AIDAProcessor::histogramFactory(this)->createProfile1D(bname+"dXvsX", 200, -10000., 10000.,   residminX, residmaxX );
      _aidaHistoMapProf1D[bname+"residualdYvsX"]= AIDAProcessor::histogramFactory(this)->createProfile1D(bname+"dXvsY", 200, -10000., 10000.,   residminX, residmaxX );
      _aidaHistoMapProf1D[bname+"residualdXvsY"]= AIDAProcessor::histogramFactory(this)->createProfile1D(bname+"dYvsX", 200, -10000., 10000.,   residminX, residmaxX );
      _aidaHistoMapProf1D[bname+"residualdYvsY"]= AIDAProcessor::histogramFactory(this)->createProfile1D(bname+"dYvsY", 200, -10000., 10000.,   residminX, residmaxX );
      _aidaHistoMapProf1D[bname+"residualdZvsX"]= AIDAProcessor::histogramFactory(this)->createProfile1D(bname+"dZvsX", 200, -10000., 10000.,  -100., 100. );
      _aidaHistoMapProf1D[bname+"residualdZvsY"]= AIDAProcessor::histogramFactory(this)->createProfile1D(bname+"dZvsY", 200, -10000., 10000.,  -100., 100. );
    
      // residuals
      _aidaHistoMap2D[bname + "residualmeasZvsmeasX"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( bname + "residualmeasZvsmeasX",20 ,-1000., 1000., 20 ,-10000., 10000.);
      _aidaHistoMap2D[bname + "residualmeasZvsmeasY"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( bname + "residualmeasZvsmeasY",20 ,-1000., 1000., 20 ,-10000., 10000.);
      _aidaHistoMap2D[bname + "residualfitZvsmeasX"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( bname + "residualfitZvsmeasX",20 ,-1000., 1000., 20 ,-10000., 10000.);
      _aidaHistoMap2D[bname + "residualfitZvsmeasY"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( bname + "residualfitZvsmeasY",20 ,-1000., 1000., 20 ,-10000., 10000.);
    
      //Angles
      _aidaHistoMap[bname + "dxdz"] = AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "dxdz", 10, -0.1, 0.1);
      _aidaHistoMap[bname + "dydz"] = AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "dydz", 10, -0.1, 0.1);
    }
}

void EUTelDafBase::bookDetailedHistos(){
  for( size_t ii = 0; ii < _system.planes.size() ; ii++) {
    daffitter::FitPlane<float>& plane = _system.planes.at(ii);
    char iden[4];
    sprintf(iden, "%d", plane.getSensorID());
    string bname = static_cast< string >("pl") + iden + "_";
    _aidaHistoMap[bname + "sigmaX"] =  AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "sigmaX", 10, 0.0f, 100);
    _aidaHistoMap[bname + "sigmaY"] =  AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "sigmaY", 10, 0.0f, 100);
    _aidaHistoMap[bname + "hitChi2"] =  AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "hitChi2", 10, 0, 100);
    _aidaHistoMap[bname + "pullX"] =  AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "pullX", 10, -2, 2);
    _aidaHistoMap[bname + "pullY"] =  AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "pullY", 10, -2, 2);
  }
}

void EUTelDafBase::end() {
  dafEnd();
  
  streamlog_out ( MESSAGE5 ) << endl;
  streamlog_out ( MESSAGE5 ) << "Number of found hit candidates: " << _nCandidates << endl;
  streamlog_out ( MESSAGE5 ) << "Tracks with ok ndof: " << n_passedNdof << endl;
  streamlog_out ( MESSAGE5 ) << "Tracks with BAD ndof: " << n_failedNdof << endl;
  streamlog_out ( MESSAGE5 ) << "Tracks with ok chi2/ndof: " << n_passedChi2OverNdof << endl;
  streamlog_out ( MESSAGE5 ) << "Tracks with BAD chi2/ndof: " << n_failedChi2OverNdof << endl;
  streamlog_out ( MESSAGE5 ) << "Tracks with no NaNs: " << n_passedIsnan<< endl;
  streamlog_out ( MESSAGE5 ) << "Tracks with NaNs: " << n_failedIsnan<< endl;
  streamlog_out ( MESSAGE5 ) << "Number of fitted tracks: " << _nTracks << endl;
  streamlog_out ( MESSAGE5 ) << "Successfully finished" << endl;
  for( size_t ii = 0; ii < _system.planes.size() ; ii++){
    daffitter::FitPlane<float>& plane = _system.planes.at(ii);
    char iden[4];
    sprintf(iden, "%d", plane.getSensorID());
    string bname = static_cast< string >("pl") + iden + "_";
    if( _aidaHistoMap[bname + "residualX"] != 0 && _aidaHistoMap[bname + "residualY"] != 0 )
      streamlog_out ( MESSAGE5 ) << "plane:" << ii <<
	"  x-stat :" <<  _aidaHistoMap[bname + "residualX"]->allEntries() <<
	"  x-mean:"  <<  _aidaHistoMap[bname + "residualX"]->mean() << 
	"  x-rms :"  <<  _aidaHistoMap[bname + "residualX"]->rms() << 
	"  y-stat :" <<  _aidaHistoMap[bname + "residualY"]->allEntries() <<
	"  y-mean:"  <<  _aidaHistoMap[bname + "residualY"]->mean() << 
	"  y-rms :"  <<  _aidaHistoMap[bname + "residualY"]->rms() << endl;
  }


}
#endif // USE_GEAR
