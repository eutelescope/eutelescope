// Version: $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
// built only if GEAR and MARLINUTIL are used 
#if defined(USE_GEAR)
// eutelescope includes ".h"
#include "EUTelAPIXKalman.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelSparseCluster2Impl.h"
#include "EUTelExceptions.h"
#include "EUTelPStream.h"
#include "EUTelAlignmentConstant.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"
#include "marlin/Exceptions.h"
#include "marlin/AIDAProcessor.h"

// marlin util includes
#include "mille/Mille.h"

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
using namespace APIXFitter;

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelAPIXKalman::_numberTracksLocalname   = "NumberTracks";
std::string EUTelAPIXKalman::_chi2Localname          = "Chi2";
std::string EUTelAPIXKalman::_ndofLocalname          = "ndof";
std::string EUTelAPIXKalman::_dxdzLocalname          = "dxdz";
std::string EUTelAPIXKalman::_dydzLocalname          = "dydz";
std::string EUTelAPIXKalman::_residualXLocalname     = "ResidualX";
std::string EUTelAPIXKalman::_residualYLocalname     = "ResidualY";

#endif

EUTelAPIXKalman::EUTelAPIXKalman () 
: Processor("EUTelAPIXKalman"),
  cutAngleX(0),
  cutAngleY(0),
  cutChi2(0),
  cutResids(0),
  cutInTime(0),
  _hitCollectionName(),
  _trackCollectionName(""),
  _fittrackvec(NULL),
  _excludePlanes(),
  _fixedPlanes(),
  _fixedTranslations(),
  _fixedZRotations(),
  _fixedScales(),
  _inTimeCheck(),
  _fixedX(false),
  _fixedY(false),
  _binaryFilename(""),
  _telescopeResolution(),
  _useResidualCuts(false),
  _residualsXMin(),
  _residualsYMin(),
  _residualsXMax(),
  _residualsYMax(),
  _shiftsX(),
  _shiftsY(),
  _scalesX(),
  _scalesY(),
  _generatePedeSteerfile(0),
  _pedeSteerfileName(""),
  _useHitResol(false),
  _runPede(false),
  _addToLCIO(false),
  _usePedeUserStartValues(false),
  _doScatter(false),
  _pedeUserStartValuesX(),
  _pedeUserStartValuesY(),
  _pedeUserStartValuesGamma(),
  _alignmentConstantLCIOFile(""),
  _alignmentConstantCollectionName(""),
  _normalizedResidualsMax(0.0),
  _eBeam(0.0),
  _nSkipMax(0),
  _minDxDz(0.0),
  _maxDxDz(0.0),
  _minDyDz(0.0),
  _maxDyDz(0.0),
  _maxChi2(0.0),
  _fitter(NULL),
  _zSort(),
  _planeHits(),
  _checkInTime(false),
  _nTracks(0),
  _expectedTracks(0),
  _iRun(0),
  _iEvt(0),
  _nMilleDataPoints(0),
  _nMilleTracks(0),
  _mille(NULL),
  _siPlanesParameters(NULL),
  _siPlanesLayerLayout(NULL),
  _aidaHistoMap(),
  _siPlaneZPosition(),
  _histogramSwitch(false)
 {
  _description = "This processor preforms track fitting. The tracks can be used for obtaining alignment using the Millepede program, or as final track fits.";

  // input collection
  std::vector<std::string > HitCollectionNameVecExample;
  HitCollectionNameVecExample.push_back("hit");
  registerInputCollections(LCIO::TRACKERHIT,"HitCollectionName", "Names of input hit collections", _hitCollectionName,HitCollectionNameVecExample);

  //Track finder options
  registerOptionalParameter("DistanceMax","Track finding: The maximum allowed normalized residual between predicted state and measurement", _normalizedResidualsMax, static_cast<float>(300.0));
  registerOptionalParameter("AllowNSkippedPlanes","How many planes is the fitter allowed to temporarily exclude? Should be 0 in case of alignment",_nSkipMax ,static_cast <int> (0));

  //Tracker system options
  registerOptionalParameter("ExcludePlanes","Exclude planes from fit/alignment, by sensor ids.",_excludePlanes ,std::vector<int>());
  registerOptionalParameter("IncludeScatter","Include scattering in the fit?",_doScatter ,static_cast <bool> (false));
  registerOptionalParameter("UseHitResolution","Use resolution found by hitmaker? If not use sigmas from TelescopeResolution",_useHitResol, static_cast<bool> (true));
  registerOptionalParameter("TelescopeResolution","Sigmas of the telescope resolution for fitter and Millepede. Ordered by z position of plane.",_telescopeResolution, std::vector<float> ());
  registerOptionalParameter("Ebeam", "Beam energy [GeV]", _eBeam,  static_cast < double > (120.0));

  // Alignment options
  registerOptionalParameter("FixedPlanes","Fix sensor planes in the alignment, by sensor ids.",_fixedPlanes ,std::vector<int>());
  registerOptionalParameter("FixedTranslations","Fix translations for sensor planes in the alignment, by sensor ids.",_fixedTranslations ,std::vector<int>());
  registerOptionalParameter("FixedZRotations","Fix z-rotation for sensor planes in the alignment, by sensor ids.",_fixedZRotations ,std::vector<int>());
  registerOptionalParameter("FixedScales","Fix scales for sensor planes in the alignment, by sensor ids.",_fixedScales ,std::vector<int>());

  registerOptionalParameter("FixedX","Fix scales and trans for all x.", _fixedX , false);
  registerOptionalParameter("FixedY","Fix scales and trans for all y.", _fixedY , false);

  //Millepede options
  registerOptionalParameter("PedeSteerfileName","Name of the steering file for the pede program.",_pedeSteerfileName, string("steer_mille.txt"));
  registerOptionalParameter("RunPede","Build steering file, binary input file, and execute the pede program.",_runPede, static_cast <bool> (false));
  registerOptionalParameter("BinaryFilename","Name of binary input file for Millepede.",_binaryFilename, string ("mille.bin"));
  registerOptionalParameter("AlignmentConstantLCIOFile","Name of LCIO db file where alignment constantds will be stored"
                            "constants (add .slcio)",_alignmentConstantLCIOFile, static_cast< string > ( "alignment.slcio" ) );
  registerOptionalParameter("AlignmentConstantCollectionName", "This is the name of the alignment collection to be saved into the slcio file",
                            _alignmentConstantCollectionName, static_cast< string > ( "alignment" ));

  //Track quality parameters
  registerOptionalParameter("UseResidualCuts","Use cuts on the residuals to reduce the combinatorial background?", _useResidualCuts, static_cast <bool> (false));
  registerOptionalParameter("ResidualsXMin","Minimal values of the hit residuals in the X direction for a track. Ordered by z position of plane.",_residualsXMin, std::vector<float>());
  registerOptionalParameter("ResidualsXMax","Maximal values of the hit residuals in the X direction for a track. Ordered by z position of plane.",_residualsXMax, std::vector<float>());
  registerOptionalParameter("ResidualsYMin","Minimal values of the hit residuals in the Y direction for a track. Ordered by z position of plane.",_residualsYMin, std::vector<float>());
  registerOptionalParameter("ResidualsYMax","Maximal values of the hit residuals in the Y direction for a track. Ordered by z position of plane.",_residualsYMax, std::vector<float>());

  registerOptionalParameter("MinDxDz", "Minimum allowed dx/dz", _minDxDz, static_cast<double> (-10.0));
  registerOptionalParameter("MaxDxDz", "Maximum allowed dx/dz", _maxDxDz, static_cast<double> ( 10.0));
  registerOptionalParameter("MinDyDz", "Minimum allowed dy/dz", _minDyDz, static_cast<double> (-10.0));
  registerOptionalParameter("MaxDyDz", "Maximum allowed dy/dz", _maxDyDz, static_cast<double> ( 10.0));
  registerOptionalParameter("MaxChi2", "Maximum allowed global chi2", _maxChi2, static_cast<double> ( 9999.0));
  registerOptionalParameter("InTimeCheck","DUT sensor id's to be checked with residual cuts after the final track fit.",_inTimeCheck,std::vector<int>());

  //Track fitter options
  registerOptionalParameter("AddToLCIO","Add track collection to LCIO event.",_addToLCIO, static_cast <bool> (false));
  registerOutputCollection(LCIO::TRACK,"TrackCollectionName", "Collection name for fitted tracks", _trackCollectionName, string ("fittracks"));

}

void EUTelAPIXKalman::init() {
  //Move Fixed plane vector to vector for spesific parameters 
  for(uint ii = 0; ii < _fixedPlanes.size() ;ii++){
    int sensid = _fixedPlanes.at(ii);
    _fixedScales.push_back(sensid);
    _fixedZRotations.push_back(sensid);
    _fixedTranslations.push_back(sensid);
  }
  printParameters ();
  _iRun = 0; _iEvt = 0;

  if ( Global::GEAR == NULL ) {
    streamlog_out ( ERROR2) << "The GearMgr is not available." << endl;
    exit(-1);
  }

  _siPlanesParameters  = const_cast<gear::SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

  //Use map to sort planes by z
  _zSort.clear();
  for(int plane = 0; plane < _siPlanesLayerLayout->getNLayers(); plane++){
    _zSort[ _siPlanesLayerLayout->getLayerPositionZ(plane) ] = plane;
  }
  //Making space for all hits
  _planeHits.resize( _zSort.size() );
  //If we do not want to use supplied resolution, make sure vector is empty
  if(_useHitResol){ _telescopeResolution.clear();}
  //The fitter is implemented in the APIXFitter namespace
  _fitter = new APIXKalman();
  //Adding planes to fitter with increasing z
  map<double, int>::iterator zit = _zSort.begin();
  int index(0), nActive(0);
  for(; zit != _zSort.end(); ++zit, index++){
    int sensorID = _siPlanesLayerLayout->getID( (*zit).second );
    double zPos  = (*zit).first; 
    if(index >= static_cast< int >(_telescopeResolution.size())){
      _telescopeResolution.push_back( 1000.0 * _siPlanesLayerLayout->getSensitiveResolution( (*zit).second ) );
    }
    double errX = _telescopeResolution.at(index);
    double errY = _telescopeResolution.at(index);

    bool excluded = (find(_excludePlanes.begin(), _excludePlanes.end(), sensorID)
		     != _excludePlanes.end());
    
    if(not excluded){nActive++;}

    double scatter = _doScatter ? getScatterCov( (*zit).second ) : 0.0;
    //The fitter holds the information about the planes in a FitPlane instance. To account
    //for scatter, excluded planes should be added as well
    FitPlane* pl = new FitPlane(index, sensorID, zPos * 1000.0, errX, errY, scatter, excluded);
    pl->print();
    _fitter->addPlane(index, pl);
  }
  if(nActive - _nSkipMax  < 2) {
    streamlog_out ( ERROR5 ) << "Too few active planes(" << nActive << ") when " << _nSkipMax << " planes can be skipped." 
			    << "Please check your configuration." << endl;
    exit(1);
  }
  _checkInTime = (!_inTimeCheck.empty());

  _nMilleDataPoints = 0;
  _nMilleTracks = 0;
  if(_runPede){
    _mille = new Mille(_binaryFilename.c_str());
    streamlog_out ( MESSAGE5 ) << "The filename for the mille binary file is: " << _binaryFilename.c_str() << endl;
    if(_nSkipMax > 0){
      streamlog_out ( ERROR5 ) << "I'm planning on doing alignment, in this case skipping planes is a very bad idea. Setting number of allowed skipped planes to 0." << endl;
      _nSkipMax = 0;
    }
  }
  _histogramSwitch = true;
  bookHistos();
}

double EUTelAPIXKalman::getScatterCov(int index){
  // x / x0
  double radLength = _siPlanesLayerLayout->getLayerThickness( index ) /  _siPlanesLayerLayout->getLayerRadLength( index );
  //From pdg live
  double scatterTheta = 0.0136/_eBeam * sqrt( radLength ) *  (1 + 0.038 * std::log(radLength) );
  //We really want tangent of angle, propagate error.
  float scatterTan = (1 + tan(scatterTheta) * (1 + scatterTheta)) * scatterTheta;
  return(scatterTan * scatterTan);
}

void EUTelAPIXKalman::processRunHeader (LCRunHeader * rdr) {
  auto_ptr<EUTelRunHeaderImpl> header ( new EUTelRunHeaderImpl (rdr) );
  header->addProcessor( type() ) ;
  ++_iRun;
}

int EUTelAPIXKalman::getPlaneIndex(double zPos){
  map<double,int>::iterator it = _zSort.begin();
  int index(0);
  bool foundIt(false);
  for(;it != _zSort.end(); index++, it++){
    if( fabs((*it).first - zPos) < 2.5){ 
      foundIt = true; break;
    }
  }
  if(not foundIt){ 
    streamlog_out ( ERROR5 ) << "Found hit at z=" << zPos << " , not able to assign to any plane!" << endl; 
    return(-1);
  }
  return(index);
}

void EUTelAPIXKalman::readHitCollection(LCEvent* event){
  //Clear _planeHits
  for(int ii = 0; ii < static_cast< int >(_planeHits.size()); ii++){ _planeHits.at(ii).clear();}
  //Extract hits from collection, add to 
  LCCollection* collection;
  for(size_t i =0;i < _hitCollectionName.size();i++){
    try {
      collection = event->getCollection(_hitCollectionName[i]);
    } catch (DataNotAvailableException& e) {
      streamlog_out ( WARNING2 ) << "No input collection " << _hitCollectionName[i] << " found for event " << event->getEventNumber()
				 << " in run " << event->getRunNumber() << endl;
      throw SkipEventException(this);
    }
    
    //Add all hits in collection to corresponding plane

    for ( int iHit = 0; iHit < collection->getNumberOfElements(); iHit++ ) {
      TrackerHitImpl* hit = static_cast<TrackerHitImpl*> ( collection->getElementAt(iHit) );
      int planeIndex = getPlaneIndex( hit->getPosition()[2] );
      //Check charge of telescope clusters, assume single hits to be noise

      if( _runPede and (hit->getType() == kEUTelAPIXClusterImpl) ){
	auto_ptr<EUTelVirtualCluster> cluster( new EUTelSparseClusterImpl< EUTelAPIXSparsePixel >
					       ( static_cast<TrackerDataImpl *> ( hit->getRawHits()[0] )));
	float xPos(0), yPos(0);
	cluster->getCenterOfGravity(xPos, yPos);
	//Discard all hits close to edge, and in the region of ganged pixels
	//Unless my logic is flawed, this discards single column clusters with all hits in col 0.
	//In cases of charge sharing, edge should behave as normal
	if( (xPos < 0.1) ) { continue; }
	if( (xPos > 16.6) ) { continue; }
	//if( (yPos < 15) or (yPos > 144) ) { continue; }
	//Drop huge clusters
	//int xSize(0), ySize(0);
	//cluster->getClusterSize(xSize, ySize);
	//if( xSize != 2 ) { continue; }
      }
      if(planeIndex >=0 ) { _planeHits.at(planeIndex).push_back( hit );}
    }
  }
}

void EUTelAPIXKalman::processEvent (LCEvent * event) {
  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);
  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG2 ) << "EORE found: nothing else to do." << endl;
    return;
  }
  _nTracks = 0;
  //Dump hit collection to collection sorted by plane
  readHitCollection(event);
  if(_addToLCIO){
    // Define output track and hit collections
    _fittrackvec = new LCCollectionVec(LCIO::TRACK);
    // Set flag for storing track hits in track collection
    LCFlagImpl flag(_fittrackvec->getFlag());
    flag.setBit( LCIO::TRBIT_HITS );
    _fittrackvec->setFlag(flag.getFlag());
  }
  // Run the fit
  TrackEstimate* estim = new TrackEstimate();
  estim->makeSeed();
  fitPermutations(0, NULL, estim, 0);
  delete estim;
  //Plot number of found tracks
  tryFill( _numberTracksLocalname, _nTracks);
  if(_addToLCIO){ event->addCollection(_fittrackvec,_trackCollectionName); }
  if(event->getEventNumber() % 100 == 0){
    streamlog_out ( MESSAGE5 ) << "Accepted " << _nMilleTracks <<" tracks at event " << event->getEventNumber() << endl;
  }
}
void EUTelAPIXKalman::fitPermutations(int plane, FitPlane* prev, TrackEstimate* est, int nSkipped){
  //Got track?
  if(_nTracks > 15) { return; }
  if(plane == static_cast< int >(_planeHits.size())){
    finalizeTrack();
    return;
  }
  FitPlane* cur = _fitter->getPlane(plane);
  //Extrapolate track to current plane
  if( plane > 0){
    _fitter->predict(prev, cur, est);
  }

  // Is plane excluded?
  if(cur->excluded){
    _fitter->update(cur, est);
    fitPermutations(plane + 1, cur, est, nSkipped);
    return;
  }

  //Does the estimate est appear to be a seed?
  bool isSeed = (gsl_vector_get(est->param, 0) == 0.0) and (gsl_vector_get(est->param, 1) == 0.0);
  //How many tracks are found before this track candidate has been fully checked?
  int tmpNtracks = _nTracks;
  for(uint hit = 0; hit < _planeHits.at(plane).size(); hit++){
    addPlaneHit(cur, _planeHits.at(plane).at(hit));
    //Check normalized residuals if estimate is a real prediction
    if( (not isSeed) ){ 
      if( fabs( gsl_vector_get(est->param, 0) - cur->hitPosX)/cur->errX > _normalizedResidualsMax ) continue;
      if( fabs( gsl_vector_get(est->param, 1) - cur->hitPosY)/cur->errY > _normalizedResidualsMax ) continue;
    }
    //Since a track candidate can branch out, we need to clone the estimate in order to not
    //contaminate the other track candidates
    TrackEstimate* clone = new TrackEstimate();
    if(plane == 0){ clone->makeSeed(); }
    else { clone->copy(est);}
    _fitter->update(cur, clone);
    fitPermutations(plane + 1, cur, clone, nSkipped);
    delete clone;
  }
  //If no measurement in this plane lead to an accepted track, try to skip this plane by
  //temporarily excluding this plane. Note, this will only work right when we expect at most
  //one track per trigger. For APIX data this is the case.
  if( (_nTracks == tmpNtracks)  and nSkipped < _nSkipMax){
    cur->excluded = true;
    _fitter->update(cur, est);
    fitPermutations(plane + 1, cur, est, nSkipped + 1);
    cur->excluded = false;
  }
}
void EUTelAPIXKalman::addPlaneHit(FitPlane* pl, TrackerHitImpl* hit){
  const double * pos = hit->getPosition();
  //Fitter uses microns, framework uses mm
  pl->hitPosX = pos[0] * 1000.0;
  pl->hitPosY = pos[1] * 1000.0;

  pl->errX = _telescopeResolution.at( pl->index);
  pl->errY = _telescopeResolution.at( pl->index);
  if(pl->sensorID > 7){
    pl->errY = _telescopeResolution.at( pl->index) * 8.0;
  }
  
  // if(not _useHitResol) { return; }
  // //If _useHitResol and non singular hit error matrix, extract errors from hit.
  // const EVENT::FloatVec cov = hit->getCovMatrix();
  // if(cov.at(0)>0.) { pl->errX = 1000.0 * sqrt(cov.at(0)); }
  // if(cov.at(2)>0.) { pl->errX = 1000.0 * sqrt(cov.at(2)); }
}
void EUTelAPIXKalman::addToMille(){
  const int nLC = 4; //number of local parameters
  const int nGL = _planeHits.size() * 5; // number of global parameters

  float *derLC = new float[nLC]; // array of derivatives for local parameters
  float *derGL = new float[nGL]; // array of derivatives for global parameters
  int *label = new int[nGL]; // array of parameter labels
  
  for(int ii = 0; ii < nGL; ii++){ 
    label[ii] = ii + 1;
    derGL[ii] = 0;
  }
  for(int ii = 0; ii < nLC; ii++){ derLC[ii] = 0; }

  int labelMe = 0;
  map<int, FitPlane*>::iterator it = _fitter->indexToPlane.begin();
  for(;it != _fitter->indexToPlane.end(); ++it){
    if((*it).second->excluded) { continue; }
    derGL[(labelMe * 5)    ] = -1; //Derivatives of residuals w.r.t. shift in x
    derGL[(labelMe * 5) + 2] = (*it).second->hitPosY; //Derivatives of residuals w.r.t. z rotations
    derGL[(labelMe * 5) + 3] = (*it).second->hitPosX; //Derivatives of residuals w.r.t. scale
						      //of x axis
    derLC[0] = 1; //Derivatives of fit pos w.r.t. x
    derLC[2] = (*it).second->posZ; //Derivatives of fit pos w.r.t. dx/dz
    _mille->mille(nLC, derLC, nGL, derGL, label, (*it).second->resX, (*it).second->errX);

    derGL[(labelMe * 5)]     = 0;
    derGL[(labelMe * 5) + 2] = 0;
    derGL[(labelMe * 5) + 3] = 0;
    derLC[0] = 0;
    derLC[2] = 0;

    derGL[(labelMe * 5) + 1] = -1; //Derivatives of residuals w.r.t. shift in y
    derGL[(labelMe * 5) + 2] = -1 * (*it).second->hitPosX; //Derivatives of residuals
							   //w.r.t. z rotations
    derGL[(labelMe * 5) + 4] = (*it).second->hitPosY;//Derivatives of residuals w.r.t. scales
						     //of y axis

    derLC[1] = 1; //Derivatives of fit pos w.r.t. y
    derLC[3] = (*it).second->posZ; //Derivatives of fit pos w.r.t. dy/dz
    
    _mille->mille(nLC, derLC, nGL, derGL, label, (*it).second->resY, (*it).second->errY);
    
    derGL[(labelMe * 5) + 1] = 0;
    derGL[(labelMe * 5) + 2] = 0;
    derGL[(labelMe * 5) + 4] = 0;
    derLC[1] = 0;
    derLC[3] = 0;

    labelMe++;
    _nMilleDataPoints++;
  }

  delete [] derLC;
  delete [] derGL;
  delete [] label;

  _mille->end();
}

bool EUTelAPIXKalman::goodResiduals(FitPlane* pl){
  //Check if residuals in plane passes cuts
  if( pl->index < static_cast< int >(_residualsXMax.size()) and pl->index < static_cast< int >(_residualsXMin.size()) ){
    if(pl->resX > _residualsXMax.at( pl->index )) return(false); 
    if(pl->resX < _residualsXMin.at( pl->index )) return(false); 
  }
  if( pl->index < static_cast< int >(_residualsYMax.size()) and pl->index < static_cast< int >(_residualsYMin.size()) ){
    if(pl->resY > _residualsYMax.at( pl->index )) return(false); 
    if(pl->resY < _residualsYMin.at( pl->index )) return(false); 
  }
  return(true);
}

bool EUTelAPIXKalman::inTimeGood(FitPlane* pl){
  //Check if track is intime with plane using residuals
  if(find (_inTimeCheck.begin(), _inTimeCheck.end(), pl->sensorID) == _inTimeCheck.end()){ return(false);}
  for(int hit = 0; hit < static_cast< int >(_planeHits.at(pl->index).size()); hit++ ){
    addPlaneHit(pl, _planeHits.at(pl->index).at(hit));
    if( goodResiduals(pl) ){ return(true); }
  }
  return(false);
}
void EUTelAPIXKalman::finalizeTrack(){
  // Run Kalman smoother to obtain optimal estimates at each plane.
  _fitter->smoothPlanes();
  double chi2(0.0);
  int ndof(-4);
  int intimeplanes(0);
  FitPlane* ppl = _fitter->getPlane(0);
  //Check angles
  if( (ppl->fitdxdz > _maxDxDz) or (ppl->fitdxdz < _minDxDz) ) { return; }
  if( (ppl->fitdydz > _maxDyDz) or (ppl->fitdydz < _minDyDz) ) { return; }

  map<int, FitPlane*>::reverse_iterator rit = _fitter->indexToPlane.rbegin();
  for(; rit != _fitter->indexToPlane.rend(); ++rit){
    FitPlane* pl = (*rit).second;
    if( pl->excluded) { 
      if( _checkInTime and inTimeGood(pl)) { intimeplanes++; }
      continue; 
    }
    // Calculate track chi2 increment
    chi2 += (pl->resX * pl->resX) / (pl->errX * pl->errX);
    chi2 += (pl->resY * pl->resY) / (pl->errY * pl->errY);
    //Kill track if chi2 is to large
    if(chi2 > _maxChi2) { return; }
    //Add measurements to ndof
    ndof += 2;
    // Check residuals
    if( _useResidualCuts and (not goodResiduals(pl))){ return; }
  }
  //Kill track if we want an intime check and no matching DUT hits are found
  if( (!_inTimeCheck.empty()) and intimeplanes == 0) { return; }
  _nTracks++;
  _nMilleTracks++;
  if( _runPede ) { addToMille(); }
  if( _histogramSwitch ){ plotResiduals(chi2, ndof); }
  if( _addToLCIO) { addToLCIO(chi2, ndof); }
}

void EUTelAPIXKalman::addToLCIO(double chi2, int ndof){

  TrackImpl * fittrack = new TrackImpl();
  // Impact parameters are useless and set to 0
  fittrack->setD0(0.);        // impact paramter of the track in (r-phi)
  fittrack->setZ0(0.);        // impact paramter of the track in (r-z)
  fittrack->setTanLambda(0.); // dip angle of the track at reference point

  TrackEstimate* est = _fitter->estimates.at(0);

  //No good way of storing the track angles, so
  fittrack->setOmega( gsl_vector_get(est->param, 2)); // Storing dxdz as Omega
  fittrack->setPhi( gsl_vector_get(est->param, 3));   // Storing dx/dy as phi

  fittrack->setChi2(chi2);
  fittrack->setNdf(ndof);
  
//  fittrack->setIsReferencePointPCA(false);
  float refpoint[3];
  
  map<int, FitPlane*>::iterator it = _fitter->indexToPlane.begin();
  for(;it != _fitter->indexToPlane.end(); ++it){
    FitPlane* pl = (*it).second;
    TrackEstimate* estim = _fitter->estimates.at( pl->index);
    TrackerHitImpl * fitpoint = new TrackerHitImpl();
    // Hit type 32 means a fitted hit
    fitpoint->setType(32);

    double pos[3];
    pos[0]= pl->fitX / 1000.0;
    pos[1]= pl->fitY / 1000.0;
    pos[2]= pl->posZ / 1000.0;
    
    fitpoint->setPosition(pos);
    // Covariance matrix of the fitted position
    // (stored as lower triangle matrix, i.e.  cov(xx),cov(y,x),cov(y,y) ).
    float cov[TRKHITNCOVMATRIX];
    cov[0]= gsl_matrix_get(estim->Cov, 0, 0);
    cov[1]= gsl_matrix_get(estim->Cov, 0, 1);
    cov[2]= gsl_matrix_get(estim->Cov, 1, 1);
    //Error of z position of fit is unclear to me, this would be a systematic alignment
    //error. Set to 0 along with all covariances.
    cov[3]=cov[4]=cov[5]=0.;
    fitpoint->setCovMatrix(cov);

    fittrack->addHit(fitpoint);
    if(pl->index == 0){
      refpoint[0] = pos[0];
      refpoint[1] = pos[1];
      refpoint[2] = pos[2];
    }
  }
  fittrack->setReferencePoint(refpoint);
  _fittrackvec->addElement(fittrack);
}
void EUTelAPIXKalman::tryFill(string name, double value){
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  AIDA::IHistogram1D* tmp = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[name]);
  if( tmp == NULL ){
    streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << name << endl;
    streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
    _histogramSwitch = false;
  }
  tmp->fill(value);
#endif
}

void EUTelAPIXKalman::plotResiduals(double chi2, int ndof){
  if(not _histogramSwitch) { return;}
  string tempHistoName;
  
  map<int, FitPlane*>::iterator it = _fitter->indexToPlane.begin();
  tryFill( _chi2Localname, chi2);
  tryFill( _ndofLocalname, ndof);
  tryFill( _dxdzLocalname, (*it).second->fitdxdz);
  tryFill( _dydzLocalname, (*it).second->fitdydz);

  for(; it != _fitter->indexToPlane.end(); ++it ){
    //Plot in-time planes and non excluded planes
    if((*it).second->excluded and
       not (find (_inTimeCheck.begin(), _inTimeCheck.end(), (*it).second->sensorID) != _inTimeCheck.end())){
      continue; 
    }
    int sensorID = (*it).second->sensorID;
    tryFill( _residualXLocalname + "_d" + to_string( sensorID), (*it).second->resX );
    tryFill( _residualYLocalname + "_d" + to_string( sensorID), (*it).second->resY );
  }
}

void EUTelAPIXKalman::tryBook(string name, string title,  int bins, double min, double max){
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  AIDA::IHistogram1D * tmpHisto =
    AIDAProcessor::histogramFactory(this)->createHistogram1D(name , bins, min, max);
  if(tmpHisto == NULL){
    streamlog_out ( ERROR2 ) << "Problem booking the " << (name) << endl;
    streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
    _histogramSwitch = false;
    return;
  }
  tmpHisto->setTitle(title);
  _aidaHistoMap.insert( make_pair( name, tmpHisto ) );
#endif
}

void EUTelAPIXKalman::bookHistos() {
  const int    tracksNBin = 20;
  const double tracksMin = -0.5;
  const double tracksMax = 19.5;
  const int    Chi2NBin =1000; 
  const double Chi2Min  = 0.;
  const double Chi2Max  =static_cast< int >(_maxChi2 + 1);
  const int    NBin = 10000;
  const double Min  = -20000.;
  const double Max  = 20000.;
  const int ddzBins = 100000;
  const double ddzMin = -1;
  const double ddzMax =  1;
    
  tryBook(_numberTracksLocalname,"number of tracks after cuts", tracksNBin, tracksMin, tracksMax );
  tryBook(_chi2Localname, "Chi2", Chi2NBin, Chi2Min, Chi2Max);
  tryBook(_ndofLocalname, "ndof", 15, 0, 15);
  tryBook(_dxdzLocalname, "dxdz", ddzBins, ddzMin, ddzMax);
  tryBook(_dydzLocalname, "dydz", ddzBins, ddzMin, ddzMax);

  map<int, FitPlane*>::iterator it = _fitter->indexToPlane.begin();
  for(; it != _fitter->indexToPlane.end(); ++it ){
    string sensorID = to_string((*it).second->sensorID);
    string resXtitle = "Residual X for " + sensorID;
    tryBook(_residualXLocalname + "_d" + sensorID, resXtitle, NBin, Min, Max);

    string resYtitle = "Residual Y for " + sensorID;
    tryBook(_residualYLocalname + "_d" + sensorID, resYtitle, NBin, Min, Max);
  }
}

vector<double> EUTelAPIXKalman::initAlignParams(){
  //Get straight line fit to fixed planes in the setup
  map<int, FitPlane*>::iterator it = _fitter->indexToPlane.begin();
  string tempHistoName;
  for(;it!= _fitter->indexToPlane.end(); ++it){
    FitPlane* pl = (*it).second;
    bool fixedTrans = (find (_fixedTranslations.begin(), _fixedTranslations.end(), 
			     pl->sensorID) != _fixedTranslations.end());
    if( fixedTrans and (not (*it).second->excluded)){
      (*it).second->excluded = false;
    } else { (*it).second->excluded = true; }    
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    tempHistoName = _residualXLocalname + "_d" + to_string( (*it).second->sensorID );
    AIDA::IHistogram1D* residx_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName.c_str()]);
    tempHistoName = _residualYLocalname + "_d" + to_string( (*it).second->sensorID );
    AIDA::IHistogram1D* residy_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName.c_str()]);
#endif
    if( not (residx_histo == NULL) ){
      pl->hitPosX = residx_histo->mean();
      pl->errX = 10.0;
    } else { pl->excluded = true; }
    if( not (residy_histo == NULL) ){
      pl->hitPosY = residy_histo->mean();
      pl->errY = 10.0;
    } else { pl->excluded = true; }
  }
  _fitter->fitPlanes();
  //Use straight line fit to initialize shifts. Rotations and scales are set to 0.
  vector<double> alignParams;
  for(it = _fitter->indexToPlane.begin(); it != _fitter->indexToPlane.end(); ++it){
    FitPlane* pl = (*it).second;
    if(find (_excludePlanes.begin(), _excludePlanes.end(), 
	     pl->sensorID) != _excludePlanes.end()) {
      pl->excluded = true;
      continue;
    }
    pl->excluded = false;
    bool fixedTrans = (find (_fixedTranslations.begin(), _fixedTranslations.end(), 
			     pl->sensorID) != _fixedTranslations.end());
    if( fixedTrans or _fixedX ){
      alignParams.push_back(0.0);
    } else {
      alignParams.push_back( pl->fitX - pl->hitPosX );
    }
    if( fixedTrans or _fixedY){
      alignParams.push_back(0.0);
    } else {
      alignParams.push_back( pl->fitY - pl->hitPosY );
    }
    alignParams.push_back( 0.0);
    alignParams.push_back( 0.0);
    alignParams.push_back( 0.0);
  }
  return(alignParams);
}

bool EUTelAPIXKalman::generatePedeSteeringFile( vector<double> &startVals , bool fixShift, bool fixRotate, bool fixScale){
  ofstream steerFile;
  steerFile.open(_pedeSteerfileName.c_str());
  if (not steerFile.is_open()) {
    streamlog_out ( ERROR2 ) << "Unable to open pede steering file." << endl;
    return(false);
  }
  steerFile << "Cfiles" << endl;
  steerFile << _binaryFilename << endl;
  steerFile << endl;
  steerFile << "Parameter" << endl;
  int counter(0);
  map<int, FitPlane*>::iterator it = _fitter->indexToPlane.begin();
  for(; it != _fitter->indexToPlane.end(); ++it){
    FitPlane* pl = (*it).second;
    if(find (_excludePlanes.begin(), _excludePlanes.end(), 
	     pl->sensorID) != _excludePlanes.end()) { continue; }
    bool fixedTrans = (find (_fixedTranslations.begin(), _fixedTranslations.end(), 
			     pl->sensorID) != _fixedTranslations.end());
    if(fixedTrans or fixShift or _fixedX){
      steerFile << (counter + 1) << " " << startVals.at(counter) << " -1.0" << endl;
      counter++;
    } else {
      steerFile << (counter + 1) << " " << startVals.at(counter) << " 0.0" << endl;
      counter++;
    }
    if(fixedTrans or fixShift or _fixedY){
      steerFile << (counter + 1) << " " << startVals.at(counter) << " -1.0" << endl;
      counter++;
    } else {
      steerFile << (counter + 1) << " " << startVals.at(counter) << " 0.0" << endl;
      counter++;
    }

    bool fixedZRot = (find (_fixedZRotations.begin(), _fixedZRotations.end(), 
			    pl->sensorID) != _fixedZRotations.end());
    if(fixedZRot or fixRotate){ 
      steerFile << (counter + 1) << " " << startVals.at(counter) << " -1.0" << endl;
      counter++;
    } else {
      steerFile << (counter + 1) << " " << startVals.at(counter) << " 0.0" << endl;
      counter++;
    }

    bool fixedScale = (find (_fixedScales.begin(), _fixedScales.end(), 
			     pl->sensorID) != _fixedScales.end());
    if(fixedScale or fixScale or _fixedX){ 
      steerFile << (counter + 1) << " " << startVals.at(counter) << " -1.0" << endl;
      counter++;
    } else {
      steerFile << (counter + 1) << " " << startVals.at(counter) << " 0.0" << endl;
      counter++;
    }
    if(fixedScale or fixScale or _fixedY){ 
      steerFile << (counter + 1) << " " << startVals.at(counter) << " -1.0" << endl;
      counter++;
    } else {
      steerFile << (counter + 1) << " " << startVals.at(counter) << " 0.0" << endl;
      counter++;
    }
  }
  steerFile << endl;
  steerFile << "! chiscut 5.0 2.5" << endl;
  steerFile << "! outlierdownweighting 4" << endl;
  steerFile << endl;
  steerFile << "method inversion 10 0.001" << endl;
  steerFile << endl;
  steerFile << "histprint" << endl;
  steerFile << endl;
  steerFile << "end" << endl;
  steerFile.close();
  streamlog_out ( MESSAGE5 ) << "File " << _pedeSteerfileName << " written." << endl;
  return(true);
}

void EUTelAPIXKalman::runPede(vector<double> &alignmentParams ){
  std::string command = "pede " + _pedeSteerfileName;
  // before starting pede, let's check if it is in the path
  bool isPedeInPath = true; 
  // create a new process
  redi::ipstream which("which pede");
  // wait for the process to finish
  which.close();
  
  // if status = 255 then the program wasn't found in the path
  isPedeInPath = !( which.rdbuf()->status() == 255 );
  
  if ( !isPedeInPath ) {
    streamlog_out( ERROR5 ) << "Cannot find pede program in path. Nothing to do." << endl;
    return;
  }
  streamlog_out ( MESSAGE5 ) << "Starting pede..." << endl;

  redi::ipstream pede( command.c_str() );
  string output;
  while ( getline( pede, output ) ) { streamlog_out( MESSAGE5 ) << output << endl; }
  // wait for the pede execution to finish
  pede.close();
  // check the exit value of pede
  if ( pede.rdbuf()->status() == 0 ) {
    streamlog_out ( MESSAGE5 ) << "Pede successfully finished" << endl;
  } else {
    streamlog_out ( ERROR2 ) << "Pede exited abnormally with exit code " << pede.rdbuf()->status() << endl;
  }
  // reading back the millepede.res file and getting the results. 
  string millepedeResFileName = "millepede.res";

  streamlog_out ( MESSAGE5 ) << "Reading back " << millepedeResFileName << endl
			     << "Saving the alignment constant into " << _alignmentConstantLCIOFile << endl;


  // open the millepede ASCII output file
  ifstream millepede( millepedeResFileName.c_str() );

  if(not millepede.is_open() ){
    streamlog_out ( MESSAGE5 ) << "Unable to open " << millepedeResFileName << " for reading" << endl; 
    return;
  }

  // reopen the LCIO file this time in append mode
  LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();

  try {
    lcWriter->open( _alignmentConstantLCIOFile, LCIO::WRITE_NEW );
  } catch ( IOException& e ) {
    streamlog_out ( ERROR4 ) << e.what() << endl;
    exit(-1);
  }

  // write an almost empty run header
  LCRunHeaderImpl * lcHeader  = new LCRunHeaderImpl;
  lcHeader->setRunNumber( 0 );
  lcWriter->writeRunHeader(lcHeader);
  delete lcHeader;

  LCEventImpl * event = new LCEventImpl;
  event->setRunNumber( 0 );
  event->setEventNumber( 0 );

  LCTime * now = new LCTime;
  event->setTimeStamp( now->timeStamp() );
  delete now;

  LCCollectionVec * constantsCollection = new LCCollectionVec( LCIO::LCGENERICOBJECT );

  if ( millepede.bad() ) {
    streamlog_out ( ERROR4 ) << "Error opening the " << millepedeResFileName << endl
			     << "The alignment slcio file cannot be saved" << endl;
    return;
  }

  vector<double > tokens;
  stringstream tokenizer;
  string line;
  double buffer;

  // get the first line and throw it away since it is a comment!
  getline( millepede, line );

  int counter = 0;
  map<int, FitPlane*>::iterator it = _fitter->indexToPlane.begin();
  while ( not millepede.eof() ) {
    EUTelAlignmentConstant * constant = new EUTelAlignmentConstant;
    bool goodLine = true;
    for ( unsigned int iParam = 0 ; iParam < 5 ; ++iParam ) {
      getline( millepede, line );
      goodLine = not line.empty();

      tokens.clear();
      tokenizer.clear();
      tokenizer.str( line );

      while ( tokenizer >> buffer ) { tokens.push_back( buffer ); }

      goodLine = ( tokens.size() == 3 ) || ( tokens.size() == 5) || (tokens.size() == 6);
      if(not goodLine) { continue; }
      bool isFixed = ( tokens.size() == 3 );
      if ( isFixed ) {
	streamlog_out ( DEBUG0 ) << "Parameter " << tokens[0] << " is at " << ( tokens[1] / 1000 )
				 << " (fixed)"  << endl;
      } else {
	streamlog_out ( DEBUG0 ) << "Parameter " << tokens[0] << " is at " << (tokens[1] / 1000 )
				 << " +/- " << ( tokens[4] / 1000 )  << endl;
      }
      alignmentParams.at(counter++) = tokens[1];
      
      if ( iParam == 0 ) {
	constant->setXOffset( tokens[1] / 1000 );
	if ( ! isFixed ) {
	  double err  = tokens[4] / 1000;
	  constant->setXOffsetError( err ) ;
	}
      }
      if ( iParam == 1 ) {
	constant->setYOffset( tokens[1] / 1000 ) ;
	if ( ! isFixed ) constant->setYOffsetError( tokens[4] / 1000 ) ;
      }
      if ( iParam == 2 ) {
	constant->setGamma( tokens[1]  ) ;
	if ( ! isFixed ) constant->setGammaError( tokens[4] ) ;
      }
      if ( iParam == 3 ) {
	constant->setAlpha( tokens[1]  ) ;
	if ( ! isFixed ) constant->setAlphaError( tokens[4] ) ;
      }
      if ( iParam == 4 ) {
	constant->setBeta( tokens[1]  ) ;
	if ( ! isFixed ) constant->setBetaError( tokens[4] ) ;
      }
    }
    // add the constant to the collection
    if ( goodLine ) {
      while ( find (_excludePlanes.begin(), _excludePlanes.end(), 
		    (*it).second->sensorID) != _excludePlanes.end()) {++it;}
      constant->setSensorID( (*it).second->sensorID );
      ++it;
      constantsCollection->push_back( constant );
      streamlog_out ( MESSAGE5 ) << (*constant) << endl;
    } else{ delete constant; }
  }
  event->addCollection( constantsCollection, _alignmentConstantCollectionName );
  lcWriter->writeEvent( event );
  delete event;
  lcWriter->close();
  millepede.close();
}

void EUTelAPIXKalman::end() {
  if(_runPede){
    delete _mille;
    vector<double> alignParams = initAlignParams();
    generatePedeSteeringFile(alignParams, false, false, false);
    runPede(alignParams);
  }
  streamlog_out ( MESSAGE5 ) << endl;
  streamlog_out ( MESSAGE5 ) << "Number of data points used: " << _nMilleDataPoints << endl;
  streamlog_out ( MESSAGE5 ) << "Number of tracks used: " << _nMilleTracks << endl;
  streamlog_out ( MESSAGE5 ) << endl;
  streamlog_out ( MESSAGE5 ) << "Successfully finished" << endl << flush << flush;
}
#endif // USE_GEAR
