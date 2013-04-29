// Version: $Id$
#include "EUTelDafTrackerSystem.h"

#include <iostream>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Array>

#include <TVector3.h>


using namespace std;
using namespace daffitter;

FitPlane::FitPlane(int sensorID, float zPos, float sigmaX, float sigmaY, float scatterThetaSqr, bool excluded){
  this->sensorID = sensorID;
  this->zPosition = zPos;
  this->sigmas(0) = sigmaX; sigmas(1) = sigmaY;
  this->variances = sigmas.cwise().square();
  this->invMeasVar(0) = 1/(sigmaX * sigmaX);
  this->invMeasVar(1) = 1/(sigmaY * sigmaY);
  this->scatterThetaSqr = scatterThetaSqr;
  this->excluded = excluded;
}

void FitPlane::print(){
  //Pretty print initialized plane info
  streamlog_out(DEBUG0) << "SensorID: " << sensorID;
  printf(" Z-position: %012.3f ", zPosition);
  streamlog_out(DEBUG0) << " Hit sigma x: " << sigmas(0)
       << " Hit sigma y: " << sigmas(1)
       << " Sigma scatter theta: " << sqrt(scatterThetaSqr)
       << " Excluded: " << excluded << endl;
}


TrackerSystem::TrackerSystem() : m_inited(false), m_maxCandidates(500), m_minClusterSize(3), m_nXdz(0.0f), m_nYdz(0.0) {;}

void TrackerSystem::setTruth(int plane, float x, float y, float xdz, float ydz){
  mcTruth.at(plane)->params(0) = x;
  mcTruth.at(plane)->params(1) = y;
  mcTruth.at(plane)->params(2) = xdz;
  mcTruth.at(plane)->params(3) = ydz;
}

void TrackerSystem::setMaxCandidates(int nCandidates){
  if( m_inited ){
    streamlog_out(ERROR5) << "ERROR: call to TrackerSystem::setMaxCandidates() must be called before call to TrackerSystem::init(). Quitting." << endl;
    exit(1);
  }
  m_maxCandidates = nCandidates;
}

void TrackerSystem::addPlane(int sensorID, float zPos, float sigmaX, float sigmaY, float scatterVariance, bool excluded){
  if(m_inited){
    streamlog_out(ERROR5) << "All planes must be added before call to TrackerSystem::init" << endl;
    exit(1);
  }

  FitPlane a(sensorID, zPos, sigmaX, sigmaY, scatterVariance, excluded);
  planes.push_back(a);
}

bool planeSort(FitPlane p1, FitPlane p2){ return( p1.getZpos() < p2.getZpos() );}

void TrackerSystem::init(){
  sort(planes.begin(), planes.end(), planeSort);
  mcTruth.resize(planes.size());
  for(int ii = 0; ii < static_cast< int >(planes.size()); ii++){
    planes.at(ii).print();
    mcTruth.at(ii) = new TrackEstimate();
  }
  m_fitter = new EigenFitter( planes.size() );
  tracks.resize( m_maxCandidates);
  for(size_t ii = 0; ii < m_maxCandidates; ii++){
    tracks.at(ii) = new TrackCandidate();
    tracks.at(ii)->init(planes.size());
  }
  m_inited = true;
}

void TrackerSystem::clear(){
  for(int ii = 0; ii < static_cast< int >(planes.size()); ii++){ planes.at(ii).clear(); }
  m_nTracks = 0;
}

void TrackerSystem::addMeasurement(size_t planeIndex, float x, float y, float z,  bool goodRegion, size_t measindex){
  planes.at(planeIndex).addMeasurement(x,y, z, goodRegion, measindex);
//  printf("addMeas: %5d at %5.2f %5.2f %5.2f good::%5d index:%5d \n", planeIndex, x, y, z, goodRegion, measindex);
}

void TrackerSystem::getChi2Kf(TrackCandidate *candidate){
  float chi2(0.0), ndof(0.0);
  Eigen::Vector2f chi2v;
  for(int plane = 0; plane < static_cast< int >(planes.size()); plane++){
    FitPlane& p = planes.at(plane);
    if(p.isExcluded()) { continue; }
    if(candidate->indexes.at(plane) < 0) { continue; }
    Measurement& m = p.meas.at( candidate->indexes.at(plane));
    ndof += 1.0;
    chi2v = (m.getM() - candidate->estimates.at(plane)->params.start<2>() ).cwise() / p.getSigmas();
    chi2v = chi2v.cwise().square();
    chi2 += chi2v.sum();
  }
  candidate->chi2 = chi2; candidate->ndof = (ndof * 2) - 4;
}

void TrackerSystem::getChi2Daf(TrackCandidate *candidate){
  float chi2(0.0), ndof(0.0);
  Vector2f chi2v;
  for(int plane = 0; plane <static_cast< int >(planes.size()); plane++){
    FitPlane& p = planes.at(plane);
    if(p.isExcluded()){continue;}
    for(size_t meas = 0; meas < p.meas.size(); meas++){
      Measurement& m = p.meas.at(meas);
      chi2v = (m.getM() - m_fitter->smoothed.at(plane)->params.start<2>()).cwise() / p.getSigmas();
      chi2v = chi2v.cwise().square();
      chi2 += p.weights(meas) * chi2v.sum();
      ndof += p.weights(meas);
    }
  }
  candidate->chi2 = chi2; candidate->ndof = (ndof * 2) - 4;
}

bool clusterSort(PlaneHit& a, PlaneHit& b) {
  //Sort by radius
  return( a.getM().squaredNorm() > b.getM().squaredNorm()  ); 
}

int TrackerSystem::addNeighbors(vector<PlaneHit> &candidate, list<PlaneHit> &hits){
  int counter(0); //How many hits are added to the cluster this iteration?
  for(list<PlaneHit>::iterator hit = hits.begin(); hit != hits.end(); hit++ ){
    for(vector<PlaneHit>::iterator it = candidate.begin(); it != candidate.end(); it++){
      Eigen::Vector2f resids = (*hit).getM() - (*it).getM();
     if(resids.squaredNorm() > m_sqrClusterRadius  ) { continue;}
//     if(counter>3)  printf("addNeigh %5d %8.3f gt/lt %8.3f \n", counter, resids.squaredNorm(), m_sqrClusterRadius );
      candidate.push_back((*hit));
      hit = hits.erase(hit);
      counter ++; break;
    }
  }
  return (counter);
}


void TrackerSystem::clusterTracker(){
  list<PlaneHit> availableHits;
  //Add all meas points to list
  for(size_t ii = 0; ii < planes.size(); ii++){
    if(planes.at(ii).isExcluded()) { continue;}
    float xShift = -1 * getNominalXdz() * planes.at(ii).getZpos();
    float yShift = -1 * getNominalYdz() * planes.at(ii).getZpos();
    for(size_t mm = 0; mm < planes.at(ii).meas.size(); mm++){
      PlaneHit a(planes.at(ii).meas.at(mm).getX() + xShift, planes.at(ii).meas.at(mm).getY() + yShift, ii, mm);
//    printf("clusterTracker:: ii=%5d mm=%5d X=%8.3f Y=%8.3f \n",ii,mm,planes.at(ii).meas.at(mm).getX() + xShift, planes.at(ii).meas.at(mm).getY() + yShift);
      availableHits.push_back( a );
    }
  }
  //Sort by radius from origin
  availableHits.sort(clusterSort);

  while( not availableHits.empty() ){
    vector<PlaneHit> candidate;
    candidate.push_back( availableHits.front() );
    availableHits.pop_front();
    while( addNeighbors( candidate , availableHits ) > 0) {;}
    //If we find enough hits, we make a candidate

    if(candidate.size() < getMinClusterSize() ){ continue; }
    if(m_nTracks >= m_maxCandidates) {
      streamlog_out(WARNING5) << "Maximum number of track candidates(" << m_maxCandidates 
		<< ") reached in DAF fitter! If this happens a lot, your configuration is probably off." 
		<< " If you are sure you config is right, see trackersystem.h on how to increase it." << std::endl;
      return;
    }
//printf("TrackerSystem::clusterTracker candidate size is OK\n");
    TrackCandidate* cnd = tracks.at(m_nTracks);
    for(size_t ii = 0; ii < planes.size(); ii++){
     if( planes.at(ii).meas.size() > 0 ) { 
       cnd->weights.at(ii).resize( planes.at(ii).meas.size());
       cnd->weights.at(ii).setZero();
     }
    }
    for(size_t ii = 0; ii < candidate.size(); ii++){
      PlaneHit& hit = candidate.at(ii);
      cnd->weights.at( hit.getPlane() )( hit.getIndex()) = 1.0;
    }
    m_nTracks++;
//printf("TrackerSystem::clusterTracker m_nTracks = %5d \n", m_nTracks);
  }
}

void TrackerSystem::fitPlanesInfo(TrackCandidate *candidate){
  //Biased fitter
  size_t nPlanes = planes.size();
  TrackEstimate* e = new TrackEstimate();
  e->cov.setZero();
  e->params.setZero();
  e->cov(2,2) = e->cov(3,3) = 1.0e-5f;
  //Forward fitter
  m_fitter->updateInfo( planes.at(0), candidate->indexes.at(0), e );
  m_fitter->forward.at(0)->copy(e);
  for(size_t ii = 1; ii < nPlanes ; ii++ ){
    m_fitter->predictInfo( planes.at( ii - 1), planes.at(ii), e );
    m_fitter->updateInfo( planes.at(ii), candidate->indexes.at(ii), e );
    m_fitter->forward.at(ii)->copy(e);
  }

  //Backward fitter, never bias
  e->cov.setZero();
  e->cov(2,2) = e->cov(3,3) = 1.0e-5f;
  e->params.setZero();
  m_fitter->backward.at( nPlanes -1 )->copy(e);
  m_fitter->updateInfo( planes.at(nPlanes -1 ), candidate->indexes.at(nPlanes -1), e );
  for(int ii = nPlanes -2; ii >= 0; ii-- ){
    m_fitter->predictInfo( planes.at( ii + 1 ), planes.at(ii), e );
    m_fitter->backward.at(ii)->copy(e);
    m_fitter->updateInfo( planes.at(ii), candidate->indexes.at(ii), e );
  }
  delete e;

  m_fitter->smoothInfo();

  for(int ii = 0; ii < static_cast< int >(planes.size()); ii++ ){
    candidate->estimates.at(ii)->copy( m_fitter->smoothed.at(ii) );
  }
  getChi2Kf(candidate);
}

void TrackerSystem::intersect(){
  for(int plane = 0; plane < static_cast< int >(planes.size()); plane++ ){
//    printf("intersect:: plane %5d of %5d \n", plane, planes.size() );
    FitPlane& pl = planes.at(plane);
    TrackEstimate* estim = m_fitter->smoothed.at(plane);
    //Line intersects with plane where
    // d = (p0 - l0) . n / ( l . n )
    // p0 is refpoint in plane
    Vector3f& refPoint = pl.getRef0();
    // n is vector normal of plane
    Vector3f& normVec = pl.getPlaneNorm();
    // l, point at line is defined by esimate x, y, and prev z
    Vector3f linePoint= Vector3f( estim->getX(), estim->getY(), pl.getMeasZ() );
//    printf(" refPoint: %5.3f %5.3f %5.3f \n", refPoint(0),refPoint(1), refPoint(2));
//    printf(" normVec: %5.3f %5.3f %5.3f \n", normVec(0), normVec(1), normVec(2));
//    printf(" linePoint: %5.3f %5.3f %5.3f \n", linePoint(0), linePoint(1), linePoint(2));
   // l, unit direction of line is normal of dx/dz, dy/dz, 1.0
    Vector3f lineDir( estim->getXdz(), estim->getYdz(), 1.0f);
    lineDir = lineDir.normalized();
//    printf(" lineDir  : %5.3f %5.3f %5.3f \n", lineDir(0), lineDir(1), lineDir(2));
    // p0 - l0
    Vector3f distance = refPoint - linePoint;
    float d = distance.dot(normVec) / lineDir.dot(normVec);
    //Propagate along track to track/plane intersection
    //estim->params(0) += d * lineDir(0);
    //estim->params(1) += d * lineDir(1);
    pl.setMeasZ( pl.getMeasZ() + d * lineDir(2));
/*
//// overwrite the line to plane intersection Z ::
        TVector3 hitInPlane( refPoint(0), refPoint(1), refPoint(2) ); // go back to mm
        TVector3 norm2Plane( normVec(0), normVec(1), normVec(2) );
        TVector3 lpoint( linePoint(0), linePoint(1), linePoint(2) );
        TVector3 lvector( lineDir(0), lineDir(1), lineDir(2) );
        TVector3 point( 1.,1.,1. );
          
        double linecoord_numenator   = norm2Plane.Dot( hitInPlane - lpoint );
        double linecoord_denumenator = norm2Plane.Dot( lvector );
//        printf("line: numenator: %8.3f denum: %8.3f [%5.3f %5.3f %5.3f::%5.3f %5.3f %5.3f]\n", linecoord_numenator, linecoord_denumenator, norm2Plane[0], norm2Plane[1], norm2Plane[2], lvector[0], lvector[1], lvector[2] );

        point = (linecoord_numenator/linecoord_denumenator)*lvector + lpoint ;
    printf("intersect: %5.2f %5.2f %5.2f \n", point(0), point(1), point(2) );
//    printf("plane mes: %5.2f %5.2f %5.2f \n", pl.getMeasX(), pl.getMeasY(), pl.getMeasZ() );
////
*/
  }
}

float TrackerSystem::runTweight(float t){
//printf("TrackerSystem::runTweight \n");
  m_fitter->setT(t);
  m_fitter->calculateWeights(planes, getDAFChi2Cut());
  float ndof = fitPlanesInfoDafInner();
  intersect();
  return( ndof );
}

void TrackerSystem::fitPlanesInfoDaf(TrackCandidate *candidate){
  float ndof = -4.0f;
//printf("TrackerSystem::fitPlanesInfoDaf\n");
  for(int plane = 0; plane < static_cast< int >(planes.size()); plane++ ){
    //Copy weights from candidate, get tot weight per plane

    planes.at(plane).weights.resize( candidate->weights.at(plane).size() );
    planes.at(plane).weights = candidate->weights.at(plane);
    if ( planes.at(plane).weights.size() > 0 ){
      planes.at(plane).setTotWeight( planes.at(plane).weights.sum() );
    } else {
      planes.at(plane).setTotWeight( 0.0f );
    }
    if( planes.at(plane).getTotWeight() > 1.0f){
      planes.at(plane).weights *= 1.0f / planes.at(plane).getTotWeight();
      planes.at(plane).setTotWeight( 1.0f );
    }
    ndof += planes.at(plane).getTotWeight() * 2.0;
//printf("plane %5d ndof:%8.3f  weight=%5.2f \n", plane, ndof, planes.at(plane).getTotWeight() );
  }
//printf("in mid of TrackerSystem::fitPlanesInfoDaf \n");
  if(ndof > 0.0f){ fitPlanesInfoDafInner();}
  // Temperatures should be given from top level. This should be fixed.
  // if(ndof > 0.0f) { ndof = runTweight(15.0); }
  // if(ndof > 0.0f) { ndof = runTweight(10.0); }
  // if(ndof > 0.0f) { ndof = runTweight(7.0); }
  // if(ndof > 0.0f) { ndof = runTweight(3.0); }
  if(ndof > 0.0f) { ndof = runTweight(1.2); }
  if(ndof > 0.0f) { ndof = runTweight(1.1); }
  if(ndof > 0.0f) { ndof = runTweight(1.0); }
  if(ndof > 0.0f) { ndof = runTweight(.1); }
  if(ndof > 0.0f) { ndof = runTweight(.1); }
//printf("and now ndof %5.3f\n",ndof);  
  if(ndof > 0.0f) {
    for(int ii = 0; ii <static_cast< int >(planes.size()); ii++ ){
      //Store estimates and weights in candidate
      candidate->estimates.at(ii)->copy( m_fitter->smoothed.at(ii) );
      candidate->weights.at(ii) = planes.at(ii).weights;
    }
    fitPlanesInfoDafBiased();
    getChi2Daf(candidate);
  } else {
//printf("and now candidate ndof %5.3f\n",ndof);  
    candidate->ndof = ndof;
    candidate->chi2 = 0;
  }
//printf("fitPlanesInfoDaf end \n");
}

void TrackerSystem::checkNan(TrackEstimate* e){
  if( isnan(e->params(0)) or
      isnan(e->params(1)) or
      isnan(e->params(2)) or
      isnan(e->params(3)) or
      isnan(e->cov(0,0)) or
      isnan(e->cov(1,1)) or
      isnan(e->cov(2,2)) or
      isnan(e->cov(3,3))){
    streamlog_out(ERROR5) << "Found nan!" << endl;
    streamlog_out(ERROR5) << e->params << endl;
    streamlog_out(ERROR5) << e->cov << endl;
    exit(1);
  }
}

float TrackerSystem::fitPlanesInfoDafInner(){
//printf("fitPlanesInfodafInner \n");

  size_t nPlanes = planes.size();
  TrackEstimate* e = new TrackEstimate();
  e->cov.setZero();
  e->params.setZero();
  //Forward fitter
  m_fitter->forward.at(0)->copy(e);
  m_fitter->updateInfoDaf( planes.at(0), e );
  float ndof( -1.0f * e->params.size());
  ndof += 2 * planes.at(0).getTotWeight();
//printf("ndof= - %8.3f +2* %8.3f \n", e->params.size(),  planes.at(0).getTotWeight() );
  for(size_t ii = 1; ii < nPlanes ; ii++ ){
    if(not planes.at(ii).isExcluded()){
      ndof += 2 * planes.at(ii).getTotWeight();
    }
//    printf("forw:[i=%5d] ndof += 2*%8.3f  => ndof=%8.3f \n", ii, planes.at(ii).getTotWeight(),ndof );
    m_fitter->predictInfo( planes.at( ii - 1), planes.at(ii), e );
    m_fitter->forward.at(ii)->copy(e);
    m_fitter->updateInfoDaf( planes.at(ii), e );
  }
//printf("ndof %5.2f <? 2.5 [return?]\n", ndof);
  //No reason to complete
  if(ndof < 2.5) { delete e; return(ndof);}
  
  //Backward fitter, never bias
  e->cov.setZero();
  e->params.setZero();
  m_fitter->backward.at( nPlanes -1 )->copy(e);
  m_fitter->updateInfoDaf( planes.at(nPlanes -1 ), e );
  for(int ii = nPlanes -2; ii >= 0; ii-- ){
//    printf("back:[i=%5d] totWeight:%8.3f   \n", ii, planes.at(ii).getTotWeight() );
    m_fitter->predictInfo( planes.at( ii + 1 ), planes.at(ii), e );
    m_fitter->backward.at(ii)->copy(e);
    m_fitter->updateInfoDaf( planes.at(ii), e );
  }
  delete e;

//  printf("returning ndof=%8.3f \n", ndof);

  m_fitter->smoothInfo();
  return(ndof);
}

float TrackerSystem::fitPlanesInfoDafBiased(){
//printf("TrackerSystem::fitPlanesInfoDafBiased\n");

  size_t nPlanes = planes.size();
  TrackEstimate* e = new TrackEstimate();
  e->cov.setZero();
  e->params.setZero();
  //Forward fitter
  m_fitter->updateInfoDaf( planes.at(0), e );
  m_fitter->forward.at(0)->copy(e);
  float ndof( -1.0f * e->params.size());
  ndof += 2 * planes.at(0).getTotWeight();
  for(size_t ii = 1; ii < nPlanes ; ii++ ){
    ndof += 2 * planes.at(ii).getTotWeight();
    m_fitter->predictInfo( planes.at( ii - 1), planes.at(ii), e );
    m_fitter->updateInfoDaf( planes.at(ii), e );
    m_fitter->forward.at(ii)->copy(e);
//    printf("forward: m_fitter: %5d  ndof=%5.2f \n", ii, ndof); 
  }
  //No reason to complete
  if(ndof < 2.5) { delete e; return(ndof);}
  
  //Backward fitter, never bias
  e->cov.setZero();
  e->params.setZero();
  m_fitter->backward.at( nPlanes -1 )->copy(e);
  m_fitter->updateInfoDaf( planes.at(nPlanes -1 ), e );
  for(int ii = nPlanes -2; ii >= 0; ii-- ){
    m_fitter->predictInfo( planes.at( ii + 1 ), planes.at(ii), e );
    m_fitter->backward.at(ii)->copy(e);
    m_fitter->updateInfoDaf( planes.at(ii), e );
//    printf("backward: m_fitter: %5d  ndof=%5.2f \n", ii, ndof); 
  }
  delete e;

  m_fitter->smoothInfo();
  return(ndof);
}


void TrackerSystem::truthTracker(){
  TrackCandidate *candidate = tracks.at(0);
  for(int ii = 0 ; ii < static_cast< int >(planes.size()); ii++){
    candidate->weights.at(ii).resize( planes.at(ii).meas.size());
    candidate->weights.at(ii).setZero();
    if( (not planes.at(ii).isExcluded()) and (planes.at(ii).meas.size() > 0)){
      candidate->weights.at(ii)(0) = 1.0;
      candidate->indexes.at(ii) = 0;
    }
    else{
      candidate->indexes.at(ii) = -1;
    }
 }
  m_nTracks++;
}

void TrackerSystem::weightToIndex(TrackCandidate* cnd){
  //Prepare a DAF fitted track for standard KF
  for(size_t plane = 0 ; plane < planes.size(); plane++){
    int index(-1);
    float maxWeight(0.5f);
    if( not planes.at(plane).isExcluded() ){
      for(size_t m = 0; m < planes.at(plane).meas.size(); m++){
	if( cnd->weights.at(plane)(m) < maxWeight){ continue; }
	maxWeight = cnd->weights.at(plane)(m);
	index = m;
      }
    }
    cnd->indexes.at(plane) = index;
  }
}
