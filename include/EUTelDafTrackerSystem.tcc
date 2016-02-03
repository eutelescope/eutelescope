#include <iostream>
#include <algorithm>
#include <Eigen/Core>

using namespace std;
using namespace daffitter;

namespace daffitter{
  template <typename T, size_t N>
  inline void TrackEstimate<T, N>::makeSeed(bool keepState){
    // Make a seed for a Kalman filter in the standard formulation

    //Initial state of the estimate,
    if(not keepState){ params.setZero(); }
    //Initial state of Cov
    cov.setZero();
    T posVar = 1e6f;
    T angleVar = 1e1f;
    cov(0,0) = posVar; cov(1,1) = posVar;;
    cov(2,2) = angleVar; cov(3,3) = angleVar;
  }

  template<typename T, size_t N>
  inline bool TrackEstimate<T, N>::isSeed(){
    //Check if the estimate is a seed.
    return( (cov(0, 0) > 1e5f ) and 
	    (cov(1, 1) > 1e5f ) ); 
  }
  
  template<typename T, size_t N>
  void TrackEstimate<T, N>::print(){
    //Print info about FitPlane
    if( cov.fullPivLu().isInvertible() ){
      std::cout << "Cov: " << "Invertible!" << std::endl;
    } else {
      std::cout << "Cov: " << "non-inverible! :(" << std::endl;
    }
    std::cout << "Params: " << std::endl << cov.inverse() * params << std::endl;
      std::cout << "Cov: " << std::endl << cov.inverse() << std::endl;
      std::cout << std::endl << params << std::endl << std::endl << cov << std::endl;
  }

  
  template <typename T, size_t N>
  inline void TrackEstimate<T, N>::makeSeedInfo(){
    // Make a seed for the Kalman filter in the information filter formulation

    //Initial state of the estimate,
    params.setZero();
    //Initial state of Cov
    cov.setZero();
  }
}

template<typename T, size_t N>
void TrackCandidate<T,N>::print(){
  //Print info on TrackCandidate
  std::cout << "Track candidate with " << indexes.size() << " planes:" << std::endl;
  for(size_t ii = 0; ii < indexes.size(); ii++){
    std::cout << "Pl: " << ii << " Index " << indexes.at(ii) << " Weights: ";
    for(int jj = 0; jj < weights.at(ii).size(); jj++){
      std::cout << weights.at(ii)(jj) << ", " ;
    }
    std::cout << std::endl;
  }
}

template<typename T, size_t N>
void TrackCandidate<T,N>::init(int nPlanes){
  //Initialize track candidate for nPlanes planes.
  indexes.resize(nPlanes);
  weights.resize(nPlanes);
  estimates.resize(nPlanes);
}

template<typename T, size_t N>
TrackCandidate<T,N>::TrackCandidate(int nPlanes){
  init(nPlanes);
}

template<typename T>
Measurement<T>::Measurement(T x, T y, T z, bool goodRegion, size_t iden):
  m_goodRegion(goodRegion), zPos(z), m_iden(iden) {
  m(0) = x; m(1) = y;
}

template<typename T>
FitPlane<T>::FitPlane(int sensorID, T zPos, T sigmaX, T sigmaY, T scatterThetaSqr, bool excluded):
  sensorID(sensorID), scatterThetaSqr(scatterThetaSqr), excluded(excluded), zPosition(zPos){
  //Constructor for a telescope or material plane.
  
  sigmas(0) = sigmaX; sigmas(1) = sigmaY;
  variances = sigmas.array().square();
  invMeasVar(0) = 1/(sigmaX * sigmaX);
  invMeasVar(1) = 1/(sigmaY * sigmaY);
  ref0(0) = 0.0f; ref0(1) = 0.0f; ref0(2) = zPos;
  norm(0) = 0.0f; norm(1) = 0.0f; norm(2) = 1.0f;
}

template<typename T>
void FitPlane<T>::print(){
  //Pretty print initialized plane info
  cout << "SensorID: " << sensorID
       << " Z-position: " << zPosition
       << " Hit sigma x: " << sigmas(0)
       << " Hit sigma y: " << sigmas(1)
       << " Sigma scatter theta: " << sqrt(scatterThetaSqr)
       << " Excluded: " << excluded << endl;
}

template<typename T>
void FitPlane<T>::scaleErrors(T scaleX, T scaleY){
  //Scale the sigmas
  sigmas(0) *= scaleX; sigmas(1) *= scaleY;
  invMeasVar(0) = 1.0f / ( sigmas(0) * sigmas(0));
  invMeasVar(1) = 1.0f / ( sigmas(1) * sigmas(1));
}

template<typename T>
void FitPlane<T>::setSigmas(T sigmaX, T sigmaY){
  //Set both measurement sigmas (std deviations)
  sigmas(0) = sigmaX; sigmas(1) = sigmaY;
  invMeasVar(0) = 1.0f / ( sigmas(0) * sigmas(0));
  invMeasVar(1) = 1.0f / ( sigmas(1) * sigmas(1));
}

template<typename T>
void FitPlane<T>::setSigmaX(T sigmaX){
  //Set measurement sigma in X direction (std deviations)
  sigmas(0) = sigmaX; invMeasVar(0) = 1.0/(sigmaX * sigmaX);
}

template<typename T>
void FitPlane<T>::setSigmaY(T sigmaY){
  //Set measurement sigma in Y direction (std deviations)
  sigmas(1) = sigmaY; invMeasVar(1) = 1.0/(sigmaY * sigmaY);
}

template<typename T>
void PlaneHit<T>::print() {
  //Print info on PlaneHit
  std::cout << "Plane " << plane << ", index " << index << ", meas " << std::endl << xy << std::endl; 
}

template <typename T, size_t N>
TrackerSystem<T, N>::TrackerSystem() : m_inited(false), m_maxCandidates(100), m_minClusterSize(3), m_nXdz(0.0f), m_nYdz(0.0),
				       m_nXdzdeviance(0.01),m_nYdzdeviance(0.01), m_skipMax(2) {
  //Constructor for the system of detector planes.
}

template <typename T, size_t N>
TrackerSystem<T, N>::TrackerSystem(const TrackerSystem<T,N>& sys) : m_inited(false), m_maxCandidates(sys.m_maxCandidates), 
								    m_minClusterSize(sys.m_minClusterSize), 
								    m_nXdz(sys.m_nXdz), m_nYdz(sys.m_nYdz),
								    m_nXdzdeviance(sys.m_nXdzdeviance), m_nYdzdeviance(sys.m_nYdzdeviance),
								    m_dafChi2(sys.m_dafChi2), m_ckfChi2(sys.m_ckfChi2), 
								    m_chi2OverNdof(sys.m_chi2OverNdof), m_sqrClusterRadius(sys.m_sqrClusterRadius),
								    m_skipMax(sys.m_skipMax){
  //Copy constructor. Copy relevant info from sys, add planes and init.
  for(size_t ii = 0; ii < sys.planes.size(); ii++){
    //const FitPlane<T>& pl = sys.planes.at(ii);
    const FitPlane<T>& pl = sys.planes.at(ii);
    this->addPlane(pl.getSensorID(), pl.getZpos(), pl.getSigmaX(), pl.getSigmaY(), pl.getScatterThetaSqr(), pl.isExcluded());
  }
  this->init(true);
}

template <typename T,size_t N>
inline Eigen::Matrix<T, 2, 1> TrackerSystem<T, N>::getBiasedResidualErrors(FitPlane<T>& pl, TrackEstimate<T, N>& est){
  //Get covariance matrix diagonal vector for biased residuals: distance from filtered/updated estimate to measurement
  return(pl.getSigmas().array().square() - est.cov.diagonal().head(2).array());
}

template <typename T,size_t N>
inline Eigen::Matrix<T, 2, 1> TrackerSystem<T, N>::getUnBiasedResidualErrors(FitPlane<T>& pl, TrackEstimate<T, N>& est){
  //Get covariance matrix diagonal vector for unbiased residuals: distance from unfiltered/predicted estimate to measurement
  return(pl.getSigmas().array().square() + est.cov.diagonal().head(2).array());
}

template <typename T,size_t N>
inline Eigen::Matrix<T, 2, 1> TrackerSystem<T, N>::getResiduals(Measurement<T>& m, TrackEstimate<T, N>& est){
  //Vector of distances from estimate to measurementx
  return(m.getM() - est.params.head(2));
}

template <typename T,size_t N>
void TrackerSystem<T, N>::setMaxCandidates(int nCandidates){
  //Set the maximum allowed number of track candidates.
  if( m_inited ){
    cerr << "ERROR: call to TrackerSystem::setMaxCandidates() must be called before call to TrackerSystem::init(). Quitting." << endl;
    exit(1);
  }
  m_maxCandidates = nCandidates;
}

template <typename T,size_t N>
void TrackerSystem<T, N>::addPlane(int sensorID, T zPos, T sigmaX, T sigmaY, T scatterVariance, bool excluded){
  //Add a detector or material plane to the tracker system
  if(m_inited){
    cerr << "All planes must be added before call to TrackerSystem::init" << endl;
    exit(1);
  }

  FitPlane<T>  a(sensorID, zPos, sigmaX, sigmaY, scatterVariance, excluded);
  planes.push_back(a);
}

template<typename T>
inline bool planeSort(FitPlane<T>  p1, FitPlane<T>  p2){ return( p1.getZpos() < p2.getZpos() );}

template <typename T,size_t N>
void TrackerSystem<T, N>::init(bool quiet){
  //Initialize the tracker system:
  //  - All planes are sorted by z position, information about the planes are printed to screen.
  //  - Memory is allocated for track candidates and MC truth
  sort(planes.begin(), planes.end(), planeSort<T>);
  if(not quiet){
    for(int ii = 0; ii < (int) planes.size(); ii++){
      planes.at(ii).print();
    }
  }
  m_fitter.init(planes.size());
  m_inited = true;
}

template <typename T,size_t N>
void TrackerSystem<T, N>::clear(){
  // Prepare tracker system for a new event.
  for(int ii = 0; ii < (int)planes.size(); ii++){ planes.at(ii).clear(); }
  tracks.clear();
  m_nTracks = 0;
}

template <typename T,size_t N>
inline void TrackerSystem<T, N>::addMeasurement(size_t planeIndex, T x, T y, T z,  bool goodRegion, size_t measiden){
  // Add a measurement to the tracker system
  planes.at(planeIndex).addMeasurement(x,y, z, goodRegion, measiden);
}

template <typename T,size_t N>
inline void TrackerSystem<T, N>::addMeasurement(Measurement<T>& meas){
  // Add a measurement to the tracker system
  for(size_t ii = 0; ii < planes.size(); ii++){
    if( meas.getIden() == planes.at(ii).getSensorID()){
      planes.at(ii).addMeasurement(meas);
      break;
    }
  }
}

template <typename T>
inline bool clusterSort(PlaneHit<T>& a, PlaneHit<T>& b) {
  //Sort by radius
  return( a.getM().squaredNorm() > b.getM().squaredNorm()  ); 
}

template <typename T,size_t N>
inline int TrackerSystem<T, N>::addNeighbors(vector<PlaneHit<T> > &candidate, list<PlaneHit<T> > &hits){
  // Part of the cluster tracker. Checks distances to measurements, and adds measurement to cluster maybe.
  int counter(0); //How many hits are added to the cluster this iteration?
  for(typename list<PlaneHit<T> >::iterator hit = hits.begin(); hit != hits.end(); hit++ ){
    for(typename vector<PlaneHit<T> >::iterator cand = candidate.begin(); cand != candidate.end(); cand++){
      Eigen::Matrix<T, 2, 1> resids = (*hit).getM() - (*cand).getM();
      if(resids.squaredNorm() > m_sqrClusterRadius ) { continue;}
      candidate.push_back((*hit));
      hit = hits.erase(hit);
      counter ++;
      break;
    }
  }
  return (counter);
}

template <typename T,size_t N>
void TrackerSystem<T, N>::index0tracker(){
  //Create a track candidate, ehere every plane has a hit with index 0. Used by EstMat.
  TrackCandidate<T, N> cnd(planes.size());
  cnd.ndof = 0;
  cnd.chi2 = 0.0f;
  for(size_t ii = 0; ii < planes.size(); ii++){
    cnd.indexes[ii] = 0;
  }
  tracks.push_back(cnd);
  m_nTracks++;
}

template <typename T,size_t N>
void TrackerSystem<T, N>::clusterTracker(){
  //A track fitter that propagates measurements into z = 0, then assumes measurement clusters are track candidates.
  list<PlaneHit<T> > availableHits;
  //Add all meas points to list
  for(size_t ii = 0; ii < planes.size(); ii++){
    if(planes.at(ii).isExcluded()) { continue;}
    T xShift = -1 * getNominalXdz() * planes.at(ii).getZpos();
    T yShift = -1 * getNominalYdz() * planes.at(ii).getZpos();
    for(size_t mm = 0; mm < planes.at(ii).meas.size(); mm++){
      PlaneHit<T> a(planes.at(ii).meas.at(mm).getX() + xShift, planes.at(ii).meas.at(mm).getY() + yShift, ii, mm);
      availableHits.push_back( a );
    }
  }
  //Sort by radius from origin
  //availableHits.sort(clusterSort<T>);
  while( not availableHits.empty() ){
    vector<PlaneHit<T> > candidate;
    candidate.push_back( availableHits.front() );
    availableHits.pop_front();
    while( addNeighbors( candidate , availableHits ) > 0) {;}
    //If we find enough hits, we make a candidate

    if(candidate.size() < getMinClusterSize() ){ continue; }
    if(m_nTracks >= m_maxCandidates) {
      std::cout << "Maximum number of track candidates(" << m_maxCandidates 
		<< ") reached in DAF fitter! If this happens a lot, your configuration is probably off." 
		<< " If you are sure you config is right, see trackersystem.h on how to increase it." << std::endl;
      return;
    }

    TrackCandidate<T,N> cnd(planes.size());

    cnd.ndof = 0;
    cnd.chi2 = 0;
    for(size_t ii = 0; ii < planes.size(); ii++){
      cnd.weights.at(ii).resize( planes.at(ii).meas.size());
      if( planes.at(ii).meas.size() > 0 ) { 
	cnd.weights.at(ii).setZero();
      }
    }
    for(size_t ii = 0; ii < candidate.size(); ii++){
      PlaneHit<T>& hit = candidate.at(ii);
      cnd.weights.at( hit.getPlane() )( hit.getIndex()) = 1.0;
    }
    tracks.push_back(cnd);
    m_nTracks++;
  }
}

template <typename T,size_t N>
void TrackerSystem<T, N>::fitInfoFWUnBiased(TrackCandidate<T, N>& candidate){
  //Fit planes in the FW direction (increasing z position). Unbiased means the predictions are saved
  size_t nPlanes = planes.size();
  TrackEstimate<T,N> e;
  e.makeSeedInfo();
  //Forward fitter, scattering not included at plane, biased
  m_fitter.forward.at(0) = e;
  m_fitter.updateInfo( planes.at(0), candidate.indexes.at(0), e );
  for(size_t ii = 1; ii < nPlanes ; ii++ ){
    m_fitter.predictInfo( planes.at( ii - 1), planes.at(ii), e );
    m_fitter.forward.at(ii) = e;
    m_fitter.updateInfo( planes.at(ii), candidate.indexes.at(ii), e );
    m_fitter.addScatteringInfo( planes.at(ii), e);
  }
}

template <typename T,size_t N>
void TrackerSystem<T, N>::fitInfoFWBiased(TrackCandidate<T, N>& candidate){
  //Fit planes in the FW direction (increasing z position). Biased means the updated estimates are saved
  size_t nPlanes = planes.size();
  TrackEstimate<T,N> e;
  e.makeSeedInfo();

  //Forward fitter, scattering not included at plane, biased
  m_fitter.updateInfo( planes.at(0), candidate.indexes.at(0), e );
  m_fitter.forward.at(0) = e;
  //m_fitter.addScatteringInfo( planes.at(0), e);
  for(size_t ii = 1; ii < nPlanes ; ii++ ){
    m_fitter.predictInfo( planes.at( ii - 1), planes.at(ii), e );
    m_fitter.updateInfo( planes.at(ii), candidate.indexes.at(ii), e );
    m_fitter.forward.at(ii) = e;
    m_fitter.addScatteringInfo( planes.at(ii), e);
  }
}

template <typename T,size_t N>
void TrackerSystem<T, N>::fitInfoBWUnBiased(TrackCandidate<T, N>& candidate){
  //Fit planes in the BW direction (decreasing z position). Unbiased means the predictions are saved
  size_t nPlanes = planes.size();
  TrackEstimate<T,N> e;
  e.makeSeedInfo();
  
  //Backward fitter, scattering included at plane, unbiased
  m_fitter.backward.at( nPlanes -1 ) = e;
  m_fitter.updateInfo( planes.at(nPlanes -1 ), candidate.indexes.at(nPlanes -1), e );
  for(int ii = nPlanes -2; ii >= 0; ii-- ){
    m_fitter.predictInfo( planes.at( ii + 1 ), planes.at(ii), e );
    m_fitter.addScatteringInfo( planes.at(ii), e);
    m_fitter.backward.at(ii) = e;
    m_fitter.updateInfo( planes.at(ii), candidate.indexes.at(ii), e );
  }
}

template <typename T,size_t N>
void TrackerSystem<T, N>::fitInfoBWBiased(TrackCandidate<T, N>& candidate){
  //Fit planes in the BW direction (decreasing z position). Biased means the predictions are saved
  size_t nPlanes = planes.size();
  TrackEstimate<T,N> e;
  e.makeSeedInfo();

  //Backward fitter, scattering included at plane, biased
  m_fitter.updateInfo( planes.at(nPlanes -1 ), candidate.indexes.at(nPlanes -1), e );
  m_fitter.backward.at( nPlanes -1 ) = e;
  for(int ii = nPlanes -2; ii >= 0; ii-- ){
    m_fitter.predictInfo( planes.at( ii + 1 ), planes.at(ii), e );
    m_fitter.addScatteringInfo( planes.at(ii), e);
    m_fitter.updateInfo( planes.at(ii), candidate.indexes.at(ii), e );
    m_fitter.backward.at(ii) = e;
  }
}

template <typename T,size_t N>
void TrackerSystem<T, N>::smoothInfo(TrackCandidate<T, N>& candidate){
  // Run a smoother for the information filter
  m_fitter.smoothInfo();
  intersect();
  for(int ii = 0; ii < (int) planes.size() ; ii++ ){
    candidate.estimates.at(ii) = m_fitter.smoothed.at(ii);
  }
}

template <typename T,size_t N>
void TrackerSystem<T, N>::fitPlanesInfoBiased(TrackCandidate<T, N>& candidate){
  // Get smoothed,biased, estimates for all planes.  Also get the chi2, ndof
  fitInfoFWBiased(candidate);
  fitInfoBWUnBiased(candidate);
  smoothInfo(candidate);
  getChi2BiasedInfo(candidate);
}

template <typename T,size_t N>
void TrackerSystem<T, N>::fitPlanesInfoUnBiased(TrackCandidate<T, N>& candidate){
  // Get smoothed, unbiased, estimates for all planes. Also get the chi2, ndof
  fitInfoFWUnBiased(candidate);
  fitInfoBWUnBiased(candidate);
  smoothInfo(candidate);
  getChi2UnBiasedInfo(candidate);
}

template <typename T,size_t N>
void TrackerSystem<T, N>::getChi2Kf(TrackCandidate<T, N>& candidate){
  //Get the chi2 for a standard formulated, unbiased KF
  T chi2(0.0), ndof(0.0);
  Eigen::Matrix<T, 2, 1> resv, errv;
  for(int plane = 0; plane < (int) planes.size(); plane++){
    FitPlane<T>& p = planes.at(plane);
    if(p.isExcluded()) { continue; }
    if(candidate.indexes.at(plane) < 0) { continue; }
    if(ndof < 1.5) { ndof+= 1.0; continue;}
    Measurement<T>& m = p.meas.at( candidate.indexes.at(plane));
    TrackEstimate<T,N>& est = m_fitter.forward.at(plane);
    resv = getResiduals(m, est).array().square();
    errv = getUnBiasedResidualErrors(p, est);
    chi2 += (resv.array() / errv.array()).sum();
    ndof += 1.0;
  }
  candidate.chi2 = chi2; 
  candidate.ndof = (ndof * 2) - 4;
}

template <typename T,size_t N>
void TrackerSystem<T, N>::getChi2BiasedInfo(TrackCandidate<T, N>& candidate){
  //Get the chi2 for a biased information filter
  Eigen::Matrix<T, N, N> tmp4x4;
  T ndof(0), chi2(0.0);
  Eigen::Matrix<T, 2, 1> errv, resv;
  TrackEstimate<T,N> estim;

  for(size_t ii = 0; ii < planes.size(); ii++){
    int index = candidate.indexes.at(ii); 
    if( index < 0) { continue;}
    if(planes.at(ii).isExcluded()) { continue;}
    //Chi2 increment of first is 0 and non invertible
    if(ndof < 1.5){ ndof += 1.0; continue;}
    ndof += 1.0;

    Measurement<T>& m = planes.at(ii).meas.at(index);
    estim = m_fitter.forward.at(ii);
    fastInvert(estim.cov);
    estim.params = estim.cov * estim.params;

    resv = getResiduals( m, estim ).array().square();
    errv = getBiasedResidualErrors( planes.at(ii), estim);
    chi2 += (resv.array() / errv.array()).sum();
  }
  candidate.chi2 = chi2;
  candidate.ndof = ndof * 2 - 4;
}

template <typename T,size_t N>
void TrackerSystem<T, N>::getChi2UnBiasedInfo(TrackCandidate<T, N>& candidate){
  //Get the chi2 for a unbiased information filter
  Eigen::Matrix<T, N, N> tmp4x4;
  T ndof(0), chi2(0.0);
  Eigen::Matrix<T, 2, 1> errv, resv;
  TrackEstimate<T,N> estim;

  for(size_t ii = 0; ii < planes.size(); ii++){
    int index = candidate.indexes.at(ii);
    if( index < 0) { continue;}
    if(planes.at(ii).isExcluded()) { continue;}
    //Chi2 increment of first is 0 and non invertible
    if(ndof < 1.5){ ndof += 1.0; continue;}
    ndof += 1.0;

    Measurement<T>& m = planes.at(ii).meas.at(index);
    estim = m_fitter.forward.at(ii);
    fastInvert(estim.cov);
    estim.params = estim.cov * estim.params;

    resv = getResiduals( m, estim ).array().square();
    errv = getUnBiasedResidualErrors( planes.at(ii), estim);
    chi2 += (resv.array() / errv.array()).sum();
  }
  candidate.chi2 = chi2;
  candidate.ndof = ndof * 2 - 4;
}

template <typename T,size_t N>
void TrackerSystem<T, N>::getChi2UnBiasedInfoDaf(TrackCandidate<T, N>& candidate){
  //Get the chi2 for a unbiased information filter with weighted measurements.
  Eigen::Matrix<T, N, N> tmp4x4;
  T ndof(0.0), chi2(0.0);
  T sumWeight(0.0);
  Eigen::Matrix<T, 2, 1> errv, resv;
  TrackEstimate<T,N> estim;

  for(size_t ii = 0; ii < planes.size(); ii++){
    if(planes.at(ii).isExcluded()) { continue;}
    //Chi2 increment of first is 0 and non invertible
    double newWeight = candidate.weights.at(ii).sum();
    sumWeight += newWeight;
    if(newWeight < 0.05) { continue; }
    if(ndof < 1.5){ ndof += 1.0; continue;}
    ndof += 1.0;
    estim = m_fitter.forward.at(ii);
    fastInvert(estim.cov);
    estim.params = estim.cov * estim.params;
    for(size_t mm = 0; mm < planes.at(ii).meas.size(); mm++){
      Measurement<T>& m = planes.at(ii).meas.at(mm);
      T weight = candidate.weights.at(ii)(mm);
      resv = getResiduals( m, estim ).array().square();
      errv = getUnBiasedResidualErrors( planes.at(ii), estim);
      chi2 += weight * (resv.array() / errv.array()).sum();
    }
  }
  candidate.chi2 = chi2;
  candidate.ndof = sumWeight * 2 - 4;
  if(isnan(candidate.chi2)){ cout << "NAN CHI2" << endl << endl;}
}

template <typename T,size_t N>
void TrackerSystem<T, N>::intersect(){
  //Check where a straight line track intersects with a measurement plane, update the measurement z position
  for(int plane = 0; plane < (int) planes.size(); plane++ ){
    FitPlane<T>& pl = planes.at(plane);
    TrackEstimate<T,N>& estim = m_fitter.smoothed.at(plane);
    //Line intersects with plane where
    // d = (p0 - l0) . n / ( l . n )
    // p0 is refpoint in plane
    Eigen::Matrix<T, 3, 1>& refPoint = pl.getRef0();
    // n is vector normal of plane
    Eigen::Matrix<T, 3, 1>& normVec = pl.getPlaneNorm();
    // l, point at line is defined by esimate x, y, and prev z
    Eigen::Matrix<T, 3, 1> linePoint( estim.getX(), estim.getY(), pl.getMeasZ() );
    // l, unit direction of line is normal of dx/dz, dy/dz, 1.0
    Eigen::Matrix<T, 3, 1> lineDir( estim.getXdz(), estim.getYdz(), 1.0f);
    lineDir = lineDir.normalized();
    // p0 - l0
    Eigen::Matrix<T, 3, 1> distance = refPoint - linePoint;
    T d = normVec.dot(distance) / normVec.dot(lineDir);
    //Propagate along track to track/plane intersection
    pl.setMeasZ( pl.getMeasZ() + d * lineDir(2));
  }
}

template <typename T,size_t N>
T TrackerSystem<T, N>::runTweight(T t, TrackCandidate<T, N>& candidate){
  //A DAF iteration with temperature t.
  m_fitter.setT(t);
  m_fitter.calculateWeights(planes, getDAFChi2Cut(), candidate.weights);
  T ndof = fitPlanesInfoDafInner(candidate);
  intersect();
  return( ndof );
}

template <typename T,size_t N>
void TrackerSystem<T, N>::fitPlanesInfoDaf(TrackCandidate<T, N>& candidate){
  // Get smoothed estimates for all planes using the unbiased DAF
  T ndof = -4.0f;
  for(int plane = 0; plane < (int) planes.size(); plane++ ){
    //set tot weight per plane
    if ( candidate.weights.at(plane).size() > 0 ){
      planes.at(plane).setTotWeight( candidate.weights.at(plane).sum() );
    } else {
      planes.at(plane).setTotWeight( 0.0f );
    }
    if( planes.at(plane).getTotWeight() > 1.0f){
      candidate.weights.at(plane) *= 1.0f / planes.at(plane).getTotWeight();
      planes.at(plane).setTotWeight( 1.0f );
    }
    ndof += planes.at(plane).getTotWeight() * 2.0;
  }
  fitPlanesInfoDafInner(candidate);
  if(isnan(ndof)) { ndof = -10.0; }

  // Running with fixed annealing schedule.
  if(ndof > -1.0f) { ndof = runTweight(25.0, candidate); }
  if(ndof > -1.0f) { ndof = runTweight(20.0, candidate); }
  if(ndof > -1.9f) { ndof = runTweight(14.0, candidate); }
  if(ndof > -1.9f) { ndof = runTweight( 8.0, candidate); }
  if(ndof > -1.9f) { ndof = runTweight( 4.0, candidate); }
  if(ndof > -1.9f) { ndof = runTweight( 1.0, candidate); }
  
  if(ndof > -1.9f) {
    for(int ii = 0; ii <(int)  planes.size() ; ii++ ){
      //Store estimates and weights in candidate
      candidate.estimates.at(ii) = m_fitter.smoothed.at(ii);
    }
    getChi2UnBiasedInfoDaf(candidate);
    weightToIndex(candidate);
  } else{
    candidate.ndof = ndof;
    candidate.chi2 = 0;
  }
}

template <typename T,size_t N>
void TrackerSystem<T, N>::checkNan(TrackEstimate<T, N>& e){
  //See if there are any nans in the estimate. For debugging numerical problems.
  if( isnan(e.params(0)) or
      isnan(e.params(1)) or
      isnan(e.params(2)) or
      isnan(e.params(3)) or
      isnan(e.cov(0,0)) or
      isnan(e.cov(1,1)) or
      isnan(e.cov(2,2)) or
      isnan(e.cov(3,3))){
    cout << "Found nan!" << endl;
    cout << e.params << endl;
    cout << e.cov << endl;
    exit(1);
  }
}

template <typename T,size_t N>
T TrackerSystem<T, N>::fitPlanesInfoDafInner(TrackCandidate<T, N>& candidate){
  // Get smoothed estimates for all planes using the weighted information filter
  size_t nPlanes = planes.size();
  TrackEstimate<T,N> e;
  e.makeSeedInfo();

  //Forward fitter
  m_fitter.forward.at(0) = e;
  m_fitter.updateInfoDaf( planes.at(0), e, candidate.weights.at(0));
  T ndof( -1.0f * e.params.size());
  ndof += 2 * planes.at(0).getTotWeight();
  for(size_t ii = 1; ii < nPlanes ; ii++ ){
    if(not planes.at(ii).isExcluded()){
      ndof += 2 * planes.at(ii).getTotWeight();
    }
    m_fitter.predictInfo( planes.at( ii - 1), planes.at(ii), e );
    m_fitter.forward.at(ii) = e;
    m_fitter.updateInfoDaf( planes.at(ii), e, candidate.weights.at(ii) );
    m_fitter.addScatteringInfo( planes.at(ii), e);
  }
  //No reason to complete unless >1 measurements are in
  if(ndof < -2.1) { return(ndof);}
  
  //Backward fitter, never bias
  e.makeSeedInfo();
  m_fitter.backward.at( nPlanes -1 ) = e;
  m_fitter.updateInfoDaf( planes.at(nPlanes -1 ), e, candidate.weights.at(nPlanes - 1) );
  for(int ii = nPlanes -2; ii >= 0; ii-- ){
    m_fitter.predictInfo( planes.at( ii + 1 ), planes.at(ii), e );
    m_fitter.addScatteringInfo( planes.at(ii), e);
    m_fitter.backward.at(ii) = e;
    m_fitter.updateInfoDaf( planes.at(ii), e, candidate.weights.at(ii) );
  }

  m_fitter.smoothInfo();
  return(ndof);
}

template <typename T,size_t N>
T TrackerSystem<T, N>::fitPlanesInfoDafBiased(TrackCandidate<T, N>& candidate){
  // Get smoothed estimates for all planes using the biased DAF

  size_t nPlanes = planes.size();
  TrackEstimate<T,N> e;
  e.makeSeedInfo();

  //Forward fitter
  m_fitter.updateInfoDaf( planes.at(0), e, candidate.weights.at(0));
  m_fitter.forward.at(0) = e;
  T ndof( -1.0f * e.params.size());
  ndof += 2 * planes.at(0).getTotWeight();
  for(size_t ii = 1; ii < nPlanes ; ii++ ){
    ndof += 2 * planes.at(ii).getTotWeight();
    m_fitter.predictInfo( planes.at( ii - 1), planes.at(ii), e );
    m_fitter.addScatteringInfo( planes.at(ii), e);
    m_fitter.updateInfoDaf( planes.at(ii), e, candidate.weights.at(ii));
    m_fitter.forward.at(ii) = e;
  }
  //No reason to complete if less than one measurement is in
  if(ndof < -2.1) { return(ndof);}
  
  //Backward fitter, never bias
  e.makeSeedInfo();
  m_fitter.backward.at( nPlanes -1 ) = e;
  m_fitter.updateInfoDaf( planes.at(nPlanes -1 ), e, candidate.weights.at(nPlanes -1));
  for(int ii = nPlanes -2; ii >= 0; ii-- ){
    m_fitter.predictInfo( planes.at( ii + 1 ), planes.at(ii), e );
    m_fitter.backward.at(ii) = e;
    m_fitter.updateInfoDaf( planes.at(ii), e, candidate.weights.at(ii));
    m_fitter.addScatteringInfo( planes.at(ii), e);
  }

  m_fitter.smoothInfo();
  return(ndof);
}


template <typename T,size_t N>
void TrackerSystem<T, N>::truthTracker(){
  //A track finmder that assumes the 0th measurement should be in the fit if 
  // the plane is not excluded, it has measurements, it is in the goodRegion.
  TrackCandidate<T,N> candidate(planes.size());

  for(size_t ii = 0 ; ii < planes.size() ; ii++){
    candidate.weights.at(ii).resize( planes.at(ii).meas.size());
    candidate.weights.at(ii).setZero();
    if( (not planes.at(ii).isExcluded() ) and 
	(planes.at(ii).meas.size() > 0) and
	(planes.at(ii).meas.at(0).goodRegion())){
      candidate.weights.at(ii)(0) = 1.0;
      candidate.indexes.at(ii) = 0;
    }
    else{
      candidate.indexes.at(ii) = -1;
    }
  }
  tracks.push_back(candidate);
  m_nTracks++;
}

template <typename T,size_t N>
void TrackerSystem<T, N>::weightToIndex(TrackCandidate<T, N>& cnd){
  //Convert weights to indexes. A DAF/cluster tracker track can now be used for a non weighter information filter
  for(size_t plane = 0 ; plane < planes.size(); plane++){
    int index(-1);
    T maxWeight(0.5f);
    if( not planes.at(plane).isExcluded() ){
      for(size_t m = 0; m < planes.at(plane).meas.size(); m++){
	if( cnd.weights.at(plane)(m) < maxWeight){ continue; }
	maxWeight = cnd.weights.at(plane)(m);
	index = m;
      }
    }
    cnd.indexes.at(plane) = index;
  }
}

template <typename T,size_t N>
void TrackerSystem<T, N>::indexToWeight(TrackCandidate<T, N>& cnd){
  //Convert indexes to weights. A CKF/KF track can now be used for a weighted information filter/DAF
  for(size_t plane = 0 ; plane < planes.size(); plane++){
    int index = cnd.indexes.at(plane);
    cnd.weights.at(plane).resize( planes.at(plane).meas.size() );
    if(cnd.weights.at(plane).size() > 0){ cnd.weights.at(plane).setZero();}
    if( index >= 0){
      cnd.weights.at(plane)(index) = 1.0;
    }
  }
}

template <typename T,size_t N>
void TrackerSystem<T, N>::fitPlanesKF(TrackCandidate<T, N>& candidate){
  //Fit planes with a standard Kalman filter. Unbiased.
  TrackEstimate<T,N> e;

  e.makeSeed(false);
  //Unbiased fit for track parameters
  //Forward fitter
  for(size_t ii = 0; ii <  planes.size() ; ii++ ){
    if(not e.isSeed() ){
      m_fitter.predict( planes.at( ii - 1), planes.at(ii), e );
    }
    m_fitter.forward.at(ii) = e;
    m_fitter.kfUpdate( planes.at(ii), candidate.indexes.at(ii), e ); 
  }
  
  //Backward fitter, never bias
  e.makeSeed(true);
  for(int ii = planes.size() -1; ii >= 0; ii-- ){
    if(not e.isSeed() ){
      m_fitter.predict( planes.at( ii + 1), planes.at(ii), e );
    }
    m_fitter.backward.at(ii) = e;
    m_fitter.kfUpdate( planes.at(ii), candidate.indexes.at(ii), e );
  }
  
  getChi2Kf(candidate);
  m_fitter.smooth();
  for(int ii = 0; ii < (int) planes.size() ; ii++ ){
    candidate.estimates.at(ii) = m_fitter.smoothed.at(ii);
  }
  indexToWeight( candidate );
}

//Combinatorial KF
template <typename T,size_t N>
void TrackerSystem<T, N>::combinatorialKF(){
  // Combinatorial Kalman filter track finder.
  vector<int> indexes(planes.size(), -1);
  TrackEstimate<T,N> e;

  //Check for tracks missing a hits in first planes plane 0
  for(size_t ii = 0; ii < m_skipMax + 1; ii++){
    if( ii > 0){ indexes.at(ii -1 ) = -1;}
    for(size_t hit = 0; hit < planes.at(ii).meas.size(); hit++){
      if( ii > 0){
	bool doContinue(false);
	//Look for accepted track including hit
	for(size_t track = 0; track < getNtracks() ; track++){
	  if(tracks.at(track).indexes.at(ii) == static_cast<int>(hit)){
	    doContinue = true; break;
	  }
	}
	if(doContinue){ continue; } //Skip if measurement is included in another track
      }
      e.makeSeedInfo();
      indexes.at(ii) = hit;
      m_fitter.updateInfo(planes.at(ii), hit, e);
      fitPermutation(ii + 1, e, ii, indexes, 1, 0.0f);
    }
  }
}

template <typename T,size_t N>
void TrackerSystem<T, N>::finalizeCKFTrack(TrackEstimate<T, N>& est, vector<int>& indexes, int nMeas, T chi2){
  //Check the fitted tracl for chi2/ndof
  fastInvert(est.cov);
  est.params = est.cov * est.params;
  if(( fabs( est.getXdz() - getNominalXdz()) > getXdzMaxDeviance() ) or
     ( fabs( est.getYdz() - getNominalYdz()) > getYdzMaxDeviance())){
    return;
  }
  // Either reject the track, or save it
  TrackCandidate<T,N> candidate(planes.size());

  candidate.ndof = nMeas * 2 - 4;
  candidate.chi2 = chi2;
  if(candidate.chi2/candidate.ndof > getChi2OverNdofCut()) { return;}
  
  //Copy indexes, assign weights
  for(int plane = 0; plane < (int) planes.size(); plane++){
    candidate.indexes.at(plane) = indexes.at(plane);
  }
  indexToWeight( candidate );
  tracks.push_back(candidate);
  m_nTracks++;
}

template <typename T,size_t N>
void TrackerSystem<T, N>::fitPermutation(int plane, TrackEstimate<T, N> &est, size_t nSkipped, vector<int> &indexes, int nMeas, T chi2){
  //Check a branch of the track tree. Either kill it or, let it live.
  if( getNtracks() >= m_maxCandidates){
    cout << "Reached maximum number of track candidates, " << m_maxCandidates << endl;
    return;
  }
  //Last plane, save and quit
  if(plane == (int) planes.size()){
    finalizeCKFTrack(est, indexes, nMeas, chi2);
    return;
  }
  //Propagate
  if(nMeas > 1) { m_fitter.addScatteringInfo( planes.at(plane - 1), est);}
  m_fitter.predictInfo(planes.at( plane - 1), planes.at(plane), est);

  //Excluded plane, propagate without looking for measurement
  if( planes.at(plane).isExcluded()){
    m_fitter.forward.at(plane) = est;
    indexes.at(plane) = -1;
    fitPermutation(plane + 1, est, nSkipped, indexes, nMeas, chi2);
    return;
  }
  //Prepare for branch generation
  size_t tmpNtracks = getNtracks();
  Eigen::Matrix<T,2,1> resv, errv;
  Eigen::Matrix<T,4,1> state;
  double chi2m = 0;
  double oldX(0.0), oldY(0.0), oldZ(0.0);
  //Get prediction explicitly if needed
  if(nMeas > 1){
    Eigen::Matrix<T, N, N> tmp4x4 = est.cov;
    fastInvert(tmp4x4);
    state = tmp4x4 * est.params;
    errv = planes.at(plane).getSigmas().array().square() + tmp4x4.diagonal().head(2).array();
  }
  //If only one measurement has been read in. prepare for checking angles
  if(nMeas == 1){
    for(int ii = 0; ii < plane; ii++){
      int index = indexes.at(ii);
      if( index >= 0){
	oldX = planes.at(ii).meas.at(index).getX();
	oldY = planes.at(ii).meas.at(index).getY();
	oldZ = planes.at(ii).getZpos();
	break;
      }
    }
  }

  for(int hit = 0; hit < (int)planes.at(plane).meas.size(); hit++){
    Measurement<T>& mm = planes.at(plane).meas.at(hit);
    bool filterMeas = false;
    if( nMeas > 1) { 
      //If more than 1 measurements, get chi2
      resv = (state.head(2) - mm.getM()).array().square();
      chi2m = (resv.array() / errv.array()).sum();
      
      if (chi2m <  getCKFChi2Cut() ) { 
	filterMeas = true;
      }
    } else if(nMeas == 1){ 
      //Check angle of second plane
      double newZ = planes.at(plane).getZpos();
      if ( (fabs((mm.getX() - oldX)/(newZ - oldZ) - getNominalXdz()) < getXdzMaxDeviance()) and
	   (fabs((mm.getY() - oldY)/(newZ - oldZ) - getNominalYdz()) < getYdzMaxDeviance())){ filterMeas = true;}
    }
    //Did the measurement pass cuts? If so propagate branch
    if ( filterMeas ){ 
      TrackEstimate<T,N> clone(est);

      m_fitter.updateInfo(planes.at(plane), hit, clone);
      indexes.at(plane)= hit;
      fitPermutation(plane + 1, clone, nSkipped, indexes, nMeas + 1, chi2 + chi2m);
    }
  }
  //Skip plane if we are allowed to skip more measurements, and including a measurement did not lead to 
  // a new track candidate being accepted.
  if( tmpNtracks == getNtracks() and nSkipped < m_skipMax){
    indexes.at(plane) = -1;
    fitPermutation(plane + 1, est, nSkipped + 1, indexes, nMeas, chi2);
  }
}
