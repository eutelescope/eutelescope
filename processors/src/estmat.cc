#include "estmat.h"
#include <gsl/gsl_multimin.h>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <TH2D.h>

//#include <thread>         // std::this_thread::sleep_for
//#include <chrono>         // std::chrono::seconds
 

inline double getScatterSigma(double eBeam, double radLength){
  radLength = fabs(radLength);
  double scatterTheta = 0.0136f/ eBeam * sqrt( radLength ) *  (1.0f + 0.038f * std::log(radLength) );
  return(scatterTheta);
}

void EstMat::addTrack( std::vector<Measurement<FITTERTYPE> > track){
  //Add a track to memory
  tracks.push_back(track);
}

void EstMat::readTrack(int track, TrackerSystem<FITTERTYPE, 4>& system){
  //Read a track into the tracker system into memory
  for(size_t meas = 0; meas < tracks.at(track).size(); meas++){
    Measurement<FITTERTYPE>& m1 = tracks.at(track).at(meas);
    for(size_t ii = 0; ii < system.planes.size(); ii++){
      if( (int) m1.getIden() == (int) system.planes.at(ii).getSensorID()){
	double x = m1.getX() * ( 1.0 + xScale.at(ii)) + m1.getY() * zRot.at(ii);
	double y = m1.getY() * ( 1.0 + yScale.at(ii)) - m1.getX() * zRot.at(ii);
	x += xShift.at(ii);
	y += yShift.at(ii); 
	
	double z = m1.getZ();
	system.addMeasurement(ii, x, y, z, true, m1.getIden());
	break;
      }
    }
  }
}

void EstMat::readTracksToArray(float** measX, float** measY, int nTracks, int nPlanes){
  if(static_cast<size_t>(nTracks) > tracks.size()){
    throw std::runtime_error("Trying to read too many tracks!");
  }
  for(int tr = 0; tr < nTracks; tr++){
    if(tracks.at(tr).size() != 9 or nPlanes != 9){
      cout << "nPlanes = " << nPlanes << endl;
      throw std::runtime_error("SDR2CL currently needs exactly nine measurements in all the tracks.");
    }
    for(int pl = 0; pl < nPlanes; pl++){
      measX[pl][tr] = tracks.at(tr).at(pl).getX();
      measY[pl][tr] = tracks.at(tr).at(pl).getY();
    }
  }
}

void EstMat::readTracksToDoubleArray(float** measX, int nTracks, int nPlanes){
  if(static_cast<size_t>(nTracks) > tracks.size()){
    throw std::runtime_error("Trying to read too many tracks!");
  }
  for(int tr = 0; tr < nTracks; tr++){
    if(tracks.at(tr).size() != 9 or nPlanes != 9){
      cout << "nPlanes = " << nPlanes << endl;
      throw std::runtime_error("SDR2CL currently needs exactly nine measurements in all the tracks.");
    }
    for(int pl = 0; pl < nPlanes; pl++){
      measX[pl][(2 * tr)]     = tracks.at(tr).at(pl).getX();
      measX[pl][(2 * tr) + 1] = tracks.at(tr).at(pl).getY();
    }
  }
}

void EstMat::getExplicitEstimate(TrackEstimate<FITTERTYPE, 4>& estim){
  //Get the explicit estimates
  fastInvert(estim.cov);
  //Sparse multiply
  FITTERTYPE x(estim.params(0)), y(estim.params(1)), dx(estim.params(2)), dy(estim.params(3));
  Eigen::Matrix<FITTERTYPE, 4, 1> offdiag, reverse(dx, dy, x, y);
  offdiag(0) = offdiag(2) = estim.cov(0,2);
  offdiag(1) = offdiag(3) = estim.cov(1,3);
  estim.params = estim.params.array() * estim.cov.diagonal().array() + offdiag.array() * reverse.array();

  // estim.params(0) = estim.cov(0,0) * x + estim.cov(0,2) * dx;
  // estim.params(1) = estim.cov(1,1) * y + estim.cov(1,3) * dy;
  // estim.params(2) = estim.cov(2,2) * dx+ estim.cov(2,0) * x;
  // estim.params(3) = estim.cov(3,3) * dy+ estim.cov(3,1) * y;
} 

void FakeChi2::init(){
  Minimizer::init();
  firstRun = true;
}

void FakeChi2::calibrate(TrackerSystem<FITTERTYPE,4>& system){
  cout << "Calculating residual errors" << endl;
  resFWErrorX.resize(system.planes.size());
  resFWErrorY.resize(system.planes.size());
  resBWErrorX.resize(system.planes.size());
  resBWErrorY.resize(system.planes.size());
  system.clear();
  mat.readTrack(0,system);
  system.index0tracker();
  TrackCandidate<FITTERTYPE,4> candidate = system.tracks.at(0);
  system.fitInfoFWUnBiased(candidate);
  //Matrix<FITTERTYPE, 2, 1> errors = 

  Eigen::Matrix<FITTERTYPE,2,1> errv;
  
  for(size_t ii = 2; ii < system.planes.size(); ii++){
    if(system.planes.at(ii).isExcluded()) { continue;}
    TrackEstimate<FITTERTYPE,4>& estim = system.m_fitter.forward.at(ii);
    mat.getExplicitEstimate(estim);
    errv = system.getUnBiasedResidualErrors( system.planes.at(ii), estim);
    resFWErrorX.at(ii) = errv[0];
    resFWErrorY.at(ii) = errv[1];
  }

  system.fitInfoBWUnBiased(candidate);
  for(size_t ii = 0; ii < system.planes.size() -2; ii++){
    if(system.planes.at(ii).isExcluded()) { continue;}
    TrackEstimate<FITTERTYPE,4>& estim = system.m_fitter.backward.at(ii);
    mat.getExplicitEstimate(estim);
    errv = system.getUnBiasedResidualErrors( system.planes.at(ii), estim);
    resBWErrorX.at(ii) = errv[0];
    resBWErrorY.at(ii) = errv[1];
  }
  firstRun = false;
}

void FakeChi2::operator() (size_t offset, size_t stride){
  //Get the global chi2 of the track sample
  TrackerSystem<FITTERTYPE, 4>& system = systems.at(offset);
  
  //Track candidate is the same for all tracks
  system.index0tracker();
  TrackCandidate<FITTERTYPE,4> candidate = system.tracks.at(0);
  
  {
#ifdef DOTHREAD
    boost::mutex::scoped_lock lock(resultGurad);
#endif
    if(firstRun){ calibrate(system); }
  }
  Eigen::Matrix<FITTERTYPE, 2, 1> resv;
  
  FITTERTYPE chi2 = 0;
  for(int track = offset; track < mat.itMax; track += stride){
    //prepare system for new track: clear system from prev go around, read track from memory, run track finder
    system.clear();
    mat.readTrack(track,system);
    system.fitInfoFWUnBiased(candidate);
    //Get explicit estimates
    for(size_t pl = 2; pl < system.planes.size(); pl++){
      TrackEstimate<FITTERTYPE,4>& fw = system.m_fitter.forward[pl];
      mat.getExplicitEstimate(fw);
      Measurement<FITTERTYPE>& meas = system.planes[pl].meas.at(0);
      resv = system.getResiduals(meas, fw);
      chi2 += resv(0) * resv(0)/resFWErrorX[pl] + resv(1) * resv(1)/resFWErrorY[pl]; 
    }

    //BW
    system.fitInfoBWUnBiased(candidate);
    //Get explicit estimates
    for(size_t pl = 0; pl < system.planes.size() - 2; pl++){
      TrackEstimate<FITTERTYPE,4>& bw = system.m_fitter.backward[pl];
      mat.getExplicitEstimate(bw);
      Measurement<FITTERTYPE>& meas = system.planes[pl].meas.at(0);
      resv = system.getResiduals(meas, bw);
      chi2 += resv(0) * resv(0)/resBWErrorX[pl] + resv(1) * resv(1)/resBWErrorY[pl]; 
    }
  }

  {
#ifdef DOTHREAD
    boost::mutex::scoped_lock lock(resultGurad);
#endif
    result += chi2;
  }
}

void FakeAbsDev::operator() (size_t offset, size_t stride){
  //Get the global chi2 of the track sample
  TrackerSystem<FITTERTYPE,4>& system = systems.at(offset);
  
  //Track candidate is the same for all tracks
  system.index0tracker();
  TrackCandidate<FITTERTYPE,4> candidate = system.tracks.at(0);
  
  {
#ifdef DOTHREAD
    boost::mutex::scoped_lock lock(resultGurad);
#endif
    if(firstRun){ calibrate(system); }
  }
  Eigen::Matrix<FITTERTYPE, 2, 1> resv;

  FITTERTYPE chi2 = 0;
  for(int track = offset; track < mat.itMax; track += stride){
    //prepare system for new track: clear system from prev go around, read track from memory, run track finder
    system.clear();
    mat.readTrack(track,system);
    system.fitInfoFWUnBiased(candidate);
    //Get explicit estimates
    for(size_t pl = 2; pl < system.planes.size(); pl++){
      TrackEstimate<FITTERTYPE,4>& fw = system.m_fitter.forward[pl];
      mat.getExplicitEstimate(fw);
      Measurement<FITTERTYPE>& meas = system.planes[pl].meas.at(0);
      resv = system.getResiduals(meas, fw);
      chi2 += fabs(resv[0]/sqrt(resFWErrorX[pl])) + fabs(resv[1]/sqrt(resFWErrorY[pl])); 
    }

    //BW
    system.fitInfoBWUnBiased(candidate);
    //Get explicit estimates
    for(size_t pl = 0; pl < system.planes.size() - 2; pl++){
      TrackEstimate<FITTERTYPE,4>& bw = system.m_fitter.backward[pl];
      mat.getExplicitEstimate(bw);
      Measurement<FITTERTYPE>& meas = system.planes[pl].meas.at(0);
      resv = system.getResiduals(meas, bw);
      chi2 += fabs(resv[0]/sqrt(resBWErrorX[pl])) + fabs(resv[1]/sqrt(resBWErrorY[pl]));
    }
  }

  {
#ifdef DOTHREAD
    boost::mutex::scoped_lock lock(resultGurad);
#endif
    result += chi2;
  }
}

void Chi2::operator() (size_t offset, size_t stride){
  //Get the global chi2 of the track sample
  TrackerSystem<FITTERTYPE,4>& system = systems.at(offset);
  
  //Track candidate is the same for all tracks
  system.index0tracker();
  TrackCandidate<FITTERTYPE,4> candidate = system.tracks.at(0);
  
  double varchi2(0.0);
  for(int track = offset; track < mat.itMax; track+= stride){
    system.clear();
    mat.readTrack(track, system);
    system.fitInfoFWBiased(candidate);
    system.getChi2BiasedInfo(candidate);
    varchi2 += candidate.chi2;
  }
  
  {
#ifdef DOTHREAD
    boost::mutex::scoped_lock lock(resultGurad);
#endif
    result += varchi2;
  }
}
  
void SDR::operator() (size_t offset, size_t stride){
  //Get the mean^2 + (1 - variance) of the standardized residuals of chi2 increments and or pull distributions
  TrackerSystem<FITTERTYPE,4>& system = systems.at(offset);
  std::vector< double > sqrPullXFW(system.planes.size() - 2, 0.0 );
  std::vector< double > sqrPullXBW(system.planes.size() - 2, 0.0 );
  std::vector< double > sqrPullYFW(system.planes.size() - 2, 0.0 );
  std::vector< double > sqrPullYBW(system.planes.size() - 2, 0.0 );
  std::vector< std::vector<double> > sqrParams(system.planes.size() - 3, std::vector<double>(4 ,0.0) );
  int nTracks = 0;
    
  //Track candidate is the same for all tracks
  system.index0tracker();
  TrackCandidate<FITTERTYPE,4> candidate = system.tracks.at(0);

  for(int track = offset; track < mat.itMax; track += stride){
    //prepare system for new track: clear system from prev go around, read track from memory, run track finder
    system.clear();
    mat.readTrack(track, system);
    
    //Only one track! Skip track finder
    //Run FW fitter, get p-values
    system.fitInfoFWBiased(candidate);
    system.fitInfoBWUnBiased(candidate);
    //Get explicit estimates
    TrackEstimate<FITTERTYPE,4>& fw1 = system.m_fitter.forward.at(1);
    mat.getExplicitEstimate(fw1);
    
    for(size_t pl = 0; pl < system.planes.size() -2; pl++){
      //We're using an information filter, we need explicit state and covariance
      TrackEstimate<FITTERTYPE,4>& fw = system.m_fitter.forward.at(pl + 2);
      TrackEstimate<FITTERTYPE,4>& bw = system.m_fitter.backward.at(pl);
      mat.getExplicitEstimate(fw);
      mat.getExplicitEstimate(bw);
    }
    //Chi2 increments FW
    for(size_t pl = 2; pl < system.planes.size(); pl++){
      TrackEstimate<FITTERTYPE,4>& result = system.m_fitter.forward.at(pl);
      Measurement<FITTERTYPE>& meas = system.planes.at(pl).meas.at(0);
      Eigen::Matrix<FITTERTYPE, 2, 1> resids = system.getResiduals(meas, result);
      Eigen::Matrix<FITTERTYPE, 2, 1> variance = system.getBiasedResidualErrors(system.planes.at(pl), result);
      //Squared pulls to calculate pull variance
      Eigen::Matrix<FITTERTYPE, 2, 1> pull2 = resids.array().square() / variance.array();
      sqrPullXFW.at(pl - 2) += pull2(0); 
      sqrPullYFW.at(pl - 2) += pull2(1); 
    }
    //Chi2 increments BW
    for(size_t pl = 0; pl < system.planes.size() - 2; pl++){
      TrackEstimate<FITTERTYPE,4>& result = system.m_fitter.backward.at(pl);
      Measurement<FITTERTYPE>& meas = system.planes.at(pl).meas.at(0);
      Eigen::Matrix<FITTERTYPE, 2, 1> resids = system.getResiduals(meas, result);
      Eigen::Matrix<FITTERTYPE, 2, 1> variance = system.getUnBiasedResidualErrors(system.planes.at(pl), result);
      //Squared pulls to calculate pull variance
      Eigen::Matrix<FITTERTYPE, 2, 1> pull2 = resids.array().square() / variance.array();
      sqrPullXBW.at(pl) += pull2(0);
      sqrPullYBW.at(pl) += pull2(1);
    }
    //Difference in parameters
    if(not cholDec){
      for(size_t pl = 1; pl < system.planes.size() -2; pl++){
	TrackEstimate<FITTERTYPE,4>& fw = system.m_fitter.forward.at(pl);
	TrackEstimate<FITTERTYPE,4>& bw = system.m_fitter.backward.at(pl);
	for(size_t param = 0; param < 4; param++){
	  double var = fw.cov(param,param) + bw.cov(param,param);
	  double res = fw.params(param) - bw.params(param);
	  sqrParams.at(pl -1).at(param) += res * res / var;
	}
      }
    } else {
      //Attempt at getting pull values from cholesky decomposed covariance.
      for(size_t pl = 1; pl < system.planes.size() -2; pl++){
	TrackEstimate<FITTERTYPE,4>& fw = system.m_fitter.forward.at(pl);
	TrackEstimate<FITTERTYPE,4>& bw = system.m_fitter.backward.at(pl);
	Eigen::Matrix<double, 4, 1> tmpDiff = (fw.params - bw.params).cast<double>();
	Eigen::Matrix<FITTERTYPE, 4, 4> tmpCov =  fw.cov + bw.cov;
	
	//Use doubles, else it fails
	Eigen::LLT<Eigen::Matrix4d> tmpChol;
	tmpChol.compute(tmpCov.cast<double>());
	tmpChol.solve(tmpDiff);
	//tmpChol.matrixL().marked<Eigen::Lower>().solveTriangularInPlace(tmpDiff);
	// tmpChol.matrixL().trinagularView<Lower>().solveInPlace(tmpDiff);
	for(size_t param = 0; param < 4; param++){
	  sqrParams.at(pl -1).at(param) += tmpDiff(param) * tmpDiff(param);
	}
      }
    }
    nTracks++;
  }
  
  double varvar(0.0);
  if(SDR2){
    for( size_t pl = 0; pl < system.planes.size() - 2; pl++){
      double resvar = 1.0f - (sqrPullXFW.at(pl)/(nTracks - 1));
      varvar += resvar * resvar;
      resvar = 1.0f - (sqrPullYFW.at(pl)/(nTracks - 1));
      varvar += resvar * resvar;
      resvar = 1.0f - (sqrPullXBW.at(pl)/(nTracks - 1));
      varvar += resvar * resvar;
      resvar = 1.0f - (sqrPullYBW.at(pl)/(nTracks - 1));
      varvar += resvar * resvar;
    }
  }
  if(SDR1){
    for( size_t pl = 1; pl < system.planes.size() - 2; pl++){
      for(int param = 0; param < 4; param++){
	double resvar = 1.0f - (sqrParams.at(pl - 1).at(param) / (nTracks - 1));
	varvar += resvar * resvar;
      }
    }
  }
  {
#ifdef DOTHREAD
    boost::mutex::scoped_lock lock(resultGurad);
#endif
    result += varvar;
  }
}

void FwBw::operator() (size_t offset, size_t stride){
  //Get the negative log likelihood of the state difference of a forward and
  //backward running Kalman filter.
  TrackerSystem<FITTERTYPE,4>& system = systems.at(offset);
  
  //Track candidate is the same for all tracks
  system.index0tracker();
  TrackCandidate<FITTERTYPE,4> candidate = system.tracks.at(0);

  double logL(0.0);
  std::vector< double > sqrPullXFW(system.planes.size() - 2, 0.0 );
  std::vector< double > sqrPullXBW(system.planes.size() - 2, 0.0 );
  std::vector< double > sqrPullYFW(system.planes.size() - 2, 0.0 );
  std::vector< double > sqrPullYBW(system.planes.size() - 2, 0.0 );
  std::vector< std::vector<double> > sqrParams(system.planes.size() - 3, std::vector<double>(4 ,0.0) ); 
  int nTracks = 0;
    
  for(int track = offset; track < mat.itMax; track += stride){
    //prepare system for new track: clear system from prev go around, read track from memory, run track finder
    system.clear();
    mat.readTrack(track,system);
    nTracks++;
    //Translate candidate from DAF to KF
    system.fitInfoFWBiased(candidate);
    system.fitInfoBWUnBiased(candidate);
    //Get explicit estimates
    TrackEstimate<FITTERTYPE,4>& fw1 = system.m_fitter.forward.at(1);
    mat.getExplicitEstimate(fw1);
    for(size_t pl = 0; pl < system.planes.size() -2; pl++){
      TrackEstimate<FITTERTYPE,4>& fw = system.m_fitter.forward.at(pl + 2);
      TrackEstimate<FITTERTYPE,4>& bw = system.m_fitter.backward.at(pl);
      mat.getExplicitEstimate(fw);
      mat.getExplicitEstimate(bw);
    }
    Eigen::Matrix<FITTERTYPE, 4, 4> cov;
    Eigen::Matrix<FITTERTYPE, 4, 1> resids;
    //Difference in parameters
    for(size_t pl = 1; pl < system.planes.size() -2; pl++){
      TrackEstimate<FITTERTYPE,4>& fw = system.m_fitter.forward.at(pl);
      TrackEstimate<FITTERTYPE,4>& bw = system.m_fitter.backward.at(pl);
      resids = fw.params - bw.params;
      cov = fw.cov + bw.cov;
      double determinant = cov.determinant();
	
      fastInvert(cov);
      double exponent = (resids.transpose() * cov * resids)(0,0);
      logL -= log( determinant ) +  exponent;
    }
    //Chi2 increments FW
    for(size_t pl = 2; pl < system.planes.size(); pl++){
      //I'm using an information filter, need explicit state and covariance
      TrackEstimate<FITTERTYPE,4>& result = system.m_fitter.forward.at(pl);
      Measurement<FITTERTYPE>& meas = system.planes.at(pl).meas.at(0);
      //Get residuals in x and y, squared
      Eigen::Matrix<FITTERTYPE, 2, 1> resids = system.getResiduals(meas, result);
      //Get variance of residuals in x and y
      Eigen::Matrix<FITTERTYPE, 2, 1> variance = system.getBiasedResidualErrors(system.planes.at(pl), result);
      Eigen::Matrix<FITTERTYPE, 2, 1> pull2 = resids.array().square() / variance.array();
      sqrPullXFW.at(pl - 2) += pull2(0); 
      sqrPullYFW.at(pl - 2) += pull2(1); 
    }
    //Chi2 increments BW
    for(size_t pl = 0; pl < system.planes.size() - 2; pl++){
      TrackEstimate<FITTERTYPE,4>& result = system.m_fitter.backward.at(pl);
      Measurement<FITTERTYPE>& meas = system.planes.at(pl).meas.at(0);
      Eigen::Matrix<FITTERTYPE, 2, 1> resids = system.getResiduals(meas, result);
      Eigen::Matrix<FITTERTYPE, 2, 1> variance = system.getUnBiasedResidualErrors(system.planes.at(pl), result);
      Eigen::Matrix<FITTERTYPE, 2, 1> pull2 = resids.array().square() / variance.array();
      sqrPullXBW.at(pl) += pull2(0); 
      sqrPullYBW.at(pl) += pull2(1); 
    }
  }
  FITTERTYPE return2 = 0.0;
  for( size_t pl = 0; pl < system.planes.size() - 2; pl++){
    double resvar = 1.0 - sqrPullXFW.at(pl)/(nTracks - 1);
    return2 += resvar * resvar;
    resvar = 1.0 - sqrPullYFW.at(pl)/(nTracks - 1);
    return2 += resvar * resvar;
    resvar = 1.0 - sqrPullXBW.at(pl)/(nTracks - 1);
    return2 += resvar * resvar;
    resvar = 1.0 - sqrPullYBW.at(pl)/(nTracks - 1);
    return2 += resvar * resvar;
  }
  {
#ifdef DOTHREAD
    boost::mutex::scoped_lock lock(resultGurad);
#endif
    result += -1.0 * logL;
    retVal2 += return2;
  }
}

void Minimizer::init(){
  //Initialize nThread threads
#ifndef DOTHREAD
  cout << "Not using threads!" << endl;
  nThreads = 1;
#endif
  if( not inited){
    systems.assign(nThreads, mat.system);
  }
  inited = true;
}

void Minimizer::prepareThreads(){
  //Copy thicknesses and resolutions, reset resturn values
  for(size_t ii = 0; ii < mat.system.planes.size(); ii++){
    FitPlane<FITTERTYPE>& plO = mat.system.planes.at(ii); //Original plane
    for(size_t thread = 0; thread < nThreads; thread++){
      FitPlane<FITTERTYPE>& plT = systems.at(thread).planes.at(ii); //Thread plane
      plT.setScatterThetaSqr( plO.getScatterThetaSqr());
      plT.setSigmas( plO.getSigmaX(), plO.getSigmaY());
    }
  }
  result = 0;
  retVal2 = 0.0f;
}

FITTERTYPE Minimizer::operator() (void){
  //Start threads if DOTHREAD, run job in main thread if not.
  prepareThreads();
#ifdef DOTHREAD
  boost::thread_group threads;
  for(int ii = 0; ii < nThreads; ii++){
    threads.create_thread( boost::bind(&Minimizer::operator(), this, ii, nThreads));
  }
  threads.join_all();
#else
  (*this)(0, 1);
#endif
  return( result );
}

// void EstMat::simulate(int nTracks){
//   // Toy simulation of a straight track with Gaussian uncertainties and scattering
//   cout << "Simulation parameters" << endl;
//   printParams( (char*) "params[\"RadiationLengths\"]", radLengths, true, "%4.6f ");
//   printParams( (char*) "params[\"ResolutionX\"]", resX, true, "%4.6f ");
//   printParams( (char*) "params[\"ResolutionY\"]", resY, true, "%4.6f ");

//   size_t nPlanes = system.planes.size();
  
//   //tracks
//   for( int track = 0; track < nTracks; track++){
//     std::vector<Measurement<FITTERTYPE> > simTrack;
//     if( track + 1 % 100000 == 0){
//       std::cout << "Simulating track: " << track << std::endl;
//     }
//     double x(0.0f), y(0.0f), dx(0.0f), dy(0.0f);
//     gaussRand(x, y);
//     x *= 10000.0f;
//     y *= 10000.0f;

//     //Shit in positive direction
//     x += 20000.0f;
//     y += 20000.0f;
    
//     double g1, g2;
//     gaussRand(g1, g2);
//     dx += g1 * 0.0001f;
//     dy += g2 * 0.0001f;
    
//     double zPos = 0;
//     for(size_t pl = 0; pl < nPlanes; pl++){
//       double zDistance = system.planes.at(pl).getZpos() - zPos;
//       zPos = system.planes.at(pl).getZpos();
//       x += dx * zDistance;
//       y += dy * zDistance;
//       gaussRand(g1, g2);
//       dx += g1 * getScatterSigma(eBeam, radLengths.at(pl));
//       dy += g2 * getScatterSigma(eBeam, radLengths.at(pl));
//       gaussRand(g1, g2);
//       Eigen::Matrix<FITTERTYPE, 2, 1> sigmas;
//       sigmas(0) = resX.at(pl);
//       sigmas(1) = resY.at(pl);

//       simTrack.push_back( Measurement<FITTERTYPE>(x + g1 * sigmas(0), y + g2 * sigmas(1), system.planes.at(pl).getZpos(), true, pl) );
//     }
//     addTrack(simTrack);
//   }
// }

// void EstMat::initSim(int nPlanes){
//   //Initialize random number generators, tracker system, ...
//   //seed with time
//   srandom ( time(NULL) );

//   system.setClusterRadius( 100000.0f);
//   system.setNominalXdz(0.0f);
//   system.setNominalYdz(0.0f);
//   system.setChi2OverNdofCut( 100.0f);
//   system.setDAFChi2Cut( 1000.0f);

//   for(int ii = 0; ii < nPlanes; ii++){
//     double scatterTheta = getScatterSigma(eBeam, 1.0);
//     system.addPlane(ii, ii, 4.3f, 4.3f, scatterTheta * scatterTheta, false);
//   }

//   system.setMaxCandidates(1);
//   system.init();
//   xShift.resize(nPlanes);
//   yShift.resize(nPlanes);
// }

void EstMat::movePlaneZ(int planeIndex, double zPos){
  //Move the plane in the z direction
  FitPlane<FITTERTYPE>& pl = system.planes.at(planeIndex);
  pl.setZpos( zPos);
  Eigen::Matrix<FITTERTYPE, 3, 1> ref = pl.getRef0();
  ref(2) += zPos;
  pl.setRef0( ref );
}

void EstMat::setPlane(int index, double sigmaX, double sigmaY, double radLength){
  //Change the state of a plane
  double scatterTheta = getScatterSigma(eBeam, radLength);
  system.planes.at(index).setScatterThetaSqr( scatterTheta * scatterTheta);
  system.planes.at(index).setSigmas(sigmaX, sigmaY);
}

void EstMat::printAllFreeParams(){
  bool resXYp = resXYIndex.size() or resXYMulti.size();
  printParams("params[\"RadiationLengths\"]", radLengths, radLengthsIndex.size() or radLengthsMulti.size(), "%4.6f " );
  printParams("params[\"ResolutionX\"]", resX, resXIndex.size() or resXMulti.size() or resXYp, "%4.6f ");
  printParams("params[\"ResolutionY\"]", resY, resYIndex.size() or resYMulti.size() or resXYp, "%4.6f ");
  printParams("params[\"XShift\"]", xShift, xShiftIndex.size(), "%4.2f ");
  printParams("params[\"YShift\"]", yShift, yShiftIndex.size(), "%4.2f ");
  printParams("params[\"XScale\"]", xScale, xScaleIndex.size(), "%4.3f ");
  printParams("params[\"YScale\"]", yScale, yScaleIndex.size(), "%4.3f ");
  printParams("params[\"ZRot\"]", zRot, zRotIndex.size(), "%4.3f ");
  printParams("params[\"ZPosition\"]", zPos, zPosIndex.size(), "%4.2f ");
}

void printHisto(TH1D* histo){
  cout << "(list" 
       <<  " :min " << histo->GetXaxis()->GetXmin()
       << " :bin-size " << histo->GetBinWidth( 0 ) << endl
       << " :data (list" << endl;
  for(int ii = 1; ii <= histo->GetNbinsX(); ii++){
    cout << " " << histo->GetBinContent( ii );
  }
  cout << endl << "))" << endl;
}

void EstMat::printParams(std::string name, std::vector<FITTERTYPE>& params, bool plot, const char* valString){
  //Print estimation parameters to screen
  if(not plot){ return;}
  std::cout << name << " = \" ";
  for( size_t pos = 0; pos < params.size(); pos++){
    printf( valString, params.at(pos));
  }
  std::cout << " \"" << std::endl;
}

void EstMat::plot(char* fname){
  //Make, fill and save some plots. The full track sample is fitted for this
  for(size_t plane = 0; plane < system.planes.size(); plane ++){
    double scatterSigma = getScatterSigma(eBeam, radLengths.at(plane));
    system.planes.at(plane).setScatterThetaSqr( scatterSigma * scatterSigma);
  }
  
  TFile* tfile = new TFile(fname, "RECREATE");
  std::vector<TH1D*> resX, resY, pullX, pullY;
  std::vector<TH1D*> xFB, yFB, dxFB, dyFB, pValFW;
  //TH1D* chi2 = new TH1D("chi2","chi2", 150, 0, 150);
  TH1D* chi2 = new TH1D("chi2","chi2", 100, 0, 100);
  TH1D* ndof = new TH1D("ndof", "ndof", 15, 0, 15);
  TH1D* pValue = new TH1D("pValue", "pValue", 100,0,1);
  TH1D* chi2OverNeod = new TH1D("chi2ndof","chi2 over ndof", 100, 0, 20);
  TH2D* corr34X = new TH2D("correlations34x", "correlations34x", 100, -4, 4, 100, -4, 4);
  TH2D* corr34Y = new TH2D("correlations34y", "correlations34y", 100, -4, 4, 100, -4, 4);
  TH2D* corr12X = new TH2D("correlations12x", "correlations12x", 100, -4, 4, 100, -4, 4);
  TH2D* corr12Y = new TH2D("correlations12y", "correlations12y", 100, -4, 4, 100, -4, 4);
  for(int ii = 0; ii < (int) system.planes.size(); ii++){
    char name[200];
    int sensorID = system.planes.at(ii).getSensorID();
    sprintf(name, "resX %i", sensorID); 
    resX.push_back( new TH1D(name,name, 800, -400, 400));
    sprintf(name, "resY %i", sensorID); 
    resY.push_back( new TH1D(name,name, 800, -400, 400));
    sprintf(name, "pullX unbiased  %i", sensorID); 
    //pullX.push_back( new TH1D(name,name, 100, -10, 10));
    pullX.push_back( new TH1D(name,name, 600, -30, 30));
    sprintf(name, "pullY biased  %i", sensorID); 
    //pullY.push_back( new TH1D(name,name, 100, -10, 10));
    pullY.push_back( new TH1D(name,name, 600, -30, 30));
  }
  
  cout << "Inited plots" << endl;
  
  //Loop over all tracks
  for(size_t track = 0; track < tracks.size(); track++){
    system.clear();
    readTrack(track, system);
    system.clusterTracker();
    for(size_t ii = 0; ii < system.getNtracks(); ii++ ){
      TrackCandidate<FITTERTYPE,4>& candidate = system.tracks.at(ii);
      system.weightToIndex(candidate);
      system.fitPlanesInfoUnBiased(candidate);
      //get p-val
      pValue->Fill( TMath::Gamma( candidate.ndof / 2, candidate.chi2 / 2.0) );
      chi2->Fill( candidate.chi2 );
      ndof->Fill( candidate.ndof );
      chi2OverNeod->Fill( candidate.chi2 / candidate.ndof );

      for(int ii = 0; ii < (int) system.planes.size() ; ii++ ){
	candidate.estimates.at(ii) = system.m_fitter.smoothed.at(ii);
      }
      for(size_t pl = 0; pl < system.planes.size(); pl++){
      	if( system.planes.at(pl).meas.size() < 1){ continue; }
      	Measurement<FITTERTYPE>& meas = system.planes.at(pl).meas.at(0);
	TrackEstimate<FITTERTYPE,4>& est = candidate.estimates.at(pl);
	Eigen::Matrix<FITTERTYPE, 2, 1> resids = system.getResiduals(meas, est);
	Eigen::Matrix<FITTERTYPE, 2, 1> variance = system.getUnBiasedResidualErrors(system.planes.at(pl), est);
      	resX.at(pl)->Fill(resids(0));
      	resY.at(pl)->Fill(resids(1));
      	pullX.at(pl)->Fill(resids(0)/sqrt(variance(0)));
      	pullY.at(pl)->Fill(resids(1)/sqrt(variance(1)));
      }
      //Plot correlations
      if( (system.planes.at(3).meas.size() > 0) and
	  (system.planes.at(4).meas.size() > 0)){
	system.planes.at(3).exclude();
	system.planes.at(4).exclude();
	system.fitPlanesInfoUnBiased(candidate);

	TrackEstimate<FITTERTYPE,4>& est3 = candidate.estimates.at(3);
      	Measurement<FITTERTYPE>& meas3 = system.planes.at(3).meas.at(0);
	Eigen::Matrix<FITTERTYPE, 2, 1> resids3 = system.getResiduals(meas3, est3);
	Eigen::Matrix<FITTERTYPE, 2, 1> variance3 = system.getUnBiasedResidualErrors(system.planes.at(3), est3);
	
	TrackEstimate<FITTERTYPE,4>& est4 = candidate.estimates.at(4);
      	Measurement<FITTERTYPE>& meas4 = system.planes.at(4).meas.at(0);
	Eigen::Matrix<FITTERTYPE, 2, 1> resids4 = system.getResiduals(meas4, est4);
	Eigen::Matrix<FITTERTYPE, 2, 1> variance4 = system.getUnBiasedResidualErrors(system.planes.at(4), est4);
	
	corr34X->Fill( resids3(0)/sqrt(variance3(0)), resids4(0)/sqrt(variance4(0)));
	corr34Y->Fill( resids3(1)/sqrt(variance3(1)), resids4(1)/sqrt(variance4(1)));
	system.planes.at(3).include();
	system.planes.at(4).include();
      }
      
      if( (system.planes.at(1).meas.size() > 0) and
	  (system.planes.at(2).meas.size() > 0)){
	
	system.planes.at(1).exclude();
	system.planes.at(2).exclude();
	system.fitPlanesInfoUnBiased(candidate);
	TrackEstimate<FITTERTYPE,4>& est6 = candidate.estimates.at(6);
      	Measurement<FITTERTYPE>& meas6 = system.planes.at(6).meas.at(0);
	Eigen::Matrix<FITTERTYPE, 2, 1> resids6 = system.getResiduals(meas6, est6);
	Eigen::Matrix<FITTERTYPE, 2, 1> variance6 = system.getUnBiasedResidualErrors(system.planes.at(6), est6);
	
	TrackEstimate<FITTERTYPE,4>& est2 = candidate.estimates.at(2);
      	Measurement<FITTERTYPE>& meas2 = system.planes.at(2).meas.at(0);
	Eigen::Matrix<FITTERTYPE, 2, 1> resids2 = system.getResiduals(meas2, est2);
	Eigen::Matrix<FITTERTYPE, 2, 1> variance2 = system.getUnBiasedResidualErrors(system.planes.at(2), est2);
	
	// corr12X->Fill( resids2(0), resids6(0));
	// corr12Y->Fill( resids2(1), resids6(1));
	corr12X->Fill( resids2(0)/sqrt(variance2(0)), resids6(0)/sqrt(variance6(0)));
	corr12Y->Fill( resids2(1)/sqrt(variance2(1)), resids6(1)/sqrt(variance6(1)));
	system.planes.at(1).include();
	system.planes.at(2).include();
      }
    }
  }
  cout << "Saving " << endl;
  tfile->cd();
  chi2->Write();
  pValue->Write();
  ndof->Write();
  chi2OverNeod->Write();
  corr34X->Write();
  corr34Y->Write();
  corr12X->Write();
  corr12Y->Write();
  for(size_t ii = 0; ii < system.planes.size(); ii++){
    resX.at(ii)->Write();
    resY.at(ii)->Write();
    pullX.at(ii)->Write();
    pullY.at(ii)->Write();
  }

  cout << "Write" << endl;
  tfile->Write();
  cout << "Done plotting to " << fname << endl;
}

//simplex stuff
size_t EstMat::getNSimplexParams(){
  //The number of parameters to be estimated
  size_t params(0);
  params += resXIndex.size();
  params += resYIndex.size();
  params += radLengthsIndex.size();
  params += resXMulti.size();
  params += resYMulti.size();
  params += resXYMulti.size();
  params += xShiftIndex.size();
  params += yShiftIndex.size();
  params += xScaleIndex.size();
  params += yScaleIndex.size();
  params += zRotIndex.size();
  params += zPosIndex.size();
  return( params );
}

void EstMat::vectToEst(gsl_vector* v, std::vector<int>& indexVector, vector<FITTERTYPE>& dataVector, size_t& param){
  for(size_t ii = 0; ii<indexVector.size();ii++){
    gsl_vector_set(v, param++, dataVector.at( indexVector.at(ii)));
  }
}
void EstMat::multiVectToEst(gsl_vector* v, std::vector<int>& indexVector, vector<FITTERTYPE>& dataVector, size_t& param){
  if(indexVector.size() > 0){
    gsl_vector_set(v,param++, dataVector.at(0));
  }
}

gsl_vector* EstMat::systemToEst(){
  //Read the state of the tracker system into a gsl vector for the simplex search
  int nParams = getNSimplexParams();
  gsl_vector* v = gsl_vector_alloc(nParams);
  size_t param = 0;
  //Material and resolutions
  vectToEst(v, resXIndex, resX, param);
  vectToEst(v, resYIndex, resY, param);
  vectToEst(v, radLengthsIndex, radLengths, param);
  multiVectToEst(v, resXMulti, resX, param);
  multiVectToEst(v, resYMulti, resY, param);
  multiVectToEst(v, resXYMulti, resX, param);

  //Alignment
  vectToEst(v, xShiftIndex, xShift, param);
  vectToEst(v, yShiftIndex, yShift, param);
  vectToEst(v, xScaleIndex, xScale, param);
  vectToEst(v, yScaleIndex, yScale, param);
  vectToEst(v, zRotIndex, zRot, param);
  vectToEst(v, zPosIndex, zPos, param);
  return(v);
}

void EstMat::estToSystem( const gsl_vector* params, TrackerSystem<FITTERTYPE,4>& system){
  //Set the system state from the gsl estimation vector
  double param = 0;
  vector<int>::iterator it;
  for(it = resXIndex.begin(); it != resXIndex.end(); it++){
    resX.at((*it)) = gsl_vector_get(params, param++);
    system.planes.at((*it)).setSigmaX( resX.at((*it)));
  }
  for(it = resYIndex.begin(); it != resYIndex.end(); it++){
    resY.at((*it)) = gsl_vector_get(params, param++);
    system.planes.at((*it)).setSigmaY( resY.at((*it)));
  }
  for(it = radLengthsIndex.begin(); it != radLengthsIndex.end(); it++){
    radLengths.at((*it)) = gsl_vector_get(params, param++);
    double scatterTheta = getScatterSigma(eBeam, radLengths.at((*it)));
    system.planes.at((*it)).setScatterThetaSqr( scatterTheta * scatterTheta);
  }
  if( resXMulti.size() > 0){
   for(size_t ii = 0; ii < resXMulti.size(); ii++){
      int index = resXMulti.at(ii);
      resX.at(index) = gsl_vector_get(params, param++);
      system.planes.at(index).setSigmaX(resX.at(index));
    }
  }
  if( resYMulti.size() > 0){
    for(size_t ii = 0; ii < resYMulti.size(); ii++){
      int index = resYMulti.at(ii);
      resY.at(index) = gsl_vector_get(params, param++);
      system.planes.at(index).setSigmaY(resX.at(index));
    }
  }
  if( resXYMulti.size() > 0){
    for(size_t ii = 0; ii < resXYMulti.size(); ii++){
      int index = resXYMulti.at(ii);
      resX.at(index) = resY.at(index) = gsl_vector_get(params, param++);
      system.planes.at(index).setSigmas(resX.at(index), resY.at(index));
    }
  }

  //Alignment
  for(it = xShiftIndex.begin(); it != xShiftIndex.end(); it++){
    xShift.at((*it)) = (gsl_vector_get(params, param++));
  }
  for(it = yShiftIndex.begin(); it != yShiftIndex.end(); it++){
    yShift.at((*it)) = (gsl_vector_get(params, param++));
  }
  for(it = xScaleIndex.begin(); it != xScaleIndex.end(); it++){
    xScale.at((*it)) = (gsl_vector_get(params, param++));
  }
  for(it = yScaleIndex.begin(); it != yScaleIndex.end(); it++){
    yScale.at((*it)) = (gsl_vector_get(params, param++));
  }
  for(it = zRotIndex.begin(); it != zRotIndex.end(); it++){
    zRot.at((*it)) = (gsl_vector_get(params, param++));
  }
  for(it = zPosIndex.begin(); it != zPosIndex.end(); it++){
    zPos.at((*it)) = gsl_vector_get(params, param++);
    movePlaneZ((*it), zPos.at((*it)));
  }
}

gsl_vector* EstMat::simplesStepSize(){
  //Initial step sizes for the simplex search
  int nParams = getNSimplexParams();
  gsl_vector* s = gsl_vector_alloc(nParams);

  double param = 0;
  for(size_t ii = 0; ii < resXIndex.size(); ii ++){
    gsl_vector_set(s, param++, 0.1 * resX.at(resXIndex.at(ii)) );
  }
  for(size_t ii = 0; ii < resYIndex.size(); ii ++){
    gsl_vector_set(s, param++, 0.1 * resX.at(resYIndex.at(ii)) );
  }
  for(size_t ii = 0; ii < radLengthsIndex.size(); ii ++){
    gsl_vector_set(s, param++, 0.1 * radLengths.at(radLengthsIndex.at(ii)) );
  }
  if( resXMulti.size() > 0){
    gsl_vector_set(s, param++, 1.0 );
  }
  if( resYMulti.size() > 0){
    gsl_vector_set(s, param++, 1.0 );
  }
  if( resXYMulti.size() > 0){
    gsl_vector_set(s, param++, 1.0 );
  }

  //Alignment
  for(size_t ii = 0; ii < xShiftIndex.size(); ii++){
    gsl_vector_set(s, param++, 50.0);
  }
  for(size_t ii = 0; ii < yShiftIndex.size(); ii++){
    gsl_vector_set(s, param++, 50.0);
  }
  for(size_t ii = 0; ii < xScaleIndex.size(); ii++){
    gsl_vector_set(s, param++, 0.2);
  }
  for(size_t ii = 0; ii < yScaleIndex.size(); ii++){
    gsl_vector_set(s, param++, 0.2);
  }
  for(size_t ii = 0; ii < zRotIndex.size(); ii++){
    gsl_vector_set(s, param++, 0.2);
  }
  for(size_t ii = 0; ii < zPosIndex.size(); ii++){
    gsl_vector_set(s, param++, 10000.0);
  }
  return(s);
}

double estimateSimplex(const gsl_vector* v, void * params){
  //A wrapper function for the simplex search for passing to the C library GSL
  Minimizer* minimize = (Minimizer*) params;
  EstMat& estMat = minimize->mat;
  estMat.estToSystem(v, minimize->mat.system);
  int zPrev = -999999999;
  for(size_t ii = 0; ii < estMat.system.planes.size(); ii++){
    if(zPrev > estMat.system.planes.at(ii).getZpos()){
      cout << "Bad z poses! " << ii << endl;
      //return( GSL_NAN );
    }
    zPrev = estMat.system.planes.at(ii).getZpos();
  }
  estMat.fitCount++;
  return((*minimize)());
}

void EstMat::simplexSearch(Minimizer* minimizeMe, size_t iterations, int restarts, size_t maxIterations){
  //Do the simplex search on the configured Minimizer object
  minimizeMe->init();
  
  cout << "Initial guesses" << endl;
  printAllFreeParams();
  
  if(tracks.size() < maxIterations){
    itMax = tracks.size();
  } else {
    itMax = maxIterations;
  }

  for(int ii = 0; ii < restarts; ii++){
    cout << "iteration " << ii << endl;
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;
    
    size_t iter = 0;
    int status;
    double size;
    size_t nParams = getNSimplexParams();
    
    /* Starting point */
    x = systemToEst();
    
    /* Set initial step sizes  */
    ss = simplesStepSize();
    
    /* Initialize method and iterate */
    minex_func.n = nParams;
    minex_func.f = estimateSimplex;
    minex_func.params = minimizeMe;
    
    
    s = gsl_multimin_fminimizer_alloc (T, nParams);

    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    
    do{
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      
      if (status){ break;}
      
      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-6);
      
      if (status == GSL_SUCCESS) {
	printf ("converged to minimum at\n");
      }

      if(iter % 100 == 0 or status == GSL_SUCCESS){
	cout << "Iteration " << iter << ", restart " << ii << ", nTracks = " << itMax << endl;
	printf ("%5d %10.3e %10.3e f() = %7.7f size = %.7f\n", 
		(int)iter,
		gsl_vector_get (s->x, 0), 
		gsl_vector_get (s->x, 1), 
		s->fval, size);
	printAllFreeParams();
      }
    }
    while (status == GSL_CONTINUE && iter < iterations);
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
    cout << "Status: " << status << endl;
  }
  cout << "The dataset has been fitted " << fitCount << " times." << endl;
}

FITTERTYPE EstMat::stepVector(gsl_vector* vc, size_t index, FITTERTYPE value, bool doMSE, Minimizer* minimize){
  //Obtain value for derivation step
  gsl_vector_set(vc, index, value);
  estToSystem(vc, system);
  double retVal = (*minimize)();
  if(doMSE){
    retVal = minimize->retVal2;
  }
  return(retVal);
}

void EstMat::quasiNewtonHomeMade(Minimizer* minimizeMe, int iterations){
  //Try to find minimum using independent steps of variables using Newtons method
  minimizeMe->init();
  // Successive steps in Newtons method to find the minimum
  cout << "Initial guesses" << endl;
  printAllFreeParams();

  size_t nParams = getNSimplexParams();
  gsl_vector* vc = systemToEst();
  itMax = tracks.size();
  double mseval = 0, fwbwval = 0;

  size_t resSize = resXIndex.size() + resYIndex.size();
  cout << "resSize " << resSize << endl;

  for(int iteration  = 0; iteration < iterations; iteration++){
    cout << "it " << iteration << endl;
    for(size_t ii = 0; ii < nParams; ii++){
      bool doMse = (ii < resSize) and minimizeMe->twoRetVals();

      double val = gsl_vector_get(vc, ii);
      double step = val * .04 + 0.0001; //1 percent fails for floats

      //Get five points along index i
      double p0h =  stepVector(vc, ii, val + 0.0 * step, doMse, minimizeMe);
      fwbwval = minimizeMe->result;
      mseval = minimizeMe->retVal2;
      double p1h =  stepVector(vc, ii, val + 1.0 * step, doMse, minimizeMe);
      double p2h =  stepVector(vc, ii, val + 2.0 * step, doMse, minimizeMe);
      double m1h =  stepVector(vc, ii, val - 1.0 * step, doMse, minimizeMe);
      double m2h =  stepVector(vc, ii, val - 2.0 * step, doMse, minimizeMe);

      //Go back to initial value
      gsl_vector_set(vc, ii, val);
      estToSystem(vc, system);

      double firstDeriv = ( - 1.0 * p2h + 8.0 * p1h - 8.0 * m1h + 1.0 * m2h) / (12.0 * step);
      double secondDeriv = ( -1.0 * p2h + 16 * p1h - 30.0 * p0h + 16.0 * m1h - 1.0 * m2h)/( 12.0 * step * step );
      double newVal = firstDeriv / secondDeriv;
      
      //Check the new value for NaNs, steps in wrong direction or large steps
      if(std::isnan(newVal) or std::isinf(newVal)){
	cout << "Nan step "<< endl;
      } else {
	if(secondDeriv < 0 ){ 
	  //iteration = 0;
	  cout << "Maximizing, " << val << " , " << val - newVal << endl;
	  if( fabs(newVal) > fabs( val ) * 0.1 ){
	    newVal = 0.15 * val;
	    if(firstDeriv < 0){
	      newVal *= -1.0;
	    }
	  } else {
	    newVal *= -1.0;
	  }
	}
	while(fabs(newVal) > 5 * fabs(val) + 0.00001){
	  cout << "Large step, damping: " << val << " , " << val - newVal << endl;
	  newVal *= 0.1;
	}
	gsl_vector_set(vc, ii, val - newVal );
	estToSystem(vc, system);
      }
    }
    printf ("fwbw() = %10.5f mse() = %10.5f\n", fwbwval, mseval);
    printAllFreeParams();
  }
  gsl_vector_free(vc);
}
