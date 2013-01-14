// Version: $Id$
#include "EUTelDafTrackerSystem.h"
#include <Eigen/Core>
#include <Eigen/Array>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <float.h>
#include <limits>

using namespace daffitter;

EigenFitter::EigenFitter(int nPlanes){
  //H maps parameter vector in to the local x-y coordinates of the plane
  H(0,0) = 1; H(1,1) =1;

  //Transportation matrix
  transM.setIdentity();
  transMtranspose = transM.transpose();

  //Storage of track estimates per plane for the forward, backward running filters, and for the final smoothed estimate
  backward.resize(nPlanes);  
  forward.resize(nPlanes);
  smoothed.resize(nPlanes);
  for(int ii = 0; ii < nPlanes; ii++){
    backward.at(ii) = new TrackEstimate();
    forward.at(ii) = new TrackEstimate();
    smoothed.at(ii) = new TrackEstimate();
  }
}

//Weight estimates
void EigenFitter::calculatePlaneWeight(FitPlane& plane, TrackEstimate* e, float chi2cutoff){
  //Calculate neasurement weights based on residuals
  size_t nMeas = plane.meas.size();
  //Nothing to do if no measurements are found
  if(nMeas < 1){
    plane.weights.resize(0);
    plane.setTotWeight(0.0f);
    return;
  }
  plane.weights.resize(nMeas);
  plane.weights.setZero();
  //Get the value exp( -chi2 / 2t) for each measurement
  for(size_t m = 0; m < nMeas ; m++){
    const Measurement &meas = plane.meas.at(m);
    resids = e->params.start<2>() - meas.getM();
//printf("%8.3f %8.3f <->%8.3f %8.3f \n", e->params(0),e->params(1), meas.getM()(0), meas.getM()(1) );

    chi2s = resids.cwise().square();
    chi2s(0) /= plane.getSigmaX() * plane.getSigmaX() + e->cov(0,0);
    chi2s(1) /= plane.getSigmaY() * plane.getSigmaY() + e->cov(1,1);
    //resids = plane.getVars() + Vector2f( e->cov(0,0), e->cov(1,1) );
    //chi2s = chi2s.cwise() / resids;
//printf("X: chi2:%8.3f  plane.sigmaX:%8.3f e-cov(00):%8.3f res(0):%8.3f\n", chi2s(0), plane.getSigmaX(), e->cov(0,0), resids(0) );
//printf("Y: chi2:%8.3f  plane.sigmaY:%8.3f e-cov(11):%8.3f res(1):%8.3f \n", chi2s(1), plane.getSigmaY(), e->cov(1,1), resids(1) );


    float chi2 = chi2s.sum();
    plane.weights(m) = exp( -1 * chi2 / (2 * tval));
//printf(" plane.weights(%4d)=%8.3f, chi2=%8.3f, tval=%8.3f \n",m,plane.weights(m),chi2, tval);
  }
  float cutWeight = exp( -1 * chi2cutoff / (2 * tval));
  plane.weights /= (cutWeight + plane.weights.sum() + FLT_MIN);
//printf("EigenFitter::calculatePlaneWeight, plane.getTotWeight=%8.3f\n", plane.getTotWeight() );

  plane.setTotWeight( plane.weights.sum() );
//printf("EigenFitter::calculatePlaneWeight, plane.getTotWeight=%8.3f\n", plane.getTotWeight() );


}

void EigenFitter::calculateWeights(std::vector<FitPlane> &planes, float chi2cut){
  //Estimate measurement weights in all planes based on smoothed track estimate
  //Track estimates should be unbiased.
  for(size_t plane = 0; plane < planes.size(); plane++){
    calculatePlaneWeight( planes.at(plane), smoothed.at(plane), chi2cut);
  }
}

//Information filter implementation of DAF
void EigenFitter::predictInfo(const FitPlane &prev, const FitPlane &cur, TrackEstimate* e){
  //float dxdz(e->params(2)), dydz(e->params(3));
  invScatterCov(0) = 1.0f / prev.getScatterThetaSqr();
  invScatterCov(1) = 1.0f /  prev.getScatterThetaSqr();

  //Add scattering to weight matrix using Woodbury matrix identity
  //inv( C + H Q H') = W - W H (inv(Q) + H' W H ) H' W
  //inv( C + H Q H')x = Wx - W H (inv(Q) + H' W H ) H' Wx
  tmp4x4.setZero();
  tmp4x4(2,2) = 1.0f / (invScatterCov(0) + e->cov(2,2));
  tmp4x4(3,3) = 1.0f / (invScatterCov(1) + e->cov(3,3));

  // //The rest
  tmp4x4_2 = e->cov * tmp4x4;
  e->cov -= tmp4x4_2 * e->cov;
  e->params -= tmp4x4_2 * e->params;
  
  //Inverse jacobian is just the oposite transformation
  transM(0,2) = transM(1,3) = transMtranspose(2,0) = transMtranspose(3,1) = prev.getMeasZ() - cur.getMeasZ();

  //New weight matrix is inv(F)' inv(C) inv(F)
  tmp4x4 = transMtranspose * e->cov;
  e->cov = tmp4x4  * transM;
  //Weigt vector bacomes
  e->params = transMtranspose * e->params;
}

void EigenFitter::updateInfo(const FitPlane &pl, const int index, TrackEstimate *e){
  if(pl.isExcluded() or (index <0)) { return;}
  const Measurement& m = pl.meas.at(index);
  //Weight matrix:
  //C = C + H'GH, with diagonal G
  e->cov.diagonal().start<2>() += pl.invMeasVar;
  //weight vector update
  // x = x + H'G M
  e->params.start<2>() += m.getM().cwise() * pl.invMeasVar;
}

void EigenFitter::updateInfoDaf(const FitPlane &pl, TrackEstimate* e){
  if(pl.isExcluded()) { return;}
  //Weight matrix:
  //C = C + H'GH
  e->cov.diagonal().start<2>() += pl.invMeasVar * pl.getTotWeight();

  //weight vector update
  // x = x + H'G M
  for(size_t ii = 0 ; ii < pl.meas.size(); ii++){
    e->params.start<2>() += pl.weights(ii) * (pl.meas.at(ii).getM().cwise() * pl.invMeasVar);
  }
}

void EigenFitter::smoothInfo(){
  for(size_t ii = 0 ; ii < smoothed.size(); ii++){
    getAvgInfo( forward.at(ii), backward.at(ii), smoothed.at(ii));
  }
}

void EigenFitter::getAvgInfo(TrackEstimate* e1, TrackEstimate* e2, TrackEstimate* result){
  tmp4x4 = e1->cov + e2->cov;
  tmp4x4.computeInverse(&result->cov);
  result->params = result->cov * ( e1->params + e2->params); 
  
  // std::cout << "Cov1:" << std::endl << e1->cov << std::endl << std::endl;
  // std::cout << "Cov2:" << std::endl << e2->cov << std::endl << std::endl;
  
  //if( isnan(result->params(0))) { std::cout << "Found nan! at " << tval << std::endl; }

  //Solve instead of invert, only gets state vector
  // result->cov = e1->cov + e2->cov;
  // result->params = e1->params + e2->params;
  // result->cov.ldlt().solveInPlace(result->params);
 
}
