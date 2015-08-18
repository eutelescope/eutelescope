#include "EUTelDafTrackerSystem.h"
#include <Eigen/Core>
#include <Eigen/Geometry> 
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <float.h>
#include <limits>

using namespace daffitter;

template <typename T, size_t N>
EigenFitter<T,N>::EigenFitter(){
  //Constructor, allocate memory for matrices and track estimates
  //H maps parameter vector in to the local x-y coordinates of the plane
  H(0,0) = 1; H(1,1) =1;
  //Transportation matrix
  transM.setIdentity();
  transMtranspose = transM.transpose();
}

template <typename T, size_t N>
void EigenFitter<T,N>::init(int nPlanes){
  //Resize storage of track estimates per plane for the forward, backward running filters, and for the final smoothed estimate
  backward.resize(nPlanes);  
  forward.resize(nPlanes);
  smoothed.resize(nPlanes);
}

//Weigh estimates
template <typename T, size_t N>
void EigenFitter<T,N>::calculatePlaneWeight(FitPlane<T>& plane, TrackEstimate<T,N>& e,
					    T chi2cutoff, Eigen::Matrix<T, Eigen::Dynamic, 1> &weights){
  //Calculate neasurement weights based on residuals
  size_t nMeas = plane.meas.size();
  weights.resize(nMeas);
  if(nMeas > 0) weights.setZero();
  //Get the value exp( -chi2 / 2t) for each measurement
  for(size_t m = 0; m < nMeas ; m++){
    const Measurement<T> &meas = plane.meas[m];
    resids = e.params.head(2) - meas.getM();
    chi2s = resids.array().square();
    chi2s(0) /= plane.getSigmaX() * plane.getSigmaX() + e.cov(0,0);
    chi2s(1) /= plane.getSigmaY() * plane.getSigmaY() + e.cov(1,1);
    T chi2 = chi2s.sum();
    weights(m) = exp( -1 * chi2 / (2 * tval));
    //if(isnan(plane.weights(m))){exit(1);}
  }
  T cutWeight = exp( -1 * chi2cutoff / (2 * tval));
  weights /= (cutWeight + weights.sum() + FLT_MIN);
  plane.setTotWeight( weights.sum() );
}

template <typename T, size_t N>
void EigenFitter<T,N>::calculateWeights(std::vector<FitPlane<T> > &planes, T chi2cut,
					std::vector< Eigen::Matrix<T, Eigen::Dynamic, 1> > &weights){
  //Estimate measurement weights in all planes based on smoothed track estimate
  //Track estimates should be unbiased.
  size_t nPlanes = planes.size();
  for(size_t plane = 0; plane < nPlanes; plane++){
    calculatePlaneWeight( planes[plane], smoothed[plane], chi2cut, weights.at(plane));
  }
}

//Information filter implementation of DAF
template <typename T, size_t N>
void EigenFitter<T,N>::addScatteringInfo(const FitPlane<T>& pl, TrackEstimate<T,N>& e){
  //Add scattering to weight matrix using Woodbury matrix identity
  //inv( C + H Q H') = W - W H inv(inv(Q) + H' W H ) H' W
  //inv( C + H Q H')x = Wx - W H inv(inv(Q) + H' W H ) H' Wx
 
  //Assuming diagonal Q, and no covariance between angles in the fit.
  //W H inv(inv(Q) + H' W H ) H' and e.cov are not fully populated,
  //the following takes advantage of this.
  T scattervar2 = 1.0f/(e.cov(2,2) + 1.0f/ pl.getScatterThetaSqr());
  T scattervar3 = 1.0f/(e.cov(3,3) + 1.0f/ pl.getScatterThetaSqr());
  T c20 = e.cov(2,0);
  T c31 = e.cov(3,1);
  T c22 = e.cov(2,2);
  T c33 = e.cov(3,3);
  e.cov(0,0) -= c20 * c20 * scattervar2; 
  e.cov(0,2) -= c22 * c20 * scattervar2; 
  e.cov(1,1) -= c31 * c31 * scattervar3; 
  e.cov(1,3) -= c31 * c33 * scattervar3; 
  e.cov(2,2) -= c22 * c22 * scattervar2; 
  e.cov(3,3) -= c33 * c33 * scattervar3; 
  e.cov(2,0) = e.cov(0,2);
  e.cov(3,1) = e.cov(1,3);

  T p2 = e.params(2);
  T p3 = e.params(3);
  e.params(0) -= scattervar2 * c20 * p2; 
  e.params(1) -= scattervar3 * c31 * p3; 
  e.params(2) -= scattervar2 * c22 * p2; 
  e.params(3) -= scattervar3 * c33 * p3; 
}

namespace daffitter{
  template <typename T, size_t N>
  inline void EigenFitter<T, N>::predictInfo(const FitPlane<T> &prev, const FitPlane<T> &cur, TrackEstimate<T, N>& e){
    //New weight matrix is inv(F)' inv(C) inv(F)
    //inv(F) is not fully populated, neither is inv(C), inv(C) is symmetric. The following takes advantage of this.
    T dz =  prev.getMeasZ() - cur.getMeasZ();
    T c02 = e.cov(0,2);
    T c13 = e.cov(1,3);

    e.cov(0,2) += dz * e.cov(0,0);
    e.cov(2,0) = e.cov(0,2);
    e.cov(1,3) += dz * e.cov(1,1);
    e.cov(3,1) = e.cov(1,3);
    //NOTE c02 is no longer equal to e.cov(0,2)! Same goes for c13.
    e.cov(2,2) += dz * c02 + dz * e.cov(0,2);
    e.cov(3,3) += dz * c13 + dz * e.cov(1,3);
    
    //Information vector becomes
    e.params(2) += dz * e.params(0);
    e.params(3) += dz * e.params(1);
  }
  
  template <typename T, size_t N>
  inline void EigenFitter<T,N>::updateInfo(const FitPlane<T> &pl, const int index, TrackEstimate<T,N>& e){
    //Read a measurement into the information filter
    if(pl.isExcluded() or (index <0)) { return;}
    const Measurement<T>& m = pl.meas[index];
    //Weight matrix:
    //C = C + H'GH, with diagonal G
    e.cov.diagonal().head(2) += pl.invMeasVar;
    //Information vector update
    // x = x + H'G M
    e.params.head(2).array() += m.getM().array() * pl.invMeasVar.array();
  }
  
  template <typename T, size_t N>
  inline void EigenFitter<T,N>::updateInfoDaf(const FitPlane<T> &pl, TrackEstimate<T,N>& e,
				       Eigen::Matrix<T, Eigen::Dynamic, 1> &weights){
    //Read a measurement into the weighted information filter
    if(pl.isExcluded()) { return;}
    //Weight matrix:
    //C = C + H'GH
    e.cov.diagonal().head(2) += pl.invMeasVar * pl.getTotWeight();
    
    //weight vector update
    // x = x + H'G M
    double nMeas = pl.meas.size();
    for(size_t ii = 0 ; ii < nMeas; ii++){
      e.params.head(2).array() += weights(ii) * (pl.meas[ii].getM().array() * pl.invMeasVar.array());
    }
  }
  
  template <typename T, size_t N>
  void EigenFitter<T,N>::smoothInfo(){
    //Get smoothed estimates for all planes
    double nPlanes = smoothed.size();
    for(size_t ii = 0 ; ii < nPlanes; ii++){
      getAvgInfo( forward[ii], backward[ii], smoothed[ii]);
    }
  }

  template <typename T, size_t N>
  inline void EigenFitter<T, N>::getAvgInfo(TrackEstimate<T, N>& e1, TrackEstimate<T, N>& e2, TrackEstimate<T, N>& result){
    //Get the weighted average of two estimates

    result.cov = e1.cov + e2.cov;
    fastInvert(result.cov);
    result.params = e1.params + e2.params; 
    T x(result.params(0)), y(result.params(1)), dx(result.params(2)), dy(result.params(3));
    result.params(0) = result.cov(0,0) * x + result.cov(0,2) * dx;
    result.params(1) = result.cov(1,1) * y + result.cov(1,3) * dy;
    result.params(2) = result.cov(2,0) * x + result.cov(2,2) * dx;
    result.params(3) = result.cov(3,1) * y + result.cov(3,3) * dy;
  }

  //Standard KF formulation
  template <typename T, size_t N>
  inline void EigenFitter<T, N>::predict(const FitPlane<T>  &prev, const FitPlane<T>  &cur, TrackEstimate<T, N>& e){
    //Get prediction at plane cur given plane prev

    //Setting up transportation matrix(straight line)
    T dz = cur.getZpos() - prev.getZpos();
    transM(0,2) = transM(1,3) = dz;
    
    //Propagating parameters
    tmpState1 = transM * e.params;
    e.params = tmpState1;
    
    //Propagating errors
    // C* = F (C + Q) F'    
    e.cov(2,2) += prev.getScatterThetaSqr();
    e.cov(3,3) += prev.getScatterThetaSqr();
    tmpNxN = transM * e.cov;
    e.cov = tmpNxN * transM.transpose();
  }
}
  
template <typename T, size_t N>
inline void EigenFitter<T,N>::getKFGain(const FitPlane<T>  &p1, TrackEstimate<T,N>& e){
  // Calculate KF gain matrix
  // Gain = K = CH' inv(V + HCH')
  //HCH'
  tmp2xN = H * e.cov;
  tmp2x2 = tmp2xN * H.transpose();
  // inc(HCH + V)
  tmp2x2(0,0) += p1.getSigmaX() * p1.getSigmaX();
  tmp2x2(1,1) += p1.getSigmaY() * p1.getSigmaY();
  tmp2x2.computeInverse(&tmp2x2_2);
  // CH'
  tmpNx2 = e.cov * H.transpose();
  //Result
  kalmanGain = tmpNx2 * tmp2x2_2;
}

template <typename T, size_t N>
inline void EigenFitter<T,N>::kfUpdate(const FitPlane<T>  &pl, int index, TrackEstimate<T,N>& e){
  //Include measurement index from plane p1, standard formulation
  if( index < 0 or pl.isExcluded()){ return;}
  const Measurement<T>& m = pl.meas[index];
  getKFGain(pl,e);
  //Update state
  // m - Hx
  resids = m.getM() - e.params.head(2);
  // K ( m - Hx)
  tmpState1 = kalmanGain * resids;
  //m + K( m - Hx)
  e.params += tmpState1;  
  //Update C
  //C = (I - KH) C
  tmpNxN = kalmanGain * H;
  tmpNxN_2.setIdentity();
  tmpNxN_2 -= tmpNxN;
  tmpNxN = tmpNxN_2 * e.cov;
  e.cov = tmpNxN;
}

template <typename T, size_t N>
void EigenFitter<T,N>::smooth(){
  //Smooth all planes
  double nPlanes = smoothed.size(); 
  for(size_t ii = 0 ; ii < nPlanes; ii++){
    getAvg( forward[ii], backward[ii], smoothed[ii]);
  }
}

template <typename T, size_t N>
inline void EigenFitter<T,N>::getAvg(TrackEstimate<T,N>& e1, TrackEstimate<T,N>& e2, TrackEstimate<T,N>& result){
  //Calculate the weighted average of two estimates
  if(e1.isSeed()){ result =e2; return;}
  if(e2.isSeed()){ result = e1; return;}
  //Covariance matrix of average
  //C3 = inv( inv(C1) + inv(C2)) 
  e1.cov.computeInverse(&tmpNxN);
  e2.cov.computeInverse(&tmpNxN_2);
  tmpNxN_3 = tmpNxN + tmpNxN_2;
  tmpNxN_3.computeInverse(&result.cov);
  //State of average
  // x3 = C3 ( inv(C1) x1 + inv(C2) x2)
  tmpState1 = tmpNxN * e1.params;
  tmpState2 = tmpNxN_2 * e2.params;
  tmpState1 += tmpState2;
  result.params = result.cov * tmpState1;
}
