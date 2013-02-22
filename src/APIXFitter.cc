// Version: $Id$
// See R.Fruhwirt: Application of Kalman Filtering to track and vertex fitting
#include "APIXFitter.h"
using namespace APIXFitter;

void TrackEstimate::printParams(bool parvec, bool covmat){
  if(parvec){
    cout << "Params: " << endl
	 << gsl_vector_get(param,0) << " , " << gsl_vector_get(param,1) << " , "
	 << gsl_vector_get(param,2) << " , " << gsl_vector_get(param,3) << endl;
  }
  if(covmat){
    cout << "Cov: " <<  endl 
	 << gsl_matrix_get(Cov,0,0) << " , " << gsl_matrix_get(Cov,0,1) << " , "
	 << gsl_matrix_get(Cov,0,2) << " , " << gsl_matrix_get(Cov,0,3) << endl
	 << gsl_matrix_get(Cov,1,0) << " , " << gsl_matrix_get(Cov,1,1) << " , "
	 << gsl_matrix_get(Cov,1,2) << " , " << gsl_matrix_get(Cov,1,3) << endl
	 << gsl_matrix_get(Cov,2,0) << " , " << gsl_matrix_get(Cov,2,1) << " , "
	 << gsl_matrix_get(Cov,2,2) << " , " << gsl_matrix_get(Cov,2,3) << endl
	 << gsl_matrix_get(Cov,3,0) << " , " << gsl_matrix_get(Cov,3,1) << " , "
	 << gsl_matrix_get(Cov,3,2) << " , " << gsl_matrix_get(Cov,3,3) << endl;
  }
}

void TrackEstimate::makeSeed(){
  //Initial state of the estimate,
  gsl_vector_set(param, 0, 0.0); gsl_vector_set(param, 1, 0.0);
  gsl_vector_set(param, 2, 0.0); gsl_vector_set(param, 3, 0.0);
  //Initial state of Cov
  gsl_matrix_set_zero(Cov);
  gsl_matrix_set(Cov,0,0,1E15); gsl_matrix_set(Cov,1,1,1E15);
  gsl_matrix_set(Cov,2,2,1E3); gsl_matrix_set(Cov,3,3,1E3);
}

APIXKalman::APIXKalman(){
  tmp4x4 = gsl_matrix_calloc(4,4); tmp4x4_2 = gsl_matrix_calloc(4,4); tmp4x4_3 = gsl_matrix_calloc(4,4);
  tmp2x2 = gsl_matrix_calloc(2,2); tmp2x2_2 = gsl_matrix_calloc(2,2);
  tmp4x2 = gsl_matrix_calloc(4,2);
  tmp2x4 = gsl_matrix_calloc(2,4);
  transM = gsl_matrix_calloc(4,4);
  kalmanGain = gsl_matrix_calloc(4,2);
  smoothGain = gsl_matrix_calloc(4,4);
  perm = gsl_permutation_alloc(4); perm2 = gsl_permutation_alloc(2);
  H = gsl_matrix_calloc(2,4);
  gsl_matrix_set(H, 0 , 0 , 1.0); gsl_matrix_set(H, 1 , 1 , 1.0);
  tmpState1 = gsl_vector_calloc(4); tmpState2 = gsl_vector_calloc(4);
  resids = gsl_vector_calloc(2);
  internalEstimate = new TrackEstimate();
}

void APIXKalman::predict(FitPlane* p1, FitPlane* p2, TrackEstimate* e) {
  //Propagating parameters
  gsl_matrix_set_identity(transM);
  double dz = p2->posZ - p1->posZ;
  gsl_matrix_set(transM, 0, 2, dz);
  gsl_matrix_set(transM, 1, 3, dz);
  gsl_blas_dgemv(CblasNoTrans,
		 1.0, transM, e->param, 
		 0.0, tmpState1);
  gsl_vector_memcpy(e->param, tmpState1);
  //Propagating errors
  // new  C  = F C F' + Q
  gsl_matrix_set(e->Cov, 2, 2, gsl_matrix_get(e->Cov, 2, 2) + p1->scatterVariance);
  gsl_matrix_set(e->Cov, 3, 3, gsl_matrix_get(e->Cov, 3, 3) + p1->scatterVariance);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, transM, e->Cov,
                  0.0, tmp4x4);
  gsl_blas_dgemm (CblasNoTrans, CblasTrans,
                  1.0, tmp4x4, transM,
                  0.0, e->Cov);
}

void APIXKalman::getKalmanGain(FitPlane* p1, TrackEstimate* e){
  //Get kalman gain matrix
  //K = CH' inv(V+ HCH')
  //HCH'
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
		 1.0, H, e->Cov, 
		 0.0, tmp2x4);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans,
		 1.0, tmp2x4, H, 
		 0.0, tmp2x2);
  //Adding V
  gsl_matrix_set(tmp2x2, 0, 0, (p1->errX * p1->errX) + gsl_matrix_get(tmp2x2, 0, 0) );
  gsl_matrix_set(tmp2x2, 1, 1, (p1->errX * p1->errX) + gsl_matrix_get(tmp2x2, 1, 1) );
  //Inverting
  int s = 0;
  gsl_linalg_LU_decomp(tmp2x2, perm2, &s);
  gsl_linalg_LU_invert(tmp2x2, perm2, tmp2x2_2);
  //CH'
  gsl_blas_dgemm(CblasNoTrans, CblasTrans,
		 1.0, e->Cov, H,
		 0.0, tmp4x2);
  //K
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
		 1.0, tmp4x2, tmp2x2_2,
		 0.0, kalmanGain);
}

void APIXKalman::update(FitPlane* p1, TrackEstimate* e){
  //If plane is excluded from fit, nothing to update
  if(p1->excluded) {  
    estimates.at(p1->index)->copy(e);
    return;
  }
  //Calculate gain matrix
  getKalmanGain(p1, e);
  //Get residuals
  gsl_vector_set(resids, 0, p1->hitPosX - gsl_vector_get(e->param, 0));
  gsl_vector_set(resids, 1, p1->hitPosY - gsl_vector_get(e->param, 1));
  //Get K(mk - Hxk|k-1)
  gsl_blas_dgemv(CblasNoTrans, 
		 1.0, kalmanGain, resids,
		 0.0, tmpState1);
  gsl_vector_add(e->param, tmpState1);
  //Update covariance matrix
  // new C = (I - KH)C 
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
		 1.0, kalmanGain, H,
		 0.0, tmp4x4);
  gsl_matrix_set_identity(tmp4x4_2);
  gsl_matrix_sub(tmp4x4_2, tmp4x4);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 
		 1.0, tmp4x4_2, e->Cov,
		 0.0, tmp4x4);
  gsl_matrix_memcpy(e->Cov, tmp4x4);
  estimates.at(p1->index)->copy(e);
}

void APIXKalman::addPlane(int index, FitPlane* plane){ 
  indexToPlane[index] = plane;
  estimates.push_back(new TrackEstimate());
}

void APIXKalman::setState(FitPlane* plane, TrackEstimate* est){
  plane->fitX = gsl_vector_get( est->param, 0);
  plane->fitY = gsl_vector_get( est->param, 1);
  plane->fitdxdz = gsl_vector_get( est->param, 2);
  plane->fitdydz = gsl_vector_get( est->param, 3);
  plane->resX = plane->fitX - plane->hitPosX;
  plane->resY = plane->fitY - plane->hitPosY;
}

void APIXKalman::getSmoothGain(FitPlane* kplus1, FitPlane* k){
  //Ak = Ck Fk' inv(Ck+1|k)
  //Ck in tmp4x4
  internalEstimate->copy( estimates.at(k->index) );
  //Ck+1|k in tmp4x4_2
  //Predict also keeps transM in tackt
  predict(k, kplus1, internalEstimate);
  gsl_matrix_memcpy(tmp4x4_2, internalEstimate->Cov);
  int s;
  gsl_linalg_LU_decomp(tmp4x4_2, perm, &s);
  gsl_linalg_LU_invert(tmp4x4_2, perm, tmp4x4_3);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 
		 1.0, transM, tmp4x4_3,
		 0.0, tmp4x4_2);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
		 1.0, estimates.at(k->index)->Cov, tmp4x4_2,
		 0.0, smoothGain);
}

void APIXKalman::getSmoothState(FitPlane* kplus1, FitPlane* k){
  //xk|n = xk|k + Ak( xk+1|n - xk+1|k)
  //xk|k
  internalEstimate->copy(estimates.at(k->index));
  //xk+1|k
  predict(k, kplus1, internalEstimate);
  //xk+1|n
  gsl_vector_memcpy(tmpState2, estimates.at(kplus1->index)->param);
  //xk+1|n - xk+1|k
  gsl_vector_sub(tmpState2, internalEstimate->param);
  gsl_blas_dgemv(CblasNoTrans, 
		 1.0, smoothGain, tmpState2, 
		 0.0, internalEstimate->param);
  gsl_vector_add(internalEstimate->param, estimates.at(k->index)->param);
  gsl_vector_memcpy(estimates.at(k->index)->param, internalEstimate->param);
}

void APIXKalman::smoothPlanes(){
  map<int, FitPlane*>::reverse_iterator it = indexToPlane.rbegin();
  FitPlane* kplus1 = (*it).second;
  setState(kplus1, estimates.at( kplus1->index));
  ++it;
  for(; it != indexToPlane.rend(); ++it){
    FitPlane* k = (*it).second;
    getSmoothGain(kplus1, k);
    getSmoothState(kplus1, k);
    setState(k, estimates.at(k->index));
    kplus1 = k;
  }
}

void APIXKalman::fitPlanes(){
  internalEstimate->makeSeed();
  map<int, FitPlane*>::iterator it = indexToPlane.begin();
  FitPlane* prev = NULL;
  for(; it != indexToPlane.end(); ++it){
    FitPlane* cur = (*it).second;
    if(not ( prev == NULL)) {
      predict(prev, cur, internalEstimate);
    }
    update(cur, internalEstimate);
    prev = cur;
  }
  smoothPlanes();
}

void APIXKalman::extrapolateBack(){
  map<int, FitPlane*>::reverse_iterator it = indexToPlane.rbegin();
  internalEstimate->copy( estimates.at(estimates.size() -1));
  setState((*it).second, internalEstimate);
  while(true){
    FitPlane* cur = (*it).second;
    FitPlane* prev = (*++it).second;
    if(it == indexToPlane.rend()){ break; }
    predict(cur, prev, internalEstimate);
    setState(prev, internalEstimate);
  }
}
