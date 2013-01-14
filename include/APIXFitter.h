// Version: $Id$
#ifndef APIXFITTER_H
#define APIXFITTER_H

#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h> 

#include <map>
#include <vector>
#include <utility>
#include <iostream>
using namespace std;

namespace APIXFitter{

  class FitPlane{
  public:
    bool excluded;
    int index, sensorID;
    double hitPosX, hitPosY, posZ;
    double errX, errY;
    double fitX, fitY, fitdxdz, fitdydz;
    double resX, resY;
    double scatterVariance;
    FitPlane(int index, int sensorID, double posZ, double errX, double errY, double scatter, bool excluded = false){
      this->index = index;
      this->sensorID = sensorID;
      this->posZ = posZ;
      this->errX = errX; this->errY = errY;
      this->scatterVariance = scatter;
      this->excluded = excluded;
      this->hitPosX = 0.0; this->hitPosY = 0.0;
    };
    void print(){
      if(excluded){ cout << "Excluded plane" ; }
      else { cout << "Active plane"; }
      cout << " with Index: " << index
	   << " sensorID: " << sensorID
	   << " posZ: " << posZ 
	   << " hitX: " << hitPosX << " hitY: " << hitPosY
	   << " errX: " << errX  << " errY: " << errY 
	   << " scatter: " << scatterVariance << endl;
    };
    void setHitpos(double x, double y){ hitPosX = x; hitPosY = y; };
  };

  class TrackEstimate{
  public:
    gsl_matrix* Cov;
    gsl_vector* param;
    TrackEstimate(){
      Cov = gsl_matrix_calloc(4,4);
      param = gsl_vector_calloc(4);
    };
    ~TrackEstimate(){
      gsl_matrix_free(Cov);
      gsl_vector_free(param);
    }
    void copy(TrackEstimate* estimate){
      gsl_matrix_memcpy ( this->Cov, estimate->Cov);
      gsl_vector_memcpy ( this->param, estimate->param);
    };
    void makeSeed();
    void printParams(bool param = true, bool cov = true);
  };

  class APIXKalman {
  private:
    //Tmp variables, avoid consing
    gsl_matrix* tmp4x4, *tmp4x4_2, *tmp4x4_3, *tmp2x2, *tmp2x2_2, *tmp4x2, *tmp2x4;
    gsl_matrix* kalmanGain, *smoothGain, *transM, *H;
    gsl_permutation *perm, *perm2;
    gsl_vector *tmpState1, *tmpState2, *resids;
    TrackEstimate* internalEstimate;

    void getKalmanGain(FitPlane* p1, TrackEstimate* e);
    void invVplusHCHT(FitPlane *p1, TrackEstimate* estimate);
    void setState(FitPlane* plane, TrackEstimate* est);
    void getSmoothGain(FitPlane* cur, FitPlane* next);
    void getSmoothState(FitPlane* cur, FitPlane* next);
  public:
    vector<TrackEstimate*> estimates;
    map<int, FitPlane* > indexToPlane;
    APIXKalman();

    void addPlane(int index, FitPlane* plane);
    FitPlane* getPlane(int index){return(indexToPlane[index]);};
    void smoothPlanes();
    void fitPlanes();
    void extrapolateBack();
    void predict(FitPlane* p1, FitPlane* p2, TrackEstimate* estimate1);
    void update(FitPlane* p1, TrackEstimate* e);
  };
}
#endif
