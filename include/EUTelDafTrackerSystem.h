// Version: $Id$
#ifndef TRACKERSYSTEM_H
#define TRACKERSYSTEM_H

// Uncomment as soon as we have a recent version of Eigen from system libraries or elsewhere
//#define EIGEN2_SUPPORT
//#include <Eigen/Eigen2Support>
#include <Eigen/Core>

#include <marlin/AIDAProcessor.h>
#include "marlin/Processor.h"
#include <list>
#include <vector>
#include <cmath>
#include <iostream>

USING_PART_OF_NAMESPACE_EIGEN

namespace daffitter{
  class TrackEstimate{
  public:
    TrackEstimate():params(),cov(){;}
    Vector4f params;
    Matrix4f cov;
    void copy(TrackEstimate* e){
      this->cov = e->cov;
      this->params = e->params;
    }
    float getX() const { return( params(0) ); }
    float getY() const { return( params(1) ); }
    float getXdz() const { return( params(2) ); }
    float getYdz() const { return( params(3) ); }
  
    float getSigmaX() const {return(  std::sqrt(cov(0,0))); }
    float getSigmaY() const {return( std::sqrt(cov(1,1))); }
    float getSigmaXdz() const {return( std::sqrt(cov(2,2))); }
    float getSigmaYdz() const {return( std::sqrt(cov(3,3))); }
  };

  class Measurement{
    Vector2f m;
    bool m_goodRegion;
    float zPos;
    size_t m_iden;
  public:
    Vector2f getM() const { return(m); }
    float getX() const { return(m(0)); }
    float getY() const { return(m(1)); }
    float getZ() const {return(zPos);}
    bool goodRegion() const { return(m_goodRegion); }
    size_t getIden() const { return(m_iden); }
    
    Measurement(float x, float y, float z, bool goodRegion, size_t iden) : m(), m_goodRegion(goodRegion), zPos(z), m_iden(iden){ 
      m(0) = x; m(1) = y;
      zPos = z;
      m_goodRegion = goodRegion;
      m_iden = iden;
    }
  };

  class TrackCandidate{
    public:
    //Information about track candidate
    //Measurement indexes for KF
    std::vector<int> indexes;
    //Weights for DAF
    std::vector< VectorXf > weights;
    //Results from fit
    float chi2, ndof;
    std::vector<TrackEstimate*> estimates;
    void print(){
      std::cout << "Track candidate with " << indexes.size() << " planes:" << std::endl;
      for(size_t ii = 0; ii < indexes.size(); ii++){
	std::cout << "Pl: " << ii << " Index " << indexes.at(ii) << " Weights: ";
	for(int jj = 0; jj < weights.at(ii).size(); jj++){
	  std::cout << weights.at(ii)(jj) << ", " ;
	}
	std::cout << std::endl;
      }
    }
    void init(int nPlanes){
      indexes.resize(nPlanes);
      weights.resize(nPlanes);
      estimates.resize(nPlanes);
      for(int ii = 0; ii < nPlanes; ii++){
	estimates.at(ii) = new TrackEstimate();
      }
    }
  };

  class FitPlane{
    // Info about a detector plane in the system
    //Unique identifier
    int sensorID;
    //Increase to angular uncertainties from scattering
    float scatterThetaSqr;
    //If excluded, plane should not be used in any fit
    bool excluded;
    float zPosition;
    float measZ;
    // Uncertainties of measurements in this plane
    Vector2f sigmas;
    Vector2f variances;
    //Sum of DAF weights for all measurements
    float sumWeights;
    //Ref point
    Vector3f ref0, ref1, ref2;
    //Norm vector
    Vector3f norm;

  public:
    //Measurements in plane
    std::vector<Measurement> meas;
    //Daf weights of measurements
    VectorXf weights;
    Vector2f invMeasVar;

    FitPlane(int sensorID, float zPos, float sigmaX, float sigmaY, float scatterVariance, bool excluded);
    int getSensorID()  const {return(this->sensorID); }
    bool isExcluded() const { return(this->excluded); }
    void exclude() { this->excluded = true; }
    void include() { this->excluded = false; }
    float getZpos()    const { return( this->zPosition); }
    float getSigmaX()  const { return(sigmas(0));}
    float getSigmaY()  const { return(sigmas(1));}
    Vector2f getSigmas() const { return(sigmas);}
    Vector2f getVars() const { return(variances);}
    void print();
    float getScatterThetaSqr() const {return(scatterThetaSqr);}
    void addMeasurement(float x, float y, float z, bool goodRegion, size_t measIden){ Measurement a(x,y, z, goodRegion, measIden); meas.push_back(a); }
    void setTotWeight(float weight){ sumWeights = weight;}
    float getTotWeight() const { return(sumWeights); };
    void clear(){ meas.clear(); measZ = zPosition;}
    float getMeasZ() const { return(measZ); }
    void setMeasZ(float z)  { measZ = z; }
    //ref points
    Vector3f& getRef0()  { return(ref0); }
    Vector3f& getRef1() { return(ref1); }
    Vector3f& getRef2() { return(ref2); }
    void setRef0(Vector3f point) { ref0 = point; }
    void setRef1(Vector3f point) { ref1 = point; }
    void setRef2(Vector3f point) { ref2 = point; }
    //Norm vector
    Vector3f& getPlaneNorm() { return(norm); }
    void setPlaneNorm(Vector3f n) { norm = n.normalized(); }
    void scaleErrors(float scaleX, float scaleY){
      sigmas(0) *= scaleX; sigmas(1) *= scaleY;
      invMeasVar(0) = 1.0f / ( sigmas(0) * sigmas(0));
      invMeasVar(1) = 1.0f / ( sigmas(1) * sigmas(1));
    }
  };
  class PlaneHit {
  private:
    Vector2f xy;
    int plane, index;
  public:
  PlaneHit(float x, float y, int plane, int index): plane(plane), index(index){ xy(0) = x; xy(1) = y; }
  PlaneHit(Vector2f xy, int plane, int index) : xy(xy), plane(plane), index(index) {}
    const Vector2f& getM() { return(xy); }
    int getPlane() const {return(plane); }
    int getIndex() const{return(index); };
    void print() {
      std::cout << "Plane " << plane << ", index " << index << ", meas " << std::endl << xy << std::endl; 
    }
  };

  class EigenFitter{
    //Eigen recommends fixed size matrixes up to 4x4
    Matrix4f transM, transMtranspose, tmp4x4, tmp4x4_2, tmp4x4_3;
    Matrix2f tmp2x2, tmp2x2_2;
    Matrix<float, 4, 2> tmp4x2, kalmanGain;
    Matrix<float, 2, 4> tmp2x4, H;
    Vector4f tmpState1, tmpState2;
    Vector2f resids, chi2s, residsum;
    Vector2f invScatterCov;
    //DAF temperature
    float tval;
  
  public:
    std::vector<TrackEstimate*> forward;
    std::vector<TrackEstimate*> backward;
    std::vector<TrackEstimate*> smoothed;
  
    EigenFitter(int nPlanes);

    //daf weights
    void setT(float tval) {this->tval = tval;};
    float getT() { return(this->tval); };
    void calculateWeights(std::vector<FitPlane> &pl, float chi2cut);
    void calculatePlaneWeight(FitPlane &pl, TrackEstimate *e, float chi2cutoff);

    //Information filter
    void predictInfo(const FitPlane &prev, const FitPlane &cur, TrackEstimate* e);
    void updateInfo(const FitPlane &pl, const int index, TrackEstimate* e);
    void updateInfoDaf(const FitPlane &pl, TrackEstimate* e);
    void getAvgInfo(TrackEstimate* e1, TrackEstimate* e2, TrackEstimate* result);
    void smoothInfo();
  };

  class TrackerSystem{
    EigenFitter* m_fitter;
    bool m_inited;
    size_t m_nTracks, m_maxCandidates, m_minClusterSize;
 
    float m_dafChi2, m_chi2OverNdof, m_sqrClusterRadius;

    float m_nXdz, m_nYdz;
    
    void getChi2Daf(daffitter::TrackCandidate *candidate);
    void getChi2Kf(daffitter::TrackCandidate *candidate);

    int addNeighbors(std::vector<PlaneHit> &candidate, std::list<PlaneHit> &hits);

    float runTweight(float t);
    float fitPlanesInfoDafInner();
    float fitPlanesInfoDafBiased();
    float getNominalXdz() const { return(m_nXdz); }
    float getNominalYdz() const { return(m_nYdz); }
    size_t getMinClusterSize() const { return(m_minClusterSize); }
    void checkNan(TrackEstimate* e);
  public:
    std::vector<daffitter::FitPlane> planes;
    std::vector<daffitter::TrackEstimate*> mcTruth;
    std::vector<daffitter::TrackCandidate*> tracks;

    TrackerSystem();
    TrackerSystem(const TrackerSystem &z);
    void addPlane(int sensorID, float zPos, float sigmaX, float sigmaY, float scatterVariance, bool excluded);
    void addMeasurement(size_t planeIndex, float x, float y, float z, bool goodRegion, size_t iden);
    void init();
    void clear();
    void setTruth(int plane, float x, float y, float xdz, float ydz);
    void setMaxCandidates(int nCandidates);
    size_t getNtracks() const { return(m_nTracks); };
    void weightToIndex(daffitter::TrackCandidate* cnd);

    //Set cut values for track finder
    void setDAFChi2Cut(float chival) { m_dafChi2 = chival;}
    float getDAFChi2Cut() const { return(m_dafChi2);};
    void setChi2OverNdofCut(float chival) { m_chi2OverNdof = chival;}
    float getChi2OverNdofCut() const { return(m_chi2OverNdof);};
    void setClusterRadius(float rad) { m_sqrClusterRadius = rad * rad;}
    void setNominalXdz(float xdz) { m_nXdz = xdz; }
    void setNominalYdz(float ydz) { m_nYdz = ydz; }
    void setMinClusterSize( size_t n) { m_minClusterSize = n; }
    void intersect();

    //Track finders
    void clusterTracker();
    void truthTracker();

    //Fitters
    void fitPlanesInfo(daffitter::TrackCandidate *candidate);
    void fitPlanesInfoDaf(daffitter::TrackCandidate*);
  };
}
#endif
