#ifndef TRACKERSYSTEM_H
#define TRACKERSYSTEM_H

#include <Eigen/Core>
#include <list>
#include <vector>
#include <cmath>
#include <iostream>

using namespace Eigen;

namespace daffitter{
  template <typename T, size_t N>
  class TrackEstimate{
  public:
    Matrix<T, N, 1> params;
    Matrix<T, N, N> cov;
    void copy(TrackEstimate<T, N>* e){
      this->cov = e->cov;
      this->params = e->params;
    }
    void print(){
      if( cov.fullPivLu().isInvertible() ){
	std::cout << "Cov: " << "Invertible!" << std::endl;
      } else {
	std::cout << "Cov: " << "non-inverible! :(" << std::endl;
      }
      std::cout << "Params: " << std::endl << cov.inverse() * params << std::endl;
      std::cout << "Cov: " << std::endl << cov.inverse() << std::endl;
      std::cout << std::endl << params << std::endl << std::endl << cov << std::endl;
    }
    void makeSeed(bool keepState, bool freeP = false);
    void makeSeedInfo(bool freeP = false);
    bool isSeed();

    T getX() const { return( params(0) ); }
    T getY() const { return( params(1) ); }
    T getXdz() const { return( params(2) ); }
    T getYdz() const { return( params(3) ); }
    T getMomentum();

    T getSigmaX() const {return(  std::sqrt(cov(0,0))); }
    T getSigmaY() const {return( std::sqrt(cov(1,1))); }
    T getSigmaXdz() const {return( std::sqrt(cov(2,2))); }
    T getSigmaYdz() const {return( std::sqrt(cov(3,3))); }
  };

  template <typename T>
  class Measurement{
    Matrix<T, 2, 1> m;
    bool m_goodRegion;
    T zPos;
    size_t m_iden;
  public:
    Matrix<T, 2, 1> getM() const { return(m); }
    T getX() const { return(m(0)); }
    T getY() const { return(m(1)); }
    T getZ() const {return(zPos);}
    bool goodRegion() const { return(m_goodRegion); }
    size_t getIden() const { return(m_iden); }
    
    Measurement(T x, T y, T z, bool goodRegion, size_t iden){ 
      m(0) = x; m(1) = y;
      zPos = z;
      m_goodRegion = goodRegion;
      m_iden = iden;
    }
    Measurement(const Measurement<T>& m ){
      this->m = m.m;
      this->m_goodRegion = m.m_goodRegion;
      this->zPos = m.zPos;
      this->m_iden = m.m_iden;
    }
  };

  template <typename T, size_t N>
  class TrackCandidate{
    public:
    //Information about track candidate
    //Measurement indexes for KF
    std::vector<int> indexes;
    //Weights for DAF
    std::vector< Matrix<T, Eigen::Dynamic, 1> > weights;
    //Results from fit
    T chi2, ndof;
    std::vector<TrackEstimate<T,N>*> estimates;
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
	estimates.at(ii) = new TrackEstimate<T,N>();
      }
    }
  };


  template<typename T>
  class FitPlane{
    // Info about a detector plane in the system
    //Unique identifier
    int sensorID;
    //Increase to angular uncertainties from scattering
    T scatterThetaSqr;
    //If excluded, plane should not be used in any fit
    bool excluded;
    T zPosition;
    T measZ;
    // Uncertainties of measurements in this plane
    Matrix<T, 2, 1> sigmas;
    Matrix<T, 2, 1> variances;
    //Sum of DAF weights for all measurements
    T sumWeights;
    //Ref point
    Matrix<T, 3, 1> ref0, ref1, ref2;
    //Norm vector
    Matrix<T, 3, 1> norm;

  public:
    //Measurements in plane
    std::vector< Measurement<T> > meas;
    //Daf weights of measurements
    //Matrix<T, Eigen::Dynamic, 1> weights;
    Matrix<T, 2, 1> invMeasVar;
    FitPlane(int sensorID, T zPos, T sigmaX, T sigmaY, T scatterVariance, bool excluded);
    int getSensorID()  const {return(this->sensorID); }
    bool isExcluded() const { return(this->excluded); }
    void exclude() { this->excluded = true; }
    void include() { this->excluded = false; }
    T getZpos()    const { return( this->zPosition); }
    void setZpos(T zPos)    { this->zPosition = zPos; }
    T getSigmaX()  const { return(sigmas(0));}
    T getSigmaY()  const { return(sigmas(1));}
    Matrix<T, 2, 1> getSigmas() const { return(sigmas);}
    Matrix<T, 2, 1> getVars() const { return(variances);}
    void print();
    T getScatterThetaSqr() const {return(scatterThetaSqr);}
    void setScatterThetaSqr(T variance) { scatterThetaSqr = variance;}
    void addMeasurement(T x, T y, T z, bool goodRegion, size_t measIden){ Measurement<T> a(x,y, z, goodRegion, measIden); meas.push_back(a); }
    void addMeasurement(Measurement<T> m) { meas.push_back(m);}
    void setTotWeight(T weight){ sumWeights = weight;}
    T getTotWeight() const { return(sumWeights); };
    void clear(){ meas.clear(); measZ = zPosition;}
    T getMeasZ() const { return(measZ); }
    void setMeasZ(T z)  { measZ = z; }
    //ref points
    Matrix<T, 3, 1>& getRef0()  { return(ref0); }
    Matrix<T, 3, 1>& getRef1() { return(ref1); }
    Matrix<T, 3, 1>& getRef2() { return(ref2); }
    void setRef0(Matrix<T, 3, 1> point) { ref0 = point; }
    void setRef1(Matrix<T, 3, 1> point) { ref1 = point; }
    void setRef2(Matrix<T, 3, 1> point) { ref2 = point; }
    //Norm vector
    Matrix<T, 3, 1>& getPlaneNorm() { return(norm); }
    void setPlaneNorm(Matrix<T, 3, 1> n) { norm = n.normalized(); }
    void scaleErrors(T scaleX, T scaleY){
      sigmas(0) *= scaleX; sigmas(1) *= scaleY;
      invMeasVar(0) = 1.0f / ( sigmas(0) * sigmas(0));
      invMeasVar(1) = 1.0f / ( sigmas(1) * sigmas(1));
    }
    void setSigmas(T sigmaX, T sigmaY){
      sigmas(0) = sigmaX; sigmas(1) = sigmaY;
      invMeasVar(0) = 1.0f / ( sigmas(0) * sigmas(0));
      invMeasVar(1) = 1.0f / ( sigmas(1) * sigmas(1));
    }
    void setSigmaX(T sigmaX){
      sigmas(0) = sigmaX; invMeasVar(0) = 1.0/(sigmaX * sigmaX);
    }
    void setSigmaY(T sigmaY){
      sigmas(1) = sigmaY; invMeasVar(1) = 1.0/(sigmaY * sigmaY);
    }
  };
  template <typename T>
  class PlaneHit {
  private:
    Matrix<T, 2, 1> xy;
    int plane, index;
  public:
  PlaneHit(T x, T y, int plane, int index): plane(plane), index(index){ xy(0) = x; xy(1) = y; }
  PlaneHit(Matrix<T, 2, 1> xy, int plane, int index) : xy(xy), plane(plane), index(index) {}
    const Matrix<T, 2, 1>& getM() { return(xy); }
    int getPlane() const {return(plane); }
    int getIndex() const{return(index); };
    void print() {
      std::cout << "Plane " << plane << ", index " << index << ", meas " << std::endl << xy << std::endl; 
    }
  };

  template <typename T, size_t N>
  class EigenFitter{
    //Eigen recommends fixed size matrixes up to 4x4
    Matrix<T, N, N> transM, transMtranspose, tmpNxN, tmpNxN_2, tmpNxN_3;
    Matrix<T, 2, 2> tmp2x2, tmp2x2_2;
    Matrix<T, N, 2> tmpNx2, kalmanGain;
    Matrix<T, 2, N> tmp2xN, H;
    Matrix<T, N, 1> tmpState1, tmpState2;
    Matrix<T, 2, 1> resids, chi2s, residsum;
    Matrix<T, 2, 1> invScatterCov;
    //DAF temperature
    T tval;
  public:
    std::vector<TrackEstimate<T,N>*> forward;
    std::vector<TrackEstimate<T,N>*> backward;
    std::vector<TrackEstimate<T,N>*> smoothed;
  
    EigenFitter(int nPlanes);

    //daf weights
    void setT(T tval) {this->tval = tval;};
    T getT() { return(this->tval); };
    void calculateWeights(std::vector<FitPlane<T> > &pl, T chi2cut, std::vector< Matrix<T, Eigen::Dynamic, 1> > &weights);
    void calculatePlaneWeight(FitPlane<T>  &pl, TrackEstimate<T,N> *e, T chi2cutoff, Matrix<T, Eigen::Dynamic, 1> &weights);

    //Information filter
    void predictInfo(const FitPlane<T>  &prev, const FitPlane<T>  &cur, TrackEstimate<T,N>* e);
    void addScatteringInfo(const FitPlane<T> & pl, TrackEstimate<T,N>* e);
    void updateInfo(const FitPlane<T>  &pl, const int index, TrackEstimate<T,N>* e);
    void updateInfoDaf(const FitPlane<T>  &pl, TrackEstimate<T,N>* e, Matrix<T, Eigen::Dynamic, 1> &weights);
    void getAvgInfo(TrackEstimate<T,N>* e1, TrackEstimate<T,N>* e2, TrackEstimate<T,N>* result);
    void smoothInfo();
    //Standard formulation
    void predict(const FitPlane<T>  &prev, const FitPlane<T>  &cur, TrackEstimate<T,N> *e);
    //smoother
    void smooth();
    void getAvg(TrackEstimate<T,N> *f, TrackEstimate<T,N> *b, TrackEstimate<T,N> *result);
    //kf
    void kfUpdate(const FitPlane<T>  &cur, int index, TrackEstimate<T,N> *e);
    void getKFGain(const FitPlane<T>  &p1, TrackEstimate<T,N> *e);
    void predictB(const FitPlane<T>  &prev, const FitPlane<T>  &cur, TrackEstimate<T,N> *e);
  };

  template <typename T, size_t N>
  class TrackerSystem{
    //EigenFitter* m_fitter;
    bool m_inited;
    size_t m_nTracks, m_maxCandidates, m_minClusterSize;

    T m_nXdz, m_nYdz, m_nXdzdeviance, m_nYdzdeviance;
    T m_dafChi2, m_ckfChi2, m_chi2OverNdof, m_sqrClusterRadius;
    size_t m_skipMax;
    
    int addNeighbors(std::vector<PlaneHit<T> > &candidate, std::list<PlaneHit<T> > &hits);
    T runTweight(T t, daffitter::TrackCandidate<T,N>* candidate);
    T fitPlanesInfoDafInner(daffitter::TrackCandidate<T,N>* candidate);
    T fitPlanesInfoDafBiased(daffitter::TrackCandidate<T,N>* candidate);
    size_t getMinClusterSize() const { return(m_minClusterSize); }
    void checkNan(TrackEstimate<T,N>* e);
    //CKF
    void finalizeCKFTrack(TrackEstimate<T,N> *est, std::vector<int>& indexes, int nMeas, T chi2);
    void fitPermutation(int plane, TrackEstimate<T,N> *est, int nSkipped, std::vector<int> &indexes, int nMeas, T chi2);
    
  public: 
    //std::vector<int> m_chi2vals;

    EigenFitter<T,N>* m_fitter;
    std::vector<daffitter::FitPlane<T> > planes;
    std::vector<daffitter::TrackCandidate<T,N>*> tracks;

    TrackerSystem();
    TrackerSystem(const TrackerSystem<T,N>& sys);
    void addPlane(int sensorID, T zPos, T sigmaX, T sigmaY, T scatterVariance, bool excluded);
    void addMeasurement(size_t planeIndex, T x, T y, T z, bool goodRegion, size_t iden);
    void addMeasurement(Measurement<T>& meas);
    void init(bool quiet = false);
    void clear();
    void setMaxCandidates(int nCandidates);
    size_t getNtracks() const { return(m_nTracks); };
    void weightToIndex(daffitter::TrackCandidate<T,N>* cnd);
    void indexToWeight(daffitter::TrackCandidate<T,N>* cnd);
    Matrix<T, 2, 1> getBiasedResidualErrors(FitPlane<T> & pl, TrackEstimate<T,N>* estim);
    Matrix<T, 2, 1> getUnBiasedResidualErrors(FitPlane<T> & pl, TrackEstimate<T,N>* estim);
    Matrix<T, 2, 1> getResiduals(Measurement<T>& meas, TrackEstimate<T,N>* estim);
      
    //Set cut values for track finder
    void setDAFChi2Cut(T chival) { m_dafChi2 = chival;}
    void setCKFChi2Cut(T chival) { m_ckfChi2 = chival;}
    T getDAFChi2Cut() const { return(m_dafChi2);};
    T getCKFChi2Cut() const { return(m_ckfChi2);};
    void setChi2OverNdofCut(T chival) { m_chi2OverNdof = chival;}
    T getChi2OverNdofCut() const { return(m_chi2OverNdof);};
    void setClusterRadius(T rad) { m_sqrClusterRadius = rad * rad;}
    void setNominalXdz(T xdz) { m_nXdz = xdz; }
    void setNominalYdz(T ydz) { m_nYdz = ydz; }
    void setXdzMaxDeviance(T xdz) { m_nXdzdeviance = xdz; }
    void setYdzMaxDeviance(T ydz) { m_nYdzdeviance = ydz; }
    void setMinClusterSize( size_t n) { m_minClusterSize = n; }
    void intersect();
    void setMaxSkippedHits(size_t n) { m_skipMax; }

    T getNominalXdz() const { return(m_nXdz); }
    T getNominalYdz() const { return(m_nYdz); }
    T getXdzMaxDeviance() const { return(m_nXdzdeviance); }
    T getYdzMaxDeviance() const { return(m_nYdzdeviance); }

    //Track finders
    void clusterTracker();
    void truthTracker();
    void combinatorialKF();
    void index0tracker();

    //Fitters
    void fitPlanesInfoBiased(daffitter::TrackCandidate<T,N> *candidate);
    void fitPlanesInfoUnBiased(daffitter::TrackCandidate<T,N> *candidate);
    void fitPlanesInfoDaf(daffitter::TrackCandidate<T,N>*);
    void fitPlanesKF(daffitter::TrackCandidate<T,N>* candidate);
    //partial fitters
    void fitInfoFWBiased(TrackCandidate<T,N> *candidate);
    void fitInfoFWUnBiased(TrackCandidate<T,N> *candidate);
    void fitInfoBWUnBiased(TrackCandidate<T,N> *candidate);
    void fitInfoBWBiased(TrackCandidate<T,N> *candidate);
    //Getting chi2
    void getChi2Kf(daffitter::TrackCandidate<T,N> *candidate);
    void getChi2BiasedInfo(TrackCandidate<T,N> *candidate);
    void getChi2UnBiasedInfo(TrackCandidate<T,N> *candidate);
    void getChi2UnBiasedInfoDaf(TrackCandidate<T,N> *candidate);

    void smoothInfo(TrackCandidate<T,N> *candidate);
  };

  // Invert sparce matrix
  template <typename T>
  inline void partialFastInvert(Matrix<T, 4, 4> &cov, size_t p1, size_t p2){
    // Matrix is sparce, if state vector was [x,dx/dz, y,dy/dz] it would be block diagonal
    //Invert as if it was

    //Invert block
    T a(cov(p1,p1)), d(cov(p2,p2)), b(cov(p1,p2));
    T det = 1.0f / (a * d - b * b);
    cov(p1,p1) = det * d;
    cov(p2,p2) = det * a;
    cov(p2,p1) = cov(p1,p2) = det * -b;
  }
  
  template <typename T>
  inline void fastInvert(Matrix<T, 4, 4> &cov){
    // Matrix is sparce, if state vector was [x,dx/dz, y,dy/dz] it would be block diagonal
    //Invert as if it was
    // "Block" [x, dx]
    partialFastInvert(cov, 0, 2);
    // "Block" [y, dy]
    partialFastInvert(cov, 1, 3);
  }
}
#include <EUTelDafTrackerSystem.tcc>
#include <EUTelDafEigenFitter.tcc>

#endif
