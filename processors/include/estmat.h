#ifndef ESTMAT_H
#define ESTMAT_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <vector>
#include <map>
#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Cholesky>

#ifdef DOTHREAD
#include <boost/thread.hpp>
#include <boost/thread/detail/thread_group.hpp>
#include <boost/bind.hpp>
#endif

#include "EUTelDafTrackerSystem.h"
//#include "simutils.h"
#include <stdexcept>

#include <gsl/gsl_vector_double.h>

using namespace daffitter;
//typedef double FITTERTYPE;
typedef float FITTERTYPE;

class Minimizer;
class FwBw;

class EstMat{
private:
  //simplex search functions
  size_t getNSimplexParams();
  gsl_vector* systemToEst();
  gsl_vector* simplesStepSize();
  void vectToEst(gsl_vector* v1, std::vector<int>& indexVector, vector<FITTERTYPE>& dataVector, size_t& param);
  void multiVectToEst(gsl_vector* v1, std::vector<int>& indexVector, vector<FITTERTYPE>& dataVector, size_t& param);
  //newtons method
  FITTERTYPE stepVector(gsl_vector* vc, size_t index, FITTERTYPE value, bool doMSE, Minimizer* minimize);
  //data
  std::vector< std::vector<Measurement<FITTERTYPE> > > tracks;
public:
  int fitCount;
  //parameters
  std::vector<FITTERTYPE> resX, resY, radLengths, zPos, xShift, yShift, xScale, yScale, zRot;
  //single parameter iteration indexes
  std::vector<int> resXIndex, resYIndex, resXYIndex, radLengthsIndex, zPosIndex, xShiftIndex, yShiftIndex, xScaleIndex, yScaleIndex, zRotIndex;
  //multi parameter iteration indexes
  std::vector<int> resXMulti, resYMulti, resXYMulti, radLengthsMulti;
  //ebeam
  double eBeam;
  TrackerSystem<FITTERTYPE, 4> system;

  //Fake constructor
  void init(double eBeam, size_t nPlanes) {
    fitCount =0;
    this->eBeam = eBeam;
    radLengths.assign(nPlanes, 0.01);
    resX.assign(nPlanes, 4.3);
    resY.assign(nPlanes, 4.3);
    xShift.assign(nPlanes, 0.0);
    yShift.assign(nPlanes, 0.0);
    xScale.assign(nPlanes, 0.0);
    yScale.assign(nPlanes, 0.0);
    zRot.assign(nPlanes, 0.0);
    zPos.assign(nPlanes, 0.0);
  }

  //Action
  void plot(char* fname);

  //simulation
  //void initSim(int nplanes);
  //void simulate(int nTracks);

  //initialization
  void setPlane(int index, double sigmaX, double sigmaY, double radLength);
  void addTrack( std::vector<Measurement<FITTERTYPE> > track);
  void movePlaneZ(int planeIndex, double deltaZ);
  
  void estToSystem( const gsl_vector* params, TrackerSystem<FITTERTYPE, 4>& system);
  void simplexSearch(Minimizer* minimizeMe, size_t iterations, int restarts, size_t itMax = 1000000);
  void quasiNewtonHomeMade(Minimizer* minimizeMe, int iterations);
  
  int itMax;
  void readTrack(int track, TrackerSystem<FITTERTYPE,4>& system);
  void readTracksToArray(float** measX, float** measY, int nTracks, int nPlanes);
  void readTracksToDoubleArray(float** measX, int nTracks, int nPlanes);
  void clear(){ tracks.clear(); }
  void getExplicitEstimate(TrackEstimate<FITTERTYPE, 4>& estim);
  void printParams( std::string name, std::vector<FITTERTYPE>& params, bool plot, const char* valString);
  void printAllFreeParams();
};

class Minimizer{
  bool inited;
public:
  EstMat& mat;
  FITTERTYPE retVal2;
  size_t nThreads;
  FITTERTYPE result;
#ifdef DOTHREAD
  boost::mutex resultGurad;
#endif
  vector<TrackerSystem<FITTERTYPE, 4> > systems;
  
  //Minimizer(EstMat& mat) : mat(mat) {;}
  Minimizer(EstMat& mat) : inited(false), mat(mat), nThreads(4) {;}
  virtual ~Minimizer(){;};

  FITTERTYPE operator() (void);
  virtual void operator() (size_t offset, size_t stride) = 0;
  void prepareThreads();
  virtual void init ();
  virtual bool twoRetVals(){ return(false); }
};

class Chi2: public Minimizer {
public:
  Chi2(EstMat& mat) : Minimizer(mat) {;}
  virtual void operator() (size_t offset, size_t stride) ;
};

class FakeChi2: public Minimizer {
protected:
  std::vector<FITTERTYPE> resFWErrorX;
  std::vector<FITTERTYPE> resFWErrorY;
  std::vector<FITTERTYPE> resBWErrorX;
  std::vector<FITTERTYPE> resBWErrorY;
  bool firstRun;
public:
  FakeChi2(EstMat& mat) : Minimizer(mat), firstRun(false) {;}
  void calibrate(TrackerSystem<FITTERTYPE,4>& system);
  virtual void init() ;
  virtual void operator() (size_t offset, size_t stride) ;
};

class FakeAbsDev: public FakeChi2 {
public:
  FakeAbsDev(EstMat& mat) : FakeChi2(mat) {;}
  virtual void operator() (size_t offset, size_t stride) ;
};

class SDR: public Minimizer {
public:
  bool SDR1, SDR2, cholDec;
  SDR(bool SDR1, bool SDR2, bool cholDec,  EstMat& mat): Minimizer(mat), SDR1(SDR1), SDR2(SDR2), cholDec(cholDec) {;}
  virtual void operator() (size_t offset, size_t stride) ;
};

class FwBw: public Minimizer {
public:
  vector <FITTERTYPE> results2;
  FwBw(EstMat& mat): Minimizer(mat), results2(vector<FITTERTYPE>(4,0.0)) {;}
  virtual void operator() (size_t offset, size_t stride) ;
  virtual bool twoRetVals() { return(true);};
};

#endif
