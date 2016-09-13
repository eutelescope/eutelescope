// based on EUTelDafBase.h and EUTelDafFitter.h

#ifndef EUTELPROCESSORTRUEHITDAFFitter_H
#define EUTELPROCESSORTRUEHITDAFFitter_H

// eutelescope includes
#include "EUTelUtility.h"
#include "EUTelDafTrackerSystem.h"
#include "EUTelAlignmentConstant.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>

// AIDA includes <.h>
#include <AIDA/IBaseHistogram.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>

// system includes <>
#include <string>
#include <fstream>
#include <vector>
#include <map>

namespace eutelescope {

  class EUTelProcessorTrueHitDafFitter : public marlin::Processor {
  public:
    // Marlin processor interface funtions
    //! Returns a new instance of EUTelDafFitter
    virtual Processor * newProcessor() {
       return new EUTelProcessorTrueHitDafFitter;
    }
    //! Default constructor
    EUTelProcessorTrueHitDafFitter ();
    //EUTelProcessorTrueHitDafFitter (std::string);
    //! Called at the job beginning.
    virtual void init ();
    //! Called for every run.
    virtual void processRunHeader (LCRunHeader * run);
    //! Called every event
    virtual void processEvent (LCEvent * evt);
    //! Called after data processing.
    virtual void end();
    bool defineSystemFromData();
    
    virtual inline bool ReferenceHitVecIsSet(){ return _referenceHitVec==0; }    

    enum DafTrackFinder { simpleCluster, combinatorialKF };

  protected:
    std::ofstream trackstream;
    //! Input hit collection name
    std::string _reconstructedHitCollectionName;
	std::string _trueHitCollectionName;
    std::vector<std::string > _alignColNames;
    std::vector<int> _nRef;
    bool _initializedSystem;
    LCCollectionVec* _trueHitCollectionVec;
	LCCollectionVec* _reconstructedHitCollectionVec;

    EVENT::StringVec		_mcCollectionStr;
    EVENT::StringVec		_mcCollectionExample;
    LCCollection*               _mcCollection;

    std::vector<int> _colMin, _colMax, _rowMin, _rowMax;
    std::map<int, std::pair<int,int> > _rowMinMax, _colMinMax;

    // List of sensor IDs identifying telescopes and duts
    std::vector<int > _telPlanes;
    std::vector<int > _dutPlanes;

    //! resolution of sensor planes
    float _telResX, _telResY, _dutResX, _dutResY;
    //! Nominal beam energy
    float _eBeam;

	//! Type of track finder used (e.g. cluster or KF)
    //DafTrackFinder _trackFinderType = combinatorialKF;
    DafTrackFinder _trackFinderType = simpleCluster;
    std::string _clusterFinderName;
   
	//! Radius for track finder finder
    /*! 
     * Track finder works by projecting all hits into plane 0, assuming a beam parallel to
     * the z-axis, then running a cluster finder on these hits. This radius determines
     * whether a hit is included or not.
     */
    float _normalizedRadius;

    //! Cutoff value for DAF
    /*!
     * This determines the maximum distance between a track and a measurement for the
     * measurement to be included in the fit.
     */
    float _chi2cutoff;
    float _nXdz, _nYdz, _nXdzMaxDeviance, _nYdzMaxDeviance;
    int _nDutHits;
   
    float _nSkipMax;
    float _ndofMin;
    //! maximum allowed chi2 /ndof for track to be accepted.
    float _maxChi2;
    float _scaleScatter;

    //virtual void dafInit(){;}
    //virtual void dafEvent(LCEvent * /*evt*/){;}  //evt commented out because it causes a warning, function doesn't seem to do anything here but is probably used in another file through inheritance
    //virtual void dafEnd(){;}
    //virtual void dafParams(){;}
    

    size_t getPlaneIndex(float zPos);
    float getScatterThetaVar(float radLength);
    void readHitCollection(LCEvent* event);
    void bookHistos();
    void bookDetailedHistos();
	void fillHitost();
    //void fillPlots(daffitter::TrackCandidate<float,4>& track);
    //void fillDetailPlots(daffitter::TrackCandidate<float,4>& track);
    bool checkTrack(daffitter::TrackCandidate<float,4>& track);
    int checkInTime(daffitter::TrackCandidate<float,4>& track, daffitter::TrackerSystem<float, 4>& const system);
    void printStats();
    //alignment stuff
    /*void gearRotate(size_t index, size_t gearIndex, daffitter::TrackerSystem<float, 4>& const system);
    Eigen::Vector3f applyAlignment(EUTelAlignmentConstant* alignment, Eigen::Vector3f point);
    void alignRotate(std::string collectionName, LCEvent* event);
    void getPlaneNorm(daffitter::FitPlane<float>& pl);*/

	daffitter::TrackerSystem<float,4> _trueHitSystem;
    daffitter::TrackerSystem<float,4> _reconstructedHitSystem;
    std::map<float, int> _zSort;
    std::map<int, int> _indexIDMap;
    std::vector<float> _radLength;
    std::vector<float> _sigmaX, _sigmaY;

    //! Counters
    int _iRun, _iEvt, _nTracks, _nCandidates, n_failedNdof, n_failedChi2OverNdof, n_failedIsnan, n_passedNdof, n_passedChi2OverNdof, n_passedIsnan;

    //! reference HitCollection name 
    /*!
     */
    std::string      _referenceHitCollectionName;
    std::string      _clusterCollectionName;
    bool             _useReferenceHitCollection;
    LCCollectionVec* _referenceHitVec;    
    LCCollectionVec* _clusterVec;    
 
    //! Silicon planes parameters as described in GEAR
    gear::SiPlanesParameters * _siPlanesParameters;
    gear::SiPlanesLayerLayout * _siPlanesLayerLayout;

    std::map<std::string, AIDA::IHistogram1D * > _aidaHistoMap;
    std::map<std::string, AIDA::IHistogram2D * > _aidaHistoMap2D;
    std::map<std::string, AIDA::IProfile1D * >   _aidaHistoMapProf1D;
 
    AIDA::IHistogram2D* _aidaZvHitX;
    AIDA::IHistogram2D* _aidaZvFitX;
    AIDA::IHistogram2D* _aidaZvHitY;
    AIDA::IHistogram2D* _aidaZvFitY;

	// stuff from EUTelDafFitter.h
    //! Output track collection name
    std::string _trueTrackCollectionName;
    std::string _reconstructedTrackCollectionName;
    //! Output track collection
    LCCollectionVec* _trueTrackVec;
    LCCollectionVec* _trueFitPointVec;
    LCCollectionVec* _reconstructedTrackVec;
    LCCollectionVec* _reconstructedFitPointVec;
    void addTrueTrackToLCIO(daffitter::TrackCandidate<float,4>& track, LCCollectionVec* lcvec);
	void addReconstructedTrackToLCIO(daffitter::TrackCandidate<float, 4>& track, LCCollectionVec* lcvec)
    //! LCIO switch
    bool _fitDuts;
  };

	EUTelProcessorTrueHitDafFitter gEUTelProcessorTrueHitDafFitter;
}
#endif
