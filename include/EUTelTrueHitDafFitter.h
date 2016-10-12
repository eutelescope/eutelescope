// based on EUTelDafBase.h

#ifndef EUTELTRUEHITDAFFITTER_H
#define EUTELTRUEHITDAFFITTER_H

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

	class EUTelTrueHitDafFitter : public marlin::Processor {

	public:

		// Marlin processor interface funtions
		// Returns a new instance of EUTelTrueHitDafFitter
		virtual Processor * newProcessor() { return new EUTelTrueHitDafFitter;}

		// Default constructor
		EUTelTrueHitDafFitter ();

		// Called at the job beginning.
		virtual void init ();
		// Called for every run.
		virtual void processRunHeader (LCRunHeader * run);
		// Called every event
		virtual void processEvent (LCEvent * evt);
		// Called after data processing.
		virtual void end();
		bool defineSystemFromData();
    
		virtual inline bool ReferenceHitVecIsSet(){ return _referenceHitVec==0; }    

		enum DafTrackFinder { simpleCluster, combinatorialKF };

	protected:

		std::ofstream trackstream;

		// Collection names
		std::string _trueHitCollectionName;
		std::string _fitpointCollectionName;
    	std::string _trackCollectionName;
		// True hit collection vector
		LCCollectionVec* _trueHitCollectionVec;

		std::vector<int> _nRef;
		bool _initializedSystem;

		std::vector<int> _colMin, _colMax, _rowMin, _rowMax;
		std::map<int, std::pair<int,int> > _rowMinMax, _colMinMax;

		// List of sensor IDs identifying telescopes and duts
		std::vector<int > _telPlanes;
		std::vector<int > _dutPlanes;

		// Resolution of sensor planes
		float _telResX, _telResY, _dutResX, _dutResY;
		// Nominal beam energy
		float _eBeam;

		// Type of track finder used (e.g. cluster or KF)
		DafTrackFinder _trackFinderType = simpleCluster;
		std::string _clusterFinderName;
   
		// Radius for track finder finder
		/* 
		*  Track finder works by projecting all hits into plane 0, assuming a beam parallel to
		*  the z-axis, then running a cluster finder on these hits. This radius determines
		*  whether a hit is included or not.
		*/
		float _normalizedRadius;

		// Cutoff value for DAF
		/*
		*  This determines the maximum distance between a track and a measurement for the
		*  measurement to be included in the fit.
		*/
		float _chi2cutoff;
		float _nXdz, _nYdz, _nXdzMaxDeviance, _nYdzMaxDeviance;
		int _nDutHits;
   
		float _nSkipMax;
		float _ndofMin;

		// Maximum allowed chi2 /ndof for track to be accepted.
		float _maxChi2;
		float _scaleScatter;

		size_t getPlaneIndex(float zPos);
		float getScatterThetaVar(float radLength);
		void readHitCollection(LCEvent* event);
		void bookHistos();
		void bookDetailedHistos();
		void fillPlots(daffitter::TrackCandidate<float,4>& track);
		void fillDetailPlots(daffitter::TrackCandidate<float,4>& track);
		bool checkTrack(daffitter::TrackCandidate<float,4>& track);
		int checkInTime(daffitter::TrackCandidate<float,4>& track);
		void printStats();

		// Alignment stuff
		void gearRotate(size_t index, size_t gearIndex);
		void getPlaneNorm(daffitter::FitPlane<float>& pl);

		daffitter::TrackerSystem<float,4> _system;
		std::map<float, int> _zSort;
		std::map<int, int> _indexIDMap;
		std::vector<float> _radLength;
		std::vector<float> _sigmaX, _sigmaY;

		// Counters
		int _iRun, _iEvt, _nTracks, _nCandidates, n_failedNdof, n_failedChi2OverNdof, n_failedIsnan, n_passedNdof, n_passedChi2OverNdof, n_passedIsnan;

		// Reference HitCollection name 
		std::string      _referenceHitCollectionName;
		std::string      _clusterCollectionName;
		bool             _useReferenceHitCollection;
		LCCollectionVec* _referenceHitVec;    
		LCCollectionVec* _clusterVec;

		// Silicon planes parameters as described in GEAR
		gear::SiPlanesParameters * _siPlanesParameters;
		gear::SiPlanesLayerLayout * _siPlanesLayerLayout;

		std::map<std::string, AIDA::IHistogram1D * > _aidaHistoMap;
		std::map<std::string, AIDA::IHistogram2D * > _aidaHistoMap2D;
		std::map<std::string, AIDA::IProfile1D * >   _aidaHistoMapProf1D;

		AIDA::IHistogram2D* _aidaZvHitX;
		AIDA::IHistogram2D* _aidaZvFitX;
		AIDA::IHistogram2D* _aidaZvHitY;
		AIDA::IHistogram2D* _aidaZvFitY;

		//2D residuals histograms
		std::map<int, AIDA::IHistogram2D*> _2DResiduals;

		// Other variables from EUTelDafFitter.h
		// Output track collection
		LCCollectionVec     * _fittrackvec;
		LCCollectionVec     * _fitpointvec;
		void addToLCIO(daffitter::TrackCandidate<float,4>& track, LCCollectionVec *lcvec);
		bool _fitDuts;
	};

	// A global instance of the processor
	EUTelTrueHitDafFitter gEUTelTrueHitDafFitter;
}
#endif
