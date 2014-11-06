//  File:   EUTelPatternRecognition.h
//  Created on July 2, 2013, 12:53 PM

#ifndef EUTelPatternRecognition_H
#define	EUTelPatternRecognition_H

// ROOT
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#endif

// system includes <>
#include <iostream>
#include <functional>
#include <algorithm>
#include <cmath>
#include <math.h>
#include <string>
#include <vector>
#include <cmath>

// EUTELESCOPE
#include "EUTelUtility.h"
#include "EUTelUtilityRungeKutta.h"
#include "EUTelEquationsOfMotion.h"
#include "EUTelTrackFitter.h"
#include "EUTelTrackStateImpl.h"
#include "EUTelTrackImpl.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelTrack.h"
#include "EUTelState.h"
//LCIO
#include "lcio.h"
#include "IMPL/TrackerHitImpl.h"
#include "IMPL/TrackImpl.h"
#include <UTIL/LCTOOLS.h>
//other
#include "streamlog/streamlog.h"
#include "gear/gearimpl/Vector3D.h"

namespace eutelescope {
	class EUTelPatternRecognition {
	private:
		DISALLOW_COPY_AND_ASSIGN(EUTelPatternRecognition) // prevent users from making (default) copies of processors
	public:
		EUTelPatternRecognition();
		~EUTelPatternRecognition();
		//GETTERS
		std::vector<EUTelTrack>& getTracks();
		inline int getEventNumber()	const {
			return _eventNumber;
		}

		inline int getAllowedMissingHits() const {
				return _allowedMissingHits;
		};
		inline double getWindowSize() const {
				return _residualsRMax;
		}
			inline double getBeamMomentum() const {
				return _beamE;
		}
		inline double getBeamCharge() const {
				return _beamQ;
		}
		inline int getNumberOfTracksAfterPruneCut(){
			return _numberOfTracksAfterPruneCut;
		}

		//SETTERS
		void setHitsVecPerPlane();

		void	setHitsVec(EVENT::TrackerHitVec& allHitsVec){ _allHitsVec = allHitsVec;}
		void setEventNumber(int eventNumber){ 
			_eventNumber = eventNumber;
		}
		inline void setAllowedSharedHitsOnTrackCandidate( int AllowedSharedHitsOnTrackCandidate) {
				this->_AllowedSharedHitsOnTrackCandidate = AllowedSharedHitsOnTrackCandidate;
		};

		inline void setAllowedMissingHits(unsigned int allowedMissingHits) {
				this->_allowedMissingHits = allowedMissingHits;
		}

		inline void setWindowSize(double window) {
				this->_residualsRMax = window;
		}

		inline void setBeamMomentum(double beam) {
				this->_beamE = beam;
		}

		inline  void setPlanesToCreateSeedsFrom(EVENT::IntVec createSeedsFromPlanes){
			this-> _createSeedsFromPlanes = createSeedsFromPlanes;
		}

		inline void setBeamCharge(double q) {
				this->_beamQ = q;
		}
		
		//Here if the user does not set a create seeds from planes x. The we set it automatically to the first plane travelling as the beam travels. 
		//This has the best of both world. No reduction on functionality. User does not even know this is here. 	
		inline 	void setAutoPlanestoCreateSeedsFrom(){
			if(_createSeedsFromPlanes.size() == 0){
				_createSeedsFromPlanes.push_back(geo::gGeometry().sensorZOrdertoIDs().at(0));
			}
		}	
		void setPlaneDimensionsVec(EVENT::IntVec);
		//COMPUTE
		TVector3 computeInitialMomentumGlobal();
		//TEST
		void testUserInput();
		void testTrackCandidates();
		void testHitsVec();

		//OTHER
		void clearEveryRun();  
		void printTrackCandidates();
		void propagateForwardFromSeedState(EUTelState&, EUTelTrack& );
		void testPlaneDimensions();
		void testHitsVecPerPlane();
		void testPositionEstimation(float position1[], float position2[]);
		void findHitsOrderVec(LCCollection* lcCollection,EVENT::TrackerHitVec& hitsOrderVec); 
		void findTracksWithEnoughHits();
		void findTrackCandidatesWithSameHitsAndRemove();
		void findTrackCandidates(); 
		void initialiseSeeds();
		void testInitialSeeds();
		void testTrackQuality();
		void clearTrackAndTrackStates();
		void clearFinalTracks();
		//VARIABLES
		int _eventNumber;
		int _totalNumberOfHits;
		int _totalNumberOfSharedHits;
		std::map< int, int > _planeDimensions;
		bool _firstExecution;
		EVENT::IntVec _createSeedsFromPlanes;
		EVENT::FloatVec _excludePlanes;         
		std::vector<EUTelTrack> _tracks;
		std::vector<EUTelTrack> _tracksAfterEnoughHitsCut;
		std::vector<EUTelTrack>	_finalTracks;
		int _numberOfTracksTotal;
		int _numberOfTracksAfterHitCut;
		int _numberOfTracksAfterPruneCut;
		void printHits();

private:


		/** Update track state and it's cov matrix */
		double updateTrackState( EUTelTrackStateImpl*, const EVENT::TrackerHit* );
	
		/** Update track propagation matrix for a given step */
const TMatrixD& getPropagationJacobianF( const EUTelTrackStateImpl*, double );
		
		/** Update Kalman gain matrix */

		/** Propagate track state */
		void propagateTrackState( EUTelState& );
		
		/** Construct LCIO track object from internal track data */
		void prepareLCIOTrack();

		/** Sort hits according to particles propagation direction */
		bool sortHitsByMeasurementLayers( const EVENT::TrackerHitVec& );
		
		// Helper functions
private:
		
		/** Calculate track momentum from track parameters */
		TVector3 getPfromCartesianParameters( const EUTelTrackStateImpl* ) const;
		
		/** Calculate position of the track in global 
		 * coordinate system for given arc length starting
		 * from track's ref. point*/
		TVector3 getXYZfromArcLength( const EUTelTrackStateImpl*, double ) const;
		TVector3 getXYZfromArcLength1( const EUTelTrackStateImpl*, double ) const;

		void setNewState(float position[],float momentum[],  EUTelState& newState);
		

		void UpdateStateUsingHitInformation(EUTelTrackStateImpl*,EVENT::TrackerHit* , const TMatrixD&, TMatrixD &, TMatrixD &);

		void UpdateTrackUsingHitInformation( EUTelTrackStateImpl* input,const EVENT::TrackerHit* hit, EUTelTrackImpl* track, const TMatrixD& jacobian, TMatrixD & KGain, TMatrixD & HMatrix);
		
		/** Calculate position of the track in global 
		 * coordinate system for given arc length starting
		 * from track's ref. point */
		TVector3 getXYZfromDzNum( const EUTelTrackStateImpl*, double ) const;

double getXYPredictionPrecision(EUTelState& ts ) const;
		
		/** Get residual vector */
		TVectorD computeResidual(  EUTelState &, const EVENT::TrackerHit* ) const;
		
		/** Convert EUTelTrackImpl to TrackImpl */
		IMPL::TrackImpl* cartesian2LCIOTrack( EUTelTrackImpl* ) const;
		
		/** Find hit closest to the track */
		const EVENT::TrackerHit* findClosestHit(EUTelState&);
		std::map<int ,EVENT::TrackerHitVec> _mapHitsVecPerPlane;
	protected:
		EVENT::TrackerHitVec _allHitsVec;//This is all the hits for a single event. 
private:       
		
		/** Final set of tracks in cartesian parameterisation */
		std::vector< EUTelTrackImpl* > _tracksCartesian;
		std::map<int, std::vector<EUTelState> > _mapSensorIDToSeedStatesVec;

		// User supplied configuration of the fitter
private:
		/** Maximum number of sensitive planes to be considered for initial seed hits */
		int _planesForPR;

		/** Maximum number of missing on a track candidate */
		int _allowedMissingHits;
		/**Allowed # of common hits on a track for a single event.*/
		int _AllowedSharedHitsOnTrackCandidate;

		/** Maximum number of track candidates to be stored */
		int _maxTrackCandidates;

		/** Maximal distance in XY plane */
		double _residualsRMax;
		
		/** Beam momentum [GeV/c] */
		double _beamE;
		
		/** Signed beam charge [e] */
		double _beamQ;

		/** Beam energy spread [%] */
		double _beamEnergyUncertainty;
		
		/** Beam angular spread (horizontal,vertical) [mr] */
		EVENT::FloatVec _beamAngularSpread;
		
private:
/** Track parameters propagation jacobian matrix */
TMatrixD _jacobianF;
		
		/** Track parameters covariance C(k,k-1) matrix */
		TMatrixD _trkParamCovCkkm1;
		
		/** Process noise matrix */
		TMatrixDSym _processNoiseQ;       
		
		/** Kalman residual covariance matrix */
		TMatrixD _residualCovR;
		
private:
		/** ODE integrator for equations of motion */
		EUTelUtilityRungeKutta* _eomIntegrator;
		
		/** ODE integrator for Kalman propagation jacobian */
		EUTelUtilityRungeKutta* _jacobianIntegrator;
		
		/** ODE for equations of motion */
		ODE* _eomODE;
		
		/** ODE for Kalman propagation jacobian */
		ODE* _jacobianODE;
		
	};

} // namespace eutelescope

#endif	/* EUTelPatternRecognition_H */

