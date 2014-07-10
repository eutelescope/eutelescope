//  File:   EUTelMagneticFieldFinder.h
//  Created on July 2, 2013, 12:53 PM

#ifndef EUTELMAGNETICFIELDFINDER_H
#define	EUTELMAGNETICFIELDFINDER_H

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
//LCIO
#include "lcio.h"
#include "IMPL/TrackerHitImpl.h"
#include "IMPL/TrackImpl.h"
//other
#include "streamlog/streamlog.h"
#include "gear/gearimpl/Vector3D.h"

namespace eutelescope {
    

    class MeasurementLayer {
    private:
        DISALLOW_COPY_AND_ASSIGN(MeasurementLayer)
        
    public:
        MeasurementLayer();
        
        explicit MeasurementLayer( int );
        
        virtual ~MeasurementLayer();

    public:
        /** Add hit */
        void addHit( EVENT::TrackerHit* );
        
        inline EVENT::TrackerHitVec& getHits() {
            return _allHits;
        }
        
        inline int sensorID() const {
            return _id;
        }
        
    private:
        /** Measurement layer id */
        int _id;
        /** Set of hit belonging to the measurement layer */
        EVENT::TrackerHitVec _allHits;
    };
    
    class EUTelKalmanFilter {
    private:
        DISALLOW_COPY_AND_ASSIGN(EUTelKalmanFilter) // prevent users from making (default) copies of processors

    public:
        EUTelKalmanFilter();

        explicit EUTelKalmanFilter(std::string name);

        virtual ~EUTelKalmanFilter();

        /** just print the list of tracks */
        void Print( std::string Name, std::vector< EUTelTrackImpl*> &_collection );

        /** fomr the list of _tracks remove _tracks_to_delete */
        void Prune(  std::vector< EUTelTrackImpl*> &_tracks , std::vector< EUTelTrackImpl*> &_tracks_to_delete );

        /** search for hit along track direction 
         *  using TGeo derived functions         */  
        void findTrackCandidates();

        /** Prune track candidates
         *  supposed to be removing track candidates which have n% hits in common      */  
        void PruneTrackCandidates();


        /** Initialise Fitter */
        void testUserInput();

				void testHitsVec();
				bool _hitsInputGood;
    public:

        inline std::vector< EUTelTrackImpl* >& getTracks() {
            return _tracksCartesian;
        }
                


				void	setHitsVec(EVENT::TrackerHitVec& allHitsVec){ _allHitsVec = allHitsVec;}
						 
        inline int getAllowedMissingHits() const {
            return _allowedMissingHits;
        };

        inline void setAllowedSharedHitsOnTrackCandidate( int AllowedSharedHitsOnTrackCandidate) {
            this->_AllowedSharedHitsOnTrackCandidate = AllowedSharedHitsOnTrackCandidate;
        };

        inline void setAllowedMissingHits(unsigned int allowedMissingHits) {
            this->_allowedMissingHits = allowedMissingHits;
        }

        inline int getMaxTrackCandidates() const {
            return _maxTrackCandidates;
        };

        inline void setMaxTrackCandidates(unsigned int maxTrackCandidates) {
            this->_maxTrackCandidates = maxTrackCandidates;
        }

        inline void setWindowSize(double window) {
            this->_residualsRMax = window;
        }

        inline double getWindowSize() const {
            return _residualsRMax;
        }
        
        inline void setBeamMomentum(double beam) {
            this->_beamE = beam;
        }

        inline double getBeamMomentum() const {
            return _beamE;
        }
        
        inline void setBeamMomentumUncertainty(double prec) {
            this->_beamEnergyUncertainty = prec;
        }
				inline  void setPlanesToCreateSeedsFrom(EVENT::IntVec createSeedsFromPlanes){
					this-> _createSeedsFromPlanes = createSeedsFromPlanes;
				}

				inline void setExcludePlanes(EVENT::FloatVec excludePlanes ){
					this->_excludePlanes = excludePlanes;
				}

        inline double getBeamMomentumUncertainty() const {
            return _beamEnergyUncertainty;
        }
        
        inline void setBeamCharge(double q) {
            this->_beamQ = q;
        }

        inline double getBeamCharge() const {
            return _beamQ;
        }
        
        inline void setBeamSpread( const EVENT::FloatVec& sp ) {
            this->_beamAngularSpread = sp;
        }

        inline EVENT::FloatVec getBeamSpread() const {
            return _beamAngularSpread;
        }

//Here if the user does not set a create seeds from planes x. The we set it automatically to the first plane travelling as the beam travels. 
//This has the best of both world. No reduction on functionality. User does not even know this is here. 	
void setAutoPlanestoCreateSeedsFrom(){
	if(_createSeedsFromPlanes.size() == 0){
		_createSeedsFromPlanes.push_back(geo::gGeometry().sensorZOrdertoIDs().at(0));
	}
}	
        
	/* type conversion:
	*
	**/
        double* toDouble(int n, const float * x){
          double *y = new double[n];
          for(int i=0;i<n;i++){
            y[i] = static_cast<double> (x[i]);
          }
          return y; 
        }


    private:
        
        /** prune seed track candidates */
        void pruneSeeds();

        /** Generate seed track candidates */

        /** update EUTelTrackState object at a new plane ID*/
        int findNextPlaneEntrance(  EUTelTrackStateImpl* , int  );
    		void propagateForwardFromSeedState(EUTelTrack, EUTelTrack& );

        /** a vector of hits found while swimming through the detector planes 
        * write down and dump into a collection in EUTelProcessorTrackerHelixSearch
        */
        EVENT::TrackerHitVec hitFittedVec;

    public:
				void setPlaneDimensionsVec();
				void testPlaneDimensions();
				void testHitsVecPerPlane();
				std::vector<int> _planeDimensions;
				void setHitsVecPerPlane();
				void findHitsOrderVec(LCCollection* lcCollection,EVENT::TrackerHitVec& hitsOrderVec); 
				void onlyRunOnce();
				bool _firstExecution=true;
				EVENT::IntVec _createSeedsFromPlanes;
				EVENT::FloatVec _excludePlanes;         

        void initialiseSeeds();
				void testInitialSeeds();
        /* need a method to get hitFittedVec
         * to be consistent with the other methods - passing the object by reference
         */     
				void printHits();

        EVENT::TrackerHitVec& getHitFittedVec() { 
          return hitFittedVec;
        }
 
    private:


        /** Update track state and it's cov matrix */
        double updateTrackState( EUTelTrackStateImpl*, const EVENT::TrackerHit* );

        /** Update track propagation matrix for a given step */
	const TMatrixD& getPropagationJacobianF( const EUTelTrackStateImpl*, double );
        
        /** Update Kalman gain matrix */
        TMatrixD updateGainK( const EUTelTrackStateImpl*, const EVENT::TrackerHit* );

        /** Propagate track state */
        void propagateTrackState( EUTelTrackStateImpl* );
        
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

				void nextStateUsingJacobianFinder(EUTelTrackStateImpl* input, EUTelTrackStateImpl* output, TMatrixD& jacobian);

				void UpdateStateUsingHitInformation(EUTelTrackStateImpl*,EVENT::TrackerHit* , const TMatrixD&, TMatrixD &, TMatrixD &);

				void UpdateTrackUsingHitInformation( EUTelTrackStateImpl* input,const EVENT::TrackerHit* hit, EUTelTrackImpl* track, const TMatrixD& jacobian, TMatrixD & KGain, TMatrixD & HMatrix);
        
        /** Calculate position of the track in global 
         * coordinate system for given arc length starting
         * from track's ref. point */
        TVector3 getXYZfromDzNum( const EUTelTrackStateImpl*, double ) const;

	double getXYPredictionPrecision( const EUTelTrackStateImpl* ts ) const;
        
        
        /** Get hit covariance matrix */
        TMatrixDSym getHitCov( const EVENT::TrackerHit* hit ) const;
        
        /** Get residual vector */
        TVectorD getResidual( const EUTelTrackStateImpl*, const EVENT::TrackerHit* ) const;
        
        /** Get residual covariance matrix */
        TMatrixDSym getResidualCov( const EUTelTrackStateImpl*, const EVENT::TrackerHit* hit );
        
        
        /** Convert EUTelTrackImpl to TrackImpl */
        IMPL::TrackImpl* cartesian2LCIOTrack( EUTelTrackImpl* ) const;
        
        /** Find hit closest to the track */
        const EVENT::TrackerHit* findClosestHit( const EUTelTrackStateImpl*, int );
				std::map<int ,EVENT::TrackerHitVec> _mapHitsVecPerPlane;
			protected:
				EVENT::TrackerHitVec _allHitsVec;//This is all the hits for a single event. 
    private:       
        
        /** Final set of tracks in cartesian parameterisation */
        std::vector< EUTelTrackImpl* > _tracksCartesian;
				std::map<int, std::vector<EUTelTrack> > _mapSensorIDToSeedStatesVec;
        /** Kalman track states */
        std::vector< EUTelTrackStateImpl* > _trackStates;

    private:
        /** Vector of hits to be processed */
        EVENT::TrackerHitVec _allHits;
        
        std::vector< MeasurementLayer* > _allMeasurements;

        // User supplied configuration of the fitter
    private:

        /** Validity of supplied hits */
        bool _isHitsOK;
        

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

#endif	/* EUTELMAGNETICFIELDFINDER_H */

