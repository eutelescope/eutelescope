/* 
 * File:   EUTelMagneticFieldFinder.h
 *
 * Created on July 2, 2013, 12:53 PM
 */

#ifndef EUTELMAGNETICFIELDFINDER_H
#define	EUTELMAGNETICFIELDFINDER_H

// system includes <>
#include <string>
#include <vector>

// EUTELESCOPE
#include "EUTelUtility.h"
#include "EUTelTrackFitter.h"
#include "EUTelTrackStateImpl.h"
#include "EUTelTrackImpl.h"

// ROOT
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TVector3.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#endif

//LCIO
#include "lcio.h"
#include "IMPL/TrackerHitImpl.h"


class IMPL::TrackImpl;
class TrackerHit;

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
        
    private:
        /** Measurement layer id */
        int _id;
        /** Set of hit belonging to the measurement layer */
        EVENT::TrackerHitVec _allHits;
    };
    
    class EUTelKalmanFilter : public EUTelTrackFitter {
    private:
        DISALLOW_COPY_AND_ASSIGN(EUTelKalmanFilter) // prevent users from making (default) copies of processors

    public:
        EUTelKalmanFilter();

        explicit EUTelKalmanFilter(std::string name);

        virtual ~EUTelKalmanFilter();

        /** Fit supplied hits */
        virtual void FitTracks();

        /** Initialise Fitter */
        bool initialise();

        // Getters and Setters
    public:

        void setHits( EVENT::TrackerHitVec& );

        inline int getAllowedMissingHits() const {
            return _allowedMissingHits;
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

        inline void setBeamMomentum(double beam) {
            this->_beamE = beam;
        }

        inline double getBeamMomentum() const {
            return _beamE;
        }
        
        inline void setBeamMomentumUncertainty(double prec) {
            this->_beamEnergyUncertainty = prec;
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
        

    private:
        /** Flush fitter data stored from previous event */
        void reset();
        
        /** Generate seed track candidates */
        void initialiseSeeds();

        /** Find intersection point of a track with geometry planes */
        std::vector< double > findIntersection( EUTelTrackStateImpl* ts ) const;
        
        /** Propagate track state by dz */
	void propagateTrackState( EUTelTrackStateImpl*, double );

        /** Update track propagation matrix for a given step */
	const TMatrixD& getPropagationJacobianF( const EUTelTrackStateImpl*, double );

        /** Propagate track state */
        void propagateTrackState( EUTelTrackStateImpl* );
        
        /** Filter track state */
        void filterTrackState( EUTelTrackStateImpl*, TVectorD& );
        
        /** Construct LCIO track object from internal track data */
        void prepareLCIOTrack();

        /** Sort hits according to particles propagation direction */
        bool sortHitsByMeasurementLayers( const EVENT::TrackerHitVec& );
        
        // Helper functions
    private:
        
        /** Calculate track momentum from track parameters */
        TVector3 getPfromCartesianParameters( const EUTelTrackStateImpl* ) const;
        
        /** Calculate position of the track in global 
         * coordinate system for given arc length calculated
         * from track's ref. point*/
        TVector3 getXYZfromArcLenght( const EUTelTrackStateImpl*, double ) const;

	double cosAlpha( const EUTelTrackStateImpl* ) const;

	double cosBeta( const EUTelTrackStateImpl* ) const;
        
        TVectorD getTrackStateVec( const EUTelTrackStateImpl* ) const;
        
        TMatrixDSym getTrackStateCov( const EUTelTrackStateImpl* ) const;
        
        TMatrixD getHitCov( const EVENT::TrackerHit* hit ) const;
        
        TVectorD getResidual( const EUTelTrackStateImpl*, const EVENT::TrackerHit* ) const;
        
        TMatrixD getResidualCov( const EUTelTrackStateImpl*, const EVENT::TrackerHit* hit ) const;
        
        TMatrixD getH( const EUTelTrackStateImpl* ) const;

        // Kalman filter states and tracks
    private:
        /** Final set of tracks */
        std::vector< EUTelTrackImpl* > _tracks;

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
        
        /** Validity of user input flag */
        bool _isReady;

        /** Maximum number of missing on a track candidate */
        int _allowedMissingHits;

        /** Maximum number of track candidates to be stored */
        int _maxTrackCandidates;

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
        
        /** Kalman filter gain matrix */
	TMatrixD _gainK;
        
        /** Kalman residual covariance matrix */
        TMatrixD _residualCovR;
    };

} // namespace eutelescope

#endif	/* EUTELMAGNETICFIELDFINDER_H */

