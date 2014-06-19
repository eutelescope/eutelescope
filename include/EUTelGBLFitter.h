/* 
 * File:   EUTelGBLFitter.h
 * Contact: denys.lontkovskyi@desy.de
 *
 * Created on January 25, 2013, 2:53 PM
 */

#ifdef USE_GBL

#ifndef EUTELGBLFITTER_H
#define	EUTELGBLFITTER_H

// eutelescope includes ".h"
#include "EUTelTrackFitter.h"
#include "EUTelUtility.h"
#include "EUTelUtilityRungeKutta.h"
#include "EUTelEquationsOfMotion.h"
#include "EUTelTrackStateImpl.h"
#include "EUTelTrackImpl.h"

// LCIO
#include <IMPL/LCCollectionVec.h>

// ROOT
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TMatrixD.h"
#else
#error *** You need ROOT to compile this code.  *** 
#endif

// GBL
#include "include/GblTrajectory.h"
#include "include/MilleBinary.h"

// system includes <>
#include <map>
#include <string>
#include <algorithm>
#include <functional>

namespace eutelescope {

    class EUTelGBLFitter : public EUTelTrackFitter {
        
    private:
        DISALLOW_COPY_AND_ASSIGN(EUTelGBLFitter)        // prevent users from making (default) copies of processors
      

    public:
        EUTelGBLFitter();
        explicit EUTelGBLFitter(std::string name);
        virtual ~EUTelGBLFitter();

      // do some clean up of internal data structures
      // will be automatically run when calling EUTelGBLFitter::FitTracks()
        void Clear();

        void SetTrackCandidates( vector <const IMPL::TrackImpl*> &);

        void SetTrackCandidates( const EVENT::TrackVec&);


        /** Fit tracks */
        // public: 

        void TrackCandidatesToGBLTrajectories();

				void FillInformationToGBLPointObject(EUTelTrackImpl* EUtrack, std::vector< gbl::GblPoint >* pointList);

				void FindHitIfThereIsOne(EUTelTrackImpl* trackimpl, EVENT::TrackerHit* hit, EUTelTrackStateImpl* state);

				void addMeasurementGBL(gbl::GblPoint& point, const double *hitPos, const double *statePos, const EVENT::FloatVec& hitCov, TMatrixD HMatrix);

				void CreateTrajectoryandFit(std::vector< gbl::GblPoint >* pointList,  gbl::GblTrajectory* traj, double* chi2, int* ndf);

				void CreateEUTelTrackFromTrajectory(gbl::GblTrajectory* traj, EUTelTrackImpl* EUtrackWithTrajInfo);

				void pushBackPointandBooleans( std::vector< gbl::GblPoint >* pointListTrack, const gbl::GblPoint pointTrack, bool originalState, bool hasHit  );

				void addSiPlaneScattererGBL(gbl::GblPoint& point, int iPlane);

        // private:
        /*
         */  
        void TrackCandidatesToGBLTrajectory( vector<const IMPL::TrackImpl*>::const_iterator&  );

        /*
         */
        void PerformFitGBLTrajectory( gbl::GblTrajectory* ,  vector<const IMPL::TrackImpl*>::const_iterator&, double );
     
        /* check that all trajectories are valid for Millepede
         * and dumpe the mille binary file
         */   
        bool PerformMille();

        void FitSingleTrackCandidate(EVENT::TrackVec::const_iterator& itTrkCand);
 
        inline void SetAlignmentMode( Utility::AlignmentMode alignmentMode) {
            this->_alignmentMode = alignmentMode;
        }

        inline Utility::AlignmentMode GetAlignmentMode() const {
            return _alignmentMode;
        }

        inline void SetBeamCharge(double beamQ) {
            this->_beamQ = beamQ;
        }

        inline double GetBeamCharge() const {
            return _beamQ;
        }

        inline void SetBeamEnergy(double beamE) {
            this->_eBeam = beamE;
        }

        inline double GetBeamEnergy() const {
            return _eBeam;
        }

        std::map<int, gbl::GblTrajectory*> GetGblTrackCandidates() const {
            return _gblTrackCandidates;
        }

      // return the fitted tracks
        IMPL::LCCollectionVec* GetFitTrackVec() const {
            return _fittrackvec;
        }

      // return the hits belonging to the fitted tracks
        IMPL::LCCollectionVec* GetFitHitsVec() const {
            return _fithitsvec;
        }

        void SetMilleBinary(gbl::MilleBinary* _mille) {
            this->_mille = _mille;
        }

        inline void SetChi2Cut(double chi2cut) {
            this->_chi2cut = chi2cut;
        }

        inline double GetChi2Cut() const {
            return _chi2cut;
        }
 
        void setParamterIdPlaneVec( const std::vector<int>& );
 
        void setParamterIdXResolutionVec( const std::vector<float>& );

        void setParamterIdYResolutionVec( const std::vector<float>& );
       
        void setParamterIdXRotationsMap( const std::map<int, int>& );
        
        void setParamterIdYRotationsMap( const std::map<int, int>& );
        
        void setParamterIdZRotationsMap( const std::map<int, int>& );
        
        void setParamterIdZShiftsMap( const std::map<int, int>& );
        
        void setParamterIdYShiftsMap( const std::map<int, int>& );
        
        void setParamterIdXShiftsMap(const std::map<int, int>& );
        
        const std::vector<  int>& getParamterIdPlaneVec() const;

        const std::vector<  float>& getParamterIdXResolutionVec() const;

        const std::vector<  float>& getParamterIdYResolutionVec() const;

        const std::map<int, int>& getParamterIdXRotationsMap() const;
        
        const std::map<int, int>& getParamterIdYRotationsMap() const;
        
        const std::map<int, int>& getParamterIdZRotationsMap() const;
        
        const std::map<int, int>& getParamterIdZShiftsMap() const;
        
        const std::map<int, int>& getParamterIdYShiftsMap() const;
        
        const std::map<int, int>& getParamterIdXShiftsMap() const;
        
        void setMEstimatorType( const std::string& _mEstimatorType );
        
        std::string getMEstimatorType() const;
        
        const std::map<long, int>& getHitId2GblPointLabel() const;
        
        void setExcludeFromFitPlanes(const std::vector<int>&);
        
        std::vector<int> getExcludeFromFitPlanes() const;

    private:
        void CalculateProjMatrix(TMatrixD&, double*);

        TMatrixD PropagatePar(  double, double, double, double, double, double, double );

	TVectorD getXYZfromDzNum( double, double, double, double, double, double, double ) const;
        
        TMatrixD propagatePar(double);

        void addMeasurementsGBL( gbl::GblPoint&, TVectorD&, TVectorD&, const double*, const double*, const EVENT::FloatVec&, TMatrixD& );
        
    //    void addSiPlaneScattererGBL( gbl::GblPoint&, TVectorD&, TVectorD&, int, double );
        
        void addGlobalParametersGBL( gbl::GblPoint&, TMatrixD&, std::vector<int>&, int, const double*, double, double );
        
        void pushBackPoint( std::vector< gbl::GblPoint >&, const gbl::GblPoint&, int );
        void pushBackPointMille( std::vector< gbl::GblPoint >&, const gbl::GblPoint&, int );
 
        void prepareLCIOTrack( gbl::GblTrajectory*, vector<const IMPL::TrackImpl*>::const_iterator&,
                                double, int, double, double, double, double, double );
      

        void prepareMilleOut( const IMPL::TrackImpl & );

// to be obsolete:
            void prepareLCIOTrack( gbl::GblTrajectory*, const EVENT::TrackerHitVec&,
                                double, int, double, double, double, double, double );
            void prepareMilleOut( gbl::GblTrajectory*, const EVENT::TrackVec::const_iterator& );
//

    private:
        vector<const IMPL::TrackImpl*> _trackCandidatesVec;

        EVENT::TrackVec _trackCandidates;

        std::map< int, gbl::GblTrajectory* > _gblTrackCandidates;

      // contains the fitted tracks, accessible through class methods
        IMPL::LCCollectionVec * _fittrackvec;
      // contains the fitted hits, accessible through class methods
        IMPL::LCCollectionVec * _fithitsvec;

    private:
        /** Parameter propagation jacobian */
        TMatrixD _parPropJac;


    private:
        /** Beam charge in [e] */
        double _beamQ;

        /** Beam energy in [GeV] */
        double _eBeam;

        /** Hit id to GBL point label lookup table */
        std::map<long,int> _hitId2GblPointLabel;
        std::map<long,int> _hitId2GblPointLabelMille;
       
        // Alignment 
    private:
        /** Alignment degrees of freedom */
        Utility::AlignmentMode _alignmentMode;

        /** Outlier downweighting option */
        std::string _mEstimatorType;
        
        /** Milipede binary file handle */
        gbl::MilleBinary* _mille;
        
        /** Parameter resolutions */
        std::vector<int> _parameterIdPlaneVec;
 
        /** Parameter resolutions */
        std::vector< float> _parameterIdXResolutionVec;
 
        /** Parameter resolutions */
        std::vector< float> _parameterIdYResolutionVec;

        /** Parameter ids */
        std::map<int,int> _parameterIdXShiftsMap;
        
        /** Parameter ids */
        std::map<int,int> _parameterIdYShiftsMap;
        
        /** Parameter ids */
        std::map<int,int> _parameterIdZShiftsMap;
        
        /** Parameter ids */
        std::map<int,int> _parameterIdXRotationsMap;
        
        /** Parameter ids */
        std::map<int,int> _parameterIdYRotationsMap;
        
        /** Parameter ids */
        std::map<int,int> _parameterIdZRotationsMap;
				//this maps the original states from the pattern recognition with true and the other scatterers to take into account volumes inbetween them as false. This is needed since you may not know the number of scatterers inbetween in future versions. 
				std::map<int,bool> _pointlabeltoBoolPatternState;

				std::map<int,bool> _pointlabeltoHasHit;
        
        /** Planes ids to be excluded from refit */
        std::vector< int > _excludeFromFit;
        
        // Track requirements for alignment step
        
        /** Maximum allowed track 2hi2 value*/
        double _chi2cut;

	/** ODE integrator for equations of motion */
        EUTelUtilityRungeKutta* _eomIntegrator;

	/** ODE for equations of motion */
        ODE* _eomODE;
        
    };

}
#endif	/* EUTELGBLFITTER_H */


#endif
