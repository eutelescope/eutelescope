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

        void SetTrackCandidates(const std::vector< EVENT::TrackerHitVec >&);

        /** Fit tracks */
        void FitTracks();

        inline void SetAlignmentMode( Utility::AlignmentMode alignmentMode) {
            this->_alignmentMode = alignmentMode;
        }

        inline Utility::AlignmentMode GetAlignmentMode() const {
            return _alignmentMode;
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

        IMPL::LCCollectionVec* GetFitTrackVec() const {
            return _fittrackvec;
        }

        std::vector<std::vector<gbl::GblPoint> >& GetGblTracksPoints() {
            return _gblTrackPoints;
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
        
        void SetXRotationsVec(std::vector<int>& );
        
        void SetYRotationsVec(std::vector<int>& );
        
        void SetZRotationsVec(std::vector<int>& );
        
        void SetZShiftsVec(std::vector<int>& );
        
        void SetYShiftsVec(std::vector<int>& );
        
        void SetXShiftsVec(std::vector<int>& );
        
        std::map<int, int> GetParamterIdXRotationsMap() const;
        
        std::map<int, int> GetParamterIdYRotationsMap() const;
        
        std::map<int, int> GetParamterIdZRotationsMap() const;
        
        std::map<int, int> GetParamterIdZShiftsMap() const;
        
        std::map<int, int> GetParamterIdYShiftsMap() const;
        
        std::map<int, int> GetParamterIdXShiftsMap() const;

    private:
        TMatrixD PropagatePar(double);

        //        double* GetTrackOffset( const Utility::HitsPVec& ) const;
        //        double* GetTrackSlope( const Utility::HitsPVec& ) const;

        double InterpolateTrackX(const EVENT::TrackerHitVec&, const double) const;
        double InterpolateTrackY(const EVENT::TrackerHitVec&, const double) const;
        
        double GetTrackSlopeX(const EVENT::TrackerHitVec&) const;
        double GetTrackSlopeY(const EVENT::TrackerHitVec&) const;

        void Reset();

        void AddMeasurementsGBL( gbl::GblPoint&, TVectorD&, TVectorD&, const double*, double, double, const EVENT::FloatVec&, TMatrixD& );
        
        void AddScattererGBL( gbl::GblPoint&, TVectorD&, TVectorD&, int, double );
        
        void AddGlobalParametersGBL( gbl::GblPoint&, TMatrixD&, std::vector<int>&, int, double, double, double, double );
        
    private:
        std::vector< EVENT::TrackerHitVec > _trackCandidates;

        std::map< int, gbl::GblTrajectory* > _gblTrackCandidates;

        IMPL::LCCollectionVec *_fittrackvec;

        std::vector< std::vector< gbl::GblPoint > > _gblTrackPoints;

    private:
        /** Parameter propagation jacobian */
        TMatrixD _parPropJac;


    private:
        /** Beam energy in [GeV] */
        double _eBeam;

        // Alignment 
    private:
        /** Alignment degrees of freedom */
        Utility::AlignmentMode _alignmentMode;

        /** Milipede binary file handle */
        gbl::MilleBinary* _mille;
        
        /** Parameter ids */
        std::map<int,int> _paramterIdXShiftsMap;
        
        /** Parameter ids */
        std::map<int,int> _paramterIdYShiftsMap;
        
        /** Parameter ids */
        std::map<int,int> _paramterIdZShiftsMap;
        
        /** Parameter ids */
        std::map<int,int> _paramterIdXRotationsMap;
        
        /** Parameter ids */
        std::map<int,int> _paramterIdYRotationsMap;
        
        /** Parameter ids */
        std::map<int,int> _paramterIdZRotationsMap;
        
        // Track requirements for alignment step
        
        /** Maximum allowed track 2hi2 value*/
        double _chi2cut;
        
    };

}
#endif	/* EUTELGBLFITTER_H */


#endif
