/* 
 * File:   EUTelGBLFitter.h
 * Contact: denys.lontkovskyi@desy.de
 *
 * Created on January 25, 2013, 2:53 PM
 */

#ifndef EUTELGBLFITTER_H
#define	EUTELGBLFITTER_H

// eutelescope includes ".h"
#include "EUTelTrackFitter.h"

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

namespace eutelescope {

    class EUTelGBLFitter : public EUTelTrackFitter {
    public:
        EUTelGBLFitter();
	EUTelGBLFitter(std::string name);
        virtual ~EUTelGBLFitter();
        
        void SetTrackCandidates( std::map< int, EVENT::TrackerHitVec >& );
        void FitTracks();

	enum AlignmentMode { kXYShift, kXYShiftZRot };

        inline void SetAlignmentMode( AlignmentMode alignmentMode ) { this->_alignmentMode = alignmentMode; }
        inline AlignmentMode GetAlignmentMode() const { return _alignmentMode; }

        std::map<int, gbl::GblTrajectory*> GetGblTrackCandidates() const {
            return _gblTrackCandidates;
        }

        IMPL::LCCollectionVec* GetFitTrackVec() const {
            return _fittrackvec;
        }

        std::vector<std::vector<gbl::GblPoint> >& GetGblTracksPoints() {
            //std::vector<std::vector<gbl::GblPoint> >& gblTracksPoints = _gblTracksPoints;
            return _gblTrackPoints;
        }

        void SetMilleBinary( gbl::MilleBinary* _mille ) {
            this->_mille = _mille;
        }

    private:
	TMatrixD PropagatePar( double );
        
//        double* GetTrackOffset( const Utility::HitsPVec& ) const;
//        double* GetTrackSlope( const Utility::HitsPVec& ) const;
        
        double InterpolateTrackX( const EVENT::TrackerHitVec& , const double  ) const;
        double InterpolateTrackY( const EVENT::TrackerHitVec& , const double  ) const;

        void Reset();
        
    private:
        std::map< int, EVENT::TrackerHitVec > _trackCandidates;
        
        std::map< int, gbl::GblTrajectory* > _gblTrackCandidates;

        IMPL::LCCollectionVec     *_fittrackvec;
        
        std::vector< std::vector< gbl::GblPoint > > _gblTrackPoints;
        
	AlignmentMode _alignmentMode;

	TMatrixD _parPropJac;
        
        gbl::MilleBinary* _mille;

    };

}
#endif	/* EUTELGBLFITTER_H */

