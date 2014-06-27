/* 
 * File:   EUTelTrackFitter.cpp
 * Contact: diont
 * 
 * Created on January 25, 2013, 6:55 PM
 */

#include "EUTelTrackFitter.h"

// system includes <>
#include <string>

namespace eutelescope {

    EUTelTrackFitter::EUTelTrackFitter() : _name("DefaultTrackFitter"), _LCIO_fittrackvec(0), _LCIO_fithitvec(0) {
 
        // reset the hit and track vectors
        Clear();

        std::cout << "default EUTelTrackFitter::EUTelTrackFitter " <<  " with " << _LCIO_fittrackvec->size() << " tracks " << std::endl;
    }
/*
    EUTelTrackFitter::EUTelTrackFitter(const EUTelTrackFitter& orig) {
    }
*/
    EUTelTrackFitter::EUTelTrackFitter( std::string name ) : _name(name), _LCIO_fittrackvec(0), _LCIO_fithitvec(0)  {
 
        // reset the hit and track vectors
        Clear();

        std::cout << " normal EUTelTrackFitter::EUTelTrackFitter " << name << " with " << _LCIO_fittrackvec->size() << " tracks " << std::endl;
    }

    EUTelTrackFitter::~EUTelTrackFitter() {
    }

    void EUTelTrackFitter::SetTrackCandidates( const EVENT::TrackVec& ) {
	return;
    }

    void EUTelTrackFitter::SetTrackCandidates( std::vector<const IMPL::TrackImpl*> & ) {
        return;
    }

    void EUTelTrackFitter::Clear() {
        
        
        if( getFitHitVec() == 0 ) 
	    try {
        	_LCIO_fithitvec  = new LCCollectionVec(LCIO::TRACKERHIT);
	        LCFlagImpl flag( _LCIO_fithitvec->getFlag());
        	flag.setBit( LCIO::TRBIT_HITS );
	        _LCIO_fithitvec->setFlag( flag.getFlag( ) );
	    } catch (...) {
        	streamlog_out(WARNING2) << "Can't allocate output _LCIO_fittracksvec collection" << std::endl;
	    }
        else
	    _LCIO_fithitvec->clear();
 
        if( getFitTrackVec() == 0 ) 
	    try {
        	_LCIO_fittrackvec  = new LCCollectionVec(LCIO::TRACK);
	        LCFlagImpl flag( _LCIO_fittrackvec->getFlag());
	       	flag.setBit( LCIO::TRBIT_HITS );
	        _LCIO_fittrackvec->setFlag( flag.getFlag( ) );
	    } catch (...) {
        	streamlog_out(WARNING2) << "Can't allocate output _LCIO_fithitvec collection" << std::endl;
	    }
        else
            _LCIO_fittrackvec->clear();

	return;
    }

    void EUTelTrackFitter::FitTracks() {
	return;
    }

    void EUTelTrackFitter::TrackCandidatesToGBLTrajectories() {
	return;
    }

    bool EUTelTrackFitter:: PerformMille() {
        return false;
    }

    void EUTelTrackFitter:: PerformFitGBLTrajectories() {
        return;
    }

    void EUTelTrackFitter::FitSingleTrackCandidate() {
	return;
    }


}
