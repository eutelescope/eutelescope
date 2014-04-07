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

    EUTelTrackFitter::EUTelTrackFitter() : _name("DefaultTrackFitter") {
    }
/*
    EUTelTrackFitter::EUTelTrackFitter(const EUTelTrackFitter& orig) {
    }
*/
    EUTelTrackFitter::EUTelTrackFitter( std::string name ) : _name(name) {
    }

    EUTelTrackFitter::~EUTelTrackFitter() {
    }

    void EUTelTrackFitter::SetTrackCandidates( const EVENT::TrackVec& ) {
	return;
    }

    void EUTelTrackFitter::Clear() {
	return;
    }

    void EUTelTrackFitter::FitTracks() {
	return;
    }

    void EUTelTrackFitter::FitSingleTrackCandidate() {
	return;
    }


}
