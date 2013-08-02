/* 
 * File:   EUTelLCObjectTrackCandidate.cpp
 */

#include "EUTelLCObjectTrackCandidate.h"

using namespace eutelescope;

EUTelLCObjectTrackCandidate::EUTelLCObjectTrackCandidate() : IMPL::LCGenericObjectImpl(),
_trackCandidates(){
}

EUTelLCObjectTrackCandidate::EUTelLCObjectTrackCandidate(const EUTelLCObjectTrackCandidate& /*orig*/) : IMPL::LCGenericObjectImpl(),
_trackCandidates(){
}

EUTelLCObjectTrackCandidate::EUTelLCObjectTrackCandidate(const std::vector< EVENT::TrackerHitVec >& cand) : IMPL::LCGenericObjectImpl(),
_trackCandidates(cand)
{ }

EUTelLCObjectTrackCandidate::~EUTelLCObjectTrackCandidate() {
}

void EUTelLCObjectTrackCandidate::setTrackCandates( const std::vector< EVENT::TrackerHitVec >& ) {
}

std::vector< EVENT::TrackerHitVec > EUTelLCObjectTrackCandidate::getTrackCandates() const {
    return this->_trackCandidates;
}
