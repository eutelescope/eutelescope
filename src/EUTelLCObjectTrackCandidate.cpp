/* 
 * File:   EUTelLCObjectTrackCandidate.cpp
 */

#include "EUTelLCObjectTrackCandidate.h"

using namespace eutelescope;

EUTelLCObjectTrackCandidate::EUTelLCObjectTrackCandidate() {
}

EUTelLCObjectTrackCandidate::EUTelLCObjectTrackCandidate(const EUTelLCObjectTrackCandidate& orig) {
}

EUTelLCObjectTrackCandidate::EUTelLCObjectTrackCandidate(const std::vector< EVENT::TrackerHitVec >& cand) :
_trackCandidates(cand)
{ }

EUTelLCObjectTrackCandidate::~EUTelLCObjectTrackCandidate() {
}

void EUTelLCObjectTrackCandidate::setTrackCandates( const std::vector< EVENT::TrackerHitVec >& ) {
    this->_trackCandidates;
}

std::vector< EVENT::TrackerHitVec > EUTelLCObjectTrackCandidate::getTrackCandates() const {
    return this->_trackCandidates;
}