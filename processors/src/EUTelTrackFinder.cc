/* 
 * File:   EUTelTrackFinder.cc
 * Contact: denys.lontkovskyi@desy.de
 * 
 * Created on January 22, 2013, 2:13 PM
 */

// eutelescope includes ".h"
#include "EUTelTrackFinder.h"
#include "EUTELESCOPE.h"
//#include "EUTelUtility.h"

// lcio includes <.h>
#include <IMPL/TrackerHitImpl.h>

namespace eutelescope {
    

    EUTelTrackFinder::EUTelTrackFinder() : _isReady(kIsNotReady),
                                           _searchResult(kFailed),
                                           _name("DefaultTrackFinder"),
                                           _allHits(),
                                           _trackCandidates() { }

    EUTelTrackFinder::EUTelTrackFinder(
            std::string name ) :
                                            _isReady(kIsNotReady),
                                            _searchResult(kFailed),
                                            _name(name),
                                            _allHits(),
                                            _trackCandidates() { }

    EUTelTrackFinder::~EUTelTrackFinder() {
    }
    
    
    std::vector< EVENT::TrackerHitVec > EUTelTrackFinder::GetTrackCandidates() const {
        
        if( !_isReady ) {
            streamlog_out( ERROR ) << "Track finder " << _name << std::endl;
            streamlog_out( ERROR ) << "Tracks were not searched. " 
                                   << "Initialise track finder by supplying hits first" << std::endl;
            streamlog_out( ERROR ) << "Returning empty list of track candidates" << std::endl;
                                   
        }
        
        if( _isReady && _searchResult == kFailed ){
            streamlog_out( MESSAGE ) << "Track finder " << _name << std::endl;
            streamlog_out( MESSAGE ) << "Failed. Returning empty list of track candidates" << std::endl;
        }
        
        return _trackCandidates;
    }

    void EUTelTrackFinder::SetAllHits( const std::vector< EVENT::TrackerHitVec >& allHits) {
        
        streamlog_out( DEBUG1 ) << "EUTelTrackFinder::SetAllHits()." << std::endl;
        streamlog_out( DEBUG1 ) << " Hits for track search ready." << std::endl;
        _isReady = kIsReady;
	this->_allHits = allHits;

	return;
    }

    std::vector< EVENT::TrackerHitVec > EUTelTrackFinder::GetAllHits() const {
        return _allHits;
    }

    void EUTelTrackFinder::Reset() {
        _trackCandidates.clear();
	_searchResult = kFailed;
    }
    
    EUTelTrackFinder::SearchResult EUTelTrackFinder::SearchTracks() {
        
        if( !_isReady ) {
            streamlog_out( ERROR ) << "Track finder " << _name << std::endl;
            streamlog_out( ERROR ) << "is not ready. ..." << std::endl;
            _searchResult = kFailed;
            return _searchResult;
        }
        
        if( _allHits.empty() ) {
            streamlog_out( ERROR ) << "Track finder " << _name << std::endl;
            streamlog_out( MESSAGE4 ) << "No hits for tracks search" << std::endl;
        }
        
        streamlog_out( DEBUG1 ) << "EUTelTrackFinder::SearchTracks()" << std::endl;
        streamlog_out( DEBUG1 ) << "Looking for tracks ..." << std::endl;
        
        _searchResult = DoTrackSearch();
        
        return _searchResult;
    }
    
    EUTelTrackFinder::SearchResult EUTelTrackFinder::DoTrackSearch() {
        return kFailed;
    }

}
