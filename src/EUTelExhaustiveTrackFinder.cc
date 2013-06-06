/*! 
 * File:   EUTelExhaustiveTrackFinder.cc
 * Author: diont
 * 
 * Created on January 27, 2013, 6:25 PM
 */

// eutelescope includes ".h"
#include "EUTelExhaustiveTrackFinder.h"
#include "EUTELESCOPE.h"
#include "EUTelTrackFinder.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

// lcio includes <.h>
#include "LCIOSTLTypes.h"
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>

// system includes <>
#include <map>
#include <utility>
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>

namespace eutelescope {

        void EUTelExhaustiveTrackFinder::PruneTrackCandidates( std::vector< EVENT::TrackerHitVec >& trackCandidates ) {
            std::vector< EVENT::TrackerHitVec >::iterator itr;
            for( itr = trackCandidates.begin(); itr != trackCandidates.end();  ) {
                if( (*itr).size() < 6 ) trackCandidates.erase( itr++ );
            }
            return;
        }
          
        EUTelTrackFinder::SearchResult EUTelExhaustiveTrackFinder::DoTrackSearch() {
            streamlog_out(DEBUG1) << "EUTelExhaustiveTrackFinder::DoTrackSearch()" << std::endl;
            streamlog_out(DEBUG1) << "Looking for tracks ..." << std::endl;
            // start with -1 because increment at the beginning make it 0
            EVENT::TrackerHitVec vec;
            int missinghits = 0;
	    _trackCandidates.clear();
            if( _mode == 1 ) FindTracks1( missinghits, _trackCandidates, vec, _allHits, 0, NULL );
            if( _mode == 2 ) FindTracks2( missinghits, _trackCandidates, vec, _allHits, 0, NULL );
            
            //PruneTrackCandidates( this->_trackCandidates );
            return kSuccess;
        }
        
        void EUTelExhaustiveTrackFinder::FindTracks2( int& missinghits,
                                                      std::vector< EVENT::TrackerHitVec >& trackCandidates,
                                                      EVENT::TrackerHitVec& vec,
                                                      std::map<int, EVENT::TrackerHitVec>& _allHitsArray,
                                                      int iPlane,
                                                      EVENT::TrackerHit* hitInPrevPlane) {

            streamlog_out(DEBUG0) << "EUTelExhaustiveTrackFinder::findtracks2()" << std::endl;
            streamlog_out(DEBUG0) << "Looking for tracks ..." << std::endl;
            streamlog_out(DEBUG0) << "Missing hits:" << missinghits << std::endl;

//            const EVENT::IntVec sensorIDVecZOrder = geo::gGeometry._sensorIDVecZOrder;
            
//            // recursive chain is dropped here;
//            if ( missinghits > GetAllowedMissingHits() ) {
//		 streamlog_out(DEBUG0) << "indexarray size:" << trackCandidates.size() << std::endl;
//		 return;
//	    }
            
            // recall hit id from the plane (i-1)
            if ( hitInPrevPlane != 0 ) { 
                vec.pop_back();
                vec.push_back(hitInPrevPlane);
            }

            // search for hits in the next plane if plane iPlane has no hits
            const size_t nPlanes = _allHitsArray.size();
            if ( _allHitsArray[ iPlane ].empty() && (iPlane < nPlanes - 1) ) {
                    ++missinghits;
                    if( missinghits > GetAllowedMissingHits() ) return;
		    FindTracks2(missinghits, trackCandidates, vec, _allHitsArray, iPlane + 1, NULL);
	    }

            // recursive chain is dropped here
            if ( missinghits > GetAllowedMissingHits() ) return;

            //assume missing hit in this plane
            if( (missinghits < GetAllowedMissingHits()) && (iPlane < nPlanes - 1) ) {
                int missinghits_copy = missinghits;
                ++missinghits;
                EVENT::TrackerHitVec vec_copy = vec;
                FindTracks2(missinghits, trackCandidates, vec, _allHitsArray, iPlane + 1, NULL );
                vec = vec_copy;
                missinghits = missinghits_copy;
            }
        
            // loop over hits in current plane
            const size_t nHitsInCurrentPlane = _allHitsArray[ iPlane ].size();
            for ( int j = 0; j < nHitsInCurrentPlane; j++) {

                EVENT::TrackerHit* iHit = _allHitsArray[iPlane][j];
		streamlog_out(DEBUG0) << "Hit " << j << " in plane " << iPlane << std::endl;
                vec.push_back(iHit); //index of the cluster in the last plane
                
                // if we are not in the last plane, call this method again
                if ( iPlane < nPlanes - 1 ) {
                    // track candidate requirements
                    bool takehit = true;
                    const int e = vec.size() - 2;
                    if (e >= 0) {
                        double residualX = (vec[ e ]->getPosition()[0] - vec[e + 1]->getPosition()[0]);
                        double residualY = (vec[ e ]->getPosition()[1] - vec[e + 1]->getPosition()[1]);
                        double residualZ = (vec[ e ]->getPosition()[2] - vec[e + 1]->getPosition()[2]);
                        streamlog_out(DEBUG0) << "residuals:" << std::endl;
                        streamlog_out(DEBUG0) << residualX << std::endl;
                        streamlog_out(DEBUG0) << residualY << std::endl;
                        streamlog_out(DEBUG0) << residualZ << std::endl;
                        
                        if( residualX < _residualsXMin[e] || residualX > _residualsXMax[e] ||
                            residualY < _residualsYMin[e] || residualY > _residualsYMax[e] ) {
                            takehit = false;
                        }
                    }
//                    vec.pop_back();
                    if (takehit) FindTracks2(missinghits, trackCandidates, vec, _allHitsArray, iPlane + 1, iHit);
                    
                } else { //we are in the last plane
                    //track candidate requirements
                    bool taketrack = true;
                    const int e = vec.size() - 2;
                    if (e >= 0) {
//                    for (size_t e = 0; e < vec.size() - 1; e++) {
                        double residualX = -999999.;
                        double residualY = -999999.;
                        double residualZ = -999999.;

                        // now loop through all hits on a track candidate "vec"
                        // start at the end, stop on the first non-zero hit
                        if (vec[e] != NULL ) { // non zero hit has a pointer vec[e] != NULL
                            residualX = (vec[ e ]->getPosition()[0] - vec[e + 1]->getPosition()[0]);
                            residualY = (vec[ e ]->getPosition()[1] - vec[e + 1]->getPosition()[1]);
                            residualZ = (vec[ e ]->getPosition()[2] - vec[e + 1]->getPosition()[2]);
                        }

                        if( residualX < _residualsXMin[e] || residualX > _residualsXMax[e] ||
                            residualY < _residualsYMin[e] || residualY > _residualsYMax[e] ) {
                            taketrack = false;
                        }
                    }

                    if( static_cast< int >(trackCandidates.size()) >= GetMaxTrackCandidates() ) taketrack = false;

                    //we are in the last plane. if the hit satisfies track candidate requirement then store this track candidate.
                    if( taketrack ) {
                        EVENT::TrackerHitVec aCandidate = vec;
//			IMPL::TrackImpl* track = new IMPL::TrackImpl();
//			for_each( aCandidate.begin(), aCandidate.end(), std::bind1st( std::mem_fun(&IMPL::TrackImpl::addHit), track ) );
//			trackCandidates.push_back( track );
//			delete track;
                        trackCandidates.push_back( aCandidate );
                        streamlog_out(DEBUG0) << "indexarray size at last plane:" << trackCandidates.size() << std::endl;
		    }
//                    vec.pop_back(); //last element must be removed because the
                    //vector is still used -> we are in a last plane hit loop!
                }
            }
        }

        void EUTelExhaustiveTrackFinder::FindTracks1( int& missinghits,
                                                      std::vector< EVENT::TrackerHitVec >& indexarray,
                                                      EVENT::TrackerHitVec& vec,
                                                      std::map< int, EVENT::TrackerHitVec>& _allHitsArray,
                                                      int iPlane,
                                                      EVENT::TrackerHit* y) {

            //streamlog_out(DEBUG0) << "EUTelExhaustiveTrackFinder::findtracks2()" << std::endl;
            //streamlog_out(DEBUG0) << "Looking for tracks ..." << std::endl;
            // 
            const unsigned int nPlanes = _allHitsArray.size();
            if (y == NULL) missinghits++;
            streamlog_out(DEBUG0) << "Missing hits:" << missinghits << std::endl;

            // recursive chain is dropped here;
            if (missinghits > GetAllowedMissingHits()) {
		 streamlog_out(DEBUG0) << "indexarray size:" << indexarray.size() << std::endl;
		 return;
	    }
            // recall hit id from the plane (i-1)
            if (iPlane > 0) vec.push_back(y);

            // search for hits in the next plane if plane i has no hits
            if ( _allHitsArray[ iPlane ].empty() && (iPlane < nPlanes - 1)) {
		    FindTracks1(missinghits, indexarray, vec, _allHitsArray, iPlane + 1, NULL);
	    }

           const unsigned int nHitsInCurrentPlane = _allHitsArray[ iPlane ].size();
           for (size_t j = 0; j < nHitsInCurrentPlane; j++) {
                //if we are not in the last plane, call this method again
                EVENT::TrackerHit* iHit = _allHitsArray[iPlane][j];
		streamlog_out(DEBUG0) << "Hit " << j << " in plane " << iPlane << std::endl;

                // if we are not in the last plane, call this method again
                if (iPlane < nPlanes - 1) {
                    vec.push_back(iHit); //index of the cluster in the last plane

                    //track candidate requirements
                    bool taketrack = true;
                    const int e = vec.size() - 2;
                    if (e >= 0) {
                        double distance = sqrt(
                                pow(vec[ e ]->getPosition()[0] - vec[ e + 1 ]->getPosition()[0], 2) +
                                pow(vec[ e ]->getPosition()[1] - vec[ e + 1 ]->getPosition()[1], 2)
                                );
                        double distance_z = vec[ e + 1 ]->getPosition()[2] - vec[ e ]->getPosition()[2];
                        const double dM = _distanceMaxVec[e];
                        double distancemax = dM * fabs(distance_z);

                        if (distance >= distancemax) taketrack = false;
//                        if (/*_onlySingleHitEvents == 1 && */(_allHitsArray[e].size() != 1 || _allHitsArray[e + 1].size() != 1)) taketrack = false;
                    }
                    vec.pop_back();
                    if (taketrack) FindTracks1(missinghits, indexarray, vec, _allHitsArray, iPlane + 1, iHit);
                    
                } else {
                    //we are in the last plane
                    vec.push_back(iHit); //index of the cluster in the last plane

                    //track candidate requirements
                    bool taketrack = true;
                    for (size_t e = 0; e < vec.size() - 1; e++) {
                        double distance = sqrt(
                                pow(vec[ e ]->getPosition()[0] - vec[ e + 1 ]->getPosition()[0], 2) +
                                pow(vec[ e ]->getPosition()[1] - vec[ e + 1 ]->getPosition()[1], 2)
                                );
                        double distance_z = vec[ e ]->getPosition()[2] - vec[ e + 1 ]->getPosition()[2];
                        const double dM = _distanceMaxVec[e];
                        double distancemax = dM * fabs(distance_z);

                        if (distance >= distancemax) taketrack = false;
//                        if (_allHitsArray[e].size() != 1 || _allHitsArray[e + 1].size() != 1) taketrack = false;    //only one hit in each planes

                    }
                    if( static_cast< int >(indexarray.size()) >= GetMaxTrackCandidates() ) taketrack = false;

                    if (taketrack) {
                            EVENT::TrackerHitVec aCandidate = vec;         //copy vec, because it will be reused
                            indexarray.push_back( aCandidate );
                            streamlog_out(DEBUG0) << "indexarray size at last plane:" << indexarray.size() << std::endl;
                    }
                    vec.pop_back(); //last element must be removed because the
                    //vector is still used
                }
            }
            vec.pop_back();
    }
        
}
