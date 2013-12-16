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
#include "EUTelUtility.h"

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
        /** Check if track candidate satisfies selection requirements.
         *  Track candidate passes if all of it hits are in a window
         *  defined by _residualsXMin, _residualsXMax and _distanceMaxVec etc.
         *  Window size changes with distance between planes. By default 150mm
         *  is assumed.
         * 
         * @param trackCandidate vector of hits to be checked
         * @return true if track candidate satisfies all requirements
         */
        bool EUTelExhaustiveTrackFinder::IsGoodCandidate( const EVENT::TrackerHitVec& trackCandidate ) {
            bool isGood = true;

  
            // iterate over candidate's hits
            EVENT::TrackerHitVec::const_iterator itrPrevHit = trackCandidate.begin();

            const double zSpacing =  10.;   // [mm] rubinskiy 30-11-2013

            EVENT::TrackerHitVec::const_iterator itrHit;
            for ( itrHit = trackCandidate.begin(); itrHit != trackCandidate.end(); ++itrHit ) {
                
                // previous hit
                const int sensorIDPrev = Utility::GuessSensorID(  static_cast< IMPL::TrackerHitImpl* >(*itrPrevHit) );
                const double* posPrevHit = (*itrPrevHit)->getPosition();
                double posPrevHitGlob[] = {0.,0.,0.};
                geo::gGeometry().local2Master( sensorIDPrev, posPrevHit, posPrevHitGlob);
                
                // current hit
                const int sensorID = Utility::GuessSensorID( static_cast< IMPL::TrackerHitImpl* >(*itrHit) );
                const double* posHit     = (*itrHit)->getPosition();

                double posHitGlob[] = {0.,0.,0.};
                geo::gGeometry().local2Master( sensorID, posHit, posHitGlob );
//                
//                const double resX = posHit[ 0 ] - posPrevHit[ 0 ];
//                const double resY = posHit[ 1 ] - posPrevHit[ 1 ];
//                const double resZ = posHit[ 2 ] - posPrevHit[ 2 ];
//                
                const double resX = posHitGlob[ 0 ] - posPrevHitGlob[ 0 ];
                const double resY = posHitGlob[ 1 ] - posPrevHitGlob[ 1 ];
                const double resZ = posHitGlob[ 2 ] - posPrevHitGlob[ 2 ];

                const double resR = resX*resX + resY*resY;
 
                const int numberAlongZ = geo::gGeometry().sensorIDtoZOrder( sensorID );
              
                streamlog_out(DEBUG1) <<  sensorIDPrev << " to " << sensorID << " resX: " << resX << 
                  "  _residualsXMax[ "<< numberAlongZ << "]= "<< _residualsXMax[ numberAlongZ ] 
                 << " resZ : " << resZ << " zSpacing : " << zSpacing <<std::endl;

                streamlog_out(DEBUG1) <<  sensorIDPrev << " to " << sensorID << " resY: " << resY << 
                  "  _residualsYMax[ "<< numberAlongZ << "]= "<< _residualsYMax[ numberAlongZ ] 
                 << " resZ : " << resZ << " zSpacing : " << zSpacing <<std::endl;

                streamlog_out(DEBUG1) <<  sensorIDPrev << " to " << sensorID << " resR: " << resR << 
                  "  _distanceMaxVec[ "<< numberAlongZ << "]= "<< _distanceMaxVec[ numberAlongZ ] 
                 << " resZ : " << resZ << " zSpacing : " << zSpacing <<std::endl;


                if(  numberAlongZ < 0 ) return false;
                if( _mode == 1 ) {
                    if( resX > _residualsXMax[ numberAlongZ ] * resZ / zSpacing ) {
                        isGood = false;
                        break;
                    }
                    if( resX < _residualsXMin[ numberAlongZ ] * resZ / zSpacing ) {
                        isGood = false;
                        break;
                    }
                    if( resY > _residualsYMax[ numberAlongZ ] * resZ / zSpacing ) {
                        isGood = false;
                        break;
                    }
                    if( resY < _residualsYMin[ numberAlongZ ] * resZ / zSpacing ) {
                        isGood = false;
                        break;
                    }
                } else {
                    if ( sqrt( resR ) > _distanceMaxVec [ numberAlongZ ] * resZ / zSpacing  ) {
                        isGood = false;
                        break;
                    }
                }
                itrPrevHit = itrHit;
            }
            
            return isGood;
        }
          
        EUTelTrackFinder::SearchResult EUTelExhaustiveTrackFinder::DoTrackSearch() {
            streamlog_out(DEBUG1) << "EUTelExhaustiveTrackFinder::DoTrackSearch()" << std::endl;
            streamlog_out(DEBUG1) << "Looking for tracks ..." << std::endl;
           
	    _trackCandidates.clear();

            FindTracks( _allowedMissingHits, _trackCandidates, _allHits );
            
            return kSuccess;
        }
        
        void EUTelExhaustiveTrackFinder::FindTracks( int allowedmissinghits,
//                                                     EVENT::TrackVec& trackCandidates,
//                                                     EVENT::TrackVec& allHitsArray ) {
            
                                                    std::vector< EVENT::TrackerHitVec >& trackCandidates,
                                                    std::vector< EVENT::TrackerHitVec>& allHitsArray ) {
            
            const int nPlanes = geo::gGeometry().nPlanes();
            
//            if ( allHitsArray.size() != nPlanes ) return;
            
            // look for full length tracks first
            EVENT::TrackerHitVec comb;
            if( _nEmptyPlanes == 0 ) {
            CombinationGenerator st(this,allHitsArray);
                do {
                    comb = st.getCurrentCombination();
                    if( IsGoodCandidate( comb ) ) 
                        trackCandidates.push_back(comb);
//                 std::cout << " comb size: " << comb.size() << " ";
//for(int i=0;i<comb.size();i++) { std::cout<< comb[i] << " ";}
//std::cout << " trackCandidates size: " << trackCandidates.size() << std::endl;
                } while (st.incrementCurrentCombination());
           }
            // if missing hits were allowed
            // sample possible combinations of planes with missing hits
            for( int missinghits = _nEmptyPlanes + 1; missinghits <= allowedmissinghits; ++missinghits ) {
                Utility::BinomialCombination com(nPlanes, missinghits);
                std::vector < std::vector < int > > dropedPlanesCombinations;
                dropedPlanesCombinations = com.sampleCombinations();

                std::vector < int > dropedPlanes;
                // iterate over combinations of planes with missing hits
                std::vector < std::vector < int > >::const_iterator itrDropPlanes;
                for (itrDropPlanes = dropedPlanesCombinations.begin();
                        itrDropPlanes != dropedPlanesCombinations.end(); ++itrDropPlanes) {
                    dropedPlanes = *itrDropPlanes;

                    // construct all possible combinations of hits from remaining planes
                    std::vector< EVENT::TrackerHitVec> RemainigHits;

                    // remove planes with assumed missing hits
                    bool isTake;
                    for ( size_t i = 0; i < allHitsArray.size(); ++i ) {
                        isTake = true;
                        for (size_t j = 0; j < dropedPlanes.size(); ++j) {
                            if ( static_cast<int>(i) == dropedPlanes[j] ) { isTake = false; break; }
                        }
                        if ( isTake ) RemainigHits.push_back( allHitsArray[i] );
                    }

                    CombinationGenerator st(this,RemainigHits);
                    do {
                        comb = st.getCurrentCombination();
                        if( IsGoodCandidate( comb ) ) 
                            trackCandidates.push_back(comb);
                    } while (st.incrementCurrentCombination());
                } // end of loop over each droped planes combination
            } // for( int missinghits = 0; missinghits <= allowedmissinghits; ++missinghits )
        } // FindTracks()        
}
