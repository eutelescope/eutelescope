/* 
 * 
 */

#include "EUTelFilterHitFilter.h"

// C++
#include <algorithm>
#include <functional>
#include <iterator>

// MARLIN
#include "marlin/VerbosityLevels.h"

using namespace std;
using namespace lcio;

namespace eutelescope {

    EUTelFilterHitFilter::EUTelFilterHitFilter( ) : 
        EUTelFilter< TrackerHit* >( "HitFilter", "Hit selector" ),
        _wantPlaneIDs() {
    }

//    EUTelHitFilter::EUTelHitFilter( const EUTelHitFilter& orig ) {
//    }

    EUTelFilterHitFilter::~EUTelFilterHitFilter( ) {
    }

    void EUTelFilterHitFilter::SetWantPlaneIDs( const IntVec& _wantPlaneIDs ) {
        this->_wantPlaneIDs = _wantPlaneIDs;
    }

    IntVec EUTelFilterHitFilter::GetWantPlaneIDs( ) const {
        return _wantPlaneIDs;
    }
    
    /**
     * This is an example routine that does nothing at the moment.
     * Just to show a principle.
     * @param hit
     * @return 
     */
    bool EUTelFilterHitFilter::Take(const TrackerHit* hit ) const {
        bool skip = true;
        
        // take hit if specific plane IDs are not specified
        if( this->_wantPlaneIDs.empty() ) return false;
        
        IntVec::const_iterator itr;
        for( itr = this->_wantPlaneIDs.begin(); itr != this->_wantPlaneIDs.end(); ++itr ) {
            if( hit->getCellID0() == *itr ) {
                skip = false;
                break;
            }
        }
        
        return skip;
    }
        
} // eutelescope
