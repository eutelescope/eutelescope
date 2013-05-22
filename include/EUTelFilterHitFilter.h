/** 
 * \class Implementation of a hit filter
 * \brief Performs selection of hits passing some requirements
 * 
 * With EUTelHitFilter one can select hits that satisfy conditions like
 * - hit belongs to certain plain
 * - hit is in fiducial sensor area
 */

#ifndef EUTELHITFILTER_H
#define	EUTELHITFILTER_H

// C++
#include <vector>

// LCIO
#include "lcio.h"

#include "LCIOSTLTypes.h"

#include "IMPL/TrackerHitImpl.h"
#include "IMPL/LCCollectionVec.h"

// EUTELESCOPE
#include "EUTelFilter.h"

using namespace std;
using namespace lcio;

namespace eutelescope {

    class EUTelFilterHitFilter : public EUTelFilter < TrackerHit* > {
    public:
        EUTelFilterHitFilter();
//        EUTelHitFilter(const EUTelHitFilter& orig);
        virtual ~EUTelFilterHitFilter();

        void SetWantPlaneIDs(const IntVec& _wantPlaneIDs);
        IntVec GetWantPlaneIDs() const;

        virtual bool Take( const TrackerHit* ) const;

    private:
        IntVec _wantPlaneIDs;


    };

} // eutelescope

#endif	/* EUTELHITFILTER_H */

