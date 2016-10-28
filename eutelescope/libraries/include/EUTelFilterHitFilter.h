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

namespace eutelescope {

  class EUTelFilterHitFilter : public EUTelFilter < EVENT::TrackerHit* > {
    public:
        EUTelFilterHitFilter();
//        EUTelHitFilter(const EUTelHitFilter& orig);
        virtual ~EUTelFilterHitFilter();

        void SetWantPlaneIDs(const lcio::IntVec& _wantPlaneIDs);
	lcio::IntVec GetWantPlaneIDs() const;

        virtual bool Take( const EVENT::TrackerHit* ) const;

    private:
	lcio::IntVec _wantPlaneIDs;


    };

} // eutelescope

#endif	/* EUTELHITFILTER_H */

