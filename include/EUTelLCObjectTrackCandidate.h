/* 
 * File:   EUTelLCObjectTrackCandidate.h
 */

#ifndef EUTELLCOBJECTTRACKCANDIDATE_H
#define	EUTELLCOBJECTTRACKCANDIDATE_H

// C++
#include <vector>

// LCIO
#include "IMPL/TrackerHitImpl.h"
#include "IMPL/LCGenericObjectImpl.h"

namespace eutelescope {

    /** @class EUTelLCObjectTrackCandidate LCGenericObjectImpl wrapper for storing 
     *  track candidates hits collection from pattern recognition step
     * 
     */
    class EUTelLCObjectTrackCandidate : public IMPL::LCGenericObjectImpl {
    public:

        /** Default constructor */
        EUTelLCObjectTrackCandidate();

        /** Copy constructor */
        EUTelLCObjectTrackCandidate(const EUTelLCObjectTrackCandidate& orig);
        
        /** Collection creation constructor */
        EUTelLCObjectTrackCandidate(const std::vector< EVENT::TrackerHitVec >&);  

        /** Destructor */
        virtual ~EUTelLCObjectTrackCandidate();


        // Setters and getters 
    public:

        void setTrackCandates(const std::vector< EVENT::TrackerHitVec >&);

        std::vector< EVENT::TrackerHitVec > getTrackCandates() const;

    private:

        /** Tracks candidates hits 
         * This is actual data structure that has to be stored
         * in lcio collections.
         * 
         * Single element of this vector contains a set of hits
         * identified as those comming from a single track.
         */
        std::vector< EVENT::TrackerHitVec > _trackCandidates;

    };

}
#endif	/* EUTELLCOBJECTTRACKCANDIDATE_H */


