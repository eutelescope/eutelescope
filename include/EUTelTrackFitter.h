/* 
 * File:   EUTelTrackFitter.h
 * Contact: denys.lontkovskyi@desy.de
 *
 * Created on January 25, 2013, 6:55 PM
 */

#ifndef EUTELTRACKFITTER_H
#define	EUTELTRACKFITTER_H

// eutelescope includes ".h"
#include <IMPL/TrackerHitImpl.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include "lcio.h"

// system includes <>
#include <string>
#include <map>

// EUTELESCOPE
#include "EUTelUtility.h"

namespace eutelescope {

    class EUTelTrackFitter {
      
    private:
        DISALLOW_COPY_AND_ASSIGN(EUTelTrackFitter)        // prevent users from making (default) copies of processors
      
    public:
        EUTelTrackFitter();
 
        explicit EUTelTrackFitter( std::string name );

        virtual ~EUTelTrackFitter();
        
        inline void SetName( std::string name) { this->_name = name; }
        inline std::string GetName() const { return _name; }

        virtual void SetTrackCandidates( const EVENT::TrackVec& );

        virtual void SetTrackCandidates( const std::vector<IMPL::TrackImpl*> & );

      // do some clean up of internal data structures
        virtual void Clear();


        virtual void SearchTrackCandidates(){};
        /** Prune track candidates
         *  supposed to be removing track candidates which have n% hits in common      */  
        virtual void PruneTrackCandidates(){};

        virtual void FitTracks();
        virtual void TrackCandidatesToGBLTrajectories();
        virtual void PerformFitGBLTrajectories();
        virtual void FitSingleTrackCandidate();

    protected:
        std::string _name;

    };

}
#endif	/* EUTELTRACKFITTER_H */

