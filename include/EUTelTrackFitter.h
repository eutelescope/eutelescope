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

namespace eutelescope {

    class EUTelTrackFitter {
    public:
        EUTelTrackFitter();
 //       EUTelTrackFitter(const EUTelTrackFitter& orig);

        EUTelTrackFitter( std::string name );

        virtual ~EUTelTrackFitter();
        
        inline void SetName( std::string name) { this->_name = name; }
        inline std::string GetName() const { return _name; }

        virtual void SetTrackCandidates( std::map< int, EVENT::TrackerHitVec >& );

        virtual void FitTracks();
    protected:
        std::string _name;

    };

}
#endif	/* EUTELTRACKFITTER_H */

