/*
 * File:   EUTelTrackFitter.h
 * Contact: denys.lontkovskyi@desy.de
 *
 * Created on January 25, 2013, 6:55 PM
 */

#ifndef EUTELTRACKFITTER_H
#define EUTELTRACKFITTER_H

// eutelescope includes ".h"

// EVENT includes
#include <EVENT/LCCollection.h>
#include <IMPL/TrackerHitImpl.h>

// LCIO includes
#include "lcio.h"
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackerHitImpl.h>

// system includes <>
#include <map>
#include <string>

// EUTELESCOPE
#include "EUTelUtility.h"

namespace eutelescope {

  class EUTelTrackFitter {

  protected:
    std::string _name;

  private:
    DISALLOW_COPY_AND_ASSIGN(EUTelTrackFitter) // prevent users from making
                                               // (default) copies of processors

    // contains the fitted tracks, accessible through class methods
    IMPL::LCCollectionVec *_LCIO_fittrackvec;

    // contains the fitted hits, accessible through class methods
    IMPL::LCCollectionVec *_LCIO_fithitvec;

  public:
    EUTelTrackFitter();

    explicit EUTelTrackFitter(std::string name);

    virtual ~EUTelTrackFitter();

    inline void SetName(std::string name) { this->_name = name; }
    inline std::string GetName() const { return _name; }

    virtual void SetTrackCandidates(const EVENT::TrackVec &);

    virtual void SetTrackCandidates(std::vector<const IMPL::TrackImpl *> &);

    // do some clean up of internal data structures
    virtual void Clear();

    virtual void SearchTrackCandidates(){};
    /** Prune track candidates
     *  supposed to be removing track candidates which have n% hits in common */
    virtual void PruneTrackCandidates(){};

    virtual void FitTracks();
    virtual void TrackCandidatesToGBLTrajectories();
    virtual void PerformFitGBLTrajectories();
    virtual bool PerformMille();

    virtual void FitSingleTrackCandidate();

    /** */
    virtual IMPL::LCCollectionVec *getFitHitVec() const {
      return _LCIO_fithitvec;
    }

    /** */
    virtual IMPL::LCCollectionVec *getFitTrackVec() const {
      return _LCIO_fittrackvec;
    }
  };
}
#endif /* EUTELTRACKFITTER_H */
