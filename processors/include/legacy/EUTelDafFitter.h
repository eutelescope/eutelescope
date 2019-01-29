/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELDAFFITTER_H
#define EUTELDAFFITTER_H

// built only if GEAR is available
#ifdef USE_GEAR

// eutelescope includes
#include "EUTelDafBase.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCRunHeader.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IHistogram1D.h>
#endif

// system includes <>
#include <map>
#include <string>
#include <vector>

namespace eutelescope {
 
  class EUTelDafFitter : EUTelDafBase {
 
  public:
    //! Returns a new instance of EUTelDafFitter
    virtual Processor *newProcessor() { return new EUTelDafFitter; }
    //! Default constructor
    EUTelDafFitter();
    //! Called at the job beginning.
    virtual void dafInit();
    //! Called every event
    virtual void dafEvent(LCEvent *event);
    //! Called after data processing.
    virtual void dafEnd();
    //! Set fitter specific params
    virtual void dafParams();

  protected:
    //! Output track collection name
    std::string _trackCollectionName;
    //! Output track collection
    LCCollectionVec *_fittrackVec;
    LCCollectionVec *_fitpointVec;
    //! function for adding track to LCIO
    void addToLCIO(daffitter::TrackCandidate<float, 4> &track,
                   LCCollectionVec *lcvec);
    //! LCIO switch
    bool _addToLCIO;
    //! switch to include DUT in track fit
    bool _fitDuts;
  };
  //! A global instance of the processor
  EUTelDafFitter gEUTelDafFitter;
}
#endif
#endif
