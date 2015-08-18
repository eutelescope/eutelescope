// Version: $Id$
//! Author Havard Gjersdal <haavagj@fys.uio.no>
/*
 *   This source code is part of the Eutelescope package of Marlin.
 */

#ifndef EUTELDAF_H
#define EUTELDAF_H

// built only if GEAR is available
#ifdef USE_GEAR

// marlin includes ".h"
#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IHistogram1D.h>
#endif

// system includes <>
#include <string>
#include <vector>
#include <map>

#include "EUTelDafBase.h"

namespace eutelescope {
  class EUTelDafFitter : EUTelDafBase{
  public:
    // Marlin processor interface funtions
    //! Returns a new instance of EUTelDafFitter
    virtual Processor * newProcessor() {
      return new EUTelDafFitter;
    }
    //! Default constructor
    EUTelDafFitter ();
    //! Called at the job beginning.
    virtual void dafInit ();
    //! Called every event
    virtual void dafEvent (LCEvent* event);
    //! Called after data processing.
    virtual void dafEnd();
    //! Set fitter specific params
    virtual void dafParams();
    //! Calculate Z coordinate for a pair X:Y on a plane define by refhit collection and "plane" id 
    /* input - plane and pos[0] and pos[1]
     * output - pos[2] to be overwritten with the right value
     */
    virtual double getZfromRefHit(int plane,int sensorID, double *pos);
 
  protected:
    //! Output track collection name
    std::string _trackCollectionName;
    //! Output track collection
    LCCollectionVec     * _fittrackvec;
    LCCollectionVec     * _fitpointvec;
    void addToLCIO(daffitter::TrackCandidate<float,4>& track, LCCollectionVec *lcvec);
    //! LCIO switch
    bool _addToLCIO, _fitDuts;

  };
  //! A global instance of the processor
  EUTelDafFitter gEUTelDafFitter;
}
#endif
#endif
