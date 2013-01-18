// Version: $Id$
//! Author Havard Gjersdal <haavagj@fys.uio.no>
/*
 *   This source code is part of the Eutelescope package of Marlin.
 */

#ifndef EUTELDAF_H
#define EUTELDAF_H

// built only if GEAR is available
#ifdef USE_GEAR

// eutelescope includes
#include "trackersystem.h"

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
#include "millehelper.h"

namespace eutelescope {
  class EUTelDafDutAlign : public marlin::Processor, EUTelDafBase{
  protected:
    std::vector<float> _dutResXMin;
    std::vector<float> _dutResXMax;
    std::vector<float> _dutResYMin;
    std::vector<float> _dutResYMax;
    std::map<size_t, std::vector<float> > _dutResids;
    int _nDutMatches;

    int checkDutResids(daffitter::TrackCandidate* cnd);
    void excludeDuts();
    
    std::string _milleOut, _steerName, _pedeExecutable;
    millehelper::MilleHelper _milleHelper;

  public:
    // Marlin processor interface funtions
    //! Returns a new instance of EUTelDafDutAlign
    virtual Processor * newProcessor() {
      return new EUTelDafDutAlign;
    }
    //! Default constructor
    EUTelDafDutAlign ();
    //! Called at the job beginning.
    virtual void init ();
    //! Called for every run.
    virtual void processRunHeader (LCRunHeader * run);
    //! Called every event
    virtual void processEvent (LCEvent * evt);
    //! Called after data processing.
    virtual void end();
  };
  //! A global instance of the processor
  EUTelDafDutAlign gEUTelDafDutAlign;
}
#endif
#endif
