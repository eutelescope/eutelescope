//! Author Havard Gjersdal <haavagj@fys.uio.no>
// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 */

#ifndef EUTELDAFMATERIAL_H
#define EUTELDAFMATERIAL_H

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
#include "estmat.h"
namespace eutelescope {
  class EUTelDafMaterial : EUTelDafBase{
  private:
    double _angleX, _angleY;
    EstMat _matest;

    std::vector<float>_resXMin, _resXMax, _resYMin, _resYMax;
    std::map<int, std::pair<float, float> > _resX, _resY;
    int checkDutResids(daffitter::TrackCandidate<float, 4>& cnd);

    //Mat res
    std::vector<int> _radLengthIndex, _resXIndex, _resYIndex;    
    //Alignment
    std::vector<int> _shiftXIndex, _shiftYIndex, _scaleXIndex, _scaleYIndex, _zRotIndex, _zPosIndex; 
    
  public:
    // Marlin processor interface funtions
    //! Returns a new instance of EUTelDafMaterial
    virtual Processor * newProcessor() {
      return new EUTelDafMaterial;
    }
    //! Default constructor
    EUTelDafMaterial ();
    //! Called at the job beginning.
    virtual void dafInit ();
    //! Called every event
    virtual void dafEvent (LCEvent* event);
    //! Called after data processing.
    virtual void dafEnd();
    //! Set fitter specific params
    virtual void dafParams();
  };
  //! A global instance of the processor
  EUTelDafMaterial gEUTelDafMaterial;
}
#endif
#endif
