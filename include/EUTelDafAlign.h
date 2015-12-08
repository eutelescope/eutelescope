// Version: $Id$
//! Author Havard Gjersdal <haavagj@fys.uio.no>
/*
 *   This source code is part of the Eutelescope package of Marlin.
 */

#ifndef EUTELDAFALIGN_H
#define EUTELDAFALIGN_H

// built only if GEAR is available
#ifdef USE_GEAR

// marlin includes ".h"
#include "marlin/Processor.h"

// marlin util includes
#include "mille/Mille.h"

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
#include <iostream>
#include <fstream>

#include "EUTelDafBase.h"

namespace eutelescope {
  class EUTelDafAlign : EUTelDafBase{
  public:
    // Marlin processor interface funtions
    //! Returns a new instance of EUTelDafAlign
    virtual Processor * newProcessor() {
      return new EUTelDafAlign;
    }
    //! Default constructor
    EUTelDafAlign ();
    //! Called at the job beginning.
    virtual void dafInit ();
    //! Called every event
    virtual void dafEvent (LCEvent* event);
    //! Called after data processing.
    virtual void dafEnd();
    //! Set fitter specific params
    virtual void dafParams();

  protected:
    //params
    bool _runPede;
    std::string _pedeSteerfileName, _binaryFilename, _alignmentConstantLCIOFile, _alignmentConstantCollectionName;
    std::vector<int> _translate, _translateX, _translateY, _zRot, _scale, _scaleX, _scaleY;
    std::vector<float>_resXMin, _resXMax, _resYMin, _resYMax;
    //Variables
    Mille * _mille;
    std::map<int, std::pair<float, float> > _resX, _resY;
    std::vector<int> _dutMatches;
    //function
    //bool checkClusterRegion(lcio::TrackerHitImpl* hit);
    int checkDutResids(daffitter::TrackCandidate<float,4>& cnd);
    void addToMille(daffitter::TrackCandidate<float,4>& track);
    void runPede();
    void generatePedeSteeringFile();
    void steerLine(std::ofstream &steerFile, int label, int iden, std::vector<int> idens);
  };
  //! A global instance of the processor
  EUTelDafAlign gEUTelDafAlign;
}
#endif
#endif
