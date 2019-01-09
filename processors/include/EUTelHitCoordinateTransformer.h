/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELHITCOORDINATETRANSFORMER_H
#define EUTELHITCOORDINATETRANSFORMER_H

// built only if GEAR is available
#ifdef USE_GEAR

// eutelescope includes ".h"
#include "EUTelEventImpl.h"
#include "EUTelUtility.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesLayerLayout.h>
#include <gear/SiPlanesParameters.h>

// lcio includes <.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCRunHeader.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#endif
#include <IMPL/LCCollectionVec.h>

// system includes <>
#include <string>
#include <vector>

namespace eutelescope {

  class EUTelHitCoordinateTransformer : public marlin::Processor {

  private:
    DISALLOW_COPY_AND_ASSIGN(EUTelHitCoordinateTransformer)

  public:
    // Returns a new instance of EUTelHitCoordinateTransformer
    virtual Processor *newProcessor() {
      return new EUTelHitCoordinateTransformer;
    }

	// Default constructor
    EUTelHitCoordinateTransformer();
    
    // Called only at the begining of a job
    virtual void init();

    // Called every run
    virtual void processRunHeader(LCRunHeader*){};

    // Called every event.
    virtual void processEvent(LCEvent *event);

    // Called at the end of the job
    virtual void end();

  private:
    // Collection names
    std::string _hitCollectionNameInput;
    std::string _hitCollectionNameOutput;
    
    //parameter
    bool _undoAlignment;
  };

  //! A global instance of the processor
  EUTelHitCoordinateTransformer gEUTelHitCoordinateTransformer;

}
#endif
#endif
