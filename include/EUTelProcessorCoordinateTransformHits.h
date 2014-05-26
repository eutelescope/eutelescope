//Writen by Alexander Morton <alexander.morton@desy.de>
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

//If header not already defined then define now
#ifndef EUTELCOORDINATETRANSFORMHITS_H
#define EUTELCOORDINATETRANSFORMHITS_H


// built only if GEAR is available
#ifdef USE_GEAR

// eutelescope includes ".h"
#include "EUTelUtility.h"
#include "EUTelEventImpl.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCCollection.h>


// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#endif

#include <IMPL/LCCollectionVec.h>

// system includes <>
#include <string>
#include <vector>
#include <map>
#include <set>

// lcio includes <.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>


namespace eutelescope {

  class EUTelProcessorCoordinateTransformHits : public marlin::Processor {

  	private:
    DISALLOW_COPY_AND_ASSIGN(EUTelProcessorCoordinateTransformHits) //This makes ensures that no other with this name can be created

  	public:

    // Returns a new instance of EUTelLocaltoGlobalHitMaker
    virtual Processor * newProcessor() {
      return new EUTelProcessorCoordinateTransformHits;
    }

		EUTelProcessorCoordinateTransformHits(); //Default constructor 
		//Called only at the begining of a job
  	virtual void init ();

		//Called every run so here we check the geometry in the gear file and run header in teh lcio file is the same
  	virtual void processRunHeader (LCRunHeader * run);
	
		//Called every event.
  	virtual void processEvent (LCEvent * event);
		//This runs every event
		virtual void check(LCEvent *event);
		//Called at the end of the job
  	virtual void end();

		protected:

		private:
    
		//Only names wit _(name) come from the steering file.
		//Collection name
		std::string _hitCollectionNameInput;
		std::string _hitCollectionNameOutput;

	};//close class declaration

  //! A global instance of the processor
 	EUTelProcessorCoordinateTransformHits gEUTelLocaltoGlobalHitMaker;


}//close eutelescope namespace scope


#endif //Gear check ifndef
#endif //EUTelLocaltoGlobalHitMaker ifndef
