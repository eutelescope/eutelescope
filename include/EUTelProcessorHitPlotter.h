//Writen by Alexander Morton <alexander.morton@desy.de>
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

//If header not already defined then define now. This is to stop multiple declarations of header files 
#ifndef EUTELPROCESSORHITPLOTTER_H
#define EUTELPROCESSORHITPLOTTER_H

// built only if GEAR is available
#ifdef USE_GEAR

// eutelescope includes ".h"
#include "EUTelUtility.h"
#include "EUTelEventImpl.h"

// marlin includes ".h"
#include "marlin/Processor.h"



namespace eutelescope {

  class EUTelProcessorHitPlotter : public marlin::Processor {

private:
    DISALLOW_COPY_AND_ASSIGN(EUTelProcessorHitPlotter) //This makes ensures that no other with this name can be created

  	public:

    // Returns a new instance of EUTelProcessorHitPlotter
    virtual Processor * newProcessor() {
      return new EUTelProcessorHitPlotter;
    }

		EUTelProcessorHitPlotter(); //Default constructor 
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
		std::string _hitCollectionNameInput;
    

	};//close class declaration

  //! A global instance of the processor
 	EUTelProcessorHitPlotter gEUTelProcessorHitPlotter;


}//close eutelescope namespace scope




#endif //Gear check ifndef
#endif //EUTelProcessorHitPlotter ifndef
