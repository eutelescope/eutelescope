/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */

#ifndef ALIBAVAPEDESTALSUBTRACTION_H
#define ALIBAVAPEDESTALSUBTRACTION_H 1

// alibava includes ".h"
#include "AlibavaBaseProcessor.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/TrackerDataImpl.h>

// ROOT includes <>
#include "TObject.h"

// system includes <>
#include <string>
#include <list>


namespace alibava {
	
	//! Example Alibava processor for Marlin.

	
	class AlibavaPedestalSubtraction:public alibava::AlibavaBaseProcessor   {
		
	public:
		
		
		//! Returns a new instance of AlibavaPedestalSubtraction
		/*! This method returns an new instance of the this processor.  It
		 *  is called by Marlin execution framework and it shouldn't be
		 *  called/used by the final user.
		 *
		 *  @return a new AlibavaPedestalSubtraction.
		 */
		virtual Processor * newProcessor () {
			return new AlibavaPedestalSubtraction;
		}
		
		//! Default constructor
		AlibavaPedestalSubtraction ();
		
		//! Called at the job beginning.
		/*! This is executed only once in the whole execution. It prints
		 *  out the processor parameters and reset all needed data
		 *  members. In principle this could also be the right place to
		 *  book histograms, but since those are also grouped according to
		 *  the detector numbers we need to have this parameter available.
		 */
		virtual void init ();
		
		//! Called for every run.
		/*! At the beginning of every run the run header is read and
		 *  processed by this method. As a first thing, the input run
		 *  header is dynamically re-casted to a EUTelRunHeaderImpl and
		 *  then important things like the number of detectors and the
		 *  pixel detector boundaries are dumped from the file. After that
		 *  the EUTelPedestalNoiseProcess::bookHistos() is called.
		 *
		 *  @param run the LCRunHeader of the this current run
		 */
		virtual void processRunHeader (LCRunHeader * run);
		
		//! Called every event
		/*! Since the behavior of the PedestalNoise processor is different
		 *  if this is the first or one of the following loop, this method
		 *  is just calling
		 *  AlibavaPedestalSubtraction::firstLoop(LCEvent*) or
		 *  AlibavaPedestalSubtraction::otherLoop(LCEvent*)
		 *
		 *  @param evt the current LCEvent event as passed by the
		 *  ProcessMgr
		 */
		virtual void processEvent (LCEvent * evt);
		
		
		//! Check event method
		/*! Does Nothing
		 *
		 */
		virtual void check (LCEvent * evt);
		
		
		//! Book histograms
		/*! Nothing is done here
		 *
		 */
		void bookHistos();
		
		//! Fill histograms
		/*! Nothing is done here
		 *
		 */
		void fillHistos(TrackerDataImpl * trkdata);
		

		//! Called after data processing.
		/*! This method is called when the loop on events is finished. It
		 *  is checking whether the calculation is properly finished or
		 *  not.
		 *  A very common error occurs when the file finished without a
		 *  EORE or when the MaxRecordNumber was set to low to loop over
		 *  all the needed record. To check this is very easy because we
		 *  just have to crosscheck if _iLoop is equal to noOfCMIterations.
		 */
		virtual void end();
				
	protected:
		
	
		
	};
	
	//! A global instance of the processor
	AlibavaPedestalSubtraction gAlibavaPedestalSubtraction;
	
}

#endif
