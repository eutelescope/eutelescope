/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */

#ifndef ALIBAVACONSTANTCOMMONMODECUTPROCESSOR_H
#define ALIBAVACONSTANTCOMMONMODECUTPROCESSOR_H 1

// alibava includes ".h"
#include "AlibavaBaseProcessor.h"
#include "AlibavaEventImpl.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/EventModifier.h"

// lcio includes <.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/TrackerRawDataImpl.h>

// ROOT includes <>
#include "TObject.h"

// system includes <>
#include <string>
#include <list>


namespace alibava {
	
	//! Example Alibava processor for Marlin.

	
	class AlibavaConstantCommonModeCutProcessor:public alibava::AlibavaBaseProcessor  {
		
	public:
		
		
		//! Returns a new instance of AlibavaConstantCommonModeCutProcessor
		/*! This method returns an new instance of the this processor.  It
		 *  is called by Marlin execution framework and it shouldn't be
		 *  called/used by the final user.
		 *
		 *  @return a new AlibavaConstantCommonModeCutProcessor.
		 */
		virtual Processor * newProcessor () {
			return new AlibavaConstantCommonModeCutProcessor;
		}
		
		//! Default constructor
		AlibavaConstantCommonModeCutProcessor ();
		
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
		 *  AlibavaConstantCommonModeCutProcessor::firstLoop(LCEvent*) or
		 *  AlibavaConstantCommonModeCutProcessor::otherLoop(LCEvent*)
		 *
		 *  @param evt the current LCEvent event as passed by the
		 *  ProcessMgr
		 */
		virtual void processEvent (LCEvent * evt);
		
		
		//! Check event method
		/*! This method is called by the Marlin execution framework as
		 *  soon as the processEvent is over. It can be used to fill check
		 *  plots. For the time being there is nothing to check and do in
		 *  this slot.
		 *
		 *  @param evt The LCEvent event as passed by the ProcessMgr
		 *
		 */
		virtual void check (LCEvent * evt);
		
		
		//! Book histograms
		/*! This method is used to prepare the needed directory structure
		 *  within the current ITree folder and books all required
		 *  histograms. Histogram pointers are stored into
		 *  EUTelPedestalNoiseProcess::_rootObjectMap so that they can be
		 *  recalled and filled from anywhere in the code.  Apart from the
		 *  histograms listed in AlibavaConstantCommonModeCutProcessor::fillHistos()
		 *  there is also a common mode histo described here below:
		 *
		 *  \li commonModeHisto: 1D histogram to store the calculated
		 *  common mode value for each event. This histogram is booked and
		 *  filled only if the loop counter is greater-equal than 1,
		 *  because for _iLoop == 0 there is no common mode suppression.
		 *  This histo is not filled with the other because it needs to be
		 *  updated every event.
		 *
		 *  @see AlibavaConstantCommonModeCutProcessor::fillHistos() for the todos
		 */
		void bookHistos();
		
		//! Fill histograms
		/*! This method is used to fill in histograms for each channel. 
		 *
		 */
		void fillHistos(AlibavaEventImpl * anAlibavaEvent);
		

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
		
		
		// getter for _timecutmin
		float getCommonModeCutMin();
		
		// getter for _timecutmax
		float getCommonModeCutMax();
		
		
	protected:

		//! Common Mode Ccollection Name
		/*! The common mode collection name
		 */
		string _commonModeColName;

		//! CommonModeCutMin
		/*! The minimum common mode noise that is acceptable to use that Event
		 */
		float _commonModeCutMin;
		
		//! CommonModeCutMax
		/*! The maximum common mode noise that is acceptable to use that Event
		 */
		float _commonModeCutMax;
		
		//! Mask If Any Chips Common Mode Is Not In Range
		/*! The maximum common mode noise that is acceptable to use that Event
		 */
		bool _maskIfAnyChipsCommonModeIsNotInRange;

		// Number of masked events by this processor
		int _numberOfMaskedEvents;

		// Number of masked events by this processor
		int _totalNumberOfMaskedEvents;
		
		// Total number of events in this data
		int _totalNumberOfEvents;

		// The name of the histogram which contains event number of masked events
		string _hMaskedEventsNumberName;
		
		
		
	};
	
	//! A global instance of the processor
	AlibavaConstantCommonModeCutProcessor gAlibavaConstantCommonModeCutProcessor;
	
}

#endif
