/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@cern.ch
 */

#ifndef ALIBAVACOMMONMODESUBTRACTION_H
#define ALIBAVACOMMONMODESUBTRACTION_H 1

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

	
	class AlibavaCommonModeSubtraction:public alibava::AlibavaBaseProcessor   {
		
	public:
		
		
		//! Returns a new instance of AlibavaCommonModeSubtraction
		/*! This method returns an new instance of the this processor.  It
		 *  is called by Marlin execution framework and it shouldn't be
		 *  called/used by the final user.
		 *
		 *  @return a new AlibavaCommonModeSubtraction.
		 */
		virtual Processor * newProcessor () {
			return new AlibavaCommonModeSubtraction;
		}
		
		//! Default constructor
		AlibavaCommonModeSubtraction ();
		
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
		 *  AlibavaCommonModeSubtraction::firstLoop(LCEvent*) or
		 *  AlibavaCommonModeSubtraction::otherLoop(LCEvent*)
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
		 *  histograms listed in AlibavaCommonModeSubtraction::fillHistos()
		 *  there is also a common mode histo described here below:
		 *
		 *  \li commonModeHisto: 1D histogram to store the calculated
		 *  common mode value for each event. This histogram is booked and
		 *  filled only if the loop counter is greater-equal than 1,
		 *  because for _iLoop == 0 there is no common mode suppression.
		 *  This histo is not filled with the other because it needs to be
		 *  updated every event.
		 *
		 *  @see AlibavaCommonModeSubtraction::fillHistos() for the todos
		 */
		void bookHistos();
		
		//! Fill histograms
		/*! This method is used to fill in histograms for each channel. 
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
		
		
		
		//! Common mode collection name.
		/*! See _commonmodeCollectionName for the detailed description
		 */
		std::string _commonmodeCollectionName;
		
		//! Common mode error collection name.
		/*! See _commonmodeCollectionName for the detailed description
		 */
		std::string _commonmodeerrorCollectionName;

		
		
	protected:

		//! The name of the histogram used to calculate common mode
		/*! For every channel a histogram will be created
		 *  and filled with the readings of that channel
		 *  
		 *  The actual name of the histogram will be 
		 *  _chanRawDataHistoName+chanNum
		 */
		std::string _chanDataHistoName;
		

		//! The function that returns name of the histogram
		/*!
		 *  returns
		 *  _chanRawDataHistoName+chanNum
		 */
		std::string getChanDataHistoName(int chipnum, int ichan);

		//! The function that returns name of the signal correction histo
		/*!
		 *  returns a name
		 */
		std::string getSignalCorrectionName();
	

	};
	
	//! A global instance of the processor
	AlibavaCommonModeSubtraction gAlibavaCommonModeSubtraction;
	
}

#endif
