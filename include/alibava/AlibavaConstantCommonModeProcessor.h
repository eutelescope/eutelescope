/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@cern.ch
 *
 *  modified by: Eda Yildirim eda.yildirim@cern.ch
 */

#ifndef ALIBAVACONSTANTCOMMONMODEPROCESSOR_H
#define ALIBAVACONSTANTCOMMONMODEPROCESSOR_H 1

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
	
	//! Common mode processor for Marlin.

	
	class AlibavaConstantCommonModeProcessor:public alibava::AlibavaBaseProcessor   {
		
	public:
		
		
		//! Returns a new instance of AlibavaConstantCommonModeProcessor
		/*! This method returns an new instance of the this processor.  It
		 *  is called by Marlin execution framework and it shouldn't be
		 *  called/used by the final user.
		 *
		 *  @return a new AlibavaConstantCommonModeProcessor.
		 */
		virtual Processor * newProcessor () {
			return new AlibavaConstantCommonModeProcessor;
		}
		
		//! Default constructor
		AlibavaConstantCommonModeProcessor ();
		
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
		 *  the EUTelCommonModeProcess::bookHistos() is called.
		 *
		 *  @param run the LCRunHeader of the this current run
		 */
		virtual void processRunHeader (LCRunHeader * run);
		
		//! Called every event
		/*! Since the behavior of the CommonMode processor is different
		 *  if this is the first or one of the following loop, this method
		 *  is just calling
		 *  AlibavaConstantCommonModeProcessor::firstLoop(LCEvent*) or
		 *  AlibavaConstantCommonModeProcessor::otherLoop(LCEvent*)
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
		 *  EUTelCommonModeProcess::_rootObjectMap so that they can be
		 *  recalled and filled from anywhere in the code.  Apart from the
		 *  histograms listed in AlibavaConstantCommonModeProcessor::fillHistos()
		 *  there is also a common mode histo described here below:
		 *
		 *  \li commonModeHisto: 1D histogram to store the calculated
		 *  common mode value for each event. This histogram is booked and
		 *  filled only if the loop counter is greater-equal than 1,
		 *  because for _iLoop == 0 there is no common mode suppression.
		 *  This histo is not filled with the other because it needs to be
		 *  updated every event.
		 *
		 *  @see AlibavaConstantCommonModeProcessor::fillHistos() for the todos
		 */
		void bookHistos();
		
		//! Fill histograms
		/*! This method is used to fill in histograms for each channel. 
		 *
		 */
		void fillHistos(TrackerDataImpl * trkdata, int event);
		

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
		
		//! Common Mode Error Calculation Iteration
		/*! The number of iteration that should be used in common mode calculation
		 */
		int _Niteration;
		
		//! Noise Deviation
		/*! The limit to the deviation of noise. The data exceeds this deviation will be considered as signal and not be included in common mode error calculation
		 */
		float _NoiseDeviation;
		
		// getter and setter for _commonmodeCollectionName
		void setCommonModeCollectionName(std::string CommonModeCollectionName);
		std::string getCommonModeCollectionName();
		
		// getter and setter for _commonmodeerrorCollectionName
		void setCommonModeErrorCollectionName(std::string CommonModeErrorCollectionName);
		std::string getCommonModeErrorCollectionName();

		// getter and setter for _commonmode
		void setCommonModeVec(EVENT::FloatVec common);
		EVENT::FloatVec getCommonModeVec();
		
		// getter and setter for _commonmodeerror
		void setCommonModeErrorVec(EVENT::FloatVec commonerror);
		EVENT::FloatVec getCommonModeErrorVec();

		
	protected:
										
		//! Name of the Common mode histogram 
		/*!
		 */
		std::string _commonmodeHistoName;
		
		//! Name of the Common mode error histogram 
		/*!
		 */
		std::string _commonmodeerrorHistoName;
		
		//! Calculates common mode values
		/*! Fills the histograms
		 */

		void calculateConstantCommonMode(TrackerDataImpl * trkdata);
		

		//! The function that returns name of the signal correction histo
		/*!
		 *  returns a name
		 */
		std::string getCommonCorrectionName();
		
		
		//! vector to store intermediate/final common mode value
		EVENT::FloatVec _commonmode;
		
		//! vector to store intermediate/final common mode error value
		EVENT::FloatVec _commonmodeerror;
		
		
	};
	
	//! A global instance of the processor
	AlibavaConstantCommonModeProcessor gAlibavaConstantCommonModeProcessor;
	
}

#endif
