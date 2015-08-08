/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */

#ifndef ALIBAVADATAHISTOGRAMMAKER_H
#define ALIBAVADATAHISTOGRAMMAKER_H 1

// alibava includes ".h"
#include "AlibavaBaseHistogramMaker.h"

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
	
	class AlibavaDataHistogramMaker : public alibava::AlibavaBaseHistogramMaker   {
		
	public:
		
		
		//! Returns a new instance of AlibavaDataHistogramMaker
		/*! This method returns an new instance of the this processor.  It
		 *  is called by Marlin execution framework and it shouldn't be
		 *  called/used by the final user.
		 *
		 *  @return a new AlibavaDataHistogramMaker.
		 */
		virtual Processor * newProcessor () {
			return new AlibavaDataHistogramMaker;
		}
		
		//! Default constructor
		AlibavaDataHistogramMaker ();
		
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
		 *  AlibavaDataHistogramMaker::firstLoop(LCEvent*) or
		 *  AlibavaDataHistogramMaker::otherLoop(LCEvent*)
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
		 *  AlibavaBaseProcessor::_rootObjectMap so that they can be
		 *  recalled and filled from anywhere in the code.  Apart from the
		 *  histograms listed in AlibavaDataHistogramMaker::fillHistos()
		 *  there is also a common mode histo described here below:
		 *
		 *  @see AlibavaDataHistogramMaker::fillHistos() for the todos
		 */
		virtual void bookHistos();
	
		// Books histogram for the event
		virtual void bookEventHisto(int eventnum);

		//! Fill histograms
		/*! This method is used to fill in histograms for each event.
		 *
		 */
		virtual void fillHistos(TrackerDataImpl * trkdata);
		

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
		
		
		// Fills the list of histogram names that has to be created
		virtual void fillListOfHistos();

    	///////////
		// Event //
		///////////
		
		// plots histogram for the event
		void fillEventHisto(int eventnum, TrackerDataImpl * trkdata);
				
		////////////
		// Others //
		////////////
		
		void fillOtherHistos(TrackerDataImpl * trkdata, float tdctime, float temperature);

		void createRulesToChangeXMLValues();

		
	protected:

    	////////////
		// Others //
		////////////

		// Signal
		//! Name of the Signal histogram
		std::string _signalHistoName;
		//! Name of the Signal vs TDCTime histogram
		std::string _signalVsTimeHistoName;
		//! Name of the Signal vs Temperature histogram
		std::string _signalVsTempHistoName;

		// SNR
		//! Name of the SNR histogram
		std::string _snrHistoName;
		//! Name of the SNR vs TDCTime histogram
		std::string _snrVsTimeHistoName;
		//! Name of the SNR vs Temperature histogram
		std::string _snrVsTempHistoName;

		// Time
		//! Name of the TDCTime histogram
		std::string _timeHistoName;
		//! Name of the TDCTime vs EventNum histogram
		std::string _timeVsEventNumHistoName;

		//Temperature
		//! Name of the Temperature histogram
		std::string _temperatureHistoName;
		//! Name of the Temperature vs EventNum histogram
		std::string _temperatureVsEventNumHistoName;

	};
	
	//! A global instance of the processor
	AlibavaDataHistogramMaker gAlibavaDataHistogramMaker;
	
}

#endif
