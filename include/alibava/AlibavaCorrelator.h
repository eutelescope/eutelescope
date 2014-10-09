/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */

#ifndef ALIBAVACORRELATOR_H
#define ALIBAVACORRELATOR_H 1

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
	
	class AlibavaCorrelator : public alibava::AlibavaBaseHistogramMaker   {
		
	public:
		
		
		//! Returns a new instance of AlibavaCorrelator
		/*! This method returns an new instance of the this processor.  It
		 *  is called by Marlin execution framework and it shouldn't be
		 *  called/used by the final user.
		 *
		 *  @return a new AlibavaCorrelator.
		 */
		virtual Processor * newProcessor () {
			return new AlibavaCorrelator;
		}
		
		//! Default constructor
		AlibavaCorrelator ();
		
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
		 *  AlibavaCorrelator::firstLoop(LCEvent*) or
		 *  AlibavaCorrelator::otherLoop(LCEvent*)
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
		 *  histograms listed in AlibavaCorrelator::fillHistos()
		 *  there is also a common mode histo described here below:
		 *
		 *  @see AlibavaCorrelator::fillHistos() for the todos
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
		
		void createRulesToChangeXMLValues();

		
	protected:


		// Hit position X
		//! Name of the histogram of hit position X
		std::string _hHitPosX;

		// Hit position Y
		//! Name of the histogram of hit position Y
		std::string _hHitPosY;

		// Hit correlation in X
		//! Name of the histogram of hit correlation X
		std::string _hCorX;

		// Hit correlation in Y
		//! Name of the histogram of hit correlation Y
		std::string _hCorY;

		// Hit residual in X vs eventnumber
		//! Name of the histogram of hit redidual in X versus event number.
		/*! This plot shows the syncronization of events
		 */
		std::string _hSyncX;

		// Hit residual in Y vs eventnumber
		//! Name of the histogram of hit redidual in Y versus event number.
		/*! This plot shows the syncronization of events
		 */
		std::string _hSyncY;

		// number of detectors
		IntVec _detectorIDs;
		
		
		// just adding "_d" and detector number to the end of string
		std::string getHistoNameForDetector(std::string name, int detID);
		
		// just adding "_d" detID1 and "_d" detID2 to the end of string
		std::string getHistoNameForDetector(std::string name, int detID1, int detID2);

		// clones TH1D Histogram for each detector
		void createClones_hHitPos(string histoName);

		// clones Correlation Histogram for each detector combination
		void createClones_hCor(string histoName);

		// clones Synchronisation Histogram for each detector combination
		void createClones_hSync(string histoName);

		// checks if detID is in _detectorIDs list
		bool isInDetectorIDsList(int detID);

	};
	
	//! A global instance of the processor
	AlibavaCorrelator gAlibavaCorrelator;
	
}

#endif
