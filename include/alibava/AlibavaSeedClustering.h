/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */

#ifndef ALIBAVASEEDCLUSTERING_H
#define ALIBAVASEEDCLUSTERING_H 1

// alibava includes ".h"
#include "AlibavaBaseProcessor.h"
#include "AlibavaCluster.h"


// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/TrackerDataImpl.h>

// ROOT includes <>
#include "TObject.h"

// system includes <>
#include <string>
#include <list>


namespace alibava {
	
	//! Example Alibava processor for Marlin.

	
	class AlibavaSeedClustering:public alibava::AlibavaBaseProcessor   {
		
	public:
		
		
		//! Returns a new instance of AlibavaSeedClustering
		/*! This method returns an new instance of the this processor.  It
		 *  is called by Marlin execution framework and it shouldn't be
		 *  called/used by the final user.
		 *
		 *  @return a new AlibavaSeedClustering.
		 */
		virtual Processor * newProcessor () {
			return new AlibavaSeedClustering;
		}
		
		//! Default constructor
		AlibavaSeedClustering ();
		
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
		 *  AlibavaSeedClustering::firstLoop(LCEvent*) or
		 *  AlibavaSeedClustering::otherLoop(LCEvent*)
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
		 *  histograms listed in AlibavaSeedClustering::fillHistos()
		 *  there is also a common mode histo described here below:
		 *
		 *  @see AlibavaSeedClustering::fillHistos() for the todos
		 */
		void bookHistos();
		
		//! Fill histograms
		/*! This method is used to fill in histograms for each channel. 
		 *
		 */
		void fillHistos(AlibavaCluster anAlibavaCluster);
		

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
		
		//////////////////////////
		// Processor Parameters //
		//////////////////////////
		
		// SNR Cuts
		// The signal/noise ratio that channels have to pass to be considered as seed channel
		float _seedCut;
		
		// The signal/noise ratio that neigbour channels have to pass to be added to the cluster
		float _neighCut;
		
		
		// Sensitive Axis
		// The sensitive axis of the strip sensor(s) according to telescope. It has to set to either "X" ot "Y" any other value will be disregarded and sensitive axis will assumed to be "X" 
		int _sensitiveAxisX;

      // Signal Polarity
		// Polarity of the signal. Set this parameter to -1 for negative signals, any other value will be disregarded and the signal will be assumed to be positive
		int _signalPolarity;
		
		string getHistoNameForChip(string histoName, int ichip);

	protected:
		
		// Finds clusters in AlibavaCluster format and fills the histograms
		// Then calls processCluster for each cluster
		vector<AlibavaCluster> findClusters(TrackerDataImpl * trkdata);

		// to calculate Eta
		float calculateEta(TrackerDataImpl * trkdata, int seedChan);
		
	//	void convertAlibavaCluster(AlibavaCluster alibavaCluster, LCCollectionVec * clusterColVec, LCCollectionVec * sparseClusterColVec);
		
		/////////////////////
		// Histogram Names //
		/////////////////////

		// Eta histogram name for ClusterSize > 1
		string _etaHistoName;
		
		// Cluster Size
		string _clusterSizeHistoName;
		
		//
		bool _isSensitiveAxisX;
		
	
	};
	
	//! A global instance of the processor
	AlibavaSeedClustering gAlibavaSeedClustering;
	
}

#endif
