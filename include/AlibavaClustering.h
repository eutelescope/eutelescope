/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@cern.ch
 */

#ifndef ALIBAVACLUSTERING_H
#define ALIBAVACLUSTERING_H 1

// alibava includes ".h"
#include "AlibavaBaseProcessor.h"
#include "AlibavaEventImpl.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/TrackerDataImpl.h>

// ROOT includes <>
#include "TObject.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TH1.h>
#include <TF1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>

// system includes <>
#include <string>
#include <list>


using namespace std;

	//! This is used for the landau gaus function and fits
	/*! Details see:
	 *	http://root.cern.ch/root/html/tutorials/fit/langaus.C.html
	 */
	Double_t langaufun(Double_t *x, Double_t *par);
	TF1 *langaufit(TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF);
	Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM);

namespace alibava {
	
	//! Example Alibava processor for Marlin.
	
	class AlibavaClustering:public alibava::AlibavaBaseProcessor   {
		
	public:
		
		//! Returns a new instance of AlibavaClustering
		/*! This method returns an new instance of the this processor.  It
		 *  is called by Marlin execution framework and it shouldn't be
		 *  called/used by the final user.
		 *
		 *  @return a new AlibavaClustering.
		 */
		virtual Processor * newProcessor () {
			return new AlibavaClustering;
		}

		//! Default constructor
		AlibavaClustering ();

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
		 *  @param run the LCRu		std::string _commonmodeCollectionName;nHeader of the this current run
		 */
		virtual void processRunHeader (LCRunHeader * run);

		//! Called every event
		/*! Since the behavior of the PedestalNoise processor is different
		 *  if this is the first or one of the following loop, this method
		 *  is just calling
		 *  AlibavaClustering::firstLoop(LCEvent*) or
		 *  AlibavaClustering::otherLoop(LCEvent*)
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
		 *  histograms listed in AlibavaClustering::fillHistos()
		 *  there is also a common mode histo described here below:
		 *
		 *  \li commonModeHisto: 1D histogram to store the calculated
		 *  common mode value for each event. This histogram is booked and
		 *  filled only if the loop counter is greater-equal than 1,
		 *  because for _iLoop == 0 there is no common mode suppression.
		 *  This histo is not filled with the other because it needs to be
		 *  updated every event.
		 *
		 *  @see AlibavaClustering::fillHistos() for the todos
		 */
		void bookHistos();

		//! Fill histograms
		/*! This method is used to fill in histograms for each channel. 
		 *
		 */
		void fillHistos();

		//! Fill the clustersize histogram
		/*! This method is used to fill in the clustersize histogram.
		 *
		 */
		void fillClusterHisto(int clusize);

		//! Fill the Eta histogram
		/*! This method is used to fill in the eta histogram.
		 */
		void fillEtaHisto(float etaratio);

		//! Fill the second Eta histogram
		/*! This method is used to fill in the second eta histogram.
		 */
		void fillEtaHisto2(float etaratio);

		//! Fill the charge distribution histogram
		/*! This method is used to fill in the charge distribution histogram.
		 */
		void fillChargeDistHisto(float a, float b, float c, float d, float e, float f, float g);

		//! Fill the SNR histogram
		/*! This method is used to fill in the SNR histogram.
		 */
		void fillSNRHisto(float signal);

		//! Fill the signal histogram
		/*! This method is used to fill in the signal histogram.
		 */
		void fillSignalHisto(float signal);

		//! Fill the hitmap histogram
		/*! This method is used to fill in the hitmap histogram.
		 */
		void fillHitmapHisto(int ichan, int negclustersize, int posclustersize);

		//! Fill the seed histogram
		/*! This method is used to fill in the seed histogram.
		 */
		void fillSeedHisto(int ichan);

		//! Fill the seed charge histogram
		/*! This method is used to fill in the seed charge histogram.
		 */
		void fillSeedChargeHisto(float signal);
		
		//! Clustering is done here.
		void findSeedClusters(TrackerDataImpl * trkdata, LCCollectionVec * clusterCollection, LCCollectionVec * sparseClusterCollectionVec, AlibavaEventImpl * alibavaEvent, float othercoordinate);

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

		//! getter and setter for _clusterCollectionName
		void setClusterCollectionName(std::string clusterCollectionName);
		std::string getClusterCollectionName();
		
		//! Cluster collection name.
		/*! See _clusterCollectionName for the detailed description
		 */
		std::string _clusterCollectionName;		
		
		//! The ratio over noise to find seeds
		float _seedcut;

		//! The ratio over noise to find clusters
		float _clustercut;

		//! The charges around a seed are summed up here
		float _clustercharge[5];

		//! The final count of clusters
		int _clustercount;

		//! The name of the sparse cluster collection
		string _sparseclusterCollectionName;

		//! The polarity of the sensor
		int _polarity;

		//! The unsensitve axis of our strip sensor
		string _nonsensitiveaxis;
		
		//! The reading instance
		LCReader* lcReader;

		//! The flag if the file is open
		bool _telescopeopen;

		//! The reading function
		LCEvent *readTelescope ();

		//! The function to get the telescope data
		float getTelescope();

		//! The telescope plane we want to get the coordinate from, -1 to deactivate
		int _telescopePlane;

		//! The array storing the coordinate
		std::vector<float> _telescopecoordinate;

		//! The telescope collection name
		string _telescopeCollectionName;

		//! The telescope file we want to read
		string _telescopeFile;

		//! The maximum size of our run, this should be below the telescope run size, since EOFs are not always caught
		int _maxcount;

		//! The call to do a landau gaussian fit
		void dolandaugausfit(string tempHistoName);

	protected:

		IMPL::LCRunHeaderImpl* _runHeader;

	};

	//! A global instance of the processor
	AlibavaClustering gAlibavaClustering;

}

#endif
