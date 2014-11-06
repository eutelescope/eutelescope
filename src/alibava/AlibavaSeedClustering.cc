/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */

// alibava includes ".h"
#include "AlibavaSeedClustering.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"
#include "AlibavaPedNoiCalIOManager.h"
#include "AlibavaCluster.h"


// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <lcio.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCEventImpl.h>


// ROOT includes ".h"
#include "TH1D.h"

// system includes <>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <memory>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;


AlibavaSeedClustering::AlibavaSeedClustering () :
AlibavaBaseProcessor("AlibavaSeedClustering"),
_seedCut(3),
_neighCut(2),
_sensitiveAxisX(1),
_signalPolarity(-1),
_etaHistoName("hEta"),
_clusterSizeHistoName("hClusterSize"),
_isSensitiveAxisX(true)
{
	
	// modify processor description
	_description =
	"AlibavaSeedClustering finds clusters using seed and neighbour cuts ";
	
	
	// first of register the input /output collection
	registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
									 "Input collection name, it should be pedestal subtracted",
									 _inputCollectionName, string("recodata") );
	
	registerOutputCollection (LCIO::TRACKERDATA, "OutputCollectionName",
									  "Output data collection name",
									  _outputCollectionName, string("alibava_clusters") );
	
	
	
	// if needed one can change these to optional parameters
	
	registerProcessorParameter ("NoiseInputFile",
										 "The filename where the pedestal and noise values stored",
										 _pedestalFile , string("pedestal.slcio"));
	
	registerProcessorParameter ("NoiseCollectionName",
										 "Noise collection name, better not to change",
										 _noiseCollectionName, string ("noise"));
	
	registerProcessorParameter ("SeedSNRCut",
										 "The signal/noise ratio that channels have to pass to be considered as seed channel",
										 _seedCut, float (3));
	
	registerProcessorParameter ("NeighbourSNRCut",
										 "The signal/noise ratio that neigbour channels have to pass to be added to the cluster",
										 _neighCut, float (2));
	
	registerProcessorParameter ("IsSensitiveAxisX",
										 "The default sensitive axis of the strip sensor(s) according to telescope is X. If sensitive axis is Y then set this parameter to zero (0). Any other value will be disregarded and sensitive axis will assumed to be \"X\" ",
										 _sensitiveAxisX, int(1));
	
	
	registerProcessorParameter ("SignalPolarity",
										 "Polarity of the signal. Set this parameter to -1 for negative signals, any other value will be disregarded and the signal will be assumed to be positive ",
										 _signalPolarity, int (-1));
	
	
}


void AlibavaSeedClustering::init () {
	streamlog_out ( MESSAGE4 ) << "Running init" << endl;
	
	// this method is called only once even when the rewind is active
	
	/* To set of channels to be used
	 ex.The format should be like $ChipNumber:StartChannel-EndChannel$
	 ex. $0:5-20$ $0:30-100$ $1:50-70$
	 means from chip 0 channels between 5-20 and 30-100, from chip 1 channels between 50-70 will be used (all numbers included). the rest will be masked and not used
	 Note that the numbers should be in ascending order and there should be no space between two $ character
	 */
	if (Global::parameters->isParameterSet(ALIBAVA::CHANNELSTOBEUSED))
		Global::parameters->getStringVals(ALIBAVA::CHANNELSTOBEUSED,_channelsToBeUsed);
	else {
		streamlog_out ( MESSAGE4 ) << "The Global Parameter "<< ALIBAVA::CHANNELSTOBEUSED <<" is not set!" << endl;
	}
	
	
	/* To choose if processor should skip masked events
	 ex. Set the value to 0 for false, to 1 for true
	 */
	if (Global::parameters->isParameterSet(ALIBAVA::SKIPMASKEDEVENTS))
		_skipMaskedEvents = bool ( Global::parameters->getIntVal(ALIBAVA::SKIPMASKEDEVENTS) );
	else {
		streamlog_out ( MESSAGE4 ) << "The Global Parameter "<< ALIBAVA::SKIPMASKEDEVENTS <<" is not set! Masked events will be used!" << endl;
	}
	
	// check signal Polarity
	if (_signalPolarity != -1)
		_signalPolarity = 1;
	
	// check sensitive axis
	if (_sensitiveAxisX == 0)
		_isSensitiveAxisX = false;
	else
		_isSensitiveAxisX = true;
	
	
	// usually a good idea to
	printParameters ();
	
}

void AlibavaSeedClustering::processRunHeader (LCRunHeader * rdr) {
	streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;
	
	// Add processor name to the runheader
	auto_ptr<AlibavaRunHeaderImpl> arunHeader ( new AlibavaRunHeaderImpl(rdr)) ;
	arunHeader->addProcessor(type());
	
	// get and set selected chips
	setChipSelection(arunHeader->getChipSelection());
	
	// set channels to be used (if it is defined)
	setChannelsToBeUsed();
	
	// set pedestal and noise values
	setPedestals();
	
	// if you want
	bookHistos();
	
	// set number of skipped events to zero (defined in AlibavaBaseProcessor)
	_numberOfSkippedEvents = 0;
	
}


void AlibavaSeedClustering::processEvent (LCEvent * anEvent) {
	
	if ( anEvent->getEventNumber() % 1000 == 0 )
		streamlog_out ( MESSAGE4 ) << "Looping events "<<anEvent->getEventNumber() << endl;
	
	AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
	
	if (_skipMaskedEvents && (alibavaEvent->isEventMasked()) ) {
		_numberOfSkippedEvents++;
		return;
	}
	
	// As an input collection we get pedestal subtracted alibava data
	// This collection contains AlibavaEvents
	LCCollectionVec * inputColVec;
	
	// The Alibava Cluster collection
	LCCollectionVec * clusterColVec = new LCCollectionVec(LCIO::TRACKERDATA);
	// cell id encode for AlibavaCluster
	CellIDEncoder<TrackerDataImpl> clusterIDEncoder(ALIBAVA::ALIBAVACLUSTER_ENCODE,clusterColVec);
	
	unsigned int noOfChip;
	try
	{
		inputColVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( getInputCollectionName() ) ) ;
		noOfChip = inputColVec->getNumberOfElements();
		
		for ( size_t i = 0; i < noOfChip; ++i ){
			// get your data from the collection and do what ever you want
			
			TrackerDataImpl * trkdata = dynamic_cast< TrackerDataImpl * > ( inputColVec->getElementAt( i ) ) ;
			vector<AlibavaCluster> clusters = findClusters(trkdata);
			
			// loop over clusters
			for (unsigned int icluster=0; icluster<clusters.size(); icluster++) {
				AlibavaCluster acluster = clusters[icluster];
				// create a TrackerDataImpl for each cluster
				TrackerDataImpl * alibavaCluster = new TrackerDataImpl();
				acluster.createTrackerData(alibavaCluster);
				
				// now store chip number, seed channel, cluster ID, cluster size, sensitive axis and signal polarity in CellIDEncode
				
				// set cluster ID
				clusterIDEncoder[ALIBAVA::ALIBAVACLUSTER_ENCODE_CLUSTERID] = acluster.getClusterID();
				// set sensitive axis
				if (acluster.getIsSensitiveAxisX())
					clusterIDEncoder[ALIBAVA::ALIBAVACLUSTER_ENCODE_ISSENSITIVEAXISX] = 1;
				else
					clusterIDEncoder[ALIBAVA::ALIBAVACLUSTER_ENCODE_ISSENSITIVEAXISX] = 0;
				// set signal polarity
				if (acluster.getSignalPolarity() == -1)
					clusterIDEncoder[ALIBAVA::ALIBAVACLUSTER_ENCODE_ISSIGNALNEGATIVE] = 1;
				else
					clusterIDEncoder[ALIBAVA::ALIBAVACLUSTER_ENCODE_ISSIGNALNEGATIVE] = 0;
				// set chip number
				clusterIDEncoder[ALIBAVA::ALIBAVACLUSTER_ENCODE_CHIPNUM] = acluster.getChipNum();
				// set cluster size
				clusterIDEncoder[ALIBAVA::ALIBAVACLUSTER_ENCODE_CLUSTERSIZE] = acluster.getClusterSize();
				// set seed channel number
				clusterIDEncoder[ALIBAVA::ALIBAVACLUSTER_ENCODE_SEED] = acluster.getSeedChanNum();
				clusterIDEncoder.setCellID(alibavaCluster);
				
				clusterColVec->push_back(alibavaCluster);
				
			} // end of loop over clusters
			
		} // end of loop ever detectors
		
		alibavaEvent->addCollection(clusterColVec, getOutputCollectionName());
		
	} catch ( lcio::DataNotAvailableException ) {
		// do nothing again
		streamlog_out( ERROR5 ) << "Collection ("<<getInputCollectionName()<<") not found! " << endl;
	}
	
}

vector<AlibavaCluster> AlibavaSeedClustering::findClusters(TrackerDataImpl * trkdata){
	
	// first get chip number
	int chipnum = getChipNum(trkdata);
	
	// then get the data vector
	FloatVec dataVec;
	dataVec = trkdata->getChargeValues();
	
	// we will need noise vector too
	FloatVec noiseVec = getNoiseOfChip(chipnum);
	
	// then check which channels we can add to a cluster
	// obviously not the ones masked
	vector<bool> channel_can_be_used;
	for (int ichan=0; ichan<int(dataVec.size()); ichan++)
		channel_can_be_used.push_back(!isMasked(chipnum,ichan));
	
	
	// Now mask channels that cannot pass NeighbourSNRCut
	// And add the channel numbers as seed candidate if it pass SeedSNRCut
	vector<int> seedCandidates;
	seedCandidates.clear();
	for (int ichan=0; ichan<int(channel_can_be_used.size()); ichan++) {
		// if this channel is already masked, continue
		if (!channel_can_be_used[ichan]) continue;
		// now calculate snr = signal/noise
		float snr = (_signalPolarity * dataVec[ichan])/noiseVec[ichan];
		// mask channels that cannot pass NeighbourSNRCut
		if (snr < _neighCut)
			channel_can_be_used[ichan]=false;
		else if (snr > _seedCut) {
			seedCandidates.push_back(ichan);
		}
	}
	
	// sort seed channels according to their SNR, highest comes first!
	for (unsigned int i=0; i+1<seedCandidates.size(); i++) {
		int ichan = seedCandidates[i];
		float isnr = (_signalPolarity*dataVec[ichan])/noiseVec[ichan];
		int next_chan = seedCandidates[i+1];
		float next_snr = (_signalPolarity*dataVec[next_chan])/noiseVec[next_chan];
		if(isnr < next_snr){ // swap
			seedCandidates[i] = next_chan;
			seedCandidates[i+1] = ichan;
			//start from beginning
			i=0;
		}
	}
	
	// now form clusters starting from the seed channel that has highest SNR
	int clusterID = 0;
	// form clusters and store them in a vector
	vector<AlibavaCluster> clusterVector;
	for (unsigned int iseed=0; iseed<seedCandidates.size(); iseed++) {
		// if this seed channel used in another cluster, skip it
		int seedChan = seedCandidates[iseed];
		if (!channel_can_be_used[seedChan]) continue;
		
		AlibavaCluster acluster;
		acluster.setChipNum(chipnum);
		acluster.setSeedChanNum(seedChan);
		acluster.setEta( calculateEta(trkdata,seedChan) );
		acluster.setIsSensitiveAxisX(_isSensitiveAxisX);
		acluster.setSignalPolarity(_signalPolarity);
		acluster.setClusterID(clusterID);
		clusterID++;
		
		// add seed channel to the cluster
		acluster.add(seedChan, dataVec[seedChan]);
		// now mask seed channel so no other cluster can use it!
		channel_can_be_used[seedChan]=false;
		
		// We will check if one of the neighbours not bonded. If it is not this is not a good cluster we will not use this
		bool thereIsNonBondedChan = false;
				
		// add channels on the left
		int ichan = seedChan-1;
		bool continueSearch = true;
		while (continueSearch) {
			// first if ichan < 0
			if (ichan < 0) {
				// stop search
				continueSearch = false;
				// this also means that left neighbour is not bonded
				thereIsNonBondedChan = true;
				// then exit while loop
				break;
			}
			
			// if chan masked
			if (isMasked(chipnum,ichan)==true) {
				// this means that it is not bonded.
				thereIsNonBondedChan = true;
				// then we are sure that there is no other channel on the left
				break;
				// We will still form the cluster but we will not use this cluster.
				// the members on the right that should be in this cluster, should not be used by some other cluster.
				// so to mark them not to be used we will continue search
			}
			
			// now check if this channel can be added to the cluster
			if ( channel_can_be_used[ichan] ){
				// if yes, add it
				acluster.add(ichan, dataVec[ichan]);
				// and mask it so that it will not be added to any other cluster
				channel_can_be_used[ichan]=false;
			}
			else{
				// if you cannot use this channel there is no other channel on the left
				break;
			}
			
			// go to the next channel on left
			ichan--;
		}
		
		continueSearch = true;
		ichan = seedChan+1;
		// add channels on right
		while (continueSearch) {
			// first if ichan < 0
			if (ichan >= int(channel_can_be_used.size()) ) {
				// stop search
				continueSearch = false;
				// this also means that right neighbour is not bonded
				thereIsNonBondedChan = true;
				// then exit while loop
				break;
			}
			
			// if chan masked
			if (isMasked(chipnum,ichan)==true) {
				// this means that it is not bonded.
				thereIsNonBondedChan = true;
				// then we are sure that there is no other channel on the left
				break;
			}
			
			// now check if this channel can be added to the cluster
			if ( channel_can_be_used[ichan] ){
				// if yes, add it
				acluster.add(ichan, dataVec[ichan]);
				// and mask it so that it will not be added to any other cluster
				channel_can_be_used[ichan]=false;
			}
			else{
				// if you cannot use this channel there is no other channel on the right
				break;
			}
			
			// go to the next channel on right
			ichan++;
		}
		
		// now if there is no neighbour not bonded
		if(thereIsNonBondedChan == false){
			// fill the histograms and add them to the cluster vector
			fillHistos(acluster);
			clusterVector.push_back(acluster);
		}
		
	}
	return clusterVector;
}

float AlibavaSeedClustering::calculateEta(TrackerDataImpl *trkdata, int seedChan){
	
	// first get chip number
	int chipnum = getChipNum(trkdata);
	// then get the data vector
	FloatVec dataVec;
	dataVec = trkdata->getChargeValues();

	// we will multiply all signal values by _signalPolarity to work on positive signal always
	float seedSignal = _signalPolarity * dataVec.at(seedChan);

	// now we make an unrealistic signal
	float unrealisticSignal = -10000; // you cannot get this signal from alibava
	
	int leftChan = seedChan - 1;
	float leftSignal = unrealisticSignal;
	// check if the channel on the left is masked
	if ( leftChan >= 0 && isMasked(chipnum,leftChan)==false ) {
		leftSignal = _signalPolarity * dataVec.at(leftChan);
	}
	
	int rightChan = seedChan+1;
	float rightSignal = unrealisticSignal;
	// check if the channel on the right is masked
	if ( rightChan < int( dataVec.size() ) && isMasked(chipnum, rightChan) == false ) {
		rightSignal = _signalPolarity * dataVec.at(rightChan);
	}
	
	// now compare right and left channel to see which one is higher

	// if both right anf left channel is masked. Simply return -1
	// this case should not be saved by clustering algorithm anyways
	if (rightSignal == unrealisticSignal && leftSignal == unrealisticSignal ) {
		streamlog_out (DEBUG1) << "Both neighbours are masked!"<<endl;
		return -1;
	}
	
	float eta = -1;
	// compare left and right signals
	// here both signal has to be positive (if not noise)
	// if one of the channel is masked, since it will be set to unrealisticSignal other signal will always be higher then the masked one.

	if ( leftSignal > rightSignal) {
		// then seed channel is on the right
		eta = leftSignal / ( leftSignal + seedSignal );
	}
	else {
		// seed channel is on the left
		eta = seedSignal / (seedSignal + rightSignal);
	}
	
	// Eta calculation: chargeOnLeftChannel / (chargeOnLeftChannel + chargeOnRightChannel)

	return eta;
	
}


void AlibavaSeedClustering::check (LCEvent * /* evt */ ) {
	// nothing to check here - could be used to fill check plots in reconstruction processor
}


void AlibavaSeedClustering::end() {
	
	if (_numberOfSkippedEvents > 0)
		streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents<<" events skipped since they are masked" << endl;
	streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
	
}

void AlibavaSeedClustering::fillHistos(AlibavaCluster anAlibavaCluster){
	
	int ichip = anAlibavaCluster.getChipNum();
	
	float eta = anAlibavaCluster.getEta();
	TH1D * hEta = dynamic_cast<TH1D*> (_rootObjectMap[getHistoNameForChip(_etaHistoName,ichip)]);
	hEta->Fill(eta);
	
	int clusterSize = anAlibavaCluster.getClusterSize();
	TH1D * hClusterSize = dynamic_cast<TH1D*> (_rootObjectMap[getHistoNameForChip(_clusterSizeHistoName,ichip)]);
	hClusterSize->Fill( clusterSize );
	
}

void AlibavaSeedClustering::bookHistos(){
	
	AIDAProcessor::tree(this)->cd(this->name());
	
	EVENT::IntVec chipSelection = getChipSelection();
	string histoName;
	string title;
	for ( unsigned int i=0; i<chipSelection.size(); i++) {
		int ichip=chipSelection[i];
		
		// Clustersize histogram
		histoName = getHistoNameForChip(_clusterSizeHistoName,ichip);
		TH1D * hClusterSize = new TH1D (histoName.c_str(),"",10, 0, 10);
		title = string("Cluster Size (chip ") +to_string(ichip)+string(");Number of Entries;Cluster Size");
		hClusterSize->SetTitle(title.c_str());
		_rootObjectMap.insert(make_pair(histoName, hClusterSize));
		
		// Eta histogram ClusterSize > 1
		histoName = getHistoNameForChip(_etaHistoName,ichip);
		TH1D * hEta = new TH1D (histoName.c_str(),"",100, -0.1, 1.1);
		title = string("Eta distribution ClusterSize > 1 (chip ") +to_string(ichip)+string(");Number of Entries;Eta");
		hEta->SetTitle(title.c_str());
		_rootObjectMap.insert(make_pair(histoName, hEta));
		
	} // end of loop over selected chips
	
	streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << endl;
}

string AlibavaSeedClustering::getHistoNameForChip(string histoName, int ichip){
	stringstream s;
	s<< histoName <<"_chip" << ichip;
	return s.str();
}




