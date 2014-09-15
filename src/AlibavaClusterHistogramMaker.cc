/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


// alibava includes ".h"
#include "AlibavaClusterHistogramMaker.h"
#include "AlibavaBaseHistogramMaker.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"
#include "AlibavaPedNoiCalIOManager.h"
#include "AlibavaCluster.h"


// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"
#include "marlin/tinyxml.h"

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
#include "TH2D.h"
#include "TProfile.h"

// system includes <>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <memory>
#include <cstdlib>
#include <ctime>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;


AlibavaClusterHistogramMaker::AlibavaClusterHistogramMaker () :
AlibavaBaseHistogramMaker("AlibavaClusterHistogramMaker"),
// List of Histogram names, initialized here.
_etaHistoName("hEta"),
_etaHistoNameCS2("hEtaCS2"),
_etaHistoNameCS3("hEtaCS3"),
_etaHistoNameCS4("hEtaCS4"),
_etaHistoNameCS5("hEtaCS5"),
_etaHistoNameCSgt5("hEtaCSgt5"),
_clusterSizeHistoName("hClusterSize"),
_etaVSCoG("hEta_vs_CoG"),
_etaVSClusterSize("hEta_vs_ClusterSize")
{
	
	// modify processor description
	_description =
	"AlibavaClusterHistogramMaker takes some type of Alibava data and produces histograms ";
	
	
	// first of register the input collection
	registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
									 "Input raw data collection name",
									 _inputCollectionName, string("rawdata") );
	// Details about HistoXMLFile
	registerProcessorParameter ("HistoXMLFile",
										 "The path of XML file where the histograms are defined",
										 _histoXMLFileName , string("AlibavaHistoList.xml"));
	
	registerProcessorParameter ("TopTagInXMLFile",
										 "The top tag in HistoXMLFile",
										 _topTag , string("AlibavaHistoList"));
	
	registerProcessorParameter ("TagToProcess",
										 "The tag in TopTagInXMLFile. This processor will only consider the histogram definitions inside this tag. This tag should be inside <TopTagInXMLFile> ... <TopTagInXMLFile/>",
										 _tagToProcess , string("myAlibavaClusterHistogramMaker"));
	
	// optional parameters, these parameters are defined in AlibavaBaseHistogramMaker, but they are optional
	registerOptionalParameter ("NoiseInputFile",
										"The filename where the pedestal and noise values stored",
										_pedestalFile , string(ALIBAVA::NOTSET));
	
	registerOptionalParameter ("CalibrationInputFile",
										"The filename where the calibration values stored",
										_calibrationFile , string(ALIBAVA::NOTSET));
	
	registerOptionalParameter ("NoiseCollectionName",
										"Noise collection name, better not to change",
										_noiseCollectionName, string (ALIBAVA::NOTSET));
	
	registerOptionalParameter ("ChargeCalibrationCollectionName",
										"Charge calibration collection name, better not to change",
										_chargeCalCollectionName, string (ALIBAVA::NOTSET));
	
	registerOptionalParameter ("PlotNoise",
										"Choose if noise should be plotted. If you want to plot only noise set this variable true and set the noise collection name.",
										_plotPedestalAndNoise, bool(false));
	
	registerOptionalParameter ("EventsToPlot",
										"Choose if pedestal and noise should be plotted. If you want to plot only noise or pedestal set this variable true and only set the noise or pedeatal collection name you want to be plotted.",
										_plotEvents, IntVec());
	
	registerOptionalParameter ("MultiplySignalBy",
										"In case this variable is set, all signals will be multipled by this value.",
										_multiplySignalby, float(1.0));
	
	registerOptionalParameter ("PlotSomePercentOfEvents",
										"In case this variable is set (say x), x percent of total events will be plotted randomly. The number should be between 0-100",
										_plotXPercentOfEvents, float(0.0));
	
	// add whatever you need for this histogram maker
	
}


void AlibavaClusterHistogramMaker::init () {
	streamlog_out ( MESSAGE4 ) << "Running init" << endl;
	
	// this method is called only once even when the rewind is active
	// usually a good idea to
	printParameters ();
	
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
	
	// if you want to change any variable defined in HistoXMLFile use this function
	// see example function below
	createRulesToChangeXMLValues();
}

void AlibavaClusterHistogramMaker::createRulesToChangeXMLValues(){
	// here we define which variables in HistoXMLFile we want to change or modify
	// You can only change or modify these variables:
	// xBin, xMin, xMax, yBin, yMin, yMax, title, labelX and labelY
	// others cannot be changed!
	
	
	changeXMLVariable(_maskedEventsHistoName, "xMax", float(_totalNumberOfEvents));
	
	
	// if MultiplySignalBy is set, it is a good idea to add it to the label of the "signal" axis
	if (_multiplySignalby != 1.0){
		// For this first create the string you want add
		string signalMultipliedby;
		signalMultipliedby = to_string(_multiplySignalby);
		signalMultipliedby = string(" (")+signalMultipliedby+string(") ");
		
		//	addToXMLTitle(_someHistoName, "labelY", "left", signalMultipliedby);
	}
	
}


void AlibavaClusterHistogramMaker::processRunHeader (LCRunHeader * rdr) {
	streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;
	
	// Add processor name to the runheader
	auto_ptr<AlibavaRunHeaderImpl> arunHeader ( new AlibavaRunHeaderImpl(rdr)) ;
	arunHeader->addProcessor(type());
	
	// set total number of events in this run
	_totalNumberOfEvents = arunHeader->getNoOfEvents();
	streamlog_out ( DEBUG1 ) << "N events "<<_totalNumberOfEvents << endl;
	
	// get and set selected chips
	setChipSelection(arunHeader->getChipSelection());
	
	// set channels to be used (if it is defined)
	setChannelsToBeUsed();
	
	// set pedestal and noise values (if it is defined)
	setPedestals();
	
	if (_plotPedestalAndNoise) plotPedestalAndNoise();
	
	// reads and creates the histograms defined in HistoXMLFile
	bookHistos();
	
	// set number of skipped events to zero (defined in AlibavaBaseProcessor)
	_numberOfSkippedEvents = 0;
	
	
}

void AlibavaClusterHistogramMaker::fillListOfHistos(){
	// Checks if all the histograms needed by this processor is defined in _histoXMLFileName
	// Unfortunately histogram names are hard coded, so we are checking if all histo names exists in the _rootObjectMap
	
	// We need these histograms
	
	//////////////////////
	// One per each chip
	addToHistoCheckList_PerChip(_etaHistoName);
	addToHistoCheckList_PerChip(_etaHistoNameCS2);
	addToHistoCheckList_PerChip(_etaHistoNameCS3);
	addToHistoCheckList_PerChip(_etaHistoNameCS4);
	addToHistoCheckList_PerChip(_etaHistoNameCS5);
	addToHistoCheckList_PerChip(_etaHistoNameCSgt5);
	addToHistoCheckList_PerChip(_clusterSizeHistoName);
	addToHistoCheckList_PerChip(_etaVSCoG);
	addToHistoCheckList_PerChip(_etaVSClusterSize);
	
	//////////////////////
	// One for all chips
	addToHistoCheckList(_maskedEventsHistoName);
	
	// here it checks
	checkListOfHistosCreatedByXMLFile();
}

void AlibavaClusterHistogramMaker::bookHistos(){
	// create histograms defined in HistoXMLFile
	processHistoXMLFile();
	
	
	// If you set _plotEvents or _plotXPercentOfEvents
	// events
	if (_plotEvents.size() != 0 || _plotXPercentOfEvents != 0 ){
		// create the directory for event histograms
		AIDAProcessor::tree(this)->cd(this->name());
		AIDAProcessor::tree(this)->mkdir(_eventHistoDir.c_str());
	}
	
	streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << endl;
}



void AlibavaClusterHistogramMaker::processEvent (LCEvent * anEvent) { // HERE look for it

	
	AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
	int eventnum = alibavaEvent->getEventNumber();

	bool printClusters = false;
	if (eventnum == 7830) {
		printClusters = true;
	}
	if ( eventnum % 1000 == 0 ){
		streamlog_out ( DEBUG1 ) << "Looping events "<<eventnum << endl;
		printClusters = true;
	}
	
	// if _skipMaskedEvents is set
	if (_skipMaskedEvents && (alibavaEvent->isEventMasked()) ) {
		_numberOfSkippedEvents++;
		
		// Fill number of masked events histogram
		TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[_maskedEventsHistoName]);
		histo->Fill(eventnum);
		
		return;
	}
	
	// if you set PlotSomePercentOfEvents or EventsToPlot paramaters
	bool plotThisEvent = false;
	if(isEventToBePlotted(eventnum)){
		// bookEventHisto will book a histogram for each chip
		bookEventHisto(eventnum);
		plotThisEvent = true;
	}
	
	
	/////////////////////////////
	// Now loop ever detectors //
	LCCollectionVec * collectionVec;
	unsigned int noOfChip;
	try{
		collectionVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( getInputCollectionName() ) ) ;
		noOfChip = collectionVec->getNumberOfElements();
		
		for ( size_t i = 0; i < noOfChip; ++i ){
			TrackerDataImpl * trkdata = dynamic_cast< TrackerDataImpl * > ( collectionVec->getElementAt( i ) ) ;
			
			// if this event selected to be plotted fill the histogram
			if(plotThisEvent){
				// here we fill the event histogram for all clusters in this event
				// how it will fill is defined in fillEventHisto
				fillEventHisto(eventnum,trkdata);
			}
			
			// for everything else
			fillHistos(trkdata);
			if (printClusters) {
				printAlibavaCluster(trkdata);
			}
			
		}
		
	} catch ( lcio::DataNotAvailableException ) {
		// do nothing again
		streamlog_out( ERROR5 ) << "Collection ("<<getInputCollectionName()<<") not found! " << endl;
	}
	
}

void AlibavaClusterHistogramMaker::printAlibavaCluster(TrackerDataImpl* trkdata){
	AlibavaCluster anAlibavaCluster(trkdata);
	anAlibavaCluster.print();
}

void AlibavaClusterHistogramMaker::fillHistos(TrackerDataImpl * trkdata ){
	
	// We will first convert his trkdata to AlibavaCluster
	// see AlibavaCluster::AlibavaCluster(TrackerDataImpl * trkdata)
	
	AlibavaCluster anAlibavaCluster(trkdata);
	int ichip = anAlibavaCluster.getChipNum();
	
	TH1D * histo;
	TProfile * profile;
	
	// Lets fill Cluster size histogram
	int clusterSize = anAlibavaCluster.getClusterSize();
	histo = dynamic_cast<TH1D*> (_rootObjectMap[getHistoNameForChip(_clusterSizeHistoName,ichip)]);
	histo->Fill( clusterSize );
	
	// Then fill eta histograms
	float eta = anAlibavaCluster.getEta();
	histo= dynamic_cast<TH1D*> (_rootObjectMap[getHistoNameForChip(_etaHistoName,ichip)]);
	histo->Fill(eta);

	// center of gravity
	float CoG = anAlibavaCluster.getCenterOfGravity();
	profile = dynamic_cast<TProfile*> (_rootObjectMap[getHistoNameForChip(_etaVSCoG,ichip)]);
	profile->Fill(CoG, eta);

	// center of gravity
	profile = dynamic_cast<TProfile*> (_rootObjectMap[getHistoNameForChip(_etaVSClusterSize,ichip)]);
	profile->Fill(clusterSize, eta);
	
	// Fill histos depending on their Cluster size

	if (clusterSize == 2) {
		histo = dynamic_cast<TH1D*> (_rootObjectMap[getHistoNameForChip(_etaHistoNameCS2,ichip)]);
		histo->Fill(eta);
	}
	else if (clusterSize == 3){
		histo = dynamic_cast<TH1D*> (_rootObjectMap[getHistoNameForChip(_etaHistoNameCS3,ichip)]);
		histo->Fill(eta);
	}
	else if (clusterSize == 4){
		histo = dynamic_cast<TH1D*> (_rootObjectMap[getHistoNameForChip(_etaHistoNameCS4,ichip)]);
		histo->Fill(eta);
	}
	else if (clusterSize == 5){
		histo = dynamic_cast<TH1D*> (_rootObjectMap[getHistoNameForChip(_etaHistoNameCS5,ichip)]);
		histo->Fill(eta);
	}
	else if (clusterSize > 5){
		histo = dynamic_cast<TH1D*> (_rootObjectMap[getHistoNameForChip(_etaHistoNameCSgt5,ichip)]);
		histo->Fill(eta);
	}
	
}

void AlibavaClusterHistogramMaker::check (LCEvent * /* evt */ ) {
	// nothing to check here
}


void AlibavaClusterHistogramMaker::end() {
	
	if (_numberOfSkippedEvents > 0)
		streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents<<" events skipped since they are masked" << endl;
	streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
	
}

// You have to define what bookEventHisto will do
void AlibavaClusterHistogramMaker::bookEventHisto(int eventnum){
	
	AIDAProcessor::tree(this)->cd(this->name());
	AIDAProcessor::tree(this)->cd(_eventHistoDir.c_str());
	
	EVENT::IntVec chipSelection = getChipSelection();
	
	for (unsigned int i=0; i<chipSelection.size(); i++) {
		int ichip=chipSelection[i];
		
		string histoName = getEventHistoName(eventnum,ichip);
		// to be sure that we won't create more than one histogram per event per chip
		if (!doesRootObjectExists(histoName)){
			TH1D * eventHisto = new TH1D (histoName.c_str(),"",ALIBAVA::NOOFCHANNELS, -0.5, ALIBAVA::NOOFCHANNELS-0.5);
			
			stringstream sp; //title string for histogram
			sp<< "Event "<< eventnum <<" (chip "<<ichip<<");Channel Number;";
			// if signal multiplied by some number put it to the label
			if(_multiplySignalby != 1)
				sp <<" ("<< _multiplySignalby << ") ";
			// in any case here is labelY
			sp <<"Signal (ADCs)";
			
			eventHisto->SetTitle((sp.str()).c_str());
			
			_rootObjectMap.insert(make_pair(histoName, eventHisto));
		}
	} // end of loop over selected chips
	
}

// You have to define what fillEventHisto will do
void AlibavaClusterHistogramMaker::fillEventHisto(int eventnum, TrackerDataImpl * trkdata){
	//noiseHisto->SetBinContent(ichan+1,noiVec[ichan]);
	
	//	streamlog_out (DEBUG1) << "Plotting Event "<< eventnum<<endl;
	AlibavaCluster anAlibavaCluster(trkdata);
	int ichip = anAlibavaCluster.getChipNum();

	
	string histoName =getEventHistoName(eventnum,ichip);
	if (!doesRootObjectExists(histoName)) {
		streamlog_out(ERROR3)<<"Root object for Event "<<eventnum<<" doesn't exists!"<<endl;
		return;
	}
	TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[histoName]);
	
	// convert the trkdata to AlibavaCluster
	int Nmembers = anAlibavaCluster.getClusterSize();
	for (int imember=0; imember<Nmembers; imember++) {
		int channelNum = anAlibavaCluster.getChanNum(imember);
		float signal = _multiplySignalby * anAlibavaCluster.getSignal(imember);
		histo->SetBinContent(channelNum+1, signal);
	}
	
	
	
}
