/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


// alibava includes ".h"
#include "ExampleAlibavaHistogramMaker.h"
#include "AlibavaBaseHistogramMaker.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"
#include "AlibavaPedNoiCalIOManager.h"


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


ExampleAlibavaHistogramMaker::ExampleAlibavaHistogramMaker () :
AlibavaBaseHistogramMaker("ExampleAlibavaHistogramMaker"),
// List of Histogram names, initialized here. As an example we put only 2
_someHistoName("hSomeHisto"),
_someOtherHistoName("hSomeOtherHisto")
{
	
	// modify processor description
	_description =
	"ExampleAlibavaHistogramMaker takes some type of Alibava data and produces histograms ";
	
	
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
										 _tagToProcess , string("myExampleAlibavaHistogramMaker"));

	// optional parameters, these parameters are defined in AlibavaBaseHistogramMaker, but they are optional
	registerOptionalParameter ("PedestalInputFile",
										 "The filename where the pedestal and noise values stored",
										 _pedestalFile , string(ALIBAVA::NOTSET));

	registerOptionalParameter ("CalibrationInputFile",
										 "The filename where the calibration values stored",
										 _calibrationFile , string(ALIBAVA::NOTSET));

	registerOptionalParameter ("PedestalCollectionName",
										"Pedestal collection name, better not to change",
										_pedestalCollectionName, string (ALIBAVA::NOTSET));
	
	registerOptionalParameter ("NoiseCollectionName",
										"Noise collection name, better not to change",
										_noiseCollectionName, string (ALIBAVA::NOTSET));
	
	registerOptionalParameter ("ChargeCalibrationCollectionName",
										 "Charge calibration collection name, better not to change",
										 _chargeCalCollectionName, string (ALIBAVA::NOTSET));

	registerOptionalParameter ("PlotPedestalAndNoise",
										"Choose if pedestal and noise should be plotted. If you want to plot only noise or pedestal set this variable true and only set the noise or pedeatal collection name you want to be plotted.",
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


void ExampleAlibavaHistogramMaker::init () {
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

void ExampleAlibavaHistogramMaker::createRulesToChangeXMLValues(){
	// here we define which variables in HistoXMLFile we want to change or modify
	// You can only change or modify these variables:
	          // xBin, xMin, xMax, yBin, yMin, yMax, title, labelX and labelY
	          // others cannot be changed!
	
	
	/*! If you want to replace one of these variables: xBin, xMin, xMax, yBin, yMin, yMax
	 *  We will use the function
	 *      		void changeXMLVariable(string histoName, string variableName, float newValue);
	 *          defined in AlibavaBaseHistogramMaker
	 *          Note that variableNames are case sensitive!
	 */
	
	 //  For example usually it is good idea to replace maximum event number for example
	 //  Lets say _someOtherHistoName histograms X axis should be total number of events
	changeXMLVariable(_someOtherHistoName, "xMax", float(_totalNumberOfEvents));
	// And _someOtherHistoName histograms Y axis should be 1000 since it is the max value it can get
	changeXMLVariable(_someHistoName, "yMax", 1000.0);
	// It is really bad idea, but for some reason if you want to hard code the binning of X axis, here how you can do. Note that float number will be changed to integer in AlibavaBaseHistogramMaker::updateXMLVariables() method but here we should give it as float
	changeXMLVariable(_someOtherHistoName, "xBin", float(100));
	
	//HERE
	/*! If you want to replace or modify one of these variables: title, labelX and labelY
	 *  We will use the function
	 *      		void addToXMLTitle(string histoName, string titleName, string whichSide, string stringToAdd);
	 *          defined in AlibavaBaseHistogramMaker
	 *          Note that titleName and whichSide are case sensitive!
	 */

	// For example if MultiplySignalBy is set, it is a good idea to add it to the label of the "signal" axis
	if (_multiplySignalby != 1.0){
		// For this first create the string you want add
		string signalMultipliedby;
		signalMultipliedby = to_string(_multiplySignalby);
		signalMultipliedby = string(" (")+signalMultipliedby+string(") ");
		
		// Say that you have the signal in X axis of _someHistoName.
		// Assuming that "Signal" is written in HistoXMLFile as this histograms labelY, to get "(_multiplySignalby) Signal" you should add signalMultipliedby string to the left
		addToXMLTitle(_someHistoName, "labelY", "left", signalMultipliedby);
		// Again it is usually bad idea but if you want to replace the X label here how you can do
		addToXMLTitle(_someOtherHistoName, "labelX", "all", string("Event Number"));
		// And here how you add string to right of title
		addToXMLTitle(_someHistoName, "title", "right", string("Some thing"));
	}
	
}


void ExampleAlibavaHistogramMaker::processRunHeader (LCRunHeader * rdr) {
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

void ExampleAlibavaHistogramMaker::fillListOfHistos(){
	// Checks if all the histograms needed by this processor is defined in _histoXMLFileName
	// Unfortunately histogram names are hard coded, so we are checking if all histo names exists in the _rootObjectMap
	
	// We need these histograms
	
	//////////////////////
	// One per each chip
	// lets say _someHistoName should be created for each chip
	addToHistoCheckList_PerChip(_someHistoName);
	
	//////////////////////
	// One for all chips
	// and we need one _someOtherHistoName for all event
	addToHistoCheckList(_someOtherHistoName);

	// here it checks
	checkListOfHistosCreatedByXMLFile();
}

void ExampleAlibavaHistogramMaker::bookHistos(){
	// create histograms defined in HistoXMLFile
	processHistoXMLFile();
	
	// If you set _plotEvents or _plotXPercentOfEvents
	/*
	 // events
	 if (_plotEvents.size() != 0 || _plotXPercentOfEvents != 0 ){
	 // create the directory for event histograms
	 AIDAProcessor::tree(this)->cd(this->name());
	 AIDAProcessor::tree(this)->mkdir(_eventHistoDir.c_str());
	 }
	 */
	
	streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << endl;
}



void ExampleAlibavaHistogramMaker::processEvent (LCEvent * anEvent) {
	
	if ( anEvent->getEventNumber() % 1000 == 0 )
		streamlog_out ( DEBUG1 ) << "Looping events "<<anEvent->getEventNumber() << endl;
	
	AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
	int eventnum = alibavaEvent->getEventNumber();

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
	
	//////////////////////////////////////////////
	// fill the histograms common for the event //
	
	// Remember _someOtherHistoName common for all event
	TH1D * histo;
	histo = dynamic_cast<TH1D*> (_rootObjectMap[_someOtherHistoName]);
	// Stupid but just fill with eventnum
	histo->Fill(eventnum);
	
	
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
				// here you fill the event histogram for this chip
				// Dont forget to define what fillEventHisto method should do
				fillEventHisto(eventnum,trkdata);
			}
			
			// for everything else
			fillHistos(trkdata);
			
		}
		
	} catch ( lcio::DataNotAvailableException ) {
		// do nothing again
		streamlog_out( ERROR5 ) << "Collection ("<<getInputCollectionName()<<") not found! " << endl;
	}
	
}


void ExampleAlibavaHistogramMaker::fillHistos(TrackerDataImpl * trkdata ){
	// you fill the histograms here
		
	FloatVec datavec;
	datavec = trkdata->getChargeValues();
	int ichip = getChipNum(trkdata);
	
	// Remember we had _someHistoName for each chip, they are booked in AlibavaBaseHistogramMaker::processHistoXMLFile()
	TH2D * histo = dynamic_cast<TH2D*> (_rootObjectMap[getHistoNameForChip(_someHistoName,ichip)]);
	// OK it is again stupid but just as an example, fill x and y with signals in trkdata
	for (int ichan=0; ichan<int( datavec.size() ); ichan++) {
		// if channel is masked, do not fill histo
		if (isMasked(ichip,ichan)) continue;
		
		float data = _multiplySignalby*datavec[ichan];
		histo->Fill(data,data);
				
	}
}

void ExampleAlibavaHistogramMaker::check (LCEvent * /* evt */ ) {
	// nothing to check here
}


void ExampleAlibavaHistogramMaker::end() {
	
	if (_numberOfSkippedEvents > 0)
		streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents<<" events skipped since they are masked" << endl;
	streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;

}

// You have to define what bookEventHisto will do
void ExampleAlibavaHistogramMaker::bookEventHisto(int eventnum){
	
	AIDAProcessor::tree(this)->cd(this->name());
	AIDAProcessor::tree(this)->cd(_eventHistoDir.c_str());
	
	EVENT::IntVec chipSelection = getChipSelection();
	
	for (unsigned int i=0; i<chipSelection.size(); i++) {
		int ichip=chipSelection[i];
		
		TH1D * eventHisto = new TH1D (getEventHistoName(eventnum,ichip).c_str(),"",ALIBAVA::NOOFCHANNELS, -0.5, ALIBAVA::NOOFCHANNELS-0.5);
		
		stringstream sp; //title string for histogram
		sp<< "Event "<< eventnum <<" (chip "<<ichip<<");Channel Number;Signal (ADCs)";
		eventHisto->SetTitle((sp.str()).c_str());
		
		_rootObjectMap.insert(make_pair(getEventHistoName(eventnum,ichip), eventHisto));
		
	} // end of loop over selected chips
	
	
}

// You have to define what fillEventHisto will do
void ExampleAlibavaHistogramMaker::fillEventHisto(int eventnum, TrackerDataImpl * trkdata){
		
//	streamlog_out (DEBUG1) << "Plotting Event "<< eventnum<<endl;
	FloatVec datavec;
	datavec = trkdata->getChargeValues();
	int ichip = getChipNum(trkdata);
	
	TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[getEventHistoName(eventnum,ichip)]);
	// say you want to fill it with what ever you have in trkdata
	// we assume that
	for (int ichan=0; ichan<int( datavec.size() ); ichan++) {
		// if channel is masked, do not fill histo
		if (isMasked(ichip,ichan)) continue;
		else histo->SetBinContent(ichan+1,datavec[ichan]);
	}
}
