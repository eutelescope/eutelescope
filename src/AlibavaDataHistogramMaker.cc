/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


// alibava includes ".h"
#include "AlibavaDataHistogramMaker.h"
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


AlibavaDataHistogramMaker::AlibavaDataHistogramMaker () :
AlibavaBaseHistogramMaker("AlibavaDataHistogramMaker"),
// Signal
_signalHistoName("hSignal"),
_signalVsTimeHistoName("hSignal_vs_Time"),
_signalVsTempHistoName("hSignal_vs_Temperature"),
// SNR
_snrHistoName("hSNR"),
_snrVsTimeHistoName("hSNR_vs_Time"),
_snrVsTempHistoName("hSNR_vs_Temperature"),
// Time
_timeHistoName("hTDCTime"),
_timeVsEventNumHistoName("hTDCTime_vs_EventNum"),
// Temperature
_temperatureHistoName("hEventTemperatures"),
_temperatureVsEventNumHistoName("hEventTemperature_vs_EventNum")
{
	
	// modify processor description
	_description =
	"AlibavaDataHistogramMaker takes TrackerData of Alibava data and produces histograms ";
	
	
	// first of register the input collection
	registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
									  "Input raw data collection name",
									  _inputCollectionName, string("rawdata") );
	
	registerProcessorParameter ("HistoXMLFile",
										 "The path of XML file where the histograms are defined",
										 _histoXMLFileName , string("AlibavaHistoList.xml"));

	registerProcessorParameter ("TopTagInXMLFile",
										 "The top tag in HistoXMLFile",
										 _topTag , string("AlibavaHistoList"));

	registerProcessorParameter ("TagToProcess",
										 "The tag in TopTagInXMLFile. This processor will only consider the histogram definitions inside this tag. This tag should be inside <TopTagInXMLFile> ... <TopTagInXMLFile/>",
										 _tagToProcess , string("myAlibavaDataHistogramMaker"));

	// optional parameters
	
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
	

}


void AlibavaDataHistogramMaker::init () {
	streamlog_out ( MESSAGE4 ) << "Running init" << endl;


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

	// this method is called only once even when the rewind is active
	// usually a good idea to
	printParameters ();

	createRulesToChangeXMLValues();
	
}

void AlibavaDataHistogramMaker::processRunHeader (LCRunHeader * rdr) {
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
	
	// set pedestal and noise values (id it is defined)
	setPedestals();
	
	if (_plotPedestalAndNoise) plotPedestalAndNoise();
		
	// if you want
	bookHistos();
	
	// set number of skipped events to zero (defined in AlibavaBaseProcessor)
	_numberOfSkippedEvents = 0;
	

}

void AlibavaDataHistogramMaker::processEvent (LCEvent * anEvent) {
	TH1D * histo;
	TProfile * profile;
	
	if ( anEvent->getEventNumber() % 1000 == 0 )
		streamlog_out ( DEBUG1 ) << "Looping events "<<anEvent->getEventNumber() << endl;
	
	AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
	int eventnum = alibavaEvent->getEventNumber();

	if (_skipMaskedEvents && (alibavaEvent->isEventMasked()) ) {
		_numberOfSkippedEvents++;

		// Fill number of masked events histogram
		histo = dynamic_cast<TH1D*> (_rootObjectMap[_maskedEventsHistoName]);
		histo->Fill(eventnum);

		return;
	}

	float tdctime = alibavaEvent->getEventTime();
	float temperature = alibavaEvent->getEventTemp();
	
	//////////////////////////////////////////////
	// fill the histograms common for the event //

	// TDC time
	histo = dynamic_cast<TH1D*> (_rootObjectMap[_timeHistoName]);
	histo->Fill(tdctime);

	// TDC time vs EventNum
	profile = dynamic_cast<TProfile*>(_rootObjectMap[_timeVsEventNumHistoName]);
	profile->Fill(eventnum, tdctime);

	// Temperature
	histo = dynamic_cast<TH1D*> (_rootObjectMap[_temperatureHistoName]);
	histo->Fill(temperature);
	
	// Temperature vs EventNum
	profile = dynamic_cast<TProfile*>(_rootObjectMap[_temperatureVsEventNumHistoName]);
	profile->Fill(eventnum, temperature);
	
	// Calibration charge values
	histo = dynamic_cast<TH1D*> (_rootObjectMap[_calChargeHistoName]);
	histo->Fill(alibavaEvent->getCalCharge());

	// Calibration delay values
	histo = dynamic_cast<TH1D*> (_rootObjectMap[_delayHistoName]);
	histo->Fill(alibavaEvent->getCalDelay());
	
	bool plotThisEvent = false;
	if(isEventToBePlotted(eventnum)){
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
			if(plotThisEvent)
				fillEventHisto(eventnum,trkdata);
			
			// for everything else
			fillOtherHistos(trkdata, tdctime, temperature);
			
		}
		
	} catch ( lcio::DataNotAvailableException ) {
		// do nothing again
		streamlog_out( ERROR5 ) << "Collection ("<<getInputCollectionName()<<") not found! " << endl;
	}
	
}


void AlibavaDataHistogramMaker::bookHistos(){
	// events
	if (_plotEvents.size() != 0 || _plotXPercentOfEvents != 0 ){
		// create the directory for event histograms
		AIDAProcessor::tree(this)->cd(this->name());
		AIDAProcessor::tree(this)->mkdir(_eventHistoDir.c_str());
	}
			
	// create histograms defined in HistoXMLFile
	processHistoXMLFile();
	
	streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << endl;
}



void AlibavaDataHistogramMaker::fillListOfHistos(){
	// Checks if all the histograms needed by this processor is defined in _histoXMLFileName
	// Unfortunately histogram names are hard coded, so we are checking if all histo names exists in the _rootObjectMap
	
	// We need these histograms
	
	//////////////////////
	// One per each chip
	// Signal
	addToHistoCheckList_PerChip(_signalHistoName);
	addToHistoCheckList_PerChip(_signalVsTimeHistoName);
	addToHistoCheckList_PerChip(_signalVsTempHistoName);
	// SNR
	addToHistoCheckList_PerChip(_snrHistoName);
	addToHistoCheckList_PerChip(_snrVsTimeHistoName);
	addToHistoCheckList_PerChip(_snrVsTempHistoName);
	
	//////////////////////
	// One for all chips
	//	Time
	addToHistoCheckList(_timeHistoName);
	addToHistoCheckList(_timeVsEventNumHistoName);
	// Temperature
	addToHistoCheckList(_temperatureHistoName);
	addToHistoCheckList(_temperatureVsEventNumHistoName);
	// Others
	addToHistoCheckList(_calChargeHistoName);
	addToHistoCheckList(_delayHistoName);
	addToHistoCheckList(_maskedEventsHistoName);
	
	// here it checks
	checkListOfHistosCreatedByXMLFile();
}

void AlibavaDataHistogramMaker::check (LCEvent * /* evt */ ) {
	// nothing to check here - could be used to fill check plots in reconstruction processor
}


void AlibavaDataHistogramMaker::end() {
	
	if (_numberOfSkippedEvents > 0)
		streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents<<" events skipped since they are masked" << endl;
	streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
	
}
void AlibavaDataHistogramMaker::fillHistos(TrackerDataImpl * /* trkdata */){
	// nothing to do here, see AlibavaDataHistogramMaker::fillOtherHistos(TrackerDataImpl * trkdata, float tdctime, float temperature)
}

void AlibavaDataHistogramMaker::fillOtherHistos(TrackerDataImpl * trkdata, float tdctime, float temperature){

	FloatVec datavec;
	datavec = trkdata->getChargeValues();
	int ichip = getChipNum(trkdata);
	
	FloatVec noiseVec;
	if (isNoiseValid())
		noiseVec= getNoiseOfChip(ichip);
	
	TH1D * histoSignal = dynamic_cast<TH1D*> (_rootObjectMap[getHistoNameForChip(_signalHistoName,ichip)]);
	TH2D * histoSignalVsTime = dynamic_cast<TH2D*> (_rootObjectMap[getHistoNameForChip(_signalVsTimeHistoName,ichip)]);
	TH2D * histoSignalVsTemp = dynamic_cast<TH2D*> (_rootObjectMap[getHistoNameForChip(_signalVsTempHistoName,ichip)]);


	TH1D * histoSNR = dynamic_cast<TH1D*> (_rootObjectMap[getHistoNameForChip(_snrHistoName,ichip)]);
	TH2D * histoSNRVsTime = dynamic_cast<TH2D*> (_rootObjectMap[getHistoNameForChip(_snrVsTimeHistoName,ichip)]);
	TH2D * histoSNRVsTemp = dynamic_cast<TH2D*> (_rootObjectMap[getHistoNameForChip(_snrVsTempHistoName,ichip)]);
	
	for (int ichan=0; ichan<ALIBAVA::NOOFCHANNELS; ichan++) {
		// if channel is masked, do not fill histo
		if (isMasked(ichip,ichan)) continue;
		
		float data = _multiplySignalby*datavec[ichan];
		float noise = noiseVec[ichan];
		
		histoSignal->Fill(data);
		histoSignalVsTime->Fill(tdctime,data);
		histoSignalVsTemp->Fill(temperature,data);
		
		if (isNoiseValid() && noise != 0){
			float snr = data/noise;
			histoSNR->Fill(snr);
			histoSNRVsTime->Fill(tdctime,snr);
			histoSNRVsTemp->Fill(temperature,snr);
		}
		
		
	}
	

}

///////////
// Event //
///////////

void AlibavaDataHistogramMaker::bookEventHisto(int eventnum){
	
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

void AlibavaDataHistogramMaker::fillEventHisto(int eventnum, TrackerDataImpl * trkdata){
		
//	streamlog_out (DEBUG1) << "Plotting Event "<< eventnum<<endl;
	FloatVec datavec;
	datavec = trkdata->getChargeValues();
	int ichip = getChipNum(trkdata);
	
	TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[getEventHistoName(eventnum,ichip)]);
	
	for (int ichan=0; ichan<ALIBAVA::NOOFCHANNELS; ichan++) {
		// if channel is masked, do not fill histo
		if (isMasked(ichip,ichan)) continue;
		else histo->SetBinContent(ichan+1,datavec[ichan]);
	}
}

void AlibavaDataHistogramMaker::createRulesToChangeXMLValues(){
	
	// replace maximum event number
	changeXMLVariable(_timeVsEventNumHistoName, "xMax", float(_totalNumberOfEvents));
	changeXMLVariable(_temperatureVsEventNumHistoName, "xMax", float(_totalNumberOfEvents));
	changeXMLVariable(_maskedEventsHistoName, "xMax", float(_totalNumberOfEvents));

	// if signal multiplied by some number (other than 1)
	if (_multiplySignalby != 1.0){
		string signalMultipliedby;
		signalMultipliedby = to_string(_multiplySignalby);
		signalMultipliedby = string(" (")+signalMultipliedby+string(") ");

		addToXMLTitle(_signalHistoName, "labelX", "left", signalMultipliedby);
		addToXMLTitle(_snrHistoName, "labelX", "left", signalMultipliedby);

		addToXMLTitle(_signalVsTimeHistoName, "labelY", "left", signalMultipliedby);
		addToXMLTitle(_snrVsTimeHistoName, "labelY", "left", signalMultipliedby);
		addToXMLTitle(_signalVsTempHistoName, "labelY", "left", signalMultipliedby);
		addToXMLTitle(_snrVsTempHistoName, "labelY", "left", signalMultipliedby);
	}

}
