/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


// alibava includes ".h"
#include "ExampleAlibavaProcessor.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"
#include "AlibavaPedNoiCalIOManager.h"


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

// system includes <>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <memory>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;


ExampleAlibavaProcessor::ExampleAlibavaProcessor () :
AlibavaBaseProcessor("ExampleAlibavaProcessor")
{
	
	// modify processor description
	_description =
	"ExampleAlibavaProcessor does whatever it does :) ";
	
	
	// first of register the input collection
	registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
									  "Input raw data collection name",
									  _inputCollectionName, string("rawdata") );

	// if needed one can change these to optional parameters
	
	registerProcessorParameter ("PedestalInputFile",
										 "The filename where the pedestal and noise values stored",
										 _pedestalFile , string("pedestal.slcio"));

	registerProcessorParameter ("CalibrationInputFile",
										 "The filename where the calibration values stored",
										 _calibrationFile , string("calibration.slcio"));


	// now the optional parameters
	registerProcessorParameter ("PedestalCollectionName",
										"Pedestal collection name, better not to change",
										_pedestalCollectionName, string ("pedestal"));
	
	registerProcessorParameter ("NoiseCollectionName",
										"Noise collection name, better not to change",
										_noiseCollectionName, string ("noise"));
	
	registerProcessorParameter ("ChargeCalibrationCollectionName",
										 "Charge calibration collection name, better not to change",
										 _chargeCalCollectionName, string ("chargeCal"));

}


void ExampleAlibavaProcessor::init () {
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

	
}
void ExampleAlibavaProcessor::processRunHeader (LCRunHeader * rdr) {
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


void ExampleAlibavaProcessor::processEvent (LCEvent * anEvent) {
	
	if ( anEvent->getEventNumber() % 1000 == 0 )
		streamlog_out ( MESSAGE4 ) << "Looping events "<<anEvent->getEventNumber() << endl;
	
	AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);

	if (_skipMaskedEvents && (alibavaEvent->isEventMasked()) ) {
		_numberOfSkippedEvents++;
		return;
	}

	LCCollectionVec * collectionVec;
	unsigned int noOfChip;
	try
	{
		collectionVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( getInputCollectionName() ) ) ;
		noOfChip = collectionVec->getNumberOfElements();
		
		for ( size_t i = 0; i < noOfChip; ++i )
		{
			// get your data from the collection and do what ever you want
			/*
			 Example:
			 TrackerDataImpl * trkdata = dynamic_cast< TrackerDataImpl * > ( collectionVec->getElementAt( i ) ) ;
			 fillHistos(trkdata);
			 */
			
		}
		
	} catch ( lcio::DataNotAvailableException ) {
		// do nothing again
		streamlog_out( ERROR5 ) << "Collection ("<<getInputCollectionName()<<") not found! " << endl;
	}
	
}

void ExampleAlibavaProcessor::check (LCEvent * /* evt */ ) {
	// nothing to check here - could be used to fill check plots in reconstruction processor
}


void ExampleAlibavaProcessor::end() {

	if (_numberOfSkippedEvents > 0)
		streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents<<" events skipped since they are masked" << endl;
	streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
	
}

void ExampleAlibavaProcessor::fillHistos(TrackerDataImpl * /* trkdata */){
	// this is an example

	/*
	 int ichip = getChipNum(trkdata);
	 FloatVec datavec;
	 datavec = trkdata->getChargeValues();
	for (size_t ichan=0; ichan<datavec.size();ichan++) {
		if(isMasked(ichip, ichan)) continue;
		
		string tempHistoName = getChanDataHistoName(ichip,ichan);
		if ( TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[tempHistoName]) )
			histo->Fill(datavec[ichan]);
	}
	 */

}


void ExampleAlibavaProcessor::bookHistos(){
	// this is an example
	/*
	AIDAProcessor::tree(this)->cd(this->name());
	int nChips = getNumberOfChips();

	TH1D * pedestalHisto = new TH1D (_pedestalHistoName.c_str(),"",nChips*ALIBAVA::NOOFCHANNELS, -0.5, nChips*ALIBAVA::NOOFCHANNELS-0.5);
	_rootObjectMap.insert(make_pair(_pedestalHistoName, pedestalHisto));
	pedestalHisto->SetTitle("Pedestal;Channel Number;Pedestal (ADCs)");
	
	AIDAProcessor::tree(this)->mkdir(getInputCollectionName().c_str());
	AIDAProcessor::tree(this)->cd(getInputCollectionName().c_str());
	
	if (!(getNumberOfChips()>0)) {
		streamlog_out ( ERROR5 )<<"Number of chips in the collection ("<<getInputCollectionName()<<") is = "<< nChips<<endl;
	}
	
	 */
	streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << endl;
}


