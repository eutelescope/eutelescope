/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@cern.ch
 */


// alibava includes ".h"
#include "AlibavaCommonModeSubtraction.h"
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
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCEventImpl.h>


// ROOT includes ".h"
#include "TH1D.h"
#include "TF1.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TSystem.h"

// system includes <>
#include <string>
#include <iostream>
#include <memory>


using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;


AlibavaCommonModeSubtraction::AlibavaCommonModeSubtraction () :
AlibavaBaseProcessor("AlibavaCommonModeSubtraction"),
_commonmodeCollectionName(ALIBAVA::NOTSET),
_commonmodeerrorCollectionName(ALIBAVA::NOTSET),
_chanDataHistoName ("Common_and_Pedestal_subtracted_data_channel")
{
	
	// modify processor description
	_description =
	"AlibavaCommonModeSubtraction subtracts the provided common mode values from the input reco (pedestal subtracted) data. ";
	
	
	// first of register the input collection
	registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
									  "Input reco data collection name",
									  _inputCollectionName, string("recodata") );

	registerOutputCollection (LCIO::TRACKERDATA, "OutputCollectionName",
									 "Output data collection name",
									 _outputCollectionName, string("recodata_cmmd") );

	
	// if needed one can change these to optional parameters

	// now the optional parameters
	registerProcessorParameter ("CommonModeCollectionName",
										"Common mode collection name, better not to change",
										_commonmodeCollectionName, string ("commonmode"));
	
	registerProcessorParameter ("CommonModeErrorCollectionName",
										"Common mode error collection name, better not to change",
										_commonmodeerrorCollectionName, string ("commonmodeerror"));
	
}


void AlibavaCommonModeSubtraction::init () {
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

}
void AlibavaCommonModeSubtraction::processRunHeader (LCRunHeader * rdr) {
	streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;

	// Add processor name to the runheader
	auto_ptr<AlibavaRunHeaderImpl> arunHeader ( new AlibavaRunHeaderImpl(rdr)) ;
	arunHeader->addProcessor(type());
	
	
	setChipSelection( arunHeader->getChipSelection() );
	// set channels to be used (if it is defined)
	setChannelsToBeUsed();
	
	// if you want
	bookHistos();

	// set number of skipped events to zero (defined in AlibavaBaseProcessor)
	_numberOfSkippedEvents = 0;

}


void AlibavaCommonModeSubtraction::processEvent (LCEvent * anEvent) {

	if ( anEvent->getEventNumber() % 1000 == 0 )
		streamlog_out ( MESSAGE4 ) << "Looping events "<<anEvent->getEventNumber() << endl;
	
	AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
	
	if (_skipMaskedEvents && (alibavaEvent->isEventMasked()) ) {
		_numberOfSkippedEvents++;
		return;
	}

	LCCollectionVec * dataColVec;
	LCCollectionVec * cmmdColVec;
	LCCollectionVec * newColVec = new LCCollectionVec(LCIO::TRACKERDATA);

	CellIDEncoder<TrackerDataImpl> chipIDEncoder(ALIBAVA::ALIBAVADATA_ENCODE,newColVec);

	unsigned int noOfChips;
	try
	{
		dataColVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( getInputCollectionName() ) ) ;
		cmmdColVec =dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( _commonmodeCollectionName ) );

		// check these collections have same number of elements
		if ( (dataColVec->getNumberOfElements()) != (cmmdColVec->getNumberOfElements()) ) {
			streamlog_out( ERROR5 ) << "Number of elements in collections are not equal!" <<endl;
			streamlog_out( ERROR5 ) << getInputCollectionName() << " has " << dataColVec->getNumberOfElements() << " elements while "<< _commonmodeCollectionName <<" has "<< cmmdColVec->getNumberOfElements() << endl;
		}
		noOfChips = dataColVec->getNumberOfElements();
		
		
		for ( size_t i = 0; i < noOfChips; ++i )
		{
			// get data from the collection
			TrackerDataImpl * dataImpl = dynamic_cast< TrackerDataImpl * > ( dataColVec->getElementAt( i ) ) ;
			TrackerDataImpl * cmmdImpl = dynamic_cast< TrackerDataImpl * > ( cmmdColVec->getElementAt( i ) ) ;
			TrackerDataImpl * newdataImpl = new TrackerDataImpl();

			// check that they belong to same chip
			if ( (getChipNum(dataImpl)) != (getChipNum(cmmdImpl)) ) {
				streamlog_out( ERROR5 ) << "The chip numbers in the collections is not same! " << endl;
			}
			int chipnum = getChipNum(dataImpl);

			// set chip number for newdataImpl
			chipIDEncoder[ALIBAVA::ALIBAVADATA_ENCODE_CHIPNUM] = chipnum;
			chipIDEncoder.setCellID(newdataImpl);
			
			
			FloatVec datavec, cmmdvec, newdatavec;
			datavec = dataImpl->getChargeValues();
			cmmdvec = cmmdImpl->getChargeValues();
			newdatavec.clear();
			
			
			// check size of data sets are equal to ALIBAVA::NOOFCHANNELS
			if ( int(datavec.size()) != ALIBAVA::NOOFCHANNELS )
				streamlog_out( ERROR5 ) << "Number of channels in input data is not equal to ALIBAVA::NOOFCHANNELS! "<< endl;
			if ( int(cmmdvec.size()) != ALIBAVA::NOOFCHANNELS )
				streamlog_out( ERROR5 ) << "Number of channels in common mode data is not equal to ALIBAVA::NOOFCHANNELS! " << endl;

			
			// now subtract common mode values from all channels
			for (size_t ichan=0; ichan<datavec.size();ichan++) {
				if(isMasked(chipnum, ichan)) {
					newdatavec.push_back(0);
					continue;
				}

				float newdata = datavec[ichan] - cmmdvec[ichan];
				
				//store it in new data vector
				newdatavec.push_back(newdata);
				
			}
			
			newdataImpl->setChargeValues(newdatavec);
			newColVec->push_back(newdataImpl);
						
			fillHistos(newdataImpl);
		}
		alibavaEvent->addCollection(newColVec, getOutputCollectionName());

		
	} catch ( lcio::DataNotAvailableException ) {
		// do nothing again
		streamlog_out( ERROR5 ) << "Collection ("<<getInputCollectionName()<<") not found! " << endl;
	}
	
}

void AlibavaCommonModeSubtraction::check (LCEvent * /* evt */ ) {
	// nothing to check here - could be used to fill check plots in reconstruction processor
}


void AlibavaCommonModeSubtraction::end() {
	
	if (_numberOfSkippedEvents > 0)
		streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents<<" events skipped since they are masked" << endl;

	streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
	
}

void AlibavaCommonModeSubtraction::fillHistos(TrackerDataImpl * trkdata){

	// Fill the histograms with the corrected data

	FloatVec datavec;
	datavec = trkdata->getChargeValues();
	int chipnum = getChipNum(trkdata);

	for ( size_t ichan = 0 ; ichan < datavec.size() ; ichan++ )
	{
		if ( isMasked(chipnum, ichan) ) continue;
		
		string tempHistoName = getChanDataHistoName(chipnum, ichan);
		if ( TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[tempHistoName]) )
			histo->Fill(datavec[ichan]);
		
		string tempHistoName1 = getSignalCorrectionName();
		if ( TH1D * histo1 = dynamic_cast<TH1D*> (_rootObjectMap[tempHistoName1]) )
			histo1->Fill(datavec[ichan]);
	}

}

string AlibavaCommonModeSubtraction::getChanDataHistoName(int chipnum, int ichan){
	stringstream s;
	s << _chanDataHistoName << "_chip_"<<chipnum<<"_chan_" << ichan;
	return s.str();
}

string AlibavaCommonModeSubtraction::getSignalCorrectionName(){
	string s;
	s = "Final Pedestal Common Mode Corrected Signal";
	return s;
}

void AlibavaCommonModeSubtraction::bookHistos(){

	string tempHistoName;

	// a histogram showing the corrected signals
	tempHistoName = getSignalCorrectionName();
	stringstream tempHistoTitle1;
	tempHistoTitle1 << tempHistoName << ";ADCs;NumberofEntries";

	TH1D * signalHisto =
	new TH1D (tempHistoName.c_str(),"",2000,-1000,1000);
	_rootObjectMap.insert(make_pair(tempHistoName, signalHisto));
	string tmp_string1 = tempHistoTitle1.str();
	signalHisto->SetTitle(tmp_string1.c_str());


	AIDAProcessor::tree(this)->cd(this->name());
		
	EVENT::IntVec chipVec = getChipSelection();
	
	// here are the histograms for each channel to show the corrected data
	
	for (unsigned int i=0; i<chipVec.size(); i++) {
		int chipnum = chipVec[i];
		
		for (int ichan=0; ichan<ALIBAVA::NOOFCHANNELS; ichan++) {
			if (isMasked(chipnum, ichan)) {
				continue;
			}
			tempHistoName = getChanDataHistoName(chipnum, ichan);
			stringstream tempHistoTitle;
			tempHistoTitle << tempHistoName << ";ADCs;NumberofEntries";
			
			TH1D * chanDataHisto =
			new TH1D (tempHistoName.c_str(),"",2000,-1000,1000);
			_rootObjectMap.insert(make_pair(tempHistoName, chanDataHisto));
			string tmp_string = tempHistoTitle.str();
			chanDataHisto->SetTitle(tmp_string.c_str());
		}
		

	}
	
		streamlog_out ( MESSAGE1 )  << "End of booking histograms. " << endl;
}

