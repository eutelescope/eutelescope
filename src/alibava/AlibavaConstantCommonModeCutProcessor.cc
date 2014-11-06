/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


// alibava includes ".h"
#include "AlibavaConstantCommonModeCutProcessor.h"
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
#include <UTIL/LCTime.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCEventImpl.h>


// ROOT includes ".h"
#include "TH1D.h"

// system includes <>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <memory>
#include <vector>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;


AlibavaConstantCommonModeCutProcessor::AlibavaConstantCommonModeCutProcessor () :
AlibavaBaseProcessor("AlibavaConstantCommonModeCutProcessor"),
_commonModeColName(ALIBAVA::NOTSET),
_commonModeCutMin(-10.0),
_commonModeCutMax(10.0),
_maskIfAnyChipsCommonModeIsNotInRange(true),
_numberOfMaskedEvents(0),
_totalNumberOfMaskedEvents(0),
_totalNumberOfEvents(0),
_hMaskedEventsNumberName("EventNumberOfMaskedEvents")
{
	
	// modify processor description
	_description =
	"AlibavaConstantCommonModeCutProcessor masks the events if their common mode correction value is not in the range specified by CommonModeCutMin and CommonModeCutMax.";

	// first of register the input collection
	registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
									 "Input collection name that will be masked if common mode noise is not in the range",
									 _inputCollectionName, string("recodata") );

	registerProcessorParameter ("CommonModeCollectionName",
										 "The common mode collection name",
										 _commonModeColName , string("commonmode"));
	
	registerProcessorParameter ("CommonModeCutMin",
										 "The minimum common mode noise that is acceptable to use that Event",
										 _commonModeCutMin , float(-10.0));

	registerProcessorParameter ("CommonModeCutMax",
										 "The maximum common mode noise that is acceptable to use that Event",
										 _commonModeCutMax , float(10.0));

	registerProcessorParameter ("CommonModeCutMax",
										 "The maximum common mode noise that is acceptable to use that Event",
										 _commonModeCutMax , float(10.0));

	registerProcessorParameter ("MaskIfAnyChipsCommonModeIsNotInRange",
										 "If MaskIfAnyChipsCommonModeIsNotInRange variable is set to \"true\", the whole event will be masked even if only one chip has common mode noise that is not in the range. In case it is set to \"false\" then the event will be masked if all chips selected has common mode noise not in the range.",
										  _maskIfAnyChipsCommonModeIsNotInRange , bool(true));


}


void AlibavaConstantCommonModeCutProcessor::init () {
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

	// this method is called only once even when the rewind is active
	// usually a good idea to
	printParameters ();

	
}
void AlibavaConstantCommonModeCutProcessor::processRunHeader (LCRunHeader * rdr) {
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
		
	// set number of skipped events to zero (defined in AlibavaBaseProcessor)
	_numberOfSkippedEvents = 0;
	
	// set number of masked events by this processor to zero
	_numberOfMaskedEvents = 0;
	// set total number of skipped events to zero
	_totalNumberOfMaskedEvents = 0;
	
	bookHistos();


}


void AlibavaConstantCommonModeCutProcessor::processEvent (LCEvent * anEvent) {
	
	if ( anEvent->getEventNumber() % 1000 == 0 )
		streamlog_out ( MESSAGE4 ) << "Looping events "<<anEvent->getEventNumber() << endl;
	
	AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
	
	// if event is already mask, no need to mask it again ;)
	if (alibavaEvent->isEventMasked()) {
		_totalNumberOfMaskedEvents++;
		return;
	}
	
	LCCollectionVec * dataColVec;
	LCCollectionVec * cmmdColVec;

	// this vector will contain common mode values from each chip
	vector<float> cmmdValues;
	cmmdValues.clear();
	
	unsigned int noOfChips;
	try
	{
		dataColVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( getInputCollectionName() ) ) ;
		cmmdColVec =dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( _commonModeColName ) );
		
		// check these collections have same number of elements
		if ( (dataColVec->getNumberOfElements()) != (cmmdColVec->getNumberOfElements()) ) {
			streamlog_out( ERROR5 ) << "Number of elements in collections are not equal!" <<endl;
			streamlog_out( ERROR5 ) << getInputCollectionName() << " has " << dataColVec->getNumberOfElements() << " elements while "<< _commonModeColName <<" has "<< cmmdColVec->getNumberOfElements() << endl;
		}
		noOfChips = dataColVec->getNumberOfElements();
		
		for ( size_t i = 0; i < noOfChips; ++i )
		{
			// get data from the collection
			TrackerDataImpl * dataImpl = dynamic_cast< TrackerDataImpl * > ( dataColVec->getElementAt( i ) ) ;
			TrackerDataImpl * cmmdImpl = dynamic_cast< TrackerDataImpl * > ( cmmdColVec->getElementAt( i ) ) ;
			
			// check that they belong to same chip
			if ( (getChipNum(dataImpl)) != (getChipNum(cmmdImpl)) ) {
				streamlog_out( ERROR5 ) << "The chip numbers in the collections is not same! " << endl;
			}
			int chipnum = getChipNum(dataImpl);
			
			FloatVec datavec, cmmdvec;
			datavec = dataImpl->getChargeValues();
			cmmdvec = cmmdImpl->getChargeValues();
			
			// check size of data sets are equal to ALIBAVA::NOOFCHANNELS
			if ( int (datavec.size()) != ALIBAVA::NOOFCHANNELS )
				streamlog_out( ERROR5 ) << "Number of channels in input data is not equal to ALIBAVA::NOOFCHANNELS! "<< endl;
			if ( int(cmmdvec.size()) != ALIBAVA::NOOFCHANNELS )
				streamlog_out( ERROR5 ) << "Number of channels in common mode data is not equal to ALIBAVA::NOOFCHANNELS! " << endl;
			
			
			// now find common mode value, it should be constant so we can use any not masked channel
			float cmmdValue_forThisChip = 0;
			for (size_t ichan=0; ichan<cmmdvec.size();ichan++) {
				if(isMasked(chipnum, ichan)) {
					continue;
				}
				cmmdValue_forThisChip = cmmdvec[ichan];
			}
			
			// contain it in cmmdValues
			cmmdValues.push_back(cmmdValue_forThisChip);
		} // end of loop ever chips
		
		// now loop over cmmdValues and count how many of those not in the range
		unsigned int Nchipsfailed = 0;
		for (unsigned int i=0; i<cmmdValues.size(); i++) {
			float cmmd = cmmdValues[i];
			if ( cmmd < getCommonModeCutMin() || cmmd > getCommonModeCutMax() ) {
				Nchipsfailed++;
			}
		}
		
		// now mask the event according to the rule specified by MaskIfAnyChipsCommonModeIsNotInRange
		if (  (_maskIfAnyChipsCommonModeIsNotInRange && Nchipsfailed > 0 ) ||
			 ( (!_maskIfAnyChipsCommonModeIsNotInRange) && Nchipsfailed == noOfChips ) ){
			alibavaEvent->maskEvent();
			_numberOfMaskedEvents++;
			_totalNumberOfMaskedEvents++;
			
			streamlog_out(DEBUG1)<< "Event masked with common mode values: ";
			for (unsigned int i=0; i<cmmdValues.size(); i++) {
				streamlog_out(DEBUG1)<< cmmdValues[i] << ", ";
			}
			streamlog_out(DEBUG1)<<endl;
		}
			
		// otherwise don't mask
		
		
	} catch ( lcio::DataNotAvailableException ) {
		// do nothing again
		streamlog_out( ERROR5 ) << "Collection ("<<getInputCollectionName()<<") not found! " << endl;
	}

}

float AlibavaConstantCommonModeCutProcessor::getCommonModeCutMin(){
	return _commonModeCutMin;
}
float AlibavaConstantCommonModeCutProcessor::getCommonModeCutMax(){
	return _commonModeCutMax;
}

void AlibavaConstantCommonModeCutProcessor::check (LCEvent * /* evt */ ) {
	// nothing to check here - could be used to fill check plots in reconstruction processor
}


void AlibavaConstantCommonModeCutProcessor::end() {

	streamlog_out ( MESSAGE5 ) << _numberOfMaskedEvents << " event masked using common mode cut "<<getCommonModeCutMin()<<" - "<<getCommonModeCutMax()<< endl;
	
	streamlog_out ( MESSAGE5 ) << "In total "<< _totalNumberOfMaskedEvents << " event masked out of "<<_totalNumberOfEvents<< endl;
	
	streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
	
}

void AlibavaConstantCommonModeCutProcessor::fillHistos(AlibavaEventImpl * anAlibavaEvent){

	int eventnum = anAlibavaEvent->getEventNumber();
	
	if (TH1D * histo = dynamic_cast<TH1D*>(_rootObjectMap[_hMaskedEventsNumberName]))
	histo->Fill(eventnum);
	
}

	
void AlibavaConstantCommonModeCutProcessor::bookHistos(){

	AIDAProcessor::tree(this)->cd(this->name());

	TH1D * _hMaskedEventsNumber = new TH1D (_hMaskedEventsNumberName.c_str(),"", 10000, 0, _totalNumberOfEvents);
	_rootObjectMap.insert(make_pair(_hMaskedEventsNumberName, _hMaskedEventsNumber));
	_hMaskedEventsNumber->SetTitle("Event Number of Masked Events By ConstantCommonModeCutProcessor; Number of Entries;Event Number");
	
	
	streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << endl;
}


