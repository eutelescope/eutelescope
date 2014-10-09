/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


// alibava includes ".h"
#include "AlibavaTimeCutProcessor.h"
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

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;


AlibavaTimeCutProcessor::AlibavaTimeCutProcessor () :
AlibavaBaseProcessor("AlibavaTimeCutProcessor"),
_timecutmin(0.0),
_timecutmax(100.0),
_numberOfMaskedEvents(0),
_hMaskedEventsNumberName("EventNumberOfMaskedEvents")
{
	
	// modify processor description
	_description =
	"AlibavaTimeCutProcessor masks the events if their TDC time value is not in the range specified by TimeCutMin and TimeCutMax ";

	registerProcessorParameter ("TimeCutMin",
										 "The minimum tdc time that is acceptable to use that Event",
										 _timecutmin , float(0.0));

	registerProcessorParameter ("TimeCutMax",
										 "The maximum tdc time that is acceptable to use that Event",
										 _timecutmax , float(100.0));

}


void AlibavaTimeCutProcessor::init () {
	streamlog_out ( MESSAGE4 ) << "Running init" << endl;

	// this method is called only once even when the rewind is active
	// usually a good idea to
	printParameters ();

	
}
void AlibavaTimeCutProcessor::processRunHeader (LCRunHeader * rdr) {
	streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;

	// Add processor name to the runheader
	auto_ptr<AlibavaRunHeaderImpl> arunHeader ( new AlibavaRunHeaderImpl(rdr)) ;
	arunHeader->addProcessor(type());
		
	_numberOfMaskedEvents = 0;

	bookHistos();

}


void AlibavaTimeCutProcessor::processEvent (LCEvent * anEvent) {
	
	//if ( anEvent->getEventNumber() % 1000 == 0 )
	//	streamlog_out ( MESSAGE4 ) << "Looping events "<<anEvent->getEventNumber() << endl;
	
	AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
	
	float eventTime = alibavaEvent->getEventTime();
	streamlog_out(DEBUG0)<<"Event Number: "<< anEvent->getEventNumber()<<" Event time: "<<eventTime<<endl;
	if( eventTime < getTimeCutMin() || eventTime > getTimeCutMax() ){
		streamlog_out(DEBUG0)<< "Event time is not in the specified range."<<endl;
		alibavaEvent->maskEvent();
		streamlog_out(DEBUG0)<<"The event masked!"<<endl;
		_numberOfMaskedEvents++;
	}

}

float AlibavaTimeCutProcessor::getTimeCutMin(){
	return _timecutmin;
}
float AlibavaTimeCutProcessor::getTimeCutMax(){
	return _timecutmax;
}

void AlibavaTimeCutProcessor::check (LCEvent * /* evt */ ) {
	// nothing to check here - could be used to fill check plots in reconstruction processor
}


void AlibavaTimeCutProcessor::end() {

	streamlog_out ( MESSAGE5 ) << _numberOfMaskedEvents << " event masked using timecut "<<getTimeCutMin()<<" - "<<getTimeCutMax()<< endl;
	
	streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
	
}

void AlibavaTimeCutProcessor::fillHistos(AlibavaEventImpl * anAlibavaEvent){

	int eventnum = anAlibavaEvent->getEventNumber();
	
	if (TH1D * histo = dynamic_cast<TH1D*>(_rootObjectMap[_hMaskedEventsNumberName]))
	histo->Fill(eventnum);
	
}

	
void AlibavaTimeCutProcessor::bookHistos(){

	AIDAProcessor::tree(this)->cd(this->name());

	TH1D * _hMaskedEventsNumber = new TH1D (_hMaskedEventsNumberName.c_str(),"", 10000, 0, 1000000);
	_rootObjectMap.insert(make_pair(_hMaskedEventsNumberName, _hMaskedEventsNumber));
	_hMaskedEventsNumber->SetTitle("Event Number of Masked Events; Number of Entries;Event Number");
	
	
	streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << endl;
}


