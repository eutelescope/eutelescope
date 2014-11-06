#include "EUTelProcessorPlotTrack.h"

using namespace eutelescope;

EUTelProcessorPlotTrack::EUTelProcessorPlotTrack() :
Processor("EUTelProcessorPlotTrack"),
_trackCandidatesInputCollectionName("Default_input"),
_tracksOutputCollectionName("Default_output"){

	// Processor description
  _description = "EUTelProcessorGBLFitCandidates this will fit gbl tracks and output them into LCIO file.";

  // TrackerHit input collection
  registerInputCollection(LCIO::TRACK, "TrackCandidatesInputCollectionName", "Input track candidate collection name",_trackCandidatesInputCollectionName,std::string("TrackCandidatesCollection"));

  // Track output collection
  registerOutputCollection(LCIO::TRACK,"TracksOutputCollectionName","Output tracks collection name",_tracksOutputCollectionName,std::string("TrackCollection"));

  registerOptionalParameter("HistogramInfoFilename", "Name of histogram info xml file", _histoInfoFileName, std::string("histoinfo.xml"));






}


void EUTelProcessorPlotTrack::init(){
	streamlog_out(DEBUG5) << "Book Histograms!" << std::endl;
	EUTelHistogram histo("GBLFit", _histoInfoFileName);


}

void EUTelProcessorPlotTrack::processRunHeader(LCRunHeader * run) {}

void EUTelProcessorPlotTrack::check(LCEvent * evt){}

void EUTelProcessorPlotTrack::processEvent(LCEvent * evt){

	EUTelEventImpl * event = static_cast<EUTelEventImpl*> (evt); ///We change the class so we can use EUTelescope functions

	
	//////////////////////////////////////////////////////////////////////// Do not process last events
	if (event->getEventType() == kEORE) {
		streamlog_out(DEBUG4) << "EORE found: nothing else to do." << endl;
		return;
  }else if (event->getEventType() == kUNKNOWN) {
  	streamlog_out(WARNING2) << "Event number " << event->getEventNumber() << " in run " << event->getRunNumber() << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }
	////////////////////////////////////////////////////////////////////////
	
	//////////////////////////////////////////////////////////////////////////Try to access collection	
	LCCollection* col = NULL;
	try{
  	col = evt->getCollection(_trackCandidatesInputCollectionName);
 		streamlog_out(DEBUG1) << "collection : " << _trackCandidatesInputCollectionName << " retrieved" << std::endl;
  } catch (DataNotAvailableException e) {
  	streamlog_out(MESSAGE0) << _trackCandidatesInputCollectionName << " collection not available" << std::endl;
  	throw marlin::SkipEventException(this);
  }
	//////////////////////////////////////////////////////////////////////

	


}

void EUTelProcessorPlotTrack::end(){}




