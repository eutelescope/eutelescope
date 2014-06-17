#ifdef USE_GBL    // Not sure where this is defined. However it is used in all the other GBL processors will use it.

#include "EUTelProcessorGBLFitCandidates.h"

using namespace eutelescope;

EUTelProcessorGBLFitCandidates::EUTelProcessorGBLFitCandidates() :
Processor("EUTelProcessorGBLFitCandidates"),
_trackCandidatesInputCollectionName("Default_input"),
_tracksOutputCollectionName("Default_output"),
_nProcessedRuns(0),
_nProcessedEvents(0),
_beamQ(-1),
_eBeam(4)
{
	// Processor description
  _description = "EUTelProcessorTrackingGBLTrajectory performs track fits using GBL optionally writing data files for MILLEPEDE II.";

  // TrackerHit input collection
  registerInputCollection(LCIO::TRACK, "TrackCandidatesInputCollectionName", "Input track candidate collection name",_trackCandidatesInputCollectionName,std::string("TrackCandidatesCollection"));

  // Track output collection
  registerOutputCollection(LCIO::TRACK,"TracksOutputCollectionName","Output tracks collection name",_tracksOutputCollectionName,std::string("TrackCollection"));

  registerOptionalParameter("BeamCharge", "Beam charge [e]", _beamQ, static_cast<double> (-1));

  // Necessary processor parameters that define fitter settings
  registerProcessorParameter("BeamEnergy", "Beam energy [GeV]", _eBeam, static_cast<double> (4.0));



}

void EUTelProcessorGBLFitCandidates::init() {

	streamlog_out(DEBUG2) << "EUTelProcessorGBLFitCandidates::init( )---------------------------------------------BEGIN" << std::endl;

	// Reset counters
	_nProcessedRuns = 0;
	_nProcessedEvents = 0;


	// Getting access to geometry description
	std::string name("test.root");
	geo::gGeometry().initializeTGeoDescription(name,false);

	// Initialize GBL fitter
	EUTelGBLFitter* Fitter = new EUTelGBLFitter("GBLFitter");
  Fitter->SetBeamCharge(_beamQ);
  Fitter->SetBeamEnergy(_eBeam);
  _trackFitter = Fitter;


	if (!_trackFitter) {
		streamlog_out(ERROR) << "Can't allocate an instance of EUTelGBLFitter. Stopping ..." << std::endl;
		throw UnknownDataTypeException("Track finder was not created");
	}

    streamlog_out(DEBUG2) << "EUTelProcessorGBLFitCandidates::init( )---------------------------------------------END" << std::endl;


}

void EUTelProcessorGBLFitCandidates::processRunHeader(LCRunHeader * run) {

	auto_ptr<EUTelRunHeaderImpl> header(new EUTelRunHeaderImpl(run));
  header->addProcessor(type());


	// this is the right place also to check the geometry ID. This is a
  // unique number identifying each different geometry used at the
  // beam test. The same number should be saved in the run header and
 	// in the xml file. If the numbers are different, warn the user.

	if (header->getGeoID() == 0)
 		streamlog_out(WARNING0) << "The geometry ID in the run header is set to zero." << endl << "This may mean that the GeoID parameter was not set" << endl;


  if (header->getGeoID() != geo::gGeometry()._siPlanesParameters->getSiPlanesID()) { 
		streamlog_out(WARNING5) << "Error during the geometry consistency check: " << endl << "The run header says the GeoID is " << header->getGeoID() << endl << "The GEAR description says is     " << geo::gGeometry()._siPlanesParameters->getSiPlanesID() << endl;
  }
    
    _nProcessedRuns++;

}

void EUTelProcessorGBLFitCandidates::check(LCEvent * evt){}

void EUTelProcessorGBLFitCandidates::processEvent(LCEvent * evt){

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
	try {
  	col = evt->getCollection(_trackCandidatesInputCollectionName);
    streamlog_out(DEBUG1) << "collection : " << _trackCandidatesInputCollectionName << " retrieved" << std::endl;
  } catch (DataNotAvailableException e) {
  	streamlog_out(MESSAGE0) << _trackCandidatesInputCollectionName << " collection not available" << std::endl;
    throw marlin::SkipEventException(this);
  }
	///////////////////////////////////////////////////////////////////////////////

	// this will only be entered if the collection is available
  if (col != NULL) {
  	streamlog_out(DEBUG2) << "Collection contains data! Continue!" << endl;

		////////////////////////////////////////////////////////////////////////////For the each event we get loop over all track candidates and fit them	
    for (int iCol = 0; iCol < col->getNumberOfElements(); iCol++) {
			
    	if (!col) {
      	streamlog_out(WARNING2) << "Track collection not found found for event " << _nProcessedEvents << " in run " << _nProcessedRuns << endl;
        throw SkipEventException(this);
      }
    	IMPL::TrackImpl* trackimpl = static_cast<IMPL::TrackImpl*> (col->getElementAt(iCol));
			EUTelTrackImpl* EUtrack;
			CreateEUTrackandStates(trackimpl,EUtrack);
      streamlog_out(DEBUG1) << "Track " << iCol << " nhits " << trackimpl->getTrackerHits().size() << endl;

			//_trackFitter->Clear(); //This is a good point to clear all things that need to be reset for each event. Why should gop here?
      _trackFitter->FillInformationToGBLPointObject(EUtrack);      

		}//END OF LOOP FOR ALL TRACKS IN AN EVENT
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	}//END OF COLLECTION IS NOT NULL LOOP	

}

void EUTelProcessorGBLFitCandidates::end() {}

#endif // USE_GBL

/////////////////////////////////////////////////////Functions
void EUTelProcessorGBLFitCandidates::CreateEUTrackandStates(TrackImpl* trackimpl, EUTelTrackImpl* EUtrack){
	
  for(int i=0;i < trackimpl->getTrackStates().size(); i++){
		EUTelTrackStateImpl *EUstate = new EUTelTrackStateImpl;

  	IMPL::TrackStateImpl* state = static_cast < IMPL::TrackStateImpl*> ( trackimpl->getTrackStates().at(i) ) ;
		EUstate->setX(state->getD0()); //x position global
		EUstate->setY(state->getPhi()); //y position global
		EUstate->setTx(state->getOmega()); //tx position global
		EUstate->setTy(state->getZ0()); //ty position global
		EUstate->setInvP(state->getTanLambda()); //invp position global
		
		EUstate->setReferencePoint(state->getReferencePoint());
		EUstate->setCovMatrix(state->getCovMatrix());	

		EUtrack->addTrackState(EUstate);	 
	}

   	// Assign hits to LCIO TRACK
    const EVENT::TrackerHitVec& trkcandhits = trackimpl->getTrackerHits();
    EVENT::TrackerHitVec::const_iterator itrHit;
    for ( itrHit = trkcandhits.begin(); itrHit != trkcandhits.end(); ++itrHit ){
    	EUtrack->addHit( *itrHit );
    }

	
}


