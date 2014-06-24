#ifdef USE_GBL    // Not sure where this is defined. However it is used in all the other GBL processors will use it.

#include "EUTelProcessorGBLFitCandidates.h"
#include "include/GblTrajectory.h"

using namespace eutelescope;

EUTelProcessorGBLFitCandidates::EUTelProcessorGBLFitCandidates() :
Processor("EUTelProcessorGBLFitCandidates"),
_trackCandidatesInputCollectionName("Default_input"),
_tracksOutputCollectionName("Default_output"),
_nProcessedRuns(0),
_nProcessedEvents(0),
_beamQ(-1),
_eBeam(4),
_mEstimatorType(),
_maxChi2Cut(1000),
_milleBinaryFilename("mille.bin"),
_alignmentMode(0)
{
	// Processor description
  _description = "EUTelProcessorGBLFitCandidates this will fit gbl tracks and output them into LCIO file.";

  // TrackerHit input collection
  registerInputCollection(LCIO::TRACK, "TrackCandidatesInputCollectionName", "Input track candidate collection name",_trackCandidatesInputCollectionName,std::string("TrackCandidatesCollection"));

  // Track output collection
  registerOutputCollection(LCIO::TRACK,"TracksOutputCollectionName","Output tracks collection name",_tracksOutputCollectionName,std::string("TrackCollection"));

  registerOptionalParameter("BeamCharge", "Beam charge [e]", _beamQ, static_cast<double> (-1));

  // Necessary processor parameters that define fitter settings
  registerProcessorParameter("BeamEnergy", "Beam energy [GeV]", _eBeam, static_cast<double> (4.0));

  // Optional processor parameters that define finder settings

  registerOptionalParameter("GBLMEstimatorType", "GBL outlier down-weighting option (t,h,c)", _mEstimatorType, string() );

  registerOptionalParameter("MilleMaxChi2Cut", "Maximum chi2 of a track candidate that goes into millepede", _maxChi2Cut, double(1000.));

  registerOptionalParameter("MilleBinaryFilename", "Name of the Millepede binary file", _milleBinaryFilename, std::string("mille.bin"));

    // MILLEPEDE specific parameters
    registerOptionalParameter("AlignmentMode", "Alignment mode specifies alignment degrees of freedom to be considered\n"
            "0 - No alignment at all. Simply fit tracks assuming that alignment is correct\n"
            "1 - Alignment of XY shifts\n"
            "2 - Alignment of XY shifts + rotations around Z\n"
            "3 - Alignment of XYZ shifts + rotations around Z\n"
            "4 - Alignment of XY shifts + rotations around X and Z\n"
            "5 - Alignment of XY shifts + rotations around Y and Z\n"
            "6 - Alignment of XY shifts + rotations around X,Y and Z\n"
            "7 - Alignment of XYZ shifts + rotations around X,Y and Z\n",
            _alignmentMode, static_cast<int> (0));



}

void EUTelProcessorGBLFitCandidates::init() {

	streamlog_out(DEBUG2) << "EUTelProcessorGBLFitCandidates::init( )---------------------------------------------BEGIN" << std::endl;

	// Reset counters
	_nProcessedRuns = 0;
	_nProcessedEvents = 0;

	// Getting access to geometry description
	std::string name("test.root");
	geo::gGeometry().initializeTGeoDescription(name,false);

	//Create the size of the jacobian and parameter list for alignment
	TMatrixD* alDer; // alignment derivatives
	std::vector<int>* globalLabels;

	// Initialize GBL fitter
	EUTelGBLFitter* Fitter = new EUTelGBLFitter("GBLFitter");
	EUTelMillepede* Mille  = new EUTelMillepede(_alignmentMode);
  Fitter->SetBeamCharge(_beamQ);
  Fitter->SetBeamEnergy(_eBeam);
	Fitter->setMEstimatorType(_mEstimatorType);
  Fitter->SetMilleBinaryName(_milleBinaryFilename);
  Fitter->SetChi2Cut(_maxChi2Cut);
	Fitter->SetMillepede(Mille);
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


  	if (header->getGeoID() != geo::gGeometry().getSiPlanesLayoutID()) {  
		streamlog_out(WARNING5) << "Error during the geometry consistency check: " << endl << "The run header says the GeoID is " << header->getGeoID() << endl << "The GEAR description says is     " << geo::gGeometry().getSiPlanesLayoutID() << endl;
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

	//This is a good point to clear all things that need to be reset for each event. Why should gop here?
        _trackFitter->Clear(); 


	// this will only be entered if the collection is available
  	if (col != NULL) {
  		streamlog_out(DEBUG2) << "Collection contains data! Continue!" << endl;

		////////////////////////////////////////////////////////////////////////////For the each event we get loop over all track candidates and fit them
		std::vector< EUTelTrackImpl* > EUtracks;	
		for (int iCol = 0; iCol < col->getNumberOfElements(); iCol++) {
			
			if (!col) {
		      		streamlog_out(WARNING2) << "Track collection not found found for event " << _nProcessedEvents << " in run " << _nProcessedRuns << endl;
		        	throw SkipEventException(this);
     			}
	   		IMPL::TrackImpl* trackimpl = static_cast<IMPL::TrackImpl*> (col->getElementAt(iCol));
			EUTelTrackImpl* EUtrack = new EUTelTrackImpl(*trackimpl);
      			streamlog_out(DEBUG1) << "Track " << iCol << " nhits " << trackimpl->getTrackerHits().size() << endl;

			//_trackFitter->Clear(); //This is a good point to clear all things that need to be reset for each event. Why should gop here?
			std::vector< gbl::GblPoint >* pointList;
      _trackFitter->FillInformationToGBLPointObject(EUtrack, pointList);

 			const gear::BField& B = geo::gGeometry().getMagneticFiled();
      const double Bmag = B.at( TVector3(0.,0.,0.) ).r2();


  		gbl::GblTrajectory* traj = 0;
      if ( Bmag < 1.E-6 ) {
      	traj = new gbl::GblTrajectory( *pointList, false ); //Must make sure this is not a memory leak
      } else {
      	traj = new gbl::GblTrajectory( *pointList, true );
      }
			double * chi2=0; 
			int* ndf=0;
			_trackFitter->CreateTrajectoryandFit(pointList,traj, chi2,ndf);

			//Update track and then state variables//////////////////////////////////////////////BEGIN
			EUtrack->setChi2(*chi2);
			EUtrack->setNdf(*ndf);
			_trackFitter->UpdateTrackFromGBLTrajectory(traj, pointList); 
			//////////////////////////////////////////////////////////////////////////////////////END
			EUtracks.push_back(EUtrack);
			
			      
		}//END OF LOOP FOR ALL TRACKS IN AN EVENT
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			outputLCIO(evt, EUtracks); // Important to note that EUtracks are still only the original states at input. The scatterers are used in the fit but are not included here.
	}//END OF COLLECTION IS NOT NULL LOOP	


}

void EUTelProcessorGBLFitCandidates::end() {}

#endif // USE_GBL

/////////////////////////////////////////////////////Functions THIS IS OLD WILL REMOVE

void EUTelProcessorGBLFitCandidates::outputLCIO(LCEvent* evt, std::vector< EUTelTrackImpl* >& trackCartesian){

        streamlog_out( DEBUG4 ) << " ---------------- EUTelProcessorGBLFitCandidates::outputLCIO ---------- BEGIN ------------- " << std::endl;

	//Create once per event    
	LCCollectionVec * trkCandCollection = new LCCollectionVec(LCIO::TRACK);

	// Prepare output collection
  	LCFlagImpl flag(trkCandCollection->getFlag());
  	flag.setBit( LCIO::TRBIT_HITS );
  	trkCandCollection->setFlag( flag.getFlag( ) );

	//Loop through all tracks
	vector< EUTelTrackImpl* >::const_iterator itTrackCartesian;
	for ( itTrackCartesian = trackCartesian.begin(); itTrackCartesian != trackCartesian.end(); itTrackCartesian++){

                if(streamlog_level(DEBUG4) ) (*itTrackCartesian)->Print();

		IMPL::TrackImpl* LCIOtrack = (*itTrackCartesian)->CreateLCIOTrack();
		

		//For every track add this to the collection
    		trkCandCollection->push_back( LCIOtrack );
	}//END TRACK LOOP

	//Now add this collection to the 
  	evt->addCollection(trkCandCollection, _tracksOutputCollectionName);

        streamlog_out( DEBUG4 ) << " ---------------- EUTelProcessorGBLFitCandidates::outputLCIO ---------- END ------------- " << std::endl;
}





void EUTelProcessorGBLFitCandidates::CreateEUTrackandStates(TrackImpl* trackimpl, EUTelTrackImpl* EUtrack){
	
  	for(int i=0;i < trackimpl->getTrackStates().size(); i++){
		EUTelTrackStateImpl *EUstate = new EUTelTrackStateImpl;

  		IMPL::TrackStateImpl* state = static_cast < IMPL::TrackStateImpl*> ( trackimpl->getTrackStates().at(i) ) ;
		EUstate->setX(state->getD0()); //x position global
		EUstate->setY(state->getPhi()); //y position global
		EUstate->setTx(state->getOmega()); //tx position global
		EUstate->setTy(state->getZ0()); //ty position global
		EUstate->setInvP(state->getTanLambda()); //invp position global

		EUstate->setbeamQ(_beamQ); //Beam charge
		
		EUstate->setReferencePoint(state->getReferencePoint());
		EUstate->setCovMatrix(state->getCovMatrix());			
		//EUstate->setHit(getHit()
		EUtrack->addTrackState(EUstate);	 
	}

    	const EVENT::TrackerHitVec& trkcandhits = trackimpl->getTrackerHits();
    	EVENT::TrackerHitVec::const_iterator itrHit;
    	for ( itrHit = trkcandhits.begin(); itrHit != trkcandhits.end(); ++itrHit ){
    		//EUtrack->addHit( *itrHit );
    	}


}



