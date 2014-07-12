#include "EUTelProcessorTrackingHelixTrackSearch.h"
/**  EUTelProcessorTrackingHelixTrackSearch
 * 
 *  If compiled with MARLIN_USE_AIDA 
 *  it creates histograms.
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Requires collection of hits from which track
 *  candidates are build.
 *  Note we do not have to consider any other possible input.
 *  Since we need only the basic functions of a hit in all cases.
 *  <h4>Output</h4> 
 *  <li> Histograms.
 *  <li> Collection of track candidates.
 */
//If you are using static variables then for the most part you are doing something wrong. 
//TO DO: Create new Histograming procedure! 
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelProcessorTrackingHelixTrackSearch::_histName::_numberTracksCandidatesHistName = "NumberTracksCandidates";
std::string EUTelProcessorTrackingHelixTrackSearch::_histName::_numberOfHitOnTrackCandidateHistName = "NumberOfHitsOnTrackCandidate";
std::string EUTelProcessorTrackingHelixTrackSearch::_histName::_HitOnTrackCandidateHistName = "HitsOnTrackCandidate";
std::string EUTelProcessorTrackingHelixTrackSearch::_histName::_chi2CandidateHistName = "chi2CandidateHistName";
#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

/** Default constructor */
EUTelProcessorTrackingHelixTrackSearch::EUTelProcessorTrackingHelixTrackSearch() :
Processor("EUTelProcessorTrackingHelixTrackSearch"),
_hitInputCollectionName("HitCollection"),
_trackCandidateHitsOutputCollectionName("TrackCandidatesCollection"),
_trackFitter(0),
_maxMissingHitsPerTrackCand(0),
_AllowedSharedHitsOnTrackCandidate(0),
_maxNTracks(10),
_eBeam(-1.),
_qBeam(-1.),
_eBeamUncertatinty(0.),
_beamSpread(2,-1.),
_histoInfoFileName("histoinfo.xml"),
_nProcessedRuns(0),
_nProcessedEvents(0),
_aidaHistoMap1D() {

	_description = "EUTelProcessorTrackingHelixTrackSearch preforms track pattern recognition.";

	registerInputCollection (LCIO::TRACKERDATA, "ZSDataCollectionName", "LCIO converted data files",_zsDataCollectionName, string("original_zsdata"));


	// TrackerHit input collection
	registerInputCollection(LCIO::TRACKERHIT,"HitInputCollectionName","Input hits collection name",_hitInputCollectionName,std::string("HitCollection"));

	// TrackerHit output collection
	registerInputCollection(LCIO::TRACKERHIT,"TrackerHitOutputCollectionName","Pattern recognition track - output hits collection name",_hitFittedOutputCollectionName,std::string("HitFittedCollection"));

	// Track candidate hits output collection
	registerOutputCollection(LCIO::TRACK,"TrackCandHitOutputCollectionName","Output track candidates hits collection name",_trackCandidateHitsOutputCollectionName,std::string("TrackCandidateHitCollection"));

	registerOptionalParameter("MaxMissingHitsPerTrack", "Maximal number of missing hits on a track candidate",
	_maxMissingHitsPerTrackCand, static_cast<int> (0)); // Search full-length tracks by default

	registerOptionalParameter("AllowedSharedHitsOnTrackCandidate", "The number of similar hit a track candidate can have to another track candidate, within the same event.",
	_AllowedSharedHitsOnTrackCandidate, static_cast<int> (0));

	registerOptionalParameter("MaxNTracksPerEvent", "Maximal number of track candidates to be found in events",_maxNTracks, static_cast<int> (100));

	registerOptionalParameter("NumberOfPlanesToProject", "Maximal number of track candidates to be found in events",_planesProject, static_cast<int> (1));

	registerOptionalParameter("ResidualsRMax", "Maximal allowed distance between hits entering the recognition step "
	"per 10 cm space between the planes. One value for each neighbor planes. "
	"DistanceMax will be used for each pair if this vector is empty. Units are mm.",
	_residualsRMax, static_cast<double> (10));

	registerOptionalParameter("BeamEnergy", "Beam energy [GeV]", _eBeam, static_cast<double> (4.0));

	registerOptionalParameter("BeamCharge", "Beam charge [e]", _qBeam, static_cast<double> (-1));

	registerOptionalParameter("BeamSpread", "Angular spread of the beam (horizontal,vertical) [mrad] (for beam constraint). ""No beam constraints if negative values are supplied.",_beamSpread, EVENT::FloatVec(2,100.) );

	registerOptionalParameter("BeamEnergyUncertainty", "Uncertainty of beam energy [%]", _eBeamUncertatinty, static_cast<double> (0.) );

	// Histogram information

	registerOptionalParameter("HistogramInfoFilename", "Name of histogram info xml file", _histoInfoFileName, std::string("histoinfo.xml"));

	registerOptionalParameter("PlanesToCreateSeedsFrom", "This is the planes you want to create seeds from", _createSeedsFromPlanes,IntVec());

	registerOptionalParameter("ExcludePlanes", "This is the planes that will not be included in analysis", _excludePlanes ,FloatVec());

}

void EUTelProcessorTrackingHelixTrackSearch::init() {

	streamlog_out(DEBUG) << "EUTelProcessorTrackingHelixTrackSearch::init( )" << std::endl;

	_nProcessedRuns = 0;
	_nProcessedEvents = 0;

	std::string name("test.root"); //This is the name outputed at the end to store geo info.
	geo::gGeometry().initializeTGeoDescription(name,false);
	geo::gGeometry().initialisePlanesNotExcluded(_excludePlanes);
	// Instantiate track finder. This is a working horse of the processor.
	streamlog_out(DEBUG) << "Initialisation of track finder" << std::endl;

	EUTelKalmanFilter* Finder = new EUTelKalmanFilter("KalmanTrackFinder");
	if (!Finder) {
		streamlog_out(ERROR) << "Can't allocate an instance of EUTelExhaustiveTrackFinder. Stopping ..." << std::endl;
		throw UnknownDataTypeException("Track finder was not created");
	}
	 
	Finder->setAllowedMissingHits( _maxMissingHitsPerTrackCand );
	Finder->setAllowedSharedHitsOnTrackCandidate( _AllowedSharedHitsOnTrackCandidate );
	Finder->setWindowSize( _residualsRMax );
	Finder->setPlanesToCreateSeedsFrom(_createSeedsFromPlanes);
	Finder->setExcludePlanes(_excludePlanes);

	Finder->setBeamMomentum( _eBeam );
	Finder->setBeamCharge( _qBeam );
	Finder->setBeamMomentumUncertainty( _eBeamUncertatinty );
	Finder->setBeamSpread( _beamSpread );

	_trackFitter = Finder;
	_trackFitter->setAutoPlanestoCreateSeedsFrom();
	_trackFitter->testUserInput();
	// Book histograms. Yet again this should be replaced. TO DO:Create better histogram method.
 	bookHistograms();
}

void EUTelProcessorTrackingHelixTrackSearch::processRunHeader(LCRunHeader* run) {

	auto_ptr<EUTelRunHeaderImpl> header(new EUTelRunHeaderImpl(run));
	header->addProcessor(type());//Add what processor has acted to collection here. 

	// this is the right place also to check the geometry ID. This is a
	// unique number identifying each different geometry used at the
	// beam test. The same number should be saved in the run header and
	// in the xml file. If the numbers are different, warn the user.

	if (header->getGeoID() == 0)
		streamlog_out(WARNING0) << "The geometry ID in the run header is set to zero." << endl
		<< "This may mean that the GeoID parameter was not set" << endl;


	if (header->getGeoID() != geo::gGeometry().getSiPlanesLayoutID()) {
		streamlog_out(WARNING5) << "Error during the geometry consistency check: " << endl
				<< "The run header says the GeoID is " << header->getGeoID() << endl
				<< "The GEAR description says is     " << geo::gGeometry().getSiPlanesLayoutID() << endl;
	}

		_nProcessedRuns++;
}

void EUTelProcessorTrackingHelixTrackSearch::processEvent(LCEvent * evt) {

	EUTelEventImpl * event = static_cast<EUTelEventImpl*> (evt); //Change the LCIO object to EUTel object. This is a simple way to extend functionality of the object.

	// Do not process last or unknown events
	if (event->getEventType() == kEORE) {
		streamlog_out(DEBUG4) << "EORE found: nothing else to do." << std::endl;
		return;
	} else if (event->getEventType() == kUNKNOWN) {
		streamlog_out(WARNING2) << "Event number " << event->getEventNumber() << " in run " << event->getRunNumber()<< " is of unknown type. Continue considering it as a normal Data Event." << std::endl;
	}

	LCCollection* hitMeasuredCollection = NULL;
	try {
		hitMeasuredCollection = evt->getCollection(_hitInputCollectionName);
		streamlog_out(DEBUG1) << "collection : " <<_hitInputCollectionName << " retrieved" << std::endl;
	} catch (DataNotAvailableException e) {
		streamlog_out(WARNING2) << _hitInputCollectionName << " collection not available" << std::endl;
		throw marlin::SkipEventException(this);
	}

	// this will only be entered if the collection is available
	if ( hitMeasuredCollection == NULL) {
		streamlog_out(DEBUG2) << "EUTelProcessorTrackingHelixTrackSearch :: processEvent() hitMeasuredCollection is void" << std::endl;
		throw marlin::SkipEventException(this);
	}
	streamlog_out(DEBUG2) << "EUTelProcessorTrackingHelixTrackSearch::processEvent()" << std::endl;

	// Prepare hits for track finder
	EVENT::TrackerHitVec allHitsVec;
	_trackFitter->findHitsOrderVec(hitMeasuredCollection,allHitsVec); 
	_trackFitter->setHitsVec(allHitsVec);
	_trackFitter->testHitsVec();//Here we simply make sure the vector contains hits.
  _trackFitter->printHits();
	_trackFitter->setHitsVecPerPlane();//set hits vectors ordered by plane(Determined by geometry)
	_trackFitter->testHitsVecPerPlane();
	_trackFitter->onlyRunOnce();//This will only execute once. It can not be placed in init since it needs hit information. It currently only determines hit dimension. However can be used for other thing latter. 
	_trackFitter->testPlaneDimensions();
	streamlog_out(DEBUG1) << "Event #" << _nProcessedEvents << std::endl;
	streamlog_out( DEBUG1 ) << "Trying to find tracks..." << endl;
	_trackFitter->initialiseSeeds();
	_trackFitter->testInitialSeeds();
// searching for hits along the expected track direction 
	_trackFitter->findTrackCandidates( );
	_trackFitter->findTracksWithEnoughHits();
// remove possible duplicates (the amount of commont hits on the split tracks is controled via the processor paraemter)
	_trackFitter->findTrackCandidatesWithSameHitsAndRemove();

	streamlog_out( DEBUG1 ) << "Retrieving track candidates..." << endl;

std::vector< EUTelTrackImpl* >& trackCartesian = static_cast < EUTelKalmanFilter* > (_trackFitter)->getTracks(); //Need to cast this since error: ‘class eutelescope::EUTelTrackFitter’ otherwise Why??

plotHistos(trackCartesian);

outputLCIO(evt,trackCartesian);

_nProcessedEvents++;

//  if (isFirstEvent()) _isFirstEvent = false; Is this needed
}

void EUTelProcessorTrackingHelixTrackSearch::outputLCIO(LCEvent* evt, std::vector< EUTelTrackImpl* >& trackCartesian){

        streamlog_out( DEBUG4 ) << " ---------------- EUTelProcessorTrackingHelixTrackSearch::outputLCIO ---------- BEGIN ------------- " << std::endl;

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
  	evt->addCollection(trkCandCollection, _trackCandidateHitsOutputCollectionName);

        streamlog_out( DEBUG4 ) << " ---------------- EUTelProcessorTrackingHelixTrackSearch::outputLCIO ---------- END ------------- " << std::endl;
}

/** 
 * Plot few histos.
 * 
 */
void EUTelProcessorTrackingHelixTrackSearch::plotHistos( vector< EUTelTrackImpl* >& trackCandidates)  {

            const int nTracks = trackCandidates.size( );
 
            static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_numberTracksCandidatesHistName ] ) -> fill( nTracks );
            streamlog_out( MESSAGE2 ) << "Event #" << _nProcessedEvents << endl;

            int nHitsOnTrack = 0;
            vector< EUTelTrackImpl* >::const_iterator itrk;
            for ( itrk = trackCandidates.begin( ) ; itrk != trackCandidates.end( ); ++itrk ) {
						//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Here we fill the chi2 of the candidate tracks
            static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_chi2CandidateHistName ] ) -> fill( (*itrk)->getChi2() );
						////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


                const EVENT::TrackerHitVec trkHits = ( *itrk )->getHitsOnTrack();
                nHitsOnTrack = trkHits.size( );
                EVENT::TrackerHitVec::const_iterator itTrkHits;
                streamlog_out( MESSAGE1 ) << "Track hits start:==============" << std::endl;
                for ( itTrkHits = trkHits.begin( ) ; itTrkHits != trkHits.end( ); ++itTrkHits ) {
                    const double* uvpos = ( *itTrkHits )->getPosition( );
                    const int sensorID = Utility::getSensorIDfromHit( static_cast<IMPL::TrackerHitImpl*> (*itTrkHits) );
                    streamlog_out( MESSAGE1 ) << "Hit (id=" << setw(3) << sensorID << ") local(u,v) coordinates: (" 
                             << setw(7) << setprecision(4) << uvpos[0] << "," <<  setw(7) << setprecision(4) << uvpos[1] << ")";
                    double globalHit[] = {0.,0.,0.};
                    geo::gGeometry().local2Master( sensorID, uvpos, globalHit);
                    streamlog_out( MESSAGE1 ) << " WorldC: " << setw(7) << globalHit[0] << setw(7) << globalHit[1] << setw(7) << globalHit[2] ;
                    streamlog_out( MESSAGE1 )  << std::endl;
                    static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_HitOnTrackCandidateHistName ] ) -> fill( sensorID );
                }
                streamlog_out( MESSAGE1 ) << "Track hits end:==============" << std::endl;
                static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_numberOfHitOnTrackCandidateHistName ] ) -> fill( nHitsOnTrack );
            }
}

/**
 * Dump track candidate fitted hits into lcio collection.
 * 
 * @param evt event pointer
 * @param trackCandidateHitFirred  vector of hits from Pattern recognition
 */
void EUTelProcessorTrackingHelixTrackSearch::addTrackCandidateHitFittedToCollection(LCEvent* evt, EVENT::TrackerHitVec& trackCandidateHitFitted ) {
    // Prepare output collection

    // Try to access Output collection
    LCCollectionVec* hitFittedCollection = NULL;
    try {
        hitFittedCollection =  static_cast<LCCollectionVec*> ( evt->getCollection(_hitFittedOutputCollectionName) ) ;
        streamlog_out(DEBUG1) << "collection : " <<_hitFittedOutputCollectionName << " retrieved" << std::endl;
    } catch (DataNotAvailableException e) {
        streamlog_out(WARNING2) << _hitFittedOutputCollectionName << " collection not available, creating one ... " << std::endl;
        hitFittedCollection = new LCCollectionVec(LCIO::TRACKERHIT);
    }

    // Fill 
    EVENT::TrackerHitVec::iterator ihit;
    for ( ihit = trackCandidateHitFitted.begin(); ihit != trackCandidateHitFitted.end(); ++ihit ) {
        streamlog_out( MESSAGE1 ) << "hit " << (*ihit)->getType() << " hits" << endl;
        hitFittedCollection->push_back( (*ihit) );

    } // for (size_t itrk = 0; itrk < trackCandidates.size(); itrk++) 


    // Write track candidates collection
    try {
        streamlog_out(MESSAGE1) << "Getting collection " << _hitFittedOutputCollectionName << endl;
        evt->getCollection(_trackCandidateHitsOutputCollectionName);
    } catch (...) {
        streamlog_out(MESSAGE1) << "Adding collection " << _hitFittedOutputCollectionName << endl;
        evt->addCollection( hitFittedCollection, _hitFittedOutputCollectionName);
    }

}


/**
 * Dump track candidate into lcio collection.
 * 
 * @param evt event pointer
 * @param trackCandidates  vectors of hits assigned to track candidates
 */
void EUTelProcessorTrackingHelixTrackSearch::addTrackCandidateToCollection1(LCEvent* evt, std::vector< IMPL::TrackImpl* >& trackCandidates ) {
    // Prepare output collection
    LCCollectionVec * trkCandCollection = 0;
    try {
        trkCandCollection = new LCCollectionVec(LCIO::TRACK);
        LCFlagImpl flag(trkCandCollection->getFlag());
        flag.setBit( LCIO::TRBIT_HITS );
        trkCandCollection->setFlag( flag.getFlag( ) );
    } catch (...) {
        streamlog_out(WARNING2) << "Can't allocate output collection" << endl;
    }

    // Fill track parameters with nonsense except of hits
    std::vector< IMPL::TrackImpl* >::iterator itrk;
    for ( itrk = trackCandidates.begin(); itrk != trackCandidates.end(); ++itrk ) {
        streamlog_out( MESSAGE1 ) << "Track has " << (*itrk)->getTrackerHits().size() << " hits" << endl;
        trkCandCollection->push_back( (*itrk) );

    } // for (size_t itrk = 0; itrk < trackCandidates.size(); itrk++) 

    // Write track candidates collection
    try {
        streamlog_out(MESSAGE1) << "Getting collection " << _trackCandidateHitsOutputCollectionName << endl;
        evt->getCollection(_trackCandidateHitsOutputCollectionName);
    } catch (...) {
        streamlog_out(MESSAGE1) << "Adding collection " << _trackCandidateHitsOutputCollectionName << endl;
        evt->addCollection(trkCandCollection, _trackCandidateHitsOutputCollectionName);
    }

}

/**
 * Dumpt track candidate into lcio collection.
 * 
 * @param evt event pointer
 * @param trackCandidates  vectors of hits assigned to track candidates
 */
void EUTelProcessorTrackingHelixTrackSearch::addTrackCandidateToCollection(LCEvent* evt, const vector< EVENT::TrackerHitVec >& trackCandidates) {
    // Prepare output collection
    LCCollectionVec * trkCandCollection = 0;
    try {
        trkCandCollection = new LCCollectionVec(LCIO::TRACK);
        LCFlagImpl flag(trkCandCollection->getFlag());
        flag.setBit( LCIO::TRBIT_HITS );
        trkCandCollection->setFlag( flag.getFlag( ) );
    } catch (...) {
        streamlog_out(WARNING2) << "Can't allocate output collection" << endl;
    }

    // Fill track parameters with nonsense except of hits
    for (size_t itrk = 0; itrk < trackCandidates.size(); itrk++) {
        IMPL::TrackImpl* trackcand = new IMPL::TrackImpl;           // Don't free it manually, because it is owned by trkCandCollection

        // Assign hits to LCIO TRACK
        EVENT::TrackerHitVec trkcandhits = trackCandidates[itrk];
        for (size_t ihit = 0; ihit < trkcandhits.size(); ihit++) {
            TrackerHitImpl* hit = static_cast<TrackerHitImpl*>(trkcandhits[ihit]);
            trackcand->addHit( hit );
        }
        streamlog_out( DEBUG1 ) << "Track has " << trackcand->getTrackerHits().size() << " hits" << endl;
        trkCandCollection->push_back(trackcand);

    } // for (size_t itrk = 0; itrk < trackCandidates.size(); itrk++) 

    // Write track candidates collection
    try {
        streamlog_out(DEBUG1) << "Getting collection " << _trackCandidateHitsOutputCollectionName << endl;
        evt->getCollection(_trackCandidateHitsOutputCollectionName);
    } catch (...) {
        streamlog_out(DEBUG1) << "Adding collection " << _trackCandidateHitsOutputCollectionName << endl;
        evt->addCollection(trkCandCollection, _trackCandidateHitsOutputCollectionName);
    }

}

void EUTelProcessorTrackingHelixTrackSearch::check(LCEvent * /*evt*/) {
    // nothing to check here
}

void EUTelProcessorTrackingHelixTrackSearch::end() {

    delete _trackFitter;
    
    streamlog_out(MESSAGE4) << "EUTelProcessorTrackingHelixTrackSearch::end()  " << name()
            << " processed " << _nProcessedEvents << " events in " << _nProcessedRuns << " runs "
            << " av.tracks : " << static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_numberTracksCandidatesHistName ] ) -> mean()
            << " track candidates : " << static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_numberOfHitOnTrackCandidateHistName ] ) -> allEntries()
            << std::endl;
}

void EUTelProcessorTrackingHelixTrackSearch::bookHistograms() {
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

    try {
        streamlog_out(DEBUG) << "Booking histograms..." << std::endl;

        auto_ptr<EUTelHistogramManager> histoMgr( new EUTelHistogramManager( _histoInfoFileName ));
        EUTelHistogramInfo    * histoInfo;
        bool                    isHistoManagerAvailable;

        try {
            isHistoManagerAvailable = histoMgr->init( );
        } catch ( ios::failure& e ) {
            streamlog_out( ERROR5 ) << "I/O problem with " << _histoInfoFileName << "\n"
                    << "Continuing without histogram manager using default settings"    << endl;
            isHistoManagerAvailable = false;
        } catch ( ParseException& e ) {
            streamlog_out( ERROR5 ) << e.what( ) << "\n"
                    << "Continuing without histogram manager using default settings" << endl;
            isHistoManagerAvailable = false;
        }
        
        const int tracksNBin = 20;    
        const double tracksXMin = -0.5;
        const double tracksXMax = 19.5;

        histoInfo = histoMgr->getHistogramInfo(_histName::_numberTracksCandidatesHistName);        
        int NBin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : tracksNBin;
        double XMin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : tracksXMin;
        double XMax = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : tracksXMax;

        // Number of track candidates from pattern recognition step
        AIDA::IHistogram1D * numberTracksCandidates =
                marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_numberTracksCandidatesHistName, NBin,
                XMin, XMax);

        if (numberTracksCandidates) {
            numberTracksCandidates->setTitle("Number of track candidates;N tracks;N events");
            _aidaHistoMap1D.insert(make_pair(_histName::_numberTracksCandidatesHistName, numberTracksCandidates));
        } else {
            streamlog_out(ERROR2) << "Problem booking the " << (_histName::_numberTracksCandidatesHistName) << endl;
            streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
	        }
			const int chiNbins = 5000;
			const double chiXmin = 0;
			const double chiXmax = 500000;												


				///////////////////////////////////////////////////////////////////////////////////////////////////////////////Chi2 create plot. Useful to determine the behaviour of the pattern recognition
        histoInfo = histoMgr->getHistogramInfo(_histName::_chi2CandidateHistName);
//      NBin =  histoInfo->_xBin; 
//      XMin =  histoInfo->_xMin; 
//      XMax = histoInfo->_xMax;
        NBin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : chiNbins;
        XMin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin :chiXmin;
        XMax = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : chiXmax;

        AIDA::IHistogram1D * chi2 =
                marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_chi2CandidateHistName, NBin,
                XMin, XMax);

            numberTracksCandidates->setTitle("Number of track candidates;N tracks;N events");
            _aidaHistoMap1D.insert(make_pair(_histName::_chi2CandidateHistName, chi2));
				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				  
				
         
        const int hitsNBin = 10;
        const double hitsMin = -0.5;
        const double hitsMax = 9.5;
        
        histoInfo = histoMgr->getHistogramInfo(_histName::_numberOfHitOnTrackCandidateHistName);        
        NBin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : hitsNBin;
        XMin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : hitsMin;
        XMax = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : hitsMax;
        
        // Number of hit per track candidate
        AIDA::IHistogram1D * numberHitsOnTrackCandidates =
                marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_numberOfHitOnTrackCandidateHistName, NBin,
                XMin, XMax);

        if (numberHitsOnTrackCandidates) {
            numberHitsOnTrackCandidates->setTitle("Number of hits on track candidates;N hits;N tracks");
            _aidaHistoMap1D.insert(make_pair(_histName::_numberOfHitOnTrackCandidateHistName, numberHitsOnTrackCandidates));
        } else {
            streamlog_out(ERROR2) << "Problem booking the " << (_histName::_numberOfHitOnTrackCandidateHistName) << endl;
            streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        }

        const int    nhitsNBin = 40;
        const double nhitsMin = -0.5;
        const double nhitsMax =39.5;
        
        histoInfo = histoMgr->getHistogramInfo(_histName::_HitOnTrackCandidateHistName);        
        NBin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : nhitsNBin;
        XMin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : nhitsMin;
        XMax = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : nhitsMax;
        
        // Number of hit per track candidate
        AIDA::IHistogram1D * HitsOnTrackCandidates =
                marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_HitOnTrackCandidateHistName, NBin,
                XMin, XMax);

        if (HitsOnTrackCandidates) {
            HitsOnTrackCandidates->setTitle("hits on track candidates;N hits;N tracks");
            _aidaHistoMap1D.insert(make_pair(_histName::_HitOnTrackCandidateHistName, HitsOnTrackCandidates));
        } else {
            streamlog_out(ERROR2) << "Problem booking the " << (_histName::_HitOnTrackCandidateHistName) << endl;
            streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        }
        
    } catch (lcio::Exception& e) {
        streamlog_out(WARNING2) << "Can't allocate histograms. Continue without histogramming" << endl;
    }
#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
}

int EUTelProcessorTrackingHelixTrackSearch::getMaxTrackCandidates() const {
    return _maxNTracks;
}

int EUTelProcessorTrackingHelixTrackSearch::getAllowedMissingHits() const {
    return _maxMissingHitsPerTrackCand;
}
