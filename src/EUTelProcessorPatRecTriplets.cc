#include "EUTelProcessorPatRecTriplets.h"
///TO DO: Create new Histograming procedure! 

using namespace eutelescope;

std::string EUTelProcessorPatRecTriplets::_histName::_numberTracksCandidatesHistName = "NumberTracksCandidates";
std::string EUTelProcessorPatRecTriplets::_histName::_HitOnTrackCandidateHistName = "HitsOnTrackCandidate";
std::string EUTelProcessorPatRecTriplets::_histName::_HitOnTrackTimeHistName = "HitsOnTrackTime";
std::string EUTelProcessorPatRecTriplets::_histName::_chi2CandidateHistName = "chi2CandidateHistName";

/** Default constructor */
EUTelProcessorPatRecTriplets::EUTelProcessorPatRecTriplets() :
Processor("EUTelProcessorPatRecTriplets"),
_aidaHistoMap1D(),
_hitInputCollectionName("HitCollection"),
_trackCandidateHitsOutputCollectionName("TrackCandidatesCollection"),
_histoInfoFileName("histoinfo.xml"),
_trackFitter(0),
_doubletDistCut(0,0),
_tripletSlopeCuts(0,0),
_nProcessedRuns(0),
_nProcessedEvents(0),
_eBeam(-1.),
_qBeam(-1.),
_dataMissNumber(0),
_DoubletXseperationHistoRight(),
_DoubletYseperationHistoRight(),
_DoubletXseperationHistoLeft(),
_DoubletYseperationHistoLeft(),
_TripletXseperationHistoRight(),
_TripletYseperationHistoRight(),
_TripletXseperationHistoLeft(),
_TripletYseperationHistoLeft(),
_TripletDistCutXHisto(),
_TripletDistCutYHisto()
{
	///The standard description that comes with every processor 
	_description = "EUTelProcessorPatRecTriplets preforms track pattern recognition.";

	/// TrackerHit input collection
	registerInputCollection(LCIO::TRACKERHIT,"HitInputCollectionName","Input hits collection name",_hitInputCollectionName,std::string("HitCollection"));

	/// Track candidate hits output collection
	registerOutputCollection(LCIO::TRACK,"TrackCandHitOutputCollectionName","Output track candidates hits collection name",_trackCandidateHitsOutputCollectionName,std::string("TrackCandidateHitCollection"));
    ///Link hits in global frame with xy dist under this cut. 
	registerOptionalParameter("DoubletDistCut", "Doublet distance cuts", _doubletDistCut, FloatVec()); 
    ///abs(Pred-Hit(Centre))<cut. 
	registerOptionalParameter("DoubletCenDistCut","Doublet hit acceptance distance from central plane ", _doubletCenDistCut,FloatVec() );
	registerOptionalParameter("localDistDUT", "The local displacement before rotations ", _localDistDUT, FloatVec());
    ///abs(Triplet.Pred-Triplet.Pred)<cut. //Both triplets prediction of where the track would propagate to must be under this cut.
	registerOptionalParameter("TripletConnectDistCut","The distance cut to allow the triplets to be associated with each other", _tripletConnectDistCut,FloatVec() );
    ///abs(Triplet.Slo-Triplet.Slo)<cut. //Both triplets prediction of slope is compared.
	registerOptionalParameter("TripletSlopeCuts", "Triplet slope difference which is allowed to create track ",
	_tripletSlopeCuts, FloatVec());
    registerOptionalParameter("minHits", "Minimum number of hits needed", _minHits ,int(6));
    registerOptionalParameter("mode", "Alignment or basic track fitting. This is either 1 for alignment and 0 for basic track fitting. ", _mode ,int(1));
    registerProcessorParameter("DUTWindow", "The distance the DUT hit must be from the track", _dutDistCut, static_cast<double> (10));

	///This is needed if we have a magnetic field to determine curvature
	registerOptionalParameter("BeamEnergy", "Beam energy [GeV]", _eBeam, static_cast<double> (4.0));

	///This is needed if we have a magnetic field to determine curvature. 
	registerOptionalParameter("BeamCharge", "Beam charge [e]", _qBeam, static_cast<double> (-1));

	/// Histogram information. This is the xml file that specifies the name and bin size etc of the create histograms.
    /// This is not used in this processor at the moment!!!!
	registerOptionalParameter("HistogramInfoFilename", "Name of histogram info xml file", _histoInfoFileName, std::string("histoinfo.xml"));

	/// This is planes that we should not look for hits or create a state. Effectively this removes this plane from the analysis. However scattering due to the material is still taken into account
	registerOptionalParameter("ExcludePlanes", "This is the planes that will not be included in analysis", _excludePlanes ,IntVec());

	/// This initial distance to the first plane is needed if we are in a magnetic field. Since the particle will begin to curve before we reach the first plane. Therefore to accurately model the particles trajectory we need to take this into account. 
	registerOptionalParameter("InitialDisplacement", "This is the initial distance the particle must travel to reach the first plane", _initialDisplacement ,float(0));
    ///This specifies if the planes are strip or pixel sensors.
      registerOptionalParameter("planeDimensions", "This is a number 1(strip sensor) or 2(pixel sensor) to identify the type of detector. Must be in z order and include all planes.", _planeDimension, IntVec());

}
//This is the inital function that Marlin will run only once when we run jobsub
void EUTelProcessorPatRecTriplets::init(){

	try{

	  _DoubletXseperationHistoRight = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("DoubletDist X  (Down stream)", 400, -0.5, 0.5);
	  _DoubletYseperationHistoRight =marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("DoubletDist Y  (Down stream)", 400, -0.5, 0.5);
	  _DoubletXseperationHistoLeft =marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("DoubletDist X  (Up stream)", 400, -0.5, 0.5);
	  _DoubletYseperationHistoLeft =marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("DoubletDist  Y  (Up stream)", 400, -0.5, 0.5);
	  _TripletXseperationHistoRight =marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("DoubletCenDist X , (Down stream)", 400, -0.5, 0.5);
	  _TripletYseperationHistoRight =marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("DoubletCenDist Y , (Down stream)", 400, -0.5, 0.5);
	  _TripletXseperationHistoLeft =marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("DoubletCenDist X (Up stream)", 400, -0.5, 0.5);
	  _TripletYseperationHistoLeft =marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("DoubletCenDist Y (Up stream)", 400, -0.5, 0.5);
	  _TripletDistCutXHisto = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("TripletConnectDist X", 400, -1,1);
	  _TripletDistCutYHisto = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("TripletConnectDist Y", 400, -1,1);
	  _TripletSlopeHistoX = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("TripletSlopeX", 400, -0.01,0.01);
	  _TripletSlopeHistoY = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("TripletSlopeY", 400, -0.01,0.01);
	  _DUTWindowHisto = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("DUTWindow", 500, 0,5);
		_nProcessedRuns = 0;
		_nProcessedEvents = 0;
		std::string name = EUTELESCOPE::GEOFILENAME;
		geo::gGeometry().initializeTGeoDescription(name,false);
		geo::gGeometry().initializeLocalDistDUT(_localDistDUT);
		geo::gGeometry().setInitialDisplacementToFirstPlane(_initialDisplacement); 
		streamlog_out(MESSAGE5) << "These are the planes you will create a state from. Mass inbetween states will be turned to scatterers in GBLTrackProcessor." << std::endl;
//		for(size_t i =0 ; i < geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size(); ++i)
//		{
//			streamlog_out(MESSAGE5) << geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(i) << "  ";
//		}
		streamlog_out(MESSAGE5) << std::endl;
		_trackFitter = new EUTelPatRecTriplets(_DoubletXseperationHistoRight, _DoubletYseperationHistoRight,_DoubletXseperationHistoLeft,
					   _DoubletYseperationHistoLeft, _TripletXseperationHistoRight, _TripletYseperationHistoRight,
					   _TripletXseperationHistoLeft, _TripletYseperationHistoLeft, _TripletDistCutXHisto,
					   _TripletDistCutYHisto,_TripletSlopeHistoX,_TripletSlopeHistoY,_DUTWindowHisto);
        _trackFitter->setMode(_mode);
		_trackFitter->setNumHits(_minHits);
		_trackFitter->setDUTCut(_dutDistCut);

        EUTelExcludedPlanes::setRelativeComplementSet(_excludePlanes);
		_trackFitter->setDoubletDistCut(_doubletDistCut);
		_trackFitter->setTripletSlopeCuts(_tripletSlopeCuts);
		_trackFitter->setDoubletCenDistCut(_doubletCenDistCut);
        _trackFitter->setTripletConnectDistCut(_tripletConnectDistCut);
		_trackFitter->setBeamMomentum(_eBeam);
        _trackFitter->setPlaneDimensionsVec(_planeDimension);
		_trackFitter->testUserInput();
		bookHistograms();		
	}
	catch(std::string &e){
		streamlog_out(MESSAGE9) << e << std::endl;
		throw marlin::StopProcessingException( this ) ;
	}
	catch(lcio::Exception& e){
		streamlog_out(MESSAGE9) << e.what() <<std::endl;
		throw marlin::StopProcessingException( this ) ;

	}
	catch(...){
		throw marlin::StopProcessingException( this ) ;
	}
}

void EUTelProcessorPatRecTriplets::processRunHeader(LCRunHeader* run) {

	std::auto_ptr<EUTelRunHeaderImpl> header(new EUTelRunHeaderImpl(run));
	header->addProcessor(type());//Add what processor has acted to collection here. 

	// this is the right place also to check the geometry ID. This is a
	// unique number identifying each different geometry used at the
	// beam test. The same number should be saved in the run header and
	// in the xml file. If the numbers are different, warn the user.

//	if (header->getGeoID() == 0)
//		streamlog_out(WARNING0) << "The geometry ID in the run header is set to zero." << std::endl
//		<< "This may mean that the GeoID parameter was not set" << std::endl;


//	if ((unsigned int)header->getGeoID() != geo::gGeometry().getSiPlanesLayoutID()) {
//		streamlog_out(WARNING5) << "Error during the geometry consistency check: " << std::endl
//				<< "The run header says the GeoID is " << header->getGeoID() << std::endl
//				<< "The GEAR description says is     " << geo::gGeometry().getSiPlanesLayoutID() << std::endl;
//	}
		_nProcessedRuns++;
}

void EUTelProcessorPatRecTriplets::processEvent(LCEvent* evt)
{
	UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder ( EUTELESCOPE::HITENCODING );

	try{
        /// Must be placed here so we record the event even if we have an exception.
		_nProcessedEvents++;
        /// Change the LCIO object to EUTel object. This is a simple way to extend functionality of the object.
		EUTelEventImpl* event = static_cast<EUTelEventImpl*> (evt); 
		_trackFitter->setEventNumber(_nProcessedEvents);//This is so we can use the event number with this class. 
		/// Do not process last event. For unknown events just a warning will do. 
		if (event->getEventType() == kEORE) {
			streamlog_out(DEBUG4) << "EORE found: nothing else to do." << std::endl;
			return;
		} else if (event->getEventType() == kUNKNOWN) {
			streamlog_out(WARNING2) << "Event number " << event->getEventNumber() << " in run " << event->getRunNumber()<< " is of unknown type. Continue considering it as a normal Data Event." << std::endl;
		}
		LCCollection* hitMeasuredCollection = NULL;
		hitMeasuredCollection = evt->getCollection(_hitInputCollectionName);
		streamlog_out(DEBUG1) << "collection : " <<_hitInputCollectionName << " retrieved" << std::endl;
		if ( hitMeasuredCollection == NULL) {
			streamlog_out(DEBUG2) << "EUTelProcessorPatRecTriplets :: processEvent() hitMeasuredCollection is void" << std::endl;
			throw marlin::SkipEventException(this);
		}
		EVENT::TrackerHitVec allHitsVec;
		for(int iHit = 0; iHit < hitMeasuredCollection->getNumberOfElements(); iHit++) 
		{
			TrackerHitImpl* hit = static_cast<TrackerHitImpl*>(hitMeasuredCollection->getElementAt(iHit));
            int sensorID = hitDecoder(static_cast<IMPL::TrackerHitImpl*>(hit))["sensorID"];
			if ( sensorID >= 0 ) allHitsVec.push_back(hit);
		}
		if(allHitsVec.empty()) throw lcio::Exception("No hits!");
		_trackFitter->setHitsVec(allHitsVec);  
		_trackFitter->printHits();
		streamlog_out( DEBUG1 ) << "Trying to find tracks..." << std::endl;
		std::vector<EUTelTrack> tracks = _trackFitter->getTracks();
		streamlog_out( DEBUG1 ) << "Trying to find tracks...We have " << tracks.size()<<" tracks"<<std::endl;

		plotHistos(tracks);
		outputLCIO(evt,tracks);

	}
	catch (DataNotAvailableException e) {
        ///We expect to get some data exceptions but should make sure we do not get under 5% of events with data.
        _dataMissNumber++;
        const float diff = _nProcessedEvents-_dataMissNumber;
        const float event1 = _nProcessedEvents+1;
        if(diff/event1 < 0.05 and _nProcessedEvents > 100){
            streamlog_out(MESSAGE0) << "The number of events with data is under 5 percent after 100 events"  <<std::endl;
            throw marlin::StopProcessingException( this ) ;
        }
		throw marlin::SkipEventException(this);
	}
	catch(std::string &e){
		streamlog_out(MESSAGE0) << e << std::endl;
		throw marlin::SkipEventException( this ) ;
	}
	catch(lcio::Exception& e){
		streamlog_out(MESSAGE0) << e.what() <<std::endl;
		throw marlin::StopProcessingException( this ) ;
	}
	catch(...){
		streamlog_out(MESSAGE0)<<"Unknown exception in process function of pattern recognition" <<std::endl;
		throw marlin::StopProcessingException( this ) ;
	}
}

void EUTelProcessorPatRecTriplets::outputLCIO(LCEvent* evt, std::vector<EUTelTrack>& tracks)
{
   if(!tracks.empty()){
        for(size_t i=0 ; i< tracks.size(); i++){
            streamlog_out(DEBUG1)<<"Found "<<tracks.size()<<" track for event " << evt->getEventNumber() <<".  Track number  " << i <<std::endl;
            tracks.at(i).print();
        }
        EUTelReaderGenericLCIO reader = EUTelReaderGenericLCIO();
        reader.getColVec(tracks, evt, _trackCandidateHitsOutputCollectionName);
	}else{
        streamlog_out(DEBUG1)<<"tracks.empty() !!!!"<<std::endl;
    }
}

//TO DO: find a more generic way to plot histograms
void EUTelProcessorPatRecTriplets::plotHistos( std::vector<EUTelTrack>& trackCandidates)  {

	const int nTracks = trackCandidates.size( );
 
	static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_numberTracksCandidatesHistName ] ) -> fill( nTracks );
	streamlog_out( MESSAGE2 ) << "Event #" << _nProcessedEvents << std::endl;
	int numberOfHits =0;
	for (size_t i = 0; i< trackCandidates.size( ) ; ++i ) {//loop over all tracks
		for(size_t j = 0; j <trackCandidates[i].getStates().size(); ++j){//loop over all states 
			if(!trackCandidates[i].getStates()[j].getStateHasHit()){//We only ever store on hit per state
				continue;
			}
			numberOfHits++;
			int sensorID = static_cast<int>(trackCandidates[i].getStates().at(j).getLocation());//since we store sensor ID in Z0
			static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_HitOnTrackCandidateHistName ] ) -> fill( sensorID );
			//static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_HitOnTrackTimeHistName ] ) -> fill(trackCandidates[i].getStates()[j].getHit().getTime()  );
		}
		streamlog_out( MESSAGE1 ) << "Track hits end:==============" << std::endl;
		}
}

void EUTelProcessorPatRecTriplets::check(LCEvent * /*evt*/) {
    // nothing to check here
}

void EUTelProcessorPatRecTriplets::end() {
    streamlog_out(MESSAGE9) << "EUTelProcessorPatRecTriplets::end()  " << name()
            << " processed " << _nProcessedEvents << " events in " << _nProcessedRuns << " runs "
            << " av.tracks : " << static_cast<float>(_trackFitter->_numberOfTracksTotal)/static_cast<float>(_nProcessedEvents)
            <<" total number of tracks: " << _trackFitter->_numberOfTracksTotal
            << std::endl;
            delete _trackFitter;

}
//TO DO: Create a better way of booking histograms.
void EUTelProcessorPatRecTriplets::bookHistograms() {
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

    try {
        streamlog_out(DEBUG) << "Booking histograms..." << std::endl;

        std::auto_ptr<EUTelHistogramManager> histoMgr( new EUTelHistogramManager( _histoInfoFileName ));
        EUTelHistogramInfo    * histoInfo;
        bool                    isHistoManagerAvailable;

        try {
            isHistoManagerAvailable = histoMgr->init( );
        } catch ( std::ios::failure& e ) {
            streamlog_out( ERROR5 ) << "I/O problem with " << _histoInfoFileName << "\n"
                    << "Continuing without histogram manager using default settings"    << std::endl;
            isHistoManagerAvailable = false;
        } catch ( marlin::ParseException& e ) {
            streamlog_out( ERROR5 ) << e.what( ) << "\n"
                    << "Continuing without histogram manager using default settings" << std::endl;
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
            streamlog_out(ERROR2) << "Problem booking the " << (_histName::_numberTracksCandidatesHistName) << std::endl;
            streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
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
            streamlog_out(ERROR2) << "Problem booking the " << (_histName::_HitOnTrackCandidateHistName) << std::endl;
            streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
        }

 const int    timeNBin = 40;
        const double timeMin = -10.;
        const double timeMax =10.;
        
        histoInfo = histoMgr->getHistogramInfo(_histName::_HitOnTrackTimeHistName);        
        NBin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : timeNBin;
        XMin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : timeMin;
        XMax = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : timeMax;
        
        // Number of hit per track candidate
        AIDA::IHistogram1D * HitsOnTrackTimes =
                marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_HitOnTrackTimeHistName, NBin,
                XMin, XMax);

        if (HitsOnTrackTimes) {
            HitsOnTrackTimes->setTitle("delta time on hit for track candidates;N hits;N tracks");
            _aidaHistoMap1D.insert(make_pair(_histName::_HitOnTrackTimeHistName, HitsOnTrackTimes));
        } else {
            streamlog_out(ERROR2) << "Problem booking the " << (_histName::_HitOnTrackTimeHistName) << std::endl;
            streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
        }
       
    } catch (lcio::Exception& e) {
        streamlog_out(WARNING2) << "Can't allocate histograms. Continue without histogramming" << std::endl;
    }
#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
}

