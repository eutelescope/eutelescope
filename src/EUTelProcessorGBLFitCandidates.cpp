//Written by Alexander Morton using code by  Denys Lontkovskyi.
//contact alexander.morton@desy.de
#ifdef USE_GBL   
#include "EUTelProcessorGBLFitCandidates.h"
using namespace eutelescope;
//TO DO:
//This way of making histograms makes no sense to me. We should have a class that when called will book any histograms in xml file automatically. So you dont have to book in every processor. It should also return a vector of names to access these histograms. I began this but have not finished. Therefore the silly way of doing the residuals
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelProcessorGBLFitCandidates::_histName::_chi2CandidateHistName = "chi2HistName";
std::string EUTelProcessorGBLFitCandidates::_histName::_fitsuccessHistName = "FitSuccessfulHistName";
std::string EUTelProcessorGBLFitCandidates::_histName::_residGblFitHistNameX0 = "Residual0X";
std::string EUTelProcessorGBLFitCandidates::_histName::_residGblFitHistNameX1 = "Residual1X";
std::string EUTelProcessorGBLFitCandidates::_histName::_residGblFitHistNameX2 = "Residual2X";
std::string EUTelProcessorGBLFitCandidates::_histName::_residGblFitHistNameX3 = "Residual3X";
std::string EUTelProcessorGBLFitCandidates::_histName::_residGblFitHistNameX4 = "Residual4X";
std::string EUTelProcessorGBLFitCandidates::_histName::_residGblFitHistNameX5 = "Residual5X";
std::string EUTelProcessorGBLFitCandidates::_histName::_residGblFitHistNameY0 = "Residual0Y";
std::string EUTelProcessorGBLFitCandidates::_histName::_residGblFitHistNameY1 = "Residual1Y";
std::string EUTelProcessorGBLFitCandidates::_histName::_residGblFitHistNameY2 = "Residual2Y";
std::string EUTelProcessorGBLFitCandidates::_histName::_residGblFitHistNameY3 = "Residual3Y";
std::string EUTelProcessorGBLFitCandidates::_histName::_residGblFitHistNameY4 = "Residual4Y";
std::string EUTelProcessorGBLFitCandidates::_histName::_residGblFitHistNameY5 = "Residual5Y";

#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

EUTelProcessorGBLFitCandidates::EUTelProcessorGBLFitCandidates() :
Processor("EUTelProcessorGBLFitCandidates"),
_trackCandidatesInputCollectionName("Default_input"),
_tracksOutputCollectionName("Default_output"),
_nProcessedRuns(0),
_nProcessedEvents(0),
_beamQ(-1),
_eBeam(4),
_mEstimatorType(), //This is used by the GBL software for outliers down weighting
_maxChi2Cut(1000)
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

  registerOptionalParameter("HistogramInfoFilename", "Name of histogram info xml file", _histoInfoFileName, std::string("histoinfo.xml"));

  registerOptionalParameter("xResolutionPlane", "x resolution of planes given in Planes", _SteeringxResolutions, FloatVec());
  registerOptionalParameter("yResolutionPlane", "y resolution of planes given in Planes", _SteeringyResolutions, FloatVec());

}

void EUTelProcessorGBLFitCandidates::init() {

	streamlog_out(DEBUG2) << "EUTelProcessorGBLFitCandidates::init( )---------------------------------------------BEGIN" << std::endl;
	streamlog_out(DEBUG2) << "Beam charge= " << _beamQ <<" Beam energy= " << _eBeam << std::endl;

	// Reset counters
	_nProcessedRuns = 0;
	_nProcessedEvents = 0;

	// Getting access to geometry description
	std::string name("test.root");
	geo::gGeometry().initializeTGeoDescription(name,false);

	// Initialize GBL fitter. This is the class that does all the work. Seems to me a good practice for the most part create a class that does the work. Since then you can use the same functions in another processor.
	EUTelGBLFitter* Fitter = new EUTelGBLFitter("GBLFitter");
  Fitter->SetBeamCharge(_beamQ);
  Fitter->SetBeamEnergy(_eBeam);
	Fitter->setMEstimatorType(_mEstimatorType);
  Fitter->setChi2Cut(_maxChi2Cut);
  Fitter->setParamterIdXResolutionVec(_SteeringxResolutions);
  Fitter->setParamterIdYResolutionVec(_SteeringyResolutions);
  _trackFitter = Fitter;


	if (!_trackFitter) {
		throw(lcio::Exception(Utility::outputColourString("Could not create instance of fitter class .", "RED")));
	}
	bookHistograms();//TO DO: Remove this and replace with generic histogram class 
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
	streamlog_out(DEBUG5) << "Start of event " << _nProcessedEvents << endl;

	EUTelEventImpl * event = static_cast<EUTelEventImpl*> (evt); ///We change the class so we can use EUTelescope functions

	if (event->getEventType() == kEORE) {
		streamlog_out(DEBUG4) << "EORE found: nothing else to do." << endl;
		return;
	}else if (event->getEventType() == kUNKNOWN) {
		streamlog_out(WARNING2) << "Event number " << event->getEventNumber() << " in run " << event->getRunNumber() << " is of unknown type. Continue considering it as a normal Data Event." << endl;
	}
	LCCollection* col = NULL;
	try {
		col = evt->getCollection(_trackCandidatesInputCollectionName);
		streamlog_out(DEBUG1) << "collection : " << _trackCandidatesInputCollectionName << " retrieved" << std::endl;
	} catch (DataNotAvailableException e) {
		streamlog_out(MESSAGE0) << _trackCandidatesInputCollectionName << " collection not available" << std::endl;
		throw marlin::SkipEventException(this);
	}


	if (col == NULL) {
		streamlog_out(MESSAGE0)<<Utility::outputColourString("The collection is NULL for this event.", "YELLOW")<<std::endl;
		throw marlin::SkipEventException(this);
	}
	std::vector<EUTelTrack> allTracksForThisEvent;
	for (int iCol = 0; iCol < col->getNumberOfElements(); iCol++) {
		//TO DO: The original states from pattern recognition are changed for some reason. We want to create a completely new collection since we may want to run over this data again. 
		//Also this is a memory leak.
		EUTelTrack* trackPointer = new  EUTelTrack(*(static_cast<EUTelTrack*> (col->getElementAt(iCol))));
		EUTelTrack track = *trackPointer;
		_trackFitter->resetPerTrack(); //Here we reset the label that connects state to GBL point to 1 again. Also we set the list of states->labels to 0
		track.print();//Print the track use for debugging
		_trackFitter->testTrack(track);//Check the track has states and hits  
		std::vector< gbl::GblPoint > pointList;
		_trackFitter->setInformationForGBLPointList(track, pointList);
		const gear::BField& B = geo::gGeometry().getMagneticFiled();
		const double Bmag = B.at( TVector3(0.,0.,0.) ).r2();
		gbl::GblTrajectory* traj = 0;
		if ( Bmag < 1.E-6 ) {
			traj = new gbl::GblTrajectory( pointList, false ); 
		}else {
			traj = new gbl::GblTrajectory( pointList, true );
		}
		_trackFitter->setListStateAndLabelAfterTrajectory(pointList,traj);//We must do this after trajectory. Since trajectory will label the points. 
		double  chi2=0; 
		int ndf=0;
		int ierr=0;
		_trackFitter->computeTrajectoryAndFit(pointList,traj, &chi2,&ndf, ierr);
		if(ierr == 0 ){
			 streamlog_out(DEBUG5) << "Ierr is: " << ierr << " Entering loop to update track information " << endl;
			static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_chi2CandidateHistName ] ) -> fill( (chi2)/(ndf));
			static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_fitsuccessHistName ] ) -> fill(1.0);
			track.setChi2(chi2);
			track.setNdf(ndf);
			_trackFitter->UpdateTrackFromGBLTrajectory(traj, pointList,track);
			map< int, map< float, float > >  SensorResidual; 
			_trackFitter->getResidualOfTrackandHits(traj, pointList,track, SensorResidual);
			plotResidual(SensorResidual, _first_time);
			_first_time = false;

		}
		else{
			 streamlog_out(DEBUG5) << "Ierr is: " << ierr << " Do not update track information " << endl;
				static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_fitsuccessHistName ] ) -> fill(0.0);
		}	
	allTracksForThisEvent.push_back(track);
	}//END OF LOOP FOR ALL TRACKS IN AN EVENT
	outputLCIO(evt, allTracksForThisEvent); 
	allTracksForThisEvent.clear();//We clear this so we don't add the same track twice
	streamlog_out(DEBUG5) << "End of event " << _nProcessedEvents << endl;
	_nProcessedEvents++;
}


//TO DO:This is a very stupid way to histogram but will add new class to do this is long run 
void EUTelProcessorGBLFitCandidates::plotResidual(map< int, map<float, float > >  & SensorResidualError, bool &first_time){

		std::map< int, map< float, float > >::iterator sensor_residual_Err_it;
		for(sensor_residual_Err_it = SensorResidualError.begin(); sensor_residual_Err_it != SensorResidualError.end(); sensor_residual_Err_it++) {

				map<float, float> map = sensor_residual_Err_it->second;
				if( !map.empty()){
					float res = map.begin()->first;	
					if( sensor_residual_Err_it->first == 0){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX0 ] ) -> fill(res);}
if( sensor_residual_Err_it->first == 1){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX1 ] ) -> fill(res);}
if( sensor_residual_Err_it->first == 2){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX2 ] ) -> fill(res);}
				if( sensor_residual_Err_it->first == 3){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX3 ] ) -> fill(res);}
				if( sensor_residual_Err_it->first == 4){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX4 ] ) -> fill(res);}
				if( sensor_residual_Err_it->first == 5){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX5 ] ) -> fill(res);}
				}else{
					streamlog_out(DEBUG5) << "The map is NULL" <<std::endl;
				}
									float res2 = map.begin()->second;

					if( sensor_residual_Err_it->first == 0){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameY0 ] ) -> fill(res2);}
if( sensor_residual_Err_it->first == 1){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameY1 ] ) -> fill(res2);}
if( sensor_residual_Err_it->first == 2){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameY2 ] ) -> fill(res2);}
			if( sensor_residual_Err_it->first == 3){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameY3 ] ) -> fill(res2);}
				if( sensor_residual_Err_it->first == 4){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameY4 ] ) -> fill(res2);}
				if( sensor_residual_Err_it->first == 5){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameY5 ] ) -> fill(res2);}
				
        }

}


void EUTelProcessorGBLFitCandidates::end() {}

#endif // USE_GBL


void EUTelProcessorGBLFitCandidates::outputLCIO(LCEvent* evt, std::vector<EUTelTrack>& tracks){

	streamlog_out( DEBUG4 ) << " ---------------- EUTelProcessorGBLFitCandidates::outputLCIO ---------- BEGIN ------------- " << std::endl;

	//Create once per event//Note that this will not cause a memory leak since lcio will delete the memory automatically for you    
	LCCollectionVec * trkCandCollection = new LCCollectionVec(LCIO::TRACK);
	LCCollectionVec * stateCandCollection = new LCCollectionVec(LCIO::TRACK);

	// Prepare output collection
	LCFlagImpl flag(trkCandCollection->getFlag());
	flag.setBit( LCIO::TRBIT_HITS );
	trkCandCollection->setFlag( flag.getFlag( ) );

	LCFlagImpl flag2(stateCandCollection->getFlag());
	flag2.setBit( LCIO::TRBIT_HITS );
	stateCandCollection->setFlag( flag2.getFlag( ) );

	//Loop through all tracks
	for ( int i = 0 ; i < tracks.size(); ++i){
		EUTelTrack* trackheap = new  EUTelTrack();
		//For every track add this to the collection
		for(int j = 0;j < tracks.at(i).getStates().size();++j){
			EUTelState* stateheap = new EUTelState(tracks.at(i).getStates().at(j));
			trackheap->addTrack(stateheap);
			stateCandCollection->push_back(static_cast<EVENT::Track*>(stateheap));
		}
		trkCandCollection->push_back(static_cast<EVENT::Track*>(trackheap));
	}//END TRACK LOOP

	//Now add this collection to the 
	evt->addCollection(trkCandCollection, _tracksOutputCollectionName);
	evt->addCollection(stateCandCollection, "States_After_GBL_Track_Fit");

	streamlog_out( DEBUG4 ) << " ---------------- EUTelProcessorGBLFitCandidates::outputLCIO ---------- END ------------- " << std::endl;
}




void EUTelProcessorGBLFitCandidates::bookHistograms() {
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

 try {
        streamlog_out(DEBUG) << "Booking histograms..." << std::endl;

        auto_ptr<EUTelHistogramManager> histoMgr( new EUTelHistogramManager( _histoInfoFileName ));
        EUTelHistogramInfo    * histoInfo;
        bool                    isHistoManagerAvailable;

////////////////////////////////////////////////////////This is for the residual//Thi si a hack must fix so can accept any number of planes. Really should be a separate processor

        int NBinX=320;
        double MinX=-0.2;
        double MaxX=0.2;


        EUTelHistogramInfo    *    histoInfo0  = histoMgr->getHistogramInfo( _histName::_residGblFitHistNameX0);
        EUTelHistogramInfo    *    histoInfo1  = histoMgr->getHistogramInfo( _histName::_residGblFitHistNameX1);
        EUTelHistogramInfo    *    histoInfo2  = histoMgr->getHistogramInfo( _histName::_residGblFitHistNameX2);
        EUTelHistogramInfo    *    histoInfo3  = histoMgr->getHistogramInfo( _histName::_residGblFitHistNameX3);
        EUTelHistogramInfo    *    histoInfo4  = histoMgr->getHistogramInfo( _histName::_residGblFitHistNameX4);
        EUTelHistogramInfo    *    histoInfo5  = histoMgr->getHistogramInfo( _histName::_residGblFitHistNameX5);
 
        AIDA::IHistogram1D * residGblFit0X = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameX0, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit1X = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameX1, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit2X = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameX2, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit3X = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameX3, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit4X = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameX4, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit5X = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameX5, NBinX, MinX, MaxX); 


              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX0, residGblFit0X));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX1, residGblFit1X));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX2, residGblFit2X));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX3, residGblFit3X));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX4, residGblFit4X));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX5, residGblFit5X));
											
        AIDA::IHistogram1D * residGblFit0Y = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameY0, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit1Y = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameY1, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit2Y = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameY2, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit3Y = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameY3, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit4Y = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameY4, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit5Y = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameY5, NBinX, MinX, MaxX); 


              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameY0, residGblFit0Y));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameY1, residGblFit1Y));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameY2, residGblFit2Y));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameY3, residGblFit3Y));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameY4, residGblFit4Y));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameY5, residGblFit5Y));
											

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////Chi2 create plot.

	const int chiNbins = 5000;
	const double chiXmin = 0;
	const double chiXmax = 100;	
  histoInfo = histoMgr->getHistogramInfo(_histName::_chi2CandidateHistName);

  int NBin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : chiNbins;
  double XMin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin :chiXmin;
  double XMax = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : chiXmax;

  AIDA::IHistogram1D * chi2 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_chi2CandidateHistName, NBin, XMin, XMax);

 	chi2->setTitle("Chi2 of tracks");
  _aidaHistoMap1D.insert(make_pair(_histName::_chi2CandidateHistName, chi2));
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Create plot if fit was successful
	const int Nbins = 2;
	const double Xmin = -0.5;
	const double Xmax = 1.5;	


  histoInfo = histoMgr->getHistogramInfo(_histName::_fitsuccessHistName);

  NBin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : Nbins;
  XMin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin :Xmin;
  XMax = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : Xmax;

  AIDA::IHistogram1D * successful = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_fitsuccessHistName, NBin, XMin, XMax);

 	successful->setTitle("Was the fit successful?");
  _aidaHistoMap1D.insert(make_pair(_histName::_fitsuccessHistName, successful));


}
catch (lcio::Exception& e) {
        streamlog_out(WARNING2) << "Can't allocate histograms. Continue without histogramming" << endl;
}
#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

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



