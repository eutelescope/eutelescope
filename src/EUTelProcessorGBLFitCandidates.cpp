//Written by Alexander Morton using code by  Denys Lontkovskyi.
//contact alexander.morton@desy.de


#ifdef USE_GBL    // Not sure where this is defined. However it is used in all the other GBL processors will use it.

#include "EUTelProcessorGBLFitCandidates.h"

using namespace eutelescope;

//They way we histogram make no sense to me. We should have a class that when called will book any histograms in xml file automatically. So you don not have to book in every processor. It should also return a vector of names to access these histograms. I begain this but have not finished. Therefore the silly way of doing the residuals
// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelProcessorGBLFitCandidates::_histName::_chi2CandidateHistName = "chi2HistName";
std::string EUTelProcessorGBLFitCandidates::_histName::_fitsuccessHistName = "FitSuccessfulHistName";
std::string EUTelProcessorGBLFitCandidates::_histName::_residGblFitHistName0 = "Residual0";
std::string EUTelProcessorGBLFitCandidates::_histName::_residGblFitHistName1 = "Residual1";
std::string EUTelProcessorGBLFitCandidates::_histName::_residGblFitHistName2 = "Residual2";
std::string EUTelProcessorGBLFitCandidates::_histName::_residGblFitHistName3 = "Residual3";
std::string EUTelProcessorGBLFitCandidates::_histName::_residGblFitHistName4 = "Residual4";
std::string EUTelProcessorGBLFitCandidates::_histName::_residGblFitHistName5 = "Residual5";
#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)



EUTelProcessorGBLFitCandidates::EUTelProcessorGBLFitCandidates() :
Processor("EUTelProcessorGBLFitCandidates"),
_trackCandidatesInputCollectionName("Default_input"),
_tracksOutputCollectionName("Default_output"),
_nProcessedRuns(0),
_nProcessedEvents(0),
_beamQ(-1),
_eBeam(4),
_mEstimatorType(), //This is used by the GBL software to outliers down weighting
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

	// Initialize GBL fitter. This is the class that does all the work. Seems to me a good practice to for the most part create a class that does the work. Since then you can use the same functions in another processor.
	EUTelGBLFitter* Fitter = new EUTelGBLFitter("GBLFitter");
  Fitter->SetBeamCharge(_beamQ);
  Fitter->SetBeamEnergy(_eBeam);
	Fitter->setMEstimatorType(_mEstimatorType);
  Fitter->SetChi2Cut(_maxChi2Cut);
  Fitter->setParamterIdXResolutionVec(_SteeringxResolutions);
  Fitter->setParamterIdYResolutionVec(_SteeringyResolutions);
  _trackFitter = Fitter;


	if (!_trackFitter) {
		streamlog_out(ERROR) << "Can't allocate an instance of EUTelGBLFitter. Stopping ..." << std::endl;
		throw UnknownDataTypeException("Track finder was not created");
	}
	bookHistograms();
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
      			streamlog_out(DEBUG1) << "Track " << iCol << " nhits " << EUtrack->getHitsOnTrack().size() << endl;

			//_trackFitter->Clear(); //This is a good point to clear all things that need to be reset for each event. Why should gop here?
			std::vector< gbl::GblPoint > pointList;
      _trackFitter->FillInformationToGBLPointObject(EUtrack, &pointList);

 			const gear::BField& B = geo::gGeometry().getMagneticFiled();
      const double Bmag = B.at( TVector3(0.,0.,0.) ).r2();


  		gbl::GblTrajectory* traj = 0;
      if ( Bmag < 1.E-6 ) {
      	traj = new gbl::GblTrajectory( pointList, false ); //Must make sure this is not a memory leak
      } else {
      	traj = new gbl::GblTrajectory( pointList, true );
      }
			double  chi2=0; 
			int ndf=0;
			int ierr=0;
			_trackFitter->CreateTrajectoryandFit(&pointList,traj, &chi2,&ndf, ierr);
			if(ierr == 0 ){
		     streamlog_out(DEBUG5) << "Ierr is: " << ierr << " Entering loop to update track information " << endl;
      	static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_chi2CandidateHistName ] ) -> fill( (chi2)/(ndf));
      	static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_fitsuccessHistName ] ) -> fill(1.0);
				//Update track and then state variables//////////////////////////////////////////////BEGIN
				EUtrack->setChi2(chi2);
				EUtrack->setNdf(ndf);
				_trackFitter->UpdateTrackFromGBLTrajectory(traj, &pointList);
				EUtracks.push_back(EUtrack); 
				//////////////////////////////////////////////////////////////////////////////////////END
				///////////////////////////////////////////////////////////////////////////////////////////Determine the residual of this new track and plot// BEGIN
				map< int, map< float, float > >  SensorResidualError; 
				_trackFitter->getResidualOfTrackandHits(traj, &pointList, SensorResidualError);
				plotResidual(SensorResidualError, _first_time);
				_first_time = false;
 
			}
			else{
		     streamlog_out(DEBUG5) << "Ierr is: " << ierr << " Do not update track information " << endl;
      		static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_fitsuccessHistName ] ) -> fill(0.0);
			}	
	_trackFitter->Clear();		
			      
		}//END OF LOOP FOR ALL TRACKS IN AN EVENT
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			outputLCIO(evt, EUtracks); // Important to note that EUtracks are still only the original states at input. The scatterers are used in the fit but are not included here.
	}//END OF COLLECTION IS NOT NULL LOOP	
_nProcessedEvents++;

}



void EUTelProcessorGBLFitCandidates::plotResidual(map< int, map<float, float > >  & SensorResidualError, bool &first_time){

		std::map< int, map< float, float > >::iterator sensor_residual_Err_it;
		for(sensor_residual_Err_it = SensorResidualError.begin(); sensor_residual_Err_it != SensorResidualError.end(); sensor_residual_Err_it++) {

				map<float, float> map = sensor_residual_Err_it->second;
				if( !map.empty()){
					float res = map.begin()->first;	
					if( sensor_residual_Err_it->first == 0){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistName0 ] ) -> fill(res);}
if( sensor_residual_Err_it->first == 1){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistName1 ] ) -> fill(res);}
if( sensor_residual_Err_it->first == 2){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistName2 ] ) -> fill(res);}
				if( sensor_residual_Err_it->first == 3){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistName3 ] ) -> fill(res);}
				if( sensor_residual_Err_it->first == 4){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistName4 ] ) -> fill(res);}
				if( sensor_residual_Err_it->first == 5){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistName5 ] ) -> fill(res);}
				}else{
					streamlog_out(DEBUG5) << "The map is NULL" <<std::endl;
				}
				
        }

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


        EUTelHistogramInfo    *    histoInfo0  = histoMgr->getHistogramInfo( _histName::_residGblFitHistName0);
        EUTelHistogramInfo    *    histoInfo1  = histoMgr->getHistogramInfo( _histName::_residGblFitHistName1);
        EUTelHistogramInfo    *    histoInfo2  = histoMgr->getHistogramInfo( _histName::_residGblFitHistName2);
        EUTelHistogramInfo    *    histoInfo3  = histoMgr->getHistogramInfo( _histName::_residGblFitHistName3);
        EUTelHistogramInfo    *    histoInfo4  = histoMgr->getHistogramInfo( _histName::_residGblFitHistName4);
        EUTelHistogramInfo    *    histoInfo5  = histoMgr->getHistogramInfo( _histName::_residGblFitHistName5);
 
        AIDA::IHistogram1D * residGblFit0 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistName0, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit1 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistName1, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit2 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistName2, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit3 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistName3, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit4 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistName4, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit5 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistName5, NBinX, MinX, MaxX); 


              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistName0, residGblFit0));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistName1, residGblFit1));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistName2, residGblFit2));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistName3, residGblFit3));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistName4, residGblFit4));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistName5, residGblFit5));
											
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



