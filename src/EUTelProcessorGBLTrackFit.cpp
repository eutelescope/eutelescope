//Written by Alexander Morton, Denys Lontkovskyi and Igor Rubinskiy.
//contact alexander.morton975@gmail.com
#ifdef USE_GBL   
#include "EUTelProcessorGBLTrackFit.h"
using namespace eutelescope;
//TO DO:
//This way of making histograms makes no sense to me. We should have a class that when called will book any histograms in xml file automatically. So you dont have to book in every processor. It should also return a vector of names to access these histograms. I began this but have not finished. Therefore the silly way of doing the residuals
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelProcessorGBLTrackFit::_histName::_chi2CandidateHistName = "chi2HistName";
std::string EUTelProcessorGBLTrackFit::_histName::_fitsuccessHistName = "FitSuccessfulHistName";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameX0 = "Residual0X";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameX1 = "Residual1X";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameX2 = "Residual2X";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameX3 = "Residual3X";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameX4 = "Residual4X";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameX5 = "Residual5X";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameY0 = "Residual0Y";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameY1 = "Residual1Y";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameY2 = "Residual2Y";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameY3 = "Residual3Y";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameY4 = "Residual4Y";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameY5 = "Residual5Y";

std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameX0p = "Residual0Xpull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameX1p = "Residual1Xpull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameX2p = "Residual2Xpull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameX3p = "Residual3Xpull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameX4p = "Residual4Xpull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameX5p = "Residual5Xpull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameY0p = "Residual0Ypull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameY1p = "Residual1Ypull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameY2p = "Residual2Ypull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameY3p = "Residual3Ypull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameY4p = "Residual4Ypull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameY5p = "Residual5Ypull";

#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

EUTelProcessorGBLTrackFit::EUTelProcessorGBLTrackFit() :
Processor("EUTelProcessorGBLTrackFit"),
_trackCandidatesInputCollectionName("Default_input"),
_tracksOutputCollectionName("Default_output"),
_nProcessedRuns(0),
_nProcessedEvents(0),
_beamQ(-1),
_eBeam(4),
_mEstimatorType(), //This is used by the GBL software for outliers down weighting
_first_time(true)
{
	// Processor description
	_description = "EUTelProcessorGBLTrackFit this will fit gbl tracks and output them into LCIO file.";
  // TrackerHit input collection
  registerInputCollection(LCIO::TRACK, "TrackCandidatesInputCollectionName", "Input track candidate collection name",_trackCandidatesInputCollectionName,std::string("TrackCandidatesCollection"));
  // Track output collection
  registerOutputCollection(LCIO::TRACK,"TracksOutputCollectionName","Output tracks collection name",_tracksOutputCollectionName,std::string("TrackCollection"));

  registerOptionalParameter("BeamCharge", "Beam charge [e]", _beamQ, static_cast<double> (-1));

  registerProcessorParameter("BeamEnergy", "Beam energy [GeV]", _eBeam, static_cast<double> (4.0));
	//This is the determines the how we down weight our outliers. This by default is set that each point will have the same weighting.
  registerOptionalParameter("GBLMEstimatorType", "GBL outlier down-weighting option (t,h,c)", _mEstimatorType, string() );

  registerOptionalParameter("HistogramInfoFilename", "Name of histogram info xml file", _histoInfoFileName, std::string("histoinfo.xml"));
	//This is the estimated resolution of the planes and DUT in x/y direction
  registerOptionalParameter("xResolutionPlane", "x resolution of planes given in Planes", _SteeringxResolutions, FloatVec());
  registerOptionalParameter("yResolutionPlane", "y resolution of planes given in Planes", _SteeringyResolutions, FloatVec());
}

void EUTelProcessorGBLTrackFit::init() {
	try{
		streamlog_out(DEBUG2) << "EUTelProcessorGBLTrackFit::init( )---------------------------------------------BEGIN" << std::endl;
		streamlog_out(DEBUG2) << "Beam charge= " << _beamQ <<" Beam energy= " << _eBeam << std::endl;
		// Reset counters
		_nProcessedRuns = 0;
		_nProcessedEvents = 0;
		//Create TGeo description from the gear.
		std::string name("test.root");
		geo::gGeometry().initializeTGeoDescription(name,false);
		// Initialize GBL fitter. This is the class that does all the work. Seems to me a good practice for the most part create a class that does the work. Since then you can use the same functions in another processor.
		EUTelGBLFitter* Fitter = new EUTelGBLFitter();
		Fitter->setBeamCharge(_beamQ);
		Fitter->setBeamEnergy(_eBeam);
		Fitter->setMEstimatorType(_mEstimatorType);//As said before this is to do with how we deal with outliers and the function we use to weight them.
		Fitter->setParamterIdXResolutionVec(_SteeringxResolutions);
		Fitter->setParamterIdYResolutionVec(_SteeringyResolutions);
		Fitter->testUserInput();
		_trackFitter = Fitter;
		if (!_trackFitter) {
			throw(lcio::Exception(Utility::outputColourString("Could not create instance of fitter class .", "RED")));
		}
		bookHistograms();//TO DO: Remove this and replace with generic histogram class 
		streamlog_out(DEBUG2) << "EUTelProcessorGBLTrackFit::init( )---------------------------------------------END" << std::endl;
	}	
	catch(std::string &e){
		streamlog_out(MESSAGE9) << e << std::endl;
		throw StopProcessingException( this ) ;
	}
	catch(lcio::Exception& e){
		streamlog_out(MESSAGE9) << e.what() <<endl;
		throw StopProcessingException( this ) ;
	}
	catch(...){
		streamlog_out(MESSAGE9)<<Utility::outputColourString("Unknown exception in init function of EUTelProcessorGBLTrackFit.","RED") <<endl;
		throw StopProcessingException( this ) ;
	}
}

void EUTelProcessorGBLTrackFit::processRunHeader(LCRunHeader * run) {
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
	_chi2NdfVec.clear();//TO DO:This is needed to determine if the track is near chi2 of one. Do we need this however?
    
	_nProcessedRuns++;
}

void EUTelProcessorGBLTrackFit::check(LCEvent * evt){}

void EUTelProcessorGBLTrackFit::processEvent(LCEvent * evt){
	try{
		streamlog_out(DEBUG5) << "Start of event " << _nProcessedEvents << endl;

		EUTelEventImpl * event = static_cast<EUTelEventImpl*> (evt); ///We change the class so we can use EUTelescope functions

		if (event->getEventType() == kEORE) {
			streamlog_out(DEBUG4) << "EORE found: nothing else to do." << endl;
			return;
		}else if (event->getEventType() == kUNKNOWN) {
			streamlog_out(WARNING2) << "Event number " << event->getEventNumber() << " in run " << event->getRunNumber() << " is of unknown type. Continue considering it as a normal Data Event." << endl;
		}
		LCCollection* col = NULL;
		col = evt->getCollection(_trackCandidatesInputCollectionName);
		streamlog_out(DEBUG1) << "collection : " << _trackCandidatesInputCollectionName << " retrieved" << std::endl;

		if (col == NULL) {
			streamlog_out(MESSAGE0)<<Utility::outputColourString("The collection is NULL for this event.", "YELLOW")<<std::endl;
			throw marlin::SkipEventException(this);
		}
		std::vector<EUTelTrack> allTracksForThisEvent;//GBL will analysis the track one at a time. However we want to save to lcio per event.
		for (int iCol = 0; iCol < col->getNumberOfElements(); iCol++) {
			//TO DO: Make sure that the original input collection is not changed by this.
			EUTelTrack track =  EUTelTrack(*(static_cast<EUTelTrack*> (col->getElementAt(iCol))));//TO DO: Is there a more elegant way to copy this?
			_trackFitter->resetPerTrack(); //Here we reset the label that connects state to GBL point to 1 again. Also we set the list of states->labels to 0
			track.print();//Print the track use for debugging
			_trackFitter->testTrack(track);//Check the track has states and hits  
			std::vector< gbl::GblPoint > pointList;
			_trackFitter->setInformationForGBLPointList(track, pointList);//Here we describe the whole setup. Geometry, scattering, data...
			const gear::BField& B = geo::gGeometry().getMagneticFiled();//We need this to determine if we should fit a curve or a straight line.
			const double Bmag = B.at( TVector3(0.,0.,0.) ).r2();
			_trackFitter->setPairMeasurementStateAndPointLabelVec(pointList);//This will create a link between the states that have a hit associated with them and the GBL label that is associated with the state.
			gbl::GblTrajectory* traj = 0;
			//Here we create the trajectory from the points created by setInformationForGBLPointList. This will take the points and propagation jacobian and split this into smaller matrices to describe the problem in terms of offsets. Here is the difference between GBL and other fitting algorithms.  
			if ( Bmag < 1.E-6 ) {
				traj = new gbl::GblTrajectory( pointList, false ); 
			}else {
				traj = new gbl::GblTrajectory( pointList, true );
			}
			_trackFitter->setPairAnyStateAndPointLabelVec(pointList,traj);//This will create a link between any state and it's GBL point label. 
			double  chi2=0; 
			int ndf=0;
			int ierr=0;
			_trackFitter->computeTrajectoryAndFit(pointList,traj, &chi2,&ndf, ierr);//This will do the minimisation of the chi2 and produce the most probable trajectory.
			if(ierr == 0 ){
				 streamlog_out(DEBUG5) << "Ierr is: " << ierr << " Entering loop to update track information " << endl;
				static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_chi2CandidateHistName ] ) -> fill( (chi2)/(ndf));
				static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_fitsuccessHistName ] ) -> fill(1.0);
				if(chi2 ==0 or ndf ==0){
					throw(lcio::Exception(Utility::outputColourString("Your fitted track has zero degrees of freedom or a chi2 of 0.", "RED"))); 	
				}

				track.setChi2(chi2);
				track.setNdf(ndf);
				_chi2NdfVec.push_back(chi2/static_cast<float>(ndf));
				map<int, vector<double> >  mapSensorIDToCorrectionVec;//This is not used now. However it maybe useful to be able to access the corrections that GBL makes to the original track. Since if this is too large then GBL may give th wrong trajectory. Since all the equations are only to first order. 
				_trackFitter->updateTrackFromGBLTrajectory(traj, pointList,track,mapSensorIDToCorrectionVec);
				map< int, map< float, float > >  SensorResidual; 
				map< int, map< float, float > >  SensorResidualError; 
				_trackFitter->getResidualOfTrackandHits(traj, pointList,track, SensorResidual, SensorResidualError);
				plotResidual(SensorResidual,SensorResidualError, _first_time);//TO DO: Need to fix how we histogram.
				_first_time = false;

			}
			else{
				 streamlog_out(DEBUG5) << "Ierr is: " << ierr << " Do not update track information " << endl;
					static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_fitsuccessHistName ] ) -> fill(0.0);
					continue;//We continue so we don't add an empty track
			}	
		allTracksForThisEvent.push_back(track);
		}//END OF LOOP FOR ALL TRACKS IN AN EVENT
		outputLCIO(evt, allTracksForThisEvent); 
		allTracksForThisEvent.clear();//We clear this so we don't add the same track twice
		streamlog_out(DEBUG5) << "End of event " << _nProcessedEvents << endl;
		_nProcessedEvents++;
	}
	catch (DataNotAvailableException e) {
		streamlog_out(MESSAGE0) << _trackCandidatesInputCollectionName << " collection not available" << std::endl;
		throw marlin::SkipEventException(this);
	}
	catch(std::string &e){
		streamlog_out(MESSAGE9) << e << std::endl;
		throw StopProcessingException( this ) ;
	}
	catch(lcio::Exception& e){
		streamlog_out(MESSAGE9) << e.what() <<endl;
		throw StopProcessingException( this ) ;
	}
	catch(...){
		streamlog_out(MESSAGE9)<<"Unknown exception in processEvent function of EUTelProcessorGBLTrackFit" <<endl;
		throw StopProcessingException( this ) ;
	}
	
}


//TO DO:This is a very stupid way to histogram but will add new class to do this is long run 
void EUTelProcessorGBLTrackFit::plotResidual(map< int, map<float, float > >  & sensorResidual, map< int, map<float, float > >  & sensorResidualError, bool &first_time){
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Residual plot
	std::map< int, map< float, float > >::iterator sensor_residual_it;
	for(sensor_residual_it = sensorResidual.begin(); sensor_residual_it != sensorResidual.end(); sensor_residual_it++) {
		map<float, float> map = sensor_residual_it->second;
		if( !map.empty()){
			float res = map.begin()->first;	
			if( sensor_residual_it->first == 0){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX0 ] ) -> fill(res);}
			if( sensor_residual_it->first == 1){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX1 ] ) -> fill(res);}
			if( sensor_residual_it->first == 2){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX2 ] ) -> fill(res);}
			if( sensor_residual_it->first == 3){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX3 ] ) -> fill(res);}
			if( sensor_residual_it->first == 4){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX4 ] ) -> fill(res);}
			if( sensor_residual_it->first == 5){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX5 ] ) -> fill(res);}
			if( sensor_residual_it->first == 20){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ "Residual6" ] ) -> fill(res);}
		}else{
			streamlog_out(DEBUG5) << "The map is NULL" <<std::endl;
		}
			float res2 = map.begin()->second;

			if( sensor_residual_it->first == 0){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameY0 ] ) -> fill(res2);}
			if( sensor_residual_it->first == 1){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameY1 ] ) -> fill(res2);}
			if( sensor_residual_it->first == 2){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameY2 ] ) -> fill(res2);}
			if( sensor_residual_it->first == 3){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameY3 ] ) -> fill(res2);}
			if( sensor_residual_it->first == 4){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameY4 ] ) -> fill(res2);}
			if( sensor_residual_it->first == 5){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameY5 ] ) -> fill(res2);}
				
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Residual Error plot
	std::map< int, map< float, float > >::iterator sensor_residualerror_it;
	for(sensor_residualerror_it = sensorResidualError.begin(); sensor_residualerror_it != sensorResidualError.end(); sensor_residualerror_it++) {
			map<float, float> maperror = sensor_residualerror_it->second;
			map<float, float> mapres = sensorResidual.at(sensor_residualerror_it->first);

			if( !maperror.empty()){
				float res = mapres.begin()->first;	
				float reserror = maperror.begin()->first;	
				if( sensor_residualerror_it->first == 0){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX0p ] ) -> fill(res/reserror);}
				if( sensor_residualerror_it->first == 1){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX1p ] ) -> fill(res/reserror);}
				if( sensor_residualerror_it->first == 2){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX2p ] ) -> fill(res/reserror);}
				if( sensor_residualerror_it->first == 3){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX3p ] ) -> fill(res/reserror);}
				if( sensor_residualerror_it->first == 4){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX4p ] ) -> fill(res/reserror);}
				if( sensor_residualerror_it->first == 5){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX5p ] ) -> fill(res/reserror);}
				}else{
					streamlog_out(DEBUG5) << "The map is NULL" <<std::endl;
				}
				float res2 = mapres.begin()->second;
				float res2error = maperror.begin()->second;

				if( sensor_residualerror_it->first == 0){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameY0p ] ) -> fill(res2/res2error);}
				if( sensor_residualerror_it->first == 1){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameY1p ] ) -> fill(res2/res2error);}
				if( sensor_residualerror_it->first == 2){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameY2p ] ) -> fill(res2/res2error);}
				if( sensor_residualerror_it->first == 3){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameY3p ] ) -> fill(res2/res2error);}
				if( sensor_residualerror_it->first == 4){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameY4p ] ) -> fill(res2/res2error);}
				if( sensor_residualerror_it->first == 5){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameY5p ] ) -> fill(res2/res2error);}
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}


void EUTelProcessorGBLTrackFit::end() {
	float total = 0;
	double sizeFittedTracks = _chi2NdfVec.size();
	cout<<"This is the chi2s"<<endl;
	for(int i=0; i<_chi2NdfVec.size(); ++i){
		total= total + _chi2NdfVec.at(i);//TO DO: This is does not seem to output the correct average chi2. Plus do we really need this to fit?
		cout<<_chi2NdfVec.at(i)<<endl;
	}
	//TO DO: We really should have a better way to look track per track	and see if the correction is too large. 
	std::vector<double> correctionTotal = _trackFitter->getCorrectionsTotal();
	streamlog_out(MESSAGE9)<<"This is the average correction for omega: " <<correctionTotal.at(0)/sizeFittedTracks<<endl;	
	streamlog_out(MESSAGE9)<<"This is the average correction for local xz inclination: " <<correctionTotal.at(1)/sizeFittedTracks<<endl;	
	streamlog_out(MESSAGE9)<<"This is the average correction for local yz inclination: " <<correctionTotal.at(2)/sizeFittedTracks<<endl;	
	streamlog_out(MESSAGE9)<<"This is the average correction for local x: " <<correctionTotal.at(3)/sizeFittedTracks<<endl;	
	streamlog_out(MESSAGE9)<<"This is the average correction for local y: " <<correctionTotal.at(4)/sizeFittedTracks<<endl;	

  float average = total/sizeFittedTracks;
	streamlog_out(MESSAGE9) << "This is the average chi2 -"<< average <<std::endl;

}

#endif // USE_GBL

//TO DO: If you have a missing hit from a state this does not seem to work. However you can exclude planes which produces no hit which is fine
void EUTelProcessorGBLTrackFit::outputLCIO(LCEvent* evt, std::vector<EUTelTrack>& tracks){

	streamlog_out( DEBUG4 ) << " ---------------- EUTelProcessorGBLTrackFit::outputLCIO ---------- BEGIN ------------- " << std::endl;

	//Create once per event//Note that this will not cause a memory leak since lcio will delete the memory automatically for you    
	//Should create a new track and state in lcio memory. So we don't change the old stuff.
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
		EUTelTrack* trackheap = new  EUTelTrack(tracks.at(i), false); //We dont want to copy contents so set to false. We only want the track object but not the states or hits. Since the states should have there own place in memory before saving.
		trackheap->print();
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
	string name = _tracksOutputCollectionName + "_GBLstates" ;
	evt->addCollection(stateCandCollection, name);

	streamlog_out( DEBUG4 ) << " ---------------- EUTelProcessorGBLTrackFit::outputLCIO ---------- END ------------- " << std::endl;
}

void EUTelProcessorGBLTrackFit::bookHistograms() {
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
 try {
        streamlog_out(DEBUG) << "Booking histograms..." << std::endl;

        auto_ptr<EUTelHistogramManager> histoMgr( new EUTelHistogramManager( _histoInfoFileName ));
        EUTelHistogramInfo    * histoInfo;
        bool                    isHistoManagerAvailable;

////////////////////////////////////////////////////////This is for the residual//Thi si a hack must fix so can accept any number of planes. Really should be a separate processor

        int NBinX=300;
        double MinX=-0.04;  //-0.2;
        double MaxX=0.04;


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

        AIDA::IHistogram1D * residGblFit6X = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Residual6" , NBinX, MinX, MaxX); 


              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX0, residGblFit0X));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX1, residGblFit1X));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX2, residGblFit2X));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX3, residGblFit3X));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX4, residGblFit4X));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX5, residGblFit5X));
              _aidaHistoMap1D.insert(std::make_pair("Residual6", residGblFit6X));

											
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

			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////??Corrections		
/*				NBinX=300;
        MinX=-0.0001;  //-0.2;
        MaxX=0.0001;

        AIDA::IHistogram1D * correction0Plane0 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction0 Plane0", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction1Plane0 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction1 Plane0", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction2Plane0 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction2 Plane0", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction3Plane0 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction3 Plane0", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction4Plane0 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction4 Plane0", NBinX, MinX, MaxX); 

				_mapSensorIDToHistogramCorrection0.insert(std::make_pair(0,correction0Plane0));
				_mapSensorIDToHistogramCorrection1.insert(std::make_pair(0,correction1Plane0));
				_mapSensorIDToHistogramCorrection2.insert(std::make_pair(0,correction2Plane0));
				_mapSensorIDToHistogramCorrection3.insert(std::make_pair(0,correction3Plane0));
				_mapSensorIDToHistogramCorrection4.insert(std::make_pair(0,correction4Plane0));

        AIDA::IHistogram1D * correction0Plane1 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction0 Plane1", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction1Plane1 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction1 Plane1", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction2Plane1 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction2 Plane1", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction3Plane1 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction3 Plane1", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction4Plane1 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction4 Plane1", NBinX, MinX, MaxX); 

				_mapSensorIDToHistogramCorrection0.insert(std::make_pair(1,correction0Plane1));
				_mapSensorIDToHistogramCorrection1.insert(std::make_pair(1,correction1Plane1));
				_mapSensorIDToHistogramCorrection2.insert(std::make_pair(1,correction2Plane1));
				_mapSensorIDToHistogramCorrection3.insert(std::make_pair(1,correction3Plane1));
				_mapSensorIDToHistogramCorrection4.insert(std::make_pair(1,correction4Plane1));

        AIDA::IHistogram1D * correction0Plane2 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction0 Plane2", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction1Plane2 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction1 Plane2", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction2Plane2 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction2 Plane2", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction3Plane2 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction3 Plane2", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction4Plane2 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction4 Plane2", NBinX, MinX, MaxX); 

				_mapSensorIDToHistogramCorrection0.insert(std::make_pair(2,correction0Plane2));
				_mapSensorIDToHistogramCorrection1.insert(std::make_pair(2,correction1Plane2));
				_mapSensorIDToHistogramCorrection2.insert(std::make_pair(2,correction2Plane2));
				_mapSensorIDToHistogramCorrection3.insert(std::make_pair(2,correction3Plane2));
				_mapSensorIDToHistogramCorrection4.insert(std::make_pair(2,correction4Plane2));

        AIDA::IHistogram1D * correction0Plane3 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction0 Plane3", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction1Plane3 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction1 Plane3", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction2Plane3 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction2 Plane3", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction3Plane3 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction3 Plane3", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction4Plane3 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction4 Plane3", NBinX, MinX, MaxX); 

				_mapSensorIDToHistogramCorrection0.insert(std::make_pair(3,correction0Plane3));
				_mapSensorIDToHistogramCorrection1.insert(std::make_pair(3,correction1Plane3));
				_mapSensorIDToHistogramCorrection2.insert(std::make_pair(3,correction2Plane3));
				_mapSensorIDToHistogramCorrection3.insert(std::make_pair(3,correction3Plane3));
				_mapSensorIDToHistogramCorrection4.insert(std::make_pair(3,correction4Plane3));

        AIDA::IHistogram1D * correction0Plane4 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction0 Plane4", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction1Plane4 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction1 Plane4", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction2Plane4 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction2 Plane4", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction3Plane4 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction3 Plane4", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction4Plane4 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction4 Plane4", NBinX, MinX, MaxX); 

				_mapSensorIDToHistogramCorrection0.insert(std::make_pair(4,correction0Plane4));
				_mapSensorIDToHistogramCorrection1.insert(std::make_pair(4,correction1Plane4));
				_mapSensorIDToHistogramCorrection2.insert(std::make_pair(4,correction2Plane4));
				_mapSensorIDToHistogramCorrection3.insert(std::make_pair(4,correction3Plane4));
				_mapSensorIDToHistogramCorrection4.insert(std::make_pair(4,correction4Plane4));

        AIDA::IHistogram1D * correction0Plane5 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction0 Plane5", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction1Plane5 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction1 Plane5", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction2Plane5 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction2 Plane5", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction3Plane5 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction3 Plane5", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction4Plane5 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction4 Plane5", NBinX, MinX, MaxX); 

				_mapSensorIDToHistogramCorrection0.insert(std::make_pair(5,correction0Plane5));
				_mapSensorIDToHistogramCorrection1.insert(std::make_pair(5,correction1Plane5));
				_mapSensorIDToHistogramCorrection2.insert(std::make_pair(5,correction2Plane5));
				_mapSensorIDToHistogramCorrection3.insert(std::make_pair(5,correction3Plane5));
				_mapSensorIDToHistogramCorrection4.insert(std::make_pair(5,correction4Plane5));

       AIDA::IHistogram1D * correction0Plane20 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction0 Plane20", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction1Plane20 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction1 Plane20", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction2Plane20 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction2 Plane20", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction3Plane20 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction3 Plane20", NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * correction4Plane20 = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("Correction4 Plane20", NBinX, MinX, MaxX); 

				_mapSensorIDToHistogramCorrection0.insert(std::make_pair(20,correction0Plane20));
				_mapSensorIDToHistogramCorrection1.insert(std::make_pair(20,correction1Plane20));
				_mapSensorIDToHistogramCorrection2.insert(std::make_pair(20,correction2Plane20));
				_mapSensorIDToHistogramCorrection3.insert(std::make_pair(20,correction3Plane20));
				_mapSensorIDToHistogramCorrection4.insert(std::make_pair(20,correction4Plane20));
*/

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////The average residual plots per plane
//			AIDA::ICloud2D* plotAverageResidual = 	marlin::AIDAProcessor::histogramFactory(this)->createCloud2D("aveRes", "The Average residual with plane number"); 				
	//		plotAverageResidual->fill(1,10,1);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Pull plots											
        int NBinp=120;
        double MinXp=-10;
        double MaxXp=10;


        EUTelHistogramInfo    *    histoInfo0p  = histoMgr->getHistogramInfo( _histName::_residGblFitHistNameX0p);
        EUTelHistogramInfo    *    histoInfo1p  = histoMgr->getHistogramInfo( _histName::_residGblFitHistNameX1p);
        EUTelHistogramInfo    *    histoInfo2p  = histoMgr->getHistogramInfo( _histName::_residGblFitHistNameX2p);
        EUTelHistogramInfo    *    histoInfo3p  = histoMgr->getHistogramInfo( _histName::_residGblFitHistNameX3p);
        EUTelHistogramInfo    *    histoInfo4p  = histoMgr->getHistogramInfo( _histName::_residGblFitHistNameX4p);
        EUTelHistogramInfo    *    histoInfo5p  = histoMgr->getHistogramInfo( _histName::_residGblFitHistNameX5p);
 
        AIDA::IHistogram1D * residGblFit0Xp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameX0p, NBinp, MinXp, MaxXp); 
        AIDA::IHistogram1D * residGblFit1Xp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameX1p, NBinp, MinXp, MaxXp); 
        AIDA::IHistogram1D * residGblFit2Xp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameX2p, NBinp, MinXp, MaxXp); 
        AIDA::IHistogram1D * residGblFit3Xp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameX3p, NBinp, MinXp, MaxXp); 
        AIDA::IHistogram1D * residGblFit4Xp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameX4p, NBinp, MinXp, MaxXp); 
        AIDA::IHistogram1D * residGblFit5Xp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameX5p, NBinp, MinXp, MaxXp); 


              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX0p, residGblFit0Xp));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX1p, residGblFit1Xp));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX2p, residGblFit2Xp));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX3p, residGblFit3Xp));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX4p, residGblFit4Xp));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX5p, residGblFit5Xp));
											
        AIDA::IHistogram1D * residGblFit0Yp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameY0p, NBinp, MinXp, MaxXp); 
        AIDA::IHistogram1D * residGblFit1Yp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameY1p, NBinp, MinXp, MaxXp); 
        AIDA::IHistogram1D * residGblFit2Yp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameY2p, NBinp, MinXp, MaxXp); 
        AIDA::IHistogram1D * residGblFit3Yp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameY3p, NBinp, MinXp, MaxXp); 
        AIDA::IHistogram1D * residGblFit4Yp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameY4p, NBinp, MinXp, MaxXp); 
        AIDA::IHistogram1D * residGblFit5Yp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameY5p, NBinp, MinXp, MaxXp); 


              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameY0p, residGblFit0Yp));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameY1p, residGblFit1Yp));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameY2p, residGblFit2Yp));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameY3p, residGblFit3Yp));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameY4p, residGblFit4Yp));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameY5p, residGblFit5Yp));
											



	///////////////////////////////////////////////////////////////////////////////////////////////////////////////Chi2 create plot.

	const int chiNbins = 100;
	const double chiXmin = 0;
	const double chiXmax = 10;	
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

