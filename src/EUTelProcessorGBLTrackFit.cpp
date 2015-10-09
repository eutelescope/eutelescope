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
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameXDut1 = "ResidualDut1X";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameXDut2 = "ResidualDut2X";

std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameY0 = "Residual0Y";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameY1 = "Residual1Y";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameY2 = "Residual2Y";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameY3 = "Residual3Y";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameY4 = "Residual4Y";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameY5 = "Residual5Y";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameYDut1 = "ResidualDut1Y";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameYDut2 = "ResidualDut2Y";

std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameX0p = "Residual0Xpull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameX1p = "Residual1Xpull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameX2p = "Residual2Xpull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameX3p = "Residual3Xpull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameX4p = "Residual4Xpull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameX5p = "Residual5Xpull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameXDut1p = "ResidualDut1Xpull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameXDut2p = "ResidualDut2Xpull";

std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameY0p = "Residual0Ypull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameY1p = "Residual1Ypull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameY2p = "Residual2Ypull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameY3p = "Residual3Ypull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameY4p = "Residual4Ypull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameY5p = "Residual5Ypull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameYDut1p = "ResidualDut1Ypull";
std::string EUTelProcessorGBLTrackFit::_histName::_residGblFitHistNameYDut2p = "ResidualDut2Ypull";

#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

EUTelProcessorGBLTrackFit::EUTelProcessorGBLTrackFit() :
Processor("EUTelProcessorGBLTrackFit"),
_nProcessedRuns(0),
_nProcessedEvents(0),
_beamQ(-1),
_eBeam(4),
_trackCandidatesInputCollectionName("Default_input"),
_tracksOutputCollectionName("Default_output"),
_mEstimatorType() //This is used by the GBL software for outliers down weighting
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
  registerOptionalParameter("GBLMEstimatorType", "GBL outlier down-weighting option (t,h,c)", _mEstimatorType, std::string() );

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
			throw(lcio::Exception("Could not create instance of fitter class."));
		}
		//Create millepede output
//		_Mille  = new EUTelMillepede(); 

		bookHistograms();//TO DO: Remove this and replace with generic histogram class 
		streamlog_out(DEBUG2) << "EUTelProcessorGBLTrackFit::init( )---------------------------------------------END" << std::endl;
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
		streamlog_out(MESSAGE9)<< "Unknown exception in init function of EUTelProcessorGBLTrackFit." << std::endl;
		throw marlin::StopProcessingException( this ) ;
	}
}

void EUTelProcessorGBLTrackFit::processRunHeader(LCRunHeader * run) {
	std::auto_ptr<EUTelRunHeaderImpl> header(new EUTelRunHeaderImpl(run));
	header->addProcessor(type());
	// this is the right place also to check the geometry ID. This is a
	// unique number identifying each different geometry used at the
	// beam test. The same number should be saved in the run header and
 	// in the xml file. If the numbers are different, warn the user.
	if (header->getGeoID() == 0)
 		streamlog_out(WARNING0) << "The geometry ID in the run header is set to zero." << std::endl << "This may mean that the GeoID parameter was not set" << std::endl;
	if ((unsigned int)header->getGeoID() != geo::gGeometry().getSiPlanesLayoutID()) {  
		streamlog_out(WARNING5) << "Error during the geometry consistency check: " << std::endl << "The run header says the GeoID is " << header->getGeoID() << std::endl << "The GEAR description says is     " << geo::gGeometry().getSiPlanesLayoutID() << std::endl;
	}
	_chi2NdfVec.clear();//TO DO:This is needed to determine if the track is near chi2 of one. Do we need this however?
    
	_nProcessedRuns++;
}

void EUTelProcessorGBLTrackFit::processEvent(LCEvent* evt){
	try{
		streamlog_out(DEBUG5) << "Start of event " << _nProcessedEvents << std::endl;

		EUTelEventImpl * event = static_cast<EUTelEventImpl*> (evt); ///We change the class so we can use EUTelescope functions

		if (event->getEventType() == kEORE) {
			streamlog_out(DEBUG4) << "EORE found: nothing else to do." << std::endl;
			return;
		}else if (event->getEventType() == kUNKNOWN) {
			streamlog_out(WARNING2) << "Event number " << event->getEventNumber() << " in run " << event->getRunNumber() << " is of unknown type. Continue considering it as a normal Data Event." << std::endl;
		}
        EUTelReaderGenericLCIO reader = EUTelReaderGenericLCIO();
        std::vector<EUTelTrack> tracks = reader.getTracks(evt, _trackCandidatesInputCollectionName );
		std::vector<EUTelTrack> allTracksForThisEvent;//GBL will analysis the track one at a time. However we want to save to lcio per event.
		for (size_t iTrack = 0; iTrack < tracks.size(); iTrack++) {
			EUTelTrack track = tracks.at(iTrack); 
            streamlog_out(DEBUG1)<<"Found "<<tracks.size()<<" tracks for event " << evt->getEventNumber() << "  This is track:  " << iTrack <<std::endl;
            track.print();
			streamlog_out(DEBUG1) << "//////////////////////////////////// " << std::endl;
			_trackFitter->resetPerTrack(); //Here we reset the label that connects state to GBL point to 1 again. Also we set the list of states->labels to 0
			_trackFitter->testTrack(track);//Check the track has states and hits  
			std::vector< gbl::GblPoint > pointList;
			_trackFitter->setInformationForGBLPointList(track, pointList);//Here we describe the whole setup. Geometry, scattering, data...
			const gear::BField& B = geo::gGeometry().getMagneticField();//We need this to determine if we should fit a curve or a straight line.
			const double Bmag = B.at( TVector3(0.,0.,0.) ).r2();
			_trackFitter->setPairMeasurementStateAndPointLabelVec(pointList);//This will create a link between the states that have a hit associated with them and the GBL label that is associated with the state.
			gbl::GblTrajectory* traj = 0;
			//Here we create the trajectory from the points created by setInformationForGBLPointList. This will take the points and propagation jacobian and split this into smaller matrices to describe the problem in terms of offsets. Here is the difference between GBL and other fitting algorithms.  
			if ( Bmag < 1.E-6 ) {
				traj = new gbl::GblTrajectory( pointList, false ); 
			}else {
				traj = new gbl::GblTrajectory( pointList, true );
			}
			_trackFitter->setPairAnyStateAndPointLabelVec(traj);//This will create a link between any state and it's GBL point label. 
			double  chi2=0; 
			int ndf=0;
			int ierr=0;
			_trackFitter->computeTrajectoryAndFit(traj, &chi2,&ndf, ierr);//This will do the minimisation of the chi2 and produce the most probable trajectory.
			if(ierr == 0 ){
				streamlog_out(DEBUG5) << "Ierr is: " << ierr << " Entering loop to update track information " << std::endl;
				//If the fit succeeded then write into the binary file.
//				traj->milleOut(*(_Mille->_milleGBL));
				static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_chi2CandidateHistName ] ) -> fill( (chi2)/(ndf));
				static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_fitsuccessHistName ] ) -> fill(1.0);
				if(chi2 ==0 or ndf ==0){
					throw(lcio::Exception("Your fitted track has zero degrees of freedom or a chi2 of 0.")); 	
					}
				track.setChi2(chi2);
				track.setNdf(ndf);
				_chi2NdfVec.push_back(chi2/static_cast<float>(ndf));
				std::map<int, std::vector<double> >  mapSensorIDToCorrectionVec;//This is not used now. However it maybe useful to be able to access the corrections that GBL makes to the original track. Since if this is too large then GBL may give th wrong trajectory. Since all the equations are only to first order. 
				_trackFitter->updateTrackFromGBLTrajectory(traj,track,mapSensorIDToCorrectionVec);
				std::map< int, std::map< float, float > >  SensorResidual; 
				std::map< int, std::map< float, float > >  SensorResidualError; 
				_trackFitter->getResidualOfTrackandHits(traj, pointList,track, SensorResidual, SensorResidualError);
				if(chi2/static_cast<float>(ndf) < 5){
				  plotResidual(SensorResidual,SensorResidualError);
				}
			}else{
				streamlog_out(DEBUG5) << "Ierr is: " << ierr << " Do not update track information " << std::endl;
				static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_fitsuccessHistName ] ) -> fill(0.0);
				continue;//We continue so we don't add an empty track
			}	
			allTracksForThisEvent.push_back(track);
			}//END OF LOOP FOR ALL TRACKS IN AN EVENT
			outputLCIO(evt, allTracksForThisEvent); 
			allTracksForThisEvent.clear();//We clear this so we don't add the same track twice
			streamlog_out(DEBUG5) << "End of event " << _nProcessedEvents << std::endl;
			_nProcessedEvents++;
	}
	catch (DataNotAvailableException e) {
		streamlog_out(MESSAGE0) << _trackCandidatesInputCollectionName << " collection not available" << std::endl;
		throw marlin::SkipEventException(this);
	}
	catch(std::string &e){
		streamlog_out(MESSAGE9) << e << std::endl;
		throw marlin::SkipEventException( this ) ;
	}
	catch(lcio::Exception& e){
		streamlog_out(MESSAGE9) << e.what() <<std::endl;
		throw marlin::StopProcessingException( this ) ;
	}
	catch(...){
		streamlog_out(MESSAGE9)<<"Unknown exception in processEvent function of EUTelProcessorGBLTrackFit" <<std::endl;
		throw marlin::StopProcessingException( this ) ;
	}
	
}


//TO DO:This is a very stupid way to histogram but will add new class to do this is long run 
void EUTelProcessorGBLTrackFit::plotResidual(std::map< int, std::map<float, float > >  & sensorResidual, std::map< int, std::map<float, float > >  & sensorResidualError){
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Residual plot
	/* Obtain DUT IDs */
        std::map< int, int>::iterator planes_it;
	int dut1 = -10;
	int dut2 = -10;
	bool flag = false;
	for ( std::vector<int>::const_iterator itPla = geo::gGeometry().sensorIDsVec().begin(); itPla != geo::gGeometry().sensorIDsVec().end(); ++itPla) {
	  if ( *itPla> 5 && flag == true) {
	    dut2 = *itPla;
	  }
	  if ( *itPla> 5) {
	    dut1 = *itPla;
	    flag = true;
	  }
	  
	}
	/* Fill histograms */
	std::map< int, std::map< float, float > >::iterator sensor_residual_it;
	for(sensor_residual_it = sensorResidual.begin(); sensor_residual_it != sensorResidual.end(); sensor_residual_it++) {
	  std::map<float, float> map = sensor_residual_it->second;
	  if( !map.empty()){
	    float res = map.begin()->first;
	    if( sensor_residual_it->first == 0 ){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX0 ] ) -> fill(res);}
	    if( sensor_residual_it->first == 1 ){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX1 ] ) -> fill(res);}
	    if( sensor_residual_it->first == 2 ){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX2 ] ) -> fill(res);}
	    if( sensor_residual_it->first == 3 ){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX3 ] ) -> fill(res);}
	    if( sensor_residual_it->first == 4 ){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX4 ] ) -> fill(res);}
	    if( sensor_residual_it->first == 5 ){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX5 ] ) -> fill(res);}
	    if( sensor_residual_it->first == dut1 ){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameXDut1 ] ) -> fill(res);}
	    if( sensor_residual_it->first == dut2 ){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameXDut2 ] ) -> fill(res);}
	    
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
			if( sensor_residual_it->first == dut1){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameYDut1 ] ) -> fill(res2);}
			if( sensor_residual_it->first == dut2){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameYDut2 ] ) -> fill(res2);}

				
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Residual Error plot
	std::map< int, std::map< float, float > >::iterator sensor_residualerror_it;
	for(sensor_residualerror_it = sensorResidualError.begin(); sensor_residualerror_it != sensorResidualError.end(); sensor_residualerror_it++) {
	  std::map<float, float> maperror = sensor_residualerror_it->second;
	  std::map<float, float> mapres = sensorResidual.at(sensor_residualerror_it->first);

	  if( !maperror.empty()){
	    float res = mapres.begin()->first;	
	    float reserror = maperror.begin()->first;
	    if( sensor_residualerror_it->first == 0){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX0p ] ) -> fill(res/reserror);}
	    if( sensor_residualerror_it->first == 1){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX1p ] ) -> fill(res/reserror);}
	    if( sensor_residualerror_it->first == 2){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX2p ] ) -> fill(res/reserror);}
	    if( sensor_residualerror_it->first == 3){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX3p ] ) -> fill(res/reserror);}
	    if( sensor_residualerror_it->first == 4){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX4p ] ) -> fill(res/reserror);}
	    if( sensor_residualerror_it->first == 5){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameX5p ] ) -> fill(res/reserror);}
	    if( sensor_residualerror_it->first == dut1){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameXDut1p ] ) -> fill(res/reserror);}
	    if( sensor_residualerror_it->first == dut2){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameXDut2p ] ) -> fill(res/reserror);}
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
	  if( sensor_residualerror_it->first == dut1){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameYDut1p ] ) -> fill(res2/res2error);}
	  if( sensor_residualerror_it->first == dut2){static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_residGblFitHistNameYDut2p ] ) -> fill(res2/res2error);}
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}


void EUTelProcessorGBLTrackFit::end() {
	float total = 0;
	double sizeFittedTracks = _chi2NdfVec.size();
	for(size_t i=0; i<_chi2NdfVec.size(); ++i)
	{
		total= total + _chi2NdfVec.at(i);//TO DO: This is does not seem to output the correct average chi2. Plus do we really need this to fit?
	}
	//TO DO: We really should have a better way to look track per track	and see if the correction is too large. 
	std::vector<double> correctionTotal = _trackFitter->getCorrectionsTotal();
	streamlog_out(MESSAGE9)<<"This is the average correction for omega: " <<correctionTotal.at(0)/sizeFittedTracks<<std::endl;	
	streamlog_out(MESSAGE9)<<"This is the average correction for local xz inclination: " <<correctionTotal.at(1)/sizeFittedTracks<<std::endl;	
	streamlog_out(MESSAGE9)<<"This is the average correction for local yz inclination: " <<correctionTotal.at(2)/sizeFittedTracks<<std::endl;	
	streamlog_out(MESSAGE9)<<"This is the average correction for local x: " <<correctionTotal.at(3)/sizeFittedTracks<<std::endl;	
	streamlog_out(MESSAGE9)<<"This is the average correction for local y: " <<correctionTotal.at(4)/sizeFittedTracks<<std::endl;	

  float average = total/sizeFittedTracks;
	streamlog_out(MESSAGE9) << "This is the average chi2 -"<< average <<std::endl;

}

#endif // USE_GBL

void EUTelProcessorGBLTrackFit::outputLCIO(LCEvent* evt, std::vector<EUTelTrack>& tracks){
    if(!tracks.empty()){
        for(unsigned int i=0 ; i< tracks.size(); i++){
            streamlog_out(DEBUG1)<<"Found "<<tracks.size()<<" track for event " << evt->getEventNumber() <<".  Track number  " << i <<std::endl;
            tracks.at(i).print();
        }
        EUTelReaderGenericLCIO reader = EUTelReaderGenericLCIO();
        reader.getColVec(tracks, evt,_tracksOutputCollectionName);
    }
}

void EUTelProcessorGBLTrackFit::bookHistograms() {
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
 try {
        streamlog_out(DEBUG) << "Booking histograms..." << std::endl;

        std::auto_ptr<EUTelHistogramManager> histoMgr( new EUTelHistogramManager( _histoInfoFileName ));
        EUTelHistogramInfo    * histoInfo;
        bool                    isHistoManagerAvailable;

////////////////////////////////////////////////////////This is for the residual//Thi si a hack must fix so can accept any number of planes. Really should be a separate processor

        int NBinX=2000;
        double MinX=-0.5;  //-0.2;
        double MaxX=0.5;


        AIDA::IHistogram1D * residGblFit0X = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameX0, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit1X = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameX1, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit2X = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameX2, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit3X = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameX3, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit4X = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameX4, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit5X = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameX5, NBinX, MinX, MaxX); 
	AIDA::IHistogram1D * residGblFitDut1X = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameXDut1, NBinX, -0.4, 0.4);
	AIDA::IHistogram1D * residGblFitDut2X = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameXDut2, NBinX, -0.4, 0.4);
        

              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX0, residGblFit0X));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX1, residGblFit1X));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX2, residGblFit2X));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX3, residGblFit3X));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX4, residGblFit4X));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX5, residGblFit5X));
	      _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameXDut1, residGblFitDut1X));
	      _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameXDut2, residGblFitDut2X));

											
        AIDA::IHistogram1D * residGblFit0Y = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameY0, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit1Y = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameY1, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit2Y = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameY2, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit3Y = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameY3, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit4Y = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameY4, NBinX, MinX, MaxX); 
        AIDA::IHistogram1D * residGblFit5Y = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameY5, NBinX, MinX, MaxX); 
	AIDA::IHistogram1D * residGblFitDut1Y = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameYDut1, NBinX, -0.4, 0.4); 
	AIDA::IHistogram1D * residGblFitDut2Y = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameYDut2, NBinX, -0.4, 0.4); 

              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameY0, residGblFit0Y));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameY1, residGblFit1Y));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameY2, residGblFit2Y));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameY3, residGblFit3Y));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameY4, residGblFit4Y));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameY5, residGblFit5Y));
	      _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameYDut1, residGblFitDut1Y));
	      _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameYDut2, residGblFitDut2Y));

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Pull plots											
        int NBinp=120;
        double MinXp=-10;
        double MaxXp=10;


        AIDA::IHistogram1D * residGblFit0Xp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameX0p, NBinp, MinXp, MaxXp); 
        AIDA::IHistogram1D * residGblFit1Xp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameX1p, NBinp, MinXp, MaxXp); 
        AIDA::IHistogram1D * residGblFit2Xp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameX2p, NBinp, MinXp, MaxXp); 
        AIDA::IHistogram1D * residGblFit3Xp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameX3p, NBinp, MinXp, MaxXp); 
        AIDA::IHistogram1D * residGblFit4Xp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameX4p, NBinp, MinXp, MaxXp); 
        AIDA::IHistogram1D * residGblFit5Xp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameX5p, NBinp, MinXp, MaxXp);
	AIDA::IHistogram1D * residGblFitDut1Xp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameXDut1p, NBinp, MinXp, MaxXp);
	AIDA::IHistogram1D * residGblFitDut2Xp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameXDut2p, NBinp, MinXp, MaxXp);


              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX0p, residGblFit0Xp));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX1p, residGblFit1Xp));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX2p, residGblFit2Xp));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX3p, residGblFit3Xp));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX4p, residGblFit4Xp));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameX5p, residGblFit5Xp));
	      _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameXDut1p, residGblFitDut1Xp));
	      _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameXDut2p, residGblFitDut2Xp));
								      
        AIDA::IHistogram1D * residGblFit0Yp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameY0p, NBinp, MinXp, MaxXp); 
        AIDA::IHistogram1D * residGblFit1Yp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameY1p, NBinp, MinXp, MaxXp); 
        AIDA::IHistogram1D * residGblFit2Yp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameY2p, NBinp, MinXp, MaxXp); 
        AIDA::IHistogram1D * residGblFit3Yp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameY3p, NBinp, MinXp, MaxXp); 
        AIDA::IHistogram1D * residGblFit4Yp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameY4p, NBinp, MinXp, MaxXp); 
        AIDA::IHistogram1D * residGblFit5Yp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameY5p, NBinp, MinXp, MaxXp); 
	AIDA::IHistogram1D * residGblFitDut1Yp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameYDut1p, NBinp, MinXp, MaxXp);
	AIDA::IHistogram1D * residGblFitDut2Yp = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_residGblFitHistNameYDut2p, NBinp, MinXp, MaxXp);
		

              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameY0p, residGblFit0Yp));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameY1p, residGblFit1Yp));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameY2p, residGblFit2Yp));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameY3p, residGblFit3Yp));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameY4p, residGblFit4Yp));
              _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameY5p, residGblFit5Yp));
	      _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameYDut1p, residGblFitDut1Yp));
	      _aidaHistoMap1D.insert(std::make_pair(_histName::_residGblFitHistNameYDut2p, residGblFitDut2Yp));



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
        streamlog_out(WARNING2) << "Can't allocate histograms. Continue without histogramming" << std::endl;
}
#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

}

