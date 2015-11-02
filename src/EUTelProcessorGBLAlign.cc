
#include "EUTelProcessorGBLAlign.h"

using namespace eutelescope;

EUTelProcessorGBLAlign::EUTelProcessorGBLAlign() :
Processor("EUTelProcessorGBLAlign"),
_milleBinaryFilename("mille.bin"),
_milleSteeringFilename("pede-steer.txt"),
_milleResultFileName("millepede.res"),
_gear_aligned_file("gear-00001-aligned.xml"),
_nProcessedRuns(0),
_nProcessedEvents(0),
_beamQ(-1),
_eBeam(4),
_createBinary(true),
_mEstimatorType()
{
  // TrackerHit input collection
  registerInputCollection(LCIO::TRACK, "TrackCandidatesInputCollectionName", "Input track candidate collection name",_trackCandidatesInputCollectionName,std::string("TrackCandidatesCollection"));

  // Track output collection
  registerOutputCollection(LCIO::TRACK,"TracksOutputCollectionName","Output tracks collection name",_tracksOutputCollectionName,std::string("TrackCollection"));

	registerOptionalParameter("MilleBinaryFilename", "Name of the Millepede binary file", _milleBinaryFilename, std::string("mille.bin"));

	registerOptionalParameter("MilleSteeringFilename", "Name of the Millepede steering file to be created", _milleSteeringFilename, std::string("pede-steer.txt"));
    
	registerOptionalParameter("MilleResultFilename", "Name of the Millepede result file", _milleResultFileName, std::string("millepede.res"));
    
	registerOptionalParameter("GearAlignedFile", "Suffix to add to the new Gear with alignment corrections", _gear_aligned_file, std::string("gear-00001-aligned.xml"));

  registerOptionalParameter("BeamCharge", "Beam charge [e]", _beamQ, static_cast<double> (-1));

  // Necessary processor parameters that define fitter settings
  registerProcessorParameter("BeamEnergy", "Beam energy [GeV]", _eBeam, static_cast<double> (4.0));
  registerOptionalParameter("ExcludePlanes", "This is the planes that will not be included in analysis", _excludePlanes ,IntVec());


  // Optional processor parameters that define finder settings

  registerOptionalParameter("GBLMEstimatorType", "GBL outlier down-weighting option (t,h,c)", _mEstimatorType, std::string() );
  registerOptionalParameter("IncMed", "Do you want to include the medium as addtional scattering", _incMed ,int(0));

  registerOptionalParameter("CreateBinary", "Should we create a binary file for millepede containing the data that millepede needs  ", _createBinary, bool(true));

  registerOptionalParameter("xResolutionPlane", "x resolution of planes given in Planes", _SteeringxResolutions, FloatVec());
  registerOptionalParameter("yResolutionPlane", "y resolution of planes given in Planes", _SteeringyResolutions, FloatVec());

    registerOptionalParameter("FixedAlignmentPlanesXshift", "Ids of planes for which X shift will be fixed during millepede call", _fixedAlignmentXShfitPlaneIds, IntVec());
    
    registerOptionalParameter("FixedAlignmentPlanesYshift", "Ids of planes for which Y shift will be fixed during millepede call", _fixedAlignmentYShfitPlaneIds, IntVec());
    
    registerOptionalParameter("FixedAlignmentPlanesZshift", "Ids of planes for which Z shift will be fixed during millepede call", _fixedAlignmentZShfitPlaneIds, IntVec());
    
    registerOptionalParameter("FixedAlignmentPlanesXrotation", "Ids of planes for which rotation around X will be fixed during millepede call", _fixedAlignmentXRotationPlaneIds, IntVec());
    
    registerOptionalParameter("FixedAlignmentPlanesYrotation", "Ids of planes for which rotation around Y will be fixed during millepede call", _fixedAlignmentYRotationPlaneIds, IntVec());
    
    registerOptionalParameter("FixedAlignmentPlanesZrotation", "Ids of planes for which rotation around Z will be fixed during millepede call", _fixedAlignmentZRotationPlaneIds, IntVec());

    registerOptionalParameter("PedeSteeringAdditionalCmds","FOR EXPERTS: List of commands that should be included in the pede steering file. Use '\\' to seperate options and introduce a line break.",_pedeSteerAddCmds, StringVec());

		registerOptionalParameter("AlignmentConstantLCIOFile","This is the name of the LCIO file name with the output alignment"
                            "constants (add .slcio)",_alignmentConstantLCIOFile, static_cast< std::string > ( "alignment.slcio" ) );
		registerOptionalParameter("ExcludePlanes", "This is the planes that will not be included in analysis", _excludePlanes ,IntVec());


}

void EUTelProcessorGBLAlign::init() {
	try{
		streamlog_out(DEBUG2) << "EUTelProcessorGBLAlign::init( )---------------------------------------------BEGIN" << std::endl;
		_nProcessedRuns = 0;
		_nProcessedEvents = 0;
		std::string name("test.root");
		geo::gGeometry().initializeTGeoDescription(name,false);
        EUTelExcludedPlanes::setRelativeComplementSet(_excludePlanes);///Only used if GBL will do the track parameterisation.

		// Initialize GBL fitter
		EUTelGBLFitter* Fitter = new EUTelGBLFitter();
		_Mille  = new EUTelMillepede(_milleResultFileName, _gear_aligned_file); 
		_Mille->setSteeringFileName(_milleSteeringFilename);// The steering file will store the labels for global variables in text file, along with errors and seeds guess.
		_Mille->setXShiftFixed(_fixedAlignmentXShfitPlaneIds);
		_Mille->setYShiftFixed(_fixedAlignmentYShfitPlaneIds);
		_Mille->setZShiftFixed(_fixedAlignmentZShfitPlaneIds);
		_Mille->setXRotationsFixed(_fixedAlignmentXRotationPlaneIds);
		_Mille->setYRotationsFixed(_fixedAlignmentYRotationPlaneIds);
		_Mille->setZRotationsFixed(_fixedAlignmentZRotationPlaneIds);
		_Mille->setBinaryFileName(_milleBinaryFilename);//The binary file holds for each state: Hold all the information needed for Millepede to work 
		_Mille->setResultsFileName(_milleResultFileName);
		_Mille->testUserInput();
		_Mille->printFixedPlanes();
		Fitter->setMEstimatorType(_mEstimatorType);//Outliers are hits that do not appear to follow errors which are Gaussian. We want to downweight the effect these hits have on the fit.
		Fitter->setParamterIdXResolutionVec(_SteeringxResolutions);//We set the accuracy of the residual information since we have no correct hit error analysis yet.
		Fitter->setParamterIdYResolutionVec(_SteeringyResolutions);
		Fitter->setMillepede(_Mille);//We need to have a connection between GBL and Millepede since GBL knows nothing about sensor orientations.
        Fitter->setIncMed(_incMed); 
		Fitter->testUserInput();
		_trackFitter = Fitter;


		if (!_trackFitter) {
			streamlog_out(ERROR) << "Can't allocate an instance of EUTelGBLFitter. Stopping ..." << std::endl;
			throw UnknownDataTypeException("Track finder was not created");
		}

		streamlog_out(DEBUG2) << "EUTelProcessorGBLAlign::init( )---------------------------------------------END" << std::endl;
	}	
	catch (DataNotAvailableException e) {
		streamlog_out(MESSAGE0) << _trackCandidatesInputCollectionName << " collection not available" << std::endl;
		throw marlin::SkipEventException(this);
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
		streamlog_out(MESSAGE9)<<"Unknown exception in processEvent function of EUTelProcessorGBLTrackFit" <<std::endl;
		throw marlin::StopProcessingException( this ) ;
	}
}

void EUTelProcessorGBLAlign::processRunHeader(LCRunHeader * run) {
	std::auto_ptr<EUTelRunHeaderImpl> header(new EUTelRunHeaderImpl(run));
	header->addProcessor(type());


	// this is the right place also to check the geometry ID. This is a
	// unique number identifying each different geometry used at the
	// beam test. The same number should be saved in the run header and
 	// in the xml file. If the numbers are different, warn the user.

	if (header->getGeoID() == 0)	streamlog_out(WARNING0) << "The geometry ID in the run header is set to zero." << std::endl << "This may mean that the GeoID parameter was not set" << std::endl;
  	if ((unsigned int)header->getGeoID() != geo::gGeometry().getSiPlanesLayoutID()) {  
			streamlog_out(WARNING5) << "Error during the geometry consistency check: " << std::endl << "The run header says the GeoID is " << header->getGeoID() << std::endl << "The GEAR description says is     " << geo::gGeometry().getSiPlanesLayoutID() << std::endl;
  	}
    
	_nProcessedRuns++;
	_totalTrackCount=0;	
}

void EUTelProcessorGBLAlign::processEvent(LCEvent * evt){
	try{
		if(_createBinary){
			EUTelEventImpl * event = static_cast<EUTelEventImpl*> (evt); ///We change the class so we can use EUTelescope functions

			if (event->getEventType() == kEORE) {
				streamlog_out(DEBUG4) << "EORE found: nothing else to do." << std::endl;
				return;
			}else if (event->getEventType() == kUNKNOWN) {
				streamlog_out(WARNING2) << "Event number " << event->getEventNumber() << " in run " << event->getRunNumber() << " is of unknown type. Continue considering it as a normal Data Event." << std::endl;
			}
            EUTelReaderGenericLCIO reader = EUTelReaderGenericLCIO();
            std::vector<EUTelTrack> tracks = reader.getTracks(evt, _trackCandidatesInputCollectionName);
            for (size_t iTrack = 0; iTrack < tracks.size(); ++iTrack) {
                _totalTrackCount++;
                EUTelTrack track = tracks.at(iTrack);
                std::vector< gbl::GblPoint > pointList;
                std::map<  unsigned int,unsigned int >  linkGL;
                std::map< unsigned int, unsigned int >  linkMeas;
                ///This will create the initial GBL trajectory
                _trackFitter->getGBLPointsFromTrack(track, pointList);
                ///NOTE: This is the only difference between a GBL track fit and alignment step with respect to GBL
                /// The rest of the work come from reading this to the gear file and making sure the transformations are correct.
                _trackFitter->getGloPar(pointList,track); 
                const gear::BField& B = geo::gGeometry().getMagneticField();
                const double Bmag = B.at( TVector3(0.,0.,0.) ).r2();
                gbl::GblTrajectory* traj = 0;
//					printPointsInformation(pointList);
                if ( Bmag < 1.E-6 ) {
                    traj = new gbl::GblTrajectory( pointList, false ); //Must make sure this is not a memory leak
                } else {
                    traj = new gbl::GblTrajectory( pointList, true );
                }
                double chi2, loss;
                int ndf2;
                traj->fit(chi2, ndf2, loss, _mEstimatorType );
                streamlog_out ( DEBUG0 ) << "This is the trajectory we are just about to fit: " << std::endl;
                streamlog_message( DEBUG0, traj->printTrajectory(10);, std::endl; );
//				std::cout<<"WRITE TO MILLEPEDE. EVENT: " << 	event->getEventNumber() << "  Total number of tracks: " << _totalTrackCount << std::endl;	
                traj->milleOut(*(_Mille->_milleGBL));
            }//END OF LOOP FOR ALL TRACKS IN AN EVENT
//			if(event->getEventNumber() == 1){
//				throw marlin::StopProcessingException( this ) ;
//			}
		}
	}
	catch (DataNotAvailableException e) {
//		streamlog_out(MESSAGE9) << "Data not avaliable skip event. " << std::endl;
		throw marlin::SkipEventException(this);
	}
	catch(std::string &e){
//		streamlog_out(MESSAGE9) << e << std::endl;
		throw marlin::SkipEventException( this ) ;
	}
	catch(lcio::Exception& e){
		streamlog_out(MESSAGE9) << e.what() <<std::endl;
		throw marlin::StopProcessingException( this ) ;
	}
	catch(...){
		streamlog_out(MESSAGE9)<<"Unknown exception in processEvent function of EUTelProcessorGBLAlign" <<std::endl;
		throw marlin::StopProcessingException( this ) ;
	}

}

void EUTelProcessorGBLAlign::end(){
	_Mille->_milleGBL->~MilleBinary();
	streamlog_out (MESSAGE9) <<"TOTAL NUMBER OF TRACKS PASSED TO ALIGNMENT: "<< _totalTrackCount << std::endl;
	if(_totalTrackCount<1000){
		streamlog_out(WARNING5)<<"You are trying to align with fewer than 1000 tracks. This could be too small a number." <<std::endl;
	}
    /// Create inital steering file
	_Mille->writeMilleSteeringFile(_pedeSteerAddCmds);
    /// This will run pede and link to the results file: _milleResultFileName
	bool tooManyRejects = 	_Mille->runPede();
	if(!tooManyRejects){
        _Mille->getNewGear();
	}else{
		streamlog_out (MESSAGE9) <<"THE NUMBER OF REJECTED TRACKS IS TOO LARGE."<< std::endl;
	}
}

void EUTelProcessorGBLAlign::printPointsInformation(std::vector<gbl::GblPoint>& pointList){
	typedef std::vector<gbl::GblPoint>::iterator IteratorType;
	streamlog_out(MESSAGE5) << "THE START OF THE TRACK POINTS///////////////" <<std::endl;
	for(IteratorType point = pointList.begin(); point != pointList.end(); point++){
		streamlog_out(MESSAGE5)<<std::endl <<"Point label: " << point->getLabel()<<std::endl;
		streamlog_out(MESSAGE5) << "GLOBAL DERIVATIVE MATRIX"<<std::endl;
		streamlog_message( MESSAGE5, point->getGlobalDerivatives().Print();, std::endl; );
		std::vector<int> label = point->getGlobalLabels();
//		streamlog_out(MESSAGE5)<<"Size  of : "<<label.size() <<std::endl;
		streamlog_out(MESSAGE5) << "GLOBAL LABELS: " <<std::endl;
		for( std::vector<int>::const_iterator i = label.begin(); i != label.end(); ++i){
			streamlog_out(MESSAGE5) << *i << ' ';
		}

	}
}	
double  EUTelProcessorGBLAlign::printSize(const std::string& address) {
	std::fstream motd(address.c_str(), std::ios::binary|std::ios::in|std::ios::ate);
	if(motd) {
		std::fstream::pos_type size = motd.tellg();
		return size;
	} else {
		perror(address.c_str());
		return 0.0;
	}
}
