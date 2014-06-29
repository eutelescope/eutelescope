
#include "EUTelProcessorMilleAlign.h"

using namespace eutelescope;

EUTelProcessorMilleAlign::EUTelProcessorMilleAlign() :
Processor("EUTelProcessorMilleAlign"),
_milleBinaryFilename("mille.bin"),
_milleSteeringFilename("pede-steer.txt"),
_milleResultFileName("millepede.res"),
_gear_aligned_file("gear-00001-aligned.xml"),
_nProcessedRuns(0),
_nProcessedEvents(0),
_beamQ(-1),
_eBeam(4),
_mEstimatorType(),
_maxChi2Cut(1000),
_alignmentMode(0){

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

  // Optional processor parameters that define finder settings

  registerOptionalParameter("GBLMEstimatorType", "GBL outlier down-weighting option (t,h,c)", _mEstimatorType, string() );

  registerOptionalParameter("MilleMaxChi2Cut", "Maximum chi2 of a track candidate that goes into millepede", _maxChi2Cut, double(1000.));

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

    registerOptionalParameter("AlignmentPlanes", "Ids of planes to be used in alignment", _alignmentPlaneIds, IntVec());

    registerOptionalParameter("FixedAlignmentPlanesXshift", "Ids of planes for which X shift will be fixed during millepede call", _fixedAlignmentXShfitPlaneIds, IntVec());
    
    registerOptionalParameter("FixedAlignmentPlanesYshift", "Ids of planes for which Y shift will be fixed during millepede call", _fixedAlignmentYShfitPlaneIds, IntVec());
    
    registerOptionalParameter("FixedAlignmentPlanesZshift", "Ids of planes for which Z shift will be fixed during millepede call", _fixedAlignmentZShfitPlaneIds, IntVec());
    
    registerOptionalParameter("FixedAlignmentPlanesXrotation", "Ids of planes for which rotation around X will be fixed during millepede call", _fixedAlignmentXRotationPlaneIds, IntVec());
    
    registerOptionalParameter("FixedAlignmentPlanesYrotation", "Ids of planes for which rotation around Y will be fixed during millepede call", _fixedAlignmentYRotationPlaneIds, IntVec());
    
    registerOptionalParameter("FixedAlignmentPlanesZrotation", "Ids of planes for which rotation around Z will be fixed during millepede call", _fixedAlignmentZRotationPlaneIds, IntVec());

    registerOptionalParameter("PedeSteeringAdditionalCmds","FOR EXPERTS: List of commands that should be included in the pede steering file. Use '\\' to seperate options and introduce a line break.",_pedeSteerAddCmds, StringVec());

   registerOptionalParameter("AlignmentConstantLCIOFile","This is the name of the LCIO file name with the output alignment"
                            "constants (add .slcio)",_alignmentConstantLCIOFile, static_cast< string > ( "alignment.slcio" ) );


}

void EUTelProcessorMilleAlign::init() {

	streamlog_out(DEBUG2) << "EUTelProcessorMilleAlign::init( )---------------------------------------------BEGIN" << std::endl;

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
	_Mille  = new EUTelMillepede(_alignmentMode);
  Fitter->SetBeamCharge(_beamQ);
  Fitter->SetBeamEnergy(_eBeam);
	Fitter->setMEstimatorType(_mEstimatorType);
  Fitter->SetMilleBinaryName(_milleBinaryFilename);
  Fitter->SetChi2Cut(_maxChi2Cut);
	Fitter->SetMillepede(_Mille);
  _trackFitter = Fitter;


	if (!_trackFitter) {
		streamlog_out(ERROR) << "Can't allocate an instance of EUTelGBLFitter. Stopping ..." << std::endl;
		throw UnknownDataTypeException("Track finder was not created");
	}

	streamlog_out(DEBUG2) << "EUTelProcessorMilleAlign::init( )---------------------------------------------END" << std::endl;




}

void EUTelProcessorMilleAlign::processRunHeader(LCRunHeader * run) {
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

void EUTelProcessorMilleAlign::processEvent(LCEvent * evt){
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
		std::vector< EUTelTrackImpl* > EUtracks;	
		for (int iCol = 0; iCol < col->getNumberOfElements(); iCol++) {
			
			if (!col) {
		      		streamlog_out(WARNING2) << "Track collection not found found for event " << _nProcessedEvents << " in run " << _nProcessedRuns << endl;
		        	throw SkipEventException(this);
     			}
	   		IMPL::TrackImpl* trackimpl = static_cast<IMPL::TrackImpl*> (col->getElementAt(iCol));
			EUTelTrackImpl* EUtrack = new EUTelTrackImpl(*trackimpl);
      			streamlog_out(DEBUG1) << "Track " << iCol << " nhits " << trackimpl->getTrackerHits().size() << endl;

			///Create points but do not fit this time
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
	streamlog_out ( DEBUG0 ) << "This is the trajectory we are just about to fit: " << endl;
	  streamlog_message( DEBUG0, traj->printTrajectory(10);, std::endl; );
			
			//////////////////////////////////////////////////////////////////////////////////////END
       _trackFitter->CreateAlignmentToMeasurementJacobian(&pointList ); //This is place in GBLFitter since millepede has not idea about states and points. Only GBLFitter know about that
			
			              _trackFitter->Clear();
		}//END OF LOOP FOR ALL TRACKS IN AN EVENT
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//	outputLCIO(evt, EUtracks); // Important to note that EUtracks are still only the original states at input. The scatterers are used in the fit but are not included here.
	}//END OF COLLECTION IS NOT NULL LOOP	


}
void EUTelProcessorMilleAlign::check(LCEvent * evt){}

void EUTelProcessorMilleAlign::end(){
_Mille->writeMilleSteeringFile(_pedeSteerAddCmds);
_Mille->runPede();
_Mille->parseMilleOutput(_alignmentConstantLCIOFile, _gear_aligned_file);



}
