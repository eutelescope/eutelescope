#ifdef USE_GBL

// C++
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <cstdio>
#include <algorithm>

// LCIO
#include <EVENT/LCCollection.h>

// MARLIN
#include "marlin/Exceptions.h"
#include "marlin/Global.h"
#include "marlin/Processor.h"
#include "marlin/VerbosityLevels.h"

// ROOT
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TMath.h"
#endif

// EUTELESCOPE
#include "EUTelProcessorAlignmentGBLMille.h"

#include "EUTelHistogramManager.h"
#include "EUTelPStream.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelUtility.h"
#include "EUTelLCObjectTrackCandidate.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

// Cluster types
class EUTelSparseClusterImpl;
class EUTelBrickedClusterImpl;
class EUTelDFFClusterImpl;
class EUTelFFClusterImpl;

using namespace lcio;
using namespace std;
using namespace marlin;
using namespace eutelescope;

/**  EUTelProcessorAlignmentGBLMille
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Requires collection of track candidates hits.
 *
 *  <h4>Output</h4> 
 *  <li> Histograms.
 *  <li> Collection of fitted tracks.
 */


/** Default constructor */
EUTelProcessorAlignmentGBLMille::EUTelProcessorAlignmentGBLMille() :
Processor("EUTelProcessorAlignmentGBLMille"),
_qBeam(-1),
_eBeam(-1.),
_alignmentMode(0),
_mEstimatorType(),
_xShiftsMap(),
_yShiftsMap(),
_zShiftsMap(),
_xRotationsMap(),
_yRotationsMap(),
_zRotationsMap(),
_milleBinaryFilename("mille.bin"),
_milleSteeringFilename("pede-steer.txt"),
_milleResultFileName("millepede.res"),
_pedeSteerAddCmds(),
_alignmentPlaneIds(),
_planeIds(),
_SteeringxResolutions(),
_SteeringyResolutions(),
_fixedAlignmentXShfitPlaneIds(),
_fixedAlignmentYShfitPlaneIds(),
_fixedAlignmentZShfitPlaneIds(),
_fixedAlignmentXRotationPlaneIds(),
_fixedAlignmentYRotationPlaneIds(),
_fixedAlignmentZRotationPlaneIds(),
_excludePlanesFromFit(),
_runPede(false),
_alignmentConstantLCIOFile("alignment.slcio"),
_maxMilleChi2Cut(1000.),
_tgeoFileName("TELESCOPE.root"),
_histoInfoFileName("histoinfo.xml"),
_trackCandidatesInputCollectionName("TrackCandidatesInputCollection"),
_trackFitter(0),
_milleGBL(0),
_seedAlignmentConstants(),
_nProcessedRuns(0),
_nProcessedEvents(0),
_flag_nohistos(false)
{

    // Processor description
    _description = "EUTelProcessorAlignmentGBLMille performs track fits using GBL optionally writing data files for MILLEPEDE II.";

    // TrackerHit input collection
    registerInputCollection(LCIO::TRACK,
            "TrackCandidatesInputCollectionName",
            "Input track candidates hits collection name",
            _trackCandidatesInputCollectionName,
            std::string("TrackCandidatesInputCollection"));


    // Necessary processor parameters that define fitter settings
    registerProcessorParameter("BeamEnergy", "Beam energy [GeV]", _eBeam, static_cast<double> (4.0));

    registerOptionalParameter("BeamCharge", "Beam charge [e]", _qBeam, static_cast<double> (-1));

    // Optional processor parameters that define finder settings

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

    registerOptionalParameter("GBLMEstimatorType", "GBL outlier down-weighting option (t,h,c)", _mEstimatorType, string() );
    
    registerOptionalParameter("MilleBinaryFilename", "Name of the Millepede binary file", _milleBinaryFilename, std::string("mille.bin"));

    registerOptionalParameter("MilleSteeringFilename", "Name of the Millepede steering file to be created", _milleSteeringFilename, std::string("pede-steer.txt"));
    
    registerOptionalParameter("MilleResultFilename", "Name of the Millepede result file", _milleResultFileName, std::string("millepede.res"));
    
    registerOptionalParameter("GearAlignedFile", "Suffix to add to the new Gear with alignment corrections", _gear_aligned_file, std::string("gear-00001-aligned.xml"));

    registerOptionalParameter("PedeSteeringAdditionalCmds","FOR EXPERTS: List of commands that should be included in the pede steering file. Use '\\' to seperate options and introduce a line break.",_pedeSteerAddCmds, StringVec());

    registerOptionalParameter("MilleMaxChi2Cut", "Maximum chi2 of a track candidate that goes into millepede", _maxMilleChi2Cut, double(1000.));

    registerOptionalParameter("AlignmentPlanes", "Ids of planes to be used in alignment", _alignmentPlaneIds, IntVec());
 
    //
    registerOptionalParameter("Planes", "Ids of planes to be used in alignment", _planeIds, IntVec());
    registerOptionalParameter("xResolutionPlane", "x resolution of planes given in AlignmentPlanes", _SteeringxResolutions, FloatVec());
    registerOptionalParameter("yResolutionPlane", "y resolution of planes given in AlignmentPlanes", _SteeringyResolutions, FloatVec());
   
    // Fixed planes parameters
    
    registerOptionalParameter("FixedAlignmentPlanesXshift", "Ids of planes for which X shift will be fixed during millepede call", _fixedAlignmentXShfitPlaneIds, IntVec());
    
    registerOptionalParameter("FixedAlignmentPlanesYshift", "Ids of planes for which Y shift will be fixed during millepede call", _fixedAlignmentYShfitPlaneIds, IntVec());
    
    registerOptionalParameter("FixedAlignmentPlanesZshift", "Ids of planes for which Z shift will be fixed during millepede call", _fixedAlignmentZShfitPlaneIds, IntVec());
    
    registerOptionalParameter("FixedAlignmentPlanesXrotation", "Ids of planes for which rotation around X will be fixed during millepede call", _fixedAlignmentXRotationPlaneIds, IntVec());
    
    registerOptionalParameter("FixedAlignmentPlanesYrotation", "Ids of planes for which rotation around Y will be fixed during millepede call", _fixedAlignmentYRotationPlaneIds, IntVec());
    
    registerOptionalParameter("FixedAlignmentPlanesZrotation", "Ids of planes for which rotation around Z will be fixed during millepede call", _fixedAlignmentZRotationPlaneIds, IntVec());
    
    registerOptionalParameter("ExcludePlanesFromFit", "Ids of planes that will be excluded from the fit", _excludePlanesFromFit, IntVec());
    
    // Pede run control
    
    registerOptionalParameter("RunPede","Execute the pede at the end of processing using the generated steering file.",_runPede, static_cast <bool> (false));
    
    registerOptionalParameter("AlignmentConstantLCIOFile","This is the name of the LCIO file name with the output alignment"
                            "constants (add .slcio)",_alignmentConstantLCIOFile, static_cast< string > ( "alignment.slcio" ) );
    
    //@TODO Implement geometry description
    // Geometry definition
    
    //registerOptionalParameter("GeometryFilename", "Name of the TGeo geometry definition file", _tgeoFileName, std::string("TELESCOPE.root"));
    
    // Histogram information

    registerOptionalParameter("HistogramInfoFilename", "Name of histogram info xml file", _histoInfoFileName, string("histoinfo.xml"));
}

void EUTelProcessorAlignmentGBLMille::init() {

    streamlog_out(DEBUG2) << "EUTelProcessorAlignmentGBLMille::init( )" << std::endl;

    // usually a good idea to
    printParameters();

    // Reset counters
    _nProcessedRuns = 0;
    _nProcessedEvents = 0;


    // Getting access to geometry description
    std::string name("test.root");
    geo::gGeometry().initializeTGeoDescription(name,false);

    // Instantiate millepede output. 
    {
        streamlog_out(DEBUG) << "Initialising Mille..." << std::endl;

        const unsigned int reserveSize = 80000;
        _milleGBL = new gbl::MilleBinary(_milleBinaryFilename, reserveSize);

        if (_milleGBL == 0) {
            streamlog_out(ERROR) << "Can't allocate an instance of mMilleBinary. Stopping ..." << std::endl;
            throw UnknownDataTypeException("MilleBinary was not created");
        }
    }

    // Instantiate track fitter. This is a working horse of the processor.
    {
        streamlog_out(DEBUG) << "Initialisation of track fitter" << std::endl;

        Utility::AlignmentMode alignmentMode = Utility::noAlignment;
        if (_alignmentMode==0) {
            alignmentMode = Utility::noAlignment;
        } else if (_alignmentMode==1) {
            alignmentMode = Utility::XYShift;
        } else if (_alignmentMode==2) {
            alignmentMode = Utility::XYShiftXYRot; 
        } else if (_alignmentMode==3) {
            alignmentMode = Utility::XYZShiftXYRot;
        } else if (_alignmentMode==4) {
            alignmentMode = Utility::XYShiftYZRotXYRot;
        } else if (_alignmentMode==5) {
            alignmentMode = Utility::XYShiftXZRotXYRot;
        } else if (_alignmentMode==6) {
            alignmentMode = Utility::XYShiftXZRotYZRotXYRot;
        } else if (_alignmentMode==7) {
            alignmentMode = Utility::XYZShiftXZRotYZRotXYRot;
        }else {
            streamlog_out(WARNING3) << "Alignment mode was not recognized:" << _alignmentMode << std::endl;
            streamlog_out(WARNING3) << "Alignment will not be performed" << std::endl;
            alignmentMode = Utility::noAlignment;
        }

        // fill MILLEPEDE alignment parameters labels
        fillMilleParametersLabels();
        
        // Initialize GBL fitter
        EUTelGBLFitter* Fitter = new EUTelGBLFitter("GBLFitter");
        Fitter->SetAlignmentMode(alignmentMode);
        Fitter->setParamterIdPlaneVec(_planeIds);
        Fitter->setParamterIdXResolutionVec(_SteeringxResolutions);
        Fitter->setParamterIdYResolutionVec(_SteeringyResolutions);

        Fitter->setExcludeFromFitPlanes( _excludePlanesFromFit );

        Fitter->setParamterIdXShiftsMap(_xShiftsMap);
        Fitter->setParamterIdYShiftsMap(_yShiftsMap);
        Fitter->setParamterIdZShiftsMap(_zShiftsMap);
        Fitter->setParamterIdXRotationsMap(_xRotationsMap);
        Fitter->setParamterIdYRotationsMap(_yRotationsMap);
        Fitter->setParamterIdZRotationsMap(_zRotationsMap);
        Fitter->SetMilleBinary(_milleGBL);
        Fitter->SetBeamEnergy(_eBeam);
        Fitter->SetBeamCharge(_qBeam);
        Fitter->setChi2Cut(_maxMilleChi2Cut);
        if (!_mEstimatorType.empty() ) Fitter->setMEstimatorType(_mEstimatorType);
        _trackFitter = Fitter;

        if (!_trackFitter) {
            streamlog_out(ERROR) << "Can't allocate an instance of EUTelGBLFitter. Stopping ..." << std::endl;
            throw UnknownDataTypeException("Track finder was not created");
        }
    }
}

void EUTelProcessorAlignmentGBLMille::processRunHeader(LCRunHeader * run) {

    auto_ptr<EUTelRunHeaderImpl> header(new EUTelRunHeaderImpl(run));
    header->addProcessor(type());


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

    // Flush seeding alignment constants
    {
        streamlog_out( DEBUG1 ) << "Flush alignment constants (X shifts)" << std::endl;
        for ( std::map<int, double >::iterator iDet = _seedAlignmentConstants._xResiduals.begin( );
                iDet != _seedAlignmentConstants._xResiduals.end( ); ++iDet ) {
            iDet->second = 0.;
        }

        for ( std::map<int, int >::iterator iDet = _seedAlignmentConstants._nxResiduals.begin( );
                iDet != _seedAlignmentConstants._nxResiduals.end( ); ++iDet ) {
            iDet->second = 1;
        }

        streamlog_out( DEBUG1 ) << "Flush alignment constants (Y shifts)" << std::endl;
        for ( std::map<int, double >::iterator iDet = _seedAlignmentConstants._yResiduals.begin( );
                iDet != _seedAlignmentConstants._yResiduals.end( ); ++iDet ) {
            iDet->second = 0.;
        }

        for ( std::map<int, int >::iterator iDet = _seedAlignmentConstants._nyResiduals.begin( );
                iDet != _seedAlignmentConstants._nyResiduals.end( ); ++iDet ) {
            iDet->second = 1;
        }
    }
    
    // Set millepede result file name according to run number
    const int runNumber = run->getRunNumber();
    
    const string defaultFileName = "millepede.res";
    if ( _milleResultFileName.compare(0,static_cast<int>(defaultFileName.size()),defaultFileName,0,static_cast<int>(defaultFileName.size())) == 0 ) {
        _milleResultFileName = "millepede-result-";
        _milleResultFileName += to_string(runNumber);
        _milleResultFileName += ".res";
    }
    _nProcessedRuns++;
}

void EUTelProcessorAlignmentGBLMille::processEvent(LCEvent * evt) {

    EUTelEventImpl * event = static_cast<EUTelEventImpl*> (evt);

    // Do not process last events
    if (event->getEventType() == kEORE) {
        streamlog_out(DEBUG4) << "EORE found: nothing else to do." << endl;
        return;
    } else if (event->getEventType() == kUNKNOWN) {
        streamlog_out(WARNING2) << "Event number " << event->getEventNumber() << " in run " << event->getRunNumber()
                << " is of unknown type. Continue considering it as a normal Data Event." << endl;
    }

    // Try to access collection

    LCCollection* col = NULL;
    try {
        col = evt->getCollection(_trackCandidatesInputCollectionName);
        streamlog_out(DEBUG1) << "collection : " << _trackCandidatesInputCollectionName << " retrieved" << std::endl;
    } catch (DataNotAvailableException e) {
        streamlog_out(MESSAGE0) << _trackCandidatesInputCollectionName << " collection not available" << std::endl;
        throw marlin::SkipEventException(this);
    }

    _trackFitter->Clear();

    // this will only be entered if the collection is available
    if (col != NULL) {
        streamlog_out(DEBUG2) << "EUTelProcessorTrackingGBLTrajectory" << endl;

        vector< const IMPL::TrackImpl* > trackCandidates;
        for (int iCol = 0; iCol < col->getNumberOfElements(); iCol++) {
            const IMPL::TrackImpl* trackimpl = static_cast<const IMPL::TrackImpl*> (col->getElementAt(iCol));

            if (!col) {
                streamlog_out(WARNING2) << "EUTelLCObjectTrackCandidate collection not found found for event " << _nProcessedEvents <<
                        " in run " << _nProcessedRuns << endl;
                throw SkipEventException(this);
            }

            streamlog_out(DEBUG1) << "Track " << iCol << " nhits " << trackimpl->getTrackerHits().size() << endl;
            trackCandidates.push_back( trackimpl );
        } //for ( int iCol = 0; iCol < col->getNumberOfElements() ; iCol++ )

        // Perform fit for all found track candidates
        // ------------------------------------------
       
        const int nTracks = trackCandidates.size();
        streamlog_out(DEBUG1) << "N tracks found " << nTracks << endl;


// read from LCIO collection into a _trackFitter vector
            _trackFitter->SetTrackCandidates(trackCandidates);

// create and fill GBL trajectory
// by the flag _alignmentMode it should either fille GBL trajectory with Measurement and Scattering points for the Fit
// or with Mesaurement and global derivatives points for Mille  
            _trackFitter->TrackCandidatesToGBLTrajectories();
 
            _trackFitter->PerformMille();
    }

    _nProcessedEvents++;

    if (isFirstEvent()) _isFirstEvent = false;
}

void EUTelProcessorAlignmentGBLMille::check(LCEvent * /* evt */) {
    // nothing to check here
}

void EUTelProcessorAlignmentGBLMille::end() {
    
    // Free file resource before running pede exe
    delete _milleGBL;
    
    writeMilleSteeringFile();

    streamlog_out(DEBUG) << "EUTelProcessorAlignmentGBLMille::end()  " << name()
            << " processed " << _nProcessedEvents << " events in " << _nProcessedRuns << " runs "
            << std::endl;
    
    if ( _runPede ) runPede();

    delete _trackFitter;
}

void EUTelProcessorAlignmentGBLMille::runPede() {
    // check if alignment was requested
    if ( _alignmentMode == Utility::noAlignment ) {
        streamlog_out( WARNING1 ) << "RunPede was required, but alignment mode is noAlignment. Stop." << endl;
        return;
    }

    std::string command = "pede " + _milleSteeringFilename;
    streamlog_out ( MESSAGE5 ) << "Starting pede...: " << command.c_str( ) << endl;

    // run pede and create a streambuf that reads its stdout and stderr
    redi::ipstream pede( command.c_str( ), redi::pstreams::pstdout | redi::pstreams::pstderr );

    if ( !pede.is_open( ) ) {
        streamlog_out( ERROR5 ) << "Pede cannot be executed: command not found in the path" << endl;
    } else {

        // Currently unused variable:
        //bool encounteredError = false;

        // output multiplexing: parse pede output in both stdout and stderr and echo messages accordingly
        char buf[1024];
        std::streamsize n;
        std::stringstream pedeoutput; // store stdout to parse later
        std::stringstream pedeerrors;
        bool finished[2] = { false, false };
        while ( !finished[0] || !finished[1] ) {
            if ( !finished[0] ) {
                while ( ( n = pede.err( ).readsome( buf, sizeof (buf ) ) ) > 0 ) {
                    streamlog_out( ERROR5 ).write( buf, n ).flush( );
                    string error ( buf, n );
                    pedeerrors << error;
                    //encounteredError = true;
                }
                if ( pede.eof( ) ) {
                    finished[0] = true;
                    if ( !finished[1] )
                        pede.clear( );
                }
            }

            if ( !finished[1] ) {
                while ( ( n = pede.out( ).readsome( buf, sizeof (buf ) ) ) > 0 ) {
                    streamlog_out( MESSAGE4 ).write( buf, n ).flush( );
                    string output ( buf, n );
                    pedeoutput << output;
                }
                if ( pede.eof( ) ) {
                    finished[1] = true;
                    if ( !finished[0] )
                        pede.clear( );
                }
            }
        }
        // wait for the pede execution to finish
        pede.close( );
        
        // Parse and rename MILLEPEDE result file
        if ( parseMilleOutput( "millepede.res" ) ) moveMilleResultFile( "millepede.res", _milleResultFileName );
    }
}

bool EUTelProcessorAlignmentGBLMille::parseMilleOutput( const string& milleResultFileName ) {
    
    bool isOK = true;
    
    // Check if the file is avaliable
    ifstream file( milleResultFileName.c_str() );
    if ( !file.good( ) ) {
        streamlog_out( WARNING2 ) << "Can't read/find " << milleResultFileName << " in current directory." << endl;
        isOK = false;
        return isOK;
    }
 
   
    const string command = "parsemilleout.sh " + _milleSteeringFilename + " " + milleResultFileName + " " + _alignmentConstantLCIOFile + 
                           " " + Global::parameters->getStringVal("GearXMLFile" ) + " " + _gear_aligned_file;
    streamlog_out ( MESSAGE5 ) << "Convering millepede results to LCIO collections... " << endl;
    streamlog_out ( MESSAGE5 ) << command << endl;

    // run pede and create a streambuf that reads its stdout and stderr
    redi::ipstream parsepede( command.c_str( ), redi::pstreams::pstdout | redi::pstreams::pstderr );

    if ( !parsepede.is_open( ) ) {
        streamlog_out( ERROR5 ) << "Pede cannot be executed: command not found in the path" << endl;
    } else {
        // Currently unused variable:
        // bool encounteredError = false;
        // output multiplexing: parse parsepede output in both stdout and stderr and echo messages accordingly
        char buf[1024];
        std::streamsize n;
        std::stringstream parsepedeoutput; // store stdout to parse later
        std::stringstream parsepedeerrors;
        bool finished[2] = { false, false };
        while ( !finished[0] || !finished[1] ) {
            if ( !finished[0] ) {
                while ( ( n = parsepede.err( ).readsome( buf, sizeof (buf ) ) ) > 0 ) {
                    streamlog_out( ERROR5 ).write( buf, n ).flush( );
                    string error ( buf, n );
                    parsepedeerrors << error;
                    //encounteredError = true;
                }
                if ( parsepede.eof( ) ) {
                    finished[0] = true;
                    if ( !finished[1] )
                        parsepede.clear( );
                }
            }

            if ( !finished[1] ) {
                while ( ( n = parsepede.out( ).readsome( buf, sizeof (buf ) ) ) > 0 ) {
                    streamlog_out( MESSAGE4 ).write( buf, n ).flush( );
                    string output ( buf, n );
                    parsepedeoutput << output;
                }
                if ( parsepede.eof( ) ) {
                    finished[1] = true;
                    if ( !finished[0] )
                        parsepede.clear( );
                }
            }
        }
        // wait for the parsepede execution to finish
        parsepede.close( );
    }
    
    return isOK;
    
}

void EUTelProcessorAlignmentGBLMille::moveMilleResultFile( const string& oldMilleResultFileName, const string& newMilleResultFileName ) {
    
    // Looking for oldMilleResultFileName in current folder   
    // Check if the file is avaliable
    ifstream infile( oldMilleResultFileName.c_str() );
    if ( !infile.good( ) ) {
        streamlog_out( WARNING2 ) << "Can't read/find " << oldMilleResultFileName << " in current directory." << endl;
        streamlog_out( WARNING2 ) << "Probably MILLEPEDE did not converged " << endl;
        return;
    }
    
    // Check if the destination file already exists
    ifstream outfile( newMilleResultFileName.c_str() );
    if ( outfile.good( ) ) {
        streamlog_out( WARNING2 ) << newMilleResultFileName << " exists in current directory." << endl;
        streamlog_out( WARNING2 ) << newMilleResultFileName << " will not be renamed." << endl;
        return;
    }
    
    // If file was found in current folder rename it
    const int result = rename( oldMilleResultFileName.c_str() , newMilleResultFileName.c_str() );
    if ( result == 0 ) {
        streamlog_out( MESSAGE4 ) << "File " << oldMilleResultFileName << " was renamed to " << newMilleResultFileName << endl;
    }
    else {
        streamlog_out( ERROR1 ) << "Error renaming file " << oldMilleResultFileName << endl;
    }
    
}

void EUTelProcessorAlignmentGBLMille::writeMilleSteeringFile() {

    streamlog_out(DEBUG2) << "writeMilleSteeringFile" << endl;

    // Prepare millepede steering files only if alignment was requested
    if (_alignmentMode == Utility::noAlignment) {
        streamlog_out(WARNING1) << "Alignment steering file will not be created" << endl;
        return;
    }

    ofstream steerFile;
    steerFile.open(_milleSteeringFilename.c_str());

    if (!steerFile.is_open()) {
        streamlog_out(ERROR2) << "Could not open steering file." << _milleSteeringFilename << endl;
        return;
    }

    steerFile << "Cfiles" << endl;
    steerFile << _milleBinaryFilename << endl;
    steerFile << endl;
    //
    steerFile << "Parameter" << endl;

    int counter = 0;

    EUTelGBLFitter* fitter = dynamic_cast < EUTelGBLFitter* > ( _trackFitter );
    std::map<int, int> XShiftsMap = fitter->getParamterIdXShiftsMap();
    std::map<int, int> YShiftsMap = fitter->getParamterIdYShiftsMap();
    std::map<int, int> ZShiftsMap = fitter->getParamterIdZShiftsMap();
    std::map<int, int> XRotationsMap = fitter->getParamterIdXRotationsMap();
    std::map<int, int> YRotationsMap = fitter->getParamterIdYRotationsMap();
    std::map<int, int> ZRotationsMap = fitter->getParamterIdZRotationsMap();
    
    // loop over all planes
    // @TODO assumes that planes have ids 0..._nplanes !generaly wrong
    for (unsigned int help = 0; help < geo::gGeometry().nPlanes(); help++) {

        const int sensorId = geo::gGeometry().sensorZOrderToID(help);
        const bool isPlaneExcluded = std::find(_alignmentPlaneIds.begin(), _alignmentPlaneIds.end(), sensorId) == _alignmentPlaneIds.end();
        
        // check if plane has to be used as fixed
        const bool isFixedXShift = std::find(_fixedAlignmentXShfitPlaneIds.begin(), _fixedAlignmentXShfitPlaneIds.end(), sensorId) != _fixedAlignmentXShfitPlaneIds.end();
        const bool isFixedYShift = std::find(_fixedAlignmentYShfitPlaneIds.begin(), _fixedAlignmentYShfitPlaneIds.end(), sensorId) != _fixedAlignmentYShfitPlaneIds.end();
        const bool isFixedZShift = std::find(_fixedAlignmentZShfitPlaneIds.begin(), _fixedAlignmentZShfitPlaneIds.end(), sensorId) != _fixedAlignmentZShfitPlaneIds.end();
        const bool isFixedXRotation = std::find(_fixedAlignmentXRotationPlaneIds.begin(), _fixedAlignmentXRotationPlaneIds.end(), sensorId) != _fixedAlignmentXRotationPlaneIds.end();
        const bool isFixedYRotation = std::find(_fixedAlignmentYRotationPlaneIds.begin(), _fixedAlignmentYRotationPlaneIds.end(), sensorId) != _fixedAlignmentYRotationPlaneIds.end();
        const bool isFixedZRotation = std::find(_fixedAlignmentZRotationPlaneIds.begin(), _fixedAlignmentZRotationPlaneIds.end(), sensorId) != _fixedAlignmentZRotationPlaneIds.end();
        
        // if plane not excluded
        if ( !isPlaneExcluded ) {

            const string initUncertaintyXShift = (isFixedXShift) ? "-1." : "0.01";
            const string initUncertaintyYShift = (isFixedYShift) ? "-1." : "0.01";
            const string initUncertaintyZShift = (isFixedZShift) ? "-1." : "0.01";
            const string initUncertaintyXRotation = (isFixedXRotation) ? "-1." : "0.01";
            const string initUncertaintyYRotation = (isFixedYRotation) ? "-1." : "0.01";
            const string initUncertaintyZRotation = (isFixedZRotation) ? "-1." : "0.01";
            
            const double initXshift = (isFixedXShift) ? 0. : _seedAlignmentConstants._xResiduals[sensorId]/_seedAlignmentConstants._nxResiduals[sensorId];
            const double initYshift = (isFixedYShift) ? 0. : _seedAlignmentConstants._yResiduals[sensorId]/_seedAlignmentConstants._nyResiduals[sensorId];
            
            if( fitter->GetAlignmentMode()==Utility::XYZShiftXYRot ) {
                steerFile << left << setw(25) << XShiftsMap[sensorId] << setw(25) << -initXshift << setw(25) << initUncertaintyXShift
                           << setw(25) << " ! X shift " << setw(25) << sensorId << endl;
                steerFile << left << setw(25) << YShiftsMap[sensorId] << setw(25) << -initYshift << setw(25) << initUncertaintyYShift
                           << setw(25) << " ! Y shift " << setw(25) << sensorId << endl;
                steerFile << left << setw(25) << ZShiftsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyZShift
                           << setw(25) << " ! Z shift " << setw(25) << sensorId << endl;
                steerFile << left << setw(25) << ZRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyZRotation
                          << setw(25) << " ! XY rotation " << sensorId << endl;
            } else if( fitter->GetAlignmentMode()==Utility::XYShiftYZRotXYRot ) {
                steerFile << left << setw(25) << XShiftsMap[sensorId] << setw(25) << -initXshift << setw(25) << initUncertaintyXShift
                          << setw(25) << " ! X shift " << sensorId << endl;
                steerFile << left << setw(25) << YShiftsMap[sensorId] << setw(25)  << -initYshift << setw(25) << initUncertaintyYShift
                          << setw(25) << " ! Y shift " << sensorId << endl;
                steerFile << left << setw(25) << XRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyXRotation
                          << setw(25) << " ! YZ rotation " << sensorId << endl;
                steerFile << left << setw(25) << ZRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyZRotation
                          << setw(25) << " ! XY rotation " << sensorId << endl;
            } else if( fitter->GetAlignmentMode()==Utility::XYShiftXZRotXYRot ) {
                steerFile << left << setw(25) << XShiftsMap[sensorId] << setw(25) << -initXshift << setw(25) << initUncertaintyXShift
                          << setw(25) << " ! X shift " << sensorId << endl;
                steerFile << left << setw(25) << YShiftsMap[sensorId] << setw(25) << -initYshift << setw(25) << initUncertaintyYShift
                          << setw(25) << " ! Y shift " << sensorId << endl;
                steerFile << left << setw(25) << YRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyYRotation
                          << setw(25) << " ! XZ rotation " << sensorId << endl;
                steerFile << left << setw(25) << ZRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyZRotation
                         << setw(25)  << " ! XY rotation " << sensorId << endl;
            } else if( fitter->GetAlignmentMode()==Utility::XYShiftXZRotYZRotXYRot ) {
                steerFile << left << setw(25) << XShiftsMap[sensorId] << setw(25) << -initXshift << setw(25) << initUncertaintyXShift
                          << setw(25) << " ! X shift " << sensorId << endl;
                steerFile << left << setw(25) << YShiftsMap[sensorId] << setw(25) << -initYshift << setw(25) << initUncertaintyYShift
                          << setw(25) << " ! Y shift " << sensorId << endl;
                steerFile << left << setw(25) << YRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyYRotation
                          << setw(25) << " ! XZ rotation " << sensorId << endl;
                steerFile << left << setw(25) << XRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyXRotation
                          << setw(25) << " ! YZ rotation " << sensorId << endl;
                steerFile << left << setw(25) << ZRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyZRotation
                         << setw(25)  << " ! XY rotation " << sensorId << endl;
            } else if( fitter->GetAlignmentMode()==Utility::XYZShiftXZRotYZRotXYRot ) {
                steerFile << left << setw(25) << XShiftsMap[sensorId] << setw(25) << -initXshift << setw(25) << initUncertaintyXShift
                          << setw(25) << " ! X shift " << sensorId << endl;
                steerFile << left << setw(25) << YShiftsMap[sensorId] << setw(25) << -initYshift << setw(25) << initUncertaintyYShift
                          << setw(25) << " ! Y shift " << sensorId << endl;
                steerFile << left << setw(25) << ZShiftsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyZShift
                          << setw(25) << " ! Z shift " << sensorId << endl;
                steerFile << left << setw(25) << YRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyYRotation
                          << setw(25) << " ! XZ rotation " << sensorId << endl;
                steerFile << left << setw(25) << XRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyXRotation
                          << setw(25) << " ! YZ rotation " << sensorId << endl;
                steerFile << left << setw(25) << ZRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyZRotation
                         << setw(25)  << " ! XY rotation " << sensorId << endl;
            } else if ( fitter->GetAlignmentMode()==Utility::XYShiftXYRot ) {
                steerFile << left << setw(25) << XShiftsMap[sensorId] << setw(25) << -initXshift << setw(25) << initUncertaintyXShift
                          << setw(25) << " ! X shift " << sensorId << endl;
                steerFile << left << setw(25) << YShiftsMap[sensorId] << setw(25) << -initYshift << setw(25) << initUncertaintyYShift
                          << setw(25) << " ! Y shift " << sensorId << endl;
                steerFile << left << setw(25) << ZRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyZRotation
                          << setw(25) << " ! XY rotation " << sensorId << endl;
            } else if ( fitter->GetAlignmentMode()==Utility::XYShift ) {
                steerFile << left << setw(25) << XShiftsMap[sensorId] << setw(25) << -initXshift << setw(25) << initUncertaintyXShift
                          << setw(25) << " ! X shift " << sensorId << endl;
                steerFile << left << setw(25) << YShiftsMap[sensorId] << setw(25) << -initYshift << setw(25) << initUncertaintyYShift
                          << setw(25) << " ! Y shift " << sensorId << endl;
                steerFile << left << setw(25) << ZRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << "-1.0"
                          << setw(25) << " ! XY rotation fixed" << sensorId << endl;
            }  

            counter++;

        } // end if plane not excluded

    } // end loop over all planes

//    steerFile << "method diagonalization 15 0.1" << endl;
//    steerFile << "hugecut 500." << endl;
//    steerFile << "!chiscut 50. 25." << endl;
//    steerFile << "outlierdownweighting 4" << endl;
//    steerFile << "dwfractioncut 0.2" << endl;

    steerFile << endl;
    for ( StringVec::iterator it = _pedeSteerAddCmds.begin( ); it != _pedeSteerAddCmds.end( ); ++it ) {
        // two backslashes will be interpreted as newline
        if ( *it == "\\\\" )
            steerFile << endl;
        else
            steerFile << *it << " ";
    }
    steerFile << endl;
    steerFile << "end" << endl;

    steerFile.close();

    if( _alignmentMode != Utility::noAlignment ) streamlog_out(MESSAGE5) << "File " << _milleSteeringFilename << " written." << endl;

}

void EUTelProcessorAlignmentGBLMille::fillMilleParametersLabels() {

    int currentLabel = 0;
    const IntVec sensorIDsVec = geo::gGeometry().sensorIDsVec();
    IntVec::const_iterator itr;
    for( itr = sensorIDsVec.begin(); itr != sensorIDsVec.end(); ++itr ) {
        _xShiftsMap.insert( make_pair(*itr, ++currentLabel) );
    }
    for( itr = sensorIDsVec.begin(); itr != sensorIDsVec.end(); ++itr ) {
        _yShiftsMap.insert( make_pair(*itr, ++currentLabel) );
    }
    for( itr = sensorIDsVec.begin(); itr != sensorIDsVec.end(); ++itr ) {
        _zShiftsMap.insert( make_pair(*itr, ++currentLabel) );
    }
    for( itr = sensorIDsVec.begin(); itr != sensorIDsVec.end(); ++itr ) {
        _xRotationsMap.insert( make_pair(*itr, ++currentLabel) );
    }
    for( itr = sensorIDsVec.begin(); itr != sensorIDsVec.end(); ++itr ) {
        _yRotationsMap.insert( make_pair(*itr, ++currentLabel) );
    }
    for( itr = sensorIDsVec.begin(); itr != sensorIDsVec.end(); ++itr ) {
        _zRotationsMap.insert( make_pair(*itr, ++currentLabel) );
    }
}

#endif // USE_GBL

