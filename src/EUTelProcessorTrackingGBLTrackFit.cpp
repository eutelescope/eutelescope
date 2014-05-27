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

// AIDA
#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IProfile2D.h>
#endif // MARLIN_USE_AIDA

// ROOT
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TMath.h"
#endif

// EUTELESCOPE
#include "EUTelProcessorTrackingGBLTrackFit.h"

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

/**  EUTelProcessorTrackingGBLTrackFit
 * 
 *  If compiled with MARLIN_USE_AIDA 
 *  it creates histograms.
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Requires collection of track candidates hits.
 *
 *  <h4>Output</h4> 
 *  <li> Histograms.
 *  <li> Collection of fitted tracks.
 */

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

std::string EUTelProcessorTrackingGBLTrackFit::_histName::_orchi2GblFitHistName    = "orchi2GblFit";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_orchi2ndfGblFitHistName = "orchi2ndfGblFit";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_orprobGblFitHistName    = "orprobGblFit";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_chi2GblFitHistName = "chi2GblFit";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_chi2ndfGblFitHistName = "chi2ndfGblFit";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_probGblFitHistName = "probGblFit";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_momentumGblFitHistName = "momentumGblFit";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_residGblFitHistName = "ResidualsGblFit";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_normResidGblFitHistName = "NormResidualsGblFit";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_residGblFitHistNameX = "ResidualsGblFit_x";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_residGblFitHistNameY = "ResidualsGblFit_y";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_resid2DGblFitHistNameX = "ResidualsGblFitXY_x";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_resid2DGblFitHistNameY = "ResidualsGblFitXY_y";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_normResidGblFitHistNameX = "NormResidualsGblFit_x";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_normResidGblFitHistNameY = "NormResidualsGblFit_y";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_resid2DGblFitHistNameXvsX = "Residuals2DGblFit_xVSx";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_resid2DGblFitHistNameYvsX = "Residuals2DGblFit_yVSx";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_resid2DGblFitHistNameXvsY = "Residuals2DGblFit_xVSy";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_resid2DGblFitHistNameYvsY = "Residuals2DGblFit_yVSy";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_normResid2DGblFitHistNameXvsX = "NormResiduals2DGblFit_xVSx";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_normResid2DGblFitHistNameYvsX = "NormResiduals2DGblFit_yVSx";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_normResid2DGblFitHistNameXvsY = "NormResiduals2DGblFit_xVSy";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_normResid2DGblFitHistNameYvsY = "NormResiduals2DGblFit_yVSy";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_kinkGblFitHistNameX = "KinksGblFit_x";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_kinkGblFitHistNameY = "KinksGblFit_y";
#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

/** Default constructor */
EUTelProcessorTrackingGBLTrackFit::EUTelProcessorTrackingGBLTrackFit() :
Processor("EUTelProcessorTrackingGBLTrackFit"),
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
_trackCandidateHitsInputCollectionName("TrackCandidateHitCollection"),
_tracksOutputCollectionName("TrackCollection"),
_trackFitter(0),
_milleGBL(0),
_seedAlignmentConstants(),
_nProcessedRuns(0),
_nProcessedEvents(0),
_flag_nohistos(false),
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
_aidaHistoMap1D(),
_aidaHistoMap2D(),
_aidaProfileMap2D()
#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
{

    // Processor description
    _description = "EUTelProcessorTrackingGBLTrackFit performs track fits using GBL optionally writing data files for MILLEPEDE II.";

    // TrackerHit input collection
    registerInputCollection(LCIO::TRACK,
            "TrackCandHitInputCollectionName",
            "Input track candidates hits collection name",
            _trackCandidateHitsInputCollectionName,
            std::string("TrackCandidateHitCollection"));

    // Track output collection
    registerOutputCollection(LCIO::TRACK,
            "TracksOutputCollectionName",
            "Output tracks collection name",
            _tracksOutputCollectionName,
            std::string("TrackCollection"));

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

void EUTelProcessorTrackingGBLTrackFit::init() {

    streamlog_out(DEBUG2) << "EUTelProcessorTrackingGBLTrackFit::init( )" << std::endl;

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
        Fitter->SetChi2Cut(_maxMilleChi2Cut);
        if (!_mEstimatorType.empty() ) Fitter->setMEstimatorType(_mEstimatorType);
        _trackFitter = Fitter;

        if (!_trackFitter) {
            streamlog_out(ERROR) << "Can't allocate an instance of EUTelGBLFitter. Stopping ..." << std::endl;
            throw UnknownDataTypeException("Track finder was not created");
        }
    }
    // Book histograms
    bookHistograms();
}

void EUTelProcessorTrackingGBLTrackFit::processRunHeader(LCRunHeader * run) {

    auto_ptr<EUTelRunHeaderImpl> header(new EUTelRunHeaderImpl(run));
    header->addProcessor(type());


    // this is the right place also to check the geometry ID. This is a
    // unique number identifying each different geometry used at the
    // beam test. The same number should be saved in the run header and
    // in the xml file. If the numbers are different, warn the user.

    if (header->getGeoID() == 0)
        streamlog_out(WARNING0) << "The geometry ID in the run header is set to zero." << endl
            << "This may mean that the GeoID parameter was not set" << endl;


    if (header->getGeoID() != geo::gGeometry()._siPlanesParameters->getSiPlanesID()) {
        streamlog_out(WARNING5) << "Error during the geometry consistency check: " << endl
                << "The run header says the GeoID is " << header->getGeoID() << endl
                << "The GEAR description says is     " << geo::gGeometry()._siPlanesParameters->getSiPlanesID() << endl;
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

void EUTelProcessorTrackingGBLTrackFit::processEvent(LCEvent * evt) {

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
        col = evt->getCollection(_trackCandidateHitsInputCollectionName);
        streamlog_out(DEBUG1) << "collection : " << _trackCandidateHitsInputCollectionName << " retrieved" << std::endl;
    } catch (DataNotAvailableException e) {
        streamlog_out(MESSAGE0) << _trackCandidateHitsInputCollectionName << " collection not available" << std::endl;
        throw marlin::SkipEventException(this);
    }

    // this will only be entered if the collection is available
    if (col != NULL) {
        streamlog_out(DEBUG2) << "EUTelProcessorTrackingGBLTrackFit" << endl;
// trackCanidadtes type !!! rubinskiy 30-11-2013
//        EVENT::TrackVec trackCandidates;
        EVENT::TrackVec trackCandidates;
//        vector < EVENT::TrackerHitVec > trackCandidates;
        for (int iCol = 0; iCol < col->getNumberOfElements(); iCol++) {
            TrackImpl* track = static_cast<TrackImpl*> (col->getElementAt(iCol));

            if (!col) {
                streamlog_out(WARNING2) << "EUTelLCObjectTrackCandidate collection not found found for event " << _nProcessedEvents <<
                        " in run " << _nProcessedRuns << endl;
                throw SkipEventException(this);
            }

            streamlog_out(DEBUG1) << "Track " << iCol << " nhits " << track->getTrackerHits().size() << endl;
            trackCandidates.push_back( track );
//            trackCandidates.push_back(track->getTrackerHits());
        } //for ( int iCol = 0; iCol < col->getNumberOfElements() ; iCol++ )

        //        return;

        // Perform fit for all found track candidates
        // ------------------------------------------
        TVectorD residual(2);
        TVectorD residualErr(2);
        
        unsigned int numData;
        TVectorD residualGBL(2);
        TVectorD measErrGBL(2);
        TVectorD residualErrGBL(2);
        TVectorD downWeightGBL(2);
        const int nTracks = trackCandidates.size();
        streamlog_out(DEBUG1) << "N tracks found " << nTracks << endl;

        // for alignment nTracks == 1, for tracking nTracks>0
        if (  nTracks > 0 ) {
//       if ( _alignmentMode ? nTracks == 1 : nTracks > 0 ) {
            _trackFitter->SetTrackCandidates(trackCandidates);
            _trackFitter->FitTracks();
//             _trackFitter->FitSingleTrackCandidate();

            IMPL::LCCollectionVec* fittrackvec;
            fittrackvec = static_cast<EUTelGBLFitter*> (_trackFitter)->GetFitTrackVec();
            
            if( fittrackvec == 0 )  {
              streamlog_out( MESSAGE0 ) << " fittrackvec is null " << std::endl;  
              return;
            }

            if( fittrackvec->size() == 0 )  {
              streamlog_out( MESSAGE0 ) << " fittrackvec is empty " << std::endl;  
              return;
            }

            //
            if( !_flag_nohistos && 1==1 ) {

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
            // build plots: // makes sense to put into a separate method, aah ?
            //
            double chi2Trk = 0.;
            int ndfTrk = 0;
            double p = 0.;

            IMPL::LCCollectionVec::const_iterator itFitTrack;


            // Loop over fitted tracks
            for (itFitTrack = fittrackvec->begin(); itFitTrack != fittrackvec->end(); ++itFitTrack) {

                chi2Trk =static_cast<TrackImpl*> (*itFitTrack)->getChi2();
                ndfTrk = static_cast<TrackImpl*> (*itFitTrack)->getNdf();
                static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ _histName::_orchi2GblFitHistName ]) -> fill(chi2Trk);
                if(ndfTrk>0) static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ _histName::_orchi2ndfGblFitHistName ]) -> fill(chi2Trk/ndfTrk);
                static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ _histName::_orprobGblFitHistName ]) -> fill( TMath::Prob(chi2Trk, ndfTrk) );
                
                p = _qBeam * ( 1./static_cast<TrackImpl*> (*itFitTrack)->getOmega() );
                static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ _histName::_momentumGblFitHistName ]) -> fill(p);


                // Retrieve original GBL information
                std::map< int, gbl::GblTrajectory* > gblTracks = static_cast < EUTelGBLFitter* > ( _trackFitter )->GetGblTrackCandidates( );
                IMPL::LCCollectionVec::const_iterator begin = fittrackvec->begin( );
                int iCounter = std::distance( begin, itFitTrack );
                gbl::GblTrajectory* gblTraj = ( *gblTracks.find( iCounter ) ).second;
                if ( gblTraj == NULL ) {
                    streamlog_out( WARNING1 ) << "Can't find GBL trajectory object #"<< iCounter << " of " << fittrackvec->size() <<". Skipping track" << std::endl;
                    continue;
                }else{
                    streamlog_out( MESSAGE1 ) << "found GBL trajectory object #"<< iCounter << " of " << fittrackvec->size() <<". chi2:" << chi2Trk << " max:" << _maxMilleChi2Cut << std::endl;
                }
 

                const std::map<long, int> gblPointLabel = static_cast < EUTelGBLFitter* > ( _trackFitter )->getHitId2GblPointLabel( );
                const EVENT::TrackerHitVec& trackHits = static_cast < TrackImpl* > ( *itFitTrack )->getTrackerHits( );

                // If this is an alignment run plot residuals for Millepede selected tracks. Plot everything otherwise.
//                if ( ( _alignmentMode > 0 && chi2Trk < _maxMilleChi2Cut ) || ( _alignmentMode == 0 ) ) {
                if ( chi2Trk < _maxMilleChi2Cut   ) {
 
                    static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ _histName::_chi2GblFitHistName ]) -> fill(chi2Trk);
                    if(ndfTrk>0) static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ _histName::_chi2ndfGblFitHistName ]) -> fill(chi2Trk/ndfTrk);
                    static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ _histName::_probGblFitHistName ]) -> fill( TMath::Prob(chi2Trk, ndfTrk) );

 
                    streamlog_out(DEBUG1) << "Event number " << event->getEventNumber() << " in run " << event->getRunNumber()
                     << "event: " 
                     << "fittrack: " << (static_cast<TrackImpl*> (*itFitTrack))->id() 
                     << " fittrackvec->end(): " << fittrackvec->size()
                     << " trackHits.end(): " << trackHits.size() << endl ; 

                    // Loop over fitted hit positions
                    EVENT::TrackerHitVec::const_iterator itrHit;
                    for ( itrHit = trackHits.begin( ); itrHit != trackHits.end( ); ++itrHit ) {

                       // Get input hit information
                        IMPL::TrackerHitImpl* originalHit =  static_cast < IMPL::TrackerHitImpl* > ( ( *itrHit )->getRawHits( ).front( ) );
                        const double* originalhitpos = originalHit->getPosition( );

                        // Get fitted hit information
                        const double* hitpos = ( *itrHit )->getPosition( );
                        const int planeID = geo::gGeometry().getSensorID( originalHit );

                        bool excludeFromFit = false;
                        if ( std::find( _excludePlanesFromFit.begin(), _excludePlanesFromFit.end(), planeID ) != _excludePlanesFromFit.end() ) excludeFromFit = true;
 
                        int xSize = -1; 
                        int ySize = -1; 
                        Utility::getClusterSize(originalHit, xSize, ySize);
                        
//                        std::cout << "plane:" << planeID << " hittype: " << originalHit->getType()  <<  " x: " << xSize << " y: " << ySize << " " << std::endl;
                

                        if ( planeID < 0 ) continue;

                        // extract hit - cluster information - along X and Y   
                        

                        // Calculate residuals between original and fitted hits
                        residual[0] = hitpos[0] - originalhitpos[0];
                        residual[1] = hitpos[1] - originalhitpos[1];

                        // Covariance includes correlations from fit parameters
                        residualErr[0] = sqrt( ( *itrHit )->getCovMatrix( )[0] );
                        residualErr[1] = sqrt( ( *itrHit )->getCovMatrix( )[2] );

                        // Spatial GBL residuals. residualGBL is the as residual and residualErrGBL is the same as residualErr.
                        int hitGblLabel = gblPointLabel.at( originalHit->id( ) );
//                        cout<<"     hitGblLabel: " << hitGblLabel << " : " << numData << " " << " gblTraj : " << gblTraj << endl;
           if( !excludeFromFit)
            {
//                        gblTraj->printTrajectory(1);  
                        gblTraj->getMeasResults( hitGblLabel, numData, residualGBL, measErrGBL, residualErrGBL, downWeightGBL );
//cout<<"nach hitGblLabel: " << hitGblLabel << " : " << numData << endl;
            }
                        std::stringstream sstr;
                        std::stringstream sstrNorm;

                        // 1D histograms
                        // SPATIAL RESIDUALS
                        const double um = 1000.;
                        sstr << _histName::_residGblFitHistNameX << planeID;
                        sstrNorm << _histName::_normResidGblFitHistNameX << planeID;
                        if ( planeID == 5 ) streamlog_out( DEBUG0 ) << planeID << " " << std::setw( 15 ) << std::setprecision( 5 ) << residual[0] << std::setw( 15 ) << std::setprecision( 5 ) << residualErr[0] << std::endl;
                        static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ sstr.str( ) ] ) -> fill( residualGBL[0] * um, downWeightGBL[0] );
                        static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ sstrNorm.str( ) ] ) -> fill( residual[0] / residualErr[0], downWeightGBL[0] );

if( event->getEventNumber()/999*999   == event->getEventNumber() )
{
  streamlog_out(MESSAGE3) <<  sstr.str( ) 
            << " " << hitpos[0] << " " << hitpos[1] << " " << hitpos[2] 
            << " " <<  static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ sstr.str( ) ] ) -> mean()  
            << " " <<  static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ sstr.str( ) ] ) -> rms() << endl;
}



                        _seedAlignmentConstants._xResiduals[planeID] += ( residual[0] );
                        _seedAlignmentConstants._nxResiduals[planeID]++;
                        sstr.str( std::string( ) );
                        sstrNorm.str( std::string( ) );
                        sstr << _histName::_residGblFitHistNameY << planeID;
                        sstrNorm << _histName::_normResidGblFitHistNameY << planeID;
                        if ( planeID == 5 ) streamlog_out( DEBUG0 ) << planeID << " " << std::setw( 15 ) << std::setprecision( 5 ) << residual[1] << std::setw( 15 ) << std::setprecision( 5 ) << residualErr[1] << std::endl;
                        static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ sstr.str( ) ] ) -> fill( residual[1] * um, downWeightGBL[1] );
                        static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ sstrNorm.str( ) ] ) -> fill( residual[1] / residualErr[1], downWeightGBL[1] );
                        _seedAlignmentConstants._yResiduals[planeID] += ( residual[1] );
                        _seedAlignmentConstants._nyResiduals[planeID]++;
                        sstr.str( std::string( ) );
                        sstrNorm.str( std::string( ) );

                        // 2D histograms
                        // SPATIAL RESIDUALS
                        sstr << _histName::_resid2DGblFitHistNameX << planeID;
                        static_cast < AIDA::IProfile2D* > ( _aidaProfileMap2D[ sstr.str( ) ] ) -> fill( hitpos[0], hitpos[1], residual[0] * um, downWeightGBL[0] );
                        sstr.str( std::string( ) );

                        sstr << _histName::_resid2DGblFitHistNameY << planeID;
                        static_cast < AIDA::IProfile2D* > ( _aidaProfileMap2D[ sstr.str( ) ] ) -> fill( hitpos[0], hitpos[1], residual[1] * um, downWeightGBL[1] );
                        sstr.str( std::string( ) );

                        sstr << _histName::_resid2DGblFitHistNameXvsX << planeID;
                        sstrNorm << _histName::_normResid2DGblFitHistNameXvsX << planeID;
                        static_cast < AIDA::IHistogram2D* > ( _aidaHistoMap2D[ sstr.str( ) ] ) -> fill( hitpos[0], residual[0] * um, downWeightGBL[0] );
                        static_cast < AIDA::IHistogram2D* > ( _aidaHistoMap2D[ sstrNorm.str( ) ] ) -> fill( hitpos[0], residual[0] / residualErr[0], downWeightGBL[0] );
                        sstr.str( std::string( ) );
                        sstrNorm.str( std::string( ) );
                        sstr << _histName::_resid2DGblFitHistNameXvsY << planeID;
                        sstrNorm << _histName::_normResid2DGblFitHistNameXvsY << planeID;
                        static_cast < AIDA::IHistogram2D* > ( _aidaHistoMap2D[ sstr.str( ) ] ) -> fill( hitpos[1], residual[0] * um, downWeightGBL[0] );
                        static_cast < AIDA::IHistogram2D* > ( _aidaHistoMap2D[ sstrNorm.str( ) ] ) -> fill( hitpos[1], residual[0] / residualErr[0], downWeightGBL[0] );
                        sstr.str( std::string( ) );
                        sstrNorm.str( std::string( ) );
                        sstr << _histName::_resid2DGblFitHistNameYvsX << planeID;
                        sstrNorm << _histName::_normResid2DGblFitHistNameYvsX << planeID;
                        static_cast < AIDA::IHistogram2D* > ( _aidaHistoMap2D[ sstr.str( ) ] ) -> fill( hitpos[0], residual[1] * um, downWeightGBL[1] );
                        static_cast < AIDA::IHistogram2D* > ( _aidaHistoMap2D[ sstrNorm.str( ) ] ) -> fill( hitpos[0], residual[1] / residualErr[1], downWeightGBL[1] );
                        sstr.str( std::string( ) );
                        sstrNorm.str( std::string( ) );
                        sstr << _histName::_resid2DGblFitHistNameYvsY << planeID;
                        sstrNorm << _histName::_normResid2DGblFitHistNameYvsY << planeID;
                        static_cast < AIDA::IHistogram2D* > ( _aidaHistoMap2D[ sstr.str( ) ] ) -> fill( hitpos[1], residual[1] * um, downWeightGBL[1] );
                        static_cast < AIDA::IHistogram2D* > ( _aidaHistoMap2D[ sstrNorm.str( ) ] ) -> fill( hitpos[1], residual[1] / residualErr[1], downWeightGBL[1] );
                        sstr.str( std::string( ) );
                        sstrNorm.str( std::string( ) );

                        // 1D histograms
                        // KINKS
                        gblTraj->getScatResults( hitGblLabel, numData, residualGBL, measErrGBL, residualErrGBL, downWeightGBL );
                        sstr << _histName::_kinkGblFitHistNameX << planeID;
                        streamlog_out( DEBUG0 ) << std::setw( 15 ) << std::setprecision( 5 ) << residualGBL[0] << std::setw( 15 ) << std::setprecision( 5 ) << residualErr[0] << std::endl;
                        static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ sstr.str( ) ] ) -> fill( residualGBL[0], downWeightGBL[0] );
                        sstr.str( std::string( ) );
                        sstr << _histName::_kinkGblFitHistNameY << planeID;
                        streamlog_out( DEBUG0 ) << std::setw( 15 ) << std::setprecision( 5 ) << residualGBL[1] << std::setw( 15 ) << std::setprecision( 5 ) << residualErr[1] << std::endl;
                        static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ sstr.str( ) ] ) -> fill( residualGBL[1], downWeightGBL[1] );
                        sstr.str( std::string( ) );
                        sstrNorm.str( std::string( ) );

                   }
               } // if ( (_alignmentMode > 0 && chi2Trk < _maxMilleChi2Cut) || ( _alignmentMode == 0 ) )

            }
            // Write track candidates collection
            try {
                streamlog_out( DEBUG1 ) << "Getting collection " << _tracksOutputCollectionName << endl;
                evt->getCollection( _tracksOutputCollectionName );
            } catch ( ... ) {
                streamlog_out( DEBUG1 ) << "Adding collection " << _tracksOutputCollectionName << endl;
                evt->addCollection( static_cast < EUTelGBLFitter* > ( _trackFitter )->GetFitHitsVec( ), _tracksOutputCollectionName+"_fittedhits" );
                evt->addCollection( static_cast < EUTelGBLFitter* > ( _trackFitter )->GetFitTrackVec( ), _tracksOutputCollectionName );
            }   
        
    
        } //if( _alignmentMode ? nTracks == 1 : nTracks > 0 )

// building histograms done.
#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
        }
 
    } //if( col != NULL )

    _trackFitter->Clear();

    _nProcessedEvents++;

    if (isFirstEvent()) _isFirstEvent = false;
}

void EUTelProcessorTrackingGBLTrackFit::check(LCEvent * /* evt */) {
    // nothing to check here
}

void EUTelProcessorTrackingGBLTrackFit::end() {
    
    // Free file resource before running pede exe
    delete _milleGBL;
    
    writeMilleSteeringFile();

    streamlog_out(DEBUG) << "EUTelProcessorTrackingGBLTrackFit::end()  " << name()
            << " processed " << _nProcessedEvents << " events in " << _nProcessedRuns << " runs "
            << std::endl;
    
    if ( _runPede ) runPede();

    delete _trackFitter;
}

void EUTelProcessorTrackingGBLTrackFit::runPede() {
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

bool EUTelProcessorTrackingGBLTrackFit::parseMilleOutput( const string& milleResultFileName ) {
    
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

void EUTelProcessorTrackingGBLTrackFit::moveMilleResultFile( const string& oldMilleResultFileName, const string& newMilleResultFileName ) {
    
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

void EUTelProcessorTrackingGBLTrackFit::writeMilleSteeringFile() {

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

void EUTelProcessorTrackingGBLTrackFit::fillMilleParametersLabels() {

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

void EUTelProcessorTrackingGBLTrackFit::bookHistograms() {
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
//
        histoInfo = histoMgr->getHistogramInfo(_histName::_orchi2GblFitHistName);
        int    orchi2NBin =          ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 1000;    
        double orchi2Min =        ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : 0.;
        double orchi2Max =        ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 5000.;

        // GBL fits
        AIDA::IHistogram1D * orchi2GblFit =
                marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_orchi2GblFitHistName, orchi2NBin, orchi2Min, orchi2Max);
        if (orchi2GblFit) {
            orchi2GblFit->setTitle("#chi^{2} of track candidates; #chi^{2};N Tracks");
            _aidaHistoMap1D.insert(std::make_pair(_histName::_orchi2GblFitHistName, orchi2GblFit));
        } else {
            streamlog_out(ERROR2) << "Problem booking the " << (_histName::_orchi2GblFitHistName) << std::endl;
            streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
        }
 
        histoInfo = histoMgr->getHistogramInfo(_histName::_orchi2ndfGblFitHistName);

        int    orchi2ndfNBin =       ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 1000;    
        double orchi2ndfMin =        ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : 0.;
        double orchi2ndfMax =        ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 5000.;

        // GBL fits
        AIDA::IHistogram1D * orchi2ndfGblFit =
                marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_orchi2ndfGblFitHistName, orchi2ndfNBin, orchi2ndfMin, orchi2ndfMax);
        if (orchi2ndfGblFit) {
            orchi2ndfGblFit->setTitle("Normilised #chi^{2} of track candidates; #chi^{2}/ndf;N Tracks");
            _aidaHistoMap1D.insert(std::make_pair(_histName::_orchi2ndfGblFitHistName, orchi2ndfGblFit));
        } else {
            streamlog_out(ERROR2) << "Problem booking the " << (_histName::_orchi2ndfGblFitHistName) << std::endl;
            streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
        }
        
        histoInfo = histoMgr->getHistogramInfo(_histName::_orprobGblFitHistName);
        int    orprobNBin =          ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 1000;
        double orprobMin =        ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : 0.;
        double orprobMax =        ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 1.;

        AIDA::IHistogram1D * orprobGblFit =
                marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_orprobGblFitHistName, orprobNBin, orprobMin, orprobMax);
        if (orprobGblFit) {
            orprobGblFit->setTitle("Probability of track fit; Prob;N Tracks");
            _aidaHistoMap1D.insert(std::make_pair(_histName::_orprobGblFitHistName, orprobGblFit));
        } else {
            streamlog_out(ERROR2) << "Problem booking the " << (_histName::_orprobGblFitHistName) << std::endl;
            streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
        }


//
        histoInfo = histoMgr->getHistogramInfo(_histName::_chi2GblFitHistName);
        int chi2NBin =          ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 1000;    
        double chi2Min =        ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : 0.;
        double chi2Max =        ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 5000.;

        // GBL fits
        AIDA::IHistogram1D * chi2GblFit =
                marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_chi2GblFitHistName, chi2NBin, chi2Min, chi2Max);
        if (chi2GblFit) {
            chi2GblFit->setTitle("#chi^{2} of track candidates; #chi^{2};N Tracks");
            _aidaHistoMap1D.insert(std::make_pair(_histName::_chi2GblFitHistName, chi2GblFit));
        } else {
            streamlog_out(ERROR2) << "Problem booking the " << (_histName::_chi2GblFitHistName) << std::endl;
            streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
        }
 
        histoInfo = histoMgr->getHistogramInfo(_histName::_chi2ndfGblFitHistName);

        int    chi2ndfNBin =       ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 1000;    
        double chi2ndfMin =        ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : 0.;
        double chi2ndfMax =        ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 5000.;

        // GBL fits
        AIDA::IHistogram1D * chi2ndfGblFit =
                marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_chi2ndfGblFitHistName, chi2ndfNBin, chi2ndfMin, chi2ndfMax);
        if (chi2ndfGblFit) {
            chi2ndfGblFit->setTitle("Normilised #chi^{2} of track candidates; #chi^{2}/ndf;N Tracks");
            _aidaHistoMap1D.insert(std::make_pair(_histName::_chi2ndfGblFitHistName, chi2ndfGblFit));
        } else {
            streamlog_out(ERROR2) << "Problem booking the " << (_histName::_chi2ndfGblFitHistName) << std::endl;
            streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
        }
        
        histoInfo = histoMgr->getHistogramInfo(_histName::_probGblFitHistName);
        int    probNBin =          ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 1000;
        double probMin =        ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : 0.;
        double probMax =        ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 1.;

        AIDA::IHistogram1D * probGblFit =
                marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_probGblFitHistName, probNBin, probMin, probMax);
        if (probGblFit) {
            probGblFit->setTitle("Probability of track fit; Prob;N Tracks");
            _aidaHistoMap1D.insert(std::make_pair(_histName::_probGblFitHistName, probGblFit));
        } else {
            streamlog_out(ERROR2) << "Problem booking the " << (_histName::_probGblFitHistName) << std::endl;
            streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
        }

//        
        histoInfo = histoMgr->getHistogramInfo(_histName::_momentumGblFitHistName);
        int momNBin =          ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 1000;
        double momMin =        ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : 0.;
        double momMax =        ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 10.;
        
        AIDA::IHistogram1D * momentumGblFit =
                marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_momentumGblFitHistName, momNBin, momMin, momMax);
        if (momentumGblFit) {
            momentumGblFit->setTitle("Momentum of track; P [GeV/c];N Tracks");
            _aidaHistoMap1D.insert(std::make_pair(_histName::_momentumGblFitHistName, momentumGblFit));
        } else {
            streamlog_out(ERROR2) << "Problem booking the " << (_histName::_momentumGblFitHistName) << std::endl;
            streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
        }

        const int    resid1dNBinX = 1000;
        const double residMinX    = -20.;
        const double residMaxX    = 20.;
        
        int NBinX;
        double MinX;
        double MaxX;
        int NBinY;
        double MinY;
        double MaxY;
        double MinZ;
        double MaxZ;
        // Residuals after fit
        std::stringstream sstm;
        std::string residGblFitHistName;
        std::string histTitle;
        // normalised residuals y
        for (size_t iPlane = 0; iPlane < geo::gGeometry().nPlanes(); iPlane++) {
            sstm << _histName::_normResidGblFitHistNameX << geo::gGeometry().sensorIDsVec().at(iPlane);
            residGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Normalised residuals. Plane " << geo::gGeometry().sensorIDsVec().at(iPlane) << "X direction; Normalised residuals; N hits";
            histTitle = sstm.str();
            sstm.str(std::string(""));
            histoInfo = histoMgr->getHistogramInfo(residGblFitHistName);
            NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : resid1dNBinX;
            MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : residMinX;
            MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : residMaxX;
            AIDA::IHistogram1D * residGblFit =
                    marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(residGblFitHistName, NBinX, MinX, MaxX);
            if (residGblFit) {
                residGblFit->setTitle(histTitle);
                _aidaHistoMap1D.insert(std::make_pair(residGblFitHistName, residGblFit));
            } else {
                streamlog_out(ERROR2) << "Problem booking the " << (residGblFitHistName) << std::endl;
                streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            }
            sstm.str(std::string(""));
        }
        
        // plane residuals y
        for (size_t iPlane = 0; iPlane < geo::gGeometry().nPlanes(); iPlane++) {
            sstm << _histName::_residGblFitHistNameX << geo::gGeometry().sensorIDsVec().at(iPlane);
            residGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Residuals. Plane " << geo::gGeometry().sensorIDsVec().at(iPlane) << "X direction; residuals [um]; N hits";
            histTitle = sstm.str();
            sstm.str(std::string(""));
            histoInfo = histoMgr->getHistogramInfo(residGblFitHistName);
            NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : resid1dNBinX;
            MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : residMinX;
            MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : residMaxX;
            AIDA::IHistogram1D * residGblFit =
                    marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(residGblFitHistName, NBinX, MinX, MaxX);
            if (residGblFit) {
                residGblFit->setTitle(histTitle);
                _aidaHistoMap1D.insert(std::make_pair(residGblFitHistName, residGblFit));
            } else {
                streamlog_out(ERROR2) << "Problem booking the " << (residGblFitHistName) << std::endl;
                streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            }
            sstm.str(std::string(""));
        }

        // normalised residuals y
        for (size_t iPlane = 0; iPlane < geo::gGeometry().nPlanes(); iPlane++) {
            sstm << _histName::_normResidGblFitHistNameY << geo::gGeometry().sensorIDsVec().at(iPlane);
            residGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Normalised residuals. Plane " << geo::gGeometry().sensorIDsVec().at(iPlane) << "Y direction; Normalised residuals; N hits";
            histTitle = sstm.str();
            sstm.str(std::string(""));
            histoInfo = histoMgr->getHistogramInfo(residGblFitHistName);
            NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : resid1dNBinX;
            MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : residMinX;
            MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : residMaxX;
            AIDA::IHistogram1D * residGblFit =
                    marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(residGblFitHistName, NBinX, MinX, MaxX);
            if (residGblFit) {
                residGblFit->setTitle(histTitle);
                _aidaHistoMap1D.insert(std::make_pair(residGblFitHistName, residGblFit));
            } else {
                streamlog_out(ERROR2) << "Problem booking the " << (residGblFitHistName) << std::endl;
                streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            }
            sstm.str(std::string(""));
        }
        
        // plane residuals y
        for (size_t iPlane = 0; iPlane < geo::gGeometry().nPlanes(); iPlane++) {
            sstm << _histName::_residGblFitHistNameY << geo::gGeometry().sensorIDsVec().at(iPlane);
            residGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Residuals. Plane " << geo::gGeometry().sensorIDsVec().at(iPlane) << "Y direction; residuals [um]; N hits";
            histTitle = sstm.str();
            sstm.str(std::string(""));
            histoInfo = histoMgr->getHistogramInfo(residGblFitHistName);
            NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : resid1dNBinX;
            MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : residMinX;
            MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : residMaxX;
            AIDA::IHistogram1D * residGblFit =
                    marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(residGblFitHistName, NBinX, MinX, MaxX);
            if (residGblFit) {
                residGblFit->setTitle(histTitle);
                _aidaHistoMap1D.insert(std::make_pair(residGblFitHistName, residGblFit));
            } else {
                streamlog_out(ERROR2) << "Problem booking the " << (residGblFitHistName) << std::endl;
                streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            }
            sstm.str(std::string(""));
        }

        // 2D histograms
        const int residNBinY = 100;
        const int residNBinX = 100;
        const double hitposMinX = -10.;
        const double hitposMaxX =  10.;
        const double hitposMinY =  -5.;
        const double hitposMaxY =   5.;
        const double residMinY = -200.;
        const double residMaxY =  200.;
        const double residMinZ = -1000.;
        const double residMaxZ =  1000.;
        
        NBinY = residNBinY;
        NBinX = residNBinX;
        MinX = hitposMinX;
        MaxX = hitposMaxX;
        MinY = residMinY;
        MaxY = residMaxY;
        MinZ = residMinZ;
        MaxZ = residMaxZ;
        std::string resid2DGblFitHistName;
        
        // normalised residuals x
        for (size_t iPlane = 0; iPlane < geo::gGeometry().nPlanes(); iPlane++) {
            sstm << _histName::_normResid2DGblFitHistNameXvsX << geo::gGeometry().sensorIDsVec().at(iPlane);
            resid2DGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Normalised residuals. Plane " << geo::gGeometry().sensorIDsVec().at(iPlane) << "; x (mm); Normalised residuals rx";
            histTitle = sstm.str();
            sstm.str(std::string(""));
            histoInfo = histoMgr->getHistogramInfo(resid2DGblFitHistName);
            NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : residNBinX;
            MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : hitposMinX;
            MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : hitposMaxX;
            NBinY = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : residNBinY;
            MinY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : residMinY;
            MaxY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax : residMaxY;
            AIDA::IHistogram2D * residGblFit1 =
                    marlin::AIDAProcessor::histogramFactory(this)->createHistogram2D(resid2DGblFitHistName, NBinX, MinX, MaxX, NBinY, MinY, MaxY);
            if (residGblFit1) {
                residGblFit1->setTitle(histTitle);
                _aidaHistoMap2D.insert(std::make_pair(resid2DGblFitHistName, residGblFit1));
            } else {
                streamlog_out(ERROR2) << "Problem booking the " << (resid2DGblFitHistName) << std::endl;
                streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            }
            sstm.str(std::string(""));

            sstm << _histName::_normResid2DGblFitHistNameXvsY << geo::gGeometry().sensorIDsVec().at(iPlane);
            resid2DGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Normalised residuals. Plane " << geo::gGeometry().sensorIDsVec().at(iPlane) << "; y (mm); Normalised residuals rx";
            histTitle = sstm.str();
            sstm.str(std::string(""));
            histoInfo = histoMgr->getHistogramInfo(resid2DGblFitHistName);
            NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : residNBinX;
            MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : hitposMinX;
            MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : hitposMaxX;
            NBinY = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : residNBinY;
            MinY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : residMinY;
            MaxY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax : residMaxY;
            AIDA::IHistogram2D * residGblFit2 =
                    marlin::AIDAProcessor::histogramFactory(this)->createHistogram2D(resid2DGblFitHistName, NBinX, MinX, MaxX, NBinY, MinY, MaxY);
            if (residGblFit2) {
                residGblFit2->setTitle(histTitle);
                _aidaHistoMap2D.insert(std::make_pair(resid2DGblFitHistName, residGblFit2));
            } else {
                streamlog_out(ERROR2) << "Problem booking the " << (resid2DGblFitHistName) << std::endl;
                streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            }
            sstm.str(std::string(""));
        }
        
        // plane residuals x
        for (size_t iPlane = 0; iPlane < geo::gGeometry().nPlanes(); iPlane++) {
            sstm << _histName::_resid2DGblFitHistNameXvsX << geo::gGeometry().sensorIDsVec().at(iPlane);
            resid2DGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Residuals. Plane " << geo::gGeometry().sensorIDsVec().at(iPlane) << "; x (mm); Residuals rx [um]";
            histTitle = sstm.str();
            sstm.str(std::string(""));
            histoInfo = histoMgr->getHistogramInfo(resid2DGblFitHistName);
            NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : residNBinX;
            MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : hitposMinX;
            MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : hitposMaxX;
            NBinY = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : residNBinY;
            MinY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : residMinY;
            MaxY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax : residMaxY;
            AIDA::IHistogram2D * residGblFit1 =
                    marlin::AIDAProcessor::histogramFactory(this)->createHistogram2D(resid2DGblFitHistName, NBinX, MinX, MaxX, NBinY, MinY, MaxY);
            if (residGblFit1) {
                residGblFit1->setTitle(histTitle);
                _aidaHistoMap2D.insert(std::make_pair(resid2DGblFitHistName, residGblFit1));
            } else {
                streamlog_out(ERROR2) << "Problem booking the " << (resid2DGblFitHistName) << std::endl;
                streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            }
            sstm.str(std::string(""));

            sstm << _histName::_resid2DGblFitHistNameXvsY << geo::gGeometry().sensorIDsVec().at(iPlane);
            resid2DGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Residuals. Plane " << geo::gGeometry().sensorIDsVec().at(iPlane) << "; y (mm); Residuals rx [um]";
            histTitle = sstm.str();
            sstm.str(std::string(""));
            histoInfo = histoMgr->getHistogramInfo(resid2DGblFitHistName);
            NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : residNBinX;
            MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : hitposMinX;
            MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : hitposMaxX;
            NBinY = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : residNBinY;
            MinY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : residMinY;
            MaxY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax : residMaxY;
            AIDA::IHistogram2D * residGblFit2 =
                    marlin::AIDAProcessor::histogramFactory(this)->createHistogram2D(resid2DGblFitHistName, NBinX, MinX, MaxX, NBinY, MinY, MaxY);
            if (residGblFit2) {
                residGblFit2->setTitle(histTitle);
                _aidaHistoMap2D.insert(std::make_pair(resid2DGblFitHistName, residGblFit2));
            } else {
                streamlog_out(ERROR2) << "Problem booking the " << (resid2DGblFitHistName) << std::endl;
                streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            }
            sstm.str(std::string(""));
        }

        // normalised residuals y
        for (size_t iPlane = 0; iPlane < geo::gGeometry().nPlanes(); iPlane++) {
            sstm << _histName::_normResid2DGblFitHistNameYvsX << geo::gGeometry().sensorIDsVec().at(iPlane);
            resid2DGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Normalised residuals. Plane " << geo::gGeometry().sensorIDsVec().at(iPlane) << "; x (mm); Normalised residuals ry";
            histTitle = sstm.str();
            sstm.str(std::string(""));
            histoInfo = histoMgr->getHistogramInfo(resid2DGblFitHistName);
            NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : residNBinX;
            MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : hitposMinX;
            MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : hitposMaxX;
            NBinY = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : residNBinY;
            MinY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : residMinY;
            MaxY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax : residMaxY;
            AIDA::IHistogram2D * residGblFit1 =
                    marlin::AIDAProcessor::histogramFactory(this)->createHistogram2D(resid2DGblFitHistName, NBinX, MinX, MaxX, NBinY, MinY, MaxY);
            if (residGblFit1) {
                residGblFit1->setTitle(histTitle);
                _aidaHistoMap2D.insert(std::make_pair(resid2DGblFitHistName, residGblFit1));
            } else {
                streamlog_out(ERROR2) << "Problem booking the " << (resid2DGblFitHistName) << std::endl;
                streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            }
            sstm.str(std::string(""));

            sstm << _histName::_normResid2DGblFitHistNameYvsY << geo::gGeometry().sensorIDsVec().at(iPlane);
            resid2DGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Normalised residuals. Plane " << geo::gGeometry().sensorIDsVec().at(iPlane) << "; y (mm); Normalised residuals ry";
            histTitle = sstm.str();
            sstm.str(std::string(""));
            histoInfo = histoMgr->getHistogramInfo(resid2DGblFitHistName);
            NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : residNBinX;
            MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : hitposMinX;
            MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : hitposMaxX;
            NBinY = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : residNBinY;
            MinY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : residMinY;
            MaxY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax : residMaxY;
            AIDA::IHistogram2D * residGblFit2 =
                    marlin::AIDAProcessor::histogramFactory(this)->createHistogram2D(resid2DGblFitHistName, NBinX, MinX, MaxX, NBinY, MinY, MaxY);
            if (residGblFit2) {
                residGblFit2->setTitle(histTitle);
                _aidaHistoMap2D.insert(std::make_pair(resid2DGblFitHistName, residGblFit2));
            } else {
                streamlog_out(ERROR2) << "Problem booking the " << (resid2DGblFitHistName) << std::endl;
                streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            }
            sstm.str(std::string(""));
        }
        
        // plane residuals y
        for (size_t iPlane = 0; iPlane < geo::gGeometry().nPlanes(); iPlane++) {
            sstm << _histName::_resid2DGblFitHistNameYvsX << geo::gGeometry().sensorIDsVec().at(iPlane);
            resid2DGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Residuals. Plane " << geo::gGeometry().sensorIDsVec().at(iPlane) << "; x (mm); Residuals ry [um]";
            histTitle = sstm.str();
            sstm.str(std::string(""));
            histoInfo = histoMgr->getHistogramInfo(resid2DGblFitHistName);
            NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : residNBinX;
            MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : hitposMinX;
            MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : hitposMaxX;
            NBinY = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : residNBinY;
            MinY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : residMinY;
            MaxY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax : residMaxY;
            AIDA::IHistogram2D * residGblFit1 =
                    marlin::AIDAProcessor::histogramFactory(this)->createHistogram2D(resid2DGblFitHistName, NBinX, MinX, MaxX, NBinY, MinY, MaxY);
            if (residGblFit1) {
                residGblFit1->setTitle(histTitle);
                _aidaHistoMap2D.insert(std::make_pair(resid2DGblFitHistName, residGblFit1));
            } else {
                streamlog_out(ERROR2) << "Problem booking the " << (resid2DGblFitHistName) << std::endl;
                streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            }
            sstm.str(std::string(""));

            sstm << _histName::_resid2DGblFitHistNameYvsY << geo::gGeometry().sensorIDsVec().at(iPlane);
            resid2DGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Residuals. Plane " << geo::gGeometry().sensorIDsVec().at(iPlane) << "; y (mm); Residuals ry [um]";
            histTitle = sstm.str();
            sstm.str(std::string(""));
            histoInfo = histoMgr->getHistogramInfo(resid2DGblFitHistName);
            NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : residNBinX;
            MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : hitposMinX;
            MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : hitposMaxX;
            NBinY = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : residNBinY;
            MinY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : residMinY;
            MaxY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax : residMaxY;
            AIDA::IHistogram2D * residGblFit2 =
                    marlin::AIDAProcessor::histogramFactory(this)->createHistogram2D(resid2DGblFitHistName, NBinX, MinX, MaxX, NBinY, MinY, MaxY);
            if (residGblFit2) {
                residGblFit2->setTitle(histTitle);
                _aidaHistoMap2D.insert(std::make_pair(resid2DGblFitHistName, residGblFit2));
            } else {
                streamlog_out(ERROR2) << "Problem booking the " << (resid2DGblFitHistName) << std::endl;
                streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            }
            sstm.str(std::string(""));
        }

        // 2D profile residuals x
        for (size_t iPlane = 0; iPlane < geo::gGeometry().nPlanes(); iPlane++) {
            sstm << _histName::_resid2DGblFitHistNameX << geo::gGeometry().sensorIDsVec().at(iPlane);
            residGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Residuals. Plane " << geo::gGeometry().sensorIDsVec().at(iPlane) << ";X direction; Y direction";
            histTitle = sstm.str();
            sstm.str(std::string(""));
            histoInfo = histoMgr->getHistogramInfo(residGblFitHistName);
            NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : residNBinY;
            MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : hitposMinX;
            MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : hitposMaxX;
            NBinY = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : residNBinY;
            MinY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : hitposMinY;
            MaxY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax : hitposMaxY;
            MinZ =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_zMin : residMinZ;
            MaxZ =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_zMax : residMaxZ;
            AIDA::IProfile2D *  residGblFit =
                    marlin::AIDAProcessor::histogramFactory(this)->createProfile2D(residGblFitHistName,  NBinX, MinX, MaxX, NBinY, MinY, MaxY, MinZ,MaxZ);
            if (residGblFit) {
                residGblFit->setTitle(histTitle);
                _aidaProfileMap2D.insert(std::make_pair(residGblFitHistName, residGblFit));
            } else {
                streamlog_out(ERROR2) << "Problem booking the " << (residGblFitHistName) << std::endl;
                streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            }
            sstm.str(std::string(""));
        }


        // 2D profile residuals y
        for (size_t iPlane = 0; iPlane < geo::gGeometry().nPlanes(); iPlane++) {
            sstm << _histName::_resid2DGblFitHistNameY << geo::gGeometry().sensorIDsVec().at(iPlane);
            residGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Residuals. Plane " << geo::gGeometry().sensorIDsVec().at(iPlane) << ";X direction; Y direction";
            histTitle = sstm.str();
            sstm.str(std::string(""));
            histoInfo = histoMgr->getHistogramInfo(residGblFitHistName);
            NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : residNBinY;
            MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : hitposMinX;
            MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : hitposMaxX;
            NBinY = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : residNBinY;
            MinY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : hitposMinY;
            MaxY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax : hitposMaxY;
            MinZ =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_zMin : residMinZ;
            MaxZ =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_zMax : residMaxZ;
            AIDA::IProfile2D *  residGblFit =
                    marlin::AIDAProcessor::histogramFactory(this)->createProfile2D(residGblFitHistName,  NBinX, MinX, MaxX, NBinY, MinY, MaxY, MinZ,MaxZ);
            if (residGblFit) {
                residGblFit->setTitle(histTitle);
                _aidaProfileMap2D.insert(std::make_pair(residGblFitHistName, residGblFit));
            } else {
                streamlog_out(ERROR2) << "Problem booking the " << (residGblFitHistName) << std::endl;
                streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            }
            sstm.str(std::string(""));
        }


        // Kink angles after fit
        const int    kinkNBinX = 100;
        const double kinkMinX = -0.001;
        const double kinkMaxX = 0.001;
        NBinX = kinkNBinX;
        MinX = kinkMinX;
        MaxX = kinkMaxX;
        std::string kinkGblFitHistName;
        for (size_t iPlane = 0; iPlane < geo::gGeometry().nPlanes(); iPlane++) {
            sstm << _histName::_kinkGblFitHistNameX << geo::gGeometry().sensorIDsVec().at(iPlane);
            kinkGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Kink angles. Plane " << geo::gGeometry().sensorIDsVec().at(iPlane) << "X direction; kink (rad); N hits";
            histTitle = sstm.str();
            sstm.str(std::string(""));
            histoInfo = histoMgr->getHistogramInfo(kinkGblFitHistName);
            NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : kinkNBinX;
            MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : kinkMinX;
            MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : kinkMaxX;
            AIDA::IHistogram1D * residGblFit =
                    marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(kinkGblFitHistName, NBinX, MinX, MaxX);
            if (residGblFit) {
                residGblFit->setTitle(histTitle);
                _aidaHistoMap1D.insert(std::make_pair(kinkGblFitHistName, residGblFit));
            } else {
                streamlog_out(ERROR2) << "Problem booking the " << (kinkGblFitHistName) << std::endl;
                streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            }
            sstm.str(std::string(""));
        }

        for (size_t iPlane = 0; iPlane < geo::gGeometry().nPlanes(); iPlane++) {
            sstm << _histName::_kinkGblFitHistNameY << geo::gGeometry().sensorIDsVec().at(iPlane);
            kinkGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Kink angles. Plane " << geo::gGeometry().sensorIDsVec().at(iPlane) << "Y direction; kink (rad); N hits";
            histTitle = sstm.str();
            sstm.str(std::string(""));
            histoInfo = histoMgr->getHistogramInfo(kinkGblFitHistName);
            NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : kinkNBinX;
            MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : kinkMinX;
            MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : kinkMaxX;
            AIDA::IHistogram1D * residGblFit =
                    marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(kinkGblFitHistName, NBinX, MinX, MaxX);
            if (residGblFit) {
                residGblFit->setTitle(histTitle);
                _aidaHistoMap1D.insert(std::make_pair(kinkGblFitHistName, residGblFit));
            } else {
                streamlog_out(ERROR2) << "Problem booking the " << (kinkGblFitHistName) << std::endl;
                streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            }
            sstm.str(std::string(""));
        }
        _flag_nohistos = false; 
   } catch (lcio::Exception& e) {
        streamlog_out(WARNING2) << "Can't allocate histograms. Continue without histogramming; exception: " << e.what() << endl;
        _flag_nohistos = true; 
    }
#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
}

#endif // USE_GBL

