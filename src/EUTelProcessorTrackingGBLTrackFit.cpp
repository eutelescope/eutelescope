#ifdef USE_GBL

// C++
#include <map>
#include <memory>
#include <string>
#include <vector>

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
#include "EUTelExhaustiveTrackFinder.h"
#include "EUTelLCObjectTrackCandidate.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

// Cluster types
class EUTelSparseCluster2Impl;
class EUTelSparseClusterImpl;
class EUTelBrickedClusterImpl;
class EUTelDFFClusterImpl;
class EUTelFFClusterImpl;

using namespace lcio;
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
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_chi2GblFitHistName = "chi2GblFit";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_probGblFitHistName = "probGblFit";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_residGblFitHistName = "ResidualsGblFit";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_residGblFitHistNameX = "ResidualsGblFit_x";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_residGblFitHistNameY = "ResidualsGblFit_y";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_resid2DGblFitHistNameXvsX = "Residuals2DGblFit_xVSx";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_resid2DGblFitHistNameYvsX = "Residuals2DGblFit_yVSx";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_resid2DGblFitHistNameXvsY = "Residuals2DGblFit_xVSy";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_resid2DGblFitHistNameYvsY = "Residuals2DGblFit_yVSy";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_kinkGblFitHistNameX = "KinksGblFit_x";
std::string EUTelProcessorTrackingGBLTrackFit::_histName::_kinkGblFitHistNameY = "KinksGblFit_y";
#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

/** Default constructor */
EUTelProcessorTrackingGBLTrackFit::EUTelProcessorTrackingGBLTrackFit() :
Processor("EUTelProcessorTrackingGBLTrackFit"),
_eBeam(-1.),
_alignmentMode(0),
_xShiftsVec(),
_yShiftsVec(),
_zShiftsVec(),
_xRotationsVec(),
_yRotationsVec(),
_zRotationsVec(),
_milleBinaryFilename("mille.bin"),
_milleSteeringFilename("pede-steer.txt"),
_alignmentPlaneIds(),
_runPede(false),
_maxChi2Cut(1000.),
_tgeoFileName("TELESCOPE.root"),
_histoInfoFileName("histoinfo.xml"),
_trackCandidateHitsInputCollectionName("TrackCandidateHitCollection"),
_tracksOutputCollectionName("TrackCollection"),
_trackFitter(0),
_milleGBL(0),
_alignmentConstants(),
_nProcessedRuns(0),
_nProcessedEvents(0),
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
_aidaHistoMap1D(),
_aidaHistoMap2D()
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

    // Optional processor parameters that define finder settings

    // MILLIPEDE specific parameters
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

    registerOptionalParameter("MilleParametersXShifts", "Plane ids and parameter ids for X shifts", _xShiftsVec, IntVec());

    registerOptionalParameter("MilleParametersYShifts", "Plane ids and parameter ids for Y shifts", _yShiftsVec, IntVec());
    
    registerOptionalParameter("MilleParametersZShifts", "Plane ids and parameter ids for Z shifts", _zShiftsVec, IntVec());

    registerOptionalParameter("MilleParametersXRotations", "Plane ids and parameter ids for X rotations", _xRotationsVec, IntVec());
    
    registerOptionalParameter("MilleParametersYRotations", "Plane ids and parameter ids for Y rotations", _yRotationsVec, IntVec());
    
    registerOptionalParameter("MilleParametersZRotations", "Plane ids and parameter ids for Z rotations", _zRotationsVec, IntVec());

    registerOptionalParameter("MilleBinaryFilename", "Name of the Millepede binary file", _milleBinaryFilename, std::string("mille.bin"));

    registerOptionalParameter("MilleSteeringFilename", "Name of the Millepede steering file to be created", _milleSteeringFilename, std::string("pede-steer.txt"));

    registerOptionalParameter("MilleMaxChi2Cut", "Maximum chi2 of a track candidate that goes into millepede", _maxChi2Cut, double(1000.));

    registerOptionalParameter("AlignmentPlanes", "Ids of planes to be used in alignment", _alignmentPlaneIds, IntVec());
    
    registerOptionalParameter("RunPede","Execute the pede at the end of processing using the generated steering file.",_runPede, static_cast <bool> (false));
    
    // Geometry definition
    
    registerOptionalParameter("GeometryFilename", "Name of the TGeo geometry definition file", _tgeoFileName, std::string("TELESCOPE.root"));
    
    // Histogram information

    registerOptionalParameter("HistogramInfoFilename", "Name of histogram info xml file", _histoInfoFileName, std::string("histoinfo.xml"));
}

void EUTelProcessorTrackingGBLTrackFit::init() {

    streamlog_out(DEBUG2) << "EUTelProcessorTrackingGBLTrackFit::init( )" << std::endl;

    // usually a good idea to
    printParameters();

    // Reset counters
    _nProcessedRuns = 0;
    _nProcessedEvents = 0;


    // Getting access to geometry description
    geo::gGeometry().initializeTGeoDescription(_tgeoFileName);

    // Instantiate millipede output. 
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

        EUTelGBLFitter* Fitter = new EUTelGBLFitter("myGBLFitter");
        Fitter->SetAlignmentMode(alignmentMode);
        Fitter->SetXShiftsVec(_xShiftsVec);
        Fitter->SetYShiftsVec(_yShiftsVec);
        Fitter->SetZShiftsVec(_zShiftsVec);
        Fitter->SetXRotationsVec(_xRotationsVec);
        Fitter->SetYRotationsVec(_yRotationsVec);
        Fitter->SetZRotationsVec(_zRotationsVec);
        Fitter->SetMilleBinary(_milleGBL);
        Fitter->SetBeamEnergy(_eBeam);
        Fitter->SetChi2Cut(_maxChi2Cut);
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

    
    streamlog_out(DEBUG1) << "Flush alignment constants (X shifts)" << std::endl;
    for (std::map<int, double >::iterator iDet = _alignmentConstants._xResiduals.begin();
            iDet != _alignmentConstants._xResiduals.end(); ++iDet) {
        iDet->second=0.;
    }
    
    for (std::map<int, int >::iterator iDet = _alignmentConstants._nxResiduals.begin();
            iDet != _alignmentConstants._nxResiduals.end(); ++iDet) {
        iDet->second=1;
    }
    
    streamlog_out(DEBUG1) << "Flush alignment constants (Y shifts)" << std::endl;
    for (std::map<int, double >::iterator iDet = _alignmentConstants._yResiduals.begin();
            iDet != _alignmentConstants._yResiduals.end(); ++iDet) {
        iDet->second=0.;
    }
    
    for (std::map<int, int >::iterator iDet = _alignmentConstants._nyResiduals.begin();
            iDet != _alignmentConstants._nyResiduals.end(); ++iDet) {
        iDet->second=1;
    }
    
    _nProcessedRuns++;
}

void EUTelProcessorTrackingGBLTrackFit::processEvent(LCEvent * evt) {

    if (isFirstEvent()) {
        ;
    }

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
    } catch (DataNotAvailableException e) {
        streamlog_out(WARNING) << _trackCandidateHitsInputCollectionName << " collection not available" << std::endl;
        throw marlin::SkipEventException(this);
    }

    // this will only be entered if the collection is available
    if (col != NULL) {
        if (_nProcessedEvents % 1000 == 1) streamlog_out(DEBUG2) << "EUTelProcessorTrackingGBLTrackFit" << endl;

        vector < EVENT::TrackerHitVec > trackCandidates;
        for (int iCol = 0; iCol < col->getNumberOfElements(); iCol++) {
            TrackImpl* track = static_cast<TrackImpl*> (col->getElementAt(iCol));

            if (!col) {
                streamlog_out(WARNING) << "EUTelLCObjectTrackCandidate collection not found found for event " << _nProcessedEvents <<
                        " in run " << _nProcessedRuns << endl;
                throw SkipEventException(this);
            }

            streamlog_out(DEBUG1) << "Track " << iCol << " nhits " << track->getTrackerHits().size() << endl;
            trackCandidates.push_back(track->getTrackerHits());
        } //for ( int iCol = 0; iCol < col->getNumberOfElements() ; iCol++ )

        //        return;

        // Perform fit for all found track candidates
        // ------------------------------------------
        unsigned int numData;
        TVectorD residual(200);
        TVectorD measErr(200);
        TVectorD residualErr(200);
        TVectorD downWeight(200);
        int nTracks = trackCandidates.size();
        streamlog_out(DEBUG1) << "N tracks found " << nTracks << endl;
        if (nTracks != 0) {
            _trackFitter->SetTrackCandidates(trackCandidates);
            _trackFitter->FitTracks();
            //
            double chi2Trk = 0.;
            int ndfTrk = 0;

            IMPL::LCCollectionVec* fittrackvec;
            fittrackvec = static_cast<EUTelGBLFitter*> (_trackFitter)->GetFitTrackVec();
            IMPL::LCCollectionVec::const_iterator itFitTrack;

            int iCounter = 0;
            for (itFitTrack = fittrackvec->begin(); itFitTrack != fittrackvec->end(); ++itFitTrack) {
                chi2Trk = static_cast<TrackImpl*> (*itFitTrack)->getChi2();
                ndfTrk = static_cast<TrackImpl*> (*itFitTrack)->getNdf();

                static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ _histName::_chi2GblFitHistName ]) -> fill(chi2Trk);
                static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ _histName::_probGblFitHistName ]) -> fill(TMath::Prob(chi2Trk, ndfTrk));


                std::map< int, gbl::GblTrajectory* > gblTracks = static_cast<EUTelGBLFitter*> (_trackFitter)->GetGblTrackCandidates();

                std::stringstream sstr;
                gbl::GblTrajectory* gblTraj = gblTracks[ iCounter ];
                //                gblTraj->printTrajectory( );
                //                gblTraj->printPoints( );
                //                gblTraj->printData( );
                std::vector< gbl::GblPoint > gblPointVec = static_cast<EUTelGBLFitter*> (_trackFitter)->GetGblTracksPoints()[iCounter];
                std::vector< gbl::GblPoint >::const_iterator itGblPoint = gblPointVec.begin();
                int iPlane = 0; // wrong in case of missing planes
                for (; itGblPoint != gblPointVec.end(); ++itGblPoint) {
                    if (iPlane > 5) continue;
                    //if ( itGblPoint->getLabel() < 1000 )
                    if (itGblPoint->getLabel() % 3 == 1) {
                        streamlog_out(DEBUG0) << std::setw(15) << itGblPoint->getLabel() << std::endl;
                        // spatial residuals
                        gblTraj->getMeasResults(itGblPoint->getLabel(), numData, residual, measErr, residualErr, downWeight);
                        sstr << _histName::_residGblFitHistNameX << iPlane;
                        streamlog_out(DEBUG0) << std::setw(15) << std::setprecision(5) << residual[0] << std::setw(15) << std::setprecision(5) << residualErr[0] << std::endl;
                        static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ sstr.str() ]) -> fill(residual[0] / residualErr[0], downWeight[0]);
                        _alignmentConstants._xResiduals[iPlane]+=(residual[0]);
                        _alignmentConstants._nxResiduals[iPlane]++;
                        sstr.str(std::string());
                        sstr << _histName::_residGblFitHistNameY << iPlane;
                        streamlog_out(DEBUG0) << std::setw(15) << std::setprecision(5) << residual[1] << std::setw(15) << std::setprecision(5) << residualErr[1] << std::endl;
                        static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ sstr.str() ]) -> fill(residual[1] / residualErr[1], downWeight[1]);
                        _alignmentConstants._yResiduals[iPlane]+=(residual[1]);
                        _alignmentConstants._nyResiduals[iPlane]++;
                        sstr.str(std::string());
                        // kinks
                        gblTraj->getScatResults(itGblPoint->getLabel(), numData, residual, measErr, residualErr, downWeight);
                        sstr << _histName::_kinkGblFitHistNameX << iPlane;
                        streamlog_out(DEBUG0) << std::setw(15) << std::setprecision(5) << residual[0] << std::setw(15) << std::setprecision(5) << residualErr[0] << std::endl;
                        static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ sstr.str() ]) -> fill(residual[0], downWeight[0]);
                        sstr.str(std::string());
                        sstr << _histName::_kinkGblFitHistNameY << iPlane;
                        streamlog_out(DEBUG0) << std::setw(15) << std::setprecision(5) << residual[1] << std::setw(15) << std::setprecision(5) << residualErr[1] << std::endl;
                        static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ sstr.str() ]) -> fill(residual[1], downWeight[1]);
                        sstr.str(std::string());

                        // 2D histograms
                        EVENT::TrackerHitVec hit_track1 = trackCandidates.front();
                        const double* hitpos = hit_track1[iPlane]->getPosition(); // wrong in case of empty planes
                        sstr << _histName::_resid2DGblFitHistNameXvsX << iPlane;
                        static_cast<AIDA::IHistogram2D*> (_aidaHistoMap2D[ sstr.str() ]) -> fill(hitpos[0], residual[0] / residualErr[0], downWeight[0]);
                        sstr.str(std::string());
                        sstr << _histName::_resid2DGblFitHistNameXvsY << iPlane;
                        static_cast<AIDA::IHistogram2D*> (_aidaHistoMap2D[ sstr.str() ]) -> fill(hitpos[1], residual[0] / residualErr[0], downWeight[0]);
                        sstr.str(std::string());
                        sstr << _histName::_resid2DGblFitHistNameYvsX << iPlane;
                        static_cast<AIDA::IHistogram2D*> (_aidaHistoMap2D[ sstr.str() ]) -> fill(hitpos[0], residual[1] / residualErr[1], downWeight[1]);
                        sstr.str(std::string());
                        sstr << _histName::_resid2DGblFitHistNameYvsY << iPlane;
                        static_cast<AIDA::IHistogram2D*> (_aidaHistoMap2D[ sstr.str() ]) -> fill(hitpos[1], residual[1] / residualErr[1], downWeight[1]);
                        sstr.str(std::string());

                        if (itGblPoint->getLabel() < 1000)++iPlane;
                    }
                }

                IMPL::LCCollectionVec::const_iterator itFitTrack;
                iCounter++;
                
                delete gblTraj;
            }
        } //if( _ntracks != 0 && _ntracks == 1)

    } //if( col != NULL )

    _nProcessedEvents++;

    if (isFirstEvent()) _isFirstEvent = false;
}

void EUTelProcessorTrackingGBLTrackFit::check(LCEvent * evt) {
    // nothing to check here
}

void EUTelProcessorTrackingGBLTrackFit::end() {

    delete _milleGBL;

    writeMilleSteeringFile();

    streamlog_out(DEBUG) << "EUTelProcessorTrackingGBLTrackFit::end()  " << name()
            << " processed " << _nProcessedEvents << " events in " << _nProcessedRuns << " runs "
            << std::endl;
    
    if ( _runPede ) {

        // check if alignment was requested
        if ( _alignmentMode == Utility::noAlignment) {
                streamlog_out(WARNING1) << "RunPede was required, but alignment mode is noAlignment. Stop." << endl;
                return;
        }

        std::string command = "pede " + _milleSteeringFilename;
        streamlog_out ( MESSAGE5 ) << "Starting pede...: " << command.c_str() << endl;
    
        // run pede and create a streambuf that reads its stdout and stderr
        redi::ipstream pede( command.c_str(), redi::pstreams::pstdout|redi::pstreams::pstderr ); 
      
        if (!pede.is_open()) {
            streamlog_out( ERROR5 ) << "Pede cannot be executed: command not found in the path" << endl; 
        } else {
            
            bool encounteredError = false;
            // output multiplexing: parse pede output in both stdout and stderr and echo messages accordingly
            char buf[1024];
            std::streamsize n;
            std::stringstream pedeoutput; // store stdout to parse later
            std::stringstream pedeerrors;
            bool finished[2] = { false, false };
            while (!finished[0] || !finished[1])
              {
                if (!finished[0])
                  {
                    while ((n = pede.err().readsome(buf, sizeof(buf))) > 0){
                      streamlog_out( ERROR5 ).write(buf, n).flush();
                      string error (buf, n);
                      pedeerrors << error;
                      encounteredError = true;
                    }
                    if (pede.eof())
                      {
                        finished[0] = true;
                        if (!finished[1])
                          pede.clear();
                      }
                  }

                if (!finished[1])
                  {
                    while ((n = pede.out().readsome(buf, sizeof(buf))) > 0){
                      streamlog_out( MESSAGE4 ).write(buf, n).flush();
                      string output (buf, n);
                      pedeoutput << output;
                    }
                    if (pede.eof())
                      {
                        finished[1] = true;
                        if (!finished[0])
                          pede.clear();
                      }
                  }
              }
            // wait for the pede execution to finish
            pede.close();
        }
    }
}

void EUTelProcessorTrackingGBLTrackFit::writeMilleSteeringFile() {

    streamlog_out(DEBUG2) << "writeMilleSteeringFile" << endl;

    // Prepare millipede steering files only if alignment was requested
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
    std::map<int, int> XShiftsMap = fitter->GetParamterIdXShiftsMap();
    std::map<int, int> YShiftsMap = fitter->GetParamterIdYShiftsMap();
    std::map<int, int> ZShiftsMap = fitter->GetParamterIdZShiftsMap();
    std::map<int, int> XRotationsMap = fitter->GetParamterIdXRotationsMap();
    std::map<int, int> YRotationsMap = fitter->GetParamterIdYRotationsMap();
    std::map<int, int> ZRotationsMap = fitter->GetParamterIdZRotationsMap();
    
    // loop over all planes
    // @TODO assumes that planes have ids 0..._nplanes !generaly wrong
    for (unsigned int help = 0; help < geo::gGeometry()._nPlanes; help++) {

        int sensorId = geo::gGeometry()._sensorIDVecZOrder[help];
        bool isPlaneExcluded = std::find(_alignmentPlaneIds.begin(), _alignmentPlaneIds.end(), sensorId) == _alignmentPlaneIds.end();
              
        // if plane not excluded
        if ( !isPlaneExcluded ) {

            if( fitter->GetAlignmentMode()==Utility::XYZShiftXYRot ) {
                steerFile << left << setw(25) << XShiftsMap[sensorId] << setw(25) << _alignmentConstants._xResiduals[sensorId]/_alignmentConstants._nxResiduals[sensorId] << setw(25) << "0.01" 
                           << setw(25) << " ! X shift " << setw(25) << sensorId << endl;
                steerFile << left << setw(25) << YShiftsMap[sensorId] << setw(25) << _alignmentConstants._yResiduals[sensorId]/_alignmentConstants._nyResiduals[sensorId] << setw(25) << "0.01" 
                           << setw(25) << " ! Y shift " << setw(25) << sensorId << endl;
                steerFile << left << setw(25) << ZShiftsMap[sensorId] << setw(25) << "0.0" << setw(25) << "0.01" 
                           << setw(25) << " ! Z shift " << setw(25) << sensorId << endl;
                steerFile << left << setw(25) << ZRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << "0.01"
                          << setw(25) << " ! XY rotation " << sensorId << endl;
            } else if( fitter->GetAlignmentMode()==Utility::XYShiftYZRotXYRot ) {
                steerFile << left << setw(25) << XShiftsMap[sensorId] << setw(25) << _alignmentConstants._xResiduals[sensorId]/_alignmentConstants._nxResiduals[sensorId] << setw(25) << "0.01" 
                          << setw(25) << " ! X shift " << sensorId << endl;
                steerFile << left << setw(25) << YShiftsMap[sensorId] << setw(25)  << _alignmentConstants._yResiduals[sensorId]/_alignmentConstants._nyResiduals[sensorId] << setw(25) << "0.01" 
                          << setw(25) << " ! Y shift " << sensorId << endl;
                steerFile << left << setw(25) << XRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << "0.01" 
                          << setw(25) << " ! YZ rotation " << sensorId << endl;
                steerFile << left << setw(25) << ZRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << "0.01"
                          << setw(25) << " ! XY rotation " << sensorId << endl;
            } else if( fitter->GetAlignmentMode()==Utility::XYShiftXZRotXYRot ) {
                steerFile << left << setw(25) << XShiftsMap[sensorId] << setw(25) << _alignmentConstants._xResiduals[sensorId]/_alignmentConstants._nxResiduals[sensorId] << setw(25) << "0.01" 
                          << setw(25) << " ! X shift " << sensorId << endl;
                steerFile << left << setw(25) << YShiftsMap[sensorId] << setw(25) << _alignmentConstants._yResiduals[sensorId]/_alignmentConstants._nyResiduals[sensorId] << setw(25) << "0.01" 
                          << setw(25) << " ! Y shift " << sensorId << endl;
                steerFile << left << setw(25) << YRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << "0.01" 
                          << setw(25) << " ! XZ rotation " << sensorId << endl;
                steerFile << left << setw(25) << ZRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << "0.01"
                         << setw(25)  << " ! XY rotation " << sensorId << endl;
            } else if( fitter->GetAlignmentMode()==Utility::XYShiftXZRotYZRotXYRot ) {
                steerFile << left << setw(25) << XShiftsMap[sensorId] << setw(25) << _alignmentConstants._xResiduals[sensorId]/_alignmentConstants._nxResiduals[sensorId] << setw(25) << "0.01" 
                          << setw(25) << " ! X shift " << sensorId << endl;
                steerFile << left << setw(25) << YShiftsMap[sensorId] << setw(25) << _alignmentConstants._yResiduals[sensorId]/_alignmentConstants._nyResiduals[sensorId] << setw(25) << "0.01" 
                          << setw(25) << " ! Y shift " << sensorId << endl;
                steerFile << left << setw(25) << YRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << "0.01" 
                          << setw(25) << " ! XZ rotation " << sensorId << endl;
                steerFile << left << setw(25) << XRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << "0.01" 
                          << setw(25) << " ! YZ rotation " << sensorId << endl;
                steerFile << left << setw(25) << ZRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << "0.01"
                         << setw(25)  << " ! XY rotation " << sensorId << endl;
            } else if( fitter->GetAlignmentMode()==Utility::XYZShiftXZRotYZRotXYRot ) {
                steerFile << left << setw(25) << XShiftsMap[sensorId] << setw(25) << _alignmentConstants._xResiduals[sensorId]/_alignmentConstants._nxResiduals[sensorId] << setw(25) << "0.01" 
                          << setw(25) << " ! X shift " << sensorId << endl;
                steerFile << left << setw(25) << YShiftsMap[sensorId] << setw(25) << _alignmentConstants._yResiduals[sensorId]/_alignmentConstants._nyResiduals[sensorId] << setw(25) << "0.01" 
                          << setw(25) << " ! Y shift " << sensorId << endl;
                steerFile << left << setw(25) << ZShiftsMap[sensorId] << setw(25) << "0.0" << setw(25) << "0.01" 
                          << setw(25) << " ! Z shift " << sensorId << endl;
                steerFile << left << setw(25) << YRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << "0.01" 
                          << setw(25) << " ! XZ rotation " << sensorId << endl;
                steerFile << left << setw(25) << XRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << "0.01" 
                          << setw(25) << " ! YZ rotation " << sensorId << endl;
                steerFile << left << setw(25) << ZRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << "0.01"
                         << setw(25)  << " ! XY rotation " << sensorId << endl;
            } else if ( fitter->GetAlignmentMode()==Utility::XYShiftXYRot ) {
                steerFile << left << setw(25) << XShiftsMap[sensorId] << setw(25) << _alignmentConstants._xResiduals[sensorId]/_alignmentConstants._nxResiduals[sensorId] << setw(25) << "0.01" 
                          << setw(25) << " ! X shift " << sensorId << endl;
                steerFile << left << setw(25) << YShiftsMap[sensorId] << setw(25) << _alignmentConstants._yResiduals[sensorId]/_alignmentConstants._nyResiduals[sensorId] << setw(25) << "0.01" 
                          << setw(25) << " ! Y shift " << sensorId << endl;
                steerFile << left << setw(25) << ZRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << "0.01" 
                          << setw(25) << " ! XY rotation " << sensorId << endl;
            } else if ( fitter->GetAlignmentMode()==Utility::XYShift ) {
                steerFile << left << setw(25) << XShiftsMap[sensorId] << setw(25) << _alignmentConstants._xResiduals[sensorId]/_alignmentConstants._nxResiduals[sensorId] << setw(25) << "0.01"
                          << setw(25) << " ! X shift " << sensorId << endl;
                steerFile << left << setw(25) << YShiftsMap[sensorId] << setw(25) << _alignmentConstants._yResiduals[sensorId]/_alignmentConstants._nyResiduals[sensorId] << setw(25) << "0.01"
                          << setw(25) << " ! Y shift " << sensorId << endl;
                steerFile << left << setw(25) << ZRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << "-1.0"
                          << setw(25) << " ! XY rotation fixed" << sensorId << endl;
            }  

            counter++;

        } // end if plane not excluded

    } // end loop over all planes

    steerFile << endl;
    steerFile << "method inversion 5 0.001" << endl;
    steerFile << "chiscut 50.0 10." << endl;
    steerFile << "! outlierdownweighting 4" << endl;
    steerFile << "! histprint" << endl;
    steerFile << endl;
    steerFile << "end" << endl;

    steerFile.close();

    if( _alignmentMode != Utility::noAlignment ) streamlog_out(MESSAGE5) << "File " << _milleSteeringFilename << " written." << endl;

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
        
        histoInfo = histoMgr->getHistogramInfo(_histName::_chi2GblFitHistName);
        int chi2NBin =          ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 1000;    
        double chi2Min =        ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : 0.;
        double chi2Max =        ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 1000.;

        histoInfo = histoMgr->getHistogramInfo(_histName::_probGblFitHistName);
        int probNBin =          ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 1000;
        double probMin =        ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : 0.;
        double probMax =        ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 1.;

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

        AIDA::IHistogram1D * probGblFit =
                marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_probGblFitHistName, probNBin, probMin, probMax);
        if (probGblFit) {
            probGblFit->setTitle("Probability of track fit; Prob;N Tracks");
            _aidaHistoMap1D.insert(std::make_pair(_histName::_probGblFitHistName, probGblFit));
        } else {
            streamlog_out(ERROR2) << "Problem booking the " << (_histName::_probGblFitHistName) << std::endl;
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
        // Residuals after fit
        std::stringstream sstm;
        std::string residGblFitHistName;
        std::string histTitle;
        for (size_t iPlane = 0; iPlane < geo::gGeometry()._nPlanes; iPlane++) {
            sstm << _histName::_residGblFitHistNameX << geo::gGeometry()._sensorIDVec.at(iPlane);
            residGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Normalised residuals. Plane " << geo::gGeometry()._sensorIDVec.at(iPlane) << "X direction; Normalised residuals; N hits";
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

        for (size_t iPlane = 0; iPlane < geo::gGeometry()._nPlanes; iPlane++) {
            sstm << _histName::_residGblFitHistNameY << geo::gGeometry()._sensorIDVec.at(iPlane);
            residGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Normalised residuals. Plane " << geo::gGeometry()._sensorIDVec.at(iPlane) << "Y direction; Normalised residuals; N hits";
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
        const double hitposMinX = -15.;
        const double hitposMaxX = 15.;
        const double residMinY = -20.;
        const double residMaxY = 20.;
        
        NBinY = residNBinY;
        NBinX = residNBinX;
        MinX = hitposMinX;
        MaxX = hitposMaxX;
        MinY = residMinY;
        MaxY = residMaxY;
        std::string resid2DGblFitHistName;
        for (size_t iPlane = 0; iPlane < geo::gGeometry()._nPlanes; iPlane++) {
            sstm << _histName::_resid2DGblFitHistNameXvsX << geo::gGeometry()._sensorIDVec.at(iPlane);
            resid2DGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Normalised residuals. Plane " << geo::gGeometry()._sensorIDVec.at(iPlane) << "; x (mm); Normalised residuals rx";
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

            sstm << _histName::_resid2DGblFitHistNameXvsY << geo::gGeometry()._sensorIDVec.at(iPlane);
            resid2DGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Normalised residuals. Plane " << geo::gGeometry()._sensorIDVec.at(iPlane) << "; y (mm); Normalised residuals rx";
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

        for (size_t iPlane = 0; iPlane < geo::gGeometry()._nPlanes; iPlane++) {
            sstm << _histName::_resid2DGblFitHistNameYvsX << geo::gGeometry()._sensorIDVec.at(iPlane);
            resid2DGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Normalised residuals. Plane " << geo::gGeometry()._sensorIDVec.at(iPlane) << "; x (mm); Normalised residuals ry";
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

            sstm << _histName::_resid2DGblFitHistNameYvsY << geo::gGeometry()._sensorIDVec.at(iPlane);
            resid2DGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Normalised residuals. Plane " << geo::gGeometry()._sensorIDVec.at(iPlane) << "; y (mm); Normalised residuals ry";
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

        // Kink angles after fit
        const int    kinkNBinX = 100;
        const double kinkMinX = -0.001;
        const double kinkMaxX = 0.001;
        NBinX = kinkNBinX;
        MinX = kinkMinX;
        MaxX = kinkMaxX;
        std::string kinkGblFitHistName;
        for (int iPlane = 0; iPlane < geo::gGeometry()._nPlanes; iPlane++) {
            sstm << _histName::_kinkGblFitHistNameX << geo::gGeometry()._sensorIDVec.at(iPlane);
            kinkGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Kink angles. Plane " << geo::gGeometry()._sensorIDVec.at(iPlane) << "X direction; kink (rad); N hits";
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

        for (size_t iPlane = 0; iPlane < geo::gGeometry()._nPlanes; iPlane++) {
            sstm << _histName::_kinkGblFitHistNameY << geo::gGeometry()._sensorIDVec.at(iPlane);
            kinkGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Kink angles. Plane " << geo::gGeometry()._sensorIDVec.at(iPlane) << "Y direction; kink (rad); N hits";
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
    } catch (lcio::Exception& e) {
        streamlog_out(WARNING2) << "Can't allocate histgrams. Continue without histogramming" << endl;
    }
#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
}

#endif // USE_GBL

