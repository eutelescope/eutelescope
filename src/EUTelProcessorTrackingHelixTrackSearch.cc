#include "EUTelProcessorTrackingHelixTrackSearch.h"

// C++
#include <map>
#include <memory>
#include <string>
#include <vector>

// LCIO
#include <EVENT/LCCollection.h>
#include <EVENT/TrackerHit.h>
#include <IMPL/TrackImpl.h>

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

// EUTELESCOPE
#include "EUTelHistogramManager.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelUtility.h"
#include "EUTelMagneticFieldFinder.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
// Cluster types
#include "EUTelSparseCluster2Impl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelFFClusterImpl.h"


using namespace lcio;
using namespace marlin;
using namespace eutelescope;

/**  EUTelProcessorTrackingHelixTrackSearch
 * 
 *  If compiled with MARLIN_USE_AIDA 
 *  it creates histograms.
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Requires collection of hits from which track
 *  candidates are build.
 *
 *  <h4>Output</h4> 
 *  <li> Histograms.
 *  <li> Collection of track candidates.
 */

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelProcessorTrackingHelixTrackSearch::_histName::_numberTracksCandidatesHistName = "NumberTracksCandidates";
std::string EUTelProcessorTrackingHelixTrackSearch::_histName::_numberOfHitOnTrackCandidateHistName = "NumberOfHitsOnTrackCandidate";
#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

/** Default constructor */
EUTelProcessorTrackingHelixTrackSearch::EUTelProcessorTrackingHelixTrackSearch() :
Processor("EUTelProcessorTrackingHelixTrackSearch"),
_hitInputCollectionName("HitCollection"),
_trackCandidateHitsOutputCollectionName("TrackCandidatesCollection"),
_hotPixelCollectionName("hotpixel"),
_hotPixelMap(),
_trackFitter(0),
_tgeoFileName("TELESCOPE.root"),
_maxMissingHitsPerTrackCand(0),
_maxNTracks(10),
_eBeam(-1.),
_qBeam(-1.),
_eBeamUncertatinty(0.),
_beamSpread(2,-1.),
_histoInfoFileName("histoinfo.xml"),
_nProcessedRuns(0),
_nProcessedEvents(0),
_aidaHistoMap1D() {

    // Processor description
    _description = "EUTelProcessorTrackingHelixTrackSearch preforms track pattern recognition.";

    // TrackerHit input collection
    registerInputCollection(LCIO::TRACKERHIT,
            "HitInputCollectionName",
            "Input hits collection name",
            _hitInputCollectionName,
            std::string("HitCollection"));

    // Track candidate hits output collection
    registerOutputCollection(LCIO::TRACK,
            "TrackCandHitOutputCollectionName",
            "Output track candidates hits collection name",
            _trackCandidateHitsOutputCollectionName,
            std::string("TrackCandidateHitCollection"));


    // Optional processor parameters that define track finder settings

    // Geometry definition
    
    registerOptionalParameter("GeometryFilename", "Name of the TGeo geometry definition file", _tgeoFileName, std::string("TELESCOPE.root"));
    
    // Fitter settings
    
    registerOptionalParameter("MaxMissingHitsPerTrack", "Maximal number of missing hits on a track candidate",
            _maxMissingHitsPerTrackCand, static_cast<int> (0)); // Search full-length tracks by default

    registerOptionalParameter("MaxNTracksPerEvent", "Maximal number of track candidates to be found in events",
            _maxNTracks, static_cast<int> (100));
    
    registerOptionalParameter("BeamEnergy", "Beam energy [GeV]", _eBeam, static_cast<double> (4.0));

    registerOptionalParameter("BeamCharge", "Beam charge [e]", _qBeam, static_cast<double> (-1));
    
    registerOptionalParameter("BeamSpread", "Angular spread of the beam (horizontal,vertical) [mrad] (for beam constraint). "
                                            "No beam constraints if negative values are supplied.",
                               _beamSpread, EVENT::FloatVec(2,0.) );
    
    registerOptionalParameter("BeamEnergyUncertainty", "Uncertainty of beam energy [%]", _eBeamUncertatinty, static_cast<double> (0.) );

    // Histogram information

    registerOptionalParameter("HistogramInfoFilename", "Name of histogram info xml file", _histoInfoFileName, std::string("histoinfo.xml"));

    /**    @TODO must be a part of a separate data structure 
     *     Do we need this at all
     */
    registerOptionalParameter("HotPixelCollectionName", "Name of the hot pixel collection ",
            _hotPixelCollectionName, static_cast<std::string> ("hotpixel"));

}

void EUTelProcessorTrackingHelixTrackSearch::init() {

    streamlog_out(DEBUG) << "EUTelProcessorTrackingHelixTrackSearch::init( )" << std::endl;

    // usually a good idea to
    printParameters();

    // Reset counters
    _nProcessedRuns = 0;
    _nProcessedEvents = 0;

    // Getting access to geometry description
    geo::gGeometry().initializeTGeoDescription(_tgeoFileName);
    
    // Instantiate track finder. This is a working horse of the processor.
    {
        streamlog_out(DEBUG) << "Initialisation of track finder" << std::endl;

        EUTelKalmanFilter* Finder = new EUTelKalmanFilter("KalmanTrackFinder");
        if (!Finder) {
            streamlog_out(ERROR) << "Can't allocate an instance of EUTelExhaustiveTrackFinder. Stopping ..." << std::endl;
            throw UnknownDataTypeException("Track finder was not created");
        }

        Finder->setAllowedMissingHits( _maxMissingHitsPerTrackCand );
        Finder->setBeamMomentum( _eBeam );
        Finder->setBeamCharge( _qBeam );
        Finder->setBeamMomentumUncertainty( _eBeamUncertatinty );
        Finder->setBeamSpread( _beamSpread );
        
        _trackFitter = Finder;
    }

    // Book histograms
    bookHistograms();
}

void EUTelProcessorTrackingHelixTrackSearch::processRunHeader(LCRunHeader* run) {

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

    _nProcessedRuns++;
}

void EUTelProcessorTrackingHelixTrackSearch::processEvent(LCEvent * evt) {

    if (isFirstEvent()) {
        // Fill hot pixel map   
        _hotPixelMap = Utility::FillHotPixelMap(evt, _hotPixelCollectionName);
    }

    EUTelEventImpl * event = static_cast<EUTelEventImpl*> (evt);

    // Do not process last or unknown events
    
    if (event->getEventType() == kEORE) {
        streamlog_out(DEBUG4) << "EORE found: nothing else to do." << std::endl;
        return;
    } else if (event->getEventType() == kUNKNOWN) {
        streamlog_out(WARNING2) << "Event number " << event->getEventNumber() << " in run " << event->getRunNumber()
                << " is of unknown type. Continue considering it as a normal Data Event." << std::endl;
    }

    // Try to access collection

    LCCollection* col = NULL;
    try {
        col = evt->getCollection(_hitInputCollectionName);
    } catch (DataNotAvailableException e) {
        streamlog_out(WARNING) << _hitInputCollectionName << " collection not available" << std::endl;
        throw marlin::SkipEventException(this);
    }

    // this will only be entered if the collection is available
    if (col != NULL) {
        streamlog_out(DEBUG2) << "EUTelProcessorTrackingHelixTrackSearch::processEvent()" << std::endl;

        // Prepare hits for track finder
        EVENT::TrackerHitVec allHitsVec;
        FillHits(evt, col, allHitsVec);

        // Search tracks
        streamlog_out(DEBUG1) << "Event #" << _nProcessedEvents << std::endl;
        streamlog_out(DEBUG1) << "Initialising hits for _theFinder..." << std::endl;
        
        static_cast<EUTelKalmanFilter*>(_trackFitter)->setHits(allHitsVec);
        static_cast<EUTelKalmanFilter*>(_trackFitter)->initialise();
        streamlog_out(DEBUG1) << "Trying to find tracks..." << endl;
        _trackFitter->FitTracks();
        streamlog_out(DEBUG1) << "Retrieving track candidates..." << endl;
        vector< EUTelTrackImpl* >& trackCandidates = static_cast<EUTelKalmanFilter*>(_trackFitter)->getTracks();
       
        const int nTracks = trackCandidates.size();
        static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ _histName::_numberTracksCandidatesHistName ]) -> fill(nTracks);
        streamlog_out(DEBUG1) << "Event #" << _nProcessedEvents << endl;
        streamlog_out(DEBUG1) << "Track finder " << _trackFitter->GetName() << " found  " << nTracks << endl;
        
        int nHitsOnTrack = 0;
        vector< EUTelTrackImpl* >::const_iterator itrk;
        for ( itrk = trackCandidates.begin() ; itrk != trackCandidates.end(); ++itrk ) {
            const EVENT::TrackerHitVec& trkHits = (*itrk)->getTrackerHits();
            nHitsOnTrack = trkHits.size();
            EVENT::TrackerHitVec::const_iterator itTrkHits;
            streamlog_out(MESSAGE0) << "Track hits start:==============" << std::endl;
            for ( itTrkHits = trkHits.begin() ; itTrkHits != trkHits.end(); ++itTrkHits ) {
                const double* uvpos = (*itTrkHits)->getPosition();
                streamlog_out(MESSAGE0) << "Hit (id=" << (*itrk)->id() << ") local(u,v) coordinates: (" << uvpos[0] << "," << uvpos[1] << ")" << std::endl;
            }
            streamlog_out(MESSAGE0) << "Track hits end:==============" << std::endl;
            static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ _histName::_numberOfHitOnTrackCandidateHistName ]) -> fill(nHitsOnTrack);
        }
        
//        std::vector< EVENT::TrackerHitVec >& trkHits = static_cast<EUTelKalmanFilter*>(_trackFitter)->getTrackCandidates();
//        int nHitsOnTrack = 0;
//        std::vector< EVENT::TrackerHitVec >::const_iterator itrk;
//        for ( itrk = trkHits.begin() ; itrk != trkHits.end(); ++itrk ) {
//            nHitsOnTrack = (*itrk).size();
//            static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ _histName::_numberOfHitOnTrackCandidateHistName ]) -> fill(nHitsOnTrack);
//        }
        
        // Write output collection
//        addTrackCandidateToCollection(evt, trackCandidates);

//        _trackFitter->Reset();

        _nProcessedEvents++;

        if (isFirstEvent()) _isFirstEvent = false;
    } // if (col != NULL)
}

/** 
 * Fill hits in a vector.
 * 
 * @TODO LCEvent* evt must be removed from the list of arguments, idealy.
 * 
 * @param evt event pointer
 * @param collection hits collection pointer
 * @param [out] allHitsVec vector of hits in event
 */
void EUTelProcessorTrackingHelixTrackSearch::FillHits(LCEvent * evt,
        LCCollection* collection,
        EVENT::TrackerHitVec& allHitsVec) const {

    allHitsVec.clear();
    allHitsVec.resize( 0 );
    
    // loop over all hits in collection
    for (int iHit = 0; iHit < collection->getNumberOfElements(); iHit++) {

        TrackerHitImpl * hit = static_cast<TrackerHitImpl*> (collection->getElementAt(iHit));
        if (Utility::HitContainsHotPixels(hit, _hotPixelMap)) {
            streamlog_out(DEBUG3) << "Hit " << iHit << " contains hot pixels; skip this one. " << std::endl;
            continue;
        }

        LCObjectVec clusterVector = hit->getRawHits();

        EUTelVirtualCluster * cluster;
        cluster = Utility::GetClusterFromHit(hit);

        if (hit->getType() == kEUTelSparseClusterImpl) {
            // ok the cluster is of sparse type, but we also need to know
            // the kind of pixel description used. This information is
            // stored in the corresponding original data collection.
            LCCollectionVec * sparseClusterCollectionVec = dynamic_cast<LCCollectionVec *> (evt->getCollection("original_zsdata"));

            TrackerDataImpl * oneCluster = dynamic_cast<TrackerDataImpl*> (sparseClusterCollectionVec->getElementAt(0));
            CellIDDecoder<TrackerDataImpl > anotherDecoder(sparseClusterCollectionVec);
            SparsePixelType pixelType = static_cast<SparsePixelType> (static_cast<int> (anotherDecoder(oneCluster)["sparsePixelType"]));

            // now we know the pixel type. So we can properly create a new
            // instance of the sparse cluster
            if (pixelType == kEUTelSimpleSparsePixel) {
                cluster = new EUTelSparseClusterImpl< EUTelSimpleSparsePixel > (static_cast<TrackerDataImpl *> (clusterVector[ 0 ]));
            } else if (pixelType == kEUTelAPIXSparsePixel) {
                cluster = new EUTelSparseClusterImpl<EUTelAPIXSparsePixel > (static_cast<TrackerDataImpl *> (clusterVector[ 0 ]));
            } else {
                streamlog_out(ERROR4) << "Unknown pixel type.  Sorry for quitting." << std::endl;
                throw UnknownDataTypeException("Pixel type unknown");
            }
        }

        const int localSensorID = Utility::GuessSensorID( hit );  // localSensorID == -1, if detector ID was not found
        if ( localSensorID >= 0 ) allHitsVec.push_back( hit );
        
        delete cluster; // <--- destroying the cluster
    } // end loop over all hits in collection
    
}

/**
 * Dumpt track candidate into lcio collection.
 * 
 * @param evt event pointer
 * @param trackCandidates  vectors of hits assigned to track candidates
 */
void EUTelProcessorTrackingHelixTrackSearch::addTrackCandidateToCollection(LCEvent* evt, const vector< EVENT::TrackerHitVec >& trackCandidates) {
    // Prepare output collection
    LCCollectionVec * trkCandCollection = 0;
    try {
        trkCandCollection = new LCCollectionVec(LCIO::TRACK);
        LCFlagImpl flag(trkCandCollection->getFlag());
        flag.setBit( LCIO::TRBIT_HITS );
        trkCandCollection->setFlag( flag.getFlag( ) );
    } catch (...) {
        streamlog_out(WARNING2) << "Can't allocate output collection" << endl;
    }

    // Fill track parameters with nonsense except of hits
    for (size_t itrk = 0; itrk < trackCandidates.size(); itrk++) {
        TrackImpl* trackcand = new TrackImpl;           // Don't free it manually, because it is owned by trkCandCollection

        // Assign hits to LCIO TRACK
        EVENT::TrackerHitVec trkcandhits = trackCandidates[itrk];
        for (size_t ihit = 0; ihit < trkcandhits.size(); ihit++) {
            TrackerHitImpl* hit = static_cast<TrackerHitImpl*>(trkcandhits[ihit]);
            trackcand->addHit( hit );
        }
        streamlog_out( DEBUG1 ) << "Track has " << trackcand->getTrackerHits().size() << " hits" << endl;
        trkCandCollection->push_back(trackcand);

    } // for (size_t itrk = 0; itrk < trackCandidates.size(); itrk++) 

    // Write track candidates collection
    try {
        streamlog_out(DEBUG1) << "Getting collection " << _trackCandidateHitsOutputCollectionName << endl;
        evt->getCollection(_trackCandidateHitsOutputCollectionName);
    } catch (...) {
        streamlog_out(DEBUG1) << "Adding collection " << _trackCandidateHitsOutputCollectionName << endl;
        evt->addCollection(trkCandCollection, _trackCandidateHitsOutputCollectionName);
    }

}

void EUTelProcessorTrackingHelixTrackSearch::check(LCEvent * /*evt*/) {
    // nothing to check here
}

void EUTelProcessorTrackingHelixTrackSearch::end() {

    delete _trackFitter;
    
    streamlog_out(MESSAGE) << "EUTelProcessorTrackingHelixTrackSearch::end()  " << name()
            << " processed " << _nProcessedEvents << " events in " << _nProcessedRuns << " runs "
            << std::endl;

}

void EUTelProcessorTrackingHelixTrackSearch::bookHistograms() {
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
            streamlog_out(ERROR2) << "Problem booking the " << (_histName::_numberTracksCandidatesHistName) << endl;
            streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        }
        
        const int hitsNBin = 10;
        const double hitsMin = -0.5;
        const double hitsMax = 9.5;
        
        histoInfo = histoMgr->getHistogramInfo(_histName::_numberOfHitOnTrackCandidateHistName);        
        NBin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : hitsNBin;
        XMin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : hitsMin;
        XMax = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : hitsMax;
        
        // Number of hit per track candidate
        AIDA::IHistogram1D * numberHitsOnTrackCandidates =
                marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_numberOfHitOnTrackCandidateHistName, NBin,
                XMin, XMax);

        if (numberHitsOnTrackCandidates) {
            numberHitsOnTrackCandidates->setTitle("Number of hits on track candidates;N hits;N tracks");
            _aidaHistoMap1D.insert(make_pair(_histName::_numberOfHitOnTrackCandidateHistName, numberHitsOnTrackCandidates));
        } else {
            streamlog_out(ERROR2) << "Problem booking the " << (_histName::_numberOfHitOnTrackCandidateHistName) << endl;
            streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        }
        
    } catch (lcio::Exception& e) {
        streamlog_out(WARNING2) << "Can't allocate histograms. Continue without histogramming" << endl;
    }
#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
}

int EUTelProcessorTrackingHelixTrackSearch::getMaxTrackCandidates() const {
    return _maxNTracks;
}

int EUTelProcessorTrackingHelixTrackSearch::getAllowedMissingHits() const {
    return _maxMissingHitsPerTrackCand;
}
