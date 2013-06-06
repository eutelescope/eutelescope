#include "EUTelProcessorTrackingExhaustiveTrackSearch.h"

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
#include "EUTelExhaustiveTrackFinder.h"
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

/**  EUTelProcessorTrackingExhaustiveTrackSearch
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
std::string EUTelProcessorTrackingExhaustiveTrackSearch::_histName::_numberTracksCandidatesHistName = "NumberTracksCandidates";
std::string EUTelProcessorTrackingExhaustiveTrackSearch::_histName::_numberOfHitOnTrackCandidateHistName = "NumberOfHitsOnTrackCandidate";
#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

/** Default constructor */
EUTelProcessorTrackingExhaustiveTrackSearch::EUTelProcessorTrackingExhaustiveTrackSearch() :
Processor("EUTelProcessorTrackingExhaustiveTrackSearch"),
_hitInputCollectionName("HitCollection"),
_trackCandidateHitsOutputCollectionName("TrackCandidatesCollection"),
_hotPixelCollectionName("hotpixel"),
_hotPixelMap(),
_trackFinder(0),
_maxMissingHitsPerTrackCand(0),
_maxNTracks(10),
_finderMode(2),
_residualsXMin(6, -300.),
_residualsYMin(6, -300.),
_residualsXMax(6, 300.),
_residualsYMax(6, 300.),
_residualsRMax(6, 300.),
_histoInfoFileName("histoinfo.xml"),
_nProcessedRuns(0),
_nProcessedEvents(0),
_aidaHistoMap1D() {

    // Processor description
    _description = "EUTelProcessorTrackingExhaustiveTrackSearch preforms track pattern recognition.";

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

    registerOptionalParameter("MaxMissingHitsPerTrack", "Maximal number of missing hits on a track candidate",
            _maxMissingHitsPerTrackCand, static_cast<int> (0)); // Search full-length tracks by default

    registerOptionalParameter("MaxNTracksPerEvent", "Maximal number of track candidates to be found in events",
            _maxNTracks, static_cast<int> (10));

    registerOptionalParameter("FinderMode", "Finder mode. Possible values are 1, 2",
            _finderMode, static_cast<int> (2));

    const FloatVec MinimalResidualsX(6, -300.); // assumes there are at most 10 planes. But this doesn't matter.
    registerOptionalParameter("ResidualsXMin",
            "Minimal values of the hit residuals in the X direction for a track. "
            "Note: these numbers are ordered according to the z position of "
            "the sensors and NOT according to the sensor id.",
            _residualsXMin, MinimalResidualsX);

    const FloatVec MinimalResidualsY(6, -300.); // assumes there are at most 10 planes. But this doesn't matter.
    registerOptionalParameter("ResidualsYMin", "Minimal values of the hit residuals in the Y direction for a track. "
            "Note: these numbers are ordered according to the z position of "
            "the sensors and NOT according to the sensor id.",
            _residualsYMin, MinimalResidualsY);

    const FloatVec MaximalResidualsX(6, 300.); // assumes there are at most 10 planes. But this doesn't matter.
    registerOptionalParameter("ResidualsXMax", "Maximal values of the hit residuals in the X direction for a track. "
            "Note: these numbers are ordered according to the z position of "
            "the sensors and NOT according to the sensor id.",
            _residualsXMax, MaximalResidualsX);

    const FloatVec MaximalResidualsY(6, 300.); // assumes there are at most 10 planes. But this doesn't matter.
    registerOptionalParameter("ResidualsYMax", "Maximal values of the hit residuals in the Y direction for a track. "
            "Note: these numbers are ordered according to the z position of "
            "the sensors and NOT according to the sensor id.",
            _residualsYMax, MaximalResidualsY);

    const FloatVec MaximalResidualsR(6, 300.);
    registerOptionalParameter("ResidualsRMax", "Maximal allowed distance between hits entering the recognition step "
            "per 10 cm space between the planes. One value for each neighbor planes. "
            "DistanceMax will be used for each pair if this vector is empty.",
            _residualsRMax, MaximalResidualsR);

    // Histogram information

    registerOptionalParameter("HistogramInfoFilename", "Name of histogram info xml file", _histoInfoFileName, std::string("histoinfo.xml"));

    /**    @TODO must be a part of a separate data structure 
     *     Do we need this at all
     */
    registerOptionalParameter("HotPixelCollectionName", "Name of the hot pixel collection ",
            _hotPixelCollectionName, static_cast<std::string> ("hotpixel"));

}

void EUTelProcessorTrackingExhaustiveTrackSearch::init() {

    streamlog_out(DEBUG) << "EUTelProcessorTrackingExhaustiveTrackSearch::init( )" << std::endl;

    // usually a good idea to
    printParameters();

    // Reset counters
    _nProcessedRuns = 0;
    _nProcessedEvents = 0;


    // Instantiate track finder. This is a working horse of the processor.
    {
        streamlog_out(DEBUG) << "Initialisation of track finder" << std::endl;

        EUTelExhaustiveTrackFinder* Finder = new EUTelExhaustiveTrackFinder("TrackFinder", getAllowedMissingHits(), getMaxTrackCandidates());
        if (!Finder) {
            streamlog_out(ERROR) << "Can't allocate an instance of EUTelExhaustiveTrackFinder. Stopping ..." << std::endl;
            throw UnknownDataTypeException("Track finder was not created");
        }
        Finder->SetMode(getFinderMode());
        Finder->SetDistanceMaxVec(_residualsRMax);
        Finder->SetResidualsYMax(_residualsYMax);
        Finder->SetResidualsYMin(_residualsYMin);
        Finder->SetResidualsXMax(_residualsXMax);
        Finder->SetResidualsXMin(_residualsXMin);

        _trackFinder = Finder;
    }

    // Book histograms
    bookHistograms();
}

void EUTelProcessorTrackingExhaustiveTrackSearch::processRunHeader(LCRunHeader* run) {

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

void EUTelProcessorTrackingExhaustiveTrackSearch::processEvent(LCEvent * evt) {

    if (isFirstEvent()) {
        // Fill hot pixel map   
        _hotPixelMap = Utility::FillHotPixelMap(evt, _hotPixelCollectionName);
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
        col = evt->getCollection(_hitInputCollectionName);
    } catch (DataNotAvailableException e) {
        streamlog_out(WARNING) << _hitInputCollectionName << " collection not available" << std::endl;
        throw marlin::SkipEventException(this);
    }

    // this will only be entered if the collection is available
    if (col != NULL) {
        streamlog_out(DEBUG2) << "EUTelProcessorTrackingExhaustiveTrackSearch::processEvent()" << endl;

        // Prepare hits for track finder
        map< int, EVENT::TrackerHitVec > allHitsArray;
        FillHits(evt, col, allHitsArray);

        // Search tracks
        int nTracks = 0;
        int nGoodTracks = 0;

        vector< EVENT::TrackerHitVec > trackCandidates;

        streamlog_out(DEBUG1) << "Event #" << _nProcessedEvents << std::endl;
        streamlog_out(DEBUG1) << "Initialising hits for _theFinder..." << std::endl;
        TrackerHitVec::iterator itHit;
        const EVENT::IntVec sensorIDVecZOrder = geo::gGeometry()._sensorIDVecZOrder;
        const size_t nplanes = geo::gGeometry()._nPlanes;
        for (size_t iLayer = 0; iLayer < nplanes; iLayer++) {
            streamlog_out(DEBUG1) << "Plane " << sensorIDVecZOrder.at(iLayer) << endl;
            streamlog_out(DEBUG1) << allHitsArray[sensorIDVecZOrder.at(iLayer)].size() << " hit(s) in plane" << endl;
            for (itHit = allHitsArray[sensorIDVecZOrder.at(iLayer)].begin(); itHit != allHitsArray[sensorIDVecZOrder.at(iLayer)].end(); ++itHit) {
                streamlog_out(DEBUG1) << *itHit << endl;
                streamlog_out(DEBUG1) << setw(10) << Utility::GuessSensorID(static_cast<TrackerHitImpl*>(*itHit)) << endl;
                streamlog_out(DEBUG1) << setw(10) << (*itHit)->id() << endl;
                streamlog_out(DEBUG1) << setw(10) << (*itHit)->getPosition()[0]
                        << setw(10) << (*itHit)->getPosition()[1]
                        << setw(10) << (*itHit)->getPosition()[2] << endl;
            }
        }

        _trackFinder->SetAllHits(allHitsArray);
        streamlog_out(DEBUG1) << "Trying to find tracks..." << endl;
        EUTelTrackFinder::SearchResult searchResult = _trackFinder->SearchTracks();
        streamlog_out(DEBUG1) << "Search results retrieved..." << endl;
        streamlog_out(DEBUG1) << "Search results = " << (int) searchResult << endl;
        streamlog_out(DEBUG1) << "Retrieving track candidates..." << endl;
        trackCandidates = _trackFinder->GetTrackCandidates();

        nTracks = (int) trackCandidates.size();
        static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ _histName::_numberTracksCandidatesHistName ]) -> fill(nTracks);
        streamlog_out(DEBUG1) << "Event #" << _nProcessedEvents << endl;
        streamlog_out(DEBUG1) << "Track finder " << _trackFinder->GetName() << " found  " << nTracks << endl;
        
        int nHitsOnTrack = 0;
        vector< EVENT::TrackerHitVec >::const_iterator itrk = trackCandidates.begin();
        for ( ; itrk != trackCandidates.end(); ++itrk ) {
            nHitsOnTrack = (int)itrk->size();
            static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ _histName::_numberOfHitOnTrackCandidateHistName ]) -> fill(nHitsOnTrack);
        }

        
        // Write output collection
        addTrackCandidateToCollection(evt, trackCandidates);

        _trackFinder->Reset(); // Breaks Fitter track candidates reference???!!!

        _nProcessedEvents++;

        if (isFirstEvent()) _isFirstEvent = false;
    }
}

// @TODO LCEvent* evt must be removed from the list of arguments

void EUTelProcessorTrackingExhaustiveTrackSearch::FillHits(LCEvent * evt,
        LCCollection* collection,
        map< int, EVENT::TrackerHitVec >& allHitsArray) const {

    int layerIndex = -1;

    // loop over all hits in collection
    for (int iHit = 0; iHit < collection->getNumberOfElements(); iHit++) {

        TrackerHitImpl * hit = static_cast<TrackerHitImpl*> (collection->getElementAt(iHit));
        ///code hoes here;
        if (Utility::HitContainsHotPixels(hit, _hotPixelMap)) {
            streamlog_out(DEBUG3) << "Hit " << iHit << " contains hot pixels; skip this one. " << std::endl;
            continue;
        }

        LCObjectVec clusterVector = hit->getRawHits();

        EUTelVirtualCluster * cluster;
        cluster = Utility::GetClusterFromHit(hit);

        //Warning! It is very bad idea to do like this. This code must ideally be implemented in GetClusterFromHit
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

        // @TODO add support of kEUTelAPIXClusterImpl type
        if (hit->getType() == kEUTelAPIXClusterImpl) {
            streamlog_out(MESSAGE4) << "kEUTelAPIXClusterImpl is not supported at the moment" << std::endl;
            continue;
            //                
            //            TrackerDataImpl * clusterFrame = static_cast<TrackerDataImpl*> (clusterVector[0]);
            //            EUTelSparseClusterImpl< EUTelAPIXSparsePixel > *apixCluster = new EUTelSparseClusterImpl< EUTelAPIXSparsePixel > (clusterFrame);
            //            int sensorID = apixCluster->getDetectorID();
            //            bool skipHit = false;
            //            for (size_t iPixel = 0; iPixel < apixCluster->size(); ++iPixel) {
            //                EUTelAPIXSparsePixel apixPixel;
            //                apixCluster->getSparsePixelAt(iPixel, &apixPixel);
            //                int pixelX = apixPixel.getXCoord();
            //                int pixelY = apixPixel.getYCoord();
            //                skipHit = !_rect.isInside(sensorID, pixelX, pixelY);
            //            }
            //
            //            if (skipHit) {
            //                streamlog_out(MESSAGE4) << "Skipping cluster due to rectangular cut!" << std::endl;
            //                continue;
            //            }
            //            if (skipHit) {
            //                streamlog_out(MESSAGE4) << "TEST: This message should never appear!" << std::endl;
            //            }
        }

        unsigned int localSensorID = Utility::GuessSensorID(hit);
        layerIndex = geo::gGeometry()._sensorIDVecMap[ localSensorID ];
        allHitsArray[layerIndex].push_back(hit);

        delete cluster; // <--- destroying the cluster
    } // end loop over all hits in collection

}

void EUTelProcessorTrackingExhaustiveTrackSearch::addTrackCandidateToCollection(LCEvent* evt, const vector< EVENT::TrackerHitVec >& trackCandidates) {
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
        TrackImpl* trackcand = new TrackImpl;           // Don't free it manually, because it is (presumably) owned by trkCandCollection
        trackcand->setD0(999.);
        trackcand->setZ0(999.);
        trackcand->setTanLambda(999.);
        trackcand->setOmega(999.);
        trackcand->setPhi(999.);
        trackcand->setChi2(999.);
        trackcand->setNdf(999);     

        // Assign hits to LCIO TRACK
        EVENT::TrackerHitVec trkcandhits = trackCandidates[itrk];
        for (size_t ihit = 0; ihit < trkcandhits.size(); ihit++) {
            TrackerHitImpl* hit = static_cast<TrackerHitImpl*>(trkcandhits[ihit]);
//            trackcand->addHit( assignCov(hit) );
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

void EUTelProcessorTrackingExhaustiveTrackSearch::check(LCEvent * evt) {
    // nothing to check here
}

void EUTelProcessorTrackingExhaustiveTrackSearch::end() {

    streamlog_out(MESSAGE) << "EUTelProcessorTrackingExhaustiveTrackSearch::end()  " << name()
            << " processed " << _nProcessedEvents << " events in " << _nProcessedRuns << " runs "
            << std::endl;

}

void EUTelProcessorTrackingExhaustiveTrackSearch::bookHistograms() {
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
        streamlog_out(WARNING2) << "Can't allocate histgrams. Continue without histogramming" << endl;
    }
#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
}

int EUTelProcessorTrackingExhaustiveTrackSearch::getMaxTrackCandidates() const {
    return _maxNTracks;
}

int EUTelProcessorTrackingExhaustiveTrackSearch::getAllowedMissingHits() const {
    return _maxMissingHitsPerTrackCand;
}

int EUTelProcessorTrackingExhaustiveTrackSearch::getFinderMode() const {
    return _finderMode;
}
