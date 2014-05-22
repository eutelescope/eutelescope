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
std::string EUTelProcessorTrackingHelixTrackSearch::_histName::_HitOnTrackCandidateHistName = "HitsOnTrackCandidate";
#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

/** Default constructor */
EUTelProcessorTrackingHelixTrackSearch::EUTelProcessorTrackingHelixTrackSearch() :
Processor("EUTelProcessorTrackingHelixTrackSearch"),
_hitInputCollectionName("HitCollection"),
_trackCandidateHitsOutputCollectionName("TrackCandidatesCollection"),
//_hotPixelCollectionName("hotpixel"),
//_hotPixelMap(),
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

    // TrackerHit output collection
    registerInputCollection(LCIO::TRACKERHIT,
            "TrackerHitOutputCollectionName",
            "Pattern recognition track - output hits collection name",
            _hitFittedOutputCollectionName,
            std::string("HitFittedCollection"));

    // Track candidate hits output collection
    registerOutputCollection(LCIO::TRACK,
            "TrackCandHitOutputCollectionName",
            "Output track candidates hits collection name",
            _trackCandidateHitsOutputCollectionName,
            std::string("TrackCandidateHitCollection"));

    // Optional processor parameters that define track finder settings
     
    // Fitter settings
    
    registerOptionalParameter("MaxMissingHitsPerTrack", "Maximal number of missing hits on a track candidate",
            _maxMissingHitsPerTrackCand, static_cast<int> (0)); // Search full-length tracks by default
    
    registerOptionalParameter("MaxNTracksPerEvent", "Maximal number of track candidates to be found in events",
            _maxNTracks, static_cast<int> (100));
    
    registerOptionalParameter("NumberOfPlanesToProject", "Maximal number of track candidates to be found in events",
            _planesProject, static_cast<int> (1));
    
    registerOptionalParameter("ResidualsRMax", "Maximal allowed distance between hits entering the recognition step "
            "per 10 cm space between the planes. One value for each neighbor planes. "
            "DistanceMax will be used for each pair if this vector is empty. Units are mm.",
            _residualsRMax, static_cast<double> (0.4));
    
    registerOptionalParameter("BeamEnergy", "Beam energy [GeV]", _eBeam, static_cast<double> (4.0));

    registerOptionalParameter("BeamCharge", "Beam charge [e]", _qBeam, static_cast<double> (-1));
    
    registerOptionalParameter("BeamSpread", "Angular spread of the beam (horizontal,vertical) [mrad] (for beam constraint). "
                                            "No beam constraints if negative values are supplied.",
                               _beamSpread, EVENT::FloatVec(2,100.) );
    
    registerOptionalParameter("BeamEnergyUncertainty", "Uncertainty of beam energy [%]", _eBeamUncertatinty, static_cast<double> (0.) );

    // Histogram information

    registerOptionalParameter("HistogramInfoFilename", "Name of histogram info xml file", _histoInfoFileName, std::string("histoinfo.xml"));

    /**    @TODO must be a part of a separate data structure 
     *     Do we need this at all
     */
//    registerOptionalParameter("HotPixelCollectionName", "Name of the hot pixel collection ",
//            _hotPixelCollectionName, static_cast<std::string> ("hotpixel"));

}

void EUTelProcessorTrackingHelixTrackSearch::init() {

    streamlog_out(DEBUG) << "EUTelProcessorTrackingHelixTrackSearch::init( )" << std::endl;

    // usually a good idea to
    printParameters();

    // Reset counters
    _nProcessedRuns = 0;
    _nProcessedEvents = 0;

    // Getting access to geometry description
//    geo::gGeometry().initializeTGeoDescription(_tgeoFileName);
    std::string name("test.root");
    geo::gGeometry().initializeTGeoDescription(name,false);
    // Instantiate track finder. This is a working horse of the processor.
    {
        streamlog_out(DEBUG) << "Initialisation of track finder" << std::endl;

        EUTelKalmanFilter* Finder = new EUTelKalmanFilter("KalmanTrackFinder");
        if (!Finder) {
            streamlog_out(ERROR) << "Can't allocate an instance of EUTelExhaustiveTrackFinder. Stopping ..." << std::endl;
            throw UnknownDataTypeException("Track finder was not created");
        }

        Finder->setAllowedMissingHits( _maxMissingHitsPerTrackCand );
        Finder->setWindowSize( _residualsRMax );
        Finder->setBeamMomentum( _eBeam );
        Finder->setBeamCharge( _qBeam );
        Finder->setPlanesProject( _planesProject );
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

    EUTelEventImpl * event = static_cast<EUTelEventImpl*> (evt);

    // Do not process last or unknown events
    
    if (event->getEventType() == kEORE) {
        streamlog_out(DEBUG4) << "EORE found: nothing else to do." << std::endl;
        return;
    } else if (event->getEventType() == kUNKNOWN) {
        streamlog_out(WARNING2) << "Event number " << event->getEventNumber() << " in run " << event->getRunNumber()
                << " is of unknown type. Continue considering it as a normal Data Event." << std::endl;
    }

    // Try to access Input collection
    LCCollection* hitMeasuredCollection = NULL;
    try {
        hitMeasuredCollection = evt->getCollection(_hitInputCollectionName);
        streamlog_out(DEBUG1) << "collection : " <<_hitInputCollectionName << " retrieved" << std::endl;
    } catch (DataNotAvailableException e) {
        streamlog_out(WARNING2) << _hitInputCollectionName << " collection not available" << std::endl;
        throw marlin::SkipEventException(this);
    }

    // this will only be entered if the collection is available
    if ( hitMeasuredCollection == NULL) {
        streamlog_out(DEBUG2) << "EUTelProcessorTrackingHelixTrackSearch :: processEvent() hitMeasuredCollection is void" << std::endl;
    }
    else if ( hitMeasuredCollection != NULL) {
        streamlog_out(DEBUG2) << "EUTelProcessorTrackingHelixTrackSearch::processEvent()" << std::endl;

        // Prepare hits for track finder
        EVENT::TrackerHitVec allHitsVec;
        FillHits(evt, hitMeasuredCollection, allHitsVec);

	streamlog_out(MESSAGE0) << "All hits in event start:==============" << std::endl;
	EVENT::TrackerHitVec::const_iterator itHits;
	for ( itHits = allHitsVec.begin() ; itHits != allHitsVec.end(); ++itHits ) {
		const double* uvpos = (*itHits)->getPosition();
                    const int sensorID = geo::gGeometry().getSensorID( static_cast<IMPL::TrackerHitImpl*> (*itHits) );
//                    const int usensorID = Utility::GuessSensorID( static_cast<EVENT::TrackerHit*> (*itHits) );
		streamlog_out(MESSAGE0) << "Hit (id=" << setw(3) << sensorID << ") local(u,v) coordinates: (" 
                       << setw(7) << setprecision(4) << uvpos[0] << "," << setw(7) << setprecision(4) << uvpos[1] << ")" << std::endl;
	}
	streamlog_out(MESSAGE0) << "All hits in event end:==============" << std::endl;

        // Search tracks
        streamlog_out(DEBUG1) << "Event #" << _nProcessedEvents << std::endl;
        streamlog_out(DEBUG1) << "Initialising hits for _theFinder..." << std::endl;
        static_cast<EUTelKalmanFilter*>(_trackFitter)->setHits(allHitsVec);
        bool isReady = static_cast<EUTelKalmanFilter*>(_trackFitter)->initialise();
        if ( isReady )  {
            streamlog_out( DEBUG1 ) << "Trying to find tracks..." << endl;
//replace this one
//            _trackFitter->FitTracks( );
// with two new ones:
            _trackFitter->SearchTrackCandidates( );
            _trackFitter->PruneTrackCandidates( );

            streamlog_out( DEBUG1 ) << "Retrieving track candidates..." << endl;

            // retrieve the result of the PR :: track candidates ::
            vector< IMPL::TrackImpl* >& trackCandidates    = static_cast < EUTelKalmanFilter* > ( _trackFitter )->getTracks( );

            // not needed any more ? ::
            EVENT::TrackerHitVec trackCandidateHitFitted = static_cast < EUTelKalmanFilter* > ( _trackFitter )->getHitFittedVec( );

            // plot only :
            plotHistos( trackCandidates );

            // Write output collection
//            addTrackCandidateHitFittedToCollection( evt, trackCandidateHitFitted );
 
           // Write output collection
            addTrackCandidateToCollection1( evt, trackCandidates );

        }
        _nProcessedEvents++;

        if (isFirstEvent()) _isFirstEvent = false;
    } // if ( hitMeasuredCollection != NULL)
}


/** 
 * Plot few histos.
 * 
 */
void EUTelProcessorTrackingHelixTrackSearch::plotHistos( vector< IMPL::TrackImpl* >& trackCandidates )  {

            const int nTracks = trackCandidates.size( );
 
            static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_numberTracksCandidatesHistName ] ) -> fill( nTracks );
            streamlog_out( MESSAGE2 ) << "Event #" << _nProcessedEvents << endl;
            streamlog_out( MESSAGE2 ) << "Track finder " << _trackFitter->GetName( ) << " found  " << nTracks << endl;

            int nHitsOnTrack = 0;
            vector< IMPL::TrackImpl* >::const_iterator itrk;
            for ( itrk = trackCandidates.begin( ) ; itrk != trackCandidates.end( ); ++itrk ) {
                const EVENT::TrackerHitVec& trkHits = ( *itrk )->getTrackerHits( );
                nHitsOnTrack = trkHits.size( );
                EVENT::TrackerHitVec::const_iterator itTrkHits;
                streamlog_out( MESSAGE1 ) << "Track hits start:==============" << std::endl;
                for ( itTrkHits = trkHits.begin( ) ; itTrkHits != trkHits.end( ); ++itTrkHits ) {
                    const double* uvpos = ( *itTrkHits )->getPosition( );
                    const int sensorID = geo::gGeometry().getSensorID( static_cast<IMPL::TrackerHitImpl*> (*itTrkHits) );
//                    const int usensorID = Utility::GuessSensorID( static_cast<EVENT::TrackerHit*> (*itTrkHits) );
                    streamlog_out( MESSAGE1 ) << "Hit (id=" << setw(3) << sensorID << ") local(u,v) coordinates: (" 
                             << setw(7) << setprecision(4) << uvpos[0] << "," <<  setw(7) << setprecision(4) << uvpos[1] << ")";
                    double globalHit[] = {0.,0.,0.};
                    geo::gGeometry().local2Master( sensorID, uvpos, globalHit);
                    streamlog_out( MESSAGE1 ) << " WorldC: " << setw(7) << globalHit[0] << setw(7) << globalHit[1] << setw(7) << globalHit[2] ;
                    streamlog_out( MESSAGE1 )  << std::endl;
                    static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_HitOnTrackCandidateHistName ] ) -> fill( sensorID );
                }
                streamlog_out( MESSAGE1 ) << "Track hits end:==============" << std::endl;
                static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_numberOfHitOnTrackCandidateHistName ] ) -> fill( nHitsOnTrack );
            }
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
//        if (Utility::HitContainsHotPixels(hit, _hotPixelMap)) {
//            streamlog_out(DEBUG3) << "Hit " << iHit << " contains hot pixels; skip this one. " << std::endl;
//            continue;
//        }

        if( hit->getType() > 31 ) continue; // how do we mark measurement hits ? those that contain cluster information ??

        LCObjectVec clusterVector = hit->getRawHits();

        EUTelVirtualCluster * cluster;

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
            delete cluster; // <--- destroying the cluster
        }

        const int localSensorID = geo::gGeometry().getSensorID( static_cast<IMPL::TrackerHitImpl*> (hit) );
//       const int localSensorID = Utility::GuessSensorID( hit );  // localSensorID == -1, if detector ID was not found
        if ( localSensorID >= 0 ) allHitsVec.push_back( hit );
        
    } // end loop over all hits in collection
    
}

/**
 * Dump track candidate fitted hits into lcio collection.
 * 
 * @param evt event pointer
 * @param trackCandidateHitFirred  vector of hits from Pattern recognition
 */
void EUTelProcessorTrackingHelixTrackSearch::addTrackCandidateHitFittedToCollection(LCEvent* evt, EVENT::TrackerHitVec& trackCandidateHitFitted ) {
    // Prepare output collection

    // Try to access Output collection
    LCCollectionVec* hitFittedCollection = NULL;
    try {
        hitFittedCollection =  static_cast<LCCollectionVec*> ( evt->getCollection(_hitFittedOutputCollectionName) ) ;
        streamlog_out(DEBUG1) << "collection : " <<_hitFittedOutputCollectionName << " retrieved" << std::endl;
    } catch (DataNotAvailableException e) {
        streamlog_out(WARNING2) << _hitFittedOutputCollectionName << " collection not available, creating one ... " << std::endl;
        hitFittedCollection = new LCCollectionVec(LCIO::TRACKERHIT);
    }

    // Fill 
    EVENT::TrackerHitVec::iterator ihit;
    for ( ihit = trackCandidateHitFitted.begin(); ihit != trackCandidateHitFitted.end(); ++ihit ) {
        streamlog_out( MESSAGE1 ) << "hit " << (*ihit)->getType() << " hits" << endl;
        hitFittedCollection->push_back( (*ihit) );

    } // for (size_t itrk = 0; itrk < trackCandidates.size(); itrk++) 


    // Write track candidates collection
    try {
        streamlog_out(MESSAGE1) << "Getting collection " << _hitFittedOutputCollectionName << endl;
        evt->getCollection(_trackCandidateHitsOutputCollectionName);
    } catch (...) {
        streamlog_out(MESSAGE1) << "Adding collection " << _hitFittedOutputCollectionName << endl;
        evt->addCollection( hitFittedCollection, _hitFittedOutputCollectionName);
    }

}


/**
 * Dump track candidate into lcio collection.
 * 
 * @param evt event pointer
 * @param trackCandidates  vectors of hits assigned to track candidates
 */
void EUTelProcessorTrackingHelixTrackSearch::addTrackCandidateToCollection1(LCEvent* evt, std::vector< IMPL::TrackImpl* >& trackCandidates ) {
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
    std::vector< IMPL::TrackImpl* >::iterator itrk;
    for ( itrk = trackCandidates.begin(); itrk != trackCandidates.end(); ++itrk ) {
        streamlog_out( MESSAGE1 ) << "Track has " << (*itrk)->getTrackerHits().size() << " hits" << endl;
        trkCandCollection->push_back( (*itrk) );

    } // for (size_t itrk = 0; itrk < trackCandidates.size(); itrk++) 

    // Write track candidates collection
    try {
        streamlog_out(MESSAGE1) << "Getting collection " << _trackCandidateHitsOutputCollectionName << endl;
        evt->getCollection(_trackCandidateHitsOutputCollectionName);
    } catch (...) {
        streamlog_out(MESSAGE1) << "Adding collection " << _trackCandidateHitsOutputCollectionName << endl;
        evt->addCollection(trkCandCollection, _trackCandidateHitsOutputCollectionName);
    }

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
        IMPL::TrackImpl* trackcand = new IMPL::TrackImpl;           // Don't free it manually, because it is owned by trkCandCollection

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
    
    streamlog_out(MESSAGE4) << "EUTelProcessorTrackingHelixTrackSearch::end()  " << name()
            << " processed " << _nProcessedEvents << " events in " << _nProcessedRuns << " runs "
            << " av.tracks : " << static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_numberTracksCandidatesHistName ] ) -> mean()
            << " track candidates : " << static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_numberOfHitOnTrackCandidateHistName ] ) -> allEntries()
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

        const int    nhitsNBin = 40;
        const double nhitsMin = -0.5;
        const double nhitsMax =39.5;
        
        histoInfo = histoMgr->getHistogramInfo(_histName::_HitOnTrackCandidateHistName);        
        NBin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : nhitsNBin;
        XMin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : nhitsMin;
        XMax = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : nhitsMax;
        
        // Number of hit per track candidate
        AIDA::IHistogram1D * HitsOnTrackCandidates =
                marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_HitOnTrackCandidateHistName, NBin,
                XMin, XMax);

        if (HitsOnTrackCandidates) {
            HitsOnTrackCandidates->setTitle("hits on track candidates;N hits;N tracks");
            _aidaHistoMap1D.insert(make_pair(_histName::_HitOnTrackCandidateHistName, HitsOnTrackCandidates));
        } else {
            streamlog_out(ERROR2) << "Problem booking the " << (_histName::_HitOnTrackCandidateHistName) << endl;
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
