/*
 * Rewritten by Thomas Eichhorn
 *  (2017 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */


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
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

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

using namespace lcio;
using namespace marlin;
using namespace eutelescope;

AIDA::IHistogram1D * trackHist;
AIDA::IHistogram1D * hitTrackHist;

EUTelProcessorTrackingHelixTrackSearch::EUTelProcessorTrackingHelixTrackSearch ( ) : Processor ( "EUTelProcessorTrackingHelixTrackSearch" )
{

    _description = "EUTelProcessorTrackingHelixTrackSearch preforms track pattern recognition in magnetic fields.";

    registerInputCollection ( LCIO::TRACKERHIT, "HitInputCollectionName", "Input hits collection name", _hitInputCollectionName, std::string ( "InputHitCollection" ) );

    registerOutputCollection ( LCIO::TRACK, "TrackCandHitOutputCollectionName", "Output track collection name", _trackCandidateHitsOutputCollectionName, std::string ( "TrackCandidateHitCollection" ) );

    registerOptionalParameter ( "MaxMissingHitsPerTrack", "Maximal number of missing hits on a track candidate", _maxMissingHitsPerTrackCand, static_cast < int > ( 0 ) );

    registerOptionalParameter ( "ResidualsRMax", "Maximal allowed distance in mm between hits in the recognition step", _residualsRMax, static_cast < double > ( 0.4 ) );

    registerOptionalParameter ( "BeamEnergy", "Beam energy in GeV", _eBeam, static_cast < double > ( 4.0 ) );

    registerOptionalParameter ( "BeamCharge", "Beam charge in e", _qBeam, static_cast < double > ( -1 ) );

    registerOptionalParameter ( "BeamSpread", "Angular spread of the beam (horizontal,vertical) in mrad for beam constraint. No beam constraints if negative values are supplied.", _beamSpread, EVENT::FloatVec ( 2.0, 2.0 ) );

    registerOptionalParameter ( "BeamEnergyUncertainty", "Uncertainty of the beam energy in %", _eBeamUncertatinty, static_cast < double > ( 0.0 ) );

}


void EUTelProcessorTrackingHelixTrackSearch::init ( )
{

    streamlog_out ( MESSAGE4 ) << "Running init" << std::endl;

    // usually a good idea to
    printParameters ( );

    // reset counters
    _nProcessedRuns = 0;
    _nProcessedEvents = 0;

    // get access to geometry description
    std::string name ( "test.root" );
    geo::gGeometry ( ) .initializeTGeoDescription ( name, false );

    // instantiate track finder
    {
	streamlog_out ( DEBUG0 ) << "Initialisation of track finder" << std::endl;

	EUTelKalmanFilter* Finder = new EUTelKalmanFilter ( "KalmanTrackFinder" );

	if ( !Finder )
	{
	    streamlog_out ( ERROR5 ) << "Can't allocate an instance of KalmanTrackFinder! Stopping ..." << std::endl;
            throw UnknownDataTypeException ( "Track finder was not created" );
	}

	Finder -> setAllowedMissingHits ( _maxMissingHitsPerTrackCand );
	Finder -> setWindowSize ( _residualsRMax );
	Finder -> setBeamMomentum ( _eBeam );
	Finder -> setBeamCharge ( _qBeam );
	Finder -> setBeamMomentumUncertainty ( _eBeamUncertatinty );
	Finder -> setBeamSpread ( _beamSpread );

	_trackFitter = Finder;
    }

    // book histograms
    bookHistograms ( );

}


void EUTelProcessorTrackingHelixTrackSearch::processRunHeader ( LCRunHeader * run )
{
    streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;

    auto arunHeader = std::make_unique < EUTelRunHeaderImpl > ( run );
    arunHeader -> addProcessor ( type ( ) );

    _nProcessedRuns++;

}


void EUTelProcessorTrackingHelixTrackSearch::processEvent ( LCEvent * evt )
{

    EUTelEventImpl * event = static_cast < EUTelEventImpl* > ( evt );

    // do not process last or unknown events
    if ( event -> getEventType ( ) == kEORE )
    {
	streamlog_out ( MESSAGE4 ) << "EORE found: nothing else to do." << std::endl;
	return;
    }
    else if ( event -> getEventType ( ) == kUNKNOWN )
    {
	streamlog_out ( WARNING2 ) << "Event number " << event -> getEventNumber ( ) << " in run " << event -> getRunNumber ( ) << " is of unknown type. Continue considering it as a normal data event." << std::endl;
    }

    // Try to access collection
    LCCollection* col = NULL;
    try
    {
	col = evt -> getCollection ( _hitInputCollectionName );
    }
    catch ( DataNotAvailableException e )
    {
	streamlog_out ( WARNING2 ) << "Collection " << _hitInputCollectionName << " not available" << std::endl;
	throw marlin::SkipEventException ( this );
    }

    // if collection is available
    if ( col != NULL )
    {
	streamlog_out ( DEBUG0 ) << "Running processEvent in event " << event -> getEventNumber ( ) << std::endl;

	// Prepare hits for track finder
	EVENT::TrackerHitVec allHitsVec;
	FillHits ( col, allHitsVec );

	streamlog_out ( DEBUG0 ) << "All hits in this event:" << std::endl;
	EVENT::TrackerHitVec::const_iterator itHits;

	for ( itHits = allHitsVec.begin ( ); itHits != allHitsVec.end ( ); ++itHits )
	{
	    const double* uvpos = ( *itHits ) -> getPosition ( );
	    streamlog_out ( DEBUG0 ) << "Hit with sensor id = " << Utility::GuessSensorID ( *itHits ) << ", local(u,v) coordinates: (" << uvpos[0] << "," << uvpos[1] << ")" << std::endl;
	}

	// search tracks
	streamlog_out ( DEBUG0 ) << "Initialising hits..." << std::endl;
	static_cast < EUTelKalmanFilter* > ( _trackFitter ) -> setHits ( allHitsVec );

	bool isReady = static_cast < EUTelKalmanFilter* > ( _trackFitter ) -> initialise ( );

	if ( isReady )
	{
	    streamlog_out ( DEBUG0 ) << "Trying to find tracks..." << std::endl;
	    LCCollectionVec * myoutvec;
	    myoutvec = new LCCollectionVec ( LCIO::TRACKERHIT );

	    LCCollectionVec * tempvec = new LCCollectionVec ( LCIO::TRACKERHIT );

	    _trackFitter -> FitTracks ( );

	    streamlog_out ( DEBUG0 ) << "Retrieving track candidates..." << std::endl;
	    vector < IMPL::TrackImpl* > & trackCandidates = static_cast < EUTelKalmanFilter* > ( _trackFitter ) -> getTracks ( );

	    const int nTracks = trackCandidates.size ( );
	    trackHist -> fill ( nTracks );
	    streamlog_out( DEBUG0 ) << "Found " << nTracks << " in event " << event -> getEventNumber ( ) << std::endl;

	    int nHitsOnTrack = 0;
	    vector < IMPL::TrackImpl* > ::const_iterator itrk;
	    for ( itrk = trackCandidates.begin ( ); itrk != trackCandidates.end ( ); ++itrk )
	    {
		const EVENT::TrackerHitVec& trkHits = ( *itrk ) -> getTrackerHits ( );
		nHitsOnTrack = trkHits.size ( );

		EVENT::TrackerHitVec::const_iterator itTrkHits;
		streamlog_out ( DEBUG0 ) << "Track hits:" << std::endl;
		for ( itTrkHits = trkHits.begin ( ); itTrkHits != trkHits.end ( ); ++itTrkHits )
		{
		    const double* uvpos = ( *itTrkHits ) -> getPosition ( );
		    streamlog_out ( DEBUG0 ) << "Hit with sensor id = " << ( *itrk )->id( ) << ", local(u,v) coordinates: (" << uvpos[0] << "," << uvpos[1] << "," << uvpos[3] << ")" << std::endl;
		}

		myoutvec = static_cast < EUTelKalmanFilter* > ( _trackFitter ) -> getFitHits ( );

		for ( size_t ihit = 0; ihit < myoutvec -> size ( ); ihit++ )
		{
		    CellIDEncoder < TrackerHitImpl > fitHitEncoder ( EUTELESCOPE::HITENCODING, tempvec );
		    CellIDDecoder < TrackerHitImpl > inputCellIDDecoder ( EUTELESCOPE::HITENCODING );
		    TrackerHitImpl * inputHit = dynamic_cast < TrackerHitImpl* > ( myoutvec -> getElementAt ( ihit ) );
		    TrackerHitImpl * newHit = new TrackerHitImpl;

		    const int sensorID = inputCellIDDecoder ( inputHit ) ["sensorID"];
		    streamlog_out ( DEBUG0 ) << "id " << sensorID << std::endl;
		    fitHitEncoder["sensorID"] =  sensorID * 2;
		    fitHitEncoder["properties"] = 99;
		    //fitHitEncoder.setCellID ( newHit );

		    // copy hit position
		    const double* hitPos = inputHit -> getPosition ( );
		    newHit -> setPosition ( &hitPos[0] );

		    // copy cov. matrix
		    newHit -> setCovMatrix ( inputHit -> getCovMatrix ( ) );

		    // copy type
		    newHit -> setType ( inputHit -> getType ( ) );

		    // copy rawhits
		    LCObjectVec clusterVec = inputHit -> getRawHits ( );
		    newHit -> rawHits ( ) = clusterVec;

		    // copy cell IDs
		    newHit -> setCellID0 ( inputHit -> getCellID0 ( ) );
		    newHit -> setCellID1 ( inputHit -> getCellID1 ( ) );

		    // copy EDep
		    newHit -> setEDep ( inputHit -> getEDep ( ) );

		    // copy EDepError
		    newHit -> setEDepError ( inputHit -> getEDepError ( ) );

		    // copy Time
		    newHit -> setTime ( inputHit -> getTime ( ) );

		    // copy Quality
		    newHit -> setQuality ( inputHit -> getQuality ( ) );

		    newHit -> setCovMatrix( inputHit -> getCovMatrix ( ) );

		    tempvec -> push_back ( newHit );

		    ( *itrk ) -> addHit ( newHit );
		}

		hitTrackHist -> fill ( nHitsOnTrack );
	    }

	    // Write output collection
	    if ( nTracks > 0 )
	    {
		addTrackCandidateToCollection ( evt, trackCandidates );
		//evt -> addCollection ( tempvec, "testing2" );

	    }

	}
	_nProcessedEvents++;

	if ( isFirstEvent ( ) )
	{
	    _isFirstEvent = false;
	}

    } // if (col != NULL)
}


void EUTelProcessorTrackingHelixTrackSearch::FillHits ( LCCollection* collection, EVENT::TrackerHitVec& allHitsVec ) const
{

    allHitsVec.clear ( );
    allHitsVec.resize ( 0 );

    // loop over all hits in collection
    for ( int iHit = 0; iHit < collection -> getNumberOfElements ( ); iHit++ )
    {

	TrackerHitImpl * hit = static_cast < TrackerHitImpl* > ( collection -> getElementAt ( iHit ) );
	const int localSensorID = Utility::GuessSensorID ( hit );
        if ( localSensorID >= 0 )
	{
	    allHitsVec.push_back ( hit );
	}

    }
}


void EUTelProcessorTrackingHelixTrackSearch::end ( )
{
    delete _trackFitter;
    streamlog_out ( DEBUG0 ) << "Processed " << _nProcessedEvents << " events in " << _nProcessedRuns << " runs " << std::endl;
}


void EUTelProcessorTrackingHelixTrackSearch::check ( LCEvent * /*evt*/ )
{
    // nothing to check here
}


void EUTelProcessorTrackingHelixTrackSearch::addTrackCandidateToCollection ( LCEvent* evt, std::vector < IMPL::TrackImpl* > & trackCandidates )
{
    // prepare output collection
    LCCollectionVec * trkCandCollection = 0;
    try
    {
	trkCandCollection = new LCCollectionVec ( LCIO::TRACK );
	LCFlagImpl flag ( trkCandCollection -> getFlag ( ) );
	flag.setBit ( LCIO::TRBIT_HITS );
	trkCandCollection -> setFlag ( flag.getFlag ( ) );
    }
    catch ( ... )
    {
	streamlog_out ( WARNING2 ) << "Can't allocate output collection!" << endl;
    }

    // fill track parameters
    std::vector < IMPL::TrackImpl* > ::iterator itrk;
    for ( itrk = trackCandidates.begin ( ); itrk != trackCandidates.end ( ); ++itrk )
    {
	streamlog_out ( DEBUG1 ) << "Track has " << ( *itrk ) -> getTrackerHits ( ) .size ( ) << " hits" << endl;
	trkCandCollection -> push_back ( ( *itrk ) );

    }

    // write collection
    try
    {
	streamlog_out ( DEBUG0 ) << "Getting collection " << _trackCandidateHitsOutputCollectionName << endl;
	evt -> getCollection ( _trackCandidateHitsOutputCollectionName );
    }
    catch ( ... )
    {
	streamlog_out ( DEBUG0 ) << "Adding collection " << _trackCandidateHitsOutputCollectionName << endl;
	evt -> addCollection ( trkCandCollection, _trackCandidateHitsOutputCollectionName );
    }

}


void EUTelProcessorTrackingHelixTrackSearch::bookHistograms ( )
{
    #if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

    trackHist = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( "trackHist", 20, -0.5, 19.5 );
    trackHist -> setTitle ( "Track Candidates per Event;Number of Tracks;Entries" );

    hitTrackHist = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( "hitTrackHist", 50, -0.5, 49.5 );
    hitTrackHist -> setTitle ( "Hits per Track Candidates;Number of Hits per Track;Entries" );

    #endif
}
