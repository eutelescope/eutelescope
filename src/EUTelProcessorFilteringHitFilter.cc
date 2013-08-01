#include "EUTelProcessorFilteringHitFilter.h"

// C++
#include <vector>

// LCIO
#include <EVENT/LCCollection.h>

// MARLIN
#include "marlin/VerbosityLevels.h"

// AIDA
#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#endif // MARLIN_USE_AIDA

// EUTELESCOPE
#include "EUTelExhaustiveTrackFinder.h"

using namespace lcio;
using namespace marlin;
using namespace eutelescope;

/**  Pattern recognition processor
 * 
 *  Processor selects hits that fulfill all specified requirements from input collection.
 * 
 *  If compiled with MARLIN_USE_AIDA 
 *  it creates histograms.
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Requires the collection of hits.
 *
 *  <h4>Output</h4> 
 *  A collection of hits classified according to presumable track candidates.
 *  Histograms.
 * 
 *  @param CollectionName Name of the --insert here-- collection
 * 
 */


EUTelProcessorFilteringHitFilter aEUTelProcessorFilteringHitFilter;

//! Default constructor
EUTelProcessorFilteringHitFilter::EUTelProcessorFilteringHitFilter( ) :
 Processor( "EUTelProcessorFilteringHitFilter" ),
 _hitInputCollectionName( "HitCollection" ),
 _nProcessedRuns( 0 ),
 _nProcessedEvents( 0 ) {

    // Processor description
    _description = "EUTelProcessorFilteringHitFilter selects hits that fulfill all specified requirements from input collection.";
    
    //--------------------------------------------------------------------------------------------------------------------------------------------------------//
    
    // Hits input collection
    registerInputCollection( LCIO::TRACKERHIT, "HitInputCollectionName", "Input hits collection name", _hitInputCollectionName, "HitCollection" );
    
    // Filtered hits output collection
    registerOutputCollection( LCIO::TRACKERHIT, "HitOutputCollectionName", "Output hits collection name", _hitOutputCollectionName, "FilteredHitCollection" );
    
    //--------------------------------------------------------------------------------------------------------------------------------------------------------//
    
    // Optional parameters definition
    
    const IntVec defaultPlaneIDs( 0 );
    registerOptionalParameter ( "WantPlaneID", "Select hits with given plane IDs (all - if empty)", _wantPlaneID, defaultPlaneIDs );
}

void EUTelProcessorFilteringHitFilter::init( ) {

    streamlog_out( DEBUG ) << "EUTelProcessorFilteringHitFilter::init( )" << std::endl;

    // usually a good idea to
    printParameters( );

    _nProcessedRuns = 0;
    _nProcessedEvents = 0;

}

void EUTelProcessorFilteringHitFilter::processRunHeader( LCRunHeader* /*run*/ ) {

    _nProcessedRuns++;
}

void EUTelProcessorFilteringHitFilter::processEvent( LCEvent * evt ) {

    // Try to access the collection
    
    LCCollection* col = NULL;
    try{
        col = evt->getCollection( _hitInputCollectionName );
    }
    catch( lcio::DataNotAvailableException e )
    {
        streamlog_out( WARNING ) << _hitInputCollectionName << " collection not available" << std::endl;
        col = NULL;
    }

    // this will only be entered if the collection is available
    if ( col != NULL ) {

    }

    

    _nProcessedEvents++;
}

void EUTelProcessorFilteringHitFilter::check( LCEvent * /*evt*/ ) {
    // nothing to check here
}

void EUTelProcessorFilteringHitFilter::end( ) {

       streamlog_out( MESSAGE )  << "EUTelProcessorTrackingExhaustiveTrackSearch::end()  " << name() 
                                 << " processed " << _nProcessedEvents << " events in " << _nProcessedRuns << " runs "
                                 << std::endl;

}

