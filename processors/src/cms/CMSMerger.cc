/*
 * Created by Thomas Eichhorn
 *  (2017 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/ITree.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#endif

// lcio includes <.h>
#include <lcio.h>
#include <Exceptions.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCGenericObjectImpl.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IO/LCWriter.h>
#include <IO/LCReader.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTime.h>
#include <UTIL/LCTOOLS.h>

// system includes <>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <memory>
#include <glob.h>
#include <vector>
#include <set>
#include <map>

// eutelescope includes ""
#include "anyoption.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelVirtualCluster.h"
#include "EUTelSparseClusterImpl.h"

#include "CMSMerger.h"


using namespace std;
using namespace lcio;
using namespace marlin;
using namespace IMPL;
using namespace eutelescope;


AIDA::IHistogram1D * multiplicityhisto;
AIDA::IHistogram2D * mergecorrelation_x;
AIDA::IHistogram2D * mergecorrelation_y;


CMSMerger::CMSMerger ( ) : Processor ( "CMSMerger" )
{

    _description = "CMSMerger merges the CBC data stream with the telescope data stream, based on events or TLU time stamps.";

    registerProcessorParameter ( "CBCDataCollectionName1", "The name of the CBC data collection we want to merge - sensor 1", _cbcDataCollectionName1, string ( "cbc_collection1d" ) );

    registerProcessorParameter ( "CBCPulseCollectionName1", "The name of the CBC pulse collection we want to merge - sensor 1", _cbcPulseCollectionName1, string ( "cbc_collection1p" ) );

    registerProcessorParameter ( "CBCDataCollectionName2", "The name of the CBC data collection we want to merge - sensor 2", _cbcDataCollectionName2, string ( "cbc_collection2d" ) );

    registerProcessorParameter ( "CBCPulseCollectionName2", "The name of the CBC pulse collection we want to merge - sensor 2", _cbcPulseCollectionName2, string ( "cbc_collection2p" ) );

    registerProcessorParameter ( "CorrelationPlane", "The sensorID of the telescope plane to use for correlation plots", _correlationPlaneID, 2 );

    registerProcessorParameter ( "EventMerge", "Merge events based on event number (true) or on event time (false)", _eventmerge, true );

    registerProcessorParameter ( "OutputCollectionName", "The name of the output collection we want to create", _outputCollectionName, string ( "output_collection1" ) );

    registerProcessorParameter ( "OutputCollectionName2", "The name of the secondary output collection we want to create", _outputCollectionName2, string ( "output_collection2" ) );

    registerProcessorParameter ( "OutputCollectionName3", "The name of the tertiary output collection we want to create", _outputCollectionName3, string ( "output_collection3" ) );

    registerProcessorParameter ( "ReadCBCAhead", "Read this number of CBC events before reading the first telescope event", _eventdifferenceCBC, 0 );

    registerProcessorParameter ( "ReadTelescopeAhead", "Read this number of telescope events before reading the first CBC event", _eventdifferenceTelescope, 0 );

    registerProcessorParameter ( "TelescopeCollectionName", "The name of the telescope collection we want to merge", _telescopeCollectionName, string ( "telescope_collection1" ) );

    registerProcessorParameter ( "TelescopeCollectionName2", "The name of the secondary telescope collection we want to merge", _telescopeCollectionName2, string ( "telescope_collection2" ) );

    registerProcessorParameter ( "TelescopeFile", "The filename where the telescope data is stored", _telescopeFile, string ( "dummy_telescope.slcio" ) );

}


void CMSMerger::init ( )
{
    streamlog_out ( MESSAGE4 ) << "Running init" << endl;
    printParameters ( );

    // set times, init with telescope smaller, so that we read the first event
    _cbceventtime = -1;
    _telescopeeventtime = -2;
    _maxevents = 2;
    _readcount = 0;


    for ( int i = 0; i < _eventdifferenceTelescope; i++ )
    {
	LCEvent *evt = readTelescope ( );
	streamlog_out ( MESSAGE4 ) << "Skipped " << i + 1 << " telescope events!" << endl;
	streamlog_out ( MESSAGE4 ) << "Event skipped was " << evt -> getEventNumber ( ) << endl;
    }

}


void CMSMerger::processRunHeader ( LCRunHeader * rdr )
{
    streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;

    auto arunHeader = std::make_unique < EUTelRunHeaderImpl > ( rdr );
    arunHeader -> addProcessor ( type ( ) );

    // so we only open the telescope file once...
    _telescopeopen = false;

    bookHistos ( );

}


LCEvent *CMSMerger::readTelescope ( )
{
    // the telescope file is read here...
    if ( _telescopeopen == false )
    {
	lcReader = LCFactory::getInstance ( ) -> createLCReader ( IO::LCReader::directAccess );
	try
	{
	    lcReader -> open ( _telescopeFile );
	    _telescopeopen = true;
	    _maxevents = lcReader -> getNumberOfEvents ( );
	}
	catch ( IOException& e )
	{
	    streamlog_out ( ERROR1 ) << "Can't open the telescope file: " << e.what ( ) << endl;
	}
    }

    try
    {
	LCEvent *evt = lcReader -> readNextEvent ( );
	if ( evt == nullptr )
	{
	    return nullptr;
	    streamlog_out ( ERROR1 ) << "FAIL! NULL Event!" << endl ;
	}

	return ( evt );

    }
    catch ( IOException& e )
    {
	streamlog_out ( ERROR1 ) << "FAIL: " << e.what ( ) << endl ;
	return nullptr;
    }
}


void CMSMerger::processEvent ( LCEvent * anEvent )
{

    // the input collections
    LCCollectionVec * cbcDataCollectionVec1;
    LCCollectionVec * cbcPulseCollectionVec1;
    LCCollectionVec * cbcDataCollectionVec2;
    LCCollectionVec * cbcPulseCollectionVec2;

    LCCollectionVec * telescopeCollectionVec;
    LCCollectionVec * telescopeCollectionVec2;
    LCCollectionVec * telescopeCollectionVec3;

    // here we have to merge two collections: the trackerdata and the trackerpulse. We also keep the zsdata for reference -> 3 collections in total

    // the output collections
    LCCollectionVec * outputVec1 = new LCCollectionVec ( LCIO::TRACKERDATA );
    LCCollectionVec * outputVec2 = new LCCollectionVec ( LCIO::TRACKERPULSE );
    LCCollectionVec * outputVec3 = new LCCollectionVec ( LCIO::TRACKERDATA );

    // these encoders produce output - the encoding string should be the same as in the clustering and fit to the telescope
    CellIDEncoder < TrackerDataImpl > outputDataEncoder ( eutelescope::EUTELESCOPE::ZSCLUSTERDEFAULTENCODING, outputVec1 );
    CellIDEncoder < TrackerPulseImpl > outputPulseEncoder ( eutelescope::EUTELESCOPE::PULSEDEFAULTENCODING, outputVec2 );
    CellIDEncoder < TrackerDataImpl > outputDataEncoder2 ( "sensorID:7,sparsePixelType:5,quality:5", outputVec3 );
    CellIDEncoder < TrackerPulseImpl > outputPulseEncoder2 ( "sensorID:7,xSeed:12,ySeed:12,xCluSize:5,yCluSize:5,type:5,quality:5", outputVec2 );

    int cbcsize1d = 0;
    int cbcsize1p = 0;
    int cbcsize2d = 0;
    int cbcsize2p = 0;
    int telescopesize = 0;
    int telescopesize2 = 0;
    int telescopesize3 = 0;

    std::vector < double > cbc_corr;
    std::vector < double > tele_corr_x;
    std::vector < double > tele_corr_y;

    if ( _eventdifferenceCBC == 0 )
    {

	try
	{

	    _cbceventtime = anEvent -> getTimeStamp ( );
	    streamlog_out ( DEBUG4 ) << "CBC time is " << _cbceventtime << endl;

	    cbcDataCollectionVec1 = dynamic_cast < LCCollectionVec * > ( anEvent -> getCollection ( _cbcDataCollectionName1 ) );
	    cbcsize1d = cbcDataCollectionVec1 -> getNumberOfElements ( );
	    streamlog_out ( DEBUG1 ) << cbcsize1d << " Elements in cbc event!" << endl;

	    // and the secondary collections
	    cbcPulseCollectionVec1 = dynamic_cast < LCCollectionVec * > ( anEvent -> getCollection ( _cbcPulseCollectionName1 ) );
	    cbcsize1p = cbcPulseCollectionVec1 -> getNumberOfElements ( );
	    streamlog_out ( DEBUG1 ) << cbcsize1p << " Elements in cbc event - collection 2!" << endl;

	    CellIDDecoder < TrackerDataImpl > inputSparseColDecoder ( cbcDataCollectionVec1 );
	    CellIDDecoder < TrackerPulseImpl > inputPulseColDecoder ( cbcPulseCollectionVec1 );

	    // Cell ID Encoders for the telescope, introduced in eutelescope::EUTELESCOPE
	    // for sparseFrame (usually called cluster collection)
	    CellIDEncoder < TrackerDataImpl > outputSparseColEncoder ( eutelescope::EUTELESCOPE::ZSCLUSTERDEFAULTENCODING, outputVec1 );
	    // for pulseFrame
	    CellIDEncoder < TrackerPulseImpl > outputPulseColEncoder ( eutelescope::EUTELESCOPE::PULSEDEFAULTENCODING, outputVec2 );

	    unsigned int nCBCClusters;
	    // go through input clusters and copy them to output cluster collection

	    nCBCClusters = cbcPulseCollectionVec1 -> getNumberOfElements ( );
	    for ( size_t i = 0; i < nCBCClusters; ++i )
	    {
		TrackerPulseImpl * outputPulseFrame = new TrackerPulseImpl ( );
		TrackerDataImpl * outputSparseFrame = new TrackerDataImpl ( );

		TrackerPulseImpl* inputPulseFrame = dynamic_cast < TrackerPulseImpl* > ( cbcPulseCollectionVec1 -> getElementAt ( i ) );
		TrackerDataImpl* inputSparseFrame = dynamic_cast < TrackerDataImpl* > ( inputPulseFrame -> getTrackerData ( ) );

		// set Cell ID for sparse collection
		outputSparseColEncoder["sensorID"] = static_cast < int > ( inputSparseColDecoder ( inputSparseFrame ) ["sensorID"] );
		outputSparseColEncoder["sparsePixelType"] = static_cast < int > ( inputSparseColDecoder ( inputSparseFrame ) ["sparsePixelType"] );
		outputSparseColEncoder["quality"] = static_cast < int > ( inputSparseColDecoder ( inputSparseFrame ) ["quality"] );
		outputSparseColEncoder.setCellID ( outputSparseFrame );

		// copy tracker data
		outputSparseFrame -> setChargeValues ( inputSparseFrame -> getChargeValues ( ) );
		// add it to the cluster collection
		outputVec1 -> push_back ( outputSparseFrame );

		// prepare a pulse for this cluster
		outputPulseColEncoder["sensorID"] = static_cast < int > ( inputPulseColDecoder ( inputPulseFrame ) ["sensorID"] );
		outputPulseColEncoder["type"] = static_cast < int > ( inputPulseColDecoder ( inputPulseFrame ) ["type"] );
		outputPulseColEncoder.setCellID ( outputPulseFrame );

		outputPulseFrame -> setCharge ( inputPulseFrame -> getCharge ( ) );
		outputPulseFrame -> setTrackerData ( outputSparseFrame );
		outputVec2 -> push_back ( outputPulseFrame );

		cbc_corr.push_back ( inputPulseColDecoder ( inputPulseFrame ) ["xSeed"] );

	    }

	}
	catch ( lcio::DataNotAvailableException& )
	{
	    streamlog_out ( DEBUG4 ) << "Collection " << _cbcDataCollectionName1 << " not found in event " << anEvent -> getEventNumber ( ) << endl;
	}

	try
	{

	    cbcDataCollectionVec2 = dynamic_cast < LCCollectionVec * > ( anEvent -> getCollection ( _cbcDataCollectionName2 ) );
	    cbcsize2d = cbcDataCollectionVec2 -> getNumberOfElements ( );
	    streamlog_out ( DEBUG1 ) << cbcsize2d << " Elements in cbc event!" << endl;

	    // and the secondary collections
	    cbcPulseCollectionVec2 = dynamic_cast < LCCollectionVec * > ( anEvent -> getCollection ( _cbcPulseCollectionName2 ) );
	    cbcsize2p = cbcPulseCollectionVec2 -> getNumberOfElements ( );
	    streamlog_out ( DEBUG1 ) << cbcsize2p << " Elements in cbc event - collection 2!" << endl;

	    CellIDDecoder < TrackerDataImpl > inputSparseColDecoder ( cbcDataCollectionVec2 );
	    CellIDDecoder < TrackerPulseImpl > inputPulseColDecoder ( cbcPulseCollectionVec2 );

	    // Cell ID Encoders for the telescope, introduced in eutelescope::EUTELESCOPE
	    // for sparseFrame (usually called cluster collection)
	    CellIDEncoder < TrackerDataImpl > outputSparseColEncoder ( eutelescope::EUTELESCOPE::ZSCLUSTERDEFAULTENCODING, outputVec1 );
	    // for pulseFrame
	    CellIDEncoder < TrackerPulseImpl > outputPulseColEncoder ( eutelescope::EUTELESCOPE::PULSEDEFAULTENCODING, outputVec2 );

	    unsigned int nCBCClusters;
	    // go through input clusters and copy them to output cluster collection

	    nCBCClusters = cbcPulseCollectionVec2 -> getNumberOfElements ( );
	    for ( size_t i = 0; i < nCBCClusters; ++i )
	    {
		TrackerPulseImpl * outputPulseFrame = new TrackerPulseImpl ( );
		TrackerDataImpl * outputSparseFrame = new TrackerDataImpl ( );

		TrackerPulseImpl* inputPulseFrame = dynamic_cast < TrackerPulseImpl* > ( cbcPulseCollectionVec2 -> getElementAt ( i ) );
		TrackerDataImpl* inputSparseFrame = dynamic_cast < TrackerDataImpl* > ( inputPulseFrame -> getTrackerData ( ) );

		// set Cell ID for sparse collection
		outputSparseColEncoder["sensorID"] = static_cast < int > ( inputSparseColDecoder ( inputSparseFrame ) ["sensorID"] );
		outputSparseColEncoder["sparsePixelType"] = static_cast < int > ( inputSparseColDecoder ( inputSparseFrame ) ["sparsePixelType"] );
		outputSparseColEncoder["quality"] = static_cast < int > ( inputSparseColDecoder ( inputSparseFrame ) ["quality"] );
		outputSparseColEncoder.setCellID ( outputSparseFrame );

		// copy tracker data
		outputSparseFrame -> setChargeValues ( inputSparseFrame -> getChargeValues ( ) );
		// add it to the cluster collection
		outputVec1 -> push_back ( outputSparseFrame );

		// prepare a pulse for this cluster
		outputPulseColEncoder["sensorID"] = static_cast < int > ( inputPulseColDecoder ( inputPulseFrame ) ["sensorID"] );
		outputPulseColEncoder["type"] = static_cast < int > ( inputPulseColDecoder ( inputPulseFrame ) ["type"] );
		outputPulseColEncoder.setCellID ( outputPulseFrame );

		outputPulseFrame -> setCharge ( inputPulseFrame -> getCharge ( ) );
		outputPulseFrame -> setTrackerData ( outputSparseFrame );
		outputVec2 -> push_back ( outputPulseFrame );

	    }

	}
	catch ( lcio::DataNotAvailableException& )
	{
	    streamlog_out ( DEBUG4 ) << "Collection " << _cbcDataCollectionName2 << " not found in event " << anEvent -> getEventNumber ( ) << endl;
	}

    }

    if ( _eventdifferenceCBC > 0 )
    {
	_eventdifferenceCBC--;
	streamlog_out ( MESSAGE4 ) << "Skipping a CBC event!" << endl;
    }

    try
    {

	// process this guy...
	LCEvent* evt;

	if ( _eventmerge == false )
	{

	    // the cbc has moved on in time, we need to read again
	    if ( _telescopeeventtime < _cbceventtime )
	    {

		multiplicityhisto -> fill ( _multiplicity );
		_multiplicity = 1;

		// the telescope is read by the function
		if ( _readcount <  _maxevents )
		{
		    evt = readTelescope ( );
		}
		else
		{
		    streamlog_out ( MESSAGE4 ) << "Reached EOF!" << endl;
		    exit ( -1 );
		}

		int tempnr = evt -> getEventNumber ( );
		streamlog_out ( DEBUG4 ) << "Reading new event nr " << tempnr << endl;
		// save for future...
		_storeevt = evt;

	    }
	    else
	    {
		_multiplicity++;
		evt = _storeevt;
		int tempnr = evt -> getEventNumber ( );

		streamlog_out ( DEBUG4 ) << "Reading old event nr " << tempnr << endl;
	    }

	}

	if ( _eventmerge == true )
	{

	    evt = readTelescope ( );

	}

	 _telescopeeventtime = evt -> getTimeStamp ( );
	streamlog_out ( DEBUG4 ) << "Telescope time is " << _telescopeeventtime << endl;

	telescopeCollectionVec = dynamic_cast < LCCollectionVec * > ( evt -> getCollection ( _telescopeCollectionName ) );
	telescopesize = telescopeCollectionVec -> getNumberOfElements ( );
	streamlog_out ( DEBUG1 ) << telescopesize << " Elements in telescope event!" << endl;

	// and the secondary collections
	telescopeCollectionVec2 = dynamic_cast < LCCollectionVec * > ( evt -> getCollection ( _telescopeCollectionName2 ) );
	telescopesize2 = telescopeCollectionVec2 -> getNumberOfElements ( );
	streamlog_out ( DEBUG1 ) << telescopesize2 << " Elements in telescope event - collection 2!" << endl;

	// this guy can have less events and is not merged, but copied
	telescopeCollectionVec3 = dynamic_cast < LCCollectionVec * > ( evt -> getCollection ( _outputCollectionName3 ) );
	telescopesize3 = telescopeCollectionVec3 -> getNumberOfElements ( );
	streamlog_out ( DEBUG1 ) << telescopesize3 << " Elements in Telescope event - collection 3!" << endl;

	CellIDDecoder < TrackerDataImpl > inputSparseColDecoder ( telescopeCollectionVec );
	CellIDDecoder < TrackerPulseImpl > inputPulseColDecoder ( telescopeCollectionVec2 );

	// Cell ID Encoders for the telescope, introduced in eutelescope::EUTELESCOPE
	// for sparseFrame (usually called cluster collection)
	CellIDEncoder < TrackerDataImpl > outputSparseColEncoder ( eutelescope::EUTELESCOPE::ZSCLUSTERDEFAULTENCODING, outputVec1 );
	// for pulseFrame
	CellIDEncoder < TrackerPulseImpl > outputPulseColEncoder ( eutelescope::EUTELESCOPE::PULSEDEFAULTENCODING, outputVec2 );

	unsigned int nTelClusters;
	// go through input clusters and copy them to output cluster collection

	nTelClusters = telescopeCollectionVec2 -> getNumberOfElements ( );
	for ( size_t i = 0; i < nTelClusters; ++i )
	{
	    TrackerPulseImpl * outputPulseFrame = new TrackerPulseImpl ( );
	    TrackerDataImpl * outputSparseFrame = new TrackerDataImpl ( );

	    TrackerPulseImpl* inputPulseFrame = dynamic_cast < TrackerPulseImpl* > ( telescopeCollectionVec2 -> getElementAt ( i ) );
	    TrackerDataImpl* inputSparseFrame = dynamic_cast < TrackerDataImpl* > ( inputPulseFrame -> getTrackerData ( ) );

	    // set Cell ID for sparse collection
	    outputSparseColEncoder["sensorID"] = static_cast < int > ( inputSparseColDecoder ( inputSparseFrame ) ["sensorID"] );
	    outputSparseColEncoder["sparsePixelType"] = static_cast < int > ( inputSparseColDecoder ( inputSparseFrame ) ["sparsePixelType"] );
	    outputSparseColEncoder["quality"] = static_cast < int > ( inputSparseColDecoder ( inputSparseFrame ) ["quality"] );
	    outputSparseColEncoder.setCellID ( outputSparseFrame );

	    // copy tracker data
	    outputSparseFrame -> setChargeValues ( inputSparseFrame -> getChargeValues ( ) );
	    // add it to the cluster collection
	    outputVec1 -> push_back ( outputSparseFrame );

	    // prepare a pulse for this cluster
	    outputPulseColEncoder["sensorID"] = static_cast < int > ( inputPulseColDecoder ( inputPulseFrame ) ["sensorID"] );
	    outputPulseColEncoder["type"] = static_cast < int > ( inputPulseColDecoder ( inputPulseFrame ) ["type"] );
	    outputPulseColEncoder.setCellID ( outputPulseFrame );

	    outputPulseFrame -> setCharge ( inputPulseFrame -> getCharge ( ) );
	    outputPulseFrame -> setTrackerData ( outputSparseFrame );
	    outputVec2 -> push_back ( outputPulseFrame );

	    int tempplane = -1;
	    tempplane = static_cast < int > ( inputPulseColDecoder ( inputPulseFrame ) ["sensorID"] );

	    if ( tempplane == _correlationPlaneID )
	    {
		// get these for the correlation plots
		EUTelVirtualCluster * tempCluster;
		tempCluster = new EUTelSparseClusterImpl < EUTelGenericSparsePixel > ( static_cast < TrackerDataImpl* > ( inputPulseFrame -> getTrackerData ( ) ) );
		float tempx = 0.0;
		float tempy = 0.0;
		tempCluster -> getCenterOfGravity ( tempx, tempy ) ;
		streamlog_out ( DEBUG0 ) << "Plot x is " << tempx << endl;
		streamlog_out ( DEBUG0 ) << "Plot y is " << tempy << endl;
		delete tempCluster;
		tele_corr_x.push_back ( tempx );
		tele_corr_y.push_back ( tempy );
	    }

	} // end of loop over input clusters

	// one more loop over the telescope, now the third collection
	for ( int j = 0; j < telescopesize3; j++ )
	{
	    streamlog_out ( DEBUG1 ) << "Reading telescope again..." << endl;
	    lcio::TrackerDataImpl * input  = dynamic_cast < lcio::TrackerDataImpl * > ( telescopeCollectionVec3 -> getElementAt ( j ) );
	    lcio::TrackerDataImpl * output = new lcio::TrackerDataImpl;
	    output -> setChargeValues ( input -> getChargeValues ( ) );
	    output -> setCellID0 ( input -> getCellID0 ( ) );
	    output -> setCellID1 ( input -> getCellID1 ( ) );
	    output -> setTime ( input -> getTime ( ) );
	    outputVec3 -> addElement ( output );
	    streamlog_out ( DEBUG1 ) << "Wrote telescope again..." << endl;
	}

    }
    catch ( lcio::DataNotAvailableException& )
    {
	streamlog_out( DEBUG4 ) << "Collection " << _telescopeCollectionName << " not found in event " << anEvent -> getEventNumber ( ) << endl;
    }

    streamlog_out ( DEBUG1 ) << "Writing out event " << anEvent -> getEventNumber ( ) << endl;
    anEvent -> addCollection ( outputVec1, _outputCollectionName );
    anEvent -> addCollection ( outputVec2, _outputCollectionName2 );
    anEvent -> addCollection ( outputVec3, _outputCollectionName3 );

}


void CMSMerger::check (LCEvent * /* evt */ )
{

}


void CMSMerger::end ( )
{
    // the telescope file is still open, we can now close it
    lcReader -> close ( );
    delete lcReader;
    streamlog_out ( MESSAGE4 ) << "Successfully finished!" << endl;
}


void CMSMerger::fillHistos ( )
{

}


void CMSMerger::bookHistos ( )
{
    multiplicityhisto = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( "Multiplicity", 10, -0.5, 9.5 );
    multiplicityhisto -> setTitle ( "Telescope Events per CBC Event;Telescope Events;Events" );

    mergecorrelation_x = AIDAProcessor::histogramFactory ( this ) -> createHistogram2D ( "Correlation in X", 1152, 0, 1151, 1016, 0, 1015 );
    mergecorrelation_x -> setTitle ( "Telescope Cluster;CBC Cluster;Entries" );

    mergecorrelation_y = AIDAProcessor::histogramFactory ( this ) -> createHistogram2D ( "Correlation in Y", 576, 0, 575, 1016, 0, 1015 );
    mergecorrelation_y -> setTitle ( "Telescope Cluster;CBC Cluster;Entries" );

}
