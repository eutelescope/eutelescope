/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

// alibava includes ".h"
#include "AlibavaMerger.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"
#include "AlibavaPedNoiCalIOManager.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <lcio.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IO/LCWriter.h>
#include <IO/LCReader.h>
#include <Exceptions.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/LCGenericObjectImpl.h>

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
#include "EUTelVirtualCluster.h"
#include "anyoption.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelSparseClusterImpl.h"

// ROOT includes ".h"
#include "TH3D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TSystem.h"

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;
using namespace IMPL;
using namespace eutelescope;

AlibavaMerger::AlibavaMerger ( ) : AlibavaBaseProcessor ( "AlibavaMerger" )
{

    _description = "AlibavaMerger merges the Alibava cluster data stream with the telescope data stream.";

    registerProcessorParameter ( "AlibavaCollectionName", "The name of the alibava collection we want to merge", _alibavaCollectionName, string ( "original_zsdata" ) );

    registerProcessorParameter ( "AlibavaCollectionName2", "The name of the secondary alibava collection we want to merge", _alibavaCollectionName2, string ( "clusters" ) );

    registerProcessorParameter ( "AlibavaFile", "The filename where the alibava data is stored", _alibavaFile, string ( "alibava.slcio" ) );

    registerProcessorParameter ( "EventdifferenceAlibava", "The event count the Alibava is behind (read: earlier than) the Telescope. 1 means alibava event 1 == telescope event 0, etc.", _eventdifferenceAlibava, 0 );

    registerProcessorParameter ( "EventdifferenceTelescope", "The event count the telescope is behind (read: earlier than) the Alibava. 1 means alibava event 0 == telescope event 1, etc.", _eventdifferenceTelescope, 0 );

    registerProcessorParameter ( "OutputCollectionName", "The name of the output collection we want to create", _outputCollectionName, string ( "original_zsdata" ) );

    registerProcessorParameter ( "OutputCollectionName2", "The name of the secondary output collection we want to create", _outputCollectionName2, string ( "combinedcluster" ) );

    registerProcessorParameter ( "OutputCollectionName3", "The name of the tertiary output collection we want to create", _outputCollectionName3, string ( "zsdata_m26" ) );

    registerProcessorParameter ( "OutputMode", "The verbosity of the merged hits: 0 to only write events where we have a hit in both systems, 1 for events where we at least have one alibava hit, 2 for events where we at least have one telescope hit and 3 for writing out all events, even if there are no hits in them", _outputmode,  0 );

    registerProcessorParameter ( "PlaneShiftTelescope", "Option to decrement/increment the telescope plane sensor id", _teleplaneshift, 0 );

    registerProcessorParameter ( "TelescopeCollectionName", "The name of the telescope collection we want to merge", _telescopeCollectionName, string ( "original_zsdata" ) );

    registerProcessorParameter ( "TelescopeCollectionName2", "The name of the secondary telescope collection we want to merge", _telescopeCollectionName2, string ( "cluster_m26" ) );

    registerProcessorParameter ( "TelescopeFile", "The filename where the telescope data is stored", _telescopeFile , string ( "telescope.slcio" ) );

    registerProcessorParameter ( "UnsensitiveAxis", "The unsensitive axis of our strip sensor", _nonsensitiveaxis, string ( "x" ) );

}


void AlibavaMerger::init ( )
{
    streamlog_out ( MESSAGE4 ) << "Running init" << endl;

    // this method is called only once even when the rewind is active
    // usually a good idea to
    printParameters ( );
}


void AlibavaMerger::processRunHeader ( LCRunHeader * rdr )
{
    streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;

    auto arunHeader = std::make_unique < AlibavaRunHeaderImpl > ( rdr );
    arunHeader -> addProcessor ( type ( ) );

    // so we only open the telescope file once...
    _telescopeopen = false;

    bookHistos ( );

    for ( int i = 0; i < _eventdifferenceTelescope; i++ )
    {
	LCEvent *evt = readTelescope ( );
	streamlog_out ( MESSAGE4 ) << "Skipped " << i + 1 << " telescope events!" << endl;
	streamlog_out ( MESSAGE4 ) << "Event skipped was " << evt -> getEventNumber ( ) << endl;
    }
}


// the telescope file is read here:
LCEvent *AlibavaMerger::readTelescope ( )
{
    if ( _telescopeopen == false )
    {
	lcReader = LCFactory::getInstance ( ) -> createLCReader ( IO::LCReader::directAccess ) ;
	try
	{
	    lcReader -> open ( _telescopeFile ) ;
	    _telescopeopen = true;
	}
	catch ( IOException& e )
	{
	    streamlog_out ( ERROR1 ) << "Can't open the telescope file: " << e.what ( ) << endl ;
	}
    }

    try
    {
	LCEvent *evt = lcReader -> readNextEvent ( );
	if ( evt == nullptr )
	{
	    return nullptr;
	    streamlog_out ( ERROR1 ) << "FAIL: " << endl ;
	}
	return ( evt );
    }
    catch ( IOException& e )
    {
	    streamlog_out ( ERROR1 ) << "FAIL: " << e.what ( ) << endl ;
	    return nullptr;
    }
}

void AlibavaMerger::addCorrelation ( float ali_x, float ali_y, float ali_z, float tele_x, float tele_y, float tele_z, int event )
{
    if ( TH2D * signalHisto = dynamic_cast < TH2D* > ( _rootObjectMap["Correlation_X"] ) )
    {
	signalHisto -> Fill ( ali_x, tele_x );
    }
    if ( TH2D * signalHisto2 = dynamic_cast < TH2D* > ( _rootObjectMap["Correlation_Y"] ) )
    {
	signalHisto2 -> Fill ( ali_y, tele_y );
    }
    if ( TH2D * signalHisto3 = dynamic_cast < TH2D* > ( _rootObjectMap["Correlation_Z"] ) )
    {
	signalHisto3 -> Fill ( ali_z, tele_z );
    }
    if ( TH2D * signalHisto4 = dynamic_cast < TH2D* > ( _rootObjectMap["Correlation_Event"] ) )
    {
	if ( _nonsensitiveaxis == "x" )
	{
	    signalHisto4 -> Fill ( event, ( tele_y / 576.0 - ali_y / 256.0 ) );
	    streamlog_out ( DEBUG4 ) << "Correlation distance in x, event : " << event << " " << tele_y / 576.0 - ali_y / 256.0 << endl;
	}
	if ( _nonsensitiveaxis == "y" )
	{
	    signalHisto4 -> Fill ( event, ( tele_x / 1152.0 - ali_x / 256.0 ) );
	    streamlog_out ( DEBUG4 ) << "Correlation distance in y, event : " << event << " " << tele_x / 1152.0 - ali_x / 256.0 << endl;
	}
    }
    if ( TH3D * multihisto = dynamic_cast < TH3D* > ( _rootObjectMap["Correlation_All"] ) )
    {
	multihisto -> Fill ( ( tele_x / 1152.0 - ali_x / 1152.0 ), ( tele_y / 576.0 - ali_y / 256.0 ), event );
    }
}

void AlibavaMerger::processEvent ( LCEvent * anEvent )
{

    // the alibava event we are processing
    AlibavaEventImpl * alibavaEvent = static_cast < AlibavaEventImpl* > ( anEvent );

    // the input collections
    LCCollectionVec * alibavaCollectionVec;
    LCCollectionVec * telescopeCollectionVec;

    // the correlation output
    float alibava_x = 0.0;
    float alibava_y = 0.0;
    float alibava_z = 0.0;
    float telescope_x = 0.0;
    float telescope_y = 0.0;
    float telescope_z = 0.0;

    std::vector < double > ali_corr;
    std::vector < double > tele_corr;

    // here we have to merge two collections: the trackerdata and the trackerpulse. We also keep the zsdata for reference -> 3 collections in total

    // out output collections
    LCCollectionVec * outputTrackerDataVec = new LCCollectionVec ( LCIO::TRACKERDATA );
    LCCollectionVec * outputTrackerPulseVec = new LCCollectionVec ( LCIO::TRACKERPULSE );
    LCCollectionVec * outputTrackerDataVec2 = new LCCollectionVec ( LCIO::TRACKERDATA );

    // the input collections
    LCCollectionVec * alibavaCollectionVec2;
    LCCollectionVec * telescopeCollectionVec2;
    LCCollectionVec * telescopeCollectionVec3;

    // these encoders produce output - the encoding string should be the same as in the clustering and fit to the telescope
    CellIDEncoder < TrackerPulseImpl > outputPulseEncoder ( eutelescope::EUTELESCOPE::PULSEDEFAULTENCODING, outputTrackerPulseVec );
    CellIDEncoder < TrackerDataImpl > outputDataEncoder ( eutelescope::EUTELESCOPE::ZSCLUSTERDEFAULTENCODING, outputTrackerDataVec );
    CellIDEncoder < TrackerDataImpl > outputDataEncoder2 ( "sensorID:7,sparsePixelType:5,quality:5", outputTrackerDataVec2 );
    CellIDEncoder < TrackerPulseImpl > outputPulseEncoder2 ( "sensorID:7,xSeed:12,ySeed:12,xCluSize:5,yCluSize:5,type:5,quality:5", outputTrackerPulseVec );

    int alibavasize = 0;
    int telescopesize = 0;
    int telescopesize2 = 0;
    int telescopesize3 = 0;

    try
    {
	// the alibava data is passed from the event
	alibavaCollectionVec = dynamic_cast < LCCollectionVec * > ( alibavaEvent -> getCollection ( _alibavaCollectionName ) ) ;
	alibavasize = alibavaCollectionVec -> getNumberOfElements ( );
	streamlog_out ( DEBUG1 ) << alibavasize << " Elements in Alibava event!" << endl;

	// so is the secondary collection of the same size (we hope )
	alibavaCollectionVec2 = dynamic_cast < LCCollectionVec * > ( alibavaEvent -> getCollection ( _alibavaCollectionName2 ) ) ;

	// this decoder produces input
	CellIDDecoder < TrackerPulseImpl > inputPulseColDecoder2 ( alibavaCollectionVec2 );

	// first, loop over the alibava event
	for ( int i = 0; i < alibavasize; i++ )
	{
	    streamlog_out ( DEBUG1 ) << "Reading alibava..." << endl;
	    lcio::TrackerDataImpl * input  = dynamic_cast< lcio::TrackerDataImpl * > ( alibavaCollectionVec -> getElementAt ( i ) ) ;
	    lcio::TrackerDataImpl * output = new lcio::TrackerDataImpl;
	    output -> setChargeValues ( input -> getChargeValues ( ) ) ;
	    output -> setCellID0 ( input -> getCellID0 ( ) ) ;
	    output -> setCellID1 ( input -> getCellID1 ( ) ) ;
	    output -> setTime ( input -> getTime ( ) ) ;
	    outputTrackerDataVec -> addElement ( output );

	    lcio::TrackerPulseImpl * inputPulseFrame  = dynamic_cast < lcio::TrackerPulseImpl * > ( alibavaCollectionVec2 -> getElementAt ( i ) ) ;
	    lcio::TrackerPulseImpl * outputPulseFrame = new lcio::TrackerPulseImpl;
	    outputPulseFrame -> setCellID0 ( inputPulseFrame -> getCellID0 ( ) ) ;
	    outputPulseFrame -> setCellID1 ( inputPulseFrame -> getCellID1 ( ) ) ;
	    outputPulseFrame -> setTime ( inputPulseFrame -> getTime ( ) ) ;
	    outputPulseFrame -> setCharge ( inputPulseFrame -> getCharge ( ) ) ;
	    outputPulseFrame -> setQuality ( inputPulseFrame -> getQuality ( ) ) ;
	    outputPulseFrame -> setTrackerData ( inputPulseFrame -> getTrackerData ( ) ) ;
	    outputTrackerPulseVec -> addElement ( outputPulseFrame );

	    streamlog_out ( DEBUG1 ) << "Wrote alibava..." << endl;

	    // output correlation into a histogram
	    alibava_x += inputPulseColDecoder2 ( inputPulseFrame ) ["xSeed"];
	    alibava_y += inputPulseColDecoder2 ( inputPulseFrame ) ["ySeed"];
	    alibava_z += 0.0;

	    if ( _nonsensitiveaxis == "x" )
	    {
		ali_corr.push_back ( inputPulseColDecoder2 ( inputPulseFrame ) ["ySeed"] );
	    }

	    if ( _nonsensitiveaxis == "y" )
	    {
		ali_corr.push_back ( inputPulseColDecoder2 ( inputPulseFrame ) ["xSeed"] );
	    }

	}
    }
    catch ( lcio::DataNotAvailableException& )
    {
	streamlog_out ( DEBUG5 ) << "Collection " <<_alibavaCollectionName << " not found in event " << anEvent -> getEventNumber ( ) << endl;
    }

    try
    {
	if ( _eventdifferenceAlibava == 0 )
	{

	    // the telescope is read by the function
	    LCEvent* evt = readTelescope ( );
	    telescopeCollectionVec = dynamic_cast < LCCollectionVec * > ( evt -> getCollection ( _telescopeCollectionName ) ) ;
	    telescopesize = telescopeCollectionVec -> getNumberOfElements ( );
	    streamlog_out ( DEBUG1 ) << telescopesize << " Elements in Telescope event!" << endl;

	    // and the secondary collections
	    telescopeCollectionVec2 = dynamic_cast < LCCollectionVec * > ( evt -> getCollection ( _telescopeCollectionName2 ) ) ;
	    telescopesize2 = telescopeCollectionVec2 -> getNumberOfElements ( );
	    streamlog_out ( DEBUG1 ) << telescopesize2 << " Elements in Telescope event - Collection 2!" << endl;

	    // this guy can have less events and is not merged, but copied
	    telescopeCollectionVec3 = dynamic_cast < LCCollectionVec * > ( evt -> getCollection ( _outputCollectionName3 ) ) ;
	    telescopesize3 = telescopeCollectionVec3 -> getNumberOfElements ( );
	    streamlog_out ( DEBUG1 ) << telescopesize3 << " Elements in Telescope event - Collection 3!" << endl;

	    CellIDDecoder < TrackerPulseImpl > inputPulseColDecoder ( telescopeCollectionVec2 );
	    CellIDDecoder < TrackerDataImpl > inputSparseColDecoder ( telescopeCollectionVec );

	    // Cell ID Encoders for the telescope, introduced in eutelescope::EUTELESCOPE
	    // for sparseFrame (usually called cluster collection)
	    CellIDEncoder < TrackerDataImpl > outputSparseColEncoder ( eutelescope::EUTELESCOPE::ZSCLUSTERDEFAULTENCODING, outputTrackerDataVec );
	    // for pulseFrame
	    CellIDEncoder < TrackerPulseImpl > outputPulseColEncoder ( eutelescope::EUTELESCOPE::PULSEDEFAULTENCODING, outputTrackerPulseVec );

	    unsigned int noOfClusters;
	    // go through input clusters and copy them to output cluster collection

	    noOfClusters = telescopeCollectionVec2 -> getNumberOfElements ( );
	    for ( size_t i = 0; i < noOfClusters; ++i )
	    {
		TrackerPulseImpl * outputPulseFrame = new TrackerPulseImpl ( );
		TrackerDataImpl * outputSparseFrame = new TrackerDataImpl ( );

		TrackerPulseImpl* inputPulseFrame = dynamic_cast < TrackerPulseImpl* > ( telescopeCollectionVec2 -> getElementAt ( i ) );
		TrackerDataImpl* inputSparseFrame = dynamic_cast < TrackerDataImpl* > ( inputPulseFrame -> getTrackerData ( ) );
		
		// set Cell ID for sparse collection
		outputSparseColEncoder["sensorID"] = static_cast < int > ( inputSparseColDecoder ( inputSparseFrame ) ["sensorID"] - _teleplaneshift );
		outputSparseColEncoder["sparsePixelType"] =static_cast < int > ( inputSparseColDecoder ( inputSparseFrame )["sparsePixelType"] );
		outputSparseColEncoder["quality"] = static_cast < int > ( inputSparseColDecoder ( inputSparseFrame ) ["quality"] );
		outputSparseColEncoder.setCellID ( outputSparseFrame );

		// copy tracker data
		outputSparseFrame->setChargeValues ( inputSparseFrame -> getChargeValues ( ) );
		// add it to the cluster collection
		outputTrackerDataVec -> push_back ( outputSparseFrame );

		// prepare a pulse for this cluster
		outputPulseColEncoder["sensorID"] = static_cast < int> ( inputPulseColDecoder ( inputPulseFrame ) ["sensorID"] - _teleplaneshift );
		outputPulseColEncoder["type"] = static_cast < int > ( inputPulseColDecoder ( inputPulseFrame ) ["type"] );
		outputPulseColEncoder.setCellID ( outputPulseFrame );

		outputPulseFrame -> setCharge ( inputPulseFrame -> getCharge ( ) );
		outputPulseFrame -> setTrackerData ( outputSparseFrame );
		outputTrackerPulseVec -> push_back ( outputPulseFrame );

		// get these for the correlation plots
		EUTelVirtualCluster * tempCluster;
		tempCluster = new EUTelSparseClusterImpl < EUTelGenericSparsePixel > ( static_cast < TrackerDataImpl *> ( inputPulseFrame -> getTrackerData ( ) ) );
		float tempx = 0.0;
		float tempy = 0.0;
		tempCluster -> getCenterOfGravity ( tempx, tempy ) ;
		telescope_x += tempx;
		telescope_y += tempy;
		streamlog_out ( DEBUG0 ) << "Plot x is " << tempx << endl;
		streamlog_out ( DEBUG0 ) << "Plot y is " << tempy << endl;
		telescope_z = 0.0;
		delete tempCluster;

		if ( _nonsensitiveaxis == "x" )
		{
		    tele_corr.push_back ( tempy );
		}

		if ( _nonsensitiveaxis == "y" )
		{
		    tele_corr.push_back ( tempx );
		}

	    } // end of loop over input clusters

	    // one more loop over the telescope, now the third collection
	    // this collection is probably not needed, but for completion, we also copy it
	    for ( int k = 0; k < telescopesize3; k++ )
	    {
		streamlog_out ( DEBUG1 ) << "Reading telescope again..." << endl;
		lcio::TrackerDataImpl * input  = dynamic_cast < lcio::TrackerDataImpl * > ( telescopeCollectionVec3 -> getElementAt ( k ) ) ;
		lcio::TrackerDataImpl * output = new lcio::TrackerDataImpl;
		output->setChargeValues ( input -> getChargeValues ( ) ) ;
		output->setCellID0 ( input -> getCellID0( ) ) ;
		output->setCellID1 ( input -> getCellID1( ) ) ;
		output->setTime ( input -> getTime ( ) ) ;
		outputTrackerDataVec2 -> addElement ( output );
		streamlog_out ( DEBUG1 ) << "Wrote telescope again..." << endl;
	    }
	}

    }
    catch ( lcio::DataNotAvailableException& )
    {
	streamlog_out ( DEBUG5 ) << "Collection "<<_telescopeCollectionName<<" not found in event " << anEvent -> getEventNumber( ) << endl;
    }

    // correlation plot. Only makes sense if we have an event in both systems
    if ( alibavasize > 0 && telescopesize > 0 )
    {
	addCorrelation ( ( alibava_x / alibavasize ), ( alibava_y / alibavasize ),( alibava_z / alibavasize ), ( telescope_x / telescopesize ), ( telescope_y / telescopesize ), ( telescope_z / telescopesize ), anEvent -> getEventNumber ( ) );
	streamlog_out ( DEBUG0 ) << "Filling histogram with: Ali_X: " << ( alibava_x / alibavasize ) << " Ali_Y: " << ( alibava_y / alibavasize ) << " Ali_Z: " << ( alibava_z / alibavasize ) << " Tele_X: " << ( telescope_x / telescopesize ) << " Tele_Y: " << ( telescope_y /telescopesize ) << " Tele_Z: " << ( telescope_z / telescopesize ) << endl;

	for ( size_t i = 0; i < ali_corr.size ( ); i ++ )
	{
	    for ( size_t j = 0; j < tele_corr.size ( ); j++ )
	    {
		if ( TH2D * histo = dynamic_cast < TH2D* > ( _rootObjectMap["Correlation"] ) )
		{
		    histo -> Fill ( ali_corr.at ( i ), tele_corr.at ( j ) );
		}
	    }
	}
    }

    // Now the output. Usually the telescope will have more events...
    if (_outputmode == 0 )
    {
	if ( alibavasize > 0 && telescopesize > 0 )
	{
	    streamlog_out ( DEBUG1 ) << "Mode 0. Writing out event " << anEvent -> getEventNumber( ) << endl;
	    alibavaEvent -> addCollection ( outputTrackerDataVec, _outputCollectionName );
	    alibavaEvent -> addCollection ( outputTrackerPulseVec, _outputCollectionName2 );
	    alibavaEvent -> addCollection ( outputTrackerDataVec2, _outputCollectionName3 );
	}
    }
    else if ( _outputmode == 1 )
    {
	if ( alibavasize > 0 )
	{
	    streamlog_out ( DEBUG1 ) << "Mode 1. Writing out event " << anEvent -> getEventNumber( ) << endl;
	    alibavaEvent -> addCollection ( outputTrackerDataVec, _outputCollectionName );
	    alibavaEvent -> addCollection ( outputTrackerPulseVec, _outputCollectionName2 );
	    alibavaEvent -> addCollection ( outputTrackerDataVec2, _outputCollectionName3 );
	}
    }
    else if ( _outputmode == 2 )
    {
	if ( telescopesize > 0 )
	{
	    streamlog_out ( DEBUG1 ) << "Mode 2. Writing out event " << anEvent -> getEventNumber( ) << endl;
	    alibavaEvent -> addCollection ( outputTrackerDataVec, _outputCollectionName );
	    alibavaEvent -> addCollection ( outputTrackerPulseVec, _outputCollectionName2 );
	    alibavaEvent -> addCollection ( outputTrackerDataVec2, _outputCollectionName3 );
	}
    }
    else if (_outputmode == 3 )
    {
	streamlog_out ( DEBUG1 ) << "Mode 3. Writing out event " << anEvent -> getEventNumber( ) << endl;
	alibavaEvent -> addCollection ( outputTrackerDataVec, _outputCollectionName );
	alibavaEvent -> addCollection ( outputTrackerPulseVec, _outputCollectionName2 );
	alibavaEvent -> addCollection ( outputTrackerDataVec2, _outputCollectionName3 );
    }
    else
    {
	streamlog_out ( ERROR1 ) << "Unknown mode! Allowed values are 0 to 3!" << endl;
    }

    if ( _eventdifferenceAlibava > 0 )
    {
	_eventdifferenceAlibava--;
    }
}

void AlibavaMerger::check ( LCEvent * /* evt */ )
{

}

void AlibavaMerger::end( )
{
    // the telescope file is still open, we can now close it
    lcReader -> close ( ) ;
    delete lcReader ;
    streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
}

void AlibavaMerger::fillHistos ( )
{

}

void AlibavaMerger::bookHistos ( )
{
    if ( _nonsensitiveaxis == "x" )
    {
	TH2D * signalHisto = new TH2D ( "Correlation_X", "", 1152, 0, 1151, 1152, 0, 1151 );
	_rootObjectMap.insert ( make_pair ( "Correlation_X", signalHisto ) );
	signalHisto -> SetTitle ( "Correlation in X;Alibava Average Channel;Telescope Average Pixel" );

	TH2D * correlation = new TH2D ( "Correlation", "", 256, 0, 255, 576, 0, 575 );
	_rootObjectMap.insert ( make_pair ( "Correlation", correlation ) );
	correlation -> SetTitle ( "Correlation in Y;Alibava Channel;Telescope Pixel" );
    }
    else
    {
	TH2D * signalHisto = new TH2D ( "Correlation_X", "", 256, 0, 255, 1152, 0, 1151 );
	_rootObjectMap.insert ( make_pair ( "Correlation_X", signalHisto ) );
	signalHisto -> SetTitle ( "Correlation in X;Alibava Average Channel;Telescope Average Pixel" );
    }

    if (_nonsensitiveaxis == "y" )
    {
	TH2D * signalHisto2 = new TH2D ( "Correlation_Y", "", 576, 0, 575, 576, 0, 575 );
	_rootObjectMap.insert ( make_pair ( "Correlation_Y", signalHisto2 ) );
	signalHisto2 -> SetTitle ( "Correlation in Y;Alibava Average Channel;Telescope Average Pixel" );

	TH2D * correlation = new TH2D ( "Correlation", "", 256, 0, 255, 1152, 0, 1151 );
	_rootObjectMap.insert ( make_pair ( "Correlation", correlation ) );
	correlation -> SetTitle ( "Correlation in X;Alibava Channel;Telescope Pixel" );
    }
    else
    {
	TH2D * signalHisto2 = new TH2D ( "Correlation_Y", "", 256, 0, 255, 576, 0, 575 );
	_rootObjectMap.insert ( make_pair ( "Correlation_Y", signalHisto2 ) );
	signalHisto2 -> SetTitle ( "Correlation in Y;Alibava Average Channel;Telescope Average Pixel" );
    }

    TH2D * signalHisto3 = new TH2D ( "Correlation_Z", "", 256, 0, 255, 1152, 0, 1151 );
    _rootObjectMap.insert ( make_pair ( "Correlation_Z", signalHisto3 ) );
    signalHisto3 -> SetTitle ( "Correlation in Z;Alibava Average Channel;Telescope Average Pixel" );

    TH2D * signalHisto4 = new TH2D ( "Correlation_Event", "", 1000, 0, 500000, 100, -1, 1 );
    _rootObjectMap.insert ( make_pair ( "Correlation_Event", signalHisto4 ) );
    signalHisto4 -> SetTitle ( "Correlation over Events;Event Nr.;Telescope Average Pixel - ALiBaVa Average Strip" );

    TH3D * multihisto = new TH3D ( "Correlation_All", "", 50, -1, 1, 50, -1, 1, 1000, 0, 500000 );
    _rootObjectMap.insert ( make_pair ( "Correlation_All", multihisto ) );
    multihisto -> SetTitle ( "Correlation over Events;Average Distance in X;Average Distance in Y;Event Nr." );

    streamlog_out ( MESSAGE4 )  << "End of Booking histograms. " << endl;
}
