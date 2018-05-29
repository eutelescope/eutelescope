/*
 * Created by Thomas Eichhorn
 *  (2017 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/DataSourceProcessor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"
#include "marlin/ProcessorMgr.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <lcio.h>
#include <Exceptions.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCGenericObjectImpl.h>
#include <IMPL/LCRunHeaderImpl.h>
//#include <IMPL/TrackerPulseImpl.h>
//#include <IMPL/TrackerRawDataImpl.h>
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

#include "CMSBuncher.h"


using namespace std;
using namespace lcio;
using namespace marlin;
using namespace IMPL;
using namespace eutelescope;


CMSBuncher::CMSBuncher ( ) : DataSourceProcessor ( "CMSBuncher" )
{

    _description = "CMSBuncher groups events together, based on LCIO time stamps. This mimics a long read-out frame.";

    registerProcessorParameter ( "TelescopeFile", "The filename where the telescope data is stored", _telescopeFile, string ( "dummy_telescope.slcio" ) );

    registerProcessorParameter ( "InputCollectionName", "The name of the telescope collection we want to read", _inputCollectionName, string ( "collection1" ) );

    registerProcessorParameter ( "OutputCollectionName", "The name of the output collection we want to create", _outputCollectionName, string ( "collection2" ) );

    registerProcessorParameter ( "SingleFrameTime", "The time of a frame. Unit is micro seconds", _singleframetime, 115 );

}


CMSBuncher * CMSBuncher::newProcessor ( )
{
    return new CMSBuncher;
}


void CMSBuncher::init ( )
{
    streamlog_out ( MESSAGE4 ) << "Running init" << endl;
    printParameters ( );

    // set time
    _frametime = 1;
    _bunchtime = 0;

    _evtcount = 0;

    _maxevents = 5;

    _telescopeopen = false;
    _usestorageevent = false;
}


LCEvent *CMSBuncher::readTelescope ( )
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
	     streamlog_out ( MESSAGE4 ) << "Running over " << _maxevents << " events!" << endl;
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


void CMSBuncher::readDataSource ( int /*numEvents*/ )
{
    // write an almost empty run header
    LCRunHeaderImpl * lcHeader  = new LCRunHeaderImpl;
    lcHeader -> setRunNumber( 0 );
    ProcessorMgr::instance ( ) -> processRunHeader ( lcHeader ) ;
    delete lcHeader;

    // the input collections
    LCCollectionVec * telescopeCollectionVec;

    int bunchcount = 0;

    while ( bunchcount < ( _maxevents - 1 ) )
    {

	_bunchtime += _singleframetime;
	streamlog_out ( DEBUG4 ) << "Output event " << _evtcount << " has a time of " << _bunchtime << endl;

	LCCollectionVec* dataCollection = new LCCollectionVec( LCIO::TRACKERDATA );
	CellIDEncoder < TrackerDataImpl > encoder ( "sensorID:7,sparsePixelType:5", dataCollection );
	LCEventImpl * event = new LCEventImpl ( );
	event -> setRunNumber ( 0 );
	event -> setEventNumber ( _evtcount );
	event -> setDetectorName ( "CBCBuncher" );
	event -> parameters ( ) .setValue ( "EventType", 2 );
	event -> setTimeStamp ( _bunchtime );

	// prepare the output vector for this event
	vector < float > outputData[10];

	try
	{
	    // process this guy...
	    LCEvent* evt;

	    bool continueloop = true;

	    while ( continueloop == true )
	    {
		if ( _usestorageevent == false )
		{
		    // the telescope is read by the function
		    bunchcount++;
		    evt = readTelescope ( );
		    streamlog_out ( DEBUG4 ) << "Reading telescope bunch " << bunchcount << endl;
		}
		else
		{
		    evt = _storeevt;
		    _usestorageevent = false;
		    streamlog_out ( DEBUG2 ) << "Using stored event!" << endl;
		}

		_frametime = evt -> getTimeStamp ( );
		// save for future...
		_storeevt = evt;

		streamlog_out ( DEBUG2 ) << "Read telescope frame has a time of " << _frametime << endl;
		
		if ( _frametime < _bunchtime )
		{
		    streamlog_out ( DEBUG2 ) << "Adding a telescope frame!" << endl;
		}
		else
		{
		    _usestorageevent = true;
		    continueloop = false;
		    streamlog_out ( DEBUG2 ) << "Over time, breaking!" << endl;
		    break;
		}

		// now read the data
		telescopeCollectionVec = dynamic_cast < LCCollectionVec * > ( evt -> getCollection ( _inputCollectionName ) );
		int readsize = telescopeCollectionVec -> getNumberOfElements ( );
		for ( int j = 0; j < readsize; j++ )
		{
		    streamlog_out ( DEBUG1 ) << "Reading telescope ..." << endl;
		    CellIDDecoder < TrackerDataImpl > inputDecoder ( telescopeCollectionVec );
		    lcio::TrackerDataImpl * input  = dynamic_cast < lcio::TrackerDataImpl * > ( telescopeCollectionVec -> getElementAt ( j ) );
		    int sensorID = inputDecoder ( input )["sensorID"];
		    streamlog_out ( DEBUG0 ) << "Reading sensorID " << sensorID << endl;
		    FloatVec inputvec;
		    inputvec = input -> getChargeValues ( );

		    // assume 10 planes tops
		    for ( int i = 0; i < 10; i++ )
		    {
			if ( i == sensorID )
			{
			    for ( size_t k = 0; k < inputvec.size ( ); k++ )
			    {
				outputData[i].push_back ( inputvec.at ( k ) );
			    }
			}
		    }
		}
	    }

	    streamlog_out ( DEBUG4 ) << "Done accumulating, writing!" << endl;

	    for ( int i = 0; i < 6; i++ )
	    {
		TrackerDataImpl* planeData = new TrackerDataImpl ( );
		encoder["sensorID"] = i;
		encoder["sparsePixelType"] = 2;
		encoder.setCellID ( planeData );
		planeData -> setChargeValues ( outputData[i] );
		dataCollection -> addElement ( planeData );
	    }
	    event -> addCollection ( dataCollection, _outputCollectionName );
    }
    catch ( lcio::DataNotAvailableException& )
    {
	streamlog_out( ERROR1 ) << "Collection " << _inputCollectionName << " not found" << endl;
    }

    ProcessorMgr::instance ( ) -> processEvent ( event ) ;

    _evtcount++;
    delete event;

    }
}


void CMSBuncher::end ( )
{
    // the telescope file is still open, we can now close it
    lcReader -> close ( );
    delete lcReader;
    streamlog_out ( MESSAGE4 ) << "Successfully finished!" << endl;
}
