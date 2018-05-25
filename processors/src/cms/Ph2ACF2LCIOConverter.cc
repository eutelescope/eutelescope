/*
 * Created by Thomas Eichhorn
 *  (2017 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

#include "Ph2ACF2LCIOConverter.h"

// marlin includes
#include "marlin/Global.h"
#include "marlin/Exceptions.h"
#include "marlin/Processor.h"
#include "marlin/DataSourceProcessor.h"
#include "marlin/ProcessorMgr.h"

// lcio includes
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCGenericObjectImpl.h>
#include <IMPL/LCEventImpl.h>
#include <UTIL/CellIDEncoder.h>

// eutelescope includes
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"

// system includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <memory>
#include <stdlib.h>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <bitset>

using namespace std;
using namespace marlin;
using namespace eutelescope;

Ph2ACF2LCIOConverter::Ph2ACF2LCIOConverter ( ) : DataSourceProcessor ( "Ph2ACF2LCIOConverter" )

{
    _description = "Reads Ph2ACF data streams and converts to LCIO";

    registerProcessorParameter ( "DataFormat", "The Ph2ACF data type, options are 'raw' and 'slink'", _dataformat, string ( "raw" ) );

    registerProcessorParameter ( "InputFileName", "This is the input file name", _fileName, string ( "runXXXXXX.dat" ) );

    registerProcessorParameter ( "MaxRecordNumber", "The maximum number of events to read", _maxRecordNumber, -1 );

    registerProcessorParameter ( "NumberOfChips", "The number of CBC chips on each front end", _nChips, 2 );

    registerProcessorParameter ( "NumberOfFrontends", "The number of front ends connected", _nFE, 1 );

    registerProcessorParameter ( "RunNumber", "Formatted run number of file", _formattedRunNumber, string ( "0" ) );

    registerOutputCollection ( LCIO::TRACKERDATA, "BottomRawDataCollectionName", "Name of the collection for the bottom sensor", _rawDataCollectionNameBottom, string ( "rawdata2" ) );

    registerOutputCollection ( LCIO::TRACKERDATA, "TopRawDataCollectionName", "Name of the collection for the top sensor", _rawDataCollectionNameTop, string ( "rawdata1" ) );

}


Ph2ACF2LCIOConverter * Ph2ACF2LCIOConverter::newProcessor ( )
{
    return new Ph2ACF2LCIOConverter;
}


void Ph2ACF2LCIOConverter::init ( )
{
    printParameters ( );

}


void Ph2ACF2LCIOConverter::readDataSource ( int /* numEvents */ )
{
    // count the converted events
    int eventCounter = 0;

    // talk to streamlog
    streamlog::out.init ( std::cout , "Ph2ACF2LCIOConverter output stream" );
    streamlog::logscope scope ( streamlog::out );
    scope.setLevel < streamlog::MESSAGE > ( );

    streamlog_out ( DEBUG4 ) << "Reading " << _fileName << " with Ph2ACF2LCIOConverter!" << endl;
    _runNumber = atoi ( _formattedRunNumber.c_str ( ) );

    // open file
    ifstream infile;
    infile.open ( _fileName.c_str ( ) );
    if ( !infile.is_open ( ) )
    {
	streamlog_out ( ERROR5 ) << "Ph2ACF2LCIOConverter could not read the file " << _fileName << " correctly. Please check the path and file names that have been input!" << endl;
	exit ( -1 );
    }
    else
    {
	streamlog_out ( DEBUG4 ) << "Input file " << _fileName << " successfully opened!" << endl;
	if ( _dataformat == "raw" )
	{
	    streamlog_out ( DEBUG4 ) << "Assuming the file is encoded in RAW file format!" << endl;
	}
	else if ( _dataformat == "slink" )
	{
	    streamlog_out ( DEBUG4 ) << "Assuming the file is encoded in SLINK file format!" << endl;
	}
	else
	{
	    streamlog_out ( ERROR5 ) << "Unknown file format set! Valid inputs are 'raw' and 'slink'!" << endl;
	    exit ( -1 );
	}
    }

    LCRunHeaderImpl * runHeader = new LCRunHeaderImpl ( );
    runHeader -> setRunNumber ( _runNumber );
    runHeader -> setDetectorName ( "CBC" );
    ProcessorMgr::instance ( ) -> processRunHeader ( runHeader ) ;
    delete runHeader;

    // the header in raw file format
    bool raw_headeropen = false;
    bool maxreached = false;

    do
    {

	// header / event parameters

	// header 1
	unsigned int header1_size = 0;
	unsigned int fe_nbr = 0;
	unsigned int block_size = 0;
	unsigned int cic_id = 0;
	unsigned int chip_id = 0;
	unsigned int data_format_ver = 0;
	unsigned int dummy_size = 0;
	unsigned int trigdata_size = 0;
	unsigned int event_nbr = 0;
	unsigned int bx_cnt = 0;
	unsigned int stubdata_size = 0;
	unsigned int tlu_trigger_id = 0;
	unsigned int tdc = 0;

	// header 2, for each FE
	std::vector < unsigned int > chip_data_mask;
	std::vector < unsigned int > header2_size;
	std::vector < unsigned int > event_size;

	// cbc trigdata, for each chip
	std::vector < std::vector < unsigned int > > lat_err;
	std::vector < std::vector < unsigned int > > buf_ovf;
	std::vector < std::vector < unsigned int > > pipeaddr;
	std::vector < std::vector < unsigned int > > l1cnt;

	// cbc stubdata, for each chip
	std::vector < std::vector < unsigned int > > stub1;
	std::vector < std::vector < unsigned int > > stub2;
	std::vector < std::vector < unsigned int > > stub3;
	std::vector < std::vector < unsigned int > > bend1;
	std::vector < std::vector < unsigned int > > bend2;
	std::vector < std::vector < unsigned int > > bend3;

	if ( eventCounter > _maxRecordNumber && _maxRecordNumber > 0 )
	{
	    maxreached = true;
	    break ;
	}

	if ( eventCounter % 1000 == 0 || eventCounter < 10 )
	{
	    streamlog_out ( DEBUG4 ) << "Processing event " << eventCounter << " in run " << _runNumber << endl;
	}

	if ( _dataformat == "raw" )
	{

	    // the first event has a header
	    if ( raw_headeropen == false )
	    {
		uint32_t cMask = 0xAAAAAAAA;
		std::vector < uint32_t > headervec;
		raw_headeropen = true;
		streamlog_out ( DEBUG0 ) << "File Header: ";
		for ( int i = 0; i < 12; i++ )
		{
		    uint32_t tempint;
		    infile.read ( reinterpret_cast < char * > ( &tempint ), sizeof ( uint32_t ) );
		    headervec.push_back ( tempint );
		    streamlog_out ( DEBUG0 ) << tempint << " ";
		}
		streamlog_out ( DEBUG0 ) << endl;
		if ( headervec.at ( 0 ) == cMask && headervec.at ( 3 ) == cMask && headervec.at ( 6 ) == cMask && headervec.at ( 9 ) == cMask && headervec.at ( 11 ) == cMask )
		{
		    char cType[8] = { 0 };
		    cType[0] = ( headervec.at ( 1 ) && 0xFF000000 ) >> 24;
		    cType[1] = ( headervec.at ( 1 ) && 0x00FF0000 ) >> 16;
		    cType[2] = ( headervec.at ( 1 ) && 0x0000FF00 ) >> 8;
		    cType[3] = ( headervec.at ( 1 ) && 0x000000FF );

		    cType[4] = ( headervec.at ( 2 ) && 0xFF000000 ) >> 24;
		    cType[5] = ( headervec.at ( 2 ) && 0x00FF0000 ) >> 16;
		    cType[6] = ( headervec.at ( 2 ) && 0x0000FF00 ) >> 8;
		    cType[7] = ( headervec.at ( 2 ) && 0x000000FF );

		    std::string cTypeString ( cType );
		    std::string fType = cTypeString;

		    uint32_t fVersionMajor = headervec.at ( 4 );
		    uint32_t fVersionMinor = headervec.at ( 5 );

		    uint32_t fBeId = headervec.at ( 7 ) & 0x000003FF;
		    uint32_t fNCbc = headervec.at ( 8 );

		    uint32_t fEventSize32 = headervec.at ( 10 );
		    streamlog_out ( DEBUG4 ) << "Board Type: " << fType << endl;
		    streamlog_out ( DEBUG4 ) << "FWMajor: " << fVersionMajor << endl;
		    streamlog_out ( DEBUG4 ) << "FWMinor: " << fVersionMinor << endl;
		    streamlog_out ( DEBUG4 ) << "BeId: " << fBeId << endl;
		    streamlog_out ( DEBUG4 ) << "NCbc: " << fNCbc << endl;
		    streamlog_out ( DEBUG4 ) << "EventSize32: " << fEventSize32 << endl;
		    streamlog_out ( DEBUG4 ) << "Valid header!" << endl;
		}
		else
		{
		    streamlog_out ( ERROR5 ) << "Error, this is not a valid header!" << endl;
		    exit ( -1 );
		}

	    }

	    // the output vectors
	    FloatVec dataoutputvec_top;
	    FloatVec dataoutputvec_bot;

	    // read event header
	    streamlog_out ( DEBUG1 ) << endl;
	    streamlog_out ( DEBUG1 ) << "CBC Header1:" << endl;
	    uint32_t tempint;
	    std::vector < uint32_t > vec_header1;
	    for ( int i = 0; i < 5; i++ )
	    {
		infile.read ( reinterpret_cast < char * > ( &tempint ), sizeof ( uint32_t ) );
		vec_header1.push_back ( tempint );
	    }

	    for ( unsigned int i = 0 ; i < vec_header1.size ( ); i++ )
	    {
		if ( i == 0 )
		{
		    streamlog_out ( DEBUG3 ) << "Part 0: " << vec_header1.at ( i ) << endl;

		    header1_size = ( vec_header1.at ( i ) >> 24 );
		    streamlog_out ( DEBUG3 ) << " header1_size " << header1_size << endl;

		    fe_nbr = ( vec_header1.at ( i ) >> 16 ) & 0xFF;
		    streamlog_out ( DEBUG3 ) << " fe_nbr " << fe_nbr << endl;

		    block_size = ( ( ( vec_header1.at ( i ) >> 8 ) & 0xFF ) + ( ( vec_header1.at ( i ) ) & 0xFF ) );
		    streamlog_out ( DEBUG3 ) << " block_size " << block_size << endl;
		}
		if ( i == 1 )
		{
		    streamlog_out ( DEBUG3 ) << "Part 1: " << vec_header1.at ( i ) << endl;

		    cic_id = ( vec_header1.at ( i ) >> 24 );
		    streamlog_out ( DEBUG3 ) << " cic_id " << cic_id << endl;

		    chip_id = ( vec_header1.at ( i ) >> 16 ) & 0xFF;
		    streamlog_out ( DEBUG3 ) << " chip_id " << chip_id << endl;

		    data_format_ver = ( vec_header1.at ( i ) >> 8 ) & 0xFF;
		    streamlog_out ( DEBUG3 ) << " data_format_ver " << data_format_ver << endl;

		    dummy_size = ( vec_header1.at ( i ) ) & 0xFF;
		    streamlog_out ( DEBUG3 ) << " dummy_size " << dummy_size << endl;
		}
		if ( i == 2 )
		{
		    streamlog_out ( DEBUG3 ) << "Part 2: " << vec_header1.at ( i ) << endl;

		    trigdata_size = ( vec_header1.at ( i ) >> 24 );
		    streamlog_out ( DEBUG3 ) << " trigdata_size " << trigdata_size << endl;

		    event_nbr = ( ( ( vec_header1.at ( i ) >> 16 ) & 0xFF ) + ( ( vec_header1.at ( i ) >> 8 ) & 0xFF ) + ( ( vec_header1.at ( i ) ) & 0xFF ) );
		    streamlog_out ( DEBUG3 ) << " event_nbr " << event_nbr << endl;
		}
		if ( i == 3 )
		{
		    streamlog_out ( DEBUG3 ) << "Part 3: " << vec_header1.at ( i ) << endl;

		    bx_cnt = ( vec_header1.at ( i ) );
		    streamlog_out ( DEBUG3 ) << " bx_cnt " << bx_cnt << endl;
		}
		if ( i == 4 )
		{
		    streamlog_out ( DEBUG3 ) << "Part 4: " << vec_header1.at ( i ) << endl;

		    stubdata_size = ( vec_header1.at ( i ) >> 24 );
		    streamlog_out ( DEBUG3 ) << " stubdata_size " << stubdata_size << endl;

		    tlu_trigger_id = ( ( ( vec_header1.at ( i ) >> 16 ) & 0xFF ) + ( ( vec_header1.at ( i ) >> 8 ) & 0xFF ) );
		    streamlog_out ( DEBUG3 ) << " tlu_trigger_id " << tlu_trigger_id << endl;

		    tdc = ( vec_header1.at ( i ) ) & 0xFF;
		    streamlog_out ( DEBUG3 ) << " tdc " << tdc << endl;
		}
	    }

	    // loop frontends
	    for ( int iFE = 0; iFE < _nFE; iFE++ )
	    {
		// read header 2
		std::vector < uint32_t > vec_header2;
		uint32_t tempint;
		infile.read ( reinterpret_cast < char * > ( &tempint ), sizeof ( uint32_t ) );
		vec_header2.push_back ( tempint );
		streamlog_out ( DEBUG2 ) << endl;
		streamlog_out ( DEBUG2 ) << "CBC Header2, FE " << iFE << ":" << endl;
		for ( unsigned int i = 0; i < vec_header2.size ( ); i++ )
		{
		    if ( i == 0 )
		    {
			chip_data_mask.push_back ( vec_header2.at ( i ) >> 24 );
			streamlog_out ( DEBUG2 ) << " chip_data_mask " << chip_data_mask.at ( iFE ) << endl;

			header2_size.push_back ( ( vec_header2.at ( i ) >> 16 ) & 0xFF );
			streamlog_out ( DEBUG2 ) << " header2_size " << header2_size.at ( iFE ) << endl;

			event_size.push_back ( ( ( vec_header2.at ( i ) >> 8 ) & 0xFF ) + ( ( vec_header2.at ( i ) ) & 0xFF ) );
			streamlog_out ( DEBUG2 ) << " event_size " << event_size.at ( iFE ) << endl;
		    }
		}
		streamlog_out ( DEBUG2 ) << endl;

		// chip loop
		// push_back for the vectors...
		lat_err.push_back ( std::vector < unsigned int > ( ) );
		buf_ovf.push_back ( std::vector < unsigned int > ( ) );
		pipeaddr.push_back ( std::vector < unsigned int > ( ) );
		l1cnt.push_back ( std::vector < unsigned int > ( ) );
		stub1.push_back ( std::vector < unsigned int > ( ) );
		stub2.push_back ( std::vector < unsigned int > ( ) );
		stub3.push_back ( std::vector < unsigned int > ( ) );
		bend1.push_back ( std::vector < unsigned int > ( ) );
		bend2.push_back ( std::vector < unsigned int > ( ) );
		bend3.push_back ( std::vector < unsigned int > ( ) );

		for ( int j = 0; j < _nChips; j++ )
		{
		    // top sensor
		    int topvec[4][32] = { 0 };
		    // bottom sensor
		    int botvec[4][32] = { 0 };

		    // 9 for trg data, 2 for stub
		    for ( int i = 0; i < 11; i++ )
		    {
			uint32_t tempint;
			infile.read ( reinterpret_cast < char * > ( &tempint ), sizeof ( uint32_t ) );

			if ( i < 4 )
			{
			    // top sensor
			    bitset < 32 > tempbit ( tempint );
			    for ( unsigned j = 0; j < tempbit.size ( ); j++ )
			    {
				if ( tempbit.test ( j ) )
				{
				    topvec[i][j] = 1;
				}
				else
				{
				    topvec[i][j] = 0;
				}
			    }
			}
			if ( i >=4 && i < 8 )
			{
			    // bottom sensor
			    bitset < 32 > tempbit ( tempint );
			    for ( unsigned j = 0; j < tempbit.size ( ); j++ )
			    {
				if ( tempbit.test ( j ) )
				{
				    botvec[i - 4][j] = 1;
				}
				else
				{
				    botvec[i - 4][j] = 0;
				}
			    }
			}
			if ( i == 8 )
			{
			    // cbc trgdata status
			    unsigned int test =  ( ( tempint & 1 ) >> 1 );
			    lat_err[iFE].push_back ( test );
			    streamlog_out ( DEBUG1 ) << " lat_err " << lat_err.at ( iFE ).at ( j ) << endl;
			    buf_ovf[iFE].push_back ( ( tempint & 2 ) >> 2 );
			    streamlog_out ( DEBUG1 ) << " buf_ovf " << buf_ovf.at ( iFE ).at ( j ) << endl;
			    pipeaddr[iFE].push_back ( ( tempint >> 4 ) & 0x09 );
			    streamlog_out ( DEBUG1 ) << " pipeaddr " << pipeaddr.at ( iFE ).at ( j ) << endl;
			    l1cnt[iFE].push_back ( ( tempint >> 16 ) & 0xFE );
			    streamlog_out ( DEBUG1 ) << " l1cnt " << l1cnt.at ( iFE ).at ( j ) << endl;
			}
			if ( i == 9 )
			{
			    // stubdata
			    stub1[iFE].push_back ( createMask ( 0, 7 ) & tempint );
			    stub2[iFE].push_back ( createMask ( 0, 7 ) & ( tempint >> 8 ) );
			    stub3[iFE].push_back ( createMask ( 0, 7 ) & ( tempint >> 16 ) );
			    streamlog_out ( DEBUG1 ) << " stub1 " << stub1.at ( iFE ).at ( j ) << " stub2 " << stub2.at ( iFE ).at ( j ) << " stub3 " << stub3.at ( iFE ).at ( j ) << endl;
			}
			if ( i == 10 )
			{
			    // stubdata
			    unsigned int sync = ( ( tempint >> 3) & 1 );
			    unsigned int or254 = ( ( tempint >> 1 ) & 1 );
			    streamlog_out ( DEBUG1 ) << " sync " << sync << " or254 " << or254 << endl;
			    bend1[iFE].push_back ( createMask ( 0, 3 ) & ( tempint >> 8 ) );
			    bend2[iFE].push_back ( createMask ( 0, 3 ) & ( tempint >> 16 ) );
			    bend3[iFE].push_back ( createMask ( 0, 3 ) & ( tempint >> 24 ) );
			    streamlog_out ( DEBUG1 ) << " bend1 " << bend1.at ( iFE ).at ( j ) << " bend2 " << bend2.at ( iFE ).at ( j ) << " bend3 " << bend3.at ( iFE ).at ( j ) << endl;

			    // check
			    if ( stub1.at ( iFE).at ( j ) == 1 )
			    {
				if ( sync != 1 || or254 != 1 )
				{
				    streamlog_out ( WARNING1 ) << "Warning! Stub found, but sync/or254 is not 1!" << endl;
				}
			    }
			}
		    }

		    int counter = 0;
		    streamlog_out ( DEBUG0 ) << "Top ";
		    for ( int i = 3; i >= 0; i-- )
		    {
			for ( int k = 0; k < 32; k++ )
			{
			    counter++;
			    if ( counter != 32 )
			    {
				streamlog_out ( DEBUG0 ) << topvec[i][k] ;
				// output to the vector
				dataoutputvec_top.push_back ( topvec[i][k] );
			    }
			}
		    }
		    streamlog_out ( DEBUG0 ) << endl;
		    counter = 0;
		    streamlog_out ( DEBUG0 ) << "Bot ";
		    for ( int i = 3; i >= 0; i-- )
		    {
			for ( int k = 0; k < 32; k++ )
			{
			    counter ++;
			    if ( counter != 32 )
			    {
				streamlog_out ( DEBUG0 ) << botvec[i][k] ;
				// output to the vector
				dataoutputvec_bot.push_back ( botvec[i][k] );
			    }
			}
		    }
		    streamlog_out ( DEBUG0 ) << endl;

		} // done chip loop

	    } // done FE loop

	    // let there be output
	    EUTelEventImpl* anEvent = new EUTelEventImpl ( );
	    const char * dummyencode = "CBCRaw:1,";

	    LCCollectionVec* rawDataCollectionTop = new LCCollectionVec ( LCIO::TRACKERDATA );
	    CellIDEncoder < TrackerDataImpl > chipIDEncoderTop ( dummyencode, rawDataCollectionTop );
	    TrackerDataImpl * rawtop = new TrackerDataImpl ( );
	    rawtop -> setChargeValues ( dataoutputvec_top );
	    chipIDEncoderTop.setCellID ( rawtop );
	    rawDataCollectionTop -> push_back ( rawtop );
	    anEvent -> addCollection ( rawDataCollectionTop, _rawDataCollectionNameTop );

	    LCCollectionVec* rawDataCollectionBot = new LCCollectionVec ( LCIO::TRACKERDATA );
	    CellIDEncoder < TrackerDataImpl > chipIDEncoderBot ( dummyencode, rawDataCollectionBot );
	    TrackerDataImpl * rawbot = new TrackerDataImpl ( );
	    rawbot -> setChargeValues ( dataoutputvec_bot );
	    chipIDEncoderBot.setCellID ( rawbot );
	    rawDataCollectionBot -> push_back ( rawbot );
	    anEvent -> addCollection ( rawDataCollectionBot, _rawDataCollectionNameBottom );

	    anEvent -> setRunNumber ( _runNumber );
	    anEvent -> setEventNumber ( eventCounter );
	    anEvent -> setDetectorName ( "CBC" );
	    anEvent -> parameters ( ).setValue ( "EventType", 2 );

	    // now we set all the header parameters
	    anEvent -> parameters ( ).setValue ( "header1_size", int ( header1_size ) );
	    anEvent -> parameters ( ).setValue ( "fe_nbr", int ( fe_nbr ) );
	    anEvent -> parameters ( ).setValue ( "block_size", int ( block_size ) );
	    anEvent -> parameters ( ).setValue ( "cic_id", int ( cic_id ) );
	    anEvent -> parameters ( ).setValue ( "chip_id", int ( chip_id ) );
	    anEvent -> parameters ( ).setValue ( "data_format_ver", int ( data_format_ver ) );
	    anEvent -> parameters ( ).setValue ( "dummy_size", int ( dummy_size ) );
	    anEvent -> parameters ( ).setValue ( "trigdata_size", int ( trigdata_size ) );
	    anEvent -> parameters ( ).setValue ( "event_nbr", int ( event_nbr ) );
	    anEvent -> parameters ( ).setValue ( "bx_cnt", int ( bx_cnt ) );
	    anEvent -> parameters ( ).setValue ( "stubdata_size", int ( stubdata_size ) );
	    anEvent -> parameters ( ).setValue ( "tlu_trigger_id", int ( tlu_trigger_id ) );
	    anEvent -> parameters ( ).setValue ( "tdc", int ( tdc ) );

	    for ( int i = 0; i < _nFE; i++ )
	    {
		anEvent -> parameters ( ).setValue ( "chip_data_mask_" + std::to_string ( i ), int ( chip_data_mask.at ( i ) ) );
		anEvent -> parameters ( ).setValue ( "header2_size_" + std::to_string ( i ), int ( header2_size.at ( i ) ) );
		anEvent -> parameters ( ).setValue ( "event_size_" + std::to_string ( i ), int ( event_size.at ( i ) ) );

		for ( int j = 0; j < _nChips; j++ )
		{
		    anEvent -> parameters ( ).setValue ( "l1cnt_" + std::to_string ( i ) + "_" + std::to_string ( j ), int ( l1cnt.at ( i ).at ( j ) ) );
		    anEvent -> parameters ( ).setValue ( "pipeaddr" + std::to_string ( i ) + "_" + std::to_string ( j ), int ( pipeaddr.at ( i ).at ( j ) ) );
		    anEvent -> parameters ( ).setValue ( "buf_ovf" + std::to_string ( i ) + "_" + std::to_string ( j ), int ( buf_ovf.at ( i ).at ( j ) ) );
		    anEvent -> parameters ( ).setValue ( "lat_err" + std::to_string ( i ) + "_" + std::to_string ( j ), int ( lat_err.at ( i ).at ( j ) ) );
		    anEvent -> parameters ( ).setValue ( "stub1_" + std::to_string ( i ) + "_" + std::to_string ( j ), int ( stub1.at ( i ).at ( j ) ) );
		    anEvent -> parameters ( ).setValue ( "stub2_" + std::to_string ( i ) + "_" + std::to_string ( j ), int ( stub2.at ( i ).at ( j ) ) );
		    anEvent -> parameters ( ).setValue ( "stub3_" + std::to_string ( i ) + "_" + std::to_string ( j ), int ( stub3.at ( i ).at ( j ) ) );
		    anEvent -> parameters ( ).setValue ( "bend1_" + std::to_string ( i ) + "_" + std::to_string ( j ), int ( bend1.at ( i ).at ( j ) ) );
		    anEvent -> parameters ( ).setValue ( "bend2_" + std::to_string ( i ) + "_" + std::to_string ( j ), int ( bend2.at ( i ).at ( j ) ) );
		    anEvent -> parameters ( ).setValue ( "bend3_" + std::to_string ( i ) + "_" + std::to_string ( j ), int ( bend3.at ( i ).at ( j ) ) );
		}

	    }

	    // FIXME this will be the TLU trigger ID
	    anEvent -> setTimeStamp ( long64 ( eventCounter * 100.0 ) );

	    ProcessorMgr::instance ( ) -> processEvent ( static_cast < LCEventImpl* > ( anEvent ) ) ;
	    eventCounter++;
	    delete anEvent;
	
	} // done _dataformat if
	else if ( _dataformat == "slink" )
	{

	    // FIXME

	    // the output vectors
	    FloatVec dataoutputvec_top;
	    FloatVec dataoutputvec_bot;

	    // FIXME
	    for ( int i = 0; i < 19; i++ )
	    {
		uint32_t tempint;
		std::vector < uint32_t > inputvec;
		infile.read ( reinterpret_cast < char * > ( &tempint ), sizeof ( uint32_t ) );
		inputvec.push_back ( tempint );
		streamlog_out ( DEBUG0 ) << tempint << " ";
	    }
	    streamlog_out ( DEBUG0 ) << endl;

	    // let there be output
	    EUTelEventImpl* anEvent = new EUTelEventImpl ( );
	    const char * dummyencode = "CBCRaw:1,";

	    LCCollectionVec* rawDataCollectionTop = new LCCollectionVec ( LCIO::TRACKERDATA );
	    CellIDEncoder < TrackerDataImpl > chipIDEncoderTop ( dummyencode, rawDataCollectionTop );
	    TrackerDataImpl * rawtop = new TrackerDataImpl ( );
	    rawtop -> setChargeValues ( dataoutputvec_top );
	    chipIDEncoderTop.setCellID ( rawtop );
	    rawDataCollectionTop -> push_back ( rawtop );
	    anEvent -> addCollection ( rawDataCollectionTop, _rawDataCollectionNameTop );

	    LCCollectionVec* rawDataCollectionBot = new LCCollectionVec ( LCIO::TRACKERDATA );
	    CellIDEncoder < TrackerDataImpl > chipIDEncoderBot ( dummyencode, rawDataCollectionBot );
	    TrackerDataImpl * rawbot = new TrackerDataImpl ( );
	    rawbot -> setChargeValues ( dataoutputvec_bot );
	    chipIDEncoderBot.setCellID ( rawbot );
	    rawDataCollectionBot -> push_back ( rawbot );
	    anEvent -> addCollection ( rawDataCollectionBot, _rawDataCollectionNameBottom );

	    anEvent -> setRunNumber ( _runNumber );
	    anEvent -> setEventNumber ( eventCounter );
	    anEvent -> setDetectorName ( "CBC" );
	    anEvent -> parameters ( ).setValue ( "EventType", 2 );

	    // FIXME this will be the TLU trigger ID
	    anEvent -> setTimeStamp ( long64 ( eventCounter * 100.0 ) );

	    ProcessorMgr::instance ( ) -> processEvent ( static_cast < LCEventImpl* > ( anEvent ) ) ;
	    eventCounter++;
	    delete anEvent;

	} // done _dataformat if

    } while ( !( infile.bad ( ) || infile.eof ( ) ) || maxreached == true );

    infile.close ( );

}


void Ph2ACF2LCIOConverter::end ( )
{
    streamlog_out ( MESSAGE5 )  << "Ph2ACF2LCIOConverter successfully finished!" << endl;
}

unsigned Ph2ACF2LCIOConverter::createMask ( unsigned a, unsigned b )
{
    unsigned r = 0;
    for ( unsigned i = a; i <= b; i++ )
    {
	r |= 1 << i;
    }
    return r;
}
