
/*
 * Created by Thomas Eichhorn
 *  (2017 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */


// ROOT includes:

// eutelescope includes ".h"
#include "CMSStubGenerator.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelGenericSparsePixel.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelGenericSparseClusterImpl.h"
#include "EUTelSimpleVirtualCluster.h"
#include "EUTelVirtualCluster.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"
#include "marlin/AIDAProcessor.h"

#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>

// lcio includes <.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

// system includes <>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>
#include <iostream>
#include <iomanip>
#include <cstdio>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace gear;
using namespace IMPL;
using namespace eutelescope;

AIDA::IHistogram1D * stubdistx;
AIDA::IHistogram1D * stubdisty;
AIDA::IHistogram1D * stubdistx_bit;
AIDA::IHistogram1D * stubdisty_bit;
AIDA::IHistogram1D * faildistx;
AIDA::IHistogram1D * faildisty;
AIDA::IHistogram2D * correx;
AIDA::IHistogram2D * correy;
AIDA::IHistogram2D * stubdt;
AIDA::IHistogram2D * stubratiodt;
AIDA::IHistogram1D * stubmap_top_x;
AIDA::IHistogram1D * stubmap_bot_x;
AIDA::IHistogram1D * stubmap_top_y;
AIDA::IHistogram1D * stubmap_bot_y;


CMSStubGenerator::CMSStubGenerator ( ) : Processor ( "CMSStubGenerator" )
{
    // modify processor description
    _description =  "CMSStubGenerator merges the two CBC sensors' hits into stubs, based on the cluster positions. In modes 1 and 2, it can also drop either sensor's hits to have only one sensor remaining.";

    registerInputCollection ( LCIO::TRACKERHIT, "InputHitCollectionName", "Input hit collection name. Hits should be in global coordinates and pre-aligned", _inputHitCollectionName, std::string ( " " ) );

    registerOutputCollection ( LCIO::TRACKERHIT, "OutputHitCollectionName", "Output hit collection name", _outputHitCollectionName, std::string ( " " ) );

    registerProcessorParameter ( "DUTPlane1", "This is the first DUT sensorID.", _dutPlane1, 6 );

    registerProcessorParameter ( "DUTPlane2", "This is the second DUT sensorID.", _dutPlane2, 7 );

    registerProcessorParameter ( "KeepDUTHits", "Keep the DUT hits in mode 0 after creating subs or discard them?", _keepDUTHits, true );

    registerProcessorParameter ( "MaxResidual", "Maximum distance in channels for two clusters to constitute a stub", _maxResidual, 0.0f );

    registerProcessorParameter ( "Mode", "0 to merge hits. 1 to drop hits from plane 1, 2 to drop hits from plane 2", _runMode, 0 );

    registerProcessorParameter ( "OutputSensorID", "SensorID of the output stub.", _outputSensorID, 8 );

    registerProcessorParameter ( "RequireStub", "Do we require an event to have the stub flag set to create an offline stub? 1 for on, 0 for off.", _requirestubflag, 1 );

}


void CMSStubGenerator::init ( )
{

    printParameters ( );

    bookHistos ( );

    _totalstubs = 0;
    _totalpl1 = 0;
    _totalpl2 = 0;

}


void CMSStubGenerator::processRunHeader ( LCRunHeader * rdr )
{
    std::unique_ptr < EUTelRunHeaderImpl > header = std::make_unique < EUTelRunHeaderImpl > ( rdr );
    header -> addProcessor ( type ( ) );

}


void CMSStubGenerator::processEvent ( LCEvent * event )
{
    int stubsinthisevent = 0;

    EUTelEventImpl * evt = static_cast < EUTelEventImpl* > ( event );

    if ( evt -> getEventType ( ) == kEORE )
    {
	streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
        return;
    }
    else if ( evt -> getEventType ( ) == kUNKNOWN )
    {
	streamlog_out ( WARNING2 ) << "Event number " << evt -> getEventNumber ( ) << " in run " << evt -> getRunNumber ( ) << " is of unknown type. Continue considering it as a normal data event." << endl;
    }

    LCCollectionVec * inputHitCollection = nullptr;
    LCCollectionVec * outputHitCollection = nullptr;

    try
    {
	inputHitCollection = static_cast < LCCollectionVec* > ( event -> getCollection ( _inputHitCollectionName ) );
    }
    catch ( DataNotAvailableException& e )
    {
	streamlog_out ( MESSAGE2 ) << "No input collection " << _inputHitCollectionName << " found in event " << event -> getEventNumber ( ) << " in run " << event -> getRunNumber ( ) << endl;
	return;
    }

    try
    {
	outputHitCollection = static_cast < LCCollectionVec* > ( event -> getCollection ( _outputHitCollectionName ) );
    }
    catch ( ... )
    {
	outputHitCollection = new LCCollectionVec ( LCIO::TRACKERHIT );
    }

    // the stub flags in the data stream
    std::string stub3 = evt -> getParameters ( ) .getStringVal ( "stub_pos_00_03_0" );
    std::string stub2 = evt -> getParameters ( ) .getStringVal ( "stub_pos_00_02_0" );
    std::string stub1 = evt -> getParameters ( ) .getStringVal ( "stub_pos_00_01_0" );

    // prepare an encoder for the hit collection
    CellIDEncoder < TrackerHitImpl > outputCellIDEncoder ( EUTELESCOPE::HITENCODING, outputHitCollection );

    CellIDDecoder < TrackerHitImpl > inputCellIDDecoder ( inputHitCollection );

    vector < int > dutPlane1Hits;
    vector < int > dutPlane2Hits;

    for ( int iInputHits = 0; iInputHits < inputHitCollection -> getNumberOfElements ( ); iInputHits++ )
    {
	TrackerHitImpl * inputHit = dynamic_cast < TrackerHitImpl* > ( inputHitCollection -> getElementAt ( iInputHits ) );
        int sensorID = inputCellIDDecoder ( inputHit ) ["sensorID"];
	
	if ( _runMode == 0 )
	{
	    bool isDUTHit = false;

	    if ( sensorID == _dutPlane1 )
	    {
		dutPlane1Hits.push_back ( iInputHits );
		isDUTHit = true;
	    }
	    if ( sensorID == _dutPlane2 )
	    {
		dutPlane2Hits.push_back ( iInputHits );
		isDUTHit = true;
	    }

	    if ( !isDUTHit || _keepDUTHits )
	    {
		outputHitCollection -> push_back ( cloneHit ( inputHit ) );
	    }
	}
	else if ( _runMode == 1 )
	{
	    if ( sensorID != _dutPlane1 )
	    {
		outputHitCollection -> push_back ( cloneHit ( inputHit ) );
	    }
	}
	else if ( _runMode == 2 )
	{
	    if ( sensorID != _dutPlane2 )
	    {
		outputHitCollection -> push_back ( cloneHit ( inputHit ) );
	    }
	}
	else
	{
	    streamlog_out ( ERROR5 ) << "Illegal setting for Mode! Only 0, 1 or 2 is a valid setting!" << endl;
	    exit ( -1 );
	}
    }

    if ( _runMode == 0)
    {
	for ( unsigned int iHitPlane1 = 0; iHitPlane1 < dutPlane1Hits.size ( ); iHitPlane1++ )
	{
	    TrackerHitImpl * Hit1 = dynamic_cast < TrackerHitImpl* > ( inputHitCollection -> getElementAt ( dutPlane1Hits[iHitPlane1] ) );
	    float x1 = -1.0;
	    float y1 = -1.0;
	    float q1 = -1.0;
	    const double* pos1 = Hit1 -> getPosition ( );
	    TrackerDataImpl* clusterVector1 = static_cast < TrackerDataImpl* > ( Hit1 -> getRawHits ( ) [0] );
	    EUTelSimpleVirtualCluster * cluster1 = nullptr;
	    cluster1 = new EUTelSparseClusterImpl < EUTelGenericSparsePixel > ( clusterVector1 );
	    if ( cluster1 != nullptr )
	    {
		cluster1 -> getCenterOfGravity ( x1, y1 );
		q1 = cluster1 -> getTotalCharge ( );
	    }

	    for ( unsigned int iHitPlane2 = 0; iHitPlane2 < dutPlane2Hits.size ( ); iHitPlane2++ )
	    {
		TrackerHitImpl * Hit2 = dynamic_cast < TrackerHitImpl* > ( inputHitCollection -> getElementAt ( dutPlane2Hits[iHitPlane2] ) );
		float x2 = -1.0;
		float y2 = -1.0;
		float q2 = -1.0;
		const double* pos2 = Hit2 -> getPosition ( );
		TrackerDataImpl* clusterVector2 = static_cast < TrackerDataImpl* > ( Hit2 -> getRawHits ( ) [0] );
		EUTelSimpleVirtualCluster * cluster2 = nullptr;
		cluster2 = new EUTelSparseClusterImpl < EUTelGenericSparsePixel > ( clusterVector2 );
		if ( cluster2 != nullptr )
		{
		    cluster2 -> getCenterOfGravity ( x2, y2 );
		    q2 = cluster2 -> getTotalCharge ( );
		}
		streamlog_out ( DEBUG0 ) << " x1 " << x1 << " y1 " << y1 << " q1 " << q1 << " x2 " << x2 << " y2 " << y2 << " q2 " << q2 << " evt " << evt -> getRunNumber ( ) << endl;

		correx -> fill ( x1, x2 );
		correy -> fill ( y1, y2 );

		float dx = -1.0;
		float dy = -1.0;
		dx = fabs ( x1 - x2 );
		dy = fabs ( y1 - y2 );

		if ( dx < _maxResidual && dy < _maxResidual )
		{
		    double newPos[3];
		    newPos[0] = ( pos1[0] + pos2[0] ) / 2.0;
		    newPos[1] = ( pos1[1] + pos2[1] ) / 2.0;
		    newPos[2] = ( pos1[2] + pos2[2] ) / 2.0;

		    const double* hitpos = newPos;
		    CellIDEncoder < TrackerHitImpl > idHitEncoder ( EUTELESCOPE::HITENCODING, outputHitCollection );
		    TrackerHitImpl* hit = new TrackerHitImpl;
		    hit -> setPosition ( &hitpos[0] );
		    float cov[TRKHITNCOVMATRIX] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
		    hit -> setCovMatrix ( cov );
		    hit -> setType ( kEUTelGenericSparseClusterImpl );
		    // assume all times are equal
		    hit -> setTime ( Hit1 -> getTime ( ) );

		    LCObjectVec clusterVec;
		    clusterVec.push_back ( clusterVector1 );
		    clusterVec.push_back ( clusterVector2 );

		    hit -> rawHits ( ) = clusterVec;

		    idHitEncoder["sensorID"] =  _outputSensorID ;
		    idHitEncoder["properties"] = 0;

		    idHitEncoder.setCellID ( hit );

		    stubdistx -> fill ( x1 - x2 );
		    stubdisty -> fill ( y1 - y2 );

		    bool bitpresent = false;
		    bool writeoutput = true;

		    if ( stub1 == "1" || stub2 == "1" || stub3 == "1" )
		    {
			stubdistx_bit -> fill ( x1 - x2 );
			stubdisty_bit -> fill ( y1 - y2 );

			bitpresent = true;


		    }

		    if ( _requirestubflag == 1 && bitpresent == false )
		    {
			writeoutput = false;
		    }

		    if ( writeoutput == true )
		    {
			outputHitCollection -> push_back ( hit );
			_totalstubs++;
			stubsinthisevent++;
			stubmap_top_x -> fill ( x1 );
			stubmap_bot_x -> fill ( x2 );
			stubmap_top_y -> fill ( y1 );
			stubmap_bot_y -> fill ( y2 );
		    }

		}
		else
		{
		    faildistx -> fill ( x1 - x2 );
		    faildisty -> fill ( y1 - y2 );
		}
	    }
	}
	_totalpl1 += dutPlane1Hits.size ();
	_totalpl2 += dutPlane2Hits.size ();
    }

    try
    {
	event -> getCollection ( _outputHitCollectionName ) ;
    }
    catch ( ... )
    {
	event -> addCollection ( outputHitCollection, _outputHitCollectionName );
    }

    double ratio = stubsinthisevent * 1.0 / ( ( dutPlane1Hits.size ( ) + dutPlane2Hits.size ( ) ) / 2.0 );
    stubratiodt -> fill ( event -> getEventNumber ( ), ratio );
    stubdt -> fill ( ( event -> getEventNumber ( ) ) * 1.0, ( ( _totalstubs * 1.0 ) / ( event -> getEventNumber ( ) * 1.0 ) ) );

}


TrackerHitImpl* CMSStubGenerator::cloneHit ( TrackerHitImpl *inputHit )
{
    TrackerHitImpl * newHit = new TrackerHitImpl;

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

    return newHit;
}


void CMSStubGenerator::end ( )
{

    streamlog_out ( MESSAGE4 ) << "Created " << _totalstubs << " stubs from on average " << ( _totalpl1 + _totalpl2 ) / 2.0 << " hits -> " << _totalstubs / ( ( _totalpl1 + _totalpl2 ) / 2.0 ) * 100.0 << "% efficiency!" << endl;

}

void CMSStubGenerator::bookHistos( )
{
    try
    {
	streamlog_out ( MESSAGE2 ) << "Booking histograms..." << endl;
	
	stubdistx = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( "Stub Distance in x", 101, -50, 50 );
	stubdistx -> setTitle ( "Stub Distance in x;x_{0} - x_{1} [channels];Entries" );

	stubdisty = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( "Stub Distance in y", 101, -50, 50 );
	stubdisty -> setTitle ( "Stub Distance in y;y_{0} - y_{1} [channels];Entries" );

	stubdistx_bit = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( "Stub Distance in x, bit required", 101, -50, 50 );
	stubdistx_bit -> setTitle ( "Stub Distance in x, bit required;x_{0} - x_{1} [channels];Entries" );

	stubdisty_bit = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( "Stub Distance in y, bit required", 101, -50, 50 );
	stubdisty_bit -> setTitle ( "Stub Distance in y, bit required;y_{0} - y_{1} [channels];Entries" );

	faildistx = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( "Fail Distance in x", 1000, -500, 500 );
	faildistx -> setTitle ( "Fail Distance in x;x_{0} - x_{1} [channels];Entries" );

	faildisty = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( "Fail Distance in y", 1000, -500, 500 );
	faildisty -> setTitle ( "Fail Distance in y;y_{0} - y_{1} [channels];Entries" );

	correx = AIDAProcessor::histogramFactory ( this ) -> createHistogram2D ( "Cluster correlation in x", 1016, 0, 1015, 1016, 0, 1015 );
	correx -> setTitle ( "Cluster Correlation in x;x_{0} [Channel];x_{1} [Channel]" );

	correy = AIDAProcessor::histogramFactory ( this ) -> createHistogram2D ( "Cluster correlation in y", 1016, 0, 1015, 1016, 0, 1015 );
	correy -> setTitle ( "Cluster Correlation in y;y_{0} [Channel];y_{1} [Channel]" );

	stubdt = AIDAProcessor::histogramFactory ( this ) -> createHistogram2D ( "Stubs per Event", 10000, 0, 500000, 100, 0, 1 );
	stubdt -> setTitle ( "Stubs per Event;Events;Stubs per Event" );

	stubratiodt = AIDAProcessor::histogramFactory ( this ) -> createHistogram2D ( "Stubs per Cluster Ratio over Events", 10000, 0, 500000, 100, 0, 1 );
	stubratiodt -> setTitle ( "Stubs per Cluster Ratio over Events;Events;Stubs per Cluster" );

	stubmap_top_x = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( "Stub map, top, x", 1016, 0, 1015 );
	stubmap_top_x -> setTitle ( "Stub Map, Top Sensor in x;Channel;Entries" );

	stubmap_bot_x = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( "Stub map, bot, x", 1016, 0, 1015 );
	stubmap_bot_x -> setTitle ( "Stub Map, Bottom Sensor in x;Channel;Entries" );

	stubmap_top_y = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( "Stub map, top, y", 1016, 0, 1015 );
	stubmap_top_y -> setTitle ( "Stub Map, Top Sensor in y;Channel;Entries" );

	stubmap_bot_y = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( "Stub map, bot, y", 1016, 0, 1015 );
	stubmap_bot_y -> setTitle ( "Stub Map, Bottom Sensor in y;Channel;Entries" );

    }
    catch ( ... )
    {

    }
}
