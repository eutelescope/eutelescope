/*
 * Created by Thomas Eichhorn
 *  (2014, 2017 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

// alibava includes ".h"
#include "AlibavaHeader.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

#if defined ( USE_AIDA ) || defined ( MARLIN_USE_AIDA )
// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCEventImpl.h>

// ROOT includes ".h"
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include "TROOT.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TProfile.h"

// system includes <>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <memory>
#include <fstream>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;


AlibavaHeader::AlibavaHeader ( ) : AlibavaBaseProcessor ( "AlibavaHeader" )
{

    _description = "AlibavaHeader does some analysis of the beetle chip header data to determine cross-talk. Input of raw header and channel data expected!";

    registerInputCollection ( LCIO::TRACKERDATA, "InputCollectionName", "Header collection name", _inputCollectionName, string ( "chipheader" ) );

    registerProcessorParameter ( "RawDataCollection", "The collection the data of channel 0 is stored in", _rawdatacollection , string ( "rawdata" ) );

    registerOptionalParameter ( "OutputFile", "The filename to write coefficients to", _filterFileName, string ( "" ) );

}

void AlibavaHeader::init ( )
{

    printParameters ( );

    if ( _filterFileName.length ( ) > 0 )
    {
	_coefffile = true;
	streamlog_out ( MESSAGE4 ) << "Writing coefficient file to " << _filterFileName << "." << endl;
    }
    else
    {
	_coefffile = false;
	streamlog_out ( MESSAGE4 ) << "Not writing coefficient file." << endl;
    }
}

void AlibavaHeader::processRunHeader ( LCRunHeader * rdr )
{

    auto arunHeader = std::make_unique < AlibavaRunHeaderImpl > ( rdr );
    arunHeader -> addProcessor ( type ( ) );

    bookHistos ( );

}

void AlibavaHeader::processEvent ( LCEvent * anEvent )
{

    if ( anEvent -> getEventNumber ( ) % 1000 == 0 )
    {
	streamlog_out ( MESSAGE4 ) << "Looping events " << anEvent -> getEventNumber ( ) << endl;
    }

    AlibavaEventImpl * alibavaEvent = static_cast < AlibavaEventImpl* > ( anEvent );

    LCCollectionVec * headColVec;
    LCCollectionVec * chanColVec;
    EVENT::IntVec chipSelection = getChipSelection ( );
    try
    {
	headColVec = dynamic_cast < LCCollectionVec * > ( alibavaEvent -> getCollection ( getInputCollectionName ( ) ) ) ;
	chanColVec = dynamic_cast < LCCollectionVec * > ( alibavaEvent -> getCollection ( _rawdatacollection ) ) ;

	for ( unsigned int j = 0; j < 2; j++ )
	{
	    TrackerDataImpl * headerdata = dynamic_cast < TrackerDataImpl * > ( headColVec -> getElementAt ( j ) ) ;
	    TrackerDataImpl * channeldata = dynamic_cast < TrackerDataImpl * > ( chanColVec -> getElementAt ( j ) ) ;
	    fillHistos ( headerdata, channeldata, j );
	    correlateLastHeader ( headerdata, channeldata, j );
	}

    }
    catch ( ... )
    {
	streamlog_out ( ERROR5 ) << "Collection ( " << getInputCollectionName ( ) << " ) not found!" << endl;
    }
}

void AlibavaHeader::check ( LCEvent * )
{

}

void AlibavaHeader::end ( )
{
    streamlog_out ( MESSAGE4 ) << "Successfully finished!" << endl;

    for ( int ichip =0; ichip < 2; ichip++ )
    {
	float x1 = 0.0;
	float x2 = 0.0;
	float x3 = 0.0;
	float x4 = 0.0;

	float y1 = 0.0;
	float y2 = 0.0;
	float y3 = 0.0;
	float y4 = 0.0;

	char tmpchar[100];

	sprintf ( tmpchar, "Chip %d, Low Low Signals", ichip );
	if ( TH2D * histo = dynamic_cast < TH2D* > ( _rootObjectMap[tmpchar] ) )
	{
	    x1 = histo -> GetMean ( 1 );
	    y1 = histo -> GetMean ( 2 );
	}
	sprintf ( tmpchar, "Chip %d, Low High Signals", ichip );
	if ( TH2D * histo = dynamic_cast < TH2D* > ( _rootObjectMap[tmpchar] ) )
	{
	    x2 = histo -> GetMean ( 1 );
	    y2 = histo -> GetMean ( 2 );
	}
	sprintf ( tmpchar, "Chip %d, High Low Signals", ichip );
	if ( TH2D * histo = dynamic_cast < TH2D* > ( _rootObjectMap[tmpchar] ) )
	{
	    x3 = histo -> GetMean ( 1 );
	    y3 = histo -> GetMean ( 2 );
	}
	sprintf ( tmpchar, "Chip %d, High High Signals", ichip );
	if ( TH2D * histo = dynamic_cast < TH2D* > ( _rootObjectMap[tmpchar] ) )
	{
	    x4 = histo -> GetMean ( 1 );
	    y4 = histo -> GetMean ( 2 );
	}

	float c1 = 0.0;
	float c2 = 0.0;

	c1 = ( y1 / x1 + y2 / x2 + y3 / x3 + y4 / x4 ) / 4.0;

	// get a value for 14:
	float tempx = ( fabs ( x1 ) + fabs ( x2 ) + fabs ( x3 ) + fabs ( x4 ) ) / 4.0;

	// a delta y:
	float tempy = ( fabs ( y2 - y4 ) + fabs ( y1 - y3 ) ) / 2.0;
	c2 = tempy / tempx;

	streamlog_out ( MESSAGE4 ) << "Coefficients chip " << ichip << ": c1: " << c1 << " c2: " << c2 << endl;

	if ( _coefffile == true )
	{
	    ofstream filterFile;
	    filterFile.open ( _filterFileName.c_str ( ) );
	    filterFile << c1 << endl;
	    filterFile << c2 << endl;
	    filterFile.close ( );
	    streamlog_out ( MESSAGE4 ) << "Coefficients written!" << endl;
	}
    }
}

// fills histograms with the 16 headers and the first channel
void AlibavaHeader::fillHistos ( TrackerDataImpl * headerdata, TrackerDataImpl * channeldata, int ichip )
{

    FloatVec headvec = headerdata -> getChargeValues ( );
    FloatVec chanvec = channeldata -> getChargeValues ( );

    for ( int ichan = 0; ichan < 16; ichan++ )
    {
	char tmpchar[100];
	sprintf ( tmpchar, "Chip %d, Header %d", ichip, ichan );
	if ( TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap[tmpchar] ) )
	{
	    histo -> Fill ( headvec[ichan] );
	}
    }
    char tmpchar[100];
    sprintf ( tmpchar, "Chip %d, Channel 0", ichip );
    if ( TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap[tmpchar] ) )
    {
	histo -> Fill ( chanvec[0] );
    }
}

// books histos
void AlibavaHeader::bookHistos ( )
{
    AIDAProcessor::tree ( this ) -> cd ( this -> name ( ) );
    for ( unsigned int ichip = 0; ichip < 2; ichip++ )
    {

	for ( int ichan = 0; ichan < 16; ichan++ )
	{
	    char tmpchar[100];
	    sprintf ( tmpchar, "Chip %d, Header %d", ichip, ichan );
	    TH1D * chanDataHisto = new TH1D ( tmpchar, "", 1000, 0, 1000 );
	    _rootObjectMap.insert ( make_pair ( tmpchar, chanDataHisto ) );
	    chanDataHisto -> SetTitle ( tmpchar );
	    chanDataHisto -> SetXTitle ( "ADCs" );
	    chanDataHisto -> SetYTitle ( "Entires" );

	    sprintf ( tmpchar, "Chip %d, Fit %d", ichip, ichan );
	    TF1 *chanDataFit = new TF1 ( tmpchar, "gaus" );
	    _rootObjectMap.insert ( make_pair ( tmpchar, chanDataFit ) );
	}

	char tmpchar[100];
	sprintf ( tmpchar, "Chip %d, Channel 0", ichip );
	TH1D * chanDataHisto = new TH1D ( tmpchar, "", 1000, 0, 1000 );
	_rootObjectMap.insert ( make_pair ( tmpchar, chanDataHisto ) );
	chanDataHisto -> SetTitle ( tmpchar );
	chanDataHisto -> SetXTitle ( "ADCs" );
	chanDataHisto -> SetXTitle ( "Entires" );

	sprintf ( tmpchar, "Chip %d, Channel Fit 0", ichip );
	TF1 * chanDataFit = new TF1 ( tmpchar, "gaus" );
	_rootObjectMap.insert ( make_pair ( "Channel1Fit", chanDataFit ) );

	sprintf ( tmpchar, "Chip %d, Correlation To Channel 0", ichip );
	TH2D * correlationHisto = new TH2D ( tmpchar, "", 1000, 0, 1000, 1000, 0, 1000 );
	_rootObjectMap.insert ( make_pair ( tmpchar, correlationHisto ) );
	correlationHisto -> SetTitle ( tmpchar );
	correlationHisto -> SetXTitle ( "Last Header [raw ADCs]" );
	correlationHisto -> SetYTitle ( "First Channel [raw ADCs]" );

	sprintf ( tmpchar, "Chip %d, Low Profile", ichip );
	TProfile * lowprofile = new TProfile ( tmpchar, "", 100, 0, 1000, 0, 1000, "s" );
	_rootObjectMap.insert ( make_pair ( tmpchar, lowprofile ) );
	lowprofile -> SetTitle ( tmpchar );
	lowprofile -> SetXTitle ( "Last Header [raw ADCs]" );
	lowprofile -> SetYTitle ( "First Channel [raw ADCs]" );

	sprintf ( tmpchar, "Chip %d, High Profile", ichip );
	TProfile * highprofile = new TProfile ( tmpchar, "", 100, 0, 1000, 0, 1000, "s" );
	_rootObjectMap.insert ( make_pair ( tmpchar, highprofile ) );
	highprofile -> SetTitle ( tmpchar );
	highprofile -> SetXTitle ( "Last Header [raw ADCs]" );
	highprofile -> SetYTitle ( "First Channel [raw ADCs]" );

	sprintf ( tmpchar, "Chip %d, Correlation Last Headers", ichip );
	TH2D * correlation2Histo = new TH2D ( tmpchar, "", 1000, 0, 1000,1000,0,1000);
	_rootObjectMap.insert ( make_pair ( tmpchar, correlation2Histo ) );
	correlation2Histo -> SetTitle ( tmpchar );
	correlation2Histo -> SetXTitle ( "Last Header [raw ADCs]" );
	correlation2Histo -> SetYTitle ( "Last but one Header [raw ADCs]" );

	sprintf ( tmpchar, "Chip %d, Low Low Signals", ichip );
	TH2D * lowlowhisto = new TH2D ( tmpchar, "", 1000, 0, 1000, 1000, 0, 1000 );
	_rootObjectMap.insert ( make_pair ( tmpchar, lowlowhisto ) );
	lowlowhisto -> SetTitle ( tmpchar );
	lowlowhisto -> SetXTitle ( "Last Header [raw ADCs]" );
	lowlowhisto -> SetYTitle ( "First Channel [raw ADCs]" );

	sprintf ( tmpchar, "Chip %d, Low High Signals", ichip );
	TH2D * lowhighhisto = new TH2D ( tmpchar, "", 1000, 0, 1000, 1000, 0, 1000 );
	_rootObjectMap.insert ( make_pair ( tmpchar, lowhighhisto ) );
	lowhighhisto -> SetTitle ( tmpchar );
	lowhighhisto -> SetXTitle ( "Last Header [raw ADCs]" );
	lowhighhisto -> SetYTitle ( "First Channel [raw ADCs]" );

	sprintf ( tmpchar, "Chip %d, High Low Signals", ichip );
	TH2D * highlowhisto = new TH2D ( tmpchar, "", 1000, 0, 1000, 1000, 0, 1000 );
	_rootObjectMap.insert ( make_pair ( tmpchar, highlowhisto ) );
	highlowhisto -> SetTitle ( tmpchar );
	highlowhisto -> SetXTitle ( "Last Header [raw ADCs]" );
	highlowhisto -> SetYTitle ( "First Channel [raw ADCs]" );

	sprintf ( tmpchar, "Chip %d, High High Signals", ichip );
	TH2D * highhighhisto = new TH2D ( tmpchar, "", 1000, 0, 1000, 1000, 0, 1000 );
	_rootObjectMap.insert ( make_pair ( tmpchar, highhighhisto ) );
	highhighhisto -> SetTitle ( tmpchar );
	highhighhisto -> SetXTitle ( "Last Header [raw ADCs]" );
	highhighhisto -> SetYTitle ( "First Channel [raw ADCs]" );

    }

    streamlog_out ( MESSAGE4 )  << "End of Booking histograms. " << endl;

}

void AlibavaHeader::correlateLastHeader ( TrackerDataImpl * headerdata, TrackerDataImpl * channeldata, unsigned int ichip )
{
    FloatVec headvec = headerdata -> getChargeValues ( );
    FloatVec chanvec = channeldata -> getChargeValues ( );

    double header1 = headvec[14];
    double header = headvec[15];
    double channel = chanvec[0];

    char tmpchar[100];

    sprintf ( tmpchar, "Chip %d, Correlation To Channel 0", ichip );
    if ( TH2D * histo = dynamic_cast < TH2D* > ( _rootObjectMap[tmpchar] ) )
    {
	histo -> Fill ( header, channel );
    }
    sprintf ( tmpchar, "Chip %d, Low Profile", ichip );
    if ( TProfile* profile = dynamic_cast < TProfile* > ( _rootObjectMap[tmpchar] ) )
    {
	if ( header < 500 )
	{
	    profile -> Fill ( header, channel, 1.0 );
	}
    }
    sprintf ( tmpchar, "Chip %d, High Profile", ichip );
    if ( TProfile* profile = dynamic_cast < TProfile* > ( _rootObjectMap[tmpchar] ) )
    {
	if ( header > 500 )
	{
	    profile -> Fill ( header, channel, 1.0 );
	}
    }
    sprintf ( tmpchar, "Chip %d, Correlation Last Headers", ichip );
    if ( TH2D * histo = dynamic_cast < TH2D* > ( _rootObjectMap[tmpchar] ) )
    {
	histo -> Fill ( header, header1 );
    }

    if ( header1 < 500 )
    {
	if ( header < 500 )
	{
	    sprintf ( tmpchar, "Chip %d, Low Low Signals", ichip );
	    if ( TH2D * histo = dynamic_cast < TH2D* > ( _rootObjectMap[tmpchar]) )
	    {
		histo -> Fill ( header, channel );
	    }
	}
	if ( header > 500 )
	{
	    sprintf ( tmpchar, "Chip %d, Low High Signals", ichip );
	    if ( TH2D * histo = dynamic_cast < TH2D* > ( _rootObjectMap[tmpchar] ) )
	    {
		histo -> Fill ( header, channel );
	    }
	}
    }

    if ( header1 > 500 )
    {
	if ( header < 500 )
	{
	    sprintf ( tmpchar, "Chip %d, High Low Signals", ichip );
	    if ( TH2D * histo = dynamic_cast < TH2D* > ( _rootObjectMap[tmpchar] ) )
	    {
		histo -> Fill ( header, channel );
	    }
	}
	if ( header > 500 )
	{
	    sprintf ( tmpchar, "Chip %d, High High Signals", ichip );
	    if ( TH2D * histo = dynamic_cast < TH2D* > ( _rootObjectMap[tmpchar] ) )
	    {
		histo -> Fill ( header, channel );
	    }
	}
    }
}
