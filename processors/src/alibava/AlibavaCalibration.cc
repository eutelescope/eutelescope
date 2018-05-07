/*
 * Created by Thomas Eichhorn
 *  (2017 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

// alibava includes ".h"
#include "AlibavaCalibration.h"
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
#include <UTIL/LCTime.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCEventImpl.h>

// ROOT includes ".h"
#include <TH1D.h>
#include <TF1.h>
#include "TROOT.h"
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


AlibavaCalibration::AlibavaCalibration ( ) : AlibavaBaseProcessor ( "AlibavaCalibration" )
{

    _description = "AlibavaCalibration analyses calibration files.";

    registerInputCollection ( LCIO::TRACKERDATA, "InputCollectionName", "Input Collection Name", _inputCollectionName, string ( "rawdata" ) );

}


void AlibavaCalibration::init ( )
{
    streamlog_out ( MESSAGE4 ) << "Running init!" << endl;

    if ( Global::parameters -> isParameterSet ( ALIBAVA::CHANNELSTOBEUSED ) )
    {
	Global::parameters -> getStringVals ( ALIBAVA::CHANNELSTOBEUSED, _channelsToBeUsed );
    }
    else
    {
	streamlog_out ( MESSAGE4 ) << "The global parameter "<< ALIBAVA::CHANNELSTOBEUSED <<" is not set!" << endl;
    }

    printParameters ( );

    _pol = -1;

    _nccalpoints = 0;

    _prevcal = 0.0;
}


void AlibavaCalibration::processRunHeader ( LCRunHeader * rdr )
{

    auto arunHeader = std::make_unique < AlibavaRunHeaderImpl > ( rdr );
    arunHeader -> addProcessor ( type ( ) );
    setChipSelection ( arunHeader -> getChipSelection ( ) );
    setChannelsToBeUsed ( );

    bookHistos ( );

}


void AlibavaCalibration::processEvent ( LCEvent * anEvent )
{

    AlibavaEventImpl * alibavaEvent = static_cast < AlibavaEventImpl* > ( anEvent );

    EVENT::IntVec chipSelection = getChipSelection ( );
    LCCollectionVec * collectionVec;
    try
    {
	collectionVec = dynamic_cast < LCCollectionVec * > ( alibavaEvent -> getCollection ( getInputCollectionName ( ) ) );

	for ( unsigned int j = 0; j < chipSelection.size ( ); j++ )
	{
	    unsigned int ichip = chipSelection[j];
	    TrackerDataImpl * trkdata = dynamic_cast < TrackerDataImpl * > ( collectionVec -> getElementAt ( j ) );

	    FloatVec datavec;
	    datavec = trkdata -> getChargeValues ( );
	
	    double calcharge = alibavaEvent -> getCalCharge ( );
	    double caldelay = alibavaEvent -> getCalDelay ( );

	    for ( size_t i = 0; i < datavec.size ( ); i++ )
	    {
		streamlog_out ( DEBUG0 ) << "Chip " << ichip << " channel " << i << " charge " << calcharge << " e, delay " << caldelay << " ns, signal " << datavec.at ( i ) << " ADCs" << endl;
		fillhisto ( ichip, i, calcharge * _pol, caldelay, datavec.at ( i ) );

		// pol catches the alternating polarity of the calibration
		_pol = _pol * -1;

		// count number of calcharge points to determine if this is a charge calibration or delay calibration run
		if ( calcharge != _prevcal )
		{
		    _nccalpoints++;
		}
		_prevcal = calcharge;
	    }
	}
    }
    catch ( ... )
    {
	streamlog_out ( ERROR5 ) << "Collection (" << getInputCollectionName ( ) << ") not found! " << endl;
    }
    _pol = _pol * -1;
}


void AlibavaCalibration::check ( LCEvent * )
{

}


void AlibavaCalibration::end ( )
{

    // arbitray guess
    if ( _nccalpoints > 10 )
    {
	EVENT::IntVec chipSelection = getChipSelection ( );

	AIDAProcessor::tree ( this ) -> cd ( this -> name ( ) );

	for ( unsigned int i = 0; i < chipSelection.size ( ); i++ )
	{
	    unsigned int ichip = chipSelection[i];

	    for ( int ichan = 0; ichan < ALIBAVA::NOOFCHANNELS; ichan++ )
	    {
		char tmpchar[100];
		if ( isMasked ( ichip, ichan ) )
		{
		    continue;
		}

		// seperate fit for positive and negative signals...
		sprintf ( tmpchar, "Charge Calibration Chip %d, Channel %d", ichip, ichan );
		TProfile * histo = dynamic_cast < TProfile* > ( _rootObjectMap[tmpchar] );
		sprintf ( tmpchar, "Charge Calibration Positive Fit Chip %d, Channel %d", ichip, ichan );
		TF1 * posfit = dynamic_cast < TF1* > ( _rootObjectMap[tmpchar] );
		posfit -> SetRange ( 0.0, 1E5 );
		sprintf ( tmpchar, "Charge Calibration Negative Fit Chip %d, Channel %d", ichip, ichan );
		TF1 * negfit = dynamic_cast < TF1* > ( _rootObjectMap[tmpchar] );
		negfit -> SetRange ( -1.0 * 1E5, 0.0 );

		histo -> Fit ( posfit, "QR+" );
		histo -> Fit ( negfit, "QR+" );

		double pos = 0.0;
		double neg = 0.0;

		pos = 1.0 / ( posfit -> GetParameter ( 1 ) );
		neg = 1.0 / ( negfit -> GetParameter ( 1 ) );

		sprintf ( tmpchar, "Charge Calibration, Chip %d, Positive", ichip );
		TH1D * poshisto = dynamic_cast < TH1D* > ( _rootObjectMap[tmpchar] );
		poshisto -> SetBinContent ( ichan, pos );

		sprintf ( tmpchar, "Charge Calibration, Chip %d, Negative", ichip );
		TH1D * neghisto = dynamic_cast < TH1D* > ( _rootObjectMap[tmpchar] );
		neghisto -> SetBinContent ( ichan, neg );

		streamlog_out ( DEBUG2 ) << "Gain chip " << ichip << " channel " << ichan << ": " << pos << " | " << neg << " e / ADC" << endl;
	    }
	}

    }
    else
    {
	streamlog_out ( MESSAGE4 ) << "Found " << _nccalpoints << " charge calibration points! This was a delay calibration run!" << endl;
    }
}


void AlibavaCalibration::fillhisto ( int chip, int chan, double calc, double cald, double q )
{

    char tmpchar[100];
    sprintf ( tmpchar, "Charge Calibration Chip %d, Channel %d", chip, chan );

    if ( TProfile * histo = dynamic_cast < TProfile* > ( _rootObjectMap[tmpchar] ) )
    {
	histo -> Fill ( calc, q );
    }
    sprintf ( tmpchar, "Delay Calibration Chip %d, Channel %d", chip, chan );
    if ( TProfile * histo = dynamic_cast < TProfile* > ( _rootObjectMap[tmpchar] ) )
    {
	histo -> Fill ( cald, q );
    }

}


void AlibavaCalibration::bookHistos ( )
{

    EVENT::IntVec chipSelection = getChipSelection ( );

    AIDAProcessor::tree ( this ) -> cd ( this -> name ( ) );
    //AIDAProcessor::tree ( this ) -> mkdir ( "ChannelData" );

    for ( unsigned int i = 0; i < chipSelection.size ( ); i++ )
    {
	unsigned int ichip = chipSelection[i];
	AIDAProcessor::tree ( this ) -> cd ( this -> name ( ) );
	//AIDAProcessor::tree ( this ) -> cd ( "ChannelData" );

	for ( int ichan = 0; ichan < ALIBAVA::NOOFCHANNELS; ichan++ )
	{
	    char tmpchar[100];
	    if ( isMasked ( ichip, ichan ) )
	    {
		continue;
	    }

	    sprintf ( tmpchar, "Charge Calibration Chip %d, Channel %d", ichip, ichan );
	    TProfile * calchisto = new TProfile ( tmpchar, tmpchar, 1000, -1.0 * 1E5, 1E5 );
	    _rootObjectMap.insert ( make_pair ( tmpchar, calchisto ) );
	    calchisto -> SetTitle ( tmpchar );
	    calchisto -> SetXTitle ( "Injected Charge in e" );
	    calchisto -> SetYTitle ( "Signal in ADCs" );

	    sprintf ( tmpchar, "Charge Calibration Positive Fit Chip %d, Channel %d", ichip, ichan );
	    TF1 *calposFit = new TF1 ( tmpchar, "pol1" );
	    _rootObjectMap.insert ( make_pair ( tmpchar, calposFit ) );

	    sprintf ( tmpchar, "Charge Calibration Negative Fit Chip %d, Channel %d", ichip, ichan );
	    TF1 *calnegFit = new TF1 ( tmpchar, "pol1" );
	    _rootObjectMap.insert ( make_pair ( tmpchar, calnegFit ) );

	}

	for ( int ichan = 0; ichan < ALIBAVA::NOOFCHANNELS; ichan++ )
	{
	    char tmpchar[100];
	    if ( isMasked ( ichip, ichan ) )
	    {
		continue;
	    }

	    sprintf ( tmpchar, "Delay Calibration Chip %d, Channel %d", ichip, ichan );
	    TProfile * caldhisto = new TProfile ( tmpchar, tmpchar, 1000, 0, 256 );
	    _rootObjectMap.insert ( make_pair ( tmpchar, caldhisto ) );
	    caldhisto -> SetTitle ( tmpchar );
	    caldhisto -> SetXTitle ( "Delay in ns" );
	    caldhisto -> SetYTitle ( "Signal in ADCs" );

	}

	AIDAProcessor::tree ( this ) -> cd ( this -> name ( ) );

	char tmpchar[100];

	sprintf ( tmpchar, "Charge Calibration, Chip %d, Positive", ichip );
	TH1D * caliposhisto = new TH1D ( tmpchar, tmpchar, ALIBAVA::NOOFCHANNELS, 0, ( ALIBAVA::NOOFCHANNELS -1 ) );
	_rootObjectMap.insert ( make_pair ( tmpchar, caliposhisto ) );
	caliposhisto -> SetTitle ( tmpchar );
	caliposhisto -> SetXTitle ( "Channel" );
	caliposhisto -> SetYTitle ( "Gain in e^{-} per ADC" );

	sprintf ( tmpchar, "Charge Calibration, Chip %d, Negative", ichip );
	TH1D * calineghisto = new TH1D ( tmpchar, tmpchar, ALIBAVA::NOOFCHANNELS, 0, ( ALIBAVA::NOOFCHANNELS -1 ) );
	_rootObjectMap.insert ( make_pair ( tmpchar, calineghisto ) );
	calineghisto -> SetTitle ( tmpchar );
	calineghisto -> SetXTitle ( "Channel" );
	calineghisto -> SetYTitle ( "Gain in e^{-} per ADC" );

    }

}
