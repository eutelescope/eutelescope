/*
 * Created by Thomas Eichhorn
 *  (2016 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

// alibava includes ".h"
#include "AlibavaSimConverter.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"
#include "marlin/ProcessorMgr.h"

#if defined ( USE_AIDA ) || defined ( MARLIN_USE_AIDA )
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
#include "TRandom.h"

// system includes <>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <memory>
#include <stdlib.h>
#include <algorithm>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;

AlibavaSimConverter::AlibavaSimConverter ( ) : AlibavaBaseProcessor ( "AlibavaSimConverter" )
{

    _description = "AlibavaSimConverter converts Allpix simulated zs lcio data into 'noisy' Alibava lcio data!";

    registerInputCollection ( LCIO::TRACKERDATA, "InputCollectionName", "Header name", _inputCollectionName, string ( "Det350" ) );

    registerProcessorParameter ( "OutputCollectionName", "The collection we output", _outputcollectionname , string ( "rawdata" ) );

    registerProcessorParameter ( "UnsensitiveAxis", "The unsensitive axis of our strip sensor", _nonsensitiveaxis, string ( "x" ) );

    registerOptionalParameter ( "ChipSelection", "Selection of chip that you want to store data from. Chip numbers start from 0. If not set, all data (i.e. chip 0 and 1) will be stored", _chipSelection, EVENT::IntVec ( ) );

    registerOptionalParameter ( "CommonmodeMean", "Mean of the commonmode noise to be added to each strip if CommonmodeSigma > 0", _commonmodemean, 0.0 );

    registerOptionalParameter ( "CommonmodeSigma", "Sigma of the commonmode noise to be added to each strip if > 0", _commonmodesigma, 20.0 );

    registerOptionalParameter ( "NoiseMean", "Mean of the Gaussian noise to be added to each strip if NoiseSigma > 0", _noisemean, 0.0 );

    registerOptionalParameter ( "NoiseSigma", "Sigma of the Gaussian noise to be added to each strip if > 0", _noisesigma, 5.0 );

    registerOptionalParameter ( "PedestalMean", "Mean of the pedestal noise to be added to each strip if PedestalSigma > 0", _pedestalmean, 500.0 );

    registerOptionalParameter ( "PedestalSigma", "Sigma of the pedestal noise to be added to each strip if > 0", _pedestalsigma, 10.0 );

    registerOptionalParameter ( "Ratio", "Scale all signals by this factor", _scalefactor, 1.0 );

    registerOptionalParameter ( "TDCMean", "Mean of the event TDC to be set", _tdcmean, 20.0 );

    registerOptionalParameter ( "TDCSigma", "Sigma of the event TDC to be set", _tdcsigma, 5.0 );

}

void AlibavaSimConverter::init ( )
{
    streamlog_out ( MESSAGE4 ) << "Running init" << endl;

    printParameters ( );

    if ( Global::parameters -> isParameterSet ( ALIBAVA::SKIPMASKEDEVENTS ) )
    {
	_skipMaskedEvents = bool ( Global::parameters -> getIntVal ( ALIBAVA::SKIPMASKEDEVENTS ) );
    }
    else
    {
	streamlog_out ( MESSAGE4 ) << "The Global Parameter " << ALIBAVA::SKIPMASKEDEVENTS << " is not set! Masked events will be used!" << endl;
    }

    checkIfChipSelectionIsValid ( );

    // generate pedestals
    if ( _pedestalsigma > 0.0 )
    {
	for ( int i = 0; i < ALIBAVA::NOOFCHANNELS; i++ )
	{
	    double pedestal = 0.0;
	    pedestal = gRandom -> Gaus ( _pedestalmean, _pedestalsigma );
	    _pedestaldb[i] = pedestal;
	    streamlog_out ( DEBUG2 ) << "Output chanel "<< i << " pedestal is " << pedestal << endl;
	}
    }
}

void AlibavaSimConverter::processRunHeader ( LCRunHeader * rdr )
{
    streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;

    // Add processor name to the runheader
    auto arunHeader = std::make_unique < AlibavaRunHeaderImpl > ( rdr );
    arunHeader -> addProcessor ( type ( ) );

    // set number of skipped events to zero (defined in AlibavaBaseProcessor)
    _numberOfSkippedEvents = 0;

    arunHeader -> setChipSelection ( _chipSelection );
}

void AlibavaSimConverter::processEvent ( LCEvent * anEvent )
{
    if ( anEvent -> getEventNumber ( ) % 1000 == 0 )
    {
	streamlog_out ( MESSAGE4 ) << "Looping events " << anEvent -> getEventNumber ( ) << endl;
    }

    AlibavaEventImpl * alibavaEvent = static_cast < AlibavaEventImpl* > ( anEvent );

    // creating LCCollection for raw data
    LCCollectionVec* rawDataCollection = new LCCollectionVec ( LCIO::TRACKERDATA );

    CellIDEncoder < TrackerDataImpl > chipIDEncoder ( ALIBAVA::ALIBAVADATA_ENCODE, rawDataCollection );

    LCCollectionVec * inputVec;
    try
    {
	inputVec = dynamic_cast < LCCollectionVec * > ( alibavaEvent -> getCollection ( getInputCollectionName ( ) ) ) ;
	TrackerDataImpl * trkdata = dynamic_cast < TrackerDataImpl * > ( inputVec -> getElementAt ( 0 ) ) ;

	// loop both chips for now...
	for ( int k = 0; k < ALIBAVA::NOOFCHIPS; k++ )
	{
	    TrackerDataImpl * newdataImpl = new TrackerDataImpl ( );
	    chipIDEncoder[ALIBAVA::ALIBAVADATA_ENCODE_CHIPNUM] = k;
	    chipIDEncoder.setCellID ( newdataImpl );
	    FloatVec datavec;
	    datavec = trkdata -> getChargeValues ( );

	    FloatVec outputvec;
	    FloatVec pixel_x;
	    FloatVec pixel_y;
	    FloatVec pixel_charge;
	    FloatVec pixel_time;

	    // only sensible if there is a charge in this event
	    streamlog_out ( DEBUG0 ) << "Simulated charges in event " << anEvent -> getEventNumber ( ) << " : " << datavec.size ( ) / 8 << endl;
	    if ( datavec.size ( ) > 0 )
	    {
		// data is x,y,q,t
		for ( size_t i = 0; i <= ( datavec.size ( ) - 4 ); i = i + 4 )
		{
		    streamlog_out ( DEBUG1 ) << "chan x " << datavec[i] << endl;
		    streamlog_out ( DEBUG1 ) << "chan y " << datavec[i + 1] << endl;
		    streamlog_out ( DEBUG1 ) << "chan q " << datavec[i + 2] << endl;
		    streamlog_out ( DEBUG1 ) << "chan t " << datavec[i + 3] << endl;
		    pixel_x.push_back ( datavec[i] );
		    pixel_y.push_back ( datavec[i + 1] );
		    double charge = datavec[i + 2];

		    // possibility to change the charge to 'sensible' numbers
		    if ( _scalefactor != 0 )
		    {
			charge = datavec[i + 2] * _scalefactor;
		    }
		    pixel_charge.push_back ( charge );
		    pixel_time.push_back ( datavec[i + 3] );
		}
	    }

	    // fill the channels with pedestals, noise and common mode (and a signal if it's there...)
	    int point = 0;

	    // common mode noise is the same for all channels in an event, so calculate beforehand
	    double commonmode = 0.0;
	    if ( _commonmodesigma > 0.0 )
	    {
		commonmode = gRandom -> Gaus ( _commonmodemean, _commonmodesigma );
		streamlog_out ( DEBUG2 ) << "Output common mode in this event is " << commonmode << endl;
	    }
	    for ( size_t j = 0; j < ALIBAVA::NOOFCHANNELS; j++ )
	    {

		double background = 0.0;
		// add pedestal
		double pedestal = gRandom -> Gaus ( _pedestaldb[j + ( k * ALIBAVA::NOOFCHANNELS ) ], _pedestalsigma );
		background += pedestal;

		// add the noise
		if ( _noisesigma > 0.0 )
		{
		    double noise = 0.0;
		    noise = gRandom -> Gaus ( _noisemean, _noisesigma );
		    background += noise;
		    streamlog_out ( DEBUG2 ) << "Output chip " << k << " chanel " << j << " noise is " << noise << endl;
		}

		// add the charge, if any
		if ( _nonsensitiveaxis == "x" && pixel_y.size ( ) > 0 )
		{
			if ( pixel_y[point] == ( j + k * ALIBAVA::NOOFCHANNELS ) )
			{
				double signal = 0.0;
				signal = pixel_charge[point];
				background += signal;
				point++;
				streamlog_out ( DEBUG2 ) << "Output chip " << k << " chanel "<< j << " signal is " << signal << endl;
			}
		}
		if ( _nonsensitiveaxis == "y" && pixel_x.size ( ) > 0 )
		{
			if ( pixel_x[point] == ( j + k * ALIBAVA::NOOFCHANNELS ) )
			{
				double signal = 0.0;
				signal = pixel_charge[point];
				background += signal;
				point++;
				streamlog_out ( DEBUG2 ) << "Output chip " << k << " chanel "<< j << " signal is " << signal << endl;
			}
		}

		// now add common mode
		background += commonmode;

		// push back this channel
		outputvec.push_back ( background );
		streamlog_out ( DEBUG3 ) << "Output chip " << k << " chanel "<< j << " final signal is " << background << endl;
	    }

	    // and let there be output
	    newdataImpl -> setChargeValues ( outputvec );
	    rawDataCollection -> push_back ( newdataImpl );
	}

	// set the time of this event
	float thetime = 0.0;
	thetime = gRandom -> Gaus ( _tdcmean, _tdcsigma );
	anEvent -> parameters ( ) .setValue ( "EventTDCTime", thetime );

	alibavaEvent -> addCollection ( rawDataCollection, _outputcollectionname );

    }
    catch ( lcio::DataNotAvailableException& )
    {
	streamlog_out( ERROR5 ) << "Collection (" << getInputCollectionName ( ) << ") not found! " << endl;
    }
}

void AlibavaSimConverter::check ( LCEvent * )
{

}

void AlibavaSimConverter::end ( )
{

}

void AlibavaSimConverter::checkIfChipSelectionIsValid ( )
{

    bool resetChipSelection = false;

    // check if there is chip selection or if there is chip selection but not valid
    if ( _chipSelection.size ( ) == 0 )
    {
	streamlog_out ( WARNING5 ) << "You didn't select any chip" << endl;
	resetChipSelection = true;
    }
    else if ( int ( _chipSelection.size ( ) ) > ALIBAVA::NOOFCHIPS )
    {
	streamlog_out ( WARNING5 ) << "The number of chips you selected ( " << _chipSelection.size ( ) << " ) is more than an alibava daughter board can have ( " << ALIBAVA::NOOFCHIPS << " )" << endl;
	resetChipSelection = true;
    }
    else
    {
	// first sort the chip numbers in ascending order
	sort( _chipSelection.begin ( ) ,_chipSelection.end ( ) );

	// check if the selected chips make sense
	for ( int ichip = 0; ichip < int ( _chipSelection.size ( ) ); ichip++ )
	{
	    bool del_this_chip = false;

	    if ( _chipSelection[ichip] < 0 )
	    {
		streamlog_out ( ERROR5 ) << "Selected chip cannot have negative value." << endl;
		del_this_chip = true;
	    }
	    else if ( _chipSelection[ichip] >= ALIBAVA::NOOFCHIPS )
	    {
		streamlog_out ( ERROR5 ) << "Chip numbering has to start from zero \"0\" and cannot be greater than " << ALIBAVA::NOOFCHIPS - 1 << endl;
		del_this_chip = true;
	    }

	    if ( del_this_chip )
	    {
		// if this chip selection is not valid, delete it.
		streamlog_out ( ERROR5 ) << "Chip " << _chipSelection[ichip] << " is deleted from the chip selection list" << endl;
		_chipSelection.erase ( _chipSelection.begin ( ) + ichip );
		ichip = ichip - 1;
	    }
	}

	// check again if there is any selected chip left
	if ( _chipSelection.size ( ) == 0 )
	{
	    resetChipSelection= true;
	}

    }

    if ( resetChipSelection )
    {
	streamlog_out ( WARNING5 ) << "I will save data from all chips" << endl;

	_chipSelection.clear ( );

	for ( int ichip = 0; ichip < ALIBAVA::NOOFCHIPS; ichip++ )
	{
	    _chipSelection.push_back ( ichip );
	}
    }

    // now there is valid chip selection!!!

    streamlog_out ( WARNING5 ) << "Final applied chip selection: ";
    for ( int ichip = 0; ichip < int ( _chipSelection.size ( ) ); ichip++ )
    {
	streamlog_out ( WARNING5 ) << _chipSelection[ichip] << " ";
    }
    streamlog_out ( WARNING5 ) << endl;
    streamlog_out ( WARNING5 ) << "Only data coming from these chips will be stored" << endl;

}
