/*
 * Created by Thomas Eichhorn
 *  (2018 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

// alibava includes ".h"
#include "AlibavaDataPlotter.h"
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
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ITree.h>
#include <AIDA/IAxis.h>
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
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

// system includes <>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <memory>
#include <cstdlib>
#include <ctime>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;

AlibavaDataPlotter::AlibavaDataPlotter ( ) : AlibavaBaseProcessor ( "AlibavaDataPlotter" )
{
    _description = "AlibavaDataPlotter reads TrackerData of Alibava data and produces histograms";

    registerInputCollection ( LCIO::TRACKERDATA, "InputCollectionName", "Input raw data collection name", _inputCollectionName, string ( "rawdata" ) );

    registerOptionalParameter ( "PedestalInputFile", "The filename where the pedestal and noise values stored", _pedestalFile , string ( ALIBAVA::NOTSET ) );

    registerOptionalParameter ( "CalibrationInputFile", "The filename where the calibration values stored", _calibrationFile , string ( ALIBAVA::NOTSET ) );

    registerOptionalParameter ( "PedestalCollectionName", "Pedestal collection name, better not to change", _pedestalCollectionName, string ( ALIBAVA::NOTSET ) );

    registerOptionalParameter ( "NoiseCollectionName", "Noise collection name, better not to change", _noiseCollectionName, string ( ALIBAVA::NOTSET ) );

    registerOptionalParameter ( "ChargeCalibrationCollectionName", "Charge calibration collection name, better not to change", _chargeCalCollectionName, string ( ALIBAVA::NOTSET ) );

    registerOptionalParameter ( "EventsToPlot", "Choose specific events to plot.", _plotEvents, IntVec ( ) );

    registerOptionalParameter ( "MultiplySignalBy", "In case this variable is set, all signals will be multipled by this value.", _multiplySignalby, float ( 1.0 ) );

    registerOptionalParameter ( "PlotSomePercentOfEvents", "In case this variable is set (say x), x percent of all events will be plotted randomly. The number should be between 0 and 100", _plotXPercentOfEvents, float ( 0.0 ) );

}

void AlibavaDataPlotter::init ( )
{
    streamlog_out ( MESSAGE4 ) << "Running init" << endl;

    printParameters ( );

    if ( Global::parameters -> isParameterSet ( ALIBAVA::CHANNELSTOBEUSED ) )
    {
	Global::parameters -> getStringVals ( ALIBAVA::CHANNELSTOBEUSED, _channelsToBeUsed );
    }
    else
    {
	streamlog_out ( MESSAGE4 ) << "The Global Parameter " << ALIBAVA::CHANNELSTOBEUSED << " is not set!" << endl;
    }

    if ( Global::parameters -> isParameterSet ( ALIBAVA::SKIPMASKEDEVENTS ) )
    {
	_skipMaskedEvents = bool ( Global::parameters -> getIntVal ( ALIBAVA::SKIPMASKEDEVENTS ) );
    }
    else
    {
	streamlog_out ( MESSAGE4 ) << "The Global Parameter " << ALIBAVA::SKIPMASKEDEVENTS << " is not set! Masked events will be used!" << endl;
    }

}

void AlibavaDataPlotter::processRunHeader ( LCRunHeader * rdr )
{
    streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;

    // Add processor name to the runheader
    auto arunHeader = std::make_unique < AlibavaRunHeaderImpl > ( rdr );
    arunHeader -> addProcessor ( type ( ) );

    // get and set selected chips
    setChipSelection ( arunHeader -> getChipSelection ( ) );

    // set channels to be used (if it is defined)
    setChannelsToBeUsed ( );

    setPedestalCollectionName ( _pedestalCollectionName );
    setNoiseCollectionName ( _noiseCollectionName );

    // set pedestal and noise values (if it is defined)
    setPedestals ( );

    bookHistos ( );

    // set number of skipped events to zero (defined in AlibavaBaseProcessor)
    _numberOfSkippedEvents = 0;

    // plot noise and pedestals

    AIDAProcessor::tree ( this ) -> cd ( this -> name ( ) );
    EVENT::IntVec chipSelection = getChipSelection ( );

    // The chipSelection should be in ascending order!

    AIDAProcessor::tree ( this ) -> mkdir ( "PedestalAndNoise" );
    AIDAProcessor::tree ( this ) -> cd ( "PedestalAndNoise" );

    for ( unsigned int i = 0; i < chipSelection.size ( ); i++ )
    {
	int ichip = chipSelection[i];

	// Pedestal
	if ( _pedestalCollectionName != string ( ALIBAVA::NOTSET ) )
	{
	    stringstream s;
	    s << "Pedestal_Chip_" << ichip;
	    TH1D * pedestalHisto = new TH1D ( s.str ( ) .c_str ( ), "", ALIBAVA::NOOFCHANNELS, -0.5, ALIBAVA::NOOFCHANNELS - 0.5 );

	    // Title string for pedestal histogram
	    stringstream sp;
	    sp << "Pedestal Chip " << ichip << ";Channel Number;Pedestal [ADCs]";
	    pedestalHisto -> SetTitle ( ( sp.str( ) ) .c_str ( ) );

	    FloatVec pedVec = getPedestalOfChip ( ichip );

	    for ( int ichan = 0; ichan < ALIBAVA::NOOFCHANNELS; ichan++ )
	    {
		// If channel is masked, do not fill histo
		if ( isMasked ( ichip, ichan ) )
		{
		    continue;
		}
		else
		{
		    pedestalHisto -> SetBinContent ( ichan + 1, pedVec[ichan] );
		}
	    }
	    _rootObjectMap.insert ( make_pair ( s.str ( ), pedestalHisto ) );
	}

	//noise
	if ( _noiseCollectionName != string ( ALIBAVA::NOTSET ) )
	{
	    stringstream s;
	    s << "Noise_Chip_" << ichip;
	    TH1D * noiseHisto = new TH1D ( s.str ( ) .c_str ( ), "", ALIBAVA::NOOFCHANNELS, -0.5, ALIBAVA::NOOFCHANNELS - 0.5 );

	    // Title string for noise histogram
	    stringstream sn;
	    sn << "Noise Chip " << ichip << ";Channel Number;Pedestal [ADCs]";
	    noiseHisto -> SetTitle ( ( sn.str ( ) ) .c_str ( ) );

	    FloatVec noiVec = getNoiseOfChip ( ichip );
	    for ( int ichan = 0; ichan < ALIBAVA::NOOFCHANNELS; ichan++ )
	    {
		// If channel is masked, do not fill histo
		if ( isMasked ( ichip, ichan ) )
		{
		    continue;
		}
		else
		{
		    noiseHisto -> SetBinContent ( ichan + 1, noiVec[ichan] );
		}
	    }

	    _rootObjectMap.insert ( make_pair ( s.str ( ), noiseHisto ) );
	}
    } // End of loop over selected chips
}



void AlibavaDataPlotter::processEvent ( LCEvent * anEvent )
{

    AlibavaEventImpl * alibavaEvent = static_cast < AlibavaEventImpl* > ( anEvent );
    int eventnum = alibavaEvent -> getEventNumber ( );

    if ( _skipMaskedEvents && ( alibavaEvent -> isEventMasked ( ) ) )
    {
	_numberOfSkippedEvents++;
	return;
    }

    float tdctime = alibavaEvent -> getEventTime ( );
    float temperature = alibavaEvent -> getEventTemp ( );

    // TDC time
    TH1D * alibavaTDCTime = dynamic_cast < TH1D* > ( _rootObjectMap["alibavaTDCTime"] );
    alibavaTDCTime -> Fill ( tdctime );
    TH2D * alibavaTDCTimeEvents = dynamic_cast < TH2D* > ( _rootObjectMap["alibavaTDCTimeEvents"] );
    alibavaTDCTimeEvents -> Fill ( eventnum, tdctime );

    // Temperature
    TH1D * alibavaEventTemperature = dynamic_cast < TH1D* > ( _rootObjectMap["alibavaEventTemperature"] );
    alibavaEventTemperature -> Fill ( temperature );
    TH2D * alibavaEventTemperatureEvents = dynamic_cast < TH2D* > ( _rootObjectMap["alibavaEventTemperatureEvents"] );
    alibavaEventTemperatureEvents -> Fill ( eventnum, temperature );

    // Calibration
    TH1D * alibavaCalCharges = dynamic_cast < TH1D* > ( _rootObjectMap["alibavaCalCharges"] );
    alibavaCalCharges -> Fill ( alibavaEvent -> getCalCharge ( ) );

    // Delay
    TH1D * alibavaDelayValues = dynamic_cast < TH1D* > ( _rootObjectMap["alibavaDelayValues"] );
    alibavaDelayValues -> Fill ( alibavaEvent -> getCalDelay ( ) );

    bool plotThisEvent = false;
    if ( isEventToBePlotted ( eventnum ) )
    {
	bookEventHisto ( eventnum );
	plotThisEvent = true;
    }

    // Now loop over detectors
    LCCollectionVec * collectionVec;
    unsigned int noOfChip;
    try
    {
	collectionVec = dynamic_cast < LCCollectionVec * > ( alibavaEvent -> getCollection ( _inputCollectionName ) );
	noOfChip = collectionVec -> getNumberOfElements ( );

	for ( size_t i = 0; i < noOfChip; ++i )
	{
	    TrackerDataImpl * trkdata = dynamic_cast < TrackerDataImpl * > ( collectionVec -> getElementAt ( i ) );

	    // if this event is selected to be plotted, fill the histogram
	    if ( plotThisEvent )
	    {
		fillEventHisto ( eventnum, trkdata );
	    }

	    // for everything else
	    fillOtherHistos ( trkdata, tdctime, temperature );
	}

    }
    catch ( lcio::DataNotAvailableException& )
    {
	// do nothing again
	streamlog_out ( ERROR5 ) << "Collection (" << _inputCollectionName << ") not found!" << endl;
    }

}

void AlibavaDataPlotter::bookHistos ( )
{

    #if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

    try
    {
	streamlog_out ( MESSAGE2 ) << "Booking histograms..." << endl;
	AIDAProcessor::tree ( this ) -> mkdir ( "Events" );

	AIDAProcessor::tree ( this ) -> mkdir ( "EventData" );
	AIDAProcessor::tree ( this ) -> cd ( "EventData" );

	TH1D * alibavaTDCTime = new TH1D ( "alibavaTDCTime", "", 100, 0, 99 );
	alibavaTDCTime -> SetTitle ( "Event TDC Time;TDC Time [nS];Entries" );
	_rootObjectMap.insert ( make_pair ( "alibavaTDCTime", alibavaTDCTime ) );

	TH2D * alibavaTDCTimeEvents = new TH2D ( "alibavaTDCTimeEvents", "", 1000, 0, 10000, 100, 0, 99 );
	alibavaTDCTimeEvents -> SetTitle ( "TDC Time over Events;Event Nr;TDC Time [ns]" );
	_rootObjectMap.insert ( make_pair ( "alibavaTDCTimeEvents", alibavaTDCTimeEvents ) );

	TH1D * alibavaEventTemperature = new TH1D ( "alibavaEventTemperature", "", 150, -50, 99 );
	alibavaEventTemperature -> SetTitle ( "Event Temperature;Temperature [#circC];Entries" );
	_rootObjectMap.insert ( make_pair ( "alibavaEventTemperature", alibavaEventTemperature ) );

	TH2D * alibavaEventTemperatureEvents = new TH2D ( "alibavaEventTemperatureEvents", "", 1000, 0, 10000, 150, -50, 99 );
	alibavaEventTemperatureEvents -> SetTitle ( "Temperature over Events;Event Nr;Temperature [#circC]" );
	_rootObjectMap.insert ( make_pair ( "alibavaEventTemperatureEvents", alibavaEventTemperatureEvents ) );

	TH1D * alibavaCalCharges = new TH1D ( "alibavaCalCharges", "", 1000, 0, 100000 );
	alibavaCalCharges -> SetTitle ( "Calibration Charge Values;Charge [e];Entries" );
	_rootObjectMap.insert ( make_pair ( "alibavaCalCharges", alibavaCalCharges ) );

	TH1D * alibavaDelayValues = new TH1D ( "alibavaDelayValues", "", 251, 0, 250 );
	alibavaDelayValues -> SetTitle ( "Calibration Delay Values;Delay [ns];Entries" );
	_rootObjectMap.insert ( make_pair ( "alibavaDelayValues", alibavaDelayValues ) );

	TH1D * alibavaSignalChip0 = new TH1D ( "alibavaSignalChip0", "", 401, -200, 200 );
	alibavaSignalChip0 -> SetTitle ( "Chip 0 Signal;Signal [ADCs];Entries" );
	_rootObjectMap.insert ( make_pair ( "alibavaSignalChip0", alibavaSignalChip0 ) );

	TH1D * alibavaSignalChip1 = new TH1D ( "alibavaSignalChip1", "", 401, -200, 200 );
	alibavaSignalChip1 -> SetTitle ( "Chip 1 Signal;Signal [ADCs];Entries" );
	_rootObjectMap.insert ( make_pair ( "alibavaSignalChip1", alibavaSignalChip1 ) );

	TH2D * alibavaSignalTDCChip0 = new TH2D ( "alibavaSignalTDCChip0", "", 100, 0, 99, 401, -200, 200 );
	alibavaSignalTDCChip0 -> SetTitle ( "Chip 0 Signal vs TDC Time;Time [ns];Signal [ADCs]" );
	_rootObjectMap.insert ( make_pair ( "alibavaSignalTDCChip0", alibavaSignalTDCChip0 ) );

	TH2D * alibavaSignalTDCChip1 = new TH2D ( "alibavaSignalTDCChip1", "", 100, 0, 99, 401, -200, 200 );
	alibavaSignalTDCChip1 -> SetTitle ( "Chip 1 Signal vs TDC Time;Time [ns];Signal [ADCs]" );
	_rootObjectMap.insert ( make_pair ( "alibavaSignalTDCChip1", alibavaSignalTDCChip1 ) );

	TH2D * alibavaSignalTempChip0 = new TH2D ( "alibavaSignalTempChip0", "", 150, -50, 99, 401, -200, 200 );
	alibavaSignalTempChip0 -> SetTitle ( "Chip 0 Signal vs Temperature;Temperature [#circC];Signal [ADCs]" );
	_rootObjectMap.insert ( make_pair ( "alibavaSignalTempChip0", alibavaSignalTempChip0 ) );

	TH2D * alibavaSignalTempChip1 = new TH2D ( "alibavaSignalTempChip1", "", 150, -50, 99, 401, -200, 200 );
	alibavaSignalTempChip1 -> SetTitle ( "Chip 1 Signal vs Temperature;Temperature [#circC];Signal [ADCs]" );
	_rootObjectMap.insert ( make_pair ( "alibavaSignalTempChip1", alibavaSignalTempChip1 ) );

	TH1D * alibavaSNRChip0 = new TH1D ( "alibavaSNRChip0", "", 401, -200, 200 );
	alibavaSNRChip0 -> SetTitle ( "Chip 0 SNR;SNR;Entries" );
	_rootObjectMap.insert ( make_pair ( "alibavaSNRChip0", alibavaSNRChip0 ) );

	TH1D * alibavaSNRChip1 = new TH1D ( "alibavaSNRChip1", "", 401, -200, 200 );
	alibavaSNRChip1 -> SetTitle ( "Chip 1 SNR;SNR;Entries" );
	_rootObjectMap.insert ( make_pair ( "alibavaSNRChip1", alibavaSNRChip1 ) );

	TH2D * alibavaSNRTDCChip0 = new TH2D ( "alibavaSNRTDCChip0", "", 100, 0, 99, 401, -200, 200 );
	alibavaSNRTDCChip0 -> SetTitle ( "Chip 0 SNR vs TDC Time;Time [ns];SNR" );
	_rootObjectMap.insert ( make_pair ( "alibavaSNRTDCChip0", alibavaSNRTDCChip0 ) );

	TH2D * alibavaSNRTDCChip1 = new TH2D ( "alibavaSNRTDCChip1", "", 100, 0, 99, 401, -200, 200 );
	alibavaSNRTDCChip1 -> SetTitle ( "Chip 1 SNR vs TDC Time;Time [ns];SNR" );
	_rootObjectMap.insert ( make_pair ( "alibavaSNRTDCChip1", alibavaSNRTDCChip1 ) );

	TH2D * alibavaSNRTempChip0 = new TH2D ( "alibavaSNRTempChip0", "", 150, -50, 99, 401, -200, 200 );
	alibavaSNRTempChip0 -> SetTitle ( "Chip 0 SNR vs Temperature;Temperature [#circC];SNR" );
	_rootObjectMap.insert ( make_pair ( "alibavaSNRTempChip0", alibavaSNRTempChip0 ) );

	TH2D * alibavaSNRTempChip1 = new TH2D ( "alibavaSNRTempChip1", "", 150, -50, 99, 401, -200, 200 );
	alibavaSNRTempChip1 -> SetTitle ( "Chip 1 SNR vs Temperature;Temperature [#circC];SNR" );
	_rootObjectMap.insert ( make_pair ( "alibavaSNRTempChip1", alibavaSNRTempChip1 ) );

    }
    catch ( ... )
    {
	
    }

    #endif
}

void AlibavaDataPlotter::check ( LCEvent * /* evt */ )
{

}

void AlibavaDataPlotter::end ( )
{
    if ( _numberOfSkippedEvents > 0 )
    {
	streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents << " events skipped since they are masked" << endl;
    }
    streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
}

void AlibavaDataPlotter::fillHistos ( TrackerDataImpl * /* trkdata */ )
{

}

void AlibavaDataPlotter::fillOtherHistos ( TrackerDataImpl * trkdata, float tdctime, float temperature )
{

    FloatVec datavec;
    datavec = trkdata -> getChargeValues ( );
    int ichip = getChipNum ( trkdata );

    FloatVec noiseVec;
    if ( isNoiseValid ( ) )
    {
	noiseVec = getNoiseOfChip ( ichip );
    }

    for ( int ichan = 0; ichan < ALIBAVA::NOOFCHANNELS; ichan++ )
    {
	// if channel is masked, do not fill histo
	if ( isMasked ( ichip, ichan ) )
	{
	    continue;
	}

	float data = _multiplySignalby * datavec[ichan];

	if ( ichip == 0 )
	{
	    TH1D * alibavaSignalChip0 = dynamic_cast < TH1D* > ( _rootObjectMap["alibavaSignalChip0"] );
	    alibavaSignalChip0 -> Fill ( data );
	    TH2D * alibavaSignalTDCChip0 = dynamic_cast < TH2D* > ( _rootObjectMap["alibavaSignalTDCChip0"] );
	    alibavaSignalTDCChip0 -> Fill ( tdctime, data );
	    TH2D * alibavaSignalTempChip0 = dynamic_cast < TH2D* > ( _rootObjectMap["alibavaSignalTempChip0"] );
	    alibavaSignalTempChip0 -> Fill ( temperature, data );

	    if ( isNoiseValid ( ) )
	    {
		float noise = noiseVec[ichan];
		if ( noise != 0 )
		{
		    TH1D * alibavaSNRChip0 = dynamic_cast < TH1D* > ( _rootObjectMap["alibavaSNRChip0"] );
		    alibavaSNRChip0 -> Fill ( data / noise );
		    TH2D * alibavaSNRTDCChip0 = dynamic_cast < TH2D* > ( _rootObjectMap["alibavaSNRTDCChip0"] );
		    alibavaSNRTDCChip0 -> Fill ( tdctime, data / noise );
		    TH2D * alibavaSNRTempChip0 = dynamic_cast < TH2D* > ( _rootObjectMap["alibavaSNRTempChip0"] );
		    alibavaSNRTempChip0 -> Fill ( temperature, data / noise );
		}
	    }
	}
	else if ( ichip == 1 )
	{
	    TH1D * alibavaSignalChip1 = dynamic_cast < TH1D* > ( _rootObjectMap["alibavaSignalChip1"] );
	    alibavaSignalChip1 -> Fill ( data );
	    TH2D * alibavaSignalTDCChip1 = dynamic_cast < TH2D* > ( _rootObjectMap["alibavaSignalTDCChip1"] );
	    alibavaSignalTDCChip1 -> Fill ( tdctime, data );
	    TH2D * alibavaSignalTempChip1 = dynamic_cast < TH2D* > ( _rootObjectMap["alibavaSignalTempChip1"] );
	    alibavaSignalTempChip1 -> Fill ( temperature, data );

	    if ( isNoiseValid ( ) )
	    {
		float noise = noiseVec[ichan];
		if ( noise != 0 )
		{
		    TH1D * alibavaSNRChip1 = dynamic_cast < TH1D* > ( _rootObjectMap["alibavaSNRChip1"] );
		    alibavaSNRChip1 -> Fill ( data / noise );
		    TH2D * alibavaSNRTDCChip1 = dynamic_cast < TH2D* > ( _rootObjectMap["alibavaSNRTDCChip1"] );
		    alibavaSNRTDCChip1 -> Fill ( tdctime, data / noise );
		    TH2D * alibavaSNRTempChip1 = dynamic_cast < TH2D* > ( _rootObjectMap["alibavaSNRTempChip1"] );
		    alibavaSNRTempChip1 -> Fill ( temperature, data / noise );
		}
	    }
	}

    }

}

void AlibavaDataPlotter::bookEventHisto ( int eventnum )
{
    AIDAProcessor::tree ( this ) -> cd ( this -> name ( ) );
    AIDAProcessor::tree ( this ) -> cd ( "Events" );

    EVENT::IntVec chipSelection = getChipSelection ( );

    for ( unsigned int i = 0; i < chipSelection.size ( ); i++ )
    {
	int ichip = chipSelection[i];
	TH1D * eventHisto = new TH1D ( getEventHistoName ( eventnum, ichip ) .c_str ( ), "", ALIBAVA::NOOFCHANNELS, -0.5, ALIBAVA::NOOFCHANNELS - 0.5 );
	stringstream sp;
	sp << "Event " << eventnum << " (chip " << ichip << ");Channel Number;Signal (ADCs)";
	eventHisto -> SetTitle ( ( sp.str ( ) ) .c_str ( ) );
	_rootObjectMap.insert ( make_pair ( getEventHistoName ( eventnum, ichip ), eventHisto ) );
    }
}

void AlibavaDataPlotter::fillEventHisto ( int eventnum, TrackerDataImpl * trkdata )
{
    streamlog_out ( DEBUG1 ) << "Plotting Event " << eventnum << endl;
    FloatVec datavec;
    datavec = trkdata -> getChargeValues ( );
    int ichip = getChipNum ( trkdata );
    TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap[getEventHistoName ( eventnum, ichip ) ] );

    for ( int ichan = 0; ichan < ALIBAVA::NOOFCHANNELS; ichan++)
    {
	// if channel is masked, do not fill histo
	if ( isMasked ( ichip, ichan ) )
	{
	    continue;
	}
	else
	{
	    histo -> SetBinContent ( ichan + 1, datavec[ichan] );
	}
    }
}

bool AlibavaDataPlotter::isEventToBePlotted ( int eventnum )
{
    for ( unsigned int iEvent = 0; iEvent < _plotEvents.size ( ); iEvent++ )
    {
	if ( _plotEvents[iEvent] == eventnum )
	{
	    return true;
	}

	// Generate a float number between 0-100 and
	float randnum = float ( rand ( ) % ( 100 * 10000 ) ) / 10000.0;
	if ( randnum <= _plotXPercentOfEvents )
	{
	    return true;
	}
    }
    // If event num is not in _plotEvents vector return false
    return false;
}

string AlibavaDataPlotter::getEventHistoName ( int eventnum,  int ichip )
{
    stringstream s;
    s << "Event_" << eventnum << "_chip" << ichip;
    return s.str ( );
}
