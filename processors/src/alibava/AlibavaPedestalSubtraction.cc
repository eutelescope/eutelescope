/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 *
 *  modified by: Thomas Eichhorn thomas.eichhorn@desy.de
 */

// alibava includes ".h"
#include "AlibavaPedestalSubtraction.h"
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
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCEventImpl.h>

// system includes <>
#include <string>
#include <iostream>
#include <memory>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;

AlibavaPedestalSubtraction::AlibavaPedestalSubtraction ( ) : AlibavaBaseProcessor ( "AlibavaPedestalSubtraction" )
{
    // modify processor description
    _description = "AlibavaPedestalSubtraction subtracts the provided pedestal values from the input raw data.";

    // first register the input collection
    registerInputCollection ( LCIO::TRACKERDATA, "InputCollectionName", "Input raw data collection name", _inputCollectionName, string ( "rawdata" ) );

    registerOutputCollection ( LCIO::TRACKERDATA, "OutputCollectionName", "Output data collection name", _outputCollectionName, string ( "recodata" ) );

    registerProcessorParameter ( "PedestalInputFile", "The filename where the pedestal and noise values are stored", _pedestalFile, string ( "pedestal.slcio" ) );

    // now the optional parameters
    registerProcessorParameter ( "PedestalCollectionName", "Pedestal collection name, better not to change", _pedestalCollectionName, string ( "pedestal" ) );

    registerProcessorParameter ( "NoiseCollectionName", "Noise collection name, better not to change", _noiseCollectionName, string ( "noise" ) );
}

void AlibavaPedestalSubtraction::init ( )
{
    streamlog_out ( MESSAGE4 ) << "Running init" << endl;

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

    printParameters ( );

}

void AlibavaPedestalSubtraction::processRunHeader ( LCRunHeader * rdr )
{
    streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;

    // Add processor name to the runheader
    auto arunHeader = std::make_unique < AlibavaRunHeaderImpl > ( rdr );
    arunHeader -> addProcessor ( type ( ) );

    // get and set selected chips
    setChipSelection ( arunHeader -> getChipSelection ( ) );

    // set channels to be used (if it is defined)
    setChannelsToBeUsed ( );

    // set pedestal and noise values
    setPedestals ( );

    // if you want
    bookHistos ( );

    // set number of skipped events to zero (defined in AlibavaBaseProcessor)
    _numberOfSkippedEvents = 0;
}

void AlibavaPedestalSubtraction::processEvent ( LCEvent * anEvent )
{
    AlibavaEventImpl * alibavaEvent = static_cast < AlibavaEventImpl* > ( anEvent );

    if ( _skipMaskedEvents && ( alibavaEvent -> isEventMasked ( ) ) )
    {
	_numberOfSkippedEvents++;
	return;
    }

    LCCollectionVec * collectionVec;
    LCCollectionVec * newDataCollection = new LCCollectionVec ( LCIO::TRACKERDATA );
    CellIDEncoder < TrackerDataImpl > chipIDEncoder ( ALIBAVA::ALIBAVADATA_ENCODE, newDataCollection );

    int noOfChip, chipnum;
    try
    {
	collectionVec = dynamic_cast < LCCollectionVec * > ( alibavaEvent -> getCollection ( getInputCollectionName ( ) ) ) ;
	noOfChip = collectionVec -> getNumberOfElements ( );
	
	for ( int i = 0; i < noOfChip; ++i )
	{
	    // get data from the collection
	    TrackerDataImpl * trkdata = dynamic_cast < TrackerDataImpl * > ( collectionVec -> getElementAt ( i ) ) ;
	    chipnum = getChipNum ( trkdata );
	    FloatVec datavec;
	    datavec = trkdata -> getChargeValues ( );
	    FloatVec newdatavec;
	    newdatavec.clear ( );
	    FloatVec pedVec = getPedestalOfChip ( chipnum );
	    // now subtract pedestal values from all channels
	    for ( size_t ichan = 0; ichan < datavec.size ( ); ichan++ )
	    {
		if ( isMasked ( chipnum, ichan ) )
		{
		    newdatavec.push_back ( 0 );
		    continue;
		}
		// now subtract pedestal
		float newdata = datavec[ichan] - pedVec[ichan];
		//store it in new data vector
		newdatavec.push_back ( newdata );
	    }
	    TrackerDataImpl * newDataImpl = new TrackerDataImpl ( );
	    newDataImpl -> setChargeValues ( newdatavec );
	    chipIDEncoder[ALIBAVA::ALIBAVADATA_ENCODE_CHIPNUM] = chipnum;
	    chipIDEncoder.setCellID ( newDataImpl );
	    newDataCollection -> push_back ( newDataImpl );
	}
	alibavaEvent -> addCollection ( newDataCollection, getOutputCollectionName ( ) );
    }
    catch ( lcio::DataNotAvailableException& )
    {
	// do nothing again
	streamlog_out ( ERROR5 ) << "Collection (" << getInputCollectionName ( ) << ") not found! " << endl;
    }
}

void AlibavaPedestalSubtraction::check ( LCEvent * /* evt */ )
{
    // nothing to check here - could be used to fill check plots in reconstruction processor
}

void AlibavaPedestalSubtraction::end ( )
{
	if ( _numberOfSkippedEvents > 0 )
	{
	    streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents << " events skipped since they are masked" << endl;
	}
	streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
}

void AlibavaPedestalSubtraction::fillHistos ( TrackerDataImpl * /* trkdata */ )
{

}

void AlibavaPedestalSubtraction::bookHistos ( )
{

}
