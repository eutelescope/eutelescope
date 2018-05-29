/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

// alibava includes ".h"
#include "AlibavaFilter.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"
#include "AlibavaPedNoiCalIOManager.h"

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
#include <UTIL/LCTOOLS.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCEventImpl.h>

// ROOT includes ".h"
#include <TH1.h>
#include <TH1D.h>

// system includes <>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <memory>
#include <fstream>
#include <numeric>

// max number of FIR filter coefficients
#define nc 5

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;

AlibavaFilter::AlibavaFilter ( ) : AlibavaBaseProcessor ( "AlibavaFilter" )
{

    // modify processor description
    _description = "AlibavaFilter filters input data to eliminate crosstalk and non-gaussian noise (Random Ghost Hits)!";

    registerInputCollection ( LCIO::TRACKERDATA, "InputCollectionName", "Input reco data collection name", _inputCollectionName, string ( "recodata_cmmd" ) );

    registerOutputCollection ( LCIO::TRACKERDATA, "OutputCollectionName", "Output filtered data collection name", _outputCollectionName, string ( "recodata_filtered" ) );

    registerOptionalParameter ( "Coefficient1", "First correction value", _initcoefficient1, 0.0373f );

    registerOptionalParameter ( "Coefficient2", "Second correction value", _initcoefficient2, 0.0162f );

    registerProcessorParameter ( "UseSimpleMethod", "Set to true to use Coefficient1 and Coefficient2 or read from file. This should be sufficient for most applications. If false, then the (define nc N) filter coefficients in the source code will be used.", _simplemethod, true );

    registerOptionalParameter ( "ReadFIRCoefficients", "FIR filter coefficients from a previous iteration can be read if this is switched on.", _readcoefficients, false );

    registerOptionalParameter ( "FIRCoefficientFile", "The filename to read/write coefficients to", _filterFileName, string ( "filtercoefficients.txt" ) );

    registerProcessorParameter ( "RGHfilter", "Set to true to use RGH filtering instead of FIR filtering. RGH filtering assumes a polarity of +1!", _rghcorrection, false );

    registerOptionalParameter ( "MinRGHnoise", "The minimum noise a channel has to have to consider RGH filtering. Set to 0 to do all good channels.", _minrghnoise, 6.0f );

    registerOptionalParameter ( "MaxSignalFactor", "The maximum signal a channel is allowed to have. This factor times a channels noise is the limit. Should be larger than the cluster seed limit, otherwise no clusters will be found.", _maxsignalfactor, 10.0f );

    registerOptionalParameter ( "MaxRGHADC", "The maximum ADC in a channel to be allowed. Higher ADCs will be set to 0 if RGH filtering is used.", _maxrghadc, 150.0f );

    registerOptionalParameter ( "SeedCut", "The used cluster seed cut, used to determine the number of cluster candidates in RGH filtering", _seedcut, 5.0f );

    registerOptionalParameter ( "MaxSuspects", "The highest number of cluster candidates in an event allowed before discarding it in RGH filtering.", _maxsuspects, 1 );

    registerOptionalParameter ( "NegativeNoise", "The ratio of noise in negative ADCs required for a neighbour to cut a seed.", _rghnegnoisecut, 2.0f );

    registerProcessorParameter ( "NoiseCollectionName", "Noise collection name", _noiseCollectionName, string ( "finalnoise" ) );

    registerProcessorParameter ( "NoiseInputFile", "The filename where the final noise is stored", _pedestalFile, string ( "finalnoise.slcio" ) );

}

void AlibavaFilter::init ( )
{
    streamlog_out ( MESSAGE4 ) << "Running init" << endl;

    printParameters ( );

    if ( Global::parameters -> isParameterSet ( ALIBAVA::CHANNELSTOBEUSED ) )
    {
	Global::parameters -> getStringVals ( ALIBAVA::CHANNELSTOBEUSED, _channelsToBeUsed );
    }
    else
    {
	streamlog_out ( MESSAGE4 ) << "The Global Parameter "<< ALIBAVA::CHANNELSTOBEUSED <<" is not set!" << endl;
    }

    if ( Global::parameters -> isParameterSet ( ALIBAVA::SKIPMASKEDEVENTS ) )
    {
	_skipMaskedEvents = bool ( Global::parameters -> getIntVal ( ALIBAVA::SKIPMASKEDEVENTS ) );
    }
    else
    {
	streamlog_out ( MESSAGE4 ) << "The Global Parameter " << ALIBAVA::SKIPMASKEDEVENTS << " is not set! Masked events will be used!" << endl;
    }

    _readcoefficient1 = 0.0;
    _readcoefficient2 = 0.0;

    // I assume this is good for about 3m of flatband ribbon between alibava daughterboard and alibava motherboard in the DESY testbeam
    // this might differ for other setups, number was found by hand :-)
    //_initcoefficient1 = 0.0373;
    //_initcoefficient2 = 0.0162;

    // load coefficients from file
    if ( _readcoefficients == true )
    {

	ifstream fileRead;
	fileRead.open ( _filterFileName.c_str ( ) );

	if ( fileRead.is_open ( ) )
	{
	    fileRead >> _readcoefficient1 >> _readcoefficient2;
	    streamlog_out ( MESSAGE4 ) << "Filter coefficients successfully loaded and set to " << _readcoefficient1 << " and " << _readcoefficient2 << " !" << endl;
	}
	else
	{
	    streamlog_out ( ERROR1 ) << "Unable to open file " << _filterFileName << " !" << endl;
	}
    }
    _dropsuspectcount = 0;
}

void AlibavaFilter::processRunHeader ( LCRunHeader * rdr )
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

void AlibavaFilter::processEvent ( LCEvent * anEvent )
{

    if ( anEvent -> getEventNumber ( ) % 1000 == 0 )
    {
	streamlog_out ( MESSAGE4 ) << "Looping events " << anEvent -> getEventNumber ( ) << endl;
    }

    AlibavaEventImpl * alibavaEvent = static_cast < AlibavaEventImpl* > ( anEvent );

    if ( _skipMaskedEvents && ( alibavaEvent -> isEventMasked ( ) )  )
    {
	_numberOfSkippedEvents++;
	return;
    }

    LCCollectionVec * newColVec = new LCCollectionVec ( LCIO::TRACKERDATA );
    CellIDEncoder < TrackerDataImpl > chipIDEncoder ( ALIBAVA::ALIBAVADATA_ENCODE, newColVec );
    LCCollectionVec * collectionVec;

    unsigned int noOfChip;
    try
    {
	collectionVec = dynamic_cast < LCCollectionVec * > ( alibavaEvent -> getCollection ( getInputCollectionName ( ) ) ) ;
	noOfChip = collectionVec -> getNumberOfElements ( );

	for ( size_t i = 0; i < noOfChip; ++i )
	{
	    TrackerDataImpl * trkdata = dynamic_cast < TrackerDataImpl * > ( collectionVec -> getElementAt ( i ) ) ;
	    TrackerDataImpl * newdataImpl = new TrackerDataImpl ( );
	    chipIDEncoder[ALIBAVA::ALIBAVADATA_ENCODE_CHIPNUM] = _chipSelection[i];
	    chipIDEncoder.setCellID ( newdataImpl );

	    FloatVec datavec, newdatavec;
	    datavec = trkdata -> getChargeValues ( );

	    // simple method with 2 coefficients and adding back the subtracted charge -> 4th order filter
	    // should be sufficient for most applications

	    if ( _simplemethod == true && _rghcorrection == false )
	    {
		double input_buffer[ALIBAVA::NOOFCHANNELS];
		double chargein = 0.0;
		double chargeout = 0.0;

		// set the input buffer
		for ( int ii = 0; ii < ALIBAVA::NOOFCHANNELS; ii++ )
		{
		    input_buffer[ii] = datavec.at ( ii );
		    streamlog_out ( DEBUG5 ) << "Reading in: Event: " << anEvent -> getEventNumber ( ) << ", chip: " << i << " , channel: " << ii << " , ADC: " << input_buffer[ii] << " !" << endl;
		    chargein += datavec.at ( ii );
		}

		// read reverse and subtract the crosstalk
		for ( int ii = ( ALIBAVA::NOOFCHANNELS - 1); ii >= 0; ii-- )
		{
		    // magic lines:
		    if ( ii == 0 || ii == 1 )
		    {
			input_buffer[ii] = input_buffer[ii];
		    }
		    if ( ii > 1 )
		    {
			// the charge we are subtracting: c1, c2 has to be added to where it crosstalked from
			// this way, the total adc count stays equal
			double c1 = ( _readcoefficient1 + _initcoefficient1 ) * input_buffer[ii - 1];
			double c2 = ( _readcoefficient2 + _initcoefficient2 ) * input_buffer[ii - 2];
			input_buffer[ii] = input_buffer[ii] - c1 - c2;
			input_buffer[ii - 1] = input_buffer[ii - 1] + c1;
			input_buffer[ii - 2] = input_buffer[ii - 2] + c2;
		    }
		}

		// write out again
		for ( int ii = 0; ii < ALIBAVA::NOOFCHANNELS; ii++ )
		{
		    newdatavec.push_back ( input_buffer[ii] );
		    chargeout += input_buffer[ii];
		}

		// we should not lose charge!
		if ( ( chargein - chargeout ) > 1 )
		{
		    streamlog_out ( ERROR1 ) << "Charge loss in filtering is : " << chargein - chargeout << " ADCs!" << endl;
		}
	    }

	    // advanced filtering with nc coefficients
	    // define nc at the top!
	    if ( _simplemethod == false && _rghcorrection == false )
	    {

		int ii = 0;
		float y;
		float b[nc];
		float input_buffer[nc];

		// filtercoefficients b[0] = b_0, ..., b[nc - 1] = b_N
		// these must be found by hand or by using some noise calculation software
		// add as many as you define in nc at the top!
		b[0] = 0.001;
		b[1] = 0.001;
		b[2] = 1;
		b[3] = -0.0373;
		b[4] = -0.0162;

		// set the input buffer to zero
		for ( ii = 0; ii < nc; ii++ )
		{
		    input_buffer[ii] = 0;
		}

		// loop over the input data
		for ( unsigned int j = 0; j < datavec.size ( ); j++ )
		{
		    // move elements in the buffer to the right
		    for ( ii = nc - 1; ii > 0; ii-- )
		    {
			input_buffer[ii] = input_buffer[ii - 1];
		    }

		    // fill the buffer with input data
		    input_buffer[0] = datavec.at ( j );

		    // calculate output
		    y = 0;

		    for ( ii = 0; ii < nc; ii++ )
		    {
			y += ( b[ii] * input_buffer[ii] );
		    }

		    streamlog_out ( DEBUG5 ) << "FIR Filter: Element: " << j << " Input was: " << datavec.at ( j ) << " , output: " << y << " !" << endl;

		    // set the output
		    newdatavec.at ( j ) = y;
		}

	    }

	    // rgh correction:

	    // the count on seed candidates in this event:
	    // discarded rghs are counted as candidates, since the "real hit" could have been in the same channel
	    int suspects = 0;

	    // this should only be done for n-type with polarity +1
	    if ( _rghcorrection == true )
	    {
		// loop the channels
		for ( int ii = 0; ii < ALIBAVA::NOOFCHANNELS; ii++ )
		{

		    // this noise of this channel
		    float noise = getNoiseAtChannel ( i, ii );

		    if ( datavec.at ( ii ) > noise )
		    {
			streamlog_out ( DEBUG5 ) << "Reading in: Event: " << anEvent -> getEventNumber ( ) << ", chip: " << i << " , channel: " << ii << " , ADC: " << datavec.at ( ii ) << " !" << endl;
		    }

		    // if it is over this limit, then there might be a rgh hit in here... if not push the original data back...
		    if ( noise > _minrghnoise )
		    {
			// replace all high adcs with zero
			// they also count as seed candidates
			float templimit = _maxrghadc;
			if ( noise * _maxsignalfactor < _maxrghadc )
			{
			    templimit = noise * _maxsignalfactor;
			}

			bool badchan = false;

			if ( datavec.at ( ii ) >= templimit )
			{
			    badchan = true;
			    suspects++;

			    // for reference we plot the amount of removed charge
			    if ( TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap["RGHfilteredchargeHighADC"] ) )
			    {
				histo -> Fill ( datavec.at ( ii ) );
			    }
			    if ( TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap["RGHPosition"] ) )
			    {
				histo -> Fill ( ii + i * ALIBAVA::NOOFCHANNELS );
			    }
			    streamlog_out ( DEBUG5 ) << "Chan " << ii << ": ADC over limit ( " << _maxrghadc << " / " << noise * _maxsignalfactor << " ) - " << suspects << " suspects in this event!" << endl;
			}

			// if a channel is below the max adc but over the seed cut
			if ( datavec.at ( ii ) > _seedcut * noise && datavec.at ( ii ) < templimit )
			{

			    // noise left/right
			    float leftnoise = getNoiseAtChannel ( i, ii - 1 );
			    float rightnoise = getNoiseAtChannel ( i, ii + 1 );

			    // if a seed candidate has high negative neighbours, it is probably a rgh too -> discard
			    if ( ii > 0 && ii < 127 )
			    {
				if ( datavec.at ( ii - 1 ) < ( -1 * _rghnegnoisecut * leftnoise ) || datavec.at ( ii + 1 ) < ( -1 * _rghnegnoisecut * rightnoise ) )
				{
				    suspects++;
				    badchan = true;

				    // for reference we plot the amount of removed charge
				    if ( TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap["RGHfilteredchargeNegNeigh"] ) )
				    {
					histo -> Fill ( datavec.at ( ii ) );
				    }
				    if ( TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap["RGHPosition"] ) )
				    {
					histo -> Fill ( ii + i * ALIBAVA::NOOFCHANNELS );
				    }
				    streamlog_out ( DEBUG5 ) << "Chan " << ii << ": Negative neighbour (l " << datavec.at ( ii - 1 ) << " /r " << datavec.at ( ii + 1 ) << " ) - " << suspects << " suspects in this event!" << endl;

				}
				else
				{
				    // candidate for a real hit
				    suspects++;
				    badchan = false;
				    streamlog_out ( DEBUG5 ) << "Chan " << ii << ": Real hit candidate!" << endl;
				}
			    }
			}

			// good or bad chan?
			if ( badchan == true )
			{
			    newdatavec.push_back ( 0.0 );
			}
			if ( badchan == false )
			{
			    newdatavec.push_back ( datavec.at ( ii ) );
			}
		    }
		    else
		    {
			// noise is small, keep adc
			newdatavec.push_back ( datavec.at ( ii ) );
		    }
		}

		// if there are too many seed candidates, we remove the event
		if ( suspects >= _maxsuspects )
		{
		    _dropsuspectcount++;
		    float droppedeventcharge = std::accumulate ( newdatavec.begin ( ), newdatavec.end ( ), 0.0 );
		    newdatavec.clear ( );
		    newdatavec.assign ( ALIBAVA::NOOFCHANNELS, 0.0 );

		    // for reference we plot the amount of removed charge
		    if ( TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap["RGHfilteredchargeCandidates"] ) )
		    {
			histo -> Fill ( droppedeventcharge );
		    }
		}
	    }

	    newdataImpl -> setChargeValues ( newdatavec );
	    newColVec -> push_back ( newdataImpl );

	}

	alibavaEvent -> addCollection ( newColVec, getOutputCollectionName ( ) );

    }
    catch ( lcio::DataNotAvailableException& )
    {
	// do nothing again
	streamlog_out ( ERROR5 ) << "Collection (" << getInputCollectionName ( ) << ") not found! " << endl;
    }
}

void AlibavaFilter::check ( LCEvent * /* evt */ )
{
    // nothing to check here - could be used to fill check plots in reconstruction processor
}

void AlibavaFilter::end ( )
{
    if ( _numberOfSkippedEvents > 0 )
    {
	streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents << " events skipped since they are masked" << endl;
    }
    streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
    streamlog_out ( MESSAGE4 ) << " " << endl;
    streamlog_out ( MESSAGE4 ) << "Dropped " << _dropsuspectcount << " events, since they failed the suspect count in RGH filtering!" << endl;

}

void AlibavaFilter::fillHistos ( TrackerDataImpl * /* trkdata */ )
{

}

void AlibavaFilter::bookHistos ( )
{

    AIDAProcessor::tree ( this ) -> cd ( this -> name ( ) );

    // a histogram showing the charge dropped due to high adcs
    string tempHistoName = "RGHfilteredchargeHighADC";
    stringstream tempHistoTitle;
    tempHistoTitle << tempHistoName << ";ADCs;Entries";
    TH1D * highadchisto = new TH1D ( tempHistoName.c_str ( ), "", 500, 0, 500 );
    _rootObjectMap.insert ( make_pair ( tempHistoName, highadchisto ) );
    string tmp_string = tempHistoTitle.str ( );
    highadchisto -> SetTitle ( tmp_string.c_str ( ) );

    // a histogram showing the charge dropped due to negative neighbours
    string tempHistoName0 = "RGHfilteredchargeNegNeigh";
    stringstream tempHistoTitle0;
    tempHistoTitle0 << tempHistoName0 << ";ADCs;Entries";

    TH1D * negneighhisto = new TH1D ( tempHistoName0.c_str ( ), "", 500, 0, 500 );
    _rootObjectMap.insert ( make_pair ( tempHistoName0, negneighhisto ) );
    string tmp_string0 = tempHistoTitle0.str ( );
    negneighhisto -> SetTitle ( tmp_string0.c_str ( ) );

    // a histogram showing the charge dropped due to number of seed candidates
    string tempHistoName1 = "RGHfilteredchargeCandidates";
    stringstream tempHistoTitle1;
    tempHistoTitle1 << tempHistoName1 << ";ADCs;Entries";

    TH1D * candidatehisto = new TH1D ( tempHistoName1.c_str ( ), "", 1000, -500, 500 );
    _rootObjectMap.insert ( make_pair ( tempHistoName1, candidatehisto ) );
    string tmp_string1 = tempHistoTitle1.str ( );
    candidatehisto -> SetTitle ( tmp_string1.c_str ( ) );

    // a histogram showing the RGH position on the sensor
    string tempHistoName2 = "RGHPosition";
    stringstream tempHistoTitle2;
    tempHistoTitle2 << tempHistoName2 << ";ADCs;Entries";

    TH1D * positionhisto = new TH1D ( tempHistoName2.c_str ( ), "", 256, 0, 255 );
    _rootObjectMap.insert ( make_pair ( tempHistoName2, positionhisto ) );
    string tmp_string2 = tempHistoTitle2.str ( );
    positionhisto -> SetTitle ( tmp_string2.c_str ( ) );

    streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << endl;
}
