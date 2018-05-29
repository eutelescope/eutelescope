/*
 * Created by Thomas Eichhorn
 *  (2017 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

#include "CBCSimConverter.h"

// marlin includes ".h"
#include "marlin/Processor.h"
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

#include <lcio.h>
#include <Exceptions.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCGenericObjectImpl.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
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
#include "EUTELESCOPE.h"
#include "EUTelGenericSparsePixel.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelGenericSparseClusterImpl.h"
#include "EUTelSimpleVirtualCluster.h"
#include "EUTelGeometricClusterImpl.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "anyoption.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace IMPL;
using namespace eutelescope;


CBCSimConverter::CBCSimConverter ( ) : Processor ( "CBCSimConverter" )
{

    _description = "CBCSimConverter converts allpix simulated data into tracker-data lcio data!";

    registerInputCollection ( LCIO::TRACKERDATA, "InputCollectionName", "Input collection name",  _inputCollectionName, string ( "Det351" ) );

    registerProcessorParameter ( "ChanCount", "The number of channels in the system", _chancount,  254 );

    registerProcessorParameter ( "NonSensitiveAxis", "The unsensitive axis of the CBC", _nonsensitiveaxis, string ( "x" ) );

    registerProcessorParameter ( "OutputCollectionName", "The collection we output", _outputCollectionName, string ( "rawdata" ) );

}


void CBCSimConverter::init ( )
{
    streamlog_out ( MESSAGE4 ) << "Running init" << endl;
    printParameters ( );

    if ( _nonsensitiveaxis == "x" || _nonsensitiveaxis == "X" )
    {
	_nonsensitiveaxis = "x";
	streamlog_out ( MESSAGE4 ) << "Non-sensitive axis is x!" << endl;
    }
    else if ( _nonsensitiveaxis == "y" || _nonsensitiveaxis == "Y" )
    {
	_nonsensitiveaxis = "y";
	streamlog_out ( MESSAGE4 ) << "Non-sensitive axis is y!" << endl;
    }
    else
    {
	streamlog_out ( ERROR5 ) << "Illegal setting for NonSensitiveAxis! Valid coordinates are x and y!" << endl;
	exit ( -1 );
    }

}

void CBCSimConverter::processRunHeader ( LCRunHeader * rdr )
{
    streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;

    // Add processor name to the runheader
    auto arunHeader = std::make_unique < EUTelRunHeaderImpl > ( rdr );
    arunHeader -> addProcessor ( type ( ) );

}

void CBCSimConverter::processEvent ( LCEvent * anEvent )
{

    if ( anEvent -> getEventNumber ( ) % 1000 == 0 )
    {
	streamlog_out ( MESSAGE4 ) << "Looping events " << anEvent -> getEventNumber ( ) << endl;
    }

    // creating LCCollection for raw data
    LCCollectionVec* rawDataCollection = new LCCollectionVec ( LCIO::TRACKERDATA );

    CellIDEncoder < TrackerDataImpl > chipIDEncoder ( eutelescope::EUTELESCOPE::ZSCLUSTERDEFAULTENCODING, rawDataCollection );

    LCCollectionVec * inputVec;
    try
    {
	inputVec = dynamic_cast < LCCollectionVec * > ( anEvent -> getCollection ( _inputCollectionName ) ) ;
	TrackerDataImpl * trkdata = dynamic_cast < TrackerDataImpl * > ( inputVec -> getElementAt ( 0 ) ) ;

	TrackerDataImpl * newdataImpl = new TrackerDataImpl ( );
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
		streamlog_out ( DEBUG1 ) << "chan y " << datavec[i+1] << endl;
		streamlog_out ( DEBUG1 ) << "chan q " << datavec[i+2] << endl;
		streamlog_out ( DEBUG1 ) << "chan t " << datavec[i+3] << endl;
		pixel_x.push_back ( datavec[i] );
		pixel_y.push_back ( datavec[i+1] );
		pixel_charge.push_back ( datavec[i+2] );
		pixel_time.push_back ( datavec[i+3] );
	    }
	}

	// fill the channels with a signal if it's there...
	int point = 0;

	for ( int j = 0; j < _chancount; j++ )
	{

	    double background = 0.0;

	    // add the charge, if any
	    if ( _nonsensitiveaxis == "x" )
	    {
		if ( pixel_y.size ( ) > 0 )
		{
		    if ( pixel_y[point] == ( j ) )
		    {
			double signal = 0.0;
			signal = pixel_charge[point];
			background += signal;
			point++;
			streamlog_out ( DEBUG2 ) << "Output chanel "<< j << " signal is " << signal << endl;
		    }
		}
	    }
	    if ( _nonsensitiveaxis == "y" )
	    {
		if ( pixel_x.size ( ) > 0 )
		{
		    if ( pixel_x[point] == ( j ) )
		    {
			double signal = 0.0;
			signal = pixel_charge[point];
			background += signal;
			point++;
			streamlog_out ( DEBUG2 ) << "Output chanel "<< j << " signal is " << signal << endl;
		    }
		}
	    }
	    // push back this channel
	    outputvec.push_back(background);
	    streamlog_out ( DEBUG3 ) << "Output chanel "<< j << " final signal is " << background << endl;
	}

	// and let there be output
	newdataImpl -> setChargeValues ( outputvec );
	rawDataCollection -> push_back ( newdataImpl );
	anEvent -> addCollection ( rawDataCollection, _outputCollectionName );
    }
    catch ( lcio::DataNotAvailableException& )
    {
	streamlog_out ( ERROR5 ) << "Collection " << _inputCollectionName << " not found! " << endl;
    }

}

void CBCSimConverter::check ( LCEvent * )
{

}

void CBCSimConverter::end ( )
{

}
