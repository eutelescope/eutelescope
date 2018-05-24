/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 *
 */

// this processor is built only if USE_AIDA and USE_GEAR are defined
#if defined ( USE_GEAR ) && ( defined ( USE_AIDA ) || defined ( MARLIN_USE_AIDA ) )

// eutelescope inlcudes
#include "EUTelDUTPosition.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelExceptions.h"
#include "EUTelAlignmentConstant.h"

// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/ITupleFactory.h>

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <Exceptions.h>

#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <memory>
#include <string>
#include <vector>
#include <cstdlib>

#include "TROOT.h"
#include "TVector3.h"

#define PI 3.14159265

using namespace std;
using namespace lcio ;
using namespace marlin ;
using namespace eutelescope;

EUTelDUTPosition::EUTelDUTPosition ( ) : Processor ( "EUTelDUTPosition" ),
_siPlanesParameters ( nullptr ),
_siPlanesLayerLayout ( nullptr )
{

    // modify processor description
    _description = "EUTelDUTPosition gets the final position of sensors and a DUT in the global coordinate system. This can then be used for futher analysis, or just for information ";

    registerOptionalParameter ( "OutputFile", "The sensors' position will be written to this file!", _outputFileName, string ( "fail.txt" ) );

    registerOptionalParameter ( "DUTPositionFile", "The DUT position will be written to this file!", _outputDUTFileName, string ( "fail2.txt" ) );

    registerOptionalParameter ( "ManualDUTPosition", "The sorted sensor position of the DUT (usually 3).", _manualDUTposition , 3 );

    registerOptionalParameter ( "ManualDUTID", "The sensor ID of the DUT.", _manualDUTid , 6 );

    registerOptionalParameter ( "FinalCollection", "The collection name after the final alignment", _finalcollectionname, string ( "AlignedHit10" ) );
}

void EUTelDUTPosition::init ( )
{
    if ( Global::GEAR == nullptr )
    {
	streamlog_out ( ERROR5 ) << "The GearMgr is not available, for an unknown reason." << endl;
	exit ( -1 );
    }

    _siPlanesParameters  = const_cast <gear::SiPlanesParameters* > ( & ( Global::GEAR -> getSiPlanesParameters ( ) ) );
    _siPlanesLayerLayout = const_cast <gear::SiPlanesLayerLayout* > ( & ( _siPlanesParameters -> getSiPlanesLayerLayout ( ) ) );

}

void EUTelDUTPosition::processRunHeader ( EVENT::LCRunHeader* )
{

}

void EUTelDUTPosition::processEvent ( LCEvent * anEvent )
{

    LCCollectionVec * telescopeCollectionVec;
    telescopeCollectionVec = dynamic_cast < LCCollectionVec * > ( anEvent -> getCollection ( _finalcollectionname ) );
    int telescopesize = telescopeCollectionVec -> getNumberOfElements ( );
    streamlog_out ( MESSAGE5 ) << "Dummy hit collection size is: " << telescopesize << endl;
    int nplanes = _siPlanesLayerLayout -> getNLayers ( );

    // writer
    LCWriter * lcWriter = LCFactory::getInstance ( ) -> createLCWriter ( );
    try 
    {
	lcWriter -> open ( _outputFileName, LCIO::WRITE_NEW );
    }
    catch ( IOException& e )
    {
	streamlog_out ( ERROR5 ) << e.what ( ) << endl;
	exit ( -1 );
    }

    // write an almost empty run header
    LCRunHeaderImpl * lcHeader  = new LCRunHeaderImpl;
    lcHeader -> setRunNumber ( 0 );
    lcWriter -> writeRunHeader ( lcHeader );
    delete lcHeader;

    // an event:
    LCEventImpl * event = new LCEventImpl;
    event -> setRunNumber ( 0 );
    event -> setEventNumber ( 0 );
    LCCollectionVec * constantsCollection = new LCCollectionVec ( LCIO::LCGENERICOBJECT );

    for ( int i = 0; i < nplanes; i++ )
    {

	lcio::TrackerHitImpl * input0  = dynamic_cast< lcio::TrackerHitImpl * > ( telescopeCollectionVec->getElementAt( i * 4 + 0 ) ) ;
	lcio::TrackerHitImpl * input1  = dynamic_cast< lcio::TrackerHitImpl * > ( telescopeCollectionVec->getElementAt( i * 4 + 1 ) ) ;
	lcio::TrackerHitImpl * input2  = dynamic_cast< lcio::TrackerHitImpl * > ( telescopeCollectionVec->getElementAt( i * 4 + 2 ) ) ;

	double zThickness = _siPlanesLayerLayout -> getSensitiveThickness ( i );
	double zposition = _siPlanesLayerLayout -> getLayerPositionZ ( i );
	double zexpect = zposition + 0.5 * zThickness;

	streamlog_out ( MESSAGE5 ) << "sensor " << i << " z position is: " << zexpect << endl;

	double * inputPosition0 = const_cast < double * > ( input0 -> getPosition ( ) );
	double * inputPosition1 = const_cast < double * > ( input1 -> getPosition ( ) );
	double * inputPosition2 = const_cast < double * > ( input2 -> getPosition ( ) );

	double xshift = ( inputPosition1[0] - inputPosition0[0] ) / 2.0 + inputPosition0[0];
	double yshift = ( inputPosition2[1] - inputPosition0[1] ) / 2.0 + inputPosition0[1];
	double zshift = ( inputPosition2[2] - inputPosition0[2] ) / 2.0;
	double alpha = atan ( ( inputPosition0[2] - inputPosition2[2] ) / ( inputPosition0[1] - inputPosition2[1]) ) * 180.0 / PI;
	double beta = atan ( ( inputPosition0[2] - inputPosition1[2] ) / ( inputPosition0[0] - inputPosition1[0]) ) * 180.0 / PI;
	double gamma = atan ( ( inputPosition0[1] - inputPosition1[1] ) / ( inputPosition0[0] - inputPosition1[0]) ) * 180.0 / PI;

	//FIXME
	int sensorid = -1;
	if ( i < 3 )
	{
	    sensorid = i;
	}
	if ( i == _manualDUTposition )
	{
	    sensorid = _manualDUTid;
	}
	if ( i > 3 )
	{
	    sensorid = i - 1;
	}

	EUTelAlignmentConstant * constant = new EUTelAlignmentConstant;
	constant -> setSensorID ( sensorid );
	constant -> setXOffset ( -1.0 * xshift );
	constant -> setYOffset ( -1.0 * yshift );
	constant -> setZOffset ( -1.0 * zshift );
	constant -> setAlpha ( -1.0 * alpha );
	constant -> setBeta ( -1.0 * beta );
	constant -> setGamma ( -1.0 * gamma );
	constantsCollection -> push_back ( constant );

	if ( i == _manualDUTposition )
	{
	    streamlog_out ( MESSAGE5 ) << " " << endl;
	    streamlog_out ( MESSAGE5 ) << " DUT: " << endl;
	    streamlog_out ( MESSAGE5 ) << " " << endl;
	    streamlog_out ( MESSAGE5 ) << " x shift     " << xshift << " mm" << endl;
	    streamlog_out ( MESSAGE5 ) << " y shift     " << yshift << " mm" << endl;
	    streamlog_out ( MESSAGE5 ) << " z shift     " << zshift << " mm" << endl;
	    streamlog_out ( MESSAGE5 ) << " a rotation: " << alpha << " °" << endl;
	    streamlog_out ( MESSAGE5 ) << " b rotation: " << beta << " °" << endl;
	    streamlog_out ( MESSAGE5 ) << " c rotation: " << gamma << " °" << endl;
	    streamlog_out ( MESSAGE5 ) << " " << endl;
	    streamlog_out ( MESSAGE5 ) << " " << endl;
	    streamlog_out ( MESSAGE5 ) << " " << endl;

	    ofstream filterFile;
	    filterFile.open ( _outputDUTFileName.c_str ( ) );
	    filterFile << xshift << endl;
	    filterFile << yshift << endl;
	    filterFile << zshift << endl;
	    filterFile << alpha << endl;
	    filterFile << beta << endl;
	    filterFile << gamma << endl;
	    filterFile.close ( );
	}
    }
    event -> addCollection ( constantsCollection, "alignment" );
    // output all this
    lcWriter -> writeEvent ( event );
    delete event;
    lcWriter -> close ( );
}

void EUTelDUTPosition::check ( LCEvent * )
{

}

void EUTelDUTPosition::end ( )
{
    streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
}

void EUTelDUTPosition::bookHistos( )
{

}

#endif
