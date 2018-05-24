/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@cern.ch
 *
 */

// this processor is built only if USE_AIDA and USE_GEAR are defined
#if defined ( USE_GEAR ) && ( defined ( USE_AIDA ) || defined ( MARLIN_USE_AIDA ) )

// eutelescope inlcudes
#include "EUTelFixDUTAlignment.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelExceptions.h"
#include "EUTelAlignmentConstant.h"
#include "EUTelSimpleVirtualCluster.h"
#include "EUTelGenericSparseClusterImpl.h"
#include "EUTelGeometricClusterImpl.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelSparseClusterImpl.h"

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
#include <UTIL/CellIDDecoder.h>
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

EUTelFixDUTAlignment::EUTelFixDUTAlignment ( ) : Processor ( "EUTelFixDUTAlignment" ), _siPlanesParameters ( nullptr ), _siPlanesLayerLayout ( nullptr )
{
    _description = "Does some very hacky stuff: Mode 0: creates a user-specified DUT alignment. Mode 1: creates 4 dummy hits for each telescope plane and the DUT. The total DUT alignment steps can then be applied to these hits to get the final alignment of the DUT sensor.";

    registerProcessorParameter ( "Mode", "Which mode: 0 or 1", _mode, 1 );

    registerProcessorParameter ( "OutputFileName", "The output file to write into (both modes)", _outputfilename, std::string ( "file.slcio" ) );

    registerProcessorParameter ( "EventsToWrite", "How many dummy events should be written?", _eventstowrite, 1 );

    registerOptionalParameter ( "DUTshiftX","Mode 0: DUT shift in X.", _shiftx , 0.0 );

    registerOptionalParameter ( "DUTshiftY","Mode 0: DUT shift in Y.", _shifty , 0.0 );

    registerOptionalParameter ( "DUTshiftZ","Mode 0: DUT shift in Z.", _shiftz , 0.0 );

    registerOptionalParameter ( "DUTrotA","Mode 0: DUT rotation in A.", _rota , 0.0 );

    registerOptionalParameter ( "DUTrotB","Mode 0: DUT rotation in B.", _rotb , 0.0 );

    registerOptionalParameter ( "DUTrotC","Mode 0: DUT rotation in C.", _rotc , 0.0 );

    registerOptionalParameter ( "ManualDUTID","Mode 0: The sensor ID of the DUT.", _manualDUTid, 6 );

    registerOptionalParameter ( "HitCollectionName", "The output hit collection name (mode 1)", _dummyhitcollectionname, std::string ( "hit" ) );

}

void EUTelFixDUTAlignment::init ( )
{
    if ( Global::GEAR == nullptr )
    {
	streamlog_out ( ERROR5 ) << "The GearMgr is not available, for an unknown reason." << endl;
	exit ( -1 );
    }

    _siPlanesParameters  = const_cast < gear::SiPlanesParameters* > ( & ( Global::GEAR -> getSiPlanesParameters ( ) ) );
    _siPlanesLayerLayout = const_cast < gear::SiPlanesLayerLayout* > ( & ( _siPlanesParameters -> getSiPlanesLayerLayout ( ) ) );

    // writer
    LCWriter * lcWriter = LCFactory::getInstance ( ) -> createLCWriter ( );
    try 
    {
	lcWriter -> open ( _outputfilename, LCIO::WRITE_NEW );
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

    for ( int j = 0; j < _eventstowrite; j++ )
    {
	// an event:
	LCEventImpl * event = new LCEventImpl;
	event -> setRunNumber ( 0 );
	event -> setEventNumber ( j );

	if ( _mode == 0 )
	{

	    streamlog_out ( MESSAGE1 ) << "Starting mode 0!" << endl;
	    // the alignment constant collection we want to write
	    // FIXME this assumes the DUT is in the middle of the telescope, with telescope planes 012 before it and 345 after it

	    LCCollectionVec * constantsCollection = new LCCollectionVec ( LCIO::LCGENERICOBJECT );

	    // upstream:
	    for ( int i = 0; i < 3; i++ )
	    {
		EUTelAlignmentConstant * constant = new EUTelAlignmentConstant;
		constant -> setSensorID ( i );
		constant -> setXOffset ( 0.0 );
		constant -> setYOffset ( 0.0 );
		constant -> setZOffset ( 0.0 );
		constant -> setAlpha ( 0.0 );
		constant -> setBeta ( 0.0 );
		constant -> setGamma ( 0.0 );
		constantsCollection -> push_back ( constant );
	    }

	    // DUT:
	    EUTelAlignmentConstant * constant = new EUTelAlignmentConstant;
	    constant -> setSensorID ( _manualDUTid );
	    constant -> setXOffset ( _shiftx );
	    constant -> setYOffset ( _shifty );
	    constant -> setZOffset ( _shiftz );
	    constant -> setAlpha ( _rota );
	    constant -> setBeta ( _rotb );
	    constant -> setGamma ( _rotc );
	    constantsCollection -> push_back ( constant );

	    // downstream:
	    for ( int i = 3; i < 6; i++ )
	    {
		EUTelAlignmentConstant * constant = new EUTelAlignmentConstant;
		constant -> setSensorID ( i );
		constant -> setXOffset ( 0.0 );
		constant -> setYOffset ( 0.0 );
		constant -> setZOffset ( 0.0 );
		constant -> setAlpha ( 0.0 );
		constant -> setBeta ( 0.0 );
		constant -> setGamma ( 0.0 );
		constantsCollection -> push_back ( constant );
	    }

	    // collection name alignment is hardcoded, as it is in EuTelMille
	    // in future versions of EUTelescope, this might change
	    event -> addCollection ( constantsCollection, "alignment" );

	    streamlog_out ( MESSAGE1 ) << "Finished mode 0!" << endl;
	}

	// make a dummy hit
	if ( _mode == 1 )
	{

	    streamlog_out ( MESSAGE1 ) << "Starting mode 1!" << endl;

	    int nplanes = _siPlanesLayerLayout -> getNLayers ( );
	    LCCollectionVec * constantsCollection = new LCCollectionVec ( LCIO::TRACKERHIT );
	    CellIDEncoder < TrackerHitImpl > idHitEncoder ( EUTELESCOPE::HITENCODING, constantsCollection );
	    for ( int j = 0; j < nplanes; j++ )
	    {

		double zThickness = _siPlanesLayerLayout -> getSensitiveThickness ( j );
		double zposition = _siPlanesLayerLayout -> getLayerPositionZ ( j );
		double xpos[4] = { -10.0, 10.0, -10.0, 10.0 };
		double ypos[4] = { -10.0, -10.0, 10.0, 10.0 };
		double zpos[4] = { 0, 0, 0, 0 };

		double gRotation[3] = { 0.0, 0.0, 0.0 };
		gRotation[0] = _siPlanesLayerLayout -> getLayerRotationXY ( j );
		gRotation[1] = _siPlanesLayerLayout -> getLayerRotationZX ( j );
		gRotation[2] = _siPlanesLayerLayout -> getLayerRotationZY ( j );
		gRotation[0] =  gRotation[0] * 3.1415926 / 180.0;
		gRotation[1] =  gRotation[1] * 3.1415926 / 180.0;
		gRotation[2] =  gRotation[2] * 3.1415926 / 180.0;

		//FIXME
		int sensorid = -1;
		if ( j < 3 )
		{
		    sensorid = j;
		}
		if ( j == 3 )
		{
		    sensorid = _manualDUTid;
		}
		if ( j > 3 )
		{
		    sensorid = j - 1;
		}

		for ( int i = 0; i < 4; i++ )
		{
		    double telpos[3];
		    telpos[0] = xpos[i];
		    telpos[1] = ypos[i];
		    telpos[2] = zpos[i];

		    _EulerRotation ( telpos, gRotation );

		    telpos[2] += zposition + 0.5 * zThickness;
		    lcio::TrackerHitImpl * constant = new lcio::TrackerHitImpl;
		    constant -> setPosition ( &telpos[0] );

		    float cov[TRKHITNCOVMATRIX] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

		    double resx = 0.0;
		    double resy = 0.0;
		    cov[0] = resx * resx;
		    cov[2] = resy * resy;
		    constant -> setCovMatrix ( cov );
		    constant -> setTime ( 0 );
		    idHitEncoder["sensorID"] =  sensorid ;
		    idHitEncoder["properties"] = kHitInGlobalCoord;
		    idHitEncoder.setCellID ( constant );
		    constantsCollection -> push_back ( constant );
		}

	    }

	    event -> addCollection ( constantsCollection, _dummyhitcollectionname );

	    streamlog_out ( MESSAGE1 ) << "Finished mode 1!" << endl;

	}

	// output all this
	lcWriter -> writeEvent ( event );
	delete event;
    }

    lcWriter -> close ( );

}

void EUTelFixDUTAlignment::check ( EVENT::LCEvent* )
{

}

void EUTelFixDUTAlignment::end ( )
{

}

void EUTelFixDUTAlignment::processEvent ( EVENT::LCEvent* )
{

}

void EUTelFixDUTAlignment::processRunHeader ( EVENT::LCRunHeader* )
{

}

void EUTelFixDUTAlignment::_EulerRotation ( double* _telPos, double* _gRotation )
{

    TVector3 _UnrotatedSensorHit ( _telPos[0], _telPos[1], 0.0 );
    TVector3 _RotatedSensorHit ( _telPos[0], _telPos[1], 0.0 );
    TVector3 _Xaxis ( 1.0, 0.0, 0.0 );
    TVector3 _Yaxis ( 0.0, 1.0, 0.0 );
    TVector3 _Zaxis ( 0.0, 0.0, 1.0 );
    if ( TMath::Abs ( _gRotation[2] ) > 1e-6 )
    {
	_RotatedSensorHit.Rotate ( _gRotation[2], _Xaxis );
    }
    if ( TMath::Abs ( _gRotation[1] ) > 1e-6 )
    {
	_RotatedSensorHit.Rotate ( _gRotation[1], _Yaxis );
    }
    if ( TMath::Abs ( _gRotation[0] ) > 1e-6 )
    {
	_RotatedSensorHit.Rotate ( _gRotation[0], _Zaxis ); // in XY
    }
    _telPos[0] = _RotatedSensorHit.X ( );
    _telPos[1] = _RotatedSensorHit.Y ( );
    _telPos[2] = _telPos[2] + _RotatedSensorHit.Z ( );

}

#endif
