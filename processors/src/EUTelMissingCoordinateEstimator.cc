/*
 * Created by Eda Yildirim
 *  (2015 DESY)
 *
 *  email:eda.yildirim@cern.ch
 *
 * Updated by Thomas Eichhorn
 *  (2017 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

// eutelescope includes ".h"
#include "EUTelMissingCoordinateEstimator.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"
#include "marlin/AIDAProcessor.h"

// aida includes <.h>
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
using namespace marlin;
using namespace gear;
using namespace eutelescope;

// definition of static members mainly used to name histograms
AIDA::IHistogram2D * hitmaphisto;
AIDA::IHistogram2D * failhitmaphisto;
AIDA::IHistogram1D * faildistancehisto;
AIDA::IHistogram2D * pointhitmaphisto;

EUTelMissingCoordinateEstimator::EUTelMissingCoordinateEstimator ( ) : Processor ( "EUTelMissingCoordinateEstimator" ),
_inputHitCollectionName ( ),
_outputHitCollectionName ( ),
_referencePlanes ( ),
_dutPlanes ( ),
_missingCoordinate ( ),
_maxResidual ( 0 ),
_missingHitPos ( 0 ),
_knownHitPos ( 0 ),
_nDutHits ( 0 ),
_nDutHitsCreated ( 0 ),
_maxExpectedCreatedHitPerDUTHit ( 10 )
{
    // modify processor description
    _description =  "EUTelMissingCoordinateEstimator:This processor estimates the missing coordinate on a strip sensor by extrapolating a straight line from two reference planes. No promises that this will work with tilted sensors and/or with magnetic fields. The merged input hits should be pre aligned for better results.";

    registerInputCollection ( LCIO::TRACKERHIT, "InputHitCollectionName", "Input hit collection name. Hits should be in global coordinates and pre-aligned", _inputHitCollectionName, std::string ( "" ) );

    registerOutputCollection ( LCIO::TRACKERHIT, "OutputHitCollectionName", "Output hit collection name", _outputHitCollectionName, std::string ( "" ) );

    registerProcessorParameter ( "ReferencePlanes","List of sensorIDs of which hits will be used to estimate the missing coordinate on the DUT. You have to give exactly 2 sensorIDs. For better results use the ones that are closest to your DUT", _referencePlanes, EVENT::IntVec ( ) );

    registerProcessorParameter ( "DUTPlanes", "List of sensorIDs with a missing coordinate to be found. Note that if the specified coordinate already exists it will be overwritten", _dutPlanes, EVENT::IntVec ( ) );

    registerProcessorParameter ( "MissingCoordinate", "The coordinate axis that needs to be estimated. You have to set this to either X or Y.", _missingCoordinate, string ( "X" ) );

    registerProcessorParameter ( "MaxResidual", "This processor will look for hits in the known coordinate to determine if the hits are correlated. The hits will be considered as correlated if the residual is smaller than MaxResidual", _maxResidual, float ( 10.0 ) );

    registerProcessorParameter ( "MultiHitMode", "Allow an individual DUT hit to be transformed into multiple hits? If false, only the closest extrapolated position will be used.", _multihitmode, true );
    

}


void EUTelMissingCoordinateEstimator::init ( )
{

    printParameters ( );

    bookHistos ( );

    // check if _referencePlanes has only and exactly two sensorIDs
    if ( _referencePlanes.size ( ) != 2 )
    {
	streamlog_out ( ERROR5 ) << "ReferencePlanes has to contain exactly 2 sensorIDs!" << std::endl;
        exit ( -1 );
    }

    // check if _missingCoordinate is valid, if it is given as lowercase convert to uppercase
    if ( _missingCoordinate == string ( "x" ) )
    {
	_missingCoordinate = string ( "X" );
    }
    if ( _missingCoordinate == string ( "y" ) )
    {
	_missingCoordinate = string ( "Y" );
    }

    if ( _missingCoordinate == string ( "X" ) || _missingCoordinate == string ( "Y" ) )
    {
	streamlog_out ( MESSAGE4 ) << "MissingCoordinate value set as: " << _missingCoordinate << std::endl;
    }
    else
    {
	streamlog_out ( ERROR5 ) << "MissingCoordinate value:" << _missingCoordinate << " is not valid!" << std::endl;
	exit ( -1 );
    }

    // now set missing and known hit position variables
    if ( _missingCoordinate == string ( "X" ) )
    {
	_missingHitPos = 0;
	_knownHitPos = 1;
    }
    if ( _missingCoordinate == string ( "Y" ) )
    {
	_missingHitPos = 1;
	_knownHitPos = 0;
    }

    // set global counters to zero
    _nDutHits = 0;
    _nDutHitsCreated = 0;
    _nNoReferenceCount = 0;
    _nResidualFailCount = 0;
    for ( unsigned int i = 0; i < 10; i++ )
    {
	_numberOfCreatedHitsPerDUTHit[i] = 0 ;
    }

}


void EUTelMissingCoordinateEstimator::processRunHeader ( LCRunHeader * rdr )
{
    std::unique_ptr < EUTelRunHeaderImpl > header = std::make_unique < EUTelRunHeaderImpl > ( rdr );
    header -> addProcessor ( type ( ) );
}


void EUTelMissingCoordinateEstimator::processEvent ( LCEvent * event )
{

    EUTelEventImpl * evt = static_cast < EUTelEventImpl* > ( event );

    if ( evt -> getEventType ( ) == kEORE )
    {
	streamlog_out ( MESSAGE4 ) << "EORE found: nothing else to do." << endl;
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
	streamlog_out ( DEBUG5 ) << "No input collection " << _inputHitCollectionName << " found on event " << event -> getEventNumber ( ) << " in run " << event -> getRunNumber ( ) << endl;
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

    // prepare an encoder for the hit collection
    CellIDEncoder < TrackerHitImpl > outputCellIDEncoder ( EUTELESCOPE::HITENCODING, outputHitCollection );

    vector < int > referencePlaneHits1;
    vector < int > referencePlaneHits2;
    vector < int > dutPlaneHits;

    // identify which hits come from the reference planes or the DUT
    for ( int i = 0; i < inputHitCollection -> getNumberOfElements ( ); i++ )
    {
	UTIL::CellIDDecoder < TrackerHitImpl > inputCellIDDecoder ( EUTELESCOPE::HITENCODING );

	TrackerHitImpl * inputHit = dynamic_cast < TrackerHitImpl* > ( inputHitCollection -> getElementAt ( i ) );
	int sensorID = inputCellIDDecoder ( inputHit ) ["sensorID"];

	bool isDUTHit = false;

	// store the reference plane hits
	if ( sensorID == _referencePlanes[0] )
	{
	    referencePlaneHits1.push_back ( i );
	}
	if ( sensorID == _referencePlanes[1] )
	{
	    referencePlaneHits2.push_back ( i );
	}

	for ( unsigned int j = 0; j < _dutPlanes.size ( ); j++)
	{
	    if ( sensorID == _dutPlanes[j] )
	    {
		dutPlaneHits.push_back ( i );
		isDUTHit = true;
		_nDutHits++;
	    }
	}

	// store all telescope hits in the new collection, we will store the DUT hits after updating their position
	if ( !isDUTHit )
	{
	    outputHitCollection -> push_back ( cloneHit ( inputHit ) );
	}
    }

    /*
     The line that passes through 2 points can be written as L(t)= P1 + V*t
     where V is the displacement vector and P1 is the starting point
     so L(t) becomes: L(t) = (x1,y1,z1) + (x2-x1, y2-y1, z2-z1)*t
     * x=x1+(x2−x1)t
     * y=y1+(y2−y1)t
     * z=z1+(z2−z1)t
     */

    // ugly way, to be optimised
    if ( _multihitmode == true )
    {
	// loop over DUT hits
	for ( unsigned int k = 0; k < dutPlaneHits.size ( ); k++ )
	{
	    TrackerHitImpl * dutHit = dynamic_cast < TrackerHitImpl* > ( inputHitCollection -> getElementAt ( dutPlaneHits[k] ) );
	    const double* dutHitPos = dutHit -> getPosition ( );
	    double newDutHitPos[3];

	    int hitsperhit = 0;

	    // loop over first reference plane hits
	    for ( unsigned int i = 0; i < referencePlaneHits1.size ( ); i++ )
	    {
		TrackerHitImpl * refHit1 = dynamic_cast < TrackerHitImpl* > ( inputHitCollection -> getElementAt ( referencePlaneHits1[i] ) );
		const double* refHit1Pos = refHit1 -> getPosition ( );

		// loop over second reference plane hits
		for ( unsigned int j = 0; j < referencePlaneHits2.size ( ); j++ )
		{
		    TrackerHitImpl * refHit2 = dynamic_cast < TrackerHitImpl* > ( inputHitCollection -> getElementAt ( referencePlaneHits2[j] ) );
		    const double* refHit2Pos = refHit2 -> getPosition ( );

		    // t = (z-z1)/(z2-z1)
		    double t = ( dutHitPos[2] - refHit1Pos[2] ) / ( refHit2Pos[2] - refHit1Pos[2] );

		    // find the known coordinate value that correcponds to that z on the line
		    double knownHitPosOnLine = refHit1Pos[_knownHitPos] + ( refHit2Pos[_knownHitPos] - refHit1Pos[_knownHitPos] ) * t;

		    pointhitmaphisto -> fill ( ( refHit1Pos[0] + ( refHit2Pos[0] - refHit1Pos[0] ) * t ), ( refHit1Pos[1] + ( refHit2Pos[1] - refHit1Pos[1] ) * t ) );

		    // if knownHitPosOnLine is close to the actual DUT hit position
		    if ( fabs ( knownHitPosOnLine - dutHitPos[_knownHitPos] ) < _maxResidual )
		    {
			// first copy old DUT hit position to the new one
			newDutHitPos[0] = dutHitPos[0];
			newDutHitPos[1] = dutHitPos[1];
			newDutHitPos[2] = dutHitPos[2];

			// then replace the unknown one with the estimated one
			double estimatedHitPos = refHit1Pos[_missingHitPos] + ( refHit2Pos[_missingHitPos] - refHit1Pos[_missingHitPos] ) * t;

			newDutHitPos[_missingHitPos] = estimatedHitPos;

			// now store new hit position in the collection
			TrackerHitImpl * newHit = cloneHit ( dutHit );
			const double* hitpos = newDutHitPos;
			newHit -> setPosition ( &hitpos[0] );
			outputHitCollection -> push_back ( newHit );
			streamlog_out ( DEBUG0 ) << "New hit: x: " << newDutHitPos[0] << ", y: " << newDutHitPos[1] << ", z: " << newDutHitPos[2] << endl;
			hitmaphisto -> fill ( newDutHitPos[0], newDutHitPos[1] );

			// count new created hits
			hitsperhit ++;
			_nDutHitsCreated++;

		    }
		    else
		    {
			streamlog_out ( DEBUG0 ) << "Failing residual cut!" << endl;
			failhitmaphisto -> fill ( dutHitPos[0], dutHitPos[1] );
			faildistancehisto -> fill ( knownHitPosOnLine - dutHitPos[_knownHitPos] );
			_nResidualFailCount++;
		    }

		} // end of loop over second reference plane hits
	    } // end of loop over first reference plane hits

	    if ( hitsperhit > 9 )
	    {
		hitsperhit = 9;
	    }
	    _numberOfCreatedHitsPerDUTHit[hitsperhit]++;
	} // end of loop over DUT hits
    } // end _multihitmode == true

    // ugly way, to be optimised
    if ( _multihitmode == false )
    {
	bool foundhit = false;
	double closestresidual = _maxResidual;
	// loop over DUT hits
	for ( unsigned int k = 0; k < dutPlaneHits.size ( ); k++ )
	{
	    TrackerHitImpl * dutHit = dynamic_cast < TrackerHitImpl* > ( inputHitCollection -> getElementAt ( dutPlaneHits[k] ) );
	    const double* dutHitPos = dutHit -> getPosition ( );
	    double newDutHitPos[3];

	    int hitsperhit = 0;

	    // loop over first reference plane hits
	    for ( unsigned int i = 0; i < referencePlaneHits1.size ( ); i++ )
	    {
		TrackerHitImpl * refHit1 = dynamic_cast < TrackerHitImpl* > ( inputHitCollection -> getElementAt ( referencePlaneHits1[i] ) );
		const double* refHit1Pos = refHit1 -> getPosition ( );

		// loop over second reference plane hits
		for ( unsigned int j = 0; j < referencePlaneHits2.size ( ); j++ )
		{
		    TrackerHitImpl * refHit2 = dynamic_cast < TrackerHitImpl* > ( inputHitCollection -> getElementAt ( referencePlaneHits2[j] ) );
		    const double* refHit2Pos = refHit2 -> getPosition ( );

		    // t = (z-z1)/(z2-z1)
		    double t = ( dutHitPos[2] - refHit1Pos[2] ) / ( refHit2Pos[2] - refHit1Pos[2] );

		    // find the known coordinate value that correcponds to that z on the line
		    double knownHitPosOnLine = refHit1Pos[_knownHitPos] + ( refHit2Pos[_knownHitPos] - refHit1Pos[_knownHitPos] ) * t;

		    pointhitmaphisto -> fill ( ( refHit1Pos[0] + ( refHit2Pos[0] - refHit1Pos[0] ) * t ), ( refHit1Pos[1] + ( refHit2Pos[1] - refHit1Pos[1] ) * t ) );

		    // if knownHitPosOnLine is close to the actual DUT hit position
		    if ( fabs ( knownHitPosOnLine - dutHitPos[_knownHitPos] ) < closestresidual )
		    {
			foundhit = true;
			closestresidual = knownHitPosOnLine - dutHitPos[_knownHitPos];
			// first copy old DUT hit position to the new one
			newDutHitPos[0] = dutHitPos[0];
			newDutHitPos[1] = dutHitPos[1];
			newDutHitPos[2] = dutHitPos[2];

			// then replace the unknown one with the estimated one
			double estimatedHitPos = refHit1Pos[_missingHitPos] + ( refHit2Pos[_missingHitPos] - refHit1Pos[_missingHitPos] ) * t;

			newDutHitPos[_missingHitPos] = estimatedHitPos;
		    }

		} // end of loop over second reference plane hits
	    } // end of loop over first reference plane hits

	    if ( foundhit == true )
	    {
		// now store new hit position in the collection
		TrackerHitImpl * newHit = cloneHit ( dutHit );
		const double* hitpos = newDutHitPos;
		newHit -> setPosition ( &hitpos[0] );
		outputHitCollection -> push_back ( newHit );
		streamlog_out ( DEBUG0 ) << "New hit: x: " << newDutHitPos[0] << ", y: " << newDutHitPos[1] << ", z: " << newDutHitPos[2] << endl;
		hitmaphisto -> fill ( newDutHitPos[0], newDutHitPos[1] );

		// count new created hits
		hitsperhit ++;
		_nDutHitsCreated++;
	    }

	    if ( hitsperhit > 9 )
	    {
		hitsperhit = 9;
	    }
	    _numberOfCreatedHitsPerDUTHit[hitsperhit]++;
	} // end of loop over DUT hits
    } // end _multihitmode == false

    if ( referencePlaneHits1.size ( ) == 0 || referencePlaneHits2.size ( ) == 0 )
    {
	streamlog_out ( DEBUG5 ) << "Couldn't create hit, no input in reference plane 1 or 2!" << endl;
	_nNoReferenceCount++;
    }

    try
    {
	event -> getCollection ( _outputHitCollectionName ) ;
    }
    catch ( ... )
    {
	event -> addCollection ( outputHitCollection, _outputHitCollectionName );
    }

}


TrackerHitImpl* EUTelMissingCoordinateEstimator::cloneHit ( TrackerHitImpl *inputHit )
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


void EUTelMissingCoordinateEstimator::end ( )
{

    streamlog_out ( MESSAGE4 )  << "Number of input DUT hits:   "<< _nDutHits << endl;
    streamlog_out ( MESSAGE4 )  << "Number of DUT hits created: "<< _nDutHitsCreated << endl;
    streamlog_out ( MESSAGE4 )  << " " << endl;
    streamlog_out ( MESSAGE4 )  << "Number of DUT hits failing the residual cut : "<< _nResidualFailCount << endl;
    streamlog_out ( MESSAGE4 )  << "Number of DUT hits with no reference hits : "<< _nNoReferenceCount << endl;
    for ( unsigned int i = 0; i < 10; i++ )
    {
	streamlog_out ( MESSAGE4 )  << "You created " << i << " hits per DUT hit "<< _numberOfCreatedHitsPerDUTHit[i] << " times" << endl;
    }
    streamlog_out ( MESSAGE4 )  << "Successfully finished" << endl;
}


void EUTelMissingCoordinateEstimator::bookHistos ( )
{

    hitmaphisto = AIDAProcessor::histogramFactory ( this ) -> createHistogram2D ( "hitmap", 1000, -20, 20, 1000, -20, 20 );
    hitmaphisto -> setTitle ( "Created Hitmap;X [mm];Y [mm]" );

    failhitmaphisto = AIDAProcessor::histogramFactory ( this ) -> createHistogram2D ( "failhitmap", 1000, -20, 20, 1000, -20, 20 );
    failhitmaphisto -> setTitle ( "Hitmap where no hit was created;X [mm];Y [mm]" );

    faildistancehisto = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( "faildistance", 1000, -50, 50 );
    faildistancehisto -> setTitle ( "Residual distance where no hit was created;Distance [mm];Entries" );

    pointhitmaphisto = AIDAProcessor::histogramFactory ( this ) -> createHistogram2D ( "pointhitmap", 1000, -20, 20, 1000, -20, 20 );
    pointhitmaphisto -> setTitle ( "Hitmap of points;X [mm];Y [mm]" );

}
