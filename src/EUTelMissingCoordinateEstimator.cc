/*
 * Created by Eda Yildirim
 *  (2015 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


// ROOT includes:

// eutelescope includes ".h"
#include "EUTelMissingCoordinateEstimator.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"


// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"

// gear includes <.h>


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

EUTelMissingCoordinateEstimator::EUTelMissingCoordinateEstimator () : Processor("EUTelMissingCoordinateEstimator"),
_inputHitCollectionName(),
_outputHitCollectionName(),
_referencePlanes(),
_dutPlanes(),
_missingCoordinate(),
_maxResidual(0),
_iRun(0),
_iEvt(0),
_missingHitPos(0),
_knownHitPos(0),
_nDutHits(0),
_nDutHitsCreated(0),
_maxExpectedCreatedHitPerDUTHit(10),
_numberOfCreatedHitPerDUTHit()
{
    // modify processor description
    _description =  "EUTelMissingCoordinateEstimator As the name suggest this processor is finds the position of the missing coordinate on your How it works is simple, it gets the hits from specified two finds the closest hit pairs, make a straight line out of it and the estimated position in one axis on your sensor you want. No promises that this will work with tilted sensors and/or with magnetic field. One needs to used this with merged hits and after pre-alignment";
    
    registerInputCollection(LCIO::TRACKERHIT,"InputHitCollectionName", "Input hit collection name. Hits should be in global coordinates and pre-aligned", _inputHitCollectionName, std::string(""));
    
    registerOutputCollection(LCIO::TRACKERHIT,"OutputHitCollectionName", "Output hit collection name", _outputHitCollectionName, std::string(""));
    
    registerProcessorParameter("ReferencePlanes","This is the list of sensorIDs that their hits will be used to estimate the missing coordinate on your DUT. You have to give exactly 2 sensorIDs. For better results use the ones that are closest to your DUT", _referencePlanes, EVENT::IntVec() );
    
    registerProcessorParameter("DUTPlanes","This is the list of sensorIDs that missing coordinate of their hits needs to be found. Notice that if the specified coordinate already exists it will be overwritten", _dutPlanes, EVENT::IntVec() );
    
    registerProcessorParameter("MissingCoordinate","The coordinate axis that needs to be estimated. You have to set this to either X or Y.", _missingCoordinate, string("X") );
    
    registerProcessorParameter("MaxResidual","This processor will look for a closest hits (in known coordinate) to determine if the hits are correlated. The hits will be considered as correlated if the residual is smaller than MaxResidual", _maxResidual, float(0) );
    
}


void EUTelMissingCoordinateEstimator::init() {
    // this method is called only once even when the rewind is active
    // usually a good idea to
    printParameters ();
    
    // set to zero the run and event counters
    _iRun = 0;
    _iEvt = 0;
    
    // check if _referencePlanes has only and exactly two sensorIDs
    if (_referencePlanes.size() != 2) {
        streamlog_out (ERROR4) << "ReferencePlanes has to contain exactly 2 sensorIDs!"<< std::endl;
        exit(-1);
    }
    
    // check if _missingCoordinate is valid
    // if it is given as lowercase make them uppercase letters
    if (_missingCoordinate == string("x")) _missingCoordinate = string("X");
    if (_missingCoordinate == string("y")) _missingCoordinate = string("Y");
    
    if (_missingCoordinate == string("X") || _missingCoordinate == string("Y") ) {
        streamlog_out (DEBUG4) << "MissingCoordinate value set as: "<< _missingCoordinate << std::endl;
    }
    else {
        streamlog_out (ERROR4) << "MissingCoordinate value ("<<_missingCoordinate<<") is not valid!"<< std::endl;
        exit(-1);
    }
    
    // now set missing and known hit position variables
    if (_missingCoordinate == string("X")) {
        _missingHitPos = 0;
        _knownHitPos = 1;
    }
    if (_missingCoordinate == string("Y")) {
        _missingHitPos = 1;
        _knownHitPos = 0;
    }
    
    
    // set counters to zero
    _nDutHits = 0;
    _nDutHitsCreated = 0;

   for (unsigned int i=0; i < _maxExpectedCreatedHitPerDUTHit+1; i++){
    	_numberOfCreatedHitPerDUTHit.push_back(0);
   }
}




void EUTelMissingCoordinateEstimator::processRunHeader (LCRunHeader * rdr) {
    std::unique_ptr<EUTelRunHeaderImpl> header = std::make_unique<EUTelRunHeaderImpl>(rdr);
    header->addProcessor(type());
    
    // this is the right place also to check the geometry ID. This is a
    // unique number identifying each different geometry used at the
    // beam test. The same number should be saved in the run header and
    // in the xml file. If the numbers are different, instead of barely
    // quitting ask the user what to do.
    
    if ( header->getGeoID() == 0 )
        streamlog_out ( WARNING0 ) <<  "The geometry ID in the run header is set to zero." << endl
        <<  "This may mean that the GeoID parameter was not set" << endl;
    
    // increment the run counter
    ++_iRun;
}


void EUTelMissingCoordinateEstimator::processEvent (LCEvent * event) {
    
    ++_iEvt;
    
    EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event) ;
    
    if ( evt->getEventType() == kEORE ) {
        streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
        return;
    } else if ( evt->getEventType() == kUNKNOWN ) {
        streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
        << " is of unknown type. Continue considering it as a normal Data Event." << endl;
    }
    

    LCCollectionVec * inputHitCollection     = 0;
    LCCollectionVec * outputHitCollection     = 0;
    
    try
    {
        inputHitCollection = static_cast<LCCollectionVec*> (event->getCollection( _inputHitCollectionName ));
    }
    catch (DataNotAvailableException& e  )
    {
        streamlog_out  ( MESSAGE2 ) <<  "No input collection " << _inputHitCollectionName << " found on event " << event->getEventNumber()
        << " in run " << event->getRunNumber() << endl;
        return ;
    }
    
    
    try
    {
        outputHitCollection  = static_cast<LCCollectionVec*> (event->getCollection( _outputHitCollectionName ));
    }
    catch(...)
    {
        outputHitCollection = new LCCollectionVec(LCIO::TRACKERHIT);
    }
    // prepare an encoder for the hit collection
    CellIDEncoder<TrackerHitImpl> outputCellIDEncoder(EUTELESCOPE::HITENCODING, outputHitCollection);
    
    CellIDDecoder<TrackerHitImpl>   inputCellIDDecoder( inputHitCollection );
  
    vector<int> referencePlaneHits1;
    vector<int> referencePlaneHits2;
    vector<int> dutPlaneHits;

    // Here identify which hits come from reference planes or DUT
    for ( int iInputHits = 0; iInputHits < inputHitCollection->getNumberOfElements(); iInputHits++ )
    {
        TrackerHitImpl * inputHit = dynamic_cast<TrackerHitImpl*> ( inputHitCollection->getElementAt( iInputHits ) );

        int sensorID    = inputCellIDDecoder(inputHit)["sensorID"];

        bool isDUTHit = false;
        // store the reference plane hits
        if (sensorID == _referencePlanes[0]) referencePlaneHits1.push_back(iInputHits);
        if (sensorID == _referencePlanes[1]) referencePlaneHits2.push_back(iInputHits);
        for (unsigned int i=0; i<_dutPlanes.size(); i++) {
            if (sensorID == _dutPlanes[i]) {
                dutPlaneHits.push_back(iInputHits);
                isDUTHit = true;
                _nDutHits++;
            }
        }
        
        // Store all telescope hits in the new collection, we will store DUT hits after updating its position
        if (!isDUTHit) {
            outputHitCollection->push_back( cloneHit(inputHit) );
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

    // with countCreatedDutHits vector we will count how many hits we create out of one DUT hit
    vector<unsigned int> countCreatedDutHits;
    countCreatedDutHits.clear();
	for (unsigned int iDutHit=0; iDutHit<dutPlaneHits.size(); iDutHit++){
		countCreatedDutHits.push_back(0);
	}
 
    // loop over first reference plane hits
    for (unsigned int iHitRefPlane1=0; iHitRefPlane1<referencePlaneHits1.size(); iHitRefPlane1++) {
        TrackerHitImpl * refHit1 = dynamic_cast<TrackerHitImpl*> ( inputHitCollection->getElementAt( referencePlaneHits1[iHitRefPlane1] ) );
        const double* refHit1Pos = refHit1->getPosition();
        
        // loop over second reference plane hits
        for (unsigned int iHitRefPlane2=0; iHitRefPlane2<referencePlaneHits2.size(); iHitRefPlane2++) {
            TrackerHitImpl * refHit2 = dynamic_cast<TrackerHitImpl*> ( inputHitCollection->getElementAt( referencePlaneHits2[iHitRefPlane2] ) );
            const double* refHit2Pos = refHit2->getPosition();
            
            // loop over second dut plane hits
            for (unsigned int iDutHit=0; iDutHit<dutPlaneHits.size(); iDutHit++) {
                TrackerHitImpl * dutHit = dynamic_cast<TrackerHitImpl*> ( inputHitCollection->getElementAt( dutPlaneHits[iDutHit] ) );
                const double* dutHitPos = dutHit->getPosition();
                double newDutHitPos[3];
                
                // find t value of the Z position of DUT
                // t = (z-z1)/(z2-z1)
                double t = ( dutHitPos[2] - refHit1Pos[2] ) / ( refHit2Pos[2] - refHit1Pos[2] );
                
                // find the known coordinate value correcponds to that z on the line
                double knownHitPosOnLine = refHit1Pos[_knownHitPos] + (refHit2Pos[_knownHitPos] - refHit1Pos[_knownHitPos]) * t;
                
                // if knownHitPosOnLine is close to the knownHitPos
                if ( fabs( knownHitPosOnLine - dutHitPos[_knownHitPos] ) < _maxResidual) {
                    // first copy old DUT hit position to the new one
                    newDutHitPos[0] = dutHitPos[0];
                    newDutHitPos[1] = dutHitPos[1];
                    newDutHitPos[2] = dutHitPos[2];
                    
                    // then replace the unknown one with the estimated one
                    double estimatedHitPos = refHit1Pos[_missingHitPos] + (refHit2Pos[_missingHitPos] - refHit1Pos[_missingHitPos]) * t;
                    
                    newDutHitPos[_missingHitPos] = estimatedHitPos;
                    
                    // now store new hit position in the TrackerHit, copy and store in the collection
                    
                    TrackerHitImpl * newHit = cloneHit(dutHit);
			const double* hitpos = newDutHitPos;
                    newHit->setPosition( &hitpos[0] );
                    outputHitCollection->push_back(newHit);

                    // count new created hits
                    _nDutHitsCreated++;

		    // increase the created DUT hits
		    countCreatedDutHits[iDutHit] ++;
                }
                
            } // end of loop over first dut plane hits
            
            
        } // end of loop over second reference plane hits
    } // end of loop over first reference plane hits
    
    for (unsigned int iDutHit=0; iDutHit<dutPlaneHits.size(); iDutHit++){
	if(_maxExpectedCreatedHitPerDUTHit < countCreatedDutHits[iDutHit]){
		_numberOfCreatedHitPerDUTHit[_maxExpectedCreatedHitPerDUTHit] ++;
	}
	else {
    		_numberOfCreatedHitPerDUTHit[countCreatedDutHits[iDutHit]] ++;
    	}
    }
    
    
    try
    {
        event->getCollection( _outputHitCollectionName ) ;
    }
    catch(...)
    {
        event->addCollection( outputHitCollection, _outputHitCollectionName );
    }
    
    if ( isFirstEvent() ) _isFirstEvent = false;
}


TrackerHitImpl* EUTelMissingCoordinateEstimator::cloneHit(TrackerHitImpl *inputHit){
    TrackerHitImpl * newHit = new TrackerHitImpl;
    
    // copy hit position
    const double* hitPos = inputHit->getPosition();
    newHit->setPosition( &hitPos[0] );
    
    // copy cov. matrix
    newHit->setCovMatrix( inputHit->getCovMatrix() );
    
    // copy type
    newHit->setType( inputHit->getType() );
    
    // copy rawhits
    LCObjectVec clusterVec = inputHit->getRawHits();
    newHit->rawHits() = clusterVec;
   
    // copy cell IDs
    newHit->setCellID0( inputHit->getCellID0() );
    newHit->setCellID1( inputHit->getCellID1() );
    
    // copy EDep
    newHit->setEDep( inputHit->getEDep() );

    // copy EDepError
    newHit->setEDepError( inputHit->getEDepError() );
    
    // copy Time
    newHit->setTime( inputHit->getTime() );
    
    // copy Quality
    newHit->setQuality( inputHit->getQuality() );

    return newHit;
}


void EUTelMissingCoordinateEstimator::end()
{
    
    streamlog_out ( MESSAGE4 )  << "Number of hits you had from all DUTs "<< _nDutHits << endl;
    streamlog_out ( MESSAGE4 )  << "Number of hits created with the estimated missing coordinate "<< _nDutHitsCreated << endl;
    for (unsigned int i=0; i<_numberOfCreatedHitPerDUTHit.size(); i++){
	
  	streamlog_out ( MESSAGE4 )  << "You created "<< i ;
	if(i==_maxExpectedCreatedHitPerDUTHit) streamlog_out ( MESSAGE4 )  << " or more";
	streamlog_out ( MESSAGE4 )  << " hits per DUT hit "<< _numberOfCreatedHitPerDUTHit[i] <<" many times"<<endl;
    }
    streamlog_out ( MESSAGE4 )  << "Successfully finished" << endl;
    
}

void EUTelMissingCoordinateEstimator::bookHistos() {
    
    // nothing to book
}
