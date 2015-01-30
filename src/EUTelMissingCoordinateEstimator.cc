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
_iRun(0),
_iEvt(0)
{
    // modify processor description
    _description =  "EUTelMissingCoordinateEstimator As the name suggest this processor is finds the position of the missing coordinate on your How it works is simple, it gets the hits from specified two finds the closest hit pairs, make a straight line out of it and the estimated position in one axis on your sensor you want. No promises that this will work with tilted sensors and/or with magnetic field. One needs to used this with merged hits and after pre-alignment";
    
    registerInputCollection(LCIO::TRACKERHIT,"InputHitCollectionName", "Input hit collection name. Hits should be in global coordinates and pre-aligned", _inputHitCollectionName, std::string(""));
    
    registerOutputCollection(LCIO::TRACKERHIT,"OutputHitCollectionName", "Output hit collection name", _outputHitCollectionName, std::string(""));
    
    registerProcessorParameter("ReferencePlanes","This is the list of sensorIDs that their hits will be used to estimate the missing coordinate on your DUT. You have to give exactly 2 sensorIDs. For better results use the ones that are closest to your DUT", _referencePlanes, EVENT::IntVec() );
    
    registerProcessorParameter("DUTPlanes","This is the list of sensorIDs that missing coordinate of their hits needs to be found. Notice that if the specified coordinate already exists it will be overwritten", _dutPlanes, EVENT::IntVec() );
    
    registerProcessorParameter("MissingCoordinate","The coordinate axis that needs to be estimated. You have to set this to either X or Y.", _missingCoordinate, string("X") );
    
    
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
    
}




void EUTelMissingCoordinateEstimator::processRunHeader (LCRunHeader * rdr) {
    
    
    auto_ptr<EUTelRunHeaderImpl> header ( new EUTelRunHeaderImpl (rdr) );
    header->addProcessor( type() );
    
    
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
            }
        }
        
        // Store all telescope hits in the new collection, we will store DUT hits after updating its position
        if (!isDUTHit) {
            outputHitCollection->push_back( cloneHit(inputHit) );
        }
      }
    
    
    
    //HERE
    
    // now using reference planes fins the missing coordinate
    
    
    
    
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
    
    // copy dEdX
    newHit->setdEdx( inputHit->getdEdx() );
    
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
    streamlog_out ( MESSAGE4 )  << "Successfully finished" << endl;
}

void EUTelMissingCoordinateEstimator::bookHistos() {
    
    // nothing to book
}
