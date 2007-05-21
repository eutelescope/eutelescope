// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelAutoPedestalNoiseProcessor.cc,v 1.2 2007-05-21 11:46:24 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h" 
#include "EUTelAutoPedestalNoiseProcessor.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"

// lcio includes <.h> 
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDEncoder.h>


// system includes <>


using namespace std;
using namespace marlin;
using namespace eutelescope;



EUTelAutoPedestalNoiseProcessor::EUTelAutoPedestalNoiseProcessor () : Processor("EUTelAutoPedestalNoiseProcessor"),
							    _pedestalCollectionVec(NULL),
							    _noiseCollectionVec(NULL),
							    _statusCollectionVec(NULL) {

  // modify processor description
  _description =
    "EUTelAutoPedestalNoiseProcessor produces initial pedestal / noise / status with user provided values";

 
  registerOutputCollection (LCIO::TRACKERDATA, "PedestalCollectionName",
			   "Pedestal local collection",
			   _pedestalCollectionName, string ("pedestal"));

  registerOutputCollection (LCIO::TRACKERDATA, "NoiseCollectionName",
			   "Noise local collection",
			   _noiseCollectionName, string("noise"));

  registerOutputCollection (LCIO::TRACKERRAWDATA, "StatusCollectionName",
			   "Pixel status collection",
			   _statusCollectionName, string("status"));

  const int nDetectorExample = 6;
  FloatVec initPedeExample(nDetectorExample, 0.);

  registerOptionalParameter("InitPedestalValue",
			    "The initial value of pedestal (one value for detector)",
			    _initPedestal, initPedeExample, initPedeExample.size() );

  FloatVec initNoiseExample(nDetectorExample, 1.);

  registerOptionalParameter("InitNoiseValue",
			    "The initial value of noise (one value for detector)",
			    _initNoise, initNoiseExample, initNoiseExample.size() );  

  
}


void EUTelAutoPedestalNoiseProcessor::init () {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

}

void EUTelAutoPedestalNoiseProcessor::processRunHeader (LCRunHeader * rdr) {

  // increment the run counter
  ++_iRun;

  EUTelRunHeaderImpl * runHeader = static_cast<EUTelRunHeaderImpl*> (rdr);
  unsigned int noOfDetector = runHeader->getNoOfDetector();

  if (noOfDetector != _initPedestal.size()) {
    message<WARNING> ( "Resizing the initial pedestal vector" );
    _initPedestal.resize(noOfDetector, _initPedestal.back());
  }

  if (noOfDetector != _initNoise.size()) {
    message<WARNING> ( "Resizing the initial noise vector");
    _initNoise.resize(noOfDetector, _initNoise.back());
  }

  // now the four vectors containing the first and the last pixel
  // along both the directions
  _minX = runHeader->getMinX();
  _maxX = runHeader->getMaxX();
  _minY = runHeader->getMinY();
  _maxY = runHeader->getMaxY();


}


void EUTelAutoPedestalNoiseProcessor::processEvent (LCEvent * event) {

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event) ;
  
  if ( evt->getEventType() == kEORE ) {
    message<DEBUG> ( "EORE found: nothing else to do." );
    return;
  }


  if (isFirstEvent()) {

    _pedestalCollectionVec = new LCCollectionVec(LCIO::TRACKERDATA);
    _noiseCollectionVec    = new LCCollectionVec(LCIO::TRACKERDATA);
    _statusCollectionVec   = new LCCollectionVec(LCIO::TRACKERRAWDATA);
    
    for (unsigned int iDetector = 0; iDetector < _initPedestal.size(); iDetector++) {

      int nPixel = ( _maxX[iDetector] - _minX[iDetector] + 1 ) * ( _maxY[iDetector] - _minY[iDetector] + 1 ) ;

      TrackerRawDataImpl * status   = new TrackerRawDataImpl;
      CellIDEncoder<TrackerRawDataImpl> statusEncoder(EUTELESCOPE::MATRIXDEFAULTENCODING, _statusCollectionVec);
      statusEncoder["sensorID"] = iDetector;
      statusEncoder["xMin"]     = _minX[iDetector];
      statusEncoder["yMin"]     = _minY[iDetector];
      statusEncoder["xMax"]     = _maxX[iDetector];
      statusEncoder["yMax"]     = _maxY[iDetector];
      statusEncoder.setCellID(status);
      ShortVec statusVec(nPixel, EUTELESCOPE::GOODPIXEL);
      status->setADCValues(statusVec);
      _statusCollectionVec->push_back(status);

      TrackerDataImpl * pedestal    = new TrackerDataImpl;
      CellIDEncoder<TrackerDataImpl>  pedestalEncoder(EUTELESCOPE::MATRIXDEFAULTENCODING, _pedestalCollectionVec);
      pedestalEncoder["sensorID"] = iDetector;
      pedestalEncoder["xMin"]     = _minX[iDetector];
      pedestalEncoder["yMin"]     = _minY[iDetector];
      pedestalEncoder["xMax"]     = _maxX[iDetector];
      pedestalEncoder["yMax"]     = _maxY[iDetector];
      pedestalEncoder.setCellID(pedestal);
      FloatVec pedestalVec(nPixel, _initPedestal[iDetector]);
      pedestal->setChargeValues(pedestalVec);
      _pedestalCollectionVec->push_back(pedestal);

      TrackerDataImpl * noise    = new TrackerDataImpl;
      CellIDEncoder<TrackerDataImpl>  noiseEncoder(EUTELESCOPE::MATRIXDEFAULTENCODING, _noiseCollectionVec);
      noiseEncoder["sensorID"] = iDetector;
      noiseEncoder["xMin"]     = _minX[iDetector];
      noiseEncoder["yMin"]     = _minY[iDetector];
      noiseEncoder["xMax"]     = _maxX[iDetector];
      noiseEncoder["yMax"]     = _maxY[iDetector];
      noiseEncoder.setCellID(noise);    
      FloatVec noiseVec(nPixel, _initNoise[iDetector]);
      noise->setChargeValues(noiseVec);
      _noiseCollectionVec->push_back(noise);
    }

    _isFirstEvent = false;
  }

  evt->addCollection(_pedestalCollectionVec, _pedestalCollectionName);
  evt->takeCollection(_pedestalCollectionName);
  evt->addCollection(_noiseCollectionVec, _noiseCollectionName);
  evt->takeCollection(_noiseCollectionName);
  evt->addCollection(_statusCollectionVec, _statusCollectionName);
  evt->takeCollection(_statusCollectionName);

  ++_iEvt;
  
}
  



void EUTelAutoPedestalNoiseProcessor::end() {

  delete _pedestalCollectionVec;
  delete _noiseCollectionVec;
  delete _statusCollectionVec;

  message<MESSAGE> ( log() << "Successfully finished" ) ;  
}

