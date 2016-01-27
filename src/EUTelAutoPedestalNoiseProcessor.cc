// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
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
#include <memory>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace marlin;
using namespace eutelescope;



EUTelAutoPedestalNoiseProcessor::EUTelAutoPedestalNoiseProcessor ()
 : Processor("EUTelAutoPedestalNoiseProcessor"),
   _pedestalCollectionName(""),
   _noiseCollectionName(""),
   _statusCollectionName(""),
   _initPedestal(),
   _initNoise(),
   _iRun(0),
   _iEvt(0),
   _pedestalCollectionVec(NULL),
   _noiseCollectionVec(NULL),
   _statusCollectionVec(NULL),
   _minX(0),
   _maxX(0),
   _minY(0),
   _maxY(0),
   _sensorIDVec()
{

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

  const size_t nDetectorExample = 6;

  // to get rid of the geometrical information in the run header, here
  // we have to provide the min and max value along x and y
  IntVec minXExample( nDetectorExample, 0 );
  registerOptionalParameter("MinXVector", "The minimum pixel along x (default 0, one value per detector)",
                            _minX, minXExample );

  IntVec maxXExample( nDetectorExample, 255 );
  registerOptionalParameter("MaxXVector", "The maximum pixel along x (default 255, one value per detector)",
                            _maxX, maxXExample );

  IntVec minYExample( nDetectorExample, 0 );
  registerOptionalParameter("MinYVector", "The minimum pixel along y (default 0, one value per detector)",
                            _minY, minYExample );

  IntVec maxYExample( nDetectorExample, 255 );
  registerOptionalParameter("MaxYVector", "The maximum pixel along y (default 255, one value per detector)",
                            _maxY, maxYExample );

  IntVec sensorIDVecExample;
  for ( size_t iDetector = 0 ; iDetector < nDetectorExample; ++iDetector ) {
    sensorIDVecExample.push_back( iDetector );
  }
  registerOptionalParameter("SensorIDVec", "The sensorID for the generated collection (one per detector)",
                            _sensorIDVec, sensorIDVecExample );


  FloatVec initPedeExample(nDetectorExample, 0.);
  registerOptionalParameter("InitPedestalValue",
                            "The initial value of pedestal (one value for detector)",
                            _initPedestal, initPedeExample );

  FloatVec initNoiseExample(nDetectorExample, 1.);

  registerOptionalParameter("InitNoiseValue",
                            "The initial value of noise (one value for detector)",
                            _initNoise, initNoiseExample );


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
  ++_iRun;
  std::unique_ptr<EUTelRunHeaderImpl> runHeader = std::make_unique<EUTelRunHeaderImpl>(rdr);
  runHeader->addProcessor( type() ) ;
}


void EUTelAutoPedestalNoiseProcessor::processEvent (LCEvent * event) {

  ++_iEvt;
  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event) ;

  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return;
  } else if ( evt->getEventType() == kUNKNOWN ) {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }


  if (isFirstEvent()) 
  {

    _pedestalCollectionVec = new LCCollectionVec(LCIO::TRACKERDATA);
    _noiseCollectionVec    = new LCCollectionVec(LCIO::TRACKERDATA);
    _statusCollectionVec   = new LCCollectionVec(LCIO::TRACKERRAWDATA);

   
    for (unsigned int iDetector = 0; iDetector < _initPedestal.size(); iDetector++) 
    {

      int nPixel = ( _maxX[iDetector] - _minX[iDetector] + 1 ) * ( _maxY[iDetector] - _minY[iDetector] + 1 ) ;

      TrackerRawDataImpl * status   = new TrackerRawDataImpl;
      CellIDEncoder<TrackerRawDataImpl> statusEncoder(EUTELESCOPE::MATRIXDEFAULTENCODING, _statusCollectionVec);
      statusEncoder["sensorID"] = _sensorIDVec[iDetector];
      statusEncoder["xMin"]     = _minX[iDetector];
      statusEncoder["yMin"]     = _minY[iDetector];
      statusEncoder["xMax"]     = _maxX[iDetector];
      statusEncoder["yMax"]     = _maxY[iDetector];
      statusEncoder.setCellID(status);
//      ShortVec statusVec(0, EUTELESCOPE::GOODPIXEL);
      ShortVec statusVec(nPixel, EUTELESCOPE::GOODPIXEL);
      status->setADCValues(statusVec);
      _statusCollectionVec->push_back(status);

      TrackerDataImpl * pedestal    = new TrackerDataImpl;
      CellIDEncoder<TrackerDataImpl>  pedestalEncoder(EUTELESCOPE::MATRIXDEFAULTENCODING, _pedestalCollectionVec);
      pedestalEncoder["sensorID"] = _sensorIDVec[iDetector];
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
      noiseEncoder["sensorID"] = _sensorIDVec[iDetector];
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

}




void EUTelAutoPedestalNoiseProcessor::end() {

  delete _pedestalCollectionVec;
  delete _noiseCollectionVec;
  delete _statusCollectionVec;

  streamlog_out ( MESSAGE2 )  << "Successfully finished" << endl;
}

