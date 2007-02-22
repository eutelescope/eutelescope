// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelUpdatePedestalNoiseProcessor.cc,v 1.1 2007-02-22 08:11:51 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h" 
#include "EUTelExceptions.h"
#include "EUTelUpdatePedestalNoiseProcessor.h"
#include "EUTELESCOPE.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h> 
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCCollectionVec.h>

// system includes <>
#include <cmath>

using namespace std;
using namespace marlin;
using namespace eutelescope;



EUTelUpdatePedestalNoiseProcessor::EUTelUpdatePedestalNoiseProcessor () : Processor("EUTelUpdatePedestalNoiseProcessor") {

  // modify processor description
  _description =
    "EUTelUpdatePedestalNoiseProcessor periodically updates the pedestal"
    "and noise values";

  registerInputCollection (LCIO::TRACKERRAWDATA, "RawDataCollectionName",
			   "Raw data collection name",
			   _rawDataCollectionName, string("rawdata"));

  registerInputCollection (LCIO::TRACKERDATA, "PedestalCollectionName",
			   "Pedestal local collection",
			   _pedestalCollectionName, string ("pedestal"));
  
  registerInputCollection (LCIO::TRACKERDATA, "NoiseCollectionName",
			   "Noise local collection",
			   _noiseCollectionName, string("noise"));
  
  registerInputCollection (LCIO::TRACKERRAWDATA, "StatusCollectionName",
			   "Pixel status collection",
			   _statusCollectionName, string("status"));

  registerProcessorParameter("UpdateAlgorithm",
			     "The algorithm to be used for pedestal update",
			     _updateAlgo, string(EUTELESCOPE::FIXEDWEIGHT));

  registerProcessorParameter("UpdateFrequency",
			     "How often the algorithm should be applied",
			     _updateFrequency, static_cast<int>(10));

  registerOptionalParameter("FixedWeightValue",
			    "The value of the fixed weight (only for fixed weight algorithm",
			    _fixedWeight, static_cast<int>(100));
  
}


void EUTelUpdatePedestalNoiseProcessor::init () {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  if ( _updateAlgo == EUTELESCOPE::FIXEDWEIGHT ) {
    if ( _fixedWeight <= 0 ) {
      throw InvalidParameterException("FixedWeightValue has to be a positive integer number");
    }
  }
  
  if ( _updateFrequency <= 0 ) {
    cerr << "[" << name() << "] The update frequency has to be greater than 0. Set it to 1." << endl;
    _updateFrequency = 1;
  }

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

}

void EUTelUpdatePedestalNoiseProcessor::processRunHeader (LCRunHeader * /*rdr*/ ) {

  // increment the run counter
  ++_iRun;

  // reset the event counter
  _iEvt = 0;

}


void EUTelUpdatePedestalNoiseProcessor::processEvent (LCEvent * evt) {

  if ( _iEvt % _updateFrequency == 0 ) {

    cout << "[" << name() << "] Updating pedestal and noise ...";

    if ( _updateAlgo == EUTELESCOPE::FIXEDWEIGHT ) 
      fixedWeightUpdate(evt);

    cout << " ok " << endl;

  }
  ++_iEvt;
  
}
  

void EUTelUpdatePedestalNoiseProcessor::fixedWeightUpdate(LCEvent * evt) {
  
  LCCollectionVec * pedestalCollection = dynamic_cast < LCCollectionVec * > (evt->getCollection(_pedestalCollectionName));
  LCCollectionVec * noiseCollection    = dynamic_cast < LCCollectionVec * > (evt->getCollection(_noiseCollectionName));
  LCCollectionVec * statusCollection   = dynamic_cast < LCCollectionVec * > (evt->getCollection(_statusCollectionName));
  LCCollectionVec * rawDataCollection  = dynamic_cast < LCCollectionVec * > (evt->getCollection(_rawDataCollectionName));


  for (int iDetector = 0; iDetector < statusCollection->getNumberOfElements(); iDetector++) {

    TrackerRawDataImpl * status   = dynamic_cast < TrackerRawDataImpl * > (statusCollection->getElementAt(iDetector));
    TrackerRawDataImpl * rawData  = dynamic_cast < TrackerRawDataImpl * > (rawDataCollection->getElementAt(iDetector));
    TrackerDataImpl    * noise    = dynamic_cast < TrackerDataImpl * >    (noiseCollection->getElementAt(iDetector));
    TrackerDataImpl    * pedestal = dynamic_cast < TrackerDataImpl * >    (pedestalCollection->getElementAt(iDetector));

    for (unsigned int iPixel = 0; iPixel < status->adcValues().size(); iPixel++) {

      if ( status->adcValues()[iPixel] == EUTELESCOPE::GOODPIXEL ) {
	pedestal->chargeValues()[iPixel] = ( (_fixedWeight - 1) * pedestal->chargeValues()[iPixel] + 
					     rawData->getADCValues()[iPixel] ) / _fixedWeight;
	noise->chargeValues()[iPixel]    = sqrt( ( (_fixedWeight - 1) * pow( noise->chargeValues()[iPixel], 2 ) +
						   pow( rawData->getADCValues()[iPixel] - pedestal->chargeValues()[iPixel], 2) ) /
						 _fixedWeight );
      }
    }
  }
}



void EUTelUpdatePedestalNoiseProcessor::end() {

  cout << "[" << name() << "] Successfully finished " << endl;
  
}

