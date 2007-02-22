// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelUpdatePedestalNoiseProcessor.cc,v 1.2 2007-02-22 13:23:39 bulgheroni Exp $
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

#ifdef MARLIN_USE_AIDA
#include "marlin/AIDAProcessor.h"
#include "AIDA/IDataPointSetFactory.h"
#include <AIDA/IDataPointSet.h>
#include <AIDA/IDataPoint.h>
#include <AIDA/IMeasurement.h>

#include <AIDA/IHistogramFactory.h>
#endif

// lcio includes <.h> 
#include <LCIOTypes.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>

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

  IntVec monitorPixelExample;
  monitorPixelExample.push_back(0);
  monitorPixelExample.push_back(10);
  monitorPixelExample.push_back(15);

  registerOptionalParameter("PixelMonitored",
			    "A pixel to be monitored (detectorID, xCoord, yCoord). Add as many line as this as you wish",
			    _monitoredPixel, monitorPixelExample );
  
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

  // reset vectors
  _monitoredPixelPedestal.clear();
  _monitoredPixelNoise.clear();

}

void EUTelUpdatePedestalNoiseProcessor::processRunHeader (LCRunHeader * /*rdr*/ ) {

  // increment the run counter
  ++_iRun;

  // reset the event counter
  _iEvt = 0;

}


void EUTelUpdatePedestalNoiseProcessor::processEvent (LCEvent * evt) {

  if (isFirstEvent()) {
    
    if ( _monitoredPixel.size() != 0 ) {
      // this means we have to fill in the vectors pedestal and noise
      // monitoring

      LCCollectionVec * pedestalCollection = dynamic_cast < LCCollectionVec * > (evt->getCollection(_pedestalCollectionName));
      LCCollectionVec * noiseCollection    = dynamic_cast < LCCollectionVec * > (evt->getCollection(_noiseCollectionName));
    
      unsigned index = 0;
      while ( index < _monitoredPixel.size() ) {
	int  iDetector = _monitoredPixel[index++];
	int  xCoord    = _monitoredPixel[index++];
	int  yCoord    = _monitoredPixel[index++];
	
	TrackerDataImpl * pedestal = dynamic_cast < TrackerDataImpl * >    (pedestalCollection->getElementAt(iDetector));
	TrackerDataImpl * noise    = dynamic_cast < TrackerDataImpl * >    (noiseCollection->getElementAt(iDetector));

	// I need to find the pixel index, for this I need the number of pixels in the x directions.
	CellIDDecoder<TrackerDataImpl> decoder(pedestalCollection);
	int xMin = decoder(pedestal)["xMin"];
	int xMax = decoder(pedestal)["xMax"];
	int yMin = decoder(pedestal)["yMin"];
	int noOfXPixel = abs( xMax - xMin ) + 1;
	if (noOfXPixel <= 0) throw InvalidParameterException("The number of pixels along has to be > 0");
	int pixelIndex = ( xCoord - xMin ) + ( yCoord - yMin ) * noOfXPixel;	

	// initialize the pedestal monitor
	FloatVec pedestalMonitor;
	pedestalMonitor.push_back(pedestal->chargeValues()[pixelIndex]);
	_monitoredPixelPedestal.push_back(pedestalMonitor);

	// and now the noise one
	FloatVec noiseMonitor;
	noiseMonitor.push_back(noise->chargeValues()[pixelIndex]);
	_monitoredPixelNoise.push_back(noiseMonitor);
	
      }  
    }
    _isFirstEvent = false;
  }


  if ( _iEvt % _updateFrequency == 0 ) {

    cout << "[" << name() << "] Updating pedestal and noise ...";

    if ( _updateAlgo == EUTELESCOPE::FIXEDWEIGHT ) 
      fixedWeightUpdate(evt);

    cout << " ok " << endl;

  }
  ++_iEvt;
  
}
  

void EUTelUpdatePedestalNoiseProcessor::pixelMonitoring(LCEvent * evt) {

  LCCollectionVec * pedestalCollection = dynamic_cast < LCCollectionVec * > (evt->getCollection(_pedestalCollectionName));
  LCCollectionVec * noiseCollection    = dynamic_cast < LCCollectionVec * > (evt->getCollection(_noiseCollectionName));
  
  unsigned int index  = 0;
  unsigned int iPixel = 0;
  while (index < _monitoredPixel.size() ) {
    int  iDetector = _monitoredPixel[index++];
    int  xCoord    = _monitoredPixel[index++];
    int  yCoord    = _monitoredPixel[index++];
    
    TrackerDataImpl * pedestal = dynamic_cast < TrackerDataImpl * >    (pedestalCollection->getElementAt(iDetector));
    TrackerDataImpl * noise    = dynamic_cast < TrackerDataImpl * >    (noiseCollection->getElementAt(iDetector));
    
    // I need to find the pixel index, for this I need the number of pixels in the x directions.
    CellIDDecoder<TrackerDataImpl> decoder(pedestalCollection);
    int xMin = decoder(pedestal)["xMin"];
    int xMax = decoder(pedestal)["xMax"];
    int yMin = decoder(pedestal)["yMin"];
    int noOfXPixel = abs( xMax - xMin ) + 1;
    if (noOfXPixel <= 0) throw InvalidParameterException("The number of pixels along has to be > 0");
    int pixelIndex = ( xCoord - xMin ) + ( yCoord - yMin ) * noOfXPixel;
    
    _monitoredPixelPedestal[iPixel].push_back(pedestal->chargeValues()[pixelIndex]);
    _monitoredPixelNoise[iPixel].push_back(noise->chargeValues()[pixelIndex]);
    ++iPixel;
  }

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

#ifdef MARLIN_USE_AIDA
  
  unsigned index  = 0;
  unsigned iPixel = 0;
  while ( index < _monitoredPixel.size() ) {
    int  iDetector = _monitoredPixel[index++];
    int  xCoord    = _monitoredPixel[index++];
    int  yCoord    = _monitoredPixel[index++];
    
    string name;
    string title;

    {
      stringstream namestream;
      namestream << "PedestalMonitor-" << iDetector << "-" << xCoord << "-" << yCoord;
      name = namestream.str();

      stringstream titlestream;
      titlestream << "Pedestal monitoring on detector " << iDetector << "(" << xCoord << "," << yCoord << ")";
      title = titlestream.str();
    }

    //    AIDA::IDataPointSet * pedestalDPS = AIDAProcessor::dataPointSetFactory(this)->create(name,title,1);

    {
      stringstream namestream;
      namestream << "NoiseMonitor-" << iDetector << "-" << xCoord << "-" << yCoord;
      name = namestream.str();

      stringstream titlestream;
      titlestream << "Noise monitoring on detector " << iDetector << "(" << xCoord << "," << yCoord << ")";
      title = titlestream.str();
    }
    
    //    AIDA::IDataPointSet * noiseDPS =  AIDAProcessor::dataPointSetFactory(this)->create();

    for (unsigned int count = 0; count < _monitoredPixelPedestal[iPixel].size(); count++) {
      cout << count << " " << _monitoredPixelPedestal[iPixel][count] << " " << _monitoredPixelNoise[iPixel][count] << endl;
//       pedestalDPS->addPoint();
//       pedestalDPS->point(count)->coordinate(0)->setValue(_monitoredPixelPedestal[iPixel][count]);
      
//       noiseDPS->addPoint();
//       noiseDPS->point(count)->coordinate(0)->setValue(_monitoredPixelNoise[iPixel][count]);
    }
    
    ++iPixel;
  }
    
#endif

  cout << "[" << name() << "] Successfully finished " << endl;
  
}

void EUTelUpdatePedestalNoiseProcessor::check( LCEvent * evt ) {

  if ( (_iEvt - 1) % _updateFrequency == 0 ) 
    pixelMonitoring(evt);

}
