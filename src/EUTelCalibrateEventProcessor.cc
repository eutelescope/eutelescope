// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelCalibrateEventProcessor.cc,v 1.7 2007-05-29 15:54:48 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h" 
#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"
#include "EUTelCalibrateEventProcessor.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"

#ifdef MARLIN_USE_AIDA
// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h> 
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

// system includes <>
#include <sstream>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;

#ifdef MARLIN_USE_AIDA
string EUTelCalibrateEventProcessor::_rawDataDistHistoName    = "RawDataDistHisto";
string EUTelCalibrateEventProcessor::_dataDistHistoName       = "DataDistHisto";
string EUTelCalibrateEventProcessor::_commonModeDistHistoName = "CommonModeDistHisto";
#endif 

EUTelCalibrateEventProcessor::EUTelCalibrateEventProcessor () :Processor("EUTelCalibrateEventProcessor") {

  // modify processor description
  _description =
    "EUTelCalibrateEventProcessor subtract the pedestal value from the input data";

  // first of all we need to register the input collection
  registerInputCollection (LCIO::TRACKERRAWDATA, "RawDataCollectionName",
			   "Input raw data collection",
			   _rawDataCollectionName, string ("rawdata"));

  registerInputCollection (LCIO::TRACKERDATA, "PedestalCollectionName",
			   "Pedestal from the condition file",
			   _pedestalCollectionName, string ("pedestal"));

  registerInputCollection (LCIO::TRACKERDATA, "NoiseCollectionName",
			   "Noise from the condition file",
			   _noiseCollectionName, string("noise"));

  registerInputCollection (LCIO::TRACKERRAWDATA, "StatusCollectionName",
			   "Pixel status from the condition file",
			   _statusCollectionName, string("status"));

  registerOutputCollection (LCIO::TRACKERDATA, "CalibratedDataCollectionName",
			    "Name of the output calibrated data collection",
			    _calibratedDataCollectionName, string("data"));

  // now the optional parameters
  registerProcessorParameter ("DebugHistoFilling",
			     "Flag to switch on (1) or off (0) the detector debug histogram filling",
			     _fillDebugHisto, static_cast<int> (0));
  
  registerProcessorParameter ("PerformCommonMode",
			      "Flag to switch on (1) or off (0) the common mode suppression algorithm",
			      _doCommonMode, static_cast<int> (1));

  registerProcessorParameter ("HitRejectionCut",
			      "Threshold of pixel SNR for hit rejection",
			      _hitRejectionCut, static_cast<float> (3.5));

  registerProcessorParameter ("MaxNoOfRejectedPixels",
			      "Maximum allowed number of rejected pixel per event",
			      _maxNoOfRejectedPixels, static_cast<int> (3000));


}


void EUTelCalibrateEventProcessor::init () {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  if ( _fillDebugHisto == 1 ) {
    message<WARNING>( "Filling debug histograms is slowing down the procedure");
  }

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

}

void EUTelCalibrateEventProcessor::processRunHeader (LCRunHeader * rdr) {

  // to make things easier re-cast the input header to the EUTelRunHeaderImpl
  EUTelRunHeaderImpl *  runHeader = static_cast<EUTelRunHeaderImpl*>(rdr);

  _noOfDetector = runHeader->getNoOfDetector();

  // increment the run counter
  ++_iRun;

}


void EUTelCalibrateEventProcessor::processEvent (LCEvent * event) {

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);

  if ( evt->getEventType() == kEORE ) {
    message<DEBUG> ( "EORE found: nothing else to do." );
    return;
  }  


  if (_iEvt % 10 == 0) 
    message<MESSAGE> ( log() << "Applying pedestal correction on event " << _iEvt );

  LCCollectionVec * inputCollectionVec    = dynamic_cast < LCCollectionVec * > (evt->getCollection(_rawDataCollectionName));
  LCCollectionVec * pedestalCollectionVec = dynamic_cast < LCCollectionVec * > (evt->getCollection(_pedestalCollectionName));
  LCCollectionVec * noiseCollectionVec    = dynamic_cast < LCCollectionVec * > (evt->getCollection(_noiseCollectionName));
  LCCollectionVec * statusCollectionVec   = dynamic_cast < LCCollectionVec * > (evt->getCollection(_statusCollectionName));
  CellIDDecoder<TrackerRawDataImpl> cellDecoder(inputCollectionVec);

  if (isFirstEvent()) {

    // this is the right place to cross check wheter the pedestal and
    // the input data are at least compatible. I mean the same number
    // of detectors and the same number of pixels in each place.
    
    if ( (inputCollectionVec->getNumberOfElements() != pedestalCollectionVec->getNumberOfElements()) ||
	 (inputCollectionVec->getNumberOfElements() != _noOfDetector) ) {
      stringstream ss;
      ss << "Input data and pedestal are incompatible" << endl
	 << "Input collection has    " << inputCollectionVec->getNumberOfElements()    << " detectors," << endl
	 << "Pedestal collection has " << pedestalCollectionVec->getNumberOfElements() << " detectors," << endl
	 << "The expected number is  " << _noOfDetector << endl;
      throw IncompatibleDataSetException(ss.str());
    }


    for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {

      TrackerRawDataImpl * rawData  = dynamic_cast < TrackerRawDataImpl * >(inputCollectionVec->getElementAt(iDetector));
      TrackerDataImpl    * pedestal = dynamic_cast < TrackerDataImpl * >   (pedestalCollectionVec->getElementAt(iDetector));
      
      if (rawData->getADCValues().size() != pedestal->getChargeValues().size()){
	stringstream ss;
	ss << "Input data and pedestal are incompatible" << endl
	   << "Detector " << iDetector << " has " <<  rawData->getADCValues().size() << " pixels in the input data " << endl
	   << "while " << pedestal->getChargeValues().size() << " in the pedestal data " << endl;
	throw IncompatibleDataSetException(ss.str());
      }

      // in the case MARLIN_USE_AIDA and the user wants to fill detector
      // debug histogram, this is the right place to book them. As usual
      // histograms are grouped in directory identifying the detector

#ifdef MARLIN_USE_AIDA
      string basePath, tempHistoName;
      {
	stringstream ss;
	ss << "detector-" << iDetector << "/";
	basePath = ss.str();
      }

      if ( ( _fillDebugHisto == 1 ) ||
	   ( _doCommonMode   == 1 ) ) {
	// it was changed from mkdir to mkdirs to be sure all
	// intermediate folders were properly created, but then I had
	// to come back to mkdir because this is the only supported in RAIDA.
	AIDAProcessor::tree(this)->mkdir(basePath.c_str());
      }

      if (_fillDebugHisto == 1) {
	// book the raw data histogram
	const int    rawDataDistHistoNBin = 4096;
	const double rawDataDistHistoMin  = -2048.5;
	const double rawDataDistHistoMax  = rawDataDistHistoMin + rawDataDistHistoNBin;
	{
	  stringstream ss;
	  ss << _rawDataDistHistoName << "-d" << iDetector;
	  tempHistoName = ss.str();
	}
	AIDA::IHistogram1D * rawDataDistHisto = 
	  AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
								    rawDataDistHistoNBin, rawDataDistHistoMin, rawDataDistHistoMax);
	_aidaHistoMap.insert(make_pair(tempHistoName, rawDataDistHisto));
	rawDataDistHisto->setTitle("Raw data distribution");

	// book the pedestal corrected data histogram
	const int    dataDistHistoNBin =  1000;
	const double dataDistHistoMin  = -100.;
	const double dataDistHistoMax  =  100.;
	{
	  stringstream ss;
	  ss << _dataDistHistoName << "-d" << iDetector;
	  tempHistoName = ss.str();
	}
	AIDA::IHistogram1D * dataDistHisto = 
	  AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
								    dataDistHistoNBin, dataDistHistoMin, dataDistHistoMax);
	_aidaHistoMap.insert(make_pair(tempHistoName, dataDistHisto));
	dataDistHisto->setTitle("Data (pedestal sub'ed) distribution");

      }

      if ( _doCommonMode == 1 ) {
	// book the common mode histo
	const int    commonModeDistHistoNBin = 100;
	const double commonModeDistHistoMin  = -10;
	const double commonModeDistHistoMax  = +10;
 	{
 	  stringstream ss;
 	  ss << _commonModeDistHistoName << "-d" << iDetector;
 	  tempHistoName = ss.str();
 	} 
 	AIDA::IHistogram1D * commonModeDistHisto = 
 	  AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), 
 								    commonModeDistHistoNBin, commonModeDistHistoMin, 
 								    commonModeDistHistoMax);
	_aidaHistoMap.insert(make_pair(tempHistoName, commonModeDistHisto));
	commonModeDistHisto->setTitle("Common mode distribution");
      }
#endif    
    }
    _isFirstEvent = false;
  }

  LCCollectionVec * correctedDataCollection = new LCCollectionVec(LCIO::TRACKERDATA);

  for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {

    // reset quantity for the common mode.
    double pixelSum      = 0.;
    double commonMode    = 0.;
    int    goodPixel     = 0;
    int    skippedPixel  = 0;


    TrackerRawDataImpl  * rawData   = dynamic_cast < TrackerRawDataImpl * >(inputCollectionVec->getElementAt(iDetector));
    TrackerDataImpl     * pedestal  = dynamic_cast < TrackerDataImpl * >   (pedestalCollectionVec->getElementAt(iDetector));
    TrackerDataImpl     * noise     = dynamic_cast < TrackerDataImpl * >   (noiseCollectionVec->getElementAt(iDetector));
    TrackerRawDataImpl  * status    = dynamic_cast < TrackerRawDataImpl * >(statusCollectionVec->getElementAt(iDetector));

    TrackerDataImpl     * corrected = new TrackerDataImpl;
    CellIDEncoder<TrackerDataImpl> idDataEncoder(EUTELESCOPE::MATRIXDEFAULTENCODING, correctedDataCollection);
    // for some unknown reason I cannot do something like
    // idDataEncoder["sensorID"] = cellDecoder(rawData)["sensorID"];
    // because it causes a compilation error, ask Frank for an
    // explanation because I'm not able to understand why it is
    // complaining
    int sensorID              = cellDecoder(rawData)["sensorID"];
    int xmin                  = cellDecoder(rawData)["xMin"];
    int xmax                  = cellDecoder(rawData)["xMax"];
    int ymin                  = cellDecoder(rawData)["yMin"];
    int ymax                  = cellDecoder(rawData)["yMax"];
    idDataEncoder["sensorID"] = sensorID;
    idDataEncoder["xMin"]     = xmin;
    idDataEncoder["xMax"]     = xmax;
    idDataEncoder["yMin"]     = ymin;
    idDataEncoder["yMax"]     = ymax;
    idDataEncoder.setCellID(corrected);

    ShortVec::const_iterator rawIter     = rawData->getADCValues().begin();
    FloatVec::const_iterator pedIter     = pedestal->getChargeValues().begin();
    FloatVec::const_iterator noiseIter   = noise->getChargeValues().begin();
    ShortVec::const_iterator statusIter  = status->getADCValues().begin();

    if ( _doCommonMode == 1 ) {
      while ( rawIter != rawData->getADCValues().end() ) {
	bool isHit   = ( ((*rawIter) - (*pedIter)) > _hitRejectionCut * (*noiseIter) );
	bool isGood  = ( (*statusIter) == EUTELESCOPE::GOODPIXEL );
	
	if ( !isHit && isGood ) {
	  pixelSum += (*rawIter) - (*pedIter);
	  ++goodPixel;
	} else if ( isHit ) {
	  ++skippedPixel;
	}
	++rawIter;
	++pedIter;
	++noiseIter;
	++statusIter;
      }
      
      if ( ( ( _maxNoOfRejectedPixels == -1 )  ||  ( skippedPixel < _maxNoOfRejectedPixels ) ) &&
	   ( goodPixel != 0 ) ) {
	
	commonMode = pixelSum / goodPixel;
#ifdef MARLIN_USE_AIDA
	string tempHistoName;
	{
	  stringstream ss;
	  ss << _commonModeDistHistoName << "-d" << iDetector;
	  tempHistoName = ss.str();
	}
	(dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName]))->fill(commonMode);
#endif

      } else {
	message<WARNING> ( log() << "Skipping event " << _iEvt << " because of common mode limit exceeded " );
	++_iEvt;
	throw SkipEventException(this);
      }
    }
  

    rawIter     = rawData->getADCValues().begin();
    pedIter     = pedestal->getChargeValues().begin();
    
    while ( rawIter != rawData->getADCValues().end() )  {
      double correctedValue = (*rawIter) - (*pedIter) - commonMode;
      corrected->chargeValues().push_back(correctedValue);
      

#ifdef MARLIN_USE_AIDA
      if (_fillDebugHisto == 1) {
	string tempHistoName;
	{
	  stringstream ss;
	  ss << _rawDataDistHistoName << "-d" << iDetector;
	  tempHistoName = ss.str();
	}
	(dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName]))->fill(*rawIter);
	{
	  stringstream ss;
	  ss << _dataDistHistoName << "-d" << iDetector;
	  tempHistoName = ss.str();
	}
	(dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName]))->fill(correctedValue);
      }
#endif
      ++rawIter;
      ++pedIter;
    }
    correctedDataCollection->push_back(corrected);
  }
  evt->addCollection(correctedDataCollection, _calibratedDataCollectionName);
  
  
  ++_iEvt;
  
}
  


void EUTelCalibrateEventProcessor::check (LCEvent * evt) {
  // nothing to check here - could be used to fill check plots in reconstruction processor
}


void EUTelCalibrateEventProcessor::end() {
  message<MESSAGE> ( "Successfully finished" );

}

