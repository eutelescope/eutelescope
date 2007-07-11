// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelCalibrateEventProcessor.cc,v 1.12 2007-07-11 06:54:16 bulgheroni Exp $
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
#include "EUTelHistogramManager.h"

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
    "EUTelCalibrateEventProcessor subtracts the pedestal value from the input data";

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
			      _fillDebugHisto, static_cast<bool> (false));
  
  registerProcessorParameter ("PerformCommonMode",
			      "Flag to switch on (1) or off (0) the common mode suppression algorithm",
			      _doCommonMode, static_cast<bool> (true));

  registerProcessorParameter ("HitRejectionCut",
			      "Threshold of pixel SNR for hit rejection",
			      _hitRejectionCut, static_cast<float> (3.5));

  registerProcessorParameter ("MaxNoOfRejectedPixels",
			      "Maximum allowed number of rejected pixel per event",
			      _maxNoOfRejectedPixels, static_cast<int> (3000));

  registerProcessorParameter("HistoInfoFileName", "This is the name of the histogram information file",
			     _histoInfoFileName, string( "histoinfo.xml" ) );
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

  try {
    
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

	// prepare the histogram manager
	auto_ptr<EUTelHistogramManager> histoMgr( new EUTelHistogramManager( _histoInfoFileName ));
	EUTelHistogramInfo    * histoInfo;
	bool                    isHistoManagerAvailable;

	try {
	  isHistoManagerAvailable = histoMgr->init();
	} catch ( ios::failure& e) {
	  message<ERROR> ( log() << "I/O problem with " << _histoInfoFileName << "\n"
			   << "Continuing without histogram manager"    );
	  isHistoManagerAvailable = false;
	} catch ( ParseException& e ) {
	  message<ERROR> ( log() << e.what() << "\n"
			   << "Continuing without histogram manager" );
	  isHistoManagerAvailable = false;
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
	  {
	    stringstream ss;
	    ss << _rawDataDistHistoName << "-d" << iDetector;
	    tempHistoName = ss.str();
	  }
	  int    rawDataDistHistoNBin = 4096;
	  double rawDataDistHistoMin  = -2048.5;
	  double rawDataDistHistoMax  = rawDataDistHistoMin + rawDataDistHistoNBin;
	  string rawDataDistTitle     = "Raw data distribution"; 
	  if ( isHistoManagerAvailable ) {
	    histoInfo = histoMgr->getHistogramInfo(_rawDataDistHistoName);
	    if ( histoInfo ) {
	      message<DEBUG> ( log() << (* histoInfo ) );
	      rawDataDistHistoNBin = histoInfo->_xBin;
	      rawDataDistHistoMin  = histoInfo->_xMin;
	      rawDataDistHistoMax  = histoInfo->_xMax;
	      if ( histoInfo->_title != "" ) rawDataDistTitle = histoInfo->_title;
	    }
	  }	  
	  AIDA::IHistogram1D * rawDataDistHisto = 
	    AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
								      rawDataDistHistoNBin, rawDataDistHistoMin, rawDataDistHistoMax);
	  if ( rawDataDistHisto ) {
	    _aidaHistoMap.insert(make_pair(tempHistoName, rawDataDistHisto));
	    rawDataDistHisto->setTitle(rawDataDistTitle.c_str());
	  } else {
	    message<ERROR> ( log() << "Problem booking the " << (basePath + tempHistoName) << ".\n"
			     << "Very likely a problem with path name. Switching off histogramming and continue w/o");
	    _fillDebugHisto = 0;
	  }

	  // book the pedestal corrected data histogram
	  {
	    stringstream ss;
	    ss << _dataDistHistoName << "-d" << iDetector;
	    tempHistoName = ss.str();
	  }
	  int    dataDistHistoNBin =  5000;
	  double dataDistHistoMin  = -500.;
	  double dataDistHistoMax  =  500.;
	  string dataDistTitle     = "Data (pedestal sub'ed) distribution";
	  if ( isHistoManagerAvailable ) {
	    histoInfo = histoMgr->getHistogramInfo(_dataDistHistoName);
	    if ( histoInfo ) {
	      message<DEBUG> ( log() << (* histoInfo ) );
	      dataDistHistoNBin = histoInfo->_xBin;
	      dataDistHistoMin  = histoInfo->_xMin;
	      dataDistHistoMax  = histoInfo->_xMax;
	      if ( histoInfo->_title != "" ) dataDistTitle = histoInfo->_title;
	    }
	  }	  
	  AIDA::IHistogram1D * dataDistHisto = 
	    AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
								      dataDistHistoNBin, dataDistHistoMin, dataDistHistoMax);
	  if ( dataDistHisto ) {
	    _aidaHistoMap.insert(make_pair(tempHistoName, dataDistHisto));
	    dataDistHisto->setTitle(dataDistTitle.c_str());
	  } else {
	    message<ERROR> ( log() << "Problem booking the " << (basePath + tempHistoName) << ".\n"
			     << "Very likely a problem with path name. Switching off histogramming and continue w/o");
	    _fillDebugHisto = 0;
	  }

	}

	if ( _doCommonMode == 1 ) {
	  // book the common mode histo
	  {
	    stringstream ss;
	    ss << _commonModeDistHistoName << "-d" << iDetector;
	    tempHistoName = ss.str();
	  } 
	  int    commonModeDistHistoNBin = 100;
	  double commonModeDistHistoMin  = -10;
	  double commonModeDistHistoMax  = +10;
	  string commonModeTitle         = "Common mode distribution";
	  if ( isHistoManagerAvailable ) {
	    histoInfo = histoMgr->getHistogramInfo(_commonModeDistHistoName);
	    if ( histoInfo ) {
	      message<DEBUG> ( log() << (* histoInfo ) );
	      commonModeDistHistoNBin = histoInfo->_xBin;
	      commonModeDistHistoMin  = histoInfo->_xMin;
	      commonModeDistHistoMax  = histoInfo->_xMax;
	      if ( histoInfo->_title != "" ) commonModeTitle = histoInfo->_title;
	    }
	  }	  
	  AIDA::IHistogram1D * commonModeDistHisto = 
	    AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), 
								      commonModeDistHistoNBin, commonModeDistHistoMin, 
								      commonModeDistHistoMax);
	  if ( commonModeDistHisto ) {
	    _aidaHistoMap.insert(make_pair(tempHistoName, commonModeDistHisto));
	    commonModeDistHisto->setTitle(commonModeTitle.c_str());
	  } else {
	    message<ERROR> ( log() << "Problem booking the " << (basePath + tempHistoName) << ".\n"
			     << "Very likely a problem with path name. Not filling common mode histograms" );
	    // setting _fillDebugHisto to 0 is not solving the problem
	    // because the common mode histogram is independent from
	    // the debug histogram. And disabling the common mode
	    // calculation is probably too much. The simplest solution
	    // would be to add another switch to the class data
	    // members. Since the fill() method is protected by a
	    // previous check on the retrieved histogram pointer,
	    // continue should be safe.
	  }
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
	  if ( AIDA::IHistogram1D* histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName]) )
	    histo->fill(commonMode);
#endif
	
	} else {
	  message<WARNING> ( log() << "Skipping event " << _iEvt 
			     << " because of common mode limit exceeded (" 
			     << skippedPixel << ")");
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
	  if ( AIDA::IHistogram1D* histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName]) )
	    histo->fill(*rawIter);
	  else {
	    message<ERROR> ( log() << "Not able to retrieve histogram pointer for " << tempHistoName 
			     << ".\nDisabling histogramming from now on " );
	    _fillDebugHisto = 0 ;
	  }

	  {
	    stringstream ss;
	    ss << _dataDistHistoName << "-d" << iDetector;
	    tempHistoName = ss.str();
	  }
	  if ( AIDA::IHistogram1D * histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName]) ) 
	    histo->fill(correctedValue);
	  else {
	    message<ERROR> ( log() << "Not able to retrieve histogram pointer for " << tempHistoName 
			     << ".\nDisabling histogramming from now on " );
	    _fillDebugHisto = 0 ;
	  }
	}
#endif
	++rawIter;
	++pedIter;
      }
      correctedDataCollection->push_back(corrected);
    }
    evt->addCollection(correctedDataCollection, _calibratedDataCollectionName);
  
  
    ++_iEvt;
  } catch (DataNotAvailableException& e) {
    message<ERROR> ( log() << e.what() << "\n"
		     << "Skipping this event " );
    throw SkipEventException(this);
  }
  
}
  


void EUTelCalibrateEventProcessor::check (LCEvent * evt) {
  // nothing to check here - could be used to fill check plots in reconstruction processor
}


void EUTelCalibrateEventProcessor::end() {
  message<MESSAGE> ( "Successfully finished" );

}

