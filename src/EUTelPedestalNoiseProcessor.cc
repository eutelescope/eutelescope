// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelPedestalNoiseProcessor.cc,v 1.9 2007-02-17 13:37:14 bulgheroni Exp $
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
#include "EUTelRunHeaderImpl.h"
#include "EUTelPedestalNoiseProcessor.h"
#include "EUTELESCOPE.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"

#ifdef MARLIN_USE_AIDA
// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile2D.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h> 
#include <lcio.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTime.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h> 

// system includes <>
#include <string>
#include <sstream>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;


// definition of static members mainly used to name histograms
#ifdef MARLIN_USE_AIDA
std::string EUTelPedestalNoiseProcessor::_pedeDistHistoName   = "PedeDist";
std::string EUTelPedestalNoiseProcessor::_noiseDistHistoName  = "NoiseDist";
std::string EUTelPedestalNoiseProcessor::_commonModeHistoName = "CommonMode";
std::string EUTelPedestalNoiseProcessor::_pedeMapHistoName    = "PedeMap";
std::string EUTelPedestalNoiseProcessor::_noiseMapHistoName   = "NoiseMap";
std::string EUTelPedestalNoiseProcessor::_statusMapHistoName  = "StatusMap";
std::string EUTelPedestalNoiseProcessor::_tempProfile2DName   = "TempProfile2D";
#endif

EUTelPedestalNoiseProcessor::EUTelPedestalNoiseProcessor () :Processor("EUTelPedestalNoiseProcessor") {

  // modify processor description
  _description =
    "EUTelPedestalNoiseProcessor computes the pedestal and noise values of a pixel detector";

  // first of all we need to register the input collection
  registerInputCollection (LCIO::TRACKERRAWDATA, "RawDataCollectionName",
			   "Input raw data collection",
			   _rawDataCollectionName, string ("rawdata"));

  // register compulsory parameters
  registerProcessorParameter ("CalculationAlgorithm",
			      "Select the algorithm for pede/noise calculation",
			      _pedestalAlgo,
			      string (EUTELESCOPE::MEANRMS));
  registerProcessorParameter("BadPixelMaskingAlgorithm",
			     "Select the algorithm for bad pixel masking",
			     _badPixelAlgo,
			     string (EUTELESCOPE::NOISEDISTRIBUTION));
  registerProcessorParameter ("NoOfCMIteration",
			      "Number of common mode suppression iterations",
			      _noOfCMIterations, static_cast < int >(1));
  registerProcessorParameter ("HitRejectionCut",
			      "Threshold for rejection of hit pixel (SNR units)",
			      _hitRejectionCut, static_cast < float >(4));
  registerProcessorParameter ("MaxNoOfRejectedPixels",
			      "Maximum allowed number of rejected pixels per event",
			      _maxNoOfRejectedPixels,
			      static_cast < int >(1000));
  registerProcessorParameter ("BadPixelMaskCut",
			      "Threshold for bad pixel identification",
			      _badPixelMaskCut, static_cast < float >(3.5));
  registerProcessorParameter ("FirstEvent", 
			      "First event for pedestal calculation",
			      _firstEvent, static_cast < int > (0));
  registerProcessorParameter ("LastEvent",
			      "Last event for pedestal calculation",
			      _lastEvent, static_cast < int > (-1));
  registerProcessorParameter ("OutputPedeFile","The filename (w/o .slcio) to store the pedestal file",
			      _outputPedeFileName , string("outputpede")); 

  // now the optional parameters
  registerOptionalParameter ("PedestalCollectionName",
			     "Pedestal collection name",
			     _pedestalCollectionName, string ("pedestal"));
  registerOptionalParameter ("NoiseCollectionName",
			     "Noise collection name", _noiseCollectionName,
			     string ("noise"));
  registerOptionalParameter ("StatusCollectionName",
			     "Status collection name",
			     _statusCollectionName, string ("status"));

}


void EUTelPedestalNoiseProcessor::init () {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

  // set the pedestal flag to true and the loop counter to zero
  _doPedestal = true;
  _iLoop = 0;

  if ( _pedestalAlgo == EUTELESCOPE::MEANRMS ) {
    // reset the temporary arrays
    _tempPede.clear ();
    _tempNoise.clear ();
    _tempEntries.clear ();
  }

#ifndef MARLIN_USE_AIDA
  if ( _pedestalAlgo == EUTELESCOPE::AIDAPROFILE ) {
    cerr << "[" << name() << "] The " << EUTELESCOPE::AIDAPROFILE 
	 << " algorithm cannot be applied since Marlin is not using AIDA" << endl
	 << " Algorithm changed to " << EUTELESCOPE::MEANRMS << endl;
    _pedestalAlgo = EUTELESCOPE::MEANRMS;
  }
#endif
  

  // reset all the final arrays
  _pedestal.clear();
  _noise.clear();
  _status.clear();

  // clean up things related to CM 
  _noOfSkippedEvent = 0;

}

void EUTelPedestalNoiseProcessor::processRunHeader (LCRunHeader * rdr) {

  _detectorName = rdr->getDetectorName();

  // to make things easier re-cast the input header to the EUTelRunHeaderImpl
  EUTelRunHeaderImpl *  runHeader = static_cast<EUTelRunHeaderImpl*>(rdr);

  // increment the run counter
  ++_iRun;

  // let me get from the run header all the available parameter
  _noOfDetector = runHeader->getNoOfDetector();

  // now the four vectors containing the first and the last pixel
  // along both the directions
  _minX = runHeader->getMinX();
  _maxX = runHeader->getMaxX();
  _minY = runHeader->getMinY();
  _maxY = runHeader->getMaxY();

  // make some test on parameters
  if ( ( _pedestalAlgo != EUTELESCOPE::MEANRMS ) &&
       ( _pedestalAlgo != EUTELESCOPE::AIDAPROFILE) 
       ) {
    throw InvalidParameterException(string("_pedestalAlgo cannot be " + _pedestalAlgo));
  }

  if ( ( _badPixelAlgo != EUTELESCOPE::NOISEDISTRIBUTION ) &&
       ( _badPixelAlgo != EUTELESCOPE::ABSOLUTENOISEVALUE )
       ) {
    throw InvalidParameterException(string("_badPixelAlgo cannot be " + _badPixelAlgo));
  }

  // get from the file header the event number
  int tempNoOfEvent = runHeader->getNoOfEvent();
  
  if (tempNoOfEvent == 0) {
    stringstream ss;
    ss << "NoOfEvent cannot be " << tempNoOfEvent;
    throw InvalidParameterException(ss.str());
  }

  if ( _lastEvent == -1) _lastEvent = tempNoOfEvent;
  if ( _lastEvent >= tempNoOfEvent) _lastEvent = tempNoOfEvent;

  if ( _iLoop == 0 ) {
    // write the current header to the output condition file
    LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();

    try {
      lcWriter->open(_outputPedeFileName, LCIO::WRITE_NEW);
    } catch (IOException& e) {
      cerr << e.what() << endl;
      return;
    }

    EUTelRunHeaderImpl * newHeader = new EUTelRunHeaderImpl;
    
    newHeader->setRunNumber(runHeader->getRunNumber());
    newHeader->setDetectorName(runHeader->getDetectorName());
    newHeader->setHeaderVersion(runHeader->getHeaderVersion());
    newHeader->setDataType(runHeader->getDataType());
    newHeader->setDateTime();
    newHeader->setDAQHWName(runHeader->getDAQHWName());
    newHeader->setDAQHWVersion(runHeader->getDAQHWVersion());
    newHeader->setDAQSWName(runHeader->getDAQSWName());
    newHeader->setDAQSWVersion(runHeader->getDAQSWVersion());  
    newHeader->setNoOfEvent(runHeader->getNoOfEvent());
    newHeader->setNoOfDetector(runHeader->getNoOfDetector());
    newHeader->setMinX(runHeader->getMinY());
    newHeader->setMaxX(runHeader->getMaxX());
    newHeader->setMinY(runHeader->getMinY());
    newHeader->setMaxY(runHeader->getMaxY());
    newHeader->addProcessor(name());
    
    lcWriter->writeRunHeader(newHeader);
    lcWriter->close();

    // also book histos
    bookHistos();
  }
}


void EUTelPedestalNoiseProcessor::processEvent (LCEvent * evt) {

  if ( _iLoop == 0 ) firstLoop(evt);
  else otherLoop(evt);

}
  


void EUTelPedestalNoiseProcessor::check (LCEvent * evt) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void EUTelPedestalNoiseProcessor::end() {
  cout << "[" << name() <<"] Successfully finished" << endl;
}

void EUTelPedestalNoiseProcessor::fillHistos() {
  
#ifdef MARLIN_USE_AIDA
  cout << "[" << name() << "] Filling final histograms " << endl;

  string tempHistoName;
  for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
    int iPixel = 0;
    for (int yPixel = _minY[iDetector]; yPixel <= _maxY[iDetector]; yPixel++) {
      for (int xPixel = _minX[iDetector]; xPixel <= _maxX[iDetector]; xPixel++) {
	{
	  stringstream ss;
	  ss << _statusMapHistoName << "-d" << iDetector << "-l" << _iLoop;
	  tempHistoName = ss.str();
	}
	(dynamic_cast<AIDA::IHistogram2D*>(_aidaHistoMap[tempHistoName]))->fill(static_cast<double>(xPixel), static_cast<double>(yPixel),
										static_cast<double> (_status[iDetector][iPixel]));

	if ( _status[iDetector][iPixel] == EUTELESCOPE::GOODPIXEL) {
	  {
	    stringstream ss;
	    ss << _pedeDistHistoName << "-d" << iDetector << "-l" << _iLoop;
	    tempHistoName = ss.str();
	  }
	  (dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(_pedestal[iDetector][iPixel]);
	  {
	    stringstream ss;
	    ss << _noiseDistHistoName << "-d" << iDetector << "-l" << _iLoop;
	    tempHistoName = ss.str();
	  }
	  (dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(_noise[iDetector][iPixel]);
	  {
	    stringstream ss;
	    ss << _pedeMapHistoName << "-d" << iDetector << "-l" << _iLoop;
	    tempHistoName = ss.str();
	  }
	  (dynamic_cast<AIDA::IHistogram2D*>(_aidaHistoMap[tempHistoName]))->fill(static_cast<double>(xPixel), static_cast<double>(yPixel),
										   _pedestal[iDetector][iPixel]);
	  {
	    stringstream ss;
	    ss << _noiseMapHistoName << "-d" << iDetector << "-l" << _iLoop;
	    tempHistoName = ss.str();
	  } 
	  (dynamic_cast<AIDA::IHistogram2D*>(_aidaHistoMap[tempHistoName]))->fill(static_cast<double>(xPixel), static_cast<double>(yPixel),
										  _noise[iDetector][iPixel]);
	}
	++iPixel;
      }
    }
  }


#else
  cout << "[" << name() << "] No histogram produced because Marlin doesn't use AIDA " << endl;
#endif

}

void EUTelPedestalNoiseProcessor::maskBadPixel() {

  double threshold;
  int    badPixelCounter = 0;
  cout <<  "[" << name() << "] Masking bad pixels ";

  if ( _badPixelAlgo == EUTELESCOPE::NOISEDISTRIBUTION ) {
    // to do it we need to know the mean value and the RMS of the noise
    // vector.
    
    for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
      
      double sumw  = 0;
      double sumw2 = 0;
      double num   = 0;
      
      // begin a first loop on all pixel to calculate the masking threshold
      for (unsigned int iPixel = 0; iPixel < _status[iDetector].size(); iPixel++) {
	if ( _status[iDetector][iPixel] == EUTELESCOPE::GOODPIXEL ) {
	  sumw  += _noise[iDetector][iPixel];
	  sumw2 += pow(_noise[iDetector][iPixel],2);
	  ++num;
	}
      } 
      double meanw      = sumw  / num;
      double meanw2     = sumw2 / num;
      double rms        = sqrt( meanw2 - pow(meanw,2));
      threshold  = meanw + (rms * _badPixelMaskCut);
      
    }
  } else if ( _badPixelAlgo == EUTELESCOPE::ABSOLUTENOISEVALUE ) {
    threshold = _badPixelMaskCut;
  }


  // scan the noise vector again and apply the cut
  for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
    for (unsigned int iPixel = 0; iPixel < _status[iDetector].size(); iPixel++) {
      if ( ( _noise[iDetector][iPixel] > threshold ) && 
	   ( _status[iDetector][iPixel] == EUTELESCOPE::GOODPIXEL ) ) {
	_status[iDetector][iPixel] = EUTELESCOPE::BADPIXEL;
	++badPixelCounter;
      }
    }
  } // end loop on detector;
  cout << "(" << badPixelCounter << ")" << endl;
}


void EUTelPedestalNoiseProcessor::firstLoop(LCEvent * evt) {

  if (  _iEvt < _firstEvent )  {
    // pedestal calculation may occuring on a subset of events
    // outside this range, just skip the current event
    if ( ( _iEvt == 0 ) || ( _iEvt ==  _firstEvent - 1) || ( _iEvt%10 == 0) )
      cout << "[" << name() << "] Skipping event " << _iEvt << endl;
    ++_iEvt;
    throw SkipEventException(this);
  }

  if ( _iEvt > _lastEvent ) {
    cout << "[" << name() <<"] Last processed event " << _iEvt - 1 << endl;
    throw StopProcessingException(this);
  }
  
  // keep the user updated
  if ( _iEvt % 10 == 0 ) {
    cout << "[" << name() << "] Performing loop " << _iLoop << " on event: " << _iEvt << endl;
  }

  // let me get the rawDataCollection. This is should contain a TrackerRawDataObject
  // for each detector plane in the telescope.
  LCCollectionVec *collectionVec = dynamic_cast < LCCollectionVec * >(evt->getCollection (_rawDataCollectionName));

  
  if ( isFirstEvent() ) {
    // the collection contains several TrackerRawData

    for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
      // _tempPedestal, _tempNoise, _tempEntries are vector of vector.
      // they have been already cleared in the init() method we are
      // already looping on detectors, so we just need to push back a
      // vector empty for each cycle
      // 
      // _tempPedestal should be initialized with the adcValues, while
      // _tempNoise and _tempEntries must be initialized to zero. Since
      // adcValues is a vector of shorts, we need to copy each
      // elements into _tempPedestal with a suitable re-casting

      // get the TrackerRawData object from the collection for this detector

      TrackerRawData *trackerRawData = dynamic_cast < TrackerRawData * >(collectionVec->getElementAt (iDetector));
      ShortVec adcValues = trackerRawData->getADCValues ();

      if ( _pedestalAlgo == EUTELESCOPE::MEANRMS ) {
	// in the case of MEANRMS we have to deal with the standard
	// vectors
	ShortVec::iterator iter = adcValues.begin();
	FloatVec tempDoubleVec;
	while ( iter != adcValues.end() ) {
	  tempDoubleVec.push_back( static_cast< double > (*iter));
	  ++iter;
	}
	_tempPede.push_back(tempDoubleVec);
	
	// initialize _tempNoise and _tempEntries with all zeros and
	// ones
	_tempNoise.push_back(FloatVec(adcValues.size(), 0.));
	_tempEntries.push_back(IntVec(adcValues.size(), 1));

      } else if ( _pedestalAlgo == EUTELESCOPE::AIDAPROFILE ) {
#ifdef MARLIN_USE_AIDA
	// in the case of AIDAPROFILE we don't need any vectors since
	// everything is done by the IProfile2D automatically
	int iPixel = 0;
	stringstream ss;
	ss << _tempProfile2DName << "-d" << iDetector;
	for (int yPixel = _minY[iDetector]; yPixel <= _maxY[iDetector]; yPixel++) {
	  for (int xPixel = _minX[iDetector]; xPixel <= _maxX[iDetector]; xPixel++) {
	    double temp = static_cast<double> (adcValues[iPixel]);
	    (dynamic_cast<AIDA::IProfile2D*> (_aidaHistoMap[ss.str()]))
	      ->fill(static_cast<double> (xPixel), static_cast<double> (yPixel), temp);
	    ++iPixel;
	  }
	}
#endif
      }

      // the status vector can be initialize as well with all
      // GOODPIXEL
      _status.push_back(ShortVec(adcValues.size(), EUTELESCOPE::GOODPIXEL));
    } // end of detector loop

    // nothing else to do in the first event
    _isFirstEvent = false;

  } else {     // end of first event

    // after the firstEvent all temp vectors and the status one have
    // the correct number of entries for both indeces
    // loop on the detectors
    for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {

      // get the TrackerRawData object from the collection for this plane
      TrackerRawData *trackerRawData = dynamic_cast < TrackerRawData * >(collectionVec->getElementAt (iDetector));
      ShortVec adcValues = trackerRawData->getADCValues ();
      

      if ( _pedestalAlgo == EUTELESCOPE::MEANRMS ) {

	// start looping on all pixels
	int iPixel = 0;
	for (int yPixel = _minY[iDetector]; yPixel <= _maxY[iDetector]; yPixel++) {
	  for (int xPixel = _minX[iDetector]; xPixel <= _maxX[iDetector]; xPixel++) {
	    _tempEntries[iDetector][iPixel] =  _tempEntries[iDetector][iPixel] + 1;
	    _tempPede[iDetector][iPixel]    = ((_tempEntries[iDetector][iPixel] - 1) * _tempPede[iDetector][iPixel]
					       + adcValues[iPixel]) / _tempEntries[iDetector][iPixel];
	    _tempNoise[iDetector][iPixel]   = sqrt(((_tempEntries[iDetector][iPixel] - 1) * pow(_tempNoise[iDetector][iPixel],2) 
						    + pow(adcValues[iPixel] - _tempPede[iDetector][iPixel], 2)) / 
						   _tempEntries[iDetector][iPixel]);
	    ++iPixel;
	  } // end loop on xPixel
	}   // end loop on yPixel	  
	
      } else if ( _pedestalAlgo == EUTELESCOPE::AIDAPROFILE ) {

#ifdef MARLIN_USE_AIDA	
	stringstream ss;
	ss << _tempProfile2DName << "-d" << iDetector;

	int iPixel = 0;
	for (int yPixel = _minY[iDetector]; yPixel <= _maxY[iDetector]; yPixel++) {
	  for (int xPixel = _minX[iDetector]; xPixel <= _maxX[iDetector]; xPixel++) {
	    (dynamic_cast<AIDA::IProfile2D*> (_aidaHistoMap[ss.str()]))
	      ->fill(static_cast<double> (xPixel), static_cast<double> (yPixel), static_cast<double> (adcValues[iPixel]));
	    ++iPixel;
	  }
	}
#endif
      }

    }     // end loop on detectors
  }

  // increment the event counter
  ++_iEvt;

  if (isLastEvent()) finalizeProcessor();

}

void EUTelPedestalNoiseProcessor::otherLoop(LCEvent * evt) {
  
  if (  _iEvt < _firstEvent  ) {
    // pedestal calculation may occuring on a subset of events
    // outside this range, just skip the current event
    if ( ( _iEvt == 0 ) || ( _iEvt ==  _firstEvent - 1) || ( _iEvt%10 == 0) )
      cout << "[" << name() << "] Skipping event " << _iEvt << endl;
    ++_iEvt;
    throw SkipEventException(this);
  }

  if ( _iEvt > _lastEvent ) {
    cout << "[" << name() <<"] Last processed event " << _iEvt - 1 << endl;
    throw StopProcessingException(this);
  }

  // keep the user updated
  if ( _iEvt % 10 == 0 ) {
    cout << "[" << name() <<"] Performing loop " << _iLoop << " on event: " << _iEvt << endl;
  }

  // let me get the rawDataCollection. This is should contain a TrackerRawDataObject
  // for each detector plane in the telescope.
  LCCollectionVec *collectionVec = dynamic_cast < LCCollectionVec * >(evt->getCollection (_rawDataCollectionName));
  
  for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
    // get the TrackerRawData object from the collection for this detector
    TrackerRawData *trackerRawData = dynamic_cast < TrackerRawData * >(collectionVec->getElementAt (iDetector));
    ShortVec adcValues = trackerRawData->getADCValues ();
    
    // prepare stuff for common mode calculation: pixelSum is the
    // sum of all good pixel signals. Pixels are identify as good if
    // their status is good and it they are not recognized as hit
    // pixel by the hit rejection threshold. goodPixel is the number
    // of good pixel in this detector used for common mode
    // calculation. commonMode is, indeed, the mean value of the
    // pixel signals pedestal sub'ed. iPixel is a pixel counter
    double pixelSum     = 0.;
    double commonMode   = 0.;
    int    goodPixel    = 0;
    int    skippedPixel = 0;
    int    iPixel       = 0;
    
    // start looping on all pixels for hit rejection
    for (int yPixel = _minY[iDetector]; yPixel <= _maxY[iDetector]; yPixel++) {
      for (int xPixel = _minX[iDetector]; xPixel <= _maxX[iDetector]; xPixel++) {
	bool isHit  = ( ( adcValues[iPixel] - _pedestal[iDetector][iPixel] ) > _hitRejectionCut * _noise[iDetector][iPixel] );
	bool isGood = ( _status[iDetector][iPixel] == EUTELESCOPE::GOODPIXEL );
	if ( !isHit && isGood ) {
	  pixelSum += adcValues[iPixel] - _pedestal[iDetector][iPixel];
	  ++goodPixel;
	} else if ( isHit ) {
	  ++skippedPixel;
	}
	++iPixel;
      } 
    }
    
    if ( ( skippedPixel < _maxNoOfRejectedPixels ) &&
	 ( goodPixel != 0 ) ) {

      commonMode = pixelSum / goodPixel;
#ifdef MARLIN_USE_AIDA      
      stringstream ss;
      ss << _commonModeHistoName << "-d" << iDetector << "-l" << _iLoop;
      (dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[ss.str()]))->fill(commonMode);
#endif
      iPixel = 0;
      for (int yPixel = _minY[iDetector]; yPixel <= _maxY[iDetector]; yPixel++) {
	for (int xPixel = _minX[iDetector]; xPixel <= _maxX[iDetector]; xPixel++) {
	  if ( _status[iDetector][iPixel] == EUTELESCOPE::GOODPIXEL ) {
	    double pedeCorrected = adcValues[iPixel] - commonMode;
	    if ( _pedestalAlgo == EUTELESCOPE::MEANRMS ) {

	      _tempEntries[iDetector][iPixel] = _tempEntries[iDetector][iPixel] + 1;
	      _tempPede[iDetector][iPixel]    = ((_tempEntries[iDetector][iPixel] - 1) * _tempPede[iDetector][iPixel]
						 + pedeCorrected) / _tempEntries[iDetector][iPixel];
	      _tempNoise[iDetector][iPixel]   = sqrt(((_tempEntries[iDetector][iPixel] - 1) * pow(_tempNoise[iDetector][iPixel],2) 
						      + pow(pedeCorrected - _tempPede[iDetector][iPixel], 2)) / 
						     _tempEntries[iDetector][iPixel]);

	    } else if ( _pedestalAlgo == EUTELESCOPE::AIDAPROFILE) {
#ifdef MARLIN_USE_AIDA	      
	      stringstream ss;
	      ss << _tempProfile2DName << "-d" << iDetector;
	      (dynamic_cast<AIDA::IProfile2D*> (_aidaHistoMap[ss.str()]))
		->fill(static_cast<double> (xPixel), static_cast<double> (yPixel), pedeCorrected);
#endif
	    }
	  }
	  ++iPixel;
	}
      }
    } else {
      cout << "Skipping event " << _iEvt << " because of max number of rejected pixels exceeded. (" << skippedPixel << ")" << endl;
      ++_noOfSkippedEvent;
    }
  }	
  ++_iEvt;
  if (isLastEvent()) finalizeProcessor();
}

void EUTelPedestalNoiseProcessor::bookHistos() {

#ifdef MARLIN_USE_AIDA
  // histograms are grouped in loops and detectors
  cout << "[" << name() << "] Booking histograms " << endl;


  string tempHistoName;

  // start looping on the number of loops. Remeber that we have one
  // loop more than the number of common mode iterations
  for (int iLoop = 0; iLoop < _noOfCMIterations + 1; iLoop++) {

    // prepare the name of the current loop directory and add it the
    // the current ITree
    string loopDirName;
    {
      stringstream ss;
      ss << "loop-" << iLoop;
      loopDirName = ss.str();
    }
    AIDAProcessor::tree(this)->mkdir(loopDirName.c_str());

    // start looping on detectors
    for (int iDetector = 0; iDetector < _noOfDetector;  iDetector++) {

      // prepare the name of the current detector and add it to the
      // current ITree inside the current loop folder
      string detectorDirName;
      {
	stringstream ss;
	ss << "detector-" << iDetector;
	detectorDirName = ss.str();
      }
      string basePath = loopDirName + "/" + detectorDirName + "/";
      AIDAProcessor::tree(this)->mkdir(basePath.c_str());
      
      // book an histogram for the edestal distribution
      const int    pedeDistHistoNBin   = 100; 
      const double pedeDistHistoMin    = -20.;
      const double pedeDistHistoMax    =  29.;
      {
	stringstream ss;
	ss << _pedeDistHistoName << "-d" << iDetector << "-l" << iLoop;
	tempHistoName = ss.str();
      } 
      AIDA::IHistogram1D * pedeDistHisto = 
	AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), 
 								  pedeDistHistoNBin, pedeDistHistoMin, pedeDistHistoMax);
      _aidaHistoMap.insert(make_pair(tempHistoName, pedeDistHisto));
      pedeDistHisto->setTitle("Pedestal distribution");

      // book an histogram for the noise distribution
      const int    noiseDistHistoNBin  =  100;
      const double noiseDistHistoMin   =  -5.;
      const double noiseDistHistoMax   =  15.;
      {
	stringstream ss;
	ss << _noiseDistHistoName << "-d" << iDetector << "-l" << iLoop;
	tempHistoName = ss.str();
      }
      AIDA::IHistogram1D * noiseDistHisto =
	AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
								  noiseDistHistoNBin, noiseDistHistoMin, noiseDistHistoMax);
      _aidaHistoMap.insert(make_pair(tempHistoName, noiseDistHisto));
      noiseDistHisto->setTitle("Noise distribution");

      // book a 1d histo for common mode only if loop >= 1
      if (iLoop >= 1) {
	const int    commonModeHistoNBin = 100;
	const double commonModeHistoMin  =  -2;
	const double commonModeHistoMax  =   2;
	{
	  stringstream ss;
	  ss << _commonModeHistoName << "-d" << iDetector << "-l" << iLoop;
	  tempHistoName = ss.str();
	}
	AIDA::IHistogram1D * commonModeHisto =
	  AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
								    commonModeHistoNBin, commonModeHistoMin, commonModeHistoMax);
	_aidaHistoMap.insert(make_pair(tempHistoName, commonModeHisto));
	commonModeHisto->setTitle("Common mode distribution");
      }

      // book a 2d histogram for pedestal map
      const int    xNoOfPixel = abs( _maxX[iDetector] - _minX[iDetector] + 1);
      const int    yNoOfPixel = abs( _maxY[iDetector] - _minY[iDetector] + 1);
      const double xMin       = _minX[iDetector] - 0.5;
      const double xMax       =  xMin + xNoOfPixel;
      const double yMin       = _minY[iDetector] - 0.5;
      const double yMax       =  yMin + yNoOfPixel;
      {
	stringstream ss;
	ss << _pedeMapHistoName << "-d" << iDetector << "-l" << iLoop;
	tempHistoName = ss.str();
      }
      AIDA::IHistogram2D * pedeMapHisto =
	AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(),
								  xNoOfPixel, xMin, xMax, yNoOfPixel, yMin, yMax);	
      _aidaHistoMap.insert(make_pair(tempHistoName, pedeMapHisto));
      pedeMapHisto->setTitle("Pedestal map");

      // book a 2d histogram for noise map
      {
	stringstream ss;
	ss << _noiseMapHistoName << "-d" << iDetector << "-l" << iLoop;
	tempHistoName = ss.str();
      }
      AIDA::IHistogram2D * noiseMapHisto =
	AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(),
								  xNoOfPixel, xMin, xMax, yNoOfPixel, yMin, yMax);	
      _aidaHistoMap.insert(make_pair(tempHistoName, noiseMapHisto));      
      noiseMapHisto->setTitle("Noise map");

      // book a 2d histogram for status map
      {
	stringstream ss;
	ss << _statusMapHistoName << "-d" << iDetector << "-l" << iLoop;
	tempHistoName = ss.str();
      }
      AIDA::IHistogram2D * statusMapHisto =
	AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(),
								  xNoOfPixel, xMin, xMax, yNoOfPixel, yMin, yMax);	
      _aidaHistoMap.insert(make_pair(tempHistoName, statusMapHisto));
      statusMapHisto->setTitle("Status map");

      if ( ( _pedestalAlgo == EUTELESCOPE::AIDAPROFILE ) &&
	   ( iLoop == 0 ) ) {
	// we just need to prepare such a 2d profile only in the case
	// we are using the AIDAPROFILE calculation algorithm and we
	// need just one copy of it for each detector.
	{
	  stringstream ss;
	  ss << _tempProfile2DName << "-d" << iDetector;
	  tempHistoName = ss.str();
	}
	AIDA::IProfile2D * tempProfile2D =
	  AIDAProcessor::histogramFactory(this)->createProfile2D( (basePath + tempHistoName).c_str(),
								  xNoOfPixel, xMin, xMax, yNoOfPixel, yMin, yMax,-1000,1000);
	_aidaHistoMap.insert(make_pair(tempHistoName, tempProfile2D));
	tempProfile2D->setTitle("Temp profile for pedestal calculation");
      }

    }
  } // end on iLoop
#endif // MARLIN_USE_AIDA
}

void EUTelPedestalNoiseProcessor::finalizeProcessor() {
  
  if ( _pedestalAlgo == EUTELESCOPE::MEANRMS ) {
    
    // the loop on events is over so we need to move temporary vectors
    // to final vectors
    _pedestal = _tempPede;
    _noise    = _tempNoise;
    
    // clear the temporary vectors
    _tempPede.clear();
    _tempNoise.clear();
    _tempEntries.clear();
    
  } else if ( _pedestalAlgo == EUTELESCOPE::AIDAPROFILE ) {
#ifdef MARLIN_USE_AIDA
    _pedestal.clear();
    _noise.clear();
    for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
      stringstream ss;
      ss << _tempProfile2DName << "-d" << iDetector;
      FloatVec tempPede;
      FloatVec tempNoise;
      for (int yPixel = _minY[iDetector]; yPixel <= _maxY[iDetector]; yPixel++) {
	for (int xPixel = _minX[iDetector]; xPixel <= _maxX[iDetector]; xPixel++) {
	  tempPede.push_back((float) (dynamic_cast<AIDA::IProfile2D*> (_aidaHistoMap[ss.str()]))->binHeight(xPixel,yPixel));
	  // WARNING: the noise part of this algorithm is still not
	  // working probably because of a bug in RAIDA implementation
	  tempNoise.push_back((float) (dynamic_cast<AIDA::IProfile2D*> (_aidaHistoMap[ss.str()]))->binRms(xPixel,yPixel));
	  //cout << xPixel << " " << yPixel << " " << tempPede.back() << " " << tempNoise.back() << endl;
	}
      }
      _pedestal.push_back(tempPede);
      _noise.push_back(tempNoise);
    }
#endif
  }
  
  // mask the bad pixels here
  maskBadPixel();
  
  // fill in the histograms
  fillHistos();
  
  // increment the loop counter
  ++_iLoop;
  
  // check if we need another loop or we can finish. Remember that we
  // have a total number of loop of _noOfCMIteration + 1
  if ( _iLoop == _noOfCMIterations + 1 ) {
    // ok this was last loop  

    cout << "[" << name() << "] Writing the output condition file" << endl;

    LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();
    
    try {
      lcWriter->open(_outputPedeFileName,LCIO::WRITE_APPEND);
    } catch (IOException& e) {
      cerr << e.what() << endl;
      return;
    }
    
    LCEventImpl * event = new LCEventImpl();
    event->setDetectorName(_detectorName);
    event->setRunNumber(_iRun);
    
    LCTime * now = new LCTime;
    event->setTimeStamp(now->timeStamp());
    delete now;


    LCCollectionVec * pedestalCollection = new LCCollectionVec(LCIO::TRACKERDATA);
    LCCollectionVec * noiseCollection    = new LCCollectionVec(LCIO::TRACKERDATA);
    LCCollectionVec * statusCollection   = new LCCollectionVec(LCIO::TRACKERRAWDATA);
    
    for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
      
      TrackerDataImpl    * pedestalMatrix = new TrackerDataImpl;
      TrackerDataImpl    * noiseMatrix    = new TrackerDataImpl;
      TrackerRawDataImpl * statusMatrix   = new TrackerRawDataImpl;
      
      CellIDEncoder<TrackerDataImpl>    idPedestalEncoder(EUTELESCOPE::MATRIXDEFAULTENCODING, pedestalCollection);
      CellIDEncoder<TrackerDataImpl>    idNoiseEncoder(EUTELESCOPE::MATRIXDEFAULTENCODING, noiseCollection);
      CellIDEncoder<TrackerRawDataImpl> idStatusEncoder(EUTELESCOPE::MATRIXDEFAULTENCODING, statusCollection);
      
      idPedestalEncoder["sensorID"] = iDetector;
      idNoiseEncoder["sensorID"]    = iDetector;
      idStatusEncoder["sensorID"]   = iDetector;
      idPedestalEncoder["xMin"]     = _minX[iDetector];
      idNoiseEncoder["xMin"]        = _minX[iDetector];
      idStatusEncoder["xMin"]       = _minX[iDetector];
      idPedestalEncoder["xMax"]     = _maxX[iDetector];
      idNoiseEncoder["xMax"]        = _maxX[iDetector];
      idStatusEncoder["xMax"]       = _maxX[iDetector];
      idPedestalEncoder["yMin"]     = _minY[iDetector];
      idNoiseEncoder["yMin"]        = _minY[iDetector];
      idStatusEncoder["yMin"]       = _minY[iDetector];
      idPedestalEncoder["yMax"]     = _maxY[iDetector];
      idNoiseEncoder["yMax"]        = _maxY[iDetector];
      idStatusEncoder["yMax"]       = _maxY[iDetector];
      idPedestalEncoder.setCellID(pedestalMatrix);
      idNoiseEncoder.setCellID(noiseMatrix);
      idStatusEncoder.setCellID(statusMatrix);
      
      pedestalMatrix->setChargeValues(_pedestal[iDetector]);
      noiseMatrix->setChargeValues(_noise[iDetector]);
      statusMatrix->setADCValues(_status[iDetector]);
      
      pedestalCollection->push_back(pedestalMatrix);
      noiseCollection->push_back(noiseMatrix);
      statusCollection->push_back(statusMatrix);
    }

    event->addCollection(pedestalCollection, _pedestalCollectionName);
    event->addCollection(noiseCollection, _noiseCollectionName);
    event->addCollection(statusCollection, _statusCollectionName);
    
    lcWriter->writeEvent(event);
    delete event;
    throw StopProcessingException(this);
    setReturnValue("IsPedestalFinished", true);
  } else {
    // now we need to loop again
    // so reset the event counter
    _iEvt = 0;

    // prepare everything for the next loop
    if ( _pedestalAlgo == EUTELESCOPE::MEANRMS ) {

      // the collection contains several TrackerRawData
      // move back the _pedestal and _noise to the _temp vector
      _tempPede  = _pedestal;
      _tempNoise = _noise;
      
      for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
	_tempEntries.push_back(IntVec( _noise[iDetector].size(), 1));
      }
      
    } else if ( _pedestalAlgo == EUTELESCOPE::AIDAPROFILE ) {
      // in case the AIDAPROFILE algorithm is used, the only thing we
      // need to do is to clean up the previous loop histograms
      // remeber to loop over all detectors
#ifdef MARLIN_USE_AIDA
      for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
	stringstream ss;
	ss << _tempProfile2DName << "-d" << iDetector;
	(dynamic_cast<AIDA::IProfile2D*> (_aidaHistoMap[ss.str()]))->reset();
      }
#endif
    }
    throw RewindDataFilesException(this);
    setReturnValue("IsPedestalFinished", false);
  }
}
