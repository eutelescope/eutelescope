// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelPedestalNoiseProcessor.cc,v 1.3 2007-02-08 09:39:32 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h" 
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
#include <AIDA/ITree.h>
#endif

// lcio includes <.h> 
#include <lcio.h>

// system includes <>
#include <string>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;

EUTelPedestalNoiseProcessor gEUTelPedestalNoiseProcessor;       // global instance of the object


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

  // reset all the final arrays
  _pedestal.clear();
  _noise.clear();
  _status.clear();

  // clean up things related to CM 
  _noOfSkippedEvent = 0;

}

void EUTelPedestalNoiseProcessor::processRunHeader (LCRunHeader * rdr) {

  // to make things easier re-cast the input header to the EUTelRunHeaderImpl
  EUTelRunHeaderImpl * runHeader = dynamic_cast<EUTelRunHeaderImpl*>(rdr);

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

  // book histograms
  bookHistos();

}


void EUTelPedestalNoiseProcessor::processEvent (LCEvent * evt) {

  if ( _iLoop == 0 ) firstLoop(evt);
  else otherLoop(evt);

}
  


void EUTelPedestalNoiseProcessor::check (LCEvent * evt) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void EUTelPedestalNoiseProcessor::end() {
  
  // the loop on events is over so we need to move temporary vectors
  // to final vectors
  _pedestal = _tempPede;
  _noise    = _tempNoise;

  // clear the temporary vectors
  _tempPede.clear();
  _tempNoise.clear();
  _tempEntries.clear();
  
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
   
    // move the pedestal, noise and status arrays to a TrackerData
    // ... //
    
    setReturnValue("IsPedestalFinished", true);
  } else {
    // now we need to loop again
    // so reset the event counter
    _iEvt = 0;
    // and the first event flag
    _isFirstEvent = true;

    // throw RewindDataFilesException(this);
    setReturnValue("IsPedestalFinished", false);
  }
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
										static_cast<double>(_status[iDetector][iPixel]));
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
    double threshold  = meanw + (rms * _badPixelMaskCut);



    // scan the noise vector again and apply the cut
    for (unsigned int iPixel = 0; iPixel < _status[iDetector].size(); iPixel++) {
      if ( ( _noise[iDetector][iPixel] > threshold ) && 
	   ( _status[iDetector][iPixel] == EUTELESCOPE::GOODPIXEL ) ) {
	_status[iDetector][iPixel] = EUTELESCOPE::BADPIXEL;
      }
    }
  } // end loop on detector;
}

void EUTelPedestalNoiseProcessor::firstLoop(LCEvent * evt) {

  if ( ( _iEvt < _firstEvent ) || ( _iEvt >= _lastEvent ) ) {
    // pedestal calculation may occuring on a subset of events
    // outside this range, just skip the current event
    if ( ( _iEvt == 0 ) || ( _iEvt ==  _firstEvent - 1) || ( _iEvt == _lastEvent ) || ( _iEvt%10 == 0) )
      cout << "[" << name() << "] Skipping event " << _iEvt << endl;
    ++_iEvt;
    return;
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
      ShortVec::iterator iter = adcValues.begin();
      int nPixel = 0;
      DoubleVec tempDoubleVec;
      while ( iter != adcValues.end() ) {
	tempDoubleVec.push_back( static_cast< double > (*iter));
	++iter;
	++nPixel;
      }
      _tempPede.push_back(tempDoubleVec);
      
      // initialize _tempNoise and _tempEntries with all zeros and
      // ones
      _tempNoise.push_back(DoubleVec(nPixel, 0.));
      _tempEntries.push_back(IntVec(nPixel, 1));
      
      // the status vector can be initialize as well with all
      // GOODPIXEL
      _status.push_back(IntVec(nPixel, EUTELESCOPE::GOODPIXEL));
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
    }     // end loop on detectors
  }

  // increment the event counter
  ++_iEvt;
}

void EUTelPedestalNoiseProcessor::otherLoop(LCEvent * evt) {
  
  if ( ( _iEvt < _firstEvent ) || ( _iEvt >= _lastEvent ) ) {
    // pedestal calculation may occuring on a subset of events
    // outside this range, just skip the current event
    if ( ( _iEvt == 0 ) || ( _iEvt ==  _firstEvent - 1) || ( _iEvt == _lastEvent ) || ( _iEvt%10 == 0) )
      cout << "[" << name() << "] Skipping event " << _iEvt << endl;
    ++_iEvt;
    return;
  }

  // keep the user updated
  if ( _iEvt % 10 == 0 ) {
    cout << "Performing " << _iLoop << " loop on event: " << _iEvt << endl;
  }

  // let me get the rawDataCollection. This is should contain a TrackerRawDataObject
  // for each detector plane in the telescope.
  LCCollectionVec *collectionVec = dynamic_cast < LCCollectionVec * >(evt->getCollection (_rawDataCollectionName));
  

  if ( isFirstEvent() ) {
    // the collection contains several TrackerRawData
    // move back the _pedestal and _noise to the _temp vector
    _tempPede  = _pedestal;
    _tempNoise = _noise;
    
    for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
      _tempEntries.push_back(IntVec( _noise[iDetector].size(), 1));
    }
    
    _isFirstEvent = false;
  }

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
      
      iPixel = 0;
      for (int yPixel = _minY[iDetector]; yPixel <= _maxY[iDetector]; yPixel++) {
	for (int xPixel = _minX[iDetector]; xPixel <= _maxX[iDetector]; xPixel++) {
	  if ( _status[iDetector][iPixel] == EUTELESCOPE::GOODPIXEL ) {
	    double pedeCorrected = adcValues[iPixel] - commonMode;
	    _tempEntries[iDetector][iPixel] = _tempEntries[iDetector][iPixel] + 1;
	    _tempPede[iDetector][iPixel]    = ((_tempEntries[iDetector][iPixel] - 1) * _tempPede[iDetector][iPixel]
					       + pedeCorrected) / _tempEntries[iDetector][iPixel];
	    _tempNoise[iDetector][iPixel]   = sqrt(((_tempEntries[iDetector][iPixel] - 1) * pow(_tempNoise[iDetector][iPixel],2) 
						    + pow(pedeCorrected - _tempPede[iDetector][iPixel], 2)) / 
						   _tempEntries[iDetector][iPixel]);
	  }
	  ++iPixel;
	}
      }
    } else {
      cout << "Skipping event " << _iEvt << " because of max number of rejected pixel exceeded." << endl;
      ++_noOfSkippedEvent;
    }
  }	
  ++_iEvt;
}


std::string EUTelPedestalNoiseProcessor::_pedeDistHistoName   = "PedeDist";
std::string EUTelPedestalNoiseProcessor::_noiseDistHistoName  = "NoiseDist";
std::string EUTelPedestalNoiseProcessor::_commonModeHistoName = "CommonMode";
std::string EUTelPedestalNoiseProcessor::_pedeMapHistoName    = "PedeMap";
std::string EUTelPedestalNoiseProcessor::_noiseMapHistoName   = "NoiseMap";
std::string EUTelPedestalNoiseProcessor::_statusMapHistoName  = "StatusMap";

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
      
      // book an histogram for the pedestal distribution
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

      // book an histogram for the noise distribution
      const int    noiseDistHistoNBin  =   30;
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

      // book a 1d histo for common mode only if loop >= 1
      if (iLoop >= 1) {
	const int    commonModeHistoNBin = 20;
	const double commonModeHistoMin  = -2;
	const double commonModeHistoMax  =  2;
	{
	  stringstream ss;
	  ss << _commonModeHistoName << "-d" << iDetector << "-l" << iLoop;
	  tempHistoName = ss.str();
	}
	AIDA::IHistogram1D * commonModeHisto =
	  AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
								    commonModeHistoNBin, commonModeHistoMin, commonModeHistoMax);
	_aidaHistoMap.insert(make_pair(tempHistoName, commonModeHisto));
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

    }
  } // end on iLoop
#endif // MARLIN_USE_AIDA
}
