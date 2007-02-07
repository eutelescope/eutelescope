// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelPedestalNoiseProcessor.cc,v 1.2 2007-02-07 16:25:13 bulgheroni Exp $
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

  // reset the temporary arrays
  _tempPede.clear ();
  _tempNoise.clear ();
  _tempEntries.clear ();

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

#ifdef MARLIN_USE_AIDA
  static StringVec loopDirName;
  static StringVec detectorDirName;
  
  for (int iLoop = 0 ; iLoop < _noOfCMIterations + 1; iLoop++) {
    loopDirName.push_back(string("loop-" + iLoop));
    AIDAProcessor::tree(this)->mkdir(loopDirName[iLoop].c_str());
    AIDAProcessor::tree(this)->cd(loopDirName[iLoop]);
    for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
      detectorDirName.push_back(string("detector-" + iDetector) + string("_loop-" + iLoop));
      AIDAProcessor::tree(this)->mkdir(detectorDirName[iDetector].c_str());
    }
    AIDAProcessor::tree(this)->cd("..");
  }
#endif

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

  // increment the loop counter
  ++_iLoop;

  // fill in the histograms
  fillHistos();

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
      for (int yPixel = _minY[iDetector]; yPixel < _maxY[iDetector]; yPixel++) {
	for (int xPixel = _minX[iDetector]; xPixel < _maxX[iDetector]; xPixel++) {
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
    for (int yPixel = _minY[iDetector]; yPixel < _maxY[iDetector]; yPixel++) {
      for (int xPixel = _minX[iDetector]; xPixel < _maxX[iDetector]; xPixel++) {
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
      for (int yPixel = _minY[iDetector]; yPixel < _maxY[iDetector]; yPixel++) {
	for (int xPixel = _minX[iDetector]; xPixel < _maxX[iDetector]; xPixel++) {
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
    }
  }	
  ++_iEvt;
}


std::string EUTelPedestalNoiseProcessor::_pedeDistHistoName   = "PedeDist-";
std::string EUTelPedestalNoiseProcessor::_noiseDistHistoName  = "NoiseDist-";
std::string EUTelPedestalNoiseProcessor::_commonModeHistoName = "CommonMode-";
std::string EUTelPedestalNoiseProcessor::_pedeMapHistoName    = "PedeMap-";
std::string EUTelPedestalNoiseProcessor::_noiseMapHistoName   = "NoiseMap-";
std::string EUTelPedestalNoiseProcessor::_statusMapHistoName  = "StatusMap-";
