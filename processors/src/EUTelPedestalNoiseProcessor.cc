// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id$
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
#include "EUTelEventImpl.h"
#include "EUTelPedestalNoiseProcessor.h"
#include "EUTelHistogramManager.h"
#include "EUTELESCOPE.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
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
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>

// system includes <>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <cstdlib>
#include <limits>
#include <algorithm>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;


// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelPedestalNoiseProcessor::_pedeDistHistoName   = "PedeDist";
std::string EUTelPedestalNoiseProcessor::_noiseDistHistoName  = "NoiseDist";
std::string EUTelPedestalNoiseProcessor::_commonModeHistoName = "CommonMode";
std::string EUTelPedestalNoiseProcessor::_pedeMapHistoName    = "PedeMap";
std::string EUTelPedestalNoiseProcessor::_noiseMapHistoName   = "NoiseMap";
std::string EUTelPedestalNoiseProcessor::_statusMapHistoName  = "StatusMap";
std::string EUTelPedestalNoiseProcessor::_tempProfile2DName   = "TempProfile2D";
std::string EUTelPedestalNoiseProcessor::_fireFreqHistoName   = "FiringFrequency";
std::string EUTelPedestalNoiseProcessor::_aPixelHistoName     = "APixelHisto";
#endif

EUTelPedestalNoiseProcessor::EUTelPedestalNoiseProcessor () :Processor("EUTelPedestalNoiseProcessor") {

  // modify processor description
  _description =
    "EUTelPedestalNoiseProcessor computes the pedestal and noise values of a pixel detector";


  vector< string > rawDataCollectionNameVecExample;
  rawDataCollectionNameVecExample.push_back( "rawdata" );

  // first of all we need to register the input collection
  registerInputCollections (LCIO::TRACKERRAWDATA, "RawDataCollectionNameVec",
                            "Input raw data collection",
                            _rawDataCollectionNameVec, rawDataCollectionNameVecExample );

  // register compulsory parameters
  registerProcessorParameter ("CalculationAlgorithm",
                              "Select the algorithm for pede/noise calculation",
                              _pedestalAlgo,
                              string (EUTELESCOPE::MEANRMS));

  registerProcessorParameter("CommonModeAlgorithm",
                             "Select the algorithm for the common mode calculation. Possible values are:\n"
                             " FullFrame: all pixels in the frame are averaged\n"
                             " RowWise: pixels are averaged line by line",
                             _commonModeAlgo, string( "FullFrame" ));


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

  registerProcessorParameter("MaxNoOfRejectedPixelPerRow",
                             "Maximum allowed number of rejected pixels per row (only with RowWise)",
                             _maxNoOfRejectedPixelPerRow,
                             static_cast < int > (25) );
  registerProcessorParameter("MaxNoOfSkippedRow",
                             "Maximum allowed number of skipped rows (only with RowWise)",
                             _maxNoOfSkippedRow,
                             static_cast< int > ( 15 ) );

  // new names for the pixel masking algorithms
  // since v00-00-09, the user can select multiple bad pixel
  // algorithms. The _badPixelAlgo has been changed into a
  // _badPixelAlgoVec.
  vector< string > badPixelAlgoVecExample;
  badPixelAlgoVecExample.push_back( "NoiseDistribution" );
  badPixelAlgoVecExample.push_back( "DeadPixel" );
  registerProcessorParameter("BadPixelMaskingAlgorithm",
                             "Select the algorithm for bad pixel masking. Possible values are:\n"
                             " NoiseDistribution: removing pixels with noise above PixelMaskUpperNoiseCut in sigma unit\n"
                             " AbsoluteNoiseValue: removing pixels with noise above PixelMaskUpperAbsNoiseCut in ADC value\n"
                             " DeadPixel: removing pixels with noise below PixelMaskLowerAbsNoiseCut in ADC value\n"
                             " AbsolutePedeValue: removing pixels having pedestal too low or high using PixelMaskUpperAbsPedeCut and PixelMaskLowerAbsPedeCut",
                             _badPixelAlgoVec, badPixelAlgoVecExample);

  registerProcessorParameter ("PixelMaskUpperNoiseCut",
                              "Upper threshold for bad pixel identification using NoiseDistribution",
                              _pixelMaskUpperNoiseCut, static_cast < float >(3.5));
  registerProcessorParameter ("PixelMaskUpperAbsNoiseCut",
                              "Upper threshold for bad pixel identification using NoiseDistribution",
                              _pixelMaskUpperAbsNoiseCut, static_cast < float >(3.5));
  registerProcessorParameter ("PixelMaskLowerAbsNoiseCut",
                              "Lower threshold for bad pixel identification using DeadPixel",
                              _pixelMaskLowerAbsNoiseCut, static_cast < float >(0.2));
  registerProcessorParameter ("PixelMaskUpperAbsPedeCut",
                              "Upper threshold for bad pixel identification using AbsolutePedeValue",
                              _pixelMaskUpperAbsPedeCut, static_cast< float > ( 15 ) );
  registerProcessorParameter ("PixelMaskLowerAbsPedeCut",
                              "Lower threshold for bad pixel identification using AbsolutePedeValue",
                              _pixelMaskLowerAbsPedeCut, static_cast< float > ( -15 ) );
  registerProcessorParameter("PixelMaskMaxFiringFrequency",
                             "This is the maximum allowed firing % frequency, being 0.1% the Gaussian limit\n"
                             "Used only during the additional masking loop",
                             _maxFiringFreq, static_cast< float > ( 0.2 ) ) ;
  registerOptionalParameter ("AdditionalMaskingLoop",
                             "Perform an additional loop for bad pixel masking",
                             _additionalMaskingLoop, static_cast<bool> ( true ) );
  registerOptionalParameter ("HitRejectionPreLoop",
                             "Perform a fast first loop to improve the efficiency of hit rejection",
                             _preLoopSwitch, static_cast< bool > ( true ) ) ;


  registerProcessorParameter ("FirstEvent",
                              "First event for pedestal calculation",
                              _firstEvent, static_cast < int > (0));
  registerProcessorParameter ("LastEvent",
                              "Last event for pedestal calculation",
                              _lastEvent, static_cast < int > (-1));
  registerProcessorParameter ("OutputPedeFile","The filename (w/o .slcio) to store the pedestal file",
                              _outputPedeFileName , string("outputpede"));

  registerProcessorParameter ("ASCIIOutputSwitch","Set to true if the pedestal should also be saved as ASCII files",
                              _asciiOutputSwitch, static_cast< bool > ( true ) );

  registerProcessorParameter("HistoInfoFileName", "This is the name of the histogram information file",
                             _histoInfoFileName, string( "histoinfo.xml" ) );

  // now the optional parameters
  registerOptionalParameter ("PedestalCollectionName",
                             "Pedestal collection name",
                             _pedestalCollectionName, string ("pedestalDB"));
  registerOptionalParameter ("NoiseCollectionName",
                             "Noise collection name", _noiseCollectionName,
                             string ("noiseDB"));
  registerOptionalParameter ("StatusCollectionName",
                             "Status collection name",
                             _statusCollectionName, string ("statusDB"));

  _histogramSwitch = true;

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

  // set the geometry ready switch to false
  _isGeometryReady = false;

  // set the loop counter
  if ( _preLoopSwitch ) _iLoop = -1;
  else _iLoop = 0;

  if ( _pedestalAlgo == EUTELESCOPE::MEANRMS ) {
    // reset the temporary arrays
    _tempPede.clear ();
    _tempNoise.clear ();
    _tempEntries.clear ();
  }

#ifndef MARLIN_USE_AIDA
  _histogramSwitch = false;
  if ( _pedestalAlgo == EUTELESCOPE::AIDAPROFILE ) {
    streamlog_out ( ERROR0 )  << "The " << EUTELESCOPE::AIDAPROFILE
                              << " algorithm cannot be applied since Marlin is not using AIDA" << endl
                              << " Algorithm changed to " << EUTELESCOPE::MEANRMS << endl;
    _pedestalAlgo = EUTELESCOPE::MEANRMS;
  }
#endif

  if ( _preLoopSwitch ) {
    _maxValuePos.clear();
    _maxValue.clear();
    _minValuePos.clear();
    _maxValue.clear();
  }


  // reset all the final arrays
  _pedestal.clear();
  _noise.clear();
  _status.clear();

  // set all the bad pixel masking switches
  setBadPixelAlgoSwitches();

  // reset the skip event list
  _skippedEventList.clear();
  _nextEventToSkip = _skippedEventList.begin();

}

void EUTelPedestalNoiseProcessor::processRunHeader (LCRunHeader * rdr) {

  _detectorName = rdr->getDetectorName();

  // to make things easier re-cast the input header to the EUTelRunHeaderImpl
  std::unique_ptr<EUTelRunHeaderImpl> runHeader = std::make_unique<EUTelRunHeaderImpl>(rdr);
  runHeader->addProcessor( type() );

  // increment the run counter
  ++_iRun;

  // make some test on parameters
  if ( ( _pedestalAlgo != EUTELESCOPE::MEANRMS ) &&
       ( _pedestalAlgo != EUTELESCOPE::AIDAPROFILE)
    ) {
    throw InvalidParameterException(string("_pedestalAlgo cannot be " + _pedestalAlgo));
  }

  // the user can decide to limit the pedestal calculation on a
  // sub range of events for many reasons. The most common one is that
  // a run is started with the beam off and some hundreds of events
  // are taken on purpose before switching the beam on. In this way
  // the same file contains both pedestal and data, with pedestal
  // events within a specific event window.
  //
  // From the Marlin steering file the best way the user has to select
  // this range is using the _firstEvent and _lastEvent parameter of
  // the processor it self.
  // There is another variable that can influence this behavior and it
  // is the global MaxRecordNumber. Being global, of course it is
  // dominant with respect to the local _lastEvent setting. Once more,
  // if the EORE is found before the _lastEvent than finalize method
  // is called anyway.
  //


  int maxRecordNumber = Global::parameters->getIntVal("MaxRecordNumber");

  streamlog_out ( DEBUG4 )  << "Event range for pedestal calculation is from " << _firstEvent << " to " << _lastEvent
                            << "\nMaxRecordNumber from the global section is   " << maxRecordNumber << endl;

  // check if the user wants an additional loop to remove the bad pixels
  int additionalLoop = 0;
  if ( _additionalMaskingLoop ) additionalLoop = 1;

  if ( _lastEvent == -1 ) {
    // the user didn't select an upper limit for the event range, so
    // we don't know on how many events the calculation should be done
    //
    // if the global MaxRecordNumber has been set, so warn the user
    // that the procedure could be wrong due to a too low number of
    // records.
    if ( maxRecordNumber != 0 ) {
      streamlog_out ( WARNING2 )  << "The MaxRecordNumber in the Global section of the steering file has been set to "
                                  << maxRecordNumber << ".\n"
                                  << "This means that in order to properly perform the pedestal calculation the maximum allowed number of events is "
                                  << maxRecordNumber / ( _noOfCMIterations + 1 + additionalLoop ) << ".\n"
                                  << "Let's hope it is correct and try to continue." << endl;
    }
  } else {
    // ok we know on how many events the calculation should be done.
    // we can compare this number with the maxRecordNumber if
    // different from 0
    if ( maxRecordNumber != 0 ) {
      if ( (_lastEvent - _firstEvent) * ( _noOfCMIterations + 1 + additionalLoop ) > maxRecordNumber ) {
        streamlog_out ( ERROR4 ) << "The pedestal calculation should be done on " << _lastEvent - _firstEvent
                                 << " times " <<  _noOfCMIterations + 1 << " iterations = "
                                 << (_lastEvent - _firstEvent) * ( _noOfCMIterations + 1 + additionalLoop ) << " records.\n"
                                 << "The global variable MarRecordNumber is limited to " << maxRecordNumber << endl;
        throw InvalidParameterException("MaxRecordNumber");
      }
    }
  }


  if ( _iLoop == 0 ) {
    // write the current header to the output condition file
    LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();

    try {
      lcWriter->open(_outputPedeFileName, LCIO::WRITE_NEW);
    } catch (IOException& e) {
      cerr << e.what() << endl;
      return;
    }

    LCRunHeaderImpl    * lcHeader  = new LCRunHeaderImpl;
    EUTelRunHeaderImpl * newHeader = new EUTelRunHeaderImpl (lcHeader);

    newHeader->lcRunHeader()->setRunNumber(runHeader->lcRunHeader()->getRunNumber());
    newHeader->lcRunHeader()->setDetectorName(runHeader->lcRunHeader()->getDetectorName());
    newHeader->setHeaderVersion(runHeader->getHeaderVersion());
    newHeader->setDataType(runHeader->getDataType());
    newHeader->setDateTime();
    newHeader->setDAQHWName(runHeader->getDAQHWName());
    newHeader->setDAQHWVersion(runHeader->getDAQHWVersion());
    newHeader->setDAQSWName(runHeader->getDAQSWName());
    newHeader->setDAQSWVersion(runHeader->getDAQSWVersion());
    newHeader->setNoOfEvent(runHeader->getNoOfEvent());

    newHeader->addProcessor(name());

    lcWriter->writeRunHeader(lcHeader);
    delete newHeader;
    delete lcHeader;

    lcWriter->close();

  }
}


void EUTelPedestalNoiseProcessor::initializeGeometry( LCEvent * evt ) {

  // starting from this version (v00-00-08-toto) the processor is
  // accepting multiple input collections.
  //
  // The number of sensors will be the sum of the elements of all the
  // input collections, while the detector boundaries are taken from
  // each rawdata cellid.

  _noOfDetector = 0;
  _noOfDetectorVec.clear();

  for ( size_t iCol = 0 ; iCol < _rawDataCollectionNameVec.size(); ++iCol ) {

    try {
      LCCollectionVec * collectionVec = dynamic_cast< LCCollectionVec * > ( evt->getCollection( _rawDataCollectionNameVec.at( iCol ) ) ) ;
      _noOfDetector  += collectionVec->size();
      _noOfDetectorVec.push_back( collectionVec->size() );


      // now to fill in the _minX, _maxX, _minY and _maxY vectors we
      // have to go through all the collection elements.
      CellIDDecoder< TrackerRawDataImpl > rawDataDecoder( collectionVec );
      for ( size_t iDetector = 0; iDetector < collectionVec->size() ; ++iDetector ) {
        TrackerRawDataImpl * rawData = dynamic_cast< TrackerRawDataImpl * > ( collectionVec->getElementAt( iDetector ) );
        _minX.push_back( rawDataDecoder( rawData ) ["xMin"] ) ;
        _maxX.push_back( rawDataDecoder( rawData ) ["xMax"] ) ;
        _minY.push_back( rawDataDecoder( rawData ) ["yMin"] ) ;
        _maxY.push_back( rawDataDecoder( rawData ) ["yMax"] ) ;
        _orderedSensorIDVec.push_back( rawDataDecoder( rawData ) ["sensorID"] );
      }

    } catch (DataNotAvailableException& e) {
      streamlog_out ( WARNING2 ) << "No input collection " << _rawDataCollectionNameVec.at( iCol ) << " is not available in the current event" << endl
                                 << "If this collection is not present in the input file, please remove it from the steering file" << endl
                                 << "Geometry initialization impossible, skipping the event" << endl;

      _isGeometryReady = false;
      throw SkipEventException(this);

    }

  }

  _isGeometryReady = true;

}

void EUTelPedestalNoiseProcessor::processEvent (LCEvent * evt) {

  EUTelEventImpl * eutelEvent = static_cast<EUTelEventImpl*> (evt);
  EventType type              = eutelEvent->getEventType();

  if ( ! _isGeometryReady ) {
    initializeGeometry(evt);
  }

  if ( type == kUNKNOWN ) {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }

  if ( _iLoop == -1 ) preLoop( evt );
  else if ( _iLoop == 0 ) firstLoop(evt);
  else if ( _additionalMaskingLoop ) {
    if ( _iLoop == _noOfCMIterations + 1 ) {
      additionalMaskingLoop(evt);
    } else {
      otherLoop(evt);
    }
  } else {
    otherLoop(evt);
  }

}



void EUTelPedestalNoiseProcessor::check (LCEvent * /* evt */ ) {
  // nothing to check here - could be used to fill check plots in reconstruction processor
}


void EUTelPedestalNoiseProcessor::end() {


  int additionalLoop = 0;
  if ( _additionalMaskingLoop ) additionalLoop = 1;
  if ( _iLoop == _noOfCMIterations + 1 + additionalLoop )  {
    streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
  }  else {
    streamlog_out ( ERROR4 ) << "Not all the iterations have been done because of a too MaxRecordNumber.\n"
                             << "Try to increase it in the global section of the steering file." << endl;
    exit(-1);
  }

}

void EUTelPedestalNoiseProcessor::fillHistos() {

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  streamlog_out ( MESSAGE2 ) << "Filling final histograms " << endl;

  string tempHistoName;
  for (size_t iDetector = 0; iDetector < _noOfDetector; iDetector++) {
    int iPixel = 0;
    for (int yPixel = _minY[iDetector]; yPixel <= _maxY[iDetector]; yPixel++) {
      for (int xPixel = _minX[iDetector]; xPixel <= _maxX[iDetector]; xPixel++) {
        if ( _histogramSwitch ) {
          tempHistoName =  _statusMapHistoName + "_d" + to_string( _orderedSensorIDVec.at( iDetector ) ) + "_l" + to_string( _iLoop );
          if ( AIDA::IHistogram2D * histo = dynamic_cast<AIDA::IHistogram2D*>(_aidaHistoMap[tempHistoName]) )
            histo->fill(static_cast<double>(xPixel), static_cast<double>(yPixel), static_cast<double> (_status[iDetector][iPixel]));
          else {
            streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                      << ".\nDisabling histogramming from now on " << endl;
            _histogramSwitch = false;
          }
        }


        if ( _status[iDetector][iPixel] == EUTELESCOPE::GOODPIXEL) {
          if ( _histogramSwitch ) {
            tempHistoName = _pedeDistHistoName + "_d" + to_string( _orderedSensorIDVec.at( iDetector ) ) + "_l" + to_string( _iLoop );
            if ( AIDA::IHistogram1D * histo = dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]) )
              histo->fill(_pedestal[iDetector][iPixel]);
            else {
              streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                        << ".\nDisabling histogramming from now on " << endl;
              _histogramSwitch = false;
            }
          }

          if ( streamlog::out.write< streamlog::DEBUG0 > () ) {
            if ( (xPixel == 10) && (yPixel == 10 )) {
              streamlog::out()  << "Detector " << _orderedSensorIDVec.at( iDetector ) << " pedestal " << (_pedestal[iDetector][iPixel]) << endl;
            }
          }

          if ( _histogramSwitch ) {
            tempHistoName = _noiseDistHistoName + "_d" + to_string( _orderedSensorIDVec.at( iDetector ) ) + "_l" + to_string( _iLoop );
            if ( AIDA::IHistogram1D * histo = dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]) )
              histo->fill(_noise[iDetector][iPixel]);
            else {
              streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                        << ".\nDisabling histogramming from now on " << endl;
              _histogramSwitch = false;
            }
          }

          if ( _histogramSwitch ) {
            tempHistoName = _pedeMapHistoName + "_d" + to_string( _orderedSensorIDVec.at( iDetector ) ) + "_l" + to_string( _iLoop );
            if ( AIDA::IHistogram2D * histo = dynamic_cast<AIDA::IHistogram2D*>(_aidaHistoMap[tempHistoName]) )
              histo-> fill(static_cast<double>(xPixel), static_cast<double>(yPixel), _pedestal[iDetector][iPixel]);
            else {
              streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                        << ".\nDisabling histogramming from now on " << endl;
              _histogramSwitch = false;
            }
          }

          if ( _histogramSwitch ) {
            tempHistoName = _noiseMapHistoName + "_d" + to_string( _orderedSensorIDVec.at( iDetector ) ) + "_l" + to_string( _iLoop );
            if ( AIDA::IHistogram2D* histo = dynamic_cast<AIDA::IHistogram2D*>(_aidaHistoMap[tempHistoName]) )
              histo->fill(static_cast<double>(xPixel), static_cast<double>(yPixel), _noise[iDetector][iPixel]);
            else {
              streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                        << ".\nDisabling histogramming from now on " << endl;
              _histogramSwitch = false;
            }
          }
        }
        ++iPixel;
      }
    }
  }


#else
  if ( _iEvt == 0 ) streamlog_out ( MESSAGE4 )  << "No histogram produced because Marlin doesn't use AIDA" << endl;
#endif

}

void EUTelPedestalNoiseProcessor::maskBadPixel() {


  // a global counter of bad pixels
  vector<int >  badPixelCounterVec( _noOfDetector, 0 );

  if ( ( !_additionalMaskingLoop ) ||
       ( _iLoop < _noOfCMIterations + 1 )) {

    // prepare vectors of thresholds
    //
    // thresholdNoiseDistVec is the vector of threshold for the
    // NOISEDISTRIBUTION algorithm
    vector<double > thresholdNoiseDistVec;

    // thresholdAbsNoiseVec is the vector of threshold for the
    // ABSOLUTENOISEVALUE algorithm
    vector<double > thresholdAbsNoiseVec( _noOfDetector,  _pixelMaskUpperAbsNoiseCut );

    // thresholdDeadPixelVec is the vector of threshold for the
    // DEADPIXEL algorithm
    vector<double > thresholdDeadPixelVec( _noOfDetector, _pixelMaskLowerAbsNoiseCut );

    // thresholdLowerAbsPedeVec is the vector of threshold for the
    // ABSOLUTEPEDEVALUE algorithm
    vector<double > thresholdLowerAbsPedeVec( _noOfDetector, _pixelMaskLowerAbsPedeCut );
    vector<double > thresholdUpperAbsPedeVec( _noOfDetector, _pixelMaskUpperAbsPedeCut );

    if ( _badPixelNoiseDistributionSwitch ) {
      // to do it we need to know the mean value and the RMS of the noise
      // vector.

      for (size_t iDetector = 0; iDetector < _noOfDetector; iDetector++) {

        double sumw  = 0;
        double sumw2 = 0;
        double num   = 0;

        // begin a first loop on all pixel to calculate the masking threshold
        for ( size_t iPixel = 0; iPixel < _status[iDetector].size(); iPixel++) {
          if ( _status[iDetector][iPixel] == EUTELESCOPE::GOODPIXEL ) {
            sumw  += _noise[iDetector][iPixel];
            sumw2 += pow(_noise[iDetector][iPixel],2);
            ++num;
          }
        }
        double meanw      = sumw  / num;
        double meanw2     = sumw2 / num;
        double rms        = sqrt( meanw2 - pow(meanw,2));
        thresholdNoiseDistVec.push_back( meanw + (rms * _pixelMaskUpperNoiseCut) );
        streamlog_out ( DEBUG4 ) << "Mean noise value is " << meanw << " ADC\n"
          "RMS of noise is " << rms << " ADC\n"
          "Masking threshold is set to " << thresholdNoiseDistVec[iDetector] << endl;

        // now reloop on pixels to mask the bad out
        for ( size_t iPixel = 0; iPixel < _status[iDetector].size(); ++iPixel ) {
          if (  ( _noise[iDetector][iPixel]  > thresholdNoiseDistVec[iDetector] ) &&
                ( _status[iDetector][iPixel] == EUTELESCOPE::GOODPIXEL ) ) {
            // first of all make the pixel bad
            _status[iDetector][iPixel] = EUTELESCOPE::BADPIXEL;
            streamlog_out ( DEBUG0 ) <<  "Masking pixel number " << iPixel
                                     << " on detector " << _orderedSensorIDVec.at( iDetector )
                                     << " (" << _noise[iDetector][iPixel]
                                     << " > " << thresholdNoiseDistVec[iDetector] << ")" << endl;
            badPixelCounterVec[iDetector]++;
          }
        }
      }
    }


    if ( _badPixelAbsNoiseSwitch ) {

      for ( size_t iDetector = 0 ; iDetector < _noOfDetector; iDetector++ )  {
        for ( size_t iPixel = 0; iPixel < _status[iDetector].size(); ++iPixel ) {
          if ( ( _noise[iDetector][iPixel] > thresholdAbsNoiseVec[iDetector] ) &&
               ( _status[iDetector][iPixel] == EUTELESCOPE::GOODPIXEL ) ) {
            _status[iDetector][iPixel] = EUTELESCOPE::BADPIXEL;
            streamlog_out ( DEBUG0 ) <<  "Masking pixel number " << iPixel
                                     << " on detector " << _orderedSensorIDVec.at( iDetector )
                                     << " (" << _noise[iDetector][iPixel]
                                     << " > " << thresholdAbsNoiseVec[iDetector] << ")" << endl;
            badPixelCounterVec[iDetector]++;
          }
        }
      }

    }


    if ( _badPixelDeadPixelSwitch ) {

      for ( size_t iDetector = 0 ; iDetector < _noOfDetector; ++iDetector ) {
        for ( size_t iPixel = 0; iPixel < _status[iDetector].size(); ++iPixel ) {
          if ( ( _noise[iDetector][iPixel]  < thresholdDeadPixelVec[iDetector] ) &&
               ( _status[iDetector][iPixel] == EUTELESCOPE::GOODPIXEL ) ) {
            _status[iDetector][iPixel] = EUTELESCOPE::BADPIXEL;
            streamlog_out ( DEBUG0 ) <<  "Masking pixel number " << iPixel
                                     << " on detector " << _orderedSensorIDVec.at( iDetector )
                                     << " (" << _noise[iDetector][iPixel]
                                     << " > " << thresholdDeadPixelVec[iDetector] << ")" << endl;
            badPixelCounterVec[iDetector]++;
          }
        }
      }

    }


    if ( _badPixelAbsPedeSwitch ) {

      for ( size_t iDetector = 0; iDetector < _noOfDetector ; ++iDetector ) {
        for ( size_t iPixel = 0; iPixel < _status[iDetector].size(); ++iPixel ) {
          if ( ( ( _pedestal[iDetector][iPixel] > thresholdUpperAbsPedeVec[iDetector] ) ||
                 ( _pedestal[iDetector][iPixel] < thresholdLowerAbsPedeVec[iDetector] ) ) &&
               _status[iDetector][iPixel] == EUTELESCOPE::GOODPIXEL ) {
            _status[iDetector][iPixel] = EUTELESCOPE::BADPIXEL;
            streamlog_out( DEBUG0 ) << "Masking pixel number " << iPixel
                                    << " on detector " << _orderedSensorIDVec.at( iDetector ) << endl;
            badPixelCounterVec[iDetector]++;
          }
        }


      }
    }

    streamlog_out ( MESSAGE4 )  << "Masking summary after loop " << _iLoop << ": " << endl;
    int totalBad = 0;
    int total    = 0;
    for ( size_t iDetector = 0; iDetector < _noOfDetector; ++iDetector ) {
      streamlog_out ( MESSAGE4 ) << "Detector ID " << setw(4) << _orderedSensorIDVec[ iDetector ] << " has "
                                 << setw(8) << badPixelCounterVec[iDetector] << " bad pixels (" << 100 * (1.0 * badPixelCounterVec[iDetector] )/ _status[iDetector].size()
                                 << "%). " << endl;
      totalBad += badPixelCounterVec[iDetector];
      total    += _status[iDetector].size();
    }
    streamlog_out ( MESSAGE4 ) << "Total masked pixels = " << totalBad << " (" << 100 * (1.0 * totalBad) / total << "%). " << endl;
  }

  if ( (  _additionalMaskingLoop ) &&
       ( _iLoop == _noOfCMIterations + 1 )) {
    // now masking relying on the additional loop
    for ( size_t iDetector = 0 ; iDetector < _noOfDetector; iDetector++ ) {

      for (unsigned int iPixel = 0; iPixel < _status[iDetector].size(); iPixel++) {
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
        if ( _histogramSwitch ) {
          string tempHistoName;
          tempHistoName = _fireFreqHistoName + "_d" + to_string( _orderedSensorIDVec.at( iDetector ) ) + "_l" + to_string( _iLoop );
          if ( AIDA::IHistogram1D * histo = dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[ tempHistoName ] ))
            histo->fill( (static_cast<double> ( _hitCounter[ iDetector ][ iPixel ] )) / _iEvt * 100. );
        }
#endif
        if ( static_cast< double > ( _hitCounter[ iDetector ][ iPixel ] ) / _iEvt * 100. > _maxFiringFreq  ) {
          _status[ iDetector ][ iPixel ] = EUTELESCOPE::BADPIXEL;
          badPixelCounterVec[iDetector]++;
        }
      }

    } // end loop on detector
    streamlog_out ( MESSAGE4 )  << "Masking summary after loop " << _iLoop << ": " << endl;
    int totalBad = 0;
    int total    = 0;
    for ( size_t iDetector = 0; iDetector < _noOfDetector; ++iDetector ) {
      streamlog_out ( MESSAGE4 ) << "Detector ID " << setw(4) << _orderedSensorIDVec[ iDetector ] << " has "
                                 << setw(8) << badPixelCounterVec[iDetector] << " bad pixels (" << 100 * (1.0 * badPixelCounterVec[iDetector] )/ _status[iDetector].size()
                                 << "%). " << endl;
      totalBad += badPixelCounterVec[iDetector];
      total    += _status[iDetector].size();
    }
    streamlog_out ( MESSAGE4 ) << "Total masked pixels = " << totalBad << " (" << 100 * (1.0 * totalBad) / total << "%). " << endl;
  }
}


void EUTelPedestalNoiseProcessor::preLoop( LCEvent * event ) {

  // first re-cast to event to EUTelEventImpl for better access
  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);

  // In case this is the last event, or it is the kEORE, then rewind
  // the data and return immediately
  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG4 ) << "EORE found: calling simpleRewind()." << endl;
    _iLoop = 0;
    simpleRewind();
    return;
  }

  if ( ( _lastEvent != -1 ) && ( _iEvt >= _lastEvent ) ) {
    streamlog_out ( DEBUG4 ) << "Looping limited by _lastEvent: calling simpleRewind()." << endl;
    _preLoopSwitch = false;
    _iLoop = 0;
    simpleRewind();
    return;
  }

  // in case this event is before the first to be considered, just
  // skip it
  if ( _iEvt < _firstEvent ) {
    ++_iEvt;
    throw SkipEventException(this);
  }

  // make some initialization (only in the first event
  if ( isFirstEvent() ) {

    for ( size_t iCol = 0;  iCol < _rawDataCollectionNameVec.size(); ++iCol ) {

      try {
        LCCollectionVec *collectionVec = dynamic_cast < LCCollectionVec * >(evt->getCollection( _rawDataCollectionNameVec.at( iCol ) ));

        for ( size_t iDetector = 0; iDetector < collectionVec->size(); iDetector++) {

          TrackerRawData * trackerRawData = dynamic_cast< TrackerRawData * > ( collectionVec->getElementAt(  iDetector ) );

          ShortVec adcValues = trackerRawData->getADCValues ();

          // we have to initialize all the vectors only if this is the
          // first collection
          ShortVec maxValue   ( adcValues.size(), numeric_limits< short >::min() );
          ShortVec maxValuePos( adcValues.size(),                             -1 );
          ShortVec minValue   ( adcValues.size(), numeric_limits< short >::max() );
          ShortVec minValuePos( adcValues.size(),                             -1 );


          _maxValue.push_back   ( maxValue    );
          _maxValuePos.push_back( maxValuePos );
          _minValue.push_back   ( minValue    );
          _minValuePos.push_back( minValuePos );
        }

      } catch (DataNotAvailableException& e) {
        streamlog_out ( WARNING2 ) << "No input collection " << _rawDataCollectionNameVec.at( iCol ) << " is not available in the current event" << endl;
      }
    }

    _isFirstEvent = false;
  }


  // here is the real begin
  for ( size_t iCol = 0 ; iCol < _rawDataCollectionNameVec.size(); ++iCol ) {

    try {

      LCCollectionVec *collectionVec = dynamic_cast < LCCollectionVec * >(evt->getCollection( _rawDataCollectionNameVec.at( iCol ) ));

      for ( size_t iDetector = 0 ; iDetector < collectionVec->size() ; iDetector++) {

        size_t detectorOffset = ( iCol == 0 ) ? 0 : _noOfDetectorVec.at( iCol - 1 );

        TrackerRawData *trackerRawData = dynamic_cast < TrackerRawData * >(collectionVec->getElementAt( iDetector ) );
        ShortVec adcValues = trackerRawData->getADCValues ();

        for ( size_t iPixel = 0; iPixel < adcValues.size(); ++iPixel ) {
          short currentVal = adcValues[ iPixel ];

          if ( currentVal > _maxValue[ iDetector + detectorOffset] [ iPixel ] ) {
            _maxValue   [ iDetector + detectorOffset ] [ iPixel ] = currentVal;
            _maxValuePos[ iDetector + detectorOffset ] [ iPixel ] = _iEvt;
          }

          if ( currentVal < _minValue[ iDetector + detectorOffset ] [ iPixel ] ) {
            _minValue   [ iDetector + detectorOffset ] [ iPixel ] = currentVal;
            _minValuePos[ iDetector + detectorOffset ] [ iPixel ] = _iEvt;
          }

        }

      }

    } catch (DataNotAvailableException& e) {
      streamlog_out ( WARNING2 ) << "No input collection " << _rawDataCollectionNameVec.at( iCol ) << " is not available in the current event ("
                                 << event->getEventNumber() << ")" << endl;
    }

  }

  ++_iEvt;
}

void EUTelPedestalNoiseProcessor::firstLoop(LCEvent * event) {

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);

  // do some checks in order to see if we have to continue or to stop
  // with the processor.
  //
  // 1. we have to go immediately to the finalize if this is a EORE event
  // 2. we have to go to the finalize if the user select an event
  // range for pedestal calculation and the current event number is
  // out of range
  // 3. we have to skip this event if _iEvt is < than the first event
  // selected for pedestal calculation


  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG4 ) << "EORE found: calling finalizeProcessor()." << endl;
    finalizeProcessor( false );
  }

  if ( ( _lastEvent != -1 ) && ( _iEvt >= _lastEvent ) ) {
    streamlog_out ( DEBUG4 ) << "Looping limited by _lastEvent: calling finalizeProcessor()." << endl;
    finalizeProcessor( false );
  }

  if ( _iEvt < _firstEvent ) {
    ++_iEvt;
    throw SkipEventException(this);
  }

  if ( isFirstEvent() ) {

    for ( size_t iCol = 0; iCol < _rawDataCollectionNameVec.size() ; ++iCol ) {

      try {

        LCCollectionVec * collectionVec = dynamic_cast < LCCollectionVec * >(evt->getCollection (_rawDataCollectionNameVec.at( iCol ) ));

        for ( size_t iDetector = 0 ; iDetector < collectionVec->size() ; ++iDetector ) {

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

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
            // in the case of AIDAPROFILE we don't need any vectors since
            // everything is done by the IProfile2D automatically
            int iPixel = 0;
            size_t detectorOffset = ( iCol == 0 ) ? 0 : _noOfDetectorVec.at( iCol - 1 );
            stringstream ss;
            ss << _tempProfile2DName << "_d" << _orderedSensorIDVec.at( iDetector + detectorOffset );
            for (int yPixel = _minY[ iDetector + detectorOffset ]; yPixel <= _maxY[ iDetector + detectorOffset ]; yPixel++) {
              for (int xPixel = _minX[ iDetector + detectorOffset ]; xPixel <= _maxX[ iDetector + detectorOffset ]; xPixel++) {
                double temp = static_cast<double> (adcValues[iPixel]);
                if ( AIDA::IProfile2D * profile = dynamic_cast<AIDA::IProfile2D*> (_aidaHistoMap[ss.str()]) ) {
                  profile ->fill(static_cast<double> (xPixel), static_cast<double> (yPixel), temp);
                } else {
                  streamlog_out ( ERROR4 )  << "Irreversible error: " << ss.str() << " is not available. Sorry for quitting." << endl;
                  exit(-1);
                }
                ++iPixel;
              }
            }
#endif
          }

          // the status vector can be initialize as well with all
          // GOODPIXEL
          _status.push_back(ShortVec(adcValues.size(), EUTELESCOPE::GOODPIXEL));

          // if the user wants to add an additional loop on events to
          // mask pixels singing too loud, so the corresponding counter
          // vector should be reset
          if ( _additionalMaskingLoop ) _hitCounter.push_back( ShortVec( adcValues.size(), 0) );

        } // end of detector loop

      } catch (DataNotAvailableException& e) {
        streamlog_out ( WARNING2 ) << "No input collection " << _rawDataCollectionNameVec.at( iCol ) << " is not available in the current event" << endl;
      }

    }

    bookHistos();

    _isFirstEvent = false;

  } else {

    // this is when it is not the first event
    for ( size_t iCol = 0 ; iCol < _rawDataCollectionNameVec.size() ; ++iCol ) {

      try {
        LCCollectionVec *collectionVec = dynamic_cast < LCCollectionVec * >(evt->getCollection (_rawDataCollectionNameVec.at( iCol ) ));

        // after the firstEvent all temp vectors and the status one have
        // the correct number of entries for both indexes
        // loop on the detectors
        for ( size_t iDetector = 0; iDetector < collectionVec->size() ; iDetector++) {

          // get the TrackerRawData object from the collection for this plane
          TrackerRawData *trackerRawData = dynamic_cast < TrackerRawData * >(collectionVec->getElementAt (iDetector));
          ShortVec adcValues = trackerRawData->getADCValues ();

          size_t detectorOffset = ( iCol == 0 ) ? 0 : _noOfDetectorVec.at( iCol - 1 );

          if ( _pedestalAlgo == EUTELESCOPE::MEANRMS ) {

            // start looping on all pixels
            int iPixel = 0;
            for (int yPixel = _minY[ iDetector + detectorOffset ]; yPixel <= _maxY[ iDetector + detectorOffset ]; yPixel++) {
              for (int xPixel = _minX[ iDetector + detectorOffset ]; xPixel <= _maxX[ iDetector + detectorOffset ]; xPixel++) {
                short currentVal = adcValues[iPixel];
                bool use = true;
                if ( _preLoopSwitch && ( ( _iEvt == _maxValuePos[ iDetector + detectorOffset ] [ iPixel ] ) ||
                                         ( _iEvt == _minValuePos[ iDetector + detectorOffset ] [ iPixel ] ) )  ) {
                  use = false;
                }

                if ( use ) {

                  _tempEntries[iDetector + detectorOffset][iPixel]    =   _tempEntries[iDetector + detectorOffset][iPixel] + 1;
                  _tempPede   [iDetector + detectorOffset][iPixel]    = ((_tempEntries[iDetector + detectorOffset][iPixel] - 1)
                                                                         * _tempPede[iDetector + detectorOffset][iPixel]
                                                                         +  currentVal) / _tempEntries[iDetector+detectorOffset][iPixel];
                  _tempNoise  [iDetector + detectorOffset][iPixel]   = sqrt(((_tempEntries[iDetector+detectorOffset][iPixel] - 1)
                                                                             * pow(_tempNoise[iDetector+detectorOffset][iPixel],2)
                                                                             + pow( currentVal - _tempPede[iDetector+detectorOffset][iPixel], 2)) /
                                                                            _tempEntries[iDetector+detectorOffset][iPixel]);
                }

                ++iPixel;
              } // end loop on xPixel
            }   // end loop on yPixel


          } else if ( _pedestalAlgo == EUTELESCOPE::AIDAPROFILE ) {

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
            stringstream ss;
            ss << _tempProfile2DName << "_d" << _orderedSensorIDVec.at( iDetector + detectorOffset );

            int iPixel = 0;
            for (int yPixel = _minY[iDetector+detectorOffset]; yPixel <= _maxY[iDetector+detectorOffset]; yPixel++) {
              for (int xPixel = _minX[iDetector+detectorOffset]; xPixel <= _maxX[iDetector+detectorOffset]; xPixel++) {
                bool use = true;
                if ( _preLoopSwitch && ( ( _iEvt == _maxValuePos[ iDetector + detectorOffset ] [ iPixel ] ) ||
                                         ( _iEvt == _minValuePos[ iDetector + detectorOffset ] [ iPixel ] ) )  ) {
                  use = false;
                }
                if ( use ) {
                  if ( AIDA::IProfile2D* profile = dynamic_cast<AIDA::IProfile2D*> (_aidaHistoMap[ss.str()]) )
                    profile->fill(static_cast<double> (xPixel), static_cast<double> (yPixel), static_cast<double> (adcValues[iPixel]));
                  else {
                    streamlog_out ( ERROR5 ) << "Irreversible error: " << ss.str() << " is not available. Sorry for quitting." << endl;
                    exit(-1);
                  }
                }
                ++iPixel;
              }
            }
#endif
          }

        }     // end loop on detectors

      } catch (DataNotAvailableException& e) {
        streamlog_out ( WARNING2 ) << "No input collection " << _rawDataCollectionNameVec.at( iCol ) << " is not available in the current event" << endl;
      }

    }
    // increment the event number
    ++_iEvt;
  } // end elif firstEvent

}

void EUTelPedestalNoiseProcessor::otherLoop(LCEvent * event) {

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);

  // do some checks in order to see if we have to continue or to stop
  // with the processor.
  //
  // 1. we have to go immediately to the finalize if this is a EORE event
  // 2. we have to go to the finalize if the user select an event
  // range for pedestal calculation and the current event number is
  // out of range
  // 3. we have to skip this event if _iEvt is < than the first event
  // selected for pedestal calculation


  if ( evt->getEventType() == kEORE ) finalizeProcessor( false );
  if ( ( _lastEvent != -1 ) && ( _iEvt >= _lastEvent ) ) finalizeProcessor( false );
  if ( _iEvt < _firstEvent ) {
    ++_iEvt;
    throw SkipEventException(this);
  }


  for ( size_t iCol = 0 ; iCol < _rawDataCollectionNameVec.size() ; ++iCol ) {

    // let me get the rawDataCollection. This is should contain a TrackerRawDataObject
    // for each detector plane in the telescope.
    try {
      LCCollectionVec *collectionVec = dynamic_cast < LCCollectionVec * >(evt->getCollection (_rawDataCollectionNameVec.at( iCol )));

      for ( size_t iDetector = 0; iDetector < collectionVec->size(); iDetector++) {

        // get the TrackerRawData object from the collection for this detector
        TrackerRawData *trackerRawData = dynamic_cast < TrackerRawData * >(collectionVec->getElementAt (iDetector));
        ShortVec adcValues = trackerRawData->getADCValues ();

        // new approach for a better common mode calculation. The idea
        // is that instead of using, as before, a single value of
        // common mode per matrix, we will have a vector of floats
        // containing the common mode correction for each pixel
        vector< float > commonModeCorVec;
        commonModeCorVec.clear();

        bool isEventValid = true;
        int    skippedPixel = 0;
        int    skippedRow   = 0;

        size_t detectorOffset = ( iCol == 0 ) ? 0 : _noOfDetectorVec.at( iCol - 1 );

        if ( _commonModeAlgo == EUTELESCOPE::FULLFRAME ) {

          double pixelSum     = 0.;
          double commonMode   = 0.;
          int    goodPixel    = 0;
          int    iPixel       = 0;

          // start looping on all pixels for hit rejection
          for (int yPixel = _minY[iDetector + detectorOffset ]; yPixel <= _maxY[iDetector + detectorOffset ]; yPixel++) {
            for (int xPixel = _minX[iDetector + detectorOffset ]; xPixel <= _maxX[iDetector + detectorOffset ]; xPixel++) {
              bool isHit  = ( ( adcValues[iPixel] - _pedestal[iDetector + detectorOffset][iPixel] ) > _hitRejectionCut * _noise[iDetector + detectorOffset][iPixel] );
              bool isGood = ( _status[iDetector + detectorOffset ][iPixel] == EUTELESCOPE::GOODPIXEL );
              if ( !isHit && isGood ) {
                pixelSum += adcValues[iPixel] - _pedestal[iDetector + detectorOffset ][iPixel];
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
            commonModeCorVec.insert( commonModeCorVec.begin(), iPixel + 1 , commonMode );
            isEventValid = true;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
              string histoname = _commonModeHistoName + "_d" + to_string( _orderedSensorIDVec.at( iDetector + detectorOffset ) )
                + "_l" + to_string( _iLoop );
              AIDA::IHistogram1D * histo = (dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[ histoname ]));
              if ( histo ) {
                histo->fill(commonMode);
              }
#endif

          } else {

            isEventValid = false;

          }

        } else if ( _commonModeAlgo == EUTELESCOPE::ROWWISE ) {

          int    iPixel       = 0;
          int    colCounter   = 0;
          int    rowLength    = _maxX[iDetector + detectorOffset] -  _minX[iDetector + detectorOffset] + 1;

          for (int yPixel = _minY[iDetector + detectorOffset]; yPixel <= _maxY[iDetector + detectorOffset]; yPixel++) {

            double pixelSum           = 0.;
            double commonMode         = 0.;
            int    goodPixel          = 0;
            int    skippedPixelPerRow = 0;

            for ( int xPixel = _minX[iDetector + detectorOffset]; xPixel <= _maxX[iDetector + detectorOffset]; xPixel++) {
              bool isHit  = ( ( adcValues[iPixel] - _pedestal[iDetector + detectorOffset][iPixel] ) > _hitRejectionCut * _noise[iDetector + detectorOffset][iPixel] );
              bool isGood = ( _status[iDetector + detectorOffset][iPixel] == EUTELESCOPE::GOODPIXEL );
              if ( !isHit && isGood ) {
                pixelSum += adcValues[iPixel] - _pedestal[iDetector + detectorOffset][iPixel];
                ++goodPixel;
              } else if ( isHit ) {
                ++skippedPixelPerRow;
                ++skippedPixel;
              }
              ++iPixel;
            }

            // we are now at the end of the row, so let's calculate the
            // common mode
            if ( ( skippedPixelPerRow < _maxNoOfRejectedPixelPerRow ) &&
                 ( goodPixel != 0 ) ) {
              commonMode = pixelSum / goodPixel ;
              commonModeCorVec.insert( commonModeCorVec.begin() + colCounter * rowLength, rowLength, commonMode );

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
              string histoname = _commonModeHistoName + "_d" + to_string( _orderedSensorIDVec.at( iDetector + detectorOffset ) )
                + "_l" + to_string( _iLoop );
              AIDA::IHistogram1D * histo = (dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[ histoname ]));
              if ( histo ) {
                histo->fill(commonMode);
              }
#endif

            } else {
              commonModeCorVec.insert( commonModeCorVec.begin() + colCounter * rowLength, rowLength, 0. );
              ++skippedRow;
            }

            ++colCounter;
          }


          if ( skippedRow < _maxNoOfSkippedRow ) {

            isEventValid = true;

          } else {

            isEventValid = false;

          }

        } else {
          streamlog_out ( ERROR4 ) << "Unknown common mode algorithm. Using flat null correction" << endl;
          commonModeCorVec.insert( commonModeCorVec.begin(),
                                   ( _maxY[iDetector + detectorOffset] - _minY[iDetector + detectorOffset] + 1 ) *
                                   ( _maxX[iDetector + detectorOffset] - _minX[iDetector + detectorOffset] + 1 ),
                                   0. );
          isEventValid = true;
        }


        if ( isEventValid ) {

          int iPixel = 0;
          for (int yPixel = _minY[iDetector + detectorOffset]; yPixel <= _maxY[iDetector + detectorOffset]; yPixel++) {
            for (int xPixel = _minX[iDetector + detectorOffset]; xPixel <= _maxX[iDetector + detectorOffset]; xPixel++) {
              if ( _status[iDetector + detectorOffset][iPixel] == EUTELESCOPE::GOODPIXEL ) {
                double pedeCorrected = adcValues[iPixel] - commonModeCorVec[iPixel];
                if ( std::abs( pedeCorrected - _pedestal[iDetector + detectorOffset][iPixel] ) < _hitRejectionCut * _noise[iDetector + detectorOffset][iPixel] ) {
                  if ( _pedestalAlgo == EUTELESCOPE::MEANRMS ) {

                    bool use = true;
                    if ( _preLoopSwitch && ( ( _iEvt == _maxValuePos[ iDetector  + detectorOffset ] [ iPixel ] ) ||
                                             ( _iEvt == _minValuePos[ iDetector  + detectorOffset ] [ iPixel ] ) )  ) {
                      use = false;
                    }
                    if ( use ) {
                      _tempEntries[iDetector + detectorOffset][iPixel] = _tempEntries[iDetector + detectorOffset][iPixel] + 1;
                      _tempPede[iDetector + detectorOffset][iPixel]    = ((_tempEntries[iDetector + detectorOffset][iPixel] - 1)
                                                                          * _tempPede[iDetector + detectorOffset][iPixel]
                                                                          + pedeCorrected) / _tempEntries[iDetector + detectorOffset][iPixel];
                      _tempNoise[iDetector + detectorOffset][iPixel]   = sqrt(((_tempEntries[iDetector + detectorOffset][iPixel] - 1)
                                                                               * pow(_tempNoise[iDetector + detectorOffset][iPixel],2)
                                                                               + pow(pedeCorrected - _tempPede[iDetector + detectorOffset][iPixel], 2)) /
                                                                              _tempEntries[iDetector + detectorOffset][iPixel]);
                    }
                  } else if ( _pedestalAlgo == EUTELESCOPE::AIDAPROFILE) {
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
                    bool use = true;
                    if ( _preLoopSwitch && ( ( _iEvt == _maxValuePos[ iDetector  + detectorOffset ] [ iPixel ] ) ||
                                             ( _iEvt == _minValuePos[ iDetector  + detectorOffset ] [ iPixel ] ) )  ) {
                      use = false;
                    }
                    if ( use ) {
                      stringstream ss;
                      ss << _tempProfile2DName << "_d" << _orderedSensorIDVec.at( iDetector  + detectorOffset );
                      (dynamic_cast<AIDA::IProfile2D*> (_aidaHistoMap[ss.str()]))
                        ->fill(static_cast<double> (xPixel), static_cast<double> (yPixel), pedeCorrected);
                    }
#endif
                  }
                }
              }
              ++iPixel;
            }
          }
        } else {
          if ( _commonModeAlgo == EUTELESCOPE::FULLFRAME ) {
            streamlog_out ( WARNING2 ) <<  "Skipping event " << _iEvt << " because of max number of rejected pixels exceeded. ("
                                       << skippedPixel << ") on detector " << _orderedSensorIDVec.at( iDetector ) << endl;
          } else if ( _commonModeAlgo == EUTELESCOPE::ROWWISE ) {
            streamlog_out ( WARNING2 ) <<  "Skipping event " << _iEvt << " because of max number of skipped rows is reached. ("
                                       << skippedRow << ") on detector " << _orderedSensorIDVec.at( iDetector ) << endl;
          }

          // the event has been skipped, so add this event number to the
          // skipped list
          _skippedEventList.push_back( _iEvt );
        }
      }
    } catch (DataNotAvailableException& e) {
      streamlog_out ( WARNING2 ) << "No input collection " << _rawDataCollectionNameVec.at( iCol ) << " is not available in the current event" << endl;
    }
  }
  ++_iEvt;

}

void EUTelPedestalNoiseProcessor::bookHistos() {

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  // histograms are grouped in loops and detectors
  streamlog_out ( MESSAGE2 ) << "Booking histograms " << endl;

  std::unique_ptr<EUTelHistogramManager> histoMgr = std::make_unique<EUTelHistogramManager>(_histoInfoFileName);
  EUTelHistogramInfo    * histoInfo;
  bool                    isHistoManagerAvailable;

  try {
    isHistoManagerAvailable = histoMgr->init();
  } catch ( ios::failure& e) {
    streamlog_out ( ERROR1 )  << "I/O problem with " << _histoInfoFileName << endl
                              << "Continuing without histogram manager"   << endl;
    isHistoManagerAvailable = false;
  } catch ( ParseException& e ) {
    streamlog_out ( ERROR1 )  << e.what() << endl  << "Continuing without histogram manager" << endl;
    isHistoManagerAvailable = false;
  }


  string tempHistoName;

  // start looping on the number of loops. Remember that we have one
  // loop more than the number of common mode iterations
  for (int iLoop = 0; iLoop < _noOfCMIterations + 1; iLoop++) {

    // prepare the name of the current loop directory and add it the
    // the current ITree
    string loopDirName = "loop_" + to_string( iLoop );
    AIDAProcessor::tree(this)->mkdir(loopDirName.c_str());

    // start looping on detectors
    for ( size_t iDetector = 0; iDetector < _noOfDetector;  iDetector++) {

      // prepare the name of the current detector and add it to the
      // current ITree inside the current loop folder
      string detectorDirName = "detector_" + to_string( _orderedSensorIDVec.at( iDetector ) );
      string basePath = loopDirName + "/" + detectorDirName + "/";
      AIDAProcessor::tree(this)->mkdir(basePath.c_str());

      // book an histogram for the pedestal distribution
      int    pedeDistHistoNBin   = 100;
      double pedeDistHistoMin    = -20.;
      double pedeDistHistoMax    =  29.;
      tempHistoName = _pedeDistHistoName + "_d" + to_string( _orderedSensorIDVec.at( iDetector ) ) + "_l" + to_string( iLoop );
      if ( isHistoManagerAvailable ) {
        histoInfo = histoMgr->getHistogramInfo( _pedeDistHistoName );
        if ( histoInfo ) {
          streamlog_out ( DEBUG1 ) << (*histoInfo) << endl;
          pedeDistHistoNBin  = histoInfo->_xBin;
          pedeDistHistoMin   = histoInfo->_xMin;
          pedeDistHistoMax   = histoInfo->_xMax;
        }
      }

      AIDA::IHistogram1D * pedeDistHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                  pedeDistHistoNBin, pedeDistHistoMin, pedeDistHistoMax);
      if ( pedeDistHisto ) {
        _aidaHistoMap.insert(make_pair(tempHistoName, pedeDistHisto));
        pedeDistHisto->setTitle("Pedestal distribution");
      } else {
        streamlog_out ( ERROR1 ) << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                 << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }

      // book an histogram for the noise distribution
      int    noiseDistHistoNBin  =  100;
      double noiseDistHistoMin   =  -5.;
      double noiseDistHistoMax   =  15.;
      tempHistoName = _noiseDistHistoName + "_d" + to_string( _orderedSensorIDVec.at( iDetector ) ) + "_l" + to_string( iLoop );

      if ( isHistoManagerAvailable ) {
        histoInfo = histoMgr->getHistogramInfo(_noiseDistHistoName);
        if ( histoInfo ) {
          streamlog_out ( DEBUG1 )  << (* histoInfo ) << endl;
          noiseDistHistoNBin = histoInfo->_xBin;
          noiseDistHistoMin  = histoInfo->_xMin;
          noiseDistHistoMax  = histoInfo->_xMax;
        }
      }
      AIDA::IHistogram1D * noiseDistHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                  noiseDistHistoNBin, noiseDistHistoMin, noiseDistHistoMax);
      if ( noiseDistHisto ) {
        _aidaHistoMap.insert(make_pair(tempHistoName, noiseDistHisto));
        noiseDistHisto->setTitle("Noise distribution");
      } else {
        streamlog_out ( ERROR1 ) << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                 << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }

      // book a 1d histo for common mode only if loop >= 1
      if (iLoop >= 1) {
        int    commonModeHistoNBin = 100;
        double commonModeHistoMin  =  -2;
        double commonModeHistoMax  =   2;
        tempHistoName = _commonModeHistoName + "_d" + to_string( _orderedSensorIDVec.at( iDetector ) ) + "_l" + to_string( iLoop );
        if ( isHistoManagerAvailable ) {
          histoInfo = histoMgr->getHistogramInfo(_commonModeHistoName);
          if ( histoInfo ) {
            streamlog_out ( DEBUG1 )  << (* histoInfo ) << endl;
            commonModeHistoNBin = histoInfo->_xBin;
            commonModeHistoMin  = histoInfo->_xMin;
            commonModeHistoMax  = histoInfo->_xMax;
          }
        }
        AIDA::IHistogram1D * commonModeHisto =
          AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                    commonModeHistoNBin, commonModeHistoMin, commonModeHistoMax);
        if ( commonModeHisto ) {
          _aidaHistoMap.insert(make_pair(tempHistoName, commonModeHisto));
          commonModeHisto->setTitle("Common mode distribution");
        } else {
          streamlog_out ( ERROR1 ) << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                   << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
          _histogramSwitch = false;
        }
      }

      // book a 2d histogram for pedestal map
      const int    xNoOfPixel = abs( _maxX[iDetector] - _minX[iDetector] + 1);
      const int    yNoOfPixel = abs( _maxY[iDetector] - _minY[iDetector] + 1);
      const double xMin       = _minX[iDetector] - 0.5;
      const double xMax       =  xMin + xNoOfPixel;
      const double yMin       = _minY[iDetector] - 0.5;
      const double yMax       =  yMin + yNoOfPixel;
      tempHistoName = _pedeMapHistoName + "_d" + to_string( _orderedSensorIDVec.at( iDetector ) ) + "_l" + to_string( iLoop );
      AIDA::IHistogram2D * pedeMapHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(),
                                                                  xNoOfPixel, xMin, xMax, yNoOfPixel, yMin, yMax);
      if ( pedeMapHisto ) {
        _aidaHistoMap.insert(make_pair(tempHistoName, pedeMapHisto));
        pedeMapHisto->setTitle("Pedestal map");
      } else {
        streamlog_out ( ERROR1 ) << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                 << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }

      // book a 2d histogram for noise map
      tempHistoName = _noiseMapHistoName + "_d" + to_string( _orderedSensorIDVec.at( iDetector ) ) + "_l" + to_string( iLoop );
      AIDA::IHistogram2D * noiseMapHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(),
                                                                  xNoOfPixel, xMin, xMax, yNoOfPixel, yMin, yMax);
      if ( noiseMapHisto ) {
        _aidaHistoMap.insert(make_pair(tempHistoName, noiseMapHisto));
        noiseMapHisto->setTitle("Noise map");
      } else {
        streamlog_out ( ERROR1 ) << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                 << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }

      // book a 2d histogram for status map
      tempHistoName = _statusMapHistoName + "_d" + to_string( _orderedSensorIDVec.at( iDetector ) ) + "_l" + to_string( iLoop );
      AIDA::IHistogram2D * statusMapHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(),
                                                                  xNoOfPixel, xMin, xMax, yNoOfPixel, yMin, yMax);
      if ( statusMapHisto ) {
        _aidaHistoMap.insert(make_pair(tempHistoName, statusMapHisto));
        statusMapHisto->setTitle("Status map");
      } else {
        streamlog_out ( ERROR1 ) << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                 << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }

      if ( ( _pedestalAlgo == EUTELESCOPE::AIDAPROFILE ) &&
           ( iLoop == 0 ) ) {
        // we just need to prepare such a 2d profile only in the case
        // we are using the AIDAPROFILE calculation algorithm and we
        // need just one copy of it for each detector.
        tempHistoName = _tempProfile2DName + "_d" + to_string( _orderedSensorIDVec.at( iDetector ) ) ;
        AIDA::IProfile2D * tempProfile2D =
          AIDAProcessor::histogramFactory(this)->createProfile2D( (basePath + tempHistoName).c_str(),
                                                                  xNoOfPixel, xMin, xMax, yNoOfPixel, yMin, yMax,-1000,1000);
        if ( tempProfile2D ) {
          _aidaHistoMap.insert(make_pair(tempHistoName, tempProfile2D));
          tempProfile2D->setTitle("Temp profile for pedestal calculation");
        } else {
          streamlog_out ( ERROR1 ) << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                   << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
          _histogramSwitch = false;
          exit(-1);
        }
      }

    }
  } // end on iLoop

  if ( _additionalMaskingLoop ) {

    int iLoop = _noOfCMIterations + 1 ;
    string loopDirName = "loop_" + to_string( iLoop ) ;
    AIDAProcessor::tree(this)->mkdir(loopDirName.c_str());

    for ( size_t iDetector = 0; iDetector < _noOfDetector; iDetector++ ) {

      // prepare the name of the current detector and add it to the
      // current ITree inside the current loop folder
      string detectorDirName = "detector_" + to_string( _orderedSensorIDVec.at( iDetector ) ) ;
      string basePath = loopDirName + "/" + detectorDirName + "/";
      AIDAProcessor::tree(this)->mkdir(basePath.c_str());

      // book an histogram for the firing frequency
      const int    fireFreqHistoNBin   = 100;
      const double fireFreqHistoMin    =  0.0;
      const double fireFreqHistoMax    =  100;
      tempHistoName = _fireFreqHistoName + "_d" + to_string( _orderedSensorIDVec.at( iDetector ) ) + "_l" + to_string( iLoop );
      AIDA::IHistogram1D * fireFreqHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                  fireFreqHistoNBin, fireFreqHistoMin, fireFreqHistoMax);
      if ( fireFreqHisto ) {
        _aidaHistoMap.insert(make_pair(tempHistoName, fireFreqHisto));
        fireFreqHisto->setTitle("Firing frequency distribution");
      } else {
        streamlog_out ( ERROR1 ) << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                 << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }

      // preparing a histogram containing the signal on all the events
      tempHistoName = _aPixelHistoName + "_d" + to_string( _orderedSensorIDVec.at( iDetector ) ) + "_l" + to_string( iLoop );
      const int    aPixelNBin =    1000;
      const double aPixelMin  =  -499.5;
      const double aPixelMax  =   500.5;
      AIDA::IHistogram1D * aPixelHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName ).c_str(),
                                                                  aPixelNBin, aPixelMin, aPixelMax );
      if ( aPixelHisto ) {
        _aidaHistoMap.insert( make_pair( tempHistoName, aPixelHisto ) );
        aPixelHisto->setTitle("Output signal of a pixel");
      } else {
        streamlog_out ( ERROR1 ) << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                 << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }

      // book a 2d histogram for status map
      const int    xNoOfPixel = abs( _maxX[iDetector] - _minX[iDetector] + 1);
      const int    yNoOfPixel = abs( _maxY[iDetector] - _minY[iDetector] + 1);
      const double xMin       = _minX[iDetector] - 0.5;
      const double xMax       =  xMin + xNoOfPixel;
      const double yMin       = _minY[iDetector] - 0.5;
      const double yMax       =  yMin + yNoOfPixel;
      tempHistoName = _statusMapHistoName + "_d" + to_string( _orderedSensorIDVec.at( iDetector ) ) + "_l" + to_string( iLoop );
      AIDA::IHistogram2D * statusMapHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(),
                                                                  xNoOfPixel, xMin, xMax, yNoOfPixel, yMin, yMax);
      if ( statusMapHisto ) {
        _aidaHistoMap.insert(make_pair(tempHistoName, statusMapHisto));
        statusMapHisto->setTitle("Status map");
      } else {
        streamlog_out ( ERROR1 ) << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                 << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }

    }
  }


#endif // MARLIN_USE_AIDA

}


void EUTelPedestalNoiseProcessor::simpleRewind() {

  _isFirstEvent = true;
  _iEvt = 0;

  setReturnValue("IsPedestalFinished", false);
  throw RewindDataFilesException(this);

}

void EUTelPedestalNoiseProcessor::finalizeProcessor(bool fromMaskingLoop) {

  if ( _iLoop > 0 ) {
    _skippedEventList.sort();
    _nextEventToSkip = _skippedEventList.begin();
    streamlog_out( MESSAGE4 ) << "Skipped " << _skippedEventList.size() << " event because of common mode ("
                              << static_cast< double > ( _skippedEventList.size() ) / _iEvt * 100
                              << "%)" << endl;
  }

  if ( ! fromMaskingLoop ) {

    // this is the case in which the function is called from a loop
    // different from the additional masking loop.

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

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
      _pedestal.clear();
      _noise.clear();
      for ( size_t iDetector = 0; iDetector < _noOfDetector; iDetector++) {
        stringstream ss;
        ss << _tempProfile2DName << "_d" << _orderedSensorIDVec.at( iDetector );
        FloatVec tempPede;
        FloatVec tempNoise;
        for (int yPixel = _minY[iDetector]; yPixel <= _maxY[iDetector]; yPixel++) {
          for (int xPixel = _minX[iDetector]; xPixel <= _maxX[iDetector]; xPixel++) {
            if ( AIDA::IProfile2D * profile = dynamic_cast<AIDA::IProfile2D*> (_aidaHistoMap[ss.str()]) ) {
              tempPede.push_back( static_cast< float >(profile->binHeight(xPixel,yPixel)));
              // WARNING: the noise part of this algorithm is still not
              // working probably because of a bug in RAIDA implementation
              tempNoise.push_back( static_cast< float >(profile->binRms(xPixel,yPixel)));
              //cout << xPixel << " " << yPixel << " " << tempPede.back() << " " << tempNoise.back() << endl;
            } else {
              streamlog_out ( ERROR4 )  << "Problem with the AIDA temporary profile.\n"
                                        << "Sorry for quitting... " << endl;
              exit(-1);
            }
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

  } else {

    // this is instead the case we where coming from the additional
    // loop. So we don't need to manipulate the pedestal and noise vectors

    // here refill the status histoMap
    maskBadPixel();

#if defined(MARLIN_USE_AIDA) || defined(USE_AIDA)
    // fill only the status map histograms
    string tempHistoName;
    for (size_t iDetector = 0; iDetector < _noOfDetector; iDetector++) {
      int iPixel = 0;
      for (int yPixel = _minY[iDetector]; yPixel <= _maxY[iDetector]; yPixel++) {
        for (int xPixel = _minX[iDetector]; xPixel <= _maxX[iDetector]; xPixel++) {
          if ( _histogramSwitch ) {
            tempHistoName =  _statusMapHistoName + "_d" + to_string( _orderedSensorIDVec.at( iDetector ) ) + "_l" + to_string( _iLoop );
            if ( AIDA::IHistogram2D * histo = dynamic_cast<AIDA::IHistogram2D*>(_aidaHistoMap[tempHistoName]) ) {
              histo->fill(static_cast<double>(xPixel), static_cast<double>(yPixel), static_cast<double> (_status[iDetector][iPixel]));
            } else {
              streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                        << ".\nDisabling histogramming from now on " << endl;
              _histogramSwitch = false;
            }
            ++iPixel;
          }
        }
      }
    }
#endif
  }


  // increment the loop counter
  ++_iLoop;

  // check if we need another loop or we can finish. Remember that we
  // have a total number of loop of _noOfCMIteration + 1 + eventually
  // the additional loop on bad pixel masking

  int additionalLoop = 0;
  if ( _additionalMaskingLoop ) additionalLoop = 1;
  if ( _iLoop == _noOfCMIterations + 1 + additionalLoop ) {

    // ok this was last loop whatever kind of loop (first, other or
    // additional) it was.

    streamlog_out ( MESSAGE4 ) << "Writing the output condition file" << endl;

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

    for ( size_t iDetector = 0; iDetector < _noOfDetector; iDetector++) {

      TrackerDataImpl    * pedestalMatrix = new TrackerDataImpl;
      TrackerDataImpl    * noiseMatrix    = new TrackerDataImpl;
      TrackerRawDataImpl * statusMatrix   = new TrackerRawDataImpl;

      CellIDEncoder<TrackerDataImpl>    idPedestalEncoder(EUTELESCOPE::MATRIXDEFAULTENCODING, pedestalCollection);
      CellIDEncoder<TrackerDataImpl>    idNoiseEncoder(EUTELESCOPE::MATRIXDEFAULTENCODING, noiseCollection);
      CellIDEncoder<TrackerRawDataImpl> idStatusEncoder(EUTELESCOPE::MATRIXDEFAULTENCODING, statusCollection);

      idPedestalEncoder["sensorID"] = _orderedSensorIDVec.at( iDetector );
      idNoiseEncoder["sensorID"]    = _orderedSensorIDVec.at( iDetector );
      idStatusEncoder["sensorID"]   = _orderedSensorIDVec.at( iDetector );
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

      if ( _asciiOutputSwitch ) {
        if ( iDetector == 0 ) streamlog_out ( MESSAGE4 ) << "Writing the ASCII pedestal files" << endl;
        stringstream ss;
        ss << _outputPedeFileName << "-b" << iDetector << ".dat";
        ofstream asciiPedeFile(ss.str().c_str());
        asciiPedeFile << "# Pedestal and noise for board number " << iDetector << endl
                      << "# calculated from run " << _outputPedeFileName << endl;

        const int subMatrixWidth = 3;
        const int xPixelWidth    = 4;
        const int yPixelWidth    = 4;
        const int pedeWidth      = 15;
        const int noiseWidth     = 15;
        const int statusWidth    = 3;
        const int precision      = 8;

        int iPixel = 0;
        for (int yPixel = _minY[iDetector]; yPixel <= _maxY[iDetector]; yPixel++) {
          for (int xPixel = _minX[iDetector]; xPixel <= _maxX[iDetector]; xPixel++) {
            asciiPedeFile << setiosflags(ios::left)
                          << setw(subMatrixWidth) << iDetector
                          << setw(xPixelWidth)    << xPixel
                          << setw(yPixelWidth)    << yPixel
                          << resetiosflags(ios::left) << setiosflags(ios::fixed) << setprecision(precision)
                          << setw(pedeWidth)      << _pedestal[iDetector][iPixel]
                          << setw(noiseWidth)     << _noise[iDetector][iPixel]
                          << resetiosflags(ios::fixed)
                          << setw(statusWidth)    << _status[iDetector][iPixel]
                          << endl;
            ++iPixel;
          }
        }
        asciiPedeFile.close();
      }
    }

    event->addCollection(pedestalCollection, _pedestalCollectionName);
    event->addCollection(noiseCollection, _noiseCollectionName);
    event->addCollection(statusCollection, _statusCollectionName);

    lcWriter->writeEvent(event);
    delete event;

    lcWriter->close();

    throw StopProcessingException(this);
    setReturnValue("IsPedestalFinished", true);

  } else if ( _iLoop < _noOfCMIterations + 1 ) {


    // now we need to loop again
    // so reset the event counter
    _iEvt = 0;

    // prepare everything for the next loop
    if ( _pedestalAlgo == EUTELESCOPE::MEANRMS ) {

      // the collection contains several TrackerRawData
      // move back the _pedestal and _noise to the _temp vector
      _tempPede  = _pedestal;
      _tempNoise = _noise;

      for ( size_t iDetector = 0; iDetector < _noOfDetector; iDetector++) {
        _tempEntries.push_back(IntVec( _noise[iDetector].size(), 1));
      }

    } else if ( _pedestalAlgo == EUTELESCOPE::AIDAPROFILE ) {
      // in case the AIDAPROFILE algorithm is used, the only thing we
      // need to do is to clean up the previous loop histograms
      // remember to loop over all detectors
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
      for ( size_t iDetector = 0; iDetector < _noOfDetector; iDetector++) {
        stringstream ss;
        ss << _tempProfile2DName << "_d" << _orderedSensorIDVec.at( iDetector );
        if ( AIDA::IProfile2D* profile = dynamic_cast<AIDA::IProfile2D*> (_aidaHistoMap[ss.str()]) )
          profile->reset();
        else {
          streamlog_out ( ERROR4 ) << "Unable to reset the AIDA temporary profile.\n"
                                   << "Sorry for quitting..." << endl;
          exit(-1);
        }
      }
#endif
    }
    setReturnValue("IsPedestalFinished", false);
    throw RewindDataFilesException(this);
  } else if ( ( _additionalMaskingLoop ) &&
              ( _iLoop ==  _noOfCMIterations + 1 ) ) {
    // additional loop!
    // now we need to loop again
    // so reset the event counter
    _iEvt = 0;
    setReturnValue("IsPedestalFinished", false);
    throw RewindDataFilesException(this);

  }
}

void EUTelPedestalNoiseProcessor::additionalMaskingLoop(LCEvent * event) {

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);

  if ( evt->getEventType() == kEORE ) finalizeProcessor( true );
  if ( ( _lastEvent != -1 ) && ( _iEvt >= _lastEvent ) ) finalizeProcessor( true );
  if ( _iEvt < _firstEvent ) {
    ++_iEvt;
    throw SkipEventException(this);
  }

  if ( *_nextEventToSkip == _iEvt ) {
    streamlog_out( MESSAGE4 ) << "Event " << _iEvt << " is skipped because labelled bad by the common mode procedure." << endl;
    ++_nextEventToSkip;
    return;
  }

  // let me get the rawDataCollection. This is should contain a TrackerRawDataObject
  // for each detector plane in the telescope.
  for ( size_t iCol = 0 ; iCol < _rawDataCollectionNameVec.size(); ++iCol ) {

    size_t detectorOffset = ( iCol == 0 ) ? 0 : _noOfDetectorVec.at( iCol - 1 );

    try {
      LCCollectionVec *collectionVec = dynamic_cast < LCCollectionVec * >(evt->getCollection (_rawDataCollectionNameVec.at( iCol )));

      for ( size_t iDetector = 0; iDetector < collectionVec->size(); iDetector++) {
        // get the TrackerRawData object from the collection for this detector
        TrackerRawData *trackerRawData = dynamic_cast < TrackerRawData * >(collectionVec->getElementAt (iDetector));
        ShortVec adcValues = trackerRawData->getADCValues ();
        for ( unsigned int iPixel = 0 ; iPixel < adcValues.size(); iPixel++ ) {
          if ( _status[iDetector + detectorOffset][iPixel] == EUTELESCOPE::GOODPIXEL ) {
            float correctedValue = adcValues[iPixel] - _pedestal[iDetector + detectorOffset][iPixel];
            float threshold      = _noise[iDetector + detectorOffset][iPixel] * 3.0 ;
#if defined(MARLIN_USE_AIDA) || defined(USE_AIDA)
            if ( _histogramSwitch  && iPixel == 1 + (adcValues.size() / 10)) {
              string tempHistoName = _aPixelHistoName + "_d" + to_string( _orderedSensorIDVec.at( iDetector  + detectorOffset ) ) + "_l" + to_string( _iLoop ) ;
              if ( AIDA::IHistogram1D * histo = dynamic_cast< AIDA::IHistogram1D*> ( _aidaHistoMap[ tempHistoName ] ) )
                histo->fill( correctedValue );
              else {
                streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                          << ".\nDisabling histogramming from now on " << endl;
                _histogramSwitch = false;
              }
            }
#endif
            if ( correctedValue > threshold ) {
              _hitCounter[iDetector + detectorOffset][iPixel]++;
            }
          }
        }
      }
    } catch (DataNotAvailableException& e) {
      streamlog_out ( WARNING2 ) << "No input collection " << _rawDataCollectionNameVec.at( iCol ) << " is not available in the current event" << endl;
    }
    ++_iEvt;
  }
}


void EUTelPedestalNoiseProcessor::setBadPixelAlgoSwitches() {

  if ( find( _badPixelAlgoVec.begin(), _badPixelAlgoVec.end(), EUTELESCOPE::NOISEDISTRIBUTION ) != _badPixelAlgoVec.end() ) {
    _badPixelNoiseDistributionSwitch = true;
    streamlog_out( DEBUG4 ) << "Using " << EUTELESCOPE::NOISEDISTRIBUTION << " algorithm" << endl;
  } else {
    _badPixelNoiseDistributionSwitch = false;
  }

  if ( find( _badPixelAlgoVec.begin(), _badPixelAlgoVec.end(), EUTELESCOPE::ABSOLUTENOISEVALUE ) != _badPixelAlgoVec.end() ) {
    _badPixelAbsNoiseSwitch = true;
    streamlog_out( DEBUG4 ) << "Using " << EUTELESCOPE::ABSOLUTENOISEVALUE << " algorithm" << endl;
  } else {
    _badPixelAbsNoiseSwitch = false;
  }

  if ( find( _badPixelAlgoVec.begin(), _badPixelAlgoVec.end(), EUTELESCOPE::DEADPIXEL ) != _badPixelAlgoVec.end() ) {
    _badPixelDeadPixelSwitch = true;
    streamlog_out( DEBUG4 ) << "Using " << EUTELESCOPE::DEADPIXEL << " algorithm" << endl;
  } else {
    _badPixelDeadPixelSwitch = false;
  }

  if ( find( _badPixelAlgoVec.begin(), _badPixelAlgoVec.end(), EUTELESCOPE::ABSOLUTEPEDEVALUE ) != _badPixelAlgoVec.end() ) {
    _badPixelAbsPedeSwitch = true;
    streamlog_out( DEBUG4 ) << "Using " << EUTELESCOPE::ABSOLUTEPEDEVALUE << " algorithm" << endl;
  } else {
    _badPixelAbsPedeSwitch = false;
  }


}

