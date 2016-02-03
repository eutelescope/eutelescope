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
#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"
#include "EUTelCalibrateEventProcessor.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelHistogramManager.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
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
#include <iostream>
#include <iomanip>
#include <memory>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
string EUTelCalibrateEventProcessor::_rawDataDistHistoName            = "RawDataDistHisto";
string EUTelCalibrateEventProcessor::_dataDistHistoName               = "DataDistHisto";
string EUTelCalibrateEventProcessor::_commonModeDistHistoName         = "CommonModeDistHisto";
string EUTelCalibrateEventProcessor::_skippedRowDistHistoName         = "SkippedRowDistHisto";
string EUTelCalibrateEventProcessor::_skippedPixelDistHistoName       = "SkippedPixelDistHisto";
string EUTelCalibrateEventProcessor::_skippedPixelPerRowDistHistoName = "SkippedPixelPerRowDistHisto";
#endif

const unsigned short EUTelCalibrateEventProcessor::_maxNoOfConsecutiveMissing = 10;

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
                              "Flag to switch on the common mode suppression algorithm. 0 -> off, 1 -> full frame,  2 -> row wise",
                              _doCommonMode, static_cast<int> (1));

  registerProcessorParameter ("HitRejectionCut",
                              "Threshold of pixel SNR for hit rejection",
                              _hitRejectionCut, static_cast<float> (3.5));

  registerProcessorParameter ("MaxNoOfRejectedPixels",
                              "Maximum allowed number of rejected pixel per event",
                              _maxNoOfRejectedPixels, static_cast<int> (3000));

  registerProcessorParameter("MaxNoOfRejectedPixelPerRow",
                             "Maximum allowed number of rejected pixels per row (only with RowWise)",
                             _maxNoOfRejectedPixelPerRow,
                             static_cast < int > (25) );

  registerProcessorParameter("MaxNoOfSkippedRow",
                             "Maximum allowed number of skipped rows (only with RowWise)",
                             _maxNoOfSkippedRow,
                             static_cast< int > ( 15 ) );

  registerProcessorParameter("HistoInfoFileName", "This is the name of the histogram information file",
                             _histoInfoFileName, string( "histoinfo.xml" ) );
}


void EUTelCalibrateEventProcessor::init () {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  if ( _fillDebugHisto == 1 ) {
    streamlog_out( WARNING2 ) << "Filling debug histograms is slowing down the procedure" << endl;
  }

  // set to zero the run counter
  _iRun = 0;

  _isGeometryReady = false;

  // reset the number of consecutive missing events
  _noOfConsecutiveMissing = 0;
}

void EUTelCalibrateEventProcessor::processRunHeader (LCRunHeader * rdr) {
  std::unique_ptr<EUTelRunHeaderImpl> runHeader = std::make_unique<EUTelRunHeaderImpl>(rdr);
  runHeader->addProcessor( type());
  ++_iRun;
}

void EUTelCalibrateEventProcessor::initializeGeometry(LCEvent * event) throw ( marlin::SkipEventException ) {

  // now I need the _ancillaryIndexMap. For this I need to take the
  // pedestal input collection
  try {
    LCCollectionVec * pedestalCol = dynamic_cast< LCCollectionVec * > ( event->getCollection( _pedestalCollectionName ) );
    CellIDDecoder< TrackerDataImpl > pedestalDecoder( pedestalCol ) ;
    for ( size_t iDetector = 0; iDetector < pedestalCol->size(); ++iDetector ) {
      TrackerDataImpl * pedestal = dynamic_cast< TrackerDataImpl * > ( pedestalCol->getElementAt( iDetector ) ) ;

      int sensorID = pedestalDecoder( pedestal )["sensorID"];

      _ancillaryIndexMap.insert( make_pair( sensorID, iDetector ) );

      unsigned int noOfPixel = ( pedestalDecoder( pedestal )["xMax"] - pedestalDecoder( pedestal )["xMin"] + 1 ) *
        ( pedestalDecoder( pedestal )["yMax"] - pedestalDecoder( pedestal )["yMin"] + 1 );

      _noOfPixelMap.insert( make_pair( sensorID, noOfPixel ) );

      unsigned int noOfPixelPerRow = ( pedestalDecoder( pedestal )["xMax"] - pedestalDecoder( pedestal )["xMin"] + 1 );
      _noOfPixelPerRowMap.insert( make_pair ( sensorID, noOfPixelPerRow ));

      unsigned int noOfRow =  ( pedestalDecoder( pedestal )["yMax"] - pedestalDecoder( pedestal )["yMin"] + 1 );
      _noOfRowMap.insert( make_pair( sensorID, noOfRow ) );

    }
  } catch ( lcio::DataNotAvailableException ) {
    streamlog_out( WARNING2 ) << "Unable to initialize the geometry with the current event. Trying with the next one" << endl;
    _isGeometryReady = false;
    throw SkipEventException( this ) ;
  }

  _isGeometryReady = true;
}

void EUTelCalibrateEventProcessor::processEvent (LCEvent * event) {


  if ( !_isGeometryReady ) {
    initializeGeometry( event ) ;
  }

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);

  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return;
  } else if ( evt->getEventType() == kUNKNOWN ) {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }

  try {

    LCCollectionVec * inputCollectionVec    = dynamic_cast < LCCollectionVec * > (evt->getCollection(_rawDataCollectionName));
    LCCollectionVec * pedestalCollectionVec = dynamic_cast < LCCollectionVec * > (evt->getCollection(_pedestalCollectionName));
    LCCollectionVec * noiseCollectionVec    = dynamic_cast < LCCollectionVec * > (evt->getCollection(_noiseCollectionName));
    LCCollectionVec * statusCollectionVec   = dynamic_cast < LCCollectionVec * > (evt->getCollection(_statusCollectionName));
    CellIDDecoder<TrackerRawDataImpl> cellDecoder( inputCollectionVec );

    // reset the number of consecutive missing
    _noOfConsecutiveMissing = 0;

    if (isFirstEvent()) {

      // since v00-00-09 the consistency check between input
      // collection size and ancillary one has been removed.

      for (unsigned int iDetector = 0; iDetector < inputCollectionVec->size(); iDetector++) {

        TrackerRawDataImpl * rawData  = dynamic_cast < TrackerRawDataImpl * >(inputCollectionVec->getElementAt(iDetector));
        int sensorID = cellDecoder(rawData)["sensorID"];

        // To get the proper pedestal collection I need to go through
        // the _ancillaryIndexMap. This takes the sensorID as input and
        // gets the position in the collection as output
        TrackerDataImpl    * pedestal = dynamic_cast < TrackerDataImpl * >   (pedestalCollectionVec->getElementAt( _ancillaryIndexMap[ sensorID ] ));

        if (rawData->getADCValues().size() != pedestal->getChargeValues().size()){
          stringstream ss;
          ss << "Input data and pedestal are incompatible\n"
             << "Detector " << sensorID << " has " <<  rawData->getADCValues().size() << " pixels in the input data \n"
             << "while " << pedestal->getChargeValues().size() << " in the pedestal data " << endl;
          throw IncompatibleDataSetException(ss.str());
        }

        // in the case MARLIN_USE_AIDA and the user wants to fill detector
        // debug histogram, this is the right place to book them. As usual
        // histograms are grouped in directory identifying the detector

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
        string basePath, tempHistoName;
        basePath = "detector_" + to_string( sensorID ) + "/";

        // prepare the histogram manager
        std::unique_ptr<EUTelHistogramManager> histoMgr = std::make_unique<EUTelHistogramManager>(_histoInfoFileName);
        EUTelHistogramInfo    * histoInfo;
        bool                    isHistoManagerAvailable;

        try {
          isHistoManagerAvailable = histoMgr->init();
        } catch ( ios::failure& e) {
          streamlog_out ( ERROR1 )  << "I/O problem with " << _histoInfoFileName << "\n"
                                    << "Continuing without histogram manager"    << endl;
          isHistoManagerAvailable = false;
        } catch ( ParseException& e ) {
          streamlog_out ( ERROR1 )  << e.what() << "\n"
                                    << "Continuing without histogram manager" << endl;
          isHistoManagerAvailable = false;
        }

        if ( ( _fillDebugHisto == 1 ) ||
             ( _doCommonMode   == 1 ) ||
             ( _doCommonMode   == 2 ) ) {
          // it was changed from mkdir to mkdirs to be sure all
          // intermediate folders were properly created, but then I had
          // to come back to mkdir because this is the only supported in RAIDA.
          AIDAProcessor::tree(this)->mkdir(basePath.c_str());
        }

        if (_fillDebugHisto == 1) {
          // book the raw data histogram
          tempHistoName = _rawDataDistHistoName + "_d" + to_string( sensorID ) ;

          int    rawDataDistHistoNBin = 4096;
          double rawDataDistHistoMin  = -2048.5;
          double rawDataDistHistoMax  = rawDataDistHistoMin + rawDataDistHistoNBin;
          string rawDataDistTitle     = "Raw data distribution";
          if ( isHistoManagerAvailable ) {
            histoInfo = histoMgr->getHistogramInfo(_rawDataDistHistoName);
            if ( histoInfo ) {
              streamlog_out ( DEBUG2 ) << (* histoInfo ) << endl;
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
            streamlog_out ( ERROR1 ) << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                     << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
            _fillDebugHisto = 0;
          }

          // book the pedestal corrected data histogram
          tempHistoName = _dataDistHistoName + "_d" + to_string( sensorID );

          int    dataDistHistoNBin =  5000;
          double dataDistHistoMin  = -500.;
          double dataDistHistoMax  =  500.;
          string dataDistTitle     = "Data (pedestal sub'ed) distribution";
          if ( isHistoManagerAvailable ) {
            histoInfo = histoMgr->getHistogramInfo(_dataDistHistoName);
            if ( histoInfo ) {
              streamlog_out ( DEBUG2 ) << (* histoInfo ) << endl;
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
            streamlog_out ( ERROR1 ) << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                     << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
            _fillDebugHisto = 0;
          }

        }

        if ( _doCommonMode == 1 || _doCommonMode == 2) {
          // book the common mode histo
          tempHistoName = _commonModeDistHistoName + "_d" + to_string( sensorID );

          int    commonModeDistHistoNBin = 100;
          double commonModeDistHistoMin  = -10;
          double commonModeDistHistoMax  = +10;
          string commonModeTitle         = "Common mode distribution";
          if ( isHistoManagerAvailable ) {
            histoInfo = histoMgr->getHistogramInfo(_commonModeDistHistoName);
            if ( histoInfo ) {
              streamlog_out ( DEBUG2 ) << (* histoInfo ) << endl;
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
            streamlog_out ( ERROR1 ) << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                     << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
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


        if ( _doCommonMode == 1 ) {

          // book full frame common mode histograms
          tempHistoName  = _skippedPixelDistHistoName + "_d" + to_string( sensorID );

          int     skippedPixelDistHistoNBin = static_cast< int > ( _noOfPixelMap[ sensorID ] / 100. );
          double  skippedPixelDistHistoMin  = 0.;
          double  skippedPixelDistHistoMax  = static_cast< double > ( _noOfPixelMap[ sensorID ] );
          string  skippedPixelDistTitle     = "Pixel exceeding the hit rejection cut";

          AIDA::IHistogram1D * skippedPixelDistHisto =
            AIDAProcessor::histogramFactory( this )->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                        skippedPixelDistHistoNBin, skippedPixelDistHistoMin,skippedPixelDistHistoMax) ;
          if ( skippedPixelDistHisto ) {
            _aidaHistoMap.insert( make_pair( tempHistoName, skippedPixelDistHisto ) );
            skippedPixelDistHisto->setTitle( skippedPixelDistTitle );
          } else {
            streamlog_out ( ERROR1 ) << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                     << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
          }
        }


        if ( _doCommonMode == 2 ) {

          // book row wise common mode histograms
          tempHistoName = _skippedPixelPerRowDistHistoName + "_d" + to_string( sensorID );

          int     skippedPixelPerRowDistHistoNBin = _noOfPixelPerRowMap[ sensorID ] + 1 ;
          double  skippedPixelPerRowDistHistoMin  = -0.5;
          double  skippedPixelPerRowDistHistoMax  = _noOfPixelPerRowMap[ sensorID ] + 0.5 ;
          string  skippedPixelPerRowDistTitle     = "Pixel in a row exceeding the hit rejection cut";

          AIDA::IHistogram1D * skippedPixelPerRowDistHisto =
            AIDAProcessor::histogramFactory( this )->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                        skippedPixelPerRowDistHistoNBin, skippedPixelPerRowDistHistoMin, skippedPixelPerRowDistHistoMax );
          if ( skippedPixelPerRowDistHisto ) {
            _aidaHistoMap.insert( make_pair( tempHistoName, skippedPixelPerRowDistHisto ) );
            skippedPixelPerRowDistHisto->setTitle( skippedPixelPerRowDistTitle );
          } else {
            streamlog_out ( ERROR1 ) << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                     << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
          }


          tempHistoName = _skippedRowDistHistoName + "_d" + to_string( sensorID );

          int     skippedRowDistHistoNBin = _noOfRowMap[ sensorID ] + 1 ;
          double  skippedRowDistHistoMin  = -0.5;
          double  skippedRowDistHistoMax  = _noOfRowMap[ sensorID ] + 0.5 ;
          string  skippedRowDistTitle     = "Row rejected in common mode calculation";

          AIDA::IHistogram1D * skippedRowDistHisto = 
            AIDAProcessor::histogramFactory( this )->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                        skippedRowDistHistoNBin, skippedRowDistHistoMin, skippedRowDistHistoMax );
          if ( skippedRowDistHisto ) {
            _aidaHistoMap.insert( make_pair( tempHistoName, skippedRowDistHisto ));
            skippedRowDistHisto->setTitle( skippedRowDistTitle ) ;
          } else {
            streamlog_out ( ERROR1 ) << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                     << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
          }
        }

#endif
      }
      _isFirstEvent = false;
    }

    LCCollectionVec * correctedDataCollection = new LCCollectionVec(LCIO::TRACKERDATA);

    _minX.clear();
    _maxX.clear();
    _minY.clear();
    _maxY.clear();

    for (unsigned int iDetector = 0; iDetector < inputCollectionVec->size(); iDetector++) {
      vector< float > commonModeCorVec;
      commonModeCorVec.clear();

      // reset quantity for the common mode.
      double pixelSum      = 0.;
      double commonMode    = 0.;
      int    goodPixel     = 0;
      int    skippedPixel  = 0;
      int    skippedRow   = 0;


      TrackerRawDataImpl  * rawData   = dynamic_cast < TrackerRawDataImpl * >(inputCollectionVec->getElementAt(iDetector));
      int sensorID                    = cellDecoder(rawData)["sensorID"];

      // this is the corresponding element in the ancillary collections
      size_t ancillaryPos = _ancillaryIndexMap[ sensorID ];
      TrackerDataImpl     * pedestal  = dynamic_cast < TrackerDataImpl * >   (pedestalCollectionVec->getElementAt( ancillaryPos ));
      TrackerDataImpl     * noise     = dynamic_cast < TrackerDataImpl * >   (noiseCollectionVec->getElementAt( ancillaryPos ));
      TrackerRawDataImpl  * status    = dynamic_cast < TrackerRawDataImpl * >(statusCollectionVec->getElementAt( ancillaryPos ));

      TrackerDataImpl     * corrected = new TrackerDataImpl;
      CellIDEncoder<TrackerDataImpl> idDataEncoder(EUTELESCOPE::MATRIXDEFAULTENCODING, correctedDataCollection);
      idDataEncoder["sensorID"] = sensorID;
      idDataEncoder["xMin"]     = static_cast<int > (cellDecoder(rawData)["xMin"]);
      idDataEncoder["xMax"]     = static_cast<int > (cellDecoder(rawData)["xMax"]);
      idDataEncoder["yMin"]     = static_cast<int > (cellDecoder(rawData)["yMin"]);
      idDataEncoder["yMax"]     = static_cast<int > (cellDecoder(rawData)["yMax"]);
      _minX.push_back( cellDecoder( rawData ) ["xMin"] ) ;
      _maxX.push_back( cellDecoder( rawData ) ["xMax"] ) ;
      _minY.push_back( cellDecoder( rawData ) ["yMin"] ) ;
      _maxY.push_back( cellDecoder( rawData ) ["yMax"] ) ;

      idDataEncoder.setCellID(corrected);

      ShortVec::const_iterator rawIter     = rawData->getADCValues().begin();
      FloatVec::const_iterator pedIter     = pedestal->getChargeValues().begin();
      FloatVec::const_iterator noiseIter   = noise->getChargeValues().begin();
      ShortVec::const_iterator statusIter  = status->getADCValues().begin();


      bool isEventValid = true;
      if ( _doCommonMode == 1 ) {

        // FULLFRAME common mode
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
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
          string tempHistoName = _commonModeDistHistoName + "_d" + to_string( sensorID );
          if ( AIDA::IHistogram1D* histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName]) )
            histo->fill(commonMode);
#endif

        } else {
          isEventValid = false;
        }

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
        string tempHistoName = _skippedPixelDistHistoName + "_d" + to_string( sensorID ) ;
        if ( AIDA::IHistogram1D* histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName]) )
          histo->fill( skippedPixel );
#endif

      } else if ( _doCommonMode == 2 ) {

        // ROWWISE common mode
        ShortVec adcValues = rawData->getADCValues ();
        FloatVec _pedestal     = pedestal->getChargeValues();
        ShortVec _status  = status->getADCValues();
        FloatVec  _noise   = noise->getChargeValues();
        int    iPixel       = 0;
        int    colCounter   = 0;
        int    rowLength    = _maxX[iDetector] -  _minX[iDetector] + 1;

        for (int yPixel = _minY[iDetector]; yPixel <= _maxY[iDetector]; yPixel++) {

          double pixelSum           = 0.;
          double commonMode         = 0.;
          int    goodPixel          = 0;
          int    skippedPixelPerRow = 0;

          for ( int xPixel = _minX[iDetector]; xPixel <= _maxX[iDetector]; xPixel++) {
            bool isHit  = ( ( adcValues[iPixel] - _pedestal[iPixel] ) > _hitRejectionCut * _noise[iPixel] );
            bool isGood = ( _status[iPixel] == EUTELESCOPE::GOODPIXEL );
            if ( !isHit && isGood ) {
              pixelSum += adcValues[iPixel] - _pedestal[iPixel];
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
            string tempHistoName = _commonModeDistHistoName + "_d" + to_string( sensorID );
            if ( AIDA::IHistogram1D* histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName]) )
              histo->fill(commonMode);
#endif
          } else {
            commonModeCorVec.insert( commonModeCorVec.begin() + colCounter * rowLength, rowLength, 0. );
            ++skippedRow;
          }

          ++colCounter;
        }
        if ( skippedRow > _maxNoOfSkippedRow ) {
          isEventValid = false;
        }

      } // end if on _doCommonMode

      if(isEventValid) {
        if(_doCommonMode == 2) {
          ShortVec adcValues = rawData->getADCValues ();
          FloatVec _pedestal = pedestal->getChargeValues();

          int iPixel = 0;

          for (int yPixel = _minY[iDetector]; yPixel <= _maxY[iDetector]; yPixel++) {
            for ( int xPixel = _minX[iDetector]; xPixel <= _maxX[iDetector]; xPixel++) {

              double correctedValue = adcValues[iPixel] - _pedestal[iPixel] - commonModeCorVec[iPixel];
              corrected->chargeValues().push_back(correctedValue);

              ++iPixel;
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
              if (_fillDebugHisto == 1) {
                string tempHistoName = _rawDataDistHistoName + "_d" + to_string( sensorID );
                if ( AIDA::IHistogram1D* histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName]) ) {
                  histo->fill(adcValues[iPixel]);
                } else {
                  streamlog_out ( ERROR1 ) << "Not able to retrieve histogram pointer for " << tempHistoName
                                           << ".\nDisabling histogramming from now on " << endl;
                  _fillDebugHisto = 0 ;
                }

                tempHistoName = _dataDistHistoName + "_d" + to_string( sensorID );
                if ( AIDA::IHistogram1D * histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName]) ) {
                  histo->fill(correctedValue);
                } else {
                  streamlog_out ( ERROR1 ) << "Not able to retrieve histogram pointer for " << tempHistoName
                                           << ".\nDisabling histogramming from now on " << endl;
                  _fillDebugHisto = 0 ;
                }
              }
#endif
            }
          }
        } else {

          // that's the case the user wants to apply the FullFrame
          // common mode or doesn't want to apply any correction at
          // all. In this last case the value of the commonMode
          // variable is taken directly from the initialization ( = 0 ).

          rawIter     = rawData->getADCValues().begin();
          pedIter     = pedestal->getChargeValues().begin();

          while ( rawIter != rawData->getADCValues().end() )  {

            double correctedValue = (*rawIter) - (*pedIter) - commonMode;
            corrected->chargeValues().push_back(correctedValue);


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
            if (_fillDebugHisto == 1) {
              string tempHistoName = _rawDataDistHistoName + "_d" + to_string( sensorID );
              if ( AIDA::IHistogram1D* histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName]) ) {
                histo->fill(*rawIter);
              } else {
                streamlog_out ( ERROR1 ) << "Not able to retrieve histogram pointer for " << tempHistoName
                                         << ".\nDisabling histogramming from now on " << endl;
                _fillDebugHisto = 0 ;
              }

              tempHistoName = _dataDistHistoName + "_d" + to_string( sensorID );
              if ( AIDA::IHistogram1D * histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName]) ) {
                histo->fill(correctedValue);
              } else {
                streamlog_out ( ERROR1 ) << "Not able to retrieve histogram pointer for " << tempHistoName
                                         << ".\nDisabling histogramming from now on " << endl;
                _fillDebugHisto = 0 ;
              }
            }
#endif
            ++rawIter;
            ++pedIter;
          }
        }

      } else {
        // this is the case the event is not valid because of common
        // mode. This is the right place to throw a SkipEventException
        // possibly motivating the reason.

        if ( _doCommonMode == 1 ) {
          streamlog_out ( WARNING4 ) << "Skipping event " << evt->getEventNumber() << " because of maximum number of pixel exceeded (" << skippedPixel << ")" << endl;
        } else if ( _doCommonMode == 2 ) {
          streamlog_out ( WARNING4 ) << "Skipping event " << evt->getEventNumber() << " because of maximum number of skipped row exceeded (" << skippedRow << ")" << endl;
        } else {
          streamlog_out ( WARNING4 ) << "Skipping event " << evt->getEventNumber() << " for an unknown reason " << endl;
        }

        throw SkipEventException( this );

      }



      correctedDataCollection->push_back(corrected);
    }
    evt->addCollection(correctedDataCollection, _calibratedDataCollectionName);


  } catch (DataNotAvailableException& e) {
    if ( _noOfConsecutiveMissing <= _maxNoOfConsecutiveMissing ) {
      streamlog_out  ( WARNING2 ) <<  "No input collection found on event " << event->getEventNumber()
                                  << " in run " << event->getRunNumber() << endl;
      if ( _noOfConsecutiveMissing == _maxNoOfConsecutiveMissing ) {
        streamlog_out ( MESSAGE2 ) << "Assuming the run was taken in ZS. Not issuing any other warning" << endl;
      }
      ++_noOfConsecutiveMissing;
    }
  }

}



void EUTelCalibrateEventProcessor::check (LCEvent * /* evt */ ) {
  // nothing to check here - could be used to fill check plots in reconstruction processor
}


void EUTelCalibrateEventProcessor::end() {
  streamlog_out ( MESSAGE2 ) <<  "Successfully finished" << endl;

}

