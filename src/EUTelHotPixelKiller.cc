// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelHotPixelKiller.cc,v 1.5 2009-07-15 17:21:28 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelHotPixelKiller.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelMatrixDecoder.h"

// marlin includes ".h"
#include "marlin/Processor.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include "marlin/AIDAProcessor.h"
#include <AIDA/ITree.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IHistogramFactory.h>
#endif

// lcio includes <.h>
#include <LCIOTypes.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>

// system includes <>
#include <map>
#include <memory>

using namespace std;
using namespace marlin;
using namespace eutelescope;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
string EUTelHotPixelKiller::_firing2DHistoName = "Firing2D";
string EUTelHotPixelKiller::_firing1DHistoName = "Firing1D";
#endif

EUTelHotPixelKiller::EUTelHotPixelKiller () : Processor("EUTelHotPixelKiller") {

  // modify processor description
  _description =
    "EUTelHotPixelKiller periodically check for pixel singing loud too often and remove them from the analysis";

  registerInputCollection (LCIO::TRACKERRAWDATA, "StatusCollectionName",
                           "Pixel status collection",
                           _statusCollectionName, string("status"));

  registerProcessorParameter("NoOfEventPerCycle",
                             "The number of events to be considered for each update cycle",
                             _noOfEventPerCycle, static_cast<int>( 100 ) );

  registerProcessorParameter("MaxAllowedFiringFreq",
                             "This float number [0,1] represents the maximum allowed firing frequency\n"
                             "within the selected number of event per cycle",
                             _maxAllowedFiringFreq, static_cast<float> (0.2) );

  registerProcessorParameter("TotalNoOfCycle",
                             "The total number of hot pixel cycle",
                             _totalNoOfCycle, static_cast<int> ( 10 ) );
}


void EUTelHotPixelKiller::init () {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

  // reset the cycle number
  _iCycle = 0;

  // reset the vector with killed pixels
  _killedPixelVec.clear();

  // reset the vector with the firing frequency
  _firingFreqVec.clear();

}

void EUTelHotPixelKiller::processRunHeader (LCRunHeader * rdr ) {

  auto_ptr<EUTelRunHeaderImpl> runHeader ( new EUTelRunHeaderImpl ( rdr ) );
  runHeader->addProcessor( type() );

  // increment the run counter
  ++_iRun;

  // reset the event counter
  _iEvt = 0;

}

void EUTelHotPixelKiller::initializeGeometry( LCEvent * event ) {

  LCCollectionVec * collection = dynamic_cast< LCCollectionVec *> ( event->getCollection( _statusCollectionName ) );

  _noOfDetectors = collection->size();

  CellIDDecoder<TrackerRawDataImpl > decoder( collection );

  for ( size_t iDetector = 0 ; iDetector < collection->size() ; ++iDetector ) {

    TrackerRawDataImpl * status = dynamic_cast< TrackerRawDataImpl * > ( collection->getElementAt( iDetector ) ) ;
    int sensorID = decoder( status ) [ "sensorID" ] ;

    _minX[ sensorID ] = decoder( status ) [ "xMin" ];
    _minY[ sensorID ] = decoder( status ) [ "yMin" ];
    _maxX[ sensorID ] = decoder( status ) [ "xMax" ];
    _maxY[ sensorID ] = decoder( status ) [ "yMax" ];

    _sensorIDVec.push_back( sensorID );

  }
}

void EUTelHotPixelKiller::processEvent (LCEvent * event) {

  if ( _iCycle > (unsigned short) _totalNoOfCycle )  return;

  if (_iEvt % 10 == 0)
    streamlog_out( MESSAGE4 ) << "Processing event "
                              << setw(6) << setiosflags(ios::right) << event->getEventNumber() << " in run "
                              << setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber() << setfill(' ')
                              << " (Total = " << setw(10) << (_iCycle * _noOfEventPerCycle) + _iEvt << ")"
                              << resetiosflags(ios::left) << endl;

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);
  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG4 ) <<  "EORE found: nothing else to do." <<  endl;
    return;
  } else if ( evt->getEventType() == kUNKNOWN ) {
    streamlog_out ( WARNING2 ) << "Event number " << event->getEventNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }

  try {

    LCCollectionVec * statusCollectionVec = dynamic_cast< LCCollectionVec * > ( event->getCollection( _statusCollectionName ) );
    CellIDDecoder<TrackerRawDataImpl>      statusCellDecoder( statusCollectionVec );


    if ( isFirstEvent() ) {

      initializeGeometry( event );
      _isFirstEvent = false;

    }

    if ( _iEvt == 0 ) {

      _firingFreqVec.clear();

      for ( int iDetector = 0; iDetector < statusCollectionVec->getNumberOfElements() ; iDetector++) {
        TrackerRawDataImpl * status = dynamic_cast< TrackerRawDataImpl * > ( statusCollectionVec->getElementAt( iDetector ) );
        vector< unsigned short > dummyVector( status->getADCValues().size(), 0 );
        _firingFreqVec.push_back( dummyVector );
      }

    }

    for ( int iDetector = 0; iDetector < statusCollectionVec->getNumberOfElements() ; iDetector++) {

      TrackerRawDataImpl * status = dynamic_cast< TrackerRawDataImpl * > ( statusCollectionVec->getElementAt( iDetector ) );
      vector< short > statusVec = status->adcValues();

      for ( unsigned int iPixel = 0; iPixel < statusVec.size(); iPixel++ ) {

        if ( statusVec[ iPixel ] == EUTELESCOPE::HITPIXEL ) {
          _firingFreqVec[ iDetector ][ iPixel ]++;
        }

      }

    }
    ++_iEvt;
  } catch (lcio::DataNotAvailableException& e ) {
    streamlog_out ( WARNING2 )  << "Input collection not found in the current event. Skipping..." << endl;
    return;
  }

}


void EUTelHotPixelKiller::end() {

  streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
  streamlog_out ( MESSAGE4 ) << printSummary() << endl;

}


string EUTelHotPixelKiller::printSummary() const {

  stringstream ss;
  int bigSpacer   = 10;
  int tinySpacer  =  5;
  int smallSpacer = 12;
  int fullLineWidth = bigSpacer + tinySpacer + _noOfDetectors * smallSpacer + 1;
  stringstream doubleLine;
  stringstream singleLine;
  for ( int i = 0; i < fullLineWidth ; i++ ) {
    doubleLine << "=";
    singleLine << "-";
  }

  if ( _killedPixelVec.size() == 0 ) {
    return "" ;
  }

  ss << doubleLine.str() << endl
     << " Hot pixel killer summary " << endl
     << doubleLine.str() << endl;


  for ( unsigned int iCycle = 0 ; iCycle < _killedPixelVec[0].size(); iCycle++ ) {
    ss << " " << setiosflags( ios::left ) << setw(bigSpacer) << "Cycle num:"
       << setiosflags( ios::right ) << setw(tinySpacer) << iCycle << resetiosflags( ios::left) ;
    for ( unsigned int iDetector = 0; iDetector < _killedPixelVec.size(); iDetector++ ) {
      ss << setw(smallSpacer) << _killedPixelVec[iDetector][iCycle] ;
    }
    ss << "\n" << singleLine.str() << endl;
  }
  ss << "\n" << doubleLine.str() << endl;

  return ss.str();
}

void EUTelHotPixelKiller::check( LCEvent * event ) {

  if ( _iCycle > (unsigned short) _totalNoOfCycle )  return;

  if ( _iEvt == _noOfEventPerCycle -1 ) {


    try {

      LCCollectionVec * statusCollectionVec = dynamic_cast< LCCollectionVec * > ( event->getCollection( _statusCollectionName ) );
      CellIDDecoder<TrackerRawDataImpl>      statusCellDecoder( statusCollectionVec );

      for ( unsigned int iDetector = 0; iDetector < _firingFreqVec.size(); iDetector++ ) {

        if ( _iCycle == 0 ) {
          _killedPixelVec.push_back( vector< unsigned short > (0, 0) );
        }

        TrackerRawDataImpl * status = dynamic_cast< TrackerRawDataImpl * > ( statusCollectionVec->getElementAt( iDetector ) );
        vector< short > statusVec = status->adcValues();
        unsigned short killerCounter = 0;

        for ( unsigned int iPixel = 0; iPixel < _firingFreqVec[iDetector].size(); iPixel++ ) {

          if ( _firingFreqVec[iDetector][iPixel] / ( (double) _iEvt ) > _maxAllowedFiringFreq ) {
            streamlog_out ( DEBUG3 ) << " Pixel " << iPixel << " on detector " << _sensorIDVec.at( iDetector )
                                     << " is firing too often (" << _firingFreqVec[iDetector][iPixel] / ((double) _iEvt )
                                     << "). Masking it now on! " << endl;
            status->adcValues()[iPixel] = EUTELESCOPE::FIRINGPIXEL;
            ++killerCounter;
          }

        }

        _killedPixelVec[ iDetector ].push_back( killerCounter );
      }


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
      bookAndFillHistos();
#endif

      // increment the cycle number
      ++_iCycle;

      // reset the _iEvt counter
      _iEvt = 0;


    } catch (lcio::DataNotAvailableException& e ) {
      streamlog_out ( WARNING2 )  << "Input collection not found in the current event. Skipping..." << endl;
      return;

    }
  }
}

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

void EUTelHotPixelKiller::bookAndFillHistos() {

  streamlog_out ( MESSAGE0 ) << "Booking and filling histograms for cycle " << _iCycle << endl;

  string tempHistoName, basePath;
  for ( int iDetector = 0; iDetector < _noOfDetectors; iDetector++ ) {
    basePath = "detector_" + to_string( _sensorIDVec.at( iDetector ) ) ;
    AIDAProcessor::tree(this)->mkdir(basePath.c_str());

    basePath += "/cycle_" + to_string( _iCycle );
    AIDAProcessor::tree(this)->mkdir(basePath.c_str());
    basePath.append("/");



    tempHistoName = _firing2DHistoName + "_d" + to_string( _sensorIDVec.at(iDetector)) + "_c" + to_string( _iCycle ) ;
    int     xBin = _maxX[_sensorIDVec.at( iDetector )] - _minX[_sensorIDVec.at( iDetector )] + 1;
    double  xMin = static_cast<double >(_minX[_sensorIDVec.at( iDetector )]) - 0.5;
    double  xMax = static_cast<double >(_maxX[_sensorIDVec.at( iDetector )]) + 0.5;
    int     yBin = _maxY[_sensorIDVec.at( iDetector )] - _minY[_sensorIDVec.at( iDetector )] + 1;
    double  yMin = static_cast<double >(_minY[_sensorIDVec.at( iDetector )]) - 0.5;
    double  yMax = static_cast<double >(_maxY[_sensorIDVec.at( iDetector )]) + 0.5;

    AIDA::IHistogram2D * firing2DHisto =
      AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(),
                                                                xBin, xMin, xMax,yBin, yMin, yMax);
    firing2DHisto->setTitle("Firing frequency map");

    tempHistoName = _firing1DHistoName + "_d" + to_string( _sensorIDVec.at( iDetector) ) + "_c" + to_string( _iCycle );

    int nBin = 100;
    double min = 0;
    double max = 1;
    AIDA::IHistogram1D * firing1DHisto =
      AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                nBin, min, max );
    firing1DHisto->setTitle("Firing frequency distribution");

    int iPixel = 0;

    for (int yPixel = _minY[_sensorIDVec.at( iDetector)]; yPixel <= _maxY[_sensorIDVec.at( iDetector)]; yPixel++) {
      for (int xPixel = _minX[_sensorIDVec.at( iDetector)]; xPixel <= _maxX[_sensorIDVec.at( iDetector)]; xPixel++) {
        firing2DHisto->fill(xPixel, yPixel, _firingFreqVec[ iDetector ][ iPixel ] );
        firing1DHisto->fill( _firingFreqVec[ iDetector ][ iPixel ] / ( (double)  _noOfEventPerCycle ));
        ++iPixel;
      }
    }
  }
}

#endif
