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

#ifdef EXPERIMENTAL
// personal includes
#include "EUTelEUDRBReader.h"
#include "EUTELESCOPE.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"

// marlin includes
#include "marlin/Processor.h"
#include "marlin/DataSourceProcessor.h"
#include "marlin/ProcessorMgr.h"

// lcio includes
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTime.h>
// #include <UTIL/LCTOOLS.h>

// system includes 
#include <fstream>

using namespace std;
using namespace marlin;

namespace eutelescope {
  
  //! A global instance of the processor
  EUTelEUDRBReader gEUTelEUDRBReader;

}

using namespace eutelescope;


EUTelEUDRBReader::EUTelEUDRBReader ():DataSourceProcessor  ("EUTelEUDRBReader"), _fileHeader(NULL) {
  
  _description =
    "Reads data files and creates LCEvent with TrackerRawData collection.\n"
    "Make sure to not specify any LCIOInputFiles in the steering in order to read SUCIMA Imager files.";
  
  registerProcessorParameter ("FileName", "Input file",
			      _fileName, std::string ("input.dat"));

  registerProcessorParameter ("CalculationAlgorithm", "Select if you want CDS or LF",
			      _algo, std::string("CDS"));
  
}

EUTelEUDRBReader * EUTelEUDRBReader::newProcessor () {
  return new EUTelEUDRBReader;
}

void EUTelEUDRBReader::init () {
  printParameters ();

}


void EUTelEUDRBReader::readDataSource (int numEvents) {

  ifstream inputFile;
  inputFile.exceptions (ifstream::failbit | ifstream::badbit );
  
  // try to open the input file....
  try {
    inputFile.open (_fileName.c_str (), ios::in | ios::binary);
  } catch (exception & e) {
    message<ERROR5> ( log() << "Problem opening file " << _fileName << ". Exiting." );
    exit (-1);
  }
  
  // read the file header
  _fileHeader = new EUDRBFileHeader;

  try {
    inputFile.read(reinterpret_cast<char*>(_fileHeader), sizeof(EUDRBFileHeader));
  } catch (exception & e) {
    message<ERROR5> ( log() << "Problem reading the file header" );
    exit(-1);
  }

  if (isFirstEvent() ) {

    IMPL::LCRunHeaderImpl rdr      = new IMPL::LCRunHeaderImpl;
    EUTelRunHeaderImpl * runHeader = new EUTelRunHeaderImpl(rdr)
    runHeader->setDAQHWName( "EUDRB" );
    runHeader->setNoOfEvent( _fileHeader->numberOfEvent + 1);
    runHeader->setNoOfDetector( _fileHeader->numberOfDetector * 4);
    IntVec minX, minY, maxX, maxY;
    for (int iDetector = 0; iDetector < _fileHeader->numberOfDetector * 4; iDetector++) {
      minX.push_back( ( _fileHeader->nXPixel ) * iDetector );
      maxX.push_back( ( _fileHeader->nXPixel ) * iDetector +  ( _fileHeader->nXPixel - 1 ) );
      minY.push_back( 0 );
      maxY.push_back( _fileHeader->nYPixel - 1 );
    }
    runHeader->setMinX( minX );
    runHeader->setMaxX( maxX );
    runHeader->setMinY( minY );
    runHeader->setMaxY( maxY );

    ProcessorMgr::instance()->processRunHeader( rdr ) ;

    _isFirstEvent = false;

    delete runHeader;
    delete rdr;

    _buffer = new int[ _fileHeader->dataSize / sizeof(int) ];

  }

  int iEvent;
  for ( iEvent = 0; iEvent < _fileHeader->numberOfEvent; iEvent++ ) {

    EUTelEventImpl     * event = new EUTelEventImpl;
    event->setDetectorName("debug_detector");
    event->setRunNumber(0);
    event->setEventNumber(iEvent);
    event->setEventType(kDE);
    
    LCTime * now = new LCTime;
    event->setTimeStamp(now->timeStamp());
    delete now;

    LCCollectionVec * rawData = new LCCollectionVec (LCIO::TRACKERRAWDATA);
    CellIDEncoder < TrackerRawDataImpl > idEncoder (EUTELESCOPE::MATRIXDEFAULTENCODING, rawData);
    
    
    EUDRBEventHeader eventHeader;
    try {
      inputFile.read(reinterpret_cast<char*>(&eventHeader), sizeof(eventHeader));
    } catch (exception & e) {
      message<ERROR5> ( log() << "Problem reading the event header for event " );
      exit(-1);
    }
    
    // check the event number consistency
    if ( iEvent != eventHeader.eventNumber ) {
      message<WARNING> ( log() << "Event number not corresponding " << eventHeader.eventNumber );
    }
    event->setRunNumber(0);
    event->setEventNumber(iEvent);
    
    
    try {
      inputFile.read(reinterpret_cast<char*>(_buffer), _fileHeader->dataSize );
    } catch (exception& e) {
      message<ERROR5> ( log() << "Problem reading the data block for event " << iEvent );
      exit(-1);
    }
    
      
    // this is made between frame 3 and frame 2
    int frameRecordSize = _fileHeader->nXPixel * _fileHeader->nYPixel * 4 /*frame*/ / 2 /*pixel per record*/;
    
    TrackerRawDataImpl * channelA = new TrackerRawDataImpl;
    idEncoder["sensorID"] = 0;
    idEncoder["xMin"]     = 0;
    idEncoder["xMax"]     = _fileHeader->nXPixel - 1;
    idEncoder["yMin"]     = 0;
    idEncoder["yMax"]     = _fileHeader->nYPixel - 1;
    idEncoder.setCellID( channelA );
    
    TrackerRawDataImpl * channelB = new TrackerRawDataImpl;
    idEncoder["sensorID"] = 1;
    idEncoder["xMin"]     = _fileHeader->nXPixel;
    idEncoder["xMax"]     = 2 * _fileHeader->nXPixel - 1;
    idEncoder["yMin"]     = 0;
    idEncoder["yMax"]     = _fileHeader->nYPixel - 1;
    idEncoder.setCellID( channelB );
      
    TrackerRawDataImpl * channelC = new TrackerRawDataImpl;
    idEncoder["sensorID"] = 2;
    idEncoder["xMin"]     = 2 * _fileHeader->nXPixel;
    idEncoder["xMax"]     = 3 * _fileHeader->nXPixel - 1;
    idEncoder["yMin"]     = 0;
    idEncoder["yMax"]     = _fileHeader->nYPixel - 1;
    idEncoder.setCellID( channelC );
    
    TrackerRawDataImpl * channelD = new TrackerRawDataImpl;
    idEncoder["sensorID"] = 3;
    idEncoder["xMin"]     = 3 * _fileHeader->nXPixel;
    idEncoder["xMax"]     = 4 * _fileHeader->nXPixel - 1;
    idEncoder["yMin"]     = 0;
    idEncoder["yMax"]     = _fileHeader->nYPixel - 1;
    idEncoder.setCellID( channelD );
    
    int firstFrame = -1;
    int secondFrame = -1;
    if ( ( _algo == "CDS32" ) || ( _algo == "LF2") ) {
      firstFrame  = 1;
      secondFrame = 2;
    } else if ( ( _algo == "CDS21" ) || ( _algo == "LF1" ) ) {
      firstFrame  = 0;
      secondFrame = 1;
    } else if ( _algo == "LF3" ) {
      firstFrame = 2;
      secondFrame = 3;
    } 
    
    for (int iRecord = firstFrame * frameRecordSize; iRecord < secondFrame * frameRecordSize; iRecord++ ) {
      
      short pixelA1 = static_cast< short > ( ( _buffer[iRecord] & _fileHeader->chACBitMask ) >> _fileHeader->chACRightShift ) ;
      short pixelB1 = static_cast< short > ( ( _buffer[iRecord] & _fileHeader->chBDBitMask ) >> _fileHeader->chBDRightShift ) ;
      if ( _algo.compare(0, 3, "CDS") == 0 ) {
	short pixelA2 = static_cast< short > ( ( _buffer[iRecord + frameRecordSize] & _fileHeader->chACBitMask ) >> 
					       _fileHeader->chACRightShift  );
	short pixelB2 = static_cast< short > ( ( _buffer[iRecord + frameRecordSize] & _fileHeader->chBDBitMask ) >> 
					       _fileHeader->chBDRightShift  );
	channelA->adcValues().push_back( pixelA2 - pixelA1 );
	channelB->adcValues().push_back( pixelB2 - pixelB1 );
      } else if ( _algo.compare(0, 2, "LF") == 0 ) {
	channelA->adcValues().push_back( pixelA1 );
	channelB->adcValues().push_back( pixelB1 );
      }


      ++iRecord;      
      short pixelC1 = static_cast< short > ( ( _buffer[iRecord] & _fileHeader->chACBitMask ) >> _fileHeader->chACRightShift ) ;
      short pixelD1 = static_cast< short > ( ( _buffer[iRecord] & _fileHeader->chBDBitMask ) >> _fileHeader->chBDRightShift ) ;
      if ( _algo.compare(0, 3, "CDS") == 0 ) {
	short pixelC2 = static_cast< short > ( ( _buffer[iRecord + frameRecordSize] & _fileHeader->chACBitMask ) >> 
					       _fileHeader->chACRightShift  );
	short pixelD2 = static_cast< short > ( ( _buffer[iRecord + frameRecordSize] & _fileHeader->chBDBitMask ) >> 
					       _fileHeader->chBDRightShift  );
	channelC->adcValues().push_back( pixelC2 - pixelC1 );
	channelD->adcValues().push_back( pixelD2 - pixelD1 );
      } else if ( _algo.compare(0, 2, "LF" ) == 0 ) {
	channelC->adcValues().push_back( pixelC1 );
	channelD->adcValues().push_back( pixelD1 );
      }

    }
    
    rawData->push_back(channelA);
    rawData->push_back(channelB);
    rawData->push_back(channelC);
    rawData->push_back(channelD);
    
    EUDRBTrailer eventTrailer;
    try {
      inputFile.read(reinterpret_cast<char*>(&eventTrailer), sizeof(eventTrailer));
    } catch (exception & e) {
      message<ERROR5> ( log() << "Problem reading the event trailer on event " << iEvent );
      exit(-1);
    }
    // crosscheck the trailer
    if (eventTrailer.trailer != 0x89abcdef ) {
      message<WARNING> ( log() << "The trailer is not correct on event " << iEvent ) ;
    }
    
    event->addCollection(rawData, "rawdata");

    ProcessorMgr::instance()->processEvent(static_cast<LCEventImpl*> (event) );
    delete event;

    if ( inputFile.eof() ) break;
  }

  // add the EORE event
  EUTelEventImpl     * event = new EUTelEventImpl;
  event->setDetectorName("debug_detector");
  event->setRunNumber(0);
  event->setEventNumber(iEvent);
  event->setEventType(kEORE);
  LCTime * now = new LCTime;
  event->setTimeStamp(now->timeStamp());
  delete now;
  event->setRunNumber(0);
  event->setEventNumber(iEvent + 1);

  ProcessorMgr::instance()->processEvent(static_cast<LCEventImpl*> (event) );
  delete event;

  inputFile.close();
}


void EUTelEUDRBReader::end () {

  delete [] _buffer;
  message<MESSAGE5> ( "Successfully finished" );

}
#endif
