// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelEUDRBReader.cc,v 1.1 2007-03-04 18:23:23 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// personal includes
#include "EUTelEUDRBReader.h"
#include "EUTELESCOPE.h"
#include "EUTelRunHeaderImpl.h"

// marlin includes
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
  inputFile.exceptions (ifstream::failbit | ifstream::badbit);
  
  // try to open the input file....
  try {
    inputFile.open (_fileName.c_str (), ios::in | ios::binary);
  } catch (exception & e) {
    cerr << "Problem opening file " << _fileName << ". Exiting." << endl;
    exit (-1);
  }
  
  // read the file header
  _fileHeader = new EUDRBFileHeader;

  try {
    inputFile.read(reinterpret_cast<char*>(_fileHeader), sizeof(EUDRBFileHeader));
  } catch (exception & e) {
    cerr << "Problem reading the file header" << endl;
    exit(-1);
  }

  if (isFirstEvent() ) {

    EUTelRunHeaderImpl * runHeader = new EUTelRunHeaderImpl;
    runHeader->setDAQHWName( "EUDRB" );
    runHeader->setNoOfEvent( _fileHeader->numberOfEvent );
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

    ProcessorMgr::instance()->processRunHeader( runHeader ) ;

    _isFirstEvent = false;

    delete runHeader;

    _buffer = new int[ _fileHeader->dataSize / sizeof(int) ];

  }

  for ( int iEvent = 0; iEvent < _fileHeader->numberOfEvent; iEvent++ ) {

    LCEventImpl     * event = new LCEventImpl;
    event->setDetectorName("debug_detector");
    event->setRunNumber(0);
    event->setEventNumber(iEvent);
    
    LCTime * now = new LCTime;
    event->setTimeStamp(now->timeStamp());
    delete now;

    LCCollectionVec * rawData = new LCCollectionVec (LCIO::TRACKERRAWDATA);
    CellIDEncoder < TrackerRawDataImpl > idEncoder (EUTELESCOPE::MATRIXDEFAULTENCODING, rawData);
    
    
    EUDRBEventHeader eventHeader;
    try {
      inputFile.read(reinterpret_cast<char*>(&eventHeader), sizeof(eventHeader));
    } catch (exception & e) {
      cerr << "Problem reading the event header for event " << iEvent << endl;
      exit(-1);
    }
    
    // check the event number consistency
    if ( iEvent != eventHeader.eventNumber ) {
      cerr << "Event number not corresponding " << eventHeader.eventNumber  << endl;
    }
    event->setRunNumber(0);
    event->setEventNumber(iEvent);
    
    
    try {
      inputFile.read(reinterpret_cast<char*>(_buffer), _fileHeader->dataSize );
    } catch (exception& e) {
      cerr << "Problem reading the data block for event " << iEvent << endl;
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
    
    int firstFrame, secondFrame;
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
      cerr << "Problem reading the event trailer on event " << iEvent << endl;
      exit(-1);
    }
    // crosscheck the trailer
    if (eventTrailer.trailer != 0x89abcdef ) {
      cerr << "The trailer is not correct " << endl;
    }
    
    event->addCollection(rawData, "rawdata");
    ProcessorMgr::instance()->processEvent(event);
    delete event;
  }

  inputFile.close();
}


void EUTelEUDRBReader::end () {

  delete [] _buffer;

}
