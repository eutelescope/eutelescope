// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelStrasMimoTelReader.cc,v 1.1 2007-06-17 22:42:25 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// personal includes
#include "EUTelStrasMimoTelReader.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"

// marlin includes
#include "marlin/DataSourceProcessor.h"
#include "marlin/ProcessorMgr.h"
#include "marlin/Exceptions.h"

// lcio includes
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTime.h>
// #include <UTIL/LCTOOLS.h>

// system includes 
#include <fstream>
#include <string>
#include <sstream>

using namespace std;
using namespace marlin;
using namespace eutelescope;

string EUTelStrasMimoTelReader::_dataFileBaseName   = "RUN_";
const string EUTelStrasMimoTelReader::_fileNameExt  = ".rz" ;

EUTelStrasMimoTelReader::EUTelStrasMimoTelReader ():DataSourceProcessor  ("EUTelStrasMimoTelReader") {
  
  _description =
    "Reads the data output file produced by the Strasbourg DAQ.\n"
    "In principle it should work for each kind of detector, but since there are some\n"
    "hardcoded numbers it is working for the MimoTel only (for the time being)";
  
  registerProcessorParameter("LEPSIRunNumber","The run to be converted (the number only)",
			     _runNumber, static_cast<int> (500) );
}

EUTelStrasMimoTelReader * EUTelStrasMimoTelReader::newProcessor () {
  return new EUTelStrasMimoTelReader;
}

void EUTelStrasMimoTelReader::init () {
  printParameters ();

  // complete the _dataFileBaseName with the run number
  stringstream ss;
  ss << _dataFileBaseName << _runNumber << "_";
  _dataFileBaseName = ss.str();

}


void EUTelStrasMimoTelReader::readDataSource (int numEvents) {

  // initialize here the event counter to be compare to numEvents
  int eventCounter = 0;

  // as a first think we need to open the run header file. This is
  // named after the run number and using the _dataFileBaseName.
  string path;
  {
    stringstream ss;
    ss << _runNumber << "/";
    path=ss.str();
  }
  
  string runHeaderFileName = path + _dataFileBaseName + "i" + _fileNameExt;
  ifstream runHeaderFile;
  runHeaderFile.exceptions( ifstream::failbit | ifstream::badbit ) ;
  
  // try to open the run header file and read out the content
  try {
    runHeaderFile.open( runHeaderFileName.c_str(), ios::in | ios::binary );
    runHeaderFile.read( reinterpret_cast<char *> (&_runHeader), sizeof(StrasRunHeader) );
    
    if ( _runNumber != _runHeader.RunNo ) {
      message<WARNING> ( log() << "The run header file contains a different run number (" 
			 << _runNumber << "!=" << _runHeader.RunNo );
    }

    message<MESSAGE> ( log() << "Reading the run header\n" << _runHeader );
    runHeaderFile.close();

    EUTelRunHeaderImpl * eutelRunHeader = new EUTelRunHeaderImpl;
    eutelRunHeader->setRunNumber( _runNumber );
    eutelRunHeader->setDataType( EUTELESCOPE::CONVDATA );
    eutelRunHeader->setDateTime();
    eutelRunHeader->setDAQHWName( EUTELESCOPE::IPHCIMAGER );
    eutelRunHeader->setDAQSWName( EUTELESCOPE::IPHCIMAGER );
    eutelRunHeader->setNoOfDetector( _runHeader.VFasPresentNb + 1 );
    eutelRunHeader->setMinX(IntVec( _runHeader.VFasPresentNb + 1, 0 ) );
    eutelRunHeader->setMaxX(IntVec( _runHeader.VFasPresentNb + 1, 263 ) );
    eutelRunHeader->setMinY(IntVec( _runHeader.VFasPresentNb + 1, 0 ) );
    eutelRunHeader->setMaxY(IntVec( _runHeader.VFasPresentNb + 1, 255 ) );
    eutelRunHeader->setDetectorName("MimoTel");
    
    ProcessorMgr::instance()->processRunHeader(eutelRunHeader);
    delete eutelRunHeader;
    
  } catch (exception& e) {
    message<ERROR> ( log() << "Unable to open file " << runHeaderFileName << ". Exiting." );
    exit(-1);
  }
  
  int nFile = _runHeader.TotEvNb / _runHeader.FileEvNb;
  _dataBuffer = new int[_runHeader.DataSz / sizeof(int)];

  
  for ( int iFile = 0; iFile < nFile; iFile++  ) {
    
    string dataFileName;
    { 
      stringstream ss;
      ss << path << _dataFileBaseName << iFile << _fileNameExt;
      dataFileName = ss.str();
    }
    
    ifstream dataFile;
    dataFile.exceptions (ifstream::failbit | ifstream::badbit );
    
    // try to open the data file
    try {
      
      message<DEBUG> ( log() << "Opening file " << dataFileName );
      dataFile.open( dataFileName.c_str(), ios::in | ios::binary );
      
      while ( !dataFile.eof() ) {
	if ( (eventCounter % 10 == 0 ) )
	  message<MESSAGE> ( log() << "Converting event " << eventCounter << " (File = " << iFile << ")" );
	    
	// get the full record
	dataFile.read( reinterpret_cast<char *> (&_eventHeader), sizeof(StrasEventHeader) );
	dataFile.read( reinterpret_cast<char *> (_dataBuffer), _runHeader.DataSz );
	dataFile.read( reinterpret_cast<char *> (&_eventTrailer), sizeof(StrasEventTrailer) );
	
	// make some checks
	if ( (unsigned) _eventTrailer.Eor != 0x89ABCDEF ) {
	  message<ERROR> ( log() << "Event trailer not found on event " << _eventHeader.EvNo << ". Exiting ");
	  exit(-1);
	}
	
	if ( _eventHeader.EvNo != (unsigned) eventCounter ) {
	  message<WARNING> ( log() << "Event number mismatch: expected " << eventCounter << " read " << _eventHeader.EvNo );
	}
	
	EUTelEventImpl * event = new EUTelEventImpl;
	event->setRunNumber( _runNumber );
	event->setEventNumber( _eventHeader.EvNo );
	LCTime * now = new LCTime;
	event->setTimeStamp(now->timeStamp());
	delete now;
	event->setEventType( kDE );

	LCCollectionVec * cdsColl = new LCCollectionVec( LCIO::TRACKERRAWDATA );
	CellIDEncoder< TrackerRawDataImpl > idEncoderCDS( EUTELESCOPE::MATRIXDEFAULTENCODING, cdsColl );
	idEncoderCDS["xMin"] = 0;
	idEncoderCDS["xMax"] = 263;
	idEncoderCDS["yMin"] = 0;
	idEncoderCDS["yMax"] = 255;
	
	for ( int iDetector = 0; iDetector < _runHeader.VFasPresentNb + 1; iDetector++ ) {
	  
	  TrackerRawDataImpl * cds = new TrackerRawDataImpl;
	  idEncoderCDS["sensorID"] = iDetector;
	  
	  // that's the place to insert data

	  cdsColl->push_back( cds );
	}
	
	event->addCollection( cdsColl, "rawdata");
	ProcessorMgr::instance()->processEvent( event ) ;
	delete event;

	++eventCounter;
	if ( eventCounter > numEvents ) {
	  dataFile.close();
	  break;
	}
      }

    } catch (ifstream::failure e) {
      if ( !dataFile.eof() ) {
	message<ERROR> ( log() << "Unable to open file " << dataFileName << ". Exiting." );
	addEORE();
	exit(-1);
      } else {
	dataFile.close();
	message<DEBUG> ( log() << "Closing file " << dataFileName ) ;
      }
    } 
  }
 
  delete [] _dataBuffer;
  addEORE();


  
}

void EUTelStrasMimoTelReader::addEORE() {
  
  EUTelEventImpl * event = new EUTelEventImpl;
  event->setEventType(kEORE);
  ProcessorMgr::instance()->processEvent( event );
  delete event;
}


void EUTelStrasMimoTelReader::end () {

  message<MESSAGE> ("Successfully finished") ;
}

