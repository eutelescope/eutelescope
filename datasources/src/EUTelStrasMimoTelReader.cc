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


// personal includes
#include "EUTelStrasMimoTelReader.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"

// marlin includes
#include "marlin/Processor.h"
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
#include <algorithm>
#include <memory>

using namespace std;
using namespace marlin;
using namespace eutelescope;

string EUTelStrasMimoTelReader::_dataFileBaseName   = "RUN_";
const string EUTelStrasMimoTelReader::_fileNameExt  = ".rz" ;
int EUTelStrasMimoTelReader::_noOfSubMatrix =   1;

EUTelStrasMimoTelReader::EUTelStrasMimoTelReader ():DataSourceProcessor  ("EUTelStrasMimoTelReader") {
  
  _description =
    "Reads the data output file produced by the Strasbourg DAQ.\n"
    "In principle it should work for each kind of detector, but since there are some\n"
    "hardcoded numbers it is working for the MimoTel only (for the time being)";
  
  registerProcessorParameter("LEPSIRunNumber","The run to be converted (the number only)",
			     _runNumber, static_cast<int> (500) );

  registerProcessorParameter("XNoOfPixel","Number of pixels along the x axis",
			     _noOfXPixel, static_cast<int> ( 66 ) );
  
  registerProcessorParameter("YNoOfPixel","Number of pixels along the x axis",
			     _noOfYPixel, static_cast<int> ( 256 ) );

  registerOutputCollection(LCIO::TRACKERRAWDATA, "Frame0CollectionName",
			   "Name of the Frame0 collection",
			   _frame0CollectionName, string( "frame0" ));

  registerOutputCollection(LCIO::TRACKERRAWDATA, "Frame1CollectionName",
			   "Name of the Frame1 collection",
			   _frame1CollectionName, string( "frame1" ));

  registerOutputCollection(LCIO::TRACKERRAWDATA, "CDSCollectionName",
			   "Name of the CDS collection",
			   _cdsCollectionName, string( "cds" ));
  
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

    message<MESSAGE5> ( log() << "Reading the run header\n" << _runHeader );
    runHeaderFile.close();

    auto lcHeader = std::make_unique<IMPL::LCRunHeaderImpl>();
    auto eutelRunHeader = std::make_unique<EUTelRunHeaderImpl>(lcHeader.get());
    eutelRunHeader->addProcessor( type() );
    eutelRunHeader->lcRunHeader()->setRunNumber( _runNumber );
    eutelRunHeader->setDataType( EUTELESCOPE::CONVDATA );
    eutelRunHeader->setDateTime();
    eutelRunHeader->setDAQHWName( EUTELESCOPE::IPHCIMAGER );
    eutelRunHeader->setDAQSWName( EUTELESCOPE::IPHCIMAGER );
    eutelRunHeader->setNoOfDetector( _runHeader.VFasPresentNb + 1 );
    eutelRunHeader->setMinX(IntVec( _runHeader.VFasPresentNb + 1, 0 ) );
    eutelRunHeader->setMaxX(IntVec( _runHeader.VFasPresentNb + 1, _noOfXPixel - 1 ) );
    eutelRunHeader->setMinY(IntVec( _runHeader.VFasPresentNb + 1, 0 ) );
    eutelRunHeader->setMaxY(IntVec( _runHeader.VFasPresentNb + 1, _noOfYPixel - 1 ) );
    eutelRunHeader->lcRunHeader()->setDetectorName("MimoTel");
    
    ProcessorMgr::instance()->processRunHeader( lcHeader.release() );
    
  } catch (exception& e) {
    message<ERROR5> ( log() << "Unable to open file " << runHeaderFileName << ". Exiting." );
  }
  
  int nFile = _runHeader.TotEvNb / _runHeader.FileEvNb;
  _dataBuffer   = new int[_runHeader.DataSz / sizeof(int) ];
  int matrixSize   = _noOfXPixel * _noOfYPixel;
  
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
      
      message<DEBUG5> ( log() << "Opening file " << dataFileName );
      dataFile.open( dataFileName.c_str(), ios::in | ios::binary );
      
      while ( !dataFile.eof() ) {
	if ( (eventCounter % 10 == 0 ) )
	  message<MESSAGE5> ( log() << "Converting event " << eventCounter << " (File = " << iFile << ")" );
	    
	// get the full record
	dataFile.read( reinterpret_cast<char *> (&_eventHeader), sizeof(StrasEventHeader) );
	dataFile.read( reinterpret_cast<char *> (_dataBuffer), _runHeader.DataSz  );
	dataFile.read( reinterpret_cast<char *> (&_eventTrailer), sizeof(StrasEventTrailer) );
	
	// make some checks
	if ( static_cast< unsigned >(_eventTrailer.Eor) != 0x89ABCDEF ) {
	  message<ERROR5> ( log() << "Event trailer not found on event " << _eventHeader.EvNo << ". Exiting ");
	  exit(-1);
	}
	
	if ( _eventHeader.EvNo != static_cast< unsigned >(eventCounter) ) {
	  message<WARNING> ( log() << "Event number mismatch: expected " << eventCounter << " read " << _eventHeader.EvNo );
	}

	if ( _eventHeader.VFasCnt < 0 ) {
	  // the trigger is not accepted. Skip the event
	  message<WARNING> ( log() << "Trigger not accepted on event " << eventCounter ) ;
	  
	} else {

	  EUTelEventImpl * event = new EUTelEventImpl;
	  event->setRunNumber( _runNumber );
	  event->setEventNumber( _eventHeader.EvNo );
	  LCTime * now = new LCTime;
	  event->setTimeStamp(now->timeStamp());
	  delete now;
	  event->setEventType( kDE );
	  
	  LCCollectionVec * cdsColl    = new LCCollectionVec( LCIO::TRACKERRAWDATA );
	  LCCollectionVec * frame0Coll = new LCCollectionVec( LCIO::TRACKERRAWDATA );
	  LCCollectionVec * frame1Coll = new LCCollectionVec( LCIO::TRACKERRAWDATA );
	  CellIDEncoder< TrackerRawDataImpl > idEncoderCDS( EUTELESCOPE::MATRIXDEFAULTENCODING, cdsColl );
	  idEncoderCDS["xMin"] = 0;
	  idEncoderCDS["xMax"] = _noOfXPixel - 1;
	  idEncoderCDS["yMin"] = 0;
	  idEncoderCDS["yMax"] = _noOfYPixel - 1;
	  CellIDEncoder< TrackerRawDataImpl > idEncoderFrame0( EUTELESCOPE::MATRIXDEFAULTENCODING, frame0Coll );
	  idEncoderFrame0["xMin"] = 0;
	  idEncoderFrame0["xMax"] = _noOfXPixel - 1;
	  idEncoderFrame0["yMin"] = 0;
	  idEncoderFrame0["yMax"] = _noOfYPixel - 1;
	  CellIDEncoder< TrackerRawDataImpl > idEncoderFrame1( EUTELESCOPE::MATRIXDEFAULTENCODING, frame1Coll );
	  idEncoderFrame1["xMin"] = 0;
	  idEncoderFrame1["xMax"] = _noOfXPixel - 1;
	  idEncoderFrame1["yMin"] = 0;
	  idEncoderFrame1["yMax"] = _noOfYPixel - 1;
	  
	  // this is  because the first matrix contains only rubbish! 
	  int offset = matrixSize ;
	  unsigned int frame0Mask  = 0xFFF;
	  unsigned int frame0Shift = 0;
	  unsigned int frame1Mask  = 0xFFF000;
	  unsigned int frame1Shift = 12;

	  for ( int iDetector = 0; iDetector < _runHeader.VFasPresentNb + 1; iDetector++ ) {
	    	    
	    TrackerRawDataImpl * frame1 = new TrackerRawDataImpl;
	    TrackerRawDataImpl * frame0 = new TrackerRawDataImpl;
	    idEncoderFrame1["sensorID"] = iDetector;
	    idEncoderFrame1.setCellID(frame1);
	    idEncoderFrame0["sensorID"] = iDetector;
	    idEncoderFrame0.setCellID(frame0);
	    TrackerRawDataImpl * cds    = new TrackerRawDataImpl;
	    idEncoderCDS["sensorID"] = iDetector;
	    idEncoderCDS.setCellID(cds);
	    
	    int iPixel =  offset + iDetector * matrixSize;
	    for ( int y = 0; y < _noOfYPixel; y++ ) {  
	      for ( int x = 0; x < _noOfXPixel; x++ ) {
		short f0 = static_cast<short>(( _dataBuffer[iPixel] & frame0Mask ) >> frame0Shift);
		short f1 = static_cast<short>(( _dataBuffer[iPixel] & frame1Mask ) >> frame1Shift);
		frame0->adcValues().push_back( f0 );
		frame1->adcValues().push_back( f1 );
		cds->adcValues().push_back( f1 - f0 );
		++iPixel;
	      }
	    }
 	    vector<short > cdsVec = cds->adcValues();
	    vector<short >::iterator begin, end;

	    //	    cout << iDetector << " before min = " << (*min_element(cdsVec.begin(), cdsVec.end())) 
	    //		 << " max = " << (*max_element(cdsVec.begin(), cdsVec.end())) << endl;

 	    // correct for the CDS sign
	    if ( _eventHeader.VFasCnt < matrixSize ) {
	      begin = cdsVec.begin() + _eventHeader.VFasCnt ;
	      end   = cdsVec.end();
	    } else { 
	      begin = cdsVec.begin();
	      end   = cdsVec.begin() + _eventHeader.VFasCnt %  matrixSize;
	    }
	    transform( begin, end, begin, negate<int>() );
	    // 	    cout << iDetector << " after min = " << (*min_element(cdsVec.begin(), cdsVec.end())) 
	    //		 << " max = " << (*max_element(cdsVec.begin(), cdsVec.end())) << endl;

	    frame0Coll->push_back( frame0 ) ;
	    frame1Coll->push_back( frame1 ) ;
	    cdsColl->push_back( cds );
	      
	  }
	  
	  event->addCollection( frame0Coll, _frame0CollectionName );
	  event->addCollection( frame1Coll, _frame1CollectionName );
	  event->addCollection( cdsColl,    _cdsCollectionName    );
	  
	  ProcessorMgr::instance()->processEvent( event ) ;
	  delete event;
	  
	}  
	++eventCounter;
	if ( eventCounter > numEvents ) {
	  dataFile.close();
	  break;
	}
	
      }
    } catch (ifstream::failure e) {
      if ( !dataFile.eof() ) {
	message<ERROR5> ( log() << "Unable to open file " << dataFileName << ". Exiting." );
	addEORE();
	exit(-1);
      } else {
	dataFile.close();
	message<DEBUG5> ( log() << "Closing file " << dataFileName ) ;
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

  message<MESSAGE5> ("Successfully finished") ;
}

