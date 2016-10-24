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
#include "EUTelSucimaImagerReader.h"
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
#include <memory>
#include <cstdlib>

using namespace std;
using namespace marlin;
using namespace eutelescope;

EUTelSucimaImagerReader::EUTelSucimaImagerReader ():DataSourceProcessor  ("EUTelSucimaImagerReader") {
  
  _description =
    "Reads SUCIMA Imager ASCII data files and creates LCEvent with TrackerRawData collection.\n"
    "Make sure to not specify any LCIOInputFiles in the steering in order to read SUCIMA Imager files.";
  
  registerProcessorParameter ("SUCIMAImagerFileName", "Input file",
			      _fileName, std::string ("input.dat"));
  registerProcessorParameter ("NoOfXPixel", "Number of pixels along X",
			      _noOfXPixel, static_cast < int >(512));
  registerProcessorParameter ("NoOfYPixel", "Number of pixels along Y",
			      _noOfYPixel, static_cast < int >(512));
  
}

EUTelSucimaImagerReader * EUTelSucimaImagerReader::newProcessor () {
  return new EUTelSucimaImagerReader;
}

void EUTelSucimaImagerReader::init () {
  printParameters ();
}


void EUTelSucimaImagerReader::readDataSource (int numEvents) {

  ifstream inputFile;
  inputFile.exceptions (ifstream::failbit | ifstream::badbit);
  
  // try to open the input file....
  try {
    inputFile.open (_fileName.c_str (), ios::in);
  } catch (exception & e) {
    message<ERROR5> ( log() << "Problem opening file " << _fileName << ". Exiting." );
    exit (-1);
  }
  
  int eventNumber = 0;
  int runNumber = 0;
  
  
  while (true)  {
    
    if (numEvents > 0 && eventNumber + 1 > numEvents) {
      break;
    }
    
    EUTelEventImpl *event = new EUTelEventImpl;
    event->setDetectorName("MIMOSA");
    event->setEventType(kDE);
    LCTime * now = new LCTime;
    event->setTimeStamp(now->timeStamp());
    delete now;

    LCCollectionVec *rawData = new LCCollectionVec (LCIO::TRACKERRAWDATA);
    
    if (eventNumber % 10 == 0)  message<MESSAGE5> ( log() << "Converting event " << eventNumber );

    
    if (isFirstEvent ()) {
      
      // in the case it is the first run so we need to process and
      // write out the run header.
      auto lcHeader = std::make_unique<IMPL::LCRunHeaderImpl>();
      auto runHeader = std::make_unique<EUTelRunHeaderImpl>(lcHeader.get());
      runHeader->addProcessor( type() );
      runHeader->lcRunHeader()->setDescription(" Events read from SUCIMA Imager ASCII input file: " + _fileName);
      runHeader->lcRunHeader()->setRunNumber (runNumber);
      runHeader->setHeaderVersion (0.0001);
      runHeader->setDataType (EUTELESCOPE::CONVDATA);
      runHeader->setDateTime ();
      runHeader->setDAQHWName (EUTELESCOPE::SUCIMAIMAGER);
      runHeader->setDAQHWVersion (0.0001);
      runHeader->setDAQSWName (EUTELESCOPE::SUCIMAIMAGER);
      runHeader->setDAQSWVersion (0.0001);
      runHeader->addIntermediateFile (_fileName);
      runHeader->addProcessor (_processorName);
      // this is a mistake here only for testing....
      runHeader->setNoOfEvent(100);
      //////////////////////////////////////////////
      runHeader->setNoOfDetector(1);
      runHeader->setMinX(IntVec(1, 0));
      runHeader->setMaxX(IntVec(1, _noOfXPixel - 1));
      runHeader->setMinY(IntVec(1, 0));
      runHeader->setMaxY(IntVec(1, _noOfYPixel - 1));
      runHeader->lcRunHeader()->setDetectorName("MIMOSA");
      // UTIL::LCTOOLS::dumpRunHeader(runHeader);
      
      // process the run header
      ProcessorMgr::instance ()->processRunHeader ( static_cast<lcio::LCRunHeader*> ( lcHeader.release()) );
      
      // prepare the short buffet
      _buffer = new short[_noOfXPixel * _noOfYPixel];
      
      // end of first event
      _isFirstEvent = false;
    }
    
    int nVal = 0;
    int xPixel, yPixel;
    try  {
      for (yPixel = 0; yPixel < _noOfYPixel; yPixel++) {
	for (xPixel = 0; xPixel < _noOfXPixel; xPixel++) {
	  inputFile >> _buffer[nVal++];
	}
      }
    } catch (exception & e) {
      if (!inputFile.eof ())
	cerr << "A read exception occurred : " << e.what () << endl;
      if ((xPixel == 0) && (yPixel == 0)) {
	// that's normal
	// we are reading the last empty line.
	// break here
	break;
      } else {
	// that's strange. it might be that the last event was not complete
	// break here
	message<ERROR5> ( log() << "Event " << eventNumber << " finished un-expectedly. " );
	message<ERROR5> ( log() << "Consider to check the input file, or limit the conversion to " 	    
			 << eventNumber - 1 << " events" << endl << "Sorry to quit!" ) ;
	exit(-1);
      }
    }
    
    event->setRunNumber (runNumber);
    event->setEventNumber (eventNumber++);
    
    // prepare a collection to store the raw-data
    
    TrackerRawDataImpl *rawMatrix = new TrackerRawDataImpl;
    CellIDEncoder < TrackerRawDataImpl > idEncoder (EUTELESCOPE::MATRIXDEFAULTENCODING, rawData);
    idEncoder["sensorID"] = 0;
    idEncoder["xMin"] = 0;
    idEncoder["xMax"] = _noOfXPixel - 1;
    idEncoder["yMin"] = 0;
    idEncoder["yMax"] = _noOfYPixel - 1;
    idEncoder.setCellID (rawMatrix);
    
    nVal = 0;
    for (int yPixel = 0; yPixel < _noOfYPixel; yPixel++) {
      for (int xPixel = 0; xPixel < _noOfXPixel; xPixel++) {
	rawMatrix->adcValues ().push_back (_buffer[nVal++]);
      }
    }
    rawData->push_back (rawMatrix);
    
    event->addCollection (rawData, "rawdata");
    ProcessorMgr::instance ()->processEvent (static_cast<LCEventImpl*> (event));
    delete event;
  }

  EUTelEventImpl *event = new EUTelEventImpl;
  event->setDetectorName("MIMOSA");
  LCTime * now = new LCTime;
  event->setTimeStamp(now->timeStamp());
  delete now;
  event->setRunNumber (runNumber);
  event->setEventNumber (eventNumber++);
  event->setEventType(kEORE);
  ProcessorMgr::instance ()->processEvent (static_cast<LCEventImpl*> (event));
  delete event;
  
}



void EUTelSucimaImagerReader::end () {
  delete[]_buffer;
  message<MESSAGE5> ("Successfully finished") ;
}
