// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelSucimaImagerReader.cc,v 1.4 2007-02-11 08:46:01 bulgheroni Exp $
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

// marlin includes
#include "marlin/DataSourceProcessor.h"
#include "marlin/ProcessorMgr.h"

// lcio includes
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <UTIL/CellIDEncoder.h>
// #include <UTIL/LCTOOLS.h>

// system includes 
#include <fstream>


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
    cerr << "Problem opening file " << _fileName << ". Exiting." << endl;
    exit (-1);
  }
  
  int eventNumber = 0;
  int runNumber = 0;
  
  
  while (true)  {
    
    if (numEvents > 0 && eventNumber + 1 > numEvents) {
      break;
    }
    
    LCEventImpl *event = new LCEventImpl;
    LCCollectionVec *rawData = new LCCollectionVec (LCIO::TRACKERRAWDATA);
    
    if (eventNumber % 10 == 0)
      cout << "[" << name() << "] Converting event " << eventNumber << endl;
    
    if (isFirstEvent ()) {
      
      // in the case it is the first run so we need to process and
      // write out the run header.
      EUTelRunHeaderImpl *runHeader = new EUTelRunHeaderImpl;
      runHeader->setDescription(" Events read from SUCIMA Imager ASCII input file: " + _fileName);
      runHeader->setRunNumber (runNumber);
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
      runHeader->setNoOfEvent( 99);
      //////////////////////////////////////////////
      runHeader->setNoOfDetector(1);
      runHeader->setMinX(IntVec(1, 0));
      runHeader->setMaxX(IntVec(1, _noOfXPixel - 1));
      runHeader->setMinY(IntVec(1, 0));
      runHeader->setMaxY(IntVec(1, _noOfYPixel - 1));

      // UTIL::LCTOOLS::dumpRunHeader(runHeader);
      
      // process the run header
      ProcessorMgr::instance ()->processRunHeader (runHeader);
      
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
	cerr << "A read exception occured : " << e.what () << endl;
      if ((xPixel == 0) && (yPixel == 0)) {
	// that's normal
	// we are reading the last empty line.
	// break here
	break;
      } else {
	// that's strange. it might be that the last event was not complete
	// break here
	cerr << "Event " << eventNumber << " finished un-expectedly. "  << endl 
	     << "Consider to check the input file, or limit the converion to "
	     << eventNumber - 1 << " events" << endl << "Sorry to quit!" << endl;
	break;
      }
    }
    
    event->setRunNumber (runNumber);
    event->setEventNumber (eventNumber++);
    
    // prepare a collection to store the raw-data
    
    TrackerRawDataImpl *rawMatrix = new TrackerRawDataImpl;
    CellIDEncoder < TrackerRawDataImpl > idEncoder ("sensorID:5,xMin:12,xMax:12,yMin:12,yMax:12", rawData);
    idEncoder["sensorID"] = 0;
    idEncoder["xMin"] = 0;
    idEncoder["xMax"] = _noOfXPixel;
    idEncoder["yMin"] = 0;
    idEncoder["yMax"] = _noOfYPixel;
    idEncoder.setCellID (rawMatrix);
    
    nVal = 0;
    for (int yPixel = 0; yPixel < _noOfYPixel; yPixel++) {
      for (int xPixel = 0; xPixel < _noOfXPixel; xPixel++) {
	rawMatrix->adcValues ().push_back (_buffer[nVal++]);
      }
    }
    rawData->push_back (rawMatrix);
    
    event->addCollection (rawData, "rawdata");
    ProcessorMgr::instance ()->processEvent (event);
    delete event;
  }
}



void EUTelSucimaImagerReader::end () {
  delete[]_buffer;
}
