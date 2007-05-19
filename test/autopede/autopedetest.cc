// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: autopedetest.cc,v 1.2 2007-05-19 09:55:25 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#include "lcio.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "IMPL/LCEventImpl.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/TrackerRawDataImpl.h"
#include "UTIL/CellIDEncoder.h"
#include "UTIL/LCTime.h"

#include <vector>
#include <fstream>
#include <cmath>

using namespace std;
using namespace lcio;
using namespace eutelescope;

const int nDetector  = 5;
const int nDataEvent = 100;

const int maxClusterPerEvent = 10;

const int   xCluSize  = 5;
const int   yCluSize  = 5;
const short clusterSignal[xCluSize * yCluSize] =  {0, 1, 2, 1, 0,
						   1, 5, 8, 6, 2,
						   3, 9, 13,7, 2,
						   2, 5, 8, 4, 1,
						   0, 1, 2, 1, 0};


const int xNPixel = 512;
const int yNPixel = 512;


void getXYFromIndex(int index, int& x, int& y);
int  getIndexFromXY(int x, int y);
void usage();

ofstream logfile;

int main(int argc, char ** argv) {
  
  float basePede = 10., baseNoise = 2.;


  if (argc == 1) {
    usage();
    return 0;
  } else {
  
    if ( argc > 1 ) basePede  = atof(argv[1]);
    if ( argc > 2 ) baseNoise = atof(argv[2]);

  }

  vector<int> minX(nDetector,0);
  vector<int> minY(nDetector,0);
  vector<int> maxX(nDetector,xNPixel - 1);
  vector<int> maxY(nDetector,yNPixel - 1);


  LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();
  try {
    lcWriter->open("data_input.slcio", LCIO::WRITE_NEW);
  } catch (IOException& e) {
      cerr << e.what() << endl;
      return 0;
    }

  EUTelRunHeaderImpl * runHeader = new EUTelRunHeaderImpl(); 
  runHeader->setRunNumber(0);
  runHeader->setDetectorName("test");
  runHeader->setHeaderVersion(0.0011);
  runHeader->setDataType(EUTELESCOPE::CONVDATA);
  runHeader->setDateTime();
  runHeader->setDAQHWName(EUTELESCOPE::SUCIMAIMAGER);
  runHeader->setDAQHWVersion(0.0001);
  runHeader->setDAQSWName(EUTELESCOPE::SUCIMAIMAGER);
  runHeader->setDAQSWVersion(0.0001);  
  runHeader->setNoOfEvent(nDataEvent);
  runHeader->setNoOfDetector(nDetector);
  runHeader->setMinX(minX);
  runHeader->setMaxX(maxX);
  runHeader->setMinY(minY);
  runHeader->setMaxY(maxY);
  
  lcWriter->writeRunHeader(runHeader);
  delete runHeader;
  
  short * matrix = new short[xNPixel * yNPixel];
  
  logfile.open("clustering.log");
  
  for (int iEvent = 0; iEvent <= nDataEvent; iEvent++) {
    if ( iEvent % 10 == 0) 
      cout << "Data on event " << iEvent << endl;
    
    logfile << "Event " << iEvent << endl;
    EUTelEventImpl * event = new EUTelEventImpl;
    event->setEventNumber(iEvent);
    event->setDetectorName("test");
    
    LCTime * now = new LCTime;
    event->setTimeStamp(now->timeStamp());
    delete now;
    

    if ( iEvent < nDataEvent) {
      event->setEventType(kDE);
      
      LCCollectionVec * rawData = new LCCollectionVec(LCIO::TRACKERRAWDATA);
      
      for (int iDetector = 0; iDetector < nDetector; iDetector++) {
	logfile << "Working on detector " << iDetector << endl;
	
	int  baseSignal = static_cast<int>( basePede  ); 
	
	TrackerRawDataImpl * rawMatrix = new TrackerRawDataImpl;
	CellIDEncoder<TrackerRawDataImpl> idEncoder("sensorID:5,xMin:12,xMax:12,yMin:12,yMax:12", rawData);
	idEncoder["sensorID"] = iDetector;
	idEncoder["xMin"]     = 0;
	idEncoder["xMax"]     = xNPixel - 1;
	idEncoder["yMin"]     = 0;
	idEncoder["yMax"]     = yNPixel - 1;
	idEncoder.setCellID(rawMatrix);
	
	int noiseDivider = static_cast<int> (sqrt(12) * baseNoise);
	
	int iPixel = 0;
	for (int yPixel = 0; yPixel < yNPixel; yPixel++) {
	  for (int xPixel = 0; xPixel < xNPixel; xPixel++) {
	    matrix[iPixel] = baseSignal +  (short) ((rand()/(RAND_MAX / noiseDivider)) - baseNoise) ;
	    ++iPixel;
	  } 
	}
	
	// set the number of cluster 
	int clusterPerEvent = rand() / (RAND_MAX / maxClusterPerEvent);
	logfile << "  injected " << clusterPerEvent << " clusters " << endl; 
	for (int iSeed = 0; iSeed < clusterPerEvent; iSeed++) {
	  int index = rand() / (RAND_MAX / (xNPixel * yNPixel));
	  int xSeed, ySeed;
	  getXYFromIndex(index, xSeed, ySeed);
	  logfile << "iSeed " << iSeed << " xSeed " << xSeed << " ySeed " << ySeed << endl;
	  int iCluPos = 0;
	  for (int yPixel = ySeed - (yCluSize / 2); yPixel <=  ySeed + (yCluSize / 2); yPixel++) {
	    for (int xPixel = xSeed - (xCluSize / 2); xPixel <=  xSeed + (xCluSize / 2); xPixel++) {
	      if ( ( xPixel >= 0 ) && ( xPixel < xNPixel ) &&
		   ( yPixel >= 0 ) && ( yPixel < yNPixel ) ) {
		index = getIndexFromXY(xPixel, yPixel);
		matrix[index] += clusterSignal[iCluPos];
		//		logfile << "x = " << xPixel << " y = " << yPixel << " s = " << matrix[index] << endl;
	      }
	      ++iCluPos;
	    }
	  }
	}
	
	iPixel = 0;
	for (int yPixel = 0; yPixel < yNPixel; yPixel++) {
	  for (int xPixel = 0; xPixel < xNPixel; xPixel++) {
	    rawMatrix->adcValues().push_back(matrix[iPixel]);
	    ++iPixel;
	  }
	}
	rawData->push_back(rawMatrix);
      }
      event->addCollection(rawData,"rawdata");
    } else event->setEventType(kEORE);

    LCEventImpl * lcevent = static_cast<LCEventImpl*> (event) ;
    lcWriter->writeEvent(lcevent);

    delete event;
  }
  
  logfile.close();
  lcWriter->close();
  delete [] matrix;
  return 0;
}

void getXYFromIndex(int index, int& x, int& y) {

  y = (index / xNPixel);
  x = index - (y * xNPixel);
}

int getIndexFromXY(int x, int y) {
  return x + y * xNPixel;
}

void usage() {
  cout << "./autopedetest basepede basenoise " << endl;

}
