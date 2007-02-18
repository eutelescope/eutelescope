#include "lcio.h"
#include "EUTelRunHeaderImpl.h"
#include "IMPL/LCEventImpl.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/TrackerRawDataImpl.h"
#include "UTIL/CellIDEncoder.h"

#include <TList.h>
#include <TH1.h>
#include <TFile.h>
#include <TProfile2D.h>
#include <iostream>

using namespace lcio;
using namespace eutelescope;
using namespace std;

int main(int argc , char ** argv) {


  const int nEvent   = 1000;
  const int nXPixel  = 10;
  const int nYPixel  = 10;

  short * matrix     = new short[nXPixel * nYPixel];
  TH1D  **histoArray = new TH1D*[nXPixel * nYPixel];
  TList * histoList  = new TList();
  
  TFile * outputFile   = TFile::Open("file.root","RECREATE");
  TProfile2D * profile = new TProfile2D("profile","profile", nXPixel, -0.5, nXPixel - 0.5,
					nYPixel, -0.5,  nYPixel - 0.5);
  
  LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();
  try {
    lcWriter->open("pedetest.slcio",LCIO::WRITE_NEW);
  } catch (IOException& e) {
    cerr << e.what() << endl;
    return 0;
  }
  
  EUTelRunHeaderImpl * runHeader = new EUTelRunHeaderImpl();
  runHeader->setRunNumber(0);
  runHeader->setDetectorName("pedetest");
  runHeader->setHeaderVersion(0.0011);
  runHeader->setDataType(EUTELESCOPE::CONVDATA);
  runHeader->setDateTime();
  runHeader->setDAQHWName(EUTELESCOPE::SUCIMAIMAGER);
  runHeader->setDAQHWVersion(0.0001);
  runHeader->setDAQSWName(EUTELESCOPE::SUCIMAIMAGER);
  runHeader->setDAQSWVersion(0.0001);  
  runHeader->setNoOfEvent(nEvent);
  runHeader->setNoOfDetector(1);
  vector<int> minX; minX.push_back(0);
  vector<int> maxX; maxX.push_back(nXPixel - 1);
  vector<int> minY; minY.push_back(0);
  vector<int> maxY; maxY.push_back(nYPixel - 1);
  runHeader->setMinX(minX);
  runHeader->setMaxX(maxX);
  runHeader->setMinY(minY);
  runHeader->setMaxY(maxY);
  
  
  lcWriter->writeRunHeader(runHeader);
  


  for (int iEvent = 0; iEvent < nEvent; iEvent++) {

    LCEventImpl * event = new LCEventImpl;
    event->setEventNumber(iEvent);
    LCCollectionVec * rawData = new LCCollectionVec(LCIO::TRACKERRAWDATA);

    int    iPixel     = 0;
    short  baseSignal = 1;
    float  noise      = 3;
    
    TrackerRawDataImpl * rawMatrix = new TrackerRawDataImpl;
    CellIDEncoder<TrackerRawDataImpl> idEncoder("sensorID:5,xMin:12,xMax:12,yMin:12,yMax:12", rawData);
    idEncoder["sensorID"] = 0;
    idEncoder["xMin"]     = 0;
    idEncoder["xMax"]     = nXPixel - 1;
    idEncoder["yMin"]     = 0;
    idEncoder["yMax"]     = nYPixel - 1;
    idEncoder.setCellID(rawMatrix);
    
    for (int yPixel = 0; yPixel < nYPixel; yPixel++) {
      for (int xPixel = 0; xPixel < nXPixel; xPixel++) {
	matrix[iPixel] = baseSignal + (short) (rand()/(RAND_MAX / noise));
	profile->Fill(xPixel, yPixel, matrix[iPixel]);
	rawMatrix->adcValues().push_back(matrix[iPixel]);

	if (iEvent == 0 ) {
	  stringstream ss;
	  ss << "h-" << xPixel << "-" << yPixel << "-" << (noise/sqrt(12));
	  histoArray[iPixel] = new TH1D(ss.str().c_str(), ss.str().c_str(), 100, 0., 0.);
	  histoArray[iPixel]->SetBit(TH1::kCanRebin);
	  histoList->Add(histoArray[iPixel]);
	}

	histoArray[iPixel]->Fill(matrix[iPixel]); 
	++iPixel; 
	baseSignal += 5;
	noise += 3;
      }
    }
    rawData->push_back(rawMatrix);
    event->addCollection(rawData,"rawdata");
    lcWriter->writeEvent(event);
    delete event;
  }
  
  

  profile->Write();
  histoList->Write();
  lcWriter->close();
  outputFile->Close();

  delete [] histoArray;
  delete [] matrix;
  
  return 0;
}
