#ifndef AnalysisNoise_h
#define AnalysisNoise_h 1

#ifdef USE_GEAR

#include "marlin/Processor.h"
#include "lcio.h"
#include <algorithm>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include "EUTelSparseDataImpl.h"

#include "TH1.h"
#include "TH2.h"

#include <vector>
#include <fstream> 

using namespace std;
using namespace lcio ;
using namespace marlin;

class AnalysisNoise : public Processor { 
public: 
  virtual Processor* newProcessor() {return new AnalysisNoise;} 
  AnalysisNoise();
  virtual void init() ;
  virtual void processEvent( LCEvent * evt ) ;
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  void bookHistos();
#endif

  virtual void end();

protected:
  std::string _zsDataCollectionName;
  LCCollectionVec *zsInputDataCollectionVec;
  bool _fillHistos;
private:
  bool _isFirstEvent;
  int _nEvent;
  vector<vector<int> > _nFiredPixel;
  int _nLayer;
  vector<int> _xPixel;
  vector<int> _yPixel;
  gear::SiPlanesParameters * _siPlanesParameters;
  gear::SiPlanesLayerLayout * _siPlanesLayerLayout;
  double _energy;
  EVENT::StringVec _chipID;
  EVENT::StringVec _irradiation;
  EVENT::StringVec _dutIDs;
  string _rate;
  string _outputSettingsFolderName;
  ofstream settingsFile[3];
  std::map<int,TH2I*> noiseMap;
  std::map<int,TH1F*> noiseOccupancy;
  TH1I* timeStampHisto;
};
#endif
#endif
