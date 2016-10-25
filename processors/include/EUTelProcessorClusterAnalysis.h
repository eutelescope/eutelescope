#ifndef EUTelProcessorClusterAnalysis_h
#define EUTelProcessorClusterAnalysis_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <UTIL/CellIDDecoder.h>
#include <IMPL/LCCollectionVec.h>

#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>

#include <IMPL/TrackerDataImpl.h>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile2D.h"
#include "cluster.h"
#include "CrossSection.hpp"

class EUTelProcessorClusterAnalysis : public marlin::Processor {
public:
  virtual Processor* newProcessor() {return new EUTelProcessorClusterAnalysis;}
  EUTelProcessorClusterAnalysis();
  virtual void init() ;
  virtual void processEvent( LCEvent * evt ) ;
  int AddressToColumn(int ARegion, int ADoubleCol, int AAddress);
  int AddressToRow(int AAddress);
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  void bookHistos();
#endif

  virtual void end();
  bool emptyMiddle(std::vector<std::vector<int> > pixVector);

protected:
  std::string _zsDataCollectionName;
  LCCollectionVec *zsInputDataCollectionVec;
  std::string _histoInfoFileName; 
  std::string _clusterAnalysisFile;
  std::ofstream clusterAnalysisOutput;
  std::string _outputSettingsFolderName;  
  std::ofstream settingsFile;
  int _nEvents;
  int _dutID;
  int _maxNumberOfPixels;
  int _layerIndex;
  int _nSectors;
  int _nTouchingBorderSectorClusters;
  int _nTouchingBorderYClusters;
  int _nOverlappingClusters;
  int _nHotpixelClusters;
  int _nNoiseMaskClusters;
  int _nDeadColumnClusters;
  int _sectorWidth;
  double _energy;
  EVENT::StringVec _chipID;
  EVENT::StringVec _irradiation;
  std::string _rate;
  bool _fillHistos;
  bool _clusterAvailable;
  bool _hotpixelAvailable;
  bool _noiseMaskAvailable;
  bool _deadColumnAvailable;
  std::string _hotPixelCollectionName;
  std::string _deadColumnCollectionName;
  std::string _noiseMaskFileName;
  int _sectorSafetyPixels;
  std::vector<Cluster> clusterVec;
  std::map<int,int> xPairs;
  std::map<int,int> yPairs;
  std::vector< std::vector<int> > symmetryGroups;
private:
  bool _isFirstEvent;
  std::vector<int> noiseMaskX;
  std::vector<int> noiseMaskY;
  int _nLayer;
  int _xPixel;
  int _yPixel;
  std::map<int,TH1I*> clusterWidthXHisto;
  std::map<int,TH1I*> clusterWidthYHisto;
  std::map<int,TH1I*> clusterSizeHisto;
  TH1I* timeStampHisto;
  TH2I* hotpixelHisto;
  TH2I* deadColumnHisto;
  TH2I* circularClusterHistos;
  TH2I* largeClusterHistos;
  TH3I* clusterShapeMap;
  std::map<int,TH1I*> clusterShapeHistoSector;
  std::map<int,TH1I*> clusterShapeHistoGroupedSector;
  LCCollectionVec * hotPixelCollectionVec;
  LCCollectionVec * deadColumnCollectionVec;
  TrackerDataImpl * hotData;
  TrackerDataImpl * deadColumn;
};
#endif
