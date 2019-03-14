#ifndef EUTelProcessorClusterAnalysis_h
#define EUTelProcessorClusterAnalysis_h 1

#include "lcio.h"
#include "marlin/Processor.h"
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

#include <IMPL/TrackerDataImpl.h>

#include "CrossSection.hpp"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile2D.h"
#include "cluster.h"

class EUTelProcessorClusterAnalysis : public marlin::Processor {
public:
  virtual Processor *newProcessor() {
    return new EUTelProcessorClusterAnalysis;
  }
  EUTelProcessorClusterAnalysis();
  virtual void init();
  virtual void processEvent(LCEvent *evt);
  int AddressToColumn(int ARegion, int ADoubleCol, int AAddress);
  int AddressToRow(int AAddress);
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  void bookHistos();
#endif

  virtual void end();
  bool emptyMiddle(std::vector<std::vector<int>> pixVector);

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
  int _chipVersion;
  int _numberofGeneratedInterestingCluster;
  int _numberofMissingInterestingCluster;
  int _number_emptyMiddle;
  int _cuttingSize;
  int _number_firing_event;
  std::vector<std::vector<int>> _before_event_memory;
  double _energy;
  EVENT::StringVec _chipID;
  EVENT::StringVec _irradiation;
  std::string _rate;
  bool _fillHistos;
  bool _clusterAvailable;
  bool _hotpixelAvailable;
  bool _noiseMaskAvailable;
  bool _deadColumnAvailable;
  bool _samecluster;
  bool _isDistanceSquareAnalysis;
  bool _isEmptyMiddleAnalysis;
  bool _isPlotSizeCutHitmap;
  bool _isPlotExampleEvents;
  bool _isDoubleFiringAnalysis;
  std::string _hotPixelCollectionName;
  std::string _deadColumnCollectionName;
  std::string _noiseMaskFileName;
  int _sectorSafetyPixels;
  std::vector<Cluster> clusterVec;
  std::map<int, int> xPairs;
  std::map<int, int> yPairs;
  std::vector<std::vector<int>> symmetryGroups;

private:
  bool _isFirstEvent;
  std::vector<int> noiseMaskX;
  std::vector<int> noiseMaskY;
  int _nLayer;
  int _xPixel;
  int _yPixel;
  int _sparseMinDistanceSquaredComparison;
  int _savedRandomEvents;
  bool _writeClustersToFile;
  std::map<int, TH1I *> clusterWidthXHisto;
  std::map<int, TH1I *> clusterWidthYHisto;
  std::map<int, TH1I *> clusterSizeHisto;
  std::map<int, TH2I *> GeneratedInterestingCluster;
  std::map<int, TH2I *> MissingInterestingCluster;
  std::map<int, TH2I *> emptyMiddleClusters;
  std::map<int, TH2I *> RandomEvent;
  std::map<int, TH2I *> Double_Firing_Events_Hitmap;
  TH1I *timeStampHisto;
  std::map<int, TH1I *> GeneratedClustersHisto;
  std::map<int, TH1I *> MissingClusterHisto;
  std::map<int, TH1I *> HowManyClusterGeneratedFromOneCluster;
  std::map<int, TH1I *> GeneratedClusterShapeHisto;
  std::map<int, TH1I *> MissingClusterShapeHisto;
  std::map<int, TH1I *> emptyMiddleClustersHisto;
  TH1I *NumberOfHits;
  TH1I *TypeOfTheEvent;
  TH1I *CLUSTER_SIZE;
  TH2I *hotpixelHisto;
  TH2I *deadColumnHisto;
  TH2I *circularClusterHistos;
  TH2I *largeClusterHistos;
  TH2I *smallerClustersHitmap;
  TH2I *biggerClustersHitmap;
  TH2I *doubleFiringPixels;
  TH2I *HIT_MAP;
  TH3I *clusterShapeMap;
  std::map<int, TH1I *> clusterShapeHistoSector;
  std::map<int, TH1I *> clusterShapeHistoGroupedSector;
  LCCollectionVec *hotPixelCollectionVec;
  LCCollectionVec *deadColumnCollectionVec;
  TrackerDataImpl *hotData;
  TrackerDataImpl *deadColumn;
};
#endif
