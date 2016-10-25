#ifndef EUTelProcessorAnalysisPALPIDEfs_h
#define EUTelProcessorAnalysisPALPIDEfs_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <UTIL/CellIDDecoder.h>
#include <IMPL/LCCollectionVec.h>
//#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
//#include <AIDA/IBaseHistogram.h>
//#include <AIDA/IHistogram1D.h>
//#include <AIDA/IHistogram2D.h>
//#endif

#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>

#include <IMPL/TrackerDataImpl.h>

#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"
#include "cluster.h"
#include "CrossSection.hpp"

class EUTelProcessorAnalysisPALPIDEfs : public marlin::Processor {
public:
  virtual Processor* newProcessor() {return new EUTelProcessorAnalysisPALPIDEfs;}
  EUTelProcessorAnalysisPALPIDEfs();
  virtual void init() ;
  virtual void processEvent( LCEvent * evt ) ;
  void _EulerRotationBack(double* _telPos, double* _gRotation);
  void _LayerRotationBack(double* pos, double& outputX, double& outputY);
  int AddressToColumn(int ARegion, int ADoubleCol, int AAddress);
  int AddressToRow(int AAddress);
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  //! Book histograms
  /*! This method is used to prepare the needed directory structure
   *  within the current ITree folder and books all required
   *  histograms. Histogram pointers are stored into
   *  vectors for class-wide access
   */
  void bookHistos();
#endif
  virtual void end();
  bool emptyMiddle(std::vector<std::vector<int> > pixVector);
  bool RemoveAlign(LCCollectionVec * preAlignmentCollectionVec, LCCollectionVec * alignmentCollectionVec, LCCollectionVec * alignmentPAlpideCollectionVec, double* fitpos, double& xposfit, double& yposfit);
protected:
  //! Fill histogram switch
  /*! This boolean is used to switch on and off the filling of
   *  histograms.
   */
  bool _fillHistos;
  std::string _inputFittedHitName;
  std::string _inputColName;
  std::string _trackCollectionName;
  std::string _alignmentPAlpideCollectionName;
  std::string _alignmentCollectionName;
  std::string _preAlignmentCollectionName;
  std::string _zsDataCollectionName;
  //! The histogram information file
  /*! This string contain the name of the histogram information
   *  file. This is selected by the user in the steering file.
   *
   *  @see eutelescope::EUTelHistogramManager
   *  @see eutelescope::EUTelHistogramInfo
   */
  LCCollectionVec *zsInputDataCollectionVec;
  std::string _histoInfoFileName;
  std::string _hotPixelCollectionName;
  std::string _deadColumnCollectionName;
  std::string _noiseMaskFileName;
  double limit;
  int _dutID;
  int _maxNumberOfPixels;
  int _nPlanesWithMoreHits;
  bool _moreTracks;
  double _energy;
  bool _writeShapes;
  std::string _shapeOutputFileName;
  std::string _outputSettingsFolderName;
  std::ofstream settingsFile;
  EVENT::StringVec _chipID;
  EVENT::StringVec _irradiation;
  std::vector<float> _holesizeX;
  std::vector<float> _holesizeY;
  std::string _rate;
  bool _oneAlignmentCollection;
  bool _clusterAvailable;
  bool _hotpixelAvailable;
  bool _noiseMaskAvailable;
  bool _deadColumnAvailable;
  int layerIndex;
  double dutZ;
  FloatVec chi2Max;
  std::vector<Cluster> clusterVec;
  std::map<int,int> xPairs;
  std::map<int,int> yPairs;
  std::vector< std::vector<int> > symmetryGroups;
  double zDistance; // TODO: which units?
  int _nEvents;
  int _nEventsFake;
  int _nEventsWithTrack;
  int _nEventsWithEfficiency;
  double _minTimeStamp;
  int _nSectors;
  int _chipVersion;
  bool _showFake;
  bool _realAssociation;
private:
  bool _isFirstEvent;
  IntVec nTracks;
  IntVec nTracksFake;
  IntVec nTracksPAlpide;
  IntVec nTracksPAlpideFake;
  IntVec nTracksAssociation;
  IntVec nTracksPAlpideAssociation;
  std::vector<int> nFakeWithTrack;
  std::vector<int> nFakeWithoutTrack;
  std::vector<int> nFake;
  std::vector<int> nFakeWithTrackCorrected;
  int nDUThits;
  int nNoPAlpideHit;
  int nWrongPAlpideHit;
  int nPlanesWithTooManyHits;
  std::vector<int> noiseMaskX;
  std::vector<int> noiseMaskY;
  double xZero;
  double yZero;
  double xPitch;
  double yPitch;
  double xSize;
  double ySize;
  int xPixel;
  int yPixel;
  double gRotation[3];
  double xPointing[2];
  double yPointing[2];
  TH2I* tracksHisto;
  TH2* tracksPAlpideHisto;
  TH2* efficiencyHisto;
  TH2I* hitmapHisto;
  TH2I* hitmapNoHitHisto;
  TH2I* hitmapWrongHitHisto;
  TH2I* nFakeHitmapHisto;
  TH2I* nFakeWithTrackHitmapHisto;
  TH2I* nFakeWithoutTrackHitmapHisto;
  TH2I* nFakeWithTrackHitmapCorrectedHisto;
  TProfile2D* scatteringAngleHisto;
  TProfile2D* chi22DHisto;
  TH2I* tmpHist;
  TH1I* nHitsPerEventHisto;
  TH1I* nHitsPerEventHistoTime;
  std::map<float, std::map<int,TH1I*> > residualXPAlpide;
  std::map<float, std::map<int,TH1I*> > residualYPAlpide;
  std::map<float, std::map<int,TH1I*> > residualZPAlpide;
  std::map<float, std::map<int,TH1I*> > residualXPCBPAlpide;
  std::map<float, std::map<int,TH1I*> > residualYPCBPAlpide;
  std::map<float, std::map<int,TH1I*> > residualZPCBPAlpide;
  std::map<float, std::map<int,TH2F*> > residualXPixel;
  std::map<float, std::map<int,TH2F*> > residualYPixel;
  std::map<float, std::map<int,TH2F*> > residualXAveragePixel;
  std::map<float, std::map<int,TH2F*> > residualYAveragePixel;
  std::map<float, std::map<int,TH2F*> > nResidualXPixel;
  std::map<float, std::map<int,TH2F*> > nResidualYPixel;
  std::map<float, std::map<int,TH2F*> > residualXPixel2by2;
  std::map<float, std::map<int,TH2F*> > residualYPixel2by2;
  std::map<float, std::map<int,TH2F*> > residualXAveragePixel2by2;
  std::map<float, std::map<int,TH2F*> > residualYAveragePixel2by2;
  std::map<float, std::map<int,TH2F*> > nResidualXPixel2by2;
  std::map<float, std::map<int,TH2F*> > nResidualYPixel2by2;
  std::map<int,TH1I*> clusterShapeX;
  std::map<int,TH1I*> clusterShapeY;
  std::map<int,TH1I*> clusterShapeDiffX;
  std::map<int,TH1I*> clusterShapeDiffY;
  std::map<int,TH1I*> chi2X;
  std::map<int,TH1I*> chi2Y;
  TH1I* chi2Histo;
  TH1I* chi2HistoNoHit;
  std::map<int,TH1I*> clusterWidthXHisto;
  std::map<int,TH1I*> clusterWidthYHisto;
  std::map<int,TH1I*> clusterSizeHisto;
  TH1I* clusterShapeHisto;
  TH1I* clusterShapeHistoGrouped;
  std::map<int,TH1I*> clusterShapeHistoSector;
  std::map<int,TH1I*> clusterShapeHistoGroupedSector;
  std::map<int,TH1F*> clusterWidthXVsXHisto;
  std::map<int,TH1F*> clusterWidthYVsYHisto;
  std::map<int,TH1F*> clusterWidthXVsXAverageHisto;
  std::map<int,TH1F*> clusterWidthYVsYAverageHisto;
  std::map<int,TH2F*> clusterSize2DHisto;
  std::map<int,TH2F*> clusterSize2D2by2Histo;
  std::map<int,TH2F*> clusterSize2DAverageHisto;
  std::map<int,TH2F*> clusterSize2DAverage2by2Histo;
  std::map<int,TH1F*> nClusterVsXHisto;
  std::map<int,TH1F*> nClusterVsYHisto;
  std::map<int,TH2F*> nClusterSizeHisto;
  std::map<int,TH2F*> nClusterSize2by2Histo;
  std::map<int,TH2F*> efficiencyPixelHisto;
  std::map<int,TH2F*> efficiencyPixel2by2Histo;
  std::map<int,TH2F*> tracksPixelHisto;
  std::map<int,TH2F*> tracksPixel2by2Histo;
  std::map<int,TH2F*> tracksPAlpidePixelHisto;
  std::map<int,TH2F*> tracksPAlpidePixel2by2Histo;
  std::map<int, std::vector<TH1F*> > efficiencyPixelCrossSection;
  std::map<int, std::vector<TH1F*> > clusterSizeCrossSection;
  std::map<int,TH1D*> tracksProjection;
  std::map<int,TH1D*> tracksPAlpideProjection;
  std::map<int,TH1D*> efficiencyProjection;
  TH1I* timeStampHisto;
  TH1F* differenceX;
  TH1F* differenceY;
  TH1F* nFakeHisto;
  TH1F* nFakeWithTrackHisto;
  TH1F* nFakeWithTrackCorrectedHisto;
  TH1F* nFakeWithoutTrackHisto;
  TH2I* hotpixelHisto;
  TH2I* deadColumnHisto;
  TH2I* circularClusterHistos;
  TH2I* largeClusterHistos;
  LCCollectionVec * hotPixelCollectionVec;
  LCCollectionVec * deadColumnCollectionVec;
  TrackerDataImpl * hotData;
  TrackerDataImpl * deadColumn;
  std::map<int,TH2I*> clusterShape2D2by2;
  std::map<int,TH2I*> clusterShape2DGrouped2by2;
  TH1I* nTrackPerEventHisto;
  TH1I* nClusterAssociatedToTrackPerEventHisto;
  TH1I* nClusterPerEventHisto;
  TH1I* nAssociatedhitsHisto;
  TH1I* nAssociatedtracksHisto;
  TH1I* stats;
  std::vector< std::vector< std::vector<double> > > posFake;
  std::vector< std::vector< std::vector<double> > > posFit;
  std::vector< std::vector<double> > posFakeTemp;
  std::vector< int > nHitsPerEvent;

  enum statsEntry : int { kAll=0, kNoLargeClusters, kGoodTimeStamp, kDataAvailable,
      kFittedHitsAvailable, kTracksAvailable, kAlignAvailable, kOnlyOneAlignAvailable,
      kClustersAvailable, kSingularRotationMatrix, kFittedHit, kTwoCloseTracks,
      kTwoCloseTracksRejected, kZposCorrect, kZposWrong, kAlignmentRemoved,
      kAlignmentRemovalFailed, kFittedHitOnChip, kFittedHitNotOnChip, kFittedHitInBorderRegion,
      kUnknownSector, kHotPixel, kMaskedPixel, kDeadColumn, kHitsOnTheSamePlane, kSingleHitsOnly,
      kRejectedMultipleHits, kDenominator, kHitInDUT, kAssociatedHitInDUT, kNoHitInDUT,
      kHitNotInDUT, kHitInDUTNotAssociated, kHitNULL
      };
};

#endif
