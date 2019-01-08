#include "EUTelProcessorClusterAnalysis.h"
#include "EUTelAlignmentConstant.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelHistogramManager.h"
#include "EUTelProcessorAnalysisPALPIDEfs.h"
#include "EUTelTrackerDataInterfacerImpl.h"

#include "marlin/AIDAProcessor.h"
#include "marlin/Global.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IAxis.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/ITree.h>
#endif

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackerHitImpl.h>

#include "TFitResult.h"
#include "TVector3.h"

#include <algorithm>
#include <cmath>
#include <memory>

using namespace lcio;
using namespace marlin;
using namespace std;
using namespace eutelescope;
using namespace gear;

EUTelProcessorClusterAnalysis aEUTelProcessorClusterAnalysis;

EUTelProcessorClusterAnalysis::EUTelProcessorClusterAnalysis()
    : Processor("EUTelProcessorClusterAnalysis"), _zsDataCollectionName(""),
      _clusterAnalysisFile(""), _outputSettingsFolderName("./"), _nEvents(0),
      _dutID(), _maxNumberOfPixels(3), _layerIndex(-1), _nSectors(),
      _nTouchingBorderSectorClusters(0), _nTouchingBorderYClusters(0),
      _nOverlappingClusters(0), _nHotpixelClusters(0), _nNoiseMaskClusters(0),
      _nDeadColumnClusters(0), _sectorWidth(), _chipVersion(4),
      _numberofGeneratedInterestingCluster(0),
      _numberofMissingInterestingCluster(0), _number_emptyMiddle(0),
      _cuttingSize(5), _number_firing_event(0), _energy(6.0), _chipID(),
      _irradiation(), _rate(""), _fillHistos(false),
      _isDistanceSquareAnalysis(false),
      _isEmptyMiddleAnalysis(false), _isPlotSizeCutHitmap(false),
      _isPlotExampleEvents(false), _isDoubleFiringAnalysis(false),
      _hotPixelCollectionName(""), _nLayer(0), _xPixel(), _yPixel(),
      _sparseMinDistanceSquaredComparison(1),
      _savedRandomEvents(0),
      _writeClustersToFile(false)

{
  _description = "Analysing cluster properties such as cluster shape and "
                 "average cluster size.";
  registerInputCollection(LCIO::TRACKERDATA, "ZSDataCollectionName",
                          "Input of Zero Suppressed data",
                          _zsDataCollectionName, string("zsdata"));
  registerProcessorParameter("HistogramFilling",
                             "Switch on or off the histogram filling",
                             _fillHistos, true);
  registerProcessorParameter(
      "HistoInfoFileName", "This is the name of the histogram information file",
      _histoInfoFileName, string("histoinfo.xml"));
  registerOptionalParameter(
      "ClusterAnalysisFileName", "This is the name of file to which all the "
                                 "information on the pixels belonging to the "
                                 "clusters is saved to.",
      _clusterAnalysisFile, static_cast<string>("clusterAnalysis.txt"));
  registerOptionalParameter(
      "HotPixelCollectionName",
      "This is the name of the hotpixel collection of the pALPIDE",
      _hotPixelCollectionName, static_cast<string>(""));
  registerOptionalParameter("DeadColumnCollectionName",
                            "This is the name of the collection containing the "
                            "pixels belonging to a dead column",
                            _deadColumnCollectionName,
                            static_cast<string>("deadColumn"));
  registerOptionalParameter("NoiseMaskFileName",
                            "This is the name of the file which contains the "
                            "pixels which were masked during datataking",
                            _noiseMaskFileName, static_cast<string>(""));
  registerOptionalParameter(
      "OutputSettingsFolderName",
      "Folder name where all the settings of each run will be saved",
      _outputSettingsFolderName, static_cast<string>("./"));
  registerProcessorParameter("dutID", "This is the ID of the DUT", _dutID,
                             6);
  registerProcessorParameter("MaxNumberOfPixels",
                             "This is the maximum number of pixels in one "
                             "cluster for the clustershape analysis",
                             _maxNumberOfPixels, 3);
  registerProcessorParameter("nSectors",
                             "This is the maximum amount of sectors", _nSectors,
                             8);
  registerOptionalParameter("SectorSafetyPixels",
                            "Safety distance (in pixel) of clusters being "
                            "associated to a sector and to the boundaries of "
                            "the chip.",
                            _sectorSafetyPixels, 2);
  registerOptionalParameter("Energy", "Particle energy [GeV]", _energy,
                            6.0);
  EVENT::StringVec _stringVecExample;
  _stringVecExample.push_back(" ");
  registerOptionalParameter("ChipID", "Chip IDs", _chipID, _stringVecExample);
  registerOptionalParameter("Irradiation", "Irradiation level", _irradiation,
                            _stringVecExample);
  registerOptionalParameter("Rate", "Data taking rate", _rate,
                            static_cast<string>(""));
  registerOptionalParameter("ChipVersion", "Chip Version", _chipVersion,
                            4);
  registerOptionalParameter(
      "IsDistanceSquareAnalysis", "Bool to analyze the clusters with a smaller "
                                  "distance allowed between the pixels",
      _isDistanceSquareAnalysis, false);
  registerOptionalParameter("IsEmptyMiddleAnalysis",
                            "Bool to analyze the clusters with empty middle",
                            _isEmptyMiddleAnalysis, false);
  registerOptionalParameter(
      "IsPlotSizeCutHitmap",
      "Bool to create hit maps with a cut on the cluster size",
      _isPlotSizeCutHitmap, false);
  registerOptionalParameter("IsPlotExampleEvents",
                            "Bool to plot the first 100 cluster examples",
                            _isPlotExampleEvents, false);
  registerOptionalParameter(
      "IsDoubleFiringAnalysis",
      "Bool to do the analysis of clusters appearing in two consecutive events",
      _isDoubleFiringAnalysis, false);
  registerOptionalParameter("SparseMinDistanceSquaredComparison",
                            "Sparse Min Distance Squared Comparison",
                            _sparseMinDistanceSquaredComparison,
                            1);
  // This cluster size is the border between the small and the big clusters:
  registerOptionalParameter(
      "CuttingSize",
      "This cluster size is the cut between the small and the big clusters",
      _cuttingSize, 5);
  registerOptionalParameter("WriteClustersToFile",
                            "Bool to write clusters to text file",
                            _writeClustersToFile, false);
  _isFirstEvent = true;
}

void EUTelProcessorClusterAnalysis::init() {
  _nLayer = geo::gGeometry().nPlanes();
  const std::vector<int> &_planeID = geo::gGeometry().sensorIDsVec();

  for (int iz = 0; iz < _nLayer; iz++)
    if (_planeID[iz] == _dutID)
      _layerIndex = iz;
  if (_layerIndex == -1) {
    cerr << "Wrong DUT ID. Exiting." << endl;
    return;
  }

  if (_chipVersion < 3)
    _nSectors = 4;
  else if (_chipVersion == 3)
    _nSectors = 8;
  else if (_chipVersion == 5)
    _nSectors = 4;
  else
    _nSectors = 1;

  // beware, sometimes dutID is 3, sometimes it is 6
  int iLayer = _dutID;
  _xPixel = geo::gGeometry().getPlaneNumberOfPixelsX(iLayer);
  _yPixel = geo::gGeometry().getPlaneNumberOfPixelsY(iLayer);
  _sectorWidth = _xPixel / _nSectors;
  //
  Cluster cluster;
  cluster.FindReferenceClusters(clusterVec, _maxNumberOfPixels);
  xPairs = cluster.SymmetryPairs(clusterVec, "x");
  yPairs = cluster.SymmetryPairs(clusterVec, "y");
  symmetryGroups = cluster.sameShape(clusterVec);
  //
  // Open Files
  // Analysis Folder
  if (_writeClustersToFile)
    clusterAnalysisOutput.open(_clusterAnalysisFile);
  // Settings folder
  bool newFile = false;
  string _outputSettingsFileName =
      _outputSettingsFolderName + Form("settings_DUT%d", _dutID) + ".txt";
  if (!std::ifstream(_outputSettingsFileName.c_str()))
    newFile = true;
  settingsFile.open(_outputSettingsFileName.c_str(), ios::out | ios::app);
  if (newFile)
    settingsFile << "Run number;Energy;Chip ID;Chip Version;Irradiation "
                    "level(0-nonIrradiated,1-2.5e12,2-1e13,3-700krad,4-"
                    "combined:1e13+700krad);Rate;BB;Ithr;Idb;Vcasn;Vcasn2;"
                    "Vclip;Vcasp;VresetP;VresetD;Threshold and their RMS for "
                    "all eight sectors;Noise and their RMS for all eight "
                    "sectors;Readout delay;Trigger delay;Strobe length;StrobeB "
                    "length;Data (1) or noise (0);Number of "
                    "events;Efficiency,Number of tracks,Number of tracks with "
                    "associated hit for all sectors"
                 << endl;
}

void EUTelProcessorClusterAnalysis::processEvent(LCEvent *evt) {
  // INIT, DEAD COLOUMN AND HOT PIXEL CHECKS
  // -----------------------------------------------------------------------------------------------------------------
  //
  //
  if (_layerIndex == -1) {
    return;
  }
  int nClusterPerEvent = 0;
  if (_isFirstEvent) {
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    if (_fillHistos) {
      bookHistos();
    }
#endif
    hotPixelCollectionVec = nullptr;
    try {
      hotPixelCollectionVec = static_cast<LCCollectionVec *>(
          evt->getCollection(_hotPixelCollectionName));
      streamlog_out(DEBUG5)
          << "hotPixelCollectionName: " << _hotPixelCollectionName.c_str()
          << " found " << endl;
    } catch (lcio::DataNotAvailableException &e) {
      streamlog_out(WARNING5)
          << "hotPixelCollectionName: " << _hotPixelCollectionName.c_str()
          << " not found " << endl;
      _hotpixelAvailable = false;
    }
    ifstream noiseMaskFile(_noiseMaskFileName.c_str());
    if (noiseMaskFile.is_open()) {
      streamlog_out(MESSAGE4)
          << "Running with noise mask: " << _noiseMaskFileName.c_str() << endl;
      int region, doubleColumn, address;
      while (noiseMaskFile >> region >> doubleColumn >> address) {
        int x = AddressToColumn(region, doubleColumn, address);
        int y = AddressToRow(address);
        noiseMaskX.push_back(x);
        noiseMaskY.push_back(y);
      }
    } else
      _noiseMaskAvailable = false;

    deadColumnCollectionVec = nullptr;
    try {
      deadColumnCollectionVec = static_cast<LCCollectionVec *>(
          evt->getCollection(_deadColumnCollectionName));
    } catch (lcio::DataNotAvailableException &e) {
      streamlog_out(WARNING5)
          << "deadPixelCollectionName: " << _deadColumnCollectionName.c_str()
          << " not found " << endl;
      _deadColumnAvailable = false;
    }
    // pALPIDE 3 settings
    settingsFile
        << evt->getRunNumber() << ";" << _energy << ";" << _chipID[_layerIndex]
        << ";3;" << _irradiation[_layerIndex] << ";" << _rate << ";"
        << evt->getParameters().getFloatVal("BackBiasVoltage") << ";"
        << evt->getParameters().getIntVal(Form("Ithr_%d", _layerIndex)) << ";"
        << evt->getParameters().getIntVal(Form("Idb_%d", _layerIndex)) << ";"
        << evt->getParameters().getIntVal(Form("Vcasn_%d", _layerIndex)) << ";"
        << evt->getParameters().getIntVal(Form("Vcasn2_%d", _layerIndex)) << ";"
        << evt->getParameters().getIntVal(Form("Vclip_%d", _layerIndex)) << ";"
        << evt->getParameters().getIntVal(Form("Vcasp_%d", _layerIndex)) << ";"
        << evt->getParameters().getIntVal(Form("VresetP_%d", _layerIndex))
        << ";"
        << evt->getParameters().getIntVal(Form("VresetD_%d", _layerIndex))
        << ";";

    for (int iSector = 0; iSector < _nSectors; iSector++)
      settingsFile << evt->getParameters().getFloatVal(
                          Form("Thr_%d_%d", _layerIndex, iSector))
                   << ";"
                   << evt->getParameters().getFloatVal(
                          Form("ThrRMS_%d_%d", _layerIndex, iSector))
                   << ";";
    for (int iSector = 0; iSector < _nSectors; iSector++)
      settingsFile << evt->getParameters().getFloatVal(
                          Form("Noise_%d_%d", _layerIndex, iSector))
                   << ";"
                   << evt->getParameters().getFloatVal(
                          Form("NoiseRMS_%d_%d", _layerIndex, iSector))
                   << ";";
    settingsFile << evt->getParameters().getIntVal(
                        Form("m_readout_delay_%d", _layerIndex))
                 << ";"
                 << evt->getParameters().getIntVal(
                        Form("m_trigger_delay_%d", _layerIndex))
                 << ";"
                 << evt->getParameters().getIntVal(
                        Form("m_strobe_length_%d", _layerIndex))
                 << ";"
                 << evt->getParameters().getIntVal(
                        Form("m_strobeb_length_%d", _layerIndex))
                 << ";1;";
    _isFirstEvent = false;
  }
  _clusterAvailable = true;
  try {
    zsInputDataCollectionVec = dynamic_cast<LCCollectionVec *>(
        evt->getCollection(_zsDataCollectionName));
    streamlog_out(DEBUG5) << "zsInputDataCollectionVec: "
                          << _zsDataCollectionName.c_str() << " found " << endl;
  } catch (lcio::DataNotAvailableException& ) {
    streamlog_out(DEBUG5) << "zsInputDataCollectionVec: "
                          << _zsDataCollectionName.c_str() << " not found "
                          << endl;
    _clusterAvailable = false;
  }
  if (_clusterAvailable) {
    CellIDDecoder<TrackerDataImpl> cellDecoder(zsInputDataCollectionVec);
    for (size_t iCluster = 0; iCluster < zsInputDataCollectionVec->size();
         iCluster++) {
      TrackerDataImpl *zsData = dynamic_cast<TrackerDataImpl *>(
          zsInputDataCollectionVec->getElementAt(iCluster));
      if (cellDecoder(zsData)["sensorID"] == _dutID)
        nClusterPerEvent++;
    }
  }

  // Write DEBUG OUTPUTS in this way//streamlog_out ( DEBUG )  <<
  // nAssociatedhits << " points for one track in DUT in event " <<
  // evt->getEventNumber() << "\t" << xposPrev << "\t" << yposPrev << "\t" <<
  // xpos << "\t" << ypos << " Fit: " << xposfit << "\t" << yposfit << " Number
  // of planes with more than one hit: " << nPlanesWithMoreHits << endl;

  //
  //
  // END INIT, DEAD COLUMN AND HOTPIXELS
  // -------------------------------------------------------------------------------------------------------------------------

  if (_clusterAvailable) {
    int numberOfHitsInAnEvent = 0;
    int numberOfSmallClusters = 0;
    int numberOfBigClusters = 0;
    std::vector<std::vector<int>> event_memory;
    for (size_t idetector = 0; idetector < zsInputDataCollectionVec->size();
         idetector++) {
      CellIDDecoder<TrackerDataImpl> cellDecoder(zsInputDataCollectionVec);
      TrackerDataImpl *zsData = dynamic_cast<TrackerDataImpl *>(
          zsInputDataCollectionVec->getElementAt(idetector));
      SparsePixelType type = static_cast<SparsePixelType>(
          static_cast<int>(cellDecoder(zsData)["sparsePixelType"]));

      // Check whether the data is the one from the DUT or not
      if (cellDecoder(zsData)["sensorID"] == _dutID) {
        int clusterSize = zsData->getChargeValues().size() / 4;
        vector<int> X(clusterSize);
        vector<int> Y(clusterSize);
        Cluster cluster;
        if (type == kEUTelGenericSparsePixel) {
          // starting actual cluster analysis
          vector<vector<int>> pixVector;
          auto sparseData =
              EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>(zsData);

          for (size_t iPixel = 0; iPixel < sparseData.size(); iPixel++) {
            auto &pixel = sparseData.at(iPixel);
            X[iPixel] = pixel.getXCoord();
            Y[iPixel] = pixel.getYCoord();

            // Check, whether a cluster is at the boundary of two sectors
            for (int safetyPixels = 0; safetyPixels <= _sectorSafetyPixels;
                 safetyPixels++) {
              // the pixels may not be on the boudary of a sector, in particular
              // not for the outer most pixels
              if ((X[iPixel] + safetyPixels % _sectorWidth == 0) ||
                  (X[iPixel] - safetyPixels % _sectorWidth == 0)) {
                _nTouchingBorderSectorClusters++;
                streamlog_out(DEBUG5)
                    << "Another cluster which touches the borders of a sector."
                    << endl;
                // if it is, discard the cluster
                goto nextCluster;
              }
            }
            // check, whether pixel touches the y outside
            if (Y[iPixel] < _sectorSafetyPixels ||
                Y[iPixel] >= (_yPixel - _sectorSafetyPixels)) {
              _nTouchingBorderYClusters++;
              streamlog_out(DEBUG5)
                  << "Another cluster which touches the y-borders of the chip."
                  << endl;
              // if it is, discard the cluster
              goto nextCluster;
            }

            // check, whether all pixels are from the same sector. If not,
            // discard the cluster and take note of it.
            if (X[iPixel] / _sectorWidth != X[0] / _sectorWidth) {
              _nOverlappingClusters++;
              streamlog_out(DEBUG5) << "Another cluster which overlaps sectors."
                                    << endl;
              // Cluster discarded since it overlapps sectors; if there is only
              // one sector, this can be used as a debug output, nextCluster
              // should always be 0 then
              goto nextCluster;
            }

            // check, if cluster includes a hotpixel, noiseMask or deadcolumn
            //
            //
            if (_hotpixelAvailable) {
              auto sparseData =
                  EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>(
                      hotData);
              for (size_t iHotPixel = 0; iHotPixel < sparseData.size();
                   iHotPixel++) {
                auto &sparsePixel = sparseData.at(iHotPixel);
                if (abs(X[iPixel] - (sparsePixel.getXCoord())) <
                        _sectorSafetyPixels &&
                    abs(Y[iPixel] - (sparsePixel.getYCoord())) <
                        _sectorSafetyPixels) {
                  _nHotpixelClusters++;
                  goto nextCluster;
                }
              }
            }
            if (_noiseMaskAvailable) {
              for (size_t iNoise = 0; iNoise < noiseMaskX.size(); iNoise++) {
                if (abs(X[iPixel] - (noiseMaskX[iNoise])) <
                        _sectorSafetyPixels &&
                    abs(Y[iPixel] - (noiseMaskY[iNoise])) <
                        _sectorSafetyPixels) {
                  _nNoiseMaskClusters++;
                  goto nextCluster;
                }
              }
            }
            if (_deadColumnAvailable) {
              auto sparseData =
                  EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>(
                      deadColumn);
              for (size_t iDeadPixel = 0; iDeadPixel < sparseData.size();
                   iDeadPixel++) {
                auto &sparsePixel = sparseData.at(iDeadPixel);
                if (abs(X[iPixel] - (sparsePixel.getXCoord())) <
                    _sectorSafetyPixels) {
                  _nDeadColumnClusters++;
                  goto nextCluster;
                }
              }
            }
            //
            //
            // End of Hotpixel, NoiseMask and DeadColumn search

            vector<int> pix;
            pix.push_back(X[iPixel]);
            pix.push_back(Y[iPixel]);
            pixVector.push_back(pix);
          }

          streamlog_out(DEBUG5) << "This is a DEBUG output to see whether the "
                                   "program gets here. The number X[0] is "
                                << X[0] << " and _sectorWidth is "
                                << _sectorWidth << endl;
          // now, since all pixels are from the same sector, the sector number
          // can be set.
          int index = X[0] / _sectorWidth;

          // This part is to analysis the effect of the distance square between
          // the pixels in one cluste
          if (_isDistanceSquareAnalysis) {
            _samecluster = true;
            int howmanyclustergeneratedfromonecluster(0);
            int AllGeneratedPixel(0);
            int AllMissingPixel(0);

            std::vector<EUTelGenericSparsePixel> hitPixelVec =
                sparseData.getPixels();
            std::vector<EUTelGenericSparsePixel> hitPixelVec2 =
                sparseData.getPixels();

            std::vector<EUTelGenericSparsePixel> newlyAdded;

            unsigned int firstclustersize = hitPixelVec.size();
            // We now cluster those hits together
            while (!hitPixelVec.empty()) {
              std::vector<EUTelGenericSparsePixel> cluCandidate;

              // First we need to take any pixel, so let's take the first one
              // Add it to the cluster as well as the newly added pixels
              newlyAdded.push_back(hitPixelVec.front());
              // sparseCluster->push_back( &(hitPixelVec.front()) );
              cluCandidate.push_back(hitPixelVec.front());
              // And remove it from the original collection
              hitPixelVec.erase(hitPixelVec.begin());

              // Now process all newly added pixels, initially this is the just
              // previously added one
              // but in the process of neighbour finding we continue to add new
              // pixels
              while (!newlyAdded.empty()) {
                bool newlyDone = true;
                int x1, x2, y1, y2, dX, dY;

                // check against all pixels in the hitPixelVec
                for (std::vector<EUTelGenericSparsePixel>::iterator hitVec =
                         hitPixelVec.begin();
                     hitVec != hitPixelVec.end(); ++hitVec) {
                  // get the relevant infos from the newly added pixel
                  x1 = newlyAdded.front().getXCoord();
                  y1 = newlyAdded.front().getYCoord();

                  // and the pixel we test against
                  x2 = hitVec->getXCoord();
                  y2 = hitVec->getYCoord();

                  dX = x1 - x2;
                  dY = y1 - y2;
                  int distance = dX * dX + dY * dY;
                  // if they pass the spatial and temporal cuts, we add them

                  if (distance <= _sparseMinDistanceSquaredComparison) {
                    // add them to the cluster as well as to the newly added
                    // ones
                    newlyAdded.push_back(*hitVec);
                    cluCandidate.push_back(*hitVec);
                    //	sparseCluster->push_back( &(*hitVec) );
                    hitPixelVec.erase(hitVec);
                    // for the pixel we test there might be other neighbours, we
                    // still have to check
                    newlyDone = false;
                    break;
                  }
                }

                // if no neighbours are found, we can delete the pixel from the
                // newly added
                // we tested against _ALL_ non cluster pixels, there are no
                // other pixels
                // which could be neighbours
                if (newlyDone) {
                  newlyAdded.erase(newlyAdded.begin());
                }
              }

              if (firstclustersize != cluCandidate.size()) {
                _samecluster = false;
                howmanyclustergeneratedfromonecluster++;
                AllGeneratedPixel += cluCandidate.size();
                GeneratedClustersHisto[index]->Fill(cluCandidate.size());
                // cout<<"I filled GeneratedClustersHisto with:
                // "<<cluCandidate.size()<<endl;

                int intrestingClusterSize = cluCandidate.size();
                Cluster interestingCluster;
                vector<int> X(intrestingClusterSize);
                vector<int> Y(intrestingClusterSize);

                int iforX = 0, Xmax = 0, Ymax = 0, Xmin = 1000000,
                    Ymin = 1000000, Xshift = 0, Yshift = 0;
                while (!cluCandidate.empty()) {
                  X[iforX] = cluCandidate.front().getXCoord();
                  Y[iforX] = cluCandidate.front().getYCoord();
                  cluCandidate.erase(cluCandidate.begin());
                  if (X[iforX] < Xmin)
                    Xmin = X[iforX];
                  if (Y[iforX] < Ymin)
                    Ymin = Y[iforX];
                  if (X[iforX] > Xmax)
                    Xmax = X[iforX];
                  if (Y[iforX] > Ymax)
                    Ymax = Y[iforX];

                  iforX++;
                }

                interestingCluster.set_values(intrestingClusterSize, X, Y);
                GeneratedClusterShapeHisto[index]->Fill(
                    interestingCluster.WhichClusterShape(interestingCluster,
                                                         clusterVec));

                Xshift = (Xmax + Xmin) / 2 - 50 / 2;
                Yshift = (Ymax + Ymin) / 2 - 50 / 2;
                for (size_t iforY = 0; iforY < Y.size() &&
                                    _numberofGeneratedInterestingCluster < 100;
                     iforY++) {
                  GeneratedInterestingCluster
                      [_numberofGeneratedInterestingCluster]
                          ->Fill(X[iforY] - Xshift, Y[iforY] - Yshift);
                }
                _numberofGeneratedInterestingCluster++;
              }
            }

            if (!_samecluster) {
              MissingClusterHisto[index]->Fill(firstclustersize);
              HowManyClusterGeneratedFromOneCluster[index]->Fill(
                  howmanyclustergeneratedfromonecluster);
              howmanyclustergeneratedfromonecluster = 0;
              AllMissingPixel = firstclustersize;

              // cout<<"I filled MissingClusterHisto with:
              // "<<firstclustersize<<endl;

              int Xmax = 0, Ymax = 0, Xmin = 1000000, Ymin = 1000000,
                  Xshift = 0, Yshift = 0;

              for (std::vector<EUTelGenericSparsePixel>::iterator
                       hitVecSparseData = hitPixelVec2.begin();
                   hitVecSparseData != hitPixelVec2.end() &&
                   _numberofMissingInterestingCluster < 100;
                   ++hitVecSparseData) {
                if (hitVecSparseData->getXCoord() < Xmin)
                  Xmin = hitVecSparseData->getXCoord();
                if (hitVecSparseData->getYCoord() < Ymin)
                  Ymin = hitVecSparseData->getYCoord();
                if (hitVecSparseData->getXCoord() > Xmax)
                  Xmax = hitVecSparseData->getXCoord();
                if (hitVecSparseData->getYCoord() > Ymax)
                  Ymax = hitVecSparseData->getYCoord();
              }
              Xshift = (Xmax + Xmin) / 2 - 50 / 2;
              Yshift = (Ymax + Ymin) / 2 - 50 / 2;

              for (std::vector<EUTelGenericSparsePixel>::iterator
                       hitVecSparseData = hitPixelVec2.begin();
                   hitVecSparseData != hitPixelVec2.end() &&
                   _numberofMissingInterestingCluster < 100;
                   ++hitVecSparseData) {
                MissingInterestingCluster[_numberofMissingInterestingCluster]
                    ->Fill(hitVecSparseData->getXCoord() - Xshift,
                           hitVecSparseData->getYCoord() - Yshift);
              }
              _numberofMissingInterestingCluster++;
            }

            if (AllGeneratedPixel != AllMissingPixel) {
              cerr << "AllMissingPixel!=AllMissingPixel" << endl;
              cerr << "AllMissingPixel: " << AllMissingPixel << endl;
              cerr << "AllGeneratedPixel: " << AllGeneratedPixel << endl;
            }
            AllMissingPixel = 0;
            AllGeneratedPixel = 0;
          }

          // The end of the part for distance analysis

          // This par looking for empty middle clusters

          if (_isEmptyMiddleAnalysis) {

            // This line is needed, because EUTelProcessorAnalysisPALPIDEfs will
            // check if the cluster is empty middle, or not.
            EUTelProcessorAnalysisPALPIDEfs *mypalpide =
                new EUTelProcessorAnalysisPALPIDEfs();

            // The next line check, if the cluster empty middled
            if (mypalpide->emptyMiddle(pixVector)) {
              // It fill, the holey clusters histo
              emptyMiddleClustersHisto[index]->Fill(pixVector.size());
              // It select holey clusters, to see them.
              int xMin = *min_element(X.begin(), X.end());
              int xMax = *max_element(X.begin(), X.end());
              int yMin = *min_element(Y.begin(), Y.end());
              int yMax = *max_element(Y.begin(), Y.end());
              int Xshift = (xMin + xMax) / 2 - 50 / 2;
              int Yshift = (yMin + yMax) / 2 - 50 / 2;
              for (size_t i_emptyMiddle = 0; i_emptyMiddle < pixVector.size() &&
                                          _number_emptyMiddle < 100;
                   i_emptyMiddle++) {
                emptyMiddleClusters[_number_emptyMiddle]->Fill(
                    pixVector[i_emptyMiddle][0] - Xshift,
                    pixVector[i_emptyMiddle][1] - Yshift);
              }
              _number_emptyMiddle++;
            }
            delete mypalpide;
          }

          // The end of the empty middle clusters part

          // This plots the clusters withe a cut in size

          if (_isPlotSizeCutHitmap) {
            if (pixVector.size() < static_cast<unsigned int> (_cuttingSize))
              numberOfSmallClusters++;
            if (pixVector.size() >= static_cast<unsigned int> (_cuttingSize))
              numberOfBigClusters++;
            for (size_t i_cuttingSize = 0; pixVector.size() < static_cast<unsigned int> (_cuttingSize) &&
                                        i_cuttingSize < pixVector.size();
                 i_cuttingSize++) {
              smallerClustersHitmap->Fill(pixVector[i_cuttingSize][0],
                                          pixVector[i_cuttingSize][1]);
            }

            for (size_t i_cuttingSize = 0; pixVector.size() >= static_cast<unsigned int> (_cuttingSize) &&
                                        i_cuttingSize < pixVector.size();
                 i_cuttingSize++) {
              biggerClustersHitmap->Fill(pixVector[i_cuttingSize][0],
                                         pixVector[i_cuttingSize][1]);
            }
          }

          // This part is for the hitmap

          for (size_t i_HITMAP = 0; i_HITMAP < pixVector.size(); i_HITMAP++) {
            HIT_MAP->Fill(pixVector[i_HITMAP][0], pixVector[i_HITMAP][1]);
          }

          // This part plots random events to see.

          if (_isPlotExampleEvents) {
            for (size_t i_random = 0;
                 i_random < pixVector.size() && _savedRandomEvents < 100;
                 i_random++) {
              RandomEvent[_savedRandomEvents]->Fill(pixVector[i_random][0],
                                                    pixVector[i_random][1]);
            }
          }

          // This part is to analysis if a pixel fires, will it fire with higher
          // probability the next event?

          if (_isDoubleFiringAnalysis) {
            for (size_t i_event_memory = 0; i_event_memory < pixVector.size();
                 i_event_memory++) {
              event_memory.push_back(pixVector[i_event_memory]);
            }
          }

          // This line is to check how many pixel fired in an event.
          numberOfHitsInAnEvent += pixVector.size();

          // set the cluster
          cluster.set_values(clusterSize, X, Y);

          // fill the file with the pixel information of a cluster, after it has
          // been checked, whether the cluster is ok
          if (_writeClustersToFile) {
            for (size_t iPixel = 0; iPixel < sparseData.size(); iPixel++) {
              // Output Pixel information to file, Comma (,) is coordinate
              // seperator, Whitespace ( ) is Pixel seperator
              if (_writeClustersToFile)
                clusterAnalysisOutput << X[iPixel] << "," << Y[iPixel] << " ";
            }
            // the next cluster will be the next thing written to the file, so
            // include the cluster seperator	semicolon
            clusterAnalysisOutput << ";";
          }
          // start the plotting
          //
          //
          clusterSizeHisto[index]->Fill(clusterSize);
          int xMin = *min_element(X.begin(), X.end());
          int xMax = *max_element(X.begin(), X.end());
          int yMin = *min_element(Y.begin(), Y.end());
          int yMax = *max_element(Y.begin(), Y.end());
          int clusterWidthX = xMax - xMin + 1;
          int clusterWidthY = yMax - yMin + 1;

          clusterWidthXHisto[index]->Fill(clusterWidthX);
          clusterWidthYHisto[index]->Fill(clusterWidthY);

          int clusterShape = cluster.WhichClusterShape(cluster, clusterVec);
          if (clusterShape >= 0) {
            clusterShapeHistoSector[index]->Fill(clusterShape);
          } else
            clusterShapeHistoSector[index]->Fill(clusterVec.size());
        }
      }
    nextCluster:;
      // End cluster for loop
    }

    // How many pixel fired in one event
    NumberOfHits->Fill(numberOfHitsInAnEvent);

    // What type of event
    if (_isPlotSizeCutHitmap) {
      if (numberOfSmallClusters == 0 && numberOfBigClusters == 0)
        TypeOfTheEvent->Fill(0);
      if (numberOfSmallClusters > 0 && numberOfBigClusters == 0)
        TypeOfTheEvent->Fill(1);
      if (numberOfSmallClusters > 0 && numberOfBigClusters == 1)
        TypeOfTheEvent->Fill(2);
      if (numberOfSmallClusters > 0 && numberOfBigClusters > 1)
        TypeOfTheEvent->Fill(3);
      if (numberOfSmallClusters == 0 && numberOfBigClusters > 0)
        TypeOfTheEvent->Fill(4);
      // cerr<<"BigClusters: "<<numberOfBigClusters<<", SmallClusters:
      // "<<numberOfSmallClusters<<endl;
    }

    // It is for generating some hitmapa for an example event
    if (_isPlotExampleEvents) {
      _savedRandomEvents++;
    }

    // This part is to analyze if a pixel fires, will it fire with higher
    // probability the next event?
    if (_isDoubleFiringAnalysis) {
      bool interestin_event = false;
      for (size_t i_event_memory = 0; i_event_memory < event_memory.size();
           i_event_memory++) {
        bool there_is_a_double_firing = false;
        for (size_t j_event_memory = 0;
             j_event_memory < _before_event_memory.size(); j_event_memory++) {
          if (event_memory[i_event_memory][0] ==
                  _before_event_memory[j_event_memory][0] &&
              event_memory[i_event_memory][1] ==
                  _before_event_memory[j_event_memory][1]) {
            there_is_a_double_firing = true;
            interestin_event = true;
            break;
          }
        }
        if (there_is_a_double_firing) {
          doubleFiringPixels->Fill(event_memory[i_event_memory][0],
                                   event_memory[i_event_memory][1]);
        }
      }
      // In the next step, it is going to create some hitmaps from the double
      // firing events. If you see a pixel in the hitmap, which has:
      // 1 entry, it just fired in the first event.
      // 2 entries, it just fired in the second event.
      // 3 entries, it fired in both of the events.

      if (interestin_event) {

        for (size_t i_event_memory = 0;
             i_event_memory < event_memory.size() && _number_firing_event < 100;
             i_event_memory++) {
          Double_Firing_Events_Hitmap[_number_firing_event]->Fill(
              event_memory[i_event_memory][0], event_memory[i_event_memory][1]);
          Double_Firing_Events_Hitmap[_number_firing_event]->Fill(
              event_memory[i_event_memory][0], event_memory[i_event_memory][1]);
        }
        for (size_t i_event_memory = 0;
             i_event_memory < _before_event_memory.size() &&
             _number_firing_event < 100;
             i_event_memory++) {
          Double_Firing_Events_Hitmap[_number_firing_event]->Fill(
              _before_event_memory[i_event_memory][0],
              _before_event_memory[i_event_memory][1]);
        }
        _number_firing_event++;
      }
      _before_event_memory = event_memory;
    }
  }

  // write the end event expression to the file, which is a linebreak
  if (_writeClustersToFile)
    clusterAnalysisOutput << endl;
  // Increment number of events
  _nEvents++;
}

void EUTelProcessorClusterAnalysis::bookHistos() {
  streamlog_out(DEBUG1) << "Booking histograms " << endl;
  auto histoMgr = std::make_unique<EUTelHistogramManager>(_histoInfoFileName);

  try {
    histoMgr->init();
  } catch (ios::failure &e) {
    streamlog_out(WARNING2) << "I/O problem with " << _histoInfoFileName << "\n"
                            << "Continuing without histogram manager" << endl;
  } catch (marlin::ParseException &e) {
    streamlog_out(WARNING2) << e.what() << "\n"
                            << "Continuing without histogram manager" << endl;
  }
  AIDAProcessor::tree(this)->cd("Analysis");
  clusterShapeMap = new TH3I(
      "ClusterShapeMap", "ClusterShapeMap;Pixel X;Pixel Y;Cluster Shape ID",
      _maxNumberOfPixels + 1, -0.5, _maxNumberOfPixels + 0.5,
      _maxNumberOfPixels + 1, -0.5, _maxNumberOfPixels + 0.5,
      clusterVec.size() + 1, -0.5, clusterVec.size() + 0.5);
  for (int iSector = 0; iSector < _nSectors; iSector++) {
    AIDAProcessor::tree(this)->mkdir(Form("Sector_%d", iSector));
    AIDAProcessor::tree(this)->cd(Form("Sector_%d", iSector));

    clusterWidthXHisto[iSector] = new TH1I(
        Form("clusterWidthXHisto_%d", iSector),
        Form("Cluster width in X in sector %d;Cluster width X (pixel);a.u.",
             iSector),
        50, 0.5, 50.5);
    clusterWidthYHisto[iSector] = new TH1I(
        Form("clusterWidthYHisto_%d", iSector),
        Form("Cluster width in Y in sector %d;Cluster width Y (pixel);a.u.",
             iSector),
        50, 0.5, 50.5);
    clusterSizeHisto[iSector] =
        new TH1I(Form("clusterSizeHisto_%d", iSector),
                 Form("Cluster size_%d;Cluster size (pixel);a.u.", iSector),
                 200, 0.5, 200.5);
    clusterShapeHistoSector[iSector] =
        new TH1I(Form("clusterShapeHisto_%d", iSector),
                 Form("Cluster shape (all rotations separately) Sector "
                      "%d;Cluster shape ID;a.u.",
                      iSector),
                 clusterVec.size() + 1, -0.5, clusterVec.size() + 0.5);
    clusterShapeHistoGroupedSector[iSector] =
        new TH1I(Form("clusterShapeHistoGrouped_%d", iSector),
                 Form("Cluster shape (all rotations treated together) Sector "
                      "%d;Cluster shape ID;a.u.",
                      iSector),
                 symmetryGroups.size(), -0.5, symmetryGroups.size() - 0.5);
    if (_isDistanceSquareAnalysis) {
      AIDAProcessor::tree(this)->mkdir(Form(
          "Sector_%d/CompareHistogramsIfReduced_sparseMinDistanceSquared=%d",
          iSector, _sparseMinDistanceSquaredComparison));
      AIDAProcessor::tree(this)->cd(Form(
          "Sector_%d/CompareHistogramsIfReduced_sparseMinDistanceSquared=%d",
          iSector, _sparseMinDistanceSquaredComparison));
      GeneratedClustersHisto[iSector] = new TH1I(
          Form("GeneratedClustersFromAnother_%d", iSector),
          Form("Generated clusters, if reduced _sparseMinDistanceSquared=%d, "
               "Sector %d;Cluster size (pixel);Number of clusters",
               _sparseMinDistanceSquaredComparison, iSector),
          200, 0.5, 200.5);
      MissingClusterHisto[iSector] =
          new TH1I(Form("DisintegratingClusters_%d", iSector),
                   Form("Disintegrating clusters, if reduced "
                        "_sparseMinDistanceSquared=%d, Sector %d;Cluster size "
                        "(pixel);Number of clusters",
                        _sparseMinDistanceSquaredComparison, iSector),
                   200, 0.5, 200.5);
      HowManyClusterGeneratedFromOneCluster[iSector] =
          new TH1I(Form("HowManyClusterGeneratedFromOneCluster_%d", iSector),
                   Form("How many cluster generated from one cluster if reduce "
                        "_sparseMinDistanceSquared=%d, Sector %d;Number of "
                        "generated clusters from one disintegrating cluster; "
                        "number of missing clusters",
                        _sparseMinDistanceSquaredComparison, iSector),
                   20, 0.5, 20.5);
      GeneratedClusterShapeHisto[iSector] =
          new TH1I(Form("GeneratedClustersShape_%d", iSector),
                   Form("These clusters generated, if reduced "
                        "_sparseMinDistanceSquared=%d, Sector %d;Cluster size "
                        "(pixel);Number of clusters",
                        _sparseMinDistanceSquaredComparison, iSector),
                   clusterVec.size() + 1, -0.5, clusterVec.size() + 0.5);
      // The next line will use, if the Cluster Shape finder process can work
      // with clusters, which have pixels touch with corner.
      // MissingClusterShapeHisto[iSector] = new
      // TH1I(Form("MissingClusterShapeHisto_%d",iSector),Form("MissingClusterShapeHisto;Type
      // of the cluster shape;Number of
      // clusters"),clusterVec.size()+1,-0.5,clusterVec.size()+0.5);
    }
    AIDAProcessor::tree(this)->cd(Form("Sector_%d", iSector));
    if (_isEmptyMiddleAnalysis)
      emptyMiddleClustersHisto[iSector] =
          new TH1I(Form("emptyMiddleClustersHisto_%d", iSector),
                   Form("Empty middle clusters histo in sector %d;Cluster size "
                        "(pixel);Number of Clusters",
                        iSector),
                   200, 0.5, 200.5);
  }
  AIDAProcessor::tree(this)->cd("");
  if (_isPlotSizeCutHitmap)
    smallerClustersHitmap =
        new TH2I(Form("smallerClustersHitmap"),
                 Form("This hitmap is filled with clusters, which smaller than "
                      "%d;X (pixel);Y (pixel)",
                      _cuttingSize),
                 1024, 0, 1024, 512, 0, 512);
  if (_isPlotSizeCutHitmap)
    biggerClustersHitmap =
        new TH2I(Form("biggerClustersHitmap"),
                 Form("This hitmap is filled with clusters, which equal or "
                      "bigger than %d;X (pixel);Y (pixel)",
                      _cuttingSize),
                 1024, 0, 1024, 512, 0, 512);
  NumberOfHits =
      new TH1I(Form("NumberOfHits"),
               Form("Number of hits in an event;n_Hits;Number of events"), 500,
               0.5, 500.5);
  // 0: There was no hit.
  // 1: There was just smaller clusters than _cuttingSize.
  // 2: There was one cluster bigger or equal than _cuttingSize and there was
  // smaller clusters than _cuttingSize.
  // 3: There was more than one cluster bigger or equal than _cuttingSize and
  // there was smaller clusters than _cuttingSize.
  // 4: There was just bigger or equal clusters than _cuttingSize.
  if (_isPlotSizeCutHitmap)
    TypeOfTheEvent =
        new TH1I(Form("TypeOfTheEvent"),
                 Form("Type of the event;Type;Number of events"), 5, -0.5, 4.5);
  // if(_isPlotSizeCutHitmap) TypeOfTheEvent->SetMarkerStyle(21);
  if (_isDoubleFiringAnalysis)
    doubleFiringPixels =
        new TH2I(Form("doubleFiringPixels"),
                 Form("Double firing pixels;X (pixel);Y (pixel)"), 1024, 0,
                 1024, 512, 0, 512);
  HIT_MAP = new TH2I("hitMap", Form("HitMap;X (pixel); Y(pixel)"), 1024, 0,
                     1024, 512, 0, 512);

  AIDAProcessor::tree(this)->mkdir(Form("First100Example"));
  for (int nInterestingCluster = 0; nInterestingCluster < 100;
       nInterestingCluster++) {
    if (_isDistanceSquareAnalysis) {
      AIDAProcessor::tree(this)->mkdir(
          Form("First100Example/GeneratedClustersFromADisintegratingCluster"));
      AIDAProcessor::tree(this)->cd(
          Form("First100Example/GeneratedClustersFromADisintegratingCluster"));
      GeneratedInterestingCluster[nInterestingCluster] =
          new TH2I(Form("GeneratedClustersFromADisintegratingCluster_%d",
                        nInterestingCluster),
                   Form(" Generated cluster, example %d;Cluster width X "
                        "(pixel);Cluster width Y (pixel)",
                        nInterestingCluster),
                   50, 0, 50, 50, 0, 50);
      AIDAProcessor::tree(this)->mkdir(
          Form("First100Example/DisintegratingCluster"));
      AIDAProcessor::tree(this)->cd(
          Form("First100Example/DisintegratingCluster"));
      MissingInterestingCluster[nInterestingCluster] =
          new TH2I(Form("MissingInterestingCluster_%d", nInterestingCluster),
                   Form(" Missing cluster, example %d;Cluster width X "
                        "(pixel);Cluster width Y (pixel)",
                        nInterestingCluster),
                   50, 0, 50, 50, 0, 50);
    }
    if (_isEmptyMiddleAnalysis) {
      AIDAProcessor::tree(this)->mkdir(
          Form("First100Example/emptyMiddleClusters"));
      AIDAProcessor::tree(this)->cd(
          Form("First100Example/emptyMiddleClusters"));
      emptyMiddleClusters[nInterestingCluster] =
          new TH2I(Form("emptyMiddleClusters_%d", nInterestingCluster),
                   Form(" Holey cluster, example %d;Cluster width X "
                        "(pixel);Cluster width Y (pixel)",
                        nInterestingCluster),
                   50, 0, 50, 50, 0, 50);
    }
    if (_isPlotExampleEvents) {
      AIDAProcessor::tree(this)->mkdir(Form("First100Example/RandomEvent"));
      AIDAProcessor::tree(this)->cd(Form("First100Example/RandomEvent"));
      RandomEvent[nInterestingCluster] =
          new TH2I(Form("RandomEvent_%d", nInterestingCluster),
                   Form(" An Event for Demonstration %d;Cluster width X "
                        "(pixel);Cluster width Y (pixel)",
                        nInterestingCluster),
                   1024, 0, 1024, 512, 0, 512);
    }
    // In the next step, it is going to create some hitmap from the double
    // firing event couples. If you see a pixel in the hitmap, which has:
    // 1 entrY, it just fired in the first event.
    // 2 entries, it just fired in the second event.
    // 3 entries, it fired in both of the events.
    if (_isDoubleFiringAnalysis) {
      AIDAProcessor::tree(this)->mkdir(
          Form("First100Example/Double_Firing_Events"));
      AIDAProcessor::tree(this)->cd(
          Form("First100Example/Double_Firing_Events"));
      Double_Firing_Events_Hitmap[nInterestingCluster] =
          new TH2I(Form("Double_Firing_Events_Hitmap_%d", nInterestingCluster),
                   Form(" Double Firing Events Hitmap %d;Cluster width X "
                        "(pixel);Cluster width Y (pixel)",
                        nInterestingCluster),
                   1024, 0, 1024, 512, 0, 512);
    }
  }
  streamlog_out(DEBUG5) << "end of Booking histograms " << endl;
}

void EUTelProcessorClusterAnalysis::end() {
  if (_layerIndex == -1) {
    return;
  }
  for (int iSector = 0; iSector < _nSectors; iSector++) {
    for (size_t i = 0; i < symmetryGroups.size(); i++) {
      string binName;
      for (size_t j = 0; j < symmetryGroups[i].size(); j++) {
        clusterShapeHistoGroupedSector[iSector]->Fill(
            i, clusterShapeHistoSector[iSector]->GetBinContent(
                   symmetryGroups[i][j] + 1));
        if (j < symmetryGroups[i].size() - 1)
          binName += Form("%d,", symmetryGroups[i][j]);
        else
          binName += Form("%d", symmetryGroups[i][j]);
      }
      clusterShapeHistoGroupedSector[iSector]->GetXaxis()->SetBinLabel(
          i + 1, binName.c_str());
    }
  }
  //
  for (unsigned int iCluster = 0; iCluster < clusterVec.size(); iCluster++) {
    Cluster clustermap = clusterVec[iCluster];
    vector<int> Xmap = clustermap.getX();
    vector<int> Ymap = clustermap.getY();
    for (unsigned int ix = 0; ix < Xmap.size(); ix++) {
      if (Xmap.size() == Ymap.size())
        clusterShapeMap->SetBinContent(Xmap[ix] + 1, Ymap[ix] + 1, iCluster + 1,
                                       1);
    }
  }
  //
  if (_isDistanceSquareAnalysis) {
    bool DistanceSquareWarning = true;
    for (int i_sector = 0; i_sector < _nSectors; i_sector++) {
      if (GeneratedClustersHisto[i_sector]->GetEntries() != 0)
        DistanceSquareWarning = false;
      if (MissingClusterHisto[i_sector]->GetEntries() != 0)
        DistanceSquareWarning = false;
    }
    if (DistanceSquareWarning) {
      streamlog_out(WARNING)
          << "The histograms of the analysis of clusters with a smaller "
             "distance allowed between the pixels are empty! It could mean you "
             "runned the clusterAnalysis with same, or higher "
             "_sparseMinDistanceSquaredComparison than the "
             "_sparseMinDistanceSquared was in clustering."
          << endl;
    }
  }
  streamlog_out(MESSAGE4) << "The amount of processed events was " << _nEvents
                          << endl;
  streamlog_out(MESSAGE4) << "The amount of ignored clusters, because they "
                             "were at an x-border of a sector were: "
                          << _nTouchingBorderSectorClusters << endl;
  streamlog_out(MESSAGE4) << "The amount of ignored clusters, because they "
                             "were at a y-border of the chip were: "
                          << _nTouchingBorderYClusters << endl;
  streamlog_out(MESSAGE4) << "The amount of ignored clusters, because it "
                             "overlapped two sectors was: "
                          << _nOverlappingClusters << endl;
  streamlog_out(MESSAGE4)
      << "The amount of ignored clusters, because of hotPixels was: "
      << _nHotpixelClusters << endl;
  streamlog_out(MESSAGE4)
      << "The amount of ignored clusters, because of noiseMasks was: "
      << _nNoiseMaskClusters << endl;
  streamlog_out(MESSAGE4)
      << "The amount of ignored clusters, because of deadColumns was: "
      << _nDeadColumnClusters << endl;
  streamlog_out(MESSAGE4) << "For your information: Number of sectors: "
                          << _nSectors
                          << " ,columns per sector: " << _sectorWidth << endl;
  if (_writeClustersToFile)
    clusterAnalysisOutput.close();

  settingsFile << _nEvents << ";";
  for (int iSector = 0; iSector < _nSectors; iSector++)
    settingsFile << "0;"; // efficiency
  for (int iSector = 0; iSector < _nSectors; iSector++)
    settingsFile << "0;"; // number of tracks
  for (int iSector = 0; iSector < _nSectors; iSector++)
    settingsFile << "0;"; // number of tracks with assoc. hits
  settingsFile << endl;

  streamlog_out(MESSAGE4) << "ClusterAnalysis finished." << endl;
}

int EUTelProcessorClusterAnalysis::AddressToColumn(int ARegion, int ADoubleCol,
                                                   int AAddress) {
  int Column =
      ARegion * 32 + ADoubleCol * 2; // Double columns before ADoubleCol
  int LeftRight = 0;
  LeftRight = ((((AAddress % 4) == 1) || ((AAddress % 4) == 2))
                   ? 1
                   : 0); // Left or right column within the double column
  Column += LeftRight;
  return Column;
}

int EUTelProcessorClusterAnalysis::AddressToRow(int AAddress) {
  // Ok, this will get ugly
  int Row = AAddress / 2; // This is OK for the top-right and the bottom-left
                          // pixel within a group of 4
  return Row;
}
