#include "EUTelProcessorClusterAnalysis.h"
#include "EUTelAlignmentConstant.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelHistogramManager.h"
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
      _dutID(), _maxNumberOfPixels(3), _nSectors(),
      _nTouchingBorderSectorClusters(0), _nTouchingBorderYClusters(0),
      _nOverlappingClusters(0), _nHotpixelClusters(0), _nNoiseMaskClusters(0),
      _nDeadColumnClusters(0), _sectorWidth(), _energy(6.0), _chipID(),
      _irradiation(), _rate(""), _fillHistos(false),
      _hotPixelCollectionName(""), _nLayer(0), _xPixel(), _yPixel()

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
  _isFirstEvent = true;
}

void EUTelProcessorClusterAnalysis::init() {
  _nLayer = geo::gGeometry().nPlanes();
  const std::vector<int> &_planeID = geo::gGeometry().sensorIDsVec();
  for (int iz = 0; iz < _nLayer; iz++)
    if (_planeID[iz] == _dutID)
      _layerIndex = iz;

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

          // set the cluster
          cluster.set_values(clusterSize, X, Y);

          // fill the file with the pixel information of a cluster, after it has
          // been checked, whether the cluster is ok
          for (size_t iPixel = 0; iPixel < sparseData.size(); iPixel++) {
            // Output Pixel information to file, Comma (,) is coordinate
            // seperator, Whitespace ( ) is Pixel seperator
            clusterAnalysisOutput << X[iPixel] << "," << Y[iPixel] << " ";
          }
          // the next cluster will be the next thing written to the file, so
          // include the cluster seperator	semicolon
          clusterAnalysisOutput << ";";

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
  }

  // write the end event expression to the file, which is a linebreak
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
        15, 0.5, 15.5);
    clusterWidthYHisto[iSector] = new TH1I(
        Form("clusterWidthYHisto_%d", iSector),
        Form("Cluster width in Y in sector %d;Cluster width Y (pixel);a.u.",
             iSector),
        15, 0.5, 15.5);
    clusterSizeHisto[iSector] =
        new TH1I(Form("clusterSizeHisto_%d", iSector),
                 Form("Cluster size_%d;Cluster size (pixel);a.u.", iSector), 20,
                 0.5, 20.5);
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
  }
  streamlog_out(DEBUG5) << "end of Booking histograms " << endl;
}

void EUTelProcessorClusterAnalysis::end() {
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
