#include "EUTelProcessorAnalysisPALPIDEfs.h"
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

EUTelProcessorAnalysisPALPIDEfs aEUTelProcessorAnalysisPALPIDEfs;

EUTelProcessorAnalysisPALPIDEfs::EUTelProcessorAnalysisPALPIDEfs()
    : Processor("EUTelProcessorAnalysisPALPIDEfs"), _fillHistos(false),
      _inputFittedHitName(""), _inputColName(""), _trackCollectionName(""),
      _alignmentPAlpideCollectionName("alignmentPAlpide"),
      _alignmentCollectionName("alignment"),
      _preAlignmentCollectionName("prealign"), _zsDataCollectionName(""),
      zsInputDataCollectionVec(nullptr), _hotPixelCollectionName(""), limit(0.05),
      _dutID(6), _maxNumberOfPixels(3), _nPlanesWithMoreHits(4),
      _moreTracks(false), _energy(6.0), _writeShapes(false),
      _shapeOutputFileName("./shapeDistribution.txt"),
      _outputSettingsFolderName("./"), _chipID(), _irradiation(), _holesizeX(),
      _holesizeY(), _rate(""), _oneAlignmentCollection(false),
      _clusterAvailable(true), _hotpixelAvailable(true),
      _noiseMaskAvailable(true), _deadColumnAvailable(true), chi2Max(1),
      _nEvents(0), _nEventsFake(5), _nEventsWithTrack(0), _minTimeStamp(0),
      _nSectors(8), _chipVersion(3), _showFake(true), _realAssociation(false),
      nTracks(8), nTracksFake(8), nTracksPAlpide(8), nTracksPAlpideFake(8),
      nTracksAssociation(8), nTracksPAlpideAssociation(8), nFakeWithTrack(8, 0),
      nFakeWithoutTrack(8, 0), nFake(8, 0), nFakeWithTrackCorrected(8, 0),
      nDUThits(0), nNoPAlpideHit(0), nWrongPAlpideHit(0),
      nPlanesWithTooManyHits(0), xZero(0), yZero(0), xPitch(0), ySize(0),
      xPixel(0), yPixel(0), hotPixelCollectionVec(nullptr)

{
  _description = "Analysis of the fitted tracks";
  registerInputCollection(LCIO::TRACKERHIT, "InputFittedHitName",
                          "Name of the input fitted TrackerHit collection",
                          _inputFittedHitName, std::string("fithit"));
  registerInputCollection(LCIO::TRACKERHIT, "InputCollectionName",
                          "Name of the input TrackerHit collection",
                          _inputColName, std::string("colhit"));
  registerInputCollection(LCIO::TRACK, "TrackCollectionName",
                          "Input track collection name", _trackCollectionName,
                          std::string("track"));
  registerInputCollection(LCIO::LCGENERICOBJECT, "AlignmentConstantName",
                          "Alignment constant from the condition file",
                          _alignmentCollectionName, string("alignment"));
  registerInputCollection(LCIO::LCGENERICOBJECT, "AlignmentPAlpideConstantName",
                          "Alignment constant from the condition file",
                          _alignmentPAlpideCollectionName,
                          string("alignmentPAlpide"));
  registerInputCollection(LCIO::LCGENERICOBJECT, "PreAlignmentConstantName",
                          "PreAlignment constant from the condition file",
                          _preAlignmentCollectionName, string("prealign"));
  registerInputCollection(LCIO::TRACKERDATA, "ZSDataCollectionName",
                          "Input of Zero Suppressed data",
                          _zsDataCollectionName, string("zsdata"));
  registerProcessorParameter("HistogramFilling",
                             "Switch on or off the histogram filling",
                             _fillHistos, true);
  registerProcessorParameter(
      "HistoInfoFileName", "This is the name of the histogram information file",
      _histoInfoFileName, string("histoinfo.xml"));
  registerProcessorParameter(
      "Limit", "This is allowed distance between the track and the hit", limit,
      0.05);
  registerProcessorParameter("dutID", "This is the ID of the DUT", _dutID,
                             6);
  registerProcessorParameter("MaxNumberOfPixels",
                             "This is the maximum number of pixels in one "
                             "cluster for the clustershape analysis",
                             _maxNumberOfPixels, 3);
  registerProcessorParameter(
      "nPlanesWithMoreHits",
      "This is the maximum number of planes that can have more than one hit",
      _nPlanesWithMoreHits, 4);
  registerProcessorParameter("MoreTracks",
                             "More tracks are allowed in one event",
                             _moreTracks, false);
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
  registerOptionalParameter("Energy", "Particle energy [GeV]", _energy,
                            6.0);
  registerProcessorParameter("WriteShapes", "Write cluster shapes to file?",
                             _writeShapes, false);
  registerOptionalParameter(
      "ShapeOutputFileName", "This is the name of the file where the IDs of "
                             "the cluster shapes will be saved",
      _shapeOutputFileName, static_cast<string>("./shapeDistribution.txt"));
  registerOptionalParameter(
      "OutputSettingsFolderName",
      "Folder name where all the settings of each run will be saved",
      _outputSettingsFolderName, static_cast<string>("./"));
  EVENT::StringVec _stringVecExample;
  _stringVecExample.push_back(" ");
  registerOptionalParameter("ChipID", "Chip IDs", _chipID, _stringVecExample);
  registerOptionalParameter("Irradiation", "Irradiation level", _irradiation,
                            _stringVecExample);
  registerOptionalParameter("Rate", "Data taking rate", _rate,
                            static_cast<string>(""));
  registerProcessorParameter(
      "MinTimeStamp", "This is minimum timestamp required to consider an event",
      _minTimeStamp, 0.);
  registerOptionalParameter("ChipVersion", "Chip Version", _chipVersion,
                            3);

  //  float defaultHoleSizeX[2] = {1, 29}; // Need to be changed if eutelescope
  //  has -std=c++11 flag
  //  float defaultHoleSizeY[2] = {9, 12.5}; // Need to be changed if
  //  eutelescope has -std=c++11 flag

  registerOptionalParameter("HoleSizeX", "Size of the hole in X axis (mm)",
                            _holesizeX, std::vector<float>{1, 29});
  registerOptionalParameter("HoleSizeY", "Size of the hole in Y axis (mm)",
                            _holesizeY, std::vector<float>{9, 12.5});
  _isFirstEvent = true;
  registerProcessorParameter("ShowFake", "Show fake efficiency", _showFake,
                             true);
  registerProcessorParameter("RealAssociation", "Calculate track to hit "
                                                "association without allowing "
                                                "the tracks to share hits",
                             _realAssociation, false);
}

void EUTelProcessorAnalysisPALPIDEfs::init() {
  printParameters();
  int _nTelPlanes = geo::gGeometry().nPlanes();
  const std::vector<int> &_planeID = geo::gGeometry().sensorIDsVec();

  if (_chipVersion < 3)
    _nSectors = 4;
  else if (_chipVersion == 3)
    _nSectors = 8;
  else if (_chipVersion == 5)
    _nSectors = 4;
  else
    _nSectors = 1;

  for (int iz = 0; iz < _nTelPlanes; iz++)
    if (_planeID[iz] == _dutID) {
      dutZ = geo::gGeometry().getPlaneZPosition(iz);
      layerIndex = iz;
      xSize = geo::gGeometry().getPlaneXSize(layerIndex);
      ySize = geo::gGeometry().getPlaneYSize(layerIndex);
      xZero = geo::gGeometry().getPlaneXPosition(layerIndex);        // mm
      yZero = geo::gGeometry().getPlaneYPosition(layerIndex);        // mm
      xSize = geo::gGeometry().getPlaneXSize(layerIndex);            // mm
      ySize = geo::gGeometry().getPlaneYSize(layerIndex);            // mm
      xPitch = geo::gGeometry().getPlaneXPitch(layerIndex);          // mm
      yPitch = geo::gGeometry().getPlaneYPitch(layerIndex);          // mm
      xPointing[0] = geo::gGeometry().planeFlip1(layerIndex); // was -1 ;
      xPointing[1] = geo::gGeometry().planeFlip2(layerIndex); // was  0 ;
      yPointing[0] = geo::gGeometry().planeFlip3(layerIndex); // was  0 ;
      yPointing[1] = geo::gGeometry().planeFlip4(layerIndex); // was -1 ;
      xPixel = geo::gGeometry().getPlaneNumberOfPixelsX(layerIndex);
      yPixel = geo::gGeometry().getPlaneNumberOfPixelsY(layerIndex);
      try {
        gRotation[0] =
            geo::gGeometry().getPlaneZRotationDegrees(layerIndex); // Euler gamma ;
        gRotation[1] =
            geo::gGeometry().getPlaneYRotationDegrees(layerIndex); // Euler beta  ;
        gRotation[2] =
            geo::gGeometry().getPlaneXRotationDegrees(layerIndex); // Euler alpha ;
      } catch (...) {
        streamlog_out(MESSAGE5) << " no sensor rotation is given in the GEAR "
                                   "steering file, assume NONE "
                                << endl;
      }
      if ((gRotation[1] != 0 && gRotation[1] != 180) ||
          (gRotation[2] != 0 && gRotation[2] != 180))
        zDistance = sqrt(xSize * xSize + ySize * ySize);
      else
        zDistance = 0.1;
      gRotation[0] = gRotation[0] * 3.1415926 / 180.; //
      gRotation[1] = gRotation[1] * 3.1415926 / 180.; //
      gRotation[2] = gRotation[2] * 3.1415926 / 180.; //
    }
  float chi2MaxTemp[1] = {30};
  for (size_t i = 0; i < chi2Max.size(); i++)
    chi2Max[i] = chi2MaxTemp[i];
  Cluster cluster;
  cluster.FindReferenceClusters(clusterVec, _maxNumberOfPixels);
  xPairs = cluster.SymmetryPairs(clusterVec, "x");
  yPairs = cluster.SymmetryPairs(clusterVec, "y");
  symmetryGroups = cluster.sameShape(clusterVec);
  for (int iSector = 0; iSector < _nSectors; iSector++) {
    nTracks[iSector] = 0;
    nTracksPAlpide[iSector] = 0;
    nTracksAssociation[iSector] = 0;
    nTracksPAlpideAssociation[iSector] = 0;
    if (_showFake) {
      nTracksFake[iSector] = 0;
      nTracksPAlpideFake[iSector] = 0;
    }
  }
  bool newFile = false;
  string _outputSettingsFileName =
      _outputSettingsFolderName + Form("settings_DUT%d", _dutID) + ".txt";
  if (!std::ifstream(_outputSettingsFileName.c_str()))
    newFile = true;
  settingsFile.open(_outputSettingsFileName.c_str(), ios::out | ios::app);
  if (newFile && _chipVersion >= 3)
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
  else if (newFile && (_chipVersion == 2 || _chipVersion == 1)) {
    settingsFile << "Run number;Energy;Chip ID;Irradiation "
                    "level(0-nonIrradiated,1-2.5e12,2-1e13,3-700krad,4-"
                    "combined:1e13+700krad);Rate;BB;Ithr;Idb;Vcasn;Vaux;Vcasp;"
                    "Vreset;Threshold and their RMS for all four sectors;Noise "
                    "and their RMS for all four sectors;Readout delay;Trigger "
                    "delay;Strobe length;StrobeB length;Data (1) or noise "
                    "(0);Number of events;Efficiency,Number of tracks,Number "
                    "of tracks with associated hit for all sectors"
                 << endl;
  }
}

void EUTelProcessorAnalysisPALPIDEfs::processEvent(LCEvent *evt) {
  // HISTOGRAM INIT
  // =============================================================================
  if (_isFirstEvent) {
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    if (_fillHistos) {
      bookHistos();
    }
#endif
  }

  int nTrackPerEvent = 0;
  int nClusterAssociatedToTrackPerEvent = 0;
  int nClusterPerEvent = 0;

  stats->Fill(kAll);

  // EVENT CUTS
  // ===================================================================================
  if (evt->getParameters().getIntVal("FLAG") == 100)
    return; // Excluding events with too large clusters
  stats->Fill(kNoLargeClusters);

  if (evt->getTimeStamp() < _minTimeStamp)
    return; // Excluding events with too small time distance
  stats->Fill(kGoodTimeStamp);

  // FIRST EVENT
  // ==================================================================================
  if (_isFirstEvent) {
    // Hot pixel collection
    // -----------------------------------------------------------------------
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
          << "not found " << endl;
      _hotpixelAvailable = false;
    }
    if (_hotpixelAvailable) {
      hotData = dynamic_cast<TrackerDataImpl *>(
          hotPixelCollectionVec->getElementAt(layerIndex));
      auto sparseData = std::make_unique<
          EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>>(hotData);
      auto &pixelVec = sparseData->getPixels();
      for (auto &sparsePixel : pixelVec) {
        hotpixelHisto->Fill(sparsePixel.getXCoord() * xPitch + xPitch / 2.,
                            sparsePixel.getYCoord() * yPitch + yPitch / 2.);
      }
    }

    // Noise mask
    // ---------------------------------------------------------------------------------
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
        hotpixelHisto->Fill(x * xPitch + xPitch / 2., y * yPitch + yPitch / 2.);
      }
    } else
      _noiseMaskAvailable = false;

    // Dead column collection
    // ---------------------------------------------------------------------
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
    if (_deadColumnAvailable) {
      deadColumn = dynamic_cast<TrackerDataImpl *>(
          deadColumnCollectionVec->getElementAt(layerIndex));
      auto sparseData = std::make_unique<
          EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>>(deadColumn);
      auto &pixelVec = sparseData->getPixels();
      for (auto &sparsePixel : pixelVec) {
        deadColumnHisto->Fill(sparsePixel.getXCoord() * xPitch + xPitch / 2.,
                              sparsePixel.getYCoord() * yPitch + yPitch / 2.);
      }
    }

    // Writing output file
    // ------------------------------------------------------------------------
    if (_chipVersion >= 3) {
      settingsFile
          << evt->getRunNumber() << ";" << _energy << ";" << _chipID[layerIndex]
          << ";" << _chipVersion << ";" << _irradiation[layerIndex] << ";"
          << _rate << ";" << evt->getParameters().getFloatVal("BackBiasVoltage")
          << ";" << evt->getParameters().getIntVal(Form("Ithr_%d", layerIndex))
          << ";" << evt->getParameters().getIntVal(Form("Idb_%d", layerIndex))
          << ";" << evt->getParameters().getIntVal(Form("Vcasn_%d", layerIndex))
          << ";"
          << evt->getParameters().getIntVal(Form("Vcasn2_%d", layerIndex))
          << ";" << evt->getParameters().getIntVal(Form("Vclip_%d", layerIndex))
          << ";" << evt->getParameters().getIntVal(Form("Vcasp_%d", layerIndex))
          << ";"
          << evt->getParameters().getIntVal(Form("VresetP_%d", layerIndex))
          << ";"
          << evt->getParameters().getIntVal(Form("VresetD_%d", layerIndex))
          << ";";
    } else if (_chipVersion == 2 || _chipVersion == 1) {
      settingsFile
          << evt->getRunNumber() << ";" << _energy << ";" << _chipID[layerIndex]
          << ";" << _irradiation[layerIndex] << ";" << _rate << ";"
          << evt->getParameters().getFloatVal("BackBiasVoltage") << ";"
          << evt->getParameters().getIntVal(Form("Ithr_%d", layerIndex)) << ";"
          << evt->getParameters().getIntVal(Form("Idb_%d", layerIndex)) << ";"
          << evt->getParameters().getIntVal(Form("Vcasn_%d", layerIndex)) << ";"
          << evt->getParameters().getIntVal(Form("Vaux_%d", layerIndex)) << ";"
          << evt->getParameters().getIntVal(Form("Vcasp_%d", layerIndex)) << ";"
          << evt->getParameters().getIntVal(Form("VresetP_%d", layerIndex))
          << ";";
    }

    for (int iSector = 0; iSector < _nSectors; iSector++)
      settingsFile << evt->getParameters().getFloatVal(
                          Form("Thr_%d_%d", layerIndex, iSector))
                   << ";"
                   << evt->getParameters().getFloatVal(
                          Form("ThrRMS_%d_%d", layerIndex, iSector))
                   << ";";

    for (int iSector = 0; iSector < _nSectors; iSector++)
      settingsFile << evt->getParameters().getFloatVal(
                          Form("Noise_%d_%d", layerIndex, iSector))
                   << ";"
                   << evt->getParameters().getFloatVal(
                          Form("NoiseRMS_%d_%d", layerIndex, iSector))
                   << ";";

    settingsFile << evt->getParameters().getIntVal(
                        Form("m_readout_delay_%d", layerIndex))
                 << ";"
                 << evt->getParameters().getIntVal(
                        Form("m_trigger_delay_%d", layerIndex))
                 << ";"
                 << evt->getParameters().getIntVal(
                        Form("m_strobe_length_%d", layerIndex))
                 << ";"
                 << evt->getParameters().getIntVal(
                        Form("m_strobeb_length_%d", layerIndex))
                 << ";1;";

    // Done
    // ---------------------------------------------------------------------------------------
    _isFirstEvent = false;
  }

  // GENERAL EVENT PROCESSING
  // =====================================================================
  timeStampHisto->Fill(evt->getTimeStamp());

  // CHECKING INPUT DATA
  // =========================================================================

  // Input collection
  // -----------------------------------------------------------------------------
  LCCollection *col = nullptr;
  try {
    col = evt->getCollection(_inputColName);
  } catch (lcio::DataNotAvailableException &e) {
    streamlog_out(DEBUG5) << "Not able to get collection " << _inputColName
                          << "\nfrom event " << evt->getEventNumber()
                          << " in run " << evt->getRunNumber() << endl;
    return;
  }
  stats->Fill(kDataAvailable);

  // Fitted hits
  // ----------------------------------------------------------------------------------
  LCCollection *colFit = nullptr;
  bool fitHitAvailable = true;
  try {
    colFit = evt->getCollection(_inputFittedHitName);
  } catch (lcio::DataNotAvailableException &e) {
    streamlog_out(DEBUG5) << "Not able to get fit collection "
                          << _inputFittedHitName << "\nfrom event "
                          << evt->getEventNumber() << " in run "
                          << evt->getRunNumber() << endl;
    fitHitAvailable = false;
  }
  if (fitHitAvailable)
    stats->Fill(kFittedHitsAvailable);

  // Tracks
  // ---------------------------------------------------------------------------------------
  LCCollection *colTrack = nullptr;
  try {
    colTrack = evt->getCollection(_trackCollectionName);
  } catch (DataNotAvailableException& e) {
    fitHitAvailable = false;
  }
  if (fitHitAvailable)
    stats->Fill(kTracksAvailable);

  // Alignment
  // ------------------------------------------------------------------------------------
  LCCollectionVec *alignmentPAlpideCollectionVec = nullptr;
  LCCollectionVec *alignmentCollectionVec = nullptr;
  LCCollectionVec *preAlignmentCollectionVec = nullptr;
  try {
    alignmentCollectionVec = dynamic_cast<LCCollectionVec *>(
        evt->getCollection(_alignmentCollectionName));
    preAlignmentCollectionVec = dynamic_cast<LCCollectionVec *>(
        evt->getCollection(_preAlignmentCollectionName));
  } catch (...) {
    streamlog_out(WARNING2) << "Alignment collection not available" << endl;
    return;
  }
  stats->Fill(kAlignAvailable);

  try {
    alignmentPAlpideCollectionVec = dynamic_cast<LCCollectionVec *>(
        evt->getCollection(_alignmentPAlpideCollectionName));
  } catch (...) {
    if (evt->getEventNumber() == 0)
      streamlog_out(WARNING2) << "Only one alignment collection for pAlpide"
                              << endl;
    _oneAlignmentCollection = true;
  }
  if (_oneAlignmentCollection)
    stats->Fill(kOnlyOneAlignAvailable);

  // Clusters
  // -------------------------------------------------------------------------------------
  _clusterAvailable = true;
  try {
    zsInputDataCollectionVec = dynamic_cast<LCCollectionVec *>(
        evt->getCollection(_zsDataCollectionName));
    streamlog_out(DEBUG4) << "zsInputDataCollectionVec: "
                          << _zsDataCollectionName.c_str() << " found " << endl;
  } catch (lcio::DataNotAvailableException&) {
    streamlog_out(DEBUG4) << "zsInputDataCollectionVec: "
                          << _zsDataCollectionName.c_str() << " not found "
                          << endl;
    _clusterAvailable = false;
  }

  if (_clusterAvailable) {
    stats->Fill(kClustersAvailable);
    CellIDDecoder<TrackerDataImpl> cellDecoder(zsInputDataCollectionVec);
    for (size_t iCluster = 0; iCluster < zsInputDataCollectionVec->size();
         iCluster++) {
      TrackerDataImpl *zsData = dynamic_cast<TrackerDataImpl *>(
          zsInputDataCollectionVec->getElementAt(iCluster));
      if (static_cast<int>(cellDecoder(zsData)["sensorID"]) == _dutID)
        nClusterPerEvent++;
    }
  }

  // CHECKING ROTATION MATRIX
  // =====================================================================
  // TODO: Why is nothing done? Why is this not checked during ::init(...) ?
  if ((xPointing[0] == xPointing[1]) && (xPointing[0] == 0)) {
    streamlog_out(ERROR4)
        << "DUT has a singular rotation matrix. Sorry for quitting" << endl;
    stats->Fill(kSingularRotationMatrix);
  }

  if ((yPointing[0] == yPointing[1]) && (yPointing[0] == 0)) {
    streamlog_out(ERROR4)
        << "Detector DUT has a singular rotation matrix. Sorry for quitting"
        << endl;
    stats->Fill(kSingularRotationMatrix);
  }

  // VARIABLES
  // ====================================================================================
  int nFitHit = 0;
  if (fitHitAvailable)
    nFitHit = colFit->getNumberOfElements();
  bool hitmapFilled = false;
  vector<int> clusterAssosiatedToTrack;
  std::vector<std::vector<double>> pT;
  std::vector<std::vector<double>> pH;

  // Evaluating fake efficiency
  // ===================================================================
  if (_showFake) {
    std::vector<std::vector<double>> posFakeEvent;
    int nHit = col->getNumberOfElements();
    for (int ihit = 0; ihit < nHit; ihit++) {
      TrackerHit *hit = dynamic_cast<TrackerHit *>(col->getElementAt(ihit));
      if (hit != nullptr) {
        const double *pos = hit->getPosition();
        std::vector<double> posFakeHit;
        posFakeHit.push_back(pos[0]);
        posFakeHit.push_back(pos[1]);
        posFakeHit.push_back(pos[2]);
        posFakeEvent.push_back(posFakeHit);
      }
    }
    posFake.insert(posFake.begin(),
                   posFakeEvent); // TODO: why insertion at the beginning?
    if (posFake.size() > _nEventsFake) {
      posFakeTemp = posFake.back();
      posFake.pop_back();
    }
  }

  // FITTED HIT LOOP
  // ===============================================================================
  bool firstTrack = true;

  for (int ifit = 0; ifit < nFitHit; ifit++) {
    TrackerHit *fithit = dynamic_cast<TrackerHit *>(colFit->getElementAt(ifit));
    double fitpos[3] = {0., 0., 0.};
    const double *fitpos0 = fithit->getPosition();
    fitpos[0] = fitpos0[0];
    fitpos[1] = fitpos0[1];
    fitpos[2] = fitpos0[2];

    stats->Fill(kFittedHit);

    bool twoTracks = false;
    double yposfitPrev = 0.;
    double xposfitPrev = 0.;
    double maxDistInPlane = 0.1;

    // check for other close tracks
    // ---------------------------------------------------------------
    for (int i = ifit + 1; i < nFitHit; i++) {
      TrackerHit *fithitcheck =
          dynamic_cast<TrackerHit *>(colFit->getElementAt(i));
      const double *fitposcheck0 = fithitcheck->getPosition();
      if (abs(fitpos[2] - fitposcheck0[2]) < maxDistInPlane) {
        stats->Fill(kTwoCloseTracks);
        twoTracks = true;
        break;
      }
    }
    if (twoTracks && !_moreTracks) {
      streamlog_out(MESSAGE5) << "2 tracks in event " << evt->getEventNumber()
                              << ", skipping the event!" << endl;
      stats->Fill(kTwoCloseTracksRejected);
      break;
    }

    // at expected z position?
    if (fitpos[2] >= dutZ - zDistance && fitpos[2] <= dutZ + zDistance) {
      stats->Fill(kZposCorrect);
      double xposfit = 0;
      double yposfit = 0;

      // Applying alignment
      bool alignSuccess =
          RemoveAlign(preAlignmentCollectionVec, alignmentCollectionVec,
                      alignmentPAlpideCollectionVec, fitpos, xposfit, yposfit);
      // streamlog_out( MESSAGE4 ) << "Removed align position: x=" << xposfit <<
      // ", y=" << yposfit << ", z="  << fitpos[2]  << " for event " <<
      // _nEvents+1 << "." << endl; // Debug output
      if (alignSuccess)
        stats->Fill(kAlignmentRemoved);
      else
        stats->Fill(kAlignmentRemovalFailed);
      if (!alignSuccess && _isFirstEvent)
        cerr << "Removing alignment did not work!" << endl;

      // check x and y position
      if (xposfit > 0 && yposfit > 0 && xposfit < xSize && yposfit < ySize) {
        stats->Fill(kFittedHitOnChip);

        // reject hits in the border region
        if (xposfit < limit || xposfit > xSize - limit || yposfit < limit ||
            yposfit > ySize - limit) {
          stats->Fill(kFittedHitInBorderRegion);
          continue;
        }

        // sector determination
        int index = -1;
        for (int iSector = 0; iSector < _nSectors; iSector++) {
          if (xposfit > xSize / static_cast<double>(_nSectors) * iSector +
                            (iSector == 0 ? 0 : 1) * (2. * xPitch + limit) &&
              xposfit < xSize / static_cast<double>(_nSectors) * (iSector + 1) -
                            (iSector == _nSectors - 1 ? 0 : 1) *
                                (2. * xPitch + limit)) {
            index = iSector;
            break;
          }
        }
        if (index == -1) {
          stats->Fill(kUnknownSector);
          continue;
        }

        // reject tracks too close to hot pixels
        if (_hotpixelAvailable) {
          auto sparseData = std::make_unique<
              EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>>(hotData);
          bool hotpixel = false;
          auto &pixelVec = sparseData->getPixels();
          for (auto &sparsePixel : pixelVec) {
            if (abs(xposfit -
                    (sparsePixel.getXCoord() * xPitch + xPitch / 2.)) < limit &&
                abs(yposfit -
                    (sparsePixel.getYCoord() * yPitch + yPitch / 2.)) < limit) {
              hotpixel = true;
              break;
            }
          }
          if (hotpixel) {
            stats->Fill(kHotPixel);
            continue;
          }
        }

        // reject tracks too close to masked pixels
        if (_noiseMaskAvailable) {
          bool noisePixel = false;
          for (size_t iNoise = 0; iNoise < noiseMaskX.size(); iNoise++) {
            if (abs(xposfit - (noiseMaskX[iNoise] * xPitch + xPitch / 2.)) <
                    limit &&
                abs(yposfit - (noiseMaskY[iNoise] * yPitch + yPitch / 2.)) <
                    limit) {
              noisePixel = true;
              break;
            }
          }
          if (noisePixel) {
            stats->Fill(kMaskedPixel);
            continue;
          }
        }

        // reject tracks too close to dead columns
        if (_deadColumnAvailable) {
          auto sparseData = std::make_unique<
              EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>>(
              deadColumn);
          bool dead = false;
          auto &pixelVec = sparseData->getPixels();
          for (auto &sparsePixel : pixelVec) {
            if (abs(xposfit -
                    (sparsePixel.getXCoord() * xPitch + xPitch / 2.)) < limit) {
              dead = true;
              break;
            }
          }
          if (dead) {
            stats->Fill(kDeadColumn);
            continue;
          }
        }

        nTrackPerEvent++;
        int nAssociatedhits = 0;
        int nDUThitsEvent = 0;
        int nHit = col->getNumberOfElements();
        double yposPrev = 0.;
        double xposPrev = 0.;
        int nPlanesWithMoreHits = 0;
        bool unfoundTrack = false;
        bool unfoundTrackFake = false;
        bool pAlpideHit = false;
        bool firstHit = true;
        int nPAlpideHits = 0;

        // HIT LOOP
        // ===============================================================================
        // Find planes with multiple hits
        // ---------------------------------------------------------
        // ToDo: why was this loop inside the other hit loop below?
        bool hitOnSamePlane = false;
        for (int j = 0; j < nHit && firstHit &&
                        nPlanesWithMoreHits <= _nPlanesWithMoreHits;
             j++) {
          TrackerHit *hitcheck1 =
              dynamic_cast<TrackerHit *>(col->getElementAt(j));
          if (!hitcheck1)
            continue;
          const double *poscheck1 = hitcheck1->getPosition();
          if (poscheck1[2] <= dutZ + zDistance &&
              poscheck1[2] >= dutZ - zDistance)
            continue;
          for (int k = j + 1;
               k < nHit && nPlanesWithMoreHits <= _nPlanesWithMoreHits; k++) {
            TrackerHit *hitcheck2 =
                dynamic_cast<TrackerHit *>(col->getElementAt(k));
            if (!hitcheck2)
              continue;
            const double *poscheck2 = hitcheck2->getPosition();
            if ((abs(poscheck1[2] - poscheck2[2]) < maxDistInPlane) &&
                !hitOnSamePlane) {
              hitOnSamePlane = true;
              nPlanesWithMoreHits++;
              break;
            } else if (k == j + 1 && poscheck1[2] != poscheck2[2])
              hitOnSamePlane = false;
          }
        }
        if (nPlanesWithMoreHits > 0)
          stats->Fill(kHitsOnTheSamePlane);
        else
          stats->Fill(kSingleHitsOnly);

        firstHit = false;

        // Apply cut on the number of planes with multiple hits
        if (nPlanesWithMoreHits > _nPlanesWithMoreHits) {
          unfoundTrack = true;
          stats->Fill(kRejectedMultipleHits);
          break;
        }

        stats->Fill(kDenominator);

        std::vector<double> posFitHit;
        posFitHit.push_back(xposfit);
        posFitHit.push_back(yposfit);
        pT.push_back(posFitHit);

        // HIT LOOP
        // ===============================================================================
        // Find hit in DUT
        // ------------------------------------------------------------------------
        for (int ihit = 0; ihit < nHit; ihit++) {
          TrackerHit *hit = dynamic_cast<TrackerHit *>(col->getElementAt(ihit));
          double pos[3] = {0., 0., 0.};
          if (hit != nullptr) {
            const double *pos0 = hit->getPosition();
            pos[0] = pos0[0];
            pos[1] = pos0[1];
            pos[2] = pos0[2];

            // Hit in the DUT?
            if (pos[2] >= dutZ - zDistance && pos[2] <= dutZ + zDistance) {
              pAlpideHit = true;
              nPAlpideHits++;
              stats->Fill(kHitInDUT);

              // Determine position of the hit on the DUT
              pos[0] -= xZero;
              pos[1] -= yZero;
              _EulerRotationBack(pos, gRotation);

              double xpos = 0.;
              double ypos = 0.;
              _LayerRotationBack(pos, xpos, ypos);

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
              if (!hitmapFilled)
                hitmapHisto->Fill(xpos, ypos);
#endif
              nDUThitsEvent++;
              if (nDUThitsEvent > 1)
                nDUThits++;

              if (firstTrack) {
                std::vector<double> posHitHit;
                posHitHit.push_back(xpos);
                posHitHit.push_back(ypos);
                pH.push_back(posHitHit);
              }

              TrackImpl *track =
                  static_cast<TrackImpl *>(colTrack->getElementAt(ifit / 7));
              float chi2 = track->getChi2();

              if (abs(xpos - xposfit) < limit && abs(ypos - yposfit) < limit) {
                nAssociatedhits++;
                stats->Fill(kAssociatedHitInDUT);
                for (int jhit = ihit + 1; jhit < nHit; jhit++) {
                  TrackerHit *hitNext =
                      dynamic_cast<TrackerHit *>(col->getElementAt(jhit));
                  double posNext[3] = {0., 0., 0.};
                  if (hitNext != nullptr) {
                    const double *pos0Next = hitNext->getPosition();
                    posNext[0] = pos0Next[0];
                    posNext[1] = pos0Next[1];
                    posNext[2] = pos0Next[2];
                    if (posNext[2] < dutZ - zDistance ||
                        posNext[2] > dutZ + zDistance)
                      continue;
                    posNext[0] -= xZero;
                    posNext[1] -= yZero;

                    _EulerRotationBack(posNext, gRotation);

                    double xposNext, yposNext;
                    _LayerRotationBack(posNext, xposNext, yposNext);

                    if (abs(xposNext - xposfit) > limit ||
                        abs(yposNext - yposfit) > limit)
                      continue;
                    nAssociatedhits++;
                    if ((xpos - xposfit) * (xpos - xposfit) +
                            (ypos - yposfit) * (ypos - yposfit) <
                        (xposNext - xposfit) * (xposNext - xposfit) +
                            (yposNext - yposfit) * (yposNext - yposfit)) {
                      ihit = jhit;
                    } else {
                      xpos = xposNext;
                      ypos = yposNext;
                      pos[2] = posNext[2];
                      hit = hitNext;
                      ihit = jhit;
                    }
                  }
                }
                if (nDUThitsEvent > 1 && nAssociatedhits == 1) {
                  tmpHist->Fill(xposfitPrev, yposfitPrev);
                  nWrongPAlpideHit--;
                }
                if (nAssociatedhits > 1)
                  streamlog_out(DEBUG)
                      << nAssociatedhits
                      << " points for one track in DUT in event "
                      << evt->getEventNumber() << "\t" << xposPrev << "\t"
                      << yposPrev << "\t" << xpos << "\t" << ypos
                      << " Fit: " << xposfit << "\t" << yposfit
                      << " Number of planes with more than one hit: "
                      << nPlanesWithMoreHits << endl;
                if (_clusterAvailable) {
                  for (size_t idetector = 0;
                       idetector < zsInputDataCollectionVec->size();
                       idetector++) {
                    CellIDDecoder<TrackerDataImpl> cellDecoder(
                        zsInputDataCollectionVec);
                    TrackerDataImpl *zsData = dynamic_cast<TrackerDataImpl *>(
                        zsInputDataCollectionVec->getElementAt(idetector));
                    SparsePixelType type =
                        static_cast<SparsePixelType>(static_cast<int>(
                            cellDecoder(zsData)["sparsePixelType"]));
                    if (hit->getTime() == zsData->getTime()) {
                      nClusterAssociatedToTrackPerEvent++;
                      clusterAssosiatedToTrack.push_back(zsData->getTime());
                      int clusterSize = zsData->getChargeValues().size() / 4;
                      vector<int> X(clusterSize);
                      vector<int> Y(clusterSize);
                      Cluster cluster;
                      if (type == kEUTelGenericSparsePixel) {
                        vector<vector<int>> pixVector;
                        auto sparseData = EUTelTrackerDataInterfacerImpl<
                            EUTelGenericSparsePixel>(zsData);
                        for (size_t iPixel = 0; iPixel < sparseData.size();
                             iPixel++) {
                          auto &pixel = sparseData.at(iPixel);
                          X[iPixel] = pixel.getXCoord();
                          Y[iPixel] = pixel.getYCoord();
                          vector<int> pix;
                          pix.push_back(X[iPixel]);
                          pix.push_back(Y[iPixel]);
                          pixVector.push_back(pix);
                        }
                        cluster.set_values(clusterSize, X, Y);
                        clusterSizeHisto[index]->Fill(clusterSize);
                        int xMin = *min_element(X.begin(), X.end());
                        int xMax = *max_element(X.begin(), X.end());
                        int yMin = *min_element(Y.begin(), Y.end());
                        int yMax = *max_element(Y.begin(), Y.end());
                        int clusterWidthX = xMax - xMin + 1;
                        int clusterWidthY = yMax - yMin + 1;

                        if ((clusterWidthX > 3 || clusterWidthY > 3) &&
                            !emptyMiddle(pixVector))
                          for (size_t iPixel = 0; iPixel < pixVector.size();
                               iPixel++)
                            largeClusterHistos->Fill(pixVector[iPixel][0],
                                                     pixVector[iPixel][1]);
                        if (emptyMiddle(pixVector)) {
                          for (size_t iPixel = 0; iPixel < pixVector.size();
                               iPixel++)
                            circularClusterHistos->Fill(pixVector[iPixel][0],
                                                        pixVector[iPixel][1]);
                        }

                        clusterWidthXHisto[index]->Fill(clusterWidthX);
                        clusterWidthYHisto[index]->Fill(clusterWidthY);
                        clusterWidthXVsXHisto[index]->Fill(
                            fmod(xposfit, xPitch), clusterWidthX);
                        clusterWidthXVsXAverageHisto[index]->Fill(
                            fmod(xposfit, xPitch), clusterWidthX);
                        clusterWidthYVsYHisto[index]->Fill(
                            fmod(yposfit, yPitch), clusterWidthY);
                        clusterWidthYVsYAverageHisto[index]->Fill(
                            fmod(yposfit, yPitch), clusterWidthY);
                        // if (yposfit > _holesizeY[0] && yposfit <
                        // _holesizeY[1] && xposfit < _holesizeX[1] && xposfit >
                        // _holesizeX[0])
                        clusterSize2DHisto[index]->Fill(fmod(xposfit, xPitch),
                                                        fmod(yposfit, yPitch),
                                                        clusterSize);
                        // if (yposfit > _holesizeY[0] && yposfit <
                        // _holesizeY[1] && xposfit < _holesizeX[1] && xposfit >
                        // _holesizeX[0])
                        clusterSize2D2by2Histo[index]->Fill(
                            fmod(xposfit, 2 * xPitch),
                            fmod(yposfit, 2 * yPitch), clusterSize);
                        // if (yposfit > _holesizeY[0] && yposfit <
                        // _holesizeY[1] && xposfit < _holesizeX[1] && xposfit >
                        // _holesizeX[0])
                        clusterSize2DAverageHisto[index]->Fill(
                            fmod(xposfit, xPitch), fmod(yposfit, yPitch),
                            clusterSize);
                        // if (yposfit > _holesizeY[0] && yposfit <
                        // _holesizeY[1] && xposfit < _holesizeX[1] && xposfit >
                        // _holesizeX[0])
                        clusterSize2DAverage2by2Histo[index]->Fill(
                            fmod(xposfit, 2 * xPitch),
                            fmod(yposfit, 2 * yPitch), clusterSize);
                        nClusterVsXHisto[index]->Fill(fmod(xposfit, xPitch));
                        nClusterVsYHisto[index]->Fill(fmod(yposfit, yPitch));
                        // if (yposfit > _holesizeY[0] && yposfit <
                        // _holesizeY[1] && xposfit < _holesizeX[1] && xposfit >
                        // _holesizeX[0])
                        nClusterSizeHisto[index]->Fill(fmod(xposfit, xPitch),
                                                       fmod(yposfit, yPitch));
                        // if (yposfit > _holesizeY[0] && yposfit <
                        // _holesizeY[1] && xposfit < _holesizeX[1] && xposfit >
                        // _holesizeX[0])
                        nClusterSize2by2Histo[index]->Fill(
                            fmod(xposfit, 2 * xPitch),
                            fmod(yposfit, 2 * yPitch));
                        int clusterShape =
                            cluster.WhichClusterShape(cluster, clusterVec);
                        if (clusterShape >= 0) {
                          clusterShapeHisto->Fill(clusterShape);
                          clusterShapeHistoSector[index]->Fill(clusterShape);
                          clusterShapeX[clusterShape]->Fill(xMin);
                          clusterShapeY[clusterShape]->Fill(yMin);
                          clusterShape2D2by2[clusterShape]->Fill(
                              fmod(xposfit, 2 * xPitch),
                              fmod(yposfit, 2 * yPitch));
                          for (size_t iGroup = 0;
                               iGroup < symmetryGroups.size(); iGroup++)
                            for (size_t iMember = 0;
                                 iMember < symmetryGroups[iGroup].size();
                                 iMember++)
                              if (symmetryGroups[iGroup][iMember] ==
                                  clusterShape)
                                clusterShape2DGrouped2by2[iGroup]->Fill(
                                    fmod(xposfit, 2 * xPitch),
                                    fmod(yposfit, 2 * yPitch));
                        } else {
                          clusterShapeHisto->Fill(clusterVec.size());
                          clusterShapeHistoSector[index]->Fill(clusterShape);
                        }
                      }
                      break;
                    }
                  }
                }
                nTracksPAlpide[index]++;

                chi22DHisto->Fill(xposfit, yposfit, chi2);
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
                chi2Histo->Fill(chi2);
                scatteringAngleHisto->Fill(xposfit, yposfit,
                                           abs(xpos - xposfit));
                for (size_t i = 0; i < chi2Max.size(); i++) {
                  if (chi2 < chi2Max[i] && yposfit < _holesizeY[1] &&
                      yposfit > _holesizeY[0] && xposfit < _holesizeX[1] &&
                      xposfit > _holesizeX[0]) {
                    residualXPAlpide[chi2Max[i]][index]->Fill(xpos - xposfit);
                    residualYPAlpide[chi2Max[i]][index]->Fill(ypos - yposfit);
                    residualZPAlpide[chi2Max[i]][index]->Fill(pos[2] -
                                                              fitpos[2]);
                    residualXPixel[chi2Max[i]][index]->Fill(
                        fmod(xposfit, xPitch), fmod(yposfit, yPitch),
                        abs(xpos - xposfit));
                    residualYPixel[chi2Max[i]][index]->Fill(
                        fmod(xposfit, xPitch), fmod(yposfit, yPitch),
                        abs(ypos - yposfit));
                    residualXAveragePixel[chi2Max[i]][index]->Fill(
                        fmod(xposfit, xPitch), fmod(yposfit, yPitch),
                        abs(xpos - xposfit));
                    residualYAveragePixel[chi2Max[i]][index]->Fill(
                        fmod(xposfit, xPitch), fmod(yposfit, yPitch),
                        abs(ypos - yposfit));
                    nResidualXPixel[chi2Max[i]][index]->Fill(
                        fmod(xposfit, xPitch), fmod(yposfit, yPitch));
                    nResidualYPixel[chi2Max[i]][index]->Fill(
                        fmod(xposfit, xPitch), fmod(yposfit, yPitch));
                    residualXPixel2by2[chi2Max[i]][index]->Fill(
                        fmod(xposfit, 2 * xPitch), fmod(yposfit, 2 * yPitch),
                        abs(xpos - xposfit));
                    residualYPixel2by2[chi2Max[i]][index]->Fill(
                        fmod(xposfit, 2 * xPitch), fmod(yposfit, 2 * yPitch),
                        abs(ypos - yposfit));
                    residualXAveragePixel2by2[chi2Max[i]][index]->Fill(
                        fmod(xposfit, 2 * xPitch), fmod(yposfit, 2 * yPitch),
                        abs(xpos - xposfit));
                    residualYAveragePixel2by2[chi2Max[i]][index]->Fill(
                        fmod(xposfit, 2 * xPitch), fmod(yposfit, 2 * yPitch),
                        abs(ypos - yposfit));
                    nResidualXPixel2by2[chi2Max[i]][index]->Fill(
                        fmod(xposfit, 2 * xPitch), fmod(yposfit, 2 * yPitch));
                    nResidualYPixel2by2[chi2Max[i]][index]->Fill(
                        fmod(xposfit, 2 * xPitch), fmod(yposfit, 2 * yPitch));
                  } else if (chi2 < chi2Max[i] && (yposfit < _holesizeY[0] ||
                                                   yposfit > _holesizeY[1] ||
                                                   xposfit > _holesizeX[1] ||
                                                   xposfit < _holesizeX[0])) {
                    residualXPCBPAlpide[chi2Max[i]][index]->Fill(xpos -
                                                                 xposfit);
                    residualYPCBPAlpide[chi2Max[i]][index]->Fill(ypos -
                                                                 yposfit);
                    residualZPCBPAlpide[chi2Max[i]][index]->Fill(pos[2] -
                                                                 fitpos[2]);
                  }
                }
                tracksPAlpideHisto->Fill(xposfit, yposfit);
                efficiencyHisto->Fill(xposfit, yposfit);
                tracksPAlpidePixelHisto[index]->Fill(fmod(xposfit, xPitch),
                                                     fmod(yposfit, yPitch));
                tracksPAlpidePixel2by2Histo[index]->Fill(
                    fmod(xposfit, 2 * xPitch), fmod(yposfit, 2 * yPitch));
                efficiencyPixelHisto[index]->Fill(fmod(xposfit, xPitch),
                                                  fmod(yposfit, yPitch));
                efficiencyPixel2by2Histo[index]->Fill(
                    fmod(xposfit, 2 * xPitch), fmod(yposfit, 2 * yPitch));
#endif
                yposPrev = ypos;
                xposPrev = xpos;

                // Move here (XYZ)
              } else {
                stats->Fill(kHitNotInDUT);
                continue;
              }
              if (nAssociatedhits < 1 &&
                  nDUThitsEvent == 1) { // Move above to XYZ
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
                hitmapWrongHitHisto->Fill(xposfit, yposfit);
                chi2HistoNoHit->Fill(chi2);
                stats->Fill(kHitInDUTNotAssociated);
#endif
                nWrongPAlpideHit++;
                xposfitPrev = xposfit;
                yposfitPrev = yposfit;
              }
              nAssociatedhitsHisto->Fill(nAssociatedhits);
            } else
              continue;
          } else {
            stats->Fill(kHitNULL);
          }
        }
        // END OF HIT LOOP

        if (nAssociatedhits < 1) {
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
          stats->Fill(kNoHitInDUT);
#endif
        }

        firstTrack = false;
        nHitsPerEvent.push_back(nPAlpideHits);

        // FAKE EFFICIENCY DETERMINATION
        // ==========================================================
        if (_showFake) {
          if (posFake.size() >= _nEventsFake) {
            int nHitFake = posFakeTemp.size();
            int nPlanesWithMoreHitsFake = 0;
            bool firstHitFake = true;

            for (int ihit = 0; ihit < nHitFake; ihit++) {
              bool hitOnSamePlane = false;
              double pos[3] = {0., 0., 0.};
              for (int j = ihit;
                   j < nHitFake && firstHitFake &&
                   nPlanesWithMoreHitsFake <= _nPlanesWithMoreHits;
                   j++) {
                double poscheck1[3] = {0., 0., 0.};
                poscheck1[0] = posFakeTemp.at(j).at(0);
                poscheck1[1] = posFakeTemp.at(j).at(1);
                poscheck1[2] = posFakeTemp.at(j).at(2);
                if (poscheck1[2] <= dutZ + zDistance &&
                    poscheck1[2] >= dutZ - zDistance)
                  continue;
                for (int k = j + 1;
                     k < nHitFake &&
                     nPlanesWithMoreHitsFake <= _nPlanesWithMoreHits;
                     k++) {
                  double poscheck2[3] = {0., 0., 0.};
                  poscheck2[0] = posFakeTemp.at(k).at(0);
                  poscheck2[1] = posFakeTemp.at(k).at(1);
                  poscheck2[2] = posFakeTemp.at(k).at(2);
                  if ((abs(poscheck1[2] - poscheck2[2]) < maxDistInPlane) &&
                      !hitOnSamePlane) {
                    hitOnSamePlane = true;
                    nPlanesWithMoreHitsFake++;
                    break;
                  } else if (k == j + 1 && poscheck1[2] != poscheck2[2])
                    hitOnSamePlane = false;
                }
              }
              firstHitFake = false;
              if (nPlanesWithMoreHitsFake > _nPlanesWithMoreHits) {
                unfoundTrackFake = true;
                break;
              }

              double pos0[3] = {0., 0., 0.};
              pos0[0] = posFakeTemp.at(ihit).at(0);
              pos0[1] = posFakeTemp.at(ihit).at(1);
              pos0[2] = posFakeTemp.at(ihit).at(2);
              pos[0] = pos0[0];
              pos[1] = pos0[1];
              pos[2] = pos0[2];

              if (pos[2] >= dutZ - zDistance && pos[2] <= dutZ + zDistance) {
                pos[0] -= xZero;
                pos[1] -= yZero;
                _EulerRotationBack(pos, gRotation);

                double xpos, ypos;
                _LayerRotationBack(pos, xpos, ypos);

                if (abs(xpos - xposfit) < limit &&
                    abs(ypos - yposfit) < limit) {
                  for (int jhit = ihit + 1; jhit < nHitFake; jhit++) {
                    double posNext[3] = {0., 0., 0.};
                    double pos0Next[3] = {0., 0., 0.};
                    pos0Next[0] = posFakeTemp.at(jhit).at(0);
                    pos0Next[1] = posFakeTemp.at(jhit).at(1);
                    pos0Next[2] = posFakeTemp.at(jhit).at(2);
                    posNext[0] = pos0Next[0];
                    posNext[1] = pos0Next[1];
                    posNext[2] = pos0Next[2];
                    if (posNext[2] < dutZ - zDistance ||
                        posNext[2] > dutZ + zDistance)
                      continue;
                    posNext[0] -= xZero;
                    posNext[1] -= yZero;

                    _EulerRotationBack(posNext, gRotation);

                    double xposNext, yposNext;
                    _LayerRotationBack(posNext, xposNext, yposNext);

                    if (abs(xposNext - xposfit) > limit ||
                        abs(yposNext - yposfit) > limit)
                      continue;
                    if ((xpos - xposfit) * (xpos - xposfit) +
                            (ypos - yposfit) * (ypos - yposfit) <
                        (xposNext - xposfit) * (xposNext - xposfit) +
                            (yposNext - yposfit) * (yposNext - yposfit)) {
                      ihit = jhit;
                    } else {
                      xpos = xposNext;
                      ypos = yposNext;
                      pos[2] = posNext[2];
                      pos0[0] = pos0Next[0];
                      pos0[1] = pos0Next[1];
                      pos0[2] = pos0Next[2];
                      ihit = jhit;
                    }
                  }
                  nTracksPAlpideFake[index]++;
                }
              } else
                continue;
            }
          }
        }
        hitmapFilled = true;
        if (unfoundTrack) {
          nPlanesWithTooManyHits++;
          break;
        } else {
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
          tracksHisto->Fill(xposfit, yposfit);
          tracksPixelHisto[index]->Fill(fmod(xposfit, xPitch),
                                        fmod(yposfit, yPitch));
          tracksPixel2by2Histo[index]->Fill(fmod(xposfit, 2 * xPitch),
                                            fmod(yposfit, 2 * yPitch));
#endif
          nTracks[index]++;
          if (!pAlpideHit) {
            hitmapNoHitHisto->Fill(xposfit, yposfit);
            nNoPAlpideHit++;
            // chi2HistoNoHit->Fill(chi2);
            streamlog_out(MESSAGE5)
                << "Tracks without a hit in the DUT in event "
                << evt->getEventNumber() << endl;
          }
        }
        if (_showFake) {
          if (unfoundTrackFake) {
            break;
          } else if (posFake.size() >= _nEventsFake)
            nTracksFake[index]++;
        }
      }
    } else {
      stats->Fill(kZposWrong);
      continue;
    }
  }

  posFit.push_back(pT);

  _nEvents++;
  if (fitHitAvailable)
    _nEventsWithTrack++;

  if (_clusterAvailable) {
    for (size_t i = 0; i < zsInputDataCollectionVec->size(); i++) {
      CellIDDecoder<TrackerDataImpl> cellDecoder(zsInputDataCollectionVec);
      TrackerDataImpl *zsData = dynamic_cast<TrackerDataImpl *>(
          zsInputDataCollectionVec->getElementAt(i));
      int sensorID = static_cast<int>(cellDecoder(zsData)["sensorID"]);
      if (sensorID == _dutID) {
        bool isAssosiated = false;
        for (size_t iAssociatedCluster = 0;
             iAssociatedCluster < clusterAssosiatedToTrack.size();
             iAssociatedCluster++)
          if (zsData->getTime() ==
              clusterAssosiatedToTrack[iAssociatedCluster]) {
            isAssosiated = true;
            break;
          }
        if (!isAssosiated) {
          int index = -1;
          SparsePixelType type = static_cast<SparsePixelType>(
              static_cast<int>(cellDecoder(zsData)["sparsePixelType"]));
          if (type == kEUTelGenericSparsePixel) {
            Cluster cluster;
            int clusterSize = zsData->getChargeValues().size() / 4;
            vector<int> X(clusterSize);
            vector<int> Y(clusterSize);
            auto sparseData =
                EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>(zsData);
            for (size_t iPixel = 0; iPixel < sparseData.size(); iPixel++) {
              auto &pixel = sparseData.at(iPixel);
              X[iPixel] = pixel.getXCoord();
              Y[iPixel] = pixel.getYCoord();
              if (!fitHitAvailable)
                nFakeWithoutTrackHitmapHisto->Fill(X[iPixel], Y[iPixel]);
              else
                nFakeWithTrackHitmapHisto->Fill(X[iPixel], Y[iPixel]);
              nFakeHitmapHisto->Fill(X[iPixel], Y[iPixel]);
            }
            cluster.set_values(clusterSize, X, Y);
            float xCenter, yCenter;
            cluster.getCenterOfGravity(xCenter, yCenter);
            for (int iSector = 0; iSector < _nSectors; iSector++) {
              if (xCenter > xPixel / static_cast<double>(_nSectors) * iSector &&
                  xCenter < xPixel / static_cast<double>(_nSectors) * (iSector + 1)) {
                index = iSector;
                break;
              }
            }
            if (index != -1) {
              if (!fitHitAvailable)
                nFakeWithoutTrack[index] += clusterSize;
              else
                nFakeWithTrack[index] += clusterSize;
              nFake[index] += clusterSize;
            }
          }
        }
      }
    }
  }
  if (_clusterAvailable && fitHitAvailable) {
    for (size_t i = 0; i < zsInputDataCollectionVec->size(); i++) {
      CellIDDecoder<TrackerDataImpl> cellDecoder(zsInputDataCollectionVec);
      TrackerDataImpl *zsData = dynamic_cast<TrackerDataImpl *>(
          zsInputDataCollectionVec->getElementAt(i));
      int sensorID = static_cast<int>(cellDecoder(zsData)["sensorID"]);
      if (sensorID == _dutID) {
        SparsePixelType type = static_cast<SparsePixelType>(
            static_cast<int>(cellDecoder(zsData)["sparsePixelType"]));
        if (type == kEUTelGenericSparsePixel) {
          Cluster cluster;
          int clusterSize = zsData->getChargeValues().size() / 4;
          vector<int> X(clusterSize);
          vector<int> Y(clusterSize);
          vector<vector<int>> pixVector;
          auto sparseData =
              EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>(zsData);
          for (size_t iPixel = 0; iPixel < sparseData.size(); iPixel++) {
            auto &pixel = sparseData.at(iPixel);
            X[iPixel] = pixel.getXCoord();
            Y[iPixel] = pixel.getYCoord();
          }
          cluster.set_values(clusterSize, X, Y);
          float xCenter, yCenter;
          cluster.getCenterOfGravity(xCenter, yCenter);
          bool isAssosiated = false;
          for (int ifit = 0; ifit < nFitHit; ifit++) {
            TrackerHit *fithit =
                dynamic_cast<TrackerHit *>(colFit->getElementAt(ifit));
            double fitpos[3] = {0., 0., 0.};
            const double *fitpos0 = fithit->getPosition();
            fitpos[0] = fitpos0[0];
            fitpos[1] = fitpos0[1];
            fitpos[2] = fitpos0[2];
            if (fitpos[2] >= dutZ - zDistance &&
                fitpos[2] <= dutZ + zDistance) {
              double xposfit = 0, yposfit = 0;
              RemoveAlign(preAlignmentCollectionVec, alignmentCollectionVec,
                          alignmentPAlpideCollectionVec, fitpos, xposfit,
                          yposfit);

              if (abs(xposfit - static_cast<double>(xCenter) / xPixel * xSize) < limit &&
                  abs(yposfit - static_cast<double>(yCenter) / yPixel * ySize) < limit) {
                isAssosiated = true;
                break;
              }
            }
          }
          if (isAssosiated)
            continue;
          for (size_t iPixel = 0; iPixel < X.size(); iPixel++)
            nFakeWithTrackHitmapCorrectedHisto->Fill(X[iPixel], Y[iPixel]);
          int index = -1;
          for (int iSector = 0; iSector < _nSectors; iSector++) {
            if (xCenter > xPixel / static_cast<double>(_nSectors) * iSector &&
                xCenter < xPixel / static_cast<double>(_nSectors) * (iSector + 1)) {
              index = iSector;
              break;
            }
          }
          if (index != -1) {
            nFakeWithTrackCorrected[index] += clusterSize;
          }
        }
      }
    }
  }
  nTrackPerEventHisto->Fill(nTrackPerEvent);
  nClusterAssociatedToTrackPerEventHisto->Fill(
      nClusterAssociatedToTrackPerEvent);
  nClusterPerEventHisto->Fill(nClusterPerEvent);
  if (_realAssociation) // Different version of track to hit association written
                        // by Martijn Dietze. It maximizes the number of
                        // association while doesn't allow hits to be shared
                        // between tracks.
  {
    int nH = pH.size();
    int nT = pT.size();
    std::vector<int> aH(nH);
    std::vector<int> aT(nT);
    for (int i = 0; i < nH; i++)
      aH[i] = -1;
    for (int i = 0; i < nT; i++)
      aT[i] = -1;
    int temp;
    bool changed = true;
    while (changed) {
      changed = false;
      for (int iT = 0; iT < nT; iT++) {
        if (aT[iT] == -1) {
          int associations = 0;
          for (int iH = 0; iH < nH; iH++) {
            if (aH[iH] == -1) {
              if (abs(pH.at(iH).at(0) - pT.at(iT).at(0)) < limit &&
                  abs(pH.at(iH).at(1) - pT.at(iT).at(1)) < limit) {
                associations++;
                temp = iH;
              }
            }
          }
          if (associations == 1) {
            changed = true;
            aH[temp] = iT;
            aT[iT] = temp;
          }
          if (associations == 0) {
            aT[iT] = -2;
          }
        }
      }
    }
    std::vector<int> aTFinal(nT);
    std::vector<int> aHFinal(nH);
    std::vector<int> order;
    for (int iT = 0; iT < nT; iT++) {
      if (aT[iT] == -1)
        order.push_back(iT);
    }
    if (order.size() > 0) {
      unsigned int maxAssociations = 0;
      double minDistance = 1e10;
      bool run1 = true;
      do {
        double totaldistance = 0;
        unsigned int associations = 0;
        std::vector<int> aTTemp(nT);
        std::vector<int> aHTemp(nH);
        for (int iT = 0; iT < nT; iT++)
          aTTemp[iT] = aT[iT];
        for (int iH = 0; iH < nH; iH++)
          aHTemp[iH] = aH[iH];
        bool run2 = true;
        for (unsigned int i = 0; i < order.size() && run2; i++) {
          int iT = order[i];
          double min = 1e10;
          bool associated = false;
          for (int iH = 0; iH < nH; iH++) {
            if (aHTemp[iH] == -1) {
              if (abs(pH.at(iH).at(0) - pT.at(iT).at(0)) < limit &&
                  abs(pH.at(iH).at(1) - pT.at(iT).at(1)) < limit) {
                double dist =
                    sqrt(pow(abs(pH.at(iH).at(0) - pT.at(iT).at(0)), 2) +
                         pow(abs(pH.at(iH).at(1) - pT.at(iT).at(1)), 2));
                if (dist < min) {
                  min = dist;
                  temp = iH;
                  associated = true;
                }
              }
            }
          }
          if (associated) {
            associations++;
            aHTemp[temp] = iT;
            aTTemp[iT] = temp;
            totaldistance += min;
          }
          if (totaldistance > minDistance)
            run2 = false;
          if ((i - associations + 1) > (order.size() - maxAssociations))
            run2 = false;
        }
        if (associations == order.size())
          run1 = false;
        if (associations >= maxAssociations) {
          if (totaldistance < minDistance) {
            for (int iT = 0; iT < nT; iT++)
              aTFinal[iT] = aTTemp[iT];
            for (int iH = 0; iH < nH; iH++)
              aHFinal[iH] = aHTemp[iH];
            maxAssociations = associations;
            minDistance = totaldistance;
          }
        }
      } while (std::next_permutation(order.begin(), order.end()) && run1);
    } else {
      for (int iT = 0; iT < nT; iT++)
        aTFinal[iT] = aT[iT];
      for (int iH = 0; iH < nH; iH++)
        aHFinal[iH] = aH[iH];
    }
    for (int iT = 0; iT < nT; iT++) {
      int index = -1;
      for (int iSector = 0; iSector < _nSectors; iSector++) {
        if (pT.at(iT).at(0) >
                xSize / static_cast<double>(_nSectors) * iSector +
                    (iSector == 0 ? 0 : 1) * (2. * xPitch + limit) &&
            pT.at(iT).at(0) <
                xSize / static_cast<double>(_nSectors) * (iSector + 1) -
                    (iSector == _nSectors ? 0 : 1) * (2. * xPitch + limit)) {
          index = iSector;
          break;
        }
      }
      if (index == -1)
        continue;
      if (aTFinal[iT] >= 0)
        nTracksPAlpideAssociation[index]++;
      nTracksAssociation[index]++;
    }
  }
}

#ifdef MARLIN_USE_AIDA
void EUTelProcessorAnalysisPALPIDEfs::bookHistos() {
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
  tracksHisto = new TH2I("tracksHisto", "Tracks;X (mm);Y (mm)", xPixel, 0,
                         xSize, yPixel, 0, ySize);
  tracksPAlpideHisto =
      new TH2I("tracksPAlpideHisto", "Tracks found in pALPIDEfs;X (mm);Y (mm)",
               xPixel, 0, xSize, yPixel, 0, ySize);
  efficiencyHisto =
      new TH2F("efficiencyHisto", "Efficiency of the pALPIDEfs;X (mm);Y (mm)",
               xPixel, 0, xSize, yPixel, 0, ySize);
  hitmapHisto = new TH2I("hitmapHisto", "Hitmap of pALPIDEfs;X (mm);Y (mm)",
                         xPixel, 0, xSize, yPixel, 0, ySize);
  hitmapNoHitHisto = new TH2I(
      "hitmapNoHitHisto",
      "Hitmap of tracks which didn't have a hit in the pALPIDEfs;X (mm);Y (mm)",
      xPixel, 0, xSize, yPixel, 0, ySize);
  hitmapWrongHitHisto = new TH2I("hitmapWrongHitHisto",
                                 "Hitmap of tracks which didn't have a hit in "
                                 "the pALPIDEfs close to them;X (mm);Y (mm)",
                                 xPixel, 0, xSize, yPixel, 0, ySize);
  nFakeHitmapHisto =
      new TH2I("nFakeHitmapHisto",
               "Hit map of hits not associated to tracks;X (pixel);Y (pixel)",
               xPixel, 0, xPixel, yPixel, 0, yPixel);
  nFakeWithTrackHitmapHisto = new TH2I(
      "nFakeWithTrackHitmapHisto", "Hit map of hits not associated to tracks "
                                   "with track in event;X (pixel);Y (pixel)",
      xPixel, 0, xPixel, yPixel, 0, yPixel);
  nFakeWithTrackHitmapCorrectedHisto =
      new TH2I("nFakeWithTrackHitmapCorrectedHisto",
               "Corrected hit map of hits not associated to tracks with track "
               "in event;X (pixel);Y (pixel)",
               xPixel, 0, xPixel, yPixel, 0, yPixel);
  nFakeWithoutTrackHitmapHisto =
      new TH2I("nFakeWithoutTrackHitmapHisto", "Hit map of hits not associated "
                                               "to trackswithout track in "
                                               "event;X (pixel);Y (pixel)",
               xPixel, 0, xPixel, yPixel, 0, yPixel);
  scatteringAngleHisto =
      new TProfile2D("scatteringAngleHisto", "Scattering angles;X (mm);Y (mm)",
                     xPixel, 0, xSize, yPixel, 0, ySize);
  chi22DHisto = new TProfile2D("chi22DHisto", "#chi^{2};X (mm);Y (mm)", xPixel,
                               0, xSize, yPixel, 0, ySize);
  tmpHist = new TH2I("tmpHist", "", xPixel, 0, xSize, yPixel, 0, ySize);
  chi2Histo = new TH1I("chi2Histo",
                       "#chi^{2} of tracks used for the analysis;#chi^{2};a.u.",
                       100, 0, 50);
  chi2HistoNoHit = new TH1I(
      "chi2HistoNoHit",
      "#chi^{2} of tracks with out a matched hit;#chi^{2};a.u.", 100, 0, 50);
  hotpixelHisto =
      new TH2I("hotpixelHisto", "Hot pixels in the DUT;X (mm);Y (mm)", xPixel,
               0, xSize, yPixel, 0, ySize);
  deadColumnHisto =
      new TH2I("deadColumnHisto", "Dead columns in the DUT;X (mm);Y (mm)",
               xPixel, 0, xSize, yPixel, 0, ySize);
  largeClusterHistos = new TH2I(
      "largeClusterHisto", "Clusters with more than 3 pixels;X (mm);Y (mm)",
      xPixel, 0, xPixel, yPixel, 0, yPixel);
  circularClusterHistos = new TH2I(
      "circularClusterHisto",
      "Circular clusters (with missing hits in the middle);X (mm);Y (mm)",
      xPixel, 0, xPixel, yPixel, 0, yPixel);
  timeStampHisto =
      new TH1I("timeStampHisto", "Distribution of the time stamp of the "
                                 "events; Time stamp (in 12.5 ns units)",
               1000, 0, 50000);
  nFakeHisto =
      new TH1F("nFakeHisto", "Noise occupancy per sector for all "
                             "events;Sector;Noise occupancy (/event/pixel)",
               _nSectors, 0, _nSectors);
  nFakeWithTrackHisto = new TH1F("nFakeWithTrackHisto",
                                 "Noise occupancy per sector for events with "
                                 "track;Sector;Noise occupancy (/event/pixel)",
                                 _nSectors, 0, _nSectors);
  nFakeWithTrackCorrectedHisto =
      new TH1F("nFakeWithTrackCorrectedHisto",
               "Corrected noise occupancy per sector for events with "
               "track;Sector;Noise occupancy (/event/pixel)",
               _nSectors, 0, _nSectors);
  nFakeWithoutTrackHisto = new TH1F(
      "nFakeWithoutTrackHisto", "Noise occupancy per sector for events without "
                                "track;Sector;Noise occupancy (/event/pixel)",
      _nSectors, 0, _nSectors);
  nTrackPerEventHisto =
      new TH1I("nTrackPerEventHisto",
               "Number of tracks per event;Number of tracks;a.u.", 30, 0, 30);
  nClusterAssociatedToTrackPerEventHisto =
      new TH1I("nClusterAssociatedToTrackPerEventHisto",
               "Number of clusters associated to tracks per event;Number of "
               "clusters;a.u.",
               30, 0, 30);
  nClusterPerEventHisto = new TH1I(
      "nClusterPerEventHisto",
      "Number of clusters per event;Number of clusters;a.u.", 30, 0, 30);
  nAssociatedhitsHisto =
      new TH1I("nAssociatedhitsHisto", "Number of hits in search region of "
                                       "track per event;Number of hits in "
                                       "search region of track;a.u.",
               30, 0, 30);
  nAssociatedtracksHisto =
      new TH1I("nAssociatedtracksHisto", "Number of tracks in search region of "
                                         "track per event;Number of tracks in "
                                         "search region of track;a.u.",
               30, 0, 30);

  stats = new TH1I("stats", "statistics of events, cuts and properties", 100,
                   0., 100.);
  stats->GetXaxis()->SetBinLabel(1 + kAll, "All");
  stats->GetXaxis()->SetBinLabel(1 + kNoLargeClusters, "No large clusters");
  stats->GetXaxis()->SetBinLabel(1 + kGoodTimeStamp, "Good timestamp");
  stats->GetXaxis()->SetBinLabel(1 + kDataAvailable, "Data available");
  stats->GetXaxis()->SetBinLabel(1 + kFittedHitsAvailable,
                                 "Fitted hits available");
  stats->GetXaxis()->SetBinLabel(1 + kTracksAvailable, "Tracks available");
  stats->GetXaxis()->SetBinLabel(1 + kAlignAvailable, "Aligment available");
  stats->GetXaxis()->SetBinLabel(1 + kOnlyOneAlignAvailable,
                                 "Only one aligment available");
  stats->GetXaxis()->SetBinLabel(1 + kClustersAvailable, "Clusters available");
  stats->GetXaxis()->SetBinLabel(1 + kSingularRotationMatrix,
                                 "Singular rotation matrix");
  stats->GetXaxis()->SetBinLabel(1 + kFittedHit, "Fitted hit found");
  stats->GetXaxis()->SetBinLabel(1 + kTwoCloseTracks, "Found two close tracks");
  stats->GetXaxis()->SetBinLabel(1 + kTwoCloseTracksRejected,
                                 "Event rejected, two close tracks");
  stats->GetXaxis()->SetBinLabel(1 + kZposCorrect, "z-position matches DUT");
  stats->GetXaxis()->SetBinLabel(1 + kZposWrong,
                                 "z-position does not match DUT");
  stats->GetXaxis()->SetBinLabel(1 + kAlignmentRemoved, "Aligment removed");
  stats->GetXaxis()->SetBinLabel(1 + kAlignmentRemovalFailed,
                                 "Aligment removal failed");
  stats->GetXaxis()->SetBinLabel(1 + kFittedHitOnChip, "Fitted hit on chip");
  stats->GetXaxis()->SetBinLabel(1 + kFittedHitNotOnChip,
                                 "Fitted hit not on chip");
  stats->GetXaxis()->SetBinLabel(1 + kFittedHitInBorderRegion,
                                 "Fited hit in border region");
  stats->GetXaxis()->SetBinLabel(1 + kUnknownSector,
                                 "Sector determination failed");
  stats->GetXaxis()->SetBinLabel(1 + kHotPixel, "Hit next to a hot pixel");
  stats->GetXaxis()->SetBinLabel(1 + kMaskedPixel,
                                 "Hit next to a masked pixel");
  stats->GetXaxis()->SetBinLabel(1 + kDeadColumn, "Hit next to a dead columns");
  stats->GetXaxis()->SetBinLabel(1 + kHitsOnTheSamePlane,
                                 "Found multiple hits on the same plane");
  stats->GetXaxis()->SetBinLabel(1 + kSingleHitsOnly, "Single hits only");
  stats->GetXaxis()->SetBinLabel(1 + kRejectedMultipleHits,
                                 "Event rejected, multiple hits");
  stats->GetXaxis()->SetBinLabel(1 + kDenominator, "Denominator");
  stats->GetXaxis()->SetBinLabel(1 + kHitInDUT, "Hit in DUT");
  stats->GetXaxis()->SetBinLabel(1 + kAssociatedHitInDUT,
                                 "Associated hit in DUT");
  stats->GetXaxis()->SetBinLabel(1 + kNoHitInDUT, "No Hit in DUT");
  stats->GetXaxis()->SetBinLabel(1 + kHitNotInDUT, "Hit not in DUT");
  stats->GetXaxis()->SetBinLabel(1 + kHitInDUTNotAssociated,
                                 "Hit in DUT not Associated");
  stats->GetXaxis()->SetBinLabel(1 + kHitNULL, "Hit object NULL");

  for (int iSector = 0; iSector < _nSectors; iSector++) {
    AIDAProcessor::tree(this)->mkdir(Form("Sector_%d", iSector));
    AIDAProcessor::tree(this)->cd(Form("Sector_%d", iSector));
    for (size_t i = 0; i < chi2Max.size(); i++) {
      residualXPAlpide[chi2Max[i]][iSector] =
          new TH1I(Form("residualXPAlpide_%.1f_%d", chi2Max[i], iSector),
                   Form("Residual X (Max chi2 = %.1f), sector %d;X (mm);a.u.",
                        chi2Max[i], iSector),
                   300, -0.2, 0.2);
      residualYPAlpide[chi2Max[i]][iSector] =
          new TH1I(Form("residualYPAlpide_%.1f_%d", chi2Max[i], iSector),
                   Form("Residual Y (Max chi2 = %.1f), sector %d;Y (mm);a.u.",
                        chi2Max[i], iSector),
                   300, -0.2, 0.2);
      residualZPAlpide[chi2Max[i]][iSector] =
          new TH1I(Form("residualZPAlpide_%.1f_%d", chi2Max[i], iSector),
                   Form("Residual Z (Max chi2 = %.1f), sector %d;Z (mm);a.u.",
                        chi2Max[i], iSector),
                   100, -0.3, 0.3);
      residualXPCBPAlpide[chi2Max[i]][iSector] =
          new TH1I(Form("residualXPCBPAlpide_%.1f_%d", chi2Max[i], iSector),
                   Form("Residual X not uner the hole (Max chi2 = %.1f), "
                        "sector %d;X (mm);a.u.",
                        chi2Max[i], iSector),
                   300, -0.2, 0.2);
      residualYPCBPAlpide[chi2Max[i]][iSector] =
          new TH1I(Form("residualYPCBPAlpide_%.1f_%d", chi2Max[i], iSector),
                   Form("Residual Y not uner the hole (Max chi2 = %.1f), "
                        "sector %d;Y (mm);a.u.",
                        chi2Max[i], iSector),
                   300, -0.2, 0.2);
      residualZPCBPAlpide[chi2Max[i]][iSector] =
          new TH1I(Form("residualZPCBPAlpide_%.1f_%d", chi2Max[i], iSector),
                   Form("Residual Z not uner the hole (Max chi2 = %.1f), "
                        "sector %d;Z (mm);a.u.",
                        chi2Max[i], iSector),
                   100, -0.3, 0.3);
      residualXPixel[chi2Max[i]][iSector] = new TH2F(
          Form("residualXPixel_%.1f_%d", chi2Max[i], iSector),
          "Residual X added up;X (mm);Y (mm)", static_cast<int>(xSize / xPixel * 1000), 0,
          xSize / xPixel, static_cast<int>(ySize / yPixel * 1000), 0, ySize / yPixel);
      residualYPixel[chi2Max[i]][iSector] = new TH2F(
          Form("residualYPixel_%.1f_%d", chi2Max[i], iSector),
          "Residual Y added up;X (mm);Y (mm)", static_cast<int>(xSize / xPixel * 1000), 0,
          xSize / xPixel, static_cast<int>(ySize / yPixel * 1000), 0, ySize / yPixel);
      residualXAveragePixel[chi2Max[i]][iSector] = new TH2F(
          Form("residualAverageXPixel_%.1f_%d", chi2Max[i], iSector),
          "Average residual X;X (mm);Y (mm)", static_cast<int>(xSize / xPixel * 1000), 0,
          xSize / xPixel, static_cast<int>(ySize / yPixel * 1000), 0, ySize / yPixel);
      residualYAveragePixel[chi2Max[i]][iSector] = new TH2F(
          Form("residualAverageYPixel_%.1f_%d", chi2Max[i], iSector),
          "Average residual Y;X (mm);Y (mm)", static_cast<int>(xSize / xPixel * 1000), 0,
          xSize / xPixel, static_cast<int>(ySize / yPixel * 1000), 0, ySize / yPixel);
      nResidualXPixel[chi2Max[i]][iSector] =
          new TH2F(Form("nResidualXPixel_%.1f_%d", chi2Max[i], iSector),
                   "Number of tracks used for residual X plot;X (mm);Y (mm)",
                   static_cast<int>(xSize / xPixel * 1000), 0, xSize / xPixel,
                   static_cast<int>(ySize / yPixel * 1000), 0, ySize / yPixel);
      nResidualYPixel[chi2Max[i]][iSector] =
          new TH2F(Form("nResidualYPixel_%.1f_%d", chi2Max[i], iSector),
                   "Number of tracks used for residual Y plot;X (mm);Y (mm)",
                   static_cast<int>(xSize / xPixel * 1000), 0, xSize / xPixel,
                   static_cast<int>(ySize / yPixel * 1000), 0, ySize / yPixel);
      residualXPixel2by2[chi2Max[i]][iSector] =
          new TH2F(Form("residualXPixel2by2_%.1f_%d", chi2Max[i], iSector),
                   "Residual X added up 2by2;X (mm);Y (mm)",
                   static_cast<int>(xSize / xPixel * 2000), 0, 2 * xSize / xPixel,
                   static_cast<int>(ySize / yPixel * 2000), 0, 2 * ySize / yPixel);
      residualYPixel2by2[chi2Max[i]][iSector] =
          new TH2F(Form("residualYPixel2by2_%.1f_%d", chi2Max[i], iSector),
                   "Residual Y added up 2by2;X (mm);Y (mm)",
                   static_cast<int>(xSize / xPixel * 2000), 0, 2 * xSize / xPixel,
                   static_cast<int>(ySize / yPixel * 2000), 0, 2 * ySize / yPixel);
      residualXAveragePixel2by2[chi2Max[i]][iSector] = new TH2F(
          Form("residualAverageXPixel2by2_%.1f_%d", chi2Max[i], iSector),
          "Average residual X 2by2;X (mm);Y (mm)", static_cast<int>(xSize / xPixel * 2000),
          0, 2 * xSize / xPixel, static_cast<int>(ySize / yPixel * 2000), 0,
          2 * ySize / yPixel);
      residualYAveragePixel2by2[chi2Max[i]][iSector] = new TH2F(
          Form("residualAverageYPixel2by2_%.1f_%d", chi2Max[i], iSector),
          "Average residual Y 2by2;X (mm);Y (mm)", static_cast<int>(xSize / xPixel * 2000),
          0, 2 * xSize / xPixel, static_cast<int>(ySize / yPixel * 2000), 0,
          2 * ySize / yPixel);
      nResidualXPixel2by2[chi2Max[i]][iSector] = new TH2F(
          Form("nResidualXPixel2by2_%.1f_%d", chi2Max[i], iSector),
          "Number of tracks used for residual X 2by2 plot;X (mm);Y (mm)",
          static_cast<int>(xSize / xPixel * 2000), 0, 2 * xSize / xPixel,
          static_cast<int>(ySize / yPixel * 2000), 0, 2 * ySize / yPixel);
      nResidualYPixel2by2[chi2Max[i]][iSector] = new TH2F(
          Form("nResidualYPixel2by2_%.1f_%d", chi2Max[i], iSector),
          "Number of tracks used for residual Y 2by2 plot;X (mm);Y (mm)",
          static_cast<int>(xSize / xPixel * 2000), 0, 2 * xSize / xPixel,
          static_cast<int>(ySize / yPixel * 2000), 0, 2 * ySize / yPixel);
    }
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
    efficiencyPixelHisto[iSector] =
        new TH2F(Form("efficiencyPixelHisto_%d", iSector),
                 Form("Efficiency as function of place in the pixel of sector "
                      "%d;X (mm);Y (mm)",
                      iSector),
                 static_cast<int>(xSize / xPixel * 1000), 0, xSize / xPixel,
                 static_cast<int>(ySize / yPixel * 1000), 0, ySize / yPixel);
    efficiencyPixel2by2Histo[iSector] =
        new TH2F(Form("efficiencyPixel2by2Histo_%d", iSector),
                 Form("Efficiency as function of place in four pixels of "
                      "sector %d;X (mm);Y (mm)",
                      iSector),
                 static_cast<int>(xSize / xPixel * 2000), 0, 2 * xSize / xPixel,
                 static_cast<int>(ySize / yPixel * 2000), 0, 2 * ySize / yPixel);
    tracksPAlpidePixelHisto[iSector] =
        new TH2F(Form("tracksPAlpidePixelHisto_%d", iSector),
                 Form("Number of found tracks as function of place in the "
                      "pixel of sector %d;X (mm);Y (mm)",
                      iSector),
                 static_cast<int>(xSize / xPixel * 1000), 0, xSize / xPixel,
                 static_cast<int>(ySize / yPixel * 1000), 0, ySize / yPixel);
    tracksPAlpidePixel2by2Histo[iSector] =
        new TH2F(Form("tracksPAlpidePixel2by2Histo_%d", iSector),
                 Form("Number of found tracks as function of place in four "
                      "pixels of sector %d;X (mm);Y (mm)",
                      iSector),
                 static_cast<int>(xSize / xPixel * 2000), 0, 2 * xSize / xPixel,
                 static_cast<int>(ySize / yPixel * 2000), 0, 2 * ySize / yPixel);
    tracksPixelHisto[iSector] =
        new TH2F(Form("tracksPixelHisto_%d", iSector),
                 Form("Number of tracks as function of place in the pixel of "
                      "sector %d;X (mm);Y (mm)",
                      iSector),
                 static_cast<int>(xSize / xPixel * 1000), 0, xSize / xPixel,
                 static_cast<int>(ySize / yPixel * 1000), 0, ySize / yPixel);
    tracksPixel2by2Histo[iSector] =
        new TH2F(Form("tracksPixel2by2Histo_%d", iSector),
                 Form("Number of tracks as function of place in four pixels of "
                      "sector %d;X (mm);Y (mm)",
                      iSector),
                 static_cast<int>(xSize / xPixel * 2000), 0, 2 * xSize / xPixel,
                 static_cast<int>(ySize / yPixel * 2000), 0, 2 * ySize / yPixel);
    clusterSize2DHisto[iSector] =
        new TH2F(Form("clusterSize2D_%d", iSector),
                 Form("Added cluster size vs the track position of sector %d;X "
                      "(mm);Y (mm)",
                      iSector),
                 static_cast<int>(xSize / xPixel * 1000), 0, xSize / xPixel,
                 static_cast<int>(ySize / yPixel * 1000), 0, ySize / yPixel);
    clusterSize2D2by2Histo[iSector] =
        new TH2F(Form("clusterSize2D2by2_%d", iSector),
                 Form("Added cluster size vs the track position in four pixels "
                      "of sector %d;X (mm);Y (mm)",
                      iSector),
                 static_cast<int>(xSize / xPixel * 2000), 0, 2 * xSize / xPixel,
                 static_cast<int>(ySize / yPixel * 2000), 0, 2 * ySize / yPixel);
    clusterSize2DAverageHisto[iSector] =
        new TH2F(Form("clusterSize2DAverage_%d", iSector),
                 Form("Average cluster size vs the track position of sector "
                      "%d;X (mm);Y (mm)",
                      iSector),
                 static_cast<int>(xSize / xPixel * 1000), 0, xSize / xPixel,
                 static_cast<int>(ySize / yPixel * 1000), 0, ySize / yPixel);
    clusterSize2DAverage2by2Histo[iSector] =
        new TH2F(Form("clusterSize2DAverage2by2_%d", iSector),
                 Form("Average cluster size vs the track position in four "
                      "pixels of sector %d;X (mm);Y (mm)",
                      iSector),
                 static_cast<int>(xSize / xPixel * 2000), 0, 2 * xSize / xPixel,
                 static_cast<int>(ySize / yPixel * 2000), 0, 2 * ySize / yPixel);
    nClusterSizeHisto[iSector] = new TH2F(
        Form("nClusterSizeHisto_%d", iSector),
        "Number of tracks used for cluster size analysis;X (mm);Y (mm)",
        static_cast<int>(xSize / xPixel * 1000), 0, xSize / xPixel,
        static_cast<int>(ySize / yPixel * 1000), 0, ySize / yPixel);
    nClusterSize2by2Histo[iSector] = new TH2F(
        Form("nClusterSize2by2Histo_%d", iSector),
        "Number of tracks used for cluster size analysis 2 by 2;X (mm);Y (mm)",
        static_cast<int>(xSize / xPixel * 2000), 0, 2 * xSize / xPixel,
        static_cast<int>(ySize / yPixel * 2000), 0, 2 * ySize / yPixel);
    clusterWidthXVsXHisto[iSector] =
        new TH1F(Form("clusterWidthXVsXHisto_%d", iSector),
                 Form("Added cluster width in X vs the track position in X of "
                      "sector %d;X (mm);a.u.",
                      iSector),
                 static_cast<int>(xSize / xPixel * 1000), 0, xSize / xPixel);
    clusterWidthYVsYHisto[iSector] =
        new TH1F(Form("clusterWidthYVsYHisto_%d", iSector),
                 Form("Added cluster width in Y vs the track position in Y of "
                      "sector %d;Y (mm);a.u",
                      iSector),
                 static_cast<int>(ySize / yPixel * 1000), 0, ySize / yPixel);
    clusterWidthXVsXAverageHisto[iSector] =
        new TH1F(Form("clusterWidthXVsXAverageHisto_%d", iSector),
                 Form("Average cluster width in X vs the track position in X "
                      "of sector %d;X (mm);a.u.",
                      iSector),
                 static_cast<int>(xSize / xPixel * 1000), 0, xSize / xPixel);
    clusterWidthYVsYAverageHisto[iSector] =
        new TH1F(Form("clusterWidthYVsYAverageHisto_%d", iSector),
                 Form("Average cluster width in Y vs the track position in Y "
                      "of sector %d;Y (mm);a.u",
                      iSector),
                 static_cast<int>(ySize / yPixel * 1000), 0, ySize / yPixel);
    nClusterVsXHisto[iSector] = new TH1F(
        Form("nClusterVsXHisto_%d", iSector),
        "Number of tracks used for the cluster width analysis in X;X (mm);a.u.",
        static_cast<int>(xSize / xPixel * 1000), 0, xSize / xPixel);
    nClusterVsYHisto[iSector] = new TH1F(
        Form("nClusterVsYHisto_%d", iSector),
        "Number of tracks used for the cluster width analysis in Y;Y (mm);a.u.",
        static_cast<int>(xSize / xPixel * 1000), 0, ySize / yPixel);
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
  AIDAProcessor::tree(this)->mkdir("ClusterShape");
  AIDAProcessor::tree(this)->cd("ClusterShape");
  clusterShapeHisto =
      new TH1I("clusterShapeHisto",
               "Cluster shape (all rotations separately);Cluster shape ID;a.u.",
               clusterVec.size() + 1, -0.5, clusterVec.size() + 0.5);
  clusterShapeHistoGrouped = new TH1I(
      "clusterShapeHistoGrouped",
      "Cluster shape (all rotations treated together);Cluster shape ID;a.u.",
      symmetryGroups.size(), -0.5, symmetryGroups.size() - 0.5);
  differenceX =
      new TH1F("differenceX", "Difference in the number of clusters which are "
                              "mirorred around x;Cluster shape ID;a.u.",
               xPairs.size(), -0.5, xPairs.size() - 0.5);
  differenceY =
      new TH1F("differenceY", "Difference in the number of clusters which are "
                              "mirorred around y;Cluster shape ID;a.u.",
               yPairs.size(), -0.5, yPairs.size() - 0.5);
  for (unsigned int i = 0; i < clusterVec.size(); i++) {
    clusterShapeX[i] = new TH1I(
        Form("clusterShapeX_%d", i),
        Form(
            "Place of cluster shape with ID %d as function of X;X (pixel);a.u.",
            i),
        xPixel, 0, xPixel);
    clusterShapeY[i] = new TH1I(
        Form("clusterShapeY_%d", i),
        Form(
            "Place of cluster shape with ID %d as function of Y;Y (pixel);a.u.",
            i),
        yPixel, 0, yPixel);
    clusterShape2D2by2[i] =
        new TH2I(Form("clusterShape2D2by2_%d", i),
                 Form("Distribuition of clusters with shape ID %d within four "
                      "pixels;X (mm);Y (mm)",
                      i),
                 static_cast<int>(xSize / xPixel * 2000), 0, 2 * xSize / xPixel,
                 static_cast<int>(ySize / yPixel * 2000), 0, 2 * ySize / yPixel);
  }
  int tmp = 0;
  for (map<int, int>::iterator it = xPairs.begin(); it != xPairs.end(); ++it) {
    clusterShapeDiffX.insert(make_pair(
        tmp,
        new TH1I(Form("clusterShapeDiffX_%d_%d", it->first, xPairs[it->first]),
                 Form("Difference of the distributions between %d and %d;X "
                      "(pixel);a.u.",
                      it->first, xPairs[it->first]),
                 xPixel, 0, xPixel)));
    clusterShapeDiffY.insert(make_pair(
        tmp,
        new TH1I(Form("clusterShapeDiffY_%d_%d", it->first, xPairs[it->first]),
                 Form("Difference of the distributions between %d and %d;Y "
                      "(pixel);a.u.",
                      it->first, xPairs[it->first]),
                 yPixel, 0, yPixel)));
    tmp++;
  }
  for (map<int, int>::iterator it = yPairs.begin(); it != yPairs.end(); ++it) {
    clusterShapeDiffX.insert(make_pair(
        tmp,
        new TH1I(Form("clusterShapeDiffX_%d_%d", it->first, yPairs[it->first]),
                 Form("Difference of the distributions between %d and %d;X "
                      "(pixel);a.u.",
                      it->first, yPairs[it->first]),
                 xPixel, 0, xPixel)));
    clusterShapeDiffY.insert(make_pair(
        tmp,
        new TH1I(Form("clusterShapeDiffY_%d_%d", it->first, yPairs[it->first]),
                 Form("Difference of the distributions between %d and %d;Y "
                      "(pixel);a.u.",
                      it->first, yPairs[it->first]),
                 yPixel, 0, yPixel)));
    tmp++;
  }
  for (unsigned int i = 0; i < symmetryGroups.size(); i++) {
    string binName;
    for (size_t j = 0; j < symmetryGroups[i].size(); j++) {
      if (j < symmetryGroups[i].size() - 1)
        binName += Form("%d,", symmetryGroups[i][j]);
      else
        binName += Form("%d", symmetryGroups[i][j]);
    }
    string title = "Distribuition of clusters with IDs " + binName +
                   " within four pixels;X (mm);Y (mm)";
    clusterShape2DGrouped2by2[i] =
        new TH2I(Form("clusterShapeGrouped2D2by2_%d", i), title.c_str(),
                 static_cast<int>(xSize / xPixel * 2000), 0, 2 * xSize / xPixel,
                 static_cast<int>(ySize / yPixel * 2000), 0, 2 * ySize / yPixel);
  }
  streamlog_out(DEBUG5) << "end of Booking histograms " << endl;
}
#endif

void EUTelProcessorAnalysisPALPIDEfs::end() {
#ifdef MARLIN_USE_AIDA
  AIDAProcessor::tree(this)->cd("Analysis");
  nHitsPerEventHistoTime = new TH1I(
      "nHitsPerEventHistoTime", "The amount of hits as a function of time",
      nHitsPerEvent.size(), 0, nHitsPerEvent.size());
  nHitsPerEventHisto =
      new TH1I("nHitsPerEventHisto",
               "The amount of hits as a function of events", 200, 0, 200);
  for (size_t i = 0; i < nHitsPerEvent.size(); i++) {
    nHitsPerEventHistoTime->SetBinContent(i, nHitsPerEvent.at(i));
    nHitsPerEventHisto->Fill(nHitsPerEvent.at(i));
  }

  for (size_t event = 0; event < posFit.size(); event++) {
    for (size_t ifit = 0; ifit < posFit.at(event).size(); ifit++) {
      int nAssociatedTracks = 0;
      for (int unsigned i = ifit + 1; i < posFit.at(event).size(); i++) {
        if (abs(posFit.at(event).at(ifit).at(0) -
                posFit.at(event).at(i).at(0)) < limit &&
            abs(posFit.at(event).at(ifit).at(1) -
                posFit.at(event).at(i).at(1)) < limit) {
          nAssociatedTracks++;
        }
      }
      nAssociatedtracksHisto->Fill(nAssociatedTracks);
    }
  }
#endif

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  efficiencyHisto->Divide(tracksHisto);
  for (int i = 0; i < xPixel; i++)
    for (int j = 0; j < yPixel; j++)
      if (tracksPAlpideHisto->GetBinContent(i, j) >
          tracksHisto->GetBinContent(i, j))
        cerr << "More hits in the pAlpide than number of traks" << endl;
  for (int iSector = 0; iSector < _nSectors; iSector++) {
    efficiencyPixelHisto[iSector]->Divide(tracksPixelHisto[iSector]);
    efficiencyPixel2by2Histo[iSector]->Divide(tracksPixel2by2Histo[iSector]);
    clusterSize2DAverageHisto[iSector]->Divide(nClusterSizeHisto[iSector]);
    clusterSize2DAverage2by2Histo[iSector]->Divide(
        nClusterSize2by2Histo[iSector]);
    clusterWidthXVsXAverageHisto[iSector]->Divide(nClusterVsXHisto[iSector]);
    clusterWidthYVsYAverageHisto[iSector]->Divide(nClusterVsYHisto[iSector]);
    AIDAProcessor::tree(this)->cd(Form("Sector_%d", iSector));
    CrossSection *csClusterSize =
        new CrossSection(clusterSize2DAverageHisto[iSector]);
    clusterSizeCrossSection[iSector] = csClusterSize->GetCrossSections();
    CrossSection *csEfficiency =
        new CrossSection(efficiencyPixelHisto[iSector]);
    efficiencyPixelCrossSection[iSector] = csEfficiency->GetCrossSections();
    for (size_t i = 0; i < chi2Max.size(); i++) {
      residualXAveragePixel[chi2Max[i]][iSector]->Divide(
          nResidualXPixel[chi2Max[i]][iSector]);
      residualYAveragePixel[chi2Max[i]][iSector]->Divide(
          nResidualYPixel[chi2Max[i]][iSector]);
      residualXAveragePixel2by2[chi2Max[i]][iSector]->Divide(
          nResidualXPixel2by2[chi2Max[i]][iSector]);
      residualYAveragePixel2by2[chi2Max[i]][iSector]->Divide(
          nResidualYPixel2by2[chi2Max[i]][iSector]);
    }
    nFakeHisto->SetBinContent(iSector + 1, static_cast<double>(nFake[iSector]) / _nEvents /
                                               (xPixel / _nSectors * yPixel));
    nFakeHisto->SetBinError(iSector + 1, sqrt(static_cast<double>(nFake[iSector])) /
                                             _nEvents /
                                             (xPixel / _nSectors * yPixel));
    nFakeWithTrackHisto->SetBinContent(
        iSector + 1, static_cast<double>(nFakeWithTrack[iSector]) / _nEventsWithTrack /
                         (xPixel / _nSectors * yPixel));
    nFakeWithTrackHisto->SetBinError(
        iSector + 1, sqrt(static_cast<double>(nFakeWithTrack[iSector])) / _nEventsWithTrack /
                         (xPixel / _nSectors * yPixel));
    nFakeWithTrackCorrectedHisto->SetBinContent(
        iSector + 1, static_cast<double>(nFakeWithTrackCorrected[iSector]) /
                         _nEventsWithTrack / (xPixel / _nSectors * yPixel));
    nFakeWithTrackCorrectedHisto->SetBinError(
        iSector + 1, sqrt(static_cast<double>(nFakeWithTrackCorrected[iSector])) /
                         _nEventsWithTrack / (xPixel / _nSectors * yPixel));
    nFakeWithoutTrackHisto->SetBinContent(iSector + 1,
                                          static_cast<double>(nFakeWithoutTrack[iSector]) /
                                              (_nEvents - _nEventsWithTrack) /
                                              (xPixel / _nSectors * yPixel));
    nFakeWithoutTrackHisto->SetBinError(
        iSector + 1, sqrt(static_cast<double>(nFakeWithoutTrack[iSector])) /
                         (_nEvents - _nEventsWithTrack) /
                         (xPixel / _nSectors * yPixel));
    tracksProjection[iSector] = tracksHisto->ProjectionY(
        Form("tracksProjection_%d", iSector), xPixel / _nSectors * iSector,
        xPixel / _nSectors * (iSector + 1));
    tracksPAlpideProjection[iSector] = tracksPAlpideHisto->ProjectionY(
        Form("tracksPAlpideProjection_%d", iSector),
        xPixel / _nSectors * iSector, xPixel / _nSectors * (iSector + 1));
    efficiencyProjection[iSector] =
        static_cast<TH1D*>(tracksPAlpideProjection[iSector]->Clone(
            Form("efficiencyProjection_%d", iSector)));
    efficiencyProjection[iSector]->Divide(tracksProjection[iSector]);
    efficiencyProjection[iSector]->SetTitle("Efficiency projection in Y");
  }
  streamlog_out(MESSAGE4) << "nEvents: " << _nEvents << endl;
  streamlog_out(MESSAGE4) << "nEvents with tracks: " << _nEventsWithTrack
                          << endl;
  streamlog_out(MESSAGE4) << "nEvents without tracks: "
                          << _nEvents - _nEventsWithTrack << endl;
  settingsFile << _nEvents << ";";

  hitmapWrongHitHisto->Add(tmpHist, -1.);
  if (_writeShapes) {
    ofstream shapeOutputFile;
    shapeOutputFile.open(_shapeOutputFileName.c_str(), ios::out);
    shapeOutputFile << "ID|Number of clusters with that shape" << endl;
    for (int i = 1; i <= clusterShapeHisto->GetNbinsX(); i++)
      shapeOutputFile << i - 1 << "|" << clusterShapeHisto->GetBinContent(i)
                      << endl;
    shapeOutputFile.close();
  }
  int tmp = 1;
  for (map<int, int>::iterator it = xPairs.begin(); it != xPairs.end(); ++it) {
    differenceX->SetBinContent(
        tmp,
        TMath::Abs(clusterShapeHisto->GetBinContent(it->first + 1) -
                   clusterShapeHisto->GetBinContent(xPairs[it->first] + 1)));
    differenceX->GetXaxis()->SetBinLabel(
        tmp, Form("%d-%d", it->first, xPairs[it->first]));
    tmp++;
  }
  tmp = 1;
  for (map<int, int>::iterator it = yPairs.begin(); it != yPairs.end(); ++it) {
    differenceY->SetBinContent(
        tmp,
        TMath::Abs(clusterShapeHisto->GetBinContent(it->first + 1) -
                   clusterShapeHisto->GetBinContent(yPairs[it->first] + 1)));
    differenceY->GetXaxis()->SetBinLabel(
        tmp, Form("%d-%d", it->first, yPairs[it->first]));
    tmp++;
  }
  for (size_t i = 0; i < symmetryGroups.size(); i++) {
    string binName;
    for (size_t j = 0; j < symmetryGroups[i].size(); j++) {
      clusterShapeHistoGrouped->Fill(
          i, clusterShapeHisto->GetBinContent(symmetryGroups[i][j] + 1));
      if (j < symmetryGroups[i].size() - 1)
        binName += Form("%d,", symmetryGroups[i][j]);
      else
        binName += Form("%d", symmetryGroups[i][j]);
    }
    clusterShapeHistoGrouped->GetXaxis()->SetBinLabel(i + 1, binName.c_str());
  }
  for (int iSector = 0; iSector < _nSectors; iSector++) {
    for (unsigned int i = 0; i < symmetryGroups.size(); i++) {
      string binName;
      for (unsigned int j = 0; j < symmetryGroups[i].size(); j++) {
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
  tmp = 0;
  for (map<int, int>::iterator it = xPairs.begin(); it != xPairs.end(); ++it) {
    for (int i = 1; i <= clusterShapeDiffX[tmp]->GetNbinsX(); i++)
      clusterShapeDiffX[tmp]->SetBinContent(
          i, TMath::Abs(clusterShapeX[it->first]->GetBinContent(i) -
                        clusterShapeX[xPairs[it->first]]->GetBinContent(i)));
    for (int i = 1; i <= clusterShapeDiffY[tmp]->GetNbinsX(); i++)
      clusterShapeDiffY[tmp]->SetBinContent(
          i, TMath::Abs(clusterShapeY[it->first]->GetBinContent(i) -
                        clusterShapeY[xPairs[it->first]]->GetBinContent(i)));
    tmp++;
  }
  for (map<int, int>::iterator it = yPairs.begin(); it != yPairs.end(); ++it) {
    for (int i = 1; i <= clusterShapeDiffX[tmp]->GetNbinsX(); i++)
      clusterShapeDiffX[tmp]->SetBinContent(
          i, TMath::Abs(clusterShapeX[it->first]->GetBinContent(i) -
                        clusterShapeX[yPairs[it->first]]->GetBinContent(i)));
    for (int i = 1; i <= clusterShapeDiffY[tmp]->GetNbinsX(); i++)
      clusterShapeDiffY[tmp]->SetBinContent(
          i, TMath::Abs(clusterShapeY[it->first]->GetBinContent(i) -
                        clusterShapeY[yPairs[it->first]]->GetBinContent(i)));
    tmp++;
  }
/*  TH1I* nX = new TH1I("nX","n",xPixel,0,xPixel);
    TH1I* nY = new TH1I("nY","n",yPixel,0,yPixel);
    for (size_t i=0; i<symmetryGroups.size(); i++)
    {
    nX->Reset();
    nY->Reset();
    string binName;
    double p = 1.0/symmetryGroups[i].size();
    for (size_t j=0; j<symmetryGroups[i].size(); j++)
    {
    if (j<symmetryGroups[i].size()-1) binName +=
   Form("%d-",symmetryGroups[i][j]);
    else binName += Form("%d",symmetryGroups[i][j]);
    for (int k=1; k<=clusterShapeX[j]->GetNbinsX(); k++)
    nX->Fill(k-1,clusterShapeX[symmetryGroups[i][j]]->GetBinContent(k));
    for (int k=1; k<=clusterShapeY[j]->GetNbinsX(); k++)
    nY->Fill(k-1,clusterShapeY[symmetryGroups[i][j]]->GetBinContent(k));
    }
    chi2X.insert(make_pair(i, new
   TH1I(Form("chi2X_%s",(char*)binName.c_str()),"",xPixel,0,xPixel)));
    chi2Y.insert(make_pair(i, new
   TH1I(Form("chi2Y_%s",(char*)binName.c_str()),"",yPixel,0,yPixel)));
    for (size_t j=0; j<symmetryGroups[i].size(); j++)
    {
    for (int k=1; k<=clusterShapeX[j]->GetNbinsX(); k++)
    {
    if (nX->GetBinContent(k) == 0) continue;
    chi2X[i]->Fill(k-1,(clusterShapeX[symmetryGroups[i][j]]->GetBinContent(k)-nX->GetBinContent(k)*p)*(clusterShapeX[symmetryGroups[i][j]]->GetBinContent(k)-nX->GetBinContent(k)*p)/nX->GetBinContent(k)/p);
    }
    for (int k=1; k<=clusterShapeY[j]->GetNbinsX(); k++)
    {
    if (nY->GetBinContent(k) == 0) continue;
    chi2Y[i]->Fill(k-1,(clusterShapeY[symmetryGroups[i][j]]->GetBinContent(k)-nY->GetBinContent(k)*p)*(clusterShapeY[symmetryGroups[i][j]]->GetBinContent(k)-nY->GetBinContent(k)*p)/nY->GetBinContent(k)/p);
    }
    }
    }
*/

#endif
  streamlog_out(MESSAGE4)
      << nPlanesWithTooManyHits
      << " tracks had too many planes with more than one hit" << endl;
  streamlog_out(MESSAGE4) << nNoPAlpideHit
                          << " tracks didn't have a hit in the pALPIDEfs"
                          << endl;
  streamlog_out(MESSAGE4)
      << "For " << nWrongPAlpideHit
      << " tracks the pALPIDEfs had a hit far from the track" << endl;
  streamlog_out(MESSAGE4) << nDUThits
                          << " hits in the DUT weren't associated to a track"
                          << endl;
  streamlog_out(MESSAGE4) << "Overall efficiency of pALPIDEfs sectors: "
                          << endl;
  for (int iSector = 0; iSector < _nSectors; iSector++) {
    streamlog_out(MESSAGE4)
        << static_cast<double>(nTracksPAlpide[iSector]) / nTracks[iSector] << "\t"
        << nTracks[iSector] << "\t" << nTracksPAlpide[iSector] << endl;
    settingsFile << static_cast<double>(nTracksPAlpide[iSector]) / nTracks[iSector] << ";"
                 << nTracks[iSector] << ";" << nTracksPAlpide[iSector] << ";";
  }
  if (_showFake) {
    streamlog_out(MESSAGE4) << "Overall fake efficiency of pALPIDEfs sectors: "
                            << endl;
    for (int iSector = 0; iSector < _nSectors; iSector++) {
      streamlog_out(MESSAGE4)
          << static_cast<double>(nTracksPAlpideFake[iSector]) / nTracksFake[iSector] << "\t"
          << nTracksFake[iSector] << "\t" << nTracksPAlpideFake[iSector]
          << endl;
      settingsFile << static_cast<double>(nTracksPAlpideFake[iSector]) / nTracksFake[iSector]
                   << ";" << nTracksFake[iSector] << ";"
                   << nTracksPAlpideFake[iSector] << ";";
    }
  }
  if (_realAssociation) // Different version of track to hit association written
                        // by Martijn Dietze. It maximizes the number of
                        // association while doesn't allow hits to be shared
                        // between tracks.
  {
    streamlog_out(MESSAGE4) << "Overall efficiency of pALPIDEfs sectors "
                               "without tracks sharing a hit: "
                            << endl;
    for (int iSector = 0; iSector < _nSectors; iSector++) {
      streamlog_out(MESSAGE4) << static_cast<double>(nTracksPAlpideAssociation[iSector]) /
                                     nTracksAssociation[iSector]
                              << "\t" << nTracksAssociation[iSector] << "\t"
                              << nTracksPAlpideAssociation[iSector] << endl;
      settingsFile << static_cast<double>(nTracksPAlpideAssociation[iSector]) /
                          nTracksAssociation[iSector]
                   << ";" << nTracksAssociation[iSector] << ";"
                   << nTracksPAlpideAssociation[iSector] << ";";
    }
  }
  settingsFile << endl;
  streamlog_out(MESSAGE4) << "Successfully finished" << endl;
}

void EUTelProcessorAnalysisPALPIDEfs::_EulerRotationBack(double *_telPos,
                                                         double *_gRotation) {
  double z = _telPos[2] - dutZ;
  TVector3 _RotatedSensorHit(_telPos[0], _telPos[1], z);
  TVector3 _Xaxis(1.0, 0.0, 0.0);
  TVector3 _Yaxis(0.0, 1.0, 0.0);
  TVector3 _Zaxis(0.0, 0.0, 1.0);
  // rotation order: X,Y,Z; rotation back order: Z,Y,X
  if (TMath::Abs(_gRotation[1]) > 1e-6) {
    _RotatedSensorHit.Rotate(-1. * _gRotation[1], _Yaxis); // in ZX
  }
  if (TMath::Abs(_gRotation[2]) > 1e-6) {
    _RotatedSensorHit.Rotate(-1. * _gRotation[2], _Xaxis); // in ZY
  }
  if (TMath::Abs(_gRotation[0]) > 1e-6) {
    _RotatedSensorHit.Rotate(-1. * _gRotation[0], _Zaxis); // in XY
  }

  _telPos[0] = _RotatedSensorHit.X();
  _telPos[1] = _RotatedSensorHit.Y();
  _telPos[2] = _RotatedSensorHit.Z() + dutZ;
}

void EUTelProcessorAnalysisPALPIDEfs::_LayerRotationBack(double *pos,
                                                         double &outputX,
                                                         double &outputY) {
  double sign = 0;
  if (xPointing[0] < -0.7)
    sign = -1;
  else if (xPointing[0] > 0.7)
    sign = 1;
  else {
    if (xPointing[1] < -0.7)
      sign = -1;
    else if (xPointing[1] > 0.7)
      sign = 1;
  }
  pos[0] -= (-1) * sign * xSize / 2;
  if (yPointing[0] < -0.7)
    sign = -1;
  else if (yPointing[0] > 0.7)
    sign = 1;
  else {
    if (yPointing[1] < -0.7)
      sign = -1;
    else if (yPointing[1] > 0.7)
      sign = 1;
  }
  pos[1] -= (-1) * sign * ySize / 2;

  outputY = (xPointing[0] * pos[1] - yPointing[0] * pos[0]) /
            (yPointing[1] * xPointing[0] - yPointing[0] * xPointing[1]);
  outputX = pos[0] / xPointing[0] - xPointing[1] / xPointing[0] * outputY;
}

int EUTelProcessorAnalysisPALPIDEfs::AddressToColumn(int ARegion,
                                                     int ADoubleCol,
                                                     int AAddress) {
  int Column =
      ARegion * 32 + ADoubleCol * 2; // Double columns before ADoubleCol
  int LeftRight = 0;
  if (_chipVersion >= 3) {
    LeftRight = ((((AAddress % 4) == 1) || ((AAddress % 4) == 2))
                     ? 1
                     : 0); // Left or right column within the double column
  } else {
    LeftRight = ((AAddress % 4) < 2
                     ? 1
                     : 0); // Left or right column within the double column
  }
  Column += LeftRight;
  return Column;
}

int EUTelProcessorAnalysisPALPIDEfs::AddressToRow(int AAddress) {
  // Ok, this will get ugly
  int Row = AAddress / 2; // This is OK for the top-right and the bottom-left
                          // pixel within a group of 4
  if (_chipVersion < 3) {
    if ((AAddress % 4) == 3)
      Row -= 1; // adjust the top-left pixel
    if ((AAddress % 4) == 0)
      Row += 1; // adjust the bottom-right pixel
  }
  return Row;
}

bool EUTelProcessorAnalysisPALPIDEfs::emptyMiddle(
    vector<vector<int>> pixVector) {

  // In this step we use an invers clustering process:
  // - First, we create a rectangle, which is bigger than the cluster in every
  // direction.
  // - Second, we fill the hitPixelVec with pixels in this rectangle, which did
  // not fired.
  // - Third, we run the clustering process with this hitPixelVec.
  // - Fourth, when this process finds one whole cluster and hitPixelVec is
  // still not empty, the original cluster has to have an empty middle, since at
  // least two cluster could be made from the empty part (one on the outside and
  // at least one in the middle of the cluster).

  std::vector<EUTelGenericSparsePixel> hitPixelVec;

  std::vector<EUTelGenericSparsePixel> newlyAdded;

  int xMax = 0, yMax = 0, xMin = 1000000, yMin = 1000000;

  for (size_t i = 0; i < pixVector.size(); i++) {
    if (pixVector[i][0] > xMax)
      xMax = pixVector[i][0];
    if (pixVector[i][0] < xMin)
      xMin = pixVector[i][0];
    if (pixVector[i][1] > yMax)
      yMax = pixVector[i][1];
    if (pixVector[i][1] < yMin)
      yMin = pixVector[i][1];
  }
  for (int n = xMin - 1; n <= xMax + 1; n++) {
    for (int m = yMin - 1; m <= yMax + 1; m++) {
      bool empty_pixel = true;
      for (size_t i = 0; i < pixVector.size(); i++) {
        if (n == pixVector[i][0] && m == pixVector[i][1]) {
          empty_pixel = false;
          break;
        }
      }
      if (empty_pixel) {
        EUTelGenericSparsePixel pixel;
        pixel.setXCoord(n);
        pixel.setYCoord(m);
        hitPixelVec.push_back(pixel);
      }
    }
  }

  // We now cluster those hits together
  // while( !hitPixelVec.empty() )
  {

    std::vector<EUTelGenericSparsePixel> cluCandidate;

    // First we need to take any pixel, so let's take the first one
    // Add it to the cluster as well as the newly added pixels
    newlyAdded.push_back(hitPixelVec.front());
    // sparseCluster->push_back( &(hitPixelVec.front()) );
    cluCandidate.push_back(hitPixelVec.front());
    // And remove it from the original collection
    hitPixelVec.erase(hitPixelVec.begin());

    // Now process all newly added pixels, initially this is the just previously
    // added one
    // but in the process of neighbour finding we continue to add new pixels
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

        int _sparseMinDistanceSquaredComparison = 1;
        if (distance <= _sparseMinDistanceSquaredComparison) {
          // add them to the cluster as well as to the newly added ones
          newlyAdded.push_back(*hitVec);
          cluCandidate.push_back(*hitVec);
          //	sparseCluster->push_back( &(*hitVec) );
          hitPixelVec.erase(hitVec);
          // for the pixel we test there might be other neighbours, we still
          // have to check
          newlyDone = false;
          break;
        }
      }
      // if no neighbours are found, we can delete the pixel from the newly
      // added
      // we tested against _ALL_ non cluster pixels, there are no other pixels
      // which could be neighbours
      if (newlyDone) {
        newlyAdded.erase(newlyAdded.begin());
      }
    }
  }
  // If the hitPixelVec is not empty there is a empty middle cluster, because
  // there must be an area inside the cluster, which is not connected to the
  // area outside the cluster.
  if (!hitPixelVec.empty())
    return true;
  else
    return false;
}

bool EUTelProcessorAnalysisPALPIDEfs::RemoveAlign(
    LCCollectionVec *preAlignmentCollectionVec,
    LCCollectionVec *alignmentCollectionVec,
    LCCollectionVec *alignmentPAlpideCollectionVec, double *fitpos,
    double &xposfit, double &yposfit) {
  // Remove the alignment in the same way it was applied by the
  // EUTelProcessorApplyAlignment.cc
  double xPlaneCenter = geo::gGeometry().getPlaneXPosition(_dutID);
  double yPlaneCenter = geo::gGeometry().getPlaneYPosition(_dutID);
  double zPlaneThickness = geo::gGeometry().getPlaneZSize(_dutID);
  double zPlaneCenter =
      geo::gGeometry().getPlaneZPosition(_dutID) + zPlaneThickness / 2.;
  TVector3 inputVec(fitpos[0] - xPlaneCenter, fitpos[1] - yPlaneCenter,
                    fitpos[2] - zPlaneCenter);
  // cerr << zPlaneThickness << "\t" << zPlaneCenter << endl;
  EUTelAlignmentConstant *alignment = nullptr;
  bool alignExist = false;
  EUTelAlignmentConstant *preAlignment = nullptr;
  bool prealignExist = false;
  EUTelAlignmentConstant *alignmentPAlpide = nullptr;
  bool alignPAlpideExist = false;
  if (!_oneAlignmentCollection) {
    for (int iAlign = 0;
         iAlign < alignmentPAlpideCollectionVec->getNumberOfElements();
         iAlign++) {
      alignmentPAlpide = static_cast<EUTelAlignmentConstant *>(
          alignmentPAlpideCollectionVec->getElementAt(iAlign));
      if (alignmentPAlpide->getSensorID() == _dutID) {
        alignPAlpideExist = true;
        inputVec[0] += alignmentPAlpide->getXOffset();
        inputVec[1] += alignmentPAlpide->getYOffset();
        inputVec[2] += alignmentPAlpide->getZOffset();
        inputVec.RotateZ(alignmentPAlpide->getGamma());
        inputVec.RotateY(alignmentPAlpide->getBeta());
        inputVec.RotateX(alignmentPAlpide->getAlpha());
        break;
      }
    }
    if (!alignPAlpideExist && _isFirstEvent)
      cerr << "No second alignment correction applied to the pAlpide!" << endl;
  }
  for (int iAlign = 0; iAlign < alignmentCollectionVec->getNumberOfElements();
       iAlign++) {
    alignment = static_cast<EUTelAlignmentConstant *>(
        alignmentCollectionVec->getElementAt(iAlign));
    if (alignment->getSensorID() == _dutID) {
      alignExist = true;
      inputVec[0] += alignment->getXOffset();
      inputVec[1] += alignment->getYOffset();
      inputVec[2] += alignment->getZOffset();
      inputVec.RotateZ(alignment->getGamma());
      inputVec.RotateY(alignment->getBeta());
      inputVec.RotateX(alignment->getAlpha());
      break;
    }
  }
  if (!alignExist)
    return false; // cerr << "No alignment correction applied!" << endl;
  for (int iAlign = 0;
       iAlign < preAlignmentCollectionVec->getNumberOfElements(); iAlign++) {
    preAlignment = static_cast<EUTelAlignmentConstant *>(
        preAlignmentCollectionVec->getElementAt(iAlign));
    if (preAlignment->getSensorID() == _dutID) {
      prealignExist = true;
      inputVec[0] += preAlignment->getXOffset();
      inputVec[1] += preAlignment->getYOffset();
      inputVec[2] += preAlignment->getZOffset();
      break;
    }
  }
  if (!prealignExist)
    return false; // cerr << "No prealignment correction applied!" << endl;

  fitpos[0] = inputVec.X();
  fitpos[1] = inputVec.Y();
  fitpos[2] = inputVec.Z() + zPlaneCenter; // since _EulerRotationBack removes
                                           // the zPlaneCenter again (and adds
                                           // it at the end)
  // cerr << fitpos[0] << "\t"  << fitpos[1] << "\t" << fitpos[2] << endl;

  _EulerRotationBack(fitpos, gRotation);
  _LayerRotationBack(fitpos, xposfit, yposfit);

  return 1;
}
