#include "EUTelProcessorAnalysisPALPIDEfs.h"
#include "EUTelHistogramManager.h"
#include "EUTelAlignmentConstant.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelTrackerDataInterfacerImpl.h"

#include "marlin/Global.h"
#include "marlin/AIDAProcessor.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ITree.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/IAxis.h>
#endif

#include <EVENT/LCCollection.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCCollectionVec.h>

#include "TVector3.h"
#include "TFitResult.h"

#include <memory>
#include <algorithm>
#include <cmath>

using namespace lcio ;
using namespace marlin ;
using namespace std;
using namespace eutelescope;
using namespace gear;

EUTelProcessorAnalysisPALPIDEfs aEUTelProcessorAnalysisPALPIDEfs;

EUTelProcessorAnalysisPALPIDEfs::EUTelProcessorAnalysisPALPIDEfs()
: Processor("EUTelProcessorAnalysisPALPIDEfs"),
  _fillHistos(false),
  _inputFittedHitName(""),
  _inputColName(""),
  _trackCollectionName(""),
  _alignmentPAlpideCollectionName("alignmentPAlpide"),
  _alignmentCollectionName("alignment"),
  _preAlignmentCollectionName("prealign"),
  _zsDataCollectionName(""),
  zsInputDataCollectionVec(NULL),
  _hotPixelCollectionName(""),
  limit(0.05),
  _dutID(6),
  _maxNumberOfPixels(3),
  _nPlanesWithMoreHits(4),
  _moreTracks(false),
  _energy(6.0),
  _writeShapes(false),
  _shapeOutputFileName("./shapeDistribution.txt"),
  _outputSettingsFolderName("./"),
  _chipID(),
  _irradiation(),
  _rate(""),
  _oneAlignmentCollection(false),
  _clusterAvailable(true),
  _hotpixelAvailable(true),
  _noiseMaskAvailable(true),
  _deadColumnAvailable(true),
  chi2Max(1),
//  chi2Max(8),
  _nEvents(0),
  _nEventsWithTrack(0),
  _minTimeStamp(0),
  nTracks(4),
  nTracksPAlpide(4),
  nFakeWithTrack(4,0),
  nFakeWithoutTrack(4,0),
  nFake(4,0),
  nFakeWithTrackCorrected(4,0),
  nDUThits(0),
  nNoPAlpideHit(0),
  nWrongPAlpideHit(0),
  nPlanesWithTooManyHits(0),
  xZero(0),
  yZero(0),
  xPitch(0),
  ySize(0),
  xPixel(0),
  yPixel(0),
  hotPixelCollectionVec(NULL)
{
  _description="Analysis of the fitted tracks";
  registerInputCollection( LCIO::TRACKERHIT,
                           "InputFittedHitName" ,
                           "Name of the input fitted TrackerHit collection"  ,
                           _inputFittedHitName ,
                           std::string("fithit") ) ;
  registerInputCollection( LCIO::TRACKERHIT,
                           "InputCollectionName" ,
                           "Name of the input TrackerHit collection"  ,
                           _inputColName ,
                           std::string("colhit") ) ;
  registerInputCollection( LCIO::TRACK,
                           "TrackCollectionName",
                           "Input track collection name",
                           _trackCollectionName,
                           std::string("track"));
  registerInputCollection (LCIO::LCGENERICOBJECT, "AlignmentConstantName",
                           "Alignment constant from the condition file",
                           _alignmentCollectionName, string ("alignment"));
  registerInputCollection (LCIO::LCGENERICOBJECT, "AlignmentPAlpideConstantName",
                           "Alignment constant from the condition file",
                           _alignmentPAlpideCollectionName, string ("alignmentPAlpide"));
  registerInputCollection (LCIO::LCGENERICOBJECT, "PreAlignmentConstantName",
                           "PreAlignment constant from the condition file",
                           _preAlignmentCollectionName, string ("prealign"));
  registerInputCollection (LCIO::TRACKERDATA, "ZSDataCollectionName",
                           "Input of Zero Suppressed data",
                           _zsDataCollectionName, string ("zsdata") );
  registerProcessorParameter("HistogramFilling","Switch on or off the histogram filling",
                             _fillHistos, static_cast< bool > ( true ) );
  registerProcessorParameter("HistoInfoFileName", "This is the name of the histogram information file",
                             _histoInfoFileName, string( "histoinfo.xml" ) );
  registerProcessorParameter("Limit", "This is allowed distance between the track and the hit",
                             limit, static_cast<double>( 0.05 ) );
  registerProcessorParameter("dutID", "This is the ID of the DUT",
                             _dutID, static_cast<int>( 6 ) );
  registerProcessorParameter("MaxNumberOfPixels", "This is the maximum number of pixels in one cluster for the clustershape analysis",
                             _maxNumberOfPixels, static_cast<int>( 3 ) );
  registerProcessorParameter("nPlanesWithMoreHits", "This is the maximum number of planes that can have more than one hit",
                             _nPlanesWithMoreHits, static_cast<int>( 4 ) );
  registerProcessorParameter("MoreTracks","More tracks are allowed in one event",
                             _moreTracks, static_cast< bool > ( false ) );
  registerOptionalParameter("HotPixelCollectionName","This is the name of the hotpixel collection of the pALPIDE",
                             _hotPixelCollectionName, static_cast< string > ( "" ) );
  registerOptionalParameter("DeadColumnCollectionName","This is the name of the collection containing the pixels belonging to a dead column",
                             _deadColumnCollectionName, static_cast< string > ( "deadColumn" ) );
  registerOptionalParameter("NoiseMaskFileName","This is the name of the file which contains the pixels which were masked during datataking",
                             _noiseMaskFileName, static_cast< string > ( "" ) );
  registerOptionalParameter("Energy","Particle energy",
                             _energy, static_cast< double > ( 6.0 ) );
  registerProcessorParameter("WriteShapes","Write cluster shapes to file?",
                             _writeShapes, static_cast< bool > ( false ) );
  registerOptionalParameter("ShapeOutputFileName","This is the name of the file where the IDs of the cluster shapes will be saved",
                             _shapeOutputFileName, static_cast< string > ( "./shapeDistribution.txt" ) );
  registerOptionalParameter("OutputSettingsFolderName","Folder name where all the settings of each run will be saved",
                             _outputSettingsFolderName, static_cast< string > ( "./" ) );
  EVENT::StringVec _stringVecExample;
  _stringVecExample.push_back(" ");
  registerOptionalParameter("ChipID","Chip IDs",
                             _chipID, _stringVecExample );
  registerOptionalParameter("Irradiation","Irradiation level",
                             _irradiation, _stringVecExample );
  registerOptionalParameter("Rate","Data taking rate",
                             _rate, static_cast< string > ( "" ) );
  registerProcessorParameter("MinTimeStamp", "This is minimum timestamp required to consider an event",
                             _minTimeStamp, static_cast<double>( 0 ) );

  _isFirstEvent = true;
}

void EUTelProcessorAnalysisPALPIDEfs::init() {
  printParameters();
  int _nTelPlanes = geo::gGeometry().nPlanes();
  const std::vector<int>& _planeID = geo::gGeometry().sensorIDsVec();
  for(int iz=0; iz < _nTelPlanes ; iz++)
    if(_planeID[iz]==_dutID)
    {
      dutZ = geo::gGeometry().siPlaneZPosition(iz);
      layerIndex = iz;
      xSize = geo::gGeometry().siPlaneXSize(layerIndex);
      ySize = geo::gGeometry().siPlaneYSize(layerIndex);
      xZero        = geo::gGeometry().siPlaneXPosition(layerIndex); // mm
      yZero        = geo::gGeometry().siPlaneYPosition(layerIndex); // mm
      xSize        = geo::gGeometry().siPlaneXSize(layerIndex);     // mm
      ySize        = geo::gGeometry().siPlaneYSize(layerIndex);     // mm
      xPitch       = geo::gGeometry().siPlaneXPitch(layerIndex);    // mm
      yPitch       = geo::gGeometry().siPlaneYPitch(layerIndex);    // mm
      xPointing[0] = geo::gGeometry().siPlaneRotation1(layerIndex); // was -1 ;
      xPointing[1] = geo::gGeometry().siPlaneRotation2(layerIndex); // was  0 ;
      yPointing[0] = geo::gGeometry().siPlaneRotation3(layerIndex); // was  0 ;
      yPointing[1] = geo::gGeometry().siPlaneRotation4(layerIndex); // was -1 ;
      xPixel       = geo::gGeometry().siPlaneXNpixels(layerIndex);
      yPixel       = geo::gGeometry().siPlaneYNpixels(layerIndex);
      try
      {
        gRotation[0] = geo::gGeometry().siPlaneZRotation(layerIndex); // Euler gamma ;
        gRotation[1] = geo::gGeometry().siPlaneYRotation(layerIndex); // Euler beta  ;
        gRotation[2] = geo::gGeometry().siPlaneXRotation(layerIndex); // Euler alpha ;
      }
      catch(...)
      {
        streamlog_out ( MESSAGE5 ) << " no sensor rotation is given in the GEAR steering file, assume NONE " << endl;
      }
      if ((gRotation[1] != 0 && gRotation[1] != 180) || (gRotation[2] != 0 && gRotation[2] != 180)) zDistance = sqrt(xSize*xSize+ySize*ySize);
      else zDistance  = 0.1;
      gRotation[0] =  gRotation[0]*3.1415926/180.; //
      gRotation[1] =  gRotation[1]*3.1415926/180.; //
      gRotation[2] =  gRotation[2]*3.1415926/180.; //
    }
//  float chi2MaxTemp[8] = {4,6,8,10,15,20,25,30};
  float chi2MaxTemp[1] = {30};
  for (unsigned int i=0; i<chi2Max.size(); i++)
    chi2Max[i] = chi2MaxTemp[i];
  Cluster cluster;
  cluster.FindReferenceClusters(clusterVec,_maxNumberOfPixels);
  xPairs = cluster.SymmetryPairs(clusterVec,"x");
  yPairs = cluster.SymmetryPairs(clusterVec,"y");
  symmetryGroups = cluster.sameShape(clusterVec);
  for (int iSector=0; iSector<4; iSector++)
  {
    nTracks[iSector] = 0;
    nTracksPAlpide[iSector] = 0;
  }
  bool newFile = false;
  string _outputSettingsFileName = _outputSettingsFolderName + Form("settings_DUT%d",_dutID) + ".txt";
  if (!std::ifstream(_outputSettingsFileName.c_str()))
    newFile = true;
  settingsFile.open (_outputSettingsFileName.c_str(), ios::out | ios::app );
    if (newFile) settingsFile << "Run number;Energy;Chip ID;Irradiation level(0-nonIrradiated,1-2.5e12,2-1e13,3-700krad,4-combined:1e13+700krad);Rate;BB;Ithr;Idb;Vcasn;Vaux;Vcasp;Vreset;Threshold and their RMS for all four sectors;Noise and their RMS for all four sectors;Readout delay;Trigger delay;Strobe length;StrobeB length;Data (1) or noise (0);Number of events;Efficiency,Number of tracks,Number of tracks with associated hit for all sectors" << endl;
}

void EUTelProcessorAnalysisPALPIDEfs::processEvent(LCEvent *evt)
{
  if (evt->getParameters().getIntVal("FLAG") == 100) return; //Excluding events with too large clusters
  int nTrackPerEvent = 0, nClusterAssociatedToTrackPerEvent = 0, nClusterPerEvent = 0;
  if (evt->getTimeStamp() < _minTimeStamp) return;
  if (_isFirstEvent)
  {
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    if ( _fillHistos)
   {
      bookHistos();
    }
#endif
    hotPixelCollectionVec = 0;
    try
    {
      hotPixelCollectionVec = static_cast< LCCollectionVec* >  (evt->getCollection( _hotPixelCollectionName )) ;
      streamlog_out ( DEBUG5 ) << "hotPixelCollectionName: " << _hotPixelCollectionName.c_str() << " found " << endl;
    }
    catch (lcio::DataNotAvailableException& e )
    {
      streamlog_out ( WARNING5 ) << "hotPixelCollectionName: " << _hotPixelCollectionName.c_str() << "not found " << endl;
      _hotpixelAvailable = false;
    }
    if (_hotpixelAvailable)
    {
      hotData = dynamic_cast< TrackerDataImpl * > ( hotPixelCollectionVec->getElementAt( layerIndex ) );
      auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel > >  sparseData(new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> ( hotData ));
      for ( unsigned int iPixel = 0; iPixel < sparseData->size(); iPixel++ )
      {
        EUTelGenericSparsePixel *sparsePixel =  new EUTelGenericSparsePixel() ;
        sparseData->getSparsePixelAt( iPixel, sparsePixel );
        hotpixelHisto->Fill(sparsePixel->getXCoord()*xPitch+xPitch/2.,sparsePixel->getYCoord()*yPitch+yPitch/2.);
        delete sparsePixel;
      }
    }
    ifstream noiseMaskFile(_noiseMaskFileName.c_str());
    if (noiseMaskFile.is_open())
    {
      streamlog_out ( MESSAGE4 ) << "Running with noise mask: " << _noiseMaskFileName.c_str() << endl;
      int region, doubleColumn, address;
      while (noiseMaskFile >> region >> doubleColumn >> address)
      {
        int x = AddressToColumn(region,doubleColumn,address);
        int y = AddressToRow(address);
        noiseMaskX.push_back(x);
        noiseMaskY.push_back(y);
        hotpixelHisto->Fill(x*xPitch+xPitch/2.,y*yPitch+yPitch/2.);
      }
    }
    else _noiseMaskAvailable = false;

    deadColumnCollectionVec = 0;
    try
    {
      deadColumnCollectionVec =  static_cast< LCCollectionVec* >  (evt->getCollection( _deadColumnCollectionName));
    }
    catch (lcio::DataNotAvailableException& e )
    {
      streamlog_out ( WARNING5 ) << "deadPixelCollectionName: " << _deadColumnCollectionName.c_str() << " not found " << endl;
      _deadColumnAvailable = false;
    }
    if (_deadColumnAvailable)
    {
      deadColumn = dynamic_cast< TrackerDataImpl * > ( deadColumnCollectionVec->getElementAt( layerIndex ) );
      auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel > >  sparseData(new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> (deadColumn));
      for ( unsigned int iPixel = 0; iPixel < sparseData->size(); iPixel++ )
      {
        EUTelGenericSparsePixel *sparsePixel =  new EUTelGenericSparsePixel() ;
        sparseData->getSparsePixelAt( iPixel, sparsePixel );
        deadColumnHisto->Fill(sparsePixel->getXCoord()*xPitch+xPitch/2.,sparsePixel->getYCoord()*yPitch+yPitch/2.);
        delete sparsePixel;
      }
    }
    settingsFile << evt->getRunNumber() << ";" << _energy << ";" << _chipID[layerIndex] << ";" << _irradiation[layerIndex] << ";" << _rate << ";" << evt->getParameters().getFloatVal("BackBiasVoltage") << ";" << evt->getParameters().getIntVal(Form("Ithr_%d",layerIndex)) << ";" << evt->getParameters().getIntVal(Form("Idb_%d",layerIndex)) << ";" << evt->getParameters().getIntVal(Form("Vcasn_%d",layerIndex)) << ";" << evt->getParameters().getIntVal(Form("Vaux_%d",layerIndex)) << ";" << evt->getParameters().getIntVal(Form("Vcasp_%d",layerIndex)) << ";" << evt->getParameters().getIntVal(Form("Vreset_%d",layerIndex)) << ";";
    for (int iSector=0; iSector<4; iSector++)
      settingsFile << evt->getParameters().getFloatVal(Form("Thr_%d_%d",layerIndex,iSector)) << ";" << evt->getParameters().getFloatVal(Form("ThrRMS_%d_%d",layerIndex,iSector)) << ";";
    for (int iSector=0; iSector<4; iSector++)
      settingsFile << evt->getParameters().getFloatVal(Form("Noise_%d_%d",layerIndex,iSector)) << ";" << evt->getParameters().getFloatVal(Form("NoiseRMS_%d_%d",layerIndex,iSector)) << ";";
    settingsFile << evt->getParameters().getIntVal(Form("m_readout_delay_%d",layerIndex)) << ";" << evt->getParameters().getIntVal(Form("m_trigger_delay_%d",layerIndex)) << ";" << evt->getParameters().getIntVal(Form("m_strobe_length_%d",layerIndex)) << ";" << evt->getParameters().getIntVal(Form("m_strobeb_length_%d",layerIndex)) << ";1;";
    _isFirstEvent = false;
  }
  timeStampHisto->Fill(evt->getTimeStamp());
  LCCollection* col;
  try
  {
    col = evt->getCollection( _inputColName ) ;
  }
  catch (lcio::DataNotAvailableException& e)
  {
    streamlog_out ( DEBUG5 ) << "Not able to get collection "
                            << _inputColName
                            << "\nfrom event " << evt->getEventNumber()
                            << " in run " << evt->getRunNumber()  << endl;
    return;
  }
  LCCollection* colFit = 0;
  bool fitHitAvailable = true;
  try
  {
    colFit = evt->getCollection( _inputFittedHitName ) ;
  }
  catch (lcio::DataNotAvailableException& e)
  {
    streamlog_out ( DEBUG5 ) << "Not able to get fit collection "
                            << _inputFittedHitName
                            << "\nfrom event " << evt->getEventNumber()
                            << " in run " << evt->getRunNumber()  << endl;
    fitHitAvailable = false;
  }
  LCCollection* colTrack = NULL;
  try
  {
    colTrack = evt->getCollection(_trackCollectionName);
  } catch (DataNotAvailableException e)
  {
    fitHitAvailable = false;
  }
  LCCollectionVec * alignmentPAlpideCollectionVec = 0;
  LCCollectionVec * alignmentCollectionVec = 0;
  LCCollectionVec * preAlignmentCollectionVec = 0;
  try {
    alignmentCollectionVec        = dynamic_cast < LCCollectionVec * > (evt->getCollection(_alignmentCollectionName));
    preAlignmentCollectionVec     = dynamic_cast < LCCollectionVec * > (evt->getCollection(_preAlignmentCollectionName));
  } catch (...) {
    streamlog_out  ( WARNING2 ) << "Alignment collection not available" << endl;
    return;
  }
  try {
    alignmentPAlpideCollectionVec        = dynamic_cast < LCCollectionVec * > (evt->getCollection(_alignmentPAlpideCollectionName));
  } catch (...) {
    if (evt->getEventNumber() == 0) streamlog_out  ( WARNING2 ) << "Only one alignment collection for pAlpide" << endl;
    _oneAlignmentCollection = true;
  }
  _clusterAvailable = true;
  try {
    zsInputDataCollectionVec = dynamic_cast< LCCollectionVec * > ( evt->getCollection( _zsDataCollectionName ) ) ;
    streamlog_out ( DEBUG4 ) << "zsInputDataCollectionVec: " << _zsDataCollectionName.c_str() << " found " << endl;
  } catch ( lcio::DataNotAvailableException ) {
    streamlog_out ( DEBUG4 ) << "zsInputDataCollectionVec: " << _zsDataCollectionName.c_str() << " not found " << endl;
    _clusterAvailable = false;
  }
  if (_clusterAvailable)
  {
    CellIDDecoder<TrackerDataImpl > cellDecoder( zsInputDataCollectionVec );
    for ( unsigned int iCluster=0; iCluster<zsInputDataCollectionVec->size(); iCluster++)
    {
      TrackerDataImpl * zsData = dynamic_cast< TrackerDataImpl * > ( zsInputDataCollectionVec->getElementAt(iCluster) );
      if((int)cellDecoder(zsData)["sensorID"] == _dutID) nClusterPerEvent++;
    }
  }
  if (  ( xPointing[0] == xPointing[1] ) && ( xPointing[0] == 0 ) ) {
    streamlog_out ( ERROR4 ) << "DUT has a singular rotation matrix. Sorry for quitting" << endl;
  }

  if (  ( yPointing[0] == yPointing[1] ) && ( yPointing[0] == 0 ) ) {
    streamlog_out ( ERROR4 ) << "Detector DUT has a singular rotation matrix. Sorry for quitting" << endl;
  }


  int nFitHit = 0;
  if (fitHitAvailable)  nFitHit = colFit->getNumberOfElements();
  bool hitmapFilled = false;
  vector<int> clusterAssosiatedToTrack;
  for(int ifit=0; ifit< nFitHit ; ifit++)
  {
    TrackerHit * fithit = dynamic_cast<TrackerHit*>( colFit->getElementAt(ifit) ) ;
    double fitpos[3]={0.,0.,0.};
    const double *fitpos0 = fithit->getPosition();
    fitpos[0] = fitpos0[0];
    fitpos[1] = fitpos0[1];
    fitpos[2] = fitpos0[2];
    bool twoTracks = false;
    double yposfitPrev = 0.;
    double xposfitPrev = 0.;
    double maxDistInPlane = 0.1;
    for(int i=ifit+1; i< nFitHit ; i++)
    {
      TrackerHit * fithitcheck = dynamic_cast<TrackerHit*>( colFit->getElementAt(i) ) ;
      const double *fitposcheck0 = fithitcheck->getPosition();
      if (abs(fitpos[2] - fitposcheck0[2]) < maxDistInPlane) {twoTracks = true; break;}
    }
    if (twoTracks && !_moreTracks) {streamlog_out ( MESSAGE1 ) << "2 tracks in event " << evt->getEventNumber() << ", skipping the event!" << endl; break;}
    if (fitpos[2] >= dutZ-zDistance && fitpos[2] <= dutZ+zDistance )
    {
      double xposfit=0, yposfit=0;
      bool alignSuccess = RemoveAlign(preAlignmentCollectionVec,alignmentCollectionVec,alignmentPAlpideCollectionVec,fitpos,xposfit,yposfit);
      if (!alignSuccess && _isFirstEvent) cerr << "Removing alignment did not work!" << endl;
      if (xposfit > 0 && yposfit > 0 && xposfit < xSize && yposfit < ySize)
      {
        if (xposfit < limit || xposfit > xSize-limit || yposfit < limit || yposfit > ySize-limit) continue;
        int index = -1;
        for (int iSector=0; iSector<4; iSector++)
        {
          if (xposfit>xSize/4.*iSector+(iSector==0?0:1)*(2.*xPitch+limit) && xposfit<xSize/4.*(iSector+1)-(iSector==3?0:1)*(2.*xPitch+limit))
          {
            index = iSector;
            break;
          }
        }
        if (index == -1) continue;
        if (_hotpixelAvailable)
        {
          auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel > >  sparseData(new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> ( hotData ));
          bool hotpixel = false;
          for ( unsigned int iPixel = 0; iPixel < sparseData->size(); iPixel++ )
          {
            EUTelGenericSparsePixel *sparsePixel =  new EUTelGenericSparsePixel() ;
            sparseData->getSparsePixelAt( iPixel, sparsePixel );
            if (abs(xposfit-(sparsePixel->getXCoord()*xPitch+xPitch/2.)) < limit && abs(yposfit-(sparsePixel->getYCoord()*yPitch+yPitch/2.)) < limit)
            {
              hotpixel = true;
              delete sparsePixel;
              break;
            }
            delete sparsePixel;
          }
          if (hotpixel) continue;
        }
        if (_noiseMaskAvailable)
        {
          bool noisePixel = false;
          for (unsigned int iNoise=0; iNoise<noiseMaskX.size(); iNoise++)
          {
            if (abs(xposfit-(noiseMaskX[iNoise]*xPitch+xPitch/2.)) < limit && abs(yposfit-(noiseMaskY[iNoise]*yPitch+yPitch/2.)) < limit)
            {
              noisePixel = true;
              break;
            }
          }
          if (noisePixel) continue;
        }
        if (_deadColumnAvailable)
        {
          auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel > >  sparseData(new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> ( deadColumn ));
          bool dead = false;
          for ( unsigned int iPixel = 0; iPixel < sparseData->size(); iPixel++ )
          {
            EUTelGenericSparsePixel *sparsePixel =  new EUTelGenericSparsePixel() ;
            sparseData->getSparsePixelAt( iPixel, sparsePixel );
            if (abs(xposfit-(sparsePixel->getXCoord()*xPitch+xPitch/2.)) < limit)
            {
              dead = true;
              delete sparsePixel;
              break;
            }
            delete sparsePixel;
          }
          if (dead) continue;

        }
        nTrackPerEvent++;
        int nAssociatedhits = 0;
        int nDUThitsEvent = 0;
        int nHit = col->getNumberOfElements();
        double yposPrev = 0.;
        double xposPrev = 0.;
        int nPlanesWithMoreHits = 0;
        bool unfoundTrack = false;
        bool pAlpideHit = false;
        bool firstHit = true;
        for(int ihit=0; ihit< nHit ; ihit++)
        {
          TrackerHit        * hit = dynamic_cast<TrackerHit*>( col->getElementAt(ihit) ) ;
          bool hitOnSamePlane = false;
          double pos[3]={0.,0.,0.};
          if( hit != 0 )
          {
            for(int j=ihit; j<nHit && firstHit && nPlanesWithMoreHits <= _nPlanesWithMoreHits; j++)
            {
              TrackerHit * hitcheck1 = dynamic_cast<TrackerHit*>( col->getElementAt(j) ) ;
              const double *poscheck1 = hitcheck1->getPosition();
              if (poscheck1[2]<=dutZ+zDistance && poscheck1[2]>=dutZ-zDistance) continue;
              for(int k=j+1; k<nHit && nPlanesWithMoreHits <= _nPlanesWithMoreHits; k++)
              {
                TrackerHit * hitcheck2 = dynamic_cast<TrackerHit*>( col->getElementAt(k) ) ;
                const double *poscheck2 = hitcheck2->getPosition();
                if ((abs(poscheck1[2] - poscheck2[2]) < maxDistInPlane) && !hitOnSamePlane)
                {
                  hitOnSamePlane = true;
                  nPlanesWithMoreHits++;
                  break;
		}
                else if (k == j+1 && poscheck1[2] != poscheck2[2]) hitOnSamePlane = false;
              }
            }
            firstHit = false;
            if (nPlanesWithMoreHits > _nPlanesWithMoreHits) {unfoundTrack = true; break;}
            const double *pos0 = hit->getPosition();
            pos[0] = pos0[0];
            pos[1] = pos0[1];
            pos[2] = pos0[2];
            if (pos[2] >= dutZ-zDistance && pos[2] <= dutZ+zDistance )
            {
              pAlpideHit = true;
              pos[0]    -= xZero;
              pos[1]    -= yZero;
              _EulerRotationBack( pos, gRotation );

              double sign = 0;
              if      ( xPointing[0] < -0.7 )       sign = -1 ;
              else if ( xPointing[0] > 0.7 )       sign =  1 ;
              else {
                if       ( xPointing[1] < -0.7 )    sign = -1 ;
                else if  ( xPointing[1] > 0.7 )    sign =  1 ;
              }
              pos[0]    -=  ( -1 ) * sign * xSize / 2;
              if      ( yPointing[0] < -0.7 )       sign = -1 ;
              else if ( yPointing[0] > 0.7 )       sign =  1 ;
              else {
                if       ( yPointing[1] < -0.7 )    sign = -1 ;
                else if  ( yPointing[1] > 0.7 )    sign =  1 ;
              }
              pos[1]    -= ( -1 ) * sign * ySize / 2;

              double ypos = (xPointing[0]*pos[1]-yPointing[0]*pos[0])/(yPointing[1]*xPointing[0]-yPointing[0]*xPointing[1]);
              double xpos = pos[0]/xPointing[0] - xPointing[1]/xPointing[0]*ypos;
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
              if (!hitmapFilled) hitmapHisto->Fill(xpos,ypos);
#endif

              nDUThitsEvent++;
              if (nDUThitsEvent > 1) nDUThits++;
              if (abs(xpos-xposfit) < limit && abs(ypos-yposfit) < limit )
              {
                nAssociatedhits++;
                for (int jhit=ihit+1; jhit< nHit ; jhit++)
                {
                  TrackerHit * hitNext = dynamic_cast<TrackerHit*>( col->getElementAt(jhit) ) ;
                  double posNext[3]={0.,0.,0.};
                  if( hitNext != 0 )
                  {
                    const double *pos0Next = hitNext->getPosition();
                    posNext[0] = pos0Next[0];
                    posNext[1] = pos0Next[1];
                    posNext[2] = pos0Next[2];
                    if (posNext[2] < dutZ-zDistance || posNext[2] > dutZ+zDistance ) continue;
                    posNext[0]    -= xZero;
                    posNext[1]    -= yZero;

                    _EulerRotationBack( posNext, gRotation );

                    if      ( xPointing[0] < -0.7 )       sign = -1 ;
                    else if ( xPointing[0] > 0.7 )       sign =  1 ;
                    else {
                      if       ( xPointing[1] < -0.7 )    sign = -1 ;
                      else if  ( xPointing[1] > 0.7 )    sign =  1 ;
                    }
                    posNext[0]    -=  ( -1 ) * sign * xSize / 2;
                    if      ( yPointing[0] < -0.7 )       sign = -1 ;
                    else if ( yPointing[0] > 0.7 )       sign =  1 ;
                    else {
                      if       ( yPointing[1] < -0.7 )    sign = -1 ;
                      else if  ( yPointing[1] > 0.7 )    sign =  1 ;
                    }
                    posNext[1]    -= ( -1 ) * sign * ySize / 2;

                    double yposNext = (xPointing[0]*posNext[1]-yPointing[0]*posNext[0])/(yPointing[1]*xPointing[0]-yPointing[0]*xPointing[1]);
                    double xposNext = posNext[0]/xPointing[0] - xPointing[1]/xPointing[0]*yposNext;
                    if (abs(xposNext-xposfit) > limit || abs(yposNext-yposfit) > limit) continue;
                    nAssociatedhits++;
                    if ((xpos-xposfit)*(xpos-xposfit)+(ypos-yposfit)*(ypos-yposfit)<(xposNext-xposfit)*(xposNext-xposfit)+(yposNext-yposfit)*(yposNext-yposfit))
                    {
                      ihit = jhit;
                    }
                    else
                    {
                      xpos = xposNext;
                      ypos = yposNext;
                      pos[2] = posNext[2];
                      hit = hitNext;
                      ihit = jhit;
                    }
                  }
                }
                if (nDUThitsEvent > 1 && nAssociatedhits == 1) {tmpHist->Fill(xposfitPrev,yposfitPrev); nWrongPAlpideHit--;}
                if (nAssociatedhits > 1)
                  streamlog_out ( DEBUG )  << nAssociatedhits << " points for one track in DUT in event " << evt->getEventNumber() << "\t" << xposPrev << "\t" << yposPrev << "\t" << xpos << "\t" << ypos << " Fit: " << xposfit << "\t" << yposfit << " Number of planes with more than one hit: " << nPlanesWithMoreHits << endl;
                if (_clusterAvailable)
                {
                  for ( unsigned int idetector=0 ; idetector<zsInputDataCollectionVec->size(); idetector++)
                  {
                    CellIDDecoder<TrackerDataImpl> cellDecoder( zsInputDataCollectionVec );
                    TrackerDataImpl * zsData = dynamic_cast< TrackerDataImpl * > ( zsInputDataCollectionVec->getElementAt(idetector) );
                    SparsePixelType   type   = static_cast<SparsePixelType> ( static_cast<int> (cellDecoder( zsData )["sparsePixelType"]) );
//                    int iCluster = 0;
                    if (hit->getTime() == zsData->getTime())
                    {
                      nClusterAssociatedToTrackPerEvent++;
                      clusterAssosiatedToTrack.push_back(zsData->getTime());
                      int clusterSize = zsData->getChargeValues().size()/4;
                      vector<int> X(clusterSize);
                      vector<int> Y(clusterSize);
                      Cluster cluster;
                      if ( type == kEUTelGenericSparsePixel )
                      {
                        vector<vector<int> > pixVector;
                        auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel > > sparseData( new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> ( zsData ) );
                        EUTelGenericSparsePixel* pixel = new EUTelGenericSparsePixel;
                        for(unsigned int iPixel = 0; iPixel < sparseData->size(); iPixel++ )
                        {
                          sparseData->getSparsePixelAt( iPixel, pixel );
                          X[iPixel] = pixel->getXCoord();
                          Y[iPixel] = pixel->getYCoord();
                          vector<int> pix;
                          pix.push_back(X[iPixel]);
                          pix.push_back(Y[iPixel]);
                          pixVector.push_back(pix);
                        }
                        delete pixel;
                        cluster.set_values(clusterSize,X,Y);
                        clusterSizeHisto[index]->Fill(clusterSize);
                        int xMin = *min_element(X.begin(), X.end());
                        int xMax = *max_element(X.begin(), X.end());
                        int yMin = *min_element(Y.begin(), Y.end());
                        int yMax = *max_element(Y.begin(), Y.end());
                        int clusterWidthX = xMax - xMin + 1;
                        int clusterWidthY = yMax - yMin + 1;

                        if ((clusterWidthX > 3 || clusterWidthY > 3) && !emptyMiddle(pixVector))
                          for (unsigned int iPixel=0; iPixel<pixVector.size(); iPixel++)
                            largeClusterHistos->Fill(pixVector[iPixel][0],pixVector[iPixel][1]);
                        if (emptyMiddle(pixVector))
                        {
                          for (unsigned int iPixel=0; iPixel<pixVector.size(); iPixel++)
                            circularClusterHistos->Fill(pixVector[iPixel][0],pixVector[iPixel][1]);
                        }

                        clusterWidthXHisto[index]->Fill(clusterWidthX);
                        clusterWidthYHisto[index]->Fill(clusterWidthY);
                        clusterWidthXVsXHisto[index]->Fill(fmod(xposfit,xPitch),clusterWidthX);
                        clusterWidthXVsXAverageHisto[index]->Fill(fmod(xposfit,xPitch),clusterWidthX);
                        clusterWidthYVsYHisto[index]->Fill(fmod(yposfit,yPitch),clusterWidthY);
                        clusterWidthYVsYAverageHisto[index]->Fill(fmod(yposfit,yPitch),clusterWidthY);
                        clusterSize2DHisto[index]->Fill(fmod(xposfit,xPitch),fmod(yposfit,yPitch),clusterSize);
                        clusterSize2D2by2Histo[index]->Fill(fmod(xposfit,2*xPitch),fmod(yposfit,2*yPitch),clusterSize);
                        clusterSize2DAverageHisto[index]->Fill(fmod(xposfit,xPitch),fmod(yposfit,yPitch),clusterSize);
                        clusterSize2DAverage2by2Histo[index]->Fill(fmod(xposfit,2*xPitch),fmod(yposfit,2*yPitch),clusterSize);
                        nClusterVsXHisto[index]->Fill(fmod(xposfit,xPitch));
                        nClusterVsYHisto[index]->Fill(fmod(yposfit,yPitch));
                        nClusterSizeHisto[index]->Fill(fmod(xposfit,xPitch),fmod(yposfit,yPitch));
                        nClusterSize2by2Histo[index]->Fill(fmod(xposfit,2*xPitch),fmod(yposfit,2*yPitch));
                        int clusterShape = cluster.WhichClusterShape(cluster, clusterVec);
                        if (clusterShape>=0)
                        {
                          clusterShapeHisto->Fill(clusterShape);
                          clusterShapeX[clusterShape]->Fill(xMin);
                          clusterShapeY[clusterShape]->Fill(yMin);
                          clusterShape2D2by2[clusterShape]->Fill(fmod(xposfit,2*xPitch),fmod(yposfit,2*yPitch));
                          for (unsigned int iGroup=0; iGroup<symmetryGroups.size(); iGroup++)
                            for (unsigned int iMember=0; iMember<symmetryGroups[iGroup].size(); iMember++)
                              if (symmetryGroups[iGroup][iMember] == clusterShape) clusterShape2DGrouped2by2[iGroup]->Fill(fmod(xposfit,2*xPitch),fmod(yposfit,2*yPitch));
                        }
                        else clusterShapeHisto->Fill(clusterVec.size());
                      }
                      break;
                    }
                  }
                }
                nTracksPAlpide[index]++;
                TrackImpl* track = static_cast<TrackImpl*> (colTrack->getElementAt(ifit/7));
                float chi2 = track->getChi2();
                chi22DHisto->Fill(xposfit,yposfit,chi2);
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
                chi2Histo->Fill(chi2);
                scatteringAngleHisto->Fill(xposfit,yposfit,abs(xpos-xposfit));
                for (unsigned int i=0; i<chi2Max.size(); i++)
                {
                  if (chi2 < chi2Max[i] && yposfit < 12.2 && yposfit > 9.7 && xposfit < 26.2 && xposfit > 2.5)
                  {
                    residualXPAlpide[chi2Max[i]][index]->Fill(xpos-xposfit);
                    residualYPAlpide[chi2Max[i]][index]->Fill(ypos-yposfit);
                    residualZPAlpide[chi2Max[i]][index]->Fill(pos[2]-fitpos[2]);
                    residualXPixel[chi2Max[i]][index]->Fill(fmod(xposfit,xPitch),fmod(yposfit,yPitch),abs(xpos-xposfit));
                    residualYPixel[chi2Max[i]][index]->Fill(fmod(xposfit,xPitch),fmod(yposfit,yPitch),abs(ypos-yposfit));
                    residualXAveragePixel[chi2Max[i]][index]->Fill(fmod(xposfit,xPitch),fmod(yposfit,yPitch),abs(xpos-xposfit));
                    residualYAveragePixel[chi2Max[i]][index]->Fill(fmod(xposfit,xPitch),fmod(yposfit,yPitch),abs(ypos-yposfit));
                    nResidualXPixel[chi2Max[i]][index]->Fill(fmod(xposfit,xPitch),fmod(yposfit,yPitch));
                    nResidualYPixel[chi2Max[i]][index]->Fill(fmod(xposfit,xPitch),fmod(yposfit,yPitch));
                    residualXPixel2by2[chi2Max[i]][index]->Fill(fmod(xposfit,2*xPitch),fmod(yposfit,2*yPitch),abs(xpos-xposfit));
                    residualYPixel2by2[chi2Max[i]][index]->Fill(fmod(xposfit,2*xPitch),fmod(yposfit,2*yPitch),abs(ypos-yposfit));
                    residualXAveragePixel2by2[chi2Max[i]][index]->Fill(fmod(xposfit,2*xPitch),fmod(yposfit,2*yPitch),abs(xpos-xposfit));
                    residualYAveragePixel2by2[chi2Max[i]][index]->Fill(fmod(xposfit,2*xPitch),fmod(yposfit,2*yPitch),abs(ypos-yposfit));
                    nResidualXPixel2by2[chi2Max[i]][index]->Fill(fmod(xposfit,2*xPitch),fmod(yposfit,2*yPitch));
                    nResidualYPixel2by2[chi2Max[i]][index]->Fill(fmod(xposfit,2*xPitch),fmod(yposfit,2*yPitch));
                  }
                  else if (chi2 < chi2Max[i] && (yposfit < 9.2 || yposfit > 12.7 || xposfit > 27.4 || xposfit < 1.2))
                  {
                    residualXPCBPAlpide[chi2Max[i]][index]->Fill(xpos-xposfit);
                    residualYPCBPAlpide[chi2Max[i]][index]->Fill(ypos-yposfit);
                    residualZPCBPAlpide[chi2Max[i]][index]->Fill(pos[2]-fitpos[2]);

                  }
                }
                tracksPAlpideHisto->Fill(xposfit,yposfit);
                efficiencyHisto->Fill(xposfit,yposfit);
                tracksPAlpidePixelHisto[index]->Fill(fmod(xposfit,xPitch),fmod(yposfit,yPitch));
                tracksPAlpidePixel2by2Histo[index]->Fill(fmod(xposfit,2*xPitch),fmod(yposfit,2*yPitch));
                efficiencyPixelHisto[index]->Fill(fmod(xposfit,xPitch),fmod(yposfit,yPitch));
                efficiencyPixel2by2Histo[index]->Fill(fmod(xposfit,2*xPitch),fmod(yposfit,2*yPitch));
#endif
                yposPrev = ypos;
                xposPrev = xpos;
              }
              else if (nAssociatedhits < 1 && nDUThitsEvent == 1)
              {
                hitmapWrongHitHisto->Fill(xposfit,yposfit);
                nWrongPAlpideHit++;
                xposfitPrev = xposfit;
                yposfitPrev = yposfit;
              }
            }
            else continue;
          }
        }
        hitmapFilled = true;
        if (unfoundTrack) {nPlanesWithTooManyHits++; break;}
        else
        {
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
          tracksHisto->Fill(xposfit,yposfit);
          tracksPixelHisto[index]->Fill(fmod(xposfit,xPitch),fmod(yposfit,yPitch));
          tracksPixel2by2Histo[index]->Fill(fmod(xposfit,2*xPitch),fmod(yposfit,2*yPitch));
#endif
          nTracks[index]++;
        }
        if(!pAlpideHit) {hitmapNoHitHisto->Fill(xposfit,yposfit); nNoPAlpideHit++;}
      }
    }
    else continue;
  }
  _nEvents++;
  if (fitHitAvailable) _nEventsWithTrack++;
  if (_clusterAvailable)
  {
    for ( unsigned int i=0 ; i<zsInputDataCollectionVec->size(); i++)
    {
      CellIDDecoder<TrackerDataImpl> cellDecoder( zsInputDataCollectionVec );
      TrackerDataImpl * zsData = dynamic_cast< TrackerDataImpl * > ( zsInputDataCollectionVec->getElementAt(i) );
      int sensorID = static_cast<int> (cellDecoder( zsData )["sensorID"]);
      if (sensorID == _dutID)
      {
        bool isAssosiated = false;
        for (unsigned int iAssociatedCluster=0; iAssociatedCluster<clusterAssosiatedToTrack.size(); iAssociatedCluster++)
          if (zsData->getTime() == clusterAssosiatedToTrack[iAssociatedCluster])
          {
            isAssosiated = true;
            break;
          }
        if (!isAssosiated)
        {
          int index = -1;
          SparsePixelType   type   = static_cast<SparsePixelType> ( static_cast<int> (cellDecoder( zsData )["sparsePixelType"]) );
          if ( type == kEUTelGenericSparsePixel )
          {
            Cluster cluster;
            int clusterSize = zsData->getChargeValues().size()/4;
            vector<int> X(clusterSize);
            vector<int> Y(clusterSize);
            auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel > > sparseData( new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> ( zsData ) );
            EUTelGenericSparsePixel* pixel = new EUTelGenericSparsePixel;
            for(unsigned int iPixel = 0; iPixel < sparseData->size(); iPixel++ )
            {
              sparseData->getSparsePixelAt( iPixel, pixel );
              X[iPixel] = pixel->getXCoord();
              Y[iPixel] = pixel->getYCoord();
              if (!fitHitAvailable) nFakeWithoutTrackHitmapHisto->Fill(X[iPixel],Y[iPixel]);
              else nFakeWithTrackHitmapHisto->Fill(X[iPixel],Y[iPixel]);
              nFakeHitmapHisto->Fill(X[iPixel],Y[iPixel]);
            }
            delete pixel;
            cluster.set_values(clusterSize,X,Y);
            float xCenter, yCenter;
            cluster.getCenterOfGravity(xCenter,yCenter);
            for (int iSector=0; iSector<4; iSector++)
            {
              if (xCenter>xPixel/4.*iSector && xCenter<xPixel/4.*(iSector+1))
              {
                index = iSector;
                break;
              }
            }
            if (index != -1)
            {
              if (!fitHitAvailable) nFakeWithoutTrack[index] += clusterSize;
              else nFakeWithTrack[index] += clusterSize;
              nFake[index] += clusterSize;
            }
          }
        }
      }
    }
  }
  if (_clusterAvailable && fitHitAvailable)
  {
    for ( unsigned int i=0 ; i<zsInputDataCollectionVec->size(); i++)
    {
      CellIDDecoder<TrackerDataImpl> cellDecoder( zsInputDataCollectionVec );
      TrackerDataImpl * zsData = dynamic_cast< TrackerDataImpl * > ( zsInputDataCollectionVec->getElementAt(i) );
      int sensorID = static_cast<int> (cellDecoder( zsData )["sensorID"]);
      if (sensorID == _dutID)
      {
        SparsePixelType   type   = static_cast<SparsePixelType> ( static_cast<int> (cellDecoder( zsData )["sparsePixelType"]) );
        if ( type == kEUTelGenericSparsePixel )
        {
          Cluster cluster;
          int clusterSize = zsData->getChargeValues().size()/4;
          vector<int> X(clusterSize);
          vector<int> Y(clusterSize);
          vector<vector<int> > pixVector;
          auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel > > sparseData( new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> ( zsData ) );
          EUTelGenericSparsePixel* pixel = new EUTelGenericSparsePixel;
          for(unsigned int iPixel = 0; iPixel < sparseData->size(); iPixel++ )
          {
            sparseData->getSparsePixelAt( iPixel, pixel );
            X[iPixel] = pixel->getXCoord();
            Y[iPixel] = pixel->getYCoord();
          }
          delete pixel;
          cluster.set_values(clusterSize,X,Y);
          float xCenter, yCenter;
          cluster.getCenterOfGravity(xCenter,yCenter);
          bool isAssosiated = false;
          for(int ifit=0; ifit< nFitHit ; ifit++)
          {
            TrackerHit * fithit = dynamic_cast<TrackerHit*>( colFit->getElementAt(ifit) ) ;
            double fitpos[3]={0.,0.,0.};
            const double *fitpos0 = fithit->getPosition();
            fitpos[0] = fitpos0[0];
            fitpos[1] = fitpos0[1];
            fitpos[2] = fitpos0[2];
            if (fitpos[2] >= dutZ-zDistance && fitpos[2] <= dutZ+zDistance )
            {
              double xposfit=0, yposfit=0;
              RemoveAlign(preAlignmentCollectionVec,alignmentCollectionVec,alignmentPAlpideCollectionVec,fitpos,xposfit,yposfit);

              if (abs(xposfit-(double)xCenter/xPixel*xSize)<limit && abs(yposfit-(double)yCenter/yPixel*ySize)<limit )
              {
                isAssosiated = true;
                break;
              }
            }
          }
          if (isAssosiated) continue;
          for (unsigned int iPixel=0; iPixel<X.size(); iPixel++)
            nFakeWithTrackHitmapCorrectedHisto->Fill(X[iPixel],Y[iPixel]);
          int index = -1;
          for (int iSector=0; iSector<4; iSector++)
          {
            if (xCenter>xPixel/4.*iSector && xCenter<xPixel/4.*(iSector+1))
            {
              index = iSector;
              break;
            }
          }
          if (index != -1)
          {
            nFakeWithTrackCorrected[index] += clusterSize;
          }
        }
      }
    }
  }
  nTrackPerEventHisto->Fill(nTrackPerEvent);
  nClusterAssociatedToTrackPerEventHisto->Fill(nClusterAssociatedToTrackPerEvent);
  nClusterPerEventHisto->Fill(nClusterPerEvent);
}


#ifdef MARLIN_USE_AIDA
void EUTelProcessorAnalysisPALPIDEfs::bookHistos()
{
  streamlog_out ( DEBUG1 )  << "Booking histograms " << endl;
  auto_ptr<EUTelHistogramManager> histoMgr( new EUTelHistogramManager( _histoInfoFileName ));

  try {
    histoMgr->init();
  } catch ( ios::failure& e) {
    streamlog_out ( WARNING2 ) << "I/O problem with " << _histoInfoFileName << "\n"
                               << "Continuing without histogram manager"  << endl;
  } catch ( marlin::ParseException& e ) {
    streamlog_out ( WARNING2 ) << e.what() << "\n"
                               << "Continuing without histogram manager" << endl;
  }
  AIDAProcessor::tree(this)->cd("Analysis");
  tracksHisto = new TH2I("tracksHisto","Tracks;X (mm);Y (mm)",xPixel,0,xSize,yPixel,0,ySize);
  tracksPAlpideHisto = new TH2I("tracksPAlpideHisto","Tracks found in pALPIDEfs;X (mm);Y (mm)",xPixel,0,xSize,yPixel,0,ySize);
  efficiencyHisto = new TH2F("efficiencyHisto","Efficiency of the pALPIDEfs;X (mm);Y (mm)",xPixel,0,xSize,yPixel,0,ySize);
  hitmapHisto = new TH2I("hitmapHisto","Hitmap of pALPIDEfs;X (mm);Y (mm)",xPixel,0,xSize,yPixel,0,ySize);
  hitmapNoHitHisto = new TH2I("hitmapNoHitHisto","Hitmap of tracks which didn't have a hit in the pALPIDEfs;X (mm);Y (mm)",xPixel,0,xSize,yPixel,0,ySize);
  hitmapWrongHitHisto = new TH2I("hitmapWrongHitHisto","Hitmap of tracks which didn't have a hit in the pALPIDEfs close to them;X (mm);Y (mm)",xPixel,0,xSize,yPixel,0,ySize);
  nFakeHitmapHisto = new TH2I("nFakeHitmapHisto","Hit map of hits not associated to tracks;X (pixel);Y (pixel)",xPixel,0,xPixel,yPixel,0,yPixel);
  nFakeWithTrackHitmapHisto = new TH2I("nFakeWithTrackHitmapHisto","Hit map of hits not associated to tracks with track in event;X (pixel);Y (pixel)",xPixel,0,xPixel,yPixel,0,yPixel);
  nFakeWithTrackHitmapCorrectedHisto = new TH2I("nFakeWithTrackHitmapCorrectedHisto","Corrected hit map of hits not associated to tracks with track in event;X (pixel);Y (pixel)",xPixel,0,xPixel,yPixel,0,yPixel);
  nFakeWithoutTrackHitmapHisto = new TH2I("nFakeWithoutTrackHitmapHisto","Hit map of hits not associated to trackswithout track in event;X (pixel);Y (pixel)",xPixel,0,xPixel,yPixel,0,yPixel);
  scatteringAngleHisto = new TProfile2D( "scatteringAngleHisto","Scattering angles;X (mm);Y (mm)",xPixel,0,xSize,yPixel,0,ySize);
  chi22DHisto = new TProfile2D( "chi22DHisto","#chi^{2};X (mm);Y (mm)",xPixel,0,xSize,yPixel,0,ySize);
  tmpHist = new TH2I("tmpHist","",xPixel,0,xSize,yPixel,0,ySize);
  chi2Histo = new TH1I("chi2Histo","#chi^{2} of tracks used for the analysis;#chi^{2};a.u.",100,0,50);
  hotpixelHisto = new TH2I("hotpixelHisto","Hot pixels in the DUT;X (mm);Y (mm)",xPixel,0,xSize,yPixel,0,ySize);
  deadColumnHisto = new TH2I("deadColumnHisto","Dead columns in the DUT;X (mm);Y (mm)",xPixel,0,xSize,yPixel,0,ySize);
  largeClusterHistos = new TH2I("largeClusterHisto","Clusters with more than 3 pixels;X (mm);Y (mm)",xPixel,0,xPixel,yPixel,0,yPixel);
  circularClusterHistos = new TH2I("circularClusterHisto","Circular clusters (with missing hits in the middle);X (mm);Y (mm)",xPixel,0,xPixel,yPixel,0,yPixel);
  timeStampHisto = new TH1I("timeStampHisto","Distribution of the time stamp of the events; Time stamp (in 12.5 ns units)",1000,0,50000);
  nFakeHisto = new TH1F("nFakeHisto","Noise occupancy per sector for all events;Sector;Noise occupancy (/event/pixel)",4,0,4);
  nFakeWithTrackHisto = new TH1F("nFakeWithTrackHisto","Noise occupancy per sector for events with track;Sector;Noise occupancy (/event/pixel)",4,0,4);
  nFakeWithTrackCorrectedHisto = new TH1F("nFakeWithTrackCorrectedHisto","Corrected noise occupancy per sector for events with track;Sector;Noise occupancy (/event/pixel)",4,0,4);
  nFakeWithoutTrackHisto = new TH1F("nFakeWithoutTrackHisto","Noise occupancy per sector for events without track;Sector;Noise occupancy (/event/pixel)",4,0,4);
  nTrackPerEventHisto = new TH1I("nTrackPerEventHisto","Number of tracks per event;Number of tracks;a.u.",30,0,30);
  nClusterAssociatedToTrackPerEventHisto = new TH1I("nClusterAssociatedToTrackPerEventHisto","Number of clusters associated to tracks per event;Number of clusters;a.u.",30,0,30);
  nClusterPerEventHisto = new TH1I("nClusterPerEventHisto","Number of clusters per event;Number of clusters;a.u.",30,0,30);
  for (int iSector=0; iSector<4; iSector++)
  {
    AIDAProcessor::tree(this)->mkdir(Form("Sector_%d",iSector));
    AIDAProcessor::tree(this)->cd(Form("Sector_%d",iSector));
    for (unsigned int i=0; i<chi2Max.size();i++)
    {
      residualXPAlpide[chi2Max[i]][iSector] = new TH1I(Form("residualXPAlpide_%.1f_%d",chi2Max[i],iSector),Form("Residual X (Max chi2 = %.1f), sector %d;X (mm);a.u.",chi2Max[i],iSector),300,-0.2,0.2);
      residualYPAlpide[chi2Max[i]][iSector] = new TH1I(Form("residualYPAlpide_%.1f_%d",chi2Max[i],iSector),Form("Residual Y (Max chi2 = %.1f), sector %d;Y (mm);a.u.",chi2Max[i],iSector),300,-0.2,0.2);
      residualZPAlpide[chi2Max[i]][iSector] = new TH1I(Form("residualZPAlpide_%.1f_%d",chi2Max[i],iSector),Form("Residual Z (Max chi2 = %.1f), sector %d;Z (mm);a.u.",chi2Max[i],iSector),100,-0.3,0.3);
      residualXPCBPAlpide[chi2Max[i]][iSector] = new TH1I(Form("residualXPCBPAlpide_%.1f_%d",chi2Max[i],iSector),Form("Residual X not uner the hole (Max chi2 = %.1f), sector %d;X (mm);a.u.",chi2Max[i],iSector),300,-0.2,0.2);
      residualYPCBPAlpide[chi2Max[i]][iSector] = new TH1I(Form("residualYPCBPAlpide_%.1f_%d",chi2Max[i],iSector),Form("Residual Y not uner the hole (Max chi2 = %.1f), sector %d;Y (mm);a.u.",chi2Max[i],iSector),300,-0.2,0.2);
      residualZPCBPAlpide[chi2Max[i]][iSector] = new TH1I(Form("residualZPCBPAlpide_%.1f_%d",chi2Max[i],iSector),Form("Residual Z not uner the hole (Max chi2 = %.1f), sector %d;Z (mm);a.u.",chi2Max[i],iSector),100,-0.3,0.3);
      residualXPixel[chi2Max[i]][iSector] = new TH2F(Form("residualXPixel_%.1f_%d",chi2Max[i],iSector),"Residual X added up;X (mm);Y (mm)",(int)(xSize/xPixel*1000),0,xSize/xPixel,(int)(ySize/yPixel*1000),0,ySize/yPixel);
      residualYPixel[chi2Max[i]][iSector] = new TH2F(Form("residualYPixel_%.1f_%d",chi2Max[i],iSector),"Residual Y added up;X (mm);Y (mm)",(int)(xSize/xPixel*1000),0,xSize/xPixel,(int)(ySize/yPixel*1000),0,ySize/yPixel);
      residualXAveragePixel[chi2Max[i]][iSector] = new TH2F(Form("residualAverageXPixel_%.1f_%d",chi2Max[i],iSector),"Average residual X;X (mm);Y (mm)",(int)(xSize/xPixel*1000),0,xSize/xPixel,(int)(ySize/yPixel*1000),0,ySize/yPixel);
      residualYAveragePixel[chi2Max[i]][iSector] = new TH2F(Form("residualAverageYPixel_%.1f_%d",chi2Max[i],iSector),"Average residual Y;X (mm);Y (mm)",(int)(xSize/xPixel*1000),0,xSize/xPixel,(int)(ySize/yPixel*1000),0,ySize/yPixel);
      nResidualXPixel[chi2Max[i]][iSector] = new TH2F(Form("nResidualXPixel_%.1f_%d",chi2Max[i],iSector),"Number of tracks used for residual X plot;X (mm);Y (mm)",(int)(xSize/xPixel*1000),0,xSize/xPixel,(int)(ySize/yPixel*1000),0,ySize/yPixel);
      nResidualYPixel[chi2Max[i]][iSector] = new TH2F(Form("nResidualYPixel_%.1f_%d",chi2Max[i],iSector),"Number of tracks used for residual Y plot;X (mm);Y (mm)",(int)(xSize/xPixel*1000),0,xSize/xPixel,(int)(ySize/yPixel*1000),0,ySize/yPixel);
      residualXPixel2by2[chi2Max[i]][iSector] = new TH2F(Form("residualXPixel2by2_%.1f_%d",chi2Max[i],iSector),"Residual X added up 2by2;X (mm);Y (mm)",(int)(xSize/xPixel*2000),0,2*xSize/xPixel,(int)(ySize/yPixel*2000),0,2*ySize/yPixel);
      residualYPixel2by2[chi2Max[i]][iSector] = new TH2F(Form("residualYPixel2by2_%.1f_%d",chi2Max[i],iSector),"Residual Y added up 2by2;X (mm);Y (mm)",(int)(xSize/xPixel*2000),0,2*xSize/xPixel,(int)(ySize/yPixel*2000),0,2*ySize/yPixel);
      residualXAveragePixel2by2[chi2Max[i]][iSector] = new TH2F(Form("residualAverageXPixel2by2_%.1f_%d",chi2Max[i],iSector),"Average residual X 2by2;X (mm);Y (mm)",(int)(xSize/xPixel*2000),0,2*xSize/xPixel,(int)(ySize/yPixel*2000),0,2*ySize/yPixel);
      residualYAveragePixel2by2[chi2Max[i]][iSector] = new TH2F(Form("residualAverageYPixel2by2_%.1f_%d",chi2Max[i],iSector),"Average residual Y 2by2;X (mm);Y (mm)",(int)(xSize/xPixel*2000),0,2*xSize/xPixel,(int)(ySize/yPixel*2000),0,2*ySize/yPixel);
      nResidualXPixel2by2[chi2Max[i]][iSector] = new TH2F(Form("nResidualXPixel2by2_%.1f_%d",chi2Max[i],iSector),"Number of tracks used for residual X 2by2 plot;X (mm);Y (mm)",(int)(xSize/xPixel*2000),0,2*xSize/xPixel,(int)(ySize/yPixel*2000),0,2*ySize/yPixel);
      nResidualYPixel2by2[chi2Max[i]][iSector] = new TH2F(Form("nResidualYPixel2by2_%.1f_%d",chi2Max[i],iSector),"Number of tracks used for residual Y 2by2 plot;X (mm);Y (mm)",(int)(xSize/xPixel*2000),0,2*xSize/xPixel,(int)(ySize/yPixel*2000),0,2*ySize/yPixel);
    }
    clusterWidthXHisto[iSector]  = new TH1I(Form("clusterWidthXHisto_%d",iSector),Form("Cluster width in X in sector %d;Cluster width X (pixel);a.u.",iSector),15,0.5,15.5);
    clusterWidthYHisto[iSector]  = new TH1I(Form("clusterWidthYHisto_%d",iSector),Form("Cluster width in Y in sector %d;Cluster width Y (pixel);a.u.",iSector),15,0.5,15.5);
    clusterSizeHisto[iSector]  = new TH1I(Form("clusterSizeHisto_%d",iSector),Form("Cluster size_%d;Cluster size (pixel);a.u.",iSector),20,0.5,20.5);
    efficiencyPixelHisto[iSector] = new TH2F(Form("efficiencyPixelHisto_%d",iSector),Form("Efficiency as function of place in the pixel of sector %d;X (mm);Y (mm)",iSector),(int)(xSize/xPixel*1000),0,xSize/xPixel,(int)(ySize/yPixel*1000),0,ySize/yPixel);
    efficiencyPixel2by2Histo[iSector] = new TH2F(Form("efficiencyPixel2by2Histo_%d",iSector),Form("Efficiency as function of place in four pixels of sector %d;X (mm);Y (mm)",iSector),(int)(xSize/xPixel*2000),0,2*xSize/xPixel,(int)(ySize/yPixel*2000),0,2*ySize/yPixel);
    tracksPAlpidePixelHisto[iSector] = new TH2F(Form("tracksPAlpidePixelHisto_%d",iSector),Form("Number of found tracks as function of place in the pixel of sector %d;X (mm);Y (mm)",iSector),(int)(xSize/xPixel*1000),0,xSize/xPixel,(int)(ySize/yPixel*1000),0,ySize/yPixel);
    tracksPAlpidePixel2by2Histo[iSector] = new TH2F(Form("tracksPAlpidePixel2by2Histo_%d",iSector),Form("Number of found tracks as function of place in four pixels of sector %d;X (mm);Y (mm)",iSector),(int)(xSize/xPixel*2000),0,2*xSize/xPixel,(int)(ySize/yPixel*2000),0,2*ySize/yPixel);
    tracksPixelHisto[iSector] = new TH2F(Form("tracksPixelHisto_%d",iSector),Form("Number of tracks as function of place in the pixel of sector %d;X (mm);Y (mm)",iSector),(int)(xSize/xPixel*1000),0,xSize/xPixel,(int)(ySize/yPixel*1000),0,ySize/yPixel);
    tracksPixel2by2Histo[iSector] = new TH2F(Form("tracksPixel2by2Histo_%d",iSector),Form("Number of tracks as function of place in four pixels of sector %d;X (mm);Y (mm)",iSector),(int)(xSize/xPixel*2000),0,2*xSize/xPixel,(int)(ySize/yPixel*2000),0,2*ySize/yPixel);
    clusterSize2DHisto[iSector] = new TH2F(Form("clusterSize2D_%d",iSector),Form("Added cluster size vs the track position of sector %d;X (mm);Y (mm)",iSector),(int)(xSize/xPixel*1000),0,xSize/xPixel,(int)(ySize/yPixel*1000),0,ySize/yPixel);
    clusterSize2D2by2Histo[iSector] = new TH2F(Form("clusterSize2D2by2_%d",iSector),Form("Added cluster size vs the track position in four pixels of sector %d;X (mm);Y (mm)",iSector),(int)(xSize/xPixel*2000),0,2*xSize/xPixel,(int)(ySize/yPixel*2000),0,2*ySize/yPixel);
    clusterSize2DAverageHisto[iSector] = new TH2F(Form("clusterSize2DAverage_%d",iSector),Form("Average cluster size vs the track position of sector %d;X (mm);Y (mm)",iSector),(int)(xSize/xPixel*1000),0,xSize/xPixel,(int)(ySize/yPixel*1000),0,ySize/yPixel);
    clusterSize2DAverage2by2Histo[iSector] = new TH2F(Form("clusterSize2DAverage2by2_%d",iSector),Form("Average cluster size vs the track position in four pixels of sector %d;X (mm);Y (mm)",iSector),(int)(xSize/xPixel*2000),0,2*xSize/xPixel,(int)(ySize/yPixel*2000),0,2*ySize/yPixel);
    nClusterSizeHisto[iSector] = new TH2F(Form("nClusterSizeHisto_%d",iSector),"Number of tracks used for cluster size analysis;X (mm);Y (mm)",(int)(xSize/xPixel*1000),0,xSize/xPixel,(int)(ySize/yPixel*1000),0,ySize/yPixel);
    nClusterSize2by2Histo[iSector] = new TH2F(Form("nClusterSize2by2Histo_%d",iSector),"Number of tracks used for cluster size analysis 2 by 2;X (mm);Y (mm)",(int)(xSize/xPixel*2000),0,2*xSize/xPixel,(int)(ySize/yPixel*2000),0,2*ySize/yPixel);
    clusterWidthXVsXHisto[iSector] = new TH1F(Form("clusterWidthXVsXHisto_%d",iSector),Form("Added cluster width in X vs the track position in X of sector %d;X (mm);a.u.",iSector),(int)(xSize/xPixel*1000),0,xSize/xPixel);
    clusterWidthYVsYHisto[iSector] = new TH1F(Form("clusterWidthYVsYHisto_%d",iSector),Form("Added cluster width in Y vs the track position in Y of sector %d;Y (mm);a.u",iSector),(int)(ySize/yPixel*1000),0,ySize/yPixel);
    clusterWidthXVsXAverageHisto[iSector] = new TH1F(Form("clusterWidthXVsXAverageHisto_%d",iSector),Form("Average cluster width in X vs the track position in X of sector %d;X (mm);a.u.",iSector),(int)(xSize/xPixel*1000),0,xSize/xPixel);
    clusterWidthYVsYAverageHisto[iSector] = new TH1F(Form("clusterWidthYVsYAverageHisto_%d",iSector),Form("Average cluster width in Y vs the track position in Y of sector %d;Y (mm);a.u",iSector),(int)(ySize/yPixel*1000),0,ySize/yPixel);
    nClusterVsXHisto[iSector] = new TH1F(Form("nClusterVsXHisto_%d",iSector),"Number of tracks used for the cluster width analysis in X;X (mm);a.u.",(int)(xSize/xPixel*1000),0,xSize/xPixel);
    nClusterVsYHisto[iSector] = new TH1F(Form("nClusterVsYHisto_%d",iSector),"Number of tracks used for the cluster width analysis in Y;Y (mm);a.u.",(int)(xSize/xPixel*1000),0,ySize/yPixel);
  }
  AIDAProcessor::tree(this)->mkdir("ClusterShape");
  AIDAProcessor::tree(this)->cd("ClusterShape");
  clusterShapeHisto = new TH1I("clusterShapeHisto","Cluster shape (all rotations separately);Cluster shape ID;a.u.",clusterVec.size()+1,-0.5,clusterVec.size()+0.5);
  clusterShapeHistoGrouped = new TH1I("clusterShapeHistoGrouped","Cluster shape (all rotations treated together);Cluster shape ID;a.u.",symmetryGroups.size(),-0.5,symmetryGroups.size()-0.5);
  differenceX = new TH1F("differenceX","Difference in the number of clusters which are mirorred around x;Cluster shape ID;a.u.",xPairs.size(),-0.5,xPairs.size()-0.5);
  differenceY = new TH1F("differenceY","Difference in the number of clusters which are mirorred around y;Cluster shape ID;a.u.",yPairs.size(),-0.5,yPairs.size()-0.5);
  for (unsigned int i=0; i<clusterVec.size();i++)
  {
    clusterShapeX[i] = new TH1I(Form("clusterShapeX_%d",i),Form("Place of cluster shape with ID %d as function of X;X (pixel);a.u.",i),xPixel,0,xPixel);
    clusterShapeY[i] = new TH1I(Form("clusterShapeY_%d",i),Form("Place of cluster shape with ID %d as function of Y;Y (pixel);a.u.",i),yPixel,0,yPixel);
   clusterShape2D2by2[i] = new TH2I(Form("clusterShape2D2by2_%d",i),Form("Distribuition of clusters with shape ID %d within four pixels;X (mm);Y (mm)",i),(int)(xSize/xPixel*2000),0,2*xSize/xPixel,(int)(ySize/yPixel*2000),0,2*ySize/yPixel);
  }
  int tmp = 0;
  for(map<int,int>::iterator it = xPairs.begin(); it != xPairs.end(); ++it)
  {
    clusterShapeDiffX.insert(make_pair(tmp,new TH1I(Form("clusterShapeDiffX_%d_%d",it->first,xPairs[it->first]),Form("Difference of the distributions between %d and %d;X (pixel);a.u.",it->first,xPairs[it->first]),xPixel,0,xPixel)));
    clusterShapeDiffY.insert(make_pair(tmp,new TH1I(Form("clusterShapeDiffY_%d_%d",it->first,xPairs[it->first]),Form("Difference of the distributions between %d and %d;Y (pixel);a.u.",it->first,xPairs[it->first]),yPixel,0,yPixel)));
    tmp++;
  }
  for(map<int,int>::iterator it = yPairs.begin(); it != yPairs.end(); ++it)
  {
    clusterShapeDiffX.insert(make_pair(tmp,new TH1I(Form("clusterShapeDiffX_%d_%d",it->first,yPairs[it->first]),Form("Difference of the distributions between %d and %d;X (pixel);a.u.",it->first,yPairs[it->first]),xPixel,0,xPixel)));
    clusterShapeDiffY.insert(make_pair(tmp,new TH1I(Form("clusterShapeDiffY_%d_%d",it->first,yPairs[it->first]),Form("Difference of the distributions between %d and %d;Y (pixel);a.u.",it->first,yPairs[it->first]),yPixel,0,yPixel)));
    tmp++;
  }
  for (unsigned int i=0; i<symmetryGroups.size(); i++)
  {
    string binName;
    for (unsigned int j=0; j<symmetryGroups[i].size(); j++)
    {
      if (j<symmetryGroups[i].size()-1) binName += Form("%d,",symmetryGroups[i][j]);
      else binName += Form("%d",symmetryGroups[i][j]);
    }
    string title = "Distribuition of clusters with IDs " + binName + " within four pixels;X (mm);Y (mm)";
    clusterShape2DGrouped2by2[i] = new TH2I(Form("clusterShapeGrouped2D2by2_%d",i),title.c_str(),(int)(xSize/xPixel*2000),0,2*xSize/xPixel,(int)(ySize/yPixel*2000),0,2*ySize/yPixel);
  }
  streamlog_out ( DEBUG5 )  << "end of Booking histograms " << endl;
}
#endif

void EUTelProcessorAnalysisPALPIDEfs::end()
{
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  efficiencyHisto->Divide(tracksHisto);
  for (int i=0; i<xPixel; i++)
    for (int j=0; j<yPixel; j++)
      if (tracksPAlpideHisto->GetBinContent(i,j) > tracksHisto->GetBinContent(i,j)) cerr << "More hits in the pAlpide than number of traks" << endl;
  for (int iSector=0; iSector<4; iSector++)
  {
    efficiencyPixelHisto[iSector]->Divide(tracksPixelHisto[iSector]);
    efficiencyPixel2by2Histo[iSector]->Divide(tracksPixel2by2Histo[iSector]);
    clusterSize2DAverageHisto[iSector]->Divide(nClusterSizeHisto[iSector]);
    clusterSize2DAverage2by2Histo[iSector]->Divide(nClusterSize2by2Histo[iSector]);
    clusterWidthXVsXAverageHisto[iSector]->Divide(nClusterVsXHisto[iSector]);
    clusterWidthYVsYAverageHisto[iSector]->Divide(nClusterVsYHisto[iSector]);
    AIDAProcessor::tree(this)->cd(Form("Sector_%d",iSector));
    CrossSection * csClusterSize = new CrossSection(clusterSize2DAverageHisto[iSector]);
    clusterSizeCrossSection[iSector] = csClusterSize->GetCrossSections();
    CrossSection * csEfficiency = new CrossSection(efficiencyPixelHisto[iSector]);
    efficiencyPixelCrossSection[iSector] = csEfficiency->GetCrossSections();
    for (unsigned int i=0; i<chi2Max.size(); i++)
    {
      residualXAveragePixel[chi2Max[i]][iSector]->Divide(nResidualXPixel[chi2Max[i]][iSector]);
      residualYAveragePixel[chi2Max[i]][iSector]->Divide(nResidualYPixel[chi2Max[i]][iSector]);
      residualXAveragePixel2by2[chi2Max[i]][iSector]->Divide(nResidualXPixel2by2[chi2Max[i]][iSector]);
      residualYAveragePixel2by2[chi2Max[i]][iSector]->Divide(nResidualYPixel2by2[chi2Max[i]][iSector]);
    }
    nFakeHisto->SetBinContent(iSector+1,(double)nFake[iSector]/_nEvents/(xPixel/4*yPixel));
    nFakeHisto->SetBinError(iSector+1,sqrt((double)nFake[iSector])/_nEvents/(xPixel/4*yPixel));
    nFakeWithTrackHisto->SetBinContent(iSector+1,(double)nFakeWithTrack[iSector]/_nEventsWithTrack/(xPixel/4*yPixel));
    nFakeWithTrackHisto->SetBinError(iSector+1,sqrt((double)nFakeWithTrack[iSector])/_nEventsWithTrack/(xPixel/4*yPixel));
    nFakeWithTrackCorrectedHisto->SetBinContent(iSector+1,(double)nFakeWithTrackCorrected[iSector]/_nEventsWithTrack/(xPixel/4*yPixel));
    nFakeWithTrackCorrectedHisto->SetBinError(iSector+1,sqrt((double)nFakeWithTrackCorrected[iSector])/_nEventsWithTrack/(xPixel/4*yPixel));
    nFakeWithoutTrackHisto->SetBinContent(iSector+1,(double)nFakeWithoutTrack[iSector]/(_nEvents-_nEventsWithTrack)/(xPixel/4*yPixel));
    nFakeWithoutTrackHisto->SetBinError(iSector+1,sqrt((double)nFakeWithoutTrack[iSector])/(_nEvents-_nEventsWithTrack)/(xPixel/4*yPixel));
    tracksProjection[iSector] = tracksHisto->ProjectionY(Form("tracksProjection_%d",iSector),xPixel/4*iSector,xPixel/4*(iSector+1));
    tracksPAlpideProjection[iSector] = tracksPAlpideHisto->ProjectionY(Form("tracksPAlpideProjection_%d",iSector),xPixel/4*iSector,xPixel/4*(iSector+1));
    efficiencyProjection[iSector] = (TH1D*)tracksPAlpideProjection[iSector]->Clone(Form("efficiencyProjection_%d",iSector));
    efficiencyProjection[iSector]->Divide(tracksProjection[iSector]);
    efficiencyProjection[iSector]->SetTitle("Efficiency projection in Y");
  }
  streamlog_out ( MESSAGE4 ) << "nEvents: " << _nEvents << endl;
  streamlog_out ( MESSAGE4 ) << "nEvents with tracks: " << _nEventsWithTrack << endl;
  streamlog_out ( MESSAGE4 ) << "nEvents without tracks: " << _nEvents-_nEventsWithTrack << endl;
  settingsFile << _nEvents << ";";

  hitmapWrongHitHisto->Add(tmpHist,-1.);
  if (_writeShapes)
  {
    ofstream shapeOutputFile;
    shapeOutputFile.open (_shapeOutputFileName.c_str(), ios::out );
    shapeOutputFile << "ID|Number of clusters with that shape" << endl;
    for (int i=1; i<=clusterShapeHisto->GetNbinsX(); i++)
      shapeOutputFile << i-1 << "|" << clusterShapeHisto->GetBinContent(i) << endl;
    shapeOutputFile.close();
  }
  int tmp=1;
  for(map<int,int>::iterator it = xPairs.begin(); it != xPairs.end(); ++it)
  {
    differenceX->SetBinContent(tmp,TMath::Abs(clusterShapeHisto->GetBinContent(it->first+1)-clusterShapeHisto->GetBinContent(xPairs[it->first]+1)));
    differenceX->GetXaxis()->SetBinLabel(tmp,Form("%d-%d",it->first,xPairs[it->first]));
    tmp++;
  }
  tmp=1;
  for(map<int,int>::iterator it = yPairs.begin(); it != yPairs.end(); ++it)
  {
    differenceY->SetBinContent(tmp,TMath::Abs(clusterShapeHisto->GetBinContent(it->first+1)-clusterShapeHisto->GetBinContent(yPairs[it->first]+1)));
    differenceY->GetXaxis()->SetBinLabel(tmp,Form("%d-%d",it->first,yPairs[it->first]));
    tmp++;
  }
  for (unsigned int i=0; i<symmetryGroups.size(); i++)
  {
    string binName;
    for (unsigned int j=0; j<symmetryGroups[i].size(); j++)
    {
      clusterShapeHistoGrouped->Fill(i,clusterShapeHisto->GetBinContent(symmetryGroups[i][j]+1));
      if (j<symmetryGroups[i].size()-1) binName += Form("%d,",symmetryGroups[i][j]);
      else binName += Form("%d",symmetryGroups[i][j]);
    }
    clusterShapeHistoGrouped->GetXaxis()->SetBinLabel(i+1,(char*)binName.c_str());
  }
  tmp = 0;
  for(map<int,int>::iterator it = xPairs.begin(); it != xPairs.end(); ++it)
  {
    for (int i=1; i<=clusterShapeDiffX[tmp]->GetNbinsX(); i++)
      clusterShapeDiffX[tmp]->SetBinContent(i,TMath::Abs(clusterShapeX[it->first]->GetBinContent(i)-clusterShapeX[xPairs[it->first]]->GetBinContent(i)));
    for (int i=1; i<=clusterShapeDiffY[tmp]->GetNbinsX(); i++)
      clusterShapeDiffY[tmp]->SetBinContent(i,TMath::Abs(clusterShapeY[it->first]->GetBinContent(i)-clusterShapeY[xPairs[it->first]]->GetBinContent(i)));
    tmp++;
  }
  for(map<int,int>::iterator it = yPairs.begin(); it != yPairs.end(); ++it)
  {
    for (int i=1; i<=clusterShapeDiffX[tmp]->GetNbinsX(); i++)
      clusterShapeDiffX[tmp]->SetBinContent(i,TMath::Abs(clusterShapeX[it->first]->GetBinContent(i)-clusterShapeX[yPairs[it->first]]->GetBinContent(i)));
    for (int i=1; i<=clusterShapeDiffY[tmp]->GetNbinsX(); i++)
      clusterShapeDiffY[tmp]->SetBinContent(i,TMath::Abs(clusterShapeY[it->first]->GetBinContent(i)-clusterShapeY[yPairs[it->first]]->GetBinContent(i)));
    tmp++;
  }
/*  TH1I* nX = new TH1I("nX","n",xPixel,0,xPixel);
  TH1I* nY = new TH1I("nY","n",yPixel,0,yPixel);
  for (unsigned int i=0; i<symmetryGroups.size(); i++)
  {
    nX->Reset();
    nY->Reset();
    string binName;
    double p = 1.0/symmetryGroups[i].size();
    for (unsigned int j=0; j<symmetryGroups[i].size(); j++)
    {
      if (j<symmetryGroups[i].size()-1) binName += Form("%d-",symmetryGroups[i][j]);
      else binName += Form("%d",symmetryGroups[i][j]);
      for (int k=1; k<=clusterShapeX[j]->GetNbinsX(); k++)
        nX->Fill(k-1,clusterShapeX[symmetryGroups[i][j]]->GetBinContent(k));
      for (int k=1; k<=clusterShapeY[j]->GetNbinsX(); k++)
        nY->Fill(k-1,clusterShapeY[symmetryGroups[i][j]]->GetBinContent(k));
    }
    chi2X.insert(make_pair(i, new TH1I(Form("chi2X_%s",(char*)binName.c_str()),"",xPixel,0,xPixel)));
    chi2Y.insert(make_pair(i, new TH1I(Form("chi2Y_%s",(char*)binName.c_str()),"",yPixel,0,yPixel)));
    for (unsigned int j=0; j<symmetryGroups[i].size(); j++)
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
//    for (int iBin=1; iBin<=chi2->GetNbinsX(); iBin++)
//      if (chi2->GetBinContent(iBin) > 7.815) cerr << "Symmetry group: " << binName << " Bin: " << iBin<< " is assymetric, with chi2 = " << chi2->GetBinContent(iBin) << endl;
  }
*/
#endif
  streamlog_out ( MESSAGE4 ) << nPlanesWithTooManyHits << " tracks had too many planes with more than one hit" << endl;
  streamlog_out ( MESSAGE4 ) << nNoPAlpideHit << " tracks didn't have a hit in the pALPIDEfs" << endl;
  streamlog_out ( MESSAGE4 ) << "For " << nWrongPAlpideHit << " tracks the pALPIDEfs had a hit far from the track" << endl;
  streamlog_out ( MESSAGE4 ) << nDUThits << " hits in the DUT weren't associated to a track" << endl;
  streamlog_out ( MESSAGE4 ) << "Overall efficiency of pALPIDEfs sectors: " << endl;
  for (int iSector=0; iSector<4; iSector++)
  {
    streamlog_out ( MESSAGE4 ) << (double)nTracksPAlpide[iSector]/nTracks[iSector] << "\t" << nTracks[iSector] << "\t" << nTracksPAlpide[iSector] << endl;
    settingsFile << (double)nTracksPAlpide[iSector]/nTracks[iSector] << ";" << nTracks[iSector] << ";" << nTracksPAlpide[iSector] << ";";
  }
  settingsFile << endl;
  streamlog_out ( MESSAGE4 )  << "Successfully finished" << endl;
}

void EUTelProcessorAnalysisPALPIDEfs::_EulerRotationBack(double* _telPos, double* _gRotation) {

    double z = _telPos[2] - dutZ;
    TVector3 _RotatedSensorHit( _telPos[0], _telPos[1], z );
    TVector3 _Xaxis( 1.0, 0.0, 0.0 );
    TVector3 _Yaxis( 0.0, 1.0, 0.0 );
    TVector3 _Zaxis( 0.0, 0.0, 1.0 );
    if( TMath::Abs(_gRotation[1]) > 1e-6 )
    {
        _RotatedSensorHit.Rotate( -1.*_gRotation[1], _Yaxis ); // in ZX
    }
    if( TMath::Abs(_gRotation[2]) > 1e-6 )
    {
        _RotatedSensorHit.Rotate( -1.*_gRotation[2], _Xaxis ); // in ZY
    }
    if( TMath::Abs(_gRotation[0]) > 1e-6 )
    {
        _RotatedSensorHit.Rotate( -1.*_gRotation[0], _Zaxis ); // in XY
    }

    _telPos[0] = _RotatedSensorHit.X();
    _telPos[1] = _RotatedSensorHit.Y();
    _telPos[2] = _RotatedSensorHit.Z() + dutZ;
}

int EUTelProcessorAnalysisPALPIDEfs::AddressToColumn(int ARegion, int ADoubleCol, int AAddress)
{
  int Column    = ARegion * 32 + ADoubleCol * 2;    // Double columns before ADoubleCol
  int LeftRight = ((AAddress % 4) < 2 ? 1:0);       // Left or right column within the double column
  Column += LeftRight;
  return Column;
}

int EUTelProcessorAnalysisPALPIDEfs::AddressToRow(int AAddress)
{
  // Ok, this will get ugly
  int Row = AAddress / 2;                // This is OK for the top-right and the bottom-left pixel within a group of 4
  if ((AAddress % 4) == 3) Row -= 1;      // adjust the top-left pixel
  if ((AAddress % 4) == 0) Row += 1;      // adjust the bottom-right pixel
  return Row;
}

bool EUTelProcessorAnalysisPALPIDEfs::emptyMiddle(vector<vector<int> > pixVector)
{
  bool holeX = false;
  bool holeY = false;
  for (unsigned int i=0; i<pixVector.size(); i++)
  {
    bool touchingX = false;
    bool lastX = true;
    for (unsigned int j=0; j<pixVector.size(); j++)
    {
      if (i==j) continue;
      if (pixVector[i][1] != pixVector[j][1]) continue;
      if (pixVector[i][0]+1 == pixVector[j][0]) {/*cerr << "Touching in x" << endl;*/ touchingX = true; break;}
      if (pixVector[i][0] <  pixVector[j][0]) {/*cerr << "Smaller in x"  << endl;*/ lastX  = false;}
    }
    if (!touchingX && !lastX) {/*cerr << "Hole in X" << endl;*/ holeX = true; break;}
  }
  for (unsigned int i=0; i<pixVector.size(); i++)
  {
    bool touchingY = false;
    bool lastY = true;
    for (unsigned int j=0; j<pixVector.size(); j++)
    {
      if (i==j) continue;
      if (pixVector[i][0] != pixVector[j][0]) continue;
      if (pixVector[i][1]+1 == pixVector[j][1]) {/*cerr << "Touching in y" << endl;*/ touchingY = true; break;}
      if (pixVector[i][1] <  pixVector[j][1]) {/*cerr << "Smaller in y"  << endl;*/ lastY  = false;}
    }
    if (!touchingY && !lastY) {/*cerr << "Hole in Y" << endl;*/ holeY = true; break;}
  }
  if (holeX && holeY) return true;
  else return false;
}

bool EUTelProcessorAnalysisPALPIDEfs::RemoveAlign(LCCollectionVec * preAlignmentCollectionVec, LCCollectionVec * alignmentCollectionVec, LCCollectionVec * alignmentPAlpideCollectionVec, double* fitpos, double& xposfit, double& yposfit)
{
      double xPlaneCenter    = geo::gGeometry().siPlaneXPosition(_dutID);
      double yPlaneCenter    = geo::gGeometry().siPlaneYPosition(_dutID);
      double zPlaneCenter    = geo::gGeometry().siPlaneZPosition(_dutID);
      TVector3 inputVec( fitpos[0] - xPlaneCenter, fitpos[1] - yPlaneCenter, fitpos[2] - zPlaneCenter);
      EUTelAlignmentConstant * alignment = 0;
      bool alignExist = false;
      EUTelAlignmentConstant * preAlignment = 0;
      bool prealignExist = false;
      EUTelAlignmentConstant * alignmentPAlpide = 0;
      bool alignPAlpideExist = false;
      for (int iAlign=0; iAlign<preAlignmentCollectionVec->getNumberOfElements(); iAlign++)
      {
        preAlignment = static_cast< EUTelAlignmentConstant * > (preAlignmentCollectionVec->getElementAt(iAlign));
        if (preAlignment->getSensorID() == _dutID)
        {
          prealignExist = true;
          inputVec[0] += preAlignment->getXOffset();
          inputVec[1] += preAlignment->getYOffset();
          inputVec[2] += preAlignment->getZOffset();
          break;
        }
      }
      if (!prealignExist) return false;//cerr << "No prealignment correction applied!" << endl;
      for (int iAlign=0; iAlign<alignmentCollectionVec->getNumberOfElements(); iAlign++)
      {
        alignment = static_cast< EUTelAlignmentConstant* > (alignmentCollectionVec->getElementAt(iAlign));
        if (alignment->getSensorID() == _dutID)
        {
          alignExist = true;
          inputVec[0] += alignment->getXOffset();
          inputVec[1] += alignment->getYOffset();
          inputVec[2] += alignment->getZOffset();
          break;
        }
      }
      if (!alignExist) return false;//cerr << "No alignment correction applied!" << endl;
      if (!_oneAlignmentCollection)
      {
        for (int iAlign=0; iAlign<alignmentPAlpideCollectionVec->getNumberOfElements(); iAlign++)
        {
          alignmentPAlpide = static_cast< EUTelAlignmentConstant* > (alignmentPAlpideCollectionVec->getElementAt(iAlign));
          if (alignmentPAlpide->getSensorID() == _dutID)
          {
            alignPAlpideExist = true;
            inputVec[0] += alignmentPAlpide->getXOffset();
            inputVec[1] += alignmentPAlpide->getYOffset();
            inputVec[2] += alignmentPAlpide->getZOffset();
            break;
          }
        }
        if (!alignPAlpideExist && _isFirstEvent) cerr << "No second alignment correction applied to the pAlpide!" << endl;
        if  (alignPAlpideExist)
        {
          inputVec.RotateX( alignmentPAlpide->getAlpha() );
          inputVec.RotateY( alignmentPAlpide->getBeta() );
          inputVec.RotateZ( alignmentPAlpide->getGamma() );
        }
      }
      if (alignExist)
      {
        inputVec.RotateX( alignment->getAlpha() );
        inputVec.RotateY( alignment->getBeta() );
        inputVec.RotateZ( alignment->getGamma() );
      }

      fitpos[0] = inputVec.X() + xPlaneCenter;
      fitpos[1] = inputVec.Y() + yPlaneCenter;

      fitpos[0] -= xZero;
      fitpos[1] -= yZero;

      _EulerRotationBack( fitpos, gRotation );

      double sign = 0;
      if      ( xPointing[0] < -0.7 )       sign = -1 ;
      else if ( xPointing[0] > 0.7 )       sign =  1 ;
      else {
        if       ( xPointing[1] < -0.7 )    sign = -1 ;
        else if  ( xPointing[1] > 0.7 )    sign =  1 ;
      }
      fitpos[0] -=  ( -1 ) * sign * xSize / 2;

      if      ( yPointing[0] < -0.7 )       sign = -1 ;
      else if ( yPointing[0] > 0.7 )       sign =  1 ;
      else {
        if       ( yPointing[1] < -0.7 )    sign = -1 ;
        else if  ( yPointing[1] > 0.7 )    sign =  1 ;
      }
      fitpos[1] -= ( -1 ) * sign * ySize / 2;
      yposfit = (xPointing[0]*fitpos[1]-yPointing[0]*fitpos[0])/(yPointing[1]*xPointing[0]-yPointing[0]*xPointing[1]);
      xposfit = fitpos[0]/xPointing[0] - xPointing[1]/xPointing[0]*yposfit;
  return 1;
}
