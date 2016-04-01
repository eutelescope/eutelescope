#include "EUTelProcessorAnalysisPALPIDEfsNoise.h"
#include "EUTELESCOPE.h"
#include "EUTelTrackerDataInterfacerImpl.h"
#include "EUTelGenericSparsePixel.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

#include "marlin/Global.h"

using namespace lcio;
using namespace marlin;
using namespace std;
using namespace eutelescope;

EUTelProcessorAnalysisPALPIDEfsNoise aEUTelProcessorAnalysisPALPIDEfsNoise;

EUTelProcessorAnalysisPALPIDEfsNoise::EUTelProcessorAnalysisPALPIDEfsNoise()
: Processor("EUTelProcessorAnalysisPALPIDEfsNoise"),
  _zsDataCollectionName(""),
  _fillHistos(false),
  _nEvent(0),
  _nFiredPixel(),
  _nLayer(0),
  _xPixel(),
  _yPixel(),
  _energy(6.0),
  _chipID(),
  _irradiation(),
  _dutIDs(),
  _rate(""),
  _outputSettingsFolderName("./")
  {
    _description="Ananlysis of noise runs";
    registerInputCollection (LCIO::TRACKERDATA, "ZSDataCollectionName",
                           "Input of Zero Suppressed data",
                           _zsDataCollectionName, string ("zsdata") );
    registerProcessorParameter("HistogramFilling","Switch on or off the histogram filling",                     _fillHistos, static_cast< bool > ( true ) );
    registerOptionalParameter("Energy","Particle energy",
                            _energy, static_cast< double > ( 6.0 ) );
    EVENT::StringVec _stringVecExample;
   _stringVecExample.push_back(" ");
    registerOptionalParameter("ChipID","Chip IDs",
                            _chipID, _stringVecExample );
    registerOptionalParameter("Irradiation","Irradiation level",
                            _irradiation, _stringVecExample );
    registerOptionalParameter("Rate","Data taking rate",
                            _rate, static_cast< string > ( "" ) );
    registerOptionalParameter("dutIDs","DUT IDs",
                            _dutIDs, _stringVecExample );
    registerOptionalParameter("OutputSettingsFolderName","Folder name where all the settings of each run will be saved",
                            _outputSettingsFolderName, static_cast< string > ( "./" ) );
    _isFirstEvent = true;
  }

void EUTelProcessorAnalysisPALPIDEfsNoise::init() {

  _nLayer = geo::gGeometry().nPlanes();
  vector<int> tmp(4,0);
  for (int iLayer=0; iLayer<_nLayer; iLayer++ )
  {
    _xPixel.push_back(geo::gGeometry().siPlaneXNpixels(iLayer));
    _yPixel.push_back(geo::gGeometry().siPlaneYNpixels(iLayer));
    _nFiredPixel.push_back(tmp);
//    cerr << iLayer << "\t" << nFiredPixel[0][iLayer] << endl;
  }
  for (unsigned int i=0; i<_dutIDs.size(); i++)
  {
    bool newFile = false;
    string _outputSettingsFileName = _outputSettingsFolderName + "settings_DUT" + _dutIDs[i] + ".txt";
    if (!std::ifstream(_outputSettingsFileName.c_str()))
      newFile = true;
    settingsFile[i].open (_outputSettingsFileName.c_str(), ios::out | ios::app );
    if (newFile) settingsFile[i] << "Run number;Energy;Chip ID;Irradiation level(0-nonIrradiated,1-2.5e12,2-1e13,3-700krad,4-combined:1e13+700krad);Rate;BB;Ithr;Idb;Vcasn;Vaux;Vcasp;Vreset;Threshold and their RMS for all four sectors;Noise and their RMS for all four sectors;Readout delay;Trigger delay;Strobe length;StrobeB length;Data (1) or noise (0);Number of events;Efficiency,Number of tracks,Number of tracks with associated hit for all sectors" << endl;
  }
}

void EUTelProcessorAnalysisPALPIDEfsNoise::processEvent(LCEvent *evt)
{
//  cerr << evt->getEventNumber() << endl;
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  if (_isFirstEvent )
  {
    for (unsigned int i=0; i< _dutIDs.size(); i++)
    {
      int dutID = atoi(_dutIDs[i].c_str());
      settingsFile[i] << evt->getRunNumber() << ";" << _energy << ";" << _chipID[dutID] << ";" << _irradiation[dutID] << ";" << _rate << ";" << evt->getParameters().getFloatVal("BackBiasVoltage") << ";" << evt->getParameters().getIntVal(Form("Ithr_%d",dutID)) << ";" << evt->getParameters().getIntVal(Form("Idb_%d",dutID)) << ";" << evt->getParameters().getIntVal(Form("Vcasn_%d",dutID)) << ";" << evt->getParameters().getIntVal(Form("Vaux_%d",dutID)) << ";" << evt->getParameters().getIntVal(Form("Vcasp_%d",dutID)) << ";" << evt->getParameters().getIntVal(Form("Vreset_%d",dutID)) << ";";
      for (int iSector=0; iSector<4; iSector++)
        settingsFile[i] << evt->getParameters().getFloatVal(Form("Thr_%d_%d",dutID,iSector)) << ";" << evt->getParameters().getFloatVal(Form("ThrRMS_%d_%d",dutID,iSector)) << ";";
      for (int iSector=0; iSector<4; iSector++)
        settingsFile[i] << evt->getParameters().getFloatVal(Form("Noise_%d_%d",dutID,iSector)) << ";" << evt->getParameters().getFloatVal(Form("NoiseRMS_%d_%d",dutID,iSector)) << ";";
      settingsFile[i] << evt->getParameters().getIntVal(Form("m_readout_delay_%d",dutID)) << ";" << evt->getParameters().getIntVal(Form("m_trigger_delay_%d",dutID)) << ";" << evt->getParameters().getIntVal(Form("m_strobe_length_%d",dutID)) << ";" << evt->getParameters().getIntVal(Form("m_strobeb_length_%d",dutID)) << ";0;";
    }
    if (_fillHistos)
      bookHistos();
    _isFirstEvent = false;
  }
  timeStampHisto->Fill(evt->getTimeStamp());
#endif
  try
  {
    zsInputDataCollectionVec = dynamic_cast< LCCollectionVec * > ( evt->getCollection( _zsDataCollectionName ) ) ;
  } catch ( lcio::DataNotAvailableException )
  {
    cerr << "In event " << evt->getEventNumber() << "_zsDataCollectionName " << _zsDataCollectionName.c_str() << " not found " << endl;
    return;
  }
  _nEvent++;
  for ( unsigned int iDetector = 0 ; iDetector < zsInputDataCollectionVec->size(); iDetector++ )
  {
    TrackerDataImpl * zsData = dynamic_cast< TrackerDataImpl * > ( zsInputDataCollectionVec->getElementAt( iDetector ) );
    auto sparseData = std::make_unique<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>>(zsData);
    auto& pixelVec = sparseData->getPixels();

    for( auto& sparsePixel: pixelVec ) {
      noiseMap[iDetector]->Fill(sparsePixel.getXCoord(),sparsePixel.getYCoord());
      for (int iSector=0; iSector<4; iSector++)
//      {
//        cerr << iSector*_xPixel[iDetector]/4 << "\t" <<(iSector+1)*_xPixel[iDetector]/4 << endl;
        if (sparsePixel.getXCoord() >= iSector*_xPixel[iDetector]/4 && sparsePixel.getXCoord() < (iSector+1)*_xPixel[iDetector]/4 )
          _nFiredPixel[iDetector][iSector]++;
//      }
//      cerr << evt->getEventNumber() << "\t" << iDetector << "\t" << sparsePixel->getXCoord() << "\t" << sparsePixel->getYCoord() << endl;
    }
  }
}

void EUTelProcessorAnalysisPALPIDEfsNoise::bookHistos()
{
  timeStampHisto = new TH1I("timeStampHisto","Distribution of the time stamp of the events; Time stamp (in 12.5 ns units)",1000,0,50000);
  for (int iLayer=0; iLayer<_nLayer; iLayer++ )
  {
    noiseMap[iLayer] = new TH2I(Form("noiseMap_%d",iLayer),Form("Noise map of layer %d",iLayer),_xPixel[iLayer],0,_xPixel[iLayer],_yPixel[iLayer],0,_yPixel[iLayer]);
    noiseOccupancy[iLayer] = new TH1F(Form("noiseOccupancy_%d",iLayer),Form("Noise occupancy in layer %d",iLayer),4,0,4);
  }
}

void EUTelProcessorAnalysisPALPIDEfsNoise::end()
{
  cout << "Total number of events: " << _nEvent << endl;
  for (unsigned int i=0; i< _dutIDs.size(); i++)
    settingsFile[i] << _nEvent << ";0;0;0;0;0;0;0;0;0;0;0;0" << endl;
  for (int iLayer=0; iLayer<_nLayer; iLayer++ )
  {
    for (int iSector=0; iSector<4; iSector++)
    {
//      cerr << "Total number of fired pixels in layer " << iLayer << ", sector " << iSector << " is " << _nFiredPixel[iLayer][iSector] << endl;
//      cerr << "Noise occupancy in layer " << iLayer << ", sector " << iSector << " is " << (double)_nFiredPixel[iLayer][iSector]/_nEvent/(_xPixel[iLayer]/4*_yPixel[iLayer]) << endl;
      noiseOccupancy[iLayer]->SetBinContent(iSector+1,(double)_nFiredPixel[iLayer][iSector]/_nEvent/(_xPixel[iLayer]/4*_yPixel[iLayer]));
      noiseOccupancy[iLayer]->SetBinError(iSector+1,sqrt((double)_nFiredPixel[iLayer][iSector])/_nEvent/(_xPixel[iLayer]/4*_yPixel[iLayer]));
    }
  }
}
