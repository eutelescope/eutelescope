// Version: $Id$
#include "EUTelX0Processor.h"
#include <cmath>
#include <cstdlib>
using namespace marlin;
using namespace eutelescope;
using namespace std;

std::string EUTelX0Processor::_histoResidualX = "ResidualX";
std::string EUTelX0Processor::_histoResidualXPlane1 = "ResidualXPlane1";
std::string EUTelX0Processor::_histoResidualXPlane2 = "ResidualXPlane2";
std::string EUTelX0Processor::_histoResidualXPlane3 = "ResidualXPlane3";
std::string EUTelX0Processor::_histoResidualXPlane4 = "ResidualXPlane4";
std::string EUTelX0Processor::_histoResidualYPlane1 = "ResidualYPlane1";
std::string EUTelX0Processor::_histoResidualYPlane2 = "ResidualYPlane2";
std::string EUTelX0Processor::_histoResidualYPlane3 = "ResidualYPlane3";
std::string EUTelX0Processor::_histoResidualYPlane4 = "ResidualYPlane4";
std::string EUTelX0Processor::_histoResidualY = "ResidualY";
std::string EUTelX0Processor::_histoResidualXZ = "ResidualXZ";
std::string EUTelX0Processor::_histoResidualYZ = "ResidualYZ";
std::string EUTelX0Processor::_histoResidualXY = "ResidualYZ";

EUTelX0Processor::EUTelX0Processor()
  :Processor("EUTelX0Processor"),
  _trackColName(""),
  _cutValue1(0.0),
  _cutValue2(0.0),
  _debug(false),
  _debugCount(0),
  _eventNumber(0),
  _finalEvent(false),
  _histoData(),
  _histoFile(""),
  _histoThing(),
  _hitInfo(),
  _inputHitCollectionVec(NULL), 
  _inputTrackCollectionVec(NULL), 
  _inputHitColName(""),
  _inputHitCollectionName(""),
  _inputTrackColName(""),
  _maxRecords(0),
  _projectedHits(),
  _referenceHitCollectionName(""),
  _referenceHitVec(NULL),
  _residual(),
  _residualAngle(),
  _residualCut(0.0),
  _residualProfile(),
  _runNumber(0),
  _trackCollectionName("")
{
  streamlog_out(DEBUG1) << "Constructing the EUTelX0Processor, setting all values to zero or NULL" << std::endl;
  _inputHitCollectionVec = new LCCollectionVec(LCIO::TRACKERHIT);//Used to store the values of the hit events
  _inputTrackCollectionVec = new LCCollectionVec(LCIO::TRACK);//Used to store the values of the hit events
  registerInputCollection(LCIO::TRACKERHIT,"InputTrackCollectionName",
                           "Collection name for corrected particle positions",
                           _trackColName, string ("alignedHit"));

  
  registerInputCollection(LCIO::TRACK,"OutputTrackCollectionName",
                           "Collection name for fitted tracks",
                           _trackCollectionName, string ("AlignedTrack"));

  registerOptionalParameter("ReferenceCollection","This is the name of the reference it collection (init at 0,0,0)", _referenceHitCollectionName, static_cast< string > ( "referenceHit" ) );//Necessary for working out which layer the particle is detected in
  registerProcessorParameter("ResidualCutValue","Used to determine cuts in the system, measured in XXX", _residualCut, static_cast< double > (50000.0));
  registerProcessorParameter("MaxRecords","Will be used to determine the final event if the final event must come before EOF", _maxRecords, static_cast< int > (0));
  registerProcessorParameter("HistoFile","Will be used to add the gaussian to the kink angle at the end of the run", _histoFile, static_cast< std::string > (""));
}

void EUTelX0Processor::init()
{
  streamlog_out(DEBUG5) << "Running EUTelX0Processor::init()" << std::endl;
  int nobins = 10000, nobinsangle = 100;//Number of bins in the histograms
  double minbin = -0.1, maxbin = 0.1;//Maximum and minimum bin values
  double minbinangle = -0.01, maxbinangle = 0.01, minbinalpha = -0.01, maxbinalpha = 0.01;
  std::vector<double> empty;  
   
  AIDA::IHistogram1D * SinglePointResidualXPlane0 = AIDAProcessor::histogramFactory(this)->createHistogram1D("SinglePointResidualXPlane0",nobins,minbin,maxbin);
  SinglePointResidualXPlane0->setTitle("SinglePointResidualXPlane0");
  _histoThing.insert(make_pair("SinglePointResidualXPlane0",SinglePointResidualXPlane0));
  _histoData["SinglePointResidualXPlane0"] = empty;
  
  AIDA::IHistogram1D * SinglePointResidualXPlane1 = AIDAProcessor::histogramFactory(this)->createHistogram1D("SinglePointResidualXPlane1",nobins,minbin,maxbin);
  SinglePointResidualXPlane1->setTitle("SinglePointResidualXPlane1");
  _histoThing.insert(make_pair("SinglePointResidualXPlane1",SinglePointResidualXPlane1));
  _histoData["SinglePointResidualXPlane1"] = empty;
 
  AIDA::IHistogram1D * SinglePointResidualXPlane2 = AIDAProcessor::histogramFactory(this)->createHistogram1D("SinglePointResidualXPlane2",nobins,minbin,maxbin);
  SinglePointResidualXPlane2->setTitle("SinglePointResidualXPlane2");
  _histoThing.insert(make_pair("SinglePointResidualXPlane2",SinglePointResidualXPlane2));
  _histoData["SinglePointResidualXPlane2"] = empty;
 
  AIDA::IHistogram1D * SinglePointResidualXPlane3 = AIDAProcessor::histogramFactory(this)->createHistogram1D("SinglePointResidualXPlane3",nobins,minbin,maxbin);
  SinglePointResidualXPlane3->setTitle("SinglePointResidualXPlane3");
  _histoThing.insert(make_pair("SinglePointResidualXPlane3",SinglePointResidualXPlane3));
  _histoData["SinglePointResidualXPlane3"] = empty;
 
  AIDA::IHistogram1D * SinglePointResidualXPlane4 = AIDAProcessor::histogramFactory(this)->createHistogram1D("SinglePointResidualXPlane4",nobins,minbin,maxbin);
  SinglePointResidualXPlane4->setTitle("SinglePointResidualXPlane4");
  _histoThing.insert(make_pair("SinglePointResidualXPlane4",SinglePointResidualXPlane4));
  _histoData["SinglePointResidualXPlane4"] = empty;
   
  AIDA::IHistogram1D * SinglePointResidualXPlane5 = AIDAProcessor::histogramFactory(this)->createHistogram1D("SinglePointResidualXPlane5",nobins,minbin,maxbin);
  SinglePointResidualXPlane5->setTitle("SinglePointResidualXPlane5");
  _histoThing.insert(make_pair("SinglePointResidualXPlane5",SinglePointResidualXPlane5));
  _histoData["SinglePointResidualXPlane5"] = empty;
   
  AIDA::IHistogram1D * SinglePointResidualYPlane0 = AIDAProcessor::histogramFactory(this)->createHistogram1D("SinglePointResidualXPlane0",nobins,minbin,maxbin);
  SinglePointResidualXPlane0->setTitle("SinglePointResidualXPlane0");
  _histoThing.insert(make_pair("SinglePointResidualXPlane0",SinglePointResidualXPlane0));
  _histoData["SinglePointResidualXPlane0"] = empty;
 
  AIDA::IHistogram1D * SinglePointResidualYPlane1 = AIDAProcessor::histogramFactory(this)->createHistogram1D("SinglePointResidualYPlane1",nobins,minbin,maxbin);
  SinglePointResidualYPlane1->setTitle("SinglePointResidualYPlane1");
  _histoThing.insert(make_pair("SinglePointResidualYPlane1",SinglePointResidualYPlane1));
  _histoData["SinglePointResidualYPlane1"] = empty;
 
  AIDA::IHistogram1D * SinglePointResidualYPlane2 = AIDAProcessor::histogramFactory(this)->createHistogram1D("SinglePointResidualYPlane2",nobins,minbin,maxbin);
  SinglePointResidualYPlane2->setTitle("SinglePointResidualYPlane2");
  _histoThing.insert(make_pair("SinglePointResidualYPlane2",SinglePointResidualYPlane2));
  _histoData["SinglePointResidualYPlane2"] = empty;
 
  AIDA::IHistogram1D * SinglePointResidualYPlane3 = AIDAProcessor::histogramFactory(this)->createHistogram1D("SinglePointResidualYPlane3",nobins,minbin,maxbin);
  SinglePointResidualYPlane3->setTitle("SinglePointResidualYPlane3");
  _histoThing.insert(make_pair("SinglePointResidualYPlane3",SinglePointResidualYPlane3));
  _histoData["SinglePointResidualYPlane3"] = empty;
 
  AIDA::IHistogram1D * SinglePointResidualYPlane4 = AIDAProcessor::histogramFactory(this)->createHistogram1D("SinglePointResidualYPlane4",nobins,minbin,maxbin);
  SinglePointResidualYPlane4->setTitle("SinglePointResidualYPlane4");
  _histoThing.insert(make_pair("SinglePointResidualYPlane4",SinglePointResidualYPlane4));
  _histoData["SinglePointResidualYPlane4"] = empty;
   
  AIDA::IHistogram1D * SinglePointResidualYPlane5 = AIDAProcessor::histogramFactory(this)->createHistogram1D("SinglePointResidualXPlane5",nobins,minbin,maxbin);
  SinglePointResidualXPlane5->setTitle("SinglePointResidualXPlane5");
  _histoThing.insert(make_pair("SinglePointResidualXPlane5",SinglePointResidualXPlane5));
  _histoData["SinglePointResidualXPlane5"] = empty;
 
  AIDA::IHistogram1D * ThreePointResidualXPlane1 = AIDAProcessor::histogramFactory(this)->createHistogram1D("ThreePointResidualXPlane1",nobins,minbin,maxbin);
  ThreePointResidualXPlane1->setTitle("ThreePointResidualXPlane1");
  _histoThing.insert(make_pair("ThreePointResidualXPlane1",ThreePointResidualXPlane1));
  _histoData["ThreePointResidualXPlane1"] = empty;
 
  AIDA::IHistogram1D * ThreePointResidualXPlane2 = AIDAProcessor::histogramFactory(this)->createHistogram1D("ThreePointResidualXPlane2",nobins,minbin,maxbin);
  ThreePointResidualXPlane2->setTitle("ThreePointResidualXPlane2");
  _histoThing.insert(make_pair("ThreePointResidualXPlane2",ThreePointResidualXPlane2));
  _histoData["ThreePointResidualXPlane2"] = empty;
 
  AIDA::IHistogram1D * ThreePointResidualXPlane3 = AIDAProcessor::histogramFactory(this)->createHistogram1D("ThreePointResidualXPlane3",nobins,minbin,maxbin);
  ThreePointResidualXPlane3->setTitle("ThreePointResidualXPlane3");
  _histoThing.insert(make_pair("ThreePointResidualXPlane3",ThreePointResidualXPlane3));
  _histoData["ThreePointResidualXPlane3"] = empty;
 
  AIDA::IHistogram1D * ThreePointResidualXPlane4 = AIDAProcessor::histogramFactory(this)->createHistogram1D("ThreePointResidualXPlane4",nobins,minbin,maxbin);
  ThreePointResidualXPlane4->setTitle("ThreePointResidualXPlane4");
  _histoThing.insert(make_pair("ThreePointResidualXPlane4",ThreePointResidualXPlane4));
  _histoData["ThreePointResidualXPlane4"] = empty;
 
  AIDA::IHistogram1D * ThreePointResidualYPlane1 = AIDAProcessor::histogramFactory(this)->createHistogram1D("ThreePointResidualYPlane1",nobins,minbin,maxbin);
  ThreePointResidualYPlane1->setTitle("ThreePointResidualYPlane1");
  _histoThing.insert(make_pair("ThreePointResidualYPlane1",ThreePointResidualYPlane1));
  _histoData["ThreePointResidualYPlane1"] = empty;
 
  AIDA::IHistogram1D * ThreePointResidualYPlane2 = AIDAProcessor::histogramFactory(this)->createHistogram1D("ThreePointResidualYPlane2",nobins,minbin,maxbin);
  ThreePointResidualYPlane2->setTitle("ThreePointResidualYPlane2");
  _histoThing.insert(make_pair("ThreePointResidualYPlane2",ThreePointResidualYPlane2));
  _histoData["ThreePointResidualYPlane2"] = empty;
 
  AIDA::IHistogram1D * ThreePointResidualYPlane3 = AIDAProcessor::histogramFactory(this)->createHistogram1D("ThreePointResidualYPlane3",nobins,minbin,maxbin);
  ThreePointResidualYPlane3->setTitle("ThreePointResidualYPlane3");
  _histoThing.insert(make_pair("ThreePointResidualYPlane3",ThreePointResidualYPlane3));
  _histoData["ThreePointResidualYPlane3"] = empty;
 
  AIDA::IHistogram1D * ThreePointResidualYPlane4 = AIDAProcessor::histogramFactory(this)->createHistogram1D("ThreePointResidualYPlane4",nobins,minbin,maxbin);
  ThreePointResidualYPlane4->setTitle("ThreePointResidualYPlane4");
  _histoThing.insert(make_pair("ThreePointResidualYPlane4",ThreePointResidualYPlane4));
  _histoData["ThreePointResidualYPlane4"] = empty;
 
  AIDA::IHistogram1D * AngleXForwardPlane0 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Angle X Forward Plane 0",nobinsangle,minbinangle,maxbinangle);
  AngleXForwardPlane0->setTitle("Angle X Forward Plane 0");
  _histoThing.insert(make_pair("Angle X Forward Plane 0",AngleXForwardPlane0));
  _histoData["Angle X Forward Plane 0"] = empty;

  AIDA::IHistogram1D * AngleXForwardPlane1 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Angle X Forward Plane 1",nobinsangle,minbinangle,maxbinangle);
  AngleXForwardPlane1->setTitle("Angle X Forward Plane 1");
  _histoThing.insert(make_pair("Angle X Forward Plane 1",AngleXForwardPlane1));
  _histoData["Angle X Forward Plane 1"] = empty;

  AIDA::IHistogram1D * AngleXForwardPlane2 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Angle X Forward Plane 2",nobinsangle,minbinangle,maxbinangle);
  AngleXForwardPlane2->setTitle("Angle X Forward Plane 2");
  _histoThing.insert(make_pair("Angle X Forward Plane 2",AngleXForwardPlane2));
  _histoData["Angle X Forward Plane 2"] = empty;

  AIDA::IHistogram1D * AngleXForwardPlane3 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Angle X Forward Plane 3",nobinsangle,minbinangle,maxbinangle);
  AngleXForwardPlane3->setTitle("Angle X Forward Plane 3");
  _histoThing.insert(make_pair("Angle X Forward Plane 3",AngleXForwardPlane3));
  _histoData["Angle X Forward Plane 3"] = empty;
 
  AIDA::IHistogram1D * AngleXForwardPlane4 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Angle X Forward Plane 4",nobinsangle,minbinangle,maxbinangle);
  AngleXForwardPlane4->setTitle("Angle X Forward Plane 4");
  _histoThing.insert(make_pair("Angle X Forward Plane 4",AngleXForwardPlane4));
  _histoData["Angle X Forward Plane 4"] = empty;

  AIDA::IHistogram1D * AngleYForwardPlane0 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Angle Y Forward Plane 0",nobinsangle,minbinangle,maxbinangle);
  AngleYForwardPlane0->setTitle("Angle Y Forward Plane 0");
  _histoThing.insert(make_pair("Angle Y Forward Plane 0",AngleYForwardPlane0));
  _histoData["Angle Y Forward Plane 0"] = empty;

  AIDA::IHistogram1D * AngleYForwardPlane1 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Angle Y Forward Plane 1",nobinsangle,minbinangle,maxbinangle);
  AngleYForwardPlane1->setTitle("Angle Y Forward Plane 1");
  _histoThing.insert(make_pair("Angle Y Forward Plane 1",AngleYForwardPlane1));
  _histoData["Angle Y Forward Plane 1"] = empty;

  AIDA::IHistogram1D * AngleYForwardPlane2 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Angle Y Forward Plane 2",nobinsangle,minbinangle,maxbinangle);
  AngleYForwardPlane2->setTitle("Angle Y Forward Plane 2");
  _histoThing.insert(make_pair("Angle Y Forward Plane 2",AngleYForwardPlane2));
  _histoData["Angle Y Forward Plane 2"] = empty;

  AIDA::IHistogram1D * AngleYForwardPlane3 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Angle Y Forward Plane 3",nobinsangle,minbinangle,maxbinangle);
  AngleYForwardPlane3->setTitle("Angle Y Forward Plane 3");
  _histoThing.insert(make_pair("Angle Y Forward Plane 3",AngleYForwardPlane3));
  _histoData["Angle Y Forward Plane 3"] = empty;
 
  AIDA::IHistogram1D * AngleYForwardPlane4 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Angle Y Forward Plane 4",nobinsangle,minbinangle,maxbinangle);
  AngleYForwardPlane4->setTitle("Angle Y Forward Plane 4");
  _histoThing.insert(make_pair("Angle Y Forward Plane 4",AngleYForwardPlane4));
  _histoData["Angle Y Forward Plane 4"] = empty;

  AIDA::IHistogram1D * AngleXBackwardPlane1 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Angle X Backward Plane 1",nobinsangle,minbinangle,maxbinangle);
  AngleXBackwardPlane1->setTitle("Angle X Backward Plane 1");
  _histoThing.insert(make_pair("Angle X Backward Plane 1",AngleXBackwardPlane1));
  _histoData["Angle X Backward Plane 1"] = empty;

  AIDA::IHistogram1D * AngleXBackwardPlane2 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Angle X Backward Plane 2",nobinsangle,minbinangle,maxbinangle);
  AngleXBackwardPlane2->setTitle("Angle X Backward Plane 2");
  _histoThing.insert(make_pair("Angle X Backward Plane 2",AngleXBackwardPlane2));
  _histoData["Angle X Backward Plane 2"] = empty;

  AIDA::IHistogram1D * AngleXBackwardPlane3 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Angle X Backward Plane 3",nobinsangle,minbinangle,maxbinangle);
  AngleXBackwardPlane3->setTitle("Angle X Backward Plane 3");
  _histoThing.insert(make_pair("Angle X Backward Plane 3",AngleXBackwardPlane3));
  _histoData["Angle X Backward Plane 3"] = empty;
 
  AIDA::IHistogram1D * AngleXBackwardPlane4 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Angle X Backward Plane 4",nobinsangle,minbinangle,maxbinangle);
  AngleXBackwardPlane4->setTitle("Angle X Backward Plane 4");
  _histoThing.insert(make_pair("Angle X Backward Plane 4",AngleXBackwardPlane4));
  _histoData["Angle X Backward Plane 4"] = empty;

  AIDA::IHistogram1D * AngleXBackwardPlane5 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Angle X Backward Plane 5",nobinsangle,minbinangle,maxbinangle);
  AngleXBackwardPlane5->setTitle("Angle X Backward Plane 5");
  _histoThing.insert(make_pair("Angle X Backward Plane 5",AngleXBackwardPlane5));
  _histoData["Angle X Backward Plane 5"] = empty;

  AIDA::IHistogram1D * AngleYBackwardPlane1 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Angle Y Backward Plane 1",nobinsangle,minbinangle,maxbinangle);
  AngleYBackwardPlane1->setTitle("Angle Y Backward Plane 1");
  _histoThing.insert(make_pair("Angle Y Backward Plane 1",AngleYBackwardPlane1));
  _histoData["Angle Y Backward Plane 1"] = empty;

  AIDA::IHistogram1D * AngleYBackwardPlane2 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Angle Y Backward Plane 2",nobinsangle,minbinangle,maxbinangle);
  AngleYBackwardPlane2->setTitle("Angle Y Backward Plane 2");
  _histoThing.insert(make_pair("Angle Y Backward Plane 2",AngleYBackwardPlane2));
  _histoData["Angle Y Backward Plane 2"] = empty;

  AIDA::IHistogram1D * AngleYBackwardPlane3 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Angle Y Backward Plane 3",nobinsangle,minbinangle,maxbinangle);
  AngleYBackwardPlane3->setTitle("Angle Y Backward Plane 3");
  _histoThing.insert(make_pair("Angle Y Backward Plane 3",AngleYBackwardPlane3));
  _histoData["Angle Y Backward Plane 3"] = empty;
 
  AIDA::IHistogram1D * AngleYBackwardPlane4 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Angle Y Backward Plane 4",nobinsangle,minbinangle,maxbinangle);
  AngleYBackwardPlane4->setTitle("Angle Y Backward Plane 4");
  _histoThing.insert(make_pair("Angle Y Backward Plane 4",AngleYBackwardPlane4));
  _histoData["Angle Y Backward Plane 4"] = empty;

  AIDA::IHistogram1D * AngleYBackwardPlane5 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Angle Y Backward Plane 5",nobinsangle,minbinangle,maxbinangle);
  AngleYBackwardPlane5->setTitle("Angle Y Backward Plane 5");
  _histoThing.insert(make_pair("Angle Y Backward Plane 5",AngleYBackwardPlane5));
  _histoData["Angle Y Backward Plane 5"] = empty;

  AIDA::IHistogram1D * ScatteringAngleXPlane1 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Scattering Angle X Plane 1",nobinsangle,minbinangle,maxbinangle);
  ScatteringAngleXPlane1->setTitle("Scattering Angle X Plane 1");
  _histoThing.insert(make_pair("Scattering Angle X Plane 1",ScatteringAngleXPlane1));
  _histoData["Scattering Angle X Plane 1"] = empty;

  AIDA::IHistogram1D * ScatteringAngleXPlane2 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Scattering Angle X Plane 2",nobinsangle,minbinangle,maxbinangle);
  ScatteringAngleXPlane2->setTitle("Scattering Angle X Plane 2");
  _histoThing.insert(make_pair("Scattering Angle X Plane 2",ScatteringAngleXPlane2));
  _histoData["Scattering Angle X Plane 2"] = empty;

  AIDA::IHistogram1D * ScatteringAngleXPlane3 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Scattering Angle X Plane 3",nobinsangle,minbinangle,maxbinangle);
  ScatteringAngleXPlane3->setTitle("Scattering Angle X Plane 3");
  _histoThing.insert(make_pair("Scattering Angle X Plane 3",ScatteringAngleXPlane3));
  _histoData["Scattering Angle X Plane 3"] = empty;
 
  AIDA::IHistogram1D * ScatteringAngleXPlane4 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Scattering Angle X Plane 4",nobinsangle,minbinangle,maxbinangle);
  ScatteringAngleXPlane4->setTitle("Scattering Angle X Plane 4");
  _histoThing.insert(make_pair("Scattering Angle X Plane 4",ScatteringAngleXPlane4));
  _histoData["Scattering Angle X Plane 4"] = empty;

  AIDA::IHistogram1D * ScatteringAngleYPlane1 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Scattering Angle Y Plane 1",nobinsangle,minbinangle,maxbinangle);
  ScatteringAngleYPlane1->setTitle("Scattering Angle Y Plane 1");
  _histoThing.insert(make_pair("Scattering Angle Y Plane 1",ScatteringAngleYPlane1));
  _histoData["Scattering Angle Y Plane 1"] = empty;

  AIDA::IHistogram1D * ScatteringAngleYPlane2 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Scattering Angle Y Plane 2",nobinsangle,minbinangle,maxbinangle);
  ScatteringAngleYPlane2->setTitle("Scattering Angle Y Plane 2");
  _histoThing.insert(make_pair("Scattering Angle Y Plane 2",ScatteringAngleYPlane2));
  _histoData["Scattering Angle Y Plane 2"] = empty;

  AIDA::IHistogram1D * ScatteringAngleYPlane3 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Scattering Angle Y Plane 3",nobinsangle,minbinangle,maxbinangle);
  ScatteringAngleYPlane3->setTitle("Scattering Angle Y Plane 3");
  _histoThing.insert(make_pair("Scattering Angle Y Plane 3",ScatteringAngleYPlane3));
  _histoData["Scattering Angle Y Plane 3"] = empty;
 
  AIDA::IHistogram1D * ScatteringAngleYPlane4 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Scattering Angle Y Plane 4",nobinsangle,minbinangle,maxbinangle);
  ScatteringAngleYPlane4->setTitle("Scattering Angle Y Plane 4");
  _histoThing.insert(make_pair("Scattering Angle Y Plane 4",ScatteringAngleYPlane4));
  _histoData["Scattering Angle Y Plane 4"] = empty;

  AIDA::IHistogram2D * KinkAnglePlane1 = AIDAProcessor::histogramFactory(this)->createHistogram2D("Kink Angle Plane 1",nobinsangle,minbinalpha,maxbinalpha,nobinsangle,minbinalpha,maxbinalpha);
  KinkAnglePlane1->setTitle("Kink Angle Plane 1");
  _histoThing.insert(make_pair("Kink Angle Plane 1",KinkAnglePlane1));
  _histoData["Kink Angle Plane 1"] = empty;

  AIDA::IHistogram2D * KinkAnglePlane2 = AIDAProcessor::histogramFactory(this)->createHistogram2D("Kink Angle Plane 2",nobinsangle,minbinalpha,maxbinalpha,nobinsangle,minbinalpha,maxbinalpha);
  KinkAnglePlane2->setTitle("Kink Angle Plane 2");
  _histoThing.insert(make_pair("Kink Angle Plane 2",KinkAnglePlane2));
  _histoData["Kink Angle Plane 2"] = empty;

  AIDA::IHistogram2D * KinkAnglePlane3 = AIDAProcessor::histogramFactory(this)->createHistogram2D("Kink Angle Plane 3",nobinsangle,minbinalpha,maxbinalpha,nobinsangle,minbinalpha,maxbinalpha);
  KinkAnglePlane3->setTitle("Kink Angle Plane 3");
  _histoThing.insert(make_pair("Kink Angle Plane 3",KinkAnglePlane3));
  _histoData["Kink Angle Plane 3"] = empty;

  AIDA::IHistogram2D * KinkAnglePlane4 = AIDAProcessor::histogramFactory(this)->createHistogram2D("Kink Angle Plane 4",nobinsangle,minbinalpha,maxbinalpha,nobinsangle,minbinalpha,maxbinalpha);
  KinkAnglePlane4->setTitle("Kink Angle Plane 4");
  _histoThing.insert(make_pair("Kink Angle Plane 4",KinkAnglePlane4));
  _histoData["Kink Angle Plane 4"] = empty;

  
}

void EUTelX0Processor::processRunHeader(LCRunHeader *run)
{
  streamlog_out(DEBUG1) << "Running EUTelX0Processor::processRunHeader(LCRunHeader *run) with run = " << run << std::endl;
  _eventNumber = 0;
  _runNumber++;
}

void EUTelX0Processor::printTrackParameters( Track* eventtrack ){
    streamlog_out(DEBUG0) << "Type: " << eventtrack->getType() << std::endl << std::endl
    << "***TRACK PARAMETERS***" << std::endl << "(Note: In EUTelTestFitter the track parameters are not yet instantiated, all currently set to zero)" << std::endl
    << "D0: " << eventtrack->getD0() << std::endl
    << "Omega: " << eventtrack->getOmega() << std::endl
    << "Phi: " << eventtrack->getPhi() << std::endl
    << "Z0: " << eventtrack->getZ0() << std::endl
    << "Tan Lambda: " << eventtrack->getTanLambda() << std::endl 
    << "Covariance Matrix of Track Parameters: (This is a square matrix with parameter order: D0, Phi, Omega, Z0, TanLambda)" << std::endl;
    std::vector< float > covmatrix = eventtrack->getCovMatrix();
    streamlog_out(DEBUG0) << covmatrix[0] << ", " << covmatrix[1] << ", " << covmatrix[3] << ", " << covmatrix[6] << ", " << covmatrix[10] << std::endl
                          << covmatrix[1] << ", " << covmatrix[2] << ", " << covmatrix[4] << ", " << covmatrix[7] << ", " << covmatrix[11] << std::endl
                          << covmatrix[3] << ", " << covmatrix[4] << ", " << covmatrix[5] << ", " << covmatrix[8] << ", " << covmatrix[12] << std::endl
                          << covmatrix[6] << ", " << covmatrix[7] << ", " << covmatrix[8] << ", " << covmatrix[9] << ", " << covmatrix[13] << std::endl
                          << covmatrix[10] << ", " << covmatrix[11] << ", " << covmatrix[12] << ", " << covmatrix[13] << ", " << covmatrix[14] << std::endl;    
    streamlog_out(DEBUG0) << std::endl
    << "***ADDITIONAL INFORMATION***" << std::endl
    << "Reference Point: " << eventtrack->getReferencePoint()[0] << ", " <<  eventtrack->getReferencePoint()[1] << ", " << eventtrack->getReferencePoint()[2] << std::endl
    << "Chi2: " << eventtrack->getChi2() << std::endl
    << "Degrees of Freedom: " << eventtrack->getNdf() << std::endl
    << "dE/dx: " << eventtrack->getdEdx() << std::endl
    << "dE/dx error: " << eventtrack->getdEdxError() << std::endl << std::endl;
  streamlog_out(DEBUG0) << "***TRACK HITS***" << std::endl;
  std::vector< TrackerHit* > trackhits = eventtrack->getTrackerHits();
  int tracknumber(0);
  for(std::vector< TrackerHit* >::iterator it = trackhits.begin(); it != trackhits.end(); ++it){
    streamlog_out(DEBUG0) << "Track Hit " << tracknumber << ":" << std::endl;
    tracknumber++;
    printHitParameters(*it);
  }
}

void EUTelX0Processor::printHitParameters( TrackerHit* trackhit ){
  streamlog_out(DEBUG0) << "Type: " << trackhit->getType() << std::endl
    << "Time of hit (in nanoseconds): " << trackhit->getTime() << std::endl
    << "dE/dx (in GeV): " << trackhit->getEDep() << std::endl
    << "Position (x, y, z and measured in mm): " << trackhit->getPosition()[0] << ", "
    << trackhit->getPosition()[1] << ", " << trackhit->getPosition()[2] << std::endl
    << "Covariance of the position (x, y, z): " << std::endl;
  std::vector< float > hitcov = trackhit->getCovMatrix(); //Size 6
  streamlog_out(DEBUG0) << hitcov[0] << ", " << hitcov[1] << ", " << hitcov[3] << std::endl
                        << hitcov[1] << ", " << hitcov[2] << ", " << hitcov[4] << std::endl
                        << hitcov[3] << ", " << hitcov[4] << ", " << hitcov[5] << std::endl << std::endl;
}

void EUTelX0Processor::kinkEstimate(Track* track){
  //This function works out an angle based on a straight line fitted from plane 0 to 2 and plane 5 to 3
  //It will also store all other angles in histograms too
  streamlog_out(DEBUG0) << "Running function kinkEstimate(Track* " << &track << ")" << std::endl;
  
  //First we extract the relevant hits from the track
  std::vector< TVector3* > hits = getHitsFromTrack(track);
  streamlog_out(DEBUG3) << "Successfully got hits from track" << std::endl;
  //Then we find the position where those lines would have hit the DUT
  //THIS IS TO BE DECIDED IF IT IS NEEDED LATER

  //Then we work out the angles of these lines with respect to XZ and YZ, plot results in histograms
  int hitsize = static_cast< int >(hits.size());
  std::vector< double > scatterx, scattery;
  for(int i = 0; i < hitsize-1; ++i){
    double x0 = hits[i]->x();
    double y0 = hits[i]->y();
    double z0 = hits[i]->z();
    double x1 = hits[i+1]->x();
    double y1 = hits[i+1]->y();
    double z1 = hits[i+1]->z();
    double deltaz = z1-z0;
    double frontanglex = atan2(x1-x0,deltaz);
    double frontangley = atan2(y1-y0,deltaz);
    double backanglex = atan2(x0-x1,deltaz);
    double backangley = atan2(y0-y1,deltaz);
    scatterx.push_back(frontanglex);
    scattery.push_back(frontangley);
    std::stringstream xforward,yforward,xbackward,ybackward;
    xforward << "Angle X Forward Plane " << i;
    yforward << "Angle Y Forward Plane " << i;
    xbackward << "Angle X Backward Plane " << i+1;
    ybackward << "Angle Y Backward Plane " << i+1;
    try{
      dynamic_cast< AIDA::IHistogram1D* >(_histoThing[xforward.str().c_str()])->fill(frontanglex);
    } catch(std::bad_cast &bc){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << xforward.str().c_str() << endl;
    }
    try{
      dynamic_cast< AIDA::IHistogram1D* >(_histoThing[yforward.str().c_str()])->fill(frontangley);
    } catch(std::bad_cast &bc){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << yforward.str().c_str() << endl;
    }
    try{
      dynamic_cast< AIDA::IHistogram1D* >(_histoThing[xbackward.str().c_str()])->fill(backanglex);
    } catch(std::bad_cast &bc){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << xbackward.str().c_str() << endl;
    }
    try{
      dynamic_cast< AIDA::IHistogram1D* >(_histoThing[ybackward.str().c_str()])->fill(backangley);
    } catch(std::bad_cast &bc){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << ybackward.str().c_str() << endl;
    }
  }
 
  int scatterxsize = static_cast< int >(scatterx.size());
  for(int i = 0; i < scatterxsize-1; ++i){
    double scatteringanglex = scatterx[i+1] - scatterx[i];
    double scatteringangley = scattery[i+1] - scattery[i];
    std::stringstream ssscatterx,ssscattery,ssscatterxy;
    ssscatterx << "Scattering Angle X Plane " << i+1;
    ssscattery << "Scattering Angle Y Plane " << i+1;
    ssscatterxy << "Kink Angle Plane " << i+1;
    try{
      dynamic_cast< AIDA::IHistogram1D* >(_histoThing[ssscatterx.str().c_str()])->fill(scatteringanglex);
    } catch(std::bad_cast &bc){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << ssscatterx.str().c_str() << endl;
    }
    try{
      dynamic_cast< AIDA::IHistogram1D* >(_histoThing[ssscattery.str().c_str()])->fill(scatteringangley);
    } catch(std::bad_cast &bc){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << ssscattery.str().c_str() << endl;
    }
    try{
      dynamic_cast< AIDA::IHistogram2D* >(_histoThing[ssscatterxy.str().c_str()])->fill(scatteringanglex,scatteringangley);
    } catch(std::bad_cast &bc){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << ssscatterxy.str().c_str() << endl;
    }
  }

  streamlog_out(DEBUG3) << "Made it to the end of kinkEstimate()" << endl;
}

void EUTelX0Processor::kinkGaussian(){
/*
  gStyle->SetOptFit(111);
  std::string _histoFile = "/scratch/hamnett/TestBeam/2013/data_X0MeasurementsOneWeek20012013/histo/000001-track-histo.root";
  TFile *file = new TFile(_histoFile.c_str(),"UPDATE");
  TDirectory *directory = (TDirectory*)file->GetDirectory("X0Scanner");
  directory->ls();
  TH1D *histogram = (TH1D*)directory->Get("KinkAngle");
  histogram->Draw();
  TF1 *fit = new TF1("fit","gaus",0,5);
  histogram->Fit("fit","R");
  file->Write();
*/
  
//  AIDA::IFitter * fit = AIDAProcessor::getIAnalysisFactory(this)->createFitFactory()->createFitter("Chi2")->fit(_histoThing["Kink Angle"],fitfunction);

}

void EUTelX0Processor::processEvent(LCEvent *evt)
{
  streamlog_out(MESSAGE0) << "Running EUTelX0Processor::processEvent(LCEvent *evt) with evt = " << evt << std::endl;
  //Take track from input parameter
  //Work out kink angle from track
  //Put kink angle into histogram
  //Check for last event
    //Fit gaussian to histogram
    //Extract sigma value from histogram
    //Deduce radiation length from sigma value - See Nadler and Fruwurth paper

  try{
    if(_eventNumber == _maxRecords - 2){ //Not sure about the -2
      _finalEvent = true;
      kinkGaussian();
    } else{
      //_referenceHitVec = evt->getCollection(_referenceHitCollectionName);
      LCCollection* trackcollection = evt->getCollection(_trackCollectionName);
      int elementnumber = trackcollection->getNumberOfElements();
      for(int i = 0; i < elementnumber; ++i){
        Track* eventtrack = dynamic_cast< Track* >(trackcollection->getElementAt(i));
        streamlog_out(DEBUG0) << "Here is all the information about the track in run " << _runNumber << ", event " << _eventNumber << ", element " << i << std::endl << std::endl;
//        printTrackParameters( eventtrack );
        kinkEstimate( eventtrack );
//        threePointResolution( eventtrack );
      }
    }
  } catch(DataNotAvailableException &datanotavailable){
    streamlog_out(WARNING4) << "Exception occured: " << datanotavailable.what() << std::endl
      << "Could not get collection '" << _trackCollectionName << "' from event " << _eventNumber << ", Skipping Event" << std::endl;
  } catch(...){
    streamlog_out(ERROR9) << "Unknown exception occured in EUTelX0Processor, here is some information which might help:" << std::endl
    << "Run Number: " << _runNumber << std::endl
    << "Event Number: " << _eventNumber << std::endl
    << "Throw occured somewhere within the try block which starts with the line: 'LCCollection* trackcollection = evt->getCollection(_trackCollectionName);'" << std::endl
    << "I hope that helps, if not then please submit a bug report or add more exception handling as appropriate, program will now exit" << std::endl;
    exit(1);
  }
  _eventNumber++;
}

int EUTelX0Processor::guessSensorID(const double * hit )
{
  int sensorID = -1;
  double minDistance = numeric_limits< double >::max() ;
  streamlog_out( DEBUG5 ) << "referencehit collection: " << _referenceHitCollectionName << " at "<< _referenceHitVec << std::endl;
  if( _referenceHitVec == 0)
  {
    streamlog_out( MESSAGE5 ) << "_referenceHitVec is empty" << std::endl;
    return 0;
  }

  if( isFirstEvent() )
  {
    // message<DEBUG > ( log() << "number of elements : " << _referenceHitVec->getNumberOfElements() << std::endl );
  }

  for(int ii = 0 ; ii < _referenceHitVec->getNumberOfElements(); ii++)
  {
    EUTelReferenceHit* refhit = static_cast< EUTelReferenceHit*> ( _referenceHitVec->getElementAt(ii) ) ;
    if(refhit == 0 ) continue;
     TVector3 hit3d( hit[0], hit[1], hit[2] );
    TVector3 hitInPlane( refhit->getXOffset(), refhit->getYOffset(), refhit->getZOffset());
    TVector3 norm2Plane( refhit->getAlpha(), refhit->getBeta(), refhit->getGamma() );
    double distance = abs( norm2Plane.Dot(hit3d-hitInPlane) );
    streamlog_out( DEBUG5 ) << " hit " << hit[0] << " "<< hit[1] << " " << hit[2] << std::endl;
    streamlog_out( DEBUG5 ) << " " << refhit->getXOffset() << " " << refhit->getYOffset() << " " << refhit->getZOffset() << std::endl;
    streamlog_out( DEBUG5 ) << " " << refhit->getAlpha() << " " << refhit->getBeta() << " " << refhit->getGamma() << std::endl;
    streamlog_out( DEBUG5 ) << " distance " << distance << std::endl;
    if ( distance < minDistance )
    {
      minDistance = distance;
      sensorID = refhit->getSensorID();
    }
  }
  //Some useful debug printouts:
  bool debug = ( _debugCount>0 );

  if(debug)
  {
    for(int ii = 0 ; ii < _referenceHitVec->getNumberOfElements(); ii++)
    {
      EUTelReferenceHit* refhit = static_cast< EUTelReferenceHit*> ( _referenceHitVec->getElementAt(ii) ) ;
      if(refhit == 0 ) continue;
      // if( sensorID != refhit->getSensorID() ) continue;
      // streamlog_out( DEBUG5 ) << " _referenceHitVec " << _referenceHitVec << " " << _referenceHitCollectionName.c_str() << " " << refhit << " at "
      // << refhit->getXOffset() << " " << refhit->getYOffset() << " " << refhit->getZOffset() << " "
      //<< refhit->getAlpha() << " " << refhit->getBeta() << " " << refhit->getGamma() << std::endl ;
      //message<DEBUG5> ( log() << "iPlane " << refhit->getSensorID() << " hitPos: [" << hit[0] << " " << hit[1] << " " << hit[2] << "] distance: " << minDistance << std::endl );
    }
  }
  return sensorID;
}
 
void EUTelX0Processor::end()
{
  //Clean up memory
  //Set all values to zero or NULL
  //


  //calculateX0();
  _hitInfo.clear();  //This stores the hit position in a TVector3. If there are multiple hits then they are all stored in the vector of TVector3's. The int key refers to the layer number
  _projectedHits.clear();  //int refers to 1 or 2, with 1 being the projection from the 01 plane and 2 being the projection from the 43 plane
  _histoThing.clear();  //This map contains all the histograms that are going to be plotted in this processor
  _residual.clear();  //The pair of doubles refers to the x and y coordinates in layer 2. The vector of TVector3's refers to the positions of the projected hits
  _residualAngle.clear();  //As above but the TVector3 doesn't contain xyz but instead theta, phi and alpha
  _residualProfile.clear(); //TODO(Phillip Hamnett): Can this be joined with _residual? //Used as above but for created a profile histogram
  _inputHitCollectionVec = NULL;  //Stores the information being brought in from the Marlin process which contains information about post-aligned hits
  _inputTrackCollectionVec = NULL;  //Stores the information being brought in from the Marlin process which contains information about post-aligned hits
  _referenceHitVec = NULL;  
}

double EUTelX0Processor::calculateX0()
{
  double X0(0);
  //Two methods for figuring out the material budget, X0, are used. The loss of energy and the scattering angle due to multiple scattering.
  //Potentially might use both a Forward-backward Kalmin filter and a Global Linear Estimator to work out the two above values, might give a more accurate meaurement of X0.
  double sigma_msi(0);  //Highland formula
  double sigma_measi(0);  //Sum of the variance of histogram alpha012 and alpha432
  double E = 5;  //The total particle energy (GeV) TODO: (Phillip Hamnett) is there some way to make this information (and the mass) into the steering file?
  double m = 0.000511;  //Electron Mass (GeV)
  double p = sqrt(E*E - m*m);  //The particle momentum (GeV)
  double theta012(0);  //The polar angle between the momentum vector and the z-axis
  double beta012(0);  //The difference between the azimuthal angle of the momentum vector and the azimuthal angle of a point on the track
  double theta432(0);
  double beta432(0);
  X0 = 100;  //An estimate of X0
  double thickness = 0.00005;  //Thickness of mimosa layer is 50 micrometers
  double X = X0/thickness;
  double BetterX = DBL_MAX;
  double currentlikelihood = DBL_MIN;
  double lastX0 = DBL_MIN;
  while(X0 != BetterX){
    sigma_measi = 5;  //HACK - This is actually the sum of the variances of theta in the forward and backward direction. 
    size_t size = _histoData["Theta012"].size();
    double loglikelihoodX0(0);
    for(size_t i = 0; i < size; ++i){
      theta012 = _histoData["Theta012"][i]*3.1415/180.0;
      theta432 = _histoData["Theta432"][i]*3.1415/180.0;
      beta012 = _histoData["Phi012"][i]*3.1415/180.0;
      beta432 = _histoData["Phi432"][i]*3.1415/180.0;
      double deltatheta = sqrt((theta432 - theta012)*(theta432 - theta012));  //Take the positive difference between the two angles in theta and beta
      double deltabeta = sqrt((beta432 - beta012)*(beta432 - beta012));
      sigma_msi = (0.015*E/(p*p))*sqrt(X/(sin(deltatheta)*cos(deltabeta)))*(1 + 0.038*std::log(X/(sin(deltatheta)*cos(deltabeta))));
      if(sigma_msi != sigma_msi){
        streamlog_out( ERROR5 ) << "Code failed with following values: " << std::endl << "E = " << E << std::endl << "p = " << p << std::endl << "X = " << X << std::endl << "theta012 = " << theta012 << std::endl << "beta012 = " << beta012 << std::endl << "theta432 = " << theta432 << std::endl << "beta432 = " << beta432 << std::endl << "deltatheta = " << deltatheta << std::endl << "deltabeta = " << deltabeta << std::endl << "sigma_msi = " << sigma_msi << std::endl << "sigma_measi = " << sigma_measi << std::endl << "loglikelihoodX0 = " << loglikelihoodX0 << std::endl << std::endl;
      }
      //Maximum log-likelihood:
      loglikelihoodX0 += -std::log(sigma_msi*sigma_msi + sigma_measi*sigma_measi) - (theta012*theta012)/(sigma_msi*sigma_msi + sigma_measi*sigma_measi); 
    }
    if(loglikelihoodX0 > currentlikelihood){
      BetterX = X0;
      currentlikelihood = loglikelihoodX0;
      if(lastX0 < X0){
        X++;
      }
      else{
        X--;
      }
    }
    else{
      if(lastX0 < X0){
        X--;
      }
      else{
        X++;
      }
    }
    X = X0/thickness;
    loglikelihoodX0 = 0;
  }
  return X0;
}

std::vector< TVector3* > EUTelX0Processor::getHitsFromTrack(Track *track){
 std::vector< TrackerHit* > trackhits = track->getTrackerHits();
 std::vector< TVector3* > hits;
 for(std::vector< TrackerHit* >::iterator it = trackhits.begin(); it != trackhits.end(); ++it){
    if((*it)->getType() == 32){  //Check if the hit type is appropriate
      Double_t x = (*it)->getPosition()[0];
      Double_t y = (*it)->getPosition()[1];
      Double_t z = (*it)->getPosition()[2];
      TVector3 *tempvec = new TVector3(x,y,z);
      hits.push_back(tempvec);
    } //End of if query for type check
  } //End of for loop running through trackhits
  return hits;
}

void EUTelX0Processor::singlePointResolution(Track *track){
//This function draws a line between two hits on adjacent planes and then projects in the forward and backward direction to the hit on the next plane.
//It then compares that hit with the actual hit on that next plane and the differencce is plotted
  streamlog_out(DEBUG1) << "Function EUTeoX0Processor::singlePointResolution(Track *" << &track << ") called" << std::endl;
  std::vector< TVector3* > hits = getHitsFromTrack(track);
  int i = 0;
  for(std::vector< TVector3* >::iterator it = hits.begin(); it != hits.end(); ++it){
    stringstream histogramnamex,histogramnamey;
    histogramnamex << "SinglePointResidualXPlane" << i;
    histogramnamey << "SinglePointResidualYPlane" << i;
    double x0, y0, x1, y1, x2, y2, z0, z1, z2;
    x0 = (*it)->x();
    y0 = (*it)->y();
    z0 = (*it)->z();
    x1 = (*it + 1)->x();
    y1 = (*it + 1)->y();
    z1 = (*it + 1)->z();
    x2 = (*it + 2)->x();
    y2 = (*it + 2)->y();
    z2 = (*it + 2)->z();
    double predictx, predicty;
    if(it == hits.begin() || it == hits.begin() + 1){
      predictx = (x2-x1)*(z2-z0)/(z2-z1);
      predicty = (y2-y1)*(z2-z0)/(z2-z1);
      try{
        dynamic_cast< AIDA::IHistogram1D* > (_histoThing[histogramnamex.str().c_str()])->fill(x0 - predictx);
      } catch(std::bad_cast &bc){
        streamlog_out(ERROR3) << "Unable to fill histogram: " << histogramnamex.str().c_str() << endl;
      }
      try{
        dynamic_cast< AIDA::IHistogram1D* > (_histoThing[histogramnamey.str().c_str()])->fill(y0 - predicty);
      } catch(std::bad_cast &bc){
        streamlog_out(ERROR3) << "Unable to fill histogram: " << histogramnamey.str().c_str() << endl;
      }
    } else {
      predictx = (x1-x0)*(z2-z0)/(z1-z0);
      predicty = (y1-y0)*(z2-z0)/(z1-z0);
       try{
        dynamic_cast< AIDA::IHistogram1D* > (_histoThing[histogramnamex.str().c_str()])->fill(x2 - predictx);
      } catch(std::bad_cast &bc){
        streamlog_out(ERROR3) << "Unable to fill histogram: " << histogramnamex.str().c_str() << endl;
      }
      try{
        dynamic_cast< AIDA::IHistogram1D* > (_histoThing[histogramnamey.str().c_str()])->fill(y2 - predicty);
      } catch(std::bad_cast &bc){
        streamlog_out(ERROR3) << "Unable to fill histogram: " << histogramnamey.str().c_str() << endl;
      }     
    }
  }

}

void EUTelX0Processor::threePointResolution(Track *track){
//This function draws a line between the hits in plane i and i+2
//then it compares where the hit in plane i+1 is with the average of the other two planes
//this is then plotted as a residual for each plane.
  streamlog_out(ERROR0) << "Function EUTelX0Processor::threePointResolution(Track *" << &track << ") called" << std::endl;

  std::vector< TrackerHit* > trackhits = track->getTrackerHits();
  std::vector< TVector3* > hits = getHitsFromTrack(track);
  int i = 1;
  for(std::vector< TVector3* >::iterator it = hits.begin(); it != hits.end() - 2; ++it){
    //This determines the guess of the position of the particle as it should hit the middle sensor
    streamlog_out(ERROR0) << "In the for loop of the three point resolution function" << endl;
    double averagingfactor = ((*it + 1)->z() - (*it)->z())/((*it + 2)->z() - (*it)->z());
    if(averagingfactor == 0){
      streamlog_out(ERROR0) << "Averaging factor == 0)" << endl;
      continue;
    }
    double averagex = ((*it)->x() + (*it + 2)->x())/averagingfactor;
    double averagey = ((*it)->y() + (*it + 2)->x())/averagingfactor;
    double middlex = (*it + 1)->x();
    double middley = (*it + 1)->y();
    std::stringstream ResidualX, ResidualY;
    ResidualX << "ThreePointResidualXPlane" << i;
    ResidualX << "ThreePointResidualYPlane" << i;
    try{
      dynamic_cast< AIDA::IHistogram1D* > (_histoThing[ResidualX.str().c_str()])->fill(averagex - middlex);
    } catch(std::bad_cast& bc){
      streamlog_out(ERROR5) << "Unable to fill histogram " << ResidualX.str() << " due to bad cast" << std::endl;
    }
    try{
      dynamic_cast< AIDA::IHistogram1D* > (_histoThing[ResidualY.str().c_str()])->fill(averagey - middley);
    } catch(std::bad_cast& bc){
      streamlog_out(ERROR5) << "Unable to fill histogram " << ResidualY.str() << " due to bad cast" << std::endl;
    }
    streamlog_out(ERROR0) << "Made it to the end of the for loop" << endl;
    i++;
  }
}

