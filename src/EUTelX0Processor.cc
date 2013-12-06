// Version: $Id$
#include "EUTelX0Processor.h"
#include <cmath>
#include "TCanvas.h"
#include "TStyle.h"
#include <cstdlib>
using namespace marlin;
using namespace eutelescope;
using namespace std;

EUTelX0Processor::EUTelX0Processor()
  :Processor("EUTelX0Processor"),
  _beamEnergy(0.0), 
  _dutPosition(0.0),
  _eventNumber(0),
  _histoThing(),
  _histoThing2D(),
  _inputTrackCollectionVec(NULL), 
  _inputTrackColName(""),
  _runNumber(0),
  _trackCollectionName(""),
  nobins(0),
  nobinsangle(0),
  minbin(0),
  maxbin(0),
  minbinangle(0),
  maxbinangle(0),
  binsx(0),
  minx(0.0),//(mm)
  maxx(0.0),
  binsy(0),
  miny(0.0),
  maxy(0.0),
  X0ProcessorDirectory(NULL),
  AngleXForwardTripleFirstThreePlanes(NULL),
  AngleXForwardTripleLastThreePlanes(NULL),
  AngleYForwardTripleFirstThreePlanes(NULL),
  AngleYForwardTripleLastThreePlanes(NULL),
  AngleXYForwardTripleFirstThreePlanes(NULL),
  AngleXYForwardTripleLastThreePlanes(NULL),
  ScatteringAngleXTriple(NULL),
  ScatteringAngleYTriple(NULL),
  ScatteringAngleXYTriple(NULL),
  ScatteringAngleXTripleMap(NULL),
  ScatteringAngleYTripleMap(NULL),
  RadiationLengthTripleMap(NULL),
  SinglePointResidualXPlane0(NULL),
  SinglePointResidualXPlane1(NULL),
  SinglePointResidualXPlane2(NULL),
  SinglePointResidualXPlane3(NULL),
  SinglePointResidualXPlane4(NULL),
  SinglePointResidualXPlane5(NULL),
  SinglePointResidualYPlane0(NULL),
  SinglePointResidualYPlane1(NULL),
  SinglePointResidualYPlane2(NULL),
  SinglePointResidualYPlane3(NULL),
  SinglePointResidualYPlane4(NULL),
  SinglePointResidualYPlane5(NULL),
  ThreePointResidualXPlane1(NULL),
  ThreePointResidualXPlane2(NULL),
  ThreePointResidualXPlane3(NULL),
  ThreePointResidualXPlane4(NULL),
  ThreePointResidualYPlane1(NULL),
  ThreePointResidualYPlane2(NULL),
  ThreePointResidualYPlane3(NULL),
  ThreePointResidualYPlane4(NULL),
  AngleXForwardPlane0(NULL),
  AngleXForwardPlane1(NULL),
  AngleXForwardPlane2(NULL),
  AngleXForwardPlane3(NULL),
  AngleXForwardPlane4(NULL),
  AngleYForwardPlane0(NULL),
  AngleYForwardPlane1(NULL),
  AngleYForwardPlane2(NULL),
  AngleYForwardPlane3(NULL),
  AngleYForwardPlane4(NULL),
  AngleXYForwardPlane0(NULL),
  AngleXYForwardPlane1(NULL),
  AngleXYForwardPlane2(NULL),
  AngleXYForwardPlane3(NULL),
  AngleXYForwardPlane4(NULL),
  ScatteringAngleXPlane1(NULL),
  ScatteringAngleXPlane2(NULL),
  ScatteringAngleXPlane3(NULL),
  ScatteringAngleXPlane4(NULL),
  ScatteringAngleYPlane1(NULL),
  ScatteringAngleYPlane2(NULL),
  ScatteringAngleYPlane3(NULL),
  ScatteringAngleYPlane4(NULL),
  KinkAnglePlane1(NULL),
  KinkAnglePlane2(NULL),
  KinkAnglePlane3(NULL),
  KinkAnglePlane4(NULL),
  ScatteringAngleXSingleMap(),
  ScatteringAngleYSingleMap(),
  RadiationLengthSingleMap(),
  ScatteringAngleXSingleMapData(),
  ScatteringAngleYSingleMapData(),
  totaleventnumber(0),
  totaltracknumber(0),
  totaltracknumberdoubledaf(0)
//  ScatteringAngleXPlane1Map(NULL),
//  ScatteringAngleXPlane2Map(NULL),
//  ScatteringAngleXPlane3Map(NULL),
//  ScatteringAngleXPlane4Map(NULL),
//  ScatteringAngleYPlane1Map(NULL),
//  ScatteringAngleYPlane2Map(NULL),
//  ScatteringAngleYPlane3Map(NULL),
//  ScatteringAngleYPlane4Map(NULL),
//  RadiationLengthPlane1Map(NULL),
//  RadiationLengthPlane2Map(NULL),
//  RadiationLengthPlane3Map(NULL),
//  RadiationLengthPlane4Map(NULL),
//  ScatteringAngleXTripleMapData(),
//  ScatteringAngleYTripleMapData(),
//  ScatteringAngleXPlane1MapData(),
//  ScatteringAngleXPlane2MapData(),
//  ScatteringAngleXPlane3MapData(),
//  ScatteringAngleXPlane4MapData(),
//  ScatteringAngleYPlane1MapData(),
//  ScatteringAngleYPlane2MapData(),
//  ScatteringAngleYPlane3MapData(),
//  ScatteringAngleYPlane4MapData()
{
  streamlog_out(DEBUG1) << "Constructing the EUTelX0Processor, setting all values to zero or NULL" << std::endl;
  _inputTrackCollectionVec = new LCCollectionVec(LCIO::TRACK);//Used to store the values of the hit events
  
  registerInputCollection(LCIO::TRACK,"InputTrackCollectionName",
                           "Collection name for fitted tracks",
                           _trackCollectionName, string ("AlignedTrack"));

  registerInputCollection(LCIO::TRACK,"InputTrackCollectionName1",
                           "Collection name for fitted tracks",
                           _trackCollectionName1, string ("AlignedTrack1"));

  registerInputCollection(LCIO::TRACK,"InputTrackCollectionName2",
                           "Collection name for fitted tracks",
                           _trackCollectionName2, string ("AlignedTrack2"));

  registerProcessorParameter("BeamEnergy","Works out the energy of the beam for radiation lengths", _beamEnergy, static_cast< double > (0.0));
  registerProcessorParameter("DoubleDafFitted","Are two tracks from seperate daf fitters being used (front and back three planes only) to make a scattering angle? If yes then this should be true", _doubleDafFitted, static_cast< bool > (false));
  registerProcessorParameter("DUTPosition","The exact position of the DUT in the z direction relative to Plane 2 of the telescope", _dutPosition, static_cast< double > (0.0));
  registerProcessorParameter("RadiationLengthMapMinX","Used to determine the minimum X for the radiation length map, measured in XXX", minx, static_cast< double > (-11.0));
  registerProcessorParameter("RadiationLengthMapMaxX","Used to determine the maximum X for the radiation length map, measured in XXX", maxx, static_cast< double > (11.0));
  registerProcessorParameter("RadiationLengthMapMinY","Used to determine the minimum Y for the radiation length map, measured in XXX", miny, static_cast< double > (-6.0));
  registerProcessorParameter("RadiationLengthMapMaxY","Used to determine the maximum Y for the radiation length map, measured in XXX", maxy, static_cast< double > (6.0));
  registerProcessorParameter("RadiationLengthMapBinSizeX","Used to determine the spatial resolution in X for the radiation length map, measured in XXX", binsizex, static_cast< double > (1.0));
  registerProcessorParameter("RadiationLengthMapBinSizeY","Used to determine the spatial resolution in Y for the radiation length map, measured in XXX", binsizey, static_cast< double > (1.0));
  registerProcessorParameter("Cut","This is the maximum allowed distance in X and y between two tracks and where they meet on the DUT", _cut, static_cast< double > (0.0));
}

void EUTelX0Processor::init()
{
  gStyle->SetOptStat("neMRuo");
  streamlog_out(DEBUG5) << "Running EUTelX0Processor::init()" << std::endl;
  nobins = 10000;
  nobinsangle = 1000;//Number of bins in the histograms
  minbin = -0.1;
  maxbin = 0.1;//Maximum and minimum bin values
  minbinangle = -0.02;
  maxbinangle = 0.02;
  std::vector<double> empty;  
  binsx = static_cast< int >((maxx-minx)/binsizex);
  binsy = static_cast< int >((maxy-miny)/binsizey);

  AngleXFrontThreePlanesDoubleDaf = new TH1D("AngleXFrontThreePlanesDoubleDaf", "Angles of Tracks in X Direction relative to the Z Axis for  Three Planes from Double Daf Fitting; \\theta_{x} (rads); Count", nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleXFrontThreePlanesDoubleDaf"] = AngleXFrontThreePlanesDoubleDaf;

  AngleYFrontThreePlanesDoubleDaf = new TH1D("AngleYFrontThreePlanesDoubleDaf", "Angles of Tracks in Y Direction relative to the Z Axis for  Three Planes from Double Daf Fitting; \\theta_{x} (rads); Count", nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleYFrontThreePlanesDoubleDaf"] = AngleYFrontThreePlanesDoubleDaf;

  AngleXBackThreePlanesDoubleDaf = new TH1D("AngleXBackThreePlanesDoubleDaf", "Angles of Tracks in X Direction relative to the Z Axis for  Three Planes from Double Daf Fitting; \\theta_{x} (rads); Count", nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleXBackThreePlanesDoubleDaf"] = AngleXBackThreePlanesDoubleDaf;

  AngleYBackThreePlanesDoubleDaf = new TH1D("AngleYBackThreePlanesDoubleDaf", "Angles of Tracks in Y Direction relative to the Z Axis for  Three Planes from Double Daf Fitting; \\theta_{x} (rads); Count", nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleYBackThreePlanesDoubleDaf"] = AngleYBackThreePlanesDoubleDaf;

  AngleXYFrontThreePlanesDoubleDaf = new TH2D("AngleXYFrontThreePlanesDoubleDaf", "Angles of Tracks in Y Direction relative to the Z Axis for  Three Planes from Double Daf Fitting; \\theta_{x} (rads); \\theta_{y} (rads)", nobinsangle,minbinangle,maxbinangle,nobinsangle,minbinangle,maxbinangle);
  _histoThing2D["AngleXYFrontThreePlanesDoubleDaf"] = AngleXYFrontThreePlanesDoubleDaf;

  AngleXYBackThreePlanesDoubleDaf = new TH2D("AngleXYBackThreePlanesDoubleDaf", "Angles of Tracks in Y Direction relative to the Z Axis for  Three Planes from Double Daf Fitting; \\theta_{x} (rads); \\theta_{y} (rads)", nobinsangle,minbinangle,maxbinangle,nobinsangle,minbinangle,maxbinangle);
  _histoThing2D["AngleXYBackThreePlanesDoubleDaf"] = AngleXYBackThreePlanesDoubleDaf;

  ScatteringAngleXDoubleDaf = new TH1D("ScatteringAngleXDoubleDaf", "Scattering Angles of Tracks in X Direction relative to the Z Axis from Double Daf Fitting; \\theta_{x} (rads); Count", nobinsangle,minbinangle,maxbinangle);
  _histoThing["ScatteringAngleXDoubleDaf"] = ScatteringAngleXDoubleDaf;

  ScatteringAngleYDoubleDaf = new TH1D("ScatteringAngleYDoubleDaf", "Scattering Angles of Tracks in Y Direction relative to the Z Axis from Double Daf Fitting; \\theta_{x} (rads); Count", nobinsangle,minbinangle,maxbinangle);
  _histoThing["ScatteringAngleYDoubleDaf"] = ScatteringAngleYDoubleDaf;

  ScatteringAngleXYDoubleDaf = new TH2D("ScatteringAngleXYDoubleDaf", "Scattering Angles of Tracks in X and Y Direction relative to the Z Axis from Double Daf Fitting; \\theta_{x} (rads); \\theta_{y} (rads)", nobinsangle,minbinangle,maxbinangle,nobinsangle,minbinangle,maxbinangle);
  _histoThing2D["ScatteringAngleXYDoubleDaf"] = ScatteringAngleXYDoubleDaf;

  ScatteringAngleXTripleMapDoubleDaf = new TH2D("ScatteringAngleXTripleMapDoubleDaf",
                                       "Scattering Angle Map in X on DUT;x (mm);y (mm); \\theta_x (rads)"
				       ,binsx,minx,maxx,binsy,miny,maxy);
  _histoThing2D["ScatteringAngleXTripleMapDoubleDaf"] = ScatteringAngleXTripleMapDoubleDaf;
 
  ScatteringAngleYTripleMapDoubleDaf = new TH2D("ScatteringAngleYTripleMapDoubleDaf",
                                       "Scattering Angle Map in Y on DUT;x (mm);y (mm); \\theta_x (rads)"
				       ,binsx,minx,maxx,binsy,miny,maxy);
  _histoThing2D["ScatteringAngleYTripleMapDoubleDaf"] = ScatteringAngleYTripleMapDoubleDaf;
 
  RadiationLengthMapDoubleDaf = new TH2D("RadiationLengthMapDoubleDaf",
                                       "Radiation Length Map on DUT;x (mm);y (mm); \\X_{0} (Percent)"
				       ,binsx,minx,maxx,binsy,miny,maxy);
  _histoThing2D["RadiationLengthMapDoubleDaf"] = RadiationLengthMapDoubleDaf;
 

  AngleXForwardTripleFirstThreePlanes = new TH1D("AngleXForwardTripleFirstThreePlanes",
                                 "Angle of Tracks in X Direction Relative to the Z Axis for First Three Planes;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleXForwardTripleFirstThreePlanes"] = AngleXForwardTripleFirstThreePlanes;
 
  AngleXForwardTripleLastThreePlanes = new TH1D("AngleXForwardTripleLastThreePlanes",
                                 "Angle of Tracks in X Direction Relative to the Z Axis for Last Three Planes;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleXForwardTripleLastThreePlanes"] = AngleXForwardTripleLastThreePlanes;
 
  AngleYForwardTripleFirstThreePlanes = new TH1D("AngleYForwardTripleFirstThreePlanes",
                                 "Angle of Tracks in Y Direction Relative to the Z Axis for First Three Planes;\\theta_y (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleYForwardTripleFirstThreePlanes"] = AngleYForwardTripleFirstThreePlanes;
 
  AngleYForwardTripleLastThreePlanes = new TH1D("AngleYForwardTripleLastThreePlanes",
                                 "Angle of Tracks in Y Direction Relative to the Z Axis for Last Three Planes;\\theta_y (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleYForwardTripleLastThreePlanes"] = AngleYForwardTripleLastThreePlanes;
 
  AngleXYForwardTripleFirstThreePlanes = new TH2D("AngleXYForwardTripleFirstThreePlanes",
                                 "Angle of Tracks in XY Direction Relative to the Z Axis for First Three Planes;\\theta_x (rads);\\theta_y (rads);Count",
				 nobinsangle,minbinangle,maxbinangle,nobinsangle,minbinangle,maxbinangle);
  _histoThing2D["AngleXYForwardTripleFirstThreePlanes"] = AngleXYForwardTripleFirstThreePlanes;
 
  AngleXYForwardTripleLastThreePlanes = new TH2D("AngleXYForwardTripleLastThreePlanes",
                                 "Angle of Tracks in XY Direction Relative to the Z Axis for Last Three Planes;\\theta_y (rads);Count",
				 nobinsangle,minbinangle,maxbinangle,nobinsangle,minbinangle,maxbinangle);
  _histoThing2D["AngleXYForwardTripleLastThreePlanes"] = AngleXYForwardTripleLastThreePlanes;
 
  ScatteringAngleXTriple = new TH1D("ScatteringAngleXTriple",
                                    "Scattering Angle in X Direction Relative to the Z Axis for Triplet Tracks;\\theta_x (rads);Count",
                                    nobinsangle,minbinangle,maxbinangle);
  _histoThing["ScatteringAngleXTriple"] = ScatteringAngleXTriple;
 
  ScatteringAngleYTriple = new TH1D("ScatteringAngleYTriple",
                                    "Scattering Angle in Y Direction Relative to the Z Axis for Triplet Tracks;\\theta_y (rads);Count",
                                    nobinsangle,minbinangle,maxbinangle);
  _histoThing["ScatteringAngleYTriple"] = ScatteringAngleYTriple;
  
  ScatteringAngleXYTriple = new TH2D("ScatteringAngleXYTriple",
                                    "Scattering Angle in XY Direction Relative to the Z Axis for Triplet Tracks;\\theta_x (rads);\\theta_y (rads);Count",
                                    nobinsangle,minbinangle,maxbinangle,nobinsangle,minbinangle,maxbinangle);
  _histoThing2D["ScatteringAngleXYTriple"] = ScatteringAngleXYTriple;
  
  ScatteringAngleXTripleMap = new TH2D("ScatteringAngleXTripleMap",
                                       "Scattering Angle Map in X on DUT;x (mm);y (mm); \\theta_x (rads)"
				       ,binsx,minx,maxx,binsy,miny,maxy);
  _histoThing2D["ScatteringAngleXTripleMap"] = ScatteringAngleXTripleMap;
 
  ScatteringAngleYTripleMap = new TH2D("ScatteringAngleYTripleMap",
                                       "Scattering Angle Map in Y on DUT;x (mm);y (mm); \\theta_y (rads)"
				       ,binsx,minx,maxx,binsy,miny,maxy);
  _histoThing2D["ScatteringAngleYTripleMap"] = ScatteringAngleYTripleMap;
 
  RadiationLengthTripleMap = new TH2D("RadiationLengthTripleMap",
                                       "Radiation Length Map on DUT;x (mm);y (mm); \\X_0 (Percent)"
				       ,binsx,minx,maxx,binsy,miny,maxy);
  _histoThing2D["RadiationLengthTripleMap"] = RadiationLengthTripleMap;
 
  SinglePointResidualXPlane1 = new TH1D("SinglePointResidualXPlane1",
                                         "Single Point Residual in X on Plane 1;\\Delta X (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["SinglePointResidualXPlane1"] = SinglePointResidualXPlane1;
 
  SinglePointResidualXPlane2 = new TH1D("SinglePointResidualXPlane2",
                                         "Single Point Residual in X on Plane 2;\\Delta X (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["SinglePointResidualXPlane2"] = SinglePointResidualXPlane2;
  
  SinglePointResidualXPlane3 = new TH1D("SinglePointResidualXPlane3",
                                         "Single Point Residual in X on Plane 3;\\Delta X (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["SinglePointResidualXPlane3"] = SinglePointResidualXPlane3;
  
  SinglePointResidualXPlane4 = new TH1D("SinglePointResidualXPlane4",
                                         "Single Point Residual in X on Plane 4;\\Delta X (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["SinglePointResidualXPlane4"] = SinglePointResidualXPlane4;
  
  SinglePointResidualXPlane5 = new TH1D("SinglePointResidualXPlane5",
                                         "Single Point Residual in X on Plane 5;\\Delta X (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["SinglePointResidualXPlane5"] = SinglePointResidualXPlane5;
 
  SinglePointResidualYPlane0 = new TH1D("SinglePointResidualYPlane0",
                                         "Single Point Residual in Y on Plane 0;\\Delta Y (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["SinglePointResidualYPlane0"] = SinglePointResidualYPlane0;
  
  SinglePointResidualYPlane1 = new TH1D("SinglePointResidualYPlane1",
                                         "Single Point Residual in Y on Plane 1;\\Delta Y (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["SinglePointResidualYPlane1"] = SinglePointResidualYPlane1;
  
  SinglePointResidualYPlane2 = new TH1D("SinglePointResidualYPlane2",
                                         "Single Point Residual in Y on Plane 2;\\Delta Y (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["SinglePointResidualYPlane2"] = SinglePointResidualYPlane2;
  
  SinglePointResidualYPlane3 = new TH1D("SinglePointResidualYPlane3",
                                         "Single Point Residual in Y on Plane 3;\\Delta Y (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["SinglePointResidualYPlane3"] = SinglePointResidualYPlane3;
  
  SinglePointResidualYPlane4 = new TH1D("SinglePointResidualYPlane4",
                                         "Single Point Residual in Y on Plane 4;\\Delta Y (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["SinglePointResidualYPlane4"] = SinglePointResidualYPlane4;
  
  SinglePointResidualYPlane5 = new TH1D("SinglePointResidualYPlane5",
                                         "Single Point Residual in Y on Plane 5;\\Delta Y (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["SinglePointResidualYPlane5"] = SinglePointResidualYPlane5;
 
  ThreePointResidualXPlane1 = new TH1D("ThreePointResidualXPlane1",
                                         "Three Point Residual in X on Plane 1;\\Delta X (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["ThreePointResidualXPlane1"] = ThreePointResidualXPlane1;
  
  ThreePointResidualXPlane2 = new TH1D("ThreePointResidualXPlane2",
                                         "Three Point Residual in X on Plane 2;\\Delta X (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["ThreePointResidualXPlane2"] = ThreePointResidualXPlane2;
  
  ThreePointResidualXPlane3 = new TH1D("ThreePointResidualXPlane3",
                                         "Three Point Residual in X on Plane 3;\\Delta X (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["ThreePointResidualXPlane3"] = ThreePointResidualXPlane3;
  
  ThreePointResidualXPlane4 = new TH1D("ThreePointResidualXPlane4",
                                         "Three Point Residual in X on Plane 4;\\Delta X (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["ThreePointResidualXPlane4"] = ThreePointResidualXPlane4;
  
  ThreePointResidualYPlane1 = new TH1D("ThreePointResidualYPlane1",
                                         "Three Point Residual in Y on Plane 1;\\Delta Y (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["ThreePointResidualYPlane1"] = ThreePointResidualYPlane1;
  
  ThreePointResidualYPlane2 = new TH1D("ThreePointResidualYPlane2",
                                         "Three Point Residual in Y on Plane 2;\\Delta Y (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["ThreePointResidualYPlane2"] = ThreePointResidualYPlane2;
  
  ThreePointResidualYPlane3 = new TH1D("ThreePointResidualYPlane3",
                                         "Three Point Residual in Y on Plane 3;\\Delta Y (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["ThreePointResidualYPlane3"] = ThreePointResidualYPlane3;
  
  ThreePointResidualYPlane4 = new TH1D("ThreePointResidualYPlane4",
                                         "Three Point Residual in Y on Plane 4;\\Delta Y (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["ThreePointResidualYPlane4"] = ThreePointResidualYPlane4;
  
  AngleXForwardPlane0 = new TH1D("AngleXForwardPlane0",
                                 "Angle of Tracks in X Direction Relative to the Z Axis on Plane 0;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleXForwardPlane0"] = AngleXForwardPlane0;
 
  AngleXForwardPlane1 = new TH1D("AngleXForwardPlane1",
                                 "Angle of Tracks in X Direction Relative to the Z Axis on Plane 1;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleXForwardPlane1"] = AngleXForwardPlane1;
 
  AngleXForwardPlane2 = new TH1D("AngleXForwardPlane2",
                                 "Angle of Tracks in X Direction Relative to the Z Axis on Plane 2;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleXForwardPlane2"] = AngleXForwardPlane2;
 
  AngleXForwardPlane3 = new TH1D("AngleXForwardPlane3",
                                 "Angle of Tracks in X Direction Relative to the Z Axis on Plane 3;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleXForwardPlane3"] = AngleXForwardPlane3;
 
  AngleXForwardPlane4 = new TH1D("AngleXForwardPlane4",
                                 "Angle of Tracks in X Direction Relative to the Z Axis on Plane 4;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleXForwardPlane4"] = AngleXForwardPlane4;
  
  AngleYForwardPlane0 = new TH1D("AngleYForwardPlane0",
                                 "Angle of Tracks in Y Direction Relative to the Z Axis on Plane 0;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleYForwardPlane0"] = AngleYForwardPlane0;
 
  AngleYForwardPlane1 = new TH1D("AngleYForwardPlane1",
                                 "Angle of Tracks in Y Direction Relative to the Z Axis on Plane 1;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleYForwardPlane1"] = AngleYForwardPlane1;
 
  AngleYForwardPlane2 = new TH1D("AngleYForwardPlane2",
                                 "Angle of Tracks in Y Direction Relative to the Z Axis on Plane 2;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleYForwardPlane2"] = AngleYForwardPlane2;
 
  AngleYForwardPlane3 = new TH1D("AngleYForwardPlane3",
                                 "Angle of Tracks in Y Direction Relative to the Z Axis on Plane 3;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleYForwardPlane3"] = AngleYForwardPlane3;
 
  AngleYForwardPlane4 = new TH1D("AngleYForwardPlane4",
                                 "Angle of Tracks in Y Direction Relative to the Z Axis on Plane 4;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleYForwardPlane4"] = AngleYForwardPlane4;
  
  AngleXYForwardPlane0 = new TH2D("AngleXYForwardPlane0",
                                 "Angle of Tracks in XY Direction Relative to the Z Axis on Plane 0;\\theta_x (rads);\\theta_y (rads)",
				 nobinsangle,minbinangle,maxbinangle,nobinsangle,minbinangle,maxbinangle);
  _histoThing2D["AngleXYForwardPlane0"] = AngleXYForwardPlane0;
  
  AngleXYForwardPlane1 = new TH2D("AngleXYForwardPlane1",
                                 "Angle of Tracks in XY Direction Relative to the Z Axis on Plane 1;\\theta_x (rads);\\theta_y (rads)",
				 nobinsangle,minbinangle,maxbinangle,nobinsangle,minbinangle,maxbinangle);
  _histoThing2D["AngleXYForwardPlane1"] = AngleXYForwardPlane1;
  
  AngleXYForwardPlane2 = new TH2D("AngleXYForwardPlane2",
                                 "Angle of Tracks in XY Direction Relative to the Z Axis on Plane 2;\\theta_x (rads);\\theta_y (rads)",
				 nobinsangle,minbinangle,maxbinangle,nobinsangle,minbinangle,maxbinangle);
  _histoThing2D["AngleXYForwardPlane2"] = AngleXYForwardPlane2;
  
  AngleXYForwardPlane3 = new TH2D("AngleXYForwardPlane3",
                                 "Angle of Tracks in XY Direction Relative to the Z Axis on Plane 3;\\theta_x (rads);\\theta_y (rads)",
				 nobinsangle,minbinangle,maxbinangle,nobinsangle,minbinangle,maxbinangle);
  _histoThing2D["AngleXYForwardPlane3"] = AngleXYForwardPlane3;
  
  AngleXYForwardPlane4 = new TH2D("AngleXYForwardPlane4",
                                 "Angle of Tracks in XY Direction Relative to the Z Axis on Plane 4;\\theta_x (rads);\\theta_y (rads)",
				 nobinsangle,minbinangle,maxbinangle,nobinsangle,minbinangle,maxbinangle);
  _histoThing2D["AngleXYForwardPlane4"] = AngleXYForwardPlane4;
  
  ScatteringAngleXPlane1 = new TH1D("ScatteringAngleXPlane1",
                                    "Scattering Angle in X Direction Relative to the Z Axis on Plane 1;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["ScatteringAngleXPlane1"] = ScatteringAngleXPlane1;
 
  ScatteringAngleXPlane2 = new TH1D("ScatteringAngleXPlane2",
                                    "Scattering Angle in X Direction Relative to the Z Axis on Plane 2;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["ScatteringAngleXPlane2"] = ScatteringAngleXPlane2;
 
  ScatteringAngleXPlane3 = new TH1D("ScatteringAngleXPlane3",
                                    "Scattering Angle in X Direction Relative to the Z Axis on Plane 3;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["ScatteringAngleXPlane3"] = ScatteringAngleXPlane3;
 
  ScatteringAngleXPlane4 = new TH1D("ScatteringAngleXPlane4",
                                    "Scattering Angle in X Direction Relative to the Z Axis on Plane 4;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["ScatteringAngleXPlane4"] = ScatteringAngleXPlane4;
  
  ScatteringAngleYPlane1 = new TH1D("ScatteringAngleYPlane1",
                                    "Scattering Angle in Y Direction Relative to the Z Axis on Plane 1;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["ScatteringAngleYPlane1"] = ScatteringAngleYPlane1;
 
  ScatteringAngleYPlane2 = new TH1D("ScatteringAngleYPlane2",
                                    "Scattering Angle in Y Direction Relative to the Z Axis on Plane 2;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["ScatteringAngleYPlane2"] = ScatteringAngleYPlane2;
 
  ScatteringAngleYPlane3 = new TH1D("ScatteringAngleYPlane3",
                                    "Scattering Angle in Y Direction Relative to the Z Axis on Plane 3;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["ScatteringAngleYPlane3"] = ScatteringAngleYPlane3;
 
  ScatteringAngleYPlane4 = new TH1D("ScatteringAngleYPlane4",
                                    "Scattering Angle in Y Direction Relative to the Z Axis on Plane 4;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["ScatteringAngleYPlane4"] = ScatteringAngleYPlane4;
     
  KinkAnglePlane1 = new TH2D("KinkAnglePlane1",
                             "Kink Angle on Plane 1;\\theta_x (rads);\\theta_y (rads)",
			     nobinsangle,minbinangle,maxbinangle,nobinsangle,minbinangle,maxbinangle);
  _histoThing2D["KinkAnglePlane1"] = KinkAnglePlane1;
  
  KinkAnglePlane2 = new TH2D("KinkAnglePlane2",
                             "Kink Angle on Plane 2;\\theta_x (rads);\\theta_y (rads)",
			     nobinsangle,minbinangle,maxbinangle,nobinsangle,minbinangle,maxbinangle);
  _histoThing2D["KinkAnglePlane2"] = KinkAnglePlane2;
  
  KinkAnglePlane3 = new TH2D("KinkAnglePlane3",
                             "Kink Angle on Plane 3;\\theta_x (rads);\\theta_y (rads)",
			     nobinsangle,minbinangle,maxbinangle,nobinsangle,minbinangle,maxbinangle);
  _histoThing2D["KinkAnglePlane3"] = KinkAnglePlane3;
  
  KinkAnglePlane4 = new TH2D("KinkAnglePlane4",
                             "Kink Angle on Plane 4;\\theta_x (rads);\\theta_y (rads)",
			     nobinsangle,minbinangle,maxbinangle,nobinsangle,minbinangle,maxbinangle);
  _histoThing2D["KinkAnglePlane4"] = KinkAnglePlane4;

  for(int i = 1; i < 5; ++i){ //This covers the inner 4 planes
    stringstream inttostring;
    inttostring << i;
    string titlex = "ScatteringAngleXPlane";
    string titlexdetail = "Scattering Angle in X on Plane";
    string titley = "ScatteringAngleYPlane";
    string titleydetail = "Scattering Angle in Y on Plane";
    string titlerad = "RadiationLengthPlane";
    string titleraddetail = "Radiation Length on Plane";
    string map = "Map";
    string axis = ";x (mm);y(mm)";
    string plane = inttostring.str();
    ScatteringAngleXSingleMap[i] = new TH2D((titlex + plane + map).c_str(),
                                       (titlexdetail + plane + axis).c_str()
                                       ,binsx,minx,maxx,binsy,miny,maxy);
    _histoThing2D[(titlex + plane + map).c_str()] = ScatteringAngleXSingleMap[i];
    ScatteringAngleYSingleMap[i] = new TH2D((titley + plane + map).c_str(),
                                       (titleydetail + plane + axis).c_str()
                                       ,binsy,miny,maxy,binsy,miny,maxy);
    _histoThing2D[(titley + plane + map).c_str()] = ScatteringAngleYSingleMap[i];
    RadiationLengthSingleMap[i] = new TH2D((titlerad + plane + map).c_str(),
                                       (titleraddetail + plane + axis).c_str()
                                       ,binsy,miny,maxy,binsy,miny,maxy);
    _histoThing2D[(titlerad + plane + map).c_str()] = RadiationLengthSingleMap[i];
  }
  gStyle->SetOptStat("neMRuo");
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
    streamlog_out(DEBUG0) << covmatrix.at(0) << ", " << covmatrix.at(1) << ", " << covmatrix.at(3) << ", " << covmatrix.at(6) << ", " << covmatrix.at(10) << std::endl
                          << covmatrix.at(1) << ", " << covmatrix.at(2) << ", " << covmatrix.at(4) << ", " << covmatrix.at(7) << ", " << covmatrix.at(11) << std::endl
                          << covmatrix.at(3) << ", " << covmatrix.at(4) << ", " << covmatrix.at(5) << ", " << covmatrix.at(8) << ", " << covmatrix.at(12) << std::endl
                          << covmatrix.at(6) << ", " << covmatrix.at(7) << ", " << covmatrix.at(8) << ", " << covmatrix.at(9) << ", " << covmatrix.at(13) << std::endl
                          << covmatrix.at(10) << ", " << covmatrix.at(11) << ", " << covmatrix.at(12) << ", " << covmatrix.at(13) << ", " << covmatrix.at(14) << std::endl;    
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

pair< vector< double >, vector< double > > EUTelX0Processor::GetSingleTrackAngles(vector< TVector3 > hits){
  gStyle->SetOptStat("neMRuo");
  streamlog_out(DEBUG1) << "Begin pair< vector< double >, vector< double > > EUTelX0Processor::GetSingleTrackAngles(vector< TVector3 > hits)" << endl;
  const size_t hitsize = hits.size();
  streamlog_out(DEBUG5) << "This track has " << hitsize << " hits" << endl;
  //Then we find the position where those lines would have hit the DUT
  //THIS IS TO BE DECIDED IF IT IS NEEDED LATER
  const double dutposition(hits[2].z());
  streamlog_out(DEBUG3) << "For now DUT position is set to the position of the 2nd plane, that is at " << dutposition << endl;
  //Then we work out the angles of these lines with respect to XZ and YZ, plot results in histograms
  if(hitsize <= 1 || hitsize > 6)
  {
    streamlog_out(DEBUG5) << "This track has one hit or less or more than 6 hits, aborting track" << endl;
    throw "Not enough hits";
  } 

  std::vector< double > scatterx, scattery;
  streamlog_out(DEBUG3) << "hitsize = " << hitsize << std::endl;
  for(size_t i = 0; i < hitsize-1; ++i){ //Fill histograms with the angles of the track at each plane and push back those angles into a vector
    const double x0 = hits.at(i).x();
    const double y0 = hits.at(i).y();
    const double z0 = hits.at(i).z();
    const double x1 = hits.at(i+1).x();
    const double y1 = hits.at(i+1).y();
    const double z1 = hits.at(i+1).z();
    const double deltaz = z1-z0;
    streamlog_out(DEBUG3) << x0 << endl << y0 << endl << z0 << endl << x1 << endl << y1 << endl << z1 << endl;
    const double frontanglex = atan2(x1-x0,deltaz);
    const double frontangley = atan2(y1-y0,deltaz);
    streamlog_out(DEBUG0) << "x0 = " << x0 << endl;
    streamlog_out(DEBUG0) << "y0 = " << y0 << endl;
    streamlog_out(DEBUG0) << "z0 = " << z0 << endl;
    streamlog_out(DEBUG0) << "x1 = " << x1 << endl;
    streamlog_out(DEBUG0) << "y1 = " << y1 << endl;
    streamlog_out(DEBUG0) << "z1 = " << z1 << endl;
    streamlog_out(DEBUG0) << "anglex = " << frontanglex << endl;
    streamlog_out(DEBUG0) << "angley = " << frontangley << endl;
    //Push back angles into a vector
    scatterx.push_back(frontanglex);
    scattery.push_back(frontangley);
    streamlog_out(DEBUG0) << "Line " << __LINE__ << endl;
    //Name the histograms
    std::stringstream xforward,yforward,xyforward;
    xforward << "AngleXForwardPlane" << i;
    yforward << "AngleYForwardPlane" << i;
    xyforward << "AngleXYForwardPlane" << i;
    //Fill the histograms with the XZ and YZ angles
    streamlog_out(DEBUG0) << "Line " << __LINE__ << endl;
    try{
    streamlog_out(DEBUG0) << "xforward = " << xforward.str().c_str() << endl;
    streamlog_out(DEBUG0) << "Line " << __LINE__ << endl;
      dynamic_cast< TH1D* >(_histoThing[xforward.str().c_str()])->Fill(frontanglex);
    streamlog_out(DEBUG0) << "Line " << __LINE__ << endl;
    } catch(std::bad_cast &bc){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << xforward.str().c_str() << endl;
    } catch(...){
      streamlog_out(ERROR9) << "Unknown Error when trying to fill histogram" << endl;
      streamlog_out(ERROR9) << "Unknown Error when trying to fill histogram: " << xforward.str().c_str() << " with the value " << frontanglex << endl;
    }
    try{
    streamlog_out(DEBUG0) << "Line " << __LINE__ << endl;
      dynamic_cast< TH1D* >(_histoThing[yforward.str().c_str()])->Fill(frontangley);
    streamlog_out(DEBUG0) << "Line " << __LINE__ << endl;
    } catch(std::bad_cast &bc){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << yforward.str().c_str() << endl;
    } catch(...){
      streamlog_out(ERROR9) << "Unknown Error when trying to fill histogram: " << yforward.str().c_str() << " with the value " << frontangley << endl;
    }
    try{
    streamlog_out(DEBUG0) << "Line " << __LINE__ << endl;
      dynamic_cast< TH2D* >(_histoThing2D[xyforward.str().c_str()])->Fill(frontanglex,frontangley);
    streamlog_out(DEBUG0) << "Line " << __LINE__ << endl;
    } catch(std::bad_cast &bc){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << xyforward.str().c_str() << endl;
    } catch(...){
      streamlog_out(ERROR9) << "Unknown Error when trying to fill histogram: " << xyforward.str().c_str() << " with the values " << frontanglex << ", " << frontangley << endl;
    }
  }
    streamlog_out(DEBUG0) << "Line " << __LINE__ << endl;
  pair< vector< double >, vector< double > > bothangles(scatterx,scattery);
    streamlog_out(DEBUG0) << "Line " << __LINE__ << endl;
  return bothangles;
}

void EUTelX0Processor::PlotTripleTrackAngleDoubleDafFitted(Track *track, bool front){
  streamlog_out(DEBUG1) << "Begin pair< vector< double >, vector< double > > EUTelX0Processor::GetTripleTrackAngles(vector< TVector3 > hits)" << endl;
  //Then we find the position where those lines would have hit the DUT
  //THIS IS TO BE DECIDED IF IT IS NEEDED LATER
  //Then we work out the angles of these lines with respect to XZ and YZ, plot results in histograms
  std::vector< TrackerHit* > trackhitsfront = track->getTrackerHits();
  std::vector< TVector3 > hitsfront;
  std::vector< TrackerHit* > trackhitsback = track->getTrackerHits();
  std::vector< TVector3 > hitsback;
  if(front == true){
    TVector3 *tempvec = new TVector3();
    for(std::vector< TrackerHit* >::iterator it = trackhitsfront.begin(); it != trackhitsfront.end(); ++it){
      if((*it)->getType() == 32){  //Check if the hit type is the aligned hits or the fitted points (we want fitted points)
        Double_t x = (*it)->getPosition()[0];
        Double_t y = (*it)->getPosition()[1];
        Double_t z = (*it)->getPosition()[2];
        tempvec->SetXYZ(x,y,z);
        hitsfront.push_back(*tempvec);
      } //End of if query for type check
    } //End of for loop running through trackhits
    delete tempvec;
  } else{
    TVector3 *tempvec = new TVector3();
    for(std::vector< TrackerHit* >::iterator it = trackhitsback.begin(); it != trackhitsback.end(); ++it){
      if((*it)->getType() == 32){  //Check if the hit type is the aligned hits or the fitted points (we want fitted points)
        Double_t x = (*it)->getPosition()[0];
        Double_t y = (*it)->getPosition()[1];
        Double_t z = (*it)->getPosition()[2];
        tempvec->SetXYZ(x,y,z);
        hitsback.push_back(*tempvec);
      } //End of if query for type check
    } //End of for loop running through trackhits
    delete tempvec;
  }
  std::vector< double > scatterx, scattery;
  double x1;
  double y1;
  double z1;
  double x2;
  double y2;
  double z2;
  double x3;
  double y3;
  double z3;
  double x4;
  double y4;
  double z4;
  double frontanglex;
  double frontangley;
  double backanglex;
  double backangley;
  //Name the histograms
  std::stringstream xforward,yforward,xyforward;
  streamlog_out(DEBUG2) << "Successfully pushed back the scatterx and y vectors and deleted the TGraphs and TFits:" << endl;
  stringstream xfront,yfront,xback,yback,xyfront,xyback;
  if(front == true){
    x1 = hitsfront.at(1).x();
    y1 = hitsfront.at(1).y();
    z1 = hitsfront.at(1).z();
    x2 = hitsfront.at(2).x();
    y2 = hitsfront.at(2).y();
    z2 = hitsfront.at(2).z();
    frontanglex = atan2((x2-x1),(z2-z1));
    frontangley = atan2((y2-y1),(z2-z1));
    xfront << "AngleXFrontThreePlanesDoubleDaf";
    yfront << "AngleYFrontThreePlanesDoubleDaf";
    streamlog_out(DEBUG2) << "Successfully named the histograms:" << endl;
    try{
      dynamic_cast< TH1D* >(_histoThing[xfront.str().c_str()])->Fill(frontanglex);
    } catch(std::bad_cast &bc){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << xfront.str().c_str() << endl;
    } catch(...){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << xfront.str().c_str() << ". Unknown Error" << endl;
    }
    try{
      dynamic_cast< TH1D* >(_histoThing[yfront.str().c_str()])->Fill(frontangley);
    } catch(std::bad_cast &bc){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << yfront.str().c_str() << endl;
    } catch(...){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << yfront.str().c_str() << ". Unknown Error" << endl;
    }
  } else{
    x3 = hitsback.at(0).x();
    y3 = hitsback.at(0).y();
    z3 = hitsback.at(0).z();
    x4 = hitsback.at(1).x();
    y4 = hitsback.at(1).y();
    z4 = hitsback.at(1).z();
    backanglex = atan2((x4-x3),(z4-z3));
    backangley = atan2((y4-y3),(z4-z3));
    xback << "AngleXBackThreePlanesDoubleDaf";
    yback << "AngleYBackThreePlanesDoubleDaf";
    streamlog_out(DEBUG2) << "Successfully named the histograms:" << endl;
    try{
     dynamic_cast< TH1D* >(_histoThing[xback.str().c_str()])->Fill(backanglex);
    } catch(std::bad_cast &bc){
    streamlog_out(ERROR3) << "Unable to fill histogram: " << xback.str().c_str() << endl;
    } catch(...){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << xback.str().c_str() << ". Unknown Error" << endl;
    }
    try{
      dynamic_cast< TH1D* >(_histoThing[yback.str().c_str()])->Fill(backangley);
    } catch(std::bad_cast &bc){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << yback.str().c_str() << endl;
    } catch(...){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << yback.str().c_str() << ". Unknown Error" << endl;
    }
  }
  streamlog_out(DEBUG2) << "Successfully filled the histograms:" << endl;
 
}

pair< vector< double >, vector< double > > EUTelX0Processor::GetTripleTrackAnglesDoubleDafFitted(Track *frontthree, Track *backthree){
  streamlog_out(DEBUG1) << "Begin pair< vector< double >, vector< double > > EUTelX0Processor::GetTripleTrackAngles(vector< TVector3 > hits)" << endl;
  //Then we find the position where those lines would have hit the DUT
  //THIS IS TO BE DECIDED IF IT IS NEEDED LATER
  //Then we work out the angles of these lines with respect to XZ and YZ, plot results in histograms
  std::vector< TrackerHit* > trackhitsfront = frontthree->getTrackerHits();
  std::vector< TVector3 > hitsfront;
  TVector3 *tempvec = new TVector3();
  for(std::vector< TrackerHit* >::iterator it = trackhitsfront.begin(); it != trackhitsfront.end(); ++it){
    if((*it)->getType() == 32){  //Check if the hit type is the aligned hits or the fitted points (we want fitted points)
      Double_t x = (*it)->getPosition()[0];
      Double_t y = (*it)->getPosition()[1];
      Double_t z = (*it)->getPosition()[2];
      tempvec->SetXYZ(x,y,z);
      hitsfront.push_back(*tempvec);
    } //End of if query for type check
  } //End of for loop running through trackhits
  std::vector< TrackerHit* > trackhitsback = backthree->getTrackerHits();
  std::vector< TVector3 > hitsback;
  for(std::vector< TrackerHit* >::iterator it = trackhitsback.begin(); it != trackhitsback.end(); ++it){
    if((*it)->getType() == 32){  //Check if the hit type is the aligned hits or the fitted points (we want fitted points)
      Double_t x = (*it)->getPosition()[0];
      Double_t y = (*it)->getPosition()[1];
      Double_t z = (*it)->getPosition()[2];
      tempvec->SetXYZ(x,y,z);
      hitsback.push_back(*tempvec);
    } //End of if query for type check
  } //End of for loop running through trackhits
  delete tempvec;
  std::vector< double > scatterx, scattery;
  const double x1 = hitsfront.at(1).x();
  const double y1 = hitsfront.at(1).y();
  const double z1 = hitsfront.at(1).z();
  const double x2 = hitsfront.at(2).x();
  const double y2 = hitsfront.at(2).y();
  const double z2 = hitsfront.at(2).z();
  const double x3 = hitsback.at(0).x();
  const double y3 = hitsback.at(0).y();
  const double z3 = hitsback.at(0).z();
  const double x4 = hitsback.at(1).x();
  const double y4 = hitsback.at(1).y();
  const double z4 = hitsback.at(1).z();
  const double frontanglex = atan2((x2-x1),(z2-z1));
  const double frontangley = atan2((y2-y1),(z2-z1));
  const double backanglex = atan2((x4-x3),(z4-z3));
  const double backangley = atan2((y4-y3),(z4-z3));
  //Push back angles into a vector
  scatterx.push_back(frontanglex);
  scattery.push_back(frontangley);
  scatterx.push_back(backanglex);
  scattery.push_back(backangley);
  //Name the histograms
  std::stringstream xforward,yforward,xyforward;
  streamlog_out(DEBUG2) << "Successfully pushed back the scatterx and y vectors and deleted the TGraphs and TFits:" << endl;
  stringstream xfront,yfront,xback,yback,xyfront,xyback;
  xfront << "AngleXFrontThreePlanesDoubleDaf";
  yfront << "AngleYFrontThreePlanesDoubleDaf";
  xback << "AngleXBackThreePlanesDoubleDaf";
  yback << "AngleYBackThreePlanesDoubleDaf";
  xyfront << "AngleXYFrontThreePlanesDoubleDaf";
  xyback << "AngleXYBackThreePlanesDoubleDaf";
  streamlog_out(DEBUG2) << "Successfully named the histograms:" << endl;
  //Fill the histograms with the XZ and YZ angles
/*
  try{
    dynamic_cast< TH1D* >(_histoThing[xfront.str().c_str()])->Fill(frontanglex);
  } catch(std::bad_cast &bc){
    streamlog_out(ERROR3) << "Unable to fill histogram: " << xfront.str().c_str() << endl;
  } catch(...){
    streamlog_out(ERROR3) << "Unable to fill histogram: " << xfront.str().c_str() << ". Unknown Error" << endl;
  }
  try{
    dynamic_cast< TH1D* >(_histoThing[yfront.str().c_str()])->Fill(frontangley);
  } catch(std::bad_cast &bc){
    streamlog_out(ERROR3) << "Unable to fill histogram: " << yfront.str().c_str() << endl;
  } catch(...){
    streamlog_out(ERROR3) << "Unable to fill histogram: " << yfront.str().c_str() << ". Unknown Error" << endl;
  }
*/
  try{
    dynamic_cast< TH2D* >(_histoThing2D[xyfront.str().c_str()])->Fill(frontanglex,frontangley);
  } catch(std::bad_cast &bc){
    streamlog_out(ERROR3) << "Unable to fill histogram: " << xyfront.str().c_str() << ". Due to bad cast." << endl;
  }
/*
  try{
    dynamic_cast< TH1D* >(_histoThing[xback.str().c_str()])->Fill(backanglex);
  } catch(std::bad_cast &bc){
    streamlog_out(ERROR3) << "Unable to fill histogram: " << xback.str().c_str() << endl;
  } catch(...){
    streamlog_out(ERROR3) << "Unable to fill histogram: " << xback.str().c_str() << ". Unknown Error" << endl;
  }
  try{
    dynamic_cast< TH1D* >(_histoThing[yback.str().c_str()])->Fill(backangley);
  } catch(std::bad_cast &bc){
    streamlog_out(ERROR3) << "Unable to fill histogram: " << yback.str().c_str() << endl;
  } catch(...){
    streamlog_out(ERROR3) << "Unable to fill histogram: " << yback.str().c_str() << ". Unknown Error" << endl;
  }
*/
  try{
    dynamic_cast< TH2D* >(_histoThing2D[xyback.str().c_str()])->Fill(backanglex,backangley);
  } catch(std::bad_cast &bc){
    streamlog_out(ERROR3) << "Unable to fill histogram: " << xyback.str().c_str() << ". Due to bad cast." << endl;
  }
  streamlog_out(DEBUG2) << "Successfully filled the histograms:" << endl;
  pair< vector< double >, vector< double > > bothangles(scatterx,scattery);
  return bothangles;
}

pair< vector< double >, vector< double > > EUTelX0Processor::GetTripleTrackAnglesStraightLines(vector< TVector3 > hits){
  streamlog_out(DEBUG1) << "Begin pair< vector< double >, vector< double > > EUTelX0Processor::GetTripleTrackAngles(vector< TVector3 > hits)" << endl;
  //Then we find the position where those lines would have hit the DUT
  //THIS IS TO BE DECIDED IF IT IS NEEDED LATER
  //Then we work out the angles of these lines with respect to XZ and YZ, plot results in histograms
  const size_t hitsize = hits.size();
  if(hitsize < 6)
  {
    streamlog_out(DEBUG5) << "This track has missing hits, aborting track" << endl;
    string error("In function: EUTelX0Processor::GetTripleTrackAngles(vector< TVector3 > hits), a track has missing hits, move on to the next track");
    throw error;
  }

  std::vector< double > scatterx, scattery;
  streamlog_out(DEBUG3) << "hitsize = " << hitsize << std::endl;
  for(size_t i = 0; i < 2; ++i){ //Fill histograms with the angles of the track at the front and back three planes and put those angles into a vector
    if(i != 0) i = 3;//Move to the next set of planes
    streamlog_out(DEBUG0) << "i = " << i << endl;
    const double x0 = hits.at(i).x();
    const double y0 = hits.at(i).y();
    const double z0 = hits.at(i).z();
    const double x1 = hits.at(i+1).x();
    const double y1 = hits.at(i+1).y();
    const double z1 = hits.at(i+1).z();
    const double x2 = hits.at(i+2).x();
    const double y2 = hits.at(i+2).y();
    const double z2 = hits.at(i+2).z();
    streamlog_out(DEBUG9) << "x0 = " << x0 << endl << "y0 = " << y0 << endl << "z0 = " << z0 << endl << "x2 = " << x2 << endl << "y1 = " << y1 << endl << "z1 = " << z1 << endl << "x2 = " << x2 << endl << "y2 = " << y2 << endl << "z2 = " << z2 << endl;
    int n = 3;
    double arrayx[3] = {x0,x1,x2};
    double arrayy[3] = {y0,y1,y2};
    double arrayz[3] = {z0,z1,z2};
    streamlog_out(DEBUG3) << "Successfully filled the arrays for the track fitting" << endl;
    TGraph *gx = new TGraph(n,arrayz,arrayx);
    TGraph *gy = new TGraph(n,arrayz,arrayy);
    streamlog_out(DEBUG2) << "Successfully created the TGraphs" << endl;
    gErrorIgnoreLevel=kError;
    TF1 *fx = new TF1("fx","[0]*x + [1]",0.0,z2);
    TF1 *fy = new TF1("fy","[0]*x + [1]",0.0,z2);
    streamlog_out(DEBUG2) << "Successfully created the TF1s" << endl;
    gx->Fit(fx,"Q");
    gy->Fit(fy,"Q");
    streamlog_out(DEBUG2) << "Successfully applied the TF1s to the TGraphs" << endl;
    const double frontanglex = atan(fx->GetParameter(0));
    const double frontangley = atan(fy->GetParameter(0));
    streamlog_out(DEBUG2) << "Successfully extracted the paratemers from the Fits:" << endl;
    //Push back angles into a vector
    scatterx.push_back(frontanglex);
    scattery.push_back(frontangley);

    delete fx;
    delete fy;
    delete gx;
    delete gy;
    //Name the histograms
    std::stringstream xforward,yforward,xyforward;
    streamlog_out(DEBUG2) << "Successfully pushed back the scatterx and y vectors and deleted the TGraphs and TFits:" << endl;

    if(i == 0){
      xforward << "AngleXForwardTripleFirstThreePlanes";
      yforward << "AngleYForwardTripleFirstThreePlanes";
      xyforward << "AngleXYForwardTripleFirstThreePlanes";
    } else{
      xforward << "AngleXForwardTripleLastThreePlanes";
      yforward << "AngleYForwardTripleLastThreePlanes";
      xyforward << "AngleXYForwardTripleLastThreePlanes";
    }
    streamlog_out(DEBUG2) << "Successfully named the histograms:" << endl;
    //Fill the histograms with the XZ and YZ angles
    try{
      dynamic_cast< TH1D* >(_histoThing[xforward.str().c_str()])->Fill(frontanglex);
    } catch(std::bad_cast &bc){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << xforward.str().c_str() << endl;
    } catch(...){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << xforward.str().c_str() << ". Unknown Error" << endl;
    }
    try{
      dynamic_cast< TH1D* >(_histoThing[yforward.str().c_str()])->Fill(frontangley);
    } catch(std::bad_cast &bc){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << yforward.str().c_str() << endl;
    } catch(...){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << yforward.str().c_str() << ". Unknown Error" << endl;
    }
    try{
      dynamic_cast< TH2D* >(_histoThing2D[xyforward.str().c_str()])->Fill(frontanglex,frontangley);
    } catch(std::bad_cast &bc){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << xyforward.str().c_str() << ". Due to bad cast." << endl;
    }
    streamlog_out(DEBUG2) << "Successfully filled the histograms:" << endl;
  }
  pair< vector< double >, vector< double > > bothangles(scatterx,scattery);
  return bothangles;
}

void EUTelX0Processor::TriplePlaneTrackScatteringAngles(vector< double > scatterx, vector< double > scattery, vector< TVector3 > hits,  bool doubledaf){
  gStyle->SetOptStat("neMRuo");
  const size_t scatterxsize = scatterx.size();
  streamlog_out(DEBUG3) << "scatterxsize = " << scatterxsize << std::endl;
  //Scattering angle is the forward angle between planes i to i+1 and plane i+1 to i+2
  const double scatteringanglex = scatterx.at(1) - scatterx.at(0);
  const double scatteringangley = scattery.at(1) - scattery.at(0);

  double x1 = hits.at(2).x() + (_dutPosition - hits.at(2).z())*tan(scatterx.at(0));
  double y1 = hits.at(2).y() + (_dutPosition - hits.at(2).z())*tan(scattery.at(0));
  //Name the histograms
  std::stringstream ssscatterx,ssscattery,ssscatterxy;
  if(doubledaf == true){
    ssscatterx << "ScatteringAngleXDoubleDaf";
    ssscattery << "ScatteringAngleYDoubleDaf";
    ssscatterxy << "ScatteringAngleXYDoubleDaf";
    _dutPosition = hits.at(2).z() + 5.1; //TODO HACK just put in the value for now from run 8876 whilst I test, need to add this variable to my runs.csv list
    double x2 = hits.at(9).x() + (_dutPosition - hits.at(9).z())*tan(scatterx.at(1));
    double y2 = hits.at(9).y() + (_dutPosition - hits.at(9).z())*tan(scattery.at(1));

    streamlog_out(DEBUG5) << "x1 = " << x1 << endl;
    streamlog_out(DEBUG5) << "x2 = " << x2 << endl;
    streamlog_out(DEBUG5) << "y1 = " << y1 << endl;
    streamlog_out(DEBUG5) << "y2 = " << y2 << endl;
    streamlog_out(DEBUG5) << "scatterx(0) = " << scatterx.at(0) << endl;
    streamlog_out(DEBUG5) << "scatterx(1) = " << scatterx.at(1) << endl;
    streamlog_out(DEBUG5) << "dutPosition = " << _dutPosition << endl;
    streamlog_out(DEBUG5) << "_cut = " << _cut << endl << endl;

    if(!(x1 < x2 + _cut && x1 > x2 - _cut && y1 < y2 + _cut && y1 > y2 - _cut)){
      //These two tracks don't come from the same real track
      return;
    }
    totaltracknumberdoubledaf++;
  } else{
    ssscatterx << "ScatteringAngleXTriple";
    ssscattery << "ScatteringAngleYTriple";
    ssscatterxy << "ScatteringAngleXYTriple";
    ++totaltracknumber;
  }

  pair< int, int > position(ConversionHitmapToX0map(x1,y1));
  if(doubledaf == true){
    ScatteringAngleXTripleMapDataDaf[position].push_back(scatteringanglex);
    ScatteringAngleYTripleMapDataDaf[position].push_back(scatteringangley);
  } else{
    cout << "filling the map data" << endl;
    ScatteringAngleXTripleMapData[position].push_back(scatteringanglex);
    ScatteringAngleYTripleMapData[position].push_back(scatteringangley);
  }
  try{
    dynamic_cast< TH1D* >(_histoThing[ssscatterx.str().c_str()])->Fill(scatteringanglex);
  } catch(std::bad_cast &bc){
    streamlog_out(ERROR3) << "Unable to fill histogram: " << ssscatterx.str().c_str() << endl;
  }
  try{
    dynamic_cast< TH1D* >(_histoThing[ssscattery.str().c_str()])->Fill(scatteringangley);
  } catch(std::bad_cast &bc){
    streamlog_out(ERROR3) << "Unable to fill histogram: " << ssscattery.str().c_str() << endl;
  }
  try{
    dynamic_cast< TH2D* >(_histoThing2D[ssscatterxy.str().c_str()])->Fill(scatteringanglex,scatteringangley);
  } catch(std::bad_cast &bc){
    streamlog_out(ERROR3) << "Unable to fill histogram: " << ssscatterxy.str().c_str() << endl;
  }
}

void EUTelX0Processor::SinglePlaneTrackScatteringAngles(vector< double > scatterx, vector< double > scattery, std::vector< TVector3 > hits){
  gStyle->SetOptStat("neMRuo");
  const size_t scatterxsize = scatterx.size();
  streamlog_out(DEBUG3) << "scatterxsize = " << scatterxsize << std::endl;
  for(size_t i = 0; i < scatterxsize-1; ++i){ //This fills the scattering angle histograms plane by plane
    //Scattering angle is the forward angle between planes i to i+1 and plane i+1 to i+2
    const double scatteringanglex = scatterx[i+1] - scatterx[i];
    const double scatteringangley = scattery[i+1] - scattery[i];

    //Name the histograms
    std::stringstream ssscatterx,ssscattery,ssscatterxy;
    ssscatterx << "ScatteringAngleXPlane" << i+1;
    ssscattery << "ScatteringAngleYPlane" << i+1;
    ssscatterxy << "KinkAnglePlane" << i+1;
    
    double x = hits[i+1].x();
    double y = hits[i+1].y();
    pair< int, int > position(ConversionHitmapToX0map(x,y));
    ScatteringAngleXSingleMapData[i+1][position].push_back(scatteringanglex);
    ScatteringAngleYSingleMapData[i+1][position].push_back(scatteringangley);
    try{
      dynamic_cast< TH1D* >(_histoThing[ssscatterx.str().c_str()])->Fill(scatteringanglex);
    } catch(std::bad_cast &bc){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << ssscatterx.str().c_str() << endl;
    }
    try{
      dynamic_cast< TH1D* >(_histoThing[ssscattery.str().c_str()])->Fill(scatteringangley);
    } catch(std::bad_cast &bc){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << ssscattery.str().c_str() << endl;
    }
    try{
      dynamic_cast< TH2D* >(_histoThing2D[ssscatterxy.str().c_str()])->Fill(scatteringanglex,scatteringangley);
    } catch(std::bad_cast &bc){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << ssscatterxy.str().c_str() << endl;
    }
  }
}

void EUTelX0Processor::kinkEstimate(Track* track){
  //This function works out an angle based on a straight line fitted from plane 0 to 2 and plane 5 to 3
  //It will also store all other angles in histograms too
  streamlog_out(DEBUG0) << "Running function kinkEstimate(Track* " << &track << ")" << std::endl;
  
  //First we extract the relevant hits from the track
  try{
    std::vector< TVector3 > hits = getHitsFromTrack(track);
  } catch(...)
  { 
    streamlog_out(ERROR3) << "Could not get hits from track" << endl;
    return;
  }
  std::vector< TVector3 > hits = getHitsFromTrack(track);

  streamlog_out(DEBUG3) << "Successfully got hits from track" << std::endl;
  try{
    pair< vector< double >, vector< double > > scatterpairsingle(GetSingleTrackAngles(hits));
    pair< vector< double >, vector< double > > scatterpairtriple(GetTripleTrackAnglesStraightLines(hits));
    std::vector< double > scatterxsingle = scatterpairsingle.first, scatterysingle = scatterpairsingle.second;
    std::vector< double > scatterxtriple = scatterpairtriple.first, scatterytriple = scatterpairtriple.second;
    SinglePlaneTrackScatteringAngles(scatterxsingle, scatterysingle, hits);
    TriplePlaneTrackScatteringAngles(scatterxtriple, scatterytriple, hits, false);
  } catch(...){
    streamlog_out(MESSAGE0) << "Caught an exception on line " << __LINE__ << " in file " << __FILE__ << endl;;
    return;
  }
  streamlog_out(DEBUG3) << "Made it to the end of kinkEstimate()" << endl;
}

void EUTelX0Processor::processEvent(LCEvent *evt)
{
  totaleventnumber++;
  if(totaleventnumber%1000 == 0){
    cout << "Total Number of Events processed: " << totaleventnumber << endl;
    cout << "Total Number of Tracks processed: " << totaltracknumber << endl;
    cout << "Total Number of Double Daf Tracks processed: " << totaltracknumberdoubledaf << endl;
  }
  streamlog_out(DEBUG0) << "Running EUTelX0Processor::processEvent(LCEvent *evt) with evt = " << evt << std::endl;
  //Take track from input parameter
  //Work out kink angle from track
  //Put kink angle into histogram
  //Check for last event
    //Fit gaussian to histogram
    //Extract sigma value from histogram
    try{
      LCCollection* trackcollection1 = evt->getCollection(_trackCollectionName1);
      int elementnumber1 = trackcollection1->getNumberOfElements();
      streamlog_out(DEBUG2) << "Created trackcollection1, it has " << elementnumber1 << " elements" << endl;
      LCCollection* trackcollection2 = evt->getCollection(_trackCollectionName2);
      int elementnumber2 = trackcollection2->getNumberOfElements();
      streamlog_out(DEBUG2) << "Created trackcollection2, it has " << elementnumber2 << " elements" << endl;
      for(int i = 0; i < elementnumber1; ++i){
        for(int j = 0; j < elementnumber2; ++j){
          streamlog_out(DEBUG0) << "i = " << i << ", j = " << j << endl;
          Track *eventtrack1 = dynamic_cast< Track* >(trackcollection1->getElementAt(i));
          Track *eventtrack2 = dynamic_cast< Track* >(trackcollection2->getElementAt(j));
          streamlog_out(DEBUG0) << "Dynamically casted tracks: " << &eventtrack1 << " and " << &eventtrack2 << endl;
          std::vector< TrackerHit* > trackhits1 = eventtrack1->getTrackerHits();
          streamlog_out(DEBUG0) << "Got the tracker hit 1" << endl;
          std::vector< TrackerHit* > trackhits2 = eventtrack2->getTrackerHits();
          streamlog_out(DEBUG0) << "Got the tracker hit 2" << endl;
          std::vector< TVector3 > hits;
          int k(0);
          TVector3 *tempvec = new TVector3();
          for(std::vector< TrackerHit* >::iterator it = trackhits1.begin(); it != trackhits1.end(); ++it){
            if((*it)->getType() == 32){  //Check if the hit type is the aligned hits or the fitted points (we want fitted points)
              streamlog_out(DEBUG0) << "k = " << k << endl;
              Double_t x = (*it)->getPosition()[0];
              Double_t y = (*it)->getPosition()[1];
              Double_t z = (*it)->getPosition()[2];
              k++;
	      tempvec->SetXYZ(x,y,z);
              hits.push_back(*tempvec);
            } //End of if query for type check
          } //End of for loop running through trackhits
          streamlog_out(DEBUG0) << "Filled hits from first track" << endl;
          for(std::vector< TrackerHit* >::iterator it = trackhits2.begin(); it != trackhits2.end(); ++it){
            if((*it)->getType() == 32){  //Check if the hit type is the aligned hits or the fitted points (we want fitted points)
              streamlog_out(DEBUG0) << "k = " << k << endl;
              Double_t x = (*it)->getPosition()[0];
              Double_t y = (*it)->getPosition()[1];
              Double_t z = (*it)->getPosition()[2];
              k++;
	      tempvec->SetXYZ(x,y,z);
              hits.push_back(*tempvec);
            } //End of if query for type check
          } //End of for loop running through trackhits
          streamlog_out(DEBUG0) << "Filled hits from second track" << endl;
          delete tempvec;
          pair< vector< double >, vector< double > > scatterpairtripledaf(GetTripleTrackAnglesDoubleDafFitted(eventtrack1, eventtrack2));
          streamlog_out(DEBUG0) << "Worked out track angles" << endl;
          TriplePlaneTrackScatteringAngles(scatterpairtripledaf.first, scatterpairtripledaf.second, hits, true);
          streamlog_out(DEBUG0) << "Worked out scattering angles" << endl;
        }
      }
      for(int i = 0; i < elementnumber1; ++i){
        Track *eventtrack1 = dynamic_cast< Track* >(trackcollection1->getElementAt(i));
        PlotTripleTrackAngleDoubleDafFitted(eventtrack1,true);
      }
      for(int j = 0; j < elementnumber2; ++j){
        Track *eventtrack2 = dynamic_cast< Track* >(trackcollection2->getElementAt(j));
        PlotTripleTrackAngleDoubleDafFitted(eventtrack2,false);
      }
    } catch(...){
      streamlog_out(MESSAGE3) << "Double Daf fitting didn't work..." << endl;
    }
    try{
      //_referenceHitVec = evt->getCollection(_referenceHitCollectionName);
      LCCollection* trackcollection = evt->getCollection(_trackCollectionName);
      int elementnumber = trackcollection->getNumberOfElements();
      for(int i = 0; i < elementnumber; ++i){
        Track* eventtrack = dynamic_cast< Track* >(trackcollection->getElementAt(i));
        streamlog_out(DEBUG0) << "Here is all the information about the track in run " << _runNumber << ", event " << _eventNumber << ", element " << i << std::endl << std::endl;
        printTrackParameters( eventtrack );
        kinkEstimate( eventtrack );
  //      pair< vector< double >, vector< double > > scatterpairtripledaf(GetTripleTrackAnglesDoubleDafFitted(track1, track2));
  //      TriplePlaneTrackScatteringAngles(scatterxtripledaf, scatterytripledaf, hits, true);
  //      singlePointResolution( eventtrack );
  //      threePointResolution( eventtrack );
      }
    } catch(DataNotAvailableException &datanotavailable){
      streamlog_out(DEBUG4) << "Exception occured: " << datanotavailable.what() << std::endl
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

std::pair<double,double> EUTelX0Processor::GetLowerAndUpperBounds(TH1D *temphisto){
  double lowrange(-9999), highrange(9999);
  double maxx = -999999;
  for(int i = 2; i < nobinsangle-2; ++i){
    if(maxx < temphisto->GetBinContent(i)){
      maxx = temphisto->GetBinContent(i);
    }
  }
  cout << "maxx = " << maxx << endl;
  for(int i = 2; i < nobinsangle-2; ++i){
    if(temphisto->GetBinContent(i) > maxx/2.0){
      lowrange = temphisto->GetBinCenter(i);
      break;
    }
  }
  if(lowrange == -9999){
    streamlog_out(WARNING3) << "Did not find low range" << endl;
  }
  for(int ii = nobinsangle-2; ii >= 2; --ii){
    if(temphisto->GetBinContent(ii) > maxx/2.0){
      highrange = temphisto->GetBinCenter(ii);
      break;
    }
  }
  if(highrange == -9999){
    streamlog_out(WARNING3) << "Did not find high range" << endl;
  }
  std::pair<double,double> ranges(lowrange,highrange);
  return ranges;
}

double EUTelX0Processor::getSigma(vector<double> vec){
  double mean(0);
  TH1D *temphisto = new TH1D("temphisto","temphisto",nobinsangle,minbinangle,maxbinangle);
  for(vector< double >::iterator it = vec.begin(); it != vec.end(); ++it){
    temphisto->Fill(*it);
//    mean += *it;
  }
  mean /= vec.size();
  double sigma(0);
//  for(vector< double >::iterator it = vec.begin(); it != vec.end(); ++it){
//    sigma += (mean - *it)*(mean - *it);
//  }
//  sigma /= static_cast<double>(vec.size());
//  sigma = sqrt(sigma);
  gErrorIgnoreLevel=kError;
  std::pair<double,double> range(GetLowerAndUpperBounds(temphisto));
  TF1 *fit = new TF1("fit","1/(sqrt(2*3.1415)*[0])*exp((-x*x)/(2*[0]*[0]))*[1]",range.first,range.second);
  int goodfit = temphisto->Fit(fit, "QRME");
  if(goodfit != 0){
    return 0;
  }
  double sigmafit = fit->GetParameter(2);
  delete fit;
  delete temphisto;
//  if(sigma != sigma){
//    return 0;
//  } else{
  return sigmafit;
//  }
}

void EUTelX0Processor::end()
{
  //Clean up memory
  //Set all values to zero or NULL
  //
  //calculateX0();
  int numberofplanes = 6;
  gErrorIgnoreLevel=kError;
  for(double i = minx; i <= maxx; i += binsizex)
  {
    streamlog_out(MESSAGE6) << "The Radiation Length Map has made it to " << i << endl;
    for(double j = miny; j <= maxy; j += binsizey)
    {
      pair< int, int > position(ConversionHitmapToX0map(i,j));
      for(int k = 1; k < numberofplanes - 1; ++k){
        double sigmasinglex = getSigma(ScatteringAngleXSingleMapData[k][position]);
        ScatteringAngleXSingleMap[k]->Fill(i,j,sigmasinglex);
        double sigmasingley = getSigma(ScatteringAngleYSingleMapData[k][position]);;
        ScatteringAngleYSingleMap[k]->Fill(i,j,sigmasingley);

        double x0single = pow(sqrt(pow(sigmasinglex,2)+pow(sigmasingley,2))*1000.0*_beamEnergy/13.6,1.0/0.555);
        RadiationLengthSingleMap[k]->Fill(i,j,x0single);

      }
      double sigmatriplex = getSigma(ScatteringAngleXTripleMapData[position]);
      ScatteringAngleXTripleMap->Fill(i,j,sigmatriplex);
      double sigmatripley = getSigma(ScatteringAngleYTripleMapData[position]);
      ScatteringAngleYTripleMap->Fill(i,j,sigmatripley);
      
      double x0triple = pow(sqrt(pow(sigmatriplex,2)+pow(sigmatripley,2))*1000.0*_beamEnergy/13.6,1.0/0.555);
      RadiationLengthTripleMap->Fill(i,j,x0triple);

      if(_doubleDafFitted){
        double sigmatriplexdaf = getSigma(ScatteringAngleXTripleMapDataDaf[position]);
        ScatteringAngleXTripleMapDoubleDaf->Fill(i,j,sigmatriplexdaf);
        double sigmatripleydaf = getSigma(ScatteringAngleYTripleMapDataDaf[position]);
        ScatteringAngleYTripleMapDoubleDaf->Fill(i,j,sigmatripleydaf);
  
        double x0tripledaf = pow(sqrt(pow(sigmatriplexdaf,2)+pow(sigmatripleydaf,2))*1000.0*_beamEnergy/13.6,1.0/0.555);
        RadiationLengthMapDoubleDaf->Fill(i,j,x0tripledaf);
      }
    }
  }
  cout << "Total Number of Events processed: " << totaleventnumber << endl;
  cout << "Total Number of Tracks processed: " << totaltracknumber << endl;
  cout << "Total Number of Double Daf Tracks processed: " << totaltracknumberdoubledaf << endl;
  _histoThing.clear(); 
  _inputTrackCollectionVec = NULL; 
}

pair< double, double > EUTelX0Processor::ConversionX0mapToHitmap(int x, int y){
  double newx = minx + binsizex*x;
  double newy = miny + binsizey*y;
  pair< double, double > position(newx,newy);
  return position;
}

pair< int, int > EUTelX0Processor::ConversionHitmapToX0map(double x, double y){
  pair< int, int > position;
  double newminbinx(minx);
  double newminbiny(miny);
  int actualbinx(0);
  int actualbiny(0);
  bool notfoundx = true;
  bool notfoundy = true;
  
  while(notfoundx){ //Begin search for which bin of the radiation length map the track is going through
    if(x > maxx){
      string xistoobig = "The x coordinate of the hit is larger than the size of the sensor. Skipping this track";
      throw xistoobig;
    } //end of: if(x > maxx)
    else if(x < newminbinx){
      string xistoosmall = "The x coordinate of the hit is either smaller than the size of the sensor, or should have been caught in the previous bin. Skipping this track";
      throw xistoosmall;
    } //end of: else if(x < newminbinx)
    else if(x < newminbinx + binsizex){
      notfoundx = false;
 
      while(notfoundy){ //Begin search for which bin of the radiation length map the track is going through
        if(y > maxy){
          string yistoobig = "The y coordinate of the hit is larger than the size of the sensor. Skipping this track";
          throw yistoobig;
        } //end of: if(y > mayy)
        else if(y < newminbiny){
          string yistoosmall = "The y coordinate of the hit is either smaller than the size of the sensor, or should have been caught in the previous bin. Skipping this track";
          throw yistoosmall;
        } //end of: else if(y < newminbiny)
        else if(y < newminbiny + binsizey){
          notfoundy = false;
          position = pair<int, int>(actualbinx,actualbiny);
        } //end of: else if(y < newminbiny + binsizey)
        else{
          newminbiny += binsizey;
          actualbiny++;
        } //end of: else [else if(y < newminbiny + binsizey)]
      } // end of: while(notfoundy)
    } //end of: else if(x < newminbinx + binsizex)
    else{
      newminbinx += binsizex;
      actualbinx++;
    } //end of: else [else if(x < newminbinx + binsizex)]
  } // end of: while(notfoundx)
  return position;
}

std::vector< TVector3 > EUTelX0Processor::getHitsFromTrack(Track *track){
  streamlog_out(DEBUG0) << "getHitsFromTrack function called" << endl;
 std::vector< TrackerHit* > trackhits = track->getTrackerHits();
 std::vector< TVector3 > hits;
 TVector3 *tempvec = new TVector3();
 for(std::vector< TrackerHit* >::iterator it = trackhits.begin(); it != trackhits.end(); ++it){
   streamlog_out(DEBUG0) << "Getting the hits" << endl;
    if((*it)->getType() < 32){  //Check if the hit type is the aligned hits or the fitted points (we want aligned hits)
      Double_t x = (*it)->getPosition()[0];
      streamlog_out(DEBUG0) << "x=" << x << endl;
      Double_t y = (*it)->getPosition()[1];
      streamlog_out(DEBUG0) << "y=" << y << endl;
      Double_t z = (*it)->getPosition()[2];
      streamlog_out(DEBUG0) << "z=" << z << endl;
      tempvec->SetXYZ(x,y,z);
      hits.push_back(*tempvec);
    } //End of if query for type check
  } //End of for loop running through trackhits
  return hits;
}

void EUTelX0Processor::singlePointResolution(Track *track){
//This function draws a line between two hits on adjacent planes and then projects in the forward and backward direction to the hit on the next plane.
//It then compares that hit with the actual hit on that next plane and the differencce is plotted
  streamlog_out(DEBUG1) << "Function EUTeoX0Processor::singlePointResolution(Track *" << &track << ") called" << std::endl;
  std::vector< TVector3 > hits = getHitsFromTrack(track);
  int i = 0;
  for(std::vector< TVector3 >::iterator it = hits.begin(); it != hits.end(); ++it){
    streamlog_out(DEBUG0) << "Entering for loop for Plane " << i << endl;
    stringstream histogramnamex,histogramnamey;
    histogramnamex << "SinglePointResidualXPlane" << i;
    histogramnamey << "SinglePointResidualYPlane" << i;
    double x0, y0, x1, y1, x2, y2, z0, z1, z2;
    x0 = (*it).x();
    y0 = (*it).y();
    z0 = (*it).z();
    x1 = (*it + 1).x();
    y1 = (*it + 1).y();
    z1 = (*it + 1).z();
    x2 = (*it + 2).x();
    y2 = (*it + 2).y();
    z2 = (*it + 2).z();
    double predictx, predicty;
    if(it == hits.begin() || it == hits.begin() + 1){
      predictx = (x2-x1)*(z2-z0)/(z2-z1);
      predicty = (y2-y1)*(z2-z0)/(z2-z1);
      try{
        dynamic_cast< TH1D* > (_histoThing[histogramnamex.str().c_str()])->Fill(x0 - predictx);
      } catch(std::bad_cast &bc){
        streamlog_out(ERROR3) << "Unable to fill histogram: " << histogramnamex.str().c_str() << endl;
      }
      try{
        dynamic_cast< TH1D* > (_histoThing[histogramnamey.str().c_str()])->Fill(y0 - predicty);
      } catch(std::bad_cast &bc){
        streamlog_out(ERROR3) << "Unable to fill histogram: " << histogramnamey.str().c_str() << endl;
      }
    } else {
      predictx = (x1-x0)*(z2-z0)/(z1-z0);
      predicty = (y1-y0)*(z2-z0)/(z1-z0);
       try{
        dynamic_cast< TH1D* > (_histoThing[histogramnamex.str().c_str()])->Fill(x2 - predictx);
      } catch(std::bad_cast &bc){
        streamlog_out(ERROR3) << "Unable to fill histogram: " << histogramnamex.str().c_str() << endl;
      }
      try{
        dynamic_cast< TH1D* > (_histoThing[histogramnamey.str().c_str()])->Fill(y2 - predicty);
      } catch(std::bad_cast &bc){
        streamlog_out(ERROR3) << "Unable to fill histogram: " << histogramnamey.str().c_str() << endl;
      }     
    }
    ++i;
  }
}

void EUTelX0Processor::threePointResolution(Track *track){
//This function draws a line between the hits in plane i and i+2
//then it compares where the hit in plane i+1 is with the average of the other two planes
//this is then plotted as a residual for each plane.
  streamlog_out(DEBUG1) << "Function EUTelX0Processor::threePointResolution(Track *" << &track << ") called" << std::endl;

  std::vector< TrackerHit* > trackhits = track->getTrackerHits();
  std::vector< TVector3 > hits = getHitsFromTrack(track);
  int i = 1;
  for(std::vector< TVector3 >::iterator it = hits.begin(); it != hits.end() - 2; ++it){
    //This determines the guess of the position of the particle as it should hit the middle sensor
    streamlog_out(DEBUG1) << "In the for loop of the three point resolution function" << endl;
    double averagingfactor = ((*it + 1).z() - (*it).z())/((*it + 2).z() - (*it).z());
    if(averagingfactor == 0){
      streamlog_out(ERROR5) << "Averaging factor == 0)" << endl;
      continue;
    }
    double averagex = ((*it).x() + (*it + 2).x())/averagingfactor;
    double averagey = ((*it).y() + (*it + 2).x())/averagingfactor;
    double middlex = (*it + 1).x();
    double middley = (*it + 1).y();
    std::stringstream ResidualX, ResidualY;
    ResidualX << "ThreePointResidualXPlane" << i;
    ResidualX << "ThreePointResidualYPlane" << i;
    try{
      dynamic_cast< TH1D* > (_histoThing[ResidualX.str().c_str()])->Fill(averagex - middlex);
    } catch(std::bad_cast& bc){
      streamlog_out(ERROR5) << "Unable to fill histogram " << ResidualX.str() << " due to bad cast" << std::endl;
    }
    try{
      dynamic_cast< TH1D* > (_histoThing[ResidualY.str().c_str()])->Fill(averagey - middley);
    } catch(std::bad_cast& bc){
      streamlog_out(ERROR5) << "Unable to fill histogram " << ResidualY.str() << " due to bad cast" << std::endl;
    }
    streamlog_out(DEBUG1) << "Made it to the end of the for loop" << endl;
    i++;
  }
}

