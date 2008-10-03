// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Philipp Roloff, DESY <mailto:philipp.roloff@desy.de>
// Version: $Id: EUTelMille.cc,v 1.30 2008-10-03 16:09:07 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
// built only if GEAR and MARLINUTIL are used
#if defined(USE_GEAR) && defined(USE_MARLINUTIL)

// eutelescope includes ".h"
#include "EUTelMille.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelSparseCluster2Impl.h"
#include "EUTelExceptions.h"
#include "EUTelPStream.h"
#include "EUTelAlignmentConstant.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"
#include "marlin/Exceptions.h"
#include "marlin/AIDAProcessor.h"

// marlin util includes
#include "mille/Mille.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <IO/LCWriter.h>
#include <UTIL/LCTime.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>

// ROOT includes
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
# include <TRandom.h>
#endif

// system includes <>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;
using namespace lcio;
using namespace marlin;
// using namespace gear;
using namespace eutelescope;

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelMille::_numberTracksLocalname   = "NumberTracks";
std::string EUTelMille::_chi2XLocalname          = "Chi2X";
std::string EUTelMille::_chi2YLocalname          = "Chi2Y";
std::string EUTelMille::_residualXLocalname      = "ResidualX";
std::string EUTelMille::_residualYLocalname      = "ResidualY";
#endif

EUTelMille::EUTelMille () : Processor("EUTelMille") {

  // modify processor description
  _description =
    "EUTelMille uses the MILLE program to write data files for MILLEPEDE II.";

  // choose input mode
  registerOptionalParameter("InputMode","Selects the source of input hits. 0 - hits read from hitfile with simple trackfinding. 1 - hits read from output of tracking processor. 2 - Test mode. Simple internal simulation and simple trackfinding.",_inputMode, static_cast <int> (0));

  // input collections

  registerInputCollection(LCIO::TRACKERHIT,"HitCollectionName",
                          "Hit collection name",
                          _hitCollectionName,std::string("hit"));

  registerInputCollection(LCIO::TRACK,"TrackCollectionName",
                          "Track collection name",
                          _trackCollectionName,std::string("fittracks"));

  // parameters

  registerOptionalParameter("DistanceMax","Maximal allowed distance between hits entering the fit per 10 cm space between the planes.",
                            _distanceMax, static_cast <float> (2000.0));

  std::vector<int> exclPlanes;

  registerOptionalParameter("ExcludePlanes","Exclude planes from fit."
                            ,_excludePlanes , exclPlanes);

  registerOptionalParameter("MaxTrackCandidates","Maximal number of track candidates."
                            ,_maxTrackCandidates, static_cast <int> (2000));

  registerOptionalParameter("BinaryFilename","Name of the Millepede binary file."
                            ,_binaryFilename, string ("mille.bin"));

  registerOptionalParameter("TelescopeResolution","Resolution of the telescope for Millepede."
                            ,_telescopeResolution, static_cast <float> (3.0));

  registerOptionalParameter("OnlySingleHitEvents","Use only events with one hit in every plane."
                            ,_onlySingleHitEvents, static_cast <int> (0));

  registerOptionalParameter("OnlySingleTrackEvents","Use only events with one track candidate."
                            ,_onlySingleTrackEvents, static_cast <int> (0));

  registerOptionalParameter("AlignMode","Number of alignment constants used. Available mode are: 1 - shifts in the X and Y directions and a rotation around the Z axis, 2 - only shifts in the X and Y directions",
                            _alignMode, static_cast <int> (1));

  registerOptionalParameter("UseResidualCuts","Use cuts on the residuals to reduce the combinatorial background. 0 for off (default), 1 for on"
                            ,_useResidualCuts, static_cast <int> (0));

  registerOptionalParameter("AlignmentConstantLCIOFile","This is the name of the LCIO file name with the output alignment"
                            "constants (add .slcio)",
                            _alignmentConstantLCIOFile, static_cast< string > ( "alignment.slcio" ) );

  FloatVec MinimalResidualsX;
  MinimalResidualsX.push_back(0.0);
  MinimalResidualsX.push_back(0.0);
  MinimalResidualsX.push_back(0.0);
  MinimalResidualsX.push_back(0.0);
  MinimalResidualsX.push_back(0.0);
  MinimalResidualsX.push_back(0.0);

  FloatVec MinimalResidualsY;
  MinimalResidualsY.push_back(0.0);
  MinimalResidualsY.push_back(0.0);
  MinimalResidualsY.push_back(0.0);
  MinimalResidualsY.push_back(0.0);
  MinimalResidualsY.push_back(0.0);
  MinimalResidualsY.push_back(0.0);

  FloatVec MaximalResidualsX;
  MaximalResidualsX.push_back(0.0);
  MaximalResidualsX.push_back(0.0);
  MaximalResidualsX.push_back(0.0);
  MaximalResidualsX.push_back(0.0);
  MaximalResidualsX.push_back(0.0);
  MaximalResidualsX.push_back(0.0);

  FloatVec MaximalResidualsY;
  MaximalResidualsY.push_back(0.0);
  MaximalResidualsY.push_back(0.0);
  MaximalResidualsY.push_back(0.0);
  MaximalResidualsY.push_back(0.0);
  MaximalResidualsY.push_back(0.0);
  MaximalResidualsY.push_back(0.0);

  registerOptionalParameter("ResidualsXMin","Minimal values of the hit residuals in the X direction for a track"
                            ,_residualsXMin,MinimalResidualsX);

  registerOptionalParameter("ResidualsYMin","Minimal values of the hit residuals in the Y direction for a track"
                            ,_residualsYMin,MinimalResidualsY);

  registerOptionalParameter("ResidualsXMax","Maximal values of the hit residuals in the X direction for a track"
                            ,_residualsXMax,MaximalResidualsX);

  registerOptionalParameter("ResidualsYMax","Maximal values of the hit residuals in the Y direction for a track"
                            ,_residualsYMax,MaximalResidualsY);

  registerOptionalParameter("GeneratePedeSteerfile","Generate a steering file for the pede program."
                            ,_generatePedeSteerfile, static_cast <int> (0));

  registerOptionalParameter("PedeSteerfileName","Name of the steering file for the pede program."
                            ,_pedeSteerfileName, string("steer_mille.txt"));

  registerOptionalParameter("RunPede","Execute the pede program using the generated steering file."
                            ,_runPede, static_cast <int> (0));

  registerOptionalParameter("UsePedeUserStartValues","Give start values for pede by hand (0 - automatic calculation of start values, 1 - start values defined by user)."
                            ,_usePedeUserStartValues, static_cast <int> (0));

  FloatVec PedeUserStartValuesX;
  PedeUserStartValuesX.push_back(0.0);
  PedeUserStartValuesX.push_back(0.0);
  PedeUserStartValuesX.push_back(0.0);
  PedeUserStartValuesX.push_back(0.0);
  PedeUserStartValuesX.push_back(0.0);
  PedeUserStartValuesX.push_back(0.0);

  registerOptionalParameter("PedeUserStartValuesX","Start values for the alignment for shifts in the X direction."
                            ,_pedeUserStartValuesX,PedeUserStartValuesX);

  FloatVec PedeUserStartValuesY;
  PedeUserStartValuesY.push_back(0.0);
  PedeUserStartValuesY.push_back(0.0);
  PedeUserStartValuesY.push_back(0.0);
  PedeUserStartValuesY.push_back(0.0);
  PedeUserStartValuesY.push_back(0.0);
  PedeUserStartValuesY.push_back(0.0);

  registerOptionalParameter("PedeUserStartValuesY","Start values for the alignment for shifts in the Y direction."
                            ,_pedeUserStartValuesY,PedeUserStartValuesY);

  FloatVec PedeUserStartValuesGamma;
  PedeUserStartValuesGamma.push_back(0.0);
  PedeUserStartValuesGamma.push_back(0.0);
  PedeUserStartValuesGamma.push_back(0.0);
  PedeUserStartValuesGamma.push_back(0.0);
  PedeUserStartValuesGamma.push_back(0.0);
  PedeUserStartValuesGamma.push_back(0.0);

  registerOptionalParameter("PedeUserStartValuesGamma","Start values for the alignment for the angle gamma."
                            ,_pedeUserStartValuesGamma,PedeUserStartValuesGamma);

  registerOptionalParameter("TestModeSensorResolution","Resolution assumed for the sensors in test mode."
                            ,_testModeSensorResolution, static_cast <float> (3.0));

  registerOptionalParameter("TestModeXTrackSlope","Width of the track slope distribution in the x direction"
                            ,_testModeXTrackSlope, static_cast <float> (0.0005));

  registerOptionalParameter("TestModeYTrackSlope","Width of the track slope distribution in the y direction"
                            ,_testModeYTrackSlope, static_cast <float> (0.0005));

  FloatVec SensorZPositions;
  SensorZPositions.push_back(20000.0);
  SensorZPositions.push_back(40000.0);
  SensorZPositions.push_back(60000.0);
  SensorZPositions.push_back(80000.0);
  SensorZPositions.push_back(100000.0);
  SensorZPositions.push_back(120000.0);

  registerOptionalParameter("TestModeSensorZPositions","Z positions of the sensors in test mode."
                            ,_testModeSensorZPositions,SensorZPositions);

  FloatVec SensorXShifts;
  SensorXShifts.push_back(0.0);
  SensorXShifts.push_back(0.0);
  SensorXShifts.push_back(0.0);
  SensorXShifts.push_back(0.0);
  SensorXShifts.push_back(0.0);
  SensorXShifts.push_back(0.0);

  FloatVec SensorYShifts;
  SensorYShifts.push_back(0.0);
  SensorYShifts.push_back(0.0);
  SensorYShifts.push_back(0.0);
  SensorYShifts.push_back(0.0);
  SensorYShifts.push_back(0.0);
  SensorYShifts.push_back(0.0);

  registerOptionalParameter("TestModeSensorXShifts","X shifts of the sensors in test mode (to be determined by the alignment)."
                            ,_testModeSensorXShifts,SensorXShifts);

  registerOptionalParameter("TestModeSensorYShifts","Y shifts of the sensors in test mode (to be determined by the alignment)."
                            ,_testModeSensorYShifts,SensorYShifts);

  FloatVec SensorGamma;
  SensorGamma.push_back(0.0);
  SensorGamma.push_back(0.0);
  SensorGamma.push_back(0.0);
  SensorGamma.push_back(0.0);
  SensorGamma.push_back(0.0);
  SensorGamma.push_back(0.0);

  registerOptionalParameter("TestModeSensorGamma","Rotation around the z axis of the sensors in test mode (to be determined by the alignment)."
                            ,_testModeSensorGamma,SensorGamma);

  FloatVec SensorAlpha;
  SensorAlpha.push_back(0.0);
  SensorAlpha.push_back(0.0);
  SensorAlpha.push_back(0.0);
  SensorAlpha.push_back(0.0);
  SensorAlpha.push_back(0.0);
  SensorAlpha.push_back(0.0);

  registerOptionalParameter("TestModeSensorAlpha","Rotation around the x axis of the sensors in test mode (to be determined by the alignment)."
                            ,_testModeSensorAlpha,SensorAlpha);

  FloatVec SensorBeta;
  SensorBeta.push_back(0.0);
  SensorBeta.push_back(0.0);
  SensorBeta.push_back(0.0);
  SensorBeta.push_back(0.0);
  SensorBeta.push_back(0.0);
  SensorBeta.push_back(0.0);

  registerOptionalParameter("TestModeSensorBeta","Rotation around the y axis of the sensors in test mode (to be determined by the alignment)."
                            ,_testModeSensorBeta,SensorBeta);

}

void EUTelMille::init() {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

  // Initialize number of excluded planes
  _nExcludePlanes = _excludePlanes.size();

  streamlog_out ( MESSAGE2 ) << "Number of planes excluded from the alignment fit: " << _nExcludePlanes << endl;

  // Initialise Mille statistics
  _nMilleDataPoints = 0;
  _nMilleTracks = 0;

  // check if Marlin was built with GEAR support or not
#ifndef USE_GEAR

  streamlog_out ( ERROR2 ) << "Marlin was not built with GEAR support." << endl;
  streamlog_out ( ERROR2 ) << "You need to install GEAR and recompile Marlin with -DUSE_GEAR before continue." << endl;

  // I'm thinking if this is the case of throwing an exception or
  // not. This is a really error and not something that can
  // exceptionally happens. Still not sure what to do
  exit(-1);

#else

  // check if the GEAR manager pointer is not null!
  if ( Global::GEAR == 0x0 ) {
    streamlog_out ( ERROR2) << "The GearMgr is not available, for an unknown reason." << endl;
    exit(-1);
  }

  _siPlanesParameters  = const_cast<gear::SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

  _histogramSwitch = true;

#endif

  _nPlanes = _siPlanesParameters->getSiPlanesNumber();

  _waferResidX = new double[_nPlanes];
  _waferResidY = new double[_nPlanes];
  _xFitPos = new double[_nPlanes];
  _yFitPos = new double[_nPlanes];

  _telescopeResolX = new double[_nPlanes];
  _telescopeResolY = new double[_nPlanes];

  // booking histograms
  bookHistos();

  streamlog_out ( MESSAGE2 ) << "Initialising Mille..." << endl;
  _mille = new Mille(_binaryFilename.c_str());
  streamlog_out ( MESSAGE2 ) << "The filename for the binary file is: " << _binaryFilename.c_str() << endl;

}

void EUTelMille::processRunHeader (LCRunHeader * rdr) {

  auto_ptr<EUTelRunHeaderImpl> header ( new EUTelRunHeaderImpl (rdr) );
  header->addProcessor( type() ) ;

  // the run header contains the number of detectors. This number
  // should be in principle the same as the number of layers in the
  // geometry description
  if ( header->getNoOfDetector() != _siPlanesParameters->getSiPlanesNumber() ) {
    streamlog_out ( ERROR2 ) << "Error during the geometry consistency check: " << endl;
    streamlog_out ( ERROR2 ) << "The run header says there are " << header->getNoOfDetector() << " silicon detectors " << endl;
    streamlog_out ( ERROR2 ) << "The GEAR description says     " << _siPlanesParameters->getSiPlanesNumber()
                             << " silicon planes" << endl;
    string answer;
    while (true) {
      streamlog_out ( ERROR2 ) << "Type Q to quit now or C to continue using the actual GEAR description anyway [Q/C]" << endl;
      cin >> answer;
      // put the answer in lower case before making the comparison.
      transform( answer.begin(), answer.end(), answer.begin(), ::tolower );
      if ( answer == "q" ) {
        exit(-1);
      } else if ( answer == "c" ) {
        break;
      }
    }
  }

  // this is the right place also to check the geometry ID. This is a
  // unique number identifying each different geometry used at the
  // beam test. The same number should be saved in the run header and
  // in the xml file. If the numbers are different, instead of barely
  // quitting ask the user what to do.

  if ( header->getGeoID() != _siPlanesParameters->getSiPlanesID() ) {
    streamlog_out ( ERROR2 ) << "Error during the geometry consistency check: " << endl;
    streamlog_out ( ERROR2 ) << "The run header says the GeoID is " << header->getGeoID() << endl;
    streamlog_out ( ERROR2 ) << "The GEAR description says is     " << _siPlanesParameters->getSiPlanesNumber() << endl;
    string answer;
    while (true) {
      streamlog_out ( ERROR2 ) << "Type Q to quit now or C to continue using the actual GEAR description anyway [Q/C]" << endl;
      cin >> answer;
      // put the answer in lower case before making the comparison.
      transform( answer.begin(), answer.end(), answer.begin(), ::tolower );
      if ( answer == "q" ) {
        exit(-1);
      } else if ( answer == "c" ) {
        break;
      }
    }
  }



  // increment the run counter
  ++_iRun;

}

void EUTelMille::FitTrack(int nPlanesFitter, double xPosFitter[], double yPosFitter[], double zPosFitter[], double xResFitter[], double yResFitter[], double chi2Fit[2], double residXFit[], double residYFit[], double angleFit[2]) {

  int sizearray;

  if (_nExcludePlanes > 0) {
    sizearray = nPlanesFitter - _nExcludePlanes;
  } else {
    sizearray = nPlanesFitter;
  }

  double * xPosFit = new double[sizearray];
  double * yPosFit = new double[sizearray];
  double * zPosFit = new double[sizearray];
  double * xResFit = new double[sizearray];
  double * yResFit = new double[sizearray];

  int nPlanesFit = 0;

  for (int help = 0; help < nPlanesFitter; help++) {

    int excluded = 0;

    // check if actual plane is excluded
    if (_nExcludePlanes > 0) {
      for (int helphelp = 0; helphelp < _nExcludePlanes; helphelp++) {
        if (help == _excludePlanes[helphelp]) {
          excluded = 1;
        }
      }
    }

    if (excluded == 1) {
      // do noting
    } else {
      xPosFit[nPlanesFit] = xPosFitter[help];
      yPosFit[nPlanesFit] = yPosFitter[help];
      zPosFit[nPlanesFit] = zPosFitter[help];
      xResFit[nPlanesFit] = xResFitter[help];
      yResFit[nPlanesFit] = yResFitter[help];
      nPlanesFit++;
    }
  }

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++ See Blobel Page 226 !!! +++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++

  int counter;

  float S1[2]   = {0,0};
  float Sx[2]   = {0,0};
  float Xbar[2] = {0,0};

  float * Zbar_X = new float[nPlanesFit];
  float * Zbar_Y = new float[nPlanesFit];
  for (counter = 0; counter < nPlanesFit; counter++){
    Zbar_X[counter] = 0.;
    Zbar_Y[counter] = 0.;
  }

  float Sy[2]     = {0,0};
  float Ybar[2]   = {0,0};
  float Sxybar[2] = {0,0};
  float Sxxbar[2] = {0,0};
  float A2[2]     = {0,0};

  // define S1
  for( counter = 0; counter < nPlanesFit; counter++ ){
    S1[0] = S1[0] + 1/pow(xResFit[counter],2);
    S1[1] = S1[1] + 1/pow(yResFit[counter],2);
  }

  // define Sx
  for( counter = 0; counter < nPlanesFit; counter++ ){
    Sx[0] = Sx[0] + zPosFit[counter]/pow(xResFit[counter],2);
    Sx[1] = Sx[1] + zPosFit[counter]/pow(yResFit[counter],2);
  }

  // define Xbar
  Xbar[0]=Sx[0]/S1[0];
  Xbar[1]=Sx[1]/S1[1];

  // coordinate transformation !! -> bar
  for( counter = 0; counter < nPlanesFit; counter++ ){
    Zbar_X[counter] = zPosFit[counter]-Xbar[0];
    Zbar_Y[counter] = zPosFit[counter]-Xbar[1];
  }

  // define Sy
  for( counter = 0; counter < nPlanesFit; counter++ ){
    Sy[0] = Sy[0] + xPosFit[counter]/pow(xResFit[counter],2);
    Sy[1] = Sy[1] + yPosFit[counter]/pow(yResFit[counter],2);
  }

  // define Ybar
  Ybar[0]=Sy[0]/S1[0];
  Ybar[1]=Sy[1]/S1[1];

  // define Sxybar
  for( counter = 0; counter < nPlanesFit; counter++ ){
    Sxybar[0] = Sxybar[0] + Zbar_X[counter] * xPosFit[counter]/pow(xResFit[counter],2);
    Sxybar[1] = Sxybar[1] + Zbar_Y[counter] * yPosFit[counter]/pow(yResFit[counter],2);
  }

  // define Sxxbar
  for( counter = 0; counter < nPlanesFit; counter++ ){
    Sxxbar[0] = Sxxbar[0] + Zbar_X[counter] * Zbar_X[counter]/pow(xResFit[counter],2);
    Sxxbar[1] = Sxxbar[1] + Zbar_Y[counter] * Zbar_Y[counter]/pow(yResFit[counter],2);
  }

  // define A2
  A2[0]=Sxybar[0]/Sxxbar[0];
  A2[1]=Sxybar[1]/Sxxbar[1];

  // Calculate chi sqaured
  // Chi^2 for X and Y coordinate for hits in all planes
  for( counter = 0; counter < nPlanesFit; counter++ ){
    chi2Fit[0] += pow(-zPosFit[counter]*A2[0]
                      +xPosFit[counter]-Ybar[0]+Xbar[0]*A2[0],2)/pow(xResFit[counter],2);
    chi2Fit[1] += pow(-zPosFit[counter]*A2[1]
                      +yPosFit[counter]-Ybar[1]+Xbar[1]*A2[1],2)/pow(yResFit[counter],2);
  }

  for( counter = 0; counter < nPlanesFitter; counter++ ) {

    residXFit[counter] = (Ybar[0]-Xbar[0]*A2[0]+zPosFitter[counter]*A2[0])-xPosFitter[counter];
    residYFit[counter] = (Ybar[1]-Xbar[1]*A2[1]+zPosFitter[counter]*A2[1])-yPosFitter[counter];

    // residXFit[counter] = xPosFitter[counter] - (Ybar[0]-Xbar[0]*A2[0]+zPosFitter[counter]*A2[0]);
    // residYFit[counter] = yPosFitter[counter] - (Ybar[1]-Xbar[1]*A2[1]+zPosFitter[counter]*A2[1]);

  }

  // define angle
  angleFit[0] = atan(A2[0]);
  angleFit[1] = atan(A2[1]);

  // clean up
  delete [] zPosFit;
  delete [] yPosFit;
  delete [] xPosFit;
  delete [] yResFit;
  delete [] xResFit;

  delete [] Zbar_X;
  delete [] Zbar_Y;

}

void EUTelMille::processEvent (LCEvent * event) {

  if (_iEvt % 10 == 0) {
    streamlog_out( MESSAGE2 ) << "Processing event "
                              << setw(6) << setiosflags(ios::right) << event->getEventNumber() << " in run "
                              << setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber() << setfill(' ')
                              << " (Total = " << setw(10) << _iEvt << ")" << resetiosflags(ios::left) << endl;
    streamlog_out( MESSAGE2 ) << "Currently having " << _nMilleDataPoints << " data points in " 
			      << _nMilleTracks << " tracks " << endl;
  }


  // fill resolution arrays
  for (int help = 0; help < _nPlanes; help++) {
    _telescopeResolX[help] = _telescopeResolution;
    _telescopeResolY[help] = _telescopeResolution;
  }

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event) ;

  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG2 ) << "EORE found: nothing else to do." << endl;
    return;
  }

  LCCollection* collection;
  try {
    if (_inputMode == 1) {
      collection = event->getCollection(_trackCollectionName);
    } else {
      collection = event->getCollection(_hitCollectionName);
    }
  } catch (DataNotAvailableException& e) {
    streamlog_out ( WARNING2 ) << "No input collection found for event " << event->getEventNumber()
                               << " in run " << event->getRunNumber() << endl;
    throw SkipEventException(this);
  }

  int detectorID    = -99; // it's a non sense
  int oldDetectorID = -100;
  int layerIndex;

  vector<EUTelMille::HitsInPlane > _hitsFirstPlane;
  vector<EUTelMille::HitsInPlane > _hitsSecondPlane;
  vector<EUTelMille::HitsInPlane > _hitsThirdPlane;
  vector<EUTelMille::HitsInPlane > _hitsFourthPlane;
  vector<EUTelMille::HitsInPlane > _hitsFifthPlane;
  vector<EUTelMille::HitsInPlane > _hitsSixthPlane;

  HitsInPlane hitsInPlane;

  // check if running in input mode 0 or 2
  if (_inputMode == 0) {

    // loop over all hits in collection
    for ( int iHit = 0; iHit < collection->getNumberOfElements(); iHit++ ) {

      TrackerHitImpl * hit = static_cast<TrackerHitImpl*> ( collection->getElementAt(iHit) );

      LCObjectVec clusterVector = hit->getRawHits();

      EUTelVirtualCluster * cluster;
      if ( hit->getType() == kEUTelFFClusterImpl ) {

        // fixed cluster implementation. Remember it can come from
        // both RAW and ZS data
        cluster = new EUTelFFClusterImpl( static_cast<TrackerDataImpl *> ( clusterVector[0] ) );
      } else if ( hit->getType() == kEUTelSparseClusterImpl ) {

        // ok the cluster is of sparse type, but we also need to know
        // the kind of pixel description used. This information is
        // stored in the corresponding original data collection.

        LCCollectionVec * sparseClusterCollectionVec = dynamic_cast < LCCollectionVec * > (evt->getCollection("original_zsdata"));
        TrackerDataImpl * oneCluster = dynamic_cast<TrackerDataImpl*> (sparseClusterCollectionVec->getElementAt( 0 ));
        CellIDDecoder<TrackerDataImpl > anotherDecoder(sparseClusterCollectionVec);
        SparsePixelType pixelType = static_cast<SparsePixelType> ( static_cast<int> ( anotherDecoder( oneCluster )["sparsePixelType"] ));

        // now we know the pixel type. So we can properly create a new
        // instance of the sparse cluster
        if ( pixelType == kEUTelSimpleSparsePixel ) {
          cluster = new EUTelSparseClusterImpl< EUTelSimpleSparsePixel >
            ( static_cast<TrackerDataImpl *> ( clusterVector[ 0 ]  ) );
        } else {
          streamlog_out ( ERROR4 ) << "Unknown pixel type.  Sorry for quitting." << endl;
          throw UnknownDataTypeException("Pixel type unknown");
        }

      } else {
        throw UnknownDataTypeException("Unknown cluster type");
      }

      detectorID = cluster->getDetectorID();

      if ( detectorID != oldDetectorID ) {
        oldDetectorID = detectorID;

        if ( _conversionIdMap.size() != (unsigned) _siPlanesParameters->getSiPlanesNumber() ) {
          // first of all try to see if this detectorID already belong to
          if ( _conversionIdMap.find( detectorID ) == _conversionIdMap.end() ) {
            // this means that this detector ID was not already inserted,
            // so this is the right place to do that
            for ( int iLayer = 0; iLayer < _siPlanesLayerLayout->getNLayers(); iLayer++ ) {
              if ( _siPlanesLayerLayout->getID(iLayer) == detectorID ) {
                _conversionIdMap.insert( make_pair( detectorID, iLayer ) );
                break;
              }
            }
          }
        }

      }

      layerIndex   = _conversionIdMap[detectorID];

      // Getting positions of the hits.
      // ------------------------------

      hitsInPlane.measuredX = 1000 * hit->getPosition()[0];
      hitsInPlane.measuredY = 1000 * hit->getPosition()[1];
      hitsInPlane.measuredZ = 1000 * hit->getPosition()[2];

      delete cluster; // <--- destroying the cluster

      // Add Hits to vector

      if (layerIndex == 0) {
        _hitsFirstPlane.push_back(hitsInPlane);
      } else if (layerIndex == 1) {
        _hitsSecondPlane.push_back(hitsInPlane);
      } else if (layerIndex == 2) {
        _hitsThirdPlane.push_back(hitsInPlane);
      } else if (layerIndex == 3) {
        _hitsFourthPlane.push_back(hitsInPlane);
      } else if (layerIndex == 4) {
        _hitsFifthPlane.push_back(hitsInPlane);
      } else if (layerIndex == 5) {
        _hitsSixthPlane.push_back(hitsInPlane);
      }

    } // end loop over all hits in collection

  } else if (_inputMode == 2) {

#if defined( USE_ROOT ) || defined(MARLIN_USE_ROOT)

    const float resolX = _testModeSensorResolution;
    const float resolY = _testModeSensorResolution;

    const float xhitpos = gRandom->Uniform(-3500.0,3500.0);
    const float yhitpos = gRandom->Uniform(-3500.0,3500.0);

    const float xslope = gRandom->Gaus(0.0,_testModeXTrackSlope);
    const float yslope = gRandom->Gaus(0.0,_testModeYTrackSlope);

    // loop over all planes
    for (int help = 0; help < _nPlanes; help++) {

      // The x and y positions are given by the sums of the measured
      // hit positions, the detector resolution, the shifts of the
      // planes and the effect due to the track slopes.

      if (help == 0) {

        hitsInPlane.measuredX = xhitpos + gRandom->Gaus(0.0,resolX) + _testModeSensorXShifts[0] + _testModeSensorZPositions[0] * tan(xslope) - _testModeSensorGamma[0] * yhitpos - _testModeSensorBeta[0] * _testModeSensorZPositions[0];
        hitsInPlane.measuredY = yhitpos + gRandom->Gaus(0.0,resolY) + _testModeSensorYShifts[0] + _testModeSensorZPositions[0] * tan(yslope) + _testModeSensorGamma[0] * xhitpos - _testModeSensorAlpha[0] * _testModeSensorZPositions[0];
        hitsInPlane.measuredZ = _testModeSensorZPositions[0];
        _hitsFirstPlane.push_back(hitsInPlane);
        _telescopeResolX[help] = resolX;
        _telescopeResolY[help] = resolY;

      } else if (help == 1) {

        hitsInPlane.measuredX = xhitpos + gRandom->Gaus(0.0,resolX) + _testModeSensorXShifts[1] + _testModeSensorZPositions[1] * tan(xslope) - _testModeSensorGamma[1] * yhitpos - _testModeSensorBeta[1] * _testModeSensorZPositions[1];
        hitsInPlane.measuredY = yhitpos + gRandom->Gaus(0.0,resolY) + _testModeSensorYShifts[1] + _testModeSensorZPositions[1] * tan(yslope) + _testModeSensorGamma[1] * xhitpos - _testModeSensorAlpha[1] * _testModeSensorZPositions[1];
        hitsInPlane.measuredZ = _testModeSensorZPositions[1];
        _hitsSecondPlane.push_back(hitsInPlane);
        _telescopeResolX[help] = resolX;
        _telescopeResolY[help] = resolY;

      } else if (help == 2) {

        hitsInPlane.measuredX = xhitpos + gRandom->Gaus(0.0,resolX) + _testModeSensorXShifts[2] + _testModeSensorZPositions[2] * tan(xslope) - _testModeSensorGamma[2] * yhitpos - _testModeSensorBeta[2] * _testModeSensorZPositions[2];
        hitsInPlane.measuredY = yhitpos + gRandom->Gaus(0.0,resolY) + _testModeSensorYShifts[2] + _testModeSensorZPositions[2] * tan(yslope) + _testModeSensorGamma[2] * xhitpos - _testModeSensorAlpha[2] * _testModeSensorZPositions[2];
        hitsInPlane.measuredZ = _testModeSensorZPositions[2];
        _hitsThirdPlane.push_back(hitsInPlane);
        _telescopeResolX[help] = resolX;
        _telescopeResolY[help] = resolY;

      } else if (help == 3) {

        hitsInPlane.measuredX = xhitpos + gRandom->Gaus(0.0,resolX) + _testModeSensorXShifts[3] + _testModeSensorZPositions[3] * tan(xslope) - _testModeSensorGamma[3] * yhitpos - _testModeSensorBeta[3] * _testModeSensorZPositions[3];
        hitsInPlane.measuredY = yhitpos + gRandom->Gaus(0.0,resolY) + _testModeSensorYShifts[3] + _testModeSensorZPositions[3] * tan(yslope) + _testModeSensorGamma[3] * xhitpos - _testModeSensorAlpha[3] * _testModeSensorZPositions[3];
        hitsInPlane.measuredZ = _testModeSensorZPositions[3];
        _hitsFourthPlane.push_back(hitsInPlane);
        _telescopeResolX[help] = resolX;
        _telescopeResolY[help] = resolY;

      } else if (help == 4) {

        hitsInPlane.measuredX = xhitpos + gRandom->Gaus(0.0,resolX) + _testModeSensorXShifts[4] + _testModeSensorZPositions[4] * tan(xslope) - _testModeSensorGamma[4] * yhitpos - _testModeSensorBeta[4] * _testModeSensorZPositions[4];
        hitsInPlane.measuredY = yhitpos + gRandom->Gaus(0.0,resolY) + _testModeSensorYShifts[4] + _testModeSensorZPositions[4] * tan(yslope) + _testModeSensorGamma[4] * xhitpos - _testModeSensorAlpha[4] * _testModeSensorZPositions[4];
        hitsInPlane.measuredZ = _testModeSensorZPositions[4];
        _hitsFifthPlane.push_back(hitsInPlane);
        _telescopeResolX[help] = resolX;
        _telescopeResolY[help] = resolY;

      } else if (help == 5) {

        hitsInPlane.measuredX = xhitpos + gRandom->Gaus(0.0,resolX) + _testModeSensorXShifts[5] + _testModeSensorZPositions[5] * tan(xslope) - _testModeSensorGamma[5] * yhitpos - _testModeSensorBeta[5] * _testModeSensorZPositions[5];
        hitsInPlane.measuredY = yhitpos + gRandom->Gaus(0.0,resolY) + _testModeSensorYShifts[5] + _testModeSensorZPositions[5] * tan(yslope) + _testModeSensorGamma[5] * xhitpos - _testModeSensorAlpha[5] * _testModeSensorZPositions[5];
        hitsInPlane.measuredZ = _testModeSensorZPositions[5];
        _hitsSixthPlane.push_back(hitsInPlane);
        _telescopeResolX[help] = resolX;
        _telescopeResolY[help] = resolY;

      }

    } // end loop over all planes

#else // USE_ROOT

    throw MissingLibraryException( this, "ROOT" );

#endif


  } // end if check running in input mode 0 or 2

  // Distances between hits in individual planes

  double distance12 = 0.0;
  double distance23 = 0.0;
  double distance34 = 0.0;
  double distance45 = 0.0;
  double distance56 = 0.0;

  double distance_plane12 = 0.0;
  double distance_plane23 = 0.0;
  double distance_plane34 = 0.0;
  double distance_plane45 = 0.0;
  double distance_plane56 = 0.0;

  double distanceMax12 = 0.0;
  double distanceMax23 = 0.0;
  double distanceMax34 = 0.0;
  double distanceMax45 = 0.0;
  double distanceMax56 = 0.0;

  _xPos = new double *[_maxTrackCandidates];
  _yPos = new double *[_maxTrackCandidates];
  _zPos = new double *[_maxTrackCandidates];

  for (int help = 0; help < _maxTrackCandidates; help++) {
    _xPos[help] = new double[_nPlanes];
    _yPos[help] = new double[_nPlanes];
    _zPos[help] = new double[_nPlanes];
  }

  int fitplane[6] = {0, 0, 0, 0, 0, 0};

  for (int help = 0; help < _nPlanes; help++) {
    fitplane[help] = 1;
  }

  int _nTracks = 0;

  int _nGoodTracks = 0;

  // check if running in input mode 0 or 2 => perform simple track finding
  if (_inputMode == 0 || _inputMode == 2) {

    // Find track candidates using the distance cuts
    // ---------------------------------------------
    //
    // This is done separately for different numbers of planes.

    // loop over all hits in first plane
    for (int firsthit = 0; uint(firsthit) < _hitsFirstPlane.size(); firsthit++) {

      // loop over all hits in second plane
      for (int secondhit = 0; uint(secondhit) < _hitsSecondPlane.size(); secondhit++) {

        distance12 = sqrt(pow(_hitsFirstPlane[firsthit].measuredX - _hitsSecondPlane[secondhit].measuredX,2) + pow(_hitsFirstPlane[firsthit].measuredY - _hitsSecondPlane[secondhit].measuredY,2));

        distance_plane12 = _hitsSecondPlane[secondhit].measuredZ - _hitsFirstPlane[firsthit].measuredZ;

        distanceMax12 = _distanceMax * (distance_plane12 / 100000.0);

        if (_nPlanes == 2 && distance12 < distanceMax12 && _nTracks < _maxTrackCandidates) {

          if (_onlySingleHitEvents == 0 || (_hitsFirstPlane.size() == 1 && _hitsSecondPlane.size() == 1)) {

            _xPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredX;
            _yPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredY;
            _zPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredZ;

            _xPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredX;
            _yPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredY;
            _zPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredZ;

            _nTracks++;

          }

        }

        // more than two planes
        if (_nPlanes > 2) {

          // loop over all hits in third plane
          for (int thirdhit = 0; uint(thirdhit) < _hitsThirdPlane.size(); thirdhit++) {

            distance23 = sqrt(pow(_hitsSecondPlane[secondhit].measuredX - _hitsThirdPlane[thirdhit].measuredX,2) + pow(_hitsSecondPlane[secondhit].measuredY - _hitsThirdPlane[thirdhit].measuredY,2));

            distance_plane23 = _hitsThirdPlane[thirdhit].measuredZ - _hitsSecondPlane[secondhit].measuredZ;

            distanceMax23 = _distanceMax * (distance_plane23 / 100000.0);

            if (_nPlanes == 3 && distance12 < distanceMax12 && distance23 < distanceMax23 && _nTracks < _maxTrackCandidates) {

              if (_onlySingleHitEvents == 0 || (_hitsFirstPlane.size() == 1 && _hitsSecondPlane.size() == 1 && _hitsThirdPlane.size() == 1)) {

                _xPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredX;
                _yPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredY;
                _zPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredZ;

                _xPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredX;
                _yPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredY;
                _zPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredZ;

                _xPos[_nTracks][2] = _hitsThirdPlane[thirdhit].measuredX;
                _yPos[_nTracks][2] = _hitsThirdPlane[thirdhit].measuredY;
                _zPos[_nTracks][2] = _hitsThirdPlane[thirdhit].measuredZ;

                _nTracks++;

              }

            }

            // more than three planes
            if (_nPlanes > 3) {

              // loop over all hits in fourth plane
              for (int fourthhit = 0; uint(fourthhit) < _hitsFourthPlane.size(); fourthhit++) {

                distance34 = sqrt(pow(_hitsThirdPlane[thirdhit].measuredX - _hitsFourthPlane[fourthhit].measuredX,2) + pow(_hitsThirdPlane[thirdhit].measuredY - _hitsFourthPlane[fourthhit].measuredY,2));

                distance_plane34 = _hitsFourthPlane[fourthhit].measuredZ - _hitsThirdPlane[thirdhit].measuredZ;

                distanceMax34 = _distanceMax * (distance_plane34 / 100000.0);

                if (_nPlanes == 4 && distance12 < distanceMax12 && distance23 < distanceMax23 && distance34 < distanceMax34 && _nTracks < _maxTrackCandidates) {

                  if (_onlySingleHitEvents == 0 || (_hitsFirstPlane.size() == 1 && _hitsSecondPlane.size() == 1 && _hitsThirdPlane.size() == 1 && _hitsFourthPlane.size() == 1)) {

                    _xPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredX;
                    _yPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredY;
                    _zPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredZ;

                    _xPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredX;
                    _yPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredY;
                    _zPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredZ;

                    _xPos[_nTracks][2] = _hitsThirdPlane[thirdhit].measuredX;
                    _yPos[_nTracks][2] = _hitsThirdPlane[thirdhit].measuredY;
                    _zPos[_nTracks][2] = _hitsThirdPlane[thirdhit].measuredZ;

                    _xPos[_nTracks][3] = _hitsFourthPlane[fourthhit].measuredX;
                    _yPos[_nTracks][3] = _hitsFourthPlane[fourthhit].measuredY;
                    _zPos[_nTracks][3] = _hitsFourthPlane[fourthhit].measuredZ;

                    _nTracks++;

                  }

                }

                // more than four planes
                if (_nPlanes > 4) {

                  // loop over all hits in fifth plane
                  for (int fifthhit = 0; uint(fifthhit) < _hitsFifthPlane.size(); fifthhit++) {

                    distance45 = sqrt(pow(_hitsFourthPlane[fourthhit].measuredX - _hitsFifthPlane[fifthhit].measuredX,2) + pow(_hitsFourthPlane[fourthhit].measuredY - _hitsFifthPlane[fifthhit].measuredY,2));

                    distance_plane45 = _hitsFifthPlane[fifthhit].measuredZ - _hitsFourthPlane[fourthhit].measuredZ;

                    distanceMax45 = _distanceMax * (distance_plane45 / 100000.0);

                    if (_nPlanes == 5 && distance12 < distanceMax12 && distance23 < distanceMax23 && distance34 < distanceMax34 && distance45 < distanceMax45 && _nTracks < _maxTrackCandidates) {

                      if (_onlySingleHitEvents == 0 || (_hitsFirstPlane.size() == 1 && _hitsSecondPlane.size() == 1 && _hitsThirdPlane.size() == 1 && _hitsFourthPlane.size() == 1 && _hitsFifthPlane.size() == 1)) {

                        _xPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredX;
                        _yPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredY;
                        _zPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredZ;

                        _xPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredX;
                        _yPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredY;
                        _zPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredZ;

                        _xPos[_nTracks][2] = _hitsThirdPlane[thirdhit].measuredX;
                        _yPos[_nTracks][2] = _hitsThirdPlane[thirdhit].measuredY;
                        _zPos[_nTracks][2] = _hitsThirdPlane[thirdhit].measuredZ;

                        _xPos[_nTracks][3] = _hitsFourthPlane[fourthhit].measuredX;
                        _yPos[_nTracks][3] = _hitsFourthPlane[fourthhit].measuredY;
                        _zPos[_nTracks][3] = _hitsFourthPlane[fourthhit].measuredZ;

                        _xPos[_nTracks][4] = _hitsFifthPlane[fifthhit].measuredX;
                        _yPos[_nTracks][4] = _hitsFifthPlane[fifthhit].measuredY;
                        _zPos[_nTracks][4] = _hitsFifthPlane[fifthhit].measuredZ;

                        _nTracks++;

                      }

                    }

                    // more than five planes
                    if (_nPlanes > 5) {

                      // loop over all hits in sixth plane
                      for (int sixthhit = 0; uint(sixthhit) < _hitsSixthPlane.size(); sixthhit++) {

                        distance56 = sqrt(pow(_hitsFifthPlane[fifthhit].measuredX - _hitsSixthPlane[sixthhit].measuredX,2) + pow(_hitsFifthPlane[fifthhit].measuredY - _hitsSixthPlane[sixthhit].measuredY,2));

                        distance_plane56 = _hitsSixthPlane[sixthhit].measuredZ - _hitsFifthPlane[fifthhit].measuredZ;

                        distanceMax56 = _distanceMax * (distance_plane56 / 100000.0);

                        if (_nPlanes == 6 && distance12 < distanceMax12 && distance23 < distanceMax23 && distance34 < distanceMax34 && distance45 < distanceMax45 && distance56 < distanceMax56 && _nTracks < _maxTrackCandidates) {

                          if (_onlySingleHitEvents == 0 || (_hitsFirstPlane.size() == 1 && _hitsSecondPlane.size() == 1 && _hitsThirdPlane.size() == 1 && _hitsFourthPlane.size() == 1 && _hitsFifthPlane.size() == 1 && _hitsSixthPlane.size() == 1)) {

                            _xPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredX;
                            _yPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredY;
                            _zPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredZ;

                            _xPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredX;
                            _yPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredY;
                            _zPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredZ;

                            _xPos[_nTracks][2] = _hitsThirdPlane[thirdhit].measuredX;
                            _yPos[_nTracks][2] = _hitsThirdPlane[thirdhit].measuredY;
                            _zPos[_nTracks][2] = _hitsThirdPlane[thirdhit].measuredZ;

                            _xPos[_nTracks][3] = _hitsFourthPlane[fourthhit].measuredX;
                            _yPos[_nTracks][3] = _hitsFourthPlane[fourthhit].measuredY;
                            _zPos[_nTracks][3] = _hitsFourthPlane[fourthhit].measuredZ;

                            _xPos[_nTracks][4] = _hitsFifthPlane[fifthhit].measuredX;
                            _yPos[_nTracks][4] = _hitsFifthPlane[fifthhit].measuredY;
                            _zPos[_nTracks][4] = _hitsFifthPlane[fifthhit].measuredZ;

                            _xPos[_nTracks][5] = _hitsSixthPlane[sixthhit].measuredX;
                            _yPos[_nTracks][5] = _hitsSixthPlane[sixthhit].measuredY;
                            _zPos[_nTracks][5] = _hitsSixthPlane[sixthhit].measuredZ;

                            _nTracks++;

                          }

                        }

                      } // end loop over all hits in sixth plane

                    } // end if more than five planes

                  } // end loop over all hits in fifth plane

                } // end if more than four planes

              } // end loop over all hits in fourth plane

            } // end if more than three planes

          } // end loop over all hits in third plane

        } // end if more than two planes

      } // end loop over all hits in second plane

    } // end loop over all hits in first plane

    // end check if running in input mode 0 or 2 => perform simple track finding
  } else if (_inputMode == 1) {

    const int nTracksHere = collection->getNumberOfElements();

    streamlog_out ( MILLEMESSAGE ) << "Number of tracks available in track collection: " << nTracksHere << endl;

    // loop over all tracks
    for (int nTracksEvent = 0; nTracksEvent < nTracksHere && nTracksEvent < _maxTrackCandidates; nTracksEvent++) {

      Track *TrackHere = dynamic_cast<Track*> (collection->getElementAt(nTracksEvent));

      // hit list assigned to track

      std::vector<EVENT::TrackerHit*> TrackHitsHere = TrackHere->getTrackerHits();

      // check for a hit in every plane
      if (_nPlanes == int(TrackHitsHere.size() / 2)) {

        // assume hits are ordered in z! start counting from 0
        int nPlaneHere = 0;

        // loop over all hits and fill arrays
        for (int nHits = 0; nHits < int(TrackHitsHere.size()); nHits++) {

          TrackerHit *HitHere = TrackHitsHere.at(nHits);

          // hit positions
          const double *PositionsHere = HitHere->getPosition();

          // assume fitted hits have type 32
          if ( HitHere->getType() == 32 ) {

            // fill hits to arrays
            _xPos[nTracksEvent][nPlaneHere] = PositionsHere[0] * 1000;
            _yPos[nTracksEvent][nPlaneHere] = PositionsHere[1] * 1000;
            _zPos[nTracksEvent][nPlaneHere] = PositionsHere[2] * 1000;

            nPlaneHere++;

          } // end assume fitted hits have type 32

        } // end loop over all hits and fill arrays

        _nTracks++;

      } else {

        streamlog_out ( MILLEMESSAGE ) << "Dropping track " << nTracksEvent << " because there is not a hit in every plane assigned to it." << endl;

      }

    } // end loop over all tracks

  }

  if (_nTracks == _maxTrackCandidates) {
    streamlog_out ( WARNING2 ) << "Maximum number of track candidates reached. Maybe further tracks were skipped" << endl;
  }

  streamlog_out ( MILLEMESSAGE ) << "Number of hits in the individual planes: " << _hitsFirstPlane.size() << " " << _hitsSecondPlane.size() << " " << _hitsThirdPlane.size() << " " << _hitsFourthPlane.size() << " " << _hitsFifthPlane.size() << " " << _hitsSixthPlane.size() << endl;

  streamlog_out ( MILLEMESSAGE ) << "Number of track candidates found: " << _iEvt << ": " << _nTracks << endl;

  // Perform fit for all found track candidates
  // ------------------------------------------

  // only one track or no single track event
  if (_nTracks == 1 || _onlySingleTrackEvents == 0) {

    double Chiquare[2] = {0,0};
    double angle[2] = {0,0};

    // loop over all track candidates
    for (int track = 0; track < _nTracks; track++) {

      _xPosHere = new double[_nPlanes];
      _yPosHere = new double[_nPlanes];
      _zPosHere = new double[_nPlanes];

      for (int help = 0; help < _nPlanes; help++) {
        _xPosHere[help] = _xPos[track][help];
        _yPosHere[help] = _yPos[track][help];
        _zPosHere[help] = _zPos[track][help];
      }

      Chiquare[0] = 0.0;
      Chiquare[1] = 0.0;

      streamlog_out ( MILLEMESSAGE ) << "Adding track using the following coordinates: ";

      // loop over all planes
      for (int help = 0; help < _nPlanes; help++) {

        int excluded = 0;

        // check if actual plane is excluded
        if (_nExcludePlanes > 0) {
          for (int helphelp = 0; helphelp < _nExcludePlanes; helphelp++) {
            if (help == _excludePlanes[helphelp]) {
              excluded = 1;
            }
          }
        }

        if (excluded == 0) {
          streamlog_out ( MILLEMESSAGE ) << _xPosHere[help] << " " << _yPosHere[help] << " " << _zPosHere[help] << "   ";
        }

      } // end loop over all planes

      streamlog_out ( MILLEMESSAGE ) << endl;

      // Calculate residuals
      FitTrack(int(_nPlanes), _xPosHere, _yPosHere, _zPosHere, _telescopeResolX, _telescopeResolY, Chiquare, _waferResidX, _waferResidY, angle);

      streamlog_out ( MILLEMESSAGE ) << "Residuals X: ";

      for (int help = 0; help < _nPlanes; help++) {
        streamlog_out ( MILLEMESSAGE ) << _waferResidX[help] << " ";
      }

      streamlog_out ( MILLEMESSAGE ) << endl;

      streamlog_out ( MILLEMESSAGE ) << "Residuals Y: ";

      for (int help = 0; help < _nPlanes; help++) {
        streamlog_out ( MILLEMESSAGE ) << _waferResidY[help] << " ";
      }

      streamlog_out ( MILLEMESSAGE ) << endl;

      int residualsXOkay = 1;
      int residualsYOkay = 1;

      // check if residal cuts are used
      if (_useResidualCuts != 0) {

        // loop over all sensors
        for (int help = 0; help < _nPlanes; help++) {

          if (_waferResidX[help] < _residualsXMin[help] || _waferResidX[help] > _residualsXMax[help]) {
            residualsXOkay = 0;
          }
          if (_waferResidY[help] < _residualsYMin[help] || _waferResidY[help] > _residualsYMax[help]) {
            residualsYOkay = 0;
          }

        } // end loop over all sensors

      } // end check if residual cuts are used

      if (_useResidualCuts != 0 && (residualsXOkay == 0 || residualsYOkay == 0)) {
        streamlog_out ( MILLEMESSAGE ) << "Track did not pass the residual cuts." << endl;
      }

      // apply track cuts (at the moment only residuals)
      if (_useResidualCuts == 0 || (residualsXOkay == 1 && residualsYOkay == 1)) {

        // Add track to Millepede
        // ---------------------------

        // Easy case: consider only shifts
        if (_alignMode == 2) {

          const int nLC = 4; // number of local parameters
          const int nGL = (_nPlanes - _nExcludePlanes) * 2; // number of global parameters

          float sigma = _telescopeResolution;

          float *derLC = new float[nLC]; // array of derivatives for local parameters
          float *derGL = new float[nGL]; // array of derivatives for global parameters

          int *label = new int[nGL]; // array of labels

          float residual;

          // create labels
          for (int help = 0; help < nGL; help++) {
            label[help] = help + 1;
          }

          for (int help = 0; help < nGL; help++) {
            derGL[help] = 0;
          }

          for (int help = 0; help < nLC; help++) {
            derLC[help] = 0;
          }

          int nExcluded = 0;

          // loop over all planes
          for (int help = 0; help < _nPlanes; help++) {

            int excluded = 0;

            // check if actual plane is excluded
            if (_nExcludePlanes > 0) {
              for (int helphelp = 0; helphelp < _nExcludePlanes; helphelp++) {
                if (help == _excludePlanes[helphelp]) {
                  excluded = 1;
                  nExcluded++;
                }
              }
            }

            // if plane is not excluded
            if (excluded == 0) {

              int helphelp = help - nExcluded; // index of plane after
                                               // excluded planes have
                                               // been removed

              derGL[(helphelp * 2)] = -1;
              derLC[0] = 1;
              derLC[2] = _zPosHere[help];
              residual = _waferResidX[help];

              _mille->mille(nLC,derLC,nGL,derGL,label,residual,sigma);

              derGL[(helphelp * 2)] = 0;
              derLC[0] = 0;
              derLC[2] = 0;

              derGL[((helphelp * 2) + 1)] = -1;
              derLC[1] = 1;
              derLC[3] = _zPosHere[help];
              residual = _waferResidY[help];

              _mille->mille(nLC,derLC,nGL,derGL,label,residual,sigma);

              derGL[((helphelp * 2) + 1)] = 0;
              derLC[1] = 0;
              derLC[3] = 0;

              _nMilleDataPoints++;

            } // end if plane is not excluded

          } // end loop over all planes

          // clean up

          delete [] derLC;
          delete [] derGL;
          delete [] label;

          // Slightly more complicated: add rotation around the z axis
        } else if (_alignMode == 1) {

          const int nLC = 4; // number of local parameters
          const int nGL = _nPlanes * 3; // number of global parameters

          float sigma = _telescopeResolution;

          float *derLC = new float[nLC]; // array of derivatives for local parameters
          float *derGL = new float[nGL]; // array of derivatives for global parameters

          int *label = new int[nGL]; // array of labels

          float residual;

          // create labels
          for (int help = 0; help < nGL; help++) {
            label[help] = help + 1;
          }

          for (int help = 0; help < nGL; help++) {
            derGL[help] = 0;
          }

          for (int help = 0; help < nLC; help++) {
            derLC[help] = 0;
          }

          int nExcluded = 0;

          // loop over all planes
          for (int help = 0; help < _nPlanes; help++) {

            int excluded = 0;

            // check if actual plane is excluded
            if (_nExcludePlanes > 0) {
              for (int helphelp = 0; helphelp < _nExcludePlanes; helphelp++) {
                if (help == _excludePlanes[helphelp]) {
                  excluded = 1;
                  nExcluded++;
                }
              }
            }

            // if plane is not excluded
            if (excluded == 0) {

              int helphelp = help - nExcluded; // index of plane after
                                               // excluded planes have
                                               // been removed

              derGL[(helphelp * 3)] = -1;
              derGL[((helphelp * 3) + 2)] = _yPosHere[help];
              derLC[0] = 1;
              derLC[2] = _zPosHere[help];
              residual = _waferResidX[help];

              _mille->mille(nLC,derLC,nGL,derGL,label,residual,sigma);

              derGL[(helphelp * 3)] = 0;
              derGL[((helphelp * 3) + 2)] = 0;
              derLC[0] = 0;
              derLC[2] = 0;

              derGL[((helphelp * 3) + 1)] = -1;
              derGL[((helphelp * 3) + 2)] = -1 * _xPosHere[help];
              derLC[1] = 1;
              derLC[3] = _zPosHere[help];
              residual = _waferResidY[help];

              _mille->mille(nLC,derLC,nGL,derGL,label,residual,sigma);

              derGL[((helphelp * 3) + 1)] = 0;
              derGL[((helphelp * 3) + 2)] = 0;
              derLC[1] = 0;
              derLC[3] = 0;

              _nMilleDataPoints++;

            } // end if plane is not excluded

          } // end loop over all planes

          // clean up

          delete [] derLC;
          delete [] derGL;
          delete [] label;

        } else {

          streamlog_out ( ERROR2 ) << _alignMode << " is not a valid mode. Please choose 1 or 2." << endl;

        }

        _nGoodTracks++;

        // end local fit
        _mille->end();

        _nMilleTracks++;

        // Fill histograms for individual tracks
        // -------------------------------------

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

        string tempHistoName;

        if ( _histogramSwitch ) {
          {
            stringstream ss;
            ss << _chi2XLocalname << endl;
          }
          if ( AIDA::IHistogram1D* chi2x_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[_chi2XLocalname]) )
            chi2x_histo->fill(Chiquare[0]);
          else {
            streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _chi2XLocalname << endl;
            streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
            _histogramSwitch = false;
          }
        }

        if ( _histogramSwitch ) {
          {
            stringstream ss;
            ss << _chi2YLocalname << endl;
          }
          if ( AIDA::IHistogram1D* chi2y_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[_chi2YLocalname]) )
            chi2y_histo->fill(Chiquare[1]);
          else {
            streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _chi2YLocalname << endl;
            streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
            _histogramSwitch = false;
          }
        }

        // loop over all detector planes
        for( int iDetector = 0; iDetector < _nPlanes; iDetector++ ) {

          if ( _histogramSwitch ) {
            {
              stringstream ss;
              ss << _residualXLocalname << "_d" << iDetector;
              tempHistoName=ss.str();
            }
            if ( AIDA::IHistogram1D* residx_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName.c_str()]) )
              residx_histo->fill(_waferResidX[iDetector]);
            else {
              streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _residualXLocalname << endl;
              streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
              _histogramSwitch = false;
            }
          }

          if ( _histogramSwitch ) {
            {
              stringstream ss;
              ss << _residualYLocalname << "_d" << iDetector;
              tempHistoName=ss.str();
            }
            if ( AIDA::IHistogram1D* residy_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName.c_str()]) )
              residy_histo->fill(_waferResidY[iDetector]);
            else {
              streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _residualYLocalname << endl;
              streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
              _histogramSwitch = false;
            }
          }

        } // end loop over all detector planes

#endif

      } // end if apply track cuts

      // clean up
      delete [] _zPosHere;
      delete [] _yPosHere;
      delete [] _xPosHere;

    } // end loop over all track candidates

  } // end if only one track or no single track event

  streamlog_out ( MILLEMESSAGE ) << "Finished fitting tracks in event " << _iEvt << endl;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  string tempHistoName;

  if ( _histogramSwitch ) {
    {
      stringstream ss;
      ss << _numberTracksLocalname << endl;
    }
    if ( AIDA::IHistogram1D* number_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[_numberTracksLocalname]) )
      number_histo->fill(_nGoodTracks);
    else {
      streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _numberTracksLocalname << endl;
      streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
      _histogramSwitch = false;
    }
  }

#endif

  // clean up
  for (int help = 0; help < _maxTrackCandidates; help++) {
    delete [] _zPos[help];
    delete [] _yPos[help];
    delete [] _xPos[help];
  }

  delete [] _zPos;
  delete [] _yPos;
  delete [] _xPos;

  // count events
  ++_iEvt;

  if ( isFirstEvent() ) _isFirstEvent = false;

}

void EUTelMille::end() {

  delete [] _telescopeResolY;
  delete [] _telescopeResolX;
  delete [] _yFitPos;
  delete [] _xFitPos;
  delete [] _waferResidY;
  delete [] _waferResidX;

  // close the output file
  delete _mille;

  // if write the pede steering file
  if (_generatePedeSteerfile) {

    streamlog_out ( MESSAGE4 ) << endl << "Generating the steering file for the pede program..." << endl;

    string tempHistoName;
    double *meanX = new double[_nPlanes];
    double *meanY = new double[_nPlanes];

    // loop over all detector planes
    for( int iDetector = 0; iDetector < _nPlanes; iDetector++ ) {

      if ( _histogramSwitch ) {
        {
          stringstream ss;
          ss << _residualXLocalname << "_d" << iDetector;
          tempHistoName=ss.str();
        }
        if ( AIDA::IHistogram1D* residx_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName.c_str()]) )
          meanX[iDetector] = residx_histo->mean();
        else {
          streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _residualXLocalname << endl;
          streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
          _histogramSwitch = false;
        }
      }

      if ( _histogramSwitch ) {
        {
          stringstream ss;
          ss << _residualYLocalname << "_d" << iDetector;
          tempHistoName=ss.str();
        }
        if ( AIDA::IHistogram1D* residy_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName.c_str()]) )
          meanY[iDetector] = residy_histo->mean();
        else {
          streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _residualYLocalname << endl;
          streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
          _histogramSwitch = false;
        }
      }

    } // end loop over all detector planes

    ofstream steerFile;
    steerFile.open(_pedeSteerfileName.c_str());

    if (steerFile.is_open()) {

      // find first and last excluded plane
      int firstnotexcl = _nPlanes;
      int lastnotexcl = 0;

      // loop over all planes
      for (int help = 0; help < _nPlanes; help++) {

        int excluded = 0;

        // loop over all excluded planes
        for (int helphelp = 0; helphelp < _nExcludePlanes; helphelp++) {
          if (help == _excludePlanes[helphelp]) {
            excluded = 1;
          }
        } // end loop over all excluded planes

        if (excluded == 0 && firstnotexcl > help) {
          firstnotexcl = help;
        }

        if (excluded == 0 && lastnotexcl < help) {
          lastnotexcl = help;
        }
      } // end loop over all planes

      // calculate average
      double averageX = (meanX[firstnotexcl] + meanX[lastnotexcl]) / 2;
      double averageY = (meanY[firstnotexcl] + meanY[lastnotexcl]) / 2;

      steerFile << "Cfiles" << endl;
      steerFile << _binaryFilename << endl;
      steerFile << endl;

      steerFile << "Parameter" << endl;

      int counter = 0;

      // loop over all planes
      for (int help = 0; help < _nPlanes; help++) {

        int excluded = 0; // flag for excluded planes

        // loop over all excluded planes
        for (int helphelp = 0; helphelp < _nExcludePlanes; helphelp++) {

          if (help == _excludePlanes[helphelp]) {
            excluded = 1;
          }

        } // end loop over all excluded planes

        // if plane not excluded
        if (excluded == 0) {

          // if fixed planes
          if (help == firstnotexcl || help == lastnotexcl) {

            if (_alignMode == 1) {
              steerFile << (counter * 3 + 1) << " 0.0 -1.0" << endl;
              steerFile << (counter * 3 + 2) << " 0.0 -1.0" << endl;
              steerFile << (counter * 3 + 3) << " 0.0 -1.0" << endl;
            } else if (_alignMode == 2) {
              steerFile << (counter * 2 + 1) << " 0.0 -1.0" << endl;
              steerFile << (counter * 2 + 2) << " 0.0 -1.0" << endl;
            }

          } else {

            if (_alignMode == 1) {

              if (_usePedeUserStartValues == 0) {
                steerFile << (counter * 3 + 1) << " " << (averageX - meanX[help]) << " 0.0" << endl;
                steerFile << (counter * 3 + 2) << " " << (averageY - meanY[help]) << " 0.0" << endl;
                steerFile << (counter * 3 + 3) << " 0.0 0.0" << endl;
              } else {
                steerFile << (counter * 3 + 1) << " " << _pedeUserStartValuesX[help] << " 0.0" << endl;
                steerFile << (counter * 3 + 2) << " " << _pedeUserStartValuesY[help] << " 0.0" << endl;
                steerFile << (counter * 3 + 3) << " " << _pedeUserStartValuesGamma[help] << " 0.0" << endl;
              }

            } else if (_alignMode == 2) {

              if (_usePedeUserStartValues == 0) {
                steerFile << (counter * 2 + 1) << " " << (averageX - meanX[help]) << " 0.0" << endl;
                steerFile << (counter * 2 + 2) << " " << (averageY - meanY[help]) << " 0.0" << endl;
              } else {
                steerFile << (counter * 2 + 1) << " " << _pedeUserStartValuesX[help] << " 0.0" << endl;
                steerFile << (counter * 2 + 2) << " " << _pedeUserStartValuesY[help] << " 0.0" << endl;
              }

            }

          }

          counter++;

        } // end if plane not excluded

      } // end loop over all planes

      steerFile << endl;
      steerFile << "! chiscut 5.0 2.5" << endl;
      steerFile << "! outlierdownweighting 4" << endl;
      steerFile << endl;
      steerFile << "method inversion 10 0.001" << endl;
      steerFile << endl;
      steerFile << "histprint" << endl;
      steerFile << endl;
      steerFile << "end" << endl;

      steerFile.close();

      streamlog_out ( MESSAGE2 ) << "File " << _pedeSteerfileName << " written." << endl;

    } else {

      streamlog_out ( ERROR2 ) << "Could not open steering file." << endl;

    }


  } // end if write the pede steering file

  streamlog_out ( MESSAGE2 ) << endl;
  streamlog_out ( MESSAGE2 ) << "Number of data points used: " << _nMilleDataPoints << endl;
  streamlog_out ( MESSAGE2 ) << "Number of tracks used: " << _nMilleTracks << endl;

  // if running pede using the generated steering file
  if (_runPede == 1) {

    // check if steering file exists
    if (_generatePedeSteerfile == 1) {

      std::string command = "pede " + _pedeSteerfileName;

      // before starting pede, let's check if it is in the path
      bool isPedeInPath = true;

      // create a new process
      redi::ipstream which("which pede");

      // wait for the process to finish
      which.close();

      // get the status
      // if it 255 then the program wasn't found in the path
      isPedeInPath = !( which.rdbuf()->status() == 255 );

      if ( !isPedeInPath ) {
        streamlog_out( ERROR ) << "Pede cannot be executed because not found in the path" << endl;
      } else {

        streamlog_out ( MESSAGE2 ) << endl;
        streamlog_out ( MESSAGE2 ) << "Starting pede..." << endl;

        redi::ipstream pede( command.c_str() );
        string output;
        while ( getline( pede, output ) ) {
          streamlog_out( MESSAGE2 ) << output << endl;
        }

        // wait for the pede execution to finish
        pede.close();

        // check the exit value of pede
        if ( pede.rdbuf()->status() == 0 ) {
          streamlog_out ( MESSAGE2 ) << "Pede successfully finished" << endl;
        }

        // reading back the millepede.res file and getting the
        // results.
        string millepedeResFileName = "millepede.res";

        streamlog_out ( MESSAGE2 ) << "Reading back the " << millepedeResFileName << endl
                                   << "Saving the alignment constant into " << _alignmentConstantLCIOFile << endl;

        // open the millepede ASCII output file
        ifstream millepede( millepedeResFileName.c_str() );

        // reopen the LCIO file this time in append mode
        LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();

        try {
          lcWriter->open( _alignmentConstantLCIOFile, LCIO::WRITE_NEW );
        } catch ( IOException& e ) {
          streamlog_out ( ERROR4 ) << e.what() << endl
                                   << "Sorry for quitting. " << endl;
          exit(-1);
        }

        // write an almost empty run header
        LCRunHeaderImpl * lcHeader  = new LCRunHeaderImpl;
        lcHeader->setRunNumber( 0 );


        lcWriter->writeRunHeader(lcHeader);

        delete lcHeader;

        LCEventImpl * event = new LCEventImpl;
        event->setRunNumber( 0 );
        event->setEventNumber( 0 );

        LCTime * now = new LCTime;
        event->setTimeStamp( now->timeStamp() );
        delete now;

        LCCollectionVec * constantsCollection = new LCCollectionVec( LCIO::LCGENERICOBJECT );


        if ( millepede.bad() ) {
          streamlog_out ( ERROR4 ) << "Error opening the " << millepedeResFileName << endl
                                   << "The alignment slcio file cannot be saved" << endl;
        } else {
          vector<double > tokens;
          stringstream tokenizer;
          string line;
          double buffer;

          // get the first line and throw it away since it is a
          // comment!
          getline( millepede, line );

	  // I need to create a map to convert plane position in
	  // sensorID
	  std::vector<int > lookUp;
	  for ( int iPlane = 0; iPlane < _nPlanes ; iPlane++ ) {
	    // if the plane in the exclude list
	    if ( find( _excludePlanes.begin(), _excludePlanes.end(), iPlane ) ==  _excludePlanes.end() ) {
	      lookUp.push_back( _siPlanesLayerLayout->getID(iPlane) );
	    }
	  }

          int counter = 0;

          while ( ! millepede.eof() ) {

            EUTelAlignmentConstant * constant = new EUTelAlignmentConstant;
 
	    bool goodLine = true;


	    constant->setSensorID( lookUp[counter] );
	    ++ counter;
	      
	    for ( unsigned int iParam = 0 ; iParam < 3 ; ++iParam ) {
		
	      getline( millepede, line );
		
	      if ( line.empty() ) {
		goodLine = false;
	      }
		
	      tokens.clear();
	      tokenizer.clear();
	      tokenizer.str( line );
		
	      while ( tokenizer >> buffer ) {
		tokens.push_back( buffer ) ;
	      }
		
	      if ( ( tokens.size() == 3 ) || ( tokens.size() == 6 ) ) {
		goodLine = true;
	      } else goodLine = false;

	      bool isFixed = ( tokens.size() == 3 );
	      if ( isFixed ) {
		streamlog_out ( DEBUG0 ) << "Parameter " << tokens[0] << " is at " << ( tokens[1] / 1000 )
                                         << " (fixed)"  << endl;
	      } else {
                streamlog_out ( DEBUG0 ) << "Parameter " << tokens[0] << " is at " << (tokens[1] / 1000 )
                                         << " +/- " << ( tokens[4] / 1000 )  << endl;
	      }
		
	      if ( iParam == 0 ) {
		constant->setXOffset( tokens[1] / 1000 );
		if ( ! isFixed ) {
		  double err  = tokens[4] / 1000;
		  constant->setXOffsetError( err ) ;
		}
	      }
	      if ( iParam == 1 ) {
		constant->setYOffset( tokens[1] / 1000 ) ;
		if ( ! isFixed ) constant->setYOffsetError( tokens[4] / 1000 ) ;
	      }
	      if ( iParam == 2 ) {
		constant->setGamma( tokens[1]  ) ;
		if ( ! isFixed ) constant->setGammaError( tokens[4] ) ;
	      }
		
	    }
	    



            // right place to add the constant to the collection
            if ( goodLine ) {
              constantsCollection->push_back( constant );
              streamlog_out ( MESSAGE0 ) << (*constant) << endl;
            }
            else delete constant;
          }

        }
        event->addCollection( constantsCollection, "alignment" );
        lcWriter->writeEvent( event );
        delete event;

        lcWriter->close();

        millepede.close();

      }
    } else {

      streamlog_out ( ERROR2 ) << "Unable to run pede. No steering file has been generated." << endl;

    }

  } // end if running pede using the generated steering file

  streamlog_out ( MESSAGE2 ) << endl;
  streamlog_out ( MESSAGE2 ) << "Successfully finished" << endl;

}

void EUTelMille::bookHistos() {


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  try {
    streamlog_out ( MESSAGE2 ) << "Booking histograms..." << endl;

    const int    tracksNBin = 20;
    const double tracksMin = -0.5;
    const double tracksMax = 19.5;
    const int    Chi2NBin = 10000000;
    const double Chi2Min  = 0.;
    const double Chi2Max  = 10000000.;
    const int    NBin = 10000;
    const double Min  = -5000.;
    const double Max  = 5000.;

    AIDA::IHistogram1D * numberTracksLocal =
      AIDAProcessor::histogramFactory(this)->createHistogram1D(_numberTracksLocalname,tracksNBin,tracksMin,tracksMax);
    if ( numberTracksLocal ) {
      numberTracksLocal->setTitle("Number of tracks after #chi^{2} cut");
      _aidaHistoMap.insert( make_pair( _numberTracksLocalname, numberTracksLocal ) );
    } else {
      streamlog_out ( ERROR2 ) << "Problem booking the " << (_numberTracksLocalname) << endl;
      streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
      _histogramSwitch = false;
    }

    AIDA::IHistogram1D * chi2XLocal =
      AIDAProcessor::histogramFactory(this)->createHistogram1D(_chi2XLocalname,Chi2NBin,Chi2Min,Chi2Max);
    if ( chi2XLocal ) {
      chi2XLocal->setTitle("Chi2 X");
      _aidaHistoMap.insert( make_pair( _chi2XLocalname, chi2XLocal ) );
    } else {
      streamlog_out ( ERROR2 ) << "Problem booking the " << (_chi2XLocalname) << endl;
      streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
      _histogramSwitch = false;
    }

    AIDA::IHistogram1D * chi2YLocal =
      AIDAProcessor::histogramFactory(this)->createHistogram1D(_chi2YLocalname,Chi2NBin,Chi2Min,Chi2Max);
    if ( chi2YLocal ) {
      chi2YLocal->setTitle("Chi2 Y");
      _aidaHistoMap.insert( make_pair( _chi2YLocalname, chi2YLocal ) );
    } else {
      streamlog_out ( ERROR2 ) << "Problem booking the " << (_chi2YLocalname) << endl;
      streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
      _histogramSwitch = false;
    }

    string tempHisto;
    string tempHistoName;
    string histoTitleXResid;
    string histoTitleYResid;

    for( int iDetector = 0; iDetector < _nPlanes; iDetector++ ){

      {
        stringstream ss;
        stringstream pp;
        stringstream tt;

        pp << "ResidualXLocal_d" << iDetector;
        tempHisto=pp.str();
        ss << _residualXLocalname << "_d" << iDetector;
        tempHistoName=ss.str();
        tt << "XResidual" << "_d" << iDetector;
        histoTitleXResid=tt.str();

      }

      AIDA::IHistogram1D *  tempXHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram1D(tempHistoName,NBin, Min,Max);
      if ( tempXHisto ) {
        tempXHisto->setTitle(histoTitleXResid);
        _aidaHistoMap.insert( make_pair( tempHistoName, tempXHisto ) );
      } else {
        streamlog_out ( ERROR2 ) << "Problem booking the " << (tempHistoName) << endl;
        streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }

      {
        stringstream ss;
        stringstream pp;
        stringstream tt;

        pp << "ResidualYLocal_d" << iDetector;
        tempHisto=pp.str();
        ss << _residualYLocalname << "_d" << iDetector;
        tempHistoName=ss.str();
        tt << "YResidual" << "_d" << iDetector;
        histoTitleYResid=tt.str();
      }

      AIDA::IHistogram1D *  tempYHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram1D(tempHistoName,NBin, Min,Max);
      if ( tempYHisto ) {
        tempYHisto->setTitle(histoTitleYResid);
        _aidaHistoMap.insert( make_pair( tempHistoName, tempYHisto ) );
      } else {
        streamlog_out ( ERROR2 ) << "Problem booking the " << (tempHistoName) << endl;
        streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }

    }

  } catch (lcio::Exception& e ) {

    streamlog_out ( ERROR2 ) << "No AIDAProcessor initialized. Type q to exit or c to continue without histogramming" << endl;
    string answer;
    while ( true ) {
      streamlog_out ( ERROR2 ) << "[q]/[c]" << endl;
      cin >> answer;
      transform( answer.begin(), answer.end(), answer.begin(), ::tolower );
      if ( answer == "q" ) {
        exit(-1);
      } else if ( answer == "c" )
        _histogramSwitch = false;
      break;
    }

  }
#endif

}

#endif // USE_GEAR


