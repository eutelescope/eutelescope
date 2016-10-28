// Contact: philipp.roloff@desy.de
// Version: $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
 
#ifdef OBSOLETE

// build only if ROOT is used
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)

// built only if GEAR is used
#ifdef USE_GEAR

// eutelescope includes ".h"
#include "EUTelAlign.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"

#include "EUTelSparseClusterImpl.h"
#include "EUTelExceptions.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"

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
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>

// ROOT includes
#include <TMinuit.h>
#include <TSystem.h>

// system includes <>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <memory>

using namespace std;
using namespace marlin;
using namespace gear;
using namespace eutelescope;

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelAlign::_distanceLocalname             = "Distance";
std::string EUTelAlign::_residualXSimpleLocalname      = "ResidualXSimple";
std::string EUTelAlign::_residualYSimpleLocalname      = "ResidualYSimple";
std::string EUTelAlign::_residualXLocalname            = "ResidualX";
std::string EUTelAlign::_residualYLocalname            = "ResidualY";
#endif

void fitFunctionWrapper(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
  EUTelAlign::Chi2Function(npar,gin,f,par,iflag);
}

vector<EUTelAlign::HitsForFit > EUTelAlign::_hitsForFit;

EUTelAlign::EUTelAlign () : Processor("EUTelAlign") {

  // modify processor description
  _description =
    "EUTelAlign makes alignment of 2 planes using predicted positions from straight line fit";

  // input collection

  registerInputCollection(LCIO::TRACKERHIT,"MeasuredHitCollectionName",
                          "Hit collection name",
                          _measHitCollectionName, string ( "hit" ));


  // input parameters

  registerProcessorParameter("AlignedPlane",
                             "Aligned plane",
                             _alignedPlane, static_cast < int > (2));

  registerProcessorParameter("AlignedBox",
                             "Aligned box",
                             _alignedBox, static_cast < int > (1));

  registerOptionalParameter("NumberPlanesFirstBox",
                            "Number of planes in the first telescope box",
                            _nPlanesFirstBox, static_cast < int > (3));

  registerOptionalParameter("ReferencePlane",
                            "Reference plane",
                            _referencePlane, static_cast < int > (1));

  registerOptionalParameter("Resolution",
                            "Resolution of aligned plane",
                            _resolution,static_cast < double > (10.0));


  registerOptionalParameter("Chi2Cut",
                            "Start value for Chi2 cut in fit",
                            _chi2Cut, static_cast < double > (1000.0));

  FloatVec startValues;
  startValues.push_back(0.0);
  startValues.push_back(0.0);
  startValues.push_back(0.0);
  startValues.push_back(0.0);
  startValues.push_back(0.0);
  startValues.push_back(0.0);

  registerOptionalParameter("StartValuesForAlignment","Start values used for the alignment:\n off_x, off_y, theta_x, theta_y, theta_z1, theta_z2",
                            _startValuesForAlignment, startValues);

  registerOptionalParameter("DistanceMin","Minimal allowed distance between hits before alignment.",
                            _distanceMin, static_cast <double> (0.0));

  registerOptionalParameter("DistanceMax","Maximal allowed distance between hits before alignment.",
                            _distanceMax, static_cast <double> (2000.0));

  registerOptionalParameter("NHitsMax","Maximal number of Hits per plane.",
                            _nHitsMax, static_cast <int> (100));

  FloatVec constantsSecondLayer;
  constantsSecondLayer.push_back(0.0);
  constantsSecondLayer.push_back(0.0);
  constantsSecondLayer.push_back(0.0);
  constantsSecondLayer.push_back(0.0);
  constantsSecondLayer.push_back(0.0);

  FloatVec constantsThirdLayer;
  constantsThirdLayer.push_back(0.0);
  constantsThirdLayer.push_back(0.0);
  constantsThirdLayer.push_back(0.0);
  constantsThirdLayer.push_back(0.0);
  constantsThirdLayer.push_back(0.0);

  registerOptionalParameter("AlignmentConstantsSecondLayer","Alignment Constants for second Telescope Layer:\n off_x, off_y, theta_x, theta_y, theta_z",
                            _alignmentConstantsSecondLayer, constantsSecondLayer);
  registerOptionalParameter("AlignmentConstantsThirdLayer","Alignment Constants for third Telescope Layer:\n off_x, off_y, theta_x, theta_y, theta_z",
                            _alignmentConstantsThirdLayer, constantsThirdLayer);

}

void EUTelAlign::init() {
  // this method is called only once even when the rewind is active
  // usually a good idea to

  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

  // check if Marlin was built with GEAR support or not
#ifndef USE_GEAR

  message<ERROR5> ( "Marlin was not built with GEAR support. You need to install GEAR and recompile Marlin with -DUSE_GEAR before continuing.");

  // I'm thinking if this is the case of throwing an exception or
  // not. This is a really error and not something that can
  // exceptionally happens. Still not sure what to do
  exit(-1);

#else

  // check if the GEAR manager pointer is not null!
  if ( Global::GEAR == 0x0 ) {
    message<ERROR5> ( "The GearMgr is not available, for an unknown reason." );
    exit(-1);
  }

  _siPlanesParameters  = const_cast<SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

#endif

  _nPlanes = _siPlanesParameters->getSiPlanesNumber();
  _intrResol = new double[_nPlanes];
  _intrResolX = new double[3];
  _intrResolY = new double[3];
  _xMeasPos = new double[_nPlanes];
  _yMeasPos = new double[_nPlanes];
  _zMeasPos = new double[_nPlanes];

}

void EUTelAlign::processRunHeader (LCRunHeader * rdr) {

  auto header = std::make_unique<EUTelRunHeaderImpl>(rdr);

  // The run header contains the number of detectors. This number
  // should be in principle the same as the number of layers in the
  // geometry description
  if ( header->getNoOfDetector() != _siPlanesParameters->getSiPlanesNumber() ) {
    message<ERROR5> ( "Error during the geometry consistency check: " );
    message<ERROR5> ( log() << "The run header says there are " << header->getNoOfDetector() << " silicon detectors " );
    message<ERROR5> ( log() << "The GEAR description says     " << _siPlanesParameters->getSiPlanesNumber() << " silicon planes" );
    exit(-1);
  }

  // This is the right place also to check the geometry ID. This is a
  // unique number identifying each different geometry used at the
  // beam test. The same number should be saved in the run header and
  // in the xml file. If the numbers are different, instead of just
  // quitting ask the user what to do.

  if ( header->getGeoID() != _siPlanesParameters->getSiPlanesID() ) {
    message<ERROR5> ( "Error during the geometry consistency check: " );
    message<ERROR5> ( log() << "The run header says the GeoID is " << header->getGeoID() );
    message<ERROR5> ( log() << "The GEAR description says is     " << _siPlanesParameters->getSiPlanesNumber() );

#ifdef EUTEL_INTERACTIVE
    string answer;
    while (true) {
      message<ERROR5> ( "Type Q to quit now or C to continue using the actual GEAR description anyway [Q/C]" );
      cin >> answer;
      // Put the answer in lower case before making the comparison.
      transform( answer.begin(), answer.end(), answer.begin(), ::tolower );
      if ( answer == "q" ) {
        exit(-1);
      } else if ( answer == "c" ) {
        break;
      }
    }
#endif

  }

  // Now book histograms
  if ( isFirstEvent() )  bookHistos();

  // Increment the run counter
  ++_iRun;
}

void EUTelAlign::FitTrack(int nPlanesFit, double xPosFit[], double yPosFit[], double zPosFit[], double xResFit[], double yResFit[], double chi2Fit[2], double &_predX, double &_predY, double _predZ) {

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

  _predX = Ybar[0]-Xbar[0]*A2[0]+_predZ*A2[0];
  _predY = Ybar[1]-Xbar[1]*A2[1]+_predZ*A2[1];

  delete [] Zbar_X;
  delete [] Zbar_Y;
}

void EUTelAlign::processEvent (LCEvent * event) {

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event) ;

  if ( evt->getEventType() == kEORE ) {
    message<DEBUG5> ( "EORE found: nothing else to do." );
    return;
  }

  int nHitsFirstPlane = 0;
  int nHitsSecondPlane = 0;

  int nHitsSecondBox = 0;

  try {

    LCCollectionVec* measHitCollection = static_cast<LCCollectionVec*>(event->getCollection( _measHitCollectionName ));
    UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder(measHitCollection);

    int detectorID    = -99; 
    int oldDetectorID = -100;
    int layerIndex;

    HitsForFit hitsForFit;

    double allHitsFirstLayerMeasuredX[20];
    double allHitsFirstLayerMeasuredY[20];
    double allHitsFirstLayerMeasuredZ[20];
    double allHitsSecondLayerMeasuredX[20];
    double allHitsSecondLayerMeasuredY[20];
    double allHitsSecondLayerMeasuredZ[20];
    double allHitsFirstLayerResolution[20];
    double allHitsSecondLayerResolution[20];

    double allHitsSecondBoxX[20];
    double allHitsSecondBoxY[20];
    double allHitsSecondBoxZ[20];

    vector<EUTelAlign::HitsInFirstBox > _hitsFirstPlane;
    vector<EUTelAlign::HitsInFirstBox > _hitsSecondPlane;
    vector<EUTelAlign::HitsInFirstBox > _hitsThirdPlane;

    HitsInFirstBox hitsInFirstBox;

    // Loop over all hits
    for ( int iHit = 0; iHit < measHitCollection->getNumberOfElements(); iHit++ ) {

      TrackerHitImpl* measHit = static_cast<TrackerHitImpl*> ( measHitCollection->getElementAt(iHit) );
      
      detectorID = hitDecoder(measHit)["sensorID"];

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

        // here we take intrinsic resolution from geometry database

        layerIndex   = _conversionIdMap[detectorID];
        _intrResol[layerIndex] = 1000*_siPlanesLayerLayout->getSensitiveResolution(layerIndex); //um

        if (layerIndex < 3) {
          _intrResolX[layerIndex] = 1000*_siPlanesLayerLayout->getSensitiveResolution(layerIndex); //um
          _intrResolY[layerIndex] = 1000*_siPlanesLayerLayout->getSensitiveResolution(layerIndex); //um
        }

      }

      layerIndex = _conversionIdMap[detectorID];

      if (_alignedBox == 2) {

        // The other layers were aligned with respect to the first one.

        double off_x, off_y, theta_x, theta_y, theta_z;

        if (layerIndex == 0) {

          hitsInFirstBox.measuredX = 1000 * measHit->getPosition()[0]; 
          hitsInFirstBox.measuredY = 1000 * measHit->getPosition()[1]; 
          hitsInFirstBox.measuredZ = 1000 * measHit->getPosition()[2];

        } else {

          if (layerIndex == 1) {

            off_x = _alignmentConstantsSecondLayer[0];
            off_y = _alignmentConstantsSecondLayer[1];
            theta_x = _alignmentConstantsSecondLayer[2];
            theta_y = _alignmentConstantsSecondLayer[3];
            theta_z = _alignmentConstantsSecondLayer[4];

          } else if (layerIndex == 2) {

            off_x = _alignmentConstantsThirdLayer[0];
            off_y = _alignmentConstantsThirdLayer[1];
            theta_x = _alignmentConstantsThirdLayer[2];
            theta_y = _alignmentConstantsThirdLayer[3];
            theta_z = _alignmentConstantsThirdLayer[4];

          } else {

            off_x = 0.0;
            off_y = 0.0;
            theta_x = 0.0;
            theta_y = 0.0;
            theta_z = 0.0;

          }

          hitsInFirstBox.measuredX = (cos(theta_y)*cos(theta_z)) * measHit->getPosition()[0] * 1000 + ((-1)*sin(theta_x)*sin(theta_y)*cos(theta_z) + cos(theta_x)*sin(theta_z)) * measHit->getPosition()[1] * 1000 + off_x;
          hitsInFirstBox.measuredY = ((-1)*cos(theta_y)*sin(theta_z)) * measHit->getPosition()[0] * 1000 + (sin(theta_x)*sin(theta_y)*sin(theta_z) + cos(theta_x)*cos(theta_z)) * measHit->getPosition()[1] * 1000 + off_y;
          hitsInFirstBox.measuredZ = 1000 * measHit->getPosition()[2];

        }

        if (layerIndex == 0) {
          _hitsFirstPlane.push_back(hitsInFirstBox);
        } else if (layerIndex == 1) {
          _hitsSecondPlane.push_back(hitsInFirstBox);
        } else if (layerIndex == 2) {
          _hitsThirdPlane.push_back(hitsInFirstBox);
        }

      } // end if _alignedBox == 2

      // If first plane: calculate extrapolation to second plane

      if (layerIndex == (_referencePlane - 1) && nHitsFirstPlane < 20) {

        allHitsFirstLayerMeasuredX[nHitsFirstPlane] = 1000 * measHit->getPosition()[0];
        allHitsFirstLayerMeasuredY[nHitsFirstPlane] = 1000 * measHit->getPosition()[1];
        allHitsFirstLayerMeasuredZ[nHitsFirstPlane] = 1000 * measHit->getPosition()[2];


        allHitsFirstLayerResolution[nHitsFirstPlane] = 1000 * _siPlanesLayerLayout->getSensitiveResolution(layerIndex); // Add multiple scattering later!

        // hitsForFit.secondLayerPredictedX = hitsForFit.firstLayerMeasuredX;
        // hitsForFit.secondLayerPredictedY = hitsForFit.firstLayerMeasuredY;

        nHitsFirstPlane++;

      } else if (layerIndex == (_alignedPlane - 1) && nHitsSecondPlane < 20) {

        allHitsSecondLayerMeasuredX[nHitsSecondPlane] = 1000 * measHit->getPosition()[0];
        allHitsSecondLayerMeasuredY[nHitsSecondPlane] = 1000 * measHit->getPosition()[1];
        allHitsSecondLayerMeasuredZ[nHitsSecondPlane] = 1000 * measHit->getPosition()[2];


        // hitsForFit.secondLayerPredictedZ = 1000 * measHit->getPosition()[2];

        allHitsSecondLayerResolution[nHitsSecondPlane] = 1000 * _siPlanesLayerLayout->getSensitiveResolution(layerIndex); // Add multiple scattering later!

        nHitsSecondPlane++;

      }

      delete cluster; // <--- destroying the cluster

    } // End loop over all hits

    if (_alignedBox == 2 && ((_hitsSecondPlane.size() < 20 && _nPlanesFirstBox == 2) || (_hitsSecondPlane.size() < 20 && _hitsThirdPlane.size() < 20 && _nPlanesFirstBox == 3))) {

      double distance12 = 0.0;
      double distance23 = 0.0;

      _xPos = new double *[500];
      _yPos = new double *[500];
      _zPos = new double *[500];

      for (int help = 0; help < 500; help++) {
        _xPos[help] = new double[3];
        _yPos[help] = new double[3];
        _zPos[help] = new double[3];
      }

      int _nTracks = 0;

      // loop over all hits in first plane
      for (int firsthit = 0; size_t(firsthit) < _hitsFirstPlane.size(); firsthit++) {

        // loop over all hits in second plane
        for (int secondhit = 0; size_t(secondhit) < _hitsSecondPlane.size(); secondhit++) {

          distance12 = sqrt(pow(_hitsFirstPlane[firsthit].measuredX - _hitsSecondPlane[secondhit].measuredX,2) + pow(_hitsFirstPlane[firsthit].measuredY - _hitsSecondPlane[secondhit].measuredY,2));

          if (_nPlanesFirstBox == 2) {

            if (distance12 < 100) {

              _xPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredX;
              _yPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredY;
              _zPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredZ;

              _xPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredX;
              _yPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredY;
              _zPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredZ;

              _nTracks++;

            }

          } else if (_nPlanesFirstBox == 3) {

            // loop over all hits in third plane
            for (int thirdhit = 0; size_t(thirdhit) < _hitsThirdPlane.size(); thirdhit++) {

              distance23 = sqrt(pow(_hitsSecondPlane[secondhit].measuredX - _hitsThirdPlane[thirdhit].measuredX,2) + pow(_hitsSecondPlane[secondhit].measuredY - _hitsThirdPlane[thirdhit].measuredY,2));

              if (distance12 < 100 && distance23 < 100) {

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

            } // end loop over all hits in third plane

          }

        } // end loop over all hits in second plane

      } // end loop over all hits in first plane

      if (nHitsSecondPlane > 0) {

        double Chiquare[2] = {0,0};

        // loop over all track candidates
        for (int track = 0; track < _nTracks; track++) {

          double _predictedX = 0.0;
          double _predictedY = 0.0;
          double _predictedZ = allHitsSecondLayerMeasuredZ[0]; // dangerous !!!

          _xPosHere = new double[_nPlanes];
          _yPosHere = new double[_nPlanes];
          _zPosHere = new double[_nPlanes];

          for (int help = 0; help < _nPlanesFirstBox; help++) {
            _xPosHere[help] = _xPos[track][help];
            _yPosHere[help] = _yPos[track][help];
            _zPosHere[help] = _zPos[track][help];
          }

          streamlog_out ( MESSAGE2 ) << "Fitting track using the following coordinates: ";

          for (int help = 0; help < 3; help++) {
            streamlog_out ( MESSAGE2 ) << _xPosHere[help] << " " << _yPosHere[help] << " " << _zPosHere[help] << "   ";
          }

          streamlog_out ( MESSAGE2 ) << endl;

          FitTrack(_nPlanesFirstBox, _xPosHere, _yPosHere, _zPosHere, _intrResolX, _intrResolY, Chiquare, _predictedX, _predictedY, _predictedZ);

          streamlog_out ( MESSAGE2 ) << "Fit Result: " << _predictedX << " " << _predictedY << " " << _predictedZ << " Chi^2(x): " << Chiquare[0] << " Chi^2(y): " << Chiquare[1] << endl;

          if (Chiquare[0] <= 20.0 && Chiquare[1] <= 20) {
            allHitsSecondBoxX[nHitsSecondBox] = _predictedX;
            allHitsSecondBoxY[nHitsSecondBox] = _predictedY;
            allHitsSecondBoxZ[nHitsSecondBox] = _predictedZ;
            nHitsSecondBox++;
          }

          // clean up
          delete [] _zPosHere;
          delete [] _yPosHere;
          delete [] _xPosHere;


        } // end if loop over all track candidates

      } // end if nHits SecondPlane > 0

      // clean up
      for (int help = 0; help < 500; help++) {
        delete [] _zPos[help];
        delete [] _yPos[help];
        delete [] _xPos[help];
      }

      delete [] _zPos;
      delete [] _yPos;
      delete [] _xPos;

    } // end if _aligendBox == 2

    // check number of hits
    if (nHitsFirstPlane <= _nHitsMax && nHitsSecondPlane <= _nHitsMax) {

      if (_alignedBox == 1) {

        // loop over all hits in first plane
        for (int firsthit = 0; firsthit < nHitsFirstPlane; firsthit++) {

          int take = -1000;
          int veto = -1000;

          // loop over all hits in second plane
          for (int secondhit = 0; secondhit < nHitsSecondPlane; secondhit++) {

            // calculate distance between hits
            double distance = sqrt(pow(allHitsFirstLayerMeasuredX[firsthit] - allHitsSecondLayerMeasuredX[secondhit],2) + pow(allHitsFirstLayerMeasuredY[firsthit] - allHitsSecondLayerMeasuredY[secondhit],2));

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

            if ( AIDA::IHistogram1D* distance_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[_distanceLocalname]) )
              distance_histo->fill(distance);
            else {
              streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " <<  _distanceLocalname << endl;
            }

#endif

            if (distance >= _distanceMin && distance <= _distanceMax) {
              if (take != -1000) {
                veto = 1;
              }
              take = secondhit;
            }

          } // end loop over hits in second plane

          if (take != -1000 && veto == -1000) {

            hitsForFit.firstLayerMeasuredX = allHitsFirstLayerMeasuredX[firsthit];
            hitsForFit.firstLayerMeasuredY = allHitsFirstLayerMeasuredY[firsthit];
            hitsForFit.firstLayerMeasuredZ = allHitsFirstLayerMeasuredZ[firsthit];

            hitsForFit.secondLayerMeasuredX = allHitsSecondLayerMeasuredX[take];
            hitsForFit.secondLayerMeasuredY = allHitsSecondLayerMeasuredY[take];
            hitsForFit.secondLayerMeasuredZ = allHitsSecondLayerMeasuredZ[take];

            hitsForFit.secondLayerPredictedX = allHitsFirstLayerMeasuredX[firsthit];
            hitsForFit.secondLayerPredictedY = allHitsFirstLayerMeasuredY[firsthit];
            hitsForFit.secondLayerPredictedZ = allHitsSecondLayerMeasuredZ[take];

            hitsForFit.firstLayerResolution = allHitsFirstLayerResolution[firsthit];
            hitsForFit.secondLayerResolution = allHitsSecondLayerResolution[take];

            _hitsForFit.push_back(hitsForFit);

          }

        } // end loop over hits in first plane

      } else if (_alignedBox == 2) {

        // loop over all hits predicted from first box
        for (int firsthit = 0; firsthit < nHitsSecondBox; firsthit++) {

          int take = -1000;
          int veto = -1000;

          // loop over all hits in second plane
          for (int secondhit = 0; secondhit < nHitsSecondPlane; secondhit++) {

            // calculate distance between hits
            double distance = sqrt(pow(allHitsSecondBoxX[firsthit] - allHitsSecondLayerMeasuredX[secondhit],2) + pow(allHitsSecondBoxY[firsthit] - allHitsSecondLayerMeasuredY[secondhit],2));

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

            if ( AIDA::IHistogram1D* distance_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[_distanceLocalname]) )
              distance_histo->fill(distance);
            else {
              streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " <<  _distanceLocalname << endl;
            }

#endif

            if (distance >= _distanceMin && distance <= _distanceMax) {
              if (take != -1000) {
                veto = 1;
              }
              take = secondhit;
            }

          } // end loop over hits in second plane

          if (take != -1000 && veto == -1000) {

            hitsForFit.secondLayerMeasuredX = allHitsSecondLayerMeasuredX[take];
            hitsForFit.secondLayerMeasuredY = allHitsSecondLayerMeasuredY[take];
            hitsForFit.secondLayerMeasuredZ = allHitsSecondLayerMeasuredZ[take];

            hitsForFit.secondLayerPredictedX = allHitsSecondBoxX[firsthit];
            hitsForFit.secondLayerPredictedY = allHitsSecondBoxY[firsthit];
            hitsForFit.secondLayerPredictedZ = allHitsSecondBoxZ[take];

            hitsForFit.firstLayerResolution = allHitsFirstLayerResolution[firsthit];
            hitsForFit.secondLayerResolution = allHitsSecondLayerResolution[take];

            _hitsForFit.push_back(hitsForFit);

          }

        } // end loop over hits in first plane

      } // end if _alignedBox == 2

    } // end if check number of hits

      // _hitsForFit.push_back(hitsForFit);

  } catch (DataNotAvailableException& e) {
    streamlog_out  ( WARNING2 ) <<  "No input collection found on event " << event->getEventNumber() << " in run " << event->getRunNumber() << endl;
  }

  // EUTelMultiLineFit::FitTrack();

  ++_iEvt;

  if ( isFirstEvent() ) _isFirstEvent = false;

  streamlog_out ( MESSAGE2 ) << "Read event: " << _iEvt << endl;
  streamlog_out ( MESSAGE2 ) << "Number of hits in first plane: " << nHitsFirstPlane << endl;
  streamlog_out ( MESSAGE2 ) << "Number of hits in the last plane: " << nHitsSecondPlane << endl;
  streamlog_out ( MESSAGE2 ) << "Hit pairs found so far: " << _hitsForFit.size() << endl;

}

void EUTelAlign::Chi2Function(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {

  // Chi2 function for Minuit
  // ------------------------
  //
  // Parameters:
  // par[0]:       off_x
  // par[1]:       off_y
  // par[2]:       theta_x
  // par[3]:       theta_y
  // par[4]:       theta_z

  double x, y, distance;
  double chi2 = 0.0;
  int usedevents = 0;

  // loop over all events
  for (size_t i = 0; i < _hitsForFit.size(); i++) {

    // just to be sure
    distance = 0.0;

    x = (cos(par[3])*cos(par[4])) * _hitsForFit[i].secondLayerMeasuredX + ((-1)*sin(par[2])*sin(par[3])*cos(par[4]) + cos(par[2])*sin(par[4])) * _hitsForFit[i].secondLayerMeasuredY + par[0];
    y = ((-1)*cos(par[3])*sin(par[4])) * _hitsForFit[i].secondLayerMeasuredX + (sin(par[2])*sin(par[3])*sin(par[4]) + cos(par[2])*cos(par[4])) * _hitsForFit[i].secondLayerMeasuredY + par[1];

    distance = ((x - _hitsForFit[i].secondLayerPredictedX) * (x - _hitsForFit[i].secondLayerPredictedX) + (y - _hitsForFit[i].secondLayerPredictedY) * (y - _hitsForFit[i].secondLayerPredictedY)) / 100;

    // add chi2 cut here later?
    if (par[5] == 0.0) {

      chi2 = chi2 + distance;
      usedevents++;

    } else if (distance < par[5]) {

      chi2 = chi2 + distance;
      usedevents++;

    }

  } // end loop over all events

  // streamlog_out ( MESSAGE2) << usedevents << " ";

  f = chi2;

}

void EUTelAlign::end() {

  streamlog_out ( MESSAGE2 ) << "Number of Events used in the fit: " << _hitsForFit.size() << endl;

  streamlog_out ( MESSAGE2 ) << "Minuit will soon be started" << endl;

  // run MINUIT
  // ----------

  gSystem->Load("libMinuit");

  // init Minuit for 6 parameters
  TMinuit *gMinuit = new TMinuit(6);

  // set print level (-1 = quiet, 0 = normal, 1 = verbose)
  gMinuit->SetPrintLevel(0);

  // set function
  gMinuit->SetFCN(fitFunctionWrapper);

  double arglist[10];
  int ierflag = 0;

  // minimization strategy (1 = standard, 2 = slower)
  arglist[0] = 1;
  gMinuit->mnexcm("SET STR",arglist,2,ierflag);

  // set error definition (1 = for chi square)
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR",arglist,1,ierflag);

  double start_off_x = _startValuesForAlignment[0];
  double start_off_y = _startValuesForAlignment[1];
  double start_theta_x = _startValuesForAlignment[2];
  double start_theta_y = _startValuesForAlignment[3];
  double start_theta_z = _startValuesForAlignment[4];
  double start_chi2 = _chi2Cut;

  // set starting values and step sizes
  gMinuit->mnparm(0,"off_x",start_off_x,1,0,0,ierflag);
  gMinuit->mnparm(1,"off_y",start_off_y,1,0,0,ierflag);
  gMinuit->mnparm(2,"theta_x",start_theta_x,0.001,0,0,ierflag);
  gMinuit->mnparm(3,"theta_y",start_theta_y,0.001,0,0,ierflag);
  gMinuit->mnparm(4,"theta_z",start_theta_z,0.001,0,0,ierflag);
  gMinuit->mnparm(5,"chi2",0.0,1,0,0,ierflag);

  gMinuit->FixParameter(2);
  gMinuit->FixParameter(3);
  gMinuit->FixParameter(4);
  gMinuit->FixParameter(5);

  streamlog_out ( MESSAGE2 ) << endl << "First iteration of alignment: only offsets" << endl;
  streamlog_out ( MESSAGE2 ) << "------------------------------------------" << endl << endl;

  // call migrad (2000 iterations, 0.1 = tolerance)
  arglist[0] = 2000;
  arglist[1] = 0.1;
  gMinuit->mnexcm("MIGRAD",arglist,1,ierflag);

  // calculate errors using MINOS
  arglist[0] = 2000;
  arglist[1] = 0.1;
  gMinuit->mnexcm("MINOS",arglist,1,ierflag);

  double off_x_simple = 0.0;
  double off_y_simple = 0.0;

  double off_x_simple_error = 0.0;
  double off_y_simple_error = 0.0;

  // get results from migrad
  gMinuit->GetParameter(0,off_x_simple,off_x_simple_error);
  gMinuit->GetParameter(1,off_y_simple,off_y_simple_error);

  // fill histograms
  double residual_x_simple = 1000.0;
  double residual_y_simple = 1000.0;

  // loop over all events
  for (size_t i = 0; i < _hitsForFit.size(); i++) {

    residual_x_simple = off_x_simple + _hitsForFit[i].secondLayerMeasuredX - _hitsForFit[i].secondLayerPredictedX;
    residual_y_simple = off_y_simple + _hitsForFit[i].secondLayerMeasuredY - _hitsForFit[i].secondLayerPredictedY;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

    if ( AIDA::IHistogram1D* residx_simple_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[_residualXSimpleLocalname]) )
      residx_simple_histo->fill(residual_x_simple);
    else {
      streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " <<  _residualXSimpleLocalname << endl;
    }

    if ( AIDA::IHistogram1D* residy_simple_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[_residualYSimpleLocalname]) )
      residy_simple_histo->fill(residual_y_simple);
    else {
      streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " <<  _residualYSimpleLocalname << endl;
    }

#endif

  } // end loop over all events

  // get results from migrad
  gMinuit->GetParameter(0,off_x_simple,off_x_simple_error);
  gMinuit->GetParameter(1,off_y_simple,off_y_simple_error);

  // release angles
  gMinuit->Release(2);
  gMinuit->Release(3);
  gMinuit->Release(4);

  streamlog_out ( MESSAGE2 ) << endl << "Second iteration of alignment: include angles" << endl;
  streamlog_out ( MESSAGE2 ) << "---------------------------------------------" << endl << endl;

  // call migrad (2000 iterations, 0.1 = tolerance)
  arglist[0] = 2000;
  arglist[1] = 0.1;
  gMinuit->mnexcm("MIGRAD",arglist,1,ierflag);

  // calculate errors using MINOS
  arglist[0] = 2000;
  arglist[1] = 0.1;
  gMinuit->mnexcm("MINOS",arglist,1,ierflag);

  streamlog_out ( MESSAGE2 ) << endl << "Third iteration of alignment: include chi^2 cut" << endl;
  streamlog_out ( MESSAGE2 ) << "-----------------------------------------------" << endl << endl;

  // release chi2
  gMinuit->mnparm(5,"chi2",start_chi2,1,0,0,ierflag);

  // call migrad (2000 iterations, 0.1 = tolerance)
  arglist[0] = 2000;
  arglist[1] = 0.1;
  gMinuit->mnexcm("MIGRAD",arglist,1,ierflag);

  // calculate errors using MINOS
  arglist[0] = 2000;
  arglist[1] = 0.1;
  gMinuit->mnexcm("MINOS",arglist,1,ierflag);

  streamlog_out ( MESSAGE2) << endl;

  double off_x = 0.0;
  double off_y = 0.0;
  double theta_x = 0.0;
  double theta_y = 0.0;
  double theta_z = 0.0;

  double off_x_error = 0.0;
  double off_y_error = 0.0;
  double theta_x_error = 0.0;
  double theta_y_error = 0.0;
  double theta_z_error = 0.0;

  // get results from migrad
  gMinuit->GetParameter(0,off_x,off_x_error);
  gMinuit->GetParameter(1,off_y,off_y_error);
  gMinuit->GetParameter(2,theta_x,theta_x_error);
  gMinuit->GetParameter(3,theta_y,theta_y_error);
  gMinuit->GetParameter(4,theta_z,theta_z_error);

  streamlog_out ( MESSAGE2 ) << endl << "Alignment constants from the fit:" << endl;
  streamlog_out ( MESSAGE2 ) << "---------------------------------" << endl;
  streamlog_out ( MESSAGE2 ) << "off_x: " << off_x << " +/- " << off_x_error << endl;
  streamlog_out ( MESSAGE2 ) << "off_y: " << off_y << " +/- " << off_y_error << endl;
  streamlog_out ( MESSAGE2 ) << "theta_x: " << theta_x << " +/- " << theta_x_error << endl;
  streamlog_out ( MESSAGE2 ) << "theta_y: " << theta_y << " +/- " << theta_y_error << endl;
  streamlog_out ( MESSAGE2 ) << "theta_z: " << theta_z << " +/- " << theta_z_error << endl;
  streamlog_out ( MESSAGE2 ) << "For copy and paste to line fit xml-file: " << off_x << " " << off_y << " " << theta_x << " " << theta_y << " " << theta_z << endl;

  // fill histograms
  // ---------------

  double x,y;
  double residual_x = 1000.0;
  double residual_y = 1000.0;

  // loop over all events
  for (size_t i = 0; i < _hitsForFit.size(); i++) {

    x = (cos(theta_y)*cos(theta_z)) * _hitsForFit[i].secondLayerMeasuredX + ((-1)*sin(theta_x)*sin(theta_y)*cos(theta_z) + cos(theta_x)*sin(theta_z)) * _hitsForFit[i].secondLayerMeasuredY + off_x;
    y = ((-1)*cos(theta_y)*sin(theta_z)) * _hitsForFit[i].secondLayerMeasuredX + (sin(theta_x)*sin(theta_y)*sin(theta_z) + cos(theta_x)*cos(theta_z)) * _hitsForFit[i].secondLayerMeasuredY + off_y;

    residual_x = x - _hitsForFit[i].secondLayerPredictedX;
    residual_y = y - _hitsForFit[i].secondLayerPredictedY;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

    if ( AIDA::IHistogram1D* residx_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[_residualXLocalname]) )
      residx_histo->fill(residual_x);
    else {
      streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " <<  _residualXLocalname << endl;
    }

    if ( AIDA::IHistogram1D* residy_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[_residualYLocalname]) )
      residy_histo->fill(residual_y);
    else {
      streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " <<  _residualYLocalname << endl;
    }

#endif

  } // end loop over all events

  delete [] _intrResolY;
  delete [] _intrResolX;
  delete [] _xMeasPos;
  delete [] _yMeasPos;
  delete [] _zMeasPos;

  streamlog_out ( MESSAGE2 ) << "Successfully finished" << endl;

}

void EUTelAlign::bookHistos() {

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  streamlog_out ( MESSAGE2 ) << "Booking histograms" << endl;

  AIDA::IHistogram1D * distanceLocal =
    AIDAProcessor::histogramFactory(this)->createHistogram1D(_distanceLocalname,1000,0.0,5000.0);
  if ( distanceLocal ) {
    distanceLocal->setTitle("Distance");
    _aidaHistoMap.insert( make_pair( _distanceLocalname, distanceLocal ) );
  } else {
    streamlog_out ( ERROR2 ) << "Problem booking the " << (_distanceLocalname) << endl;
  }

  const int    NBin = 2000;
  const double Min  = -1000.0;
  const double Max  = 1000.0;

  AIDA::IHistogram1D * residualXSimpleLocal =
    AIDAProcessor::histogramFactory(this)->createHistogram1D(_residualXSimpleLocalname,NBin,Min,Max);
  if ( residualXSimpleLocal ) {
    residualXSimpleLocal->setTitle("Residual X - only offsets");
    _aidaHistoMap.insert( make_pair( _residualXSimpleLocalname, residualXSimpleLocal ) );
  } else {
    streamlog_out ( ERROR2 ) << "Problem booking the " << (_residualXSimpleLocalname) << endl;
  }

  AIDA::IHistogram1D * residualYSimpleLocal =
    AIDAProcessor::histogramFactory(this)->createHistogram1D(_residualYSimpleLocalname,NBin,Min,Max);
  if ( residualYSimpleLocal ) {
    residualYSimpleLocal->setTitle("Residual Y - only offsets");
    _aidaHistoMap.insert( make_pair( _residualYSimpleLocalname, residualYSimpleLocal ) );
  } else {
    streamlog_out ( ERROR2 ) << "Problem booking the " << (_residualYSimpleLocalname) << endl;
  }

  AIDA::IHistogram1D * residualXLocal =
    AIDAProcessor::histogramFactory(this)->createHistogram1D(_residualXLocalname,NBin,Min,Max);
  if ( residualXLocal ) {
    residualXLocal->setTitle("Residual X");
    _aidaHistoMap.insert( make_pair( _residualXLocalname, residualXLocal ) );
  } else {
    streamlog_out ( ERROR2 ) << "Problem booking the " << (_residualXLocalname) << endl;
  }

  AIDA::IHistogram1D * residualYLocal =
    AIDAProcessor::histogramFactory(this)->createHistogram1D(_residualYLocalname,NBin,Min,Max);
  if ( residualYLocal ) {
    residualYLocal->setTitle("Residual Y");
    _aidaHistoMap.insert( make_pair( _residualYLocalname, residualYLocal ) );
  } else {
    streamlog_out ( ERROR2 ) << "Problem booking the " << (_residualYLocalname) << endl;
  }

#endif

}

#endif // USE_GEAR

#endif // USE_ROOT

#endif // OBSOLETE
