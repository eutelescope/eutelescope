// Author Philipp Roloff, DESY <mailto:philipp.roloff@desy.de>
// Version: $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// built only if GEAR is used
#ifdef USE_GEAR

// eutelescope includes ".h"
#include "EUTelMultiLineFit.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelSparseCluster2Impl.h"
#include "EUTelExceptions.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"
#include "marlin/Exceptions.h"
#include "marlin/AIDAProcessor.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>

// system includes <>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>

using namespace std;
using namespace lcio;
using namespace marlin;
// using namespace gear;
using namespace eutelescope;

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelMultiLineFit::_numberTracksLocalname   = "NumberTracks";
std::string EUTelMultiLineFit::_chi2XLocalname          = "Chi2X";
std::string EUTelMultiLineFit::_chi2YLocalname          = "Chi2Y";
std::string EUTelMultiLineFit::_angleXLocalname         = "AngleX";
std::string EUTelMultiLineFit::_angleYLocalname         = "AngleY";
std::string EUTelMultiLineFit::_residualXLocalname      = "ResidualX";
std::string EUTelMultiLineFit::_residualYLocalname      = "ResidualY";
std::string EUTelMultiLineFit::_positionXYLocalname     = "PositionXY";
std::string EUTelMultiLineFit::_seedChargeLocalname     = "SeedCharge";
std::string EUTelMultiLineFit::_clusterChargeLocalname  = "ClusterCharge";
std::string EUTelMultiLineFit::_hitDistanceLocalname    = "HitDistance";
#endif

EUTelMultiLineFit::EUTelMultiLineFit () : Processor("EUTelMultiLineFit") {

  // modify processor description
  _description =
    "EUTelMultiLineFit will fit several straight lines";

  // input collection

  registerInputCollection(LCIO::TRACKERHIT,"HitCollectionName",
                          "Hit collection name",
                          _hitCollectionName, string ( "hit" ));

  // output collection

  registerOutputCollection(LCIO::TRACK,"OutputTrackCollectionName",
                           "Collection name for fitted tracks",
                           _outputTrackColName, string ("fittracks"));

  registerOutputCollection(LCIO::TRACKERHIT,"CorrectedHitCollectionName",
                           "Collection name for corrected particle positions",
                           _correctedHitColName, string ("corrfithits"));

  registerOutputCollection(LCIO::TRACKERHIT,"OutputHitCollectionName",
                           "Collection name for fitted particle hits (positions)",
                           _outputHitColName, string ("fithits"));

  registerOptionalParameter("AlignmentMode","1 for constants from EUTelAlign (default), 2 for Millepede."
                            ,_alignmentMode, static_cast <int> (1));

  // input parameters: take these from database later

  FloatVec constantsFirstLayer;
  constantsFirstLayer.push_back(0.0);
  constantsFirstLayer.push_back(0.0);
  constantsFirstLayer.push_back(0.0);
  constantsFirstLayer.push_back(0.0);
  constantsFirstLayer.push_back(0.0);

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

  FloatVec constantsFourthLayer;
  constantsFourthLayer.push_back(0.0);
  constantsFourthLayer.push_back(0.0);
  constantsFourthLayer.push_back(0.0);
  constantsFourthLayer.push_back(0.0);
  constantsFourthLayer.push_back(0.0);

  FloatVec constantsFifthLayer;
  constantsFifthLayer.push_back(0.0);
  constantsFifthLayer.push_back(0.0);
  constantsFifthLayer.push_back(0.0);
  constantsFifthLayer.push_back(0.0);
  constantsFifthLayer.push_back(0.0);

  FloatVec constantsSixthLayer;
  constantsSixthLayer.push_back(0.0);
  constantsSixthLayer.push_back(0.0);
  constantsSixthLayer.push_back(0.0);
  constantsSixthLayer.push_back(0.0);
  constantsSixthLayer.push_back(0.0);

  registerOptionalParameter("AlignmentConstantsSecondLayer","Alignment Constants for second Telescope Layer:\n off_x, off_y, theta_x, theta_y, theta_z",
                            _alignmentConstantsSecondLayer, constantsSecondLayer);
  registerOptionalParameter("AlignmentConstantsThirdLayer","Alignment Constants for third Telescope Layer:\n off_x, off_y, theta_x, theta_y, theta_z",
                            _alignmentConstantsThirdLayer, constantsThirdLayer);
  registerOptionalParameter("AlignmentConstantsFourthLayer","Alignment Constants for fourth Telescope Layer:\n off_x, off_y, theta_x, theta_y, theta_z",
                            _alignmentConstantsFourthLayer, constantsFourthLayer);
  registerOptionalParameter("AlignmentConstantsFifthLayer","Alignment Constants for fifth Telescope Layer:\n off_x, off_y, theta_x, theta_y, theta_z"
                            ,_alignmentConstantsFifthLayer, constantsFifthLayer);
  registerOptionalParameter("AlignmentConstantsSixthLayer","Alignment Constants for sixth Telescope Layer:\n off_x, off_y, theta_x, theta_y, theta_z"
                            ,_alignmentConstantsSixthLayer, constantsSixthLayer);

  registerOptionalParameter("MillepedeConstantsFirstLayer","Millepede Constants for first Telescope Layer:\n x0, y0, alpha, beta, gamma",
                            _millepedeConstantsFirstLayer, constantsFirstLayer);
  registerOptionalParameter("MillepedeConstantsSecondLayer","Millepede Constants for second Telescope Layer:\n x0, y0, alpha, beta, gamme",
                            _millepedeConstantsSecondLayer, constantsSecondLayer);
  registerOptionalParameter("MillepedeConstantsThirdLayer","Millepede Constants for third Telescope Layer:\n x0, y0, alpha, beta, gamma",
                            _millepedeConstantsThirdLayer, constantsThirdLayer);
  registerOptionalParameter("MillepedeConstantsFourthLayer","Millepede Constants for fourth Telescope Layer:\n x0, y0, alpha, beta, gamma",
                            _millepedeConstantsFourthLayer, constantsFourthLayer);
  registerOptionalParameter("MillepedeConstantsFifthLayer","Millepede Constants for fifth Telescope Layer:\n x0, y0, alpha, beta, gamma"
                            ,_millepedeConstantsFifthLayer, constantsFifthLayer);
  registerOptionalParameter("MillepedeConstantsSixthLayer","Millepede Constants for sixth Telescope Layer:\n x0, y0, alpha, beta, gamma"
                            ,_millepedeConstantsSixthLayer, constantsSixthLayer);

  registerOptionalParameter("DistanceMax","Maximal allowed distance between hits entering the fit per 10 cm space between the planes.",
                            _distanceMax, static_cast <float> (2000.0));

  registerOptionalParameter("DistanceDUTMax","Distance used for DUT hit matching.",
			    _distanceDUTMax, static_cast <float> (50.0));

  registerOptionalParameter("Chi2XMax","Maximal chi2 for fit of x coordinate."
                            ,_chi2XMax, static_cast <float> (10000.0));

  registerOptionalParameter("Chi2YMax","Maximal chi2 for fit of y coordinate."                             ,_chi2YMax, static_cast <float> (10000.0));

  registerOptionalParameter("ExcludePlane","Exclude plane from fit."
                            ,_excludePlane, static_cast <int> (0));

  registerOptionalParameter("MaxTrackCandidates","Maximal number of track candidates"
                            ,_maxTrackCandidates, static_cast <int> (2000));

  registerOptionalParameter("MaxHitsPlane","Maximal number of hits per plane"
                            ,_maxHitsPlane, static_cast <int> (100));

  registerOptionalParameter("HitDistanceXMax","Maximal |x| for calculation of efficiencies",
                            _hitDistanceXMax, static_cast <float> (2500.0));

  registerOptionalParameter("HitDistanceYMax","Maximal |y| for calculation of efficiencies",
                            _hitDistanceYMax, static_cast <float> (2500.0));

}

void EUTelMultiLineFit::init() {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

  // set to zero hit efficiency counters
  _iHitDUT = 0;
  for ( int help = 0; help < 10; help++) {
    _iHitDUTDistanceCut[help] = 0;
  }

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

  _intrResolX = new double[_nPlanes];
  _intrResolY = new double[_nPlanes];

}

void EUTelMultiLineFit::processRunHeader (LCRunHeader * rdr) {


  auto_ptr<EUTelRunHeaderImpl> header ( new EUTelRunHeaderImpl (rdr) );
  header->addProcessor( type() ) ;

  // the run header contains the number of detectors. This number
  // should be in principle the same as the number of layers in the
  // geometry description
  if ( header->getNoOfDetector() != _siPlanesParameters->getSiPlanesNumber() ) {
    streamlog_out ( ERROR2 ) << "Error during the geometry consistency check: " << endl;
    streamlog_out ( ERROR2 ) << "The run header says there are " << header->getNoOfDetector() << " silicon detectors " << endl;
    streamlog_out ( ERROR2 ) << "The GEAR description says     " << _siPlanesParameters->getSiPlanesNumber() << " silicon planes" << endl;
    string answer;
    while (true) {
      streamlog_out ( ERROR2 ) << "Type Q to quit now or C to continue [Q/C]" << endl;
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
    streamlog_out ( ERROR2 ) << "The GEAR description says is     " << _siPlanesParameters->getSiPlanesID() << endl;

#ifdef EUTEL_INTERACTIVE
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
#endif
}

  // now book histograms plz...
  if ( isFirstEvent() )  bookHistos();

  // increment the run counter
  ++_iRun;
}

void EUTelMultiLineFit::FitTrack(int nPlanesFitter, double xPosFitter[], double yPosFitter[], double zPosFitter[], double xResFitter[], double yResFitter[], double chi2Fit[2], double residXFit[], double residYFit[], double angleFit[2]) {

  int sizearray;

  if (_excludePlane == 0) {
    sizearray = nPlanesFitter;
  } else {
    sizearray = nPlanesFitter - 1;
  }

  double * xPosFit = new double[sizearray];
  double * yPosFit = new double[sizearray];
  double * zPosFit = new double[sizearray];
  double * xResFit = new double[sizearray];
  double * yResFit = new double[sizearray];

  int nPlanesFit = 0;

  for (int help = 0; help < nPlanesFitter; help++) {
    if ((_excludePlane - 1) == help) {
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

void EUTelMultiLineFit::processEvent (LCEvent * event) {
  
  if ( _iEvt % 10  == 0 ) {
    streamlog_out( MESSAGE2 ) << "Processing event " 	 
                              << setw(6) << setiosflags(ios::right) << event->getEventNumber() << " in run " 	 
                              << setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber() << setfill(' ') 	 
                              << " (Total = " << setw(10) << _iEvt << ")" << resetiosflags(ios::left) << endl; 	 
    
 }

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event) ;

  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG2 ) << "EORE found: nothing else to do." << endl;
    return;
  }

  try {

    LCCollectionVec * hitCollection     = static_cast<LCCollectionVec*> (event->getCollection( _hitCollectionName ));

    int detectorID    = -99; // it's a non sense
    int oldDetectorID = -100;
    int layerIndex;

    double seedCharge = -1000.0;
    double clusterCharge = -1000.0;

    vector<EUTelMultiLineFit::HitsInPlane > _hitsFirstPlane;
    vector<EUTelMultiLineFit::HitsInPlane > _hitsSecondPlane;
    vector<EUTelMultiLineFit::HitsInPlane > _hitsThirdPlane;
    vector<EUTelMultiLineFit::HitsInPlane > _hitsFourthPlane;
    vector<EUTelMultiLineFit::HitsInPlane > _hitsFifthPlane;
    vector<EUTelMultiLineFit::HitsInPlane > _hitsSixthPlane;

    HitsInPlane hitsInPlane;

    for ( int iHit = 0; iHit < hitCollection->getNumberOfElements(); iHit++ ) {

      TrackerHitImpl * hit = static_cast<TrackerHitImpl*> ( hitCollection->getElementAt(iHit) );

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

      seedCharge = cluster->getSeedCharge();
      clusterCharge = cluster->getTotalCharge();

      detectorID = cluster->getDetectorID();

      if ( detectorID != oldDetectorID ) {
        oldDetectorID = detectorID;

        if ( _conversionIdMap.size() != static_cast< unsigned >(_siPlanesParameters->getSiPlanesNumber()) ) {
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
        _intrResolX[layerIndex] = 1000*_siPlanesLayerLayout->getSensitiveResolution(layerIndex); //um
        _intrResolY[layerIndex] = 1000*_siPlanesLayerLayout->getSensitiveResolution(layerIndex); //um

      }

      layerIndex   = _conversionIdMap[detectorID];

      // Getting positions of the hits.
      // ------------------------------
      //
      // Here the alignment constants are used to correct the positions.
      // The other layers were aligned with respect to the first one.
      // All distances are given in \mu m.

      // "Simple" alignment constants
      double off_x, off_y, theta_x, theta_y, theta_z;

      // Millepede alignment constants
      double x0, y0, alpha, beta, gamma;

      if (layerIndex == 0) {

        off_x = 0.0;
        off_y = 0.0;
        theta_x = 0.0;
        theta_y = 0.0;
        theta_z = 0.0;

        x0 = _millepedeConstantsFirstLayer[0];
        y0 = _millepedeConstantsFirstLayer[1];
        alpha = _millepedeConstantsFirstLayer[2];
        beta = _millepedeConstantsFirstLayer[3];
        gamma = _millepedeConstantsFirstLayer[4];

      } else if (layerIndex == 1) {

        off_x = _alignmentConstantsSecondLayer[0];
        off_y = _alignmentConstantsSecondLayer[1];
        theta_x = _alignmentConstantsSecondLayer[2];
        theta_y = _alignmentConstantsSecondLayer[3];
        theta_z = _alignmentConstantsSecondLayer[4];

        x0 = _millepedeConstantsSecondLayer[0];
        y0 = _millepedeConstantsSecondLayer[1];
        alpha = _millepedeConstantsSecondLayer[2];
        beta = _millepedeConstantsSecondLayer[3];
        gamma = _millepedeConstantsSecondLayer[4];

      } else if (layerIndex == 2) {

        off_x = _alignmentConstantsThirdLayer[0];
        off_y = _alignmentConstantsThirdLayer[1];
        theta_x = _alignmentConstantsThirdLayer[2];
        theta_y = _alignmentConstantsThirdLayer[3];
        theta_z = _alignmentConstantsThirdLayer[4];

        x0 = _millepedeConstantsThirdLayer[0];
        y0 = _millepedeConstantsThirdLayer[1];
        alpha = _millepedeConstantsThirdLayer[2];
        beta = _millepedeConstantsThirdLayer[3];
        gamma = _millepedeConstantsThirdLayer[4];

      } else if (layerIndex == 3) {

        off_x = _alignmentConstantsFourthLayer[0];
        off_y = _alignmentConstantsFourthLayer[1];
        theta_x = _alignmentConstantsFourthLayer[2];
        theta_y = _alignmentConstantsFourthLayer[3];
        theta_z = _alignmentConstantsFourthLayer[4];

        x0 = _millepedeConstantsFourthLayer[0];
        y0 = _millepedeConstantsFourthLayer[1];
        alpha = _millepedeConstantsFourthLayer[2];
        beta = _millepedeConstantsFourthLayer[3];
        gamma = _millepedeConstantsFourthLayer[4];

      } else if (layerIndex == 4) {

        off_x = _alignmentConstantsFifthLayer[0];
        off_y = _alignmentConstantsFifthLayer[1];
        theta_x = _alignmentConstantsFifthLayer[2];
        theta_y = _alignmentConstantsFifthLayer[3];
        theta_z = _alignmentConstantsFifthLayer[4];

        x0 = _millepedeConstantsFifthLayer[0];
        y0 = _millepedeConstantsFifthLayer[1];
        alpha = _millepedeConstantsFifthLayer[2];
        beta = _millepedeConstantsFifthLayer[3];
        gamma = _millepedeConstantsFifthLayer[4];

      } else if (layerIndex == 5) {

        off_x = _alignmentConstantsSixthLayer[0];
        off_y = _alignmentConstantsSixthLayer[1];
        theta_x = _alignmentConstantsSixthLayer[2];
        theta_y = _alignmentConstantsSixthLayer[3];
        theta_z = _alignmentConstantsSixthLayer[4];

        x0 = _millepedeConstantsSixthLayer[0];
        y0 = _millepedeConstantsSixthLayer[1];
        alpha = _millepedeConstantsSixthLayer[2];
        beta = _millepedeConstantsSixthLayer[3];
        gamma = _millepedeConstantsSixthLayer[4];

      } else {

        off_x = 0.0;
        off_y = 0.0;
        theta_x = 0.0;
        theta_y = 0.0;
        theta_z = 0.0;

        x0 = 0.0;
        y0 = 0.0;
        alpha = 0.0;
        beta = 0.0;
        gamma = 0.0;

      }

      if (_alignmentMode != 2) {

        // For documentation of these formulas look at EUDET-Memo-2007-20

        hitsInPlane.measuredX = (cos(theta_y)*cos(theta_z)) * hit->getPosition()[0] * 1000 + ((-1)*sin(theta_x)*sin(theta_y)*cos(theta_z) + cos(theta_x)*sin(theta_z)) * hit->getPosition()[1] * 1000 + off_x;
        hitsInPlane.measuredY = ((-1)*cos(theta_y)*sin(theta_z)) * hit->getPosition()[0] * 1000 + (sin(theta_x)*sin(theta_y)*sin(theta_z) + cos(theta_x)*cos(theta_z)) * hit->getPosition()[1] * 1000 + off_y;
        hitsInPlane.measuredZ = 1000 * hit->getPosition()[2];
        hitsInPlane.seedCharge = seedCharge;
        hitsInPlane.clusterCharge = clusterCharge;

      } else {

        hitsInPlane.measuredX = 1000 * hit->getPosition()[0] + gamma * 1000 * hit->getPosition()[1] + beta * 1000 * hit->getPosition()[2] - x0;
        hitsInPlane.measuredY = (-1) * gamma * 1000 * hit->getPosition()[0] + 1000 * hit->getPosition()[1] + alpha * 1000 * hit->getPosition()[2] - y0;
        hitsInPlane.measuredZ = 1000 * hit->getPosition()[2];
        hitsInPlane.seedCharge = seedCharge;
        hitsInPlane.clusterCharge = clusterCharge;

      }

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

    }

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
    _seedCharge = new double *[_maxTrackCandidates];
    _clusterCharge = new double *[_maxTrackCandidates];

    for (int help = 0; help < _maxTrackCandidates; help++) {
      _xPos[help] = new double[_nPlanes];
      _yPos[help] = new double[_nPlanes];
      _zPos[help] = new double[_nPlanes];
      _seedCharge[help] = new double[_nPlanes];
      _clusterCharge[help] = new double[_nPlanes];
    }

    // Currently unused variable:
    /* int fitplane[6] = {0, 0, 0, 0, 0, 0};

    for (int help = 0; help < _nPlanes; help++) {
      fitplane[help] = 1;
    }
    */

    int _nTracks = 0;

    int _nGoodTracks = 0;

    // Find track candidates using the distance cuts
    // ---------------------------------------------
    //
    // This is done separately for different numbers of planes.

    if (_hitsFirstPlane.size() <= size_t(_maxHitsPlane) && _hitsSecondPlane.size() <= size_t(_maxHitsPlane) && _hitsThirdPlane.size() <= size_t(_maxHitsPlane) && _hitsFourthPlane.size() <= size_t(_maxHitsPlane) && _hitsFifthPlane.size() <= size_t(_maxHitsPlane) && _hitsSixthPlane.size() <= size_t(_maxHitsPlane)) {

      // loop over all hits in first plane
      for (int firsthit = 0; size_t(firsthit) < _hitsFirstPlane.size(); firsthit++) {

        // loop over all hits in second plane
        for (int secondhit = 0; size_t(secondhit) < _hitsSecondPlane.size(); secondhit++) {

          distance12 = sqrt(pow(_hitsFirstPlane[firsthit].measuredX - _hitsSecondPlane[secondhit].measuredX,2) + pow(_hitsFirstPlane[firsthit].measuredY - _hitsSecondPlane[secondhit].measuredY,2));

          distance_plane12 = _hitsSecondPlane[secondhit].measuredZ - _hitsFirstPlane[firsthit].measuredZ;

          distanceMax12 = _distanceMax * (distance_plane12 / 100000.0);

          if (_nPlanes == 2 && distance12 < distanceMax12 && _nTracks < _maxTrackCandidates) {

            _xPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredX;
            _yPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredY;
            _zPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredZ;
            _seedCharge[_nTracks][0] = _hitsFirstPlane[firsthit].seedCharge;
            _clusterCharge[_nTracks][0] = _hitsFirstPlane[firsthit].clusterCharge;

            _xPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredX;
            _yPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredY;
            _zPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredZ;
            _seedCharge[_nTracks][1] = _hitsSecondPlane[secondhit].seedCharge;
            _clusterCharge[_nTracks][1] = _hitsSecondPlane[secondhit].clusterCharge;

            _nTracks++;

          }

          // more than two planes
          if (_nPlanes > 2) {

            // loop over all hits in third plane
            for (int thirdhit = 0; size_t(thirdhit) < _hitsThirdPlane.size(); thirdhit++) {

              distance23 = sqrt(pow(_hitsSecondPlane[secondhit].measuredX - _hitsThirdPlane[thirdhit].measuredX,2) + pow(_hitsSecondPlane[secondhit].measuredY - _hitsThirdPlane[thirdhit].measuredY,2));

              distance_plane23 = _hitsThirdPlane[thirdhit].measuredZ - _hitsSecondPlane[secondhit].measuredZ;

              distanceMax23 = _distanceMax * (distance_plane23 / 100000.0);

              if (_nPlanes == 3 && distance12 < distanceMax12 && distance23 < distanceMax23 && _nTracks < _maxTrackCandidates) {

                _xPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredX;
                _yPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredY;
                _zPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredZ;
                _seedCharge[_nTracks][0] = _hitsFirstPlane[firsthit].seedCharge;
                _clusterCharge[_nTracks][0] = _hitsFirstPlane[firsthit].clusterCharge;

                _xPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredX;
                _yPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredY;
                _zPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredZ;
                _seedCharge[_nTracks][1] = _hitsSecondPlane[secondhit].seedCharge;
                _clusterCharge[_nTracks][1] = _hitsSecondPlane[secondhit].clusterCharge;

                _xPos[_nTracks][2] = _hitsThirdPlane[thirdhit].measuredX;
                _yPos[_nTracks][2] = _hitsThirdPlane[thirdhit].measuredY;
                _zPos[_nTracks][2] = _hitsThirdPlane[thirdhit].measuredZ;
                _seedCharge[_nTracks][2] = _hitsThirdPlane[thirdhit].seedCharge;
                _clusterCharge[_nTracks][2] = _hitsThirdPlane[thirdhit].clusterCharge;

                _nTracks++;

              }

              // more than three planes
              if (_nPlanes > 3) {

                // loop over all hits in fourth plane
                for (int fourthhit = 0; size_t(fourthhit) < _hitsFourthPlane.size(); fourthhit++) {

                  distance34 = sqrt(pow(_hitsThirdPlane[thirdhit].measuredX - _hitsFourthPlane[fourthhit].measuredX,2) + pow(_hitsThirdPlane[thirdhit].measuredY - _hitsFourthPlane[fourthhit].measuredY,2));

                  distance_plane34 = _hitsFourthPlane[fourthhit].measuredZ - _hitsThirdPlane[thirdhit].measuredZ;

                  distanceMax34 = _distanceMax * (distance_plane34 / 100000.0);

                  if (_nPlanes == 4 && distance12 < distanceMax12 && distance23 < distanceMax23 && distance34 < distanceMax34 && _nTracks < _maxTrackCandidates) {

                    _xPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredX;
                    _yPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredY;
                    _zPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredZ;
                    _seedCharge[_nTracks][0] = _hitsFirstPlane[firsthit].seedCharge;
                    _clusterCharge[_nTracks][0] = _hitsFirstPlane[firsthit].clusterCharge;

                    _xPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredX;
                    _yPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredY;
                    _zPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredZ;
                    _seedCharge[_nTracks][1] = _hitsSecondPlane[secondhit].seedCharge;
                    _clusterCharge[_nTracks][1] = _hitsSecondPlane[secondhit].clusterCharge;

                    _xPos[_nTracks][2] = _hitsThirdPlane[thirdhit].measuredX;
                    _yPos[_nTracks][2] = _hitsThirdPlane[thirdhit].measuredY;
                    _zPos[_nTracks][2] = _hitsThirdPlane[thirdhit].measuredZ;
                    _seedCharge[_nTracks][2] = _hitsThirdPlane[thirdhit].seedCharge;
                    _clusterCharge[_nTracks][2] = _hitsThirdPlane[thirdhit].clusterCharge;

                    _xPos[_nTracks][3] = _hitsFourthPlane[fourthhit].measuredX;
                    _yPos[_nTracks][3] = _hitsFourthPlane[fourthhit].measuredY;
                    _zPos[_nTracks][3] = _hitsFourthPlane[fourthhit].measuredZ;
                    _seedCharge[_nTracks][3] = _hitsFourthPlane[fourthhit].seedCharge;
                    _clusterCharge[_nTracks][3] = _hitsFourthPlane[fourthhit].clusterCharge;

                    _nTracks++;

                  }

                  // more than four planes
                  if (_nPlanes > 4) {

                    // loop over all hits in fifth plane
                    for (int fifthhit = 0; size_t(fifthhit) < _hitsFifthPlane.size(); fifthhit++) {

                      distance45 = sqrt(pow(_hitsFourthPlane[fourthhit].measuredX - _hitsFifthPlane[fifthhit].measuredX,2) + pow(_hitsFourthPlane[fourthhit].measuredY - _hitsFifthPlane[fifthhit].measuredY,2));

                      distance_plane45 = _hitsFifthPlane[fifthhit].measuredZ - _hitsFourthPlane[fourthhit].measuredZ;

                      distanceMax45 = _distanceMax * (distance_plane45 / 100000.0);

                      if (_nPlanes == 5 && distance12 < distanceMax12 && distance23 < distanceMax23 && distance34 < distanceMax34 && distance45 < distanceMax45 && _nTracks < _maxTrackCandidates) {

                        _xPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredX;
                        _yPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredY;
                        _zPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredZ;
                        _seedCharge[_nTracks][0] = _hitsFirstPlane[firsthit].seedCharge;
                        _clusterCharge[_nTracks][0] = _hitsFirstPlane[firsthit].clusterCharge;

                        _xPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredX;
                        _yPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredY;
                        _zPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredZ;
                        _seedCharge[_nTracks][1] = _hitsSecondPlane[secondhit].seedCharge;
                        _clusterCharge[_nTracks][1] = _hitsSecondPlane[secondhit].clusterCharge;

                        _xPos[_nTracks][2] = _hitsThirdPlane[thirdhit].measuredX;
                        _yPos[_nTracks][2] = _hitsThirdPlane[thirdhit].measuredY;
                        _zPos[_nTracks][2] = _hitsThirdPlane[thirdhit].measuredZ;
                        _seedCharge[_nTracks][2] = _hitsThirdPlane[thirdhit].seedCharge;
                        _clusterCharge[_nTracks][2] = _hitsThirdPlane[thirdhit].clusterCharge;

                        _xPos[_nTracks][3] = _hitsFourthPlane[fourthhit].measuredX;
                        _yPos[_nTracks][3] = _hitsFourthPlane[fourthhit].measuredY;
                        _zPos[_nTracks][3] = _hitsFourthPlane[fourthhit].measuredZ;
                        _seedCharge[_nTracks][3] = _hitsFourthPlane[fourthhit].seedCharge;
                        _clusterCharge[_nTracks][3] = _hitsFourthPlane[fourthhit].clusterCharge;

                        _xPos[_nTracks][4] = _hitsFifthPlane[fifthhit].measuredX;
                        _yPos[_nTracks][4] = _hitsFifthPlane[fifthhit].measuredY;
                        _zPos[_nTracks][4] = _hitsFifthPlane[fifthhit].measuredZ;
                        _seedCharge[_nTracks][4] = _hitsFifthPlane[fifthhit].seedCharge;
                        _clusterCharge[_nTracks][4] = _hitsFifthPlane[fifthhit].clusterCharge;

                        _nTracks++;

                      }

                      // more than five planes
                      if (_nPlanes > 5) {

                        // loop over all hits in sixth plane
                        for (int sixthhit = 0; size_t(sixthhit) < _hitsSixthPlane.size(); sixthhit++) {

                          distance56 = sqrt(pow(_hitsFifthPlane[fifthhit].measuredX - _hitsSixthPlane[sixthhit].measuredX,2) + pow(_hitsFifthPlane[fifthhit].measuredY - _hitsSixthPlane[sixthhit].measuredY,2));

                          distance_plane56 = _hitsSixthPlane[sixthhit].measuredZ - _hitsFifthPlane[fifthhit].measuredZ;

                          distanceMax56 = _distanceMax * (distance_plane56 / 100000.0);

                          if (_nPlanes == 6 && distance12 < distanceMax12 && distance23 < distanceMax23 && distance34 < distanceMax34 && distance45 < distanceMax45 && distance56 < distanceMax56 && _nTracks < _maxTrackCandidates) {

                            _xPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredX;
                            _yPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredY;
                            _zPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredZ;
                            _seedCharge[_nTracks][0] = _hitsFirstPlane[firsthit].seedCharge;
                            _clusterCharge[_nTracks][0] = _hitsFirstPlane[firsthit].clusterCharge;

                            _xPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredX;
                            _yPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredY;
                            _zPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredZ;
                            _seedCharge[_nTracks][1] = _hitsSecondPlane[secondhit].seedCharge;
                            _clusterCharge[_nTracks][1] = _hitsSecondPlane[secondhit].clusterCharge;

                            _xPos[_nTracks][2] = _hitsThirdPlane[thirdhit].measuredX;
                            _yPos[_nTracks][2] = _hitsThirdPlane[thirdhit].measuredY;
                            _zPos[_nTracks][2] = _hitsThirdPlane[thirdhit].measuredZ;
                            _seedCharge[_nTracks][2] = _hitsThirdPlane[thirdhit].seedCharge;
                            _clusterCharge[_nTracks][2] = _hitsThirdPlane[thirdhit].clusterCharge;

                            _xPos[_nTracks][3] = _hitsFourthPlane[fourthhit].measuredX;
                            _yPos[_nTracks][3] = _hitsFourthPlane[fourthhit].measuredY;
                            _zPos[_nTracks][3] = _hitsFourthPlane[fourthhit].measuredZ;
                            _seedCharge[_nTracks][3] = _hitsFourthPlane[fourthhit].seedCharge;
                            _clusterCharge[_nTracks][3] = _hitsFourthPlane[fourthhit].clusterCharge;

                            _xPos[_nTracks][4] = _hitsFifthPlane[fifthhit].measuredX;
                            _yPos[_nTracks][4] = _hitsFifthPlane[fifthhit].measuredY;
                            _zPos[_nTracks][4] = _hitsFifthPlane[fifthhit].measuredZ;
                            _seedCharge[_nTracks][4] = _hitsFifthPlane[fifthhit].seedCharge;
                            _clusterCharge[_nTracks][4] = _hitsFifthPlane[fifthhit].clusterCharge;

                            _xPos[_nTracks][5] = _hitsSixthPlane[sixthhit].measuredX;
                            _yPos[_nTracks][5] = _hitsSixthPlane[sixthhit].measuredY;
                            _zPos[_nTracks][5] = _hitsSixthPlane[sixthhit].measuredZ;
                            _seedCharge[_nTracks][5] = _hitsSixthPlane[sixthhit].seedCharge;
                            _clusterCharge[_nTracks][5] = _hitsSixthPlane[sixthhit].clusterCharge;

                            _nTracks++;

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

    } else {

      streamlog_out ( WARNING2 ) << "Too many hits. The event will be skipped." << endl;

    }

    if (_nTracks == _maxTrackCandidates) {
      streamlog_out ( WARNING2 ) << "Maximum number of track candidates reached. Maybe further tracks were skipped" << endl;
    }

    streamlog_out ( MULTIFITTERMESSAGE ) << "Number of hits in the individual planes: " << _hitsFirstPlane.size() << " " << _hitsSecondPlane.size() << " " << _hitsThirdPlane.size() << " " << _hitsFourthPlane.size() << " " << _hitsFifthPlane.size() << " " << _hitsSixthPlane.size() << endl;

    streamlog_out ( MULTIFITTERMESSAGE ) << "Number of track candidates found: " << _iEvt << ": " << _nTracks << endl;

    // Perform fit for all found track candidates
    // ------------------------------------------

    // Define output track and hit collections
    LCCollectionVec     * fittrackvec  = new LCCollectionVec(LCIO::TRACK);
    LCCollectionVec     * corrpointvec = new LCCollectionVec(LCIO::TRACKERHIT);
    LCCollectionVec     * fitpointvec  = new LCCollectionVec(LCIO::TRACKERHIT);

    // Set flag for storing track hits in track collection
    LCFlagImpl flag(fittrackvec->getFlag());
    flag.setBit( LCIO::TRBIT_HITS );
    fittrackvec->setFlag(flag.getFlag());

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

      streamlog_out ( MULTIFITTERMESSAGE ) << "Fitting track using the following coordinates: ";

      for (int help = 0; help < _nPlanes; help++) {
        streamlog_out ( MULTIFITTERMESSAGE ) << _xPosHere[help] << " " << _yPosHere[help] << " " << _zPosHere[help] << "   ";
      }

      streamlog_out ( MULTIFITTERMESSAGE ) << endl;

      FitTrack(int(_nPlanes), _xPosHere, _yPosHere, _zPosHere, _intrResolX, _intrResolY, Chiquare, _waferResidX, _waferResidY, angle);

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

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
#endif

      // chi2 cut
      if (Chiquare[0] <= _chi2XMax && Chiquare[1] <= _chi2YMax) {

        _nGoodTracks++;

        // Save track to output file
        // -------------------------

        // Write fit result out

        TrackImpl * fittrack = new TrackImpl();

        // Following parameters are not used for Telescope
        // and are set to zero (just in case)
        fittrack->setOmega(0.);     // curvature of the track
        fittrack->setD0(0.);        // impact paramter of the track in (r-phi)
        fittrack->setZ0(0.);        // impact paramter of the track in (r-z)
        fittrack->setPhi(0.);       // phi of the track at reference point
        fittrack->setTanLambda(0.); // dip angle of the track at reference point

        // Used class members

        fittrack->setChi2(Chiquare[0]);  // x Chi2 of the fit
        fittrack->setNdf(_nPlanes);
        //  fittrack->setNdf(nBestFired); // Number of planes fired (!)

//        fittrack->setIsReferencePointPCA(false);

        // Calculate positions of fitted track in every plane

        int counter;

        for( counter = 0; counter < _nPlanes; counter++ ){

          TrackerHitImpl * corrpoint = new TrackerHitImpl;
          TrackerHitImpl * fitpoint  = new TrackerHitImpl;

          // Plane number stored as hit type
          // corrpoint->setType(counter+1);
          corrpoint->setType(counter+1);

          // Use hit type 32 to be compatible with analytic fitter
          fitpoint->setType(32);

	  // positions in the track are saved in mm as everywhere else
	  // in the code.
          double corrpos[3];
          corrpos[0] = _xPosHere[counter] / 1000.;
          corrpos[1] = _yPosHere[counter] / 1000.;
          corrpos[2] = _zPosHere[counter] / 1000.;

          double fitpos[3];
          fitpos[0] = (_waferResidX[counter] + _xPosHere[counter] ) / 1000.;
          fitpos[1] = (_waferResidY[counter] + _yPosHere[counter] ) / 1000.;
          fitpos[2] = (_zPosHere[counter] ) / 1000.;

          corrpoint->setPosition(corrpos);
          fitpoint->setPosition(fitpos);

          // store corr and fit point

          corrpointvec->push_back(corrpoint);
          fitpointvec->push_back(fitpoint);

          // add points to track

          fittrack->addHit(corrpoint);
          fittrack->addHit(fitpoint);

        }

        fittrackvec->addElement(fittrack);

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

        string tempHistoName;

        // Fill x angle histogram
        if ( _histogramSwitch ) {
          {
            stringstream ss;
            ss << _angleXLocalname << endl;
          }
          if ( AIDA::IHistogram1D* anglex_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[_angleXLocalname]) )
            anglex_histo->fill(angle[0]*180/M_PI);
          else {
            streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _angleXLocalname << endl;
            streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
            _histogramSwitch = false;
          }
        }

        // Fill y angle histogram
        if ( _histogramSwitch ) {
          {
            stringstream ss;
            ss << _angleYLocalname << endl;
          }
          if ( AIDA::IHistogram1D* angley_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[_angleYLocalname]) )
            angley_histo->fill(angle[1]*180/M_PI);
          else {
            streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _angleYLocalname << endl;
            streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
            _histogramSwitch = false;
          }
        }

        // Fill hit distance histogram
        if ( _excludePlane > 0) {

          // exclude border region
          if (fabs(_xPosHere[(_excludePlane - 1)]) < _hitDistanceXMax && fabs(_yPosHere[(_excludePlane - 1)]) < _hitDistanceYMax) {

            double hitDistance = sqrt(_waferResidX[(_excludePlane - 1)] * _waferResidX[(_excludePlane - 1)] + _waferResidY[(_excludePlane - 1)] * _waferResidY[(_excludePlane - 1)]);

            _iHitDUT = _iHitDUT + 1;

            // loop over all distance cuts
            for ( int help = 0; help < 10; help++ ) {
              // distance cut
              if ( hitDistance < (help * 20 + 20)) {
                _iHitDUTDistanceCut[help] = _iHitDUTDistanceCut[help] + 1;
              } // end if distance cut
            } // end loop over all distance cuts

            if ( _histogramSwitch ) {
              {
                stringstream ss;
                ss << _hitDistanceLocalname << endl;
              }
              if ( AIDA::IHistogram1D* hitdistance_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[_hitDistanceLocalname]) )
                hitdistance_histo->fill(hitDistance);
              else {
                streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _hitDistanceLocalname << endl;
                streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
                _histogramSwitch = false;
              }
            }

          } // end if exclude border region

        }

        // loop over all detector planes
        for( int iDetector = 0; iDetector < _nPlanes; iDetector++ ){

	  // distance matching for DUT plane
	  if ((iDetector == (_excludePlane - 1) && sqrt(_waferResidX[iDetector] * _waferResidX[iDetector] + _waferResidY[iDetector] * _waferResidY[iDetector]) < _distanceDUTMax) || iDetector != (_excludePlane - 1)) { 

	    // fill x residuals histograms
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

	    // fill y residuals histograms
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

	  } // end if distance matching for DUT plane

          // fill seed charge histograms
          if ( _histogramSwitch ) {
            {
              stringstream ss;
              ss << _seedChargeLocalname << "_d" << iDetector;
              tempHistoName=ss.str();
            }
            if ( AIDA::IHistogram1D* seedsnr_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName.c_str()]) )
              seedsnr_histo->fill(_seedCharge[track][iDetector]);
            else {
              streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _residualXLocalname << endl;
              streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
              _histogramSwitch = false;
            }
          }

          // fill cluster charge histograms
          if ( _histogramSwitch ) {
            {
              stringstream ss;
              ss << _clusterChargeLocalname << "_d" << iDetector;
              tempHistoName=ss.str();
            }
            if ( AIDA::IHistogram1D* clustersnr_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName.c_str()]) )
              clustersnr_histo->fill(_clusterCharge[track][iDetector]);
            else {
              streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _residualXLocalname << endl;
              streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
              _histogramSwitch = false;
            }
          }

          // fill xy positions histograms
          if ( _histogramSwitch ) {
            {
              stringstream ss;
              ss << _positionXYLocalname << "_d" << iDetector;
              tempHistoName=ss.str();
            }
            if ( AIDA::IHistogram2D* positionxy_histo = dynamic_cast<AIDA::IHistogram2D*>(_aidaHistoMap[tempHistoName.c_str()]) )
              positionxy_histo->fill(_waferResidX[iDetector] + _xPosHere[iDetector],_waferResidY[iDetector] + _yPosHere[iDetector]);
            else {
              streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _positionXYLocalname << endl;
              streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
              _histogramSwitch = false;
            }
          }

        } // end loop over all detector planes

#endif

      } // end if chi2 cut

      // clean up
      delete [] _zPosHere;
      delete [] _yPosHere;
      delete [] _xPosHere;

    } // end loop over all track candidates

    if (_nGoodTracks > 0) {
      event->addCollection(fittrackvec,_outputTrackColName);
      event->addCollection(corrpointvec,_correctedHitColName);
      event->addCollection(fitpointvec,_outputHitColName);
    } else {
      delete fittrackvec;
      delete corrpointvec;
      delete fitpointvec;
    }

    streamlog_out ( MULTIFITTERMESSAGE ) << "Finished fitting tracks in event " << _iEvt << endl;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

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
      delete [] _seedCharge[help];
      delete [] _clusterCharge[help];
    }

    delete [] _zPos;
    delete [] _yPos;
    delete [] _xPos;
    delete [] _seedCharge;
    delete [] _clusterCharge;

  } catch (DataNotAvailableException& e) {

    streamlog_out  ( WARNING2 ) <<  "No input collection found on event " << event->getEventNumber()
                                << " in run " << event->getRunNumber() << endl;

  }

  ++_iEvt;

  if ( isFirstEvent() ) _isFirstEvent = false;

}

void EUTelMultiLineFit::end() {

  delete [] _intrResolY;
  delete [] _intrResolX;
  delete [] _yFitPos;
  delete [] _xFitPos;
  delete [] _waferResidY;
  delete [] _waferResidX;

  streamlog_out ( MESSAGE2 ) << "Successfully finished" << endl;

  // information on DUT efficiency
  if (_excludePlane > 0) {
    // loop over all distance cuts
    for ( int help = 0; help < 10; help++) {
      streamlog_out ( MESSAGE2 ) << "DUT efficiency (distance < " << (help * 20 + 20) << " mu m :" << (double(_iHitDUTDistanceCut[help]) / _iHitDUT) << endl;
    } // end loop over all distance cuts
  }

}

void EUTelMultiLineFit::bookHistos() {


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  try {
    streamlog_out ( MESSAGE2 ) << "Booking histograms" << endl;

    const int    NBin = 20000;
    const double Chi2Min  = 0.0;
    const double Chi2Max  = 10000.0;
    const double Min  = -10000.0;
    const double Max  = 10000.0;
    const double angleMin  = -3.0;
    const double angleMax  = 3.0;
    const int    NBinHitDistance = 3000;
    const double hitDistanceMin = 0.0;
    const double hitDistanceMax = 3000.0;
    const double tracksMin = -0.5;
    const double tracksMax = 19.5;
    const double snrMin = 0.0;
    const double snrMax = 10000.0;

    const int NBin2D = 100;
    const double positionXMin = -10000.0;
    const double positionXMax = 10000.0;
    const double positionYMin = -10000.0;
    const double positionYMax = 10000.0;

    AIDA::IHistogram1D * numberTracksLocal =
      AIDAProcessor::histogramFactory(this)->createHistogram1D(_numberTracksLocalname,20,tracksMin,tracksMax);
    if ( numberTracksLocal ) {
      numberTracksLocal->setTitle("Number of tracks after #chi^{2} cut");
      _aidaHistoMap.insert( make_pair( _numberTracksLocalname, numberTracksLocal ) );
    } else {
      streamlog_out ( ERROR2 ) << "Problem booking the " << (_numberTracksLocalname) << endl;
      streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
      _histogramSwitch = false;
    }

    AIDA::IHistogram1D * chi2XLocal =
      AIDAProcessor::histogramFactory(this)->createHistogram1D(_chi2XLocalname,NBin,Chi2Min,Chi2Max);
    if ( chi2XLocal ) {
      chi2XLocal->setTitle("Chi2 X");
      _aidaHistoMap.insert( make_pair( _chi2XLocalname, chi2XLocal ) );
    } else {
      streamlog_out ( ERROR2 ) << "Problem booking the " << (_chi2XLocalname) << endl;
      streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
      _histogramSwitch = false;
    }

    AIDA::IHistogram1D * chi2YLocal =
      AIDAProcessor::histogramFactory(this)->createHistogram1D(_chi2YLocalname,NBin,Chi2Min,Chi2Max);
    if ( chi2YLocal ) {
      chi2YLocal->setTitle("Chi2 Y");
      _aidaHistoMap.insert( make_pair( _chi2YLocalname, chi2YLocal ) );
    } else {
      streamlog_out ( ERROR2 ) << "Problem booking the " << (_chi2YLocalname) << endl;
      streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
      _histogramSwitch = false;
    }

    AIDA::IHistogram1D * angleXLocal =
      AIDAProcessor::histogramFactory(this)->createHistogram1D(_angleXLocalname,NBin,angleMin,angleMax);
    if ( angleXLocal ) {
      angleXLocal->setTitle("X Angle");
      _aidaHistoMap.insert( make_pair( _angleXLocalname, angleXLocal ) );
    } else {
      streamlog_out ( ERROR2 ) << "Problem booking the " << (_angleXLocalname) << endl;
      streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
      _histogramSwitch = false;
    }

    AIDA::IHistogram1D * angleYLocal =
      AIDAProcessor::histogramFactory(this)->createHistogram1D(_angleYLocalname,NBin,angleMin,angleMax);
    if ( angleYLocal ) {
      angleYLocal->setTitle("Y Angle");
      _aidaHistoMap.insert( make_pair( _angleYLocalname, angleYLocal ) );
    } else {
      streamlog_out ( ERROR2 ) << "Problem booking the " << (_angleYLocalname) << endl;
      streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
      _histogramSwitch = false;
    }

    AIDA::IHistogram1D * hitDistanceLocal =
      AIDAProcessor::histogramFactory(this)->createHistogram1D(_hitDistanceLocalname,NBinHitDistance,hitDistanceMin,hitDistanceMax);
    if ( hitDistanceLocal ) {
      angleYLocal->setTitle("Distance of hit to predicted position");
      _aidaHistoMap.insert( make_pair( _hitDistanceLocalname, hitDistanceLocal ) );
    } else {
      streamlog_out ( ERROR2 ) << "Problem booking the " << (_hitDistanceLocalname) << endl;
      streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
      _histogramSwitch = false;
    }

    string tempHisto;
    string tempHistoName;
    string histoTitleXResid;
    string histoTitleYResid;
    string histoTitleSeedCharge;
    string histoTitleClusterCharge;
    string histoTitleXYPosition;

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
        AIDAProcessor::histogramFactory(this)->createHistogram1D(tempHistoName,NBin , Min,Max);
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

      {
        stringstream ss;
        stringstream pp;
        stringstream tt;

        pp << "SeedChargeLocal_d" << iDetector;
        tempHisto=pp.str();
        ss << _seedChargeLocalname << "_d" << iDetector;
        tempHistoName=ss.str();
        tt << "SeedCharge" << "_d" << iDetector;
        histoTitleSeedCharge=tt.str();

      }

      AIDA::IHistogram1D *  tempSeedChargeHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram1D(tempHistoName,NBin,snrMin,snrMax);
      if ( tempSeedChargeHisto ) {
        tempSeedChargeHisto->setTitle(histoTitleSeedCharge);
        _aidaHistoMap.insert( make_pair( tempHistoName, tempSeedChargeHisto ) );
      } else {
        streamlog_out ( ERROR2 ) << "Problem booking the " << (tempHistoName) << endl;
        streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }

      {
        stringstream ss;
        stringstream pp;
        stringstream tt;

        pp << "ClusterChargeLocal_d" << iDetector;
        tempHisto=pp.str();
        ss << _clusterChargeLocalname << "_d" << iDetector;
        tempHistoName=ss.str();
        tt << "ClusterCharge" << "_d" << iDetector;
        histoTitleClusterCharge=tt.str();

      }

      AIDA::IHistogram1D *  tempClusterChargeHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram1D(tempHistoName,NBin,snrMin,snrMax);
      if ( tempClusterChargeHisto ) {
        tempClusterChargeHisto->setTitle(histoTitleClusterCharge);
        _aidaHistoMap.insert( make_pair( tempHistoName, tempClusterChargeHisto ) );
      } else {
        streamlog_out ( ERROR2 ) << "Problem booking the " << (tempHistoName) << endl;
        streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }

      {
        stringstream ss;
        stringstream pp;
        stringstream tt;

        pp << "PositionXYLocal_d" << iDetector;
        tempHisto=pp.str();
        ss << _positionXYLocalname << "_d" << iDetector;
        tempHistoName=ss.str();
        tt << "PositionXY" << "_d" << iDetector;
        histoTitleXYPosition=tt.str();

      }

      AIDA::IHistogram2D * tempXYHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram2D(tempHistoName,NBin2D,positionXMin,positionXMax,NBin2D,positionYMin,positionYMax);
      if ( tempXYHisto ) {
        tempXYHisto->setTitle(histoTitleXYPosition);
        _aidaHistoMap.insert( make_pair( tempHistoName, tempXYHisto ) );
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


#endif
