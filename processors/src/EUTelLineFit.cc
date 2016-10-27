// Author Tatsiana Klimkovich, DESY <mailto:tklimk@mail.desy.de>
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
#include "EUTelLineFit.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
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
std::string EUTelLineFit::_chi2XLocalname          = "Chi2X";
std::string EUTelLineFit::_chi2YLocalname          = "Chi2Y";
std::string EUTelLineFit::_angleXLocalname          = "AngleX";
std::string EUTelLineFit::_angleYLocalname          = "AngleY";
std::string EUTelLineFit::_residualXLocalname      = "ResidualX";
std::string EUTelLineFit::_residualYLocalname      = "ResidualY";
#endif


EUTelLineFit::EUTelLineFit () : Processor("EUTelLineFit") {

  // modify processor description
  _description =
    "EUTelLineFit will fit a straight line";

  // input collection

  registerInputCollection(LCIO::TRACKERHIT,"HitCollectionName",
                          "Hit collection name",
                          _hitCollectionName, string ( "hit" ));

  // output collection

  registerOutputCollection(LCIO::TRACK,"OutputTrackCollectionName",
                           "Collection name for fitted tracks",
                           _outputTrackColName, string ("fittracks"));

  registerOutputCollection(LCIO::TRACKERHIT,"OutputHitCollectionName",
                           "Collection name for fitted particle hits (positions)",
                           _outputHitColName, string ("fithits"));

  // input parameters: take these from database later

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

}

void EUTelLineFit::init() {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

  // check if Marlin was built with GEAR support or not
#ifndef USE_GEAR

  message<ERROR5> ( "Marlin was not built with GEAR support." );
  message<ERROR5> ( "You need to install GEAR and recompile Marlin with -DUSE_GEAR before continue.");

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

  _siPlanesParameters  = const_cast<gear::SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

  _histogramSwitch = true;

#endif

  _nPlanes = _siPlanesParameters->getSiPlanesNumber();

  _xPos = new double[_nPlanes];
  _yPos = new double[_nPlanes];
  _zPos = new double[_nPlanes];
  _waferResidX = new double[_nPlanes];
  _waferResidY = new double[_nPlanes];
  _xFitPos = new double[_nPlanes];
  _yFitPos = new double[_nPlanes];

  _intrResolX = new double[_nPlanes];
  _intrResolY = new double[_nPlanes];

}

void EUTelLineFit::processRunHeader (LCRunHeader * rdr) {
  std::unique_ptr<EUTelRunHeaderImpl> header = std::make_unique<EUTelRunHeaderImpl>(rdr);
  header->addProcessor(type());

  // the run header contains the number of detectors. This number
  // should be in principle the same as the number of layers in the
  // geometry description
  if ( header->getNoOfDetector() != _siPlanesParameters->getSiPlanesNumber() ) {
    message<ERROR5> ( "Error during the geometry consistency check: " );
    message<ERROR5> ( log() << "The run header says there are " << header->getNoOfDetector() << " silicon detectors " );
    message<ERROR5> ( log() << "The GEAR description says     " << _siPlanesParameters->getSiPlanesNumber() << " silicon planes" );
    exit(-1);
  }

  // this is the right place also to check the geometry ID. This is a
  // unique number identifying each different geometry used at the
  // beam test. The same number should be saved in the run header and
  // in the xml file. If the numbers are different, instead of barely
  // quitting ask the user what to do.

  if ( header->getGeoID() != _siPlanesParameters->getSiPlanesID() ) {
    message<ERROR5> ( "Error during the geometry consistency check: " );
    message<ERROR5> ( log() << "The run header says the GeoID is " << header->getGeoID() );
    message<ERROR5> ( log() << "The GEAR description says is     " << _siPlanesParameters->getSiPlanesNumber() );
    string answer;
    while (true) {
      message<ERROR5> ( "Type Q to quit now or C to continue using the actual GEAR description anyway [Q/C]" );
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

  // now book histograms plz...
  if ( isFirstEvent() )  bookHistos();

  // increment the run counter
  ++_iRun;
}


void EUTelLineFit::processEvent (LCEvent * event) {


  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event) ;

  if ( evt->getEventType() == kEORE ) {
    message<DEBUG5> ( "EORE found: nothing else to do." );
    return;
  }


  try {

    LCCollectionVec * hitCollection     = static_cast<LCCollectionVec*> (event->getCollection( _hitCollectionName ));

    int detectorID    = -99; // it's a non sense
    int oldDetectorID = -100;
    int layerIndex;

    for ( int iHit = 0; iHit < hitCollection->getNumberOfElements(); iHit++ ) {

      TrackerHitImpl* hit = static_cast<TrackerHitImpl*> ( hitCollection->getElementAt(iHit) );
      UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder (hitCollection);
      detectorID = hitDecoder(hit)["sensorID"];

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

      // Getting positions of the hits.
      // Here the alignment constants are used to correct the positions.

      layerIndex   = _conversionIdMap[detectorID];

      // The other layers were aligned with respect to the first one.

      double off_x, off_y, theta_x, theta_y, theta_z;

      if (layerIndex == 0) {

        _xPos[iHit] = 1000 * hit->getPosition()[0]; // in um
        _yPos[iHit] = 1000 * hit->getPosition()[1]; // in um
        _zPos[iHit] = 1000 * hit->getPosition()[2]; // in um

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

        } else if (layerIndex == 3) {

          off_x = _alignmentConstantsFourthLayer[0];
          off_y = _alignmentConstantsFourthLayer[1];
          theta_x = _alignmentConstantsFourthLayer[2];
          theta_y = _alignmentConstantsFourthLayer[3];
          theta_z = _alignmentConstantsFourthLayer[4];

        } else if (layerIndex == 4) {

          off_x = _alignmentConstantsFifthLayer[0];
          off_y = _alignmentConstantsFifthLayer[1];
          theta_x = _alignmentConstantsFifthLayer[2];
          theta_y = _alignmentConstantsFifthLayer[3];
          theta_z = _alignmentConstantsFifthLayer[4];

        } else if (layerIndex == 5) {

          off_x = _alignmentConstantsSixthLayer[0];
          off_y = _alignmentConstantsSixthLayer[1];
          theta_x = _alignmentConstantsSixthLayer[2];
          theta_y = _alignmentConstantsSixthLayer[3];
          theta_z = _alignmentConstantsSixthLayer[4];

        } else {

          off_x = 0.0;
          off_y = 0.0;
          theta_x = 0.0;
          theta_y = 0.0;
          theta_z = 0.0;

        }

        // For documentation of these formulas look at EUTelAlign

        _xPos[iHit] = (cos(theta_y)*cos(theta_z)) * hit->getPosition()[0] * 1000 + ((-1)*sin(theta_x)*sin(theta_y)*cos(theta_z) + cos(theta_x)*sin(theta_z)) * hit->getPosition()[1] * 1000 + off_x;
        _yPos[iHit] = ((-1)*cos(theta_y)*sin(theta_z)) * hit->getPosition()[0] * 1000 + (sin(theta_x)*sin(theta_y)*sin(theta_z) + cos(theta_x)*cos(theta_z)) * hit->getPosition()[1] * 1000 + off_y;
        _zPos[iHit] = 1000 * hit->getPosition()[2];

      }
    }

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++ See Blobel Page 226 !!! +++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++

    int counter;

    float S1[2]   = {0,0};
    float Sx[2]   = {0,0};
    float Xbar[2] = {0,0};

    float * Zbar_X = new float[_nPlanes];
    float * Zbar_Y = new float[_nPlanes];
    for (counter = 0; counter < _nPlanes; counter++){
      Zbar_X[counter] = 0.;
      Zbar_Y[counter] = 0.;
    }

    float Sy[2]     = {0,0};
    float Ybar[2]   = {0,0};
    float Sxybar[2] = {0,0};
    float Sxxbar[2] = {0,0};
    float A2[2]     = {0,0};
    float Chiquare[2] = {0,0};
    float angle[2] = {0,0};

    // define S1
    for( counter = 0; counter < _nPlanes; counter++ ){
      S1[0] = S1[0] + 1/pow(_intrResolX[counter],2);
      S1[1] = S1[1] + 1/pow(_intrResolY[counter],2);
    }

    // define Sx
    for( counter = 0; counter < _nPlanes; counter++ ){
      Sx[0] = Sx[0] + _zPos[counter]/pow(_intrResolX[counter],2);
      Sx[1] = Sx[1] + _zPos[counter]/pow(_intrResolY[counter],2);
    }

    // define Xbar
    Xbar[0]=Sx[0]/S1[0];
    Xbar[1]=Sx[1]/S1[1];

    // coordinate transformation !! -> bar
    for( counter=0; counter < _nPlanes; counter++ ){
      Zbar_X[counter] = _zPos[counter]-Xbar[0];
      Zbar_Y[counter] = _zPos[counter]-Xbar[1];
    }

    // define Sy
    for( counter = 0; counter < _nPlanes; counter++ ){
      Sy[0] = Sy[0] + _xPos[counter]/pow(_intrResolX[counter],2);
      Sy[1] = Sy[1] + _yPos[counter]/pow(_intrResolY[counter],2);
    }

    // define Ybar
    Ybar[0]=Sy[0]/S1[0];
    Ybar[1]=Sy[1]/S1[1];

    // define Sxybar
    for( counter = 0; counter < _nPlanes; counter++ ){
      Sxybar[0] = Sxybar[0] + Zbar_X[counter] * _xPos[counter]/pow(_intrResolX[counter],2);
      Sxybar[1] = Sxybar[1] + Zbar_Y[counter] * _yPos[counter]/pow(_intrResolY[counter],2);
    }

    // define Sxxbar
    for( counter = 0; counter < _nPlanes; counter++ ){
      Sxxbar[0] = Sxxbar[0] + Zbar_X[counter] * Zbar_X[counter]/pow(_intrResolX[counter],2);
      Sxxbar[1] = Sxxbar[1] + Zbar_Y[counter] * Zbar_Y[counter]/pow(_intrResolX[counter],2);
    }

    // define A2

    A2[0]=Sxybar[0]/Sxxbar[0];
    A2[1]=Sxybar[1]/Sxxbar[1];

    // Calculate chi sqaured
    // Chi^2 for X and Y coordinate for hits in all planes

    for( counter = 0; counter < _nPlanes; counter++ ){
      Chiquare[0] += pow(-_zPos[counter]*A2[0]
                         +_xPos[counter]-Ybar[0]+Xbar[0]*A2[0],2)/pow(_intrResolX[counter],2);
      Chiquare[1] += pow(-_zPos[counter]*A2[1]
                         +_yPos[counter]-Ybar[1]+Xbar[1]*A2[1],2)/pow(_intrResolY[counter],2);
    }

    for( counter = 0; counter < _nPlanes; counter++ ){
      _waferResidX[counter] = (Ybar[0]-Xbar[0]*A2[0]+_zPos[counter]*A2[0])-_xPos[counter];
      _waferResidY[counter] = (Ybar[1]-Xbar[1]*A2[1]+_zPos[counter]*A2[1])-_yPos[counter];
    }

    // define angle


    angle[0] = atan(A2[0]);
    angle[1] = atan(A2[1]);


    // Define output track and hit collections
    LCCollectionVec     * fittrackvec = new LCCollectionVec(LCIO::TRACK);
    LCCollectionVec     * fitpointvec = new LCCollectionVec(LCIO::TRACKERHIT);

    // Set flag for storing track hits in track collection

    LCFlagImpl flag(fittrackvec->getFlag());
    flag.setBit( LCIO::TRBIT_HITS );
    fittrackvec->setFlag(flag.getFlag());


    // fill output collections...

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
    //  fittrack->setNdf(nBestFired); // Number of planes fired (!)

//    fittrack->setIsReferencePointPCA(false);

    // Calculate positions of fitted track in every plane

    for( counter = 0; counter < _nPlanes; counter++ ){

      _xFitPos[counter] = Ybar[0]-Xbar[0]*A2[0]+_zPos[counter]*A2[0];
      _yFitPos[counter] = Ybar[1]-Xbar[1]*A2[1]+_zPos[counter]*A2[1];

      TrackerHitImpl * fitpoint = new TrackerHitImpl;

      // Plane number stored as hit type
      fitpoint->setType(counter+1);
      double pos[3];
      pos[0] = _xFitPos[counter];
      pos[1] = _yFitPos[counter];
      pos[2] = _zPos[counter];

      fitpoint->setPosition(pos);

      // store fit point

      fitpointvec->push_back(fitpoint);

      //   add point to track

      fittrack->addHit(fitpoint);

    }

    fittrackvec->addElement(fittrack);

    event->addCollection(fittrackvec,_outputTrackColName);
    event->addCollection(fitpointvec,_outputHitColName);

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
        message<ERROR5> ( log() << "Not able to retrieve histogram pointer for " <<  _chi2XLocalname
                         << ".\nDisabling histogramming from now on " );
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
        message<ERROR5> ( log() << "Not able to retrieve histogram pointer for " <<  _chi2YLocalname
                         << ".\nDisabling histogramming from now on " );
        _histogramSwitch = false;
      }
    }

    if ( _histogramSwitch ) {
      {
        stringstream ss;
        ss << _angleXLocalname << endl;
      }
      if ( AIDA::IHistogram1D* anglex_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[_angleXLocalname]) )
        anglex_histo->fill(angle[0]*180/M_PI);
      else {
        message<ERROR5> ( log() << "Not able to retrieve histogram pointer for " <<  _angleXLocalname
                         << ".\nDisabling histogramming from now on " );
        _histogramSwitch = false;
      }
    }

    if ( _histogramSwitch ) {
      {
        stringstream ss;
        ss << _angleYLocalname << endl;
      }
      if ( AIDA::IHistogram1D* angley_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[_angleYLocalname]) )
        angley_histo->fill(angle[1]*180/M_PI);
      else {
        message<ERROR5> ( log() << "Not able to retrieve histogram pointer for " <<  _angleYLocalname
                         << ".\nDisabling histogramming from now on " );
        _histogramSwitch = false;
      }
    }


    for( int iDetector = 0; iDetector < _nPlanes; iDetector++ ){

      if ( _histogramSwitch ) {
        {
          stringstream ss;
          ss << _residualXLocalname << "-d" << iDetector;
          tempHistoName=ss.str();
        }
        if ( AIDA::IHistogram1D* residx_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName.c_str()]) )
          residx_histo->fill(_waferResidX[iDetector]);
        else {
          message<ERROR5> ( log() << "Not able to retrieve histogram pointer for " <<  _residualXLocalname
                           << ".\nDisabling histogramming from now on " );
          _histogramSwitch = false;
        }
      }

      if ( _histogramSwitch ) {
        {
          stringstream ss;
          ss << _residualYLocalname << "-d" << iDetector;
          tempHistoName=ss.str();
        }
        if ( AIDA::IHistogram1D* residy_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName.c_str()]) )
          residy_histo->fill(_waferResidY[iDetector]);
        else {
          message<ERROR5> ( log() << "Not able to retrieve histogram pointer for " <<  _residualYLocalname
                           << ".\nDisabling histogramming from now on " );
          _histogramSwitch = false;
        }
      }
    }

#endif

    delete [] Zbar_X;
    delete [] Zbar_Y;


  } catch (DataNotAvailableException& e) {

    streamlog_out  ( WARNING2 ) <<  "No input collection found on event " << event->getEventNumber()
                                << " in run " << event->getRunNumber() << endl;

  }

  ++_iEvt;

  if ( isFirstEvent() ) _isFirstEvent = false;

}

void EUTelLineFit::end() {

  delete [] _intrResolY;
  delete [] _intrResolX;
  delete [] _yFitPos;
  delete [] _xFitPos;
  delete [] _waferResidY;
  delete [] _waferResidX;
  delete [] _zPos;
  delete [] _yPos;
  delete [] _xPos;

  message<MESSAGE5> ( log() << "Successfully finished" ) ;

}

void EUTelLineFit::bookHistos() {


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  try {
    message<MESSAGE5> ( "Booking histograms" );

    const int    NBin = 100;
    const double Chi2Min  = 0.;
    const double Chi2Max  = 10.;
    const double Min  = -50.;
    const double Max  = 50.;
    const double angleMin  = -3.;
    const double angleMax  = 3.;


    AIDA::IHistogram1D * chi2XLocal =
      AIDAProcessor::histogramFactory(this)->createHistogram1D(_chi2XLocalname,NBin,Chi2Min,Chi2Max);
    if ( chi2XLocal ) {
      chi2XLocal->setTitle("Chi2 X");
      _aidaHistoMap.insert( make_pair( _chi2XLocalname, chi2XLocal ) );
    } else {
      message<ERROR5> ( log() << "Problem booking the " << (_chi2XLocalname) << ".\n"
                       << "Very likely a problem with path name. Switching off histogramming and continue w/o");
      _histogramSwitch = false;
    }

    AIDA::IHistogram1D * chi2YLocal =
      AIDAProcessor::histogramFactory(this)->createHistogram1D(_chi2YLocalname,NBin,Chi2Min,Chi2Max);
    if ( chi2YLocal ) {
      chi2YLocal->setTitle("Chi2 Y");
      _aidaHistoMap.insert( make_pair( _chi2YLocalname, chi2YLocal ) );
    } else {
      message<ERROR5> ( log() << "Problem booking the " << (_chi2YLocalname) << ".\n"
                       << "Very likely a problem with path name. Switching off histogramming and continue w/o");
      _histogramSwitch = false;
    }

    AIDA::IHistogram1D * angleXLocal =
      AIDAProcessor::histogramFactory(this)->createHistogram1D(_angleXLocalname,NBin,angleMin,angleMax);
    if ( angleXLocal ) {
      angleXLocal->setTitle("X Angle");
      _aidaHistoMap.insert( make_pair( _angleXLocalname, angleXLocal ) );
    } else {
      message<ERROR5> ( log() << "Problem booking the " << (_angleXLocalname) << ".\n"
                       << "Very likely a problem with path name. Switching off histogramming and continue w/o");
      _histogramSwitch = false;
    }

    AIDA::IHistogram1D * angleYLocal =
      AIDAProcessor::histogramFactory(this)->createHistogram1D(_angleYLocalname,NBin,angleMin,angleMax);
    if ( angleYLocal ) {
      angleYLocal->setTitle("Y Angle");
      _aidaHistoMap.insert( make_pair( _angleYLocalname, angleYLocal ) );
    } else {
      message<ERROR5> ( log() << "Problem booking the " << (_angleYLocalname) << ".\n"
                       << "Very likely a problem with path name. Switching off histogramming and continue w/o");
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

        pp << "ResidualXLocal-d" << iDetector;
        tempHisto=pp.str();
        ss << _residualXLocalname << "-d" << iDetector;
        tempHistoName=ss.str();
        tt << "XResidual" << "-d" << iDetector;
        histoTitleXResid=tt.str();

      }

      AIDA::IHistogram1D *  tempXHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram1D(tempHistoName,NBin, Min,Max);
      if ( tempXHisto ) {
        tempXHisto->setTitle(histoTitleXResid);
        _aidaHistoMap.insert( make_pair( tempHistoName, tempXHisto ) );
      } else {
        message<ERROR5> ( log() << "Problem booking the " << (tempHistoName) << ".\n"
                         << "Very likely a problem with path name. Switching off histogramming and continue w/o");
        _histogramSwitch = false;
      }

      {
        stringstream ss;
        stringstream pp;
        stringstream tt;

        pp << "ResidualYLocal-d" << iDetector;
        tempHisto=pp.str();
        ss << _residualYLocalname << "-d" << iDetector;
        tempHistoName=ss.str();
        tt << "YResidual" << "-d" << iDetector;
        histoTitleYResid=tt.str();
      }

      AIDA::IHistogram1D *  tempYHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram1D(tempHistoName,NBin, Min,Max);
      if ( tempYHisto ) {
        tempYHisto->setTitle(histoTitleYResid);
        _aidaHistoMap.insert( make_pair( tempHistoName, tempYHisto ) );
      } else {
        message<ERROR5> ( log() << "Problem booking the " << (tempHistoName) << ".\n"
                         << "Very likely a problem with path name. Switching off histogramming and continue w/o");
        _histogramSwitch = false;
      }

    }

  } catch (lcio::Exception& e ) {

    message<ERROR5> ( log() << "No AIDAProcessor initialized. Type q to exit or c to continue without histogramming" );
    string answer;
    while ( true ) {
      message<ERROR5> ( "[q]/[c]" );
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
