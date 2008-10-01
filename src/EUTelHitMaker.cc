// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelHitMaker.cc,v 1.24 2008-10-01 10:21:47 bulgheroni Exp $
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
#include "EUTelHitMaker.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelSparseCluster2Impl.h"
#include "EUTelExceptions.h"
#include "EUTelEtaFunctionImpl.h"

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
#include <AIDA/IHistogram2D.h>
#include <AIDA/IHistogram3D.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <UTIL/CellIDDecoder.h>

// system includes <>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace marlin;
using namespace gear;
using namespace eutelescope;

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelHitMaker::_hitHistoLocalName          = "HitHistoLocal";
std::string EUTelHitMaker::_hitHistoTelescopeName      = "HitHistoTelescope";
std::string EUTelHitMaker::_densityPlotName            = "DensityPlot";
std::string EUTelHitMaker::_clusterCenterHistoName     = "ClusterCenter";
std::string EUTelHitMaker::_clusterCenterXHistoName    = "ClusterCenterX";
std::string EUTelHitMaker::_clusterCenterYHistoName    = "ClusterCenterY";
std::string EUTelHitMaker::_clusterCenterEtaHistoName  = "ClusterCenterEta";
std::string EUTelHitMaker::_clusterCenterEtaXHistoName = "ClusterCenterEtaX";
std::string EUTelHitMaker::_clusterCenterEtaYHistoName = "ClusterCenterEtaY";
#endif

EUTelHitMaker::EUTelHitMaker () : Processor("EUTelHitMaker") {

  // modify processor description
  _description =
    "EUTelHitMaker is responsible to translate cluster centers from the local frame of reference"
    " to the external frame of reference using the GEAR geometry description";

  registerInputCollection(LCIO::TRACKERPULSE,"PulseCollectionName",
                          "Cluster (pulse) collection name",
                          _pulseCollectionName, string ( "cluster" ));

  registerOutputCollection(LCIO::TRACKERHIT,"HitCollectionName",
                           "Hit collection name",
                           _hitCollectionName, string ( "hit" ));


  registerProcessorParameter("CoGAlgorithm", "Select here how the center of gravity should be calculated.\n"
                             "FULL: using the full cluster\n"
                             "NPixel: using only the first N most significant pixels (set NPixel too)\n"
                             "NxMPixel: using a subframe of the cluster N x M pixels wide (set NxMPixel too).",
                             _cogAlgorithm, string( "NxMPixel" ) );

  registerOptionalParameter("NPixel", "The number of most significant pixels to be used if CoGAlgorithm is \"NPixel\"",
                            _nPixel, static_cast<int>( 9 ) );

  registerOptionalParameter("Enable3DHisto","If true a 3D histo will be filled. It may require large memory",
			    _3DHistoSwitch, static_cast<bool> ( true ) );


  vector<int > xyCluSizeExample(2,3);
  registerOptionalParameter("NxMPixel", "The submatrix size to be used for CoGAlgorithm = \"NxMPixel\"",
                            _xyCluSize, xyCluSizeExample, xyCluSizeExample.size());

  vector<string > etaNames;
  etaNames.push_back("xEta");
  etaNames.push_back("yEta");

  registerOptionalParameter("EtaCollectionName",
                            "The name of the collections containing the eta function (x and y respectively)",
                            _etaCollectionNames, etaNames, etaNames.size());


  registerOptionalParameter("EtaSwitch","Enable or disable eta correction",
                            _etaCorrection, static_cast< bool > ( 1 ) );

}


void EUTelHitMaker::init() {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

  // transform the algorithm string to small letters
  transform(_cogAlgorithm.begin(), _cogAlgorithm.end(), _cogAlgorithm.begin(), ::tolower);

  // check if Marlin was built with GEAR support or not
#ifndef USE_GEAR

  streamlog_out ( ERROR4 ) <<  "Marlin was not built with GEAR support." << endl
                           <<  "You need to install GEAR and recompile Marlin with -DUSE_GEAR before continue." << endl;

  // I'm thinking if this is the case of throwing an exception or
  // not. This is a really error and not something that can
  // exceptionally happens. Still not sure what to do
  exit(-1);

#else

  // check if the GEAR manager pointer is not null!
  if ( Global::GEAR == 0x0 ) {
    streamlog_out ( ERROR4 ) <<  "The GearMgr is not available, for an unknown reason." << endl;
    exit(-1);
  }

  _siPlanesParameters  = const_cast<SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

  _histogramSwitch = true;

#endif

}

void EUTelHitMaker::processRunHeader (LCRunHeader * rdr) {


  auto_ptr<EUTelRunHeaderImpl> header ( new EUTelRunHeaderImpl (rdr) );
  header->addProcessor( type() );

  // the run header contains the number of detectors. This number
  // should be in principle be the same as the number of layers in the
  // geometry description
  if ( header->getNoOfDetector() != _siPlanesParameters->getSiPlanesNumber() ) {
    streamlog_out ( ERROR4 ) << "Error during the geometry consistency check: " << endl
                             << "The run header says there are " << header->getNoOfDetector() << " silicon detectors " << endl
                             << "The GEAR description says     " << _siPlanesParameters->getSiPlanesNumber() << " silicon planes" << endl;
    exit(-1);
  }

  // this is the right place also to check the geometry ID. This is a
  // unique number identifying each different geometry used at the
  // beam test. The same number should be saved in the run header and
  // in the xml file. If the numbers are different, instead of barely
  // quitting ask the user what to do.

  if ( header->getGeoID() == 0 )
    streamlog_out ( WARNING0 ) <<  "The geometry ID in the run header is set to zero." << endl
                               <<  "This may mean that the GeoID parameter was not set" << endl;


  if ( header->getGeoID() != _siPlanesParameters->getSiPlanesID() ) {
    streamlog_out ( ERROR1 ) <<  "Error during the geometry consistency check: " << endl
                             << "The run header says the GeoID is " << header->getGeoID() << endl
                             << "The GEAR description says is     " << _siPlanesParameters->getSiPlanesID() << endl;
    string answer;
    while (true) {
      streamlog_out ( ERROR1 ) << "Type Q to quit now or C to continue using the actual GEAR description anyway [Q/C]" << endl;
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


void EUTelHitMaker::processEvent (LCEvent * event) {

  if (_iEvt % 10 == 0)
    streamlog_out( MESSAGE4 ) << "Processing event "
                              << setw(6) << setiosflags(ios::right) << event->getEventNumber() << " in run "
                              << setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber() << setfill(' ')
                              << " (Total = " << setw(10) << _iEvt << ")" << resetiosflags(ios::left) << endl;
  ++_iEvt;


  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event) ;

  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return;
  } else if ( evt->getEventType() == kUNKNOWN ) {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }


  try {

    LCCollectionVec * pulseCollection   = static_cast<LCCollectionVec*> (event->getCollection( _pulseCollectionName ));
    LCCollectionVec * hitCollection     = new LCCollectionVec(LCIO::TRACKERHIT);
    LCCollectionVec * xEtaCollection = 0x0, * yEtaCollection = 0x0;

    CellIDDecoder<TrackerPulseImpl>  pulseCellDecoder(pulseCollection);

    if ( _etaCorrection == 1 ) {
      // this means that the user wants to apply the eta correction to
      // the center of gravity, so we need to have the two Eta function
      // collections

      try {
        xEtaCollection = static_cast<LCCollectionVec*> (event->getCollection( _etaCollectionNames[0] )) ;
      } catch (DataNotAvailableException& e) {
        streamlog_out ( ERROR1 )  << "The eta collection " << _etaCollectionNames[0] << " is not available" << endl
                                  << "Continuing without eta correction " << endl;
        _etaCorrection = 0;
      }

      try {
        yEtaCollection = static_cast<LCCollectionVec*> (event->getCollection( _etaCollectionNames[1] )) ;
      } catch (DataNotAvailableException& e) {
        streamlog_out ( ERROR1 ) << "The eta collection " << _etaCollectionNames[1] << " is not available" << endl
                                 << "Continuing without eta correction " << endl;
        _etaCorrection = 0;
      }
    }

    if ( isFirstEvent() && _etaCorrection == 1) {
      if ( ( xEtaCollection->getNumberOfElements() != _siPlanesParameters->getSiPlanesNumber() ) ||
           ( yEtaCollection->getNumberOfElements() != _siPlanesParameters->getSiPlanesNumber() ) ) {
        streamlog_out ( ERROR1 ) <<  "The eta collections contain a different number of elements wrt to "
                                 << _siPlanesParameters->getSiPlanesNumber() << endl
                                 << "Continuing without eta correction " << endl;
        _etaCorrection = 0;
      } else {
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
        // this is also the right place to book the eta specific
        // histograms.
        if ( _histogramSwitch ) {
          for ( int iDet = 0 ; iDet < _siPlanesParameters->getSiPlanesNumber(); iDet++) {
            string basePath;
            {
              stringstream ss ;
              ss << "plane-" << iDet << "/";
              basePath = ss.str();
            }

            EUTelEtaFunctionImpl * xEtaFunc = static_cast<EUTelEtaFunctionImpl*> ( xEtaCollection->getElementAt(iDet) );
            EUTelEtaFunctionImpl * yEtaFunc = static_cast<EUTelEtaFunctionImpl*> ( yEtaCollection->getElementAt(iDet) );
            int xNoOfBin = xEtaFunc->getNoOfBin();
            int yNoOfBin = yEtaFunc->getNoOfBin();
            string tempHistoName;
            {
              stringstream ss;
              ss << _clusterCenterEtaHistoName << "-" << iDet;
              tempHistoName = ss.str();
            }
            AIDA::IHistogram2D * clusterCenterEta =
              AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName ).c_str(),
                                                                        1* xNoOfBin, -.5, +.5, 1 * yNoOfBin, -.5, +.5);
            if ( clusterCenterEta ) {
              clusterCenterEta->setTitle("Position of the cluster center (Eta corrected)");
              _aidaHistoMap.insert( make_pair( tempHistoName, clusterCenterEta ) );
            } else {
              streamlog_out ( ERROR1 )  << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                        << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
              _histogramSwitch = false;
            }


            {
              stringstream ss;
              ss << _clusterCenterHistoName << "-" << iDet;
              tempHistoName = ss.str();
            }
            AIDA::IHistogram2D * clusterCenter =
              AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName ).c_str(),
                                                                        1* xNoOfBin, -.5, +.5, 1 * yNoOfBin, -.5, +.5);
            if ( clusterCenter ) {
              clusterCenterEta->setTitle("Position of the cluster center");
              _aidaHistoMap.insert( make_pair( tempHistoName, clusterCenter ) );
            } else {
              streamlog_out ( ERROR1 )  << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                        << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
              _histogramSwitch = false;
            }


            {
              stringstream ss;
              ss << _clusterCenterEtaXHistoName << "-" << iDet;
              tempHistoName = ss.str();
            }
            AIDA::IHistogram1D * clusterCenterEtaX =
              AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                        1 * xNoOfBin, -0.5, +0.5 );
            if ( clusterCenterEtaX ) {
              clusterCenterEtaX->setTitle("Projection along X of the cluster center (Eta corrected)");
              _aidaHistoMap.insert( make_pair( tempHistoName, clusterCenterEtaX ) );
            } else {
              streamlog_out ( ERROR1 )  << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                        << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
              _histogramSwitch = false;
            }

            {
              stringstream ss;
              ss << _clusterCenterXHistoName << "-" << iDet;
              tempHistoName = ss.str();
            }
            AIDA::IHistogram1D * clusterCenterX =
              AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                        1 * xNoOfBin, -0.5, +0.5 );
            if ( clusterCenterX ) {
              clusterCenterX->setTitle("Projection along X of the cluster center");
              _aidaHistoMap.insert( make_pair( tempHistoName, clusterCenterX ) );
            } else {
              streamlog_out ( ERROR1 )  << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                        << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
              _histogramSwitch = false;
            }

            {
              stringstream ss;
              ss << _clusterCenterEtaYHistoName << "-" << iDet;
              tempHistoName = ss.str();
            }
            AIDA::IHistogram1D * clusterCenterEtaY =
              AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                        1 * xNoOfBin, -0.5, +0.5 );
            if ( clusterCenterEtaY ) {
              clusterCenterEtaY->setTitle("Projection along Y of the cluster center (Eta corrected)");
              _aidaHistoMap.insert( make_pair( tempHistoName, clusterCenterEtaY ) );
            } else {
              streamlog_out ( ERROR1 )  << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                        << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
              _histogramSwitch = false;
            }

            {
              stringstream ss;
              ss << _clusterCenterYHistoName << "-" << iDet;
              tempHistoName = ss.str();
            }
            AIDA::IHistogram1D * clusterCenterY =
              AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                        1 * xNoOfBin, -0.5, +0.5 );
            if ( clusterCenterY ) {
              clusterCenterY->setTitle("Projection along Y of the cluster center");
              _aidaHistoMap.insert( make_pair( tempHistoName, clusterCenterY ) );
            } else {
              streamlog_out ( ERROR1 )  << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                        << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
              _histogramSwitch = false;
            }

          }
        }
#endif

      }
    }

    int detectorID    = -99; // it's a non sense
    int oldDetectorID = -100;

    int    layerIndex = -99;
    double xZero = 0., yZero = 0., zZero = 0. ;
    double xSize = 0., ySize = 0.;
    double zThickness = 0.;
    double xPitch = 0., yPitch = 0.;
    double xPointing[2] = { 1., 0. }, yPointing[2] = { 1., 0. };

    for ( int iPulse = 0; iPulse < pulseCollection->getNumberOfElements(); iPulse++ ) {

      TrackerPulseImpl     * pulse   = static_cast<TrackerPulseImpl*> ( pulseCollection->getElementAt(iPulse) );
      EUTelVirtualCluster  * cluster;
      ClusterType type = static_cast<ClusterType>(static_cast<int>((pulseCellDecoder(pulse)["type"])));

      if ( type == kEUTelFFClusterImpl ) {

        // fixed cluster implementation. Remember it can come from
        // both RAW and ZS data
        cluster = new EUTelFFClusterImpl( static_cast<TrackerDataImpl*> (pulse->getTrackerData()) );

      } else if ( type == kEUTelSparseClusterImpl ) {

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
            ( static_cast<TrackerDataImpl *> ( pulse->getTrackerData()  ) );
        } else {
          streamlog_out ( ERROR4 ) << "Unknown pixel type.  Sorry for quitting." << endl;
          throw UnknownDataTypeException("Pixel type unknown");
        }

      } else {
        streamlog_out ( ERROR4 ) <<  "Unknown cluster type. Sorry for quitting" << endl;
        throw UnknownDataTypeException("Cluster type unknown");
      }

      // there could be several clusters belonging to the same
      // detector. So update the geometry information only if this new
      // cluster belongs to a different detector.
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

        // perfect! The full geometry description is now coming from the
        // GEAR interface. Let's keep the finger xed!
        layerIndex   = _conversionIdMap[detectorID];
        xZero        = _siPlanesLayerLayout->getSensitivePositionX(layerIndex); // mm
        yZero        = _siPlanesLayerLayout->getSensitivePositionY(layerIndex); // mm
        zZero        = _siPlanesLayerLayout->getSensitivePositionZ(layerIndex); // mm
        zThickness   = _siPlanesLayerLayout->getSensitiveThickness(layerIndex); // mm
        xPitch       = _siPlanesLayerLayout->getSensitivePitchX(layerIndex);    // mm
        yPitch       = _siPlanesLayerLayout->getSensitivePitchY(layerIndex);    // mm
        xSize        = _siPlanesLayerLayout->getSensitiveSizeX(layerIndex);     // mm
        ySize        = _siPlanesLayerLayout->getSensitiveSizeY(layerIndex);     // mm
        xPointing[0] = _siPlanesLayerLayout->getSensitiveRotation1(layerIndex); // was -1 ;
        xPointing[1] = _siPlanesLayerLayout->getSensitiveRotation2(layerIndex); // was  0 ;
        yPointing[0] = _siPlanesLayerLayout->getSensitiveRotation3(layerIndex); // was  0 ;
        yPointing[1] = _siPlanesLayerLayout->getSensitiveRotation4(layerIndex); // was -1 ;

        if (  ( xPointing[0] == xPointing[1] ) && ( xPointing[0] == 0 ) ) {
          streamlog_out ( ERROR4 ) << "Detector " << detectorID << " has a singular rotation matrix. Sorry for quitting" << endl;
        }

        if (  ( yPointing[0] == yPointing[1] ) && ( yPointing[0] == 0 ) ) {
          streamlog_out ( ERROR4 ) << "Detector " << detectorID << " has a singular rotation matrix. Sorry for quitting" << endl;
        }

      }

      // get the position of the seed pixel. This is in pixel number.
      int xCluCenter, yCluCenter;
      cluster->getCenterCoord(xCluCenter, yCluCenter);

      // with the charge center of gravity calculation, we get a shift
      // from the seed pixel center due to the charge distribution. Those
      // two numbers are the correction values in the case the Eta
      // correction is not applied.
      float xShift, yShift;
      if ( _cogAlgorithm == "full" ) {
        cluster->getCenterOfGravityShift( xShift, yShift );
      } else if ( _cogAlgorithm == "npixel" ) {
        cluster->getCenterOfGravityShift( xShift, yShift, _nPixel );
      } else if ( _cogAlgorithm == "nxmpixel") {
        cluster->getCenterOfGravityShift( xShift, yShift, _xyCluSize[0], _xyCluSize[1]);
      }

      double xCorrection = static_cast<double> (xShift) ;
      double yCorrection = static_cast<double> (yShift) ;


      if ( _etaCorrection == 1 ) {

        EUTelEtaFunctionImpl * xEtaFunc = static_cast<EUTelEtaFunctionImpl*> ( xEtaCollection->getElementAt(detectorID) );
        EUTelEtaFunctionImpl * yEtaFunc = static_cast<EUTelEtaFunctionImpl*> ( yEtaCollection->getElementAt(detectorID) );

        bool anomalous = false;
        if ( ( xShift >= -0.5 ) && ( xShift <= 0.5 ) ) {
          xCorrection = xEtaFunc->getEtaFromCoG( xShift );
        } else {
          anomalous = true;
        }
        if ( ( yShift >= -0.5 ) && ( yShift <= 0.5 ) ) {
          yCorrection = yEtaFunc->getEtaFromCoG( yShift );
        } else {
          anomalous = true;
        }

        if ( anomalous )  {
          streamlog_out ( DEBUG2 ) << "Found anomalous cluster\n" << ( * cluster ) << endl;
        }

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
        string tempHistoName;
        if ( _histogramSwitch ) {
          {
            stringstream ss;
            ss  << _clusterCenterEtaHistoName << "-" << detectorID ;
            tempHistoName = ss.str();
          }
          if ( AIDA::IHistogram2D * histo = dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[ tempHistoName ] )) {
            histo->fill( xCorrection, yCorrection );
          }

          {
            stringstream ss;
            ss  << _clusterCenterHistoName << "-" << detectorID ;
            tempHistoName = ss.str();
          }
          if ( AIDA::IHistogram2D * histo = dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[ tempHistoName ] )) {
            histo->fill( xShift, yShift );
          }

          {
            stringstream ss;
            ss << _clusterCenterEtaXHistoName << "-" << detectorID ;
            tempHistoName = ss.str();
          }
          if ( AIDA::IHistogram1D * histo = dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[ tempHistoName ] )) {
            histo->fill( xCorrection );
          }

          {
            stringstream ss;
            ss << _clusterCenterXHistoName << "-" << detectorID ;
            tempHistoName = ss.str();
          }
          if ( AIDA::IHistogram1D * histo = dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[ tempHistoName ] )) {
            histo->fill( xShift );
          }

          {
            stringstream ss;
            ss << _clusterCenterEtaYHistoName << "-" << detectorID ;
            tempHistoName = ss.str();
          }
          if ( AIDA::IHistogram1D * histo = dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[ tempHistoName ] )) {
            histo->fill( yCorrection );
          }

          {
            stringstream ss;
            ss << _clusterCenterYHistoName << "-" << detectorID ;
            tempHistoName = ss.str();
          }
          if ( AIDA::IHistogram1D * histo = dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[ tempHistoName ] )) {
            histo->fill( yShift );
          }
        }
#endif

      }

      // rescale the pixel number in millimeter
      double xDet = ( static_cast<double> (xCluCenter) + xCorrection + 0.5 ) * xPitch ;
      double yDet = ( static_cast<double> (yCluCenter) + yCorrection + 0.5 ) * yPitch ;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
      string tempHistoName;
      if ( _histogramSwitch ) {
        {
          stringstream ss;
          ss << _hitHistoLocalName << "-" << detectorID ;
          tempHistoName = ss.str();
        }
        if ( AIDA::IHistogram2D* histo = dynamic_cast<AIDA::IHistogram2D*>(_aidaHistoMap[ tempHistoName ]) )
          histo->fill(xDet, yDet);
        else {
          streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                    << ".\nDisabling histogramming from now on " << endl;
          _histogramSwitch = false;
        }

      }
#endif


      // now perform the rotation of the frame of references and put the
      // results already into a 3D array of double to be ready for the
      // setPosition method of TrackerHit
      double telPos[3];
      telPos[0] = xPointing[0] * xDet + xPointing[1] * yDet;
      telPos[1] = yPointing[0] * xDet + yPointing[1] * yDet;

      // now the translation
      // not sure about the sign. At least it is working for the current
      // configuration but we need to double check it
      double sign = 0;
      if      ( xPointing[0] < 0 )       sign = -1 ;
      else if ( xPointing[0] > 0 )       sign =  1 ;
      else {
        if       ( xPointing[1] < 0 )    sign = -1 ;
        else if  ( xPointing[1] > 0 )    sign =  1 ;
      }
      telPos[0] -= sign * ( xZero + xSize/2 );

      if      ( yPointing[0] < 0 )       sign = -1 ;
      else if ( yPointing[0] > 0 )       sign =  1 ;
      else {
        if       ( yPointing[1] < 0 )    sign = -1 ;
        else if  ( yPointing[1] > 0 )    sign =  1 ;
      }
      telPos[1] -= sign * ( yZero + ySize/2 );
      telPos[2] = zZero + 0.5 * zThickness;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
      if ( _histogramSwitch ) {
        {
          stringstream ss;
          ss << _hitHistoTelescopeName << "-" << detectorID ;
          tempHistoName = ss.str();
        }
        AIDA::IHistogram2D * histo2D = dynamic_cast<AIDA::IHistogram2D*> (_aidaHistoMap[ tempHistoName ] );
        if ( histo2D ) histo2D->fill( telPos[0], telPos[1] );
        else {
          streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                    << ".\nDisabling histogramming from now on " << endl;
          _histogramSwitch = false;
        }

	if ( _3DHistoSwitch ) {
	  AIDA::IHistogram3D * histo3D = dynamic_cast<AIDA::IHistogram3D*> (_aidaHistoMap[ _densityPlotName ] );
	  if ( histo3D ) histo3D->fill( telPos[0], telPos[1], telPos[2] );
	  else {
	    streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
				      << ".\nDisabling histogramming from now on " << endl;
	    _histogramSwitch = false;
	  }
	}
      }
#endif

      // create the new hit
      TrackerHitImpl * hit = new TrackerHitImpl;
      hit->setPosition( &telPos[0] );
      hit->setType( pulseCellDecoder(pulse)["type"] );

      // prepare a LCObjectVec to store the current cluster
      LCObjectVec clusterVec;
      clusterVec.push_back( pulse->getTrackerData() );

      // add the clusterVec to the hit
      hit->rawHits() = clusterVec;

      // add the new hit to the hit collection
      hitCollection->push_back( hit );

      // delete the eutel cluster
      delete cluster;

    }
    evt->addCollection( hitCollection, _hitCollectionName );

    if ( isFirstEvent() ) _isFirstEvent = false;
  } catch (DataNotAvailableException& e  ) {
    streamlog_out  ( WARNING2 ) <<  "No input collection found on event " << event->getEventNumber()
                                << " in run " << event->getRunNumber() << endl;
  }

}

void EUTelHitMaker::end() {

  streamlog_out ( MESSAGE4 )  << "Successfully finished" << endl;
}

void EUTelHitMaker::bookHistos() {

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  try {
    streamlog_out ( MESSAGE4 ) <<  "Booking histograms" << endl;

    string tempHistoName;

    // histograms are grouped into folders named after the
    // detector. This requires to loop on detector now.
    for (int iDet = 0 ; iDet < _siPlanesParameters->getSiPlanesNumber(); iDet++) {

      string basePath;
      {
        stringstream ss ;
        ss << "plane-" << iDet;
        basePath = ss.str();
      }
      AIDAProcessor::tree(this)->mkdir(basePath.c_str());
      basePath = basePath + "/";

      {
        stringstream ss ;
        ss <<  _hitHistoLocalName << "-" << iDet ;
        tempHistoName = ss.str();
      }



      double xMin =  0;
      double xMax =  _siPlanesLayerLayout->getSensitiveSizeX ( iDet );

      double yMin =  0;
      double yMax =  _siPlanesLayerLayout->getSensitiveSizeY ( iDet );

      int xNBin =  _siPlanesLayerLayout->getSensitiveNpixelX( iDet );
      int yNBin =  _siPlanesLayerLayout->getSensitiveNpixelY( iDet );


      AIDA::IHistogram2D * hitHistoLocal = AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(),
                                                                                                     xNBin, xMin, xMax, yNBin, yMin, yMax );
      if ( hitHistoLocal ) {
        hitHistoLocal->setTitle("Hit map in the detector local frame of reference");
        _aidaHistoMap.insert( make_pair( tempHistoName, hitHistoLocal ) );
      } else {
        streamlog_out ( ERROR1 )  << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                  << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }

      // 2 should be enough because it
      // means that the sensor is wrong
      // by all its size.
      double safetyFactor = 2.0;

      xMin = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionX( iDet ) -
                                     ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeX ( iDet ) ));
      xMax = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionX( iDet ) +
                                     ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeX ( iDet )));

      yMin = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionY( iDet ) -
                                     ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeY ( iDet )));
      yMax = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionY( iDet ) +
                                     ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeY ( iDet )) );

      xNBin = static_cast< int > ( safetyFactor ) * _siPlanesLayerLayout->getSensitiveNpixelX( iDet );
      yNBin = static_cast< int > ( safetyFactor ) * _siPlanesLayerLayout->getSensitiveNpixelY( iDet );

      {
        stringstream ss ;
        ss <<  _hitHistoTelescopeName << "-" << iDet ;
        tempHistoName = ss.str();
      }
      AIDA::IHistogram2D * hitHistoTelescope =
        AIDAProcessor::histogramFactory(this)->createHistogram2D( ( basePath + tempHistoName ).c_str(),
                                                                  xNBin, xMin, xMax, yNBin, yMin, yMax );

      if ( hitHistoTelescope ) {
        hitHistoTelescope->setTitle("Hit map in the telescope frame of reference");
        _aidaHistoMap.insert( make_pair ( tempHistoName, hitHistoTelescope ) );
      } else {
        streamlog_out ( ERROR1 )  << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                  << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }
    }

    // we have to found the boundaries of this histograms. Let's take
    // the outer positions in all directions
    double xMin  =      numeric_limits< double >::max();
    double xMax  = -1 * numeric_limits< double >::max();
    int    xNBin = numeric_limits< int >::min();

    double yMin  =      numeric_limits< double >::max();
    double yMax  = -1 * numeric_limits< double >::max();
    int    yNBin = numeric_limits< int >::min();

    for ( int iPlane = 0 ; iPlane < _siPlanesParameters->getSiPlanesNumber(); ++iPlane ) {

      // x axis
      xMin  = min( _siPlanesLayerLayout->getSensitivePositionX( iPlane ) - ( 0.5  * _siPlanesLayerLayout->getSensitiveSizeX( iPlane )), xMin);
      xMax  = max( _siPlanesLayerLayout->getSensitivePositionX( iPlane ) + ( 0.5  * _siPlanesLayerLayout->getSensitiveSizeX( iPlane )), xMax);
      xNBin = max( _siPlanesLayerLayout->getSensitiveNpixelX( iPlane ), xNBin );

      // y axis
      yMin  = min( _siPlanesLayerLayout->getSensitivePositionY( iPlane ) - ( 0.5  * _siPlanesLayerLayout->getSensitiveSizeY( iPlane )), yMin);
      yMax  = max( _siPlanesLayerLayout->getSensitivePositionY( iPlane ) + ( 0.5  * _siPlanesLayerLayout->getSensitiveSizeY( iPlane )), yMax);
      yNBin = max( _siPlanesLayerLayout->getSensitiveNpixelY( iPlane ), yNBin );

    }

    if ( _3DHistoSwitch ) {
      // since we may still have alignment problem, we have to take a
      // safety factor on the x and y direction especially.
      // here I take something less than 2 because otherwise I will have
      // a 200MB histogram.
      double safetyFactor = 1.2;
      
      double xDistance = std::abs( xMax - xMin ) ;
      double xCenter   = ( xMax + xMin ) / 2 ;
      xMin  = xCenter - safetyFactor * ( xDistance / 2 );
      xMax  = xCenter + safetyFactor * ( xDistance / 2 );
      xNBin = static_cast< int > ( xNBin * safetyFactor );
      
      // generate the x axis binning 
      vector< double > xAxis;
      double step = xDistance / xNBin;
      for ( int i = 0 ; i < xNBin ; ++i ) {
	xAxis.push_back ( xMin + i * step );
      }
      
      double yDistance = std::abs( yMax - yMin ) ;
      double yCenter   = ( yMax + yMin ) / 2 ;
      yMin  = yCenter - safetyFactor * ( yDistance / 2 );
      yMax  = yCenter + safetyFactor * ( yDistance / 2 );
      yNBin = static_cast< int > ( yNBin * safetyFactor );
      
      // generate the y axis binning
      vector< double > yAxis;
      step = yDistance / yNBin;
      for ( int i = 0 ; i < yNBin ; ++i ) {
	yAxis.push_back( yMin + i * step ) ;
      }
      
      
      // generate the z axis but not equally spaced!
      double safetyMargin = 10; // this is mm
      vector< double > zAxis;
      for ( int i = 0 ; i < 2 * _siPlanesParameters->getSiPlanesNumber(); ++i ) {
	double zPos =  _siPlanesLayerLayout->getSensitivePositionZ( i/2 );
	zAxis.push_back( zPos - safetyMargin) ;
	++i;
	zAxis.push_back( zPos + safetyMargin );
      }
      
      
      AIDA::IHistogram3D * densityPlot = AIDAProcessor::histogramFactory(this)->createHistogram3D( _densityPlotName , 
												   "Hit position in the telescope frame of reference",
												   xAxis, yAxis, zAxis, "");
      
      if ( densityPlot ) {
	_aidaHistoMap.insert( make_pair ( _densityPlotName, densityPlot ) ) ;
      } else {
	streamlog_out ( ERROR1 )  << "Problem booking the " << (_densityPlotName) << ".\n"
				  << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
	_histogramSwitch = false;
      }
    }

  } catch (lcio::Exception& e ) {

    streamlog_out ( ERROR1 ) << "No AIDAProcessor initialized. Type q to exit or c to continue without histogramming" << endl;
    string answer;
    while ( true ) {
      streamlog_out ( ERROR1 ) <<  "[q]/[c]" << endl;
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
