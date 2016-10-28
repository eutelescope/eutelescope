// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifdef USE_GEAR

// eutelescope includes ".h"
#include "EUTelExampleProcessorCorrelator.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"

// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ITree.h>
#endif

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <UTIL/CellIDDecoder.h>

// system includes <>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <cstdlib>

using namespace std;
using namespace marlin;
using namespace eutelescope;

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelExampleProcessorCorrelator::_hitXCorrelationHistoName       = "HitXCorrelatioHisto";
std::string EUTelExampleProcessorCorrelator::_hitYCorrelationHistoName       = "HitYCorrelationHisto";
#endif

EUTelExampleProcessorCorrelator::EUTelExampleProcessorCorrelator () : Processor("EUTelExampleProcessorCorrelator") {

  // modify processor description
  _description =
    "EUTelExampleProcessorCorrelator fills histograms with correlation plots between the telescope hits and the DUT ones";

  registerInputCollection(LCIO::TRACKERHIT,"TelescopeHitCollection",
                          "Telescope hit collection name",
                          _inputTelescopeCollectionName, string ( "hit" ));

  registerInputCollection(LCIO::TRACKERHIT,"DUTHitCollection",
                          "DUT hit collection name",
                          _inputDUTCollectionName, string ( "hit_DUT" ));

}


void EUTelExampleProcessorCorrelator::init() {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run counter
  _iRun = 0;

  // check if the GEAR manager pointer is not null!
  if ( Global::GEAR == 0x0 ) {
    streamlog_out ( ERROR4 ) <<  "The GearMgr is not available, for an unknown reason." << endl;
    exit(-1);
  }

  _siPlanesParameters  = const_cast<gear::SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> (&(_siPlanesParameters->getSiPlanesLayerLayout() ));

  _siPlaneZPosition = new double[ _siPlanesLayerLayout->getNLayers() ];
  for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++ ) {
    _siPlaneZPosition[ iPlane ] = _siPlanesLayerLayout->getLayerPositionZ(iPlane);
  }

  // verify that the telescope has a DUT in the steering file!
  if ( _siPlanesParameters->getSiPlanesType() != gear::SiPlanesParameters::TelescopeWithDUT ) {
    throw InvalidGeometryException( this->name() + " needs to have a telescope with DUT!" );
  }

}

void EUTelExampleProcessorCorrelator::processRunHeader (LCRunHeader * rdr) {


  EUTelRunHeaderImpl * runHeader = new EUTelRunHeaderImpl( rdr ) ;
  runHeader->addProcessor( type() );

  // this is the right place also to check the geometry ID. This is a
  // unique number identifying each different geometry used at the
  // beam test. The same number should be saved in the run header and
  // in the xml file. If the numbers are different, instead of barely
  // quitting ask the user what to do.

  if ( runHeader->getGeoID() == 0 )
    streamlog_out ( WARNING0 ) <<  "The geometry ID in the run header is set to zero." << endl
                               <<  "This may mean that the GeoID parameter was not set" << endl;


  if ( runHeader->getGeoID() != _siPlanesParameters->getSiPlanesID() ) {
    streamlog_out ( ERROR1 ) <<  "Error during the geometry consistency check: " << endl
                             << "The run header says the GeoID is " << runHeader->getGeoID() << endl
                             << "The GEAR description says is     " << _siPlanesParameters->getSiPlanesID()
                             << endl;

#ifdef EUTEL_INTERACTIVE
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
#endif
  }

  delete runHeader;

  // increment the run counter
  ++_iRun;
}


void EUTelExampleProcessorCorrelator::processEvent (LCEvent * event) {

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event) ;

  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return;
  } else if ( evt->getEventType() == kUNKNOWN ) {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event."
                               << endl;
  }

  if ( isFirstEvent() ) {

    // need to book histograms
    bookHistos();

    _isFirstEvent = false;

  }

  try {
    LCCollectionVec * inputTelescopeCollection = static_cast< LCCollectionVec *>
      ( event->getCollection( _inputTelescopeCollectionName )) ;

    LCCollectionVec * inputDUTCollection = static_cast< LCCollectionVec *>
      ( event->getCollection( _inputDUTCollectionName ));

    // loop over all the DUT hit first
    for ( size_t iDUT = 0 ; iDUT < inputDUTCollection->size(); ++iDUT ) {

      // get the current DUT hit
      TrackerHitImpl * dutHit = static_cast< TrackerHitImpl * >
        ( inputDUTCollection->getElementAt ( iDUT ) );

      // get the x,y,z position of the DUT hit
      const double * dutHitPosition =  dutHit->getPosition() ;

      // now loop over all the telescope hit
      for ( size_t iTel = 0; iTel < inputTelescopeCollection->size(); ++iTel ) {

        // get the current telescope hit
        TrackerHitImpl * telHit = static_cast< TrackerHitImpl * >
          ( inputTelescopeCollection->getElementAt( iTel ) ) ;

        // guess the sensor ID 
        int sensorID = guessSensorID( telHit );

        // get the position of the telescope plane
        const double * telHitPosition = telHit->getPosition();

        // fill the correlation histograms
        _hitXCorrelationMatrix[ sensorID ]->fill ( telHitPosition[0], dutHitPosition[0] ) ;
        _hitYCorrelationMatrix[ sensorID ]->fill ( telHitPosition[1], dutHitPosition[1] );

      } // loop over of the tel hits


    } // loop over the DUT hits


  } catch (DataNotAvailableException& e  ) {

    streamlog_out  ( WARNING2 ) <<  "No input collection found on event " << event->getEventNumber()
                                << " in run " << event->getRunNumber() << endl;
  }

#endif

}

void EUTelExampleProcessorCorrelator::end() {

  streamlog_out ( MESSAGE4 )  << "Successfully finished" << endl;
  delete [] _siPlaneZPosition;

}

void EUTelExampleProcessorCorrelator::bookHistos() {


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  try {

    streamlog_out ( MESSAGE4 ) <<  "Booking histograms" << endl;

    // create all the directories first
    vector< string > dirNames;
    dirNames.push_back ("HitX");
    dirNames.push_back ("HitY");

    for ( size_t iPos = 0 ; iPos < dirNames.size() ; iPos++ ) {

      AIDAProcessor::tree(this)->mkdir( dirNames[iPos].c_str() ) ;

    }

    string tempHistoName;
    string tempHistoTitle ;

    // get the number of detectors from GEAR.
    int nTelDetector = _siPlanesLayerLayout->getNLayers();

    // the boundaries of the histos can be read from the GEAR
    // description and for safety multiplied by a safety factor
    // to take into account possible misalignment.
    //
    // 2 should be enough because it
    // means that the sensor is wrong
    // by all its size.
    double safetyFactor = 2.0;

    // get the boundaries of the DUT detector
    double dutXMin  = safetyFactor * ( _siPlanesLayerLayout->getDUTSensitivePositionX() -
                                       ( 0.5 * _siPlanesLayerLayout->getDUTSensitiveSizeX() ));
    double dutXMax  = safetyFactor * ( _siPlanesLayerLayout->getDUTSensitivePositionX() +
                                       ( 0.5 * _siPlanesLayerLayout->getDUTSensitiveSizeX() ));
    int    dutXNBin = static_cast< int > ( safetyFactor ) * _siPlanesLayerLayout->getDUTSensitiveNpixelX();


    double dutYMin = safetyFactor * ( _siPlanesLayerLayout->getDUTSensitivePositionY() -
                                      ( 0.5 * _siPlanesLayerLayout->getDUTSensitiveSizeY() ));
    double dutYMax = safetyFactor * ( _siPlanesLayerLayout->getDUTSensitivePositionY() +
                                      ( 0.5 * _siPlanesLayerLayout->getDUTSensitiveSizeY() ));
    int    dutYNBin = static_cast< int > ( safetyFactor ) * _siPlanesLayerLayout->getDUTSensitiveNpixelY();

    // loop over all telescope planes
    for ( int iTel = 0 ; iTel < nTelDetector; ++iTel ) {

      // as declared in the gear description
      int sensorID = _siPlanesLayerLayout->getID( iTel );

      double telXMin = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionX( iTel ) -
                                        ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeX ( iTel ) ));
      double telXMax = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionX( iTel ) +
                                        ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeX ( iTel )));
      int    telXNBin = static_cast< int > ( safetyFactor ) * _siPlanesLayerLayout->getSensitiveNpixelX( iTel );


      double telYMin = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionY( iTel ) -
                                        ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeY ( iTel ) ));
      double telYMax = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionY( iTel ) +
                                        ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeY ( iTel )));
      int    telYNBin = static_cast< int > ( safetyFactor ) * _siPlanesLayerLayout->getSensitiveNpixelY( iTel );

      tempHistoName = dirNames[0] + "/" + _hitXCorrelationHistoName + to_string( sensorID );
      streamlog_out( DEBUG5 ) << "Booking histo " << tempHistoName ;

      tempHistoTitle = "Correlation of the DUT and the X detector " + to_string( sensorID );

      // book the histogram
      AIDA::IHistogram2D * histo2D =
        AIDAProcessor::histogramFactory( this )->createHistogram2D( tempHistoName.c_str(), telXNBin, telXMin, telXMax,
                                                                    dutXNBin, dutXMin, dutXMax);
      histo2D->setTitle( tempHistoTitle.c_str() );

      // add it to the associative map. I'm using iTel instead of
      // sensorID because it is way too practical
      _hitXCorrelationMatrix[ sensorID ] = histo2D;

      // repeat for the y direction
      tempHistoName = dirNames[1] + "/" + _hitYCorrelationHistoName + to_string( sensorID );
      streamlog_out( DEBUG5 ) << "Booking histo " << tempHistoName ;

      tempHistoTitle = "Correlation of the DUT and the Y detector " + to_string( sensorID );

      // book the histogram
      histo2D = AIDAProcessor::histogramFactory( this )->createHistogram2D( tempHistoName.c_str(), telYNBin, telYMin, telYMax,
                                                                            dutYNBin, dutYMin, dutYMax);
      histo2D->setTitle( tempHistoTitle.c_str() );

      // add it to the associative map. I'm using iTel instead of
      // sensorID because it is way too practical
      _hitYCorrelationMatrix[ sensorID ] = histo2D;

    }

  } catch (lcio::Exception& e ) {

    streamlog_out ( ERROR1 ) << "No AIDAProcessor initialized. Sorry for quitting..." << endl;
    exit( -1 );

  }
#endif
}

int EUTelExampleProcessorCorrelator::guessSensorID( TrackerHitImpl * hit ) {

  int sensorID = -1;
  double minDistance =  numeric_limits< double >::max() ;
  double * hitPosition = const_cast<double * > (hit->getPosition());
  double zPos = 0;

  for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); ++iPlane ) {
    double distance = std::abs( hitPosition[2] - _siPlaneZPosition[ iPlane ] );
    if ( distance < minDistance ) {
      minDistance = distance;
      sensorID = _siPlanesLayerLayout->getID( iPlane );
      zPos = _siPlaneZPosition[ iPlane ];
    }
  }
  if ( minDistance > 5 /* mm */ ) {
    // advice the user that the guessing wasn't successful 
    streamlog_out( WARNING3 ) << "A hit was found " << minDistance << " mm far from the nearest plane\n"
      "Please check the consistency of the data with the GEAR file" << endl;
    streamlog_out( WARNING3 ) << "SensorID = " << sensorID << " expected pos = " << zPos << " meas pos = " << hitPosition[2] << endl;
  }

  return sensorID;
}


#endif // USE_GEAR
