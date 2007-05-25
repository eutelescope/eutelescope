// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelHitMaker.cc,v 1.1 2007-05-25 05:16:55 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h" 
#include "EUTelHitMaker.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelEtaFunctionImpl.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"

// gear includes <.h>
#ifdef USE_GEAR
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>
#endif

// lcio includes <.h> 
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerPulseImpl.h>

// system includes <>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;
using namespace marlin;
using namespace gear;
using namespace eutelescope;



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
  
  vector<string > etaNames;
  etaNames.push_back("xEta");
  etaNames.push_back("yEta");
   
  registerOptionalParameter("EtaCollectionName", 
			    "The name of the collections containing the eta function (x and y respectively)",
			    _etaCollectionNames, etaNames, etaNames.size());


  registerOptionalParameter("EtaSwitch","Enable (==1) or disable eta correction",
			    _etaCorrection, static_cast< int > ( 1 ) ); 
  
}


void EUTelHitMaker::init () {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

  // check if Marlin was built with GEAR support or not
#ifndef USE_GEAR

  message<ERROR> ( "Marlin was not built with GEAR support." );
  message<ERROR> ( "You need to install GEAR and recompile Marlin with -DUSE_GEAR before continue.");
  
  // I'm thinking if this is the case of throwing an exception or
  // not. This is a really error and not something that can
  // exceptionally happens. Still not sure what to do
  exit(-1);

#else

  // check if the GEAR manager pointer is not null!
  if ( Global::GEAR == 0x0 ) {
    message<ERROR> ( "The GearMgr is not available, for an unknown reason." );
    exit(-1);
  }

  _siPlanesParameters  = const_cast<SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

  
#endif
}

void EUTelHitMaker::processRunHeader (LCRunHeader * rdr) {


#ifdef USE_GEAR
  EUTelRunHeaderImpl * header = static_cast<EUTelRunHeaderImpl*> (rdr);

  // the run header contains the number of detectors. This number
  // should be in principle be the same as the number of layers in the
  // geometry description
  if ( header->getNoOfDetector() != _siPlanesParameters->getSiPlanesNumber() ) {
    message<ERROR> ( "Error during the geometry consistency check: " );
    message<ERROR> ( log() << "The run header says there are " << header->getNoOfDetector() << " silicon detectors " );
    message<ERROR> ( log() << "The GEAR description says     " << _siPlanesParameters->getSiPlanesNumber() << " silicon planes" );
    exit(-1);
  }
  
  // this is the right place also to check the geometry ID. This is a
  // unique number identifying each different geometry used at the
  // beam test. The same number should be saved in the run header and
  // in the xml file. If the numbers are different, instead of barely
  // quitting ask the user what to do.

  // UNCOMMENT THE FOLLOWING PART WHEN THE GEO ID WILL BE IMPLEMENTED
  // ALSO IN THE XML GEAR DESCRIPTION
  //
  //   if ( header->getGeoID() != _siPlanesParameters->getSiPlanesNumber() ) {
  //     message<ERROR> ( "Error during the geometry consistency check: " );
  //     message<ERROR> ( log() << "The run header says the GeoID is " << header->getGeoID() );
  //     message<ERROR> ( log() << "The GEAR description says is     " << _siPlanesParameters->getSiPlanesNumber() );
  //     string answer;
  //     while (true) {
  //       message<ERROR> ( "Type Q to quit now or C to continue using the actual GEAR description anyway [Q/C]" );
  //       cin >> answer;
  //       // put the answer in lower case before making the comparison.
  //       transform( answer.begin(), answer.end(), answer.begin(), ::tolower );
  //       if ( answer == "q" ) {
  // 	exit(-1);
  //       } else if ( answer == "c" ) {
  // 	break;
  //       }
  //     }
  //   }
	
  
  
#endif 

  // increment the run counter
  ++_iRun;

  

}


void EUTelHitMaker::processEvent (LCEvent * event) {

#ifdef USE_GEAR

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event) ;
  
  if ( evt->getEventType() == kEORE ) {
    message<DEBUG> ( "EORE found: nothing else to do." );
    return;
  }

  
  LCCollectionVec * pulseCollection = static_cast<LCCollectionVec*> (event->getCollection( _pulseCollectionName ));
  LCCollectionVec * xEtaCollection, * yEtaCollection;

  if ( _etaCorrection == 1 ) {
    // this means that the user wants to apply the eta correction to
    // the center of gravity, so we need to have the two Eta function
    // collections
    
    try {  
      xEtaCollection = static_cast<LCCollectionVec*> (event->getCollection( _etaCollectionNames[0] )) ;
    } catch (DataNotAvailableException& e) {
      message<ERROR> ( log() << "The eta collection " << _etaCollectionNames[0] << " is not available" );
      message<ERROR> ( log() << "Continuing without eta correction " ) ;
      _etaCorrection = 0;
    }

    try {  
      yEtaCollection = static_cast<LCCollectionVec*> (event->getCollection( _etaCollectionNames[1] )) ;
    } catch (DataNotAvailableException& e) {
      message<ERROR> ( log() << "The eta collection " << _etaCollectionNames[1] << " is not available" );
      message<ERROR> ( log() << "Continuing without eta correction " ) ;
      _etaCorrection = 0;
    }


    if ( isFirstEvent() ) {
      if ( ( xEtaCollection->getNumberOfElements() != pulseCollection->getNumberOfElements() ) ||
	   ( yEtaCollection->getNumberOfElements() != pulseCollection->getNumberOfElements() ) ) {
	message<ERROR> ( log() 
			 <<  "The eta collections contain a different number of elements wrt to "
			 << _siPlanesParameters->getSiPlanesNumber() );
	message<ERROR> ( log() << "Continuing without eta correction " );
	_etaCorrection = 0;
      }
    }

			 

  }

 
  for ( int iPulse = 0; iPulse < pulseCollection->getNumberOfElements(); iPulse++ ) {
    
    TrackerPulseImpl     * pulse   = static_cast<TrackerPulseImpl*> ( pulseCollection->getElementAt(iPulse) );
    EUTelVirtualCluster  * cluster = static_cast<EUTelVirtualCluster*> ( pulse->getTrackerData() );
    
    int detectorID = cluster->getDetectorID();
    
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
      
    int layerIndex = _conversionIdMap[detectorID];
    double xZero   = _siPlanesLayerLayout->getSensitivePositionX(layerIndex); // mm
    double yZero   = _siPlanesLayerLayout->getSensitivePositionY(layerIndex); // mm
    double zZero   = _siPlanesLayerLayout->getSensitivePositionZ(layerIndex); // mm
    
    double xPitch  = 20 /* um */ * 0.001 ; // mm
    double yPitch  = 20 /* um */ * 0.001 ; // mm
      
    int xCluCenter, yCluCenter;
    cluster->getSeedCoord( xCluCenter, yCluCenter );
    
    float xShift, yShift;
    cluster->getCenterOfGravityShift( xShift, yShift );

    double xCorrection = 0, yCorrection = 0;

    if ( _etaCorrection == 1 ) {
      
      EUTelEtaFunctionImpl * xEtaFunc = static_cast<EUTelEtaFunctionImpl*> ( xEtaCollection->getElementAt(detectorID) );
      EUTelEtaFunctionImpl * yEtaFunc = static_cast<EUTelEtaFunctionImpl*> ( yEtaCollection->getElementAt(detectorID) );

      xCorrection = xEtaFunc->getEtaFromCoG( xShift );
      yCorrection = yEtaFunc->getEtaFromCoG( yShift );

    }

    
  }
  ++_iEvt;
  
#endif 
  
}
  



void EUTelHitMaker::end() {



  message<MESSAGE> ( log() << "Successfully finished" ) ;  
}

