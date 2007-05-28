// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelHitMaker.cc,v 1.2 2007-05-28 11:52:45 bulgheroni Exp $
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
#include "EUTelExceptions.h"
#include "EUTelEtaFunctionImpl.h"


// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// aida includes <.h>
#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud2D.h>
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

using namespace std;
using namespace marlin;
using namespace gear;
using namespace eutelescope;

// definition of static members mainly used to name histograms
#ifdef MARLIN_USE_AIDA
std::string EUTelHitMaker::_hitCloudLocalName     = "HitCloudLocal";
std::string EUTelHitMaker::_hitCloudTelescopeName = "HitCloudTelescope";
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
  
  vector<string > etaNames;
  etaNames.push_back("xEta");
  etaNames.push_back("yEta");
   
  registerOptionalParameter("EtaCollectionName", 
			    "The name of the collections containing the eta function (x and y respectively)",
			    _etaCollectionNames, etaNames, etaNames.size());


  registerOptionalParameter("EtaSwitch","Enable (==1) or disable eta correction",
			    _etaCorrection, static_cast< int > ( 1 ) ); 
  
}


void EUTelHitMaker::init() {
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

  // now book histograms plz...
  bookHistos();
}

void EUTelHitMaker::processRunHeader (LCRunHeader * rdr) {


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
	
  

  // increment the run counter
  ++_iRun;

  

}


void EUTelHitMaker::processEvent (LCEvent * event) {


  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event) ;
  
  if ( evt->getEventType() == kEORE ) {
    message<DEBUG> ( "EORE found: nothing else to do." );
    return;
  }

  
  LCCollectionVec * pulseCollection   = static_cast<LCCollectionVec*> (event->getCollection( _pulseCollectionName ));
  LCCollectionVec * clusterCollection = static_cast<LCCollectionVec*> (event->getCollection("original_data"));
  LCCollectionVec * xEtaCollection, * yEtaCollection;
  LCCollectionVec * hitCollection   = new LCCollectionVec(LCIO::TRACKERHIT);

  CellIDDecoder<TrackerPulseImpl>  pulseCellDecoder(pulseCollection);
  CellIDDecoder<TrackerDataImpl>   clusterCellDecoder(clusterCollection);

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
  }

  if ( isFirstEvent() && _etaCorrection == 1) {
    if ( ( xEtaCollection->getNumberOfElements() != pulseCollection->getNumberOfElements() ) ||
	 ( yEtaCollection->getNumberOfElements() != pulseCollection->getNumberOfElements() ) ) {
      message<ERROR> ( log() 
		       <<  "The eta collections contain a different number of elements wrt to "
		       << _siPlanesParameters->getSiPlanesNumber() );
      message<ERROR> ( log() << "Continuing without eta correction " );
      _etaCorrection = 0;
    }
  }

  int detectorID    = -99; // it's a non sense
  int oldDetectorID = -100;

  int    layerIndex;
  double xZero, yZero, zZero;
  double zThickness;
  double xPitch, yPitch;
  int    xPointing[2], yPointing[2];

  for ( int iPulse = 0; iPulse < pulseCollection->getNumberOfElements(); iPulse++ ) {
    
    TrackerPulseImpl     * pulse   = static_cast<TrackerPulseImpl*> ( pulseCollection->getElementAt(iPulse) );
    EUTelVirtualCluster  * cluster;
    ClusterType type = static_cast<ClusterType>(static_cast<int>((pulseCellDecoder(pulse)["type"])));    

    if ( type == kEUTelFFClusterImpl ) {
      cluster = static_cast<EUTelFFClusterImpl*> ( pulse->getTrackerData() );
    } else {
      message<ERROR> ( "Unknown cluster type. Sorry for quitting" );
      throw UnknownDataTypeException("Cluster type unknown");
    }
    
    // there could be several clusters belonging to the same
    // detector. So update the geometry information only if this new
    // cluster belongs to a different detector.
    detectorID = clusterCellDecoder(cluster)["sensorID"];

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
      
      layerIndex = _conversionIdMap[detectorID];
      xZero      = _siPlanesLayerLayout->getSensitivePositionX(layerIndex); // mm
      yZero      = _siPlanesLayerLayout->getSensitivePositionY(layerIndex); // mm
      zZero      = _siPlanesLayerLayout->getSensitivePositionZ(layerIndex); // mm
      zThickness = _siPlanesLayerLayout->getSensitiveThickness(layerIndex); //  ??? 
      
      xPitch  = 30 /* um */ * 0.001 ; // mm
      yPitch  = 30 /* um */ * 0.001 ; // mm
      
      xPointing[0] = -1 ;       xPointing[1] =  0 ;
      yPointing[0] =  0 ;       yPointing[1] = -1 ;

      
      if ( isFirstEvent() && (iPulse == 0) ) message<WARNING> ( "Using hardcoded values for the pitches and for the orientation." );
    }
    
    // get the position of the seed pixel. This is in pixel number.
    int xCluCenter, yCluCenter;
    cluster->getSeedCoord( xCluCenter, yCluCenter );
    
    // with the charge center of gravity calculation, we get a shift
    // from the seed pixel center due to the charge distribution. Those
    // two numbers are the correction values in the case the Eta
    // correction is not applied.
//     float xShift, yShift;
//     cluster->getCenterOfGravityShift( xShift, yShift );
//     double xCorrection = static_cast<double> (xShift) ;
//     double yCorrection = static_cast<double> (yShift) ;

    
//     if ( _etaCorrection == 1 ) {
      
//       EUTelEtaFunctionImpl * xEtaFunc = static_cast<EUTelEtaFunctionImpl*> ( xEtaCollection->getElementAt(detectorID) );
//       EUTelEtaFunctionImpl * yEtaFunc = static_cast<EUTelEtaFunctionImpl*> ( yEtaCollection->getElementAt(detectorID) );

//       xCorrection = xEtaFunc->getEtaFromCoG( xShift );
//       yCorrection = yEtaFunc->getEtaFromCoG( yShift );
      
//     }

//     // rescale the pixel number in millimeter
//     double xDet = ( static_cast<double> (xCluCenter) + xCorrection + 0.5 ) * xPitch ;
//     double yDet = ( static_cast<double> (yCluCenter) + yCorrection + 0.5 ) * yPitch ;

// #ifdef MARLIN_USE_AIDA
//     string tempHistoName;
//     {
//       stringstream ss; 
//       ss << _hitCloudLocalName << "-" << detectorID ;
//       tempHistoName = ss.str();
//     }
//     (dynamic_cast<AIDA::ICloud2D*>(_aidaHistoMap[ tempHistoName ]))->fill(xDet,yDet);

// #endif 

//     // now perform the rotation of the frame of references and put the
//     // results already into a 3D array of double to be ready for the
//     // setPosition method of TrackerHit
//     double telPos[3];
//     telPos[0] = xPointing[0] * xDet + xPointing[1] * yDet;
//     telPos[1] = yPointing[0] * xDet + yPointing[1] * yDet;

//     // now the translation
//     telPos[0] -= xZero;
//     telPos[1] -= yZero;
//     telPos[2] = zZero + 0.5 * zThickness;

// #ifdef MARLIN_USE_AIDA
//     {
//       stringstream ss;
//       ss << _hitCloudTelescopeName << "-" << detectorID ;
//       tempHistoName = ss.str();
//     }
//     (dynamic_cast<AIDA::ICloud2D*> (_aidaHistoMap[ tempHistoName ] ))->fill( telPos[0], telPos[1] );
// #endif

//     // create the new hit
//     TrackerHitImpl * hit = new TrackerHitImpl;
//     hit->setPosition( &telPos[0] );
//     hit->setType( pulseCellDecoder(pulse)["type"] );
    
//     // prepare a LCObjectVec to store the current cluster
//     LCObjectVec clusterVec;
//     clusterVec.push_back( cluster );
        
//     // add the clusterVec to the hit
//     hit->rawHits() = clusterVec;

//     // add the new hit to the hit collection
//     hitCollection->push_back( hit );

  }
  ++_iEvt;
  //   evt->addCollection( hitCollection, _hitCollectionName );
  
  if ( isFirstEvent() ) _isFirstEvent = false;

}

void EUTelHitMaker::end() {

  message<MESSAGE> ( log() << "Successfully finished" ) ;  
}

void EUTelHitMaker::bookHistos() {
  
#ifdef MARLIN_USE_AIDA
  message<MESSAGE> ( "Booking histograms" );

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
      ss <<  _hitCloudLocalName << "-" << iDet ;
      tempHistoName = ss.str();
    }
    AIDA::ICloud2D * hitCloudLocal = AIDAProcessor::histogramFactory(this)->createCloud2D( (basePath + tempHistoName).c_str() );
    _aidaHistoMap.insert( make_pair( tempHistoName, hitCloudLocal ) );

    {
      stringstream ss ;
      ss <<  _hitCloudTelescopeName << "-" << iDet ;
      tempHistoName = ss.str();
    }
    AIDA::ICloud2D * hitCloudTelescope = AIDAProcessor::histogramFactory(this)->createCloud2D( ( basePath + tempHistoName ).c_str() );
    _aidaHistoMap.insert( make_pair ( tempHistoName, hitCloudTelescope ) );
											   
  }


  
#endif
}

#endif
