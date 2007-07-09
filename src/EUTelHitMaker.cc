// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelHitMaker.cc,v 1.12 2007-07-09 13:40:52 bulgheroni Exp $
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
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ICloud.h>
#include <AIDA/ICloud2D.h>
#include <AIDA/ICloud3D.h>
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
std::string EUTelHitMaker::_hitCloudLocalName          = "HitCloudLocal";
std::string EUTelHitMaker::_hitCloudTelescopeName      = "HitCloudTelescope";
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

  _histogramSwitch = true;

#endif

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

  if ( header->getGeoID() == 0 ) 
    message<WARNING> ( "The geometry ID in the run header is set to zero.\n" 
		       "This may mean that the GeoID parameter was not set" );
  

  if ( header->getGeoID() != _siPlanesParameters->getSiPlanesID() ) {
    message<ERROR> ( "Error during the geometry consistency check: " );
    message<ERROR> ( log() << "The run header says the GeoID is " << header->getGeoID() );
    message<ERROR> ( log() << "The GEAR description says is     " << _siPlanesParameters->getSiPlanesID() );
    string answer;
    while (true) {
      message<ERROR> ( "Type Q to quit now or C to continue using the actual GEAR description anyway [Q/C]" );
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


  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event) ;
  
  if ( evt->getEventType() == kEORE ) {

    message<DEBUG> ( "EORE found: nothing else to do." );
    
    return;
  }

  if ( (_iEvt % 10) == 0 ) 
    message<MESSAGE> ( log() << "Applying geometry to the pulses of event " << _iEvt );
  
  LCCollectionVec * pulseCollection   = static_cast<LCCollectionVec*> (event->getCollection( _pulseCollectionName ));
  LCCollectionVec * clusterCollection = static_cast<LCCollectionVec*> (event->getCollection("original_data"));
  LCCollectionVec * hitCollection     = new LCCollectionVec(LCIO::TRACKERHIT);
  LCCollectionVec * xEtaCollection, * yEtaCollection;

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
    if ( ( xEtaCollection->getNumberOfElements() != _siPlanesParameters->getSiPlanesNumber() ) ||
	 ( yEtaCollection->getNumberOfElements() != _siPlanesParameters->getSiPlanesNumber() ) ) {
      message<ERROR> ( log() 
		       <<  "The eta collections contain a different number of elements wrt to "
		       << _siPlanesParameters->getSiPlanesNumber() );
      message<ERROR> ( log() << "Continuing without eta correction " );
      _etaCorrection = 0;
    } else {
#ifdef MARLIN_USE_AIDA
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
	    message<ERROR> ( log() << "Problem booking the " << (basePath + tempHistoName) << ".\n"
			     << "Very likely a problem with path name. Switching off histogramming and continue w/o");
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
	    message<ERROR> ( log() << "Problem booking the " << (basePath + tempHistoName) << ".\n"
			     << "Very likely a problem with path name. Switching off histogramming and continue w/o");
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
	    message<ERROR> ( log() << "Problem booking the " << (basePath + tempHistoName) << ".\n"
			     << "Very likely a problem with path name. Switching off histogramming and continue w/o");
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
	    message<ERROR> ( log() << "Problem booking the " << (basePath + tempHistoName) << ".\n"
			     << "Very likely a problem with path name. Switching off histogramming and continue w/o");
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
	    message<ERROR> ( log() << "Problem booking the " << (basePath + tempHistoName) << ".\n"
			     << "Very likely a problem with path name. Switching off histogramming and continue w/o");
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
	    message<ERROR> ( log() << "Problem booking the " << (basePath + tempHistoName) << ".\n"
			     << "Very likely a problem with path name. Switching off histogramming and continue w/o");
	    _histogramSwitch = false;
	  }
	  
	}
      }
#endif
      
    }
  }
  int detectorID    = -99; // it's a non sense
  int oldDetectorID = -100;

  int    layerIndex;
  double xZero, yZero, zZero;
  double xSize, ySize;
  double zThickness;
  double xPitch, yPitch;
  double xPointing[2], yPointing[2];

  for ( int iPulse = 0; iPulse < pulseCollection->getNumberOfElements(); iPulse++ ) {
    
    TrackerPulseImpl     * pulse   = static_cast<TrackerPulseImpl*> ( pulseCollection->getElementAt(iPulse) );
    EUTelVirtualCluster  * cluster;
    ClusterType type = static_cast<ClusterType>(static_cast<int>((pulseCellDecoder(pulse)["type"])));    

    if ( type == kEUTelFFClusterImpl ) {
      cluster = new EUTelFFClusterImpl( static_cast<TrackerDataImpl*> (pulse->getTrackerData()) );
    } else {
      message<ERROR> ( "Unknown cluster type. Sorry for quitting" );
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

    }
    
    // get the position of the seed pixel. This is in pixel number.
    int xCluCenter, yCluCenter;
    cluster->getSeedCoord(xCluCenter, yCluCenter);
    
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

      if ( ( xShift >= -0.5 ) && ( xShift <= 0.5 ) ) {
	xCorrection = xEtaFunc->getEtaFromCoG( xShift );
      } else {
	message<DEBUG> ( log() << "Found anomalous cluster\n" << ( * cluster ) );
      }
      if ( ( yShift >= -0.5 ) && ( yShift <= 0.5 ) ) {
	yCorrection = yEtaFunc->getEtaFromCoG( yShift );
      } else {
	message<DEBUG> ( log() << "Found anomalous cluster\n" << ( * cluster ) );	
      }

#ifdef MARLIN_USE_AIDA
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

#ifdef MARLIN_USE_AIDA
    string tempHistoName;
    if ( _histogramSwitch ) {
      {
	stringstream ss; 
	ss << _hitCloudLocalName << "-" << detectorID ;
	tempHistoName = ss.str();
      }
      if ( AIDA::ICloud2D* cloud = dynamic_cast<AIDA::ICloud2D*>(_aidaHistoMap[ tempHistoName ]) )
	cloud->fill(xDet, yDet);
      else {
	message<ERROR> ( log() << "Not able to retrieve histogram pointer for " << tempHistoName 
			 << ".\nDisabling histogramming from now on " );
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
    // configuration but we need to double checkit
    telPos[0] -= ( xZero - xSize/2 );
    telPos[1] -= ( yZero - ySize/2 );
    telPos[2] = zZero + 0.5 * zThickness;

#ifdef MARLIN_USE_AIDA
    if ( _histogramSwitch ) {
      {
	stringstream ss;
	ss << _hitCloudTelescopeName << "-" << detectorID ;
	tempHistoName = ss.str();
      }
      AIDA::ICloud2D * cloud2D = dynamic_cast<AIDA::ICloud2D*> (_aidaHistoMap[ tempHistoName ] );
      if ( cloud2D ) cloud2D->fill( telPos[0], telPos[1] );
      else {
	message<ERROR> ( log() << "Not able to retrieve histogram pointer for " << tempHistoName 
			 << ".\nDisabling histogramming from now on " );
	_histogramSwitch = false;
      }
      AIDA::ICloud3D * cloud3D = dynamic_cast<AIDA::ICloud3D*> (_aidaHistoMap[ _densityPlotName ] );
      if ( cloud3D ) cloud3D->fill( telPos[0], telPos[1], telPos[2] );
      else {
	message<ERROR> ( log() << "Not able to retrieve histogram pointer for " << tempHistoName 
			 << ".\nDisabling histogramming from now on " );
	_histogramSwitch = false;
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
  ++_iEvt;
  evt->addCollection( hitCollection, _hitCollectionName );
  
  if ( isFirstEvent() ) _isFirstEvent = false;

}

void EUTelHitMaker::end() {

#ifdef MARLIN_USE_AIDA
    if ( _histogramSwitch ) {
      
      message<DEBUG> ( log() << "Converting clouds to histograms before leaving " );
      
      map<string, AIDA::IBaseHistogram *>::iterator mapIter = _aidaHistoMap.begin();
      while ( mapIter != _aidaHistoMap.end() ) {

 	AIDA::ICloud * cloud = dynamic_cast<AIDA::ICloud*> ( mapIter->second );
 	if ( cloud )  {
	  cloud->convertToHistogram();
	}
 	++mapIter;
      }
    }
    
#endif

  message<MESSAGE> ( log() << "Successfully finished" ) ;  
}

void EUTelHitMaker::bookHistos() {
  
#ifdef MARLIN_USE_AIDA

  try {
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
      if ( hitCloudLocal ) {
	hitCloudLocal->setTitle("Hit map in the detector local frame of reference");
	_aidaHistoMap.insert( make_pair( tempHistoName, hitCloudLocal ) );
      } else {
	message<ERROR> ( log() << "Problem booking the " << (basePath + tempHistoName) << ".\n"
			 << "Very likely a problem with path name. Switching off histogramming and continue w/o");
	_histogramSwitch = false;
      }
      
      
      
      {
	stringstream ss ;
	ss <<  _hitCloudTelescopeName << "-" << iDet ;
	tempHistoName = ss.str();
      }
      AIDA::ICloud2D * hitCloudTelescope = AIDAProcessor::histogramFactory(this)->createCloud2D( ( basePath + tempHistoName ).c_str() );
      if ( hitCloudTelescope ) {
	hitCloudTelescope->setTitle("Hit map in the telescope frame of reference");
	_aidaHistoMap.insert( make_pair ( tempHistoName, hitCloudTelescope ) );
      } else {
	message<ERROR> ( log() << "Problem booking the " << (basePath + tempHistoName) << ".\n"
			 << "Very likely a problem with path name. Switching off histogramming and continue w/o");
	_histogramSwitch = false;
      }
    }

    AIDA::ICloud3D * densityPlot = AIDAProcessor::histogramFactory(this)->createCloud3D( _densityPlotName );
    if ( densityPlot ) {
      densityPlot->setTitle("Hit position in the telescope frame of reference");
      _aidaHistoMap.insert( make_pair ( _densityPlotName, densityPlot ) ) ;
    } else {
      message<ERROR> ( log() << "Problem booking the " << (_densityPlotName) << ".\n"
		       << "Very likely a problem with path name. Switching off histogramming and continue w/o");
      _histogramSwitch = false;
    }

	
      


  } catch (lcio::Exception& e ) {
    
    message<ERROR> ( log() << "No AIDAProcessor initialized. Type q to exit or c to continue without histogramming" );
    string answer;
    while ( true ) {
      message<ERROR> ( "[q]/[c]" );
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
