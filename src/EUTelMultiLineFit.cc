// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Philipp Roloff, DESY <mailto:philipp.roloff@desy.de>
// Version: $Id $
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
#ifdef MARLIN_USE_AIDA
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
#ifdef MARLIN_USE_AIDA
std::string EUTelMultiLineFit::_numberTracksLocalname   = "NumberTracks";
std::string EUTelMultiLineFit::_chi2XLocalname          = "Chi2X";
std::string EUTelMultiLineFit::_chi2YLocalname          = "Chi2Y";
std::string EUTelMultiLineFit::_angleXLocalname         = "AngleX";
std::string EUTelMultiLineFit::_angleYLocalname         = "AngleY";
std::string EUTelMultiLineFit::_residualXLocalname      = "ResidualX";
std::string EUTelMultiLineFit::_residualYLocalname      = "ResidualY";
#endif

EUTelMultiLineFit::EUTelMultiLineFit () : Processor("EUTelMultiLineFit") {
  
  // modify processor description
  _description =
    "EUTelMultiLineFit will fit a straight line";
  
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

  registerOptionalParameter("DistanceMax","Maximal allowed distance between hits entering the fit.",
                            _distanceMax, static_cast <float> (2000.0));

  registerOptionalParameter("Chi2XMax","Maximal chi2 for fit of x coordinate."
			    ,_chi2XMax, static_cast <float> (10000.0));

  registerOptionalParameter("Chi2YMax","Maximal chi2 for fit of y coordinate."                             ,_chi2YMax, static_cast <float> (10000.0));

  registerOptionalParameter("ExcludePlane","Exclude plane from fit."
			    ,_excludePlane, static_cast <int> (0));

}

void EUTelMultiLineFit::init() {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();
  
  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;
  
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
    exit(-1);
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
    // float Chiquare[2] = {0,0};
    // float angle[2] = {0,0};
    
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
	cluster = new EUTelFFClusterImpl( static_cast<TrackerDataImpl *> ( clusterVector[0] ) );
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
	
	// here we take intrinsic resolution from geometry database
	
	layerIndex   = _conversionIdMap[detectorID];     
	_intrResolX[layerIndex] = 1000*_siPlanesLayerLayout->getSensitiveResolution(layerIndex); //um
	_intrResolY[layerIndex] = 1000*_siPlanesLayerLayout->getSensitiveResolution(layerIndex); //um
	
      }
      
      // Old code. Should be removed soon.
      
      //       _xPos[iHit] = 1000 * hit->getPosition()[0]; // in um
      //       _yPos[iHit] = 1000 * hit->getPosition()[1]; // in um
      //       _zPos[iHit] = 1000 * hit->getPosition()[2]; // in um
      
      // Getting positions of the hits.
      // Here the alignment constants are used to correct the positions.
      
      layerIndex   = _conversionIdMap[detectorID];     
      
      // The other layers were aligned with respect to the first one.
      
      double off_x, off_y, theta_x, theta_y, theta_z;

      if (layerIndex == 0) {
	
	hitsInPlane.measuredX = 1000 * hit->getPosition()[0]; // in um
	hitsInPlane.measuredY = 1000 * hit->getPosition()[1]; // in um
	hitsInPlane.measuredZ = 1000 * hit->getPosition()[2]; // in um

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
	
	hitsInPlane.measuredX = (cos(theta_y)*cos(theta_z)) * hit->getPosition()[0] * 1000 + ((-1)*sin(theta_x)*sin(theta_y)*cos(theta_z) + cos(theta_x)*sin(theta_z)) * hit->getPosition()[1] * 1000 + off_x;
	hitsInPlane.measuredY = ((-1)*cos(theta_y)*sin(theta_z)) * hit->getPosition()[0] * 1000 + (sin(theta_x)*sin(theta_y)*sin(theta_z) + cos(theta_x)*cos(theta_z)) * hit->getPosition()[1] * 1000 + off_y;
	hitsInPlane.measuredZ = 1000 * hit->getPosition()[2];
	
      }
      
      delete cluster; // <--- destroying the cluster   

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

    double distance12 = 0.0;
    double distance23 = 0.0;
    double distance34 = 0.0;
    double distance45 = 0.0;
    double distance56 = 0.0;

    _xPos = new double *[500];
    _yPos = new double *[500];
    _zPos = new double *[500];

    for (int help = 0; help < 500; help++) {
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

    // loop over all hits in first plane
    for (int firsthit = 0; uint(firsthit) < _hitsFirstPlane.size(); firsthit++) {
      
      // loop over all hits in second plane
      for (int secondhit = 0; uint(secondhit) < _hitsSecondPlane.size(); secondhit++) {

	distance12 = sqrt(pow(_hitsFirstPlane[firsthit].measuredX - _hitsSecondPlane[secondhit].measuredX,2) + pow(_hitsFirstPlane[firsthit].measuredY - _hitsSecondPlane[secondhit].measuredY,2));

	if (_nPlanes == 2 && distance12 < _distanceMax) {

	  _xPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredX;
	  _yPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredY;
	  _zPos[_nTracks][0] = _hitsFirstPlane[firsthit].measuredZ;

	  _xPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredX;
	  _yPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredY;
	  _zPos[_nTracks][1] = _hitsSecondPlane[secondhit].measuredZ;

	  _nTracks++;

	}
	
	// more than two planes
	if (_nPlanes > 2) {

	  // loop over all hits in third plane
	  for (int thirdhit = 0; uint(thirdhit) < _hitsThirdPlane.size(); thirdhit++) {

	    distance23 = sqrt(pow(_hitsSecondPlane[secondhit].measuredX - _hitsThirdPlane[thirdhit].measuredX,2) + pow(_hitsSecondPlane[secondhit].measuredY - _hitsThirdPlane[thirdhit].measuredY,2));

	    if (_nPlanes == 3 && distance12 < _distanceMax && distance23 < _distanceMax) {

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

	    // more than three planes
	    if (_nPlanes > 3) {
	    
	      // loop over all hits in fourth plane
	      for (int fourthhit = 0; uint(fourthhit) < _hitsFourthPlane.size(); fourthhit++) {

		distance34 = sqrt(pow(_hitsThirdPlane[thirdhit].measuredX - _hitsFourthPlane[fourthhit].measuredX,2) + pow(_hitsThirdPlane[thirdhit].measuredY - _hitsFourthPlane[fourthhit].measuredY,2));

		if (_nPlanes == 4 && distance12 < _distanceMax && distance23 < _distanceMax && distance34 < _distanceMax) {

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
	    
		// more than four planes
		if (_nPlanes > 4) {

		  // loop over all hits in fifth plane
		  for (int fifthhit = 0; uint(fifthhit) < _hitsFifthPlane.size(); fifthhit++) {

		    distance45 = sqrt(pow(_hitsFourthPlane[fourthhit].measuredX - _hitsFifthPlane[fifthhit].measuredX,2) + pow(_hitsFourthPlane[fourthhit].measuredY - _hitsFifthPlane[fifthhit].measuredY,2));

		    if (_nPlanes == 5 && distance12 < _distanceMax && distance23 < _distanceMax && distance34 < _distanceMax && distance45 < _distanceMax) {

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

		    // more than five planes
		    if (_nPlanes > 5) {

		      // loop over all hits in sixth plane
		      for (int sixthhit = 0; uint(sixthhit) < _hitsSixthPlane.size(); sixthhit++) {
			
			distance56 = sqrt(pow(_hitsFifthPlane[fifthhit].measuredX - _hitsSixthPlane[sixthhit].measuredX,2) + pow(_hitsFifthPlane[fifthhit].measuredY - _hitsSixthPlane[sixthhit].measuredY,2));
			
			if (_nPlanes == 6 && distance12 < _distanceMax && distance23 < _distanceMax && distance34 < _distanceMax && distance45 < _distanceMax && distance56 < _distanceMax) {

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

    streamlog_out ( MESSAGE2 ) << "Number of hits in the individual planes: " << _hitsFirstPlane.size() << " " << _hitsSecondPlane.size() << " " << _hitsThirdPlane.size() << " " << _hitsFourthPlane.size() << " " << _hitsFifthPlane.size() << " " << _hitsSixthPlane.size() << endl;

    streamlog_out ( MESSAGE2 ) << "Number of tracks found in event " << _iEvt << ": " << _nTracks << endl;

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

      streamlog_out ( MESSAGE2 ) << "Fitting track using the following coordinates: ";

      for (int help = 0; help < _nPlanes; help++) {
	streamlog_out ( MESSAGE2 ) << _xPosHere[help] << " " << _yPosHere[help] << " " << _zPosHere[help] << "   ";
      }

      streamlog_out ( MESSAGE2 ) << endl;

      FitTrack(int(_nPlanes), _xPosHere, _yPosHere, _zPosHere, _intrResolX, _intrResolY, Chiquare, _waferResidX, _waferResidY, angle);

    // not needed so far
    /*

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
    
    fittrack->setIsReferencePointPCA(false);  
    
    // Calculate positions of fitted track in every plane
    
    int counter;

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

    */

#ifdef MARLIN_USE_AIDA
    
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
    
#ifdef MARLIN_USE_AIDA

	string tempHistoName;
	
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
    
    
	for( int iDetector = 0; iDetector < _nPlanes; iDetector++ ){
	  
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
	}
	
#endif 

      } // end if chi2 cut

      // clean up
      delete [] _zPosHere;
      delete [] _yPosHere;
      delete [] _xPosHere;
      
    } // end loop over all track candidates

#ifdef MARLIN_USE_AIDA

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
    for (int help = 0; help < 500; help++) {
      delete [] _zPos[help];
      delete [] _yPos[help];
      delete [] _xPos[help];
    }

    delete [] _zPos;
    delete [] _yPos;
    delete [] _xPos;
    
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

}

void EUTelMultiLineFit::bookHistos() {
  
  
#ifdef MARLIN_USE_AIDA
  
  try {
    streamlog_out ( MESSAGE2 ) << "Booking histograms" << endl;
    
    const int    NBin = 10000;
    const double Chi2Min  = 0.;      
    const double Chi2Max  = 10000.;      
    const double Min  = -5000.;
    const double Max  = 5000.;
    const double angleMin  = -3.;
    const double angleMax  = 3.;
    const double tracksMin = -0.5;
    const double tracksMax = 19.5;

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

  
#endif
