// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Philip Roloff, DESY <mailto:philipp.roloff@desy.de>
// Version: $Id: EUTelAlign.cc,v 1.14 2007-10-03 22:31:20 roloff Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// build only if ROOT is used
#ifdef MARLIN_USE_ROOT

// built only if GEAR is used
#ifdef USE_GEAR

// eutelescope includes ".h" 
#include "EUTelAlign.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelExceptions.h"

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
#ifdef MARLIN_USE_AIDA
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
			     _alignedPlane, static_cast < int > (0));
  
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

  _nPlanes = _siPlanesParameters->getSiPlanesNumber();
  _intrResol = new double[_nPlanes];
  _xMeasPos = new double[_nPlanes];
  _yMeasPos = new double[_nPlanes];
  _zMeasPos = new double[_nPlanes];
    
}

void EUTelAlign::processRunHeader (LCRunHeader * rdr) {

  auto_ptr<EUTelRunHeaderImpl> header( new EUTelRunHeaderImpl( rdr ) );


  // the run header contains the number of detectors. This number
  // should be in principle the same as the number of layers in the
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
  
  if ( header->getGeoID() != _siPlanesParameters->getSiPlanesID() ) {
    message<ERROR> ( "Error during the geometry consistency check: " );
    message<ERROR> ( log() << "The run header says the GeoID is " << header->getGeoID() );
    message<ERROR> ( log() << "The GEAR description says is     " << _siPlanesParameters->getSiPlanesNumber() );
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

void EUTelAlign::processEvent (LCEvent * event) {

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event) ;
  
  if ( evt->getEventType() == kEORE ) {
    message<DEBUG> ( "EORE found: nothing else to do." );
    return;
  }

  int nHitsFirstPlane = 0;
  int nHitsSecondPlane = 0;
  
  try {

    LCCollectionVec * measHitCollection = static_cast<LCCollectionVec*> (event->getCollection( _measHitCollectionName ));
  
    int detectorID    = -99; // it's a non sense
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

    // Loop over all hits
    for ( int iHit = 0; iHit < measHitCollection->getNumberOfElements(); iHit++ ) {
    
      TrackerHitImpl * measHit = static_cast<TrackerHitImpl*> ( measHitCollection->getElementAt(iHit) );
    
      LCObjectVec clusterVector = measHit->getRawHits();
      
      EUTelVirtualCluster * cluster;
      if ( measHit->getType() == kEUTelFFClusterImpl ) {
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
	_intrResol[layerIndex] = 1000*_siPlanesLayerLayout->getSensitiveResolution(layerIndex); //um
      
      }
      
      // If first plane: calculate extrapolation to second plane

      layerIndex = _conversionIdMap[detectorID];

      if (layerIndex == 0 && nHitsFirstPlane < 20) {
      
	allHitsFirstLayerMeasuredX[nHitsFirstPlane] = 1000 * measHit->getPosition()[0];
	allHitsFirstLayerMeasuredY[nHitsFirstPlane] = 1000 * measHit->getPosition()[1];
	allHitsFirstLayerMeasuredZ[nHitsFirstPlane] = 1000 * measHit->getPosition()[2];

	// cout << hitsForFit.firstLayerMeasuredX << " " << hitsForFit.firstLayerMeasuredY << "            ";

	allHitsFirstLayerResolution[nHitsFirstPlane] = 1000 * _siPlanesLayerLayout->getSensitiveResolution(layerIndex); // Add multiple scattering later!

	// hitsForFit.secondLayerPredictedX = hitsForFit.firstLayerMeasuredX;
	// hitsForFit.secondLayerPredictedY = hitsForFit.firstLayerMeasuredY;

	nHitsFirstPlane++;

      } else if (layerIndex == (_alignedPlane - 1) && nHitsSecondPlane < 20) {

	allHitsSecondLayerMeasuredX[nHitsSecondPlane] = 1000 * measHit->getPosition()[0];
	allHitsSecondLayerMeasuredY[nHitsSecondPlane] = 1000 * measHit->getPosition()[1];
	allHitsSecondLayerMeasuredZ[nHitsSecondPlane] = 1000 * measHit->getPosition()[2];

	// cout << hitsForFit.secondLayerMeasuredX << " " << hitsForFit.secondLayerMeasuredY << endl;
	
	// hitsForFit.secondLayerPredictedZ = 1000 * measHit->getPosition()[2];
	
	allHitsSecondLayerResolution[nHitsSecondPlane] = 1000 * _siPlanesLayerLayout->getSensitiveResolution(layerIndex); // Add multiple scattering later!

	nHitsSecondPlane++;

      }

      delete cluster; // <--- destroying the cluster   

    } // End loop over all hits

    // check number of hits
    if (nHitsFirstPlane <= _nHitsMax && nHitsSecondPlane <= _nHitsMax) {

      // loop over all hits in first plane
      for (int firsthit = 0; firsthit < nHitsFirstPlane; firsthit++) {

	int take = -1000;
	int veto = -1000;
      
	// loop over all hits in second plane
	for (int secondhit = 0; secondhit < nHitsSecondPlane; secondhit++) {
	
	  // calculate distance between hits
	  double distance = sqrt(pow(allHitsFirstLayerMeasuredX[firsthit] - allHitsSecondLayerMeasuredX[secondhit],2) + pow(allHitsFirstLayerMeasuredY[firsthit] - allHitsSecondLayerMeasuredY[secondhit],2));

#ifdef MARLIN_USE_AIDA
    
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

    } // end if check number of hits

      // _hitsForFit.push_back(hitsForFit);
  
  } catch (DataNotAvailableException& e) {
    streamlog_out  ( WARNING2 ) <<  "No input collection found on event " << event->getEventNumber() << " in run " << event->getRunNumber() << endl;
  }

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
  for (uint i = 0; i < _hitsForFit.size(); i++) {

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

  // call migrad (500 iterations, 0.1 = tolerance)
  arglist[0] = 500;
  arglist[1] = 0.1;
  gMinuit->mnexcm("MIGRAD",arglist,1,ierflag);

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
  for (uint i = 0; i < _hitsForFit.size(); i++) {
  
    residual_x_simple = off_x_simple + _hitsForFit[i].secondLayerMeasuredX - _hitsForFit[i].secondLayerPredictedX;
    residual_y_simple = off_y_simple + _hitsForFit[i].secondLayerMeasuredY - _hitsForFit[i].secondLayerPredictedY;

#ifdef MARLIN_USE_AIDA
    
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

  // call migrad (500 iterations, 0.1 = tolerance)
  arglist[0] = 500;
  arglist[1] = 0.1;
  gMinuit->mnexcm("MIGRAD",arglist,1,ierflag);

  streamlog_out ( MESSAGE2 ) << endl << "Third iteration of alignment: include chi^2 cut" << endl;
  streamlog_out ( MESSAGE2 ) << "-----------------------------------------------" << endl << endl;

  // release chi2
  gMinuit->mnparm(5,"chi2",start_chi2,1,0,0,ierflag);

  // call migrad (500 iterations, 0.1 = tolerance)
  arglist[0] = 500;
  arglist[1] = 0.1;
  gMinuit->mnexcm("MIGRAD",arglist,1,ierflag);

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
  for (uint i = 0; i < _hitsForFit.size(); i++) {

    x = (cos(theta_y)*cos(theta_z)) * _hitsForFit[i].secondLayerMeasuredX + ((-1)*sin(theta_x)*sin(theta_y)*cos(theta_z) + cos(theta_x)*sin(theta_z)) * _hitsForFit[i].secondLayerMeasuredY + off_x;
    y = ((-1)*cos(theta_y)*sin(theta_z)) * _hitsForFit[i].secondLayerMeasuredX + (sin(theta_x)*sin(theta_y)*sin(theta_z) + cos(theta_x)*cos(theta_z)) * _hitsForFit[i].secondLayerMeasuredY + off_y;

    residual_x = x - _hitsForFit[i].secondLayerPredictedX;
    residual_y = y - _hitsForFit[i].secondLayerPredictedY;

#ifdef MARLIN_USE_AIDA
    
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

//   delete [] _intrResolY;
//   delete [] _intrResolX;
  delete [] _xMeasPos;
  delete [] _yMeasPos;
  delete [] _zMeasPos;

  streamlog_out ( MESSAGE2 ) << "Successfully finished" << endl;

}

void EUTelAlign::bookHistos() {

#ifdef MARLIN_USE_AIDA
  
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

#endif

#endif
