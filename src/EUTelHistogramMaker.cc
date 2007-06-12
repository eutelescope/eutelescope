// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelHistogramMaker.cc,v 1.1 2007-06-12 13:56:08 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h" 
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelHistogramMaker.h"
#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"

// marlin includes ".h"
#include <marlin/Processor.h>
#include <marlin/AIDAProcessor.h>

#ifdef MARLIN_USE_AIDA
// aida includes <.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h> 

// system includes <>
#include <string>
#include <sstream>
#include <vector>


using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;


// definition of static members mainly used to name histograms
#ifdef MARLIN_USE_AIDA
std::string EUTelHistogramMaker::_clusterSignalHistoName   = "clusterSignal";
std::string EUTelHistogramMaker::_seedSignalHistoName      = "seedSignal";
std::string EUTelHistogramMaker::_hitMapHistoName          = "hitMap";
#endif

EUTelHistogramMaker::EUTelHistogramMaker () : Processor("EUTelHistogramMaker") {

  // modify processor description
  _description =
    "EUTelHistogramMaker fills reference and control histograms";

  // first of all we need to register the input collection
  registerInputCollection (LCIO::TRACKERPULSE, "PulseCollectionName",
			   "Input tracker pulse collection",
			   _pulseCollectionName, string("pulse"));

  IntVec clusterNxNExample;
  clusterNxNExample.push_back(3);
  clusterNxNExample.push_back(5);
  
  registerOptionalParameter("ClusterNxN", "The list of cluster NxN to be filled."
			    "For example 3 means filling the 3x3 histogram spectrum",
			    _clusterSpectraNxNVector, clusterNxNExample, clusterNxNExample.size());

  IntVec clusterNExample;
  clusterNExample.push_back(4);
  clusterNExample.push_back(9);
  clusterNExample.push_back(14);
  clusterNExample.push_back(19);
  clusterNExample.push_back(25);
  registerOptionalParameter("ClusterN", "The list of cluster N to be filled."
			    "For example 7 means filling the cluster spectra with the 7 most significant pixels",
			    _clusterSpectraNVector, clusterNExample, clusterNExample.size() );
  
  

}


void EUTelHistogramMaker::init () {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  _iEvt = 0;

}

void EUTelHistogramMaker::processRunHeader (LCRunHeader * rdr) {

  // to make things easier re-cast the input header to the EUTelRunHeaderImpl
  EUTelRunHeaderImpl *  runHeader = static_cast<EUTelRunHeaderImpl*>(rdr);

  // let me get from the run header all the available parameter
  _noOfDetector = runHeader->getNoOfDetector();

  // now the four vectors containing the first and the last pixel
  // along both the directions
  _minX = runHeader->getMinX();
  _maxX = runHeader->getMaxX();
  _minY = runHeader->getMinY();
  _maxY = runHeader->getMaxY();

#ifdef MARLIN_USE_AIDA
  bookHistos();
#endif

}


void EUTelHistogramMaker::processEvent (LCEvent * evt) {

  EUTelEventImpl * eutelEvent = static_cast<EUTelEventImpl*> (evt);
  EventType type              = eutelEvent->getEventType();
  
  if ( type == kEORE ) {
    message<DEBUG> ( "EORE found: nothing else to do.");
    return ;
  } else if ( type == kUNKNOWN ) {
    message<WARNING> ( log() << "Event number " << evt->getEventNumber() 
		       << " is of unknown type. Continue considering it as a normal Data Event."  );
  }

#ifdef MARLIN_USE_AIDA

  if ( (_iEvt % 10) == 0 ) 
    message<MESSAGE> ( log() << "Filling histogram on event " << _iEvt );
  
  LCCollectionVec * pulseCollectionVec = dynamic_cast<LCCollectionVec*> 
    (evt->getCollection(_pulseCollectionName));
  CellIDDecoder<TrackerPulseImpl> cellDecoder(pulseCollectionVec);

  for ( int iPulse = 0; iPulse < pulseCollectionVec->getNumberOfElements(); iPulse++ ) {
    TrackerPulseImpl * pulse = dynamic_cast<TrackerPulseImpl*> ( pulseCollectionVec->getElementAt(iPulse) );
    ClusterType        type  = static_cast<ClusterType> ( static_cast<int> ( cellDecoder(pulse)["type"] ));
    
    EUTelVirtualCluster * cluster;
    
    if ( type == kEUTelFFClusterImpl ) 
      cluster = new EUTelFFClusterImpl ( static_cast<TrackerDataImpl*> ( pulse->getTrackerData() ) );
    else {
      message<ERROR> ( "Unknown cluster type. Sorry for quitting" );
      throw UnknownDataTypeException("Cluster type unknown");
    }

    int detectorID = cluster->getDetectorID();
    string tempHistoName;

    {
      stringstream ss;
      ss << _clusterSignalHistoName << "-d" << detectorID;
      tempHistoName = ss.str();
    } 
    (dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(cluster->getTotalCharge());
    
    {
      stringstream ss;
      ss << _seedSignalHistoName << "-d" << detectorID;
      tempHistoName = ss.str();
    }
    (dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(cluster->getSeedCharge());

    vector<int >::iterator iter = _clusterSpectraNVector.begin();
    while ( iter != _clusterSpectraNVector.end() ) {
      {
	stringstream ss;
	ss << _clusterSignalHistoName << (*iter) << "-d" << detectorID;
	tempHistoName = ss.str();
      }
      (dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))
	->fill(cluster->getClusterCharge((*iter)));
      ++iter;
    }

    iter = _clusterSpectraNxNVector.begin();
    while ( iter != _clusterSpectraNxNVector.end() ) {
      {
	stringstream ss;
	ss << _clusterSignalHistoName << (*iter) << "x" << (*iter) << "-d" << detectorID;
	tempHistoName = ss.str();
      }
      (dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))
	->fill(cluster->getClusterCharge((*iter), (*iter)));
      ++iter;
    }

    {
      stringstream ss;
      ss << _hitMapHistoName << "-d" << detectorID;
      tempHistoName = ss.str();
    } 
    int xSeed, ySeed;
    cluster->getSeedCoord(xSeed, ySeed);
    (dynamic_cast<AIDA::IHistogram2D*> (_aidaHistoMap[tempHistoName]))
      ->fill(static_cast<double >(xSeed), static_cast<double >(ySeed), 1.);
    
    
  }

#endif
  
}

void EUTelHistogramMaker::check (LCEvent * evt) {
  // nothing to check here - could be used to fill check plots in reconstruction processor
}


void EUTelHistogramMaker::end() {

  message<MESSAGE> ( "Processor finished successfully." );

}

void EUTelHistogramMaker::bookHistos() {

#ifdef MARLIN_USE_AIDA
  // histograms are grouped in loops and detectors
  message<MESSAGE> ( log() << "Booking histograms " );


  string tempHistoName;
  string basePath;
  for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
    
    {
      stringstream ss;
      ss << "detector-" << iDetector;
      basePath = ss.str();
    }
    AIDAProcessor::tree(this)->mkdir(basePath.c_str());
    basePath.append("/");

    {
      stringstream ss;
      ss << _clusterSignalHistoName << "-d" << iDetector;
      tempHistoName = ss.str();
    } 

    const int    clusterNBin = 1000;
    const double clusterMin  = 0.;
    const double clusterMax  = 1000.;

    AIDA::IHistogram1D * clusterSignalHisto = 
      AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), 
								clusterNBin,clusterMin,clusterMax);
    _aidaHistoMap.insert(make_pair(tempHistoName, clusterSignalHisto));
    clusterSignalHisto->setTitle("Cluster spectrum with all pixels");

    
    vector<int >::iterator iter = _clusterSpectraNVector.begin();
    while ( iter != _clusterSpectraNVector.end() ) {
      {
	stringstream ss;
	ss << _clusterSignalHistoName << (*iter) << "-d" << iDetector;
	tempHistoName = ss.str();
      }
      AIDA::IHistogram1D * clusterSignalNHisto = 
	AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
								  clusterNBin, clusterMin, clusterMax);
      _aidaHistoMap.insert(make_pair(tempHistoName, clusterSignalNHisto) );
      string tempTitle = "Cluster spectrum with the " + (*iter);
      tempTitle.append(" most significant pixels ");
      clusterSignalNHisto->setTitle(tempTitle.c_str());

      ++iter;
    }

    iter = _clusterSpectraNxNVector.begin();
    while ( iter != _clusterSpectraNxNVector.end() ) {
      {
	stringstream ss;
	ss << _clusterSignalHistoName << (*iter) << "x" << (*iter) << "-d" << iDetector;
	tempHistoName = ss.str();
      }
      AIDA::IHistogram1D * clusterSignalNxNHisto = 
	AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
								  clusterNBin, clusterMin, clusterMax);
      _aidaHistoMap.insert(make_pair(tempHistoName, clusterSignalNxNHisto) );
      string tempTitle;
      {
	stringstream ss;
	ss << "Cluster spectrum with " << (*iter) << " by " << (*iter) << " pixels ";
	tempTitle = ss.str();
      }
      clusterSignalNxNHisto->setTitle(tempTitle.c_str());

      ++iter;
    }


    int    seedNBin = 500;
    double seedMin  = 0.;
    double seedMax  = 500.;

    {
      stringstream ss;
      ss << _seedSignalHistoName << "-d" << iDetector;
      tempHistoName = ss.str();
    } 
    AIDA::IHistogram1D * seedSignalHisto =
      AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
								seedNBin, seedMin, seedMax);
    _aidaHistoMap.insert(make_pair(tempHistoName, seedSignalHisto));
    seedSignalHisto->setTitle("Seed pixel spectrum");
  
    
    {
      stringstream ss;
      ss << _hitMapHistoName << "-d" << iDetector;
      tempHistoName = ss.str();
    } 
    int     xBin = _maxX[iDetector] - _minX[iDetector] + 1;
    double  xMin = static_cast<double >(_minX[iDetector]) - 0.5;
    double  xMax = static_cast<double >(_maxX[iDetector]) + 0.5;
    int     yBin = _maxY[iDetector] - _minY[iDetector] + 1;
    double  yMin = static_cast<double >(_minY[iDetector]) - 0.5;
    double  yMax = static_cast<double >(_maxY[iDetector]) + 0.5;
    AIDA::IHistogram2D * hitMapHisto = 
      AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(),
							       xBin, xMin, xMax,yBin, yMin, yMax);
    _aidaHistoMap.insert(make_pair(tempHistoName, hitMapHisto));
    hitMapHisto->setTitle("Hit map");

  }
  
  
#else
  message<MESSAGE> ( log() << "No histogram produced because Marlin doesn't use AIDA" );
#endif

  ++_iEvt;
}

