// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelHistogramMaker.cc,v 1.14 2007-07-26 06:49:20 bulgheroni Exp $
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
#include "EUTelHistogramManager.h"
#include "EUTelMatrixDecoder.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/AIDAProcessor.h"

// lcio includes <.h> 
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <Exceptions.h>

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
std::string EUTelHistogramMaker::_clusterSignalHistoName      = "clusterSignal";
std::string EUTelHistogramMaker::_seedSignalHistoName         = "seedSignal";
std::string EUTelHistogramMaker::_hitMapHistoName             = "hitMap";
std::string EUTelHistogramMaker::_seedSNRHistoName            = "seedSNR";
std::string EUTelHistogramMaker::_clusterNoiseHistoName       = "clusterNoise";
std::string EUTelHistogramMaker::_clusterSNRHistoName         = "clusterSNR";
std::string EUTelHistogramMaker::_eventMultiplicityHistoName  = "eventMultiplicity";
#endif

EUTelHistogramMaker::EUTelHistogramMaker () : Processor("EUTelHistogramMaker") {

  // modify processor description
  _description =
    "EUTelHistogramMaker fills reference and control histograms";

  // first of all we need to register the input collection
  registerInputCollection (LCIO::TRACKERPULSE, "PulseCollectionName",
			   "Input tracker pulse collection",
			   _pulseCollectionName, string("pulse"));

  registerProcessorParameter("HistoInfoFileName", "This is the name of the histogram information file",
			     _histoInfoFileName, string( "histoinfo.xml" ) );

  IntVec clusterNxNExample;
  clusterNxNExample.push_back(3);
  clusterNxNExample.push_back(5);
  
  registerOptionalParameter("ClusterNxN", "The list of cluster NxN to be filled."
			    "For example 3 means filling the 3x3 histogram spectrum",
			    _clusterSpectraNxNVector, clusterNxNExample);

  IntVec clusterNExample;
  clusterNExample.push_back(4);
  clusterNExample.push_back(9);
  clusterNExample.push_back(14);
  clusterNExample.push_back(19);
  clusterNExample.push_back(25);
  registerOptionalParameter("ClusterN", "The list of cluster N to be filled."
			    "For example 7 means filling the cluster spectra with the 7 most significant pixels",
			    _clusterSpectraNVector, clusterNExample );

  _isFirstEvent = true;

  registerOptionalParameter("StatusCollectionName","The name of the status collection.\n"
			    "Needed to fill in noise related histograms",
			    _statusCollectionName, string( "status" ) );

  registerOptionalParameter("NoiseCollectionName","The name of the noise collection.\n"
			    "Needed to fill in noise related histograms",
			    _noiseCollectionName, string( "noise" ) );

}


void EUTelHistogramMaker::init () {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  _iRun = 0;
  _iEvt = 0;

  
  // by default fill also the noise related histograms.
  _noiseHistoSwitch = true;
    

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

  ++_iRun;
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

  try {

    LCCollectionVec * pulseCollectionVec = dynamic_cast<LCCollectionVec*>  (evt->getCollection(_pulseCollectionName));
    CellIDDecoder<TrackerPulseImpl> cellDecoder(pulseCollectionVec);
    
    LCCollectionVec * noiseCollectionVec = 0x0, * statusCollectionVec = 0x0;
    
    if ( _noiseHistoSwitch ) {
      try {
	noiseCollectionVec  = dynamic_cast<LCCollectionVec *> ( evt->getCollection( _noiseCollectionName ) );
      } catch (lcio::DataNotAvailableException& e) {
	message<ERROR> ( log() << e.what() 
			 << "Switching off the noise histogram filling and continuing" );
	_noiseHistoSwitch &= false;    
      }
      try {
	statusCollectionVec = dynamic_cast<LCCollectionVec *> ( evt->getCollection( _statusCollectionName ) );
      } catch (lcio::DataNotAvailableException& e) {
	message<ERROR> ( log() << e.what() 
			 << "Switching off the noise histogram filling and continuing" );    
	_noiseHistoSwitch &= false;
      }
    }
    
    if ( isFirstEvent() ) {
      bookHistos();
      _isFirstEvent = false;
    } 
    
    // prepare and reset the hit counter
    vector<unsigned short> eventCounterVec( _noOfDetector, 0 );
    
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
      
      // increment of one unit the event counter for this plane
      eventCounterVec[detectorID]++;
      
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
      
      if ( _noiseHistoSwitch ) {
	// get the noise TrackerDataImpl corresponding to the detector
	// under analysis and the status matrix as well
	TrackerDataImpl    * noiseMatrix  = dynamic_cast<TrackerDataImpl *>    (noiseCollectionVec->getElementAt(detectorID) );
	TrackerRawDataImpl * statusMatrix = dynamic_cast<TrackerRawDataImpl *> (statusCollectionVec->getElementAt(detectorID) );
	
	EUTelMatrixDecoder noiseMatrixDecoder( _maxX[detectorID] - _minX[detectorID] + 1,
					       _maxY[detectorID] - _minY[detectorID] + 1,
					       _minX[detectorID] , _minY[detectorID] );
	
	int xClusterSize, yClusterSize;
	cluster->getClusterSize(xClusterSize, yClusterSize);
	vector<float > noiseValues;
	for ( int yPixel = ySeed - ( yClusterSize / 2 ); yPixel <= ySeed + ( yClusterSize / 2 ); yPixel++ ) {
	  for ( int xPixel = xSeed - ( xClusterSize / 2 ); xPixel <= xSeed + ( xClusterSize / 2 ); xPixel++ ) {
	    
	    // always check we are still within the sensor!!!
	    if ( ( xPixel >= _minX[detectorID] )  &&  ( xPixel <= _maxX[detectorID] ) &&
		 ( yPixel >= _minY[detectorID] )  &&  ( yPixel <= _maxY[detectorID] ) ) {
	      int index = noiseMatrixDecoder.getIndexFromXY(xPixel, yPixel);
	      
	      // the corresponding position in the status matrix has to be HITPIXEL
	      // in the EUTelClusteringProcessor, we verify also that
	      // the pixel isHit, but this cannot be done in this
	      // processor, since the status matrix could have been reset
	      // 
	      // bool isHit  = ( statusMatrix->getADCValues()[index] ==
	      // EUTELESCOPE::HITPIXEL );
	      //
	      bool isBad  = ( statusMatrix->getADCValues()[index] == EUTELESCOPE::BADPIXEL );
	      if ( !isBad ) {
		noiseValues.push_back( noiseMatrix->getChargeValues()[index] );
	      } else {
		noiseValues.push_back( 0. );
	      }
	    } else {
	      noiseValues.push_back( 0. );
	    }
	  }
	}
	try {
	  cluster->setNoiseValues( noiseValues );
	} catch ( IncompatibleDataSetException& e ) {
	  message<ERROR> ( log() << e.what() << "\n"
			   "Continuing without filling the noise histograms" );
	  _noiseHistoSwitch = false;
	}
      }
    
      if ( _noiseHistoSwitch ) {
	AIDA::IHistogram1D * histo;
	
	{
	  stringstream ss;
	  ss << _clusterNoiseHistoName << "-d" << detectorID;
	  tempHistoName = ss.str();
	}
	histo = dynamic_cast<AIDA::IHistogram1D* > ( _aidaHistoMap[tempHistoName] );
	if ( histo ) {
	  histo->fill( cluster->getClusterNoise() );
	}
	
	
	{
	  stringstream ss;
	  ss << _clusterSNRHistoName << "-d" << detectorID;
	  tempHistoName = ss.str();
	}
	histo = dynamic_cast<AIDA::IHistogram1D* > ( _aidaHistoMap[tempHistoName] );
	if ( histo ) {
	  histo->fill( cluster->getClusterSNR() );
	}
	
	{
	  stringstream ss;
	  ss << _seedSNRHistoName << "-d" << detectorID;
	  tempHistoName = ss.str();
	}
	histo = dynamic_cast<AIDA::IHistogram1D * > ( _aidaHistoMap[tempHistoName] );
	if ( histo ) {
	  histo->fill( cluster->getSeedSNR() );
	}      
	
	vector<int >::iterator iter = _clusterSpectraNxNVector.begin();
	while ( iter != _clusterSpectraNxNVector.end() ) {
	  {
	    stringstream ss;
	    ss << _clusterSNRHistoName << (*iter) << "x" << (*iter) << "-d" << detectorID;
	    tempHistoName = ss.str();
	  }
	  histo = dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName] ) ;
	  if ( histo ) {
	    histo->fill(cluster->getClusterSNR( (*iter), (*iter) ));
	  }
	  ++iter;
	}
	
	vector<float > snrs = cluster->getClusterSNR(_clusterSpectraNVector);
	for ( unsigned int i = 0; i < snrs.size() ; i++ ) {
	  {
	    stringstream ss;
	    ss << _clusterSNRHistoName << _clusterSpectraNVector[i] << "-d" << detectorID;
	    tempHistoName = ss.str();
	  }
	  histo = dynamic_cast<AIDA::IHistogram1D * > ( _aidaHistoMap[tempHistoName] ) ;
	  if ( histo ) {
	    histo->fill( snrs[i] ); 
	  }
	}
      }
      
      delete cluster;
    }

    // fill the event multiplicity here
    string tempHistoName;
    for ( int iDetector = 0; iDetector < _noOfDetector; iDetector++ ) {
      {
	stringstream ss;
	ss << _eventMultiplicityHistoName << "-d" << iDetector;
	tempHistoName = ss.str();
      }
      AIDA::IHistogram1D * histo = dynamic_cast<AIDA::IHistogram1D *> ( _aidaHistoMap[tempHistoName] );
      if ( histo ) {
	histo->fill( eventCounterVec[iDetector] );
      }
    }
    
  } catch( DataNotAvailableException& e ) {
    message<WARNING> ( log() << "No input collection found on event " << _iEvt );
  }
#endif

  ++_iEvt;
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
  
  auto_ptr<EUTelHistogramManager> histoMgr( new EUTelHistogramManager( _histoInfoFileName ));
  EUTelHistogramInfo    * histoInfo;
  bool                    isHistoManagerAvailable;

  try {
    isHistoManagerAvailable = histoMgr->init();
  } catch ( ios::failure& e) {
    message<ERROR> ( log() << "I/O problem with " << _histoInfoFileName << "\n"
		     << "Continuing without histogram manager"    );
    isHistoManagerAvailable = false;
  } catch ( ParseException& e ) {
    message<ERROR> ( log() << e.what() << "\n"
		     << "Continuing without histogram manager" );
    isHistoManagerAvailable = false;
  }

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
    
    int    clusterNBin  = 1000;
    double clusterMin   = 0.;
    double clusterMax   = 1000.;
    string clusterTitle = "Cluster spectrum with all pixels";
    if ( isHistoManagerAvailable ) {
      histoInfo = histoMgr->getHistogramInfo(_clusterSignalHistoName);
      if ( histoInfo ) {
	message<DEBUG> ( log() << (* histoInfo ) );
	clusterNBin = histoInfo->_xBin;
	clusterMin  = histoInfo->_xMin;
	clusterMax  = histoInfo->_xMax;
	if ( histoInfo->_title != "" ) clusterTitle = histoInfo->_title;
      }
    }
    AIDA::IHistogram1D * clusterSignalHisto = 
      AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), 
								clusterNBin,clusterMin,clusterMax);
    _aidaHistoMap.insert(make_pair(tempHistoName, clusterSignalHisto));
    clusterSignalHisto->setTitle(clusterTitle.c_str());

    int    clusterSNRNBin  = 300;
    double clusterSNRMin   = 0.;
    double clusterSNRMax   = 200;
    string clusterSNRTitle = "Cluster SNR";
    if ( isHistoManagerAvailable ) {
      histoInfo = histoMgr->getHistogramInfo( _clusterSNRHistoName );
      if ( histoInfo ) {
	message<DEBUG> ( log() << (* histoInfo ) );
	clusterSNRNBin = histoInfo->_xBin;
	clusterSNRMin  = histoInfo->_xMin;
	clusterSNRMax  = histoInfo->_xMax;
	if ( histoInfo->_title != "" ) clusterSNRTitle = histoInfo->_title;
      }
    }

    if ( _noiseHistoSwitch ) {
      // cluster SNR
      {
	stringstream ss;
	ss << _clusterSNRHistoName << "-d" << iDetector;
	tempHistoName = ss.str();
      }
      AIDA::IHistogram1D * clusterSNRHisto = 
	AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
								  clusterSNRNBin, clusterSNRMin, clusterSNRMax);
      _aidaHistoMap.insert( make_pair(tempHistoName, clusterSNRHisto) ) ;
      clusterSNRHisto->setTitle(clusterSNRTitle.c_str());
    }
    

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
      string tempTitle;
      {
	stringstream ss;
	ss <<  "Cluster spectrum with the " << (*iter) << " most significant pixels ";
	tempTitle=ss.str();
      }
      clusterSignalNHisto->setTitle(tempTitle.c_str());

      if ( _noiseHistoSwitch ) {
	// this is for the SNR
	{
	  stringstream ss;
	  ss << _clusterSNRHistoName << (*iter) << "-d" << iDetector;
	  tempHistoName = ss.str();
	}
	AIDA::IHistogram1D * clusterSNRNHisto = 
	  AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
								    clusterSNRNBin, clusterSNRMin, clusterSNRMax);
	{
	  stringstream ss;
	  ss << "Cluster SNR with the " << (*iter) << " most significant pixels";
	  tempTitle = ss.str();
	} 
	_aidaHistoMap.insert( make_pair( tempHistoName, clusterSNRNHisto ) );
	clusterSNRNHisto->setTitle(tempTitle.c_str());
      }
      
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

      if ( _noiseHistoSwitch ) {
	// then the SNR
	{
	  stringstream ss;
	  ss << _clusterSNRHistoName << (*iter) << "x" << (*iter) << "-d" << iDetector;
	  tempHistoName = ss.str();
	}
	AIDA::IHistogram1D * clusterSNRNxNHisto =
	  AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
								    clusterSNRNBin, clusterSNRMin, clusterSNRMax);
	_aidaHistoMap.insert(make_pair(tempHistoName, clusterSNRNxNHisto) );
	{
	  stringstream ss;
	  ss << "SNR with " << (*iter) << " by " << (*iter) << " pixels ";
	  tempTitle = ss.str();
	}
	clusterSNRNxNHisto->setTitle(tempTitle.c_str());
      }
      ++iter;
    }

    {
      stringstream ss;
      ss << _seedSignalHistoName << "-d" << iDetector;
      tempHistoName = ss.str();
    } 

    int    seedNBin  = 500;
    double seedMin   = 0.;
    double seedMax   = 500.;
    string seedTitle = "Seed pixel spectrum";
    if ( isHistoManagerAvailable ) {
      histoInfo = histoMgr->getHistogramInfo( _seedSignalHistoName );
      if ( histoInfo ) {
	message<DEBUG> ( log() << (* histoInfo ) );
	seedNBin = histoInfo->_xBin;
	seedMin  = histoInfo->_xMin;
	seedMax  = histoInfo->_xMax;
	if ( histoInfo->_title != "" ) seedTitle = histoInfo->_title;
      }
    }
    AIDA::IHistogram1D * seedSignalHisto =
      AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
								seedNBin, seedMin, seedMax);
    _aidaHistoMap.insert(make_pair(tempHistoName, seedSignalHisto));
    seedSignalHisto->setTitle(seedTitle.c_str());
  
    if ( _noiseHistoSwitch ) {
      // seed SNR
      {
	stringstream ss;
	ss << _seedSNRHistoName << "-d" << iDetector;
	tempHistoName = ss.str();
      }
      int    seedSNRNBin  =  300;
      double seedSNRMin   =    0.;
      double seedSNRMax   =  200.;
      string seedSNRTitle = "Seed SNR";
      if ( isHistoManagerAvailable ) {
	histoInfo = histoMgr->getHistogramInfo( _seedSNRHistoName );
	if ( histoInfo ) {
	  message<DEBUG> ( log() << (* histoInfo ) );
	  seedSNRNBin = histoInfo->_xBin;
	  seedSNRMin  = histoInfo->_xMin;
	  seedSNRMax  = histoInfo->_xMax;
	  if ( histoInfo->_title != "" ) seedSNRTitle = histoInfo->_title;
	}
      }
      AIDA::IHistogram1D * seedSNRHisto = 
	AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
								  seedSNRNBin, seedSNRMin, seedSNRMax);
      _aidaHistoMap.insert( make_pair(tempHistoName, seedSNRHisto ));
      seedSNRHisto->setTitle(seedSNRTitle.c_str());    
      
      {
	stringstream ss;
	ss << _clusterNoiseHistoName << "-d" << iDetector;
	tempHistoName = ss.str();
      }

      // cluster noise
      int    clusterNoiseNBin  =  300;
      double clusterNoiseMin   =    0.;
      double clusterNoiseMax   =  200.;
      string clusterNoiseTitle = "Cluster noise";
      if ( isHistoManagerAvailable ) {
	histoInfo = histoMgr->getHistogramInfo( _clusterNoiseHistoName );
	if ( histoInfo ) {
	  message<DEBUG> ( log() << (* histoInfo ) );
	  clusterNoiseNBin = histoInfo->_xBin;
	  clusterNoiseMin  = histoInfo->_xMin;
	  clusterNoiseMax  = histoInfo->_xMax;
	  if ( histoInfo->_title != "" ) clusterNoiseTitle = histoInfo->_title;
	}
      }
      AIDA::IHistogram1D * clusterNoiseHisto = 
	AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
								  clusterNoiseNBin, clusterNoiseMin, clusterNoiseMax);
      _aidaHistoMap.insert( make_pair(tempHistoName, clusterNoiseHisto ));
      clusterNoiseHisto->setTitle(clusterNoiseTitle.c_str());
      
    }

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

  
    {
      stringstream ss;
      ss << _eventMultiplicityHistoName << "-d" << iDetector;
      tempHistoName = ss.str();
    }
    int     eventMultiNBin  = 30;
    double  eventMultiMin   =  0.;
    double  eventMultiMax   = 30.;
    string  eventMultiTitle = "Event multiplicity";
    if ( isHistoManagerAvailable ) {
      histoInfo = histoMgr->getHistogramInfo(  _eventMultiplicityHistoName );
      if ( histoInfo ) {
	message<DEBUG> ( log() << (* histoInfo ) );
	eventMultiNBin  = histoInfo->_xBin;
	eventMultiMin   = histoInfo->_xMin;
	eventMultiMax   = histoInfo->_xMax;
	if ( histoInfo->_title != "" ) eventMultiTitle = histoInfo->_title;
      }
    }
    AIDA::IHistogram1D * eventMultiHisto =
      AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
								eventMultiNBin, eventMultiMin, eventMultiMax);
    _aidaHistoMap.insert( make_pair(tempHistoName, eventMultiHisto) );
    eventMultiHisto->setTitle( eventMultiTitle.c_str() );
  }

#else
  message<MESSAGE> ( log() << "No histogram produced because Marlin doesn't use AIDA" );
#endif

}

