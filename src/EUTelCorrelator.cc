// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Silvia Bonfanti, Uni. Insubria  <mailto:silviafisica@gmail.com>
// Author Loretta Negrini, Uni. Insubria  <mailto:loryneg@gmail.com>
// Version $Id: EUTelCorrelator.cc,v 1.1 2008-07-11 15:35:30 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h" 
#include "EUTelCorrelator.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelExceptions.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// aida includes <.h>
#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram2D.h>
#endif

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

using namespace std;
using namespace marlin;
using namespace eutelescope;

// definition of static members mainly used to name histograms
#ifdef MARLIN_USE_AIDA
std::string EUTelCorrelator::_clusterXCorrelationHistoName   = "ClusterXCorrelationHisto";
#endif

EUTelCorrelator::EUTelCorrelator () : Processor("EUTelCorrelator") {

  // modify processor description
  _description =
    "EUTelCorrelator fills histograms with correlation plots";

  registerInputCollection(LCIO::TRACKERPULSE,"InputCollectionName",
			  "Cluster (pulse) collection name",
			  _inputCollectionName, string ( "cluster" ));
  
  
}


void EUTelCorrelator::init() {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;


}

void EUTelCorrelator::processRunHeader (LCRunHeader * rdr) {

  EUTelRunHeaderImpl * runHeader = new EUTelRunHeaderImpl( rdr ) ; 
  
  _noOfDetectors = runHeader->getNoOfDetector(); 

  // the four vectors containing the first and the last pixel
  // along both the directions
  _minX = runHeader->getMinX();
  _maxX = runHeader->getMaxX();
  _minY = runHeader->getMinY();
  _maxY = runHeader->getMaxY();

  delete runHeader;

  // increment the run counter
  ++_iRun;
}


void EUTelCorrelator::processEvent (LCEvent * event) {

#ifdef MARLIN_USE_AIDA

  if (_iEvt % 10 == 0) 
    streamlog_out( MESSAGE4 ) << "Processing event " 
			      << setw(6) << setiosflags(ios::right) << event->getEventNumber() << " in run "
			      << setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber() 
			      << setfill(' ') << " (Total = " << setw(10) << _iEvt << ")" 
			      << resetiosflags(ios::left) << endl;
  ++_iEvt;

 
  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event) ;
  
  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return;
  } else if ( evt->getEventType() == kUNKNOWN ) {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
			       << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }

  if ( isFirstEvent() ) {
    
    bookHistos();

    _isFirstEvent = false;

  }

  
  try {

    LCCollectionVec * inputCollection   = static_cast<LCCollectionVec*> 
      (event->getCollection( _inputCollectionName ));

    CellIDDecoder<TrackerPulseImpl>  pulseCellDecoder( inputCollection );
    
    for ( size_t iExt = 0 ; iExt < inputCollection->size() ; ++iExt ) {

      TrackerPulseImpl * externalPulse = static_cast< TrackerPulseImpl * > 
	( inputCollection->getElementAt( iExt ) );

      EUTelVirtualCluster  * externalCluster;

      ClusterType type = static_cast<ClusterType>(static_cast<int>((pulseCellDecoder(externalPulse)["type"])));

      if ( type == kEUTelFFClusterImpl ) {
	externalCluster = new EUTelFFClusterImpl( static_cast<TrackerDataImpl*>
						  ( externalPulse->getTrackerData()) );
      } else {
	streamlog_out ( ERROR4 ) <<  "Unknown cluster type. Sorry for quitting" << endl;
	throw UnknownDataTypeException("Cluster type unknown");
      }

      int externalSensorID = pulseCellDecoder( externalPulse ) [ "sensorID" ] ;

      float externalXCenter;
      float externalYCenter;
      externalCluster->getCenterOfGravity( externalXCenter, externalYCenter ) ; 

      for ( size_t iInt = 0;  iInt <  inputCollection->size() ; ++iInt ) {

	TrackerPulseImpl * internalPulse = static_cast< TrackerPulseImpl * > 
	  ( inputCollection->getElementAt( iInt ) );
	EUTelVirtualCluster  * internalCluster;
	ClusterType type = static_cast<ClusterType>
	  (static_cast<int>((pulseCellDecoder(internalPulse)["type"])));
	if ( type == kEUTelFFClusterImpl ) {
	  internalCluster = new EUTelFFClusterImpl( static_cast<TrackerDataImpl*> 
						    (internalPulse->getTrackerData()) );
	} else {
	  streamlog_out ( ERROR4 ) <<  "Unknown cluster type. Sorry for quitting" << endl;
	  throw UnknownDataTypeException("Cluster type unknown");
	}
	int internalSensorID = pulseCellDecoder( internalPulse ) [ "sensorID" ] ;

	if ( internalSensorID != externalSensorID ) {

	  float internalXCenter;
	  float internalYCenter;
	  internalCluster->getCenterOfGravity( internalXCenter, internalYCenter ) ; 

	  streamlog_out ( DEBUG ) << "Filling histo " << externalSensorID << " " << internalSensorID << endl;

	  _clusterXCorrelationMatrix[ externalSensorID ][ internalSensorID ]->fill( internalXCenter, externalXCenter );
	  
	  
	}

      }
    }


  } catch (DataNotAvailableException& e  ) {
   
    streamlog_out  ( WARNING2 ) <<  "No input collection found on event " << event->getEventNumber() 
				<< " in run " << event->getRunNumber() << endl;
  }

#endif 

}

void EUTelCorrelator::end() {
  
    streamlog_out ( MESSAGE4 )  << "Successfully finished" << endl;
}

void EUTelCorrelator::bookHistos() {
  
#ifdef MARLIN_USE_AIDA

  try {

    streamlog_out ( MESSAGE4 ) <<  "Booking histograms" << endl;
    
    string tempHistoName;
    

    for ( int row = 0 ; row < _noOfDetectors; ++row ) {

      vector< AIDA::IHistogram2D * > innerVector;
      
      for ( int col = 0 ; col < _noOfDetectors; ++col ) {

	
	stringstream ss;
	ss << _clusterXCorrelationHistoName << "-d" << row 
	   << "-d" << col ;

	tempHistoName = ss.str();

	streamlog_out( MESSAGE ) << "Booking histo " << tempHistoName << endl;
	
	int     xBin = _maxX[ col ] - _minX[ col ] + 1;
	double  xMin = static_cast<double >(_minX[ col ]) - 0.5;
	double  xMax = static_cast<double >(_maxX[ col ]) + 0.5;
	int     yBin = _maxX[ row ] - _minX[ row ] + 1;
	double  yMin = static_cast<double >(_minX[ row ]) - 0.5;
	double  yMax = static_cast<double >(_maxX[ row ]) + 0.5;
	
	AIDA::IHistogram2D * histo2D = 
	  AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(), 
								    xBin, xMin, xMax, yBin, yMin, yMax );
	innerVector.push_back ( histo2D );
	
      }

      _clusterXCorrelationMatrix.push_back( innerVector ) ;

    }

  } catch (lcio::Exception& e ) {
    
    streamlog_out ( ERROR1 ) << "No AIDAProcessor initialized. Sorry for quitting..." << endl;
    exit( -1 ); 

  }
#endif
}

  

