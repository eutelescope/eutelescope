// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelApplyAlignmentProcessor.cc,v 1.4 2008-07-28 16:13:03 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h" 
#include "EUTelApplyAlignmentProcessor.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelAlignmentConstant.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"

// lcio includes <.h> 
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

// system includes <>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;



EUTelApplyAlignmentProcessor::EUTelApplyAlignmentProcessor () :Processor("EUTelApplyAlignmentProcessor") {

  // modify processor description
  _description =
    "Apply to the input hit the alignment corrections";

  // first of all we need to register the input collection
  registerInputCollection (LCIO::TRACKERHIT, "InputHitCollectionName",
			   "The name of the input hit collection",
			   _inputHitCollectionName, string ("hit"));

  registerInputCollection (LCIO::LCGENERICOBJECT, "AlignmentConstantName",
			   "Alignment constant from the condition file",
			   _alignmentCollectionName, string ("alignment"));

  registerOutputCollection (LCIO::TRACKERHIT, "OutputHitCollectionName",
			   "The name of the output hit collection",
			   _outputHitCollectionName, string("correctedHit"));


  // now the optional parameters
  registerProcessorParameter ("CorrectionMethod",
			      "Available methods are:\n"
			      " 0 --> shift only \n"
			      " 1 --> rotation first \n"
			      " 2 --> shift first ",
			      _correctionMethod, static_cast<int > (1));
  

}


void EUTelApplyAlignmentProcessor::init () {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

}

void EUTelApplyAlignmentProcessor::processRunHeader (LCRunHeader * rdr) {

  // convert the run header into something mode easy to digest
  auto_ptr<EUTelRunHeaderImpl> runHeader( new EUTelRunHeaderImpl( rdr ) );

  string eudrbGlobalMode = runHeader->getEUDRBMode() ;

  transform( eudrbGlobalMode.begin(), eudrbGlobalMode.end(), eudrbGlobalMode.begin(), ::tolower);

  _hasNZSData = false;

  if ( ( eudrbGlobalMode == "raw2" ) ||
       ( eudrbGlobalMode == "raw3" ) ||
       ( eudrbGlobalMode == "mixed" ) ) {
    _hasNZSData = true;
  }

  _hasZSData = false;

  if ( ( eudrbGlobalMode == "zs" ) || 
       ( eudrbGlobalMode == "mixed" ) ) {
    _hasZSData = true;
  }

  // increment the run counter
  ++_iRun;  

}


void EUTelApplyAlignmentProcessor::processEvent (LCEvent * event) {

  if ( _iEvt % 10 == 0 ) 
    streamlog_out ( MESSAGE4 ) << "Processing event " 
			       << setw(6) << setiosflags(ios::right) << event->getEventNumber() << " in run "
			       << setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber() << setfill(' ')
			       << " (Total = " << setw(10) << _iEvt << ")" << resetiosflags(ios::left) << endl;
  ++_iEvt;


  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);

  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return;
  } else if ( evt->getEventType() == kUNKNOWN ) {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
			       << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }



  // the cell decoder to get the sensor ID
  CellIDDecoder<TrackerDataImpl> * clusterCellDecoder = NULL; //( originalZSDataCollectionVec );
  CellIDDecoder<TrackerDataImpl> * clusterZSCellDecoder = NULL;
  
    
  try {
    
    LCCollectionVec * inputCollectionVec         = dynamic_cast < LCCollectionVec * > (evt->getCollection(_inputHitCollectionName));
    LCCollectionVec * alignmentCollectionVec     = dynamic_cast < LCCollectionVec * > (evt->getCollection(_alignmentCollectionName));

    // I also need the original data collection for ZS and NZS data
    LCCollectionVec * originalDataCollectionVec = NULL;
    LCCollectionVec * originalZSDataCollectionVec = NULL;

    
    if ( _hasNZSData ) {
      originalDataCollectionVec = dynamic_cast < LCCollectionVec * > (evt->getCollection( "original_data" ) );
      clusterCellDecoder = new CellIDDecoder<TrackerDataImpl>(  originalDataCollectionVec );
    }
    if ( _hasZSData ) {
      originalZSDataCollectionVec = dynamic_cast < LCCollectionVec * > (evt->getCollection( "original_zsdata" ) );
      clusterZSCellDecoder = new CellIDDecoder<TrackerDataImpl>(  originalZSDataCollectionVec );
    }


    

    if (isFirstEvent()) {
      
      
      
      for ( size_t iPos = 0; iPos < alignmentCollectionVec->size(); ++iPos ) {
	
	EUTelAlignmentConstant * alignment = static_cast< EUTelAlignmentConstant * > ( alignmentCollectionVec->getElementAt( iPos ) );
	_lookUpTable[ alignment->getSensorID() ] = iPos;

      }

      _isFirstEvent = false;
    }
  

    LCCollectionVec * outputCollectionVec = new LCCollectionVec(LCIO::TRACKERHIT);
    
    for (size_t iHit = 0; iHit < inputCollectionVec->size(); iHit++) {
      
      TrackerHitImpl   * inputHit   = dynamic_cast< TrackerHitImpl * >  ( inputCollectionVec->getElementAt( iHit ) ) ;
      
      // now we have to understand which layer this hit belongs to. 
      LCObjectVec        clusterVec = inputHit->getRawHits();
      TrackerDataImpl  * cluster    = dynamic_cast< TrackerDataImpl *>  ( clusterVec[0] );
    
      int sensorID;
      
      if ( _hasZSData && _hasNZSData ) {
	// it means that this is a MIXED run still I don't know what
	// to do
	streamlog_out ( ERROR ) << "This processor is unable to deal with MIXED data. Sorry for quitting..." << endl;
	exit(-01);
      }
      if ( _hasNZSData ) sensorID = (*clusterCellDecoder)( cluster ) ["sensorID"] ;
      if ( _hasZSData  ) sensorID = (*clusterZSCellDecoder)( cluster ) ["sensorID"]   ;
      
      // now that we know at which sensor the hit belongs to, we can
      // get the corresponding alignment constants
      EUTelAlignmentConstant * alignment = static_cast< EUTelAlignmentConstant * > 
	( alignmentCollectionVec->getElementAt( _lookUpTable[ sensorID ] ) );


      TrackerHitImpl   * outputHit  = new TrackerHitImpl;
      outputHit->setType( inputHit->getType() );
      outputHit->rawHits() = clusterVec;
      

      double * inputPosition      = const_cast< double * > ( inputHit->getPosition() ) ;
      double   outputPosition[3]  = { 0., 0., 0. };

      if ( _correctionMethod == 0 ) {
	
	// this is the shift only case 
	
	outputPosition[0] = inputPosition[0] - alignment->getXOffset();
	outputPosition[1] = inputPosition[1] - alignment->getYOffset();
	outputPosition[2] = inputPosition[2] - alignment->getZOffset();

      } else if ( _correctionMethod == 1 ) {
	
	// this is the rotation first
	
	// first the rotation
	outputPosition[0] = inputPosition[0] + alignment->getGamma() * inputPosition[1] + alignment->getBeta() * inputPosition[2] ;
	outputPosition[1] = -1 * alignment->getGamma() * inputPosition[0] + inputPosition[1] + alignment->getAlpha() * inputPosition[2];
	outputPosition[2] = -1 * alignment->getBeta()  * inputPosition[0] + alignment->getAlpha() * inputPosition[1] + inputPosition[2];

	// second the shift
	outputPosition[0] -= alignment->getXOffset();
	outputPosition[1] -= alignment->getYOffset();
	outputPosition[2] -= alignment->getZOffset();

      } else if ( _correctionMethod == 2 ) {
	
	// this is the translation first
	
	// first the shifts
	inputPosition[0] -= alignment->getXOffset();
	inputPosition[1] -= alignment->getYOffset();
	inputPosition[2] -= alignment->getZOffset();
	
	// second the rotation
	outputPosition[0] = inputPosition[0] + alignment->getGamma() * inputPosition[1] + alignment->getBeta() * inputPosition[2] ;
	outputPosition[1] = -1 * alignment->getGamma() * inputPosition[0] + inputPosition[1] + alignment->getAlpha() * inputPosition[2];
	outputPosition[2] = -1 * alignment->getBeta()  * inputPosition[0] + alignment->getAlpha() * inputPosition[1] + inputPosition[2];

      }
      
      outputHit->setPosition( outputPosition ) ;
      outputCollectionVec->push_back( outputHit );  

    }
    

    evt->addCollection( outputCollectionVec, _outputHitCollectionName );
  
    delete clusterCellDecoder;
    delete clusterZSCellDecoder;
  
  } catch (DataNotAvailableException& e) {
    if ( clusterCellDecoder ) delete clusterCellDecoder;
    if ( clusterZSCellDecoder ) delete clusterZSCellDecoder;
    streamlog_out  ( WARNING2 ) <<  "No input collection found on event " << event->getEventNumber() 
				<< " in run " << event->getRunNumber() << endl;
  }
  
}
  


void EUTelApplyAlignmentProcessor::check (LCEvent * evt) {
  // nothing to check here - could be used to fill check plots in reconstruction processor
}


void EUTelApplyAlignmentProcessor::end() {
  streamlog_out ( MESSAGE2 ) <<  "Successfully finished" << endl;

}

