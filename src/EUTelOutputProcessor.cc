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

// eutelescope includes ".h"
#include "EUTelOutputProcessor.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTELESCOPE.h"

// marlin includes ".h"
#include "marlin/LCIOOutputProcessor.h"

// lcio includes <.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/LCTime.h>

// system includes <>
#include <memory>

using namespace std;
using namespace marlin;
using namespace eutelescope;

 
EUTelOutputProcessor::EUTelOutputProcessor() : LCIOOutputProcessor("EUTelOutputProcessor")  {
    
  _description = "Writes the current event to the specified LCIO outputfile."
    " Eventually it adds a EORE at the of the file if it was missing"
    " Needs to be the last ActiveProcessor." ;
    
#ifndef MARLIN_VERSION_GE

  registerProcessorParameter( "LCIOOutputFile" , 
			      " name of output file "  ,
			      _lcioOutputFile ,
			      std::string("outputfile.slcio") ) ;
    
  registerProcessorParameter( "LCIOWriteMode" , 
			      "write mode for output file:  WRITE_APPEND or WRITE_NEW"  ,
			      _lcioWriteMode ,
			      std::string("None") ) ;


  StringVec dropNamesExamples ;
  dropNamesExamples.push_back("rawdata");
  dropNamesExamples.push_back("data");
  dropNamesExamples.push_back("pedestal");
  dropNamesExamples.push_back("noise");
  dropNamesExamples.push_back("status");

  registerOptionalParameter( "DropCollectionNames" , 
			     "drops the named collections from the event"  ,
			     _dropCollectionNames ,
			     dropNamesExamples ) ;
    
    
  StringVec dropTypesExample ;
  dropTypesExample.push_back("TrackerRawData");
  dropTypesExample.push_back("TrackerData");
    
  registerOptionalParameter( "DropCollectionTypes" , 
			     "drops all collections of the given type from the event"  ,
			     _dropCollectionTypes ,
			     dropTypesExample ) ;
    

  registerOptionalParameter( "SplitFileSizekB" , 
			     "will split output file if size in kB exceeds given value - doesn't work with APPEND and NEW"  ,
			     _splitFileSizekB, 
			     1992294 ) ;  // 1.9 GB in kB
#endif    

  registerProcessorParameter("SkipIntermediateEORE",
			     "Set it to true to remove intermediate EORE in merged runs",
			     _skipIntermediateEORESwitch, static_cast< bool > ( true ) );


}

void EUTelOutputProcessor::init() { 
  
  // needs to be reimplemented since it is virtual in
  // LCIOOutputProcessor
  LCIOOutputProcessor::init();

}



void EUTelOutputProcessor::processRunHeader( LCRunHeader* run) { 
  std::unique_ptr<EUTelRunHeaderImpl> runHeader = std::make_unique<EUTelRunHeaderImpl>(run);
  runHeader->addProcessor(type());
  LCIOOutputProcessor::processRunHeader(run);
} 

void EUTelOutputProcessor::processEvent( LCEvent * evt ) { 

  EUTelEventImpl * eutelEvt =  static_cast<EUTelEventImpl * > ( evt );

  if ( _skipIntermediateEORESwitch && ( eutelEvt->getEventType() == kEORE ) ) {
    // ok the user wants me to skip this because it may be an
    // intermediate EORE. But how to know this is not the last of the
    // intermediate EORE? Easy, don't care, set the _eventType to kDE
    // and in case this the last one, end() will be called soon after
    // and a kEORE will be appended!
    _eventType = kDE;
    return ;
  }

  LCIOOutputProcessor::processEvent(evt);
  _eventType = eutelEvt->getEventType();

}

void EUTelOutputProcessor::end(){ 

  if ( _eventType != kEORE ) {

    message<WARNING> ( "Adding a EORE because was missing" );

    EUTelEventImpl * event = new EUTelEventImpl;
    event->setDetectorName("kEORE fix by EUTelOutputProcessor");
    event->setEventType(kEORE);
    event->setEventNumber( _nEvt + 1 );
    
    LCTime * now = new LCTime;
    event->setTimeStamp(now->timeStamp());
    delete now;

    _lcWrt->writeEvent( static_cast<LCEventImpl*> (event) );
    delete event;

  }

  message<MESSAGE5> ( log() << "Writing the output file " << _lcioOutputFile );
  _lcWrt->close() ;

}


