// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelMimoTelReader.cc,v 1.4 2007-06-29 15:24:23 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifdef USE_EUDAQ
// personal includes
#include "EUTelMimoTelReader.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelExceptions.h"

// marlin includes
#include "marlin/Exceptions.h"
#include "marlin/DataSourceProcessor.h"
#include "marlin/ProcessorMgr.h"

// eudaq includes 
#include <eudaq/FileSerializer.hh>
#include <eudaq/EUDRBEvent.hh>
#include <eudaq/DetectorEvent.hh>
#include <eudaq/Utils.hh>

// lcio includes
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/LCEventImpl.h>
#include <UTIL/CellIDEncoder.h>

// system includes 
#include <string>
#include <vector>


using namespace std;
using namespace marlin;
using namespace eutelescope;
using namespace eudaq;

EUTelMimoTelReader::EUTelMimoTelReader (): DataSourceProcessor  ("EUTelMimoTelReader") {
  
  _description =
    "Reads EUDRB events with MimoTel data inside";
 
  registerOptionalParameter("CDS","Enable (==1) or disable (==0) the CDS calculation",
			  _cdsCalculation, static_cast< bool > ( true ) );

  registerOutputCollection(LCIO::TRACKERRAWDATA,"FirstFrameCollectionName",
			   "First frame collection name", _firstFrameCollectionName, 
			   string( "firstFrame" ));
  registerOutputCollection(LCIO::TRACKERRAWDATA,"SecondFrameCollectionName",
			   "Second frame collection name", _secondFrameCollectionName, 
			   string( "secondFrame" ));
  registerOutputCollection(LCIO::TRACKERRAWDATA,"ThirdFrameCollectionName",
			   "Third frame collection name", _thirdFrameCollectionName, 
			   string( "thirdFrame" ));
  registerOutputCollection(LCIO::TRACKERRAWDATA,"CDSCollection",
			   "CDS collection name", _cdsCollectionName, string( "rawdata" ));
  
  registerProcessorParameter("InputDataFileName", "Inpuf file",
			     _fileName, string("run012345.raw") );
  
  registerProcessorParameter("SignalPolarity", "Signal polarity (negative == -1)",
			     _polarity, float(-1));


}

EUTelMimoTelReader * EUTelMimoTelReader::newProcessor () {
  return new EUTelMimoTelReader;
}

void EUTelMimoTelReader::init () {
  printParameters ();
}


void EUTelMimoTelReader::readDataSource (int numEvents) {

  int eventNumber = 0;
  
  message<DEBUG> ( log() << "Reading " << _fileName << " with eudaq file deserializer " );

  // this is eudaq de-serialiazer for the input file 
  FileDeserializer des( _fileName );
  EUDRBDecoder * eudrbDecoder = NULL;

  while ( des.HasData() ) {
    counted_ptr<Event> ev(EventFactory::Create(des));
    DetectorEvent * dev = dynamic_cast<DetectorEvent*>(ev.get());

    if ( ev->IsBORE() ) {
      // that's the right place to process the run header
      message<DEBUG> ( log() << "Found a BORE, so processing the RunHeader");
      EUTelRunHeaderImpl * runHeader = new EUTelRunHeaderImpl;
      runHeader->setRunNumber( ev->GetRunNumber() );
      
      int noOfDetectors = 0;
      
      for (unsigned int i = 0; i < dev->NumEvents(); ++i) {
	Event * subev = dev->GetEvent(i);
	EUDRBEvent * eudev = dynamic_cast<EUDRBEvent*>(subev);
	if (eudev) {
	  noOfDetectors += from_string(eudev->GetTag("BOARDS"), 0) ;
	  string mode = eudev->GetTag("MODE");
	  string det  = eudev->GetTag("DET");
	  if ( !eudrbDecoder ) {
	    message<WARNING> ("Assuming all the EUDRB producers running with the same mode and detector");
	    eudrbDecoder = new EUDRBDecoder(*eudev);
	  } 
	  if (  mode != "RAW3" ) {
	    throw InvalidParameterException("For the time being only RAW3 is supported");
	  }
	}
      }
      runHeader->setNoOfDetector(noOfDetectors);
      runHeader->setMinX(IntVec(noOfDetectors, 0));
      runHeader->setMaxX(IntVec(noOfDetectors, 263));
      runHeader->setMinY(IntVec(noOfDetectors, 0));
      runHeader->setMaxY(IntVec(noOfDetectors, 255));
      runHeader->setDetectorName("MimoTel");
      
      runHeader->setDateTime();
      
      ProcessorMgr::instance()->processRunHeader(runHeader);
      delete runHeader;
    } else if ( ev->IsEORE() ) {
      message<DEBUG>( log() << "Found a EORE, processing a dummy empty event");
      EUTelEventImpl * event = new EUTelEventImpl;
      event->setEventType(kEORE);
      event->setTimeStamp( ev->GetTimestamp() );
      event->setRunNumber( ev->GetRunNumber() );
      event->setEventNumber( ev->GetEventNumber() );
      ProcessorMgr::instance()->processEvent( static_cast<LCEventImpl*> (event) ) ;
      delete event;
    } else if (!dev) {
      message<WARNING> ( log() << "Event number " << ev->GetEventNumber() << " does not contain any data ");
      throw SkipEventException(this);
    } else {
      // this is the last case, this is a kDE event no chance to fail!
      EUTelEventImpl * event = new EUTelEventImpl;
      event->setRunNumber( ev->GetRunNumber() );
      event->setEventNumber( ev->GetEventNumber() );
      event->setTimeStamp( ev->GetTimestamp() );
      event->setEventType( kDE );

      // prepare a collection for each frame and set up the cell ID
      // encoder for each of them. 
      LCCollectionVec * firstFrameColl  = new LCCollectionVec (LCIO::TRACKERRAWDATA);
      LCCollectionVec * secondFrameColl = new LCCollectionVec (LCIO::TRACKERRAWDATA);
      LCCollectionVec * thirdFrameColl  = new LCCollectionVec (LCIO::TRACKERRAWDATA);
      LCCollectionVec * cdsFrameColl    = new LCCollectionVec (LCIO::TRACKERRAWDATA);
      CellIDEncoder< TrackerRawDataImpl > idEncoder1   (EUTELESCOPE::MATRIXDEFAULTENCODING, firstFrameColl);
      CellIDEncoder< TrackerRawDataImpl > idEncoder2   (EUTELESCOPE::MATRIXDEFAULTENCODING, secondFrameColl);
      CellIDEncoder< TrackerRawDataImpl > idEncoder3   (EUTELESCOPE::MATRIXDEFAULTENCODING, thirdFrameColl);
      CellIDEncoder< TrackerRawDataImpl > idEncoderCDS (EUTELESCOPE::MATRIXDEFAULTENCODING, cdsFrameColl);
    
      // set here static parameters, while keeping the sensor id
      // setting for within the sensor loop
      idEncoder1["xMin"]     = 0;
      idEncoder1["xMax"]     = 263;
      idEncoder1["yMin"]     = 0;
      idEncoder1["yMax"]     = 255;
      idEncoder2["xMin"]     = 0;
      idEncoder2["xMax"]     = 263;
      idEncoder2["yMin"]     = 0;
      idEncoder2["yMax"]     = 255;
      idEncoder3["xMin"]     = 0;
      idEncoder3["xMax"]     = 263;
      idEncoder3["yMin"]     = 0;
      idEncoder3["yMax"]     = 255;
      idEncoderCDS["xMin"]   = 0;
      idEncoderCDS["xMax"]   = 263;
      idEncoderCDS["yMin"]   = 0;
      idEncoderCDS["yMax"]   = 255;
      if ( eventNumber % 10 == 0 ) 
	message<MESSAGE> ( log() << "Converting event " << eventNumber );

      // The detector event contains one sub-event for each
      // producer. It is very likely to have a sub event for the TLU
      // and another one for the EUDRB

      for (unsigned int iProducer = 0; iProducer < dev->NumEvents(); iProducer++ ) {
	// get from the producer the event into a virtual base class
	// that is to say "Event"
	Event * subev = dev->GetEvent(iProducer);
	message<DEBUG> ( log() << "Processing producer number  " << iProducer );

	// now try to see if the subev we got is a EUDRB event and
	// continue only in this case
	EUDRBEvent * eudev = dynamic_cast<EUDRBEvent*> (subev) ;
	
	if ( eudev ) {
	  // ok great this is a EUDRB event, now we need to loop on
	  // the number of boards the EUDRB producer is reading out
	  for (unsigned  int iDetector = 0; iDetector < eudev->NumBoards(); iDetector++) {
	    // EUDRBBoard is the wrapper class containing the real
	    // data we are interested in
	    EUDRBBoard & brd = eudev->GetBoard(iDetector);
	    EUDRBDecoder::arrays_t<short, short> array = eudrbDecoder->GetArrays<short,short>(brd);
	    

	    // prepare a TrackerRawDataImpl for each of the frame and
	    // the corresponding cellIDEncoder.
	    TrackerRawDataImpl * firstFrame = new TrackerRawDataImpl;
	    idEncoder1["sensorID"] = iDetector;
	    idEncoder1.setCellID( firstFrame );

	    TrackerRawDataImpl * secondFrame = new TrackerRawDataImpl;
	    idEncoder2["sensorID"] = iDetector;
	    idEncoder2.setCellID( secondFrame );

	    TrackerRawDataImpl * thirdFrame = new TrackerRawDataImpl;
	    idEncoder3["sensorID"] = iDetector;
	    idEncoder3.setCellID( thirdFrame );
	  
	    // move the adc values from the EUDRBDecoder::arrays_t to
	    // the trackerRawDataImpl
	    firstFrame->setADCValues(  array.m_adc[0] );
	    firstFrame->setTime(brd.PivotPixel());
	    secondFrame->setADCValues( array.m_adc[1] );
	    secondFrame->setTime(brd.PivotPixel());
	    thirdFrame->setADCValues(  array.m_adc[2] );
	    thirdFrame->setTime(brd.PivotPixel());

	    if ( _cdsCalculation ) {
	      TrackerRawDataImpl * cdsFrame = new TrackerRawDataImpl;
	      idEncoderCDS["sensorID"] = iDetector;
	      idEncoderCDS.setCellID( cdsFrame );
	      vector< short > cdsVector;
	      for (unsigned int iPixel = 0; iPixel < eudrbDecoder->NumPixels(brd); iPixel++) {
		cdsVector.push_back( _polarity * 
				     ( ( array.m_pivot[iPixel] - 1   ) * array.m_adc[0][iPixel] +
				       ( 1 - 2*array.m_pivot[iPixel] ) * array.m_adc[1][iPixel] +
				       ( 1*array.m_pivot[iPixel]     ) * array.m_adc[2][iPixel] ) );
	      }	      
	      cdsFrame->setADCValues( cdsVector );
	      cdsFrameColl->push_back( cdsFrame );
	    }	    


	    // add the trackerRawData to the corresponding collection
	    firstFrameColl->push_back(  firstFrame   );
	    secondFrameColl->push_back( secondFrame  );
	    thirdFrameColl->push_back(  thirdFrame   );
	  }
	} else {
	  message<DEBUG> ( log() << "Not a EUDRBEvent, very likely a TLUEvent " );
	}
      }
      event->addCollection(firstFrameColl,  _firstFrameCollectionName);
      event->addCollection(secondFrameColl, _secondFrameCollectionName);
      event->addCollection(thirdFrameColl,  _thirdFrameCollectionName);
      if ( _cdsCalculation ) event->addCollection( cdsFrameColl, _cdsCollectionName );
	
      ProcessorMgr::instance()->processEvent( static_cast<LCEventImpl*> (event) );
      delete event;
    }
    ++eventNumber;
    
    if ( eventNumber >=  numEvents ) {
      // even if there are some more events in the input file, we
      // don't want to decode them. So add (eventually another kEORE)
      // and then stop the data reading
      EUTelEventImpl * event = new EUTelEventImpl;
      event->setEventType(kEORE);
      event->setEventNumber( eventNumber );
      ProcessorMgr::instance()->processEvent( static_cast<LCEventImpl*> (event) ) ;
      delete event;
    }
  }
}    




void EUTelMimoTelReader::end () {
  message<MESSAGE> ("Successfully finished") ;
}

#endif

