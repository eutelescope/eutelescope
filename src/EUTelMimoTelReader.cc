// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelMimoTelReader.cc,v 1.6 2007-07-04 13:01:20 bulgheroni Exp $
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
#include <algorithm>

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
			     _polarity, static_cast<int> (-1));

  registerProcessorParameter("GeoID", "The geometry identification number",
			     _geoID, static_cast<int> ( 0 ));

  registerProcessorParameter("RemoveMarker","If set to true, pixels identified as markers will be removed\n"
			     "from the TrackerRawData output collections",
			     _removeMarkerSwitch, static_cast< bool > ( false ) );

  IntVec markerPositionExample;
  markerPositionExample.push_back( 0  );
  markerPositionExample.push_back( 1  );
  markerPositionExample.push_back( 66 );
  markerPositionExample.push_back( 67 );
  markerPositionExample.push_back( 132 );
  markerPositionExample.push_back( 133 );
  markerPositionExample.push_back( 198 );
  markerPositionExample.push_back( 199 );

  registerOptionalParameter("MarkerPosition","This vector of integer contains the marker positions in pixel number.\n"
			    "(Pixels are numbered starting from 0)",
			    _markerPositionVec, markerPositionExample );
			    


}

EUTelMimoTelReader * EUTelMimoTelReader::newProcessor () {
  return new EUTelMimoTelReader;
}

void EUTelMimoTelReader::init () {
  printParameters ();

  // those two values are constants!!!
  _xMax = 264;
  _yMax = 256;

  if ( _removeMarkerSwitch ) {

    message<DEBUG> ( "Data conversion with markers removal" );
    if ( _markerPositionVec.size() == 0 ) {
      message<WARNING> ( "Data conversion with markers removal selected but no markers position found\n"
			 "Disabling marker removal" );
      _removeMarkerSwitch = false;
      sort( _markerPositionVec.begin(), _markerPositionVec.end() );
    }
  } else {
    _markerPositionVec.clear();
  }
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
      runHeader->setMaxX(IntVec(noOfDetectors, _xMax - _markerPositionVec.size() -  1));
      runHeader->setMinY(IntVec(noOfDetectors, 0));
      runHeader->setMaxY(IntVec(noOfDetectors, _yMax - 1));
      runHeader->setDetectorName("MimoTel");
      runHeader->setGeoID( _geoID );

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
      idEncoder1["xMax"]     = _xMax - 1 - _markerPositionVec.size();
      idEncoder1["yMin"]     = 0;
      idEncoder1["yMax"]     = _yMax - 1;
      idEncoder2["xMin"]     = 0;
      idEncoder2["xMax"]     = _xMax - 1 - _markerPositionVec.size();
      idEncoder2["yMin"]     = 0;
      idEncoder2["yMax"]     = _yMax - 1;
      idEncoder3["xMin"]     = 0;
      idEncoder3["xMax"]     = _xMax - 1 - _markerPositionVec.size();
      idEncoder3["yMin"]     = 0;
      idEncoder3["yMax"]     = _yMax - 1;
      idEncoderCDS["xMin"]   = 0;
      idEncoderCDS["xMax"]   = _xMax - 1 - _markerPositionVec.size();
      idEncoderCDS["yMin"]   = 0;
      idEncoderCDS["yMax"]   = _yMax - 1;
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
	  
	    if ( _removeMarkerSwitch ) {
	      // the idea behind the marker removal procedure is that:
	      // markers are occurring always at the same column
	      // position, so we can repeat the procedure into a loop
	      // over all the rows. The marker positions are given by
	      // the user in the steering file and stored into the
	      // _markerPositionVec. For each row the part of the
	      // original vector in between two markers is copied into
	      // the stripped matrix.

	      vector<short > firstStrippedVec(  ( _xMax - _markerPositionVec.size() ) * _yMax );
	      vector<short > secondStrippedVec( ( _xMax - _markerPositionVec.size() ) * _yMax );
	      vector<short > thirdStrippedVec(  ( _xMax - _markerPositionVec.size() ) * _yMax );
	      vector<short >::iterator currentFirstPos  = firstStrippedVec.begin();
	      vector<short >::iterator currentSecondPos = secondStrippedVec.begin();
	      vector<short >::iterator currentThirdPos  = thirdStrippedVec.begin();
	      vector<short >::iterator firstBegin       = array.m_adc[0].begin();
	      vector<short >::iterator secondBegin      = array.m_adc[1].begin();
	      vector<short >::iterator thirdBegin       = array.m_adc[2].begin();

	      for ( int y = 0; y < _yMax ; y++ ) {
		int offset = y * _xMax;
		vector<int >::iterator marker = _markerPositionVec.begin();
		// first of all copy the part of the array going from
		// the row beginning to the first marker position
		currentFirstPos    = copy( firstBegin + offset, 
					   firstBegin + ( *(marker) + offset ),
					   currentFirstPos );
		currentSecondPos   = copy( secondBegin + offset,
					   secondBegin + ( *(marker) + offset ),
					   currentSecondPos );
		currentThirdPos    = copy( thirdBegin + offset, 
					   thirdBegin + ( *(marker) + offset ),
					   currentThirdPos );

		while ( marker != _markerPositionVec.end() ) {
		  if ( marker < _markerPositionVec.end() - 1 ) {
		    currentFirstPos   = copy( firstBegin + ( *(marker) + 1 + offset),
					      firstBegin + ( *(marker + 1 ) + offset),
					      currentFirstPos );
		    currentSecondPos  = copy( secondBegin + ( *(marker) + 1 + offset),
					      secondBegin + ( *(marker + 1 ) + offset),
					      currentSecondPos );
		    currentThirdPos   = copy( thirdBegin + ( *(marker) + 1 + offset),
					      thirdBegin + ( *(marker + 1 ) + offset),
					      currentThirdPos );
		  } else {
		    // this is the case where we have to copy the part
		    // of the original vector going from the last
		    // marker and the end of the row.
		    currentFirstPos = copy( firstBegin + ( *(marker) + 1 + offset ),
					    firstBegin + offset + _xMax, 
					    currentFirstPos );

		    currentSecondPos = copy( secondBegin + ( *(marker) + 1 + offset ),
					    secondBegin + offset + _xMax, 
					    currentSecondPos );

		    currentThirdPos = copy( thirdBegin + ( *(marker) + 1 + offset ),
					    thirdBegin + offset + _xMax, 
					    currentThirdPos );
		  }
		  ++marker;
		}
	      }

	      firstFrame->setADCValues( firstStrippedVec );
	      secondFrame->setADCValues( secondStrippedVec );
	      thirdFrame->setADCValues( thirdStrippedVec );
	      
	      
	    } else {
	      
	      // just move the adc values from the EUDRBDecoder::arrays_t to
	      // the trackerRawDataImpl


	      firstFrame->setADCValues(  array.m_adc[0] );
	      secondFrame->setADCValues( array.m_adc[1] );
	      thirdFrame->setADCValues(  array.m_adc[2] );

	    }

	    firstFrame->setTime(brd.PivotPixel());
	    secondFrame->setTime(brd.PivotPixel());	      
	    thirdFrame->setTime(brd.PivotPixel());
	    
	    if ( _cdsCalculation ) {
	      TrackerRawDataImpl * cdsFrame = new TrackerRawDataImpl;
	      idEncoderCDS["sensorID"] = iDetector;
	      idEncoderCDS.setCellID( cdsFrame );
	      vector<short > cdsVector;
	      vector<short > firstFrameVec  = array.m_adc[0];
	      vector<short > secondFrameVec = array.m_adc[1];
	      vector<short > thirdFrameVec  = array.m_adc[2]; 

	      // bad luck! the pivot pixel is obtained counting also
	      // the markers, so I prefer to calculate the cds as the
	      // marker should stay into and in case remove them afterward.
	      for (unsigned int iPixel = 0; iPixel < eudrbDecoder->NumPixels(brd); iPixel++) {
		cdsVector.push_back( static_cast<short> (_polarity )* 
				     ( ( -1 * array.m_pivot[iPixel]    ) * firstFrameVec[iPixel] +
				       ( 2 * array.m_pivot[iPixel] - 1 ) * secondFrameVec[iPixel] +
				       ( 1 - array.m_pivot[iPixel]     ) * thirdFrameVec[iPixel] ) );
	      }	      


	      if ( _removeMarkerSwitch ) {
		// strip away the markers...
		vector<short > cdsStrippedVec( ( _xMax - _markerPositionVec.size() ) * _yMax );
		vector<short >::iterator currentCDSPos = cdsStrippedVec.begin();
		vector<short >::iterator cdsBegin      = cdsVector.begin();
		for ( int y = 0; y < _yMax ; y ++ ) {
		  int offset = y * _xMax ;
		  vector<int >::iterator marker = _markerPositionVec.begin();
		  currentCDSPos = copy( cdsBegin + offset, 
					cdsBegin + ( *(marker) + offset ),
					currentCDSPos);
		  while ( marker != _markerPositionVec.end() ) {
		    if ( marker < _markerPositionVec.end() - 1 ) {
		      currentCDSPos = copy( cdsBegin + ( *(marker) + 1 + offset),
					    cdsBegin + ( *(marker + 1 ) + offset),
					    currentCDSPos );
		    } else {
		      currentCDSPos = copy( cdsBegin + ( *(marker) + 1 + offset ),
					    cdsBegin + offset + _xMax, 
					    currentCDSPos );
		    }
		    ++marker;
		  }
		}
		cdsFrame->setADCValues( cdsStrippedVec );
	      } else {
		cdsFrame->setADCValues( cdsVector );
	      }
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

