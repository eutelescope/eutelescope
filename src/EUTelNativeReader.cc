// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Author Loretta Negrini, Univ. Insubria <mailto:loryneg@gmail.com>
// Author Silvia Bonfanti, Univ. Insubria <mailto:silviafisica@gmail.com>
// Version $Id: EUTelNativeReader.cc,v 1.14 2008-11-12 12:01:37 furletova Exp $
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
#include "EUTelNativeReader.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelBaseDetector.h"
#include "EUTelPixelDetector.h"
#include "EUTelMimoTelDetector.h"
#include "EUTelMimosa18Detector.h"
#include "EUTelTLUDetector.h"
#include "EUTelDEPFETDetector.h"
#include "EUTelSparseDataImpl.h"
#include "EUTelSimpleSparsePixel.h"
#include "EUTelSetupDescription.h"

// marlin includes
#include "marlin/Exceptions.h"
#include "marlin/Processor.h"
#include "marlin/DataSourceProcessor.h"
#include "marlin/ProcessorMgr.h"

// eudaq includes
#include <eudaq/FileSerializer.hh>
#include <eudaq/EUDRBEvent.hh>
#include <eudaq/DetectorEvent.hh>
#include <eudaq/DEPFETEvent.hh>
#include <eudaq/Utils.hh>
#include <eudaq/Exception.hh>

// lcio includes
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCGenericObjectImpl.h>
#include <IMPL/LCEventImpl.h>
#include <UTIL/CellIDEncoder.h>

// system includes
#include <iostream>
#include <cassert>
#include <memory>

using namespace std;
using namespace marlin;
using namespace eutelescope;

// initialize static members
const unsigned short EUTelNativeReader::_eudrbOutOfSyncThr = 2;
const unsigned short EUTelNativeReader::_eudrbMaxConsecutiveOutOfSyncWarning = 20;

EUTelNativeReader::EUTelNativeReader (): DataSourceProcessor  ("EUTelNativeReader") {

  // initialize few variables
  _eudrbPreviousOutOfSyncEvent      = 0;
  _eudrbConsecutiveOutOfSyncWarning = 0;
  _eudrbTotalOutOfSyncEvent         = 0;

  _description =
    "Reads data streams produced by EUDAQ and produced the corresponding LCIO output";

  registerProcessorParameter("InputFileName", "This is the input file name",
                             _fileName, string("run012345.raw") );

  registerProcessorParameter("GeoID", "The geometry identification number", _geoID, static_cast<int> ( 0 ));


  // from here below only detector specific parameters.

  // ---------- //
  //  EUDRB     //
  // ---------- //

  // first the compulsory parameters

  registerOutputCollection(LCIO::TRACKERRAWDATA, "EUBRDRawModeOutputCollection",
                           "This is the eudrb producer output collection when read in RAW mode",
                           _eudrbRawModeOutputCollectionName, string("rawdata") );

  registerOutputCollection(LCIO::TRACKERDATA, "EUDRBZSModeOutputCollection",
                           "This si the mimotel output collection when read in ZS mode",
                           _eudrbZSModeOutputCollectionName, string("zsdata") );

  registerOptionalParameter("EUDRBSparsePixelType",
                            "Type of sparsified pixel data structure (use SparsePixelType enumerator)",
                            _eudrbSparsePixelType , static_cast<int> ( 1 ) );


  // to clean up the field justification
//  streamlog_out ( ERROR4 ) << setiosflags( ios::left );
//  streamlog_out ( ERROR4 ) << resetiosflags( ios::right );


}

EUTelNativeReader * EUTelNativeReader::newProcessor () {
  return new EUTelNativeReader;
}

void EUTelNativeReader::init () {
  printParameters ();


}


void EUTelNativeReader::readDataSource (int numEvents) {

  // this event counter is used to stop the processing when it is
  // greater than numEvents.
  int eventCounter = 0;

  // by definition this is the first event! 
  _isFirstEvent = true;

  // this is to make the output messages nicer
  streamlog::logscope scope(streamlog::out);
  scope.setName(name());

  streamlog_out( DEBUG4 ) << "Reading " << _fileName << " with eudaq file deserializer " << endl;

  // this is eudaq de-serialiazer for the input file
  eudaq::FileDeserializer nativeDeserializer( _fileName );

  // start the loop on all the event available into the input file
  while ( nativeDeserializer.HasData() ) {

    // get an event from the deserializer
    counted_ptr< eudaq::Event > eudaqEvent ( eudaq::EventFactory::Create( nativeDeserializer ) );

    // inform the user about the reading status
    if ( eventCounter % 10 == 0 )
      streamlog_out ( MESSAGE4 ) << "Processing event "
                                 << setw(6) << setiosflags( ios::right ) << eudaqEvent->GetEventNumber() << resetiosflags(ios::right)
                                 << " in run "  << setw(6)
                                 << setiosflags( ios::right ) << setfill('0') << eudaqEvent->GetRunNumber() << resetiosflags(ios::right)
                                 << setfill(' ') << " (Total = " << setw(10)
                                 << setiosflags( ios::right ) << eventCounter << resetiosflags(ios::right) << ")"
                                 << setiosflags( ios::left ) << endl;

    if ( eventCounter >= numEvents ) {
      // even if there are some more events in the input file, we
      // don't want to decode them. So add a kEORE
      // and then stop the data reading
      processEORE( eudaqEvent.get() );

      break; // this is breaking the loop on input events
    }

    if ( eudaqEvent->IsBORE() ) {
      // this is the case in which the eudaq event is a Begin Of Run
      // Event. This is translating into a RunHeader in the case of
      // LCIO
      processBORE( eudaqEvent.get() );

    } else if ( eudaqEvent->IsEORE() ) {
      // this is the case in which the eudaq event is a End Of Run
      // Event
      processEORE( eudaqEvent.get() );

    } else {
      // if we get to here, it means that this should be a data event,
      // but it may be empty. So before proceeding, check if there is
      // some data inside
      eudaq::DetectorEvent * eudaqDetectorEvent = dynamic_cast<eudaq::DetectorEvent *> ( eudaqEvent.get() );
      if ( !eudaqDetectorEvent ) {
        // the event is empty... strange but possible! So skip it!
        streamlog_out( WARNING2 ) << "Event number " << eudaqEvent->GetEventNumber() << " does not contain any data " << endl;

      } else {
        // finally this is the case in which we have some data to
        // process...
        auto_ptr<EUTelEventImpl> eutelEvent(new EUTelEventImpl);
        eutelEvent->setRunNumber  ( eudaqEvent->GetRunNumber() );
        eutelEvent->setEventNumber( eudaqEvent->GetEventNumber() );
        eutelEvent->setTimeStamp  ( eudaqEvent->GetTimestamp() );
        eutelEvent->setEventType  ( kDE );

        // we have to process all the different kind of subevent (one
        // for each data producer. Start here the loop...
        for ( size_t iProducer = 0; iProducer < eudaqDetectorEvent->NumEvents(); ++iProducer ) {
          // get the current sub event
          eudaq::Event * subEvent = eudaqDetectorEvent->GetEvent( iProducer );

          // once again we have to guess the event type making several
          // dynamic recasting
          //
          // let's start from the EUDRBEvents
          eudaq::EUDRBEvent * eudrbEvent = dynamic_cast< eudaq::EUDRBEvent *> ( subEvent );
          if ( eudrbEvent ) {
            // ok, this is an eudrb event, get ready to process it!
            processEUDRBDataEvent( eudrbEvent, eutelEvent.get() ) ;
          }

          eudaq::TLUEvent * tluEvent = dynamic_cast< eudaq::TLUEvent * > ( subEvent );
          if ( tluEvent ) {
            // this is a TLU event, process it!
            processTLUDataEvent( tluEvent, eutelEvent.get() );
          }

          /* Add here your event using the same procedure as above */


 
          eudaq::DEPFETEvent * depfetEvent = dynamic_cast< eudaq::DEPFETEvent * > ( subEvent );
          if ( depfetEvent ) {
            // this is a DEPFET event, process it!
	      printf("DEPFET event found!!! \n ");
            processDEPFETDataEvent( depfetEvent, eutelEvent.get() );
          }


        }

        // all the producers have been added to the LCIO event, so we
        // can process it. The processEvent wants to have a real
        // pointer, so I have to release the auto_ptr but I have to
        // keep it in mind to delete it afterwards
        EUTelEventImpl * dummyEvent = eutelEvent.release();
        ProcessorMgr::instance()->processEvent( static_cast< lcio::LCEventImpl *> ( dummyEvent ) );
        delete dummyEvent;
      }

      // in case this was the first event, toggle the boolean
      if ( _isFirstEvent ) _isFirstEvent = false;
    }

    // increment the event counter here. This number is not used to
    // set the event number but only to count how many events have
    // been processed to stop the conversion
    ++eventCounter;

  }
}


void EUTelNativeReader::processEUDRBDataEvent( eudaq::EUDRBEvent * eudrbEvent, EUTelEventImpl * eutelEvent ) {
  // here we have to process the EUDRB event
  //
  // prepare a collection for the raw data and one for the zs
  auto_ptr< lcio::LCCollectionVec > rawDataCollection ( new LCCollectionVec (LCIO::TRACKERRAWDATA) ) ;
  auto_ptr< lcio::LCCollectionVec > zsDataCollection  ( new LCCollectionVec (LCIO::TRACKERDATA) ) ;

  // set the proper cell encoder
  CellIDEncoder< TrackerRawDataImpl > rawDataEncoder ( EUTELESCOPE::MATRIXDEFAULTENCODING, rawDataCollection.get() );
  CellIDEncoder< TrackerDataImpl    > zsDataEncoder  ( EUTELESCOPE::ZSDATADEFAULTENCODING, zsDataCollection.get()  );

  // if this is the first event we also have to prepare a collection
  // of EUTelSetupDescription with the info of all the detectors
  // readout by the EUDRBs
  
  if ( isFirstEvent() ) {

    auto_ptr< lcio::LCCollectionVec > eudrbSetupCollection( new LCCollectionVec (LCIO::LCGENERICOBJECT) );

    for ( size_t iDetector = 0; iDetector < _eudrbDetectors.size() ; iDetector++) {
      EUTelSetupDescription * detector = new EUTelSetupDescription( _eudrbDetectors.at( iDetector ) );
      eudrbSetupCollection->push_back( detector );
    }
    eutelEvent->addCollection( eudrbSetupCollection.release(), "eudrbSetup" );
  }

  // to understand if we have problem with de-syncronisation, let
  // me prepare a Boolean switch and a vector of size_t to contain the
  // pivot pixel position
  bool outOfSyncFlag = false;
  vector<size_t > pivotPixelPosVec;

  // we can now loop over the boards contained into this EUDRB event
  for ( size_t iPlane = 0; iPlane < eudrbEvent->NumBoards(); ++iPlane ) {

    // I need to have a EUDRBBoard
    eudaq::EUDRBBoard& eudrbBoard = eudrbEvent->GetBoard( iPlane );

    // get from the eudrb detectors the current one
    EUTelPixelDetector * currentDetector = _eudrbDetectors.at( iPlane );

    // from now on we have to proceed in a different way depending if
    // the sensor was readout in RAW mode or in ZS
    string currentMode = currentDetector->getMode();

    if (  ( currentMode == "RAW" ) || ( currentMode == "RAW2" ) || ( currentMode == "RAW3" ) ) {

      // ----------------------------------------------------------------------------------------------------
      //
      //                                           R A W   M O D E
      //
      // ----------------------------------------------------------------------------------------------------

      rawDataEncoder["xMin"]     = currentDetector->getXMin();
      rawDataEncoder["xMax"]     = currentDetector->getXMax() - currentDetector->getMarkerPosition().size();
      rawDataEncoder["yMin"]     = currentDetector->getYMin();
      rawDataEncoder["yMax"]     = currentDetector->getYMax();
      rawDataEncoder["sensorID"] = iPlane;

      // put a try/catch box here mainly for the EUDRBDecoder::GetArrays
      try {

        // get arrays may throw a eudaq::Exception in case the array
        // sizes are different. In this case we should catch the
        // exception, inform the user and then skip the event
        eudaq::EUDRBDecoder::arrays_t<short, short > array = _eudrbDecoder->GetArrays<short, short > ( eudrbBoard );

        // this is the right point to split the processing of RAW2
        // from the one of RAW3
        if ( currentMode == "RAW3" ) {

          // ----------------------------------------------------------------------------------------------------
          //
          //                                           R A W 3   M O D E
          //
          // ----------------------------------------------------------------------------------------------------

          // let's start from the case of RAW3. I need to get the
          // three arrays of short from the arrays_t struct and to
          // create a vector of short for the CDS. The CDS calculation
          // is done on the full matrix (markers included) and
          // afterwards the markers are removed.
          vector<short > cdsValueVec;
          vector<short > firstFrameVec  = array.m_adc[0];
          vector<short > secondFrameVec = array.m_adc[1];
          vector<short > thirdFrameVec  = array.m_adc[2];

          // ready to make the calculation
          for ( size_t iPixel = 0; iPixel < (unsigned) currentDetector->getXNoOfPixel() * currentDetector->getYNoOfPixel(); ++iPixel ) {

            // in RAW3 mode the CDS calculation is done using the
            // pivot pixel
            cdsValueVec.push_back( currentDetector->getSignalPolarity() *
                                   ( ( -1 * array.m_pivot[iPixel]     ) * firstFrameVec[iPixel]  +
                                     (  2 * array.m_pivot[iPixel] - 1 ) * secondFrameVec[iPixel] +
                                     (  1 - array.m_pivot[iPixel]     ) * thirdFrameVec[iPixel]   ) );

          }

          // now we have to strip out the marker cols from the CDS
          // value. To do this I need a vector of short large enough
          // to accommodate the full matrix without the markers
          vector<short > cdsStrippedVec( currentDetector->getYNoOfPixel() * ( currentDetector->getXNoOfPixel() - currentDetector->getMarkerPosition().size() ) );

          // I need also two iterators, one for the stripped vec and
          // one for the original one.
          vector<short >::iterator currentCDSPos = cdsStrippedVec.begin();
          vector<short >::iterator cdsBegin      = cdsValueVec.begin();

          // now loop over all the pixels
          for ( size_t y = 0; y < currentDetector->getYNoOfPixel(); ++y ) {
            size_t offset = y * currentDetector->getXNoOfPixel();
            vector<size_t >::iterator marker = currentDetector->getMarkerPosition().begin();

            // first copy from the beginning of the row to the first
            // marker column
            currentCDSPos = copy( cdsBegin + offset, cdsBegin + ( *(marker) + offset ), currentCDSPos );

            // now copy from the next column to the next marker into a
            // while loop
            while ( marker != currentDetector->getMarkerPosition().end() ) {
              if ( marker < currentDetector->getMarkerPosition().end() - 1 ) {
                currentCDSPos = copy( cdsBegin + ( *(marker) + 1 + offset ), cdsBegin + ( *(marker + 1) + offset ), currentCDSPos );
              } else {
                // now from the last marker column to the end of the
                // row
                currentCDSPos = copy( cdsBegin + ( *(marker) + 1 + offset ), cdsBegin + offset + currentDetector->getXNoOfPixel(), currentCDSPos );
              }
              ++marker;
            }
          }

          // this is the right place to prepare the TrackerRawData
          // object
          auto_ptr< lcio::TrackerRawDataImpl > cdsFrame( new lcio::TrackerRawDataImpl );
          rawDataEncoder.setCellID( cdsFrame.get() );

          // add the cds stripped values to the current TrackerRawData
          cdsFrame->setADCValues( cdsStrippedVec ) ;

          // put the pivot pixel in the timestamp field of the
          // TrackerRawData. I know that is not correct, but this is
          // the only place where I can put this info
          cdsFrame->setTime( eudrbBoard.PivotPixel() );

          // this is also the right place to add the pivot pixel to
          // the pivot pixel vector for synchronization checks
          pivotPixelPosVec.push_back( eudrbBoard.PivotPixel() );

          // now append the TrackerRawData object to the corresponding
          // collection releasing the auto pointer
          rawDataCollection->push_back( cdsFrame.release() );

        } else if ( currentMode == "RAW2" ) {

          // ----------------------------------------------------------------------------------------------------
          //
          //                                           R A W 2   M O D E
          //
          // ----------------------------------------------------------------------------------------------------

          // the procedure is quite similar to the one in RAW3, simply
          // we need only 2 vectors instead of three and the CDS
          // calculation is made easier
          vector<short > cdsValueVec;
          vector<short > firstFrameVec  = array.m_adc[0];
          vector<short > secondFrameVec = array.m_adc[1];

          // ready to make the calculation
          for ( size_t iPixel = 0; iPixel < (unsigned) currentDetector->getXNoOfPixel() * currentDetector->getYNoOfPixel(); ++iPixel ) {

            // in RAW2 the two arrays are already sorted to make the
            // proper difference
            cdsValueVec.push_back( currentDetector->getSignalPolarity() * ( secondFrameVec[ iPixel ] - firstFrameVec[ iPixel ] ) );

          }

          // now we have to strip out the marker cols from the CDS
          // value. To do this I need a vector of short large enough
          // to accommodate the full matrix without the markers
          vector<short > cdsStrippedVec( currentDetector->getYNoOfPixel() * ( currentDetector->getXNoOfPixel() - currentDetector->getMarkerPosition().size() ) );

          // I need also two iterators, one for the stripped vec and
          // one for the original one.
          vector<short >::iterator currentCDSPos = cdsStrippedVec.begin();
          vector<short >::iterator cdsBegin      = cdsValueVec.begin();

          // now loop over all the pixels
          for ( size_t y = 0; y < currentDetector->getYNoOfPixel(); ++y ) {
            size_t offset = y * currentDetector->getXNoOfPixel();
            vector<size_t >::iterator marker = currentDetector->getMarkerPosition().begin();

            // first copy from the beginning of the row to the first
            // marker column
            currentCDSPos = copy( cdsBegin + offset, cdsBegin + ( *(marker) + offset ), currentCDSPos );

            // now copy from the next column to the next marker into a
            // while loop
            while ( marker != currentDetector->getMarkerPosition().end() ) {
              if ( marker < currentDetector->getMarkerPosition().end() - 1 ) {
                currentCDSPos = copy( cdsBegin + ( *(marker) + 1 + offset ), cdsBegin + ( *(marker + 1) + offset ), currentCDSPos );
              } else {
                // now from the last marker column to the end of the
                // row
                currentCDSPos = copy( cdsBegin + ( *(marker) + 1 + offset ), cdsBegin + offset + currentDetector->getXNoOfPixel(), currentCDSPos );
              }
              ++marker;
            }
          }

          // this is the right place to prepare the TrackerRawData
          // object
          auto_ptr< lcio::TrackerRawDataImpl > cdsFrame( new lcio::TrackerRawDataImpl );
          rawDataEncoder.setCellID( cdsFrame.get() );

          // add the cds stripped values to the current TrackerRawData
          cdsFrame->setADCValues( cdsStrippedVec ) ;

          // put the pivot pixel in the timestamp field of the
          // TrackerRawData. I know that is not correct, but this is
          // the only place where I can put this info
          cdsFrame->setTime( eudrbBoard.PivotPixel() );

          // this is also the right place to add the pivot pixel to
          // the pivot pixel vector for synchronization checks
          pivotPixelPosVec.push_back( eudrbBoard.PivotPixel() );

          // now append the TrackerRawData object to the corresponding
          // collection releasing the auto pointer
          rawDataCollection->push_back( cdsFrame.release() );

        } else {
          // there shouldn't be any memory leak because of
          // auto_ptr. Since we are inside a loop over the boards and
          // other boards can be operated with a different modality,
          // it is better to skip the full eudrb event returning here
          streamlog_out ( ERROR0 ) << "The current mode " << currentMode << " has been recognised as RAW mode, but unable to process it" << endl;
          return;
        }

      } catch ( eudaq::Exception& e) {
        // very luckily this is a problem with the GetArrays. It means
        // that one or more boards have been badly decoded, so it is
        // better skipping the full events returning here
        streamlog_out( ERROR ) << e.what() << endl << "Skipping the current event " << endl;
        return;
      }


    } else if ( ( currentMode == "ZS" ) ) {

      // ----------------------------------------------------------------------------------------------------
      //
      //                                            Z S   M O D E
      //
      // ----------------------------------------------------------------------------------------------------

      zsDataEncoder["sensorID"] = iPlane;
      zsDataEncoder["sparsePixelType"] = _eudrbSparsePixelType;

      // get the total number of pixel. This is written in the
      // eudrbBoard and to get it in a easy way pass through the eudrbDecoder
      size_t nPixel = _eudrbDecoder->NumPixels( eudrbBoard );

      // prepare a new TrackerData for the ZS data
      auto_ptr<lcio::TrackerDataImpl > zsFrame( new lcio::TrackerDataImpl );
      zsDataEncoder.setCellID( zsFrame.get() );

      // depending on the sparse pixel type, now you have to do
      // something different here. For the time being, only one kind
      // of sparse pixel exists, so it is easy
      if ( _eudrbSparsePixelType == kEUTelSimpleSparsePixel ) {

        // this is the structure that will host the sparse pixel
        auto_ptr<EUTelSparseDataImpl<EUTelSimpleSparsePixel > >
          sparseFrame( new EUTelSparseDataImpl<EUTelSimpleSparsePixel > ( zsFrame.get() ) );

        // ready to get the data, but put the GetArrays into a
        // try/catch box
        try {

          eudaq::EUDRBDecoder::arrays_t<short, short > array = _eudrbDecoder->GetArrays<short, short > ( eudrbBoard );

          // prepare a sparse pixel to be added to the sparse data
          auto_ptr<EUTelSimpleSparsePixel > sparsePixel( new EUTelSimpleSparsePixel );
          for ( size_t iPixel = 0; iPixel < nPixel; ++iPixel ) {

            // the data contain also the markers, so we have to strip
            // them out. First I need to have the original position
            // (with markers in) and then calculate how many pixels I
            // have to remove
            size_t originalX = array.m_x[ iPixel ] ;
	    vector<size_t > markerVec = currentDetector->getMarkerPosition();

            if ( find( markerVec.begin(), markerVec.end(), originalX ) == markerVec.end() ) {
              // the original X is not on a marker column, so I need
              // to remove a certain number of pixels depending on the
              // position

              // this counts the number of markers found on the left
              // of the original X
              short  diff = ( short ) count_if ( markerVec.begin(),markerVec.end(), bind2nd(less<short> (), originalX ) );
              sparsePixel->setXCoord( originalX - diff );

              // no problem instead with the Y coordinate
              sparsePixel->setYCoord( array.m_y[ iPixel ] );

              // last the pixel charge. The CDS is automatically
              // calculated by the EUDRB
              sparsePixel->setSignal( array.m_adc[0][ iPixel ] );

              // in case of DEBUG
              streamlog_out ( DEBUG0 ) << ( *(sparsePixel.get() ) ) << endl;

              // now add this pixel to the sparse frame
              sparseFrame->addSparsePixel( sparsePixel.get() );
            } else {
              // the original X was on a marker column, so we don't
              // need to process this pixel any further and of course
              // we don't have to add it to the sparse frame.
              streamlog_out ( DEBUG0 ) << "Found a sparse pixel ("<< iPixel
                                       <<")  on a marker column. Not adding it to the frame" << endl
                                       << (* (sparsePixel.get() ) ) << endl;
            }

          }

          // perfect! Now add the TrackerData to the collection
          zsDataCollection->push_back( zsFrame.release() );

          // for the debug of the synchronization
          pivotPixelPosVec.push_back( eudrbBoard.PivotPixel() );

        } catch ( eudaq::Exception & e ) {
          // very luckily this is a problem with the GetArrays. It means
          // that one or more boards have been badly decoded, so it is
          // better skipping the full events returning here
          streamlog_out( ERROR ) << e.what() << endl << "Skipping the current event " << endl;
          return;
        }
      }
      /*

        ADD HERE ANY NEW READOUT MODE
        } else if ( currentMode == "myMode" ) {
        ...

      */
    } else {
      // leave this final else in case there is a modality that is not
      // recognized.
      streamlog_out ( ERROR ) << "The specified readout mode (" << currentMode << ") is unknown. Skipping the current event " << endl;
      return;

    }

  } // end of the loop over all the planes into the telescope

  // check if all the boards where running in synchronous mode or
  // not. Remember that the last pivot pixel entry is the one of the
  // master board.
  vector<size_t >::iterator masterBoardPivotAddress = pivotPixelPosVec.end() - 1;
  vector<size_t >::iterator slaveBoardPivotAddress  = pivotPixelPosVec.begin();
  while ( slaveBoardPivotAddress < masterBoardPivotAddress ) {
    if ( *slaveBoardPivotAddress - *masterBoardPivotAddress >=  _eudrbOutOfSyncThr ) {
      outOfSyncFlag = true;

      // we don't need to continue looping over all boards if one of
      // them is already out of sync
      break;
    }
    ++slaveBoardPivotAddress;
  }
  if ( outOfSyncFlag ) {

    // Increase the total out of sync counter
    ++_eudrbTotalOutOfSyncEvent;

    // in case the current event is out of sync then we have to tell
    // it to the user, but not too often!
    unsigned int currentEventNumber = eutelEvent->getEventNumber();
    if ( currentEventNumber == _eudrbPreviousOutOfSyncEvent + 1 ) {
      // if the previous event out of sync was the previous event,
      // then we have to increment the counter
      ++_eudrbConsecutiveOutOfSyncWarning;
      _eudrbPreviousOutOfSyncEvent = currentEventNumber;
    } else {
      // otherwise we can reset the warning counter
      _eudrbConsecutiveOutOfSyncWarning = 0;
    }

    if ( _eudrbConsecutiveOutOfSyncWarning < _eudrbMaxConsecutiveOutOfSyncWarning ) {
      // in this case we have the responsibility to tell the user that
      // the event was out of sync
      streamlog_out ( WARNING0 ) << "Event number " << eutelEvent->getEventNumber() << " seems to be out of sync" << endl;
      vector<size_t >::iterator masterBoardPivotAddress = pivotPixelPosVec.end() - 1;
      vector<size_t >::iterator slaveBoardPivotAddress  = pivotPixelPosVec.begin();
      while ( slaveBoardPivotAddress < masterBoardPivotAddress ) {
        // print out all the slave boards first
        streamlog_out( WARNING0 ) << " --> Board (S) " <<  setw(3) << setiosflags( ios::right )
                                  << slaveBoardPivotAddress - pivotPixelPosVec.begin() << resetiosflags( ios::right )
                                  << " = " << setw(15) << setiosflags( ios::right )
                                  << (*slaveBoardPivotAddress) << resetiosflags( ios::right )
                                  << " (" << setw(15) << setiosflags( ios::right )
                                  << (signed) (*masterBoardPivotAddress) - (signed) (*slaveBoardPivotAddress) << resetiosflags(ios::right)
                                  << ")" << endl;
        ++slaveBoardPivotAddress;
      }
      // print out also the master. It is impossible that the master
      // is out of sync with respect to itself, but for completeness...
      streamlog_out( WARNING0 )  << " --> Board (M) "  <<  setw(3) << setiosflags( ios::right )
                                 << slaveBoardPivotAddress - pivotPixelPosVec.begin() << resetiosflags( ios::right )
                                 << " = " << setw(15) << setiosflags( ios::right )
                                 << (*slaveBoardPivotAddress) << resetiosflags( ios::right )
                                 << " (" << setw(15)  << setiosflags( ios::right )
                                 << (signed) (*masterBoardPivotAddress) - (signed) (*slaveBoardPivotAddress) << resetiosflags(ios::right)
                                 << ")" << endl;

    } else if ( _eudrbConsecutiveOutOfSyncWarning == _eudrbMaxConsecutiveOutOfSyncWarning ) {
      // if the number of consecutive warnings is equal to the maximum
      // allowed, don't bother the user anymore with this message,
      // because it's very luckily the run was taken unsynchronized on
      // purpose
      streamlog_out( WARNING0 ) << "The maximum number of consecutive unsychronized events has been reached." << endl
                                << "Assuming the run was taken in asynchronous mode" << endl;
    }

  }

  // before leaving we have to add the two collections (raw and zs to
  // the current event, but we do it only with not empty collections.
  if ( rawDataCollection->size() ) {
    // we have some rawdata...
    eutelEvent->addCollection( rawDataCollection.release(), _eudrbRawModeOutputCollectionName );
  }

  if ( zsDataCollection->size() ) {
    // we have some ZS data...
    eutelEvent->addCollection( zsDataCollection.release(), _eudrbZSModeOutputCollectionName );
  }

}

void EUTelNativeReader::processTLUDataEvent( eudaq::TLUEvent * tluEvent, EUTelEventImpl * eutelEvent ) {


}
void EUTelNativeReader::processDEPFETDataEvent( eudaq::DEPFETEvent * depfetEvent, EUTelEventImpl * eutelEvent ) {
    auto_ptr< lcio::LCCollectionVec > rawDataCollection ( new LCCollectionVec (LCIO::TRACKERRAWDATA) ) ;

  CellIDEncoder< TrackerRawDataImpl > rawDataEncoder ( EUTELESCOPE::MATRIXDEFAULTENCODING, rawDataCollection.get() );
  int DATA[64][128];
    printf("DEPFET DATA Event \n");
  
    if(isFirstEvent()) {
	printf("DEPFET DATA Event first 1\n");
	
	auto_ptr< lcio::LCCollectionVec > depfetSetupCollection( new LCCollectionVec (LCIO::LCGENERICOBJECT) );
	
	eutelEvent->addCollection( depfetSetupCollection.release(), "depfetSetup" );
    }
    
    printf("DEPFET DATA Event first 2=%d \n",_depfetDetectors.size());

    for ( size_t iDetector = 0; iDetector < _depfetDetectors.size() ; iDetector++) {
	printf("DEPFET DATA Event loop %d numboards=%d \n",iDetector,depfetEvent->NumBoards());     
	for(size_t iPlane=0; iPlane<depfetEvent->NumBoards(); ++iPlane){ 
	    printf("DEPFET DATA:: iPlane loop %d \n",iPlane);     
	    
	    eudaq::DEPFETBoard&depfetBoard=depfetEvent->GetBoard(iPlane);
	    eudaq::DEPFETDecoder::arrays_t<short, short > array = _depfetDecoder->GetArrays<short, short > ( depfetBoard );
	    printf("DEPFET DATA:: GetBoard\n");     
	    
	    rawDataEncoder["sensorID"] = 8;
	    rawDataEncoder["xMin"] = 0;
	    rawDataEncoder["xMax"] = 64 - 1;
	    rawDataEncoder["yMin"] = 0;
	    rawDataEncoder["yMax"] = 128 - 1;
            printf("call for get raw data\n");
	    auto_ptr< lcio::TrackerRawDataImpl > rawData( new lcio::TrackerRawDataImpl );
            printf("call for setcell id \n");
	    rawDataEncoder.setCellID(rawData.get());
	    for (int i = 0; i < 8192; i++) {
 	 
		    DATA[array.m_x[i]][array.m_y[i]]=array.m_adc[0][i];
		    //	    printf("call for pushback ipix=%d  %d %d %d \n",i,array.m_x[i],array.m_y[i],array.m_adc[0][i] );
//			if(iprint==6)printf("Det=%d ADCvalue=%d,  ipix=%d %d,\n",i,DATA[xPixel][yPixel],xPixel,yPixel);
		
	    }
          
	    for (int yPixel = 0; yPixel < 128; yPixel++) {
		for (int xPixel = 0; xPixel < 64; xPixel++) {
		    rawData->adcValues().push_back (DATA[xPixel][yPixel]);
		}
	    }
	    

	    rawDataCollection->push_back(rawData.release());
	}

    }
    printf("rawDataCollection size= %d \n",rawDataCollection->size() );
    if ( rawDataCollection->size() ) {
	// we have some rawdata...
	eutelEvent->addCollection( rawDataCollection.release(), "DEPFET" );
    }

}


void EUTelNativeReader::processEORE( eudaq::Event * eore ) {
  streamlog_out( DEBUG4 ) << "Found a EORE, so adding an EORE to the LCIO file as well" << endl;
  EUTelEventImpl * lcioEvent = new EUTelEventImpl;
  lcioEvent->setEventType(kEORE);
  lcioEvent->setTimeStamp  ( eore->GetTimestamp()   );
  lcioEvent->setRunNumber  ( eore->GetRunNumber()   );
  lcioEvent->setEventNumber( eore->GetEventNumber() );

  // sent the lcioEvent to the processor manager for further processing
  ProcessorMgr::instance()->processEvent( static_cast<LCEventImpl*> ( lcioEvent ) ) ;
  delete lcioEvent;
}

void EUTelNativeReader::processBORE( eudaq::Event * bore ) {
  streamlog_out( DEBUG4 ) << "Found a BORE, so processing the RunHeader" << endl;

  // prepare the LCIORunHeader and to make the assignment easier use
  // the EUTelRunHeader decorator pattern.
  auto_ptr<IMPL::LCRunHeaderImpl> lcHeader( new IMPL::LCRunHeaderImpl );
  auto_ptr<EUTelRunHeaderImpl>    runHeader( new EUTelRunHeaderImpl( lcHeader.get() ) );

  runHeader->lcRunHeader()->setRunNumber( bore->GetRunNumber() ) ;

  unsigned short noOfEUDRBDetectors = 0;
  unsigned short noOfTLUDetectors   = 0;
  unsigned short noOfDEPFETDetectors   = 0;

  // add here your own noOf
  // ...


  // in the BORE we have as many sub events as the number of
  // detectors into the full DAQ setup (telescope + TLUs + DUTs
  // ...), first we need to convert the general event into a
  // detector event
  eudaq::DetectorEvent * eudaqDetectorEvent = dynamic_cast< eudaq::DetectorEvent *> ( bore );
  for ( size_t iProducer = 0 ; iProducer < eudaqDetectorEvent->NumEvents(); ++iProducer ) {
    eudaq::Event * eudaqSubEvent = eudaqDetectorEvent->GetEvent( iProducer );

    // now we have to guess which type of sub event is going thourg
    // a successfull dynamic casting
    eudaq::EUDRBEvent * eudaqEUDRBEvent = dynamic_cast< eudaq::EUDRBEvent *> ( eudaqSubEvent );

    if ( eudaqEUDRBEvent ) {
      // ok this is the sub event produced by the EUDRBProducer.

      // first of all set the global detector type and mode
      string globalMode = eudaqEUDRBEvent->GetTag( "MODE" );
      string globalDet  = eudaqEUDRBEvent->GetTag( "DET" );
      runHeader->setEUDRBMode( globalMode );
      runHeader->setEUDRBDet( globalDet );

      // now loop over all the detectors in the eudrb producer
      unsigned short localNoOfEUDRBDetectors =  eudaq::from_string( eudaqEUDRBEvent->GetTag("BOARDS"), 0 );
      for ( size_t iDetector = 0 ; iDetector < localNoOfEUDRBDetectors; ++iDetector ) {
        string detectorFlag = "DET" + to_string( iDetector );
        string detectorType = eudaqEUDRBEvent->GetTag( detectorFlag.c_str() );
        string modeFlag = "MODE" + to_string( iDetector );
        string modeType = eudaqEUDRBEvent->GetTag( modeFlag.c_str() );
        if ( detectorType.compare( "MIMOTEL" ) == 0 ) {
          EUTelMimoTelDetector * detector = new EUTelMimoTelDetector;
          detector->setMode( modeType );
          _eudrbDetectors.push_back ( detector );
        } else if ( detectorType.compare( "MIMOSA18" ) == 0 ) {
          EUTelMimosa18Detector * detector = new EUTelMimosa18Detector;
          detector->setMode( modeType );
          _eudrbDetectors.push_back ( detector );
/*
  Add here your own detector type according to the following example
  } else if ( detectorType.compare("MyDetector") == 0 ) {
*/
        } else {
          throw InvalidParameterException("For the time being only the following detectors are supported: \n\n"
                                          "  1. MimoTEL \n"
                                          "  2. Mimosa18 \n");
        }
      }
      noOfEUDRBDetectors += localNoOfEUDRBDetectors;
      assert( noOfEUDRBDetectors == _eudrbDetectors.size());

      // before leaving, remember to assign the _eudrbDecoder
      _eudrbDecoder = new eudaq::EUDRBDecoder( *eudaqDetectorEvent );


    } // this is the end of the EUDRBEvent

      // try now with a TLU event
    eudaq::TLUEvent * eudaqTLUEvent = dynamic_cast< eudaq::TLUEvent * > ( eudaqSubEvent ) ;

    if ( eudaqTLUEvent ) {
      // this sub event was produced by the TLUProducer
      // nothing to do for the time being
      // ....
      EUTelTLUDetector * tluDetector = new EUTelTLUDetector;
      tluDetector->setAndMask ( from_string( eudaqTLUEvent->GetTag( "AndMask"  ), (unsigned int) 0x0) );
      tluDetector->setOrMask  ( from_string( eudaqTLUEvent->GetTag( "OrMask"   ), (unsigned int) 0x0) );
      tluDetector->setVetoMask( from_string( eudaqTLUEvent->GetTag( "VetoMask" ), (unsigned int) 0x0) );
      tluDetector->setDUTMask ( from_string( eudaqTLUEvent->GetTag( "DutMask"  ), (unsigned int) 0x0) );
      tluDetector->setFirmwareID  ( from_string( eudaqTLUEvent->GetTag( "FirmwareID" ), (unsigned int) 0x0) );
      tluDetector->setTimeInterval( from_string( eudaqTLUEvent->GetTag( "TimeInterval" ), (unsigned short) 0x0) );
      _tluDetectors.push_back( tluDetector );
      noOfTLUDetectors++;
    }

/*
  Add here your own data producer following the examples above
*/
    eudaq::DEPFETEvent * eudaqDEPFETEvent = dynamic_cast< eudaq::DEPFETEvent * > ( eudaqSubEvent ) ;

    if ( eudaqDEPFETEvent ) {
      // this sub event was produced by the TLUProducer
      // nothing to do for the time being
      // ....
	printf("JFadd====> DEPFET BOR Event found\n");
     // now loop over all the detectors in the eudrb producer
   
/*     unsigned short localNoOfDEPFETDetectors =  eudaq::from_string( eudaqDEPFETEvent->GetTag("BOARDS"), 0 );
 	printf("JFadd====> DEPFET BOR Event found det=%d\n",localNoOfDEPFETDetectors );
	for ( size_t iDetector = 0 ; iDetector < localNoOfDEPFETDetectors; ++iDetector ) {

	}
        noOfDEPFETDetectors +=localNoOfDEPFETDetectors ;
*/
	
	EUTelDEPFETDetector * detector = new EUTelDEPFETDetector;
	_depfetDetectors.push_back( detector );     
//	assert( noOfDEPFETDetectors == _depfetDetectors.size());
	
      // before leaving, remember to assign the _eudrbDecoder
        _depfetDecoder = new eudaq::DEPFETDecoder( *eudaqDetectorEvent );


  }


  }
  // print a summary of all recognized producers:
  // first the EUDRB producer
  if ( _eudrbDetectors.size() == 0 ) {
    streamlog_out( DEBUG0 ) << "No EUDRB producer found in this BORE" << endl;
  } else {
    streamlog_out ( MESSAGE2 ) << "EUDRBProducer: " << endl << endl;
    for ( size_t iDetector = 0 ; iDetector < _eudrbDetectors.size(); ++iDetector ) {
      streamlog_out ( MESSAGE2 ) << " Detector " << iDetector << endl
                                 << *(_eudrbDetectors.at( iDetector )) << endl
                                 << " -------------------------------------------- " << endl;
    }
  }

  // then the TLU producer
  if ( _tluDetectors.size() == 0 ) {
    streamlog_out( DEBUG0 ) << "No TLU producer found in this EORE" << endl;
  } else {
    streamlog_out( MESSAGE2 ) << "TLUProducer: " << endl << endl;
    for ( size_t iDetector = 0 ; iDetector < _tluDetectors.size(); ++iDetector ) {
      streamlog_out ( MESSAGE2 ) << " TLU " << iDetector << " found " << endl
                                 << *(_tluDetectors.at( iDetector ) ) << endl
                                 << " -------------------------------------------- " << endl;
    }

  }
  // then the DEPFET producer
  if ( _depfetDetectors.size() == 0 ) {
      printf("depfet detector size=0\n");
      streamlog_out( DEBUG0 ) << "No DEPFET producer found in this EORE" << endl;
  }

  // before processing the run header we have to set the important
  // information about the telescope (i.e. the EUDRBProcuder)
  vector< int > xMin, xMax, yMin, yMax;
  for ( size_t iPlane = 0; iPlane < _eudrbDetectors.size(); ++iPlane ) {
    xMin.push_back( _eudrbDetectors.at( iPlane )->getXMin() );
    xMax.push_back( _eudrbDetectors.at( iPlane )->getXMax() - _eudrbDetectors.at( iPlane )->getMarkerPosition().size() );
    yMin.push_back( _eudrbDetectors.at( iPlane )->getYMin() );
    yMax.push_back( _eudrbDetectors.at( iPlane )->getYMax() );
  }

  runHeader->setNoOfDetector( _eudrbDetectors.size() );
  runHeader->setMinX( xMin );
  runHeader->setMaxX( xMax );
  runHeader->setMinY( yMin );
  runHeader->setMaxY( yMax );
  runHeader->lcRunHeader()->setDetectorName("EUTelescope");
  runHeader->setGeoID( _geoID );
  runHeader->setDateTime();
  runHeader->addProcessor( type() );
  ProcessorMgr::instance()->processRunHeader(lcHeader.get());
  delete lcHeader.release();


}

void EUTelNativeReader::end () {

  if ( _eudrbTotalOutOfSyncEvent != 0 ) {
    streamlog_out ( MESSAGE ) << "Found " << _eudrbTotalOutOfSyncEvent << " events out of sync " << endl;
  }
  streamlog_out ( MESSAGE)  << "Successfully finished" << endl;
}

#endif


//  LocalWords:  serialiazer EUTelNativeReader MIMOTEL rawdata eudrb zsdata
//  LocalWords:  EUBRDRawModeOutputCollection EUDRBZSModeOutputCollection xMin
//  LocalWords:  EUDRBSparsePixelType sensorID xMax yMin yMax EUDRBBoard
//  LocalWords:  EUDRBDecoder TrackerRawData
