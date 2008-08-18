// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Author Loretta Negrini, Univ. Insubria <mailto:loryneg@gmail.com>
// Author Silvia Bonfanti, Univ. Insubria <mailto:silviafisica@gmail.com>
// Version $Id: EUTelNativeReader.cc,v 1.2 2008-08-18 15:04:58 bulgheroni Exp $
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
#include "EUTelMimoTelDetector.h"
#include "EUTelMimosa18Detector.h"

// marlin includes
#include "marlin/Exceptions.h"
#include "marlin/Processor.h"
#include "marlin/DataSourceProcessor.h"
#include "marlin/ProcessorMgr.h"

// eudaq includes
#include <eudaq/FileSerializer.hh>
#include <eudaq/EUDRBEvent.hh>
#include <eudaq/DetectorEvent.hh>
#include <eudaq/Utils.hh>
#include <eudaq/Exception.hh>

// lcio includes
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCEventImpl.h>
#include <UTIL/CellIDEncoder.h>

// system includes
#include <iostream>

using namespace std;
using namespace marlin;
using namespace eutelescope;

// initialize static members
const unsigned short EUTelNativeReader::_eudrbOutOfSynchThr = 2;
const unsigned short EUTelNativeReader::_eudrbMaxConsecutiveOutOfSynchWarning = 20;

EUTelNativeReader::EUTelNativeReader (): DataSourceProcessor  ("EUTelNativeReader") {

  // initialize few variables
  _eudrbPreviousOutOfSynchEvent = 0;

  _description =
    "Reads data streams produced by EUDAQ and produced the corresponding LCIO output";

  registerProcessorParameter("InputFileName", "This is the input file name",
                             _fileName, string("run012345.raw") );

  registerProcessorParameter("GeoID", "The geometry identification number", _geoID, static_cast<int> ( 0 ));


  // from here below only detector specific parameters.

  // ---------- //
  //  MIMOTEL   //
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
                                 << setw(6) << setiosflags(ios::right) << eudaqEvent->GetEventNumber() << " in run "
                                 << setw(6) << setiosflags(ios::right) << setfill('0')  << eudaqEvent->GetRunNumber() << setfill(' ')
                                 << " (Total = " << setw(10) << eventCounter << ")" << resetiosflags(ios::left) << endl;


    // increment the event counter here. This number is not used to
    // set the event number but only to count how many events have
    // been processed to stop the conversion
    ++eventCounter;

    if ( eventCounter >= numEvents ) {
      // even if there are some more events in the input file, we
      // don't want to decode them. So add a kEORE
      // and then stop the data reading
      processEORE( eudaqEvent.get() );

      break; // this is breaking the loop on input evets
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
        }

        // all the producers have been added to the LCIO event, so we
        // can process it. The processEvent wants to have a real
        // pointer, so I have to release the auto_ptr but I have to
        // keep it in mind to delete it afterwards
        EUTelEventImpl * dummyEvent = eutelEvent.release();
        ProcessorMgr::instance()->processEvent( static_cast< LCEventImpl *> ( dummyEvent ) );
        delete dummyEvent;
      }
    }
  }
}


void EUTelNativeReader::processEUDRBDataEvent( eudaq::EUDRBEvent * eudrbEvent, EUTelEventImpl * eutelEvent ) {
  // here we have to process the EUDRB event
  //
  // prepare a collection for the raw data and one for the zs
  auto_ptr< LCCollectionVec > rawDataCollection ( new LCCollectionVec (LCIO::TRACKERRAWDATA) ) ;
  auto_ptr< LCCollectionVec > zsDataCollection  ( new LCCollectionVec (LCIO::TRACKERDATA) ) ;

  // set the proper cell encoder
  CellIDEncoder< TrackerRawDataImpl > rawDataEncoder ( EUTELESCOPE::MATRIXDEFAULTENCODING, rawDataCollection.get() );
  CellIDEncoder< TrackerDataImpl    > zsDataEncoder  ( EUTELESCOPE::ZSDATADEFAULTENCODING, zsDataCollection.get()  );

  // we can now loop over the boards contained into this EUDRB event
  for ( size_t iPlane = 0; iPlane < eudrbEvent->NumBoards(); ++iPlane ) {

    // I need to have a EUDRBBoard 
    eudaq::EUDRBBoard& eudrbBoard = eudrbEvent->GetBoard( iPlane );

    // to understand if we have problem with de-synchronisation, let
    // me prepare a Boolean switch
    bool outOfSynchFlag = false;

    // get from the eudrb detectors the current one
    EUTelBaseDetector * currentDetector = _eudrbDetectors.at( iPlane );

    // from now on we have to proceed in a different way depending if
    // the sensor was readout in RAW mode or in ZS
    string currentMode = currentDetector->getMode();
    if (  ( currentMode == "RAW" ) || ( currentMode == "RAW2" ) || ( currentMode == "RAW3" ) ) {
      rawDataEncoder["xMin"]     = currentDetector->getXMin();
      rawDataEncoder["xMax"]     = currentDetector->getXMax();
      rawDataEncoder["yMin"]     = currentDetector->getYMin();
      rawDataEncoder["yMax"]     = currentDetector->getYMax();
      rawDataEncoder["sensorID"] = iPlane;

      // put a try/catch box here mainly for the EUDRBDecoder::GetArrays
      try { 
        
        // get arrays may throw a eudaq::Exception in case the array
        // sizes are different. In this case we should catch the
        // exception, inform the user and then skip the event
        eudaq::EUDRBDecoder::arrays_t<short, short > array = _eudrbDecoder->GetArrays<short, short > ( eudrbBoard );

      } catch ( eudaq::Exception& e) {
        streamlog_out( ERROR ) << e.what() << endl << "Skipping the current event " << endl;
      }


    } else if ( ( currentMode == "ZS" ) ) {


    } else {
      // unrecognised detector mode

    }
  }


}

void EUTelNativeReader::processTLUDataEvent( eudaq::TLUEvent * tluEvent, EUTelEventImpl * eutelEvent ) {


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
      noOfTLUDetectors++;
    }

/*
  Add here your own data producer following the examples above
*/

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

}

void EUTelNativeReader::end () {
  message<MESSAGE> ("Successfully finished") ;
}

#endif


//  LocalWords:  serialiazer EUTelNativeReader MIMOTEL rawdata eudrb zsdata
//  LocalWords:  EUBRDRawModeOutputCollection EUDRBZSModeOutputCollection xMin
//  LocalWords:  EUDRBSparsePixelType sensorID xMax yMin yMax EUDRBBoard
//  LocalWords:  EUDRBDecoder
