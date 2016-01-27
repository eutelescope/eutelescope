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
#include "EUTelCopyPedestalProcessor.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDEncoder.h>

// system includes <>
#include <memory>

using namespace std;
using namespace marlin;
using namespace eutelescope;



EUTelCopyPedestalProcessor::EUTelCopyPedestalProcessor () : Processor("EUTelCopyPedestalProcessor"),
                                                            _pedestalCollectionVec(NULL),
                                                            _noiseCollectionVec(NULL),
                                                            _statusCollectionVec(NULL) {

  // modify processor description
  _description =
    "EUTelCopyPedestalProcessor copies the condition data into local writable collections";

  // first of all we need to register the input collection
  registerInputCollection (LCIO::TRACKERDATA, "PedestalConditionName",
                           "Pedestal input condition",
                           _pedestalConditionName, string ("pedestalDB"));

  registerInputCollection (LCIO::TRACKERDATA, "NoiseConditionName",
                           "Noise input condition",
                           _noiseConditionName, string ("noiseDB"));

  registerInputCollection (LCIO::TRACKERRAWDATA, "StatusConditionName",
                           "Status input condition",
                           _statusConditionName, string ("statusDB"));

  registerOutputCollection (LCIO::TRACKERDATA, "PedestalCollectionName",
                            "Pedestal local collection",
                            _pedestalCollectionName, string ("pedestal"));

  registerOutputCollection (LCIO::TRACKERDATA, "NoiseCollectionName",
                            "Noise local collection",
                            _noiseCollectionName, string("noise"));

  registerOutputCollection (LCIO::TRACKERRAWDATA, "StatusCollectionName",
                            "Pixel status collection",
                            _statusCollectionName, string("status"));
}


void EUTelCopyPedestalProcessor::init () {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

}

void EUTelCopyPedestalProcessor::processRunHeader (LCRunHeader * rdr) {
  std::unique_ptr<EUTelRunHeaderImpl> runHeader = std::make_unique<EUTelRunHeaderImpl>(rdr);
  runHeader->addProcessor( type() );
  ++_iRun;
}


void EUTelCopyPedestalProcessor::processEvent (LCEvent * event) {

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);

  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return;
  } else if ( evt->getEventType() == kUNKNOWN ) {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }


  if (isFirstEvent()) {

    LCCollectionVec * pedestalCollectionVecDB = dynamic_cast< LCCollectionVec * > (evt->getCollection(_pedestalConditionName));
    LCCollectionVec * noiseCollectionVecDB    = dynamic_cast< LCCollectionVec * > (evt->getCollection(_noiseConditionName));
    LCCollectionVec * statusCollectionVecDB   = dynamic_cast< LCCollectionVec * > (evt->getCollection(_statusConditionName));

    _pedestalCollectionVec = new LCCollectionVec(LCIO::TRACKERDATA);
    _noiseCollectionVec    = new LCCollectionVec(LCIO::TRACKERDATA);
    _statusCollectionVec   = new LCCollectionVec(LCIO::TRACKERRAWDATA);
    CellIDEncoder<TrackerDataImpl>    pedeDecoder  (EUTELESCOPE::MATRIXDEFAULTENCODING,_pedestalCollectionVec);
    CellIDEncoder<TrackerDataImpl>    noiseDecoder (EUTELESCOPE::MATRIXDEFAULTENCODING,_noiseCollectionVec);
    CellIDEncoder<TrackerRawDataImpl> statusDecoder(EUTELESCOPE::MATRIXDEFAULTENCODING,_statusCollectionVec);

    for (int iDetector = 0; iDetector < pedestalCollectionVecDB->getNumberOfElements(); iDetector++) {

      TrackerRawDataImpl * statusDB = dynamic_cast< TrackerRawDataImpl * > (statusCollectionVecDB->getElementAt(iDetector));
      TrackerRawDataImpl * status   = new TrackerRawDataImpl;
      status->setADCValues(statusDB->getADCValues());
      status->setCellID0(statusDB->getCellID0());
      status->setCellID1(statusDB->getCellID1());
      _statusCollectionVec->push_back(status);


      TrackerDataImpl * pedestalDB  = dynamic_cast< TrackerDataImpl * > (pedestalCollectionVecDB->getElementAt(iDetector));
      TrackerDataImpl * pedestal    = new TrackerDataImpl;
      pedestal->setChargeValues(pedestalDB->getChargeValues());
      pedestal->setCellID0(pedestalDB->getCellID0());
      pedestal->setCellID1(pedestalDB->getCellID1());
      _pedestalCollectionVec->push_back(pedestal);


      TrackerDataImpl * noiseDB  = dynamic_cast< TrackerDataImpl * > (noiseCollectionVecDB->getElementAt(iDetector));
      TrackerDataImpl * noise    = new TrackerDataImpl;
      noise->setChargeValues(noiseDB->getChargeValues());
      noise->setCellID0(noiseDB->getCellID0());
      noise->setCellID1(noiseDB->getCellID1());
      _noiseCollectionVec->push_back(noise);

    }

    _isFirstEvent = false;
  }

  evt->addCollection(_pedestalCollectionVec, _pedestalCollectionName);
  evt->takeCollection(_pedestalCollectionName);
  evt->addCollection(_noiseCollectionVec, _noiseCollectionName);
  evt->takeCollection(_noiseCollectionName);
  evt->addCollection(_statusCollectionVec, _statusCollectionName);
  evt->takeCollection(_statusCollectionName);

  ++_iEvt;

}




void EUTelCopyPedestalProcessor::end() {

  delete _pedestalCollectionVec;
  delete _noiseCollectionVec;
  delete _statusCollectionVec;

  streamlog_out ( MESSAGE2 )  << "Successfully finished" << endl;
}

