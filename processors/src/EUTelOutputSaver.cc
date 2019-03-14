/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelOutputSaver.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"

// marlin includes ".h"
#include "marlin/LCIOOutputProcessor.h"

// lcio includes <.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/LCTime.h>

// system includes <>
#include <memory>

using namespace marlin;
using namespace eutelescope;

EUTelOutputSaver::EUTelOutputSaver()
    : LCIOOutputProcessor("EUTelOutputSaver") {

  _description = "EUTelOutputSaver writes the current event to the specified LCIO outputfile. "
                 "Eventually it adds a EORE at the of the file if it was missing. "
                 "Needs to be the last ActiveProcessor.";

  #ifndef MARLIN_VERSION_GE
  registerProcessorParameter("LCIOOutputFile", 
			     "Name of output file",
                             _lcioOutputFile,
			     std::string("outputfile.slcio"));

  registerProcessorParameter("LCIOWriteMode",
			     "Write mode for output file:  WRITE_APPEND or WRITE_NEW",
			     _lcioWriteMode,
			     std::string("None"));

  std::vector<std::string> dropNamesExamples;
  dropNamesExamples.push_back("rawdata");
  dropNamesExamples.push_back("data");
  dropNamesExamples.push_back("pedestal");
  dropNamesExamples.push_back("noise");
  dropNamesExamples.push_back("status");

  registerOptionalParameter("DropCollectionNames",
                            "Drops the named collections from the event",
                            _dropCollectionNames,
			    dropNamesExamples);

  std::vector<std::string> dropTypesExample;
  dropTypesExample.push_back("TrackerRawData");
  dropTypesExample.push_back("TrackerData");

  registerOptionalParameter("DropCollectionTypes",
			    "Drops all collections of the given type from the event",
			    _dropCollectionTypes,
			    dropTypesExample);

  registerOptionalParameter("SplitFileSizekB",
                            "Will split output file if size in kB exceeds "
                            "given value - doesn't work with APPEND and NEW",
                            _splitFileSizekB,
                            1992294); // 1.9 GB in kB
  #endif

  registerProcessorParameter("SkipIntermediateEORE",
			     "Set it to true to remove intermediate EORE in merged runs",
			     _skipIntermediateEORESwitch,
			     true);
}

void EUTelOutputSaver::init() {

  //needs to be reimplemented since it is virtual in LCIOOutputProcessor
  LCIOOutputProcessor::init();
}

void EUTelOutputSaver::processRunHeader(LCRunHeader *run) {

  std::unique_ptr<EUTelRunHeaderImpl> runHeader =
      std::make_unique<EUTelRunHeaderImpl>(run);
  runHeader->addProcessor(type());
  LCIOOutputProcessor::processRunHeader(run);
}

void EUTelOutputSaver::processEvent(LCEvent *evt) {

  EUTelEventImpl *eutelEvt = static_cast<EUTelEventImpl *>(evt);

  if(_skipIntermediateEORESwitch && (eutelEvt->getEventType() == kEORE)) {
    // ok the user wants me to skip this because it may be an
    // intermediate EORE. But how to know this is not the last of the
    // intermediate EORE? Easy, don't care, set the _eventType to kDE
    // and in case this the last one, end() will be called soon after
    // and a kEORE will be appended!
    _eventType = kDE;
    return;
  }

  LCIOOutputProcessor::processEvent(evt);
  _eventType = eutelEvt->getEventType();
}

void EUTelOutputSaver::end() {

  if(_eventType != kEORE) {
    streamlog_out(WARNING4) << "Adding a EORE because was missing" 
			    << std::endl;
			   
    EUTelEventImpl *event = new EUTelEventImpl;
    event->setDetectorName("kEORE fix by EUTelOutputSaver");
    event->setEventType(kEORE);
    event->setEventNumber(_nEvt + 1);

    LCTime *now = new LCTime;
    event->setTimeStamp(now->timeStamp());
    delete now;

    _lcWrt->writeEvent(static_cast<LCEventImpl *>(event));
    delete event;
  }

  streamlog_out(MESSAGE5) << "Writing the output file " << _lcioOutputFile 
			  << std::endl;
  _lcWrt->close();
}
