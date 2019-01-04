/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelEventNumberPrinter.h"

// marlin includes ".h"
#ifdef MARLIN_USE_AIDA
#include <AIDA/AIDA.h>
#include <marlin/AIDAProcessor.h>
#endif

// lcio includes <.h>
#include "UTIL/LCTOOLS.h"

// system includes <>
#include <iomanip>
#include <iostream>

using namespace lcio;
using namespace marlin;
using namespace eutelescope;

EUTelEventNumberPrinter aPrintEventNumber;

EUTelEventNumberPrinter::EUTelEventNumberPrinter()
    : Processor("EUTelEventNumberPrinter"), _everyNEvents(1000),
      _printTimestamp(false), totalevents(0), totalruns(0) {
  /* the constructor call first the super constructor which the name
   * of the processor
   * this must by the same as the class name.
   */

  _description = "EUTelEventNumberPrinter prints the event number to screen"
                 " depending on the verbosity level";

  /* register steering parameters: name, description, class-variable, default
   * value
   * the type will by definded by the type of the default value
   * string, double, float and int are possible
   */
  registerProcessorParameter("EveryNEvents",
                             "Print event number for every n-th event",
                             _everyNEvents, 
                             1000);
  registerOptionalParameter("printTimestamp",
                            "print the event timestamp as read from LCIO",
                            _printTimestamp, 
                            false);
}

void EUTelEventNumberPrinter::init() {
  // this method is called only once even when the rewind is active

  // usually a good idea to do
  printParameters();
}

void EUTelEventNumberPrinter::processRunHeader(LCRunHeader *run) {

  run->parameters().setValue(_processorName + "_revision", "$Rev: 699 $");

  for(ProcParamMap::iterator it = _map.begin(); it != _map.end(); it++) {
    if(!it->second->isOptional() || it->second->valueSet()) {
      run->parameters().setValue(_processorName + "_" + it->second->name(),
                                 it->second->value());
    }	
  }

  //increment the total runs counter
  totalruns++;
}

void EUTelEventNumberPrinter::processEvent(LCEvent *evt) {

  //regular output in case verbosity MESSAGE is set:
  if(evt->getEventNumber() <= 10 ||
     (evt->getEventNumber() <= 100 && evt->getEventNumber() % 10 == 0) ||
      evt->getEventNumber() % _everyNEvents == 0) {
    streamlog_out(MESSAGE5) << "Processing event " << std::setw(7)
                            << evt->getEventNumber() << " in run "
                            << std::setw(6) << std::setfill('0')
                            << evt->getRunNumber() << std::setfill(' ');
    if(totalruns > 1)
      streamlog_out(MESSAGE5) << " (" << totalevents << " events in total)";
    if(_printTimestamp)
      streamlog_out(MESSAGE5) << ", timestamp " << evt->getTimeStamp();
    streamlog_out(MESSAGE5) << std::endl;
  }
  //additional output if we have the DEBUG level: print every event.
  else {
    streamlog_out(DEBUG5) << "Processing event " << std::setw(7)
                          << evt->getEventNumber() << " in run "
                          << evt->getRunNumber();
    if(_printTimestamp)
      streamlog_out(DEBUG5) << ", timestamp " << evt->getTimeStamp();
    streamlog_out(DEBUG5) << std::endl;
  }

  //increment total event counter
  totalevents++;
}

void EUTelEventNumberPrinter::check(LCEvent * /* evt */) {
  /* Nothing to do here... */
}

void EUTelEventNumberPrinter::end() {
  streamlog_out(MESSAGE4) << "Finished. Processed " << totalevents
                          << " events in " << totalruns << " runs in total."
                          << std::endl;
}
