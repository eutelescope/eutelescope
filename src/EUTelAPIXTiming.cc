// Version: $Id$
// Author:  Georg Troska <mailto: georg.troska@uni-dortmund.de>


/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
// personal includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelAPIXTiming.h"

// system includes <>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace eutelescope;


EUTelAPIXTiming::EUTelAPIXTiming() {
  _sensorID = 0;
  _tluCounter = 0;
  _tpllCounter = 0;
  _realtimeSec = 0;
  _realtimeNs = 0;
  _type = kAPIX;
}

EUTelAPIXTiming::EUTelAPIXTiming(short sensorID, uint64_t realtime, uint32_t tlu_counter, uint32_t tpll_counter ) {
  _sensorID = sensorID;
  _tluCounter = tlu_counter;
  _tpllCounter = tpll_counter;
  _realtimeSec = realtime/1000000000;
  _realtimeNs = realtime-(_realtimeSec*1000000000);
  _type = kAPIX;
  _noOfElements = 5;

}

EUTelAPIXTiming::EUTelAPIXTiming(const EUTelAPIXTiming &orig) : EUTelTiming() {
  _sensorID = orig.getSensorID();
  _tluCounter = orig.getTLUCounter();
  _tpllCounter = orig.getTPLLCounter();
  _realtimeSec = orig.getRealtime()/1000000000;
  _realtimeNs =  orig.getRealtime()-(_realtimeSec*1000000000);
  _noOfElements = orig.getNoOfElements();
  _type = orig.getDetectorType();
}


unsigned int EUTelAPIXTiming::getNoOfElements() const {
  return _noOfElements;
}

EUTelDetectorType EUTelAPIXTiming::getDetectorType() const {
  return _type;
}

void EUTelAPIXTiming::print(std::ostream& os) const {
  int bigWidth = 50;
  for ( int i = 0 ; i < bigWidth ; ++i ) {
    os << "-";
  }
  os << endl;
  int width = 20;
  os << setw(width) << setiosflags(ios::left) << "Type: "     << static_cast< int >(_type) << endl
     << setw(width) << "Elements: " << _noOfElements << endl
     << setw(width) << "Realtime (sec): "  << _realtimeSec << endl
     << setw(width) << "Realtime (ns): "  << _realtimeNs << endl
     << setw(width) << "TLU: "  << _tluCounter << endl
     << setw(width) << "TPLL: "  << _tpllCounter << endl
     << setw(width) << "Sensor: "  << _sensorID << endl;

  for ( int i = 0 ; i < bigWidth ; ++i ) {
    os << "-";
  }
}
