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
#include "EUTelM26Timing.h"

// system includes <>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace eutelescope;


EUTelM26Timing::EUTelM26Timing() {
  _sensorID = 0;
  _tluCounter = 0;
  _type = kAPIX;
}

EUTelM26Timing::EUTelM26Timing(short sensorID, uint32_t tlu_counter ) {
  _sensorID = sensorID;
  _tluCounter = tlu_counter;
  _type = kMimosa26;
  _noOfElements = 2;

}

EUTelM26Timing::EUTelM26Timing(const EUTelM26Timing &orig) : EUTelTiming() {
  _sensorID = orig.getSensorID();
  _tluCounter = orig.getTLUCounter();
  _noOfElements = orig.getNoOfElements();
  _type = orig.getDetectorType();
}


unsigned int EUTelM26Timing::getNoOfElements() const {
  return _noOfElements;
}

EUTelDetectorType EUTelM26Timing::getDetectorType() const {
  return _type;
}

void EUTelM26Timing::print(std::ostream& os) const {
  int bigWidth = 50;
  for ( int i = 0 ; i < bigWidth ; ++i ) {
    os << "-";
  }
  os << endl;
  int width = 20;
  os << setw(width) << setiosflags(ios::left) << "Type: "     << static_cast< int >(_type) << endl
     << setw(width) << "Elements: " << _noOfElements << endl
     << setw(width) << "TLU: "  << _tluCounter << endl
     << setw(width) << "Sensor: "  << _sensorID << endl;

  for ( int i = 0 ; i < bigWidth ; ++i ) {
    os << "-";
  }
}
