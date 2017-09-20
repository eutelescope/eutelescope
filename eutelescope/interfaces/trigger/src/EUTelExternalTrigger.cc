/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// personal includes ".h"
#include "EUTelExternalTrigger.h"

using namespace eutelescope;

EUTelExternalTrigger::EUTelExternalTrigger(long long unsigned timestamp, short label){ }

long long unsigned EUTelExternalTrigger::getTimestamp() const { 
  return _timestamp; 
}

short EUTelExternalTrigger::getLabel() const { 
  return _label; 
}

void EUTelExternalTrigger::setTimestamp(long long unsigned timestamp) { 
  _timestamp = timestamp; 
}

void EUTelExternalTrigger::setLabel(short label) { 
  _label = label; 
}

unsigned int EUTelExternalTrigger::GetNoOfElements() const { 
  return _nElement; 
}
