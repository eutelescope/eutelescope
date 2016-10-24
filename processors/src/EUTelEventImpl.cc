// Author:  Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version: $Id$

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// personal includes ".h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"

// marlin includes ".h"

// lcio includes <.h>
// #include "lcio.h"

using namespace std;
using namespace eutelescope;


EUTelEventImpl::EUTelEventImpl() : IMPL::LCEventImpl() { /* NO - OP */ ; }

void EUTelEventImpl::setEventType(EventType type) {

  int iType = static_cast<EventType>(type);
  setEventType(iType);

}

void EUTelEventImpl::setEventType(int type) {

  _params.setValue(EUTELESCOPE::EVENTTYPE, type);

}
