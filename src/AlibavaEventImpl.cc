/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


// personal includes ".h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"

// marlin includes ".h"

// lcio includes <.h>

using namespace std;
using namespace alibava;


AlibavaEventImpl::AlibavaEventImpl() : IMPL::LCEventImpl() { /* NO - OP */ ; }

void AlibavaEventImpl::setEventType(EventType type) {

  int iType = static_cast<EventType>(type);
  setEventType(iType);

}

void AlibavaEventImpl::setEventType(int type) {

  _params.setValue(ALIBAVA::EVENTTYPE, type);

}

void AlibavaEventImpl::setEventValue(float avalue) {
	
	_params.setValue(ALIBAVA::EVENTVALUE, avalue);
	
}

void AlibavaEventImpl::setEventSize(int asize) {
	
	_params.setValue(ALIBAVA::EVENTSIZE, asize);
	
}

void AlibavaEventImpl::setEventCode(int acode) {
	
	_params.setValue(ALIBAVA::EVENTCODE, acode);
	
}

void AlibavaEventImpl::setEventTime(float atime) {
	
	_params.setValue(ALIBAVA::EVENTTIME, atime);
	
}

void AlibavaEventImpl::setEventTemp(float atemp) {
	
	_params.setValue(ALIBAVA::EVENTTEMP, atemp);
	
}

void AlibavaEventImpl::setCalCharge(float acharge) {
	
	_params.setValue(ALIBAVA::CALCHARGE, acharge);
	
}
void AlibavaEventImpl::setCalDelay(float adelay) {
	
	_params.setValue(ALIBAVA::CALDELAY, adelay);
	
}
void AlibavaEventImpl::maskEvent() {
	_params.setValue(ALIBAVA::EVENTMASK, 1);
	
}
void AlibavaEventImpl::unmaskEvent() {
	
	_params.setValue(ALIBAVA::EVENTMASK, 0);
	
}





