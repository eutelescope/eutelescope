/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


// personal includes ".h"
#include "AlibavaRunHeaderImpl.h"
#include "ALIBAVA.h"
// marlin includes ".h"

// lcio includes <.h>
#include "lcio.h"
#include "LCIOTypes.h"
#include "UTIL/LCTime.h"

using namespace std;
using namespace alibava;

void AlibavaRunHeaderImpl::setDetectorName (string name) {
   // sets detector name
	
   _lcHeader->setDetectorName(name);
}

void AlibavaRunHeaderImpl::setRunNumber (int runnum) {
   // sets the run number
	
   _lcHeader->setRunNumber(runnum);
}

void AlibavaRunHeaderImpl::setHeader (std::string aheader) {
   // sets the header implementation version number
	
   _lcHeader->parameters().setValue (ALIBAVA::HEADER, aheader);
}

void AlibavaRunHeaderImpl::setHeaderVersion (int ver) {
   // sets the header implementation version number

   _lcHeader->parameters().setValue (ALIBAVA::HEADERVERSION, ver);
}

void AlibavaRunHeaderImpl::setDataType (int type) {
   // sets the type of data saved in this file

   _lcHeader->parameters().setValue (ALIBAVA::DATATYPE, type);
}

void AlibavaRunHeaderImpl::setNoOfEvents ( int num ) {
  // sets the number of events in the file
  
  _lcHeader->parameters().setValue(ALIBAVA::NOOFEVENT, num);

}

void AlibavaRunHeaderImpl::setGeoID (int id) {
   // sets the current geometry identification number
   // use this number to query the geometry database

   _lcHeader->parameters().setValue (ALIBAVA::GEOID, id);
}

void AlibavaRunHeaderImpl::setDateTime (std::string atime) {
   // sets the date and time of data
   _lcHeader->parameters().setValue (ALIBAVA::DATETIME, atime);

}

void AlibavaRunHeaderImpl::setHeaderPedestal (EVENT::FloatVec v_ped) {
	// sets pedestals stored in run header
	
   _lcHeader->parameters().setValues (ALIBAVA::HEADERPEDESTAL, v_ped);
}

void AlibavaRunHeaderImpl::setHeaderNoise (EVENT::FloatVec v_noise) {
	// sets noise stored in run header
	
   _lcHeader->parameters().setValues (ALIBAVA::HEADERNOISE, v_noise);
}

void AlibavaRunHeaderImpl::setTiltAngle (float tiltAngle) {
	// sets tilt angle of the sensors
	
   _lcHeader->parameters().setValue (ALIBAVA::TILTANGLE, tiltAngle);
}

void AlibavaRunHeaderImpl::setSensorTemperature (float sensorTemperature) {
	// sets temperature of the sensors
	
   _lcHeader->parameters().setValue (ALIBAVA::SENSORTEMPERATURE, sensorTemperature);
}

void AlibavaRunHeaderImpl::setChipSelection (EVENT::IntVec achipselection) {
	// sets if there is chip selection
	
   _lcHeader->parameters().setValues (ALIBAVA::SELECTEDCHIPNUM, achipselection);
}


void AlibavaRunHeaderImpl::addProcessor (std::string processor) {
   // add a "processor" to the list of applied processors
	
   lcio::StringVec processorVec;
   _lcHeader->parameters().getStringVals (ALIBAVA::APPLIEDPROCESSOR, processorVec);
   processorVec.push_back (processor);
   _lcHeader->parameters().setValues (ALIBAVA::APPLIEDPROCESSOR, processorVec);
}


