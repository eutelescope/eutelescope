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
#include "EUTelRunHeaderImpl.h"
#include "EUTELESCOPE.h"

// marlin includes ".h"

// lcio includes <.h>
#include "lcio.h"
#include "LCIOTypes.h"
#include "UTIL/LCTime.h"

using namespace std;
using namespace eutelescope;

EUTelRunHeaderImpl::EUTelRunHeaderImpl(const EUTelRunHeaderImpl &z) : _lcHeader(NULL) {
  IMPL::LCRunHeaderImpl *zlcrunheader = const_cast< EUTelRunHeaderImpl& >(z).lcRunHeader();
  _lcHeader->setRunNumber(zlcrunheader->getRunNumber());
  _lcHeader->setDetectorName(zlcrunheader->getDetectorName());
  _lcHeader->setDescription(zlcrunheader->getDescription());
  std::vector< std::string > *activesubdetectors = const_cast< std::vector< std::string >* >(zlcrunheader->getActiveSubdetectors());
  for(std::vector< std::string >::iterator it = activesubdetectors->begin(); it != activesubdetectors->end(); ++it){
    _lcHeader->addActiveSubdetector(*it);
  }
}

EUTelRunHeaderImpl& EUTelRunHeaderImpl::operator = (const EUTelRunHeaderImpl &z){
  if (this == &z) return *this;  //This handles self assignment
  IMPL::LCRunHeaderImpl *zlcrunheader = const_cast< EUTelRunHeaderImpl& >(z).lcRunHeader();
  _lcHeader->setRunNumber(zlcrunheader->getRunNumber());
  _lcHeader->setDetectorName(zlcrunheader->getDetectorName());
  _lcHeader->setDescription(zlcrunheader->getDescription());
  std::vector< std::string > *activesubdetectors = const_cast< std::vector< std::string >* >(zlcrunheader->getActiveSubdetectors());
  for(std::vector< std::string >::iterator it = activesubdetectors->begin(); it != activesubdetectors->end(); ++it){
    _lcHeader->addActiveSubdetector(*it);
  }
  return *this;
}

void EUTelRunHeaderImpl::setHeaderVersion (float ver) {
   // sets the header implementation version number

   _lcHeader->parameters().setValue (EUTELESCOPE::HEADERVERSION, ver);
}

void EUTelRunHeaderImpl::setDataType (std::string type) {
   // sets the type of data saved in this file

   _lcHeader->parameters().setValue (EUTELESCOPE::DATATYPE, type);
}

void EUTelRunHeaderImpl::setNoOfEvent ( int num ) {
  // sets the number of events in the file
  
  _lcHeader->parameters().setValue(EUTELESCOPE::NOOFEVENT, num);

}

void EUTelRunHeaderImpl::setDateTime () {
   // sets the date and time of this file to now
   // it uses the LCTime util

   UTIL::LCTime * now = new UTIL::LCTime ();
   _lcHeader->parameters().setValue (EUTELESCOPE::DATETIME, now->getDateString ());
   delete now;
}

void EUTelRunHeaderImpl::setDAQHWName (std::string name) {
   // sets the DAQ hardware name.
   // only in the case of EUTELESCOPE::DAQDATA

   _lcHeader->parameters().setValue (EUTELESCOPE::DAQHWNAME, name);
}

void EUTelRunHeaderImpl::setDAQHWVersion (float ver) {
   // sets the DAQ hardware version number
   // only in the case of EUTELESCOPE::DAQDATA

   _lcHeader->parameters().setValue (EUTELESCOPE::DAQHWVERSION, ver);
}

void EUTelRunHeaderImpl::setDAQSWName (std::string name) {
   // sets the DAQ software name
   // only in the case of EUTELESCOPE::DAQDATA

   _lcHeader->parameters().setValue (EUTELESCOPE::DAQSWNAME, name);
}

void EUTelRunHeaderImpl::setDAQSWVersion (float ver) {
   // sets the DAQ software version number
   // only in the case of EUTELESCOPE::DAQDATA

   _lcHeader->parameters().setValue (EUTELESCOPE::DAQSWVERSION, ver);
}

void EUTelRunHeaderImpl::setSimulSWName (std::string name) {
   // sets the simulation program name used to generate this data
   // only in the case of EUTELESCOPE::SIMULDATA

   _lcHeader->parameters().setValue (EUTELESCOPE::SIMULSWNAME, name);
}

void EUTelRunHeaderImpl::setSimulSWVersion (float ver) {
   // sets the simulation program version number
   // only in the case of EUTELESCOPE::SIMULDATA

   _lcHeader->parameters().setValue (EUTELESCOPE::SIMULSWVERSION, ver);
}

void EUTelRunHeaderImpl::setGeoID (int id) {
   // sets the current geometry identification number
   // use this number to query the geometry database

   _lcHeader->parameters().setValue (EUTELESCOPE::GEOID, id);
}

void EUTelRunHeaderImpl::setBeamLocation(std::string location) {
  _lcHeader->parameters().setValue(EUTELESCOPE::BEAMLOCATION, location);
}


void EUTelRunHeaderImpl::setBeamType(std::string type) {
  _lcHeader->parameters().setValue(EUTELESCOPE::BEAMTYPE, type);
}

void EUTelRunHeaderImpl::setBeamEnergy(float energy) {
  _lcHeader->parameters().setValue(EUTELESCOPE::BEAMENERGY, energy);
}

void EUTelRunHeaderImpl::setNoOfDetector (int num) {
   // sets the number of detectors in this run.  this corresponds to
   // the number of TrackerRawData objects in each LCCollectionVec

   _lcHeader->parameters().setValue (EUTELESCOPE::NOOFDETECTOR, num);
}

void EUTelRunHeaderImpl::setMinX (lcio::IntVec xMin) {
   // sets the vector containing the minimum pixel number along X

   _lcHeader->parameters().setValues (EUTELESCOPE::MINX, xMin);
}

void EUTelRunHeaderImpl::setMaxX (lcio::IntVec xMax) {
   // sets the vector containing the maximum pixel number along X

   _lcHeader->parameters().setValues (EUTELESCOPE::MAXX, xMax);
}

void EUTelRunHeaderImpl::setMinY (lcio::IntVec yMin) {
   // sets the vector containing the minimum pixel number along Y

   _lcHeader->parameters().setValues (EUTELESCOPE::MINY, yMin);
}

void EUTelRunHeaderImpl::setMaxY (lcio::IntVec yMax) {
   // sets the vector containing the maximum pixel number along Y

   _lcHeader->parameters().setValues (EUTELESCOPE::MAXY, yMax);
}

void EUTelRunHeaderImpl::addProcessor (std::string processor) {
   // add a "processor" to the list of applied processors

   lcio::StringVec processorVec;
   _lcHeader->parameters().getStringVals (EUTELESCOPE::APPLIEDPROCESSOR, processorVec);
   processorVec.push_back (processor);
   _lcHeader->parameters().setValues (EUTELESCOPE::APPLIEDPROCESSOR, processorVec);
}

void EUTelRunHeaderImpl::addIntermediateFile (std::string file) {
   // add a "file" to the list of intermediate file

   lcio::StringVec fileVec;
   _lcHeader->parameters().getStringVals (EUTELESCOPE::INTERMEDIATEFILE, fileVec);
   fileVec.push_back (file);
   _lcHeader->parameters().setValues (EUTELESCOPE::INTERMEDIATEFILE, fileVec);
}

void EUTelRunHeaderImpl::setUserComment(std::string note) {
  _lcHeader->parameters().setValue(EUTELESCOPE::USERCOMMENT, note);
}

void EUTelRunHeaderImpl::setEUDRBMode(std::string mode) {
  _lcHeader->parameters().setValue(EUTELESCOPE::EUDRBMODE, mode);
}

void EUTelRunHeaderImpl::setEUDRBDet(std::string det) {
  _lcHeader->parameters().setValue(EUTELESCOPE::EUDRBDET, det);
}
