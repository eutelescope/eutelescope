// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelSetupDescription.cc,v 1.1 2008-09-09 08:15:46 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelSetupDescription.h"
#include "EUTELESCOPE.h"
#include "EUTelBaseDetector.h"
#include "EUTelPixelDetector.h"
#include "EUTelExceptions.h"

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/LCGenericObjectImpl.h>

// system includes <>
#include <string>

using namespace lcio;
using namespace eutelescope;
using namespace std;

EUTelSetupDescription::EUTelSetupDescription(EUTelDetectorType type, EUTelReadoutMode mode) :
  IMPL::LCGenericObjectImpl(5,0,0) {
  _typeName        = "Setup Description";
  _dataDescription = "type:i,mode:i,spare1:i,spare2:i,spare3:i";
  _isFixedSize     = true;
  setIntVal( 0, static_cast< int > ( type ) );
  setIntVal( 1, static_cast< int > ( mode ) );
}

EUTelSetupDescription::EUTelSetupDescription() :
  IMPL::LCGenericObjectImpl(5,0,0) {
  _typeName        = "Setup Description";
  _dataDescription = "type:i,mode:i,spare1:i,spare2:i,spare3:i";
  _isFixedSize     = true;
}

EUTelSetupDescription::EUTelSetupDescription(EUTelPixelDetector * detector)  :
  IMPL::LCGenericObjectImpl(5,0,0) {
  _typeName        = "Setup Description";
  _dataDescription = "type:i,mode:i,spare1:i,spare2:i,spare3:i";
  _isFixedSize     = true;

  string typeS = detector->getDetectorName();
  EUTelDetectorType typeE;

  if ( typeS == "Mimosa18" ) typeE = kMimosa18;
  else if ( typeS == "MimoTel") typeE = kMimoTel;
  else {
    throw UnknownDataTypeException( typeS + " is not a valid detector type." );
  }

  string modeS = detector->getMode();
  EUTelReadoutMode modeE;

  if ( modeS == "RAW2" ) modeE = kRAW2;
  else if ( modeS == "RAW3" ) modeE = kRAW3;
  else {
    throw UnknownDataTypeException( modeS + " is not a valid readout mode." );
  }


  setIntVal( 0, static_cast< int > ( typeE ) );
  setIntVal( 1, static_cast< int > ( modeE ) );
}

void EUTelSetupDescription::setDetectorType(EUTelDetectorType type) {
  setIntVal( 0, static_cast< int > ( type ) );
}

void EUTelSetupDescription::setReadoutMode(EUTelReadoutMode mode) {
  setIntVal( 1, static_cast< int > ( mode ) );
}

EUTelDetectorType EUTelSetupDescription::getDetectorType()  {
  return static_cast< EUTelDetectorType > ( getIntVal( 0 ) );
}

EUTelReadoutMode EUTelSetupDescription::getReadoutMode()  {
  return static_cast< EUTelReadoutMode > ( getIntVal( 1 ) );
}
