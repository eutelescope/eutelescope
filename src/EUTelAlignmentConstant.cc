// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelAlignmentConstant.cc,v 1.1 2008-07-09 13:01:26 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h" 
#include "EUTelAlignmentConstant.h"

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/LCGenericObjectImpl.h>

// system includes <>
#include <iostream>
#include <iomanip>

using namespace lcio;
using namespace eutelescope;
using namespace std;

EUTelAlignmentConstant::EUTelAlignmentConstant() : IMPL::LCGenericObjectImpl( ALIGN_CONST_MAX_SIZE, 0 , ALIGN_CONST_MAX_SIZE ) {
  _typeName        = "Alignment constant";
  _dataDescription = "sensorID, xOff, yOff, zOff, xTheta, yTheta, zTheta + 13 spare fields";
  _isFixedSize     = true;


  for ( size_t iPos = 0; iPos < ALIGN_CONST_MAX_SIZE; ++iPos ) {
    setIntVal   ( iPos, 0   );
    setDoubleVal( iPos, 0.0 );
  }


}

EUTelAlignmentConstant::EUTelAlignmentConstant( int sensorID,
					       double xOff,   double yOff,   double zOff,
					       double xTheta, double yTheta, double zTheta,
					       double xOffErr,double yOffErr,double zOffErr,
					       double xThetaErr, double yThetaErr, double zThetaErr )  : 
  
  IMPL::LCGenericObjectImpl(ALIGN_CONST_MAX_SIZE,0,ALIGN_CONST_MAX_SIZE) {
  
  _typeName        = "Alignment constant";
  _dataDescription = "sensorID,\n"
    "xOff,    yOff,    zOff,    xTheta,    yTheta,    zTheta\n"
    "xOffErr, yOffErr, zOffErr, xThetaErr, yThetaErr, zThetaErr";
  _isFixedSize     = true;
  
  for ( size_t iPos = 0; iPos < ALIGN_CONST_MAX_SIZE; ++iPos ) {
    setIntVal   ( iPos, 0   );
    setDoubleVal( iPos, 0.0 );
  }

  // set the sensor ID
  setIntVal( 0 , sensorID );

  // set the sensor shifts in mm
  setDoubleVal( 0, xOff );
  setDoubleVal( 1, yOff );
  setDoubleVal( 2, zOff );
  
  // set the sensor angles
  setDoubleVal( 3, xTheta );
  setDoubleVal( 4, yTheta );
  setDoubleVal( 5, zTheta );

  // set the sensor shifts in mm
  setDoubleVal( 6, xOffErr );
  setDoubleVal( 7, yOffErr );
  setDoubleVal( 8, zOffErr );
  
  // set the sensor angles
  setDoubleVal( 9, xThetaErr );
  setDoubleVal( 10, yThetaErr );
  setDoubleVal( 11, zThetaErr );

}

void EUTelAlignmentConstant::setSensorID( int id ) { setIntVal( 0 , id ) ; }

void EUTelAlignmentConstant::setXOffset( double off ) { setDoubleVal( 0, off ); }
void EUTelAlignmentConstant::setYOffset( double off ) { setDoubleVal( 1, off ); }
void EUTelAlignmentConstant::setZOffset( double off ) { setDoubleVal( 2, off ); }

void EUTelAlignmentConstant::setXTheta( double theta ) { setDoubleVal( 3, theta ); }
void EUTelAlignmentConstant::setYTheta( double theta ) { setDoubleVal( 4, theta ); }
void EUTelAlignmentConstant::setZTheta( double theta ) { setDoubleVal( 5, theta ); }

void EUTelAlignmentConstant::setXOffsetError( double err ) { setDoubleVal( 6, err ); }
void EUTelAlignmentConstant::setYOffsetError( double err ) { setDoubleVal( 7, err ); }
void EUTelAlignmentConstant::setZOffsetError( double err ) { setDoubleVal( 8, err ); }

void EUTelAlignmentConstant::setXThetaError( double err ) { setDoubleVal( 9, err ); }
void EUTelAlignmentConstant::setYThetaError( double err ) { setDoubleVal( 10, err ); }
void EUTelAlignmentConstant::setZThetaError( double err ) { setDoubleVal( 11, err ); }

int EUTelAlignmentConstant::getSensorID() const { return getIntVal( 0 ) ; }

double EUTelAlignmentConstant::getXOffset() const { return getDoubleVal( 0 ) ; }
double EUTelAlignmentConstant::getYOffset() const { return getDoubleVal( 1 ) ; }
double EUTelAlignmentConstant::getZOffset() const { return getDoubleVal( 2 ) ; }

double EUTelAlignmentConstant::getXTheta() const { return getDoubleVal( 3 ) ; }
double EUTelAlignmentConstant::getYTheta() const { return getDoubleVal( 4 ) ; }
double EUTelAlignmentConstant::getZTheta() const { return getDoubleVal( 5 ) ; }

double EUTelAlignmentConstant::getXOffsetError() const { return getDoubleVal( 6 ) ; }
double EUTelAlignmentConstant::getYOffsetError() const { return getDoubleVal( 7 ) ; }
double EUTelAlignmentConstant::getZOffsetError() const { return getDoubleVal( 8 ) ; }

double EUTelAlignmentConstant::getXThetaError() const { return getDoubleVal( 9 ) ; }
double EUTelAlignmentConstant::getYThetaError() const { return getDoubleVal( 10 ) ; }
double EUTelAlignmentConstant::getZThetaError() const { return getDoubleVal( 11 ) ; }

void EUTelAlignmentConstant::print(ostream & os ) const {
  
  const int maxFieldNo = 7;
  const int narrowWidth = 10;
  const int largeWidth = 14;
  const int lineWidth  = largeWidth * maxFieldNo;
  
  string dashline( lineWidth, '-' );

  os << dashline << endl
     << setw( lineWidth - narrowWidth ) << setiosflags(ios::left) << "Sensor ID " << resetiosflags(ios::left)
     << setw( narrowWidth ) << setiosflags(ios::right) << getSensorID() << resetiosflags(ios::right) << endl
     << dashline << endl
     << setw(largeWidth) << " X off [mm] "
     << setw(largeWidth) << " Y off [mm] "
     << setw(largeWidth) << " Z off [mm] "
     << setw(largeWidth) << " X angle "
     << setw(largeWidth) << " Y angle "
     << setw(largeWidth) << " Z angle "
     << endl
     << setw(largeWidth) << getXOffset() 
     << setw(largeWidth) << getYOffset()
     << setw(largeWidth) << getZOffset()
     << setw(largeWidth) << getXTheta()
     << setw(largeWidth) << getYTheta()
     << setw(largeWidth) << getZTheta()
     << endl
     << setw(largeWidth) << getXOffsetError() 
     << setw(largeWidth) << getYOffsetError()
     << setw(largeWidth) << getZOffsetError()
     << setw(largeWidth) << getXThetaError()
     << setw(largeWidth) << getYThetaError()
     << setw(largeWidth) << getZThetaError()
     << endl
     << dashline
     << endl;


  return;
}
