// Author Igor Rubinskiy, INFN <mailto:rubinky@mail.desy.de>
// Version $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelReferenceHit.h"

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/LCGenericObjectImpl.h>

// system includes <>
#include <iostream>
#include <iomanip>

using namespace lcio;
using namespace eutelescope;
using namespace std;

EUTelReferenceHit::EUTelReferenceHit() : IMPL::LCGenericObjectImpl( 1, 0 , REFHIT_CONST_MAX_SIZE ) {
  _typeName        = "Reference hit";
  _dataDescription = "sensorID, xOff, yOff, zOff, alpha, beta, gamma ";
  _isFixedSize     = true;


  setIntVal   ( 0, 0   );
  for ( size_t iPos = 0; iPos < REFHIT_CONST_MAX_SIZE; ++iPos ) {
   setDoubleVal( iPos, 0.0 );
  }


}

EUTelReferenceHit::EUTelReferenceHit( int sensorID,
                                                double xOff,   double yOff,   double zOff,
                                                double alpha, double beta, double gamma )  :

  IMPL::LCGenericObjectImpl( 1, 0, REFHIT_CONST_MAX_SIZE) {

  _typeName        = "Reference Hit";
  _dataDescription = "sensorID,\n"
    "xOff,    yOff,    zOff,    alpha,    beta,    gamma";
  _isFixedSize     = true;

  setIntVal   ( 0, 0   );
  for ( size_t iPos = 0; iPos < REFHIT_CONST_MAX_SIZE; ++iPos ) {
   setDoubleVal( iPos, 0.0 );
  }

  // set the sensor ID
  setIntVal( 0 , sensorID );

  // set the sensor shifts in mm
  setDoubleVal( 0, xOff );
  setDoubleVal( 1, yOff );
  setDoubleVal( 2, zOff );

  // set the sensor angles
  setDoubleVal( 3, alpha );
  setDoubleVal( 4, beta );
  setDoubleVal( 5, gamma );
}

void EUTelReferenceHit::setSensorID( int id ) { setIntVal( 0 , id ) ; }

void EUTelReferenceHit::setXOffset( double off ) { setDoubleVal( 0, off ); }
void EUTelReferenceHit::setYOffset( double off ) { setDoubleVal( 1, off ); }
void EUTelReferenceHit::setZOffset( double off ) { setDoubleVal( 2, off ); }

void EUTelReferenceHit::setAlpha( double theta ) { setDoubleVal( 3, theta ); }
void EUTelReferenceHit::setBeta( double theta ) { setDoubleVal( 4, theta ); }
void EUTelReferenceHit::setGamma( double theta ) { setDoubleVal( 5, theta ); }


int EUTelReferenceHit::getSensorID()   const { return getIntVal( 0 ) ; }

double EUTelReferenceHit::getXOffset() const { return getDoubleVal( 0 ) ; }
double EUTelReferenceHit::getYOffset() const { return getDoubleVal( 1 ) ; }
double EUTelReferenceHit::getZOffset() const { return getDoubleVal( 2 ) ; }

double EUTelReferenceHit::getAlpha() const { return getDoubleVal( 3 ) ; }
double EUTelReferenceHit::getBeta()  const { return getDoubleVal( 4 ) ; }
double EUTelReferenceHit::getGamma() const { return getDoubleVal( 5 ) ; }


void EUTelReferenceHit::print(ostream & os ) const {

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
     << setw(largeWidth) << getAlpha()
     << setw(largeWidth) << getBeta()
     << setw(largeWidth) << getGamma()
     << endl
     << dashline
     << endl;


  return;
}
