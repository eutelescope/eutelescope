/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with contact names in all development based on this file.
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

/**Defaut constructor, we use 1 int to store the sensorID and 12 doubles for the 3 offsets and 3 
 * angles and their corresponding uncertainties */ 
EUTelAlignmentConstant::EUTelAlignmentConstant() : IMPL::LCGenericObjectImpl(1, 0, 12)
{
  _typeName        = "Alignment constant";
  _dataDescription = "sensorID, xOff, yOff, zOff, alpha, beta, gamma, xErr, yErr, zErr, alphaErr, betaErr, gammaErr";
  _isFixedSize     = true;
}

EUTelAlignmentConstant::EUTelAlignmentConstant( int sensorID,
                                                double xOff,   double yOff,   double zOff,
                                                double alpha, double beta, double gamma,
                                                double xOffErr,double yOffErr,double zOffErr,
                                                double alphaErr, double betaErr, double gammaErr )  :

  IMPL::LCGenericObjectImpl(1, 0, 12)
{
  _typeName        = "Alignment constant";
  _dataDescription = "sensorID, xOff, yOff, zOff, alpha, beta, gamma, xErr, yErr, zErr, alphaErr, betaErr, gammaErr";
  _isFixedSize     = true;

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

  // set the sensor shifts in mm
  setDoubleVal( 6, xOffErr );
  setDoubleVal( 7, yOffErr );
  setDoubleVal( 8, zOffErr );

  // set the sensor angles
  setDoubleVal( 9, alphaErr );
  setDoubleVal( 10, betaErr );
  setDoubleVal( 11, gammaErr );

}

void EUTelAlignmentConstant::setSensorID( int id ) { setIntVal( 0 , id ) ; }

void EUTelAlignmentConstant::setXOffset( double off ) { setDoubleVal( 0, off ); }
void EUTelAlignmentConstant::setYOffset( double off ) { setDoubleVal( 1, off ); }
void EUTelAlignmentConstant::setZOffset( double off ) { setDoubleVal( 2, off ); }

void EUTelAlignmentConstant::setAlpha( double theta ) { setDoubleVal( 3, theta ); }
void EUTelAlignmentConstant::setBeta( double theta ) { setDoubleVal( 4, theta ); }
void EUTelAlignmentConstant::setGamma( double theta ) { setDoubleVal( 5, theta ); }

void EUTelAlignmentConstant::setXOffsetError( double err ) { setDoubleVal( 6, err ); }
void EUTelAlignmentConstant::setYOffsetError( double err ) { setDoubleVal( 7, err ); }
void EUTelAlignmentConstant::setZOffsetError( double err ) { setDoubleVal( 8, err ); }

void EUTelAlignmentConstant::setAlphaError( double err ) { setDoubleVal( 9, err ); }
void EUTelAlignmentConstant::setBetaError( double err ) { setDoubleVal( 10, err ); }
void EUTelAlignmentConstant::setGammaError( double err ) { setDoubleVal( 11, err ); }

int EUTelAlignmentConstant::getSensorID() const { return getIntVal( 0 ) ; }

double EUTelAlignmentConstant::getXOffset() const { return getDoubleVal( 0 ) ; }
double EUTelAlignmentConstant::getYOffset() const { return getDoubleVal( 1 ) ; }
double EUTelAlignmentConstant::getZOffset() const { return getDoubleVal( 2 ) ; }

double EUTelAlignmentConstant::getAlpha() const { return getDoubleVal( 3 ) ; }
double EUTelAlignmentConstant::getBeta() const { return getDoubleVal( 4 ) ; }
double EUTelAlignmentConstant::getGamma() const { return getDoubleVal( 5 ) ; }

double EUTelAlignmentConstant::getXOffsetError() const { return getDoubleVal( 6 ) ; }
double EUTelAlignmentConstant::getYOffsetError() const { return getDoubleVal( 7 ) ; }
double EUTelAlignmentConstant::getZOffsetError() const { return getDoubleVal( 8 ) ; }

double EUTelAlignmentConstant::getAlphaError() const { return getDoubleVal( 9 ) ; }
double EUTelAlignmentConstant::getBetaError() const { return getDoubleVal( 10 ) ; }
double EUTelAlignmentConstant::getGammaError() const { return getDoubleVal( 11 ) ; }

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
     << setw(largeWidth) << " ZY angle "
     << setw(largeWidth) << " ZX angle "
     << setw(largeWidth) << " XY angle "
     << endl
     << setw(largeWidth) << getXOffset()
     << setw(largeWidth) << getYOffset()
     << setw(largeWidth) << getZOffset()
     << setw(largeWidth) << getAlpha()
     << setw(largeWidth) << getBeta()
     << setw(largeWidth) << getGamma()
     << endl
     << setw(largeWidth) << getXOffsetError()
     << setw(largeWidth) << getYOffsetError()
     << setw(largeWidth) << getZOffsetError()
     << setw(largeWidth) << getAlphaError()
     << setw(largeWidth) << getBetaError()
     << setw(largeWidth) << getGammaError()
     << endl
     << dashline
     << endl;
  return;
}
