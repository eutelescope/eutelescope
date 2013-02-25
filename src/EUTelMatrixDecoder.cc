// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
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
#include "EUTELESCOPE.h"
#include "EUTelMatrixDecoder.h"
#include "EUTelExceptions.h"

#ifdef USE_GEAR
#include "gear/SiPlanesLayerLayout.h"
#endif

// lcio includes <.h>
#include <UTIL/CellIDDecoder.h>

// system includes <>
#include <string>
#include <iomanip>
#include <iostream>

using namespace std;
using namespace lcio;
using namespace eutelescope;

#ifdef USE_GEAR
using namespace gear;
#endif

EUTelMatrixDecoder::EUTelMatrixDecoder(int xNoOfPixel, int yNoOfPixel) throw(InvalidParameterException)
: _xNoOfPixel(xNoOfPixel),
  _yNoOfPixel(yNoOfPixel),
  _xMin(0),
  _yMin(0){

  if ( xNoOfPixel <= 0 ) throw InvalidParameterException("xNoOfPixel has to be positive");
  if ( yNoOfPixel <= 0 ) throw InvalidParameterException("yNoOfPixel has to be positive");
  _xNoOfPixel = xNoOfPixel;
  _yNoOfPixel = yNoOfPixel;
  _xMin = 0;
  _yMin = 0;

}

EUTelMatrixDecoder::EUTelMatrixDecoder(int xNoOfPixel, int yNoOfPixel, int xMin, int yMin)
: _xNoOfPixel(xNoOfPixel),
  _yNoOfPixel(yNoOfPixel),
  _xMin(xMin),
  _yMin(yMin){
  if ( xNoOfPixel <= 0 ) throw InvalidParameterException("xNoOfPixel has to be positive");
  if ( yNoOfPixel <= 0 ) throw InvalidParameterException("yNoOfPixel has to be positive");
  _xNoOfPixel = xNoOfPixel;
  _yNoOfPixel = yNoOfPixel;
  _xMin = xMin;
  _yMin = yMin;
}


#ifdef USE_GEAR
EUTelMatrixDecoder::EUTelMatrixDecoder(gear::SiPlanesLayerLayout * siPlanes, int layerIndex)
: _xNoOfPixel(0),
  _yNoOfPixel(0),
  _xMin(0),
  _yMin(0)
{
  _xNoOfPixel = siPlanes->getSensitiveNpixelX(layerIndex);
  _yNoOfPixel = siPlanes->getSensitiveNpixelY(layerIndex);
  _xMin = 0;
  _yMin = 0;

}
#endif

int EUTelMatrixDecoder::getIndexFromXY(int x, int y) const {
  
  int xCor = x - _xMin;
  int yCor = y - _yMin;
  return xCor + yCor * _xNoOfPixel;

}

void EUTelMatrixDecoder::getXYFromIndex(int index, int& x, int& y) const {
  
  y = getYFromIndex(index);
  x = getXFromIndex(index);

}

int EUTelMatrixDecoder::getXFromIndex(int index) const {
  
  return ( index % _xNoOfPixel ) + _xMin;

}

int EUTelMatrixDecoder::getYFromIndex(int index) const {
  
  return ( index / _xNoOfPixel ) + _yMin;

}

std::ostream& eutelescope::operator<< (std::ostream& os, const EUTelMatrixDecoder& decoder){
        
  int spacer = 30;
  
  os << setiosflags(ios::left) 
     << std::setw(spacer) << "Number of pixel along X " << decoder._xNoOfPixel << "\n"
     << std::setw(spacer) << "Number of pixel along Y " << decoder._yNoOfPixel << "\n"
     << std::setw(spacer) << "Origin  " << "(" << decoder._xMin << ", " << decoder._yMin  << ")"
     << resetiosflags(ios::left);
  
  return os;
}
