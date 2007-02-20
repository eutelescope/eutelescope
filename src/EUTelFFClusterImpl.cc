// Author:  Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version: $Id: EUTelFFClusterImpl.cc,v 1.1 2007-02-20 11:35:50 bulgheroni Exp $

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// personal includes ".h"
#include "EUTelFFClusterImpl.h"

// lcio includes <.h>
#include <LCIOTypes.h>
#include <IMPL/TrackerDataImpl.h>

using namespace eutelescope;

EUTelFFClusterImpl::EUTelFFClusterImpl() : IMPL::TrackerDataImpl() { ; } 


float EUTelFFClusterImpl::getTotalCharge() const {
  
  float totalCharge = 0;
  
  FloatVec::const_iterator iter = getChargeValues().begin();

  while (iter != getChargeValues().end()) {
    totalCharge += (*iter++);
  }

  return totalCharge;
}

void EUTelFFClusterImpl::getCenterOfGravity(float& xCoG, float& yCoG) const {

  int xSize, ySize;
  getClusterSize(xSize, ySize);

  int xSeed, ySeed;
  getSeedCoord(xSeed, ySeed);

  float normalization = 0;
  float tempX = 0;
  float tempY = 0;

  int iPixel = 0;
  for (int yPixel = -1 * (ySize / 2); yPixel <= (ySize / 2); yPixel++) {
    for (int xPixel = -1 * (xSize / 2); xPixel <= (xSize / 2); xPixel++) {
      normalization += getChargeValues()[iPixel];
      tempX += xPixel * getChargeValues()[iPixel];
      tempY += yPixel * getChargeValues()[iPixel];
      ++iPixel;
    }
  }

  if ( normalization != 0)  {
    xCoG = tempX / normalization + xSeed;
    yCoG = tempY / normalization + ySeed;
  } else {
    xCoG = 0;
    yCoG = 0;
  }
  
}
  
  
