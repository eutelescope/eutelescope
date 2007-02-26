// Author:  Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version: $Id: EUTelFFClusterImpl.cc,v 1.3 2007-02-26 09:25:27 bulgheroni Exp $

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

// system includes
#include <map>
#include <cmath>

using namespace eutelescope;
using namespace std;


EUTelFFClusterImpl::EUTelFFClusterImpl() : IMPL::TrackerDataImpl() { ; } 


float EUTelFFClusterImpl::getDistance(EUTelFFClusterImpl * otherCluster) const {

  int xOtherSeed, yOtherSeed;
  otherCluster->getSeedCoord(xOtherSeed, yOtherSeed);

  int xThisSeed, yThisSeed;
  this->getSeedCoord(xThisSeed, yThisSeed);

  return sqrt( pow(static_cast<double> (xThisSeed - xOtherSeed), 2) + pow(static_cast<double> (yThisSeed - yOtherSeed), 2) );

}

float EUTelFFClusterImpl::getExternalRadius() const {

  int xSize, ySize;
  getClusterSize(xSize, ySize);

  return 0.5 * sqrt( pow(static_cast<double> (xSize), 2) + pow(static_cast<double> (ySize), 2) );

}


float EUTelFFClusterImpl::getTotalCharge() const {
  
  float totalCharge = 0;
  
  FloatVec::const_iterator iter = getChargeValues().begin();

  while (iter != getChargeValues().end()) {
    totalCharge += (*iter++);
  }

  return totalCharge;
}

void EUTelFFClusterImpl::getCenterOfGravityShift(float& xCoG, float& yCoG) const {

  int xSize, ySize;
  getClusterSize(xSize, ySize);

  float normalization = 0;
  float tempX = 0;
  float tempY = 0;

  int iPixel = 0;
  for (int yPixel = -1 * (ySize / 2); yPixel <= (ySize / 2); yPixel++) {
    for (int xPixel = -1 * (xSize / 2); xPixel <= (xSize / 2); xPixel++) {
      normalization += getChargeValues()[iPixel];
      tempX         += xPixel * getChargeValues()[iPixel];
      tempY         += yPixel * getChargeValues()[iPixel];
      ++iPixel;
    }
  }

  if ( normalization != 0)  {
    xCoG = tempX / normalization;
    yCoG = tempY / normalization;
  } else {
    xCoG = 0;
    yCoG = 0;
  }

}

void EUTelFFClusterImpl::getCenterOfGravityShift(float& xCoG, float& yCoG, int xSize, int ySize) const {

  int xCluSize, yCluSize;
  getClusterSize(xCluSize, yCluSize);
  
  if ( ( xSize >= xCluSize ) && ( ySize >= yCluSize ) ) {
    getCenterOfGravityShift(xCoG, yCoG);
    return;
  } 

  float normalization = 0;
  float tempX         = 0; 
  float tempY         = 0;
  int   iPixel        = 0;

  for (int yPixel = -1 * (yCluSize / 2); yPixel <= (yCluSize / 2); yPixel++) {
    for (int xPixel = -1 * (xCluSize / 2); xPixel <= (xCluSize / 2); xPixel++) {
      if ( ( xPixel >= -1 * (xSize / 2) ) &&  ( xPixel <= (xSize / 2) ) &&
	   ( yPixel >= -1 * (ySize / 2) ) &&  ( yPixel <= (ySize / 2) ) ) {
	normalization += getChargeValues()[iPixel];
	tempX         += xPixel * getChargeValues()[iPixel];
	tempY         += yPixel * getChargeValues()[iPixel];
      } 
      ++iPixel;
    }
  }
  
  if ( normalization != 0 ) {
    xCoG = tempX / normalization;
    yCoG = tempY / normalization;
  } else {
    xCoG = 0;
    yCoG = 0;
  }
}


void EUTelFFClusterImpl::getCenterOfGravityShift(float& xCoG, float& yCoG, int nPixel) const {
  
  int xSize, ySize;
  getClusterSize(xSize, ySize);

  if ( nPixel >= xSize * ySize ) {
    getCenterOfGravityShift(xCoG, yCoG);
    return ; 
  }

  map<int, float> highSignalPixel;
  FloatVec        vectorCopy(getChargeValues());
  int             iPixel = 0;
  while ( iPixel != nPixel - 1 ) {
    
    float maxSignal = 0;
    int   maxIndex  = 0;
    int   index;
    FloatVec::iterator maxIter;
    FloatVec::iterator iter = vectorCopy.begin();
    
    while ( iter != vectorCopy.end() ) {
      if ( *iter > maxSignal ) {
	maxSignal = *(iter);
	maxIndex  = index;
	maxIter   = iter;
      }
      ++index; ++iter;
    }
    highSignalPixel.insert( make_pair(maxIndex, maxSignal) ) ;
    vectorCopy.erase(maxIter);
    ++iPixel;
  }

  iPixel = 0;
  float normalization = 0;
  float tempX         = 0;
  float tempY         = 0;
  map<int , float>::iterator mapIter;
  for (int yPixel = -1 * (ySize / 2); yPixel <= (ySize / 2); yPixel++) {
    for (int xPixel = -1 * (xSize / 2); xPixel <= (xSize / 2); xPixel++) {
      mapIter = highSignalPixel.find( iPixel );
      if ( mapIter != highSignalPixel.end() ) {
      	normalization += mapIter->second;
      	tempX         += xPixel * mapIter->second;
      	tempY         += yPixel * mapIter->second;
      }
      ++iPixel;
    }
  }

  if ( normalization != 0 ) {
    xCoG = tempX / normalization;
    yCoG = tempY / normalization;
  } else {
    xCoG = 0;
    yCoG = 0;
  }

}


void EUTelFFClusterImpl::getCenterOfGravity(float& xCoG, float& yCoG) const {

  int xSeed, ySeed;
  getSeedCoord(xSeed, ySeed);

  getCenterOfGravityShift(xCoG, yCoG);
  
  xCoG += xSeed;
  yCoG += ySeed;

  
}
  
  
