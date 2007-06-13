// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author:  Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version: $Id: EUTelFFClusterImpl.cc,v 1.11 2007-06-13 18:02:32 bulgheroni Exp $

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
#include "EUTelVirtualCluster.h"

// lcio includes <.h>
#include <LCIOTypes.h>
#include <IMPL/TrackerDataImpl.h>

// system includes
#include <map>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace eutelescope;
using namespace IMPL;
using namespace std;


EUTelFFClusterImpl::EUTelFFClusterImpl(TrackerDataImpl * data) : EUTelVirtualCluster(data) { _trackerData = data; } 


float EUTelFFClusterImpl::getDistance(EUTelVirtualCluster * otherCluster) const {

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
  
  FloatVec::const_iterator iter = _trackerData->getChargeValues().begin();

  while (iter != _trackerData->getChargeValues().end()) {
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
      normalization += _trackerData->getChargeValues()[iPixel];
      tempX         += xPixel * _trackerData->getChargeValues()[iPixel];
      tempY         += yPixel * _trackerData->getChargeValues()[iPixel];
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
	normalization += _trackerData->getChargeValues()[iPixel];
	tempX         += xPixel * _trackerData->getChargeValues()[iPixel];
	tempY         += yPixel * _trackerData->getChargeValues()[iPixel];
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
  FloatVec        vectorCopy(_trackerData->getChargeValues());
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
  
  
void EUTelFFClusterImpl::setClusterQuality(ClusterQuality quality) {

  lcio::long64 cell1 = static_cast<lcio::long64> (_trackerData->getCellID1()) ;

  int rhs = 15;
  lcio::long64  emptyMask = ~( 0x1F << rhs );
  lcio::long64  maskedQuality = ( (static_cast<int> (quality) & 0x1F ) << rhs );

  // first apply an empty mask for the quality bit ranges
  cell1 = cell1 & emptyMask;
  
  // now apply the maskedQuality
  cell1 = cell1 | maskedQuality;

  // apply the changes
  _trackerData->setCellID1(cell1);

}


float EUTelFFClusterImpl::getSeedCharge() const {
  
  return *max_element( _trackerData->getChargeValues().begin(),
		       _trackerData->getChargeValues().end() );
}

float EUTelFFClusterImpl::getClusterCharge(int nPixel) const {

  vector<float > vectorCopy(_trackerData->getChargeValues());
  sort(vectorCopy.begin(), vectorCopy.end(), greater<float>());

  vector<float >::iterator iter = vectorCopy.begin();
  float charge = 0;
  while ( iter != vectorCopy.begin() + nPixel ) {
    charge += *(iter);
    ++iter;
  }
  return charge;

}

vector<float> EUTelFFClusterImpl::getClusterCharge(vector<int > nPixels) const {
  
  vector< float > clusterSignal;

  vector<float> vectorCopy(_trackerData->getChargeValues()); 
  sort(vectorCopy.begin(), vectorCopy.end(), greater<float>());
  vector<float >::iterator iter;

  for (unsigned int i = 0; i < nPixels.size(); i++) {
    iter = vectorCopy.begin();
    float charge = 0;
    while ( iter != vectorCopy.begin() + nPixels[i] ) {
      charge += (*iter);
      ++iter;
    }
    clusterSignal.push_back(charge);
  }
  
  return clusterSignal;

}

float EUTelFFClusterImpl::getClusterCharge(int xSize, int ySize) const {
  
  int xCluSize, yCluSize;
  getClusterSize(xCluSize, yCluSize);
  
  if ( ( xSize >= xCluSize ) && ( ySize >= yCluSize ) ) {
    return getTotalCharge();
  }

  int iPixel = 0;
  float charge = 0;

  for (int yPixel = -1 * (yCluSize / 2); yPixel <= (yCluSize / 2); yPixel++) {
    for (int xPixel = -1 * (xCluSize / 2); xPixel <= (xCluSize / 2); xPixel++) {
      if ( ( xPixel >= -1 * (xSize / 2) ) &&  ( xPixel <= (xSize / 2) ) &&
	   ( yPixel >= -1 * (ySize / 2) ) &&  ( yPixel <= (ySize / 2) ) ) {
	charge += _trackerData->getChargeValues()[iPixel];
      } 
      ++iPixel;
    }
  }
  return charge;
}
