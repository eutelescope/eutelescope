
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
#include "EUTelFFClusterImpl.h"
#include "EUTelVirtualCluster.h"
#include "EUTelExceptions.h"

// lcio includes <.h>
#include <LCIOTypes.h>
#include <IMPL/TrackerDataImpl.h>
#include <Exceptions.h>

// system includes
#include <map>
#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cassert>

using namespace eutelescope;
using namespace IMPL;
using namespace std;


EUTelFFClusterImpl::EUTelFFClusterImpl(TrackerDataImpl * data) : EUTelVirtualCluster(data) {
  _trackerData = data;
  _noiseValues.clear();
  _noiseSetSwitch = false;
}


float EUTelFFClusterImpl::getDistance(EUTelVirtualCluster * otherCluster) const {

  int xOtherSeed, yOtherSeed;
  otherCluster->getCenterCoord(xOtherSeed, yOtherSeed);

  int xThisSeed, yThisSeed;
  this->getCenterCoord(xThisSeed, yThisSeed);

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

  int xSize=0, ySize=0;
  getClusterSize(xSize, ySize);

 
  float normalization = 0.;
  float tempX = 0.;
  float tempY = 0.;

  unsigned int iPixel = 0;
  int ii=0;
  for (int yPixel = -1 * (ySize-1) / 2; yPixel <= (ySize-1) / 2; yPixel++) {
    for (int xPixel = -1 * (xSize-1) / 2; xPixel <= (xSize-1) / 2; xPixel++) {
        if( _trackerData->getChargeValues().size()<= iPixel) break;
        normalization += _trackerData->getChargeValues()[iPixel];
        if( _trackerData->getChargeValues()[iPixel] > 0 ){
	  tempX         += xPixel * _trackerData->getChargeValues()[iPixel];
	  tempY         += yPixel * _trackerData->getChargeValues()[iPixel];
// 	  printf("x/y/adc,  %d, %d, %.0f\n", xPixel, yPixel,  _trackerData->getChargeValues()[iPixel]);
	  ii++;
        }      
        iPixel++;    
    }
  }

  
  if ( abs(normalization) > 1e-12 )  {
    xCoG = tempX / normalization;
    yCoG = tempY / normalization;
  } else {
    xCoG = 0.;
    yCoG = 0.;
  }
//   printf("size %d %d, nPixel %d  \n", xSize, ySize, ii);
//   printf("%5.2f %5.2f  %5.2f %5.2f  norm=%5.2f\n", xCoG, yCoG, tempX,tempY, normalization);

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

  for (int yPixel = -1 * (yCluSize-1) / 2; yPixel <= (yCluSize-1) / 2; yPixel++) {
    for (int xPixel = -1 * (xCluSize-1) / 2; xPixel <= (xCluSize-1) / 2; xPixel++) {
      if ( ( xPixel >= -1 * (xSize-1) / 2 ) &&  ( xPixel <= (xSize-1) / 2 ) &&
	   ( yPixel >= -1 * (ySize-1) / 2 )  &&  ( yPixel <= (ySize-1) / 2 ) ) {
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
  while ( iPixel != nPixel  ) {
    float maxSignal = -1 * numeric_limits<float >::max();
    int   maxIndex  = 0;
    int   index     = 0;
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
    (*maxIter) = -1 * numeric_limits<float >::max();
    ++iPixel;
  }

  iPixel = 0;
  float normalization = 0;
  float tempX         = 0;
  float tempY         = 0;
  map<int , float>::iterator mapIter;
  for (int yPixel = -1 * (ySize-1) / 2 ; yPixel <= (ySize-1) / 2 ; yPixel++) {
    for (int xPixel = -1 * (xSize-1) / 2 ; xPixel <= (xSize-1) / 2 ; xPixel++) {
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

  int xSeed=0, ySeed=0;
  getCenterCoord(xSeed, ySeed);

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

  if ( static_cast< size_t >(nPixel) >= vectorCopy.size() ) return getTotalCharge() ;

  sort(vectorCopy.begin(), vectorCopy.end(), greater<float>());

  vector<float >::iterator iter = vectorCopy.begin();
  float charge = 0;
  while ( iter != vectorCopy.begin() + nPixel ) {
    charge += *(iter);
    ++iter;
  }
  return charge;

}

std::vector<float> EUTelFFClusterImpl::getClusterCharge(std::vector<int > nPixels) const {

  vector< float > clusterSignal;

  vector<float> vectorCopy(_trackerData->getChargeValues());
  sort(vectorCopy.begin(), vectorCopy.end(), greater<float>());
  vector<float >::iterator iter;

  for (unsigned int i = 0; i < nPixels.size(); i++) {
    iter = vectorCopy.begin();
    float charge = 0;

    if ( static_cast< size_t >(nPixels[i]) >= vectorCopy.size() ) {

      clusterSignal.push_back( getTotalCharge() );

    } else {

      while ( iter != vectorCopy.begin() + nPixels[i] ) {
        charge += (*iter);
        ++iter;
      }
      clusterSignal.push_back(charge);

    }
  }

  return clusterSignal;

}

void EUTelFFClusterImpl::setNoiseValues(std::vector<float > noiseValues ) {

  // first check that the noiseValues sizes is the same of the
  // TrackerData
  if ( noiseValues.size() != _trackerData->getChargeValues().size() ) {
    _noiseSetSwitch = false;
    throw IncompatibleDataSetException("The noiseValues size is different from the TrackerData size");
  }

  _noiseValues = noiseValues;
  _noiseSetSwitch = true;

}

vector<float > EUTelFFClusterImpl::getNoiseValues() const {

  if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");
  return _noiseValues;
}

float EUTelFFClusterImpl::getClusterNoise() const {

  if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");

  float squaredSum = 0;
  vector<float >::const_iterator iter = _noiseValues.begin();
  while ( iter != _noiseValues.end() ) {
    squaredSum += pow( (*iter), 2 );
    ++iter;
  }
  return sqrt( squaredSum );

}

float EUTelFFClusterImpl::getClusterSNR() const {

  if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");
  float clusterNoise = getClusterNoise();
  if ( clusterNoise == 0 )  return 0.;
  float clusterSignal = getTotalCharge();
  return clusterSignal / clusterNoise;
}

float EUTelFFClusterImpl::getSeedSNR() const {

  if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");
  vector<float >::const_iterator chargeBegin    = _trackerData->getChargeValues().begin();
  vector<float >::const_iterator seedChargeIter = max_element( chargeBegin, _trackerData->getChargeValues().end() );
  vector<float >::const_iterator seedNoiseIter  = _noiseValues.begin() + ( seedChargeIter - chargeBegin );
  return (*seedChargeIter) / (*seedNoiseIter);
}

float EUTelFFClusterImpl::getClusterSNR(int nPixel) const {

  if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");
  int xSize, ySize;
  getClusterSize(xSize, ySize);

  if ( nPixel >= xSize * ySize )
    return getClusterSNR();

  map<int, float > highSignalPixel;
  vector<float >   vectorCopy( _trackerData->getChargeValues() );
  int              iPixel = 0;
  while ( iPixel != nPixel ) {
    float maxSignal = -1 * numeric_limits<float >::max();
    int   maxIndex  = 0;
    int   index     = 0;
    vector<float >::iterator maxIter;
    vector<float >::iterator iter = vectorCopy.begin();

    while ( iter != vectorCopy.end() ) {
      if ( *iter > maxSignal ) {
        maxSignal = (*iter);
        maxIndex  = index;
        maxIter   = iter;
      }
      ++index; ++iter;
    }
    highSignalPixel.insert( make_pair(maxIndex, maxSignal) );
    (*maxIter) = -1 * numeric_limits<float >::max();
    ++iPixel;
  }

  float signal = 0, noise2 = 0;
  map<int, float >::iterator mapIter = highSignalPixel.begin();
  while ( mapIter != highSignalPixel.end() ) {
    signal += mapIter->second;
    noise2 += pow( _noiseValues[mapIter->first], 2 );
    ++mapIter;
  }
  if ( noise2 == 0 ) return 0;
  return signal / sqrt( noise2 );

}

std::vector<float > EUTelFFClusterImpl::getClusterSNR( std::vector<int > nPixels ) const {

  if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");
  multimap<float, int > clusterSignalMap;
  vector<float >::const_iterator iter = _trackerData->getChargeValues().begin();
  int index = 0;
  while ( iter != _trackerData->getChargeValues().end() ) {
    clusterSignalMap.insert( make_pair( (*iter), index ) );
    ++index; ++iter;
  }

  vector<int >::iterator pixelIter = nPixels.begin();
  vector<float > snr;
  int xSize, ySize;
  getClusterSize(xSize, ySize);

  while ( pixelIter != nPixels.end() ) {

    if ( *pixelIter >= xSize * ySize ) {
      snr.push_back( getClusterSNR() );
    } else {

      map<float, int >::reverse_iterator mapIter = clusterSignalMap.rbegin();
      float signal = 0;
      float noise2 = 0;
      int   iPixel = 0;
      while ( (iPixel < (*pixelIter)) && (mapIter != clusterSignalMap.rend()) ) {
        signal += mapIter->first;
        noise2 += pow( _noiseValues[ mapIter->second], 2 );
        ++mapIter; ++iPixel;
      }
      if ( noise2 == 0 ) snr.push_back(0.);
      else snr.push_back( signal / sqrt( noise2 ) );
    }
    ++pixelIter;
  }
  return snr;
}



float EUTelFFClusterImpl::getClusterSNR(int xSize, int ySize) const {

  if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");

  int xCluSize, yCluSize;
  getClusterSize(xCluSize, yCluSize);

  if ( ( xSize >= xCluSize ) && ( ySize >= yCluSize ) ) {
    return getClusterSNR();
  }

  int   iPixel = 0;
  float signal = 0, noise2 = 0;

  for ( int yPixel = -1 * (yCluSize-1) / 2 ; yPixel <= (yCluSize-1) / 2 ; yPixel++ ) {
    for ( int xPixel = -1 * (xCluSize-1) / 2 ; xPixel <= (xCluSize-1) / 2 ; xPixel++ ) {
      if ( ( xPixel >= -1 * (xSize-1) / 2  ) &&  ( xPixel <= (xSize-1) / 2  ) &&
           ( yPixel >= -1 * (ySize-1) / 2  ) &&  ( yPixel <= (ySize-1) / 2  ) ) {
        signal += _trackerData->getChargeValues()[iPixel];
        noise2 += pow( _noiseValues[iPixel] , 2 );
      }
      ++iPixel;
    }
  }
  if ( noise2 == 0 ) return 0;
  return signal / sqrt( noise2 );

}

float EUTelFFClusterImpl::getClusterCharge(int xSize, int ySize) const {

  int xCluSize, yCluSize;
  getClusterSize(xCluSize, yCluSize);

  if ( ( xSize >= xCluSize ) && ( ySize >= yCluSize ) ) {
    return getTotalCharge();
  }

  int iPixel = 0;
  float charge = 0;

  for (int yPixel = -1 * (yCluSize-1) / 2; yPixel <= (yCluSize-1) / 2; yPixel++) {
    for (int xPixel = -1 * (xCluSize-1) / 2 ; xPixel <= (xCluSize-1) / 2 ; xPixel++) {
      if ( ( xPixel >= -1 * (xSize-1) / 2  ) &&  ( xPixel <= (xSize-1) / 2  ) &&
           ( yPixel >= -1 * (ySize-1) / 2  ) &&  ( yPixel <= (ySize-1) / 2  ) ) {
        charge += _trackerData->getChargeValues()[iPixel];
      }
      ++iPixel;
    }
  }
  return charge;
}

void EUTelFFClusterImpl::print(std::ostream& os ) const {

  int xSize, ySize, xSeed, ySeed;
  float xShift, yShift, xShift9, yShift9, xShift3x3, yShift3x3;
  ClusterQuality quality = getClusterQuality();
  getClusterSize(xSize,ySize);
  getSeedCoord(xSeed, ySeed);
  getCenterOfGravityShift(xShift, yShift);
  getCenterOfGravityShift(xShift9, yShift9, 9);
  getCenterOfGravityShift(xShift3x3, yShift3x3, 3, 3);

  float noise = 0., SNR = 0., SNR9 = 0., SNR3x3 = 0.;
  if ( _noiseSetSwitch ) {
    noise  = getClusterNoise();
    SNR    = getClusterSNR();
    SNR9   = getClusterSNR(9);
    SNR3x3 = getClusterSNR(3,3);
  }

  int bigspacer = 23;

  os   <<  setw(bigspacer) << setiosflags(ios::left) << "Fixed frame cluster "<< "(" << xSize << ", " << ySize << ")\n"
       <<  setw(bigspacer) <<  "Cluster ID " << getClusterID() << " on detector " << getDetectorID() << "\n"
       <<  setw(bigspacer) <<  "Cluster quality " << quality << "\n"
       <<  setw(bigspacer) <<  "Cluster total charge " << getTotalCharge() << "\n"
       <<  setw(bigspacer) <<  "Cluster charge (9) " << getClusterCharge(9) << "\n"
       <<  setw(bigspacer) <<  "Cluster charge (3x3) " << getClusterCharge(3,3) << "\n"
       <<  setw(bigspacer) <<  "Seed charge " << getSeedCharge() << " in (" << xSeed << ", " << ySeed << ")\n"
       <<  setw(bigspacer) <<  "CoG shift "<< "(" << xShift << ", " << yShift << ")\n"
       <<  setw(bigspacer) <<  "CoG(9) shift " << "(" << xShift9 << ", " << yShift9 << ")\n"
       <<  setw(bigspacer) <<  "CoG(3x3) shift " << "(" << xShift3x3 << ", " << yShift3x3 << ")\n" ;
  if ( _noiseSetSwitch ) {
    os << setw(bigspacer)  <<  "Cluster noise " << noise << "\n"
       << setw(bigspacer)  <<  "Cluster SNR " << SNR << "\n"
       << setw(bigspacer)  <<  "Cluster SNR(9) " << SNR9 << "\n"
       << setw(bigspacer)  <<  "Cluster SNR(3x3) " << SNR3x3 << "\n";
  }

  os   << resetiosflags(ios::left);

  int spacer = 14;

  os << "|";
  for ( int i = 0; i < spacer - 1; i++ ) {
    os << "-";
  }

  for (int xPixel = -1 * (xSize-1) / 2 ; xPixel <= (xSize-1) / 2 ; xPixel++) {
    os << "|";
    for ( int i = 0; i < spacer - 1; i++ ) {
      os << "-";
    }
  }
  os <<"|\n"
     << "|" << setw(spacer - 1) << " x / y ";
  for (int xPixel = -1 * (xSize-1) / 2 ; xPixel <= (xSize-1) / 2 ; xPixel++) {
    os << "|" << setw(spacer - 1 ) << xSeed - ( -1 * xPixel );
  }
  os << "|\n";
  os << "|";
  for ( int i = 0; i < spacer - 1; i++ ) {
    os << "-";
  }
  for (int xPixel = -1 * (xSize-1) / 2 ; xPixel <= (xSize-1) / 2 ; xPixel++) {
    os << "|";
    for ( int i = 0; i < spacer - 1; i++ ) {
      os << "-";
    }
  }
  os <<"|\n";
  int iPixel = 0;
  for (int yPixel = -1 * (ySize-1) / 2 ; yPixel <= (ySize-1) / 2 ; yPixel++) {
    os << "|" << setw(spacer - 1) << ySeed - ( -1 * yPixel ) << "|";
    for (int xPixel = -1 * (xSize-1) / 2 ; xPixel <= (xSize-1) / 2 ; xPixel++) {
      os <<  setw(spacer - 1) << _trackerData->getChargeValues()[iPixel] << "|" ;
      ++iPixel;
    }
    os << "\n";
    os << "|";
    for ( int i = 0; i < spacer - 1; i++ ) {
      os << "-";
    }
    for (int xPixel = -1 * (xSize-1) / 2 ; xPixel <= (xSize-1) / 2 ; xPixel++) {
      os << "|";
      for ( int i = 0; i < spacer - 1; i++ ) {
        os << "-";
      }
    }
    if ( yPixel == (ySize-1) / 2  ) os << "|";
    else os << "|\n";
  }
}

