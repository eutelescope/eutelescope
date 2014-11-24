// Version: $Id$
// Author Christian Takacs, SUS UNI HD <mailto:ctakacs@rumms.uni-mannheim.de>
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 */

// personal includes ".h"
#include "EUTelBrickedClusterImpl.h"
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


EUTelBrickedClusterImpl::EUTelBrickedClusterImpl(TrackerDataImpl * data) : EUTelVirtualCluster(data) {
  _trackerData = data;

  _noiseValues.clear();
  _noiseSetSwitch = false;
}

float EUTelBrickedClusterImpl::getDistance(EUTelVirtualCluster * otherCluster) const {

  int xOtherSeed, yOtherSeed;
  otherCluster->getCenterCoord(xOtherSeed, yOtherSeed);

  int xThisSeed, yThisSeed;
  this->getCenterCoord(xThisSeed, yThisSeed);

  return sqrt( pow(static_cast<double> (xThisSeed - xOtherSeed), 2) + pow(static_cast<double> (yThisSeed - yOtherSeed), 2) );

}

float EUTelBrickedClusterImpl::getExternalRadius() const {

    //this holds for a seed pixel with one layer of surrounding pixels:
    return sqrt(1.25f); //= sqrt (1^2 + 0.5^2)

    //NOTE fix this if you switch the implemenatation to a variable size (size aka surrounding layers):
    //NOTE for two surrounding layers this would be sqrt(2^2 + 1^2) and so on:
    /*
    int xSize, ySize;
    getClusterSize(xSize, ySize);
    if (xSize!=ySize)
        //Exception
    else
    {
        //xSize needs to be the distance |seedPixel-outermostPixelInTheSameRow| now:
        xSize /= 2; //e.g. 3 becomes 1, 5 becomes 2, ...
        return sqrt( pow( static_cast<double>(xSize) , 2) + pow( static_cast<double>(xSize*0.5) , 2) );
    }

    */
}

float EUTelBrickedClusterImpl::getTotalCharge() const {

  float totalCharge = 0;

  //!will be okay if we set zeros for the unwanted pixels
  FloatVec vectorCopy(_trackerData->getChargeValues());
  setOutsiderValuesInVectorInterpretedAsBrickedMatrix( static_cast< std::vector< float>& > (vectorCopy), 0.0f );

  FloatVec::const_iterator iter = vectorCopy.begin();

  while (iter != vectorCopy.end()) {
    totalCharge += (*iter++);
  }

  return totalCharge;
}

void EUTelBrickedClusterImpl::getCenterOfGravityShift(float& xCoG, float& yCoG) const {
    //!will be okay if we set zeros for the unwanted pixels
    FloatVec vectorCopy(_trackerData->getChargeValues());
    setOutsiderValuesInVectorInterpretedAsBrickedMatrix( static_cast< std::vector< float>& > (vectorCopy), 0.0f );

    //streamlog_out( MESSAGE4 ) << "RUNNING EUTelBrickedClusterImpl::getCenterOfGravityShift() " << endl;

    //NOTE this whole class could use a flag showing, if even rows are skewed left or right!
    //NOTE then the offset would become +/-0.5 depending on the flag, instead of fixed -0.5 here!

    int xSize, ySize;
    getClusterSize(xSize, ySize);
    if (!(xSize==3 && ySize==3))
    {
        //NOTE fix this if you switch the implemenatation to a variable size:
        streamlog_out( WARNING4 ) << " BRICKED PIXEL FIXED FRAME CLUSTER SIZE MUST BE 3x3!!!" << endl;
        streamlog_out( WARNING4 ) << " BUT IT IS " << xSize << "x" << ySize << "!!!" << endl;
    }

    int seedX, seedY;
    getSeedCoord(seedX, seedY);

    bool bSeedRowIsEven = false;
    if (seedY % 2 == 0) bSeedRowIsEven = true; //seed pixel's row is even

    float normalization = 0;
    float tempX = 0;
    float tempY = 0;

    float skewCorrectionX;
    bool  currRowIsEven;
    int   iPixel = 0;

    for (int yPixel = -1 * (ySize / 2); yPixel <= (ySize / 2); yPixel++)
    {
        //Pixel Choice Correction not needed. Unwanted pixels are set to 0.

        //Coordinate Correction:
        currRowIsEven = (bSeedRowIsEven && (yPixel % 2 == 0)) || (!bSeedRowIsEven && (yPixel % 2 != 0)); //even+even or odd+odd
        if (currRowIsEven)
            skewCorrectionX = -0.5;
        else
            skewCorrectionX = 0.0;

        //Go:
        for (int xPixel = -1 * (xSize / 2); xPixel <= (xSize / 2); xPixel++)
        {
            normalization += vectorCopy[iPixel];
            tempX         += (xPixel+skewCorrectionX) * vectorCopy[iPixel];
            tempY         +=  yPixel                  * vectorCopy[iPixel];
            ++iPixel;
        }
    }

    if ( normalization != 0)
    {
        xCoG = tempX / normalization;
        yCoG = tempY / normalization;
    }
    else
    {
        xCoG = 0;
        yCoG = 0;
    }

}

void EUTelBrickedClusterImpl::getCenterOfGravityShift(float& xCoG, float& yCoG, int , int) const {

  //!UGLY HACK
  streamlog_out( WARNING4 ) << "[getCenterOfGravityShift(float& xCoG, float& yCoG, int xSize, int ySize)] DOES NOT MAKE SENSE FOR A BRICKED PIXEL" << endl;
  streamlog_out( WARNING4 ) << "[getCenterOfGravityShift(float& xCoG, float& yCoG, int xSize, int ySize)] USING THE WHOLE FRAME TO COMPUTE THE SHIFT!" << endl;
  getCenterOfGravityShift(xCoG,yCoG);

}

void EUTelBrickedClusterImpl::getCenterOfGravityShift(float& xCoG, float& yCoG, int nPixel) const {

  //! UGLY HACK:
  //! HERE WE SET STRONG NEGATIVE VALUES OUTSIDE (IN THE UNWANTED PIXELS)
  //! SUCH THAT THEY WILL NOT BE REGARDED AS SIGNIFICANT ONES.
  //! WE JUST HAVE TO PAY ATTENTION FOR THESE PIXELS NOT TO BE TAKEN INTO
  //! ACCOUNT FOR COG SHIFT CALCULATION !!
  //! So the max. nPixel value allowed is
  //! [numberOfPixels_in_fixedFrameLyingBelow - numberOfPixelsWeWantToIgnore]
  //! = (9-2) in the 3x3 case. Adjust this to be more generic if you implement a variable cluster size!

  //! Implementation for 3x3 only!

  if ( /* (size_t) nPixel >= _trackerData->getChargeValues().size() */ (nPixel>=7) )
  {
    //!UGLY HACK! 7 fits only for a 3x3 cluster lying below!
    //! 3x3 means: one layer of pixels surrounding seed.
    //! Again: more than 7 would be bad, because 7 is all we want to
    //! take into account (exactly one layer of surrounding pixels).
    //! What we do is walking blindly across the 3x3 fixed frame below,
    //! assuming, that the charge values are set in a way that ensures
    //! the correct choice of most significant pixels!
    //! So if we went on with more than 7, then the fake strong negative
    //! values would destroy the cog shift!
    //! Why fake strong negative values?
    //! See below. They are used for these pixels to be ignored by
    //! the maximum-finding loop.
    getCenterOfGravityShift(xCoG, yCoG);
  }
  else
  {
        int xSize, ySize;
        getClusterSize(xSize, ySize);
        if (!(xSize==3 && ySize==3))
        {
            streamlog_out( WARNING4 ) << " BRICKED PIXEL FIXED FRAME CLUSTER SIZE MUST BE 3x3!!!" << endl;
            streamlog_out( WARNING4 ) << " BUT IT IS " << xSize << "x" << ySize << "!!!" << endl;
        }

        FloatVec vectorCopy(_trackerData->getChargeValues());
        //! set fake negative values such that the unwanted pixels are not among the significant ones!
        setOutsiderValuesInVectorInterpretedAsBrickedMatrix( static_cast< std::vector< float>& > (vectorCopy), (-1) * numeric_limits<float >::max() );

        map<int, float> highSignalPixel;
        int             iPixel = 0;

        //BEGIN sort the n highest pixels into the highSignalPixel map (will be sorted decending via index as the key, but that doesn't matter)
        //BEGIN  (the highest pixel will be put in first but that doesn't matter either)
        while ( iPixel != nPixel  )
        {
            float maxSignal = (-1) * numeric_limits<float >::max();
            int   maxIndex  = 0;
            int   index     = 0;
            FloatVec::iterator maxIter;
            FloatVec::iterator iter = vectorCopy.begin();

            while ( iter != vectorCopy.end() )
            {
                if ( *iter > maxSignal )
                {
                    maxSignal = *(iter);
                    maxIndex  = index;
                    maxIter   = iter;
                }
                ++index; ++iter;
            }
            highSignalPixel.insert( make_pair(maxIndex, maxSignal) ) ;
            (*maxIter) = (-1) * numeric_limits<float >::max();
            ++iPixel;
        }
        //END sort the n highest pixels into the highSignalPixel map (will be sorted decending via index as the key, but that doesn't matter))


        int seedX, seedY;
        getSeedCoord(seedX, seedY);
        bool bSeedRowIsEven = false;
        if (seedY % 2 == 0) bSeedRowIsEven = true; //seed pixel's row is even
        float skewCorrectionX;
        bool  currRowIsEven;

        iPixel = 0;
        float normalization = 0;
        float tempX         = 0;
        float tempY         = 0;
        map<int , float>::iterator mapIter;

        for (int yPixel = -1 * (ySize / 2); yPixel <= (ySize / 2); yPixel++)
        {
            //Pixel Choice Correction not needed. Unwanted pixels are set to 0.

            //Coordinate Correction:
            currRowIsEven = (bSeedRowIsEven && (yPixel % 2 == 0)) || (!bSeedRowIsEven && (yPixel % 2 != 0)); //even+even or odd+odd
            if (currRowIsEven)
                skewCorrectionX = -0.5;
            else
                skewCorrectionX = 0.0;

            //Go:
            for (int xPixel = -1 * (xSize / 2); xPixel <= (xSize / 2); xPixel++)
            {
                mapIter = highSignalPixel.find( iPixel );
                if ( mapIter != highSignalPixel.end() ) //not found in map -> mapIter == map.end()
                {
                    normalization += mapIter->second;
                    tempX         += (xPixel+skewCorrectionX) * mapIter->second;
                    tempY         += yPixel                   * mapIter->second;
                }
                ++iPixel;
            }
        }

        if ( normalization != 0 )
        {
            xCoG = tempX / normalization;
            yCoG = tempY / normalization;
        }
        else
        {
            xCoG = 0;
            yCoG = 0;
        }
  }

}

void EUTelBrickedClusterImpl::getCenterOfGravityShiftWithOutGlobalSeedCoordinateCorrection(float& xCoG, float& yCoG) const {

    // for eta we need the cog offset from the original pixel
    // -
    // but we do NOT the correction of the x coordinate
    // (seed_XCoordinateCorrection is: in case the seed pixel lies in an even row, the cog shift will reflect an additional 0.5 offset)
    //
    // So:
    //   let the normal COG Shift (with(!) seed coordinate correction) be computed and adjust the result:
    //   if the seed pixel row is even, then add the 0.5 that was substracted before!

    getCenterOfGravityShift(xCoG, yCoG);

    int seedX, seedY;
    getSeedCoord(seedX, seedY);
    if (seedY %2 == 0)
    {
        xCoG += 0.5f;
    }
}

void EUTelBrickedClusterImpl::getCenterOfGravityShiftWithOutGlobalSeedCoordinateCorrection(float& xCoG, float& yCoG, int nPixel) const {

    // for eta we need the cog offset from the original pixel
    // -
    // but we do NOT the correction of the x coordinate
    // (seed_XCoordinateCorrection is: in case the seed pixel lies in an even row, the cog shift will reflect an additional 0.5 offset)
    //
    // So:
    //   let the normal COG Shift (with(!) seed coordinate correction) be computed and adjust the result:
    //   if the seed pixel row is even, then add the 0.5 that was substracted before!

    getCenterOfGravityShift(xCoG, yCoG, nPixel);

    int seedX, seedY;
    getSeedCoord(seedX, seedY);
    if (seedY %2 == 0)
    {
        xCoG += 0.5f;
    }
}

void EUTelBrickedClusterImpl::getCenterOfGravity(float& xCoG, float& yCoG) const {

  int xSeed, ySeed;
  getCenterCoord(xSeed, ySeed);

  getCenterOfGravityShift(xCoG, yCoG);

  xCoG += xSeed;
  yCoG += ySeed;

}

float EUTelBrickedClusterImpl::getSeedCharge() const {

  return *max_element( _trackerData->getChargeValues().begin(),
                       _trackerData->getChargeValues().end() );
}

float EUTelBrickedClusterImpl::getClusterCharge(int nPixel) const {

  //!HACK TAKI, for an explanation see getCenterOfGravityShift(float& xCoG, float& yCoG, int nPixel)
  if ( /*(size_t) nPixel >= vectorCopy.size()*/ nPixel>=7 ) return getTotalCharge();

  vector<float > vectorCopy(_trackerData->getChargeValues());
  //!will only be okay if we set -inf for the unwanted pixels
  setOutsiderValuesInVectorInterpretedAsBrickedMatrix( static_cast< std::vector< float>& > (vectorCopy), (-1) * numeric_limits<float >::max() );

  sort(vectorCopy.begin(), vectorCopy.end(), greater<float>());

  vector<float >::iterator iter = vectorCopy.begin();
  float charge = 0;
  while ( iter != vectorCopy.begin() + nPixel ) {
    charge += *(iter);
    ++iter;
  }
  return charge;
}

std::vector<float> EUTelBrickedClusterImpl::getClusterCharge(std::vector<int> nPixels) const {

  //!HACK TAKI, for an explanation see getCenterOfGravityShift(float& xCoG, float& yCoG, int nPixel)

  vector<float> vectorCopy(_trackerData->getChargeValues());
  //!will only be okay if we set -inf for the unwanted pixels
  setOutsiderValuesInVectorInterpretedAsBrickedMatrix( static_cast< std::vector< float>& > (vectorCopy), (-1) * numeric_limits<float >::max() );

  vector< float > clusterSignal;

  sort(vectorCopy.begin(), vectorCopy.end(), greater<float>());
  vector<float >::iterator iter;

  for (unsigned int i = 0; i < nPixels.size(); i++)
  {

    float charge = 0;

    if ( /*(size_t)(nPixels[i]) >= vectorCopy.size()*/ nPixels[i]>=7 )
    {
      clusterSignal.push_back( getTotalCharge() );
    }
    else
    {
      iter = vectorCopy.begin();
      while ( iter != vectorCopy.begin() + nPixels[i] )
      {
        charge += (*iter);
        ++iter;
      }
      clusterSignal.push_back(charge);
    }
  }

  return clusterSignal;
}

void EUTelBrickedClusterImpl::setNoiseValues(std::vector<float > noiseValues ) {

  // first check that the noiseValues sizes is the same of the
  // TrackerData
  if ( noiseValues.size() != _trackerData->getChargeValues().size() ) {
    _noiseSetSwitch = false;
    streamlog_out( ERROR2 ) << "[EUTelBrickedClusterImpl::setNoiseValues()] The noiseValues size is different from the TrackerData size!!!" << endl;
    throw IncompatibleDataSetException("[EUTelBrickedClusterImpl::setNoiseValues()] The noiseValues size is different from the TrackerData size");
  }

  _noiseValues = noiseValues;
  setOutsiderValuesInVectorInterpretedAsBrickedMatrix( _noiseValues, 0.0f ); //!HACK TAKI
  _noiseSetSwitch = true;

}

vector<float > EUTelBrickedClusterImpl::getNoiseValues() const {

  if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");
  return _noiseValues;
}

float EUTelBrickedClusterImpl::getClusterNoise() const {

  if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");

  float squaredSum = 0; //!will be okay if we set zeros @ noise for the unwanted pixels (which we did)
  vector<float >::const_iterator iter = _noiseValues.begin();
  while ( iter != _noiseValues.end() ) {
    squaredSum += pow( (*iter), 2 );
    ++iter;
  }
  return sqrt( squaredSum );

}

float EUTelBrickedClusterImpl::getClusterSNR() const {

  if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");
  float clusterNoise = getClusterNoise();  //!will be okay if we set zeros @ noise for the unwanted pixels!!! (which we did)
  if ( clusterNoise == 0 )  return 0.;
  float clusterSignal = getTotalCharge();  //!will be okay if we set zeros @ charge for the unwanted pixels!!! (which we do in getTotalCharge() )
  return clusterSignal / clusterNoise;
}

float EUTelBrickedClusterImpl::getSeedSNR() const {

  if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");
  vector<float >::const_iterator chargeBegin    = _trackerData->getChargeValues().begin();
  vector<float >::const_iterator seedChargeIter = max_element( chargeBegin, _trackerData->getChargeValues().end() );
  vector<float >::const_iterator seedNoiseIter  = _noiseValues.begin() + ( seedChargeIter - chargeBegin );
  if (*seedNoiseIter==0.0)
  {
        streamlog_out( ERROR4 ) << "[getSeedSNR()] Just found a noise value == 0.0 !" << endl;
        return 0;
  }
  return (*seedChargeIter) / (*seedNoiseIter);
}

float EUTelBrickedClusterImpl::getClusterSNR(int nPixel) const {

  if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");

  if ( /*(size_t) nPixel >= _trackerData->getChargeValues().size()*/ nPixel >= 7 ) //!HACK TAKI, for an explanation see getCenterOfGravityShift(float& xCoG, float& yCoG, int nPixel)
    return getClusterSNR();

  vector<float>   vectorCopy( _trackerData->getChargeValues() );
  //! set fake negative values such that the unwanted pixels are not among the significant ones!
  setOutsiderValuesInVectorInterpretedAsBrickedMatrix( static_cast< std::vector< float>& > (vectorCopy), (-1) * numeric_limits<float >::max() );
  //! zeros @ noise for unwanted pixels is not even necessary here, but it is done automatically anwayway

  map<int, float > highSignalPixel;
  int              iPixel = 0;
  while ( iPixel != nPixel ) //grab the n highest pixels!
  {
    float maxSignal = (-1) * numeric_limits<float >::max();
    int   maxIndex  = 0;
    int   index     = 0;
    vector<float >::iterator maxIter;
    vector<float >::iterator iter = vectorCopy.begin();

    while ( iter != vectorCopy.end() )
    {
      if ( *iter > maxSignal )
      {
        maxSignal = (*iter);
        maxIndex  = index;
        maxIter   = iter;
      }
      ++index; ++iter;
    }
    highSignalPixel.insert( make_pair(maxIndex, maxSignal) );
    (*maxIter) = (-1) * numeric_limits<float >::max();
    ++iPixel;
  }

  float signal = 0, noise2 = 0;
  map<int, float >::iterator mapIter = highSignalPixel.begin();
  while ( mapIter != highSignalPixel.end() ) //SNR only for those highest pixels. okay.
  {
    signal += mapIter->second;
    noise2 += pow( _noiseValues[mapIter->first], 2 );
    if (_noiseValues[mapIter->first]==0.0)
    {
        streamlog_out( ERROR4 ) << "[getClusterSNR(int nPixel)] Just found a noise value == 0.0 !" << endl;
    }
    ++mapIter;
  }
  if ( noise2 == 0 )
  {
      return 0;
  }
  return signal / sqrt( noise2 );

}

std::vector<float > EUTelBrickedClusterImpl::getClusterSNR( std::vector<int> nPixels ) const {

  if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");

  //! set fake negative values such that the unwanted pixels are not among the significant ones!
  vector<float>   vectorCopy( _trackerData->getChargeValues() );
  setOutsiderValuesInVectorInterpretedAsBrickedMatrix( static_cast< std::vector< float>& > (vectorCopy), (-1) * numeric_limits<float >::max() );

  //! zeros @ noise for unwanted pixels is not even necessary here, but it is done automatically anwayway



//!DEBOUG OUTPUT ENABLE
//#define DEBUG_OUTPUT_FOR_getClusterSNR_ON_INT_VECTOR

//!USE MULTIMAP OR NORMAL MAP
#define USE_MULTIMAP



  #ifdef DEBUG_OUTPUT_FOR_getClusterSNR_ON_INT_VECTOR
    streamlog_out( MESSAGE2 ) << "[getClusterSNR( std::vector<int > nPixels )] THE CLUSTER:" << endl;
    outputBricked3x3MatrixFromVector(vectorCopy);
  #endif



  #ifdef USE_MULTIMAP
    multimap<float, int > clusterSignalMap;
  #else
    map<float, int > clusterSignalMap;
  #endif


  vector<float>::const_iterator iter = vectorCopy.begin();
  int index = 0;
  while ( iter != vectorCopy.end() )
  {
    clusterSignalMap.insert( make_pair( (*iter), index ) ); //sorted by charge value, ascending
    #ifdef DEBUG_OUTPUT_FOR_getClusterSNR_ON_INT_VECTOR
        streamlog_out( MESSAGE2 ) << "Adding a pixel to the map! Value = " << (*iter) << " Index = " << index << endl;
    #endif
    ++index; ++iter;
  }
  vector<int >::iterator pixelIter = nPixels.begin();



  #ifdef DEBUG_OUTPUT_FOR_getClusterSNR_ON_INT_VECTOR
        streamlog_out( MESSAGE2 ) << "[getClusterSNR( std::vector<int > nPixels )] SORTED CHARGES:" << endl;
        #ifdef USE_MULTIMAP
            multimap<float, int >::reverse_iterator mapIterTmp = clusterSignalMap.rbegin(); //walking down descending!
        #else
            map<float, int >::reverse_iterator mapIterTmp = clusterSignalMap.rbegin(); //walking down descending!
        #endif
        while ( mapIterTmp != clusterSignalMap.rend() )
        {
            streamlog_out( MESSAGE2 ) << "Value:" << mapIterTmp->first << ", Index:" << mapIterTmp->second << endl;
            ++mapIterTmp;
        }
  #endif


  vector<float > snr;
  pixelIter = nPixels.begin();
  while ( pixelIter != nPixels.end() )
  {
    if ( /*(size_t) (*pixelIter) >= _trackerData->getChargeValues().size()*/ (*pixelIter) >= 7  ) //! 7 because this is the maximum amount of
    {                                                                                             //! pixels for us to take care of in a bricked
      snr.push_back( getClusterSNR() );                                                           //! cluster with one surrounding layer of
    }                                                                                             //! neighbor pixels (represented by a 3x3 fixed frame)
    else
    {
#ifdef USE_MULTIMAP
      multimap<float, int >::reverse_iterator mapIter = clusterSignalMap.rbegin(); //walking down descending!
#else
      map<float, int >::reverse_iterator mapIter = clusterSignalMap.rbegin(); //walking down descending!
#endif
      float signal = 0;
      float noise2 = 0;
      int   iPixel = 0;
      while ( (iPixel < (*pixelIter)) && (mapIter != clusterSignalMap.rend()) ) //NOTE: mapIter should never come close to rend()!!
      {                                                                         //rend() = 'begin()-1'
        signal += mapIter->first;                                               //note the difference between rbegin and begin
        noise2 += pow( _noiseValues[ mapIter->second], 2 );

        //! check, if we accidentally picked the fake strong negative values! (should not happen though!)
        if ( (mapIter->first) < -10.0 )
        {
            streamlog_out( ERROR4 ) << endl << "[getClusterSNR( std::vector<int > nPixels )] Might have reached the outsider pixels!!!" << endl;
            streamlog_out( ERROR4 ) << "  Just found a signal value < -10.0 !!! VERY BAD !!!" << endl;
            streamlog_out( ERROR4 ) << "  Was looking for SNR with " << (*pixelIter) << " MSP this time." << endl;
            streamlog_out( ERROR4 ) << "  Was working on the " << (iPixel)+1 << "th MSP. With index " << (mapIter->second) << endl;
            streamlog_out( ERROR4 ) << "  Signal for this pixel is " << mapIter->first << ", Signal values are:" << endl;
            outputBricked3x3MatrixFromVector(vectorCopy);
            int xSeed, ySeed;
            getSeedCoord(xSeed, ySeed);
            streamlog_out( ERROR4 ) << "  Seed Coord: x=" << xSeed << ", y=" << ySeed << "." << endl;
            streamlog_out( ERROR4 ) << endl << endl;
        }
        if ( _noiseValues[ mapIter->second ] == 0.0 )
        {
            streamlog_out( ERROR4 ) << endl << "[getClusterSNR( std::vector<int > nPixels )] Possible error in noise map or pixel status!" << endl;
            streamlog_out( ERROR4 ) << "  Just found a noise value == 0.0 !" << endl;
            streamlog_out( ERROR4 ) << "  Was looking for SNR with " << (*pixelIter) << " MSP this time." << endl;
            streamlog_out( ERROR4 ) << "  Was working on the " << (iPixel)+1 << "th MSP. With index " << (mapIter->second) << endl;
            streamlog_out( ERROR4 ) << "  Noise for this pixel is " << _noiseValues[ mapIter->second ] << ", Noise values are:" << endl;
            outputBricked3x3MatrixFromVector(_noiseValues);
            int xSeed, ySeed;
            getSeedCoord(xSeed, ySeed);
            streamlog_out( ERROR4 ) << "  Seed Coord: x=" << xSeed << ", y=" << ySeed << "." << endl;
            streamlog_out( ERROR4 ) << endl << endl;
        }


        ++mapIter; ++iPixel;

        //! check, if we accidentally picked the fake strong negative values! (should not happen though!)
        if ( mapIter==clusterSignalMap.rend() ) //this means that the whole map has been worked on. the map should always be 9 in size, and it should only be used up to 7 pixels!
        {
            streamlog_out( ERROR4 ) << endl << "[getClusterSNR( std::vector<int > nPixels )] !!! reached the outsider pixels we wanted to ignore !!! VERY BAD !!!" << endl;
        }

      }
      if ( noise2 == 0 ) snr.push_back(0.);
      else snr.push_back( signal / sqrt( noise2 ) );
    }
    ++pixelIter;
  }
  return snr;
}

float EUTelBrickedClusterImpl::getClusterSNR(int xSize, int ySize) const {

    if (!(xSize==1 && ySize==1)) //NOTE adjust this if you implement a variable size
    {
        //throw IncompatibleDataSetException("[EUTelBrickedClusterImpl::getClusterSNR(int xSize, int ySize)] This is not really applicable within a bricked pixel structure.");
        streamlog_out( WARNING2 ) << "[getClusterSNR( int x, int y )] Not applicable for a brickedCluster. Will return whole cluster's SNR!" << endl;
        return getClusterSNR();
    }
    return getSeedSNR();
    //return getClusterSNR();
}

float EUTelBrickedClusterImpl::getClusterCharge(int xSize, int ySize) const {

    if (!(xSize==1 && ySize==1)) //NOTE adjust this if you implement a variable size
    {
        //throw IncompatibleDataSetException("[EUTelBrickedClusterImpl::getClusterCharge(int xSize, int ySize)] This is not really applicable within a bricked pixel structure.");
        streamlog_out( WARNING2 ) << "[getClusterCharge( int x, int y )] Not applicable for a brickedCluster. Will return whole cluster's charge!" << endl;
        getTotalCharge();
    }
    return getSeedCharge();
    //return getClusterCharge();
}

void EUTelBrickedClusterImpl::print(std::ostream& os ) const {

    int xSeed, ySeed, xSize, ySize;
    float xShift, yShift, xShift2, yShift2, xShift3, yShift3; //xShift3x3, yShift3x3
    ClusterQuality quality = getClusterQuality();
    getClusterSize(xSize,ySize);
    if (xSize != 3 || ySize != 3)
    {
        streamlog_out( ERROR2 ) << "[EUTelBrickedClusterImpl::print(std::ostream& os )] Wrong Cluster Size!!!" << endl;
        return; //TODO ADJUST THIS IF YOU IMPLEMENT A VARIABLE SIZE
    }

    getSeedCoord(xSeed, ySeed);
    getCenterOfGravityShift(xShift, yShift);
    getCenterOfGravityShift(xShift2, yShift2, 2);
    getCenterOfGravityShift(xShift3, yShift3, 3);
    //getCenterOfGravityShift(xShift3x3, yShift3x3, 3, 3);

    float noise = 0., SNR = 0., SNR2 = 0., SNR3 = 0. ;//SNR3x3 = 0.
    if ( _noiseSetSwitch )
    {
        noise  = getClusterNoise();
        SNR    = getClusterSNR();
        SNR2   = getClusterSNR(2);
        SNR3   = getClusterSNR(3);
        //SNR3x3 = getClusterSNR(3,3);
    }

    int bigspacer = 23;

    os  <<  setw(bigspacer) << setiosflags(ios::left) << "Bricked pixel cluster "<< "( one surrounding layer of neighbour pixels )\n"
                                                                                    //NOTE adjust this if you implement a variable size
        <<  setw(bigspacer) <<  "Cluster quality: " << quality << "\n"
        <<  setw(bigspacer) <<  "Seed coords: " << "x=" << xSeed << ", y=" << ySeed << "\n"
        <<  setw(bigspacer) <<  "Seed charge: " << getSeedCharge() << " in (" << xSeed << ", " << ySeed << ")\n"
        <<  setw(bigspacer) <<  "Total charge: " << getTotalCharge() << "\n"
        <<  setw(bigspacer) <<  "2MSP charge: " << getClusterCharge(2) << "\n"
        <<  setw(bigspacer) <<  "3MSP charge: " << getClusterCharge(3) << "\n"
        //<<  setw(bigspacer) <<  "3x3-Charge       " << getClusterCharge(3,3) << "\n"
        <<  setw(bigspacer) <<  "CoG shift FULL: " << "(" << xShift  << ", " << yShift  << ")\n"
        <<  setw(bigspacer) <<  "CoG shift 2MSP: " << "(" << xShift2 << ", " << yShift2 << ")\n"
        <<  setw(bigspacer) <<  "CoG shift 3MSP: " << "(" << xShift3 << ", " << yShift3 << ")\n"
        ;//<<  setw(bigspacer) <<  "CoG(3x3) shift " << "(" << xShift3x3 << ", " << yShift3x3 << ")\n" ;


    if ( _noiseSetSwitch )
    {
       os << setw(bigspacer)  <<  "Cluster noise: " << noise << "\n"
          << setw(bigspacer)  <<  "Cluster SNR: " << SNR  << "\n"
          << setw(bigspacer)  <<  "2MSP SNR: " << SNR2 << "\n"
          << setw(bigspacer)  <<  "3MSP SNR: " << SNR3 << "\n"
          ;//<< setw(bigspacer)  <<  "Cluster SNR(3x3) " << SNR3x3 << "\n";
    }
    else
    {
         os << setw(bigspacer)  <<  "(Cluster noise not set)\n";
    }

    os << "--- Signal Values ---" << endl;
    FloatVec vectorCopy(_trackerData->getChargeValues());
    outputBricked3x3MatrixFromVector( static_cast< std::vector< float>& > (vectorCopy) );

}


void EUTelBrickedClusterImpl::setOutsiderValuesInVectorInterpretedAsBrickedMatrix(std::vector<float>& v, float val) const
{
    //NOTE Fix this if you switch the implemenatation to a variable size:
    //     For now this is only implemented for 3x3...
    //     ...and _even rows_ being skewed to the left (aka their coordinates have to become 0.5 smaller!)
    //     If _even rows_ are skewed to the right, then this method has to remove
    //     pixels on the other sides respectively compared to what it does now.

    //     Anyway the code down there can be used for such a variable implementation.
    //     It is very well suited to be adapted a little bit more.
    //     Ee should assume a square instead of a rectangle as well because everything would become
    //     even more complicated otherwise! =/

    //!General approach would be:
    //- with each row further away from the seed pixel row consider one pixel less in x direction!
    //- only in odd rows away from seed pixel (  |yRowCurrent-yRowSeedpixel|=2n+1  ) you have
    //  to consider on which side to drop a pixel.
    //- in even rows away from seed pixel there is an uneven number of pixels to use
    //  and they are centered in xDirection vertically above/below the seed pixel
    //- ((a little drawing will help a lot!))

    //!For now:
    // This is only implemented for 3x3 size...
    // ...and even rows being skewed to the left (aka their coordinates have to become 0.5 smaller!)


    //check sizes
    int xSize, ySize;
    getClusterSize(xSize, ySize);
    if (!(xSize==3 && ySize==3))
    {
        streamlog_out( WARNING4 ) << "EUTelBrickedClusterImpl::setOutsiderValuesInVectorInterpretedAsBrickedMatrix(FloatVec& v):" << endl;
        streamlog_out( WARNING4 ) << " BRICKED PIXEL FIXED FRAME CLUSTER SIZE MUST BE 3x3!!!" << endl;
        streamlog_out( WARNING4 ) << " BUT IT IS " << xSize << "x" << ySize << "!!!" << endl;
    }
    if (! ( static_cast< size_t >(xSize*ySize) == v.size() ) )
    {
        streamlog_out( WARNING4 ) << "EUTelBrickedClusterImpl::setOutsiderValuesInVectorInterpretedAsBrickedMatrix(FloatVec& v):" << endl;
        streamlog_out( WARNING4 ) << " BRICKED PIXEL FIXED FRAME CLUSTER SIZE DOES NOT MATCH VECTOR SIZE!!!" << endl;
        streamlog_out( WARNING4 ) << " FIXED FRAME: " << xSize << "x" << ySize << ", VECTOR: " << v.size() << "." << endl;
    }

    //start removing unwanted pixels
    int seedX, seedY;
    getSeedCoord(seedX, seedY);
    bool bSeedRowIsEven = false;
    if (seedY % 2 == 0) bSeedRowIsEven = true; //seed pixel's row is even

    int startZerosX, endZerosX;

    for (int y= (-1)*ySize/2; y<=ySize/2; ++y)
    {
        if (y!=0) //only if we are not looking at the middle row: write zeros to the two pixels we don't want.
        {              //this will be a little bit more complicated in a variable approach
                       //(have to take care of rows with even and odd distances of the seed pixel row, etc.)
            if (bSeedRowIsEven) //seed pixel's row is even:
            {                   //delete one pixel in this row, on the right side
                startZerosX  = (xSize / 2);
                endZerosX    = startZerosX; //this will be a little bit more complicated in a variable approach
            }
            else //seed pixel's row is odd
            {    //delete one pixel in this row, on the left side
                startZerosX  = -1 * (xSize / 2);
                endZerosX    = startZerosX; //this will be a little bit more complicated in a variable approach
            }

            for (int x= startZerosX; x<=endZerosX; ++x)
            {
                // the pixel index is:  [ xSize*(y + ySize/2) + (x + xSize/2) ];
                v[xSize*(y + ySize/2) + (x + xSize/2)] = val;
            }
        }
    }


}


void EUTelBrickedClusterImpl::debugOutput() const
{

    streamlog_out( MESSAGE2 ) << endl;
    streamlog_out( MESSAGE2 ) << "  = Cluster Info: " << endl;

    //!GENERAL INFO
    int xSize, ySize;
    getClusterSize(xSize,ySize);
    streamlog_out( MESSAGE2 ) << "    Size: x=" << xSize << ", y= " << ySize << endl;

    int xSeed, ySeed;
    getSeedCoord(xSeed, ySeed);
    streamlog_out( MESSAGE2 ) << "    Seed: x=" << xSeed << ", y= " << ySeed << endl;


    //!CHARGE MATRIX
    FloatVec vectorCopy(_trackerData->getChargeValues());
    //setOutsiderValuesInVectorInterpretedAsBrickedMatrix( static_cast< std::vector< float>& > (vectorCopy), 777.0f );

    streamlog_out( MESSAGE2 ) << endl;
    streamlog_out( MESSAGE2 ) << "  = Charge matrix:" << endl;
    outputBricked3x3MatrixFromVector( static_cast< std::vector< float>& > (vectorCopy) );
    streamlog_out( MESSAGE2 ) << "    Total Charge: " << getTotalCharge() << endl;

    std::vector<int> numbersTwoThreeSixSeven;
    numbersTwoThreeSixSeven.push_back(2);
    numbersTwoThreeSixSeven.push_back(3);
    numbersTwoThreeSixSeven.push_back(6);
    numbersTwoThreeSixSeven.push_back(7);
    std::vector<float> fourMspCharges = getClusterCharge( numbersTwoThreeSixSeven );
    streamlog_out( MESSAGE2 ) << "    MSP2  Charge: " << getClusterCharge(2) << endl;
    streamlog_out( MESSAGE2 ) << "    MSP3  Charge: " << getClusterCharge(3) << endl;
    streamlog_out( MESSAGE2 ) << "    MSP6  Charge: " << getClusterCharge(6) << endl;
    streamlog_out( MESSAGE2 ) << "    MSP7  Charge: " << getClusterCharge(7) << endl;
    streamlog_out( MESSAGE2 ) << "    MSP2 ChargeV: " << fourMspCharges[0] << endl;
    streamlog_out( MESSAGE2 ) << "    MSP3 ChargeV: " << fourMspCharges[1] << endl;
    streamlog_out( MESSAGE2 ) << "    MSP6 ChargeV: " << fourMspCharges[2] << endl;
    streamlog_out( MESSAGE2 ) << "    MSP7 ChargeV: " << fourMspCharges[3] << endl;


    //!COG SHIFT
    streamlog_out( MESSAGE2 ) << endl;
    float xShift, yShift;
    getCenterOfGravityShift(xShift, yShift);
    streamlog_out( MESSAGE2 ) << "  = CoG Shift Global Full: x=" << xShift << ", y= " << yShift << endl;

    getCenterOfGravityShift(xShift, yShift, 2);
    streamlog_out( MESSAGE2 ) << "    CoG Shift Global 2MSP: x=" << xShift << ", y= " << yShift << endl;

    getCenterOfGravityShift(xShift, yShift, 3);
    streamlog_out( MESSAGE2 ) << "    CoG Shift Global 3MSP: x=" << xShift << ", y= " << yShift << endl;

    getCenterOfGravityShiftWithOutGlobalSeedCoordinateCorrection(xShift, yShift);
    streamlog_out( MESSAGE2 ) << "    CoG Shift ETA    FULL: x=" << xShift << ", y= " << yShift << endl;


    //!NOISE
    float noise = 0., SNR = 0., SNR2 = 0., SNR3 = 0., SNR6 = 0., SNR7 = 0.;
    if ( _noiseSetSwitch )
    {
        streamlog_out( MESSAGE2 ) << endl;
        streamlog_out( MESSAGE2 ) << "  = Noise matrix:" << endl;
        outputBricked3x3MatrixFromVector(_noiseValues);


        noise  = getClusterNoise();
        SNR    = getClusterSNR();
        SNR2   = getClusterSNR(2);
        SNR3   = getClusterSNR(3);
        SNR6   = getClusterSNR(6);
        SNR7   = getClusterSNR(7);
        std::vector<float> fourMspSNRs = getClusterSNR( numbersTwoThreeSixSeven );
        streamlog_out( MESSAGE2 ) << "  = ClusterNoise    : " << noise << endl;
        streamlog_out( MESSAGE2 ) << "    ClusterSNR FULL : " << SNR << endl;
        streamlog_out( MESSAGE2 ) << "    ClusterSNR 2MSP : " << SNR2 << endl;
        streamlog_out( MESSAGE2 ) << "    ClusterSNR 3MSP : " << SNR3 << endl;
        streamlog_out( MESSAGE2 ) << "    ClusterSNR 6MSP : " << SNR6 << endl;
        streamlog_out( MESSAGE2 ) << "    ClusterSNR 7MSP : " << SNR7 << endl;
        streamlog_out( MESSAGE2 ) << "    ClusterSNR 2MSPv: " << fourMspSNRs[0] << endl;
        streamlog_out( MESSAGE2 ) << "    ClusterSNR 3MSPv: " << fourMspSNRs[1] << endl;
        streamlog_out( MESSAGE2 ) << "    ClusterSNR 6MSPv: " << fourMspSNRs[2] << endl;
        streamlog_out( MESSAGE2 ) << "    ClusterSNR 7MSPv: " << fourMspSNRs[3] << endl;
    }
    else
    {
         streamlog_out( WARNING4 ) << "    Noise not set!" << endl;
    }
    streamlog_out( MESSAGE2 ) << endl;
}

void EUTelBrickedClusterImpl::outputBricked3x3MatrixFromVector(const std::vector<float>& v) const
{

    /**
    *
    * Short explanation on how the output works:
    * a) seedRow is even:
    *    .even rows are substracted 0.5 from their x pixel coordinate number
    *    .so odd rows are the refernce here
    *    .so the seedRow (and other even rows if the cluster is bigger than just one surroundling layer, aka 3x3)
    *     is virtually shifted left by 0.5 because that is the case in reality as well.
    *    .for this output we can only shift stuff right tough. so we shift the other (odd) rows right
    *     and leave the rest (including the pixel coordinate "header" of our small table in place)
    *
    * b) seedRow is odd:
    *    .even Rows are substracted 0.5 from their x pixel coordinate number as well
    *    .so odd rows are the refernce here as well
    *    .so the (other!) even rows next to the seed row are
    *     virtually shifted left by 0.5 because that is the case in reality as well.
    *    .for this output we can only shift stuff right tough. so we shift the seedRow right
    *     and leave the rest (including the pixel coordinate "header" of our small table in place)
    */

    if (v.size() != 9)
    {
        streamlog_out( ERROR2 ) << "    [outputBricked3x3MatrixFromVector(std::vector<float>& v)] Wrong vectorSize!!!" << endl;
        return;
    }

    int xSeed, ySeed;
    getSeedCoord(xSeed, ySeed);

    streamlog_out( ERROR4 ) << "-----------------------------------------------------" << endl;
    streamlog_out( ERROR4 ) << "| x=      |" //11 long
                            << "     " << setw(3) << xSeed-1 << "     |"
                            << "     " << setw(3) << xSeed-0 << "     |"
                            << "     " << setw(3) << xSeed+1 << "     |"
                            << endl;

    streamlog_out( ERROR4 ) << "-----------------------------------------------------" << endl;
    if (ySeed % 2 == 0)
    {
        streamlog_out( ERROR4 ) << "| " << "y= " << setw(4) << ySeed-1 << " |"  //2+3+4+2=11 long
                                << "      "  //! EXTRA SPACE HERE (THE SKEW)
                                << setw(13) << v[0] << ","
                                << setw(13) << v[1] << ","
                                << "(" << setw(11) << v[2] << ")"
                                << endl << "|" << endl;

        streamlog_out( ERROR4 ) << "| " << "y= " << setw(4) << ySeed-0 << " |"
                                << setw(13) << v[3] << ","
                                << setw(13) << v[4] << ","
                                << setw(13) << v[5]
                                << endl << "|" << endl;

        streamlog_out( ERROR4 ) << "| " << "y= " << setw(4) << ySeed+1 << " |"
                                << "      " //! EXTRA SPACE HERE (THE SKEW)
                                << setw(13) << v[6] << ","
                                << setw(13) << v[7] << ","
                                << "(" << setw(11) << v[8] << ")";
    }
    else
    {
        streamlog_out( ERROR4 ) << "| " << "y= " << setw(4) << ySeed-1 << " |"
                                << "(" << setw(11) << v[0] << "),"
                                << setw(13) << v[1] << ","
                                << setw(13) << v[2]
                                << endl << "|" << endl;

        streamlog_out( ERROR4 ) << "| " << "y= " << setw(4) << ySeed+0 << " |"
                                << "      "  //! EXTRA SPACE HERE (THE SKEW)
                                << setw(13) << v[3] << ","
                                << setw(13) << v[4] << ","
                                << setw(13) << v[5]
                                << endl << "|" << endl;

        streamlog_out( ERROR4 ) << "| " << "y= " << setw(4) << ySeed+1 << " |"
                                << "(" << setw(11) << v[6] << "),"
                                << setw(13) << v[7] << ","
                                << setw(13) << v[8];
    }
    streamlog_out( ERROR4 ) << endl;
    streamlog_out( ERROR4 ) << "-----------------------------------------------------" << endl;


    return;

}
