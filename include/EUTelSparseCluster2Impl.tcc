// Version: $Id$
// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELSPARSECLUSTERIMPL_TCC
#define EUTELSPARSECLUSTERIMPL_TCC


namespace eutelescope {

  template<class PixelType>
  EUTelSparseCluster2Impl<PixelType>::EUTelSparseCluster2Impl(IMPL::TrackerDataImpl * data) : 
    EUTelVirtualCluster(data),
    _nElement(0),
    _type(kUnknownPixelType),
    _trackerData()
  {

    std::auto_ptr<PixelType> pixel( new PixelType);
    _nElement       = pixel->getNoOfElements();
    _type           = pixel->getSparsePixelType();
    _trackerData    = data;


    // reset all the associated vectors
    _noiseValues.clear();
    _pixelVec.clear();

    // reset all the booleans
    _noiseSetSwitch    = false;
    _isPositionSorted  = false;
    _isSignalSorted    = false;
    _isOriginalOrder   = true;
  }


  template<class PixelType>
  unsigned int EUTelSparseCluster2Impl<PixelType>::size()  const {
    return _pixelVec.size();
  }

  template<class PixelType>
  void EUTelSparseCluster2Impl<PixelType>::setNoiseValues(std::vector<float > noiseValues) {
    if ( noiseValues.size() != size() ) {
      _noiseSetSwitch = false;
      throw IncompatibleDataSetException("The noiseValues size is different from the number of pixel in the cluster");
    }

    _noiseValues    = noiseValues;
    _noiseSetSwitch = true;

  }

  template<class PixelType>
  std::vector<float > EUTelSparseCluster2Impl<PixelType>::getNoiseValues() const {
    if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");
    return _noiseValues;
  }

  template<class PixelType>
  void EUTelSparseCluster2Impl<PixelType>::getSeedCoord(int& xSeed, int& ySeed) const {

    std::vector<PixelType >::const_iterator seedIterator = max_element( _pixelVec.begin(), _pixelVec.end(), SmallerSignal<PixelType>() );
    xSeed = (*seedIterator)->getXCoord();
    ySeed = (*seedIterator)->getYCoord();


  }

  template<class PixelType>
  float EUTelSparseCluster2Impl<PixelType>::getTotalCharge() const {

    std::vector<PixelType >::const_iterator iter = _pixelVec.begin();
    float charge = 0 ;
    while ( iter != _pixelVec.end() ) {
      charge += (*iter).getSignal();
      ++iter;
    }

  }

  template<class PixelType>
  float EUTelSparseCluster2Impl<PixelType>::getSeedCharge() const {
    std::vector<PixelType >::const_iterator seedIterator = max_element( _pixelVec.begin(), _pixelVec.end(), Smallersignal<PixelType>() );
    return (*seedIterator)->getSignal();
  }


  template<class PixelType>
  void EUTelSparseCluster2Impl<PixelType>::getCenterOfGravityShift(float& xCoG, float& yCoG) const {
    if ( size() == 1 ) {
      xCoG = 0;
      yCoG = 0;
      return;
    }

    float normalization = 0;
    float tempX = 0;
    float tempY = 0;

    int xSeed, ySeed;
    getSeedCoord(xSeed, ySeed);

    std::vector<PixelType >::const_iterator pixelIter = _pixelVec.begin();
    while ( pixelIter != _pixelVec.end() ) {
      tempX         += ( (*pixelIter)->getSignal() * ( (*pixelIter)->getXCoord() - xSeed ) );
      tempY         += ( (*pixelIter)->getSignal() * ( (*pixelIter)->getYCoord() - ySeed ) );
      normalzation  += ( (*pixelIter)->getSignal() );
      ++pixelIter;
    }

    if ( normalization != 0 ) {
      xCoG = tempX / normalization;
      yCoG = tempY / normalization;
    } else {
      xCoG = 0.;
      yCoG = 0.;
    }

  }

  template<class PixelType>
  void EUTelSparseCluster2Impl<PixelType>::getCenterOfGravityShift(float& xCoG, float& yCoG,
                                                                   int xSize, int ySize) const {

    if ( size() == 1 ) {
      xCoG = 0;
      yCoG = 0;
      return;
    }

    int xSeed, ySeed;
    getSeedCoord(xSeed, ySeed);

    float normalization = 0,  tempX = 0, tempY = 0;

    std::vector<float >::iterator signalIter = _signalVec.begin();
    std::vector<int   >::iterator xIter      = _xCoordVec.begin();
    std::vector<int   >::iterator yIter      = _yCoordVec.begin();

    std::vector<PixelType >::const_iterator pixelIter = _pixelVec.begin();
    while ( pixelIter != _pixelVec.end() ) {
      if ( ( abs( (*pixelIter)->getXCoord() - xSeed ) < ( xSize / 2 ) ) &&
           ( abs( (*pixelIter)->getYCoord() - ySeed ) < ( ySize / 2 ) ) ) {
        tempX         += ( (*pixelIter)->getSignal() * ( (*pixelIter)->getXCoord() - xSeed ) );
        tempY         += ( (*pixelIter)->getSignal() * ( (*pixelIter)->getYCoord() - ySeed ) );
        normalzation  += ( (*pixelIter)->getSignal() );
      }
      ++pixelIter;
    }

    if ( normalization != 0 ) {
      xCoG = tempX / normalization;
      yCoG = tempY / normalization;
    } else {
      xCoG = 0.;
      yCoG = 0.;
    }

  }

  template<class PixelType>
  void EUTelSparseCluster2Impl<PixelType>::getCenterOfGravityShift(float& xCoG, float& yCoG,  int n) const {

    if ( size() == 1 ) {
      xCoG = 0;
      yCoG = 0;
      return;
    }

    if ( ! _isSignalSorted ) sortBySignal();

    float normalization = 0,  tempX = 0, tempY = 0;
    int xSeed, ySeed;
    getSeedCoord(xSeed, ySeed);

    std::vector<PixelType >::const_iterator pixelIter = _pixelVec.begin();
    while ( pixelIter != _pixelVec.begin() + n ) {
      tempX         += ( (*pixelIter)->getSignal() * ( (*pixelIter)->getXCoord() - xSeed ) );
      tempY         += ( (*pixelIter)->getSignal() * ( (*pixelIter)->getYCoord() - ySeed ) );
      normalization += (*pixelIter)->getSignal();
      ++pixelIter;
    }

    if ( normalization != 0 ) {
      xCoG = tempX / normalization;
      yCoG = tempY / normalization;
    } else {
      xCoG = 0.;
      yCoG = 0.;
    }

  }

  template<class PixelType>
  void EUTelSparseCluster2Impl<PixelType>::getCenterOfGravity(float&  xCoG, float& yCoG) const {
    getCenterOfGravityShift(xCoG, yCoG);

    int xSeed, ySeed;
    getSeedCoord(xSeed, ySeed);

    xCoG += static_cast<float > ( xSeed );
    yCoG += static_cast<float > ( ySeed );

  }


  template<class PixelType>
  void EUTelSparseCluster2Impl<PixelType>::getClusterSize(int& xSize, int& ySize) const {

    // also the new implementation is very expensive in terms of CPU
    // cycles... so, I'm keeping the old one also because this method
    // is only seldom used

    int xMin = std::numeric_limits<int>::max(), yMin = std::numeric_limits<int>::max();
    int xMax = std::numeric_limits<int>::min(), yMax = std::numeric_limits<int>::min();
    PixelType * pixel = new PixelType;
    for ( unsigned int index = 0; index < size() ; index++ ) {
      getSparsePixelAt( index , pixel);
      short xCur = pixel->getXCoord();
      short yCur = pixel->getYCoord();
      if ( xCur < xMin ) xMin = xCur;
      if ( xCur > xMax ) xMax = xCur;
      if ( yCur < yMin ) yMin = yCur;
      if ( yCur > yMax ) yMax = yCur;
    }
    xSize = abs( xMax - xMin) + 1;
    ySize = abs( yMax - yMin) + 1;

    delete pixel;
  }


  template<class PixelType>
  void EUTelSparseCluster2Impl<PixelType>::getCenterCoord(int& xCenter, int& yCenter) const {
    // by definition, the cluster center is the pixel containing the
    // charge center of gravity.
    float xCoG, yCoG;
    getCenterOfGravity(xCoG, yCoG);

    // unfortunately I couldn't found a better way to approximate a
    // float number to the closest integer. It seems that the standard
    // libraries are missing this simple tools and this is the only
    // hack I found!
    xCenter = static_cast<int> ( floor( xCoG + 0.5 ) );
    yCenter = static_cast<int> ( floor( yCoG + 0.5 ) );

  }

  template<class PixelType>
  float EUTelSparseCluster2Impl<PixelType>::getDistance(EUTelVirtualCluster * clu) const {
    // get the cluster center for this
    int xThisCenter,  yThisCenter;
    this->getCenterCoord(xThisCenter, yThisCenter);

    // get the cluster center for the other
    int xOtherCenter, yOtherCenter;
    clu->getCenterCoord(xOtherCenter, yOtherCenter);

    return sqrt( pow( static_cast<double> ( xThisCenter - xOtherCenter), 2) +
                 pow( static_cast<double> ( yThisCenter - yOtherCenter), 2) );
  }

  template<class PixelType>
  float EUTelSparseCluster2Impl<PixelType>::getExternalRadius() const {
    int xSize, ySize;
    getClusterSize(xSize, ySize);
    return 0.5 * std::max(xSize, ySize);
  }

  template<class PixelType>
  float EUTelSparseCluster2Impl<PixelType>::getClusterCharge(int nPixel) const {

    // if the cluster is smaller than nPixel, return the total charge
    // w/o any further calculation!
    if ( nPixel >= (signed) size() ) return getTotalCharge() ;

    if ( ! _isSignalSorted ) sortBySignal();

    float charge = 0;
    std::vector<PixelType >::const_iterator pixelIter = _pixelVec.begin();
    while ( pixelIter != _pixelVec.begin() + n ) {
      charge += (*pixelIter)->getSignal();
      ++pixelIter;
    }

    return charge;
  }

  e  template<class PixelType>
  std::vector<float > EUTelSparseCluster2Impl<PixelType>::getClusterCharge(std::vector<int > nPixels) const {

    std::vector<float > clusterSignal;

    if ( ! _isSignalSorted ) sortBySignal();

    std::vector<PixelType >::const_iterator iter;

    for ( unsigned int i = 0 ; i < nPixels.size(); i++ ) {
      if ( nPixel[i] >= size() ) {
        clusterSignal.push_back( getTotalCharge() );
      } else {
        iter = _pixelVec.begin();
        float charge = 0;
        while ( iter != _pixelVec.begin() + nPixels[i] ) {
          charge += (*iter)->getSignal();
          ++iter;
        }
        clusterSignal.push_back( charge );
      }
    }
    return clusterSignal;
  }



  template<class PixelType>
  float EUTelSparseCluster2Impl<PixelType>::getClusterCharge(int xSize, int ySize) const {

    float charge = 0;

    int xSeed, ySeed;
    getSeedCoord(xSeed, ySeed);

    std::vector<PixelType >::const_iterator iter = _pixelVec.begin();
    while ( iter != _pixelVec.end() ) {
      xPixel = (int)
        yPixel = (int) pixel->getYCoord();

      if ( ( abs( xSeed - (*iter)->getXCoord() ) <= ( xSize / 2 ) ) &&
           ( abs( ySeed - (*iter)->getYCoord() ) <= ( ySize / 2 ) ) ) {
        charge += (*iter)->getSignal();
      }
    }
    return charge;
  }

  template<class PixelType>
  void EUTelSparseCluster2Impl<PixelType>::setClusterQuality(eutelescope::ClusterQuality quality)  {
    lcio::long64 cell0 = static_cast<lcio::long64> ( _trackerData->getCellID0() );
    int rhs = 18;
    lcio::long64 emptyMask     = ~( 0x1F << rhs );
    lcio::long64 maskedQuality = ( (static_cast<int> ( quality ) & 0x1F ) << rhs );

    // first apply an empty mask for the quality bit ranges, that is
    // to say reset the quality bit range but keep all the rest.
    cell0 = cell0 & emptyMask;

    // now apply the maskedQuality.
    cell0 = cell0 | maskedQuality;

    _trackerData->setCellID0(cell0);

  }

  template<class PixelType>
  float EUTelSparseCluster2Impl<PixelType>::getClusterNoise() const {

    if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");

    float squaredSum = 0;
    std::vector<float >::const_iterator iter = _noiseValues.begin();
    while ( iter != _noiseValues.end() ) {
      squaredSum += pow( (*iter), 2 );
      ++iter;
    }
    return sqrt( squaredSum );

  }

  template<class PixelType>
  float EUTelSparseCluster2Impl<PixelType>::getClusterSNR() const {
    if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");
    float clusterNoise = getClusterNoise();
    if ( clusterNoise == 0 )  return 0.;
    float clusterSignal = getTotalCharge();
    return clusterSignal / clusterNoise;
  }

  template<class PixelType>
  float EUTelSparseCluster2Impl<PixelType>::getSeedSNR() const {

    if ( ! _isOriginalOrder ) restoreOriginalOrder();

    std::vector<PixelType >::const_iterator seedIterator = max_element( _pixelVec.begin(), _pixelVec.end(), Smallersignal<PixelType>() );
    float signal = (*seedIterator)->getSignal();
    float noise  = _noiseValues[ seedIterator - _pixelVec.begin()];
    if ( noise == 0 ) return 0;
    return signal / noise;
  }


  template<class PixelType>
  float EUTelSparseCluster2Impl<PixelType>::getClusterSNR(int nPixel) const {
    // still I don't have a better idea for the clusterSNR(int
    // nPixel);

    if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");
    if ( ( unsigned ) nPixel >= size() )    return getClusterSNR();


    PixelType * pixel = new PixelType;

    // prepare a vector containing all the pixel signals
    std::vector<float > signalVec;

    for ( unsigned int iPixel = 0 ; iPixel < size() ; iPixel++ ) {
      getSparsePixelAt( iPixel, pixel );
      signalVec.push_back( pixel->getSignal() );
    }

    delete pixel;

    std::map<int , float > highSignalPixel;
    int iPixel = 0;
    while ( iPixel != nPixel ) {
      float maxSignal = -1 * std::numeric_limits<float >::max();
      int   maxIndex  = 0;
      int   index     = 0;
      std::vector<float >::iterator maxIter;
      std::vector<float >::iterator iter = signalVec.begin();

      while ( iter != signalVec.end() ) {
        if ( *iter > maxSignal ) {
          maxSignal = (*iter);
          maxIndex  = index;
          maxIter   = iter;
        }
        ++index;
        ++iter;
      }
      highSignalPixel[maxIndex] = maxSignal;
      (*maxIter) = -1 *  std::numeric_limits<float >::max();
      ++iPixel;
    }

    float signal = 0, noise2 = 0;
    std::map<int, float >::iterator mapIter = highSignalPixel.begin();
    while ( mapIter != highSignalPixel.end() ) {
      signal += mapIter->second;
      noise2 += pow( _noiseValues[mapIter->first], 2 );
      ++mapIter;
    }
    if ( noise2 == 0 ) return 0;
    return signal / sqrt( noise2 );

  }

  template<class PixelType>
  float EUTelSparseCluster2Impl<PixelType>::getClusterSNR(int xSize, int ySize) const {
    if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");

    int xSeed, ySeed;
    getSeedCoord(xSeed, ySeed) ;

    PixelType * pixel = new PixelType;
    int xPixel, yPixel;

    float charge = 0, noise2 = 0;

    for (unsigned int iPixel = 0;  iPixel < size() ; iPixel++ ) {
      getSparsePixelAt( iPixel, pixel );
      xPixel = (int) pixel->getXCoord();
      yPixel = (int) pixel->getYCoord();
      if ( ( abs( xSeed - xPixel ) <= ( xSize / 2 ) ) &&
           ( abs( ySeed - yPixel ) <= ( ySize / 2 ) ) ) {
        charge += pixel->getSignal();
        noise2 += pow( _noiseValues[iPixel] , 2 );
      }
    }
    delete pixel;
    if ( noise2 != 0 ) return charge / sqrt( noise2 );
    else return 0.;
  }


  template<class PixelType>
  std::vector<float > EUTelSparseCluster2Impl<PixelType>::getClusterSNR( std::vector<int > nPixels ) const {

    if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");

    // put all pixel values into a map
    PixelType * pixel = new PixelType;
    std::map<float, int> clusterSignalMap;

    for ( unsigned int iPixel = 0 ; iPixel < size(); iPixel++ ) {
      getSparsePixelAt( iPixel, pixel );
      clusterSignalMap.insert( make_pair( pixel->getSignal() , iPixel ) );
    }

    delete pixel;

    std::vector<int >::iterator pixelIter = nPixels.begin();
    std::vector<float > snr;

    while ( pixelIter != nPixels.end() ) {
      std::map<float, int >::reverse_iterator mapIter = clusterSignalMap.rbegin();
      float signal = 0;
      float noise2 = 0;
      int   iPixel = 0;
      while ( ( iPixel < (*pixelIter)) && (mapIter != clusterSignalMap.rend() ) ) {
        signal += mapIter->first;
        noise2 += pow( _noiseValues[ mapIter->second ], 2 );
        ++mapIter;
        ++iPixel;
      }
      if ( noise2 == 0 ) snr.push_back( 0. );
      else snr.push_back( signal / sqrt( noise2 ) );
      ++pixelIter;
    }
    return snr;
  }

  template<class PixelType>
  void EUTelSparseCluster2Impl<PixelType>::print(std::ostream& os) const {

    int   xSize, ySize, xSeed, ySeed, xCenter, yCenter;
    float xShift, yShift, xShift9, yShift9, xShift3x3, yShift3x3;
    ClusterQuality  quality = getClusterQuality();
    SparsePixelType type    = getSparsePixelType();
    getClusterSize(xSize, ySize);
    getSeedCoord(xSeed, ySeed);
    getCenterCoord(xCenter, yCenter);
    getCenterOfGravityShift(xShift, yShift);
    getCenterOfGravityShift(xShift9, yShift9, 9);
    getCenterOfGravityShift(xShift3x3, yShift3x3, 3, 3);

    float noise =0., SNR = 0., SNR9 = 0., SNR3x3 = 0. ;
    if ( _noiseSetSwitch ) {
      noise  = getClusterNoise();
      SNR    = getClusterSNR();
      SNR9   = getClusterSNR(9);
      SNR3x3 = getClusterSNR(3,3);
    }

    int bigspacer = 23;

    os   <<  std::setw(bigspacer) << std::setiosflags(std::ios::left) << "Sparse cluster made of " << type << " pixels" << std::endl
         <<  std::setw(bigspacer) << "Number of pixel " << size() << std::endl
         <<  std::setw(bigspacer) << "Cluster ID " << getClusterID() << " on detector " << getDetectorID() << std::endl
         <<  std::setw(bigspacer) << "Cluster size " << "(" << xSize << ", " << ySize << ")" << std::endl
         <<  std::setw(bigspacer) << "Cluster quality " << quality << std::endl
         <<  std::setw(bigspacer) << "Cluster total charge " << getTotalCharge() << std::endl
         <<  std::setw(bigspacer) << "Cluster charge (9) " << getClusterCharge(9) << std::endl
         <<  std::setw(bigspacer) << "Cluster charge (3x3) " << getClusterCharge(3,3) << std::endl
         <<  std::setw(bigspacer) << "Seed charge " << getSeedCharge() << " in (" << xSeed << ", " << ySeed << ")" << std::endl
         <<  std::setw(bigspacer) << "CoG shift " << "(" << xShift << ", " << yShift << ")" << std::endl
         <<  std::setw(bigspacer) <<  "CoG(9) shift " << "(" << xShift9 << ", " << yShift9 << ")" << std::endl
         <<  std::setw(bigspacer) <<  "CoG(3x3) shift " << "(" << xShift3x3 << ", " << yShift3x3 << ")" << std::endl;
    if ( _noiseSetSwitch ) {
      os << std::setw(bigspacer)  <<  "Cluster noise " << noise << std::endl
         << std::setw(bigspacer)  <<  "Cluster SNR " << SNR << std::endl
         << std::setw(bigspacer)  <<  "Cluster SNR(9) " << SNR9 << std::endl
         << std::setw(bigspacer)  <<  "Cluster SNR(3x3) " << SNR3x3 << std::endl;
    }

    int spacer = 14;
    for ( int i = 0; i < spacer - 1; i++ ) {
      os << "-";
    }
    os << std::endl;

    PixelType * pixel = new PixelType;
    for ( unsigned int iPixel = 0 ; iPixel < size() ; iPixel++ ) {
      getSparsePixelAt( iPixel, pixel );

      os << "Pixel number = " << iPixel << std::endl
         << ( * pixel ) << std::endl;
    }
    for ( int i = 0; i < spacer - 1; i++ ) {
      os << "-";
    }
    os << std::resetiosflags(std::ios::left) << std::endl;
    delete pixel;
  }
}

#endif
