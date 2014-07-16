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

#include <iostream>

namespace eutelescope {

  template<class PixelType>
  EUTelSparseClusterImpl<PixelType>::EUTelSparseClusterImpl(IMPL::TrackerDataImpl* data) : 
    EUTelVirtualCluster(data),
    _rawDataInterfacer(data),
    _nElement(0),
    _type(kUnknownPixelType),
    _noiseValues(),
    _noiseSetSwitch(false)
  {

    std::auto_ptr<PixelType> pixel( new PixelType);
    _nElement       = pixel->getNoOfElements();
    _type           = pixel->getSparsePixelType();
    _noiseValues.clear();
  }


  template<class PixelType>
  unsigned int EUTelSparseClusterImpl<PixelType>::size()  const {
    return _trackerData->getChargeValues().size() / _nElement;
  }

  template<class PixelType>
  void EUTelSparseClusterImpl<PixelType>::setNoiseValues(std::vector<float > noiseValues) {
    if ( noiseValues.size() != size() ) {
      _noiseSetSwitch = false;
      throw IncompatibleDataSetException("The noiseValues size is different from the number of pixel in the cluster");
    }
    
    _noiseValues    = noiseValues;
    _noiseSetSwitch = true;

  }

  template<class PixelType>
  std::vector<float > EUTelSparseClusterImpl<PixelType>::getNoiseValues() const {
    if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");
    return _noiseValues;
  }

  template<class PixelType>
  void EUTelSparseClusterImpl<PixelType>::getSeedCoord(int& xSeed, int& ySeed) const {
    unsigned int   maxIndex  =  0;
    float          maxSignal = -1 * std::numeric_limits<float>::max();
    PixelType * pixel = new PixelType;
    for ( unsigned int index = 0; index < size() ; index++ ) {
      getSparsePixelAt( index, pixel );
      if ( pixel->getSignal() > maxSignal ) {
 	maxSignal = pixel->getSignal();
 	maxIndex  = index;
      }
    }
    getSparsePixelAt( maxIndex, pixel ) ;
    xSeed = pixel->getXCoord();
    ySeed = pixel->getYCoord();
    delete pixel;
  }

  template<class PixelType>
  float EUTelSparseClusterImpl<PixelType>::getTotalCharge() const {
    float charge = 0;
    PixelType * pixel = new PixelType;
    for ( unsigned int index = 0; index < size() ; index++ ) {    
      getSparsePixelAt( index, pixel );
      charge += pixel->getSignal();
    }
    delete pixel;
    return charge;
  }

  template<class PixelType>
  float EUTelSparseClusterImpl<PixelType>::getSeedCharge() const {
    float          maxSignal = -1 * std::numeric_limits<float>::max();
    PixelType * pixel = new PixelType;
    for ( unsigned int index = 0; index < size() ; index++ ) {
      getSparsePixelAt(index , pixel);
      if ( pixel->getSignal() > maxSignal ) {
	maxSignal = pixel->getSignal();
      }
    }
    delete pixel;
    return maxSignal;
  }


  template<class PixelType>
  void EUTelSparseClusterImpl<PixelType>::getCenterOfGravityShift(float& xCoG, float& yCoG) const {

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
   
    PixelType * pixel = new PixelType;
    for ( unsigned int index = 0; index < size() ; index++ ) {
      getSparsePixelAt( index, pixel );
      tempX         += pixel->getSignal() * ( pixel->getXCoord() - xSeed );
      tempY         += pixel->getSignal() * ( pixel->getYCoord() - ySeed );
      normalization += pixel->getSignal() ;
//      printf(" getCoG %5d %5d %5.1f \n",  pixel->getXCoord(), pixel->getYCoord(), normalization);
    }
    if ( normalization != 0 ) {
      xCoG = tempX / normalization;
      yCoG = tempY / normalization;
    } else {
      xCoG = 0.;
      yCoG = 0.;
    }

    delete pixel;

  }

  template<class PixelType> 
  void EUTelSparseClusterImpl<PixelType>::getCenterOfGravityShift(float& xCoG, float& yCoG,
								  int xSize, int ySize) const {

    if ( size() == 1 ) {
      xCoG = 0; 
      yCoG = 0;
      return;
    }

    
    std::vector<unsigned int> goodPixelList;
    
    int xSeed, ySeed;
    getSeedCoord(xSeed, ySeed);

    PixelType * pixel = new PixelType;
    int xPixel, yPixel;

    for ( unsigned int iPixel = 0; iPixel < size() ; iPixel++ ) {
      getSparsePixelAt( iPixel, pixel );
      xPixel = static_cast<int>(pixel->getXCoord());
      yPixel = static_cast<int>(pixel->getYCoord());
      if ( ( abs( xSeed - xPixel ) <= ( xSize / 2 ) ) &&
	   ( abs( ySeed - yPixel ) <= ( ySize / 2 ) ) ) {
	goodPixelList.push_back(iPixel);
      }
    }
    
    float normalization = 0,  tempX = 0, tempY = 0;

    std::vector<unsigned int>::iterator iter = goodPixelList.begin();
    while ( iter != goodPixelList.end() ) {
      getSparsePixelAt( (*iter), pixel );
      tempX         += pixel->getSignal() * ( pixel->getXCoord() - xSeed ) ;
      tempY         += pixel->getSignal() * ( pixel->getYCoord() - ySeed ) ;
      normalization += pixel->getSignal();
      ++iter;
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
  void EUTelSparseClusterImpl<PixelType>::getCenterOfGravityShift(float& xCoG, float& yCoG,  int n) const {

    if ( size() == 1 ) {
      xCoG = 0; 
      yCoG = 0;
      return;
    }

    if ( static_cast<unsigned int>(n) >= size() ) {
      getCenterOfGravityShift( xCoG, yCoG );
      return;
    }


    std::multimap<float, unsigned int> pixelMap;
    PixelType * sparsePixel = new PixelType;
    
    for ( unsigned int iPixel = 0; iPixel < size() ; iPixel++ ) {
      getSparsePixelAt(iPixel, sparsePixel);
      pixelMap.insert( std::make_pair( sparsePixel->getSignal(), iPixel ) );
    }

    int counter = 0;
    float normalization = 0,  tempX = 0, tempY = 0;
    int xSeed, ySeed;
    getSeedCoord(xSeed, ySeed);

    xCoG = 0;
    yCoG = 0;


    std::multimap<float, unsigned int>::iterator mapIter = pixelMap.end();
    while ( counter < n ) {
      --mapIter;
      getSparsePixelAt( (*mapIter).second, sparsePixel ) ;
      tempX         += sparsePixel->getSignal() * ( sparsePixel->getXCoord() - xSeed ) ;
      tempY         += sparsePixel->getSignal() * ( sparsePixel->getYCoord() - ySeed ) ;
      normalization += sparsePixel->getSignal();
      ++counter;
    }

    if ( normalization != 0 ) {
      xCoG = tempX / normalization;
      yCoG = tempY / normalization;
    } else {
      xCoG = 0.;
      yCoG = 0.;
    }
    
    delete sparsePixel;
    return;
  }

/*
 * previous version (should be also ok, but sometimes not, why?)
  template<class PixelType> 
  void EUTelSparseClusterImpl<PixelType>::getCenterOfGravity(float&  xCoG, float& yCoG) const {
    getCenterOfGravityShift(xCoG, yCoG);
    
    int xSeed = 0;
    int ySeed = 0;
    getSeedCoord(xSeed, ySeed);
    
    xCoG += static_cast<float > ( xSeed );
    yCoG += static_cast<float > ( ySeed );

  }
*/

  //
  //direct calculation:
  //
  template<class PixelType> 
  void EUTelSparseClusterImpl<PixelType>::getCenterOfGravity(float&  xCoG, float& yCoG) const {
    
    PixelType* pixel = new PixelType;
    
    float xPos(0.0f), yPos(0.0f), totWeight(0.0f);
  
    for ( unsigned int index = 0; index < size() ; index++ ) 
    {
      getSparsePixelAt( index, pixel );
      float curSignal = pixel->getSignal(); 
      xPos += (pixel->getXCoord())*curSignal;
      yPos += (pixel->getYCoord())*curSignal;
      totWeight += curSignal;
    }

    xCoG = xPos / totWeight;
    yCoG = yPos / totWeight;
    
    delete pixel;
  }


  template<class PixelType>
  void EUTelSparseClusterImpl<PixelType>::getClusterSize(int& xSize, int& ySize) const {
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
  void EUTelSparseClusterImpl<PixelType>::getCenterCoord(int& xCenter, int& yCenter) const {
    // by definition, the cluster center is the pixel containing the
    // charge center of gravity. 
    float xCoG = 0.;
    float yCoG = 0.;
    getCenterOfGravity(xCoG, yCoG);

    // unfortunately I couldn't found a better way to approximate a
    // float number to the closest integer. It seems that the standard
    // libraries are missing this simple tools and this is the only
    // hack I found! 
    xCenter = static_cast<int> ( floor( xCoG + 0.5 ) );
    yCenter = static_cast<int> ( floor( yCoG + 0.5 ) );
    
  }
  
  template<class PixelType>
  float EUTelSparseClusterImpl<PixelType>::getDistance(EUTelVirtualCluster * clu) const {
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
  float EUTelSparseClusterImpl<PixelType>::getExternalRadius() const {
    int xSize, ySize;
    getClusterSize(xSize, ySize);
    return 0.5 * std::max(xSize, ySize);
  }

  template<class PixelType>
  float EUTelSparseClusterImpl<PixelType>::getClusterCharge(int nPixel) const {


    if ( static_cast<unsigned int> (nPixel) >= size() )
    {
      return getTotalCharge();
    }

    std::vector<float > allSignals;
    
    PixelType * sparsePixel = new PixelType;

    for ( unsigned int iPixel = 0; iPixel < size(); iPixel++ ) {
    
      getSparsePixelAt(iPixel, sparsePixel);
      allSignals.push_back( sparsePixel->getSignal() );

    }

    std::vector<float >::iterator iter = allSignals.begin();
    sort( allSignals.begin(), allSignals.end(), std::greater<float>() );

    float charge = 0;
    while ( iter != allSignals.begin() + nPixel ) {
      charge += *(iter);
      ++iter;
    }
    delete sparsePixel;
    
    return charge;
  }

  template<class PixelType>
  std::vector<float > EUTelSparseClusterImpl<PixelType>::getClusterCharge(std::vector<int > nPixels) const {
    
    std::vector<float > clusterSignal;
    
    std::vector<float > allSignals;
    PixelType * sparsePixel = new PixelType;
    for ( unsigned int iPixel = 0; iPixel < size(); iPixel++ ) {
      getSparsePixelAt(iPixel, sparsePixel);
      allSignals.push_back( sparsePixel->getSignal() );

    }
    sort( allSignals.begin(), allSignals.end(), std::greater<float>() );
    std::vector<float >::iterator iter;
    
    for (unsigned int i = 0; i < nPixels.size(); i++ ) {
      iter = allSignals.begin();
      float charge = 0;
      while ( iter != allSignals.begin() + nPixels[i]
            && iter != allSignals.end() ) {//add or for allsignals.end
	charge += (*iter);
	++iter;
      }
      clusterSignal.push_back(charge);
    }
    
    delete sparsePixel;
    return clusterSignal;
  }

  template<class PixelType> 
  float EUTelSparseClusterImpl<PixelType>::getClusterCharge(int xSize, int ySize) const {
    
    float charge = 0;
    
    int xSeed, ySeed;
    getSeedCoord(xSeed, ySeed);

    PixelType * pixel = new PixelType;
    int xPixel, yPixel;

    for ( unsigned int iPixel = 0; iPixel < size() ; iPixel++ ) {
      getSparsePixelAt( iPixel, pixel );
      xPixel = static_cast<int>( pixel->getXCoord());
      yPixel = static_cast<int>( pixel->getYCoord());
      if ( ( abs( xSeed - xPixel ) <= ( xSize / 2 ) ) &&
	   ( abs( ySeed - yPixel ) <= ( ySize / 2 ) ) ) {
	charge += pixel->getSignal();
      }
    }
    delete pixel;
    return charge;
  }

  template<class PixelType>
  void EUTelSparseClusterImpl<PixelType>::setClusterQuality(eutelescope::ClusterQuality quality)  {
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
  float EUTelSparseClusterImpl<PixelType>::getClusterNoise() const {
    
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
  float EUTelSparseClusterImpl<PixelType>::getClusterSNR() const {
    if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");
    float clusterNoise = getClusterNoise();
    if ( clusterNoise == 0 )  return 0.;
    float clusterSignal = getTotalCharge();
    return clusterSignal / clusterNoise;
  }

  template<class PixelType>
  float EUTelSparseClusterImpl<PixelType>::getSeedSNR() const {

      
    if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");
    unsigned int   maxIndex  = 0;
    float          maxSignal = -1 * std::numeric_limits<float>::max();
    PixelType * pixel = new PixelType;
    for ( unsigned int index = 0; index < size() ; index++ ) {
      getSparsePixelAt( index, pixel );
      if ( pixel->getSignal() > maxSignal ) {
 	maxSignal = pixel->getSignal();
 	maxIndex  = index;
      }
    }
    getSparsePixelAt( maxIndex, pixel ) ;
    delete pixel;
    return maxSignal / _noiseValues[maxIndex];
  }

  template<class PixelType>
  float EUTelSparseClusterImpl<PixelType>::getClusterSNR(int nPixel) const {
    if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");
    if ( static_cast<unsigned>(nPixel) >= size() ) 
      return getClusterSNR();

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
  float EUTelSparseClusterImpl<PixelType>::getClusterSNR(int xSize, int ySize) const {
    if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");
    
    int xSeed, ySeed;
    getSeedCoord(xSeed, ySeed) ;
    
    PixelType * pixel = new PixelType;
    int xPixel, yPixel;
    
    float charge = 0, noise2 = 0;

    for (unsigned int iPixel = 0;  iPixel < size() ; iPixel++ ) {
      getSparsePixelAt( iPixel, pixel );
      xPixel = static_cast<int> (pixel->getXCoord());
      yPixel = static_cast<int> (pixel->getYCoord());
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
  std::vector<float > EUTelSparseClusterImpl<PixelType>::getClusterSNR( std::vector<int > nPixels ) const {

    if ( ! _noiseSetSwitch ) throw DataNotAvailableException("No noise values set");
    
    // put all pixel values into a map
    PixelType * pixel = new PixelType;
    std::multimap<float, int> clusterSignalMap;

    for ( unsigned int iPixel = 0 ; iPixel < size(); iPixel++ ) {
      getSparsePixelAt( iPixel, pixel );
      clusterSignalMap.insert( std::make_pair( pixel->getSignal() , iPixel ) );
    }

    delete pixel;
    
    std::vector<int >::iterator pixelIter = nPixels.begin();
    std::vector<float > snr;
    
    while ( pixelIter != nPixels.end() ) {
      std::multimap<float, int >::reverse_iterator mapIter = clusterSignalMap.rbegin();
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
  void EUTelSparseClusterImpl<PixelType>::print(std::ostream& os) const {
    
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
