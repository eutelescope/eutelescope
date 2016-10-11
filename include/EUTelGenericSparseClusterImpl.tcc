/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
//TODO: Documentation
#ifndef EUTELGENERICSPARSECLUSTERIMPL_TCC
#define EUTELGENERICSPARSECLUSTERIMPL_TCC

#include <iostream>
#include <cmath>

namespace eutelescope{ 

template<class PixelType>
EUTelGenericSparseClusterImpl<PixelType>::EUTelGenericSparseClusterImpl(IMPL::TrackerDataImpl * data) : 
	EUTelSimpleVirtualCluster(data),
	EUTelClusterDataInterfacer<PixelType>(data),
	_nElement(0),
	_type(kUnknownPixelType)
{
	auto pixel = std::make_unique<PixelType>();
	_nElement = pixel->getNoOfElements();
	_type = pixel->getSparsePixelType();
}

template<class PixelType>
float EUTelGenericSparseClusterImpl<PixelType>::getTotalCharge() const 
{
	float charge = 0;
	
        auto& pixelVec = this->getPixels();
        for( auto& pixel: pixelVec ) {
		charge += pixel.getSignal();
   	 }
    return charge;
  }

template<class PixelType>
void EUTelGenericSparseClusterImpl<PixelType>::getClusterSize(int& xSize, int& ySize) const
{
	int xMin = std::numeric_limits<int>::max();	//stores the largest possible value
	int yMin = xMin;				//every pixel will be lower, so its OK for max
	int xMax = -1;					//pixel index starts at 0, so thats also ok
	int yMax = -1;

	auto& pixelVec = this->getPixels();
	for( auto& pixel: pixelVec ) {
		short xCur = pixel.getXCoord();
		short yCur = pixel.getYCoord();
		if ( xCur < xMin ) xMin = xCur;
		if ( xCur > xMax ) xMax = xCur;
		if ( yCur < yMin ) yMin = yCur;
		if ( yCur > yMax ) yMax = yCur;
	}
	xSize = xMax - xMin + 1;
	ySize = yMax - yMin + 1;
}
  
template<class PixelType>
void EUTelGenericSparseClusterImpl<PixelType>::getClusterInfo(int& xPos, int& yPos, int& xSize, int& ySize) const
{
	int xMin = std::numeric_limits<int>::max();	//stores the largest possible value
	int yMin = xMin;				//every pixel will be lower, so its OK for max
	int xMax = -1;					//pixel index starts at 0, so thats also ok
	int yMax = -1;
	
	auto& pixelVec = this->getPixels();
	for( auto& pixel: pixelVec ) {
		short xCur = pixel.getXCoord();
		short yCur = pixel.getYCoord();
		if ( xCur < xMin ) xMin = xCur;
		if ( xCur > xMax ) xMax = xCur;
		if ( yCur < yMin ) yMin = yCur;
		if ( yCur > yMax ) yMax = yCur;
	}

	xSize = xMax - xMin + 1;
	ySize = yMax - yMin + 1;
	
	xPos =  static_cast<int>( std::floor ( static_cast<float>(xMax) - 0.5 * static_cast<float>(xSize) + 0.5 ) );
	yPos =  static_cast<int>( std::floor ( static_cast<float>(yMax) - 0.5 * static_cast<float>(ySize) + 0.5 ) );
}

template<class PixelType>
void EUTelGenericSparseClusterImpl<PixelType>::getCenterOfGravity(float& xCoG, float& yCoG) const
{
	xCoG = 0;
	yCoG = 0;
	
	double totalCharge = 0;

        auto& pixelVec = this->getPixels();
        for( auto& pixel: pixelVec ) {
		double curSignal = pixel.getSignal();
		xCoG += (pixel.getXCoord())*curSignal;
		yCoG += (pixel.getYCoord())*curSignal;
		totalCharge += curSignal;
	} 
	xCoG /= totalCharge;
	yCoG /= totalCharge;
}

template<class PixelType> 
void EUTelGenericSparseClusterImpl<PixelType>::print(std::ostream& os) const {
   
    int   xSize, ySize/*, xSeed, ySeed, xCenter, yCenter*/;
    //ClusterQuality  quality = getClusterQuality();
    SparsePixelType type    = getSparsePixelType();
    getClusterSize(xSize, ySize);

    int bigspacer = 23;
  
    os   <<  std::setw(bigspacer) << std::setiosflags(std::ios::left) << "Sparse cluster made of " << type << " pixels" << std::endl
	 <<  std::setw(bigspacer) << "Number of pixel " << this->size() << std::endl
	 <<  std::setw(bigspacer) << "Cluster size " << "(" << xSize << ", " << ySize << ")" << std::endl
	 <<  std::setw(bigspacer) << "Cluster total charge " << getTotalCharge() << std::endl;

    int spacer = 14;
    for ( int i = 0; i < spacer - 1; i++ ) {
      os << "-";
    }
    os << std::endl;
    
    for( auto iPixel = this->begin(); iPixel != this->end() ; iPixel++ ) {
		os << "Pixel number = " << iPixel-(this->begin()) << std::endl
	 	   << ( *iPixel ) << std::endl;
    }
    for ( int i = 0; i < spacer - 1; i++ ) 
	{
      os << "-";
    }
    os << std::resetiosflags(std::ios::left) << std::endl;
  }
}
#endif
