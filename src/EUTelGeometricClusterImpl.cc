#include "EUTelGeometricClusterImpl.h"
#include "EUTelGeometricPixel.h"

using namespace eutelescope;

EUTelGeometricClusterImpl::EUTelGeometricClusterImpl(IMPL::TrackerDataImpl* data): EUTelGenericSparseClusterImpl<EUTelGeometricPixel>(data)
{
/*nothing else to do*/
} 

EUTelGeometricClusterImpl::~EUTelGeometricClusterImpl()
{

}

void EUTelGeometricClusterImpl::getClusterGeomInfo(float& xPos, float& yPos, float& xSize, float& ySize) const 
{
	float xMin = std::numeric_limits<float>::max(); 	//stores the largest possible value every pixel will be lower, 
	float yMin = xMin;					//so its OK for max 
	float xMax = -xMin;								
	float yMax = xMax;
	float xMinBoundary = 0;
	float xMaxBoundary = 0;
	float yMinBoundary = 0;
	float yMaxBoundary = 0;

	EUTelGeometricPixel* pixel = new EUTelGeometricPixel;
	//Loop over all the pixels in the cluster
	for( unsigned int index = 0; index < size() ; index++ ) 
	{
		//Get the pixel
		getSparsePixelAt( index , pixel);

		//And its position
		float xCur = pixel->getPosX();
		float yCur = pixel->getPosY();

		if( xCur < xMin )
		{
			xMin = xCur;
			xMinBoundary = pixel->getBoundaryX();
		}
		if ( xCur > xMax ) 
		{
			xMax = xCur;
			xMaxBoundary = pixel->getBoundaryX();
		}
		if ( yCur < yMin )
		{
			yMin = yCur;
			yMinBoundary = pixel->getBoundaryY();
		}
		if ( yCur > yMax ) 
		{
			yMax = yCur;
			yMaxBoundary = pixel->getBoundaryY();
		}
	}
	xSize = xMax + xMaxBoundary - xMin + xMinBoundary;
	ySize = yMax + yMaxBoundary - yMin + yMinBoundary;

	xPos = xMax + xMaxBoundary - 0.5 *xSize;
	yPos = yMax + yMaxBoundary - 0.5 *ySize;
delete pixel;
}

void EUTelGeometricClusterImpl::getGeometricCenterOfGravity(float& xCoG, float& yCoG) const
{
	xCoG = 0;
	yCoG = 0;
	
	double totalCharge = 0;
	
	EUTelGeometricPixel* pixel = new EUTelGeometricPixel;
	for( unsigned int index = 0; index < size() ; index++ ) 
	{
		getSparsePixelAt( index , pixel);

		double curSignal = pixel->getSignal();
		xCoG += (pixel->getPosX())*curSignal;
		yCoG += (pixel->getPosY())*curSignal;
		totalCharge += curSignal;
	} 

	xCoG /= totalCharge;
	yCoG /= totalCharge;
	delete pixel;
}


