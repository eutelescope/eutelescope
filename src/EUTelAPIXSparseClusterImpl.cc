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
#include "EUTelSparseDataImpl.h"
#include "EUTelAPIXSparseClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelVirtualCluster.h"

// lcio includes <.h>
#include <IMPL/TrackerDataImpl.h>

// system includes
#include <vector>

using namespace eutelescope;
using namespace IMPL;
using namespace std;



namespace eutelescope {

	EUTelAPIXSparseClusterImpl::EUTelAPIXSparseClusterImpl(IMPL::TrackerDataImpl * data) :
		EUTelSparseClusterImpl<EUTelAPIXSparsePixel>::EUTelSparseClusterImpl(data)	{

		for ( unsigned int k=0; k < size(); k++) {
			EUTelAPIXSparsePixel * Pixel = new EUTelAPIXSparsePixel;
			getSparsePixelAt(k, Pixel);
			int XIndex = Pixel -> getXCoord();
			int YIndex = Pixel -> getYCoord();

			float	XCoord, YCoord;

			if ( XIndex == 0 ) XCoord = 0.3;
			else if ( XIndex > 0 && XIndex < 17 ) { XCoord = 0.3 + 0.1 + 0.4 * XIndex; }
			else if (XIndex == 17) { XCoord = 0.3 + 0.1 + 17 * 0.4 + 0.1; }

			YCoord = 0.025 + 0.05 * YIndex;

			_X.push_back(XCoord);
			_Y.push_back(YCoord);
		}
	}


	void EUTelAPIXSparseClusterImpl::getSeedCoord2(float& xSeed, float& ySeed) const {
		unsigned int   maxIndex  =  0;
		float          maxSignal = -1 * std::numeric_limits<float>::max();
		EUTelAPIXSparsePixel * pixel = new EUTelAPIXSparsePixel;
		for ( unsigned int index = 0; index < size() ; index++ ) {
			getSparsePixelAt( index, pixel );
			if ( pixel->getSignal() > maxSignal ) {
				maxSignal = pixel->getSignal();
				maxIndex  = index;
			}
		}

		xSeed = _X[maxIndex];
		ySeed = _Y[maxIndex];

		delete pixel;
	}

	void EUTelAPIXSparseClusterImpl::getCenterOfGravityShift2(float& xCoG, float& yCoG) const {
		if ( size() == 1 ) {
			xCoG = 0;
			yCoG = 0;
			return;
		}

		float normalization = 0;
		float tempX = 0;
		float tempY = 0;

		float xSeed, ySeed;
		getSeedCoord2(xSeed, ySeed);

		EUTelAPIXSparsePixel * pixel = new EUTelAPIXSparsePixel;
		for ( unsigned int index = 0; index < size() ; index++ ) {
			getSparsePixelAt( index, pixel );
			tempX         += pixel->getSignal() * ( _X[index] - xSeed );
			tempY         += pixel->getSignal() * ( _Y[index] - ySeed );
			normalization += pixel->getSignal() ;
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

	EUTelAPIXSparseClusterImpl::~EUTelAPIXSparseClusterImpl(){
	}

}

