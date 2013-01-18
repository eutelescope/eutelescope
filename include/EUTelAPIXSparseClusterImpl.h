// Version: $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELAPIXSPARSECLUSTERIMPL_H
#define EUTELAPIXSPARSECLUSTERIMPL_H

// personal includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelBaseSparsePixel.h"
#include "EUTelSimpleSparsePixel.h"
#include "EUTelSparseDataImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelExceptions.h"
#include "EUTelAPIXSparsePixel.h"

#ifdef USE_MARLIN
// marling includes ".h"
#include <marlin/Exceptions.h>
#endif

// lcio includes <.h>
#include <LCIOTypes.h>
#include <IMPL/TrackerDataImpl.h>

// system includes <>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <iomanip>

// template implementation
//#include "EUTelSparseClusterImpl.hcc"
//#include "EUTelSparseClusterImpl.tcc"

namespace eutelescope {

	class EUTelAPIXSparseClusterImpl : public EUTelSparseClusterImpl<EUTelAPIXSparsePixel> {
	//class EUTelAPIXSparseClusterImpl : public EUTelVirtualCluster {
		public:
    	//! Default constructor
		EUTelAPIXSparseClusterImpl(IMPL::TrackerDataImpl * data);
		void getSeedCoord2(float& xSeed, float& ySeed) const ;
		void  getCenterOfGravityShift2(float& x, float& y) const  ;
		//! Destructor
		~EUTelAPIXSparseClusterImpl();

		//void getCenterOfGravityShift(float& xCoG, float& yCoG) const ;

		// to store coordinates of pixels in the cluster in the same order as in _trackerData
		std::vector <float> _X;
		std::vector <float> _Y;
	};



}


#endif

