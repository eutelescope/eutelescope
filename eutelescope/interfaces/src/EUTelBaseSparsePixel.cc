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
#include "EUTelBaseSparsePixel.h"

// system 
#include <cmath>

using namespace eutelescope;

float eutelescope::distance(EUTelBaseSparsePixel * first, EUTelBaseSparsePixel * second) {
  
  return std::sqrt( pow( first->getXCoord() - second->getXCoord() , 2) +
		    pow( first->getYCoord() - second->getYCoord() , 2) );

}
