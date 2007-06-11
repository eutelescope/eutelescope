// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author:  Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version: $Id: EUTelSparseDataImpl.cc,v 1.1 2007-06-11 22:18:04 bulgheroni Exp $

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifdef EXPERIMENTAL

// personal includes ".h"
#include "EUTelSparseDataImpl.h"

// lcio includes <.h>
#include <IMPL/TrackerDataImpl.h>

// system includes
#include <vector>

using namespace eutelescope;
using namespace IMPL;
using namespace std;

unsigned const int EUTelSparseDataImpl::_nElement = 8;

EUTelSparseDataImpl::EUTelSparseDataImpl(TrackerDataImpl * data) { _trackerData = data; } 

inline TrackerDataImpl * EUTelSparseDataImpl::trackerData() { return _trackerData ; }

inline unsigned int EUTelSparseDataImpl::size() { return _trackerData->adcValues().size() / _nElement ; }

void EUTelSparseDataImpl::addSparsePixel(EUTelSparsePixel * /* pixel
							     */ ) {
  return;
}

EUTelSparsePixel * EUTelSparseDataImpl::getSparsePixelAt(int /* index
							      */ ) {
  return 0x0;
}

vector<EUTelClusterImpl *> EUTelSparseDataImpl::findClusters(double /* minDistance */ ) {
  return 0;
}

#endif 
