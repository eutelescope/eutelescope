/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELGENERICSPARSECLUSTERIMPL_H
#define EUTELGENERICSPARSECLUSTERIMPL_H

// personal includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelBaseSparsePixel.h"
#include "EUTelExceptions.h"
#include "EUTelGenericSparsePixel.h"
#include "EUTelGeometricPixel.h"
#include "EUTelSimpleVirtualCluster.h"
#include "EUTelTrackerDataInterfacerImpl.h"

#ifdef USE_MARLIN
// marling includes ".h"
#include <marlin/Exceptions.h>
#endif

// lcio includes <.h>
#include <IMPL/TrackerDataImpl.h>
#include <LCIOTypes.h>

// system includes <>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <vector>

// template implementation
#include "EUTelGenericSparseClusterImpl.hcc"
#include "EUTelGenericSparseClusterImpl.tcc"

#endif
