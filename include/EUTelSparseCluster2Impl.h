// Version: $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELSPARSECLUSTER2IMPL_H
#define EUTELSPARSECLUSTER2IMPL_H

// personal includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelBaseSparsePixel.h"
#include "EUTelSimpleSparsePixel.h"
#include "EUTelSparseDataImpl.h"
#include "EUTelExceptions.h"

// marling includes ".h"
#include <marlin/Exceptions.h>

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
#include "EUTelSparseCluster2Impl.hcc"
#include "EUTelSparseCluster2Impl.tcc"

#endif

