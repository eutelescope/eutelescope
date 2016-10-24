/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELTRACKERDATAINTERFACERIMPL_H
#define EUTELTRACKERDATAINTERFACERIMPL_H

// personal includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelSimpleSparsePixel.h"
#include "EUTelGenericSparsePixel.h"
#include "EUTelGeometricPixel.h"
#include "EUTelMuPixel.h"
#include "EUTelTrackerDataInterfacer.h"

#ifdef USE_MARLIN
// marling includes ".h"
#include <marlin/Exceptions.h>
#endif

// lcio includes <.h>
#include <LCIOTypes.h>
#include <IMPL/TrackerDataImpl.h>

//system includes
#include <memory>

// template implementation
#include "EUTelTrackerDataInterfacerImpl.hcc"
#include "EUTelTrackerDataInterfacerImpl.tcc"

#endif

