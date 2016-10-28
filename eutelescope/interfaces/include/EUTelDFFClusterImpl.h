// Version: $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELDFFCLUSTERIMPL_H
#define EUTELDFFCLUSTERIMPL_H

// personal includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"

// marlin includes ".h"

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

// system includes <>
#include <vector>
#include <iostream>


namespace eutelescope {

  class EUTelDFFClusterImpl : public EUTelFFClusterImpl {
  public:
    EUTelDFFClusterImpl(TrackerDataImpl * data) : EUTelFFClusterImpl(data)
    {
      
    }

 };
 
}
#endif
