
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
#include "EUTelDFFClusterImpl.h"
#include "EUTelExceptions.h"
#include "EUTelVirtualCluster.h"

// lcio includes <.h>
#include <Exceptions.h>
#include <IMPL/TrackerDataImpl.h>
#include <LCIOTypes.h>

// system includes
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

using namespace eutelescope;
using namespace IMPL;
using namespace std;
