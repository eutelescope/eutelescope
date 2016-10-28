
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
#include "EUTelVirtualCluster.h"
#include "EUTelExceptions.h"

// lcio includes <.h>
#include <LCIOTypes.h>
#include <IMPL/TrackerDataImpl.h>
#include <Exceptions.h>

// system includes
#include <map>
#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cassert>

using namespace eutelescope;
using namespace IMPL;
using namespace std;



