// Author:  Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version: $Id: EUTelRunHeaderImpl.cc,v 1.3 2007-02-22 08:09:36 bulgheroni Exp $

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// personal includes ".h"
#include "EUTelRunHeaderImpl.h"
#include "EUTELESCOPE.h"

// marlin includes ".h"

// lcio includes <.h>
#include "lcio.h"
#include "LCIOTypes.h"
#include "IMPL/LCRunHeaderImpl.h"
#include "UTIL/LCTime.h"

using namespace std;
using namespace eutelescope;


EUTelRunHeaderImpl::EUTelRunHeaderImpl ():
IMPL::LCRunHeaderImpl ()
{ /* NO - OP */ ;
}

void
EUTelRunHeaderImpl::setHeaderVersion (float ver)
{
   // sets the header implementation version number

   _params.setValue (EUTELESCOPE::HEADERVERSION, ver);
}

void
EUTelRunHeaderImpl::setDataType (std::string type)
{
   // sets the type of data saved in this file

   _params.setValue (EUTELESCOPE::DATATYPE, type);
}

void
EUTelRunHeaderImpl::setNoOfEvent ( int num )
{
  // sets the number of events in the file
  
  _params.setValue(EUTELESCOPE::NOOFEVENT, num);

}

void
EUTelRunHeaderImpl::setDateTime ()
{
   // sets the date and time of this file to now
   // it uses the LCTime util

   UTIL::LCTime * now = new UTIL::LCTime ();
   _params.setValue (EUTELESCOPE::DATETIME, now->getDateString ());
   delete now;
}

void
EUTelRunHeaderImpl::setDAQHWName (std::string name)
{
   // sets the DAQ hardware name.
   // only in the case of EUTELESCOPE::DAQDATA

   _params.setValue (EUTELESCOPE::DAQHWNAME, name);
}

void
EUTelRunHeaderImpl::setDAQHWVersion (float ver)
{
   // sets the DAQ hardware version number
   // only in the case of EUTELESCOPE::DAQDATA

   _params.setValue (EUTELESCOPE::DAQHWVERSION, ver);
}

void
EUTelRunHeaderImpl::setDAQSWName (std::string name)
{
   // sets the DAQ software name
   // only in the case of EUTELESCOPE::DAQDATA

   _params.setValue (EUTELESCOPE::DAQSWNAME, name);
}

void
EUTelRunHeaderImpl::setDAQSWVersion (float ver)
{
   // sets the DAQ software version number
   // only in the case of EUTELESCOPE::DAQDATA

   _params.setValue (EUTELESCOPE::DAQSWVERSION, ver);
}

void
EUTelRunHeaderImpl::setSimulSWName (std::string name)
{
   // sets the simulation program name used to generate this data
   // only in the case of EUTELESCOPE::SIMULDATA

   _params.setValue (EUTELESCOPE::SIMULSWNAME, name);
}

void
EUTelRunHeaderImpl::setSimulSWVersion (float ver)
{
   // sets the simulation program version number
   // only in the case of EUTELESCOPE::SIMULDATA

   _params.setValue (EUTELESCOPE::SIMULSWVERSION, ver);
}

void
EUTelRunHeaderImpl::setGeoID (int id)
{
   // sets the current geometry identification number
   // use this number to query the geometry database

   _params.setValue (EUTELESCOPE::GEOID, id);
}

void
EUTelRunHeaderImpl::setNoOfDetector (int num)
{
   // sets the number of detectors in this run.  this corresponds to
   // the number of TrackerRawData objects in each LCCollectionVec

   _params.setValue (EUTELESCOPE::NOOFDETECTOR, num);
}

void
EUTelRunHeaderImpl::setMinX (lcio::IntVec xMin)
{
   // sets the vector containing the minimum pixel number along X

   _params.setValues (EUTELESCOPE::MINX, xMin);
}

void
EUTelRunHeaderImpl::setMaxX (lcio::IntVec xMax)
{
   // sets the vector containing the maximum pixel number along X

   _params.setValues (EUTELESCOPE::MAXX, xMax);
}

void
EUTelRunHeaderImpl::setMinY (lcio::IntVec yMin)
{
   // sets the vector containing the minimum pixel number along Y

   _params.setValues (EUTELESCOPE::MINY, yMin);
}

void
EUTelRunHeaderImpl::setMaxY (lcio::IntVec yMax)
{
   // sets the vector containing the maximum pixel number along Y

   _params.setValues (EUTELESCOPE::MAXY, yMax);
}

void
EUTelRunHeaderImpl::addProcessor (std::string processor)
{
   // add a "processor" to the list of applied processors

   lcio::StringVec processorVec;
   _params.getStringVals (EUTELESCOPE::APPLIEDPROCESSOR, processorVec);
   processorVec.push_back (processor);
   _params.setValues (EUTELESCOPE::APPLIEDPROCESSOR, processorVec);
}

void
EUTelRunHeaderImpl::addIntermediateFile (std::string file)
{
   // add a "file" to the list of intermediate file

   lcio::StringVec fileVec;
   _params.getStringVals (EUTELESCOPE::INTERMEDIATEFILE, fileVec);
   fileVec.push_back (file);
   _params.setValues (EUTELESCOPE::INTERMEDIATEFILE, fileVec);
}
