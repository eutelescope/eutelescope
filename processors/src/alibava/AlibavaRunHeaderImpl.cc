/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 *
 *  modified by: Thomas Eichhorn thomas.eichhorn@desy.de
 */

// personal includes ".h"
#include "AlibavaRunHeaderImpl.h"
#include "ALIBAVA.h"

// lcio includes <.h>
#include "lcio.h"
#include "LCIOTypes.h"
#include "UTIL/LCTime.h"

using namespace std;
using namespace alibava;

void AlibavaRunHeaderImpl::setDetectorName ( string name )
{
   _lcHeader -> setDetectorName ( name );
}

void AlibavaRunHeaderImpl::setRunNumber ( int runnum )
{
   _lcHeader -> setRunNumber ( runnum );
}

void AlibavaRunHeaderImpl::setHeader ( std::string aheader )
{
   _lcHeader -> parameters ( ) .setValue ( ALIBAVA::HEADER, aheader );
}

void AlibavaRunHeaderImpl::setHeaderVersion ( int ver )
{
    _lcHeader -> parameters ( ) .setValue ( ALIBAVA::HEADERVERSION, ver );
}

void AlibavaRunHeaderImpl::setDataType ( int type )
{
   _lcHeader -> parameters ( ) .setValue ( ALIBAVA::DATATYPE, type );
}

void AlibavaRunHeaderImpl::setNoOfEvents ( int num )
{
    _lcHeader -> parameters ( ) .setValue ( ALIBAVA::NOOFEVENT, num );
}

void AlibavaRunHeaderImpl::setGeoID ( int id )
{
    _lcHeader -> parameters ( ) .setValue ( ALIBAVA::GEOID, id );
}

void AlibavaRunHeaderImpl::setDateTime ( std::string atime )
{
   _lcHeader -> parameters ( ) .setValue ( ALIBAVA::DATETIME, atime );
}

void AlibavaRunHeaderImpl::setHeaderPedestal ( EVENT::FloatVec v_ped )
{
    _lcHeader -> parameters ( ) .setValues ( ALIBAVA::HEADERPEDESTAL, v_ped );
}

void AlibavaRunHeaderImpl::setHeaderNoise ( EVENT::FloatVec v_noise )
{
    _lcHeader -> parameters ( ) .setValues ( ALIBAVA::HEADERNOISE, v_noise );
}

void AlibavaRunHeaderImpl::setChipSelection ( EVENT::IntVec achipselection )
{
    _lcHeader -> parameters ( ) .setValues ( ALIBAVA::SELECTEDCHIPNUM, achipselection );
}

void AlibavaRunHeaderImpl::addProcessor ( std::string processor )
{
    lcio::StringVec processorVec;
    _lcHeader -> parameters ( ) .getStringVals ( ALIBAVA::APPLIEDPROCESSOR, processorVec );
    processorVec.push_back ( processor );
    _lcHeader -> parameters ( ) .setValues ( ALIBAVA::APPLIEDPROCESSOR, processorVec );
}
