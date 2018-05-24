/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 *
 *  modified by: Thomas Eichhorn thomas.eichhorn@desy.de
 */

#ifndef ALIBAVARUNHEADERIMPL_H
#define ALIBAVARUNHEADERIMPL_H

// personal includes ".h"
#include "ALIBAVA.h"
#include "EUTELESCOPE.h"

// lcio includes <.h>
#include <lcio.h>
#include <EVENT/LCRunHeader.h>
#include <IMPL/LCRunHeaderImpl.h>

namespace alibava
{

    // Implementation of the Run Header for the Alibava readout system
    /* This is used to store the run header into the LCIO files produced 
     *  both by the Alibava system. This class is using the decorator pattern 
     *  around the LCRunHeaderImpl. The following parameters have been defined:
     *
     *  \li <b>HeaderVersion</b>: a float number representing the
     *  version of this header class. Standard version number v01-23-04
     *  are converted in a float number like 1.2304.
     *
     * Author: Eda Yildirim, DESY <mailto:eda.yildirim@desy.de>
     * Based on EUTelRunHeaderImpl.cc written by Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
    */

    class AlibavaRunHeaderImpl
    {

	public:

	    AlibavaRunHeaderImpl ( lcio::LCRunHeader * lcHeader ) : _lcHeader(nullptr)
	    {
		_lcHeader = dynamic_cast < IMPL::LCRunHeaderImpl* > ( lcHeader );
	    }

	    virtual ~ AlibavaRunHeaderImpl ( )
	    {
		/* NO-OP */;
	    }

	    // Set the run number. This is an integer storing the run number of the alibava data file
	    virtual void setRunNumber ( int runnum );

	    // return the header
	    inline int getRunNumber ( ) const
	    {
		return _lcHeader -> getRunNumber ( );
	    }

	    virtual void setDetectorName ( std::string name );

	    // return the header
	    inline std::string getDetectorName ( ) const
	    {
		return _lcHeader -> getDetectorName ( );
	    }

	    // Set the header. This is a string storing the header of the alibava data file
	    virtual void setHeader ( std::string aheader );

	    // return the header
	    inline std::string getHeader ( ) const
	    {
		return _lcHeader -> parameters ( ) .getStringVal ( ALIBAVA::HEADER );
	    }

	    // Set the header version number. This is an int number representing the version of this header
	    virtual void setHeaderVersion ( int ver );

	    // return the header version
	    inline int getHeaderVersion ( ) const
	    {
		return _lcHeader -> parameters ( ) .getIntVal ( ALIBAVA::HEADERVERSION );
	    }

	    // Set the number of events in the file
	    virtual void setNoOfEvents ( int num );

	    // return the number of events
	    inline int getNoOfEvents ( ) const
	    {
		return _lcHeader -> parameters ( ) .getIntVal ( ALIBAVA::NOOFEVENT );
	    }

	    // Set the data type. Here ALIBAVA::DAQDATA and  ALIBAVA::SIMULDATA are available
	    virtual void setDataType ( int type );

	    // return the data type
	    inline int getDataType ( ) const
	    {
		return _lcHeader -> parameters ( ) .getIntVal ( ALIBAVA::DATATYPE );
	    }

	    virtual void setGeoID ( int id );

	    inline int getGeoID ( ) const
	    {
		return _lcHeader -> parameters ( ) .getIntVal ( ALIBAVA::GEOID );
	    }

	    // Set the data taking date and time 
	    virtual void setDateTime ( std::string atime );

	    inline std::string getDateTime ( ) const
	    {
		return _lcHeader -> parameters ( ) .getStringVal ( ALIBAVA::DATETIME );
	    }

	    // Set the pedestal values stored in run header
	    virtual void setHeaderPedestal ( EVENT::FloatVec v_ped );

	    // return pedestal values stored in the run header
	    inline EVENT::FloatVec getHeaderPedestal ( ) const
	    {
		EVENT::FloatVec v_ped;
		v_ped = _lcHeader -> parameters ( ) .getFloatVals ( ALIBAVA::HEADERPEDESTAL, v_ped );
		return v_ped;
	    }

	    // Set the pedestal values stored in run header
	    virtual void setHeaderNoise ( EVENT::FloatVec v_noise );

	    // return pedestal values stored in run header
	    inline EVENT::FloatVec getHeaderNoise ( ) const
	    {
		EVENT::FloatVec v_noise;
		v_noise = _lcHeader -> parameters ( ) .getFloatVals ( ALIBAVA::HEADERNOISE, v_noise );
		return v_noise;
	    }

	    // Set selected chip numbers
	    virtual void setChipSelection ( EVENT::IntVec achipselectionvec );

	    // return selected chip numbers
	    inline EVENT::IntVec getChipSelection ( ) const
	    {
		EVENT::IntVec chipSelec;
		_lcHeader -> parameters ( ) .getIntVals ( ALIBAVA::SELECTEDCHIPNUM, chipSelec );
		return chipSelec;
	    }

	    // return number of selected chips
	    inline int getNChips ( ) const
	    {
		EVENT::IntVec chipSelec = getChipSelection ( );
		return chipSelec.size ( );
	    }

	    virtual void addProcessor ( std::string processor );

	    // returns the LCRunHeaderImpl underlying object
	    inline IMPL::LCRunHeaderImpl * lcRunHeader ( )
	    {
		return  _lcHeader;
	    }

	private:

	    DISALLOW_COPY_AND_ASSIGN ( AlibavaRunHeaderImpl )
	    IMPL::LCRunHeaderImpl * _lcHeader;
    };
}

#endif
