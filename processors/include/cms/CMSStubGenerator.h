/*
 * Created by Thomas Eichhorn
 *  (2017 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

#ifndef CMSStubGenerator_H
#define CMSStubGenerator_H

// eutelescope includes ".h"
#include "EUTelUtility.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>

// system includes <>
#include <string>
#include <vector>

namespace eutelescope
{

    class CMSStubGenerator : public marlin::Processor
    {

	private:
	    DISALLOW_COPY_AND_ASSIGN ( CMSStubGenerator )

	public:

	    virtual Processor * newProcessor ( )
	    {
		return new CMSStubGenerator;
	    }

	    CMSStubGenerator ( );

	    virtual void init ( );

	    virtual void processRunHeader ( LCRunHeader * run );

	    virtual void processEvent ( LCEvent * evt );

	    virtual void end ( );

	    void bookHistos( );

	protected:

	    std::string _inputHitCollectionName;

	    std::string _outputHitCollectionName;

	    int _dutPlane1;

	    int _dutPlane2;

	    bool _keepDUTHits;

	    float _maxResidual;

	    int _outputSensorID;

	    int _requirestubflag;

	    int _runMode;

	    int _totalstubs;

	    int _totalpl1;

	    int _totalpl2;

	    TrackerHitImpl* cloneHit ( TrackerHitImpl *inputHit );

	private:

    };

    //! A global instance of the processor
    CMSStubGenerator gCMSStubGenerator;

}
#endif
