/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 *
 */

#ifndef EUTelFixDUTAlignment_h
#define EUTelFixDUTAlignment_h 1

#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include "lcio.h"
#include "EUTelEventImpl.h"
#include "EUTelAlignmentConstant.h"

// system includes <>
#include <string>
#include <vector>
#include <map>

using namespace std;

namespace eutelescope
{

    class EUTelFixDUTAlignment : public marlin::Processor
    {
	public:

	    virtual Processor*  newProcessor ( )
	    {
		return new EUTelFixDUTAlignment;
	    }

	    EUTelFixDUTAlignment ( );

	    virtual void init ( );

	    virtual void processRunHeader ( LCRunHeader* run );

	    virtual void processEvent ( LCEvent * evt );

	    virtual void check ( LCEvent * evt );

	    void bookHistos ( );

	    virtual void end ( );

	    gear::SiPlanesParameters * _siPlanesParameters;

	    gear::SiPlanesLayerLayout * _siPlanesLayerLayout;

	protected:

	    int _mode;

	    std::string _outputfilename;

	    double _shiftx;

	    double _shifty;

	    double _shiftz;

	    double _rota;

	    double _rotb;

	    double _rotc;

	    int _manualDUTid;

	    int _eventstowrite;

	    std::string _dummyhitcollectionname;

	    void _EulerRotation ( double* _telPos, double* _gRotation );

    } ;

    EUTelFixDUTAlignment aEUTelFixDUTAlignment ;

}

#endif
