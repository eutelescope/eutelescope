/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 *
 */

#ifndef EUTelDUTPosition_H
#define EUTelDUTPosition_H 1

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

    class EUTelDUTPosition : public marlin::Processor
    {

	public:

	    virtual Processor * newProcessor ( )
	    {
		return new EUTelDUTPosition;
	    }

	    EUTelDUTPosition ( );

	    virtual void init ( );

	    virtual void processRunHeader ( LCRunHeader * run );

	    virtual void processEvent ( LCEvent * evt );

	    virtual void check ( LCEvent * evt );

	    void bookHistos ( );

	    virtual void end ( );

	    //! Silicon plane parameters as described in GEAR. This object is provided by GEAR during the init() phase and stored here for local use.
	    gear::SiPlanesParameters * _siPlanesParameters;

	    //! This is the real geometry description for each layer. This object is taken from _siPlanesParameters during the init() phase and stored for local use
	    gear::SiPlanesLayerLayout * _siPlanesLayerLayout;

	protected:

	    std::string _outputFileName;

	    std::string _outputDUTFileName;

	    int _manualDUTid;

	    int _manualDUTposition;

	    std::string _finalcollectionname;

	};

    //! A global instance of the processor
    EUTelDUTPosition aEUTelDUTPosition;

}

#endif
