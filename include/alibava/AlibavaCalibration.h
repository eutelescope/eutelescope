/*
 * Created by Thomas Eichhorn
 *  (2017 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

#ifndef AlibavaCalibration_H
#define AlibavaCalibration_H 1

// alibava includes ".h"
#include "AlibavaBaseProcessor.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <IMPL/LCRunHeaderImpl.h>

// ROOT includes <>
#include "TObject.h"

// system includes <>
#include <string>


namespace alibava
{

    class AlibavaCalibration:public alibava::AlibavaBaseProcessor
    {

	public:

	    virtual Processor * newProcessor ( )
	    {
		return new AlibavaCalibration;
	    }

	    AlibavaCalibration ( );

	    virtual void check ( LCEvent * evt );

	    virtual void end ( );

	    virtual void init ( );

	    virtual void processEvent ( LCEvent * evt );

	    virtual void processRunHeader ( LCRunHeader * run );

	    void bookHistos ( );

	    void fillhisto ( int chip, int chan, double calc, double cald, double q );

	protected:

	    int _pol;

	    int _nccalpoints;

	    double _prevcal;

	    std::map < std::string , TObject * > _rootObjectMap;

	};

	AlibavaCalibration gAlibavaCalibration;

}

#endif
