/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 *
 *  modified by: Thomas Eichhorn thomas.eichhorn@desy.de
 */

#ifndef ALIBAVATIMECUTPROCESSOR_H
#define ALIBAVATIMECUTPROCESSOR_H 1

// alibava includes ".h"
#include "AlibavaBaseProcessor.h"
#include "AlibavaEventImpl.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/EventModifier.h"

// lcio includes <.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/TrackerRawDataImpl.h>

// ROOT includes <>
#include "TObject.h"

// system includes <>
#include <string>
#include <list>

namespace alibava
{

    class AlibavaTimeCutProcessor : public alibava::AlibavaBaseProcessor
    {
	public:

	    virtual Processor * newProcessor ( )
	    {
		return new AlibavaTimeCutProcessor;
	    }

	    AlibavaTimeCutProcessor ( );

	    virtual void init ( );

	    virtual void processRunHeader ( LCRunHeader * run );

	    virtual void processEvent ( LCEvent * evt );

	    virtual void check ( LCEvent * evt );

	    void bookHistos ( );

	    void fillHistos ( AlibavaEventImpl * anAlibavaEvent );

	    virtual void end ( );

	    float getTimeCutMin ( );

	    float getTimeCutMax ( );

	protected:

	    // The minimum tdc time that is acceptable to use that Event
	    float _timecutmin;

	    // The maximum tdc time that is acceptable to use that Event
	    float _timecutmax;

	    // Number of masked events
	    int _numberOfMaskedEvents;

	    // The name of the histogram which contains event number of masked events
	    std::string _hMaskedEventsNumberName;

	
    };

    AlibavaTimeCutProcessor gAlibavaTimeCutProcessor;
}

#endif
