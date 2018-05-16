/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 *
 *  modified by: Thomas Eichhorn thomas.eichhorn@desy.de
 */

#ifndef ALIBAVAPEDESTALSUBTRACTION_H
#define ALIBAVAPEDESTALSUBTRACTION_H 1

// alibava includes ".h"
#include "AlibavaBaseProcessor.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/TrackerDataImpl.h>

// ROOT includes <>
#include "TObject.h"

// system includes <>
#include <string>
#include <list>

namespace alibava
{
    class AlibavaPedestalSubtraction : public alibava::AlibavaBaseProcessor
    {
	public:

	    virtual Processor * newProcessor ( )
	    {
		return new AlibavaPedestalSubtraction;
	    }

	    AlibavaPedestalSubtraction ( );

	    virtual void init ( );

	    virtual void processRunHeader ( LCRunHeader * run );

	    virtual void processEvent ( LCEvent * evt );

	    virtual void check ( LCEvent * evt );

	    void bookHistos ( );

	    void fillHistos ( TrackerDataImpl * trkdata );

	    virtual void end ( );

	protected:

    };

    AlibavaPedestalSubtraction gAlibavaPedestalSubtraction;
}

#endif
