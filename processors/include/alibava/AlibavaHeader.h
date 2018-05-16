/*
 * Created by Thomas Eichhorn
 *  (2014, 2017 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

#ifndef AlibavaHeader_H
#define AlibavaHeader_H 1

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

    class AlibavaHeader:public alibava::AlibavaBaseProcessor
    {

	public:

	    virtual Processor * newProcessor ( )
	    {
		return new AlibavaHeader;
	    }

	    AlibavaHeader ( );

	    virtual void init ( );

	    virtual void processRunHeader ( LCRunHeader * run );

	    virtual void processEvent ( LCEvent * evt );

	    virtual void check ( LCEvent * evt );

	    void bookHistos ( );

	    virtual void end ( );

	protected:

	    bool _coefffile;

	    void fillHistos ( TrackerDataImpl * headerdata, TrackerDataImpl * channeldata, int ichip );

	    void correlateLastHeader ( TrackerDataImpl * headerdata, TrackerDataImpl * channeldata, unsigned int ichip );

	    std::string _filterFileName;

	    std::string _rawdatacollection;

	};

	AlibavaHeader gAlibavaHeader;

}

#endif
