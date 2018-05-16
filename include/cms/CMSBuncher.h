/*
 * Created by Thomas Eichhorn
 *  (2017 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

#ifndef CMSBuncher_H
#define CMSBuncher_H 1

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/DataSourceProcessor.h"

// system includes <>
#include <string>

namespace eutelescope
{

    class CMSBuncher : public marlin::DataSourceProcessor
    {

	public:

	    CMSBuncher ( );

	    virtual CMSBuncher * newProcessor ( );

	    virtual void readDataSource ( int numEvents );

	    virtual void init ( );

	    virtual void end ( );

	    bool _telescopeopen;

	    LCEvent* _storeevt;

	    std::string _telescopeFile;

	    long _frametime;

	    long _bunchtime;

	    int _evtcount;

	    bool _usestorageevent;

	    int _singleframetime;

	    int _maxevents;

	    std::string _inputCollectionName;

	    std::string _outputCollectionName;

	    LCReader* lcReader;

	    LCEvent *readTelescope ();

	protected:

    };

    //! A global instance of the processor
    CMSBuncher gCMSBuncher;

}

#endif
