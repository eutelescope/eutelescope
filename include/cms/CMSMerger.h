/*
 * Created by Thomas Eichhorn
 *  (2017 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

#ifndef CMSMERGER_H
#define CMSMERGER_H 1

// marlin includes ".h"
#include "marlin/Processor.h"

// system includes <>
#include <string>

namespace eutelescope
{

    //! CMS CBC merge processor for Marlin.

    class CMSMerger : public marlin::Processor
    {

	public:

	    virtual Processor * newProcessor ( )
	    {
		return new CMSMerger;
	    }

	    //! Default constructor
	    CMSMerger ();

	    virtual void init ();

	    virtual void processRunHeader (LCRunHeader * run);

	    virtual void processEvent (LCEvent * evt);

	    virtual void check (LCEvent * evt);

	    void bookHistos();

	    void fillHistos();

	    virtual void end();

	    std::string _telescopeFile;

	    std::string _cbcDataCollectionName1;

	    std::string _cbcPulseCollectionName1;

	    std::string _cbcDataCollectionName2;

	    std::string _cbcPulseCollectionName2;

	    std::string _telescopeCollectionName;

	    std::string _telescopeCollectionName2;

	    std::string _outputCollectionName;

	    std::string _outputCollectionName2;

	    std::string _outputCollectionName3;

	    bool _telescopeopen;

	    LCReader* lcReader;

	    long _cbceventtime;
	    long _telescopeeventtime;

	    int _readcount;
	    int _maxevents;
	    int _multiplicity;

	    LCEvent* _storeevt;

	    LCEvent *readTelescope ();

	protected:

    };

    //! A global instance of the processor
    CMSMerger gCMSMerger;

}

#endif
