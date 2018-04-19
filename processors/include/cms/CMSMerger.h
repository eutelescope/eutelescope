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
	    CMSMerger ( );

	    bool _eventmerge;

	    bool _telescopeopen;

	    int _correlationPlaneID;

	    int _eventdifferenceTelescope;

	    int _eventdifferenceCBC;

	    int _maxevents;

	    int _multiplicity;

	    int _readcount;

	    LCEvent* _storeevt;

	    LCEvent *readTelescope ( );

	    LCReader* lcReader;

	    long _cbceventtime;

	    long _telescopeeventtime;

	    std::string _telescopeFile;

	    std::string _cbcDataCollectionName1;

	    std::string _cbcPulseCollectionName1;

	    std::string _cbcDataCollectionName2;

	    std::string _cbcPulseCollectionName2;

	    std::string _outputCollectionName;

	    std::string _outputCollectionName2;

	    std::string _outputCollectionName3;

	    std::string _telescopeCollectionName;

	    std::string _telescopeCollectionName2;

	    virtual void check ( LCEvent * evt );

	    virtual void end ( );

	    virtual void init ( );

	    virtual void processRunHeader ( LCRunHeader * run );

	    virtual void processEvent ( LCEvent * evt );

	    void bookHistos ( );

	    void fillHistos ( );

	protected:

    };

    //! A global instance of the processor
    CMSMerger gCMSMerger;

}

#endif
