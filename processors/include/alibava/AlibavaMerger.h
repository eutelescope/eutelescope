/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

#ifndef ALIBAVAMERGER_H
#define ALIBAVAMERGER_H 1

// alibava includes ".h"
#include "AlibavaBaseProcessor.h"

// marlin includes ".h"
#include "marlin/Processor.h"

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
    class AlibavaMerger:public alibava::AlibavaBaseProcessor
    {
	public:

	    virtual Processor * newProcessor ( )
	    {
		return new AlibavaMerger;
	    }

	    AlibavaMerger ( );

	    virtual void init ( );

	    virtual void processRunHeader ( LCRunHeader * run );

	    virtual void processEvent ( LCEvent * evt );

	    virtual void check ( LCEvent * evt );

	    void bookHistos ( );

	    void fillHistos ( );

	    virtual void end ( );

	    //! The alibava file we want to read
	    std::string _alibavaFile;

	    //! The telescope file we want to read
	    std::string _telescopeFile;

	    //! The alibava collection name we want to read
	    std::string _alibavaCollectionName;

	    //! The telescope collection name
	    std::string _telescopeCollectionName;

	    //! The alibava collection name no2
	    std::string _alibavaCollectionName2;

	    //! The telescope collection name no2
	    std::string _telescopeCollectionName2;

	    //! The output collection names
	    std::string _outputCollectionName;

	    //! The output collection names no2
	    std::string _outputCollectionName2;

	    //! The output collection names no3
	    std::string _outputCollectionName3;

	    //! How do we want to output?
	    int _outputmode;

	    //! Move telescope sensor id?
	    int _teleplaneshift;

	    //! The reading instance
	    LCReader* lcReader;

	    //! The flag if the file is open
	    bool _telescopeopen;

	    //! The reading function
	    LCEvent *readTelescope ( );

	    void addCorrelation ( float ali_x, float ali_y, float ali_z, float tele_x, float tele_y, float tele_z, int event );

	    //! The unsensitive axis of our strip sensor
	    std::string _nonsensitiveaxis;

	    //! Difference between the systems
	    int _eventdifferenceTelescope;
	    int _eventdifferenceAlibava;

	protected:

    };

    AlibavaMerger gAlibavaMerger;
}

#endif
