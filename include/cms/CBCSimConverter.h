/*
 * Created by Thomas Eichhorn
 *  (2017 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

#ifndef CBCSimConverter_H
#define CBCSimConverter_H 1

// marlin includes ".h"
#include "marlin/Processor.h"

// system includes <>
#include <string>

namespace eutelescope
{

    //! CMS CBC clustering processor for Marlin.

    class CBCSimConverter : public marlin::Processor
    {

	public:

	    virtual Processor * newProcessor ( )
	    {
		return new CBCSimConverter;
	    }

	    //! Default constructor
	    CBCSimConverter ();

	    virtual void init ();

	    virtual void processRunHeader (LCRunHeader * run);

	    virtual void processEvent (LCEvent * evt);

	    virtual void check (LCEvent * evt);

	    void bookHistos();

	    void fillHistos();

	    virtual void end();

	    std::string _outputCollectionName;

	    std::string _inputCollectionName;

	    std::string _nonsensitiveaxis;

	    int _chancount;

	protected:

    };

    //! A global instance of the processor
    CBCSimConverter gCBCSimConverter;

}

#endif
