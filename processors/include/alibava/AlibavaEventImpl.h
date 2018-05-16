/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 *
 *  modified by: Thomas Eichhorn thomas.eichhorn@desy.de
 */

#ifndef ALIBAVAEVENTIMPL_H
#define ALIBAVAEVENTIMPL_H

// personal includes ".h"
#include "ALIBAVA.h"

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/LCEventImpl.h>

namespace alibava
{
    //! Implementation of the LCEvent for the Alibava data.
    /*!
     * @author Eda Yildirim, DESY <mailto:eda.yildirim@desy.de>
    */

    class AlibavaEventImpl : public IMPL::LCEventImpl
    {

	public:

	    AlibavaEventImpl ( );

	    virtual ~ AlibavaEventImpl ( )
	    {
		/* NO-OP */ ;
	    }

	    virtual void setEventType ( EventType type );

	    virtual void setEventType ( int type );

	    inline EventType getEventType ( ) const
	    {
		EventType type = static_cast < EventType > ( _params.getIntVal ( ALIBAVA::EVENTTYPE ) );
		return type;
	    }

	    virtual void setEventValue ( float avalue );

	    inline float getEventValue ( ) const
	    {
		float avalue = _params.getFloatVal ( ALIBAVA::EVENTVALUE );
		return avalue;
	    }

	    virtual void setEventSize ( int asize );

	    inline int getEventSize ( ) const
	    {
		int asize = _params.getIntVal ( ALIBAVA::EVENTSIZE );
		return asize;
	    }

	    virtual void setEventCode ( int acode );

	    inline int getEventCode ( ) const
	    {
		int acode = _params.getIntVal ( ALIBAVA::EVENTCODE );
		return acode;
	    }


	    virtual void setEventClock ( float aclock );

	    inline float getEventClock ( ) const
	    {
		float aclock = _params.getFloatVal ( ALIBAVA::EVENTCLOCK );
		return aclock;
	    }

	    virtual void setEventTime ( float atime );

	    inline float getEventTime ( ) const
	    {
		float atime = _params.getFloatVal ( ALIBAVA::EVENTTIME );
		return atime;
	    }

	    virtual void setEventTemp ( float atemp );

	    inline float getEventTemp ( ) const
	    {
		float atemp = _params.getFloatVal ( ALIBAVA::EVENTTEMP );
		return atemp;
	    }

	    virtual void setCalCharge ( float acharge );

	    inline float getCalCharge ( ) const
	    {
		float acharge = _params.getFloatVal ( ALIBAVA::CALCHARGE );
		return acharge;
	    }


	    virtual void setCalDelay ( float adelay );

	    inline float getCalDelay ( ) const
	    {
		float adelay = _params.getFloatVal ( ALIBAVA::CALDELAY );
		return adelay;
	    }

	    void maskEvent ( );

	    void unmaskEvent ( );

	    inline bool isEventMasked ( ) const
	    {
		bool amask = bool ( _params.getIntVal ( ALIBAVA::EVENTMASK ) );
		return amask;
	    }

    };
}

#endif // AlibavaEventImpl
