/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELEVENTIMPL_H
#define EUTELEVENTIMPL_H

// personal includes ".h"
#include "EUTELESCOPE.h"

// marlin includes ".h"

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/LCEventImpl.h>

// system includes <>


namespace eutelescope {

  //! Implementation of the LCEvent for the EUDET telescope.
  /*! This is a very trivial re-implementation of the LCEvent for the
   *  EUDET telescope and it offers the some (honestly very few for
   *  the time being) convenience methods to get useful things.
   *
   *  A typical example comes from the need during the analysis
   *  procedure to know how many events are contained into the file,
   *  or at least to know whether the current event is the last
   *  one. Being the LCIO data model based on the SIO data format,
   *  featuring a pure serial I/O access to the data on disk, it is
   *  not possible to know in advance how many events are stored into
   *  a file and if there will be another event after the current.  To
   *  workaround this limitation we are adopting the so called
   *  <b>BORE</b> (Beginning Of Run Event) and <b>EORE</b> (End Of Run
   *  Event) events used also in the Bonn DAQ software. Those events
   *  are saved into the file just with a flag identifying them among
   *  the real events. Actually there is no need to have a BORE since
   *  the first event is very well identified being the one following
   *  the RunHeader.
   *  All EUDET processors, as soon as they get the current event,
   *  they will check if this is a EORE, and if so, they simply call
   *  the <i>Finalize</i> method (if exists), otherwise it has to
   *  return immediately.
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
   * 
   */

  class EUTelEventImpl : public IMPL::LCEventImpl   {

  public:

    //! Default constructor
    EUTelEventImpl ();

    //! Destructor
    virtual ~ EUTelEventImpl() { /* NO-OP */ ; }

    //! Set the event type
    /*! This method is used to set the type of this
     *  EUTelEventImpl. For this purpose the EventType enumeration. 
     *
     *  @param type The event type in the form of a EventType enum
     *  object
     */
    virtual void setEventType(EventType type);

    //! Set the event type
    /*! This overloaded method is provided only for convenience and it
     *  can be used to set the event type using the integer number
     *  corresponding to the event type defined in the EventType enum.
     *
     *  @param type The event type as integer number
     */
    virtual void setEventType(int type);

    //! Return the event type
    /*! This method returns the event type as EventType enumeration
     *  object. 
     *
     *  @return The event type according to the EventType enum
     */
    inline EventType getEventType() const {
      EventType type = static_cast<EventType>(_params.getIntVal(EUTELESCOPE::EVENTTYPE));
      return type;
    }

  };                           // end of EUTelEventImpl
}                              // eutelescope namespace

#endif // EUTELEVENTIMPL
