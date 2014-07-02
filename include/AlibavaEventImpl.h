/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


#ifndef ALIBAVAEVENTIMPL_H
#define ALIBAVAEVENTIMPL_H

// personal includes ".h"
#include "ALIBAVA.h"

// marlin includes ".h"

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/LCEventImpl.h>

// system includes <>


namespace alibava {

  //! Implementation of the LCEvent for the Alibava data.
  /*!
	* @author Eda Yildirim, DESY <mailto:eda.yildirim@desy.de>
   */

  class AlibavaEventImpl : public IMPL::LCEventImpl   {

  public:

    //! Default constructor
    AlibavaEventImpl ();

    //! Destructor
    virtual ~ AlibavaEventImpl() { /* NO-OP */ ; }

    //! Set the event type
    /*! This method is used to set the type of this
     *  AlibavaEventImpl. For this purpose the EventType enumeration. 
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
      EventType type = static_cast<EventType>(_params.getIntVal(ALIBAVA::EVENTTYPE));
      return type;
    }

	  
	  //! Set the event value
	  /*! This method is used to set the value stored in alibava event
		*  See alibava documentation for more detail
		*
		*  @param avalue The event value
		*/
	  virtual void setEventValue(float avalue);
	  //! Return the event value
	  /*! This method returns the event value
		*
		*  @return The event value
		*/
	  inline float getEventValue() const {
		  float avalue = _params.getFloatVal(ALIBAVA::EVENTVALUE);
		  return avalue;
	  }
	  
	  //! Set the event size
	  /*! This method is used to set the size of alibava event
		*  See alibava documentation for more detail
		*
		*  @param asize The event size
		*/
	  virtual void setEventSize(int asize);
	  //! Return the event size
	  /*! This method returns the event size
		*
		*  @return The event size
		*/
	  inline int getEventSize() const {
		  int asize = _params.getIntVal(ALIBAVA::EVENTSIZE);
		  return asize;
	  }

	  //! Set the event code
	  /*! This method is used to set the code of alibava event
		*  See alibava documentation for more detail
		*
		*  @param acode The event size
		*/
	  virtual void setEventCode(int acode);
	  //! Return the event code
	  /*! This method returns the event code
		*
		*  @return The event code
		*/
	  inline int getEventCode() const {
		  int acode = _params.getIntVal(ALIBAVA::EVENTCODE);
		  return acode;
	  }

	  //! Set the TDC time
	  /*! This method is used to set the tdc time
		*  See alibava documentation for more detail
		*
		*  @param aTIME The event size
		*/
	  virtual void setEventTime(float atime);
	  //! Return the TDC time
	  /*! This method returns the event size
		*
		*  @return The event time
		*/
	  inline float getEventTime() const {
		  float atime = _params.getFloatVal(ALIBAVA::EVENTTIME);
		  return atime;
	  }

	  //! Set the temperature
	  /*! This method is used to set the temperature
		*  See alibava documentation for more detail
		*
		*  @param atemp The event size
		*/
	  virtual void setEventTemp(float atemp);
	  //! Return the temperature
	  /*! This method returns the temperature
		*
		*  @return The temperature
		*/
	  inline float getEventTemp() const {
		  float atemp = _params.getFloatVal(ALIBAVA::EVENTTEMP);
		  return atemp;
	  }

	  //! Set the injected charge
	  /*! This method is used to set the injected charge 
		*  in charge calibration run events
		*  See alibava documentation for more detail
		*
		*  @param acharge The injected charge
		*/
	  virtual void setCalCharge(float acharge);
	  //! Return the injected charge
	  /*! This method returns the injected charge
		*  in charge calibration run events
		*
		*  @return The injected charge
		*/
	  inline float getCalCharge() const {
		  float acharge = _params.getFloatVal(ALIBAVA::CALCHARGE);
		  return acharge;
	  }

	  //! Set the delay
	  /*! This method is used to set the delay
		*  in delay calibration run events
		*  See alibava documentation for more detail
		*
		*  @param adelay The delay
		*/
	  virtual void setCalDelay(float adelay);
	  //! Return the delay
	  /*! This method returns thedelay
		*  in delay calibration run events
		*
		*  @return The delay
		*/
	  inline float getCalDelay() const {
		  float adelay = _params.getFloatVal(ALIBAVA::CALDELAY);
		  return adelay;
	  }
	
	  //! Mask the Event
	  /*! This method is used to mask event
		*  You might want to mask events to apply time cut
		*  or for any other
		*/
	  void maskEvent();

	  //! Unmask the event
	  /*! This method is used to unmask event
		*/
	  void unmaskEvent();

	  
	  //! Return the mask
	  /*! This method returns true if the event is masked
		*
		*/
	  inline bool isEventMasked() const {
		  bool amask = bool( _params.getIntVal(ALIBAVA::EVENTMASK) );
		  return amask;
	  }

	  
	  
	  
	  
  };                           // end of AlibavaEventImpl
}                              // alibava namespace

#endif // AlibavaEventImpl
