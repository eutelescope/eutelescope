/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


#ifndef ALIBAVARUNHEADERIMPL_H
#define ALIBAVARUNHEADERIMPL_H

// personal includes ".h"
#include "ALIBAVA.h"

// marlin includes ".h"

// lcio includes <.h>
#include <lcio.h>
#include <EVENT/LCRunHeader.h>
#include <IMPL/LCRunHeaderImpl.h>

// system includes <>

namespace alibava {

  //! Implementation of the Run Header for the Alibava readout system
  /*! This is used to store the run header into the LCIO files produced 
	*  both by the Alibava system. This class is using the decorator pattern 
	*  around the LCRunHeaderImpl. The following parameters have been defined:
   *
   *  \li <b>HeaderVersion</b>: a float number representing the
   *  version of this header class. Standard version number v01-23-04
   *  are converted in a float number like 1.2304.
   *
	// Author: Eda Yildirim, DESY <mailto:eda.yildirim@desy.de>
	// Based on EUTelRunHeaderImpl.cc written by Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *
   */

	class AlibavaRunHeaderImpl {

  public:
    
    //! Default constructor
//		AlibavaRunHeaderImpl ( ) : _lcHeader(NULL) {};
		AlibavaRunHeaderImpl ( lcio::LCRunHeader * lcHeader ) : _lcHeader(NULL) { _lcHeader = dynamic_cast<IMPL::LCRunHeaderImpl*> (lcHeader) ; }

    //! Destructor
    virtual ~ AlibavaRunHeaderImpl ()  { /* NO-OP */ ;  }

		//! Set the run number
		/*! this is an integer storing the run number of alibava data file
		 *
		 *  @param arunnum The run number the user wants to set
		 */
		virtual void setRunNumber (int runnum);
		
		//! return the header
		inline int getRunNumber () const    {
			return _lcHeader->getRunNumber();
		}


	//! Set detector name
	/*! this is a string storing the detector name 
	*
	*  @param name The detector name
	*/
	virtual void setDetectorName (std::string name);
	
	//! return the header
	inline std::string getDetectorName() const    {
	return _lcHeader->getDetectorName();
        }
		
	  //! Set the header
	  /*! this is a string storing the header of alibava data file
		*
		*  @param aheader The version number the user wants to set
		*/
	  virtual void setHeader (std::string aheader);
     
	  //! return the header
		inline std::string getHeader () const    {
		  return _lcHeader->parameters().getStringVal (ALIBAVA::HEADER);
	  }

	  
    //! Set the header version number
    /*! this is an int number representing the version of this header
     *
     *  @param ver The version number the user wants to set
     */
    virtual void setHeaderVersion (int ver);
     
	  //! return the header version
	  inline int getHeaderVersion () const    {
		  return _lcHeader->parameters().getIntVal (ALIBAVA::HEADERVERSION);
	  }

    //! Set the number of events in the file
    /*! this is an integer number equal to the number of events in the
     *  file.
     * 
     *  @param num The number of events.
     */
    virtual void setNoOfEvents(int num);

	  //! return the number of events
	  inline int getNoOfEvents() const     {
		  return _lcHeader->parameters().getIntVal(ALIBAVA::NOOFEVENT);
	  }

    //! Set the data type
    /*! this string is used to distinguish real data run from
     *  simulation. In the ALIBAVA class few static constant
     *  string are available for this purpose: ALIBAVA::DAQDATA
     *  and ALIBAVA::SIMULDATA.
     *
     *  @param type The type of data contained into the file
     */
    virtual void setDataType (int type);
	  
	  //! return the data type
	  inline int getDataType () const     {
		  return _lcHeader->parameters().getIntVal (ALIBAVA::DATATYPE);
	  }


    //! Set the geometry identification number
    /*! this is an integer number used to link the
     *  current telescope geometrical configuration with an entry in the
     *  geometry database. This is used during the reconstruction phase
     *  into a DB query to retrieve precise information about the
     *  detector positioning and alignments.
     *
     *  @param id The identification number of the current geometry
     */
    virtual void setGeoID (int id);
	  
	  //! return the geometry identification number
	  inline int getGeoID () const    {
		  return _lcHeader->parameters().getIntVal (ALIBAVA::GEOID);
	  }

		//! Set the data taking date and time 
		/*! this string stores date and time of data
		 *
		 */
		virtual void setDateTime (std::string atime);

		
		//! return the date and time in a human readable format
		inline std::string getDateTime () const    {
			return _lcHeader->parameters().getStringVal (ALIBAVA::DATETIME);
		}

		//! Set the pedestal values stored in run header
		/*! std::vector<int> v_ped stores the values
		 *
		 */
		virtual void setHeaderPedestal (EVENT::FloatVec v_ped);
		
		
		//! return pedestal values stored in run header
		inline EVENT::FloatVec getHeaderPedestal () const    {
			EVENT::FloatVec v_ped;
			v_ped = _lcHeader->parameters().getFloatVals (ALIBAVA::HEADERPEDESTAL,v_ped);
			return v_ped;
		}

		//! Set the pedestal values stored in run header
		/*! std::vector<int> v_ped stores the values
		 *
		 */
		virtual void setHeaderNoise ( EVENT::FloatVec v_noise);
		
		
		//! return pedestal values stored in run header
		inline EVENT::FloatVec getHeaderNoise () const    {
			EVENT::FloatVec v_noise;
			v_noise = _lcHeader->parameters().getFloatVals (ALIBAVA::HEADERNOISE, v_noise);
			return v_noise;
		}

		//! Set tilt angle of the sensors
		virtual void setTiltAngle (float tiltAngle);
		
		//! return tilt angle stored in run header
		inline float getTiltAngle () const    {
			return _lcHeader->parameters().getFloatVal (ALIBAVA::TILTANGLE);
		}
		
		//! Set temperature of the sensors
		virtual void setSensorTemperature (float sensorTemperature);
		
		//! return temperature ot the sensors stored in run header
		inline float getSensorTemperature () const    {
			return _lcHeader->parameters().getFloatVal (ALIBAVA::SENSORTEMPERATURE);
		}
		
		//! Set selected chip numbers
		virtual void setChipSelection (EVENT::IntVec achipselectionvec);
		
		//! return selected chip numbers
		inline EVENT::IntVec getChipSelection () const    {
			EVENT::IntVec chipSelec;
			_lcHeader->parameters().getIntVals (ALIBAVA::SELECTEDCHIPNUM,chipSelec);
			return chipSelec;
		}

		//! return number of selected chips
		inline int getNChips () const    {
			EVENT::IntVec chipSelec = getChipSelection();
			return chipSelec.size();
		}
		
		
		//! Add a processor to the applied processor list
		/*! The analysis procedure of some input data usually requires
		 *  that many processors have been applied sequentially. Saving
		 *  this list it may be of interest when trying to reconstruct the
		 *  history of a file.
		 *
		 *  @param processor The name of the process to add to the list
		 */
		virtual void addProcessor (std::string processor);
		
		
		
    //! returns the LCRunHeaderImpl underlying object
    inline IMPL::LCRunHeaderImpl * lcRunHeader() { return  _lcHeader ; }

  private:
  #ifndef DISALLOW_COPY_AND_ASSIGN
      //Following #define stops the accidental creation of a copy or assignment operator by causing a link error. Copy and Assignment operators not allowed because they are unnecessary and the cause of many bugs
  #define DISALLOW_COPY_AND_ASSIGN(AlibavaRunHeaderImpl) \
  AlibavaRunHeaderImpl(const AlibavaRunHeaderImpl&); \
  void operator=(const AlibavaRunHeaderImpl&);

  //Private Functions
  DISALLOW_COPY_AND_ASSIGN(AlibavaRunHeaderImpl)//See #define just above
  #endif

    IMPL::LCRunHeaderImpl * _lcHeader;

  };                           // end of AlibavaRunHeaderImpl
}     //end of alibava namespace
#endif // ALIBAVARUNHEADERIMPL
