// Version: $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELAPIXTIMING_H
#define EUTELAPIXTIMING_H

// personal includes ".h"
#include "EUTelTiming.h"
#include "EUTELESCOPE.h"

// marlin includes ".h"

// lcio includes <.h>

// system includes <>
#include <stdint.h>
namespace eutelescope {

//! Helper class for simple APIX Timing informations
/*! This class contains only the timing informations (TLU and Realtime) as integer numbers.
*
*  @author Georg Troska <mailto:georg.troska@uni-dortmund.de>

*/

class EUTelAPIXTiming  :  public EUTelTiming  {

	public:

	//! Default constructor with all arguments
	EUTelAPIXTiming(short sensorID, uint64_t realtime, uint32_t tlu_counter, uint32_t tpll_counter );

	//! Default constructor with no args
	EUTelAPIXTiming();

	//! Copy constructor
	EUTelAPIXTiming(const EUTelAPIXTiming &orig);

	//! Default destructor
	virtual ~EUTelAPIXTiming() { ; }

	//! Get the number of elements in the data structure
	/*! This method returns the number of elements the sparse pixel
	*  contains.
	*
	*  @return The number of elements in the data structure
	*/
	virtual unsigned int getNoOfElements() const;

	//! Get the detector type using the enumerator
	/*! This methods returns the sparse pixel type using the
	*  enumerator defined in EUTELESCOPE.h
	*
	*  @return The detector type using the enumerator
	*/
	virtual EUTelDetectorType getDetectorType() const;

	//! Print method
	/*! This method is used to print out the contents of the sparse
	*  pixel
	*
	*  @param os The input output stream
	*/
	virtual void print(std::ostream& os) const ;

	//! Setter for the SensorID
	void setSensorID(short sensorID) { _sensorID = sensorID; }

	//! Setter for the TLU-Couter
	void setTLUCounter(uint32_t tluCounter) { _tluCounter = tluCounter; }

	//! Setter for the Realtime-Counter
	void setRealtime(uint64_t realtime) { _realtimeSec = realtime/1000000000; _realtimeNs = realtime-(_realtimeSec*1000000000);}

	//! Setter for the TPLL-Counter
	void setTPLLCounter(uint32_t tpllCounter) { _tpllCounter = tpllCounter; }

	//! Getter for the SensorID
	short getSensorID() const { return _sensorID; }

	//! Getter for the TLU-Couter
	uint32_t getTLUCounter() const { return _tluCounter; }

	//! Getter for the Realtime-Counter
	uint64_t getRealtime() const { return _realtimeSec*1000000000+_realtimeNs;}

	//! Getter for the TPLL-Counter
	uint32_t getTPLLCounter() const { return _tpllCounter; }



	private:

	//! The sensorID
	short _sensorID;

	//! The TLU Counter
	uint32_t _tluCounter;

	//! The Realtime in seconds
	uint32_t _realtimeSec;

	//!The Realtime in subsecond-ticks ... in ns
	uint32_t _realtimeNs;

	//! The TPLL Counter
	short _tpllCounter;

	EUTelDetectorType _type;

};
} // Of Namespace

#endif
