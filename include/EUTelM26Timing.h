// Version: $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELM26TIMING_H
#define EUTELM26TIMING_H

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

class EUTelM26Timing  :  public EUTelTiming  {

	public:

	//! Default constructor with all arguments
	EUTelM26Timing(short sensorID, uint32_t tlu_counter);

	//! Default constructor with no args
	EUTelM26Timing();

	//! Copy constructor
	EUTelM26Timing(const EUTelM26Timing &orig);

	//! Default destructor
	virtual ~EUTelM26Timing() { ; }

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


	//! Getter for the SensorID
	short getSensorID() const { return _sensorID; }

	//! Getter for the TLU-Couter
	uint32_t getTLUCounter() const { return _tluCounter; }




	private:

	//! The sensorID
	short _sensorID;

	//! The TLU Counter
	uint32_t _tluCounter;

	EUTelDetectorType _type;

};
} // Of Namespace

#endif
