// Version: $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELTIMING_H
#define EUTELTIMING_H

// personal includes ".h"
#include "EUTELESCOPE.h"

// marlin includes ".h"

// lcio includes <.h>

// system includes <>
#include <iostream>

namespace eutelescope {

  //! Abstract base class for timing information
  /*! This is a pure abstract class containing only the definition of
   *  methods available to all timing informations of differnet plane-types
   *  
   *  @author Georg Troska  <mailto:georg.troska@uni-dortmund.de>
   */
  class EUTelTiming {

    public:
    //! Default constructor
    EUTelTiming() { } 

    //! Default destructor
    virtual ~EUTelTiming() { ; } 

    //! Get the number of elements in the data structure
    /*! This method returns the number of elements the timing
     *  contains. 
     *
     *  @return The number of elements in the data structure
     */
    virtual unsigned int getNoOfElements() const = 0;

    //! Get the detector type using the enumerator
    /*! This methods returns the type using the
     *  enumerator defined in EUTELESCOPE.h
     *
     *  @return The detector type using the enumerator
     */
    virtual EUTelDetectorType getDetectorType() const = 0;

    //! Print method
    /*! This method is used to print out the contents of the sparse pixel
     *
     *  @param os The input output stream
     */
    virtual void print(std::ostream& os) const = 0 ;
	
    //! Overload of operator\<\<
    /*! This friend function is the overload of the operator << for
     *  the base sparse pixel class. It uses the print method that is
     *  virtually defined for all sparse pixel subclasses.
     *  
     *  @param os The input output stream as modified by the print
     *  method
     *  @param timing The base pixel timing to be streamed out
     *  @return The output stream
     */
    friend std::ostream& operator<< (std::ostream& os, const EUTelTiming& timing)  { timing.print(os); return os; }

  protected:

    //! The number of elements in the data structure
    unsigned int _noOfElements;

    //! The sparse pixel type enumerator
    EUTelDetectorType _type;


};

} // of namespace
#endif
