/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELBINARYPIXEL_H
#define EUTELBINARYPIXEL_H

// personal includes ".h"
#include "EUTelGenericSparsePixel.h"
#include "EUTELESCOPE.h"

// marlin includes ".h"

// lcio includes <.h>

// system includes <>

namespace eutelescope {

  //! Helper class decribing a binary pixel
  /*! This class inherits from EUTelGenericSparsePixel, but 
   *  additionally stores information on the timestamp of individual hits
   *  and of events.
   *  If only the generic information is required, the
   *  EUTelMuPixel can be easily cast to an
   *  EUTelGenericSparsePixel, which then can be stored in
   *  a collection.
   */ 

class EUTelMuPixel : public EUTelGenericSparsePixel  {

public:
    //! Default constructor with all arguments (individually)
  EUTelMuPixel(short xCoord, short yCoord, float signal, short time, short hit_time, long long unsigned frame_time); 

    //! Default constructor with all arguments (EUTelGenericSparsePixel, 2 additional time stamps)
  EUTelMuPixel(EUTelGenericSparsePixel& genericPixel, short hit_time, long long unsigned frame_time); 
  
  //! Default constructor with only arguments of a EUTelescopeGenericSparsePixel (geometry related values are set to 0)
  EUTelMuPixel(EUTelGenericSparsePixel& genericPixel); 
  
  //! Default constructor with no args (all values are set to 0)
  EUTelMuPixel(); 
    
  //! Destructor
  virtual ~EUTelMuPixel() {} 
  
  //! Get the number of elements in the data structure
  /*! This method returns the number of elements the sparse pixel
   *  contains. 
   *
   *  @return The number of elements in the data structure
   */
  virtual unsigned int getNoOfElements() const;
  
  //! Get the sparse pixel type using the enumerator
  /*! This methods returns the sparse pixel type using the
   *  enumerator defined in EUTELESCOPE.h
   *
   *  Overloaded for derived class, since downcast should
   *  yield the sparse pixel type of the base class.
   *
   *  @return The sparse pixel type using the enumerator
   */
  virtual SparsePixelType getSparsePixelType() const;

  //! Print method
  /*! This method is used to print out the contents of the sparse
   *  pixel
   *
   *  @param os The input output stream
   */
  virtual void print(std::ostream& os) const ;
  
  //! Setter for the hit time stamp (double precision)
  void setHitTime(short hitTime) { _hitTime = hitTime; }

  //! Setter for the frame time stamp
  void setFrameTime(long long unsigned frameTime) { _frameTime = frameTime; } 

  //! Getter for the hit time stamp
  inline short getHitTime() const { return _hitTime; }

  //! Getter for the frame time stamp
  inline long long unsigned getFrameTime() const { return _frameTime; }
  
 protected:
  //! hit time stamp
  short _hitTime;
  
  //! frame time stamp
  long long unsigned _frameTime;

  //! The number of elements in the data structure
  unsigned int _noOfElementsDerived;
  
  //! The sparse pixel type enumerator for the derived type
  /*! Required since a downcast to EUTelGenericSparsePixel
   *  needs to return the downcast type and member variable
   *  overloading is not possible in C++.
   */
  SparsePixelType _typeDerived;
};
} //namespace eutelescope

#endif
