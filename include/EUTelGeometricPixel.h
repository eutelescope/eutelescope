/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELGEOMETRICPIXEL_H
#define EUTELGEOMETRICPIXEL_H

// personal includes ".h"
#include "EUTelGenericSparsePixel.h"
#include "EUTELESCOPE.h"

// marlin includes ".h"

// lcio includes <.h>

// system includes <>

namespace eutelescope {

  //! Helper class decribing a geometric pixel
  /*! This class inherits from EUTelGenericSparsePixel, but 
   *  additionally stores information on the exact geometric 
   *  position as well as information of the bounding box 
   *  for each pixel.
   *  It thus takes roughly twice as much space as a 
   *  EUTelGenericSparsePixel and is intended to be only 
   *  used in intermediate steps where this information
   *  is essential.
   *  If only the generic information is required, the
   *  EUTelGeometricPixel can be easily cast to an
   *  EUTelGenericSparsePixel, which then can be stored in
   *  a collection.
   *  The values for the boundary are given in half the
   *  distance in x and y. Thus the position(+/-)boundary
   *  yields the box boundaries in space.
   */ 

class EUTelGeometricPixel : public EUTelGenericSparsePixel  {

public:
    //! Default constructor with all arguments (individually)
    EUTelGeometricPixel(short xCoord, short yCoord, float signal, short time, float posX, float posY, float boundX, float boundY); 

    //! Default constructor with all arguments (EUTelGenericSparsePixel, 4 additional geometry related)
    EUTelGeometricPixel(EUTelGenericSparsePixel& genericPixel, float posX, float posY, float boundX, float boundY); 

    //! Default constructor with only arguments of a EUTelescopeGenericSparsePixel (geometry related values are set to 0)
    EUTelGeometricPixel(EUTelGenericSparsePixel& genericPixel); 

    //! Default constructor with no args (all values are set to 0)
    EUTelGeometricPixel(); 
    
    //! Destructor
    virtual ~EUTelGeometricPixel() {} 

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

    //! Setter for the boundary X
    void setBoundaryX(float boundX) { _boundX = boundX; }

    //! Setter for the boundary Y
    void setBoundaryY(float boundY) { _boundY = boundY; } 

    //! Setter for the x position
    void setPosX(float posX) { _posX = posX; } 

    //! Setter for the y position
    void setPosY(float posY) { _posY = posY; }

    //! Getter for the boundary X
    inline float getBoundaryX() const { return _boundX; }

    //! Getter for the boundary Y
    inline float getBoundaryY() const { return _boundY; }

    //! Getter for the x position
    inline float getPosX() const { return _posX; }

    //! Getter for the y position
    inline float getPosY() const { return _posY; }

protected:
    //! X position on the sensor plane
    float _posX;

    //! Y position on the sensor plane
    float _posY;

    //! Half x-dimension of the boundary
    float _boundX;

    //! Half y-dimension of the boundary
    float _boundY;

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
