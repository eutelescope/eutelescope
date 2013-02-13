/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELBASESPARSEPIXEL_H
#define EUTELBASESPARSEPIXEL_H

// personal includes ".h"
#include "EUTELESCOPE.h"

// marlin includes ".h"

// lcio includes <.h>

// system includes <>
#include <iostream>

namespace eutelescope {

  //! Abstract base class for sparsified pixels
  /*! This is a pure abstract class containing only the definition of
   *  methods available to all sparsified pixel type
   *  
   *  @author Antonio Bulgheroni, INFN  <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
   */
  class EUTelBaseSparsePixel {

    public:
    //! Default constructor
    EUTelBaseSparsePixel() : _noOfElements(0), _type() { } 

    //! Default destructor
    virtual ~EUTelBaseSparsePixel() { ; } 

    //! Get the number of elements in the data structure
    /*! This method returns the number of elements the sparse pixel
     *  contains. 
     *
     *  @return The number of elements in the data structure
     */
    virtual unsigned int getNoOfElements() const = 0;

    //! Get the sparse pixel type using the enumerator
    /*! This methods returns the sparse pixel type using the
     *  enumerator defined in EUTELESCOPE.h
     *
     *  @return The sparse pixel type using the enumerator
     */
    virtual SparsePixelType getSparsePixelType() const = 0;

    //! Get the x pixel coordinates
    /*! This method returns the sparse pixel coordinate along the x
     *  direction.
     * 
     *  @return The x pixel coordinates
     */
     virtual short getXCoord() const = 0;

    //! Get the y pixel coordinates
    /*! This method returns the sparse pixel coordinate along the y
     *  direction.
     * 
     *  @return The y pixel coordinates
     */
     virtual short getYCoord() const = 0;

    //! Get the pixel signal. 
    /*! This method returns the pixel signal always as a float value
     *  even if it is an short integer.
     *
     *  @return The pixel signal
     */
    virtual float getSignal() const = 0;

    //! Print method
    /*! This method is used to print out the contents of the sparse
     *  pixel
     *
     *  @param os The input output stream
     */
    virtual void print(std::ostream& os) const = 0 ;
	
    //! Overload of operator<<
    /*! This friend function is the overload of the operator << for
     *  the base sparse pixel class. It uses the print method that is
     *  virtually defined for all sparse pixel subclasses.
     *  
     *  @param os The input output stream as modified by the print
     *  method
     *  @param pixel The base pixel to be stream out
     *  @return The output stream
     */
    friend std::ostream& operator<< (std::ostream& os, const EUTelBaseSparsePixel& pixel)  { pixel.print(os); return os; }

  protected:

    //! The number of elements in the data structure
    unsigned int _noOfElements;

    //! The sparse pixel type enumerator
    SparsePixelType _type;

  public:

    //! Helper class for signal sorting
    /*! This template class is used as a binary function for pixel
     *  sorting. 
     *
     *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
     *  @version $Id$
     */ 
    template<class T>
    class GreaterSignal :  
      public std::binary_function<T, T, bool> {
      
    public:
      //! Operator () overload
      /*! This is the method used by the sort
       *  
       *  @param pixel1 The first pixel 
       *  @param pixel2 The second pixel
       *  @return True if first pixel has greater signal than pixel2.
       */ 
      bool operator() (const T& pixel1,
		       const T& pixel2 ) const {
	return pixel1.getSignal() >= pixel2.getSignal();
      }
    };
    
    //! Helper class for signal sorting
    /*! This template class is used as a binary function for pixel
     *  sorting. 
     *
     *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
     *  @version $Id$ 
     */ 
    template<class T>
    class SmallerSignal :
      public std::binary_function<T, T, bool> {
      
    public:
      //! Operator () overload
      /*! This is the method used by the sort
       *  
       *  @param pixel1 The first pixel 
       *  @param pixel2 The second pixel
       *  @return True if first pixel has less signal than pixel2.
       */     
      bool operator() ( const T& pixel1,
			const T& pixel2 ) const {
	return pixel1.getSignal() <= pixel2.getSignal();
      }
    };
    
    //! Helper class for position sorting
    /*! This template class is used as a binary function for pixel
     *  sorting. 
     *
     *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
     *  @version $Id$ 
     */     
    template<class T>
    class IsBefore :
      public std::binary_function<T, T, bool > {
      
    public:
      //! Operator () overload
      /*! This is the method used by the sort
       *  
       *  @param pixel1 The first pixel 
       *  @param pixel2 The second pixel
       *  @return True if first pixel is positioned before the second
       *  on the matrix
       */ 
      bool operator() ( const T& pixel1, 
			const T& pixel2 ) const {
	if ( pixel1.getYCoord() < pixel2.getYCoord() ) return true;
	else if ( pixel1.getYCoord() == pixel2.getYCoord() ) return ( pixel1.getXCoord() < pixel2.getXCoord() );
	else return false;
      }
    };

    //! Helper predicate class
    /*! This template class is used as a unary function to find a
     *  pixel in the matrix
     *
     *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
     *  @version $Id$ 
     */     
    template<class T>
    class HasCoord {
    public:
      //! Default constructor
      HasCoord( int x, int y ) : _x(x), _y(y) { ; } 
      
      //! Operator () overload
      /*! This is the method used by the sort
       *  
       *  @param pixel The first pixel 
       *  @return True if the pixel has the coordinate x, y;
       */       
      bool operator() (const T& pixel) 
      { return ( ( pixel.getXCoord() == _x ) && ( pixel.getYCoord() == _y ) ) ; }
      
    private:
      //! Local copy of the x coordinate
      int _x;
      
      //! Local copy of the y coordinate
      int _y;

    };

    //! Helper predicate class
    /*! This template class is used as a unary function to find a
     *  pixel in the matrix
     *
     *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
     *  @version $Id$ 
     */ 
    template<class T>
    class HasXCoord {
    public:
      //! Default constructor
      HasXCoord( int x ) : _x(x) { ; }
      
      //! Operator () overload
      /*! This is the method used by the sort
       *  
       *  @param pixel The first pixel 
       *  @return True if the pixel has the coordinate x.
       */
      bool operator() (const T& pixel )
      { return ( pixel.getXCoord() == _x ); }
    
    private:
      //! Local copy of the x coordinate
      int _x;
    };

    //! Helper predicate class
    /*! This template class is used as a unary function to find a
     *  pixel in the matrix
     *
     *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
     *  @version $Id$ 
     */ 
    template<class T>
    class HasYCoord {
    public:
      //! Default constructor
      HasYCoord( int y ) : _y(y) { ; }

      //! Operator () overload
      /*! This is the method used by the sort
       *  
       *  @param pixel The first pixel 
       *  @return True if the pixel has the coordinate x.
       */
      bool operator() (const T& pixel )
      { return ( pixel.getYCoord() == _y ); }
    
    private:
      //! Local copy of the x coordinate
      int _y;
    };
       
  };

  //! Compute the distance between two pixels
  /*! This function is calculating the distance between two sparsified
   *  pixels (of any type) returning a floating value in pixel unit.
   *
   *  Two pixels sharing one side will have a distance equal to one,
   *  while if they are touching by a corner, they will be sqrt(2) far
   *  a part. 
   *
   */ 
  float distance(EUTelBaseSparsePixel * first, EUTelBaseSparsePixel * second) ;

}
#endif
