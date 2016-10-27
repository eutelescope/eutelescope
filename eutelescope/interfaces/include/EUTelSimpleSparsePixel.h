/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELSIMPLESPARSEPIXEL_H
#define EUTELSIMPLESPARSEPIXEL_H

// personal includes ".h"
#include "EUTelBaseSparsePixel.h"
#include "EUTELESCOPE.h"

// marlin includes ".h"

// lcio includes <.h>

// system includes <>

namespace eutelescope {

  //! Helper class for simple sparsified pixel
  /*! This class contains only the pixel coordinates and signal as
   *  integer numbers.
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
   */ 

  class EUTelSimpleSparsePixel  :  public EUTelBaseSparsePixel  {

  public:
    
    //! Default constructor with all arguments
    EUTelSimpleSparsePixel(short xCoord, short yCoord, float signal); 

    //! Default constructor with no args
    EUTelSimpleSparsePixel(); 
    
    //! Default destructor
    virtual ~EUTelSimpleSparsePixel() { ; } 

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

    //! Setter for x coordinate
    void setXCoord(short xCoord) { _xCoord = xCoord ; }

    //! Setter for y coordinate
    void setYCoord(short yCoord) { _yCoord = yCoord ; }

    //! Setter for the signal
    void setSignal(float signal) { _signal = signal ; }

    //! Getter for the x coordinate
    inline short getXCoord() const { return _xCoord ; } 

    //! Getter for the y coordinate
    inline short getYCoord() const { return _yCoord ; } 

    //! Getter for the signal
    inline float getSignal() const { return _signal ; } 


  private:
  
    //! The x coordinate
    short _xCoord;
    
    //! The y coordinate
    short _yCoord;

    //! The signal
    float _signal;

  };
}

#endif
