/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELGENERICSPARSEPIXEL_H
#define EUTELGENERICSPARSEPIXEL_H

// personal includes ".h"
#include "EUTelBaseSparsePixel.h"
#include "EUTELESCOPE.h"

namespace eutelescope {

  //! Helper class for simple sparsified pixel
  /*! This class contains the pixel coordinates, signal and
   *  time as integer numbers.
   *
   *  Based on the EUTelSimpeSparsePixel class by Antonio Bulgheroni
   */ 

class EUTelGenericSparsePixel : public EUTelBaseSparsePixel  {

public:
    //! Default constructor with all arguments
    EUTelGenericSparsePixel(short xCoord, short yCoord, short signal, short time); 

    //! Default constructor with time argument omitted
    /*! Time will be set to zero */
    EUTelGenericSparsePixel(short xCoord, short yCoord, short signal); 

    //! Default constructor with no args
    /*! Every value will be set to zero */
    EUTelGenericSparsePixel(); 
    
    //! Destructor
    virtual ~EUTelGenericSparsePixel() { ; } 

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
    void setSignal(short signal) { _signal = signal ; }

    //! Setter for the time
    void setTime(short time) { _time = time ; }

    //! Getter for the x coordinate
    inline short getXCoord() const { return _xCoord ; } 

    //! Getter for the y coordinate
    inline short getYCoord() const { return _yCoord ; } 

    //! Getter for the signal
    inline float getSignal() const { return static_cast<float> (_signal); } 

    //! Getter for the time
    inline float getTime() const { return static_cast<float> (_time); } 

protected:
     
    //! The x coordinate
    short _xCoord;
    
    //! The y coordinate
    short _yCoord;

    //! The signal
    int _signal;

	//! The time
    short _time;

};

} //namespace eutelescope
#endif
