/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELAPIXPARSEPIXEL_H
#define EUTELAPIXPARSEPIXEL_H

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
  
  class EUTelAPIXSparsePixel  :  public EUTelBaseSparsePixel  {
        
  public:
    
    //! Default constructor with all arguments
    EUTelAPIXSparsePixel(short xCoord, short yCoord, short signal, short chip, short time); 

    //! Default constructor with no args
    EUTelAPIXSparsePixel(); 
  
    //! Copy constructor
    EUTelAPIXSparsePixel(const EUTelAPIXSparsePixel &orig);
    
    //! Default destructor
    virtual ~EUTelAPIXSparsePixel() { ; } 

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

    //! Setter for the chipide
    void setChip(short chip) { _chip = chip; }

    //! Setter for the time
    void setTime(short time) { _time = time ; }

    //! Getter for the x coordinate
    short getXCoord() const { return _xCoord ; } 

    //! Getter for the y coordinate
    short getYCoord() const { return _yCoord ; } 

    //! Getter for the signal
    //should be short, kept equal to SimpleSparsePixel for compatibility
    float getSignal() const { return static_cast<float> (_signal) ; } 

    //! Getter for the chip
    short getChip() const { return _chip ; }

    //! Getter for the Time
    float getTime() const { return static_cast<float> (_time) ; } 

  private:
  
    //! The x coordinate
    short _xCoord;
    
    //! The y coordinate
    short _yCoord;

    //! The signal
    short _signal;
    
    //! The chip
    short _chip;

    //! The time
    short _time;

  };
}

#endif
