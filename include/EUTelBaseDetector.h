// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELBASEDETECTOR_H
#define EUTELBASEDETECTOR_H

// personal includes ".h"
#include "EUTELESCOPE.h"

// lcio includes <.h>

// system includes <>
#include <iostream>
#include <vector>
#include <string>

namespace eutelescope {


  //! Virtual class to describe detector in the EUTelescope framework
  /*! 
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id: EUTelBaseDetector.h,v 1.1 2008-08-06 20:37:00 bulgheroni Exp $
   */

  class EUTelBaseDetector {

  public:
    //! Default constructor
    EUTelBaseDetector() { } 

    //! Default destructor
    virtual ~EUTelBaseDetector() {;}

    //! Get the first pixel along x
    virtual unsigned short getXMin() const = 0;

    //! Get the first pixel along y
    virtual unsigned short getyMin() const = 0;

    //! Get the last pixel along x
    virtual unsigned short getXMax() const = 0;

    //! Get the last pixel along y
    virtual unsigned short getyMax() const = 0;

    //! Get the no of pixel along X
    virtual unsigned short getXNoOfPixel() const = 0;

    //! Get the no of pixel along Y
    virtual unsigned short getYNoOfPixel() const = 0;    
    
    //! Get the pixel pitch along X
    virtual float getXPitch() const = 0; 

    //! Get the pixel pitch along Y
    virtual float getYPitch() const = 0; 
 
    //! Get signal polarity
    virtual short getSignalPolarity() const = 0;

    //! Get detector name
    virtual std::string getDetectorName() const = 0;

    //! Get RO mode
    virtual std::string getMode() const = 0;

    //! Get marker position
    virtual std::vector< size_t > getMarkerPosition() const = 0;

    //! Has marker? 
    virtual bool hasMarker() const = 0 ;

    //! Print
    /*! This method is used to print out the detector
     * 
     *  @param os The input output stream
     */
    virtual void print(std::ostream& os) const                            = 0;
      
    //! Overload of operator<<
    /*! This friend function is the overload of the operator << for
     *  the detector base class. It uses the print method that is
     *  virtually defined for all cluster subclasses.
     *
     *  @param os The input output stream as modified by the print
     *  method
     *  @param clu The detector to be stream out
     *  @return The output stream
     */ 
    friend std::ostream& operator<< (std::ostream& os , const EUTelBaseDetector & clu )  { clu.print(os); return os; }

  protected:
    
    // data members
    
    //! This is the detector name!
    std::string _name;
    
    //! The first pixel along x
    unsigned short _xMin;

    //! The last pixel along x
    unsigned short _xMax;

    //! The first pixel along y
    unsigned short _yMin;

    //! The last pixel along y
    unsigned short _yMax;

    //! Picth along x in mm as usual
    float _xPitch; 

    //! Picth along y in mm as usual
    float _yPitch;

    //! The signal polarity
    short _signalPolarity;

    //! Marker position in cols number
    std::vector< size_t > _markerPos;

    //! This is the detector RO mode
    std::string _mode;

  };

}

#endif
