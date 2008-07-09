// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELALIGNMENTCONSTAN_H
#define EUTELALIGNMENTCONSTAN_H

#define ALIGN_CONST_MAX_SIZE 20

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/LCGenericObjectImpl.h>

// system includes <>
#include <string>

namespace eutelescope {

  //! Alignment constant for the EUTelescope package
  /*! The aim of this class is to store into a LCGenericObject the
   *  alignment contants obtained by the execution of EUTelMille +
   *  pede. 
   *
   *  @Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @Version $Id: EUTelAlignmentConstant.h,v 1.1 2008-07-09 13:01:26 bulgheroni Exp $
   */ 
  class EUTelAlignmentConstant : public IMPL::LCGenericObjectImpl {

  public:
    //! Default constructor
    /*! This is the default constructor of the alignment constants
     *  
     */
    EUTelAlignmentConstant(); 

    //! Constructor with all the needed parameters
    /*! This constructor is the recommended one because it set all the
     *  needed parameters in one go. 
     *
     *  @param sensorID The sensor number as from the GEAR file
     *  @param xOff Offset of the sensor along the x direction
     *  @param yOff Offset of the sensor along the y direction
     *  @param zOff Offset of the sensor along the z direction
     *  @param xTheta Angle of the sensor along x
     *  @param yTheta Angle of the sensor along y
     *  @param zTheta Angle of the sensor along z
     *  @param xThetaErr Error on the angle of the sensor along x
     *  @param yThetaErr Error on the angle of the sensor along y
     *  @param zThetaErr Error on the angle of the sensor along z
     *  @param xOffErr Error on the offset of the sensor along the x direction
     *  @param yOffErr Error on the offset of the sensor along the y direction
     *  @param zOffErr Error on the offset of the sensor along the z direction
     *
     */ 
    EUTelAlignmentConstant( int sensorID, 
			    double xOff,   double yOff,   double zOff,
			    double xTheta, double yTheta, double zTheta,
			    double xOffErr,   double yOffErr,   double zOffErr,
			    double xThetaErr, double yThetaErr, double zThetaErr );
      
    //! Default destructor
    virtual ~EUTelAlignmentConstant() { /* NO-OP */ ; }

    //! Set the sensor id
    void setSensorID( int id ) ;

    //! Set the x offset in mm 
    void setXOffset( double off ) ;

    //! Set the y offset in mm 
    void setYOffset( double off ) ;

    //! Set the z offset in mm 
    void setZOffset( double off ) ;

    //! Set the angle around x
    void setXTheta( double theta );

    //! Set the angle around y
    void setYTheta( double theta );

    //! Set the angle around z
    void setZTheta( double theta );

    //! Set the error of the offset along x
    void setXOffsetError( double err ) ;

    //! Set the error of the offset along y
    void setYOffsetError( double err ) ;

    //! Set the error of the offset along z
    void setZOffsetError( double err ) ;

    //! Set the error of the angle around x
    void setXThetaError( double err ) ;

    //! Set the error of the angle around y
    void setYThetaError( double err ) ;

    //! Set the error of the angle around z
    void setZThetaError( double err ) ;

    //! Get the sensor ID
    int getSensorID() const ;
    
    //! Get the offset along x
    double getXOffset() const;

    //! Get the offset along y
    double getYOffset() const;

    //! Get the offset along z
    double getZOffset() const;

    //! Get the theta along x
    double getXTheta() const;

    //! Get the theta along y
    double getYTheta() const;

    //! Get the theta along z
    double getZTheta() const;

    //! Get the error of the offset along x
    double getXOffsetError() const;

    //! Get the error of the offset along y
    double getYOffsetError() const;

    //! Get the error of the offset along z
    double getZOffsetError() const;

    //! Get the error of the theta along x
    double getXThetaError() const;

    //! Get the error of the theta along y
    double getYThetaError() const;

    //! Get the error of the theta along z
    double getZThetaError() const;

    //! Print the output
    /*! This method is used to print out the constant 
     *
     *  @param os The input output stream
     */
    virtual void print(std::ostream& os) const ;

    //! Overload of operator<<
    /*! This friend function is the overload of the operator << for
     *  the EUTelAlignmentConstant
     *
     *  @param os The input output stream as modified by the print
     *  method
     *  @param c The alignment constant
     *  @return The output stream
     *
     */ 
    friend std::ostream& operator<< (std::ostream& os, const EUTelAlignmentConstant & c) { c.print(os); return os; }

  };

  

}

#endif
