/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELALIGNMENTCONSTANT_H
#define EUTELALIGNMENTCONSTANT_H


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
   *  @author Contact: antonio.bulgheroni@gmail.com
   *  @version $Id$
   */ 
  class EUTelAlignmentConstant : public IMPL::LCGenericObjectImpl {

  public:
    //! Default constructor
    /*! This is the default constructor of the alignment constants
     *  
     */
    EUTelAlignmentConstant(); 

    //! Constructor with all the needed parameters
    /*! This constructor is the recommended one because it sets all the
     *  needed parameters in one go. 
     *
     *  @param sensorID The sensor number as from the GEAR file
     *  @param xOff Offset of the sensor along the x direction
     *  @param yOff Offset of the sensor along the y direction
     *  @param zOff Offset of the sensor along the z direction
     *  @param alpha Angle of the sensor along x
     *  @param beta Angle of the sensor along y
     *  @param gamma Angle of the sensor along z
     *  @param alphaErr Error on the angle of the sensor along x
     *  @param betaErr Error on the angle of the sensor along y
     *  @param gammaErr Error on the angle of the sensor along z
     *  @param xOffErr Error on the offset of the sensor along the x direction
     *  @param yOffErr Error on the offset of the sensor along the y direction
     *  @param zOffErr Error on the offset of the sensor along the z direction
     *
     */ 
    EUTelAlignmentConstant( int sensorID, 
       double xOff,   double yOff,   double zOff,
       double alpha, double beta, double gamma,
       double xOffErr,   double yOffErr,   double zOffErr,
       double alphaErr, double betaErr, double gammaErr );
      
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
    void setAlpha( double theta );

    //! Set the angle around y
    void setBeta( double theta );

    //! Set the angle around z
    void setGamma( double theta );

    //! Set the error of the offset along x
    void setXOffsetError( double err ) ;

    //! Set the error of the offset along y
    void setYOffsetError( double err ) ;

    //! Set the error of the offset along z
    void setZOffsetError( double err ) ;

    //! Set the error of the angle around x
    void setAlphaError( double err ) ;

    //! Set the error of the angle around y
    void setBetaError( double err ) ;

    //! Set the error of the angle around z
    void setGammaError( double err ) ;

    //! Get the sensor ID
    int getSensorID() const ;
    
    //! Get the offset along x
    double getXOffset() const;

    //! Get the offset along y
    double getYOffset() const;

    //! Get the offset along z
    double getZOffset() const;

    //! Get the theta along x
    double getAlpha() const;

    //! Get the theta along y
    double getBeta() const;

    //! Get the theta along z
    double getGamma() const;

    //! Get the error of the offset along x
    double getXOffsetError() const;

    //! Get the error of the offset along y
    double getYOffsetError() const;

    //! Get the error of the offset along z
    double getZOffsetError() const;

    //! Get the error of the theta along x
    double getAlphaError() const;

    //! Get the error of the theta along y
    double getBetaError() const;

    //! Get the error of the theta along z
    double getGammaError() const;

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
