/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELREFERENCEHIT_H
#define EUTELREFERENCEHIT_H

#define REFHIT_CONST_MAX_SIZE 20

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/LCGenericObjectImpl.h>

// system includes <>
#include <string>

namespace eutelescope {

  //! ReferenceHit class for the EUTelescope package
  /*! The aim of this class is to store into a LCGenericObject the
   *  hit nment contants obtained by the execution of EUTelMille +
   *  pede. 
   *
   *  @author Igor Rubinskiy, DESY <mailto:rubinsky@mail.desy.de>
   *  @version $Id$
   */ 
  class EUTelReferenceHit : public IMPL::LCGenericObjectImpl {

  public:
    //! Default constructor
    /*! This is the default constructor of the alignment constants
     *  
     */
    EUTelReferenceHit();

    //! Constructor with all the needed parameters
    /*! This constructor is the recommended one because it set all the
     *  needed parameters in one go. 
     *
     *  @param sensorID The sensor number as from the GEAR file
     *  @param xOff Offset of the sensor along the x direction
     *  @param yOff Offset of the sensor along the y direction
     *  @param zOff Offset of the sensor along the z direction
     *  @param alpha Angle of the sensor along x
     *  @param beta  Angle of the sensor along y
     *  @param gamma Angle of the sensor along z
     *
     */ 
    EUTelReferenceHit( int sensorID, 
       double xOff,   double yOff,   double zOff,
       double alpha, double beta, double gamma );
      
    //! Default destructor
    virtual ~EUTelReferenceHit() { /* NO-OP */ ; }

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

    //! Print the output
    /*! This method is used to print out the constant 
     *
     *  @param os The input output stream
     */
    virtual void print(std::ostream& os) const ;

    //! Overload of operator<<
    /*! This friend function is the overload of the operator << for
     *  the EUTelReferenceHit
     *
     *  @param os The input output stream as modified by the print
     *  method
     *  @param c The alignment constant
     *  @return The output stream
     *
     */ 
    friend std::ostream& operator<< (std::ostream& os, const EUTelReferenceHit & c) { c.print(os); return os; }

  };

  

}

#endif
