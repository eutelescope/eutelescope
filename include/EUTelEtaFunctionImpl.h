/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELETAFUNCTIONIMPL_H
#define EUTELETAFUNCTIONIMPL_H

#define ETA_VERSION 2

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/LCGenericObjectImpl.h>

// system includes <>
#include <string>

namespace eutelescope {

  //! Eta function LCIO implementation
  /*! The Eta function calculated by the EUTelCalculateEtaProcessor is
   *  saved into a Condition file using this implementation.
   *  EUTelEtaFunctionImpl is nothing more than a sub-class of
   *  LCGenericObjectImpl where the integer and float parts are masked
   *  behind the private label. The number of double entries is twice
   *  the number of bins: the first @a nBin are the bin centers (x
   *  axis) while the second @a nBin are the corresponding Eta values.
   *
   *  Some convenience methods have been introduced to allow the user
   *  to set/get the full bin center and the full eta value vectors
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
   */
  class EUTelEtaFunctionImpl : public IMPL::LCGenericObjectImpl {

  public:
    //! Default constructor
    /*! This is the default constructor of the Eta function LCIO
     *  implementation. 
     *  
     *  @param nBin This is the number of bins used in the calculation
     *  of the Eta function.
     */
    EUTelEtaFunctionImpl(int nBin); 

    //! Constructor with all the needed parameters
    /*! This convenience constructor can be used to create the object
     *  and immediately set both the bin center and the value vectors
     *
     *  @param nBin This is the number of bins used in the calculation of the Eta function
     *  @param centerVec This is a STL vector of double with the bin centers.
     *  @param valueVec This is a STL vector of double with the eta values
     */ 
    EUTelEtaFunctionImpl(int nBin, std::vector<double > centerVec, std::vector<double > valueVec);


    //! Constructor with all the needed parameters
    /*! This convenience constructor can be used to create the object
     *  and immediately set both the bin center and the value vectors
     *
     *  @param sensorID The is the sensorID of this sensor.
     *  @param nBin This is the number of bins used in the calculation of the Eta function
     *  @param centerVec This is a STL vector of double with the bin centers.
     *  @param valueVec This is a STL vector of double with the eta values
     */ 
    EUTelEtaFunctionImpl(int sensorID, int nBin, std::vector<double > centerVec, std::vector<double > valueVec);

    //! Default destructor
    virtual ~EUTelEtaFunctionImpl() { /* NO-OP */ ; }

    //! Set the bin center vector
    /*! This method can be used to set all bin centers in one call
     *  @param center A STL vector of double with all the bin centers
     */ 
    void setBinCenterVector(std::vector<double > center);

    //! Set the eta value vector
    /*! This method can be used to set all eta values in one call
     * 
     *  @param value A STL vector of double with all the eta values
     */ 
    void setEtaValueVector(std::vector<double > value);

    //! Get the number of bin
    /*! This is used to get the number of bin in the function.
     */
    int getNoOfBin() const;

    //! Get the bin center vector
    /*! @return A STL vector of double with the bin centers
     */ 
    const std::vector<double > getBinCenterVector() const;

    //! Get the eta value vector
    /*! @return A STL vector of double with the eta values
     */ 
    const std::vector<double > getEtaValueVector() const;

    //! Set the sensor ID
    /*! @param sensorID The sensor ID
     */
    void setSensorID( int sensorID ) ;

    //! Get the sensor ID 
    /*! @return the sensor ID
     */
    int getSensorID( ) const; 


    //! Get Eta for a given CoG value
    /*! When applying the Eta correction to the charge center of
     *  gravity, the user must find Eta(CoG). This convenience
     *  function is used for that purpose. This algorithm is based on
     *  a binary search of the center of gravity value @a x using the
     *  C++ implementation of lower_bound. Once the binning which @a x
     *  belongs to is found, the returned value of @a eta is linearly
     *  interpolated.
     *
     *  \li Note 1: It is possible to use the lower_bound algorithm
     *  because the CoG vector is sorted by definition.
     *
     *  \li Note 2: In the case the x axis binning during the Eta
     *  function calculation was kept constant and uniform, then the
     *  binary search algorthim can be replaced with a much faster
     *  calculation of the x closest pair of values. Because of this
     *  lack in generality, we prefer to invest some calculation power
     *  in the binary search.
     *
     *  @param x is the current CoG value
     *  @return the corresponding Eta value
     *
     */ 
    double getEtaFromCoG(double x) const ;

  protected:
    
    //! Get the begin iterator for the CoG vector
    /*! This method is used to get an iterator corresponding to the
     *  beginning of the charge center of gravity vector.
     *
     *  @return a iterator to the beginning of the CoG vector
     */ 
    std::vector<double >::const_iterator getCoGBeginConstIterator() const;

    //! Get the end iterator for the CoG vector
    /*! This method is used to get an iterator corresponding to the
     *  end of the charge center of gravity vector.
     *
     *  @return a iterator to the end of the CoG vector
     */ 
    std::vector<double >::const_iterator getCoGEndConstIterator() const;


    //! Get the begin iterator for the Eta value vector
    /*! This method is used to get an iterator corresponding to the
     *  beginning of the charge center of gravity vector.
     *
     *  @return a iterator to the beginning of the Eta value vector
     */ 
    std::vector<double >::const_iterator getEtaBeginConstIterator() const;

    //! Get the end iterator for the Eta value vector
    /*! This method is used to get an iterator corresponding to the
     *  end of the charge center of gravity vector.
     *
     *  @return a iterator to the end of the Eta value vector
     */ 
    std::vector<double >::const_iterator getEtaEndConstIterator() const;

  private:

    void getNFloat() {;}
    void getFloatVal() {;}
    void setFloatVal(unsigned int, float) {;}
    void isFixedSize() {;}

  };

}

#endif
