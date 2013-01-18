/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELEXCEPTIONS_H
#define EUTELEXCEPTIONS_H 1

// marlin includes
#ifdef USE_MARLIN
#include <marlin/Processor.h>
#endif

// lcio includes
#include <lcio.h>
#include <Exceptions.h>

// system includes
#include <string>


namespace eutelescope {

  //! Invalid parameter
  /*! This exception is thrown when a parameter is available in the
   *  run header / event or collection but it is not accepted. For *
   *  example, in EUTelPedestalNoiseProcess::_pedestalAlgo is set from
   *  a parameter given in the steering file. If the user select an
   *  invalid algorithm then this exceptions is thrown
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
   */

  class InvalidParameterException  : public lcio::Exception {

  protected:
    //! Default constructor
    InvalidParameterException() { /* NO-OP */ ;}
  public:
    //! Default destructor
    virtual ~InvalidParameterException() throw() { /* NO-OP */ ;}

    //! Default constructor with string argument
    /*! @param text Message shown by what().
     */
    InvalidParameterException( const std::string&  text ) {
      message = "eutelescope::InvalidParameterException: " + text;
    }
  };


  //! Incompatible data set
  /*! This exception is thrown all the times two data set are compared
   *  and not found to be compatible. A typical example is when you
   *  try to subtract a pedestal collection from a raw data
   *  collection acquired using two configurations that differs in
   *  the number of detectors.
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
   */
  class IncompatibleDataSetException : public lcio::Exception {

  protected:
    //! Default constructor
    IncompatibleDataSetException() { /* NO-OP */ ; }

  public:
    //! Default destructor
    virtual ~IncompatibleDataSetException() throw() { /* NO - OP */ ; }

    //! Default constructor with string argument
    /*! @param text Message shown by what().
     */
    IncompatibleDataSetException(const std::string& text) {
      message = "eutelescope::IncompatibleDataSetException: " + text;
    }
  };

  //! Unknown data type
  /*! This exception is thrown when the code is not able to understand
   *  the type of an object.
   *
   *  A typical example is the case of clusters. In the EUTelescope
   *  framework there is a virtual base class to describe a cluster
   *  inheriting from IMPL::TrackerDataImpl. There are several
   *  different way this virtual base class can be implemented
   *  according to the way the cluster was found (clustering
   *  algorithm) or the initial data were available (zero suppressed
   *  or not).
   *
   *  Clusters are stored within a TrackerPulse object, but the
   *  corresponding TrackerData needs to re-interpret according to the
   *  cluster type originally used. This is coded into the
   *  TrackerPulse cellID (ClusterType) with an enumeration
   *  (eutelescope::ClusterType).
   *
   *  Of course, if this information is missing the rest of the
   *  analysis procedure might not be able to continue and this
   *  exception is thrown.
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
   */
  class UnknownDataTypeException : public lcio::Exception {

  protected:
    //! Default constructor
    UnknownDataTypeException() { /* NO - OP */ ; }

  public:
    //! Default destructor
    virtual ~UnknownDataTypeException() throw() { /* NO - OP */ ; }

    //! Default constructor with string argument
    /*! This is the standard way to create this exception.
     *
     *  @param text Message to be shown by what()
     */
    UnknownDataTypeException(const std::string& text) {
      message = "eutelescope::UnknownDataTypeException: " + text;
    }
  };

#ifdef USE_MARLIN
  //! Missing library
  /*! This exception is thrown by a processor having optional
   *  dependency against specific library and this is missing.
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
   */
  class MissingLibraryException: public lcio::Exception {

  protected:
    //! Default constructor
    MissingLibraryException() { /* NO - OP */ ; }

  public:
    //! Default destructor
    virtual ~MissingLibraryException() throw() { /* NO - OP */ ; }

    //! Default constructor with a processor and a string argument
    /*! This is the standard way to create this exception.
     *
     *  @param proc This is the processor throwing the exception
     *  @param missingLib This is the name of missing library
     */
    MissingLibraryException( marlin::Processor * proc, std::string missingLib ) {
      message = "To execute this functionality " + proc->name() + " requires " + missingLib;
    }
  };
#endif

  //! Invalid geometry
  /*! This exception is thrown when the processor is expecting to have
   *  a different telescope geometry compared to the one described in
   *  the GEAR file
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
   */
  class InvalidGeometryException: public lcio::Exception {

  public:
    //! Default constructor
    /*! @param text The error message
     */
    InvalidGeometryException(std::string text)  {
      message = text ;
    }

    //! Default destructor
    ~InvalidGeometryException() throw()  { /* NO - OP */ ; }
  };


}

#endif
