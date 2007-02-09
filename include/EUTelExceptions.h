// -*- C++ -*-
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

#include <lcio.h>
#include <Exceptions.h>
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
   *  @version $Id: EUTelExceptions.h,v 1.1 2007-02-09 10:24:44 bulgheroni Exp $
   */

  class InvalidParameterException  : public lcio::Exception {
    
  protected:
    InvalidParameterException() { /* NO-OP */ ;}
  public:
    virtual ~InvalidParameterException() throw() { /* NO-OP */ ;}
    
    InvalidParameterException( std::string text ) {
      message = "eutelescope::InvalidParameterException: " + text;
    }
  };
}

#endif
