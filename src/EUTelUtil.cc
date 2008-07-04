// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelUtil.cc,v 1.1 2008-07-04 10:02:15 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelUtil.h"

// marlin include
#include <streamlog/streamlog.h>

// system includes
#include <string>

#ifndef MAX_LINE_SIZE
#define MAX_LINE_SIZE 4096
#endif


using namespace std;
using namespace eutelescope;


SystemCommand::SystemCommand(string command) {
  
  _command = command;

}

bool SystemCommand::execute() {

  FILE            * fpOutputFile;
  FILE            * fpErrorFile;

  string            errorFileName = "/tmp/redirect.log";

  const size_t      SIZEBUF = MAX_LINE_SIZE;
  char              buf[SIZEBUF];

  const int         fdStdError = 2;
  int               fdNewStdError;

  _stdOutputVec.clear();
  _stdErrorVec.clear();


  // taking care of the stderror first
  // copy the stderr to a new fd to allow stderr restore at the end
  if (( fdNewStdError = dup ( fdStdError ) ) == -1 ) {
    if( streamlog::out.write<streamlog::ERROR>() ){
      streamlog::out() << "SystemCommand::execute: Unable to duplicate the stderr " << _command;
      return false;
    }
  }

  // close the stderr and reopen it to our redirection file
  if (( fpErrorFile = freopen( errorFileName.c_str(),"w", stderr ) ) == NULL ) {
    if( streamlog::out.write<streamlog::ERROR>() ){
      streamlog::out() << "SystemCommand::execute: Unable to open a temporary file" << _command;
      return false;
    }
  }
  
  // from now on all messages directed to the stderr will end up in
  // our fpErrorFile.


  // now we can really process the command
  if ( ( fpOutputFile = popen ( _command.c_str(), "r" ) ) == NULL ) {
    
    // something went wrong with the creation of the new pipe: either
    // the pipe couldn't be opened or new processes couldn't be
    // created. Advise the user and return false

    if( streamlog::out.write<streamlog::ERROR>() ){
      streamlog::out() << "SystemCommand::execute: Files or processes cannot be created " << _command;
      return false;
    }

  }

  // the command is executed and now we just have to get the std
  // output and error from the opened pipe.
  string currentString;
  while ( fgets( buf, sizeof ( buf ) , fpOutputFile ) ) {
    
    // convert char to string for better and easier handling
    currentString = buf;

    // check if the current string is ending with a new line. If not,
    // this means that we are actually not getting enough character
    // for this line and the user should increase MAX_LINE_SIZE
    if ( currentString [currentString.size() - 1] != '\n' ) {
      if( streamlog::out.write<streamlog::ERROR>() ){
	streamlog::out() << "SystemCommand::execute: MAX_SILE_SIZE is too small, try recompile the code using MAX_LINE_SIZE = " 
			 << 2 * MAX_LINE_SIZE;	
	return false;
      }
    
    }
  
    // ok add the currentString to the output vector removing the last
    // character that is an end of line. This is making the future
    // manipulation of the string much easier
    _stdOutputVec.push_back( currentString.substr(0, currentString.size() - 1) );
     
  }
  
  // the std output pipe is now empty, so close it! 
  if ( pclose( fpOutputFile ) == -1 ) {
    // we have problem closing the pipe... strange but not
    // impossible...
    if( streamlog::out.write<streamlog::ERROR>() ){
      streamlog::out() << "SystemCommand::execute: Cannot close the pipe " << _command;
      return false;
    }
  }

  // now we need to restore the stderr.
  // first move fdNewStdError to the old stderror
  if ( dup2( fdNewStdError, fdStdError ) == -1 ) {
    if( streamlog::out.write<streamlog::ERROR>() ){
      streamlog::out() << "SystemCommand::execute: Unable to dup2 the new std error to the old one " << _command;
      return false;
    }
  }

  // now close the new stderror because we don't need it anymore
  if ( close( fdNewStdError ) == -1 ) {
    if( streamlog::out.write<streamlog::ERROR>() ){
      streamlog::out() << "SystemCommand::execute: Unable to close the redirected std error " << _command;
      return false;
    }
  }

  // we need to reopen the file containing the redirected std error in
  // reading mode and get all the lines into it
  if ( freopen( errorFileName.c_str(), "r", fpErrorFile ) == NULL ) {
    if( streamlog::out.write<streamlog::ERROR>() ){
      streamlog::out() << "SystemCommand::execute: Problem reopening " << errorFileName ;
      return false;
    }
  }
  
  // readback all the file lines as we did already for the stdout pipe
  while ( fgets( buf, sizeof( buf ), fpErrorFile ) ) {
    // convert char to string for better and easier handling
    currentString = buf;

    // check if the current string is ending with a new line. If not,
    // this means that we are actually not getting enough character
    // for this line and the user should increase MAX_LINE_SIZE
    if ( currentString [currentString.size() - 1] != '\n' ) {
      if( streamlog::out.write<streamlog::ERROR>() ){
	streamlog::out() << "SystemCommand::execute: MAX_SILE_SIZE is too small, try recompile the code using MAX_LINE_SIZE = " 
			 << 2 * MAX_LINE_SIZE;	
	return false;
      }
    
    }
    
    // ok add the currentString to the error vector removing the last
    // character that is an end of line. This is making the future
    // manipulation of the string much easier
    _stdErrorVec.push_back( currentString.substr(0, currentString.size() - 1) );
  }

  
  if ( fclose( fpErrorFile ) == EOF ) { 
    if( streamlog::out.write<streamlog::ERROR>() ){
      streamlog::out() << "SystemCommand::execute: Problem closing " << errorFileName ;
      return false;
    }
  }

  if ( unlink( errorFileName.c_str() ) == -1 ) {
    if( streamlog::out.write<streamlog::ERROR>() ){
      streamlog::out() << "SystemCommand::execute: Unable to remove " << errorFileName ;
      return false;
    }
  }
  
  return true;
}

inline size_t SystemCommand::getStdOutputSize() const {
  return _stdOutputVec.size();
}

inline size_t SystemCommand::getStdErrorSize() const {
  return _stdErrorVec.size();
}

inline string SystemCommand::getStdOutputAt(size_t iPos) const {
  return _stdOutputVec.at( iPos );
}

inline string SystemCommand::getStdErrorAt(size_t iPos) const {
  return _stdErrorVec.at( iPos );
}

