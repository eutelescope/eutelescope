// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELUTIL
#define EUTELUTIL 1

// eutelescope includes ".h" 



// system includes
#include <vector>
#include <string>
#include <stdexcept>

namespace eutelescope {

  //! Utility to call a command from within a processor
  /*! This utility class has been implemented in order to easily
   *  launch external command within a standard C++ code. The
   *  advantage with respect to the use of system() C call is that the
   *  user can get both the standard error and output from the
   *  command. 
   *
   *  @Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @Version $Id: EUTelUtil.h,v 1.1 2008-07-04 10:02:15 bulgheroni Exp $
   *
   */
  class SystemCommand {

  public:

    //! Default constructor
    /*! The user as to provide the command he wants to execute as a
     *  string with all the arguments.
     *  
     *  @param command String with the command and all the arguments
     */ 
    SystemCommand(std::string command);

    //! Execute the command
    /*! The command defined in the command string is actually executed
     *  only when this method is called.
     *
     *  @return True is the command execution was
     *  successful. N.B. this is referring to the command
     *  execution and not to the command output.
     */ 
    bool execute();
    
    //! Get the number of lines of stdout
    /*! During the command execution the standard output is redirected
     *  to a pipe and each line is stored into a vector of strings.
     *
     *  @return the number of lines in the std output pipe
     */ 
    size_t getStdOutputSize() const;

    //! Get the number of lines of stderr
    /*! During the command execution the standard error is redirected
     *  to a pipe and each line is stored into a vector of strings.
     *
     *  @return the number of lines in the std error pipe
     */ 
    size_t getStdErrorSize()  const;
    
    //! Get a line of the std output 
    /*! This method returns the line of the standard output at
     *  position iPos.
     *
     *  @param iPos The line number of the stdoutput
     *  @return A string with the corresponding line 
     *  @thrown std::out_of_range& if iPos is greater than the vector size
     */ 
    std::string  getStdOutputAt(size_t iPos) const;

    //! Get a line of the std error
    /*! This method returns the line of the standard error at
     *  position iPos.
     *
     *  @param iPos The line number of the stderror
     *  @return A string with the corresponding line 
     *  @thrown std::out_of_range& if iPos is greater than the vector size
     */ 
    std::string  getStdErrorAt(size_t iPos)  const;

  private:
    
    //! The command to be executed 
    std::string _command;
    
    //! The vector with all the string of the stdout
    std::vector< std::string > _stdOutputVec;

    //! The vector with all the string of the stderr
    std::vector< std::string > _stdErrorVec;
    
  };

}

#endif
