#############################################################################
# Example of a BuildSetup.cmake file for a Marlin Package
#
# For building with cmake type:
# (1) $ mkdir build
# (2) $ cd build
# (3) $ cmake -C ../BuildSetup.cmake ..
# (4) $ make install
#
# @author Jan Engels, DESY
#############################################################################


#############################################################################
# Setup path variables
#############################################################################

# ILC_HOME
SET( ILC_HOME "/home/toto/ilc"
  CACHE PATH "Path to ILC Software" FORCE )

# Path to Marlin
SET( Marlin_HOME "${ILC_HOME}/Marlin/HEAD"
  CACHE PATH "Path to Marlin" FORCE )

# Path to MarlinUtil
SET( MarlinUtil_HOME "${ILC_HOME}/MarlinUtil/HEAD"
  CACHE PATH "Path to MarlinUtil" FORCE )

# Path to LCIO
SET( LCIO_HOME "${ILC_HOME}/lcio/HEAD"
  CACHE PATH "Path to LCIO" FORCE )

## Path to GEAR
SET( GEAR_HOME "${ILC_HOME}/gear/HEAD"
  CACHE PATH "Path to GEAR" FORCE )

# Path to CLHEP
SET( CLHEP_HOME "${ILC_HOME}/CLHEP/2.0.2.2"
    CACHE PATH "Path to CLHEP" FORCE )

## Path to RAIDA
SET( RAIDA_HOME "${ILC_HOME}/RAIDA/v01-02"
  CACHE PATH "Path to RAIDA" FORCE )

# CMake Modules Path
SET( CMAKE_MODULE_PATH "${ILC_HOME}/CMakeModules/HEAD"
    CACHE PATH "Path to CMake Modules" FORCE )

# Path to EUDAQ (define this to build Eutelescope with EUDAQ)
#SET( EUDAQ_HOME "path_to_eudaq"
#  CACHE PATH "Path to EUDAQ" FORCE )

#############################################################################
# Optional packages
#############################################################################
                                                                                                                                                            
# if you want to build and link your package with one or more additional
# tools you also have to define the corresponding "home" paths above
SET( BUILD_WITH "GEAR RAIDA" CACHE STRING "Build Eutelescope with these optional packages" FORCE )
                                                                                                                                                            
#############################################################################
# Project options
#############################################################################

#SET( INSTALL_DOC OFF CACHE BOOL "Set to OFF to skip build/install Documentation" FORCE )

# set cmake build type, default value is: RelWithDebInfo
# possible options are: None Debug Release RelWithDebInfo MinSizeRel
#SET( CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build" FORCE )

#############################################################################
# Advanced options
#############################################################################

#SET( BUILD_SHARED_LIBS OFF CACHE BOOL "Set to OFF to build static libraries" FORCE )

# where do you want files (libraries, binaries, ...) to be placed for a "make install"
# if you uncomment the next line don't forget to change 'mymarlin' to your package name
#SET( CMAKE_INSTALL_PREFIX "/foo/bar" CACHE STRING "Where to install mymarlin" FORCE )
