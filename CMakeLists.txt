########################################################
# cmake file for building the Eutelescope package
#
# Original author: J. Engels, DESY
# Author: Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
# Version $Id: CMakeLists.txt,v 1.21 2009-04-21 20:49:25 bulgheroni Exp $
########################################################


########################################################
# CMake compatibility issues: don't modify this, please!
CMAKE_MINIMUM_REQUIRED( VERSION 2.4.6 )
#SET( CMAKE_BACKWARDS_COMPATIBILITY 2.4.6 )
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
# allow more human readable "if then else" constructs
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )
########################################################



####################################################################################################
# User section: modify as needed
####################################################################################################


# project name
PROJECT( Eutelescope )

# project version
SET( ${PROJECT_NAME}_MAJOR_VERSION 0 )
SET( ${PROJECT_NAME}_MINOR_VERSION 0 )
SET( ${PROJECT_NAME}_PATCH_LEVEL 8 )


### SETTINGS #################################################################

# project options
OPTION( BUILD_SHARED_LIBS "Set to OFF to build static libraries" ON )
OPTION( BUILD_32BIT_COMPATIBLE "Set to ON to build in 32 bit compatibility mode" ON )
OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" ON )

# project dependencies
# e.g. SET( ${PROJECT_NAME}_DEPENDS "Marlin MarlinUtil LCIO GEAR CLHEP GSL RAIDA ROOT" )
#SET( ${PROJECT_NAME}_DEPENDS "Marlin MarlinUtil LCIO GSL CLHEP" )
SET( ${PROJECT_NAME}_DEPENDS "Marlin LCIO" )

IF( APPLE )
 SET( ${PROJECT_NAME}_DEPENDS "${${PROJECT_NAME}_DEPENDS} streamlog" )
 SET( streamlog_HOME "${Marlin_HOME}" CACHE PATH "Path to streamlog" FORCE)
ENDIF()

# set default cmake build type to RelWithDebInfo
# possible options are: None Debug Release RelWithDebInfo MinSizeRel
IF( NOT CMAKE_BUILD_TYPE )
    SET( CMAKE_BUILD_TYPE "RelWithDebInfo" )
ENDIF()

# set default install prefix to project root directory
IF( CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" )
    SET( CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}" )
ENDIF()

#---------------- 32/64 bit issues ---------------------------------------
#INCLUDE( CheckTypeSize )
#CHECK_TYPE_SIZE( "void*" SIZE_OF_VOID_P )

IF( CMAKE_SIZEOF_VOID_P EQUAL 4 )
    MESSAGE( STATUS "32 bit architecture detected" )
ENDIF()

IF( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    MESSAGE( STATUS "64 bit architecture detected" )

    IF( BUILD_32BIT_COMPATIBLE )
        IF( COMMAND SET_PROPERTY )
            SET_PROPERTY(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS 0)
            SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32" )
        ELSE()
            MESSAGE( FATAL_ERROR "to build in 32 bit compatibility mode you need cmake >= 2.6" )
        ENDIF()
        MESSAGE( STATUS "Build in 32 bit compatibility mode" )
    ELSE()
        MESSAGE( STATUS "Build in native 64 bit mode" )
    ENDIF()
ENDIF()

#-------------------------------------------------------------------------


### INCLUDES #################################################################

#include directories
INCLUDE_DIRECTORIES( "${PROJECT_SOURCE_DIR}/include" )

# install include files
INSTALL( DIRECTORY "${PROJECT_SOURCE_DIR}/include"
        DESTINATION .
        PATTERN "*~" EXCLUDE
        PATTERN "*CVS*" EXCLUDE )


### SOURCES ##################################################################

# EUDAQ stuff
IF( DEFINED EUDAQ_HOME OR DEFINED ENV{EUDAQ} )
    IF( NOT DEFINED EUDAQ_HOME )
        SET( EUDAQ_HOME "$ENV{EUDAQ}" )
    ENDIF()
    # write EUDAQ_HOME to cache
    SET( EUDAQ_HOME "${EUDAQ_HOME}" CACHE PATH "Path to EUDAQ" FORCE )

    # find eudaq include dir
    FIND_PATH( EUDAQ_INCLUDE_DIR
        NAMES eudaq/RunControl.hh
        PATHS ${EUDAQ_HOME}
        PATH_SUFFIXES main/include
        NO_DEFAULT_PATH )
    IF( NOT EUDAQ_INCLUDE_DIR )
        MESSAGE( STATUS "EUDAQ include dir not found!!" )
    ENDIF()
   
    # find eudaq library
    FIND_LIBRARY( EUDAQ_LIB
        NAMES eudaq
        PATHS ${EUDAQ_HOME}
        PATH_SUFFIXES bin   
        NO_DEFAULT_PATH )
    IF( NOT EUDAQ_LIB )
        MESSAGE( STATUS "EUDAQ library not found!!" )
    ENDIF()
    
    IF( EUDAQ_INCLUDE_DIR AND EUDAQ_LIB)
        INCLUDE_DIRECTORIES( "${EUDAQ_INCLUDE_DIR}" )
        LINK_LIBRARIES( "${EUDAQ_LIB}" )
        ADD_DEFINITIONS( "-DUSE_EUDAQ -DEUDAQ_PLATFORM=PF_LINUX -DEUDAQ_FUNC=__PRETTY_FUNCTION__" )
    ELSE()
        MESSAGE( STATUS "WARNING: failed to configure ${PROJECT_NAME} with EUDAQ!!" )
    ENDIF()
ENDIF()

# require proper c++
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )

#----- need long long for int64 for now ------
ADD_DEFINITIONS( "-Wno-long-long" )

# add debug definitions
#IF( CMAKE_BUILD_TYPE STREQUAL "Debug" OR
#    CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
#    ADD_DEFINITIONS( "-DDEBUG" )
#ENDIF()

# get list of all source files
AUX_SOURCE_DIRECTORY( src library_sources )


### DOCUMENTATION ############################################################
FIND_PACKAGE( Doxygen )
IF( DOXYGEN_FOUND )

    ADD_CUSTOM_COMMAND(
        OUTPUT  "${PROJECT_SOURCE_DIR}/doc/html"
        COMMAND "${DOXYGEN_EXECUTABLE}"
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/doc"
        COMMENT "Building API Documentation..."
        VERBATIM )

    # add doc target
    ADD_CUSTOM_TARGET( doc DEPENDS
        "${PROJECT_SOURCE_DIR}/doc/html" )
ELSE()
    MESSAGE( STATUS "Doxygen not found in your system!!" )
    IF( INSTALL_DOC )
        MESSAGE( STATUS "INSTALL_DOC forced to OFF" )
        SET( INSTALL_DOC OFF )
    ENDIF()
ENDIF()

# install documentation
IF( INSTALL_DOC )
    # make sure doxygen is executed (make doc) before make install
    INSTALL( CODE "EXEC_PROGRAM(${CMAKE_BUILD_TOOL} ${PROJECT_BINARY_DIR} ARGS doc)" )
    # install documentation
    INSTALL( DIRECTORY "${PROJECT_SOURCE_DIR}/doc"
            DESTINATION .
            PATTERN "*CVS*" EXCLUDE )
ENDIF()


##########################################################################################
# End of User section: please try not to modify below this line
##########################################################################################


# library *nix style versioning
SET( ${PROJECT_NAME}_SOVERSION
    "${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION}" )
SET( ${PROJECT_NAME}_VERSION
    "${${PROJECT_NAME}_SOVERSION}.${${PROJECT_NAME}_PATCH_LEVEL}" )

# add library install path to the rpath list
SET( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
MARK_AS_ADVANCED( CMAKE_INSTALL_RPATH )

# append link pathes to rpath list
SET( CMAKE_INSTALL_RPATH_USE_LINK_PATH 1 )
MARK_AS_ADVANCED( CMAKE_INSTALL_RPATH_USE_LINK_PATH )

# output directories
SET( EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin" CACHE PATH
    "EXECUTABLE_OUTPUT_PATH" FORCE )
SET( LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib" CACHE PATH
    "LIBRARY_OUTPUT_PATH" FORCE )
MARK_AS_ADVANCED( EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH )

# DEPENDENCIES: this code has to be placed before adding any library or
# executable so that these are linked properly against the dependencies
IF( DEFINED ${PROJECT_NAME}_DEPENDS OR DEFINED BUILD_WITH OR DEFINED LINK_WITH )
    # load macro
    IF( NOT EXISTS "${CMAKE_MODULE_PATH}/MacroCheckDeps.cmake" )
        MESSAGE( FATAL_ERROR
            "\nSorry, could not find MacroCheckDeps.cmake...\n"
            "Please set CMAKE_MODULE_PATH correctly with: "
            "cmake -DCMAKE_MODULE_PATH=<path_to_cmake_modules>" )
    ENDIF()
    INCLUDE( "${CMAKE_MODULE_PATH}/MacroCheckDeps.cmake" )
    CHECK_DEPS()
ENDIF()

# add preprocessor flag: USE_CED
# this code can be deleted when CEDConfig.cmake exports the USE_CED # definition
IF ( CED_FOUND ) 
  ADD_DEFINITIONS("-DUSE_CED")
ENDIF()

# add preprocessor flag: USE_MARLINUTIL
# this code can be deleted when MarlinUtilConfig.cmake exports the USE_MARLINUTIL # definition
IF ( MarlinUtil_FOUND ) 
  ADD_DEFINITIONS("-DUSE_MARLINUTIL")
ENDIF()

# add preprocessor flag: USE_GSL
# this code can be deleted when GSLConfig.cmake exports the USE_GSL # definition
IF ( GSL_FOUND ) 
  ADD_DEFINITIONS("-DUSE_GSL")
ENDIF()


# LIBRARY
ADD_LIBRARY( lib_${PROJECT_NAME} ${library_sources} )
# create symbolic lib target for calling target lib_XXX
ADD_CUSTOM_TARGET( lib DEPENDS lib_${PROJECT_NAME} )
# change lib_target properties
SET_TARGET_PROPERTIES( lib_${PROJECT_NAME} PROPERTIES
    # create *nix style library versions + symbolic links
    VERSION ${${PROJECT_NAME}_VERSION}
    SOVERSION ${${PROJECT_NAME}_SOVERSION}
    # allow creating static and shared libs without conflicts
    CLEAN_DIRECT_OUTPUT 1
    # avoid conflicts between library and binary target names
    OUTPUT_NAME ${PROJECT_NAME} )

# install library
INSTALL( TARGETS lib_${PROJECT_NAME} DESTINATION lib PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE )

# EXECUTABLE
ADD_EXECUTABLE( pede2lcio src/pede2lcio.cc )
TARGET_LINK_LIBRARIES( pede2lcio lib_${PROJECT_NAME} )

# install executable
INSTALL( TARGETS pede2lcio DESTINATION bin PERMISSIONS 
  OWNER_READ OWNER_WRITE OWNER_EXECUTE     
  GROUP_READ GROUP_EXECUTE
  WORLD_READ WORLD_EXECUTE )


# create uninstall configuration file 
CONFIGURE_FILE( "${PROJECT_SOURCE_DIR}/cmake_uninstall.cmake.in"
                "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake"
                IMMEDIATE @ONLY )

# create uninstall target
ADD_CUSTOM_TARGET( uninstall
  "${CMAKE_COMMAND}" -P "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake" )

# create configuration file from .in file
CONFIGURE_FILE( "${PROJECT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake.in"
                "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake" @ONLY )

# install configuration file
INSTALL( FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake" DESTINATION . )

# display status message for important variables
MESSAGE( STATUS )
MESSAGE( STATUS "-------------------------------------------------------------------------------" )
MESSAGE( STATUS "BUILD_SHARED_LIBS = ${BUILD_SHARED_LIBS}" )
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
MESSAGE( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
MESSAGE( STATUS "CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}" )
IF( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    MESSAGE( STATUS "BUILD_32BIT_COMPATIBLE = ${BUILD_32BIT_COMPATIBLE}" )
ENDIF()
MESSAGE( STATUS "${PROJECT_NAME}_DEPENDS = \"${${PROJECT_NAME}_DEPENDS}\"" )
MESSAGE( STATUS "BUILD_WITH = \"${BUILD_WITH}\"" )
MESSAGE( STATUS "INSTALL_DOC = ${INSTALL_DOC}" )
MESSAGE( STATUS "Change a value with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "-------------------------------------------------------------------------------" )
MESSAGE( STATUS )

# force some variables that could be defined in the command line
# to be written to cache
SET( BUILD_SHARED_LIBS "${BUILD_SHARED_LIBS}" CACHE BOOL
    "Set to OFF to build static libraries" FORCE )
SET( CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH
    "Where to install ${PROJECT_NAME}" FORCE )
SET( CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE )
SET( CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" CACHE PATH
    "Path to custom CMake Modules" FORCE )
SET( INSTALL_DOC "${INSTALL_DOC}" CACHE BOOL
    "Set to OFF to skip build/install Documentation" FORCE )

# export build settings
INCLUDE( CMakeExportBuildSettings )
CMAKE_EXPORT_BUILD_SETTINGS( "${PROJECT_NAME}BuildSettings.cmake" )
INSTALL( FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}BuildSettings.cmake" DESTINATION lib/cmake )

# export library dependencies (keep this as the last line in the file)
EXPORT_LIBRARY_DEPENDENCIES( "${PROJECT_NAME}LibDeps.cmake" )
INSTALL( FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}LibDeps.cmake" DESTINATION lib/cmake )

