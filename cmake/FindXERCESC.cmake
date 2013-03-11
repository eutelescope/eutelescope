# - Find XERCESC
# This module tries to find the XERCESC library and headers.
# Once done this will define
#
#   XERCESC_FOUND - system has XERCESC headers and libraries
#   XERCESC_INCLUDE_DIRS - the include directories needed for XERCESC
#   XERCESC_LIBRARIES - the libraries needed to use XERCESC
#
# Variables used by this module, which can change the default behaviour and
# need to be set before calling find_package:
#
#   XERCESC_ROOT_DIR            Root directory to XERCESC installation. Will
#                               be used ahead of CMake default path.
#
# The following advanced variables may be used if the module has difficulty
# locating XERCESC or you need fine control over what is used.
#
#   XERCESC_INCLUDE_DIR
#
#   XERCESC_LIBRARY
#
# Copyright (c) 2009, Ben Morgan, <Ben.Morgan@warwick.ac.uk>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.


IF ( DEFINED "$ENV{XERCESC_ROOT_DIR}" )

# Look for the header - preferentially searching below XERCESC_ROOT_DIR
find_path(
    XERCESC_INCLUDE_DIR 
    NAMES xercesc/util/XercesVersion.hpp
    PATHS ${XERCESC_ROOT_DIR}
    PATH_SUFFIXES include
    NO_DEFAULT_PATH 
)

# If we didn't find it there, fall back to some standard search paths
find_path(
    XERCESC_INCLUDE_DIR
    NAMES xercesc/util/XercesVersion.hpp
)

MESSAGE ( STATUS "xerces: ${XERCESC_FOUND}" )
MESSAGE ( STATUS "xerces: ${XERCESC_ROOT_DIR}" )
MESSAGE ( STATUS "xerces: ${XERCESC_LIBRARY}" )
MESSAGE ( STATUS "xerces: ${XERCESC_INCLUDE_DIR}" )

EXECUTE_PROCESS(
    COMMAND find ${XERCESC_ROOT_DIR}/include/xercesc/ -type d # -printf "%p;" # not recognized on mac osx
    OUTPUT_VARIABLE XERCESC_INCLUDE_DIR
    RESULT_VARIABLE _exit_code
)
 

# Look for the library, preferentially searching below XERCESC_ROOT_DIR
find_library(
    XERCESC_LIBRARY
    NAMES xerces-c 
    PATHS ${XERCESC_ROOT_DIR}
    PATH_SUFFIXES lib64 lib32 lib
    NO_DEFAULT_PATH 
)


CHECK_PACKAGE_LIBS( XERCESC xerces-c  )


include(FindPackageHandleStandardArgs)


find_package_handle_standard_args(
    XERCESC
    DEFAULT_MSG
    XERCESC_LIBRARY
    XERCESC_INCLUDE_DIR
)


if (XERCESC_FOUND)
    set(XERCESC_LIBRARIES ${XERCESC_LIBRARY})
    set(XERCESC_INCLUDE_DIRS ${XERCESC_INCLUDE_DIR})
else (XERCESC_FOUND)
    set(XERCESC_LIBRARIES)
    set(XERCESC_INCLUDE_DIRS)
endif (XERCESC_FOUND)


mark_as_advanced(
    XERCESC_LIBRARY 
    XERCESC_INCLUDE_DIR
)


ENDIF()
