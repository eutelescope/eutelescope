#########################################################
# cmake module for finding EUDAQ
#
#
# returns:
#   EUDAQ_FOUND        : set to TRUE or FALSE
#   EUDAQ_INCLUDE_DIRS : paths to gsl includes
#   EUDAQ_LIBRARY_DIRS : paths to gsl libraries
#   EUDAQ_LIBRARIES    : list of gsl libraries
#
# @author Jan Engels, DESY
#########################################################

# -- fix for backwards compatibility
IF( NOT EUDAQ_DIR AND EUDAQ_HOME )
    SET( EUDAQ_DIR "${EUDAQ_HOME}" )
ENDIF()

IF( NOT EUDAQ_DIR AND DEFINED ENV{EUDAQ_HOME} )
    SET( EUDAQ_DIR "$ENV{EUDAQ_HOME}" )
ENDIF()

IF( NOT EUDAQ_DIR AND DEFINED ENV{EUDAQ} )
    SET( EUDAQ_DIR "$ENV{EUDAQ}" )
ENDIF()



# ---------- includes ---------------------------------------------------------
SET( EUDAQ_INCLUDE_DIRS EUDAQ_INCLUDE_DIRS-NOTFOUND )
MARK_AS_ADVANCED( EUDAQ_INCLUDE_DIRS )

FIND_PATH( EUDAQ_INCLUDE_DIRS NAMES eudaq/RunControl.hh
    PATHS ${EUDAQ_DIR}/main/include
    NO_DEFAULT_PATH
)

IF( NOT EUDAQ_DIR )
    FIND_PATH( EUDAQ_INCLUDE_DIRS NAMES eudaq/RunControl.hh )
ENDIF()



# ---------- libraries --------------------------------------------------------
INCLUDE( MacroCheckPackageLibs )

# EUDAQ library lives in bin directory rather than in lib
LIST( APPEND EUDAQ_LIB_SEARCH_PATH ${EUDAQ_DIR}/bin )

# only standard libraries should be passed as arguments to CHECK_PACKAGE_LIBS
# additional components are set by cmake in variable PKG_FIND_COMPONENTS
# first argument should be the package name
CHECK_PACKAGE_LIBS( EUDAQ eudaq )



# ---------- final checking ---------------------------------------------------
INCLUDE( FindPackageHandleStandardArgs )
# set EUDAQ_FOUND to TRUE if all listed variables are TRUE and not empty
FIND_PACKAGE_HANDLE_STANDARD_ARGS( EUDAQ DEFAULT_MSG EUDAQ_DIR EUDAQ_INCLUDE_DIRS EUDAQ_LIBRARIES )

