#########################################################
# cmake module for finding GBL
#
#
# returns:
#   GBL_FOUND        : set to TRUE or FALSE
#   GBL_INCLUDE_DIRS : paths to gsl includes
#   GBL_LIBRARY_DIRS : paths to gsl libraries
#   GBL_LIBRARIES    : list of gsl libraries
#
# @author Jan Engels, DESY
#########################################################


# -- fix for backwards compatibility
IF( NOT GBL_DIR AND GBL_HOME )
    SET( GBL_DIR "${GBL_HOME}" )
ENDIF()

IF( NOT GBL_DIR AND DEFINED ENV{GBL_HOME} )
    SET( GBL_DIR "$ENV{GBL_HOME}" )
ENDIF()

IF( NOT GBL_DIR AND DEFINED ENV{GBL} )
    SET( GBL_DIR "$ENV{GBL}" )
ENDIF()

IF( NOT GBL_DIR )
    SET( GBL_DIR "$ENV{EUTELESCOPE}/gbl/cpp/" )
ENDIF()

#MESSAGE (STATUS "gbl_home: ${GBL_HOME}")
#MESSAGE (STATUS "env gbl_home: $ENV{GBL_HOME}")
#MESSAGE (STATUS "env gbl: $ENV{GBL}")
#MESSAGE (STATUS "gbl_dir: ${GBL_DIR}")



# ---------- includes ---------------------------------------------------------
SET( GBL_INCLUDE_DIRS GBL_INCLUDE_DIRS-NOTFOUND )
MARK_AS_ADVANCED( GBL_INCLUDE_DIRS )

FIND_PATH( GBL_INCLUDE_DIRS NAMES  include/GblTrajectory.h
    PATHS ${GBL_DIR} 
    NO_DEFAULT_PATH
)

#INCLUDE_DIRECTORIES( ${GBL_INCLUDE_DIRS} )

MESSAGE (STATUS "${GBL_INCLUDE_DIRS}") 

EXECUTE_PROCESS( COMMAND  " cd ${GBL}  &&  make install  && cd ${EUTELESCOPE}/build" )


# ---------- libraries --------------------------------------------------------
INCLUDE( MacroCheckPackageLibs )

# GBL library lives in bin directory rather than in lib
LIST( APPEND GBL_LIB_SEARCH_PATH ${GBL_DIR}/lib  ${GBL_DIR} )

# only standard libraries should be passed as arguments to CHECK_PACKAGE_LIBS
# additional components are set by cmake in variable PKG_FIND_COMPONENTS
# first argument should be the package name
CHECK_PACKAGE_LIBS( GBL GBL )



# ---------- final checking ---------------------------------------------------
INCLUDE( FindPackageHandleStandardArgs )
# set GBL_FOUND to TRUE if all listed variables are TRUE and not empty
FIND_PACKAGE_HANDLE_STANDARD_ARGS( GBL DEFAULT_MSG GBL_DIR GBL_INCLUDE_DIRS GBL_LIBRARIES )

