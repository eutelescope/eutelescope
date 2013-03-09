#########################################################
# cmake module for finding ALLPIX
#
#
# returns:
#   ALLPIX_FOUND        : set to TRUE or FALSE
#   ALLPIX_INCLUDE_DIRS : paths to  includes
#   ALLPIX_LIBRARY_DIRS : paths to  libraries
#   ALLPIX_LIBRARIES    : list of  libraries
#
# @author Igor Rubinskiy, DESY
#########################################################

#MESSAGE ( STATUS "${ALLPIX} ALLPIX ")
#MESSAGE ( STATUS "${ALLPIX_DIR} ALLPIX_DIR ")
#MESSAGE ( STATUS "${ALLPIX_HOME} ALLPIX_HOME ")


#SET( ALLPIXLIB "$ENV{G4WORKDIR}/tmp/Linux-g++/allpix/")
#MESSAGE ( STATUS "${ALLPIXLIB}" )
 
# -- fix for backwards compatibility
IF( NOT ALLPIX_DIR AND ALLPIX_HOME )
    SET( ALLPIX_DIR "${ALLPIX_HOME}" )
ENDIF()

IF( NOT ALLPIX_DIR AND DEFINED ENV{ALLPIX_HOME} )
    SET( ALLPIX_DIR "$ENV{ALLPIX_HOME}" )
ENDIF()

IF( NOT ALLPIX_DIR AND DEFINED ENV{ALLPIX} )
    SET( ALLPIX_DIR "$ENV{ALLPIX}" )
ENDIF()

#MESSAGE ( STATUS "ALLPIX: ${ALLPIXLIB} ")



# ---------- includes ---------------------------------------------------------
SET( ALLPIX_INCLUDE_DIRS ALLPIX_INCLUDE_DIRS-NOTFOUND )
MARK_AS_ADVANCED( ALLPIX_INCLUDE_DIRS )


FIND_PATH( ALLPIX_INCLUDE_DIRS NAMES AllPixRun.hh *hh
    PATHS ${ALLPIX_DIR}/include
    NO_DEFAULT_PATH
)


# ---------- libraries --------------------------------------------------------
INCLUDE( MacroCheckPackageLibs )

LIST( APPEND ALLPIX_LIB_SEARCH_PATH $ENV{G4WORKDIR}/tmp/Linux-g++/allpix/  )

# only standard libraries should be passed as arguments to CHECK_PACKAGE_LIBS
# additional components are set by cmake in variable PKG_FIND_COMPONENTS
# first argument should be the package name
CHECK_PACKAGE_LIBS( ALLPIX allpix  )


# ---------- final checking ---------------------------------------------------
INCLUDE( FindPackageHandleStandardArgs )
# set ALLPIX_FOUND to TRUE if all listed variables are TRUE and not empty
FIND_PACKAGE_HANDLE_STANDARD_ARGS( ALLPIX DEFAULT_MSG ALLPIX_DIR ALLPIX_INCLUDE_DIRS ALLPIX_LIBRARIES )

#MESSAGE (STATUS "${ALLPIX}")
#MESSAGE (STATUS "${ALLPIX_FOUND}")
#MESSAGE (STATUS "${ALLPIX_DIR}")
#MESSAGE (STATUS "${ALLPIX_INCLUDE_DIRS}")
#MESSAGE (STATUS "${ALLPIX_LIBRARIES}")


