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


#LIST( APPEND ALLPIX_INCLUDE_DIRS /home/ilcsoft/v01-15/Eutelescope/HEAD/external/allpix/include/ ${XERCESC_INCLUDE_DIRS} )
#SET ( ALLPIX_INCLUDE_DIRS "/home/ilcsoft/v01-15/Eutelescope/HEAD/external/allpix/include/;/usr/include/xercesc/;/usr/include/xercesc/dom;/usr/include/xercesc/dom/impl;/usr/include/xercesc/parsers;/usr/include/xercesc/sax2;/usr/include/xercesc/validators;/usr/include/xercesc/validators/schema;/usr/include/xercesc/validators/schema/identity;/usr/include/xercesc/validators/DTD;/usr/include/xercesc/validators/datatype;/usr/include/xercesc/validators/common;/usr/include/xercesc/util;/usr/include/xercesc/util/MsgLoaders;/usr/include/xercesc/util/MsgLoaders/InMemory;/usr/include/xercesc/util/regx;/usr/include/xercesc/util/MutexManagers;/usr/include/xercesc/util/FileManagers;/usr/include/xercesc/util/Transcoders;/usr/include/xercesc/util/Transcoders/IconvGNU;/usr/include/xercesc/util/NetAccessors;/usr/include/xercesc/util/NetAccessors/Socket;/usr/include/xercesc/xinclude;/usr/include/xercesc/framework;/usr/include/xercesc/framework/psvi;/usr/include/xercesc/internal;/usr/include/xercesc/sax")



# ---------- libraries --------------------------------------------------------
INCLUDE( MacroCheckPackageLibs )

LIST( APPEND ALLPIX_LIB_SEARCH_PATH /home/rubinsky/geant4_workdir/tmp/Linux-g++/allpix/  )
#                                    /home/rubinsky/geant4_workdir/tmp/Linux-g++/allpix/liballpix.so
#FIND_LIBRARY(ALLPIX_LIBRARIES
#    NAMES allpix
#    DOC "The allpix library"
#    PATHS
#	${ALLPIX_LIB_SEARCH_PATH}
#    )




# only standard libraries should be passed as arguments to CHECK_PACKAGE_LIBS
# additional components are set by cmake in variable PKG_FIND_COMPONENTS
# first argument should be the package name
CHECK_PACKAGE_LIBS( ALLPIX allpix  )

SET ( ROOTLIBS "/home/ilcsoft/v01-15/root/5.32.00/lib/libXMLParser.so;/home/ilcsoft/v01-15/root/5.32.00/lib/libGeom.so;/usr/lib/libxerces-c-3.1.so;/home/ilcsoft/v01-15/root/5.32.00/lib/libGraf3d.so;/home/ilcsoft/v01-15/root/5.32.00/lib/libGraf.so")

# ---------- final checking ---------------------------------------------------
INCLUDE( FindPackageHandleStandardArgs )
# set ALLPIX_FOUND to TRUE if all listed variables are TRUE and not empty
FIND_PACKAGE_HANDLE_STANDARD_ARGS( ALLPIX DEFAULT_MSG ALLPIX_DIR ALLPIX_INCLUDE_DIRS ALLPIX_LIBRARIES ROOTLIBS )

#MESSAGE (STATUS "${ALLPIX}")
#MESSAGE (STATUS "${ALLPIX_FOUND}")
#MESSAGE (STATUS "${ALLPIX_DIR}")
#MESSAGE (STATUS "${ALLPIX_INCLUDE_DIRS}")
#MESSAGE (STATUS "${ALLPIX_LIBRARIES}")


