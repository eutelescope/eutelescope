########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################


# project name
PROJECT( Eutelescope )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 5 )
SET( ${PROJECT_NAME}_VERSION_PATCH 1 )



### DEPENDENCIES ############################################################

# for backwards compatibility
LIST( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

# FindILCSOFT_CMAKE_MODULES.cmake needs to be in CMAKE_MODULE_PATH
FIND_PACKAGE( ILCSOFT_CMAKE_MODULES 1.11 REQUIRED )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )

# ------- required packages -------------------
SET( required_packages LCIO GEAR Marlin MarlinUtil)

# ------- optional packages -------------------
SET( optional_packages EUDAQ GSL CLHEP CED AIDA ROOT )

# need to use ROOT's Minuit library
SET( ROOT_USE_COMPONENTS Minuit )

# ------- do the search -----------------------

# fix for backwards compatibility
FOREACH( pkg ${required_packages} ${optional_packages} )
    IF( NOT ${pkg}_DIR AND DEFINED ${pkg}_HOME )
        SET( ${pkg}_DIR ${${pkg}_HOME} )
    ENDIF()
ENDFOREACH()

IF( APPLE )
 SET( required_packages ${required_packages} streamlog )
 IF( NOT streamlog_DIR )
    SET( streamlog_DIR "${Marlin_DIR}" CACHE PATH "Path to streamlog" FORCE)
 ENDIF()
ENDIF()

FOREACH( pkg ${required_packages} )
    FIND_PACKAGE( ${pkg} REQUIRED )
ENDFOREACH()

FOREACH( pkg ${optional_packages} )
    FIND_PACKAGE( ${pkg} QUIET )
ENDFOREACH()

FOREACH( pkg ${required_packages} ${optional_packages} )
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
    ELSE()
        MESSAGE( STATUS "${pkg} not found" )
    ENDIF()
ENDFOREACH()



### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" ON )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( doc )
ELSE()
    MESSAGE( STATUS "INSTALL_DOC set to OFF" )
    SET( INSTALL_DOC OFF )
ENDIF()



### LIBRARY #################################################################

# ------ add preprocessor flags ---------------------------------------------
IF( EUDAQ_FOUND )
    ADD_DEFINITIONS( "-DUSE_EUDAQ -DEUDAQ_PLATFORM=PF_LINUX -DEUDAQ_FUNC=__PRETTY_FUNCTION__" )
ELSE()
    MESSAGE( STATUS "WARNING: failed to configure Eutelescope with EUDAQ!!" )
ENDIF()

IF ( CED_FOUND )
  ADD_DEFINITIONS("-DUSE_CED")
ENDIF()

IF ( Marlin_FOUND)
  ADD_DEFINITIONS("-DUSE_MARLIN")
ENDIF()

IF ( AIDA_FOUND)
  ADD_DEFINITIONS("-DUSE_AIDA")
ENDIF()

IF ( MarlinUtil_FOUND )
  ADD_DEFINITIONS("-DUSE_MARLINUTIL")
ENDIF()

IF ( GSL_FOUND )
  ADD_DEFINITIONS("-DUSE_GSL")
ENDIF()

IF ( GEAR_FOUND )
  ADD_DEFINITIONS("-DUSE_GEAR")
ENDIF()
# ---------------------------------------------------------------------------


# include directories
INCLUDE_DIRECTORIES( ./include )
INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
ADD_DEFINITIONS( "-Wno-long-long" )

# add library
SET( libname ${PROJECT_NAME} )
AUX_SOURCE_DIRECTORY( ./src library_sources )
ADD_SHARED_LIBRARY( ${libname} ${library_sources} )
INSTALL( TARGETS ${libname} DESTINATION lib )

# used for alignment if Eutelescope was build with ROOT support
IF( ROOT_FOUND AND ROOT_MINUIT_FOUND )
    TARGET_LINK_LIBRARIES( ${libname} ${ROOT_MINUIT_LIBRARY} )
ENDIF()

MACRO( ADD_EUTELESCOPE_TOOL _name )
    ADD_EXECUTABLE( ${_name} src/exec/${_name}.cxx )
    TARGET_LINK_LIBRARIES( ${_name} ${libname} )
    INSTALL( TARGETS ${_name} DESTINATION bin )
ENDMACRO()

ADD_EUTELESCOPE_TOOL( pede2lcio )
ADD_EUTELESCOPE_TOOL( pedestalmerge )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

