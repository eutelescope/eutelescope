########################################################
# cmake file for building Marlin example Package
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################


# project name
PROJECT( Eutelescope )

# project version; before changing this, read README.release_procedure
SET( ${PROJECT_NAME}_VERSION_MAJOR 1 )
SET( ${PROJECT_NAME}_VERSION_MINOR 0 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

# These switches activate/deactivate compiler switches such as the
# inclusion of debugging symbols and code optimizations. For actual
# data analysis we recommend setting CMAKE_BUILD_TYPE to "Release"
# while for developement "RelWithDebInfo" is best used. For actual
# debugging e.g. using gdb, please consider "Debug" which will turn
# off all (!) optimizations but will be "gdb-friendlier"
# !RELEASE: comment out the following line and comment in either one of the release ones
# SET(CMAKE_BUILD_TYPE Debug) # use for debugging with e.g. gdb: add. warnings, no optimizations, debug symbols
SET(CMAKE_BUILD_TYPE RelWithDebInfo) # useful for normal development, runs faster but compilation takes longer: gcc optimizations (-O2) with debug symbols and add. warnings
#SET(CMAKE_BUILD_TYPE Release) # for production: full gcc optimization (-O3), no debug symbols, no debug messages


### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED ILCTEST)

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )

# additional packages to be searched for by cmake
LIST( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake )


FIND_PACKAGE( Marlin 1.0 REQUIRED )
FIND_PACKAGE( MarlinUtil 1.2 REQUIRED )
FIND_PACKAGE( GSL )
FIND_PACKAGE( AIDA )
FIND_PACKAGE( ROOT COMPONENTS Minuit Geom )
FIND_PACKAGE( LCCD  REQUIRED )               

# search for Eigen (linear algebra) library
FIND_PACKAGE( Eigen2 REQUIRED)
INCLUDE_DIRECTORIES( ${EIGEN2_INCLUDE_DIR} )

# development mode:

# the Geant4 be compiled with SoXt and Coin3D and Xerces-C libraries
# Geant4 environment loaded as 
# $ILCSOFT/geant4/build-9.5.p02/geant4make.sh
# now the allpix lib can be loaded as well 
# allpix lib if checkout and compiled here:
# export ALLPIX="/home/ilcsoft/v01-16-02/Eutelescope/HEAD/external/allpix"
# the library is placed by default to
# GBL if checkout and compiled as shared lib:
# export GBL="/home/ilcsoft/v01-16-02/Eutelescope/HEAD/external/gbl/cpp"
# the libs is here
#
 FIND_PACKAGE( XERCESC )               
 FIND_PACKAGE( G4 )               
 FIND_PACKAGE( ALLPIX )               
# 
# General Broken Line track fitter as shared lib:
 FIND_PACKAGE( GBL )

FOREACH( pkg Marlin MarlinUtil GSL AIDA ROOT LCCD GBL XERCESC G4 ALLPIX )
    IF( ${pkg}_FOUND )
        # include as "system" libraries: gcc will be less verbose w.r.t. warnings
        INCLUDE_DIRECTORIES(SYSTEM ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
    ELSE()
        MESSAGE( STATUS "${pkg} not found" )
    ENDIF()
ENDFOREACH()

#MESSAGE (STATUS "${XERCESC_LIBRARIES}" )
#MESSAGE (STATUS "${XERCESC_INCLUDE_DIRS}" )


### LIBRARY #################################################################

IF( GBL_FOUND )
    ADD_DEFINITIONS( "-DUSE_GBL " )
    GET_FILENAME_COMPONENT( GBL_LIBRARY_FULL_PATH ${GBL_GBL_LIBRARY} REALPATH )
    INSTALL( FILES ${GBL_LIBRARY_FULL_PATH} DESTINATION lib ) # PROGRAMS is like FILES but sets executable permissions
ELSE()
    MESSAGE( STATUS "WARNING: failed to configure Eutelescope with GBL!!" )
ENDIF()

IF( ALLPIX_FOUND )
    ADD_DEFINITIONS( "-DUSE_ALLPIX " )
    GET_FILENAME_COMPONENT( ALLPIX_LIBRARY_FULL_PATH ${ALLPIX_ALLPIX_LIBRARY} REALPATH )
    INSTALL( FILES ${ALLPIX_LIBRARY_FULL_PATH} DESTINATION lib ) # PROGRAMS is like FILES but sets executable permissions

    # two ROOT libraries which are not always included in "root-config --glibs"
    find_library(
      LXMLPARSER
      NAMES XMLParser 
      PATHS  ${ROOT_LIBRARY_FULL_PATH}
      PATH_SUFFIXES lib64 lib32 lib
      NO_DEFAULT_PATH
    )
  

    find_library(
      LGEOM
      NAMES Geom 
      PATHS  ${ROOT_LIBRARY_FULL_PATH}
      PATH_SUFFIXES lib64 lib32 lib
      NO_DEFAULT_PATH
    )

    LINK_LIBRARIES( ${LGEOM} )
    LINK_LIBRARIES( ${LXMLPARSER} )
 
    MESSAGE ( STATUS "linking TGeom : ${LGEOM}" )
    MESSAGE ( STATUS "linking TXMLParser : ${LXMLPARSER}" )

ELSE()
    IF( DEFINED "$ENV{ALLPIX}" )
      MESSAGE( STATUS "WARNING: failed to configure Eutelescope with ALLPIX!!" )
    ENDIF()
ENDIF()


IF ( AIDA_FOUND)
  ADD_DEFINITIONS("-DUSE_AIDA")
ENDIF()

IF ( GSL_FOUND )
  ADD_DEFINITIONS("-DUSE_GSL")
ENDIF()

# these are needed anyway...
ADD_DEFINITIONS("-DUSE_MARLIN")
ADD_DEFINITIONS("-DUSE_GEAR")
ADD_DEFINITIONS("-DUSE_MARLINUTIL")
ADD_DEFINITIONS("-DUSE_CLHEP")
ADD_DEFINITIONS("-DUSE_CED")


# !RELEASE: switches on output parsed for automatic nightly tests
# (only useful if using CDash web frontend to CTest otherwise produces superfluous output)
ADD_DEFINITIONS("-DDO_TESTING")


# ---------------------------------------------------------------------------


# include directories
INCLUDE_DIRECTORIES( ./include )

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wall -ansi -pedantic -Wextra" )
ADD_DEFINITIONS( "-Wno-long-long" )

# debuging warnings and switches - only active if  CMAKE_BUILD_TYPE  is set to debug
SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fno-inline -fdiagnostics-show-option -Weffc++ -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wsign-promo -Wstrict-null-sentinel -Wswitch-default -Wundef"  CACHE STRING "Debug options." FORCE )
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -fdiagnostics-show-option -Weffc++ -Wcast-align -Wcast-qual -Wdisabled-optimization -Winit-self -Wmissing-include-dirs -Wredundant-decls -Wsign-promo -Wstrict-null-sentinel -Wswitch-default -Wundef"  CACHE STRING "Debug options." FORCE )
# also useful: -Wshadow (however, GCC 4.1 uses this even for system libraries, causing many warnings from Marlin & co)

# add library
SET( libname ${PROJECT_NAME} )
AUX_SOURCE_DIRECTORY( ./src library_sources )
ADD_SHARED_LIBRARY( ${libname} ${library_sources} )
INSTALL_SHARED_LIBRARY( ${libname} DESTINATION lib )
#INSTALL_SHARED_LIBRARY( GBL DESTINATION lib )

# Properly add the CMSPixelDecoder external dependency:
ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/external/CMSPixelDecoder)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/external/CMSPixelDecoder)
# ..and link it to libEUTelescope:
TARGET_LINK_LIBRARIES( ${libname} CMSPixelDecoder )


# used for alignment if Eutelescope was build with ROOT support
IF( ROOT_FOUND AND ROOT_MINUIT_FOUND )
    TARGET_LINK_LIBRARIES( ${libname} ${ROOT_MINUIT_LIBRARY} )
ENDIF()

IF( ROOT_FOUND AND ROOT_GEOM_FOUND )
    TARGET_LINK_LIBRARIES( ${libname} ${ROOT_GEOM_LIBRARY} )
ENDIF()

MACRO( ADD_EUTELESCOPE_TOOL _name )
    ADD_EXECUTABLE( ${_name} src/exec/${_name}.cxx )
    TARGET_LINK_LIBRARIES( ${_name} ${libname} )
    INSTALL( TARGETS ${_name} DESTINATION bin )
ENDMACRO()

ADD_EUTELESCOPE_TOOL( pede2lcio )
ADD_EUTELESCOPE_TOOL( pedestalmerge )


# !RELEASE: REMOVE FOR RELEASE VERSIONS
# electric fence
 FIND_LIBRARY(EFENCE_LIBRARY
     NAMES efence
     DOC "The electric fence library"
     PATHS
         /usr/lib
         /usr/local/lib
         /usr/efence/lib
         /usr/share/efence/lib
     )
 IF(EFENCE_LIBRARY)
   MESSAGE("Electric Fence library was found and will be linked against ${PROJECT_NAME}")
   target_link_libraries(${libname} ${EFENCE_LIBRARY})
 ENDIF()

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

# create some useful symbolic links in bin directory at install time
install (CODE
    "EXECUTE_PROCESS(COMMAND ln -sf ${PROJECT_SOURCE_DIR}/jobsub/jobsub.py ${PROJECT_SOURCE_DIR}/bin/jobsub)
    EXECUTE_PROCESS(COMMAND ln -sf ${PROJECT_SOURCE_DIR}/jobsub/jobsub.py ${PROJECT_SOURCE_DIR}/bin/jobsub.py)
    EXECUTE_PROCESS(COMMAND ln -sf ${PROJECT_SOURCE_DIR}/tools/pyroplot/pyroplot.py ${PROJECT_SOURCE_DIR}/bin/pyroplot.py)
    EXECUTE_PROCESS(COMMAND ln -sf ${PROJECT_SOURCE_DIR}/tools/pyroplot/pyroplot.py ${PROJECT_SOURCE_DIR}/bin/pyroplot)
    EXECUTE_PROCESS(COMMAND ln -sf ${PROJECT_SOURCE_DIR}/tools/parsepede/parsemilleout.sh ${PROJECT_SOURCE_DIR}/bin/parsemilleout.sh)"
)

#  _            _       
# | |_ ___  ___| |_ ___ 
# | __/ _ \/ __| __/ __|
# | ||  __/\__ \ |_\__ \
#  \__\___||___/\__|___/
#=============================================================================                      
# EUTelescope uses the CMake/CTest framework for daily automated,
#data-driven tests to constantly verify the validity of the output. The
#basis for the tests are the examples provided with the software
#together with a real data sample and the expected output in form of
#root files.
#
#If you are a user of EUTelescope and you have already taken data and
#performed you analysis, we would be very interested in obtaining a
#data sample with the corresponding analysis results including the
#EUTelescope configuration used so that we can include your analysis in
#the regular testing. This guarantees that your analysis will be well
#supported in future releases of EUTelescope and it might serve as a
#good starting point for others.
#
#The test can also be used by any developer who wishes to verify
#changes made to the EUTelescope code. For a short how-to see
#README.tests; for more detailed information see the CMake manual or
#web page (http://www.cmake.org/Wiki/CMake/Testing_With_CTest).

ENABLE_TESTING()
INCLUDE(CTest)


# ----------------------------------------------------------------------------
# TEST IF ENVIRONMENT IS OK FOR RUNNING TESTS
# ----------------------------------------------------------------------------

IF( NOT DEFINED ENV{EUTELESCOPE} OR NOT DEFINED ENV{MARLIN})
      MESSAGE( "trying to load build and run environment for EUTelescope" )

  IF( DEFINED ENV{ILCSOFT})
        EXECUTE_PROCESS( COMMAND sh -c " [ -f $ILCSOFT/Eutelscope/trunk/build_env.sh ] && . $ILCSOFT/Eutelescope/trunk/build_env.sh" )

	IF( NOT DEFINED ENV{EUTELESCOPE} OR NOT DEFINED ENV{MARLIN} )
        # try to load build_env.sh environment - if running 'make' from 'build' directory, it should be located in '../'
        EXECUTE_PROCESS( COMMAND sh -c " [ -f '../build_env.sh' ] && . ../build_env.sh" )
	IF( DEFINED ENV{EUTELESCOPE} AND DEFINED ENV{MARLIN} )
	    MESSAGE( "EUTELESCOPE environment automatically loaded from ../build_env.sh" )
	ENDIF()
     ENDIF()
  ENDIF()

ENDIF()

# ----------------------------------------------------------------------------
# ACTUAL TESTS BEGIN HERE
# ----------------------------------------------------------------------------

# check environment
IF( DEFINED ENV{EUTELESCOPE} AND DEFINED ENV{MARLIN})

# use valgrind for dynamic analysis
set (CTEST_MEMORYCHECK_COMMAND "/usr/bin/valgrind")
set (CTEST_MEMORYCHECK_COMMAND_OPTIONS "-q --num-callers=50 --tool=memcheck --track-origins=yes --trace-children=yes --leak-check=full --leak-resolution=low --suppressions=$ENV{ROOTSYS}/etc/valgrind-root.supp --suppressions=$ENV{EUTELESCOPE}/test/valgrind/valgrind-python.supp")

# User Options to be included in DartConfiguration.tcl
set(MEMORYCHECK_COMMAND ${CTEST_MEMORYCHECK_COMMAND} CACHE STRING "" FORCE)
set(MEMORYCHECK_COMMAND_OPTIONS ${CTEST_MEMORYCHECK_COMMAND_OPTIONS} CACHE STRING "" FORCE)


# cppcheck should be a useful check for developers to test their code against
# if cppcheck is installed, offer static code analysis test
find_program(CPPCHECK cppcheck /scratch/local/bin) # hint at path /scratch/... used on EUTelTesting server; will also parse $PATH
if(CPPCHECK)
  MESSAGE("cppcheck was found and will be made available as test.")
  ADD_TEST( TestCppcheckSources sh -c "${CPPCHECK} --enable=all -I $ENV{EUTELESCOPE}/include -q -rp=$ENV{EUTELESCOPE} --std=c99 --std=posix $ENV{EUTELESCOPE}/src" )
  SET_TESTS_PROPERTIES (TestCppcheckSources PROPERTIES
    FAIL_REGULAR_EXPRESSION "(error)"
    )
#else()
#  MESSAGE("cppcheck was not found - omitting cppcheck static code analysis test.")
endif()

# Developers: please consider using these tests to verify your code!
# to obtain the necessary data files, please check the corresponding
# README files in the example folders and/or contact the EUTelescope
# software coordinators.
# !RELEASE: the data-driven tests require DESY-AFS access and are not
# useful for typical users (comment out the following lines)
  INCLUDE(jobsub/examples/datura-noDUT/testing.cmake)
  INCLUDE(jobsub/examples/datura-alone/testing.cmake)
  INCLUDE(jobsub/examples/anemone-2FEI4/testing.cmake)

  configure_file (
	       ./CTestCustom.cmake.in   # set custom settings (e.g. tests to be skipped for valgrind)
	       ${PROJECT_BINARY_DIR}/CTestCustom.cmake
	       )
# /!RELEASE

ELSE()
    MESSAGE( "EUTELESCOPE environment variable not defined, cannot run tests - please run 'source build_env.sh' in the EUTel root directory first." )
ENDIF()

