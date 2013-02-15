########################################################
# cmake file for building Marlin example Package
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################


# project name
PROJECT( Eutelescope )

# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 7 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

# IMPORTANT FOR RELEASES:
# ======================
# 0. if new major release: create new dev. branch for this release
#    - remove the lines below that are marked with "!RELEASE": these activate debugging/extra warnings and electric fence
#    - remove DO_TESTING precompiler flag
# 1. bump release version in this file
# 2. create SVN tag for this release
# 3. update web page, update installation documentation
# 4. let ILCSoft guys know (?)

# !RELEASE: comment out the following line and comment in either one of the release ones
#SET(CMAKE_BUILD_TYPE Debug) # add. warnings, no optimizations, debug symbols
SET(CMAKE_BUILD_TYPE RelWithDebInfo) # gcc optimizations (-O2) with debug symbols and add. warnings
# SET(CMAKE_BUILD_TYPE Release) # full gcc optimization (-O3), no debug symbols, no debug messages


### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED ILCTEST)

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )

# in order to find FindEUDAQ.cmake
LIST( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake )


FIND_PACKAGE( Marlin 1.0 REQUIRED )
FIND_PACKAGE( MarlinUtil 1.2 REQUIRED )
FIND_PACKAGE( GSL )
FIND_PACKAGE( AIDA )
FIND_PACKAGE( EUDAQ )
FIND_PACKAGE( ROOT COMPONENTS Minuit )
FIND_PACKAGE( LCCD  REQUIRED )               

FOREACH( pkg Marlin MarlinUtil GSL AIDA EUDAQ ROOT LCCD )
    IF( ${pkg}_FOUND )
        # include as "system" libraries: gcc will be less verbose w.r.t. warnings
        INCLUDE_DIRECTORIES(SYSTEM ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
    ELSE()
        MESSAGE( STATUS "${pkg} not found" )
    ENDIF()
ENDFOREACH()

INCLUDE_DIRECTORIES( ./gbl/cpp )
AUX_SOURCE_DIRECTORY( ./gbl/cpp/src gbl_sources )
ADD_LIBRARY( GBL STATIC ${gbl_sources} )
IF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
 SET_TARGET_PROPERTIES(GBL PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" ) 

### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()



### LIBRARY #################################################################

# ------ add preprocessor flags ---------------------------------------------
IF( EUDAQ_FOUND )
    ADD_DEFINITIONS( "-DUSE_EUDAQ -DEUDAQ_PLATFORM=PF_LINUX -DEUDAQ_FUNC=__PRETTY_FUNCTION__" )
    GET_FILENAME_COMPONENT( EUDAQ_LIBRARY_FULL_PATH ${EUDAQ_EUDAQ_LIBRARY} REALPATH )
    INSTALL( PROGRAMS ${EUDAQ_LIBRARY_FULL_PATH} DESTINATION lib ) # PROGRAMS is like FILES but sets executable permissions
ELSE()
    MESSAGE( STATUS "WARNING: failed to configure Eutelescope with EUDAQ!!" )
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
ADD_DEFINITIONS("-DDO_TESTING")


# ---------------------------------------------------------------------------


# include directories
INCLUDE_DIRECTORIES( ./include )
INCLUDE_DIRECTORIES( ./gbl/cpp/include )

#INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wall -ansi -pedantic -Wextra" )
ADD_DEFINITIONS( "-Wno-long-long" )

# debuging warnings and switches - only active if  CMAKE_BUILD_TYPE  is set to debug
SET(CMAKE_CXX_FLAGS_DEBUG "-g -fdiagnostics-show-option -Weffc++ -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wsign-promo -Wstrict-null-sentinel -Wswitch-default -Wundef" CACHE STRING "Debug options." FORCE )
# also useful: -Wshadow (however, GCC 4.1 uses this even for system libraries, causing many warnings from Marlin & co)
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -fdiagnostics-show-option -Weffc++ -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wsign-promo -Wstrict-null-sentinel -Wswitch-default -Wundef" CACHE STRING "Debug options." FORCE )

# add library
SET( libname ${PROJECT_NAME} )
AUX_SOURCE_DIRECTORY( ./src library_sources )
ADD_SHARED_LIBRARY( ${libname} ${library_sources} )
INSTALL_SHARED_LIBRARY( ${libname} DESTINATION lib )

TARGET_LINK_LIBRARIES( ${libname} GBL)

# used for alignment if Eutelescope was build with ROOT support
IF( ROOT_FOUND AND ROOT_MINUIT_FOUND )
    TARGET_LINK_LIBRARIES( ${libname} ${ROOT_MINUIT_LIBRARY} )
ENDIF()

MACRO( ADD_EUTELESCOPE_TOOL _name )
    ADD_EXECUTABLE( ${_name} src/exec/${_name}.cxx )
    TARGET_LINK_LIBRARIES( ${_name} ${libname} )
    INSTALL( TARGETS ${_name} DESTINATION bin )
    TARGET_LINK_LIBRARIES( ${_name} GBL) 
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


#  _            _       
# | |_ ___  ___| |_ ___ 
# | __/ _ \/ __| __/ __|
# | ||  __/\__ \ |_\__ \
#  \__\___||___/\__|___/
#=============================================================================                      

# to see the output of 'make test' in case of failure, export the environmental variable CTEST_OUTPUT_ON_FAILURE=true

ENABLE_TESTING()
INCLUDE(CTest)


# ----------------------------------------------------------------------------
# TEST IF ENVIRONMENT IS OK FOR RUNNING TESTS
# ----------------------------------------------------------------------------

IF( NOT DEFINED ENV{EUTELESCOPE} OR NOT DEFINED ENV{MARLIN})
      MESSAGE( "trying to load build and run environment for EUTelescope" )

  IF( DEFINED ENV{ILCSOFT})
        EXECUTE_PROCESS( COMMAND sh -c " [ -f $ILCSOFT/Eutelscope/HEAD/build_env.sh ] && . $ILCSOFT/Eutelescope/HEAD/build_env.sh" )

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

# !RELEASE: cppcheck is not a useful test for most users (comment out the following lines)
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
# /!RELEASE

# ======================================================================
# ======================================================================
# TestPysubExampleDaturaAlone: based on config in pysub/examples/datura-alone
# ======================================================================
# ======================================================================

# !RELEASE: the data-driven tests require DESY-AFS access and are not useful for most users (comment out the following lines)
  INCLUDE(jobsub/examples/datura-alone/testing.cmake)
  INCLUDE(pysub/examples/datura-alone/testing.cmake)
  INCLUDE(pysub/examples/anemone-2FEI4/testing.cmake)

# set custom settings (e.g. tests to be skipped for valgrind)
configure_file (
	       ./CTestCustom.cmake.in 
	       ${PROJECT_BINARY_DIR}/CTestCustom.cmake
	       )
# /!RELEASE

ELSE()
    MESSAGE( "EUTELESCOPE environment variable not defined, cannot run tests - please run 'source build_env.sh' in the EUTel root directory first." )
ENDIF()

