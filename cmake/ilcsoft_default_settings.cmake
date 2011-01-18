MARK_AS_ADVANCED( CMAKE_BACKWARDS_COMPATIBILITY )
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE ) # default in cmake 2.6

# DEPRECATED
#OPTION( BUILD_SHARED_LIBS "Set to OFF to build static libraries" ON )

# include some helper macros
INCLUDE( MacroAddSharedLibrary )
INCLUDE( MacroInstallDirectory )
INCLUDE( MacroDisplayStandardVariables )
#INCLUDE( MacroGeneratePackageConfigFiles )


# set default install prefix to project root directory
# instead of the cmake default /usr/local
IF( CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" )
    SET( CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}" )
ENDIF()

# set default cmake build type to RelWithDebInfo
# possible options are: None Debug Release RelWithDebInfo MinSizeRel
IF( NOT CMAKE_BUILD_TYPE )
    SET( CMAKE_BUILD_TYPE "RelWithDebInfo" )
ENDIF()

# enable ctest
ENABLE_TESTING()
INCLUDE(CTest)
MARK_AS_ADVANCED( DART_TESTING_TIMEOUT )


# library *nix style versioning
SET( ${PROJECT_NAME}_SOVERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}" )
SET( ${PROJECT_NAME}_VERSION   "${${PROJECT_NAME}_SOVERSION}.${${PROJECT_NAME}_VERSION_PATCH}" )


# output directories
SET( EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin" )
SET( LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib" )
MARK_AS_ADVANCED( EXECUTABLE_OUTPUT_PATH )
MARK_AS_ADVANCED( LIBRARY_OUTPUT_PATH )

# add library install path to the rpath list
SET( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
MARK_AS_ADVANCED( CMAKE_INSTALL_RPATH )

# append link pathes to rpath list
SET( CMAKE_INSTALL_RPATH_USE_LINK_PATH 1 )
MARK_AS_ADVANCED( CMAKE_INSTALL_RPATH_USE_LINK_PATH )

GET_FILENAME_COMPONENT( _current_dir ${CMAKE_CURRENT_LIST_FILE} PATH )

IF( EXISTS "${_current_dir}/cmake_uninstall.cmake.in" )

    # create uninstall configuration file 
    CONFIGURE_FILE( "${_current_dir}/cmake_uninstall.cmake.in"
                    "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake" @ONLY )

    # add uninstall target
    ADD_CUSTOM_TARGET( uninstall "${CMAKE_COMMAND}" -P "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake" )

ENDIF( EXISTS "${_current_dir}/cmake_uninstall.cmake.in" )

