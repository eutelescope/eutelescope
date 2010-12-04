MARK_AS_ADVANCED( CMAKE_BACKWARDS_COMPATIBILITY )
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE ) # default in cmake 2.6

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

# library *nix style versioning
SET( ${PROJECT_NAME}_SOVERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}" )
SET( ${PROJECT_NAME}_VERSION   "${${PROJECT_NAME}_SOVERSION}.${${PROJECT_NAME}_VERSION_PATCH}" )

# add library install path to the rpath list
SET( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
MARK_AS_ADVANCED( CMAKE_INSTALL_RPATH )

# append link pathes to rpath list
SET( CMAKE_INSTALL_RPATH_USE_LINK_PATH 1 )
MARK_AS_ADVANCED( CMAKE_INSTALL_RPATH_USE_LINK_PATH )

# output directories
SET( EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin" )
SET( LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib" )
MARK_AS_ADVANCED( EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH )

# add uninstall target
IF( EXISTS "${PROJECT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in" )

    # create uninstall configuration file 
    CONFIGURE_FILE( "${PROJECT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
                    "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake" @ONLY )

    # add uninstall target
    ADD_CUSTOM_TARGET( uninstall "${CMAKE_COMMAND}" -P "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake" )

ENDIF( EXISTS "${PROJECT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in" )


# display standard variables and force write to cache
# otherwise outdated values may appear in: ccmake ..
MACRO( DISPLAY_STD_VARIABLES )
    MESSAGE( STATUS )
    MESSAGE( STATUS "-------------------------------------------------------------------------------" )
    MESSAGE( STATUS "Change values with: cmake -D<Variable>=<Value>" )

    IF( DEFINED CMAKE_INSTALL_PREFIX )
        MESSAGE( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
        SET( CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "Where to install ${PROJECT_NAME}" FORCE )
    ENDIF()


    IF( DEFINED CMAKE_BUILD_TYPE )
        MESSAGE( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
        SET( CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE )
    ENDIF()

    IF( DEFINED BUILD_SHARED_LIBS )
        MESSAGE( STATUS "BUILD_SHARED_LIBS = ${BUILD_SHARED_LIBS}" )
        SET( BUILD_SHARED_LIBS "${BUILD_SHARED_LIBS}" CACHE BOOL "Set to OFF to build static libraries" FORCE )
    ENDIF()

    IF( DEFINED CMAKE_MODULE_PATH )
        LIST( REMOVE_DUPLICATES CMAKE_MODULE_PATH )
        #MESSAGE( STATUS "CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}" )
        MESSAGE( STATUS "CMAKE_MODULE_PATH =" )
        FOREACH( _path ${CMAKE_MODULE_PATH} )
            MESSAGE( STATUS "   ${_path};" )
        ENDFOREACH( _path ${CMAKE_MODULE_PATH} )
        SET( CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" CACHE PATH "Path to custom CMake Modules" FORCE )
    ENDIF()

    MESSAGE( STATUS "-------------------------------------------------------------------------------" )
    MESSAGE( STATUS )

ENDMACRO( DISPLAY_STD_VARIABLES )

