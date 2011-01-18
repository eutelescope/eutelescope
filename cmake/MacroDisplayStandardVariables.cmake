
# helper macro to display standard cmake variables and force write to cache
# otherwise outdated values may appear in ccmake gui
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

