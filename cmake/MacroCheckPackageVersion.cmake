##############################################################################
# macro for checking a package version
#
# this macro should be called from your PKGVersion.cmake file after the
# following variables have been set:
#   PKG_VERSION_MAJOR
#   PKG_VERSION_MINOR
#   PKG_VERSION_PATCH
#   PKG_VERSION_TWEAK
#
# note that PATCH and TWEAK are not required
#
# this macro follows the following conventions:
#
#   if FIND_PACKAGE is called with EXACT argument than the version has to
#   match EXACTLY, i.e.:
#       1.5 == 1.5
#       1.5 == 1.5.0
#       1.5 == 1.5.0.0
#       1.5.2 == 1.5.2.0
#       1.5.2.1 == 1.5.2.1
#       1.5.2 != 1.5.2.1
#       1.5 != 1.5.0.1
#
#
#   if PKG_VERSION_CHECK_MINIMUM_REQUIRED is set to TRUE than any version
#   greater or equal to the requested version is compatible, i.e. the same
#   behavior as with the cmake variable CMAKE_MINIMUM_REQUIRED, e.g.:
#       searching: 1.2     --> installed: 1.5.2.2 --> compatible
#       searching: 1.5     --> installed: 1.5.2.2 --> compatible
#       searching: 1.5.2.1 --> installed: 1.5.2.2 --> compatible
#       searching: 1.5.2.2 --> installed: 1.5.2.2 --> compatible
#       searching: 1.5.2.3 --> installed: 1.5.2.2 --> unsuitable
#       searching: 1.7     --> installed: 1.5.2.2 --> unsuitable
#
#
#   otherwise a version is compatible if and only if:
#       1. MAJOR AND MINOR versions are EQUAL
#       2. remaining version parts (PATCH.TWEAK) of installed package
#             are greater or equal than the searched version i.e.:
#          searching: 1.5     --> installed: 1.5.2.2 --> compatible
#          searching: 1.5.0   --> installed: 1.5.2.2 --> compatible
#          searching: 1.5.1   --> installed: 1.5.2.2 --> compatible
#          searching: 1.5.2   --> installed: 1.5.2.2 --> compatible
#          searching: 1.5.2.1 --> installed: 1.5.2.2 --> compatible
#          searching: 1.5.2.3 --> installed: 1.5.2.2 --> unsuitable
#          searching: 1.5.3   --> installed: 1.5.2.2 --> unsuitable
#          searching: 1.4     --> installed: 1.5.2.2 --> unsuitable
#          searching: 1.4.7.7 --> installed: 1.5.2.2 --> unsuitable
#       
#
#
# following variables are returned (internally to cmake):
#   PACKAGE_VERSION_EXACT       : set to TRUE if exact version was found
#   PACKAGE_VERSION_COMPATIBLE  : set to TRUE if version is compatible
#   PACKAGE_VERSION_UNSUITABLE  : set to TRUE if version found is unsuitable
#
#
# @author Jan Engels, Desy IT
##############################################################################

SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )

MACRO( CHECK_PACKAGE_VERSION _pkgname )

    # set option PKG_VERSION_CHECK_MINIMUM_REQUIRED per default to TRUE
    IF( NOT DEFINED ${_pkgname}_VERSION_CHECK_MINIMUM_REQUIRED )
        SET( ${_pkgname}_VERSION_CHECK_MINIMUM_REQUIRED TRUE )
    ENDIF()

    # PACKAGE_FIND_NAME is not defined on FindPKG.cmake modules, therefore
    # we need to make it an argument to the macro
    IF( NOT _pkgname AND DEFINED PACKAGE_FIND_NAME )
        SET( _pkgname "${PACKAGE_FIND_NAME}" )
    ENDIF()

    SET( _installed_major ${${_pkgname}_VERSION_MAJOR} )
    SET( _installed_minor ${${_pkgname}_VERSION_MINOR} )
    SET( _installed_version ${_installed_major}.${_installed_minor}.${${_pkgname}_VERSION_PATCH}.${${_pkgname}_VERSION_TWEAK} )

    SET( _searching_major ${${_pkgname}_FIND_VERSION_MAJOR} )
    SET( _searching_minor ${${_pkgname}_FIND_VERSION_MINOR} )
    SET( _searching_version ${${_pkgname}_FIND_VERSION} )

    # only do work if FIND_PACKAGE called with a version argument
    IF( DEFINED ${_pkgname}_FIND_VERSION )

        IF( NOT ${_pkgname}_FIND_QUIETLY )
            MESSAGE( STATUS "Check for ${_pkgname}: looking for version ${_searching_version}" )
        ENDIF()

        # these variables are evaluated internally by the cmake command FIND_PACKAGE to mark this
        # package as suitable or not depending on the required version
        SET( PACKAGE_VERSION_EXACT FALSE )
        SET( PACKAGE_VERSION_COMPATIBLE FALSE )
        SET( PACKAGE_VERSION_UNSUITABLE TRUE )

        IF( _installed_version VERSION_EQUAL _searching_version ) # if version matches EXACTLY

            IF( NOT ${_pkgname}_FIND_QUIETLY )
                MESSAGE( STATUS "Check for ${_pkgname}: exact version found: ${_installed_version}" )
            ENDIF()
            
            SET( PACKAGE_VERSION_EXACT TRUE )
            SET( PACKAGE_VERSION_UNSUITABLE FALSE )
            #SET( PACKAGE_VERSION_COMPATIBLE TRUE ) # FIXME also set COMPATIBLE if version matches EXACTLY ?!

        ELSE() # if version does not match EXACTLY, check if it is at least compatible/suitable

            IF( ${_pkgname}_FIND_VERSION_EXACT ) # user required version to match EXACTLY and it failed!
                
                IF( NOT ${_pkgname}_FIND_QUIETLY )
                    MESSAGE( STATUS "Check for ${_pkgname}: version found: ${_installed_version}" )
                    MESSAGE( STATUS "Check for ${_pkgname}: could not find exact version" )
                ENDIF()

                SET( PACKAGE_VERSION_UNSUITABLE TRUE )

            ELSE() # user did not require an EXACT match, let's check if version is compatible...


                IF( ${_pkgname}_VERSION_CHECK_MINIMUM_REQUIRED )

                    IF( NOT ${_pkgname}_FIND_QUIETLY )
                        MESSAGE( STATUS "Check for ${_pkgname}_VERSION_CHECK_MINIMUM_REQUIRED" )
                    ENDIF()

                    # --- METHOD 1: MINIMUM_REQUIRED ------------------------------------------------------
                    # this method is not as strict as method 2. it only requires that installed version is
                    # greater or equal than version searched by the user, i.e. like CMAKE_MINIMUM_REQUIRED
                    # if user asks for version 1.2.5 then any version >= 1.2.5 is suitable/compatible
                    IF( _searching_version VERSION_LESS _installed_version )
                        SET( PACKAGE_VERSION_COMPATIBLE TRUE )
                    ENDIF()
                    # -------------------------------------------------------------------------------------

                ELSE()

                    # --- METHOD 2: MAJOR.MINOR MUST BE EQUAL ---------------------------------------------
                    # check for compatible version, i.e.
                    # MAJOR and MINOR versions must be EQUAL
                    IF( _installed_major VERSION_EQUAL _searching_major AND
                        _installed_minor VERSION_EQUAL _searching_minor )

                        # if given, PATCH.TWEAK must be EQUAL or GREATER, i.e.
                        # if user asks for 1.2.6 then 1.2.7 is compatible, but 1.2.5 is unsuitable
                        IF( ${_pkgname}_FIND_VERSION_COUNT LESS 3 )
                            SET( PACKAGE_VERSION_COMPATIBLE TRUE )
                        ELSE()
                            IF( _searching_version VERSION_LESS _installed_version )
                                SET( PACKAGE_VERSION_COMPATIBLE TRUE )
                            ENDIF()
                        ENDIF()
                    ENDIF()
                    # -------------------------------------------------------------------------------------

                ENDIF()

            ENDIF()

        ENDIF()


        IF( PACKAGE_VERSION_COMPATIBLE )
            SET( PACKAGE_VERSION_UNSUITABLE FALSE )
            IF( NOT ${_pkgname}_FIND_QUIETLY )
                MESSAGE( STATUS "Check for ${_pkgname}: compatible version found: ${_installed_version}" )
            ENDIF()
        ENDIF()

        IF( PACKAGE_VERSION_UNSUITABLE )
            IF( NOT ${_pkgname}_FIND_QUIETLY )
                MESSAGE( STATUS "Check for ${_pkgname}: unsuitable version found: ${_installed_version}" )
                #MESSAGE( STATUS "Check for ${_pkgname}: minimum required version not found" )
            ENDIF()
        ENDIF()

    ENDIF( DEFINED ${_pkgname}_FIND_VERSION )

ENDMACRO( CHECK_PACKAGE_VERSION )
