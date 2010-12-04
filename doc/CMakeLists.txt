########################################################
# cmake file for building Eutelescope documentation
# @author Jan Engels, Desy IT
########################################################


FIND_PACKAGE( Doxygen )

IF( DOXYGEN_FOUND )

    # FIXME should build documentation out-of-source
    SET( DOC_SRC_DIR "${PROJECT_SOURCE_DIR}/doc" )

    ADD_CUSTOM_COMMAND(
        OUTPUT  "${DOC_SRC_DIR}/html/index.html"
        COMMAND "${DOXYGEN_EXECUTABLE}"
        WORKING_DIRECTORY "${DOC_SRC_DIR}"
        COMMENT "Building API Documentation..."
    )

    # add doc target
    ADD_CUSTOM_TARGET( doc DEPENDS "${DOC_SRC_DIR}/html/index.html" )

ELSE()

    MESSAGE( STATUS "Doxygen not found -- INSTALL_DOC set to OFF" )
    SET( INSTALL_DOC OFF )

ENDIF()

IF( INSTALL_DOC )

    # make sure doxygen is executed (make doc) before make install
    INSTALL( CODE "EXECUTE_PROCESS( COMMAND ${CMAKE_BUILD_TOOL} doc)" )

    # install documentation
    INSTALL( DIRECTORY "${DOC_SRC_DIR}" DESTINATION .
        PATTERN "*/.svn*" EXCLUDE
        PATTERN "*.txt" EXCLUDE
        PATTERN "Doxyfile" EXCLUDE
    )

ENDIF()

