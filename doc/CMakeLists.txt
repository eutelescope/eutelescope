########################################################
# cmake file for building Eutelescope documentation
# @author Jan Engels, DESY
########################################################


FIND_PACKAGE( Doxygen )

IF( DOXYGEN_FOUND )

    # build documentation out-of-source
    SET( DOC_SRC_DIR "${PROJECT_SOURCE_DIR}/doc" )
    SET( DOC_BIN_DIR "${PROJECT_BINARY_DIR}/docbuild" )

    # directories to search for documentation
    #SET( DOX_INPUT introduction.html ../include ../src )


    # custom command to build documentation
    ADD_CUSTOM_COMMAND(
        OUTPUT  "${DOC_BIN_DIR}/html/index.html"
        COMMAND DOX_PROJECT_NAME=${PROJECT_NAME}
                DOX_PROJECT_NUMBER="${${PROJECT_NAME}_VERSION}"
                DOX_OUTPUT_DIRECTORY="${DOC_BIN_DIR}"
                #DOX_INPUT="${DOX_INPUT}"
                "${DOXYGEN_EXECUTABLE}"
        WORKING_DIRECTORY "${DOC_SRC_DIR}"
        COMMENT "Building API Documentation..."
        DEPENDS Doxyfile CMakeLists.txt ../include/* ../src/*
    )

    # add doc target
    ADD_CUSTOM_TARGET( doc DEPENDS "${DOC_BIN_DIR}/html/index.html" )

ELSE()

    MESSAGE( STATUS "Doxygen not found -- INSTALL_DOC set to OFF" )
    SET( INSTALL_DOC OFF )

ENDIF()



IF( INSTALL_DOC )

    # make sure doxygen is executed (make doc) before make install
    INSTALL( CODE "EXECUTE_PROCESS( COMMAND ${CMAKE_BUILD_TOOL} doc)" )

    # in-source vs out-of-source installations
    IF( CMAKE_INSTALL_PREFIX STREQUAL "${PROJECT_SOURCE_DIR}" )

        # --- in-source installations ---
        SET( DOC_INSTALL_PREFIX doc )

        # keep source tree clean
        ADD_CUSTOM_TARGET( uninstall_doc
            COMMAND "${CMAKE_COMMAND}" -E remove_directory "${CMAKE_INSTALL_PREFIX}/doc/html"
        )
        ADD_DEPENDENCIES( uninstall uninstall_doc )

    ELSE()

        # --- out-of-source installations ---
        SET( DOC_INSTALL_PREFIX doc/${PROJECT_NAME} )
    
    ENDIF()

ENDIF()

INSTALL( DIRECTORY "${DOC_BIN_DIR}/html" DESTINATION ${DOC_INSTALL_PREFIX} )
INSTALL( FILES
        "${PROJECT_SOURCE_DIR}/steering-examples/minimal.cfg"
        "${PROJECT_SOURCE_DIR}/steering-examples/complete.cfg"
    DESTINATION ${DOC_INSTALL_PREFIX}/html
)


