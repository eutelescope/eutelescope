FUNCTION(GET_TEST_REGEX INP OUTPUT_PASS OUTPUT_FAIL)
    IF(SLC_TESTS)
        FILE(STRINGS ${INP} OUTPUT_PASS_ REGEX "#PASS_SLC ")
        FILE(STRINGS ${INP} OUTPUT_FAIL_ REGEX "#FAIL_SLC ")
    ENDIF()

    IF(NOT OUTPUT_PASS_)
        FILE(STRINGS ${INP} OUTPUT_PASS_ REGEX "#PASS ")
    ENDIF()
    IF(NOT OUTPUT_FAIL_)
        FILE(STRINGS ${INP} OUTPUT_FAIL_ REGEX "#FAIL ")
    ENDIF()

    # Check for number of arguments - should only be one:
    LIST(LENGTH OUTPUT_PASS_ LISTCOUNT_PASS)
    LIST(LENGTH OUTPUT_FAIL_ LISTCOUNT_FAIL)
    IF(LISTCOUNT_PASS GREATER 1)
        MESSAGE(FATAL_ERROR "More than one PASS expressions defined in test ${INP}")
    ENDIF()
    IF(LISTCOUNT_FAIL GREATER 1)
        MESSAGE(FATAL_ERROR "More than one FAIL expressions defined in test ${INP}")
    ENDIF()

    # Escape possible regex patterns in the expected output:
    ESCAPE_REGEX("${OUTPUT_PASS_}" OUTPUT_PASS_)
    ESCAPE_REGEX("${OUTPUT_FAIL_}" OUTPUT_FAIL_)

    SET(${OUTPUT_PASS} "${OUTPUT_PASS_}" PARENT_SCOPE)
    SET(${OUTPUT_FAIL} "${OUTPUT_FAIL_}" PARENT_SCOPE)
ENDFUNCTION()

FUNCTION(ESCAPE_REGEX INP OUTPUT)
    # Escape possible regex patterns in the expected output:
    STRING(REPLACE "#PASS " "" _TMP_STR "${INP}")
    STRING(REPLACE "#FAIL " "" _TMP_STR "${_TMP_STR}")
    STRING(REPLACE "#PASS_SLC " "" _TMP_STR "${_TMP_STR}")
    STRING(REPLACE "#FAIL_SLC " "" _TMP_STR "${_TMP_STR}")
#    STRING(REPLACE "\\" "\\\\" _TMP_STR "${_TMP_STR}")
#    STRING(REPLACE "?" "\\?" _TMP_STR "${_TMP_STR}")
#    STRING(REPLACE "+" "\\+" _TMP_STR "${_TMP_STR}")
#    STRING(REPLACE "*" "\\*" _TMP_STR "${_TMP_STR}")
#    STRING(REPLACE "(" "\\\(" _TMP_STR "${_TMP_STR}")
#    STRING(REPLACE ")" "\\\)" _TMP_STR "${_TMP_STR}")
#    STRING(REPLACE "[" "\\\[" _TMP_STR "${_TMP_STR}")
#    STRING(REPLACE "]" "\\\]" _TMP_STR "${_TMP_STR}")
    SET(${OUTPUT} "${_TMP_STR}" PARENT_SCOPE)
ENDFUNCTION()

FUNCTION(ADD_EUTEL_TEST TEST SUBDIR)
    # Allow the test to specify additional CLI parameters:
    FILE(STRINGS ${TEST} OPTS REGEX "#OPTION ")
    FOREACH(OPT ${OPTS})
        STRING(REPLACE "#OPTION " "" OPT "${OPT}")
        SET(CLIOPTIONS "${CLIOPTIONS} --option ${OPT}")
    ENDFOREACH()
    
    
    FILE(STRINGS ${TEST} CONF REGEX "#CONFIG ")
    STRING(REPLACE "#CONFIG " "" CONF "${CONF}")

    FILE(STRINGS ${TEST} RUNLIST REGEX "#RUNLIST ")
    STRING(REPLACE "#RUNLIST " "" RUNLIST "${RUNLIST}")

    FILE(STRINGS ${TEST} STEP REGEX "#STEP ")
    STRING(REPLACE "#STEP " "" STEP "${STEP}")

    FILE(STRINGS ${TEST} RUNNO REGEX "#RUNNUMBER ")
    STRING(REPLACE "#RUNNUMBER " "" RUNNO "${RUNNO}")


    ADD_TEST(NAME ${TEST}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/make_output_dir_structure.sh "${SUBDIR}" "jobsub -c $ENV{EUTELESCOPE}/${CONF} -csv $ENV{EUTELESCOPE}/${RUNLIST} ${STEP} ${RUNNO} ${CLIOPTIONS}"
    )

    # Parse configuration file for pass/fail conditions:
    GET_TEST_REGEX(${TEST} EXPRESSIONS_PASS EXPRESSIONS_FAIL)
    IF(EXPRESSIONS_PASS)
        SET_TESTS_PROPERTIES(${TEST} PROPERTIES PASS_REGULAR_EXPRESSION "${EXPRESSIONS_PASS}")
    ENDIF()
    IF(EXPRESSIONS_FAIL)
        SET_TESTS_PROPERTIES(${TEST} PROPERTIES FAIL_REGULAR_EXPRESSION "${EXPRESSIONS_FAIL}")
    ENDIF()

    # Some tests might depend on others:
    FILE(STRINGS ${TEST} DEPENDENCY REGEX "#DEPENDS ")
    IF(DEPENDENCY)
        STRING(REPLACE "#DEPENDS " "" DEPENDENCY "${DEPENDENCY}")
        SET_TESTS_PROPERTIES(${TEST} PROPERTIES DEPENDS "${DEPENDENCY}")
    ENDIF()

    # Add individual timeout criteria:
    FILE(STRINGS ${TEST} TESTTIMEOUT REGEX "#TIMEOUT ")
    IF(TESTTIMEOUT)
        STRING(REPLACE "#TIMEOUT " "" TESTTIMEOUT "${TESTTIMEOUT}")
        SET_TESTS_PROPERTIES(${TEST} PROPERTIES TIMEOUT "${TESTTIMEOUT}")
    ENDIF()

ENDFUNCTION()

##############################
#  aconite-4chipLocal tests  #
##############################
OPTION(TEST_ACONITE4CHIP-LOCAL "Perform aconite-4chipLocal tests?" ON)

IF(TEST_ACONITE4CHIP-LOCAL)
    FILE(GLOB TEST_LIST_MODULES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} aconite-4chipLocal_test/test_*)
    LIST(LENGTH TEST_LIST_MODULES NUM_TEST_MODULES)
    MESSAGE(STATUS "aconite-4chipLocal tests: ${NUM_TEST_MODULES} tests")
    FOREACH(TEST ${TEST_LIST_MODULES})
        ADD_EUTEL_TEST(${TEST} aconite-4chipLocal)
    ENDFOREACH()
ELSE()
    MESSAGE(STATUS "aconite-4chipLocal tests: tests deactivated.")
ENDIF()

##############################
#  gbl_local tests           #
##############################
OPTION(TEST_GBL-LOCAL "Perform gbl_local tests?" ON)

IF(TEST_GBL-LOCAL)
    FILE(GLOB TEST_LIST_MODULES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} gbl_local_test/test_*)
    LIST(LENGTH TEST_LIST_MODULES NUM_TEST_MODULES)
    MESSAGE(STATUS "gbl_local tests: ${NUM_TEST_MODULES} tests")
    FOREACH(TEST ${TEST_LIST_MODULES})
        ADD_EUTEL_TEST(${TEST} gbl_local)
    ENDFOREACH()
ELSE()
    MESSAGE(STATUS "gbl_local tests: tests deactivated.")
ENDIF()
