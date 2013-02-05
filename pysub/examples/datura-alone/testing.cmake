#
# This file defines a number of data-driven tests based on the example configuration
# in this directory. If you have access the the files and paths defined below you
# can run the tests by running 'make test' in the build directory in the EUTelescope root
#

# ======================================================================
# ======================================================================
# TestPysubExampleDaturaAlone: based on config in pysub/examples/datura-alone
# ======================================================================
# ======================================================================

    SET( testdir "${PROJECT_BINARY_DIR}/Testing/test_datura-alone" )
    SET( pysubdir "$ENV{EUTELESCOPE}/pysub" )
    SET( exampledir "${pysubdir}/examples/datura-alone" )

    SET( RunNr "4118" )
 
   # only needed in the last step to test the results of EUTel against a set of reference files:
    SET( stattestdir "$ENV{EUTELESCOPE}/test/stattest/bin" )
    SET( referencedatadir "/afs/desy.de/group/telescopes/EutelTestData/TestExampleDaturaAlone" )


#
#  STEP 0: PREPARE TEST DIRECTORY
#
	ADD_TEST( TestPysubExampleDaturaAloneCleanup sh -c "[ -d ${testdir} ] && rm -rf ${testdir} || echo 'no cleanup needed.'" )
	ADD_TEST( TestPysubExampleDaturaAloneSetup sh -c "mkdir -p ${testdir}/output/histo && mkdir -p ${testdir}/output/results && mkdir -p ${testdir}/output/db && mkdir -p ${testdir}/output/logs && mkdir -p ${testdir}/output/lcio-raw" )
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 1: CONVERTER
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# todo:
#  - currently running only over 10000 events (see steering template); should this number be adjusted and/or tested for (using lcio_event_counter)?
#  - set property DEPENDS to make sure that failure of required steps does not cause following steps to fail as well 

    SET( executable "submit-converter.py" )

    # all this regular expressions must be matched for the test to pass
    SET( converter_pass_regex_1 "Now processing run 0*${RunNr}" )
    SET( converter_pass_regex_2 "Running Marlin" )
    SET( converter_pass_regex_3 "Marlin finished successfully" )

    SET( converter_fail_regex "Skipping to the next run" "ERROR" "CRITICAL" "segmentation violation")

    ADD_TEST( TestPysubExampleDaturaAloneConverterRun sh -c "cd ${testdir} && python ${pysubdir}/${executable} --config=${exampledir}/config.cfg --hot ${RunNr} ${RunNr}" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneConverterRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${converter_pass_regex_1}.*${converter_pass_regex_2}.*${converter_pass_regex_3}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${converter_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestPysubExampleDaturaAloneSetup
    )
    # now check if the expected output files exist and look ok
    ADD_TEST( TestPysubExampleDaturaAloneConverterLog sh -c "[ -f ${testdir}/output/logs/converter-`printf %06d ${RunNr}`.tar.gz ]" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneConverterLog PROPERTIES DEPENDS TestPysubExampleDaturaAloneConverterRun)

    ADD_TEST( TestPysubExampleDaturaAloneConverterHisto sh -c "[ -f ${testdir}/output/histo/run`printf %06d ${RunNr}`-converter-histo.root ]" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneConverterHisto PROPERTIES DEPENDS TestPysubExampleDaturaAloneConverterRun)

    ADD_TEST( TestPysubExampleDaturaAloneConverterHotpix sh -c "[ -f ${testdir}/output/db/run${RunNr}-hotpixel-db.slcio ] && lcio_check_col_elements --expelements 6  hotpixel_m26  ${testdir}/output/db/run${RunNr}-hotpixel-db.slcio" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneConverterHotpix PROPERTIES DEPENDS TestPysubExampleDaturaAloneConverterRun)

    ADD_TEST( TestPysubExampleDaturaAloneConverterOutput sh -c "[ -f ${testdir}/output/lcio-raw/run`printf %06d ${RunNr}`.slcio ] && lcio_check_col_elements --expelements 6 zsdata_m26 ${testdir}/output/lcio-raw/run`printf %06d ${RunNr}`.slcio" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneConverterOutput PROPERTIES DEPENDS TestPysubExampleDaturaAloneConverterRun)




#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 2: CLUSTERING
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    SET( executable "submit-clusearch.py" )

    # all this regular expressions must be matched for the test to pass
    SET( clusearch_pass_regex_1 "Now processing run 0*${RunNr}" )
    SET( clusearch_pass_regex_2 "Running Marlin" )
    SET( clusearch_pass_regex_3 "Marlin finished successfully" )

    SET( clusearch_fail_regex "Skipping to the next run" "ERROR" "CRITICAL" "segmentation violation")

    ADD_TEST( TestPysubExampleDaturaAloneClusearchRun sh -c "cd ${testdir} && python ${pysubdir}/${executable} --config=${exampledir}/config.cfg --hot ${RunNr} ${RunNr}" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneClusearchRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${clusearch_pass_regex_1}.*${clusearch_pass_regex_2}.*${clusearch_pass_regex_3}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${clusearch_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestPysubExampleDaturaAloneConverterOutput
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestPysubExampleDaturaAloneClusearchLog sh -c "[ -f ${testdir}/output/logs/clusearch-`printf %06d ${RunNr}`.tar.gz ]" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneClusearchLog PROPERTIES DEPENDS TestPysubExampleDaturaAloneClusearchRun)

    ADD_TEST( TestPysubExampleDaturaAloneClusearchHisto sh -c "[ -f ${testdir}/output/histo/run`printf %06d ${RunNr}`-clu-histo.root ]" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneClusearchHisto PROPERTIES DEPENDS TestPysubExampleDaturaAloneClusearchRun)

    ADD_TEST( TestPysubExampleDaturaAloneClusearchOffset sh -c "[ -f ${testdir}/output/db/run`printf %06d ${RunNr}`-offset-db.slcio ] && lcio_check_col_elements --expelements 6  preAlignment  ${testdir}/output/db/run`printf %06d ${RunNr}`-offset-db.slcio" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneClusearchOffset PROPERTIES DEPENDS TestPysubExampleDaturaAloneClusearchRun)


    ADD_TEST( TestPysubExampleDaturaAloneClusearchOutput sh -c "[ -f ${testdir}/output/results/run`printf %06d ${RunNr}`-clu-p.slcio ] && lcio_check_col_elements --expelements 6 zsdata_m26 ${testdir}/output/results/run`printf %06d ${RunNr}`-clu-p.slcio" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneClusearchOutput PROPERTIES DEPENDS TestPysubExampleDaturaAloneClusearchRun)



#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 3: HITMAKER
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    SET( executable "submit-hitmaker.py" )

    # all this regular expressions must be matched for the test to pass
    SET( hitmaker_pass_regex_1 "Running Marlin" )
    SET( hitmaker_pass_regex_2 "Successfully finished" )
    SET( hitmaker_pass_regex_3 "Marlin finished successfully" )

    SET( hitmaker_fail_regex "Skipping to the next run" "ERROR" "CRITICAL" "segmentation violation")

    ADD_TEST( TestPysubExampleDaturaAloneHitmakerRun sh -c "cd ${testdir} && python ${pysubdir}/${executable} --config=${exampledir}/config.cfg -o ${RunNr} run`printf %06d ${RunNr}`-clu-p.slcio" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneHitmakerRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${hitmaker_pass_regex_1}.*${hitmaker_pass_regex_2}.*${hitmaker_pass_regex_3}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${hitmaker_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestPysubExampleDaturaAloneClusearchOutput
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestPysubExampleDaturaAloneHitmakerLog sh -c "[ -f ${testdir}/output/logs/hitmaker-${RunNr}.tar.gz ]" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneHitmakerLog PROPERTIES DEPENDS TestPysubExampleDaturaAloneHitmakerRun)

    ADD_TEST( TestPysubExampleDaturaAloneHitmakerHisto sh -c "[ -f ${testdir}/output/histo/${RunNr}-hit-histo.root ]" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneHitmakerHisto PROPERTIES DEPENDS TestPysubExampleDaturaAloneHitmakerRun)

    ADD_TEST( TestPysubExampleDaturaAloneHitmakerPrealign sh -c "[ -f ${testdir}/output/db/${RunNr}-prealign-db.slcio ] && lcio_check_col_elements --expelements 6  alignment  ${testdir}/output/db/${RunNr}-prealign-db.slcio" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneHitmakerPrealign PROPERTIES DEPENDS TestPysubExampleDaturaAloneHitmakerRun)

    # we expect an average hit number of 35 for run 4118 using the example configuration
    ADD_TEST( TestPysubExampleDaturaAloneHitmakerOutput sh -c "[ -f ${testdir}/output/results/${RunNr}-hit.slcio ] && lcio_check_col_elements -a --expelements 35 hit ${testdir}/output/results/${RunNr}-hit.slcio" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneHitmakerOutput PROPERTIES DEPENDS TestPysubExampleDaturaAloneHitmakerRun)


#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 4: ALIGNMENT
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    SET( executable "submit-align.py" )

    # all this regular expressions must be matched for the test to pass
    SET( align_pass_regex_1 "Checking the input file ${RunNr}-hit.slcio" )
    SET( align_pass_regex_2 "Running Marlin" )
    SET( align_pass_regex_3 "Pede successfully finished" )
    SET( align_pass_regex_4 "Marlin finished successfully" )

    SET( align_fail_regex "Skipping to the next run" "ERROR" "CRITICAL" "segmentation violation")

    ADD_TEST( TestPysubExampleDaturaAloneAlignRun sh -c "cd ${testdir} && python ${pysubdir}/${executable} --config=${exampledir}/config.cfg -o ${RunNr} ${RunNr}-hit.slcio" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneAlignRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${align_pass_regex_1}.*${align_pass_regex_2}.*${align_pass_regex_3}.*${align_pass_regex_4}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${align_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestPysubExampleDaturaAloneHitmakerOutput
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestPysubExampleDaturaAloneAlignLog sh -c "[ -f ${testdir}/output/logs/align-${RunNr}.tar.gz ]" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneAlignLog PROPERTIES DEPENDS TestPysubExampleDaturaAloneAlignRun)

    ADD_TEST( TestPysubExampleDaturaAloneAlignHisto sh -c "[ -f ${testdir}/output/histo/${RunNr}-align-histo.root ]" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneAlignHisto PROPERTIES DEPENDS TestPysubExampleDaturaAloneAlignRun)

    ADD_TEST( TestPysubExampleDaturaAloneAlignDB sh -c "[ -f ${testdir}/output/db/${RunNr}-align-db.slcio ] && lcio_check_col_elements --expelements 6  alignment  ${testdir}/output/db/${RunNr}-prealign-db.slcio" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneAlignDB PROPERTIES DEPENDS TestPysubExampleDaturaAloneAlignRun)

    ADD_TEST( TestPysubExampleDaturaAloneAlignOutput sh -c "[ -f ${testdir}/output/results/${RunNr}-align-mille.bin -a -f ${testdir}/output/results/${RunNr}-pede-steer.txt ] " )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneAlignOutput PROPERTIES DEPENDS TestPysubExampleDaturaAloneAlignRun)







#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 5: FITTER
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    SET( executable "submit-fitter.py" )

    # all this regular expressions must be matched for the test to pass
    SET( fit_pass_regex_1 "Processing run header 1" )
    SET( fit_pass_regex_2 "Total number of reconstructed tracks *[0-9][0-9][0-9][0-9][0-9]+" )
    SET( fit_pass_regex_3 "Marlin finished successfully" )

    SET( fit_fail_regex "Skipping to the next run" "ERROR" "CRITICAL" "segmentation violation")

    ADD_TEST( TestPysubExampleDaturaAloneFitterRun sh -c "cd ${testdir} && python ${pysubdir}/${executable} --config=${exampledir}/config.cfg -o ${RunNr} ${RunNr}-hit.slcio -a ${RunNr}-align-db.slcio" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneFitterRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${fit_pass_regex_1}.*${fit_pass_regex_2}.*${fit_pass_regex_3}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${fit_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestPysubExampleDaturaAloneAlignRun
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestPysubExampleDaturaAloneFitterLog sh -c "[ -f ${testdir}/output/logs/fitter-${RunNr}.tar.gz ]" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneFitterLog PROPERTIES DEPENDS TestPysubExampleDaturaAloneFitterRun)

    ADD_TEST( TestPysubExampleDaturaAloneFitterHisto sh -c "[ -f ${testdir}/output/histo/${RunNr}-track-histo.root ]" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneFitterHisto PROPERTIES DEPENDS TestPysubExampleDaturaAloneFitterRun)

    # we expect to see between 1 and 7 tracks in every event 
    # but tolerate if this is not the case in 15% of the events (empty events are not counted!)
    ADD_TEST( TestPysubExampleDaturaAloneFitterOutput sh -c "[ -f ${testdir}/output/results/${RunNr}-track.slcio ] && lcio_check_col_elements --pedantic --expelements 4 --abselementerror 3 --releventerror .15 track ${testdir}/output/results/${RunNr}-track.slcio" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneFitterOutput PROPERTIES DEPENDS TestPysubExampleDaturaAloneFitterRun)

#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 6: StatTest
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

    # TODO: ask Andrea Dotti or other Geant4 people for their FindStatTest.cmake 

    SET( executable "runtests.py" )

    # all this regular expressions must be matched for the test to pass
    SET( fit_pass_regex_1 "SUCCESS" )
    SET( fit_fail_regex "FAILED" "NOT PASSED" "segmentation violation")

    # run stattest tool on output from previous step and test it against reference file; test are configured in specified config file (*.qa)

    ADD_TEST( TestPysubExampleDaturaAloneStatTestAlign sh -c "PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH python ${stattestdir}/${executable} --cdash -g ${testdir}/output/stattest_report_align.pdf ${referencedatadir}/StatTestConf_DaturaAloneAlign.qa ${testdir}/output/histo/${RunNr}-align-histo.root ${referencedatadir}/${RunNr}-align-histo.root" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneStatTestAlign PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${fit_pass_regex_1}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${fit_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestPysubExampleDaturaAloneAlignRun
	)


    ADD_TEST( TestPysubExampleDaturaAloneStatTestFitter sh -c "PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH python ${stattestdir}/${executable} --cdash  -g${testdir}/output/stattest_report_fitter.pdf ${referencedatadir}/StatTestConf_DaturaAloneFitter.qa ${testdir}/output/histo/${RunNr}-track-histo.root ${referencedatadir}/${RunNr}-track-histo.root" )
    SET_TESTS_PROPERTIES (TestPysubExampleDaturaAloneStatTestFitter PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${fit_pass_regex_1}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${fit_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestPysubExampleDaturaAloneFitterRun
	)

