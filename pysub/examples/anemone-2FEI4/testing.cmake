#
# This file defines a number of data-driven tests based on the example configuration
# in this directory. If you have access the the files and paths defined below you
# can run the tests by running 'make test' in the build directory in the EUTelescope root
#

# ======================================================================
# ======================================================================
# TestPysubExampleAnemone2FEI4: based on config in pysub/examples/anemone-2FEI4
# ======================================================================
# ======================================================================

    SET( testdir "${PROJECT_BINARY_DIR}/Testing/test_anemone-2FEI4" )
    SET( pysubdir "$ENV{EUTELESCOPE}/pysub" )
    SET( exampledir "${pysubdir}/examples/anemone-2FEI4" )

    SET( RunNr "17" )
 
   # only needed in the last step to test the results of EUTel against a set of reference files:
    SET( stattestdir "$ENV{EUTELESCOPE}/test/stattest/bin" )
    SET( referencedatadir "/afs/desy.de/group/telescopes/EutelTestData/TestPysubExampleAnemone2FEI4" )
#    SET( referencedatadir "/home/ilcsoft/TestBeam/tests/anemone2FEI4/native" )


#
#  STEP 0: PREPARE TEST DIRECTORY
#
	ADD_TEST( TestPysubExampleAnemone2FEI4Cleanup sh -c "[ -d ${testdir} ] && rm -rf ${testdir} || echo 'no cleanup needed.'" )
	ADD_TEST( TestPysubExampleAnemone2FEI4Setup sh -c "mkdir -p ${testdir}/output/histo && mkdir -p ${testdir}/output/results && mkdir -p ${testdir}/output/db && mkdir -p ${testdir}/output/logs && mkdir -p ${testdir}/output/lcio-raw" )
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

    ADD_TEST( TestPysubExampleAnemone2FEI4ConverterRun sh -c "cd ${testdir} && python ${pysubdir}/${executable} --config=${exampledir}/config.cfg --hot ${RunNr} ${RunNr}" )
    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4ConverterRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${converter_pass_regex_1}.*${converter_pass_regex_2}.*${converter_pass_regex_3}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${converter_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestPysubExampleAnemone2FEI4Setup
    )
    # now check if the expected output files exist and look ok
    ADD_TEST( TestPysubExampleAnemone2FEI4ConverterLog sh -c "[ -f ${testdir}/output/logs/converter-`printf %06d ${RunNr}`.tar.gz ]" )
    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4ConverterLog PROPERTIES DEPENDS TestPysubExampleAnemone2FEI4ConverterRun)

    ADD_TEST( TestPysubExampleAnemone2FEI4ConverterHisto sh -c "[ -f ${testdir}/output/histo/run`printf %06d ${RunNr}`-converter-histo.root ]" )
    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4ConverterHisto PROPERTIES DEPENDS TestPysubExampleAnemone2FEI4ConverterRun)

    ADD_TEST( TestPysubExampleAnemone2FEI4ConverterHotpix sh -c "[ -f ${testdir}/output/db/run${RunNr}-hotpixel-db.slcio ] && lcio_check_col_elements --expelements 6  hotpixel_m26  ${testdir}/output/db/run${RunNr}-hotpixel-db.slcio" )
    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4ConverterHotpix PROPERTIES DEPENDS TestPysubExampleAnemone2FEI4ConverterRun)

    ADD_TEST( TestPysubExampleAnemone2FEI4ConverterOutput sh -c "[ -f ${testdir}/output/lcio-raw/run`printf %06d ${RunNr}`.slcio ] && lcio_check_col_elements --expelements 6 zsdata_m26 ${testdir}/output/lcio-raw/run`printf %06d ${RunNr}`.slcio" )
    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4ConverterOutput PROPERTIES DEPENDS TestPysubExampleAnemone2FEI4ConverterRun)




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

    ADD_TEST( TestPysubExampleAnemone2FEI4ClusearchRun sh -c "cd ${testdir} && python ${pysubdir}/${executable} --config=${exampledir}/config.cfg --hot ${RunNr} ${RunNr}" )
    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4ClusearchRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${clusearch_pass_regex_1}.*${clusearch_pass_regex_2}.*${clusearch_pass_regex_3}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${clusearch_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestPysubExampleAnemone2FEI4ConverterOutput
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestPysubExampleAnemone2FEI4ClusearchLog sh -c "[ -f ${testdir}/output/logs/clusearch-`printf %06d ${RunNr}`.tar.gz ]" )
    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4ClusearchLog PROPERTIES DEPENDS TestPysubExampleAnemone2FEI4ClusearchRun)

    ADD_TEST( TestPysubExampleAnemone2FEI4ClusearchHisto sh -c "[ -f ${testdir}/output/histo/run`printf %06d ${RunNr}`-clu-histo.root ]" )
    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4ClusearchHisto PROPERTIES DEPENDS TestPysubExampleAnemone2FEI4ClusearchRun)

    # check if output file is ok and if cluster collections have expected number
    # of elements (+/- 10%) on average
    ADD_TEST( TestPysubExampleAnemone2FEI4ClusearchOutput sh -c "[ -f ${testdir}/output/results/run`printf %06d ${RunNr}`-clu-p.slcio ] && lcio_check_col_elements -a -x 35 --relelementerror 0.1 cluster_m26 ${testdir}/output/results/run`printf %06d ${RunNr}`-clu-p.slcio && lcio_check_col_elements -a -x 2 --relelementerror 0.1 cluster_apix ${testdir}/output/results/run`printf %06d ${RunNr}`-clu-p.slcio" )
    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4ClusearchOutput PROPERTIES DEPENDS TestPysubExampleAnemone2FEI4ClusearchRun)



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

    ADD_TEST( TestPysubExampleAnemone2FEI4HitmakerRun sh -c "cd ${testdir} && python ${pysubdir}/${executable} --config=${exampledir}/config.cfg -o ${RunNr} run`printf %06d ${RunNr}`-clu-p.slcio" )
    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4HitmakerRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${hitmaker_pass_regex_1}.*${hitmaker_pass_regex_2}.*${hitmaker_pass_regex_3}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${hitmaker_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestPysubExampleAnemone2FEI4ClusearchOutput
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestPysubExampleAnemone2FEI4HitmakerLog sh -c "[ -f ${testdir}/output/logs/hitmaker-${RunNr}.tar.gz ]" )
    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4HitmakerLog PROPERTIES DEPENDS TestPysubExampleAnemone2FEI4HitmakerRun)

    ADD_TEST( TestPysubExampleAnemone2FEI4HitmakerHisto sh -c "[ -f ${testdir}/output/histo/${RunNr}-hit-histo.root ]" )
    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4HitmakerHisto PROPERTIES DEPENDS TestPysubExampleAnemone2FEI4HitmakerRun)

    ADD_TEST( TestPysubExampleAnemone2FEI4HitmakerPrealign sh -c "[ -f ${testdir}/output/db/${RunNr}-prealign-db.slcio ] && lcio_check_col_elements --expelements 8  alignment  ${testdir}/output/db/${RunNr}-prealign-db.slcio" )
    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4HitmakerPrealign PROPERTIES DEPENDS TestPysubExampleAnemone2FEI4HitmakerRun)

    # we expect an average hit number of 35 for run 4118 using the example configuration
#    ADD_TEST( TestPysubExampleAnemone2FEI4HitmakerOutput sh -c "[ -f ${testdir}/output/results/${RunNr}-hit.slcio ] && lcio_check_col_elements --expelements 40 --relelementerror .5 --releventerror .4  hit ${testdir}/output/results/${RunNr}-hit.slcio" )
#    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4HitmakerOutput PROPERTIES DEPENDS TestPysubExampleAnemone2FEI4HitmakerRun)


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

    ADD_TEST( TestPysubExampleAnemone2FEI4AlignRun sh -c "cd ${testdir} && python ${pysubdir}/${executable} --config=${exampledir}/config.cfg -o ${RunNr} ${RunNr}-hit.slcio" )
    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4AlignRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${align_pass_regex_1}.*${align_pass_regex_2}.*${align_pass_regex_3}.*${align_pass_regex_4}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${align_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestPysubExampleAnemone2FEI4HitmakerOutput
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestPysubExampleAnemone2FEI4AlignLog sh -c "[ -f ${testdir}/output/logs/align-${RunNr}.tar.gz ]" )
    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4AlignLog PROPERTIES DEPENDS TestPysubExampleAnemone2FEI4AlignRun)

    ADD_TEST( TestPysubExampleAnemone2FEI4AlignHisto sh -c "[ -f ${testdir}/output/histo/${RunNr}-align-histo.root ]" )
    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4AlignHisto PROPERTIES DEPENDS TestPysubExampleAnemone2FEI4AlignRun)

    ADD_TEST( TestPysubExampleAnemone2FEI4AlignDB sh -c "[ -f ${testdir}/output/db/${RunNr}-align-db.slcio ] && lcio_check_col_elements --expelements 8  alignment  ${testdir}/output/db/${RunNr}-prealign-db.slcio" )
    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4AlignDB PROPERTIES DEPENDS TestPysubExampleAnemone2FEI4AlignRun)

    ADD_TEST( TestPysubExampleAnemone2FEI4AlignOutput sh -c "[ -f ${testdir}/output/results/${RunNr}-align-mille.bin -a -f ${testdir}/output/results/${RunNr}-pede-steer.txt ] " )
    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4AlignOutput PROPERTIES DEPENDS TestPysubExampleAnemone2FEI4AlignRun)







#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 5: FITTER
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    SET( executable "submit-fitter.py" )

    # all this regular expressions must be matched for the test to pass
    SET( fit_pass_regex_1 "Processing run header [0-9]" )
    SET( fit_pass_regex_2 "Number of fitted tracks: *[0-9][0-9][0-9][0-9][0-9]+" )
    SET( fit_pass_regex_3 "Marlin finished successfully" )

    SET( fit_fail_regex "Skipping to the next run" "ERROR" "CRITICAL" "segmentation violation")

    ADD_TEST( TestPysubExampleAnemone2FEI4FitterRun sh -c "cd ${testdir} && python ${pysubdir}/${executable} --config=${exampledir}/config.cfg -o ${RunNr} ${RunNr}-hit.slcio -a ${RunNr}-align-db.slcio" )
    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4FitterRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${fit_pass_regex_1}.*${fit_pass_regex_2}.*${fit_pass_regex_3}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${fit_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestPysubExampleAnemone2FEI4AlignRun
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestPysubExampleAnemone2FEI4FitterLog sh -c "[ -f ${testdir}/output/logs/fitter-${RunNr}.tar.gz ]" )
    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4FitterLog PROPERTIES DEPENDS TestPysubExampleAnemone2FEI4FitterRun)

    ADD_TEST( TestPysubExampleAnemone2FEI4FitterHisto sh -c "[ -f ${testdir}/output/histo/${RunNr}-track-histo.root ]" )
    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4FitterHisto PROPERTIES DEPENDS TestPysubExampleAnemone2FEI4FitterRun)

    # we expect to see between 1 and 7 tracks in every event 
    # but tolerate if this is not the case in 15% of the events (empty events are not counted!)
    ADD_TEST( TestPysubExampleAnemone2FEI4FitterOutput sh -c "[ -f ${testdir}/output/results/${RunNr}-track.slcio ] && lcio_check_col_elements --pedantic --expelements 4 --abselementerror 3 --releventerror .15 track ${testdir}/output/results/${RunNr}-track.slcio" )
    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4FitterOutput PROPERTIES DEPENDS TestPysubExampleAnemone2FEI4FitterRun)

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

    ADD_TEST( TestPysubExampleAnemone2FEI4StatTestAlign sh -c "PYTHONPATH=$ROOTSYS/lib python ${stattestdir}/${executable} -g ${testdir}/output/stattest_report_align.pdf ${referencedatadir}/StatTestConf_AnemoneFEI4Align.qa ${testdir}/output/histo/${RunNr}-align-histo.root ${referencedatadir}/${RunNr}-align-histo.root" )
    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4StatTestAlign PROPERTIES
    # test will pass if ALL of the following expressions are matched
    PASS_REGULAR_EXPRESSION "${fit_pass_regex_1}"
    # test will fail if ANY of the following expressions is matched 
    FAIL_REGULAR_EXPRESSION "${fit_fail_regex}"
    # test depends on earlier steps
    DEPENDS TestPysubExampleAnemone2FEI4AlignRun
    )


    ADD_TEST( TestPysubExampleAnemone2FEI4StatTestFitter sh -c "PYTHONPATH=$ROOTSYS/lib python ${stattestdir}/${executable} -g ${testdir}/output/stattest_report_fitter.pdf ${referencedatadir}/StatTestConf_AnemoneFEI4Fitter.qa ${testdir}/output/histo/${RunNr}-track-histo.root ${referencedatadir}/${RunNr}-track-histo.root" )
    SET_TESTS_PROPERTIES (TestPysubExampleAnemone2FEI4StatTestFitter PROPERTIES
    # test will pass if ALL of the following expressions are matched
    PASS_REGULAR_EXPRESSION "${fit_pass_regex_1}"
    # test will fail if ANY of the following expressions is matched 
    FAIL_REGULAR_EXPRESSION "${fit_fail_regex}"
    # test depends on earlier steps
    DEPENDS TestPysubExampleAnemone2FEI4FitterRun
    )

