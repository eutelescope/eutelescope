#
# This file defines a number of data-driven tests based on the example configuration
# in this directory. If you have access the the files and paths defined below you
# can run the tests by running 'make test' in the build directory in the EUTelescope root
#

# ======================================================================
# ======================================================================
# TestJobsubExampleDaturaAlone: based on config in jobsub/examples/datura-alone
# ======================================================================
# ======================================================================

    SET( testdir "${PROJECT_BINARY_DIR}/Testing/test_datura-alone" )
    SET( jobsubdir "$ENV{EUTELESCOPE}/jobsub" )
    SET( exampledir "${jobsubdir}/examples/datura-alone" )

    # the run number
    SET( RunNr "4118" )
    # run number padded with leading zeros
    execute_process(COMMAND sh -c "printf %06d ${RunNr}" OUTPUT_VARIABLE PaddedRunNr)

    SET( executable "jobsub.py" )
    SET( jobsubOptions --config=${exampledir}/config.cfg -csv ${exampledir}/table_orgmode.csv)

   # only needed in the last step to test the results of EUTel against a set of reference files:
    SET( stattestdir "$ENV{EUTELESCOPE}/test/stattest/bin" )
    SET( referencedatadir "/afs/desy.de/group/telescopes/EutelTestData/TestExampleDaturaAlone" )


#
#  STEP 0: PREPARE TEST DIRECTORY
#
	ADD_TEST( TestJobsubExampleDaturaAloneCleanup sh -c "[ -d ${testdir} ] && rm -rf ${testdir} || echo 'no cleanup needed.'" )
	ADD_TEST( TestJobsubExampleDaturaAloneSetup sh -c "mkdir -p ${testdir}/output/histo && mkdir -p ${testdir}/output/results && mkdir -p ${testdir}/output/db && mkdir -p ${testdir}/output/logs && mkdir -p ${testdir}/output/lcio-raw" )
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 1: CONVERTER
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# todo:
#  - currently running only over 10000 events (see steering template); should this number be adjusted and/or tested for (using lcio_event_counter)?
#  - set property DEPENDS to make sure that failure of required steps does not cause following steps to fail as well 

    # all this regular expressions must be matched for the test to pass
    SET( converter_pass_regex_1 "Now running Marlin" )
    SET( converter_pass_regex_2 "Processing event.*in run ${PaddedRunNr}" )
    SET( converter_pass_regex_3 "Marlin execution done" )

    SET( converter_fail_regex "ERROR" "CRITICAL" "segmentation violation")

    ADD_TEST( NAME TestJobsubExampleDaturaAloneConverterRun 
              WORKING_DIRECTORY "${testdir}"
	      COMMAND python ${jobsubdir}/${executable} ${jobsubOptions} converter ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneConverterRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${converter_pass_regex_1}.*${converter_pass_regex_2}.*${converter_pass_regex_3}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${converter_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleDaturaAloneSetup
	# converter step sometimes takes a bit longer (in s)
	TIMEOUT 2500
    )

    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleDaturaAloneConverterLog sh -c "[ -f ${testdir}/output/logs/converter-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneConverterLog PROPERTIES DEPENDS TestJobsubExampleDaturaAloneConverterRun)

    ADD_TEST( TestJobsubExampleDaturaAloneConverterHisto sh -c "[ -f ${testdir}/output/histo/run${PaddedRunNr}-converter-histo.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneConverterHisto PROPERTIES DEPENDS TestJobsubExampleDaturaAloneConverterRun)

    ADD_TEST( TestJobsubExampleDaturaAloneConverterHotpix sh -c "[ -f ${testdir}/output/db/run${PaddedRunNr}-hotpixel-db.slcio ] && lcio_check_col_elements --expelements 6  hotpixel_m26  ${testdir}/output/db/run${PaddedRunNr}-hotpixel-db.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneConverterHotpix PROPERTIES DEPENDS TestJobsubExampleDaturaAloneConverterRun)

    ADD_TEST( TestJobsubExampleDaturaAloneConverterOutput sh -c "[ -f ${testdir}/output/lcio-raw/run${PaddedRunNr}.slcio ] && lcio_check_col_elements --expelements 6 zsdata_m26 ${testdir}/output/lcio-raw/run${PaddedRunNr}.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneConverterOutput PROPERTIES DEPENDS TestJobsubExampleDaturaAloneConverterRun)




#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 2: CLUSTERING
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

    # all this regular expressions must be matched for the test to pass
    SET( clusearch_pass_regex_1 "Now running Marlin" )
    SET( clusearch_pass_regex_2 "Processing event.*in run ${PaddedRunNr}" )
    SET( clusearch_pass_regex_3 "Marlin execution done" )

    SET( clusearch_fail_regex "ERROR" "CRITICAL" "segmentation violation")

    ADD_TEST( NAME TestJobsubExampleDaturaAloneClusearchRun 
              WORKING_DIRECTORY ${testdir} 
	      COMMAND python ${jobsubdir}/${executable} ${jobsubOptions} clusearch ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneClusearchRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${clusearch_pass_regex_1}.*${clusearch_pass_regex_2}.*${clusearch_pass_regex_3}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${clusearch_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleDaturaAloneConverterOutput
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleDaturaAloneClusearchLog sh -c "[ -f ${testdir}/output/logs/clusearch-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneClusearchLog PROPERTIES DEPENDS TestJobsubExampleDaturaAloneClusearchRun)

    ADD_TEST( TestJobsubExampleDaturaAloneClusearchHisto sh -c "[ -f ${testdir}/output/histo/run${PaddedRunNr}-clu-histo.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneClusearchHisto PROPERTIES DEPENDS TestJobsubExampleDaturaAloneClusearchRun)

    ADD_TEST( TestJobsubExampleDaturaAloneClusearchOffset sh -c "[ -f ${testdir}/output/db/run${PaddedRunNr}-offset-db.slcio ] && lcio_check_col_elements --expelements 6  preAlignment  ${testdir}/output/db/run${PaddedRunNr}-offset-db.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneClusearchOffset PROPERTIES DEPENDS TestJobsubExampleDaturaAloneClusearchRun)

    ADD_TEST( TestJobsubExampleDaturaAloneClusearchOutput sh -c "[ -f ${testdir}/output/results/run${PaddedRunNr}-clu.slcio ] && lcio_check_col_elements --average --expelements 35 cluster_m26 ${testdir}/output/results/run${PaddedRunNr}-clu.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneClusearchOutput PROPERTIES DEPENDS TestJobsubExampleDaturaAloneClusearchRun)



#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 3: HITMAKER
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    # all this regular expressions must be matched for the test to pass
    SET( hitmaker_pass_regex_1 "Now running Marlin" )
    SET( hitmaker_pass_regex_2 "Processing event.*in run ${PaddedRunNr}" )
    SET( hitmaker_pass_regex_3 "Marlin execution done" )

    SET( hitmaker_fail_regex "ERROR" "CRITICAL" "segmentation violation")

    ADD_TEST( NAME TestJobsubExampleDaturaAloneHitmakerRun 
              WORKING_DIRECTORY ${testdir} 
	      COMMAND python ${jobsubdir}/${executable} ${jobsubOptions} hitmaker ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneHitmakerRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${hitmaker_pass_regex_1}.*${hitmaker_pass_regex_2}.*${hitmaker_pass_regex_3}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${hitmaker_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleDaturaAloneClusearchOutput
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleDaturaAloneHitmakerLog sh -c "[ -f ${testdir}/output/logs/hitmaker-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneHitmakerLog PROPERTIES DEPENDS TestJobsubExampleDaturaAloneHitmakerRun)

    ADD_TEST( TestJobsubExampleDaturaAloneHitmakerHisto sh -c "[ -f ${testdir}/output/histo/run${PaddedRunNr}-hit-histo.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneHitmakerHisto PROPERTIES DEPENDS TestJobsubExampleDaturaAloneHitmakerRun)

    ADD_TEST( TestJobsubExampleDaturaAloneHitmakerPrealign sh -c "[ -f ${testdir}/output/db/run${PaddedRunNr}-prealign-db.slcio ] && lcio_check_col_elements --expelements 6  alignment  ${testdir}/output/db/run${PaddedRunNr}-prealign-db.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneHitmakerPrealign PROPERTIES DEPENDS TestJobsubExampleDaturaAloneHitmakerRun)

    # we expect an average hit number of 35 for run 4118 using the example configuration
    ADD_TEST( TestJobsubExampleDaturaAloneHitmakerOutput sh -c "[ -f ${testdir}/output/results/run${PaddedRunNr}-hit.slcio ] && lcio_check_col_elements -a --expelements 35 hit ${testdir}/output/results/run${PaddedRunNr}-hit.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneHitmakerOutput PROPERTIES DEPENDS TestJobsubExampleDaturaAloneHitmakerRun)


#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 4: ALIGNMENT
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    # all this regular expressions must be matched for the test to pass
    SET( align_pass_regex_1 "Now running Marlin" )
    SET( align_pass_regex_2 "Initialising Mille" )
    SET( align_pass_regex_3 "Processing event.*in run ${PaddedRunNr}" )
    SET( align_pass_regex_4 "Pede successfully finished" )
    SET( align_pass_regex_5 "Marlin execution done" )

    SET( align_fail_regex "ERROR" "CRITICAL" "segmentation violation")

    ADD_TEST( NAME TestJobsubExampleDaturaAloneAlignRun 
              WORKING_DIRECTORY ${testdir} 
	      COMMAND python ${jobsubdir}/${executable} ${jobsubOptions} align ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneAlignRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${align_pass_regex_1}.*${align_pass_regex_2}.*${align_pass_regex_3}.*${align_pass_regex_4}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${align_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleDaturaAloneHitmakerOutput
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleDaturaAloneAlignLog sh -c "[ -f ${testdir}/output/logs/align-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneAlignLog PROPERTIES DEPENDS TestJobsubExampleDaturaAloneAlignRun)

    ADD_TEST( TestJobsubExampleDaturaAloneAlignHisto sh -c "[ -f ${testdir}/output/histo/run${PaddedRunNr}-align-histo.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneAlignHisto PROPERTIES DEPENDS TestJobsubExampleDaturaAloneAlignRun)

    ADD_TEST( TestJobsubExampleDaturaAloneAlignDB sh -c "[ -f ${testdir}/output/db/run${PaddedRunNr}-align-db.slcio ] && lcio_check_col_elements --expelements 6  alignment  ${testdir}/output/db/run${PaddedRunNr}-align-db.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneAlignDB PROPERTIES DEPENDS TestJobsubExampleDaturaAloneAlignRun)

    ADD_TEST( TestJobsubExampleDaturaAloneAlignOutput sh -c "[ -f ${testdir}/output/results/run${PaddedRunNr}-align-mille.bin -a -f ${testdir}/output/results/run${PaddedRunNr}-pede-steer.txt ] " )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneAlignOutput PROPERTIES DEPENDS TestJobsubExampleDaturaAloneAlignRun)







#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 5: FITTER
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    # all this regular expressions must be matched for the test to pass
    SET( fit_pass_regex_1 "Processing run header 1" )
    #SET( fit_pass_regex_2 "Total number of reconstructed tracks *[0-9][0-9][0-9][0-9][0-9]+" )
    SET( fit_pass_regex_3 "Marlin execution done" )

    SET( fit_fail_regex "ERROR" "CRITICAL" "segmentation violation")

    ADD_TEST( NAME TestJobsubExampleDaturaAloneFitterRun 
              WORKING_DIRECTORY ${testdir} 
	      COMMAND python ${jobsubdir}/${executable} ${jobsubOptions} fitter ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneFitterRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${fit_pass_regex_1}.*${fit_pass_regex_3}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${fit_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleDaturaAloneAlignRun
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleDaturaAloneFitterLog sh -c "[ -f ${testdir}/output/logs/fitter-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneFitterLog PROPERTIES DEPENDS TestJobsubExampleDaturaAloneFitterRun)

    ADD_TEST( TestJobsubExampleDaturaAloneFitterHisto sh -c "[ -f ${testdir}/output/histo/run${PaddedRunNr}-track-histo.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneFitterHisto PROPERTIES DEPENDS TestJobsubExampleDaturaAloneFitterRun)

    # we expect to see between 1 and 7 tracks in every event 
    # but tolerate if this is not the case in 15% of the events (empty events are not counted!)
    ADD_TEST( TestJobsubExampleDaturaAloneFitterOutput sh -c "[ -f ${testdir}/output/results/run${PaddedRunNr}-track.slcio ] && lcio_check_col_elements --pedantic --expelements 4 --abselementerror 3 --releventerror .15 track ${testdir}/output/results/run${PaddedRunNr}-track.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneFitterOutput PROPERTIES DEPENDS TestJobsubExampleDaturaAloneFitterRun)

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

    ADD_TEST( TestJobsubExampleDaturaAloneStatTestAlign sh -c "PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH python ${stattestdir}/${executable} --cdash -g ${testdir}/output/stattest_report_align.pdf ${referencedatadir}/StatTestConf_DaturaAloneAlign.qa ${testdir}/output/histo/run${PaddedRunNr}-align-histo.root ${referencedatadir}/${RunNr}-align-histo.root" )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneStatTestAlign PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${fit_pass_regex_1}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${fit_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleDaturaAloneAlignRun
	)


    ADD_TEST( TestJobsubExampleDaturaAloneStatTestFitter sh -c "PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH python ${stattestdir}/${executable} --cdash  -g${testdir}/output/stattest_report_fitter.pdf ${referencedatadir}/StatTestConf_DaturaAloneFitter.qa ${testdir}/output/histo/run${PaddedRunNr}-track-histo.root ${referencedatadir}/${RunNr}-track-histo.root" )
    SET_TESTS_PROPERTIES (TestJobsubExampleDaturaAloneStatTestFitter PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${fit_pass_regex_1}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${fit_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleDaturaAloneFitterRun
	)

