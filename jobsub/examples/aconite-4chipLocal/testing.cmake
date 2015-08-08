#
# This file defines a number of data-driven test based on the example configuration
# in this directory. If you have access the the files and paths defined below you
# can run the tests by running 'make test' in the build directory in the EUTelescope root
#

# ======================================================================
# ======================================================================
# TestJobsubExampleAconite-4chipLocal: based on config in jobsub/examples/aconite-4chip
# ======================================================================
# ======================================================================

    SET( testdir "${PROJECT_BINARY_DIR}/Testing/test_jobsub-aconite-4chipLocal" )
    SET( jobsubdir "$ENV{EUTELESCOPE}/jobsub" )
    SET( exampledir "${jobsubdir}/examples/aconite-4chipLocal" )

    # the run number
    SET( RunNr "1085" )
    # run number padded with leading zeros
    execute_process(COMMAND sh -c "printf %06d ${RunNr}" OUTPUT_VARIABLE PaddedRunNr)
 
    # only needed in the last step to test the results of EUTel against a set of reference files:
    SET( stattestdir "$ENV{EUTELESCOPE}/test/stattest/bin" )
    SET( referencedatadir "/afs/desy.de/group/telescopes/EutelTestData/TestExampleAconiteQuadFEI4" )

    SET( executable python -tt ${jobsubdir}/jobsub.py )
    # options: use config, use csv, change native path to central AFS location, reduce number of events to 200k
    SET( jobsubOptions --config=${exampledir}/config.cfg -csv ${exampledir}/runlist.csv -o NativePath=${referencedatadir} -o MaxRecordNumber=30000)


    # all this regular expressions must be matched for the tests to pass.
    # the order of the expressions must be matched in the test execution!
    # additional statements can be defined for each test individually
    SET( jobsub_pass_regex_1 "Now running Marlin" )
    SET( marlin_pass_regex_1 "Processing event.*in run ${PaddedRunNr}" )
    SET( jobsub_pass_regex_2 "Marlin execution done" )

    SET( generic_fail_regex "ERROR" "CRITICAL" "segmentation violation" "There were [0-9]* error messages reported")


#
#  STEP 0: PREPARE TEST DIRECTORY
#
    ADD_TEST( TestJobsubExampleAconite-4chipLocalCleanup sh -c "[ -d ${testdir} ] && rm -rf ${testdir} || echo 'no cleanup needed.'" )
    ADD_TEST( TestJobsubExampleAconite-4chipLocalSetup sh -c "mkdir -p ${testdir}/gear && mkdir -p ${testdir}/output/histograms && mkdir -p ${testdir}/output/database && mkdir -p ${testdir}/output/logs  && mkdir -p ${testdir}/output/lcio && mkdir -p ${testdir}/output/results" )
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 1: CONVERTER
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

    ADD_TEST( NAME TestJobsubExampleAconite-4chipLocalConverterRun 
              WORKING_DIRECTORY "${testdir}"
	      COMMAND ${executable} ${jobsubOptions} converter ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalConverterRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleAconite-4chipLocalSetup
	# converter step sometimes takes a bit longer (in s)
	TIMEOUT 2500
    )

    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleAconite-4chipLocalConverterLog sh -c "[ -f ${testdir}/output/logs/converter-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalConverterLog PROPERTIES DEPENDS TestJobsubExampleAconite-4chipLocalConverterRun)

    ADD_TEST( TestJobsubExampleAconite-4chipLocalConverterHisto sh -c "[ -f ${testdir}/output/histograms/run${PaddedRunNr}-converter-histo.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalConverterHisto PROPERTIES DEPENDS TestJobsubExampleAconite-4chipLocalConverterRun)

    ADD_TEST( TestJobsubExampleAconite-4chipLocalConverterM26Hotpix sh -c "[ -f ${testdir}/output/database/run${PaddedRunNr}-hotpixel-m26-db.slcio ] && lcio_check_col_elements --expelements 6  hotpixel_m26  ${testdir}/output/database/run${PaddedRunNr}-hotpixel-m26-db.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalConverterM26Hotpix PROPERTIES DEPENDS TestJobsubExampleAconite-4chipLocalConverterRun)

    ADD_TEST( TestJobsubExampleAconite-4chipLocalConverterAPIXHotpix sh -c "[ -f ${testdir}/output/database/run${PaddedRunNr}-hotpixel-apix-db.slcio ] && lcio_check_col_elements --expelements 2  hotpixel_apix  ${testdir}/output/database/run${PaddedRunNr}-hotpixel-apix-db.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalConverterAPIXHotpix PROPERTIES DEPENDS TestJobsubExampleAconite-4chipLocalConverterRun)

    ADD_TEST( TestJobsubExampleAconite-4chipLocalConverterOutputM26 sh -c "[ -f ${testdir}/output/lcio/run${PaddedRunNr}-converter.slcio ] && lcio_check_col_elements --expelements 6 zsdata_m26 ${testdir}/output/lcio/run${PaddedRunNr}-converter.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalConverterOutputM26 PROPERTIES DEPENDS TestJobsubExampleAconite-4chipLocalConverterRun)

    ADD_TEST( TestJobsubExampleAconite-4chipLocalConverterOutputAPIX sh -c "[ -f ${testdir}/output/lcio/run${PaddedRunNr}-converter.slcio ] && lcio_check_col_elements --expelements 2 zsdata_apix ${testdir}/output/lcio/run${PaddedRunNr}-converter.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalConverterOutputAPIX PROPERTIES DEPENDS TestJobsubExampleAconite-4chipLocalConverterRun)


#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 2: CLUSTERING
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    ADD_TEST( NAME TestJobsubExampleAconite-4chipLocalClusteringRun 
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubOptions} clustering ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalClusteringRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleAconite-4chipLocalConverterOutput
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleAconite-4chipLocalClusteringLog sh -c "[ -f ${testdir}/output/logs/clustering-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalClusteringLog PROPERTIES DEPENDS TestJobsubExampleAconite-4chipLocalClusteringRun)

    ADD_TEST( TestJobsubExampleAconite-4chipLocalClusteringHisto sh -c "[ -f ${testdir}/output/histograms/run${PaddedRunNr}-clustering-histo.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalClusteringHisto PROPERTIES DEPENDS TestJobsubExampleAconite-4chipLocalClusteringRun)

#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 3: HITMAKER
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    ADD_TEST( NAME TestJobsubExampleAconite-4chipLocalHitmakerRun 
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubOptions} hitmaker ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalHitmakerRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleAconite-4chipLocalClusteringOutput
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleAconite-4chipLocalHitmakerLog sh -c "[ -f ${testdir}/output/logs/hitmaker-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalHitmakerLog PROPERTIES DEPENDS TestJobsubExampleAconite-4chipLocalHitmakerRun)

    ADD_TEST( TestJobsubExampleAconite-4chipLocalHitmakerHisto sh -c "[ -f ${testdir}/output/histograms/run${PaddedRunNr}-hitmaker-histo.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalHitmakerHisto PROPERTIES DEPENDS TestJobsubExampleAconite-4chipLocalHitmakerRun)

    ADD_TEST( NAME TestJobsubExampleAconite-4chipLocalHitmakerPrealignmentGEAR COMMAND diff "${exampledir}/gear/gear001085_pre.xml" "/afs/desy.de/user/b/bisanz/public/gear001085_pre.xml")
    SET_TESTS_PROPERTIES(TestJobsubExampleAconite-4chipLocalHitmakerPrealignmentGEAR PROPERTIES FAIL_REGULAR_EXPRESSION "")

#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 4: ALIGNMENT
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    # all this regular expressions must be matched for the test to pass
    SET( align_pass_regex_1 "Initialising Mille" )
    SET( align_pass_regex_2 "Pede successfully finished" )

    ADD_TEST( NAME TestJobsubExampleAconite-4chipLocalAlignRun 
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubOptions} align ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalAlignRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${align_pass_regex_1}.*${marlin_pass_regex_1}.*${align_pass_regex_2}.*${jobsub_pass_regex_2}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleAconite-4chipLocalHitmakerOutput
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleAconite-4chipLocalAlignLog sh -c "[ -f ${testdir}/output/logs/align-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalAlignLog PROPERTIES DEPENDS TestJobsubExampleAconite-4chipLocalAlignRun)

    ADD_TEST( TestJobsubExampleAconite-4chipLocalAlignHisto sh -c "[ -f ${testdir}/output/histograms/run${PaddedRunNr}-aligndaf-histo.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalAlignHisto PROPERTIES DEPENDS TestJobsubExampleAconite-4chipLocalAlignRun)

    ADD_TEST( TestJobsubExampleAconite-4chipLocalAlignOutput sh -c "[ -f ${testdir}/output/database/run${PaddedRunNr}-align-mille.bin -a -f ${testdir}/output/database/run${PaddedRunNr}-pede-steer.txt ] " )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalAlignOutput PROPERTIES DEPENDS TestJobsubExampleAconite-4chipLocalAlignRun)

#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 5: FITTER
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    #SET( fit_pass_regex_1 "Total number of reconstructed tracks *[0-9][0-9][0-9][0-9][0-9]+" )

    ADD_TEST( NAME TestJobsubExampleAconite-4chipLocalFitterRun 
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubOptions} fitter ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalFitterRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleAconite-4chipLocalAlignRun
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleAconite-4chipLocalFitterLog sh -c "[ -f ${testdir}/output/logs/fitter-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalFitterLog PROPERTIES DEPENDS TestJobsubExampleAconite-4chipLocalFitterRun)

    ADD_TEST( TestJobsubExampleAconite-4chipLocalFitterHisto sh -c "[ -f ${testdir}/output/histograms/run${PaddedRunNr}-fitter-histo.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalFitterHisto PROPERTIES DEPENDS TestJobsubExampleAconite-4chipLocalFitterRun)

    # we expect to see between 1 and 3 tracks in every event 
    # but tolerate if this is not the case in 40% of the events (empty events are counted)
    ADD_TEST( TestJobsubExampleAconite-4chipLocalFitterOutput sh -c "[ -f ${testdir}/output/lcio/run${PaddedRunNr}-track.slcio ] && lcio_check_col_elements --pedantic --expelements 1+2-1 --abselementerror 2 --releventerror .40 track ${testdir}/output/lcio/run${PaddedRunNr}-track.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalFitterOutput PROPERTIES DEPENDS TestJobsubExampleAconite-4chipLocalFitterRun)

#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 6: StatTest
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#    SET( executable "python ${stattestdir}/runtests.py" )
#
#    # all this regular expressions must be matched for the test to pass
#    SET( fit_pass_regex_1 "SUCCESS" )
#    SET( fit_fail_regex "FAILED" "NOT PASSED" "segmentation violation")
#
#    # run stattest tool on output from previous step and test it against reference file; test are configured in specified config file (*.qa)
#
#    ADD_TEST( TestJobsubExampleAconite-4chipLocalStatTestClustering sh -c "PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH ${executable} --cdash -g ${testdir}/output/stattest_report_clus.pdf ${referencedatadir}/StatTestConf_Anemone2FEI4Clustering.qa ${testdir}/output/histograms/run${PaddedRunNr}-clu-histo.root ${referencedatadir}/run${PaddedRunNr}-clu-histo.root" )
#    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalStatTestClustering PROPERTIES
#        # test will pass if ALL of the following expressions are matched
#        PASS_REGULAR_EXPRESSION "${fit_pass_regex_1}"
#        # test will fail if ANY of the following expressions is matched 
#        FAIL_REGULAR_EXPRESSION "${fit_fail_regex}"
#	# test depends on earlier steps
#	DEPENDS TestJobsubExampleAconite-4chipLocalClusteringRun
#	)
#
#
#    ADD_TEST( TestJobsubExampleAconite-4chipLocalStatTestAlign sh -c "PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH ${executable} --cdash -g ${testdir}/output/stattest_report_align.pdf ${referencedatadir}/StatTestConf_Anemone2FEI4Align.qa ${testdir}/output/histograms/run${PaddedRunNr}-align-histo.root ${referencedatadir}/run${PaddedRunNr}-align-histo.root" )
#    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalStatTestAlign PROPERTIES
#        # test will pass if ALL of the following expressions are matched
#        PASS_REGULAR_EXPRESSION "${fit_pass_regex_1}"
#        # test will fail if ANY of the following expressions is matched 
#        FAIL_REGULAR_EXPRESSION "${fit_fail_regex}"
#	# test depends on earlier steps
#	DEPENDS TestJobsubExampleAconite-4chipLocalAlignRun
#	)
#
#
#    ADD_TEST( TestJobsubExampleAconite-4chipLocalStatTestFitter sh -c "PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH ${executable} --cdash  -g${testdir}/output/stattest_report_fitter.pdf ${referencedatadir}/StatTestConf_Anemone2FEI4Fitter.qa ${testdir}/output/histograms/run${PaddedRunNr}-track-histo.root ${referencedatadir}/run${PaddedRunNr}-track-histo.root" )
#    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalStatTestFitter PROPERTIES
#        # test will pass if ALL of the following expressions are matched
#        PASS_REGULAR_EXPRESSION "${fit_pass_regex_1}"
#        # test will fail if ANY of the following expressions is matched 
#        FAIL_REGULAR_EXPRESSION "${fit_fail_regex}"
#	# test depends on earlier steps
#	DEPENDS TestJobsubExampleAconite-4chipLocalFitterRun
#	)
#
#
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##  STEP 7: MemChecks
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
#  # STEP 1-5 VARIANTS USED FOR MEMCHECKS ONLY:
#    SET( executable python -tt ${jobsubdir}/jobsub.py )
#    # options for memcheck runs: reduced run range, plain output for valgrind parsing
#    SET( jobsubMemCheckOptions --config=${exampledir}/config.cfg -csv ${exampledir}/runlist.csv -o NativePath=${referencedatadir} -o MaxRecordNumber=2000 -o MemCheckFlag=MemCheck --plain)
#
#  # Converter run with reduced run range
#    ADD_TEST( NAME TestJobsubExampleAconite-4chipLocalConverterRunMemCheck
#              WORKING_DIRECTORY "${testdir}"
#	      COMMAND ${executable} ${jobsubMemCheckOptions} converter ${RunNr} )
#    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalConverterRunMemCheck PROPERTIES
#        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
#        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
#    )
#
#    ADD_TEST( NAME TestJobsubExampleAconite-4chipLocalClusteringRunMemCheck
#              WORKING_DIRECTORY ${testdir} 
#	      COMMAND ${executable} ${jobsubMemCheckOptions} clustering ${RunNr} )
#    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalClusteringRunMemCheck PROPERTIES
#        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
#        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
#	)
#
#    ADD_TEST( NAME TestJobsubExampleAconite-4chipLocalHitmakerRunMemCheck
#              WORKING_DIRECTORY ${testdir} 
#	      COMMAND ${executable} ${jobsubMemCheckOptions} hitmaker ${RunNr} )
#    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalHitmakerRunMemCheck PROPERTIES
#        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
#        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
#	)
#
#    ADD_TEST( NAME TestJobsubExampleAconite-4chipLocalAlignRunMemCheck
#              WORKING_DIRECTORY ${testdir} 
#	      COMMAND ${executable} ${jobsubMemCheckOptions} align ${RunNr} )
#    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalAlignRunMemCheck PROPERTIES
#        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
#        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
#	)
#
#    ADD_TEST( NAME TestJobsubExampleAconite-4chipLocalFitterRunMemCheck
#              WORKING_DIRECTORY ${testdir} 
#	      COMMAND ${executable} ${jobsubMemCheckOptions} fitter ${RunNr} )
#    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipLocalFitterRunMemCheck PROPERTIES
#        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
#        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
#	)
