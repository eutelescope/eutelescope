#
# This file defines a number of data-driven test based on the example configuration
# in this directory. If you have access the the files and paths defined below you
# can run the tests by running 'make test' in the build directory in the EUTelescope root
#

# ======================================================================
# ======================================================================
# TestJobsubExampleAconite-4chip: based on config in jobsub/examples/aconite-4chip
# ======================================================================
# ======================================================================

    SET( testdir "${PROJECT_BINARY_DIR}/Testing/test_jobsub-aconite-4chip" )
    SET( jobsubdir "$ENV{EUTELESCOPE}/jobsub" )
    SET( exampledir "${jobsubdir}/examples/aconite-4chip" )

    # the run number
    SET( RunNr "1085" )
    # run number padded with leading zeros
    execute_process(COMMAND sh -c "printf %06d ${RunNr}" OUTPUT_VARIABLE PaddedRunNr)
 
    # only needed in the last step to test the results of EUTel against a set of reference files:
    SET( stattestdir "$ENV{EUTELESCOPE}/test/stattest/bin" )
    SET( referencedatadir "/afs/desy.de/group/telescopes/EutelTestData/TestExampleAconiteQuadFEI4" )
    # SET( referencedatadir "/home/ilcsoft/EutelTestData/TestExampleAnemone2FEI4" )

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
    ADD_TEST( TestJobsubExampleAconite-4chipCleanup sh -c "[ -d ${testdir} ] && rm -rf ${testdir} || echo 'no cleanup needed.'" )
    ADD_TEST( TestJobsubExampleAconite-4chipSetup sh -c "mkdir -p ${testdir}/output/histograms  && mkdir -p ${testdir}/output/database && mkdir -p ${testdir}/output/logs  && mkdir -p ${testdir}/output/lcio && mkdir -p ${testdir}/output/results" )
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 1: CONVERTER
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

    ADD_TEST( NAME TestJobsubExampleAconite-4chipConverterRun 
              WORKING_DIRECTORY "${testdir}"
	      COMMAND ${executable} ${jobsubOptions} converter ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipConverterRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleAconite-4chipSetup
	# converter step sometimes takes a bit longer (in s)
	TIMEOUT 2500
    )


    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleAconite-4chipConverterLog sh -c "[ -f ${testdir}/output/logs/converter-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipConverterLog PROPERTIES DEPENDS TestJobsubExampleAconite-4chipConverterRun)

    ADD_TEST( TestJobsubExampleAconite-4chipConverterHisto sh -c "[ -f ${testdir}/output/histograms/run${PaddedRunNr}-converter-histo.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipConverterHisto PROPERTIES DEPENDS TestJobsubExampleAconite-4chipConverterRun)

    ADD_TEST( TestJobsubExampleAconite-4chipConverterM26Hotpix sh -c "[ -f ${testdir}/output/database/run${PaddedRunNr}-hotpixel-m26-db.slcio ] && lcio_check_col_elements --expelements 6  hotpixel_m26  ${testdir}/output/database/run${PaddedRunNr}-hotpixel-m26-db.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipConverterM26Hotpix PROPERTIES DEPENDS TestJobsubExampleAconite-4chipConverterRun)

    ADD_TEST( TestJobsubExampleAconite-4chipConverterAPIXHotpix sh -c "[ -f ${testdir}/output/database/run${PaddedRunNr}-hotpixel-apix-db.slcio ] && lcio_check_col_elements --expelements 2  hotpixel_apix  ${testdir}/output/database/run${PaddedRunNr}-hotpixel-apix-db.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipConverterAPIXHotpix PROPERTIES DEPENDS TestJobsubExampleAconite-4chipConverterRun)

    ADD_TEST( TestJobsubExampleAconite-4chipConverterOutputM26 sh -c "[ -f ${testdir}/output/lcio/run${PaddedRunNr}-converter.slcio ] && lcio_check_col_elements --expelements 6 zsdata_m26 ${testdir}/output/lcio/run${PaddedRunNr}-converter.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipConverterOutputM26 PROPERTIES DEPENDS TestJobsubExampleAconite-4chipConverterRun)

    ADD_TEST( TestJobsubExampleAconite-4chipConverterOutputAPIX sh -c "[ -f ${testdir}/output/lcio/run${PaddedRunNr}-converter.slcio ] && lcio_check_col_elements --expelements 2 zsdata_apix ${testdir}/output/lcio/run${PaddedRunNr}-converter.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipConverterOutputAPIX PROPERTIES DEPENDS TestJobsubExampleAconite-4chipConverterRun)


#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 2: CLUSTERING
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    ADD_TEST( NAME TestJobsubExampleAconite-4chipClusteringRun 
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubOptions} clustering ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipClusteringRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleAconite-4chipConverterOutput
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleAconite-4chipClusteringLog sh -c "[ -f ${testdir}/output/logs/clustering-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipClusteringLog PROPERTIES DEPENDS TestJobsubExampleAconite-4chipClusteringRun)

    ADD_TEST( TestJobsubExampleAconite-4chipClusteringHisto sh -c "[ -f ${testdir}/output/histograms/run${PaddedRunNr}-clustering-histo.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipClusteringHisto PROPERTIES DEPENDS TestJobsubExampleAconite-4chipClusteringRun)

    #TODO: FIXME!!!
    # we expect an average of 24.4 clusters per event
    #ADD_TEST( TestJobsubExampleAconite-4chipClusteringOutput sh -c "[ -f ${testdir}/output/results/run${PaddedRunNr}-clu.slcio ] && lcio_check_col_elements --average --expelements 38 --relelementerror 0.1 cluster_m26 ${testdir}/output/results/run${PaddedRunNr}-clu.slcio" )
    #SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipClusteringOutput PROPERTIES DEPENDS TestJobsubExampleAconite-4chipClusteringRun)



#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 3: HITMAKER
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    ADD_TEST( NAME TestJobsubExampleAconite-4chipHitmakerRun 
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubOptions} hitmaker ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipHitmakerRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleAconite-4chipClusteringOutput
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleAconite-4chipHitmakerLog sh -c "[ -f ${testdir}/output/logs/hitmaker-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipHitmakerLog PROPERTIES DEPENDS TestJobsubExampleAconite-4chipHitmakerRun)

    ADD_TEST( TestJobsubExampleAconite-4chipHitmakerHisto sh -c "[ -f ${testdir}/output/histograms/run${PaddedRunNr}-hitmaker-histo.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipHitmakerHisto PROPERTIES DEPENDS TestJobsubExampleAconite-4chipHitmakerRun)

    ADD_TEST( TestJobsubExampleAconite-4chipHitmakerPrealign sh -c "[ -f ${testdir}/output/database/run${PaddedRunNr}-prealign-db.slcio ] && lcio_check_col_elements --expelements 8  alignment  ${testdir}/output/database/run${PaddedRunNr}-prealign-db.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipHitmakerPrealign PROPERTIES DEPENDS TestJobsubExampleAconite-4chipHitmakerRun)

    #TODO: FIXME
    # we expect an average hit number of 24 for run 97 (wide geometry) using the example configuration
    #ADD_TEST( TestJobsubExampleAconite-4chipHitmakerOutput sh -c "[ -f ${testdir}/output/results/run${PaddedRunNr}-hit.slcio ] && lcio_check_col_elements -a --expelements 39 --relelementerror 0.1 hit ${testdir}/output/results/run${PaddedRunNr}-hit.slcio" )
    #SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipHitmakerOutput PROPERTIES DEPENDS TestJobsubExampleAconite-4chipHitmakerRun)


#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 4: ALIGNMENT
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    # all this regular expressions must be matched for the test to pass
    SET( align_pass_regex_1 "Initialising Mille" )
    SET( align_pass_regex_2 "Pede successfully finished" )

    ADD_TEST( NAME TestJobsubExampleAconite-4chipAlignRun 
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubOptions} align ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipAlignRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${align_pass_regex_1}.*${marlin_pass_regex_1}.*${align_pass_regex_2}.*${jobsub_pass_regex_2}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleAconite-4chipHitmakerOutput
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleAconite-4chipAlignLog sh -c "[ -f ${testdir}/output/logs/align-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipAlignLog PROPERTIES DEPENDS TestJobsubExampleAconite-4chipAlignRun)

    ADD_TEST( TestJobsubExampleAconite-4chipAlignHisto sh -c "[ -f ${testdir}/output/histograms/run${PaddedRunNr}-aligndaf-histo.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipAlignHisto PROPERTIES DEPENDS TestJobsubExampleAconite-4chipAlignRun)

    ADD_TEST( TestJobsubExampleAconite-4chipAlignDB sh -c "[ -f ${testdir}/output/database/run${PaddedRunNr}-alignment.slcio ] && lcio_check_col_elements --expelements 8  alignment  ${testdir}/output/database/run${PaddedRunNr}-alignment.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipAlignDB PROPERTIES DEPENDS TestJobsubExampleAconite-4chipAlignRun)

    ADD_TEST( TestJobsubExampleAconite-4chipAlignOutput sh -c "[ -f ${testdir}/output/database/run${PaddedRunNr}-align-mille.bin -a -f ${testdir}/output/database/run${PaddedRunNr}-pede-steer.txt ] " )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipAlignOutput PROPERTIES DEPENDS TestJobsubExampleAconite-4chipAlignRun)





#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 5: FITTER
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    #SET( fit_pass_regex_1 "Total number of reconstructed tracks *[0-9][0-9][0-9][0-9][0-9]+" )

    ADD_TEST( NAME TestJobsubExampleAconite-4chipFitterRun 
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubOptions} fitter ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipFitterRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleAconite-4chipAlignRun
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleAconite-4chipFitterLog sh -c "[ -f ${testdir}/output/logs/fitter-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipFitterLog PROPERTIES DEPENDS TestJobsubExampleAconite-4chipFitterRun)

    ADD_TEST( TestJobsubExampleAconite-4chipFitterHisto sh -c "[ -f ${testdir}/output/histograms/run${PaddedRunNr}-fitter-histo.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipFitterHisto PROPERTIES DEPENDS TestJobsubExampleAconite-4chipFitterRun)

    # we expect to see between 1 and 3 tracks in every event 
    # but tolerate if this is not the case in 40% of the events (empty events are counted)
    ADD_TEST( TestJobsubExampleAconite-4chipFitterOutput sh -c "[ -f ${testdir}/output/lcio/run${PaddedRunNr}-track.slcio ] && lcio_check_col_elements --pedantic --expelements 1+2-1 --abselementerror 2 --releventerror .40 track ${testdir}/output/lcio/run${PaddedRunNr}-track.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipFitterOutput PROPERTIES DEPENDS TestJobsubExampleAconite-4chipFitterRun)

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
#    ADD_TEST( TestJobsubExampleAconite-4chipStatTestClustering sh -c "PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH ${executable} --cdash -g ${testdir}/output/stattest_report_clus.pdf ${referencedatadir}/StatTestConf_Anemone2FEI4Clustering.qa ${testdir}/output/histograms/run${PaddedRunNr}-clu-histo.root ${referencedatadir}/run${PaddedRunNr}-clu-histo.root" )
#    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipStatTestClustering PROPERTIES
#        # test will pass if ALL of the following expressions are matched
#        PASS_REGULAR_EXPRESSION "${fit_pass_regex_1}"
#        # test will fail if ANY of the following expressions is matched 
#        FAIL_REGULAR_EXPRESSION "${fit_fail_regex}"
#	# test depends on earlier steps
#	DEPENDS TestJobsubExampleAconite-4chipClusteringRun
#	)
#
#
#    ADD_TEST( TestJobsubExampleAconite-4chipStatTestAlign sh -c "PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH ${executable} --cdash -g ${testdir}/output/stattest_report_align.pdf ${referencedatadir}/StatTestConf_Anemone2FEI4Align.qa ${testdir}/output/histograms/run${PaddedRunNr}-align-histo.root ${referencedatadir}/run${PaddedRunNr}-align-histo.root" )
#    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipStatTestAlign PROPERTIES
#        # test will pass if ALL of the following expressions are matched
#        PASS_REGULAR_EXPRESSION "${fit_pass_regex_1}"
#        # test will fail if ANY of the following expressions is matched 
#        FAIL_REGULAR_EXPRESSION "${fit_fail_regex}"
#	# test depends on earlier steps
#	DEPENDS TestJobsubExampleAconite-4chipAlignRun
#	)
#
#
#    ADD_TEST( TestJobsubExampleAconite-4chipStatTestFitter sh -c "PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH ${executable} --cdash  -g${testdir}/output/stattest_report_fitter.pdf ${referencedatadir}/StatTestConf_Anemone2FEI4Fitter.qa ${testdir}/output/histograms/run${PaddedRunNr}-track-histo.root ${referencedatadir}/run${PaddedRunNr}-track-histo.root" )
#    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipStatTestFitter PROPERTIES
#        # test will pass if ALL of the following expressions are matched
#        PASS_REGULAR_EXPRESSION "${fit_pass_regex_1}"
#        # test will fail if ANY of the following expressions is matched 
#        FAIL_REGULAR_EXPRESSION "${fit_fail_regex}"
#	# test depends on earlier steps
#	DEPENDS TestJobsubExampleAconite-4chipFitterRun
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
#    ADD_TEST( NAME TestJobsubExampleAconite-4chipConverterRunMemCheck
#              WORKING_DIRECTORY "${testdir}"
#	      COMMAND ${executable} ${jobsubMemCheckOptions} converter ${RunNr} )
#    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipConverterRunMemCheck PROPERTIES
#        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
#        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
#    )
#
#    ADD_TEST( NAME TestJobsubExampleAconite-4chipClusteringRunMemCheck
#              WORKING_DIRECTORY ${testdir} 
#	      COMMAND ${executable} ${jobsubMemCheckOptions} clustering ${RunNr} )
#    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipClusteringRunMemCheck PROPERTIES
#        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
#        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
#	)
#
#    ADD_TEST( NAME TestJobsubExampleAconite-4chipHitmakerRunMemCheck
#              WORKING_DIRECTORY ${testdir} 
#	      COMMAND ${executable} ${jobsubMemCheckOptions} hitmaker ${RunNr} )
#    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipHitmakerRunMemCheck PROPERTIES
#        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
#        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
#	)
#
#    ADD_TEST( NAME TestJobsubExampleAconite-4chipAlignRunMemCheck
#              WORKING_DIRECTORY ${testdir} 
#	      COMMAND ${executable} ${jobsubMemCheckOptions} align ${RunNr} )
#    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipAlignRunMemCheck PROPERTIES
#        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
#        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
#	)
#
#    ADD_TEST( NAME TestJobsubExampleAconite-4chipFitterRunMemCheck
#              WORKING_DIRECTORY ${testdir} 
#	      COMMAND ${executable} ${jobsubMemCheckOptions} fitter ${RunNr} )
#    SET_TESTS_PROPERTIES (TestJobsubExampleAconite-4chipFitterRunMemCheck PROPERTIES
#        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
#        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
#	)
