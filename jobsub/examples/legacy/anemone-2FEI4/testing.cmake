#
# This file defines a number of data-driven test based on the example configuration
# in this directory. If you have access the the files and paths defined below you
# can run the tests by running 'make test' in the build directory in the EUTelescope root
#

# ======================================================================
# ======================================================================
# TestJobsubExampleAnemone2FEI4: based on config in jobsub/examples/anemone-2FEI4
# ======================================================================
# ======================================================================

    SET( testdir "${PROJECT_BINARY_DIR}/Testing/test_jobsub-anemone-2FEI4" )
    SET( jobsubdir "$ENV{EUTELESCOPE}/jobsub" )
    SET( exampledir "${jobsubdir}/examples/anemone-2FEI4" )

    # the run number
    SET( RunNr "17" )
    # run number padded with leading zeros
    execute_process(COMMAND sh -c "printf %06d ${RunNr}" OUTPUT_VARIABLE PaddedRunNr)
 
   # only needed in the last step to test the results of EUTel against a set of reference files:
    SET( stattestdir "$ENV{EUTELESCOPE}/test/stattest/bin" )
    SET( referencedatadir "/afs/desy.de/group/telescopes/EutelTestData/TestExampleAnemone2FEI4" )
#    SET( referencedatadir "/home/ilcsoft/EutelTestData/TestExampleAnemone2FEI4" )

    SET( executable python -tt ${jobsubdir}/jobsub.py )
    # options: use config, use csv, change native path to central AFS location, reduce number of events to 200k
    SET( jobsubOptions --config=${exampledir}/config.cfg -csv ${exampledir}/runlist.csv -o NativePath=${referencedatadir} -o MaxRecordNumber=100000)


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
    ADD_TEST( TestJobsubExampleAnemone2FEI4Cleanup sh -c "[ -d ${testdir} ] && rm -rf ${testdir} || echo 'no cleanup needed.'" )
    ADD_TEST( TestJobsubExampleAnemone2FEI4Setup sh -c "mkdir -p ${testdir}/output/histograms  && mkdir -p ${testdir}/output/database && mkdir -p ${testdir}/output/logs  && mkdir -p ${testdir}/output/lcio && mkdir -p ${testdir}/output/results" )
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 1: CONVERTER
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

    ADD_TEST( NAME TestJobsubExampleAnemone2FEI4ConverterRun 
              WORKING_DIRECTORY "${testdir}"
	      COMMAND ${executable} ${jobsubOptions} converter ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4ConverterRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleAnemone2FEI4Setup
	# converter step sometimes takes a bit longer (in s)
	TIMEOUT 2500
    )


    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleAnemone2FEI4ConverterLog sh -c "[ -f ${testdir}/output/logs/converter-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4ConverterLog PROPERTIES DEPENDS TestJobsubExampleAnemone2FEI4ConverterRun)

    ADD_TEST( TestJobsubExampleAnemone2FEI4ConverterHisto sh -c "[ -f ${testdir}/output/histograms/run${PaddedRunNr}-converter-histo.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4ConverterHisto PROPERTIES DEPENDS TestJobsubExampleAnemone2FEI4ConverterRun)

    ADD_TEST( TestJobsubExampleAnemone2FEI4ConverterHotpixM26 sh -c "[ -f ${testdir}/output/database/run${PaddedRunNr}-hotpixel-m26-db.slcio ] && lcio_check_col_elements --expelements 6  hotpixel_m26  ${testdir}/output/database/run${PaddedRunNr}-hotpixel-m26-db.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4ConverterHotpixM26 PROPERTIES DEPENDS TestJobsubExampleAnemone2FEI4ConverterRun)

    ADD_TEST( TestJobsubExampleAnemone2FEI4ConverterOutputM26 sh -c "[ -f ${testdir}/output/lcio/run${PaddedRunNr}-converter.slcio ] && lcio_check_col_elements --expelements 6 zsdata_m26 ${testdir}/output/lcio/run${PaddedRunNr}-converter.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4ConverterOutputM26 PROPERTIES DEPENDS TestJobsubExampleAnemone2FEI4ConverterRun)

    ADD_TEST( TestJobsubExampleAnemone2FEI4ConverterHotpixAPIX sh -c "[ -f ${testdir}/output/database/run${PaddedRunNr}-hotpixel-apix-db.slcio ] && lcio_check_col_elements --expelements 2  hotpixel_apix  ${testdir}/output/database/run${PaddedRunNr}-hotpixel-apix-db.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4ConverterHotpixAPIX PROPERTIES DEPENDS TestJobsubExampleAnemone2FEI4ConverterRun)

    ADD_TEST( TestJobsubExampleAnemone2FEI4ConverterOutputAPIX sh -c "[ -f ${testdir}/output/lcio/run${PaddedRunNr}-converter.slcio ] && lcio_check_col_elements --expelements 2 zsdata_apix ${testdir}/output/lcio/run${PaddedRunNr}-converter.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4ConverterOutputAPIX PROPERTIES DEPENDS TestJobsubExampleAnemone2FEI4ConverterRun)


#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 2: CLUSTERING
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    ADD_TEST( NAME TestJobsubExampleAnemone2FEI4ClusteringRun 
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubOptions} clustering ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4ClusteringRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleAnemone2FEI4ConverterOutput
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleAnemone2FEI4ClusteringLog sh -c "[ -f ${testdir}/output/logs/clustering-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4ClusteringLog PROPERTIES DEPENDS TestJobsubExampleAnemone2FEI4ClusteringRun)

    ADD_TEST( TestJobsubExampleAnemone2FEI4ClusteringHisto sh -c "[ -f ${testdir}/output/histograms/run${PaddedRunNr}-clustering.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4ClusteringHisto PROPERTIES DEPENDS TestJobsubExampleAnemone2FEI4ClusteringRun)

    # we expect an average of 24.4 clusters per event for M26
    ADD_TEST( TestJobsubExampleAnemone2FEI4ClusteringOutputM26 sh -c "[ -f ${testdir}/output/lcio/run${PaddedRunNr}-clustering.slcio ] && lcio_check_col_elements --average --expelements 38 --relelementerror 0.1 cluster_m26_free ${testdir}/output/lcio/run${PaddedRunNr}-clustering.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4ClusteringOutputM26 PROPERTIES DEPENDS TestJobsubExampleAnemone2FEI4ClusteringRun)

    # we expect an average of 2.1 clusters per event for APIX
    ADD_TEST( TestJobsubExampleAnemone2FEI4ClusteringOutputAPIX sh -c "[ -f ${testdir}/output/lcio/run${PaddedRunNr}-clustering.slcio ] && lcio_check_col_elements --average --expelements 2 --relelementerror 0.2 cluster_apix_free ${testdir}/output/lcio/run${PaddedRunNr}-clustering.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4ClusteringOutputAPIX PROPERTIES DEPENDS TestJobsubExampleAnemone2FEI4ClusteringRun)

#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 3: HITMAKER
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    ADD_TEST( NAME TestJobsubExampleAnemone2FEI4HitmakerRun 
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubOptions} hitmaker ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4HitmakerRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleAnemone2FEI4ClusteringOutput
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleAnemone2FEI4HitmakerLog sh -c "[ -f ${testdir}/output/logs/hitmaker-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4HitmakerLog PROPERTIES DEPENDS TestJobsubExampleAnemone2FEI4HitmakerRun)

    ADD_TEST( TestJobsubExampleAnemone2FEI4HitmakerHisto sh -c "[ -f ${testdir}/output/histograms/run${PaddedRunNr}-hitmaker-histo.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4HitmakerHisto PROPERTIES DEPENDS TestJobsubExampleAnemone2FEI4HitmakerRun)

    ADD_TEST( TestJobsubExampleAnemone2FEI4HitmakerPrealign sh -c "[ -f ${testdir}/output/database/run${PaddedRunNr}-prealign-db.slcio ] && lcio_check_col_elements --expelements 8  alignment  ${testdir}/output/database/run${PaddedRunNr}-prealign-db.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4HitmakerPrealign PROPERTIES DEPENDS TestJobsubExampleAnemone2FEI4HitmakerRun)

    # we expect an average hit number of 24 for run 97 (wide geometry) using the example configuration
    ADD_TEST( TestJobsubExampleAnemone2FEI4HitmakerOutput sh -c "[ -f ${testdir}/output/lcio/run${PaddedRunNr}-hitmaker.slcio ] && lcio_check_col_elements -a --expelements 39 --relelementerror 0.1 hit ${testdir}/output/lcio/run${PaddedRunNr}-hitmaker.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4HitmakerOutput PROPERTIES DEPENDS TestJobsubExampleAnemone2FEI4HitmakerRun)


#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 4: ALIGNMENT
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    # all this regular expressions must be matched for the test to pass
    SET( align_pass_regex_1 "Initialising Mille" )
    SET( align_pass_regex_2 "Pede successfully finished" )

    ADD_TEST( NAME TestJobsubExampleAnemone2FEI4AlignRun 
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubOptions} align ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4AlignRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${align_pass_regex_1}.*${marlin_pass_regex_1}.*${align_pass_regex_2}.*${jobsub_pass_regex_2}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleAnemone2FEI4HitmakerOutput
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleAnemone2FEI4AlignLog sh -c "[ -f ${testdir}/output/logs/align-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4AlignLog PROPERTIES DEPENDS TestJobsubExampleAnemone2FEI4AlignRun)

    ADD_TEST( TestJobsubExampleAnemone2FEI4AlignHisto sh -c "[ -f ${testdir}/output/histograms/run${PaddedRunNr}-alignmentdaf.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4AlignHisto PROPERTIES DEPENDS TestJobsubExampleAnemone2FEI4AlignRun)

    ADD_TEST( TestJobsubExampleAnemone2FEI4AlignDB sh -c "[ -f ${testdir}/output/database/run${PaddedRunNr}-alignment.slcio ] && lcio_check_col_elements --expelements 8  alignment  ${testdir}/output/database/run${PaddedRunNr}-alignment.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4AlignDB PROPERTIES DEPENDS TestJobsubExampleAnemone2FEI4AlignRun)

    ADD_TEST( TestJobsubExampleAnemone2FEI4AlignOutput sh -c "[ -f ${testdir}/output/database/run${PaddedRunNr}-align-mille.bin -a -f ${testdir}/output/database/run${PaddedRunNr}-pede-steer.txt ] " )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4AlignOutput PROPERTIES DEPENDS TestJobsubExampleAnemone2FEI4AlignRun)


#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 5: FITTER
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    #SET( fit_pass_regex_1 "Total number of reconstructed tracks *[0-9][0-9][0-9][0-9][0-9]+" )

    ADD_TEST( NAME TestJobsubExampleAnemone2FEI4FitterRun 
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubOptions} fitter ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4FitterRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleAnemone2FEI4AlignRun
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleAnemone2FEI4FitterLog sh -c "[ -f ${testdir}/output/logs/fitter-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4FitterLog PROPERTIES DEPENDS TestJobsubExampleAnemone2FEI4FitterRun)

    ADD_TEST( TestJobsubExampleAnemone2FEI4FitterHisto sh -c "[ -f ${testdir}/output/histograms/run${PaddedRunNr}-fitter.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4FitterHisto PROPERTIES DEPENDS TestJobsubExampleAnemone2FEI4FitterRun)

    # we expect to see between 1 and 3 tracks in every event 
    # but tolerate if this is not the case in 40% of the events (empty events are counted)
    ADD_TEST( TestJobsubExampleAnemone2FEI4FitterOutput sh -c "[ -f ${testdir}/output/lcio/run${PaddedRunNr}-track.slcio ] && lcio_check_col_elements --pedantic --expelements 3 --abselementerror 2 --releventerror .40 track ${testdir}/output/lcio/run${PaddedRunNr}-track.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4FitterOutput PROPERTIES DEPENDS TestJobsubExampleAnemone2FEI4FitterRun)

#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 6: StatTest
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    SET( executable "python ${stattestdir}/runtests.py" )

    # all this regular expressions must be matched for the test to pass
    SET( fit_pass_regex_1 "SUCCESS" )
    SET( fit_fail_regex "FAILED" "NOT PASSED" "segmentation violation")

    # run stattest tool on output from previous step and test it against reference file; test are configured in specified config file (*.qa)

    ADD_TEST( TestJobsubExampleAnemone2FEI4StatTestClustering sh -c "PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH ${executable} --cdash -g ${testdir}/output/stattest_report_clus.pdf ${referencedatadir}/StatTestConf_Anemone2FEI4Clustering.qa ${testdir}/output/histograms/run${PaddedRunNr}-clustering.root ${referencedatadir}/run${PaddedRunNr}-clu-histo.root" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4StatTestClustering PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${fit_pass_regex_1}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${fit_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleAnemone2FEI4ClusteringRun
	)


    ADD_TEST( TestJobsubExampleAnemone2FEI4StatTestAlign sh -c "PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH ${executable} --cdash -g ${testdir}/output/stattest_report_align.pdf ${referencedatadir}/StatTestConf_Anemone2FEI4Align.qa ${testdir}/output/histograms/run${PaddedRunNr}-alignmentdaf.root ${referencedatadir}/run${PaddedRunNr}-align-histo.root" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4StatTestAlign PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${fit_pass_regex_1}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${fit_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleAnemone2FEI4AlignRun
	)


    ADD_TEST( TestJobsubExampleAnemone2FEI4StatTestFitter sh -c "PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH ${executable} --cdash  -g${testdir}/output/stattest_report_fitter.pdf ${referencedatadir}/StatTestConf_Anemone2FEI4Fitter.qa ${testdir}/output/histograms/run${PaddedRunNr}-fitter.root ${referencedatadir}/run${PaddedRunNr}-track-histo.root" )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4StatTestFitter PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${fit_pass_regex_1}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${fit_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleAnemone2FEI4FitterRun
	)


#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 7: MemChecks
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
  # STEP 1-5 VARIANTS USED FOR MEMCHECKS ONLY:
    SET( executable python -tt ${jobsubdir}/jobsub.py )
    # options for memcheck runs: reduced run range, plain output for valgrind parsing
    SET( jobsubMemCheckOptions --config=${exampledir}/config.cfg -csv ${exampledir}/runlist.csv -o NativePath=${referencedatadir} -o MaxRecordNumber=2000 -o MemCheckFlag=MemCheck --plain)

  # Converter run with reduced run range
    ADD_TEST( NAME TestJobsubExampleAnemone2FEI4ConverterRunMemCheck
              WORKING_DIRECTORY "${testdir}"
	      COMMAND ${executable} ${jobsubMemCheckOptions} converter ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4ConverterRunMemCheck PROPERTIES
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
    )

    ADD_TEST( NAME TestJobsubExampleAnemone2FEI4ClusteringRunMemCheck
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubMemCheckOptions} clustering ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4ClusteringRunMemCheck PROPERTIES
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	)

    ADD_TEST( NAME TestJobsubExampleAnemone2FEI4HitmakerRunMemCheck
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubMemCheckOptions} hitmaker ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4HitmakerRunMemCheck PROPERTIES
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	)

    ADD_TEST( NAME TestJobsubExampleAnemone2FEI4AlignRunMemCheck
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubMemCheckOptions} align ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4AlignRunMemCheck PROPERTIES
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	)

    ADD_TEST( NAME TestJobsubExampleAnemone2FEI4FitterRunMemCheck
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubMemCheckOptions} fitter ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleAnemone2FEI4FitterRunMemCheck PROPERTIES
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	)
