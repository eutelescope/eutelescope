#
# This file defines a number of data-driven tests based on the example configuration
# in this directory. If you have access the the files and paths defined below you
# can run the tests by running 'make test' in the build directory in the EUTelescope root
#

# ======================================================================
# ======================================================================
# TestJobsubExamplex0: based on config in jobsub/examples/x0
# ======================================================================
# ======================================================================

    SET( testdir "${PROJECT_BINARY_DIR}/Testing/test_x0" )
    SET( jobsubdir "$ENV{EUTELESCOPE}/jobsub" )
    SET( exampledir "${jobsubdir}/examples/X0" )

    # the run number
    SET( RunNr "5575" )
    # run number padded with leading zeros
    execute_process(COMMAND sh -c "printf %06d ${RunNr}" OUTPUT_VARIABLE PaddedRunNr)
 
   # only needed in the last step to test the results of EUTel against a set of reference files:
    SET( stattestdir "$ENV{EUTELESCOPE}/test/stattest/bin" )
    SET( datadir "/afs/desy.de/group/telescopes/EutelTestData/TestExampleX0" )
    SET( referencedatadir "/afs/desy.de/group/telescopes/EutelTestData/TestExampleX0" )

    SET( executable python -tt ${jobsubdir}/jobsub.py )
    # options: use config, use csv, change native path to central AFS location, reduce number of events to 200k
    SET( jobsubOptions --config=${exampledir}/config.cfg -csv ${exampledir}/runs.csv -o NativePath=${datadir} -o MaxRecordNumber=100000)


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
	ADD_TEST( TestJobsubExampleX0Cleanup sh -c "[ -d ${testdir} ] && rm -rf ${testdir} || echo 'no cleanup needed.'" )
	ADD_TEST( TestJobsubExampleX0Setup sh -c "mkdir -p ${testdir}/output/histograms  && mkdir -p ${testdir}/output/database && mkdir -p ${testdir}/output/logs && mkdir -p ${testdir}/output/results" )
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 1: CONVERTER
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

    ADD_TEST( NAME TestJobsubExampleX0ConverterRun 
              WORKING_DIRECTORY "${testdir}"
	      COMMAND ${executable} ${jobsubOptions} converter ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0ConverterRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleX0Setup
	# converter step sometimes takes a bit longer (in s)
	TIMEOUT 2500
    )


    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleX0ConverterLog sh -c "[ -f ${testdir}/output/logs/converter-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0ConverterLog PROPERTIES DEPENDS TestJobsubExampleX0ConverterRun)

    ADD_TEST( TestJobsubExampleX0ConverterHisto sh -c "[ -f ${testdir}/output/histograms/run${PaddedRunNr}-converter.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0ConverterHisto PROPERTIES DEPENDS TestJobsubExampleX0ConverterRun)

    ADD_TEST( TestJobsubExampleX0ConverterHotpix sh -c "[ -f ${testdir}/output/database/run${PaddedRunNr}-hotpixel.slcio ] && lcio_check_col_elements --expelements 6  hotpixel_m26  ${testdir}/output/database/run${PaddedRunNr}-hotpixel.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0ConverterHotpix PROPERTIES DEPENDS TestJobsubExampleX0ConverterRun)

    ADD_TEST( TestJobsubExampleX0ConverterOutput sh -c "[ -f ${testdir}/output/lcio/run${PaddedRunNr}-converter.slcio ] && lcio_check_col_elements --expelements 6 zsdata_m26 ${testdir}/output/lcio/run${PaddedRunNr}-converter.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0ConverterOutput PROPERTIES DEPENDS TestJobsubExampleX0ConverterRun)




#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 2: CLUSTERING
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    ADD_TEST( NAME TestJobsubExampleX0ClusteringRun 
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubOptions} clustering ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0ClusteringRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleX0ConverterOutput
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleX0ClusteringLog sh -c "[ -f ${testdir}/output/logs/clustering-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0ClusteringLog PROPERTIES DEPENDS TestJobsubExampleX0ClusteringRun)

    ADD_TEST( TestJobsubExampleX0ClusteringHisto sh -c "[ -f ${testdir}/output/histograms/run${PaddedRunNr}-clustering.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0ClusteringHisto PROPERTIES DEPENDS TestJobsubExampleX0ClusteringRun)

    # we expect an average of 24.4 clusters per event
    ADD_TEST( TestJobsubExampleX0ClusteringOutput sh -c "[ -f ${testdir}/output/lcio/run${PaddedRunNr}-clustering.slcio ] && lcio_check_col_elements --average --expelements 24 cluster_m26 ${testdir}/output/lcio/run${PaddedRunNr}-clustering.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0ClusteringOutput PROPERTIES DEPENDS TestJobsubExampleX0ClusteringRun)


#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 2B: CLUSTER FILTERING
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    ADD_TEST( NAME TestJobsubExampleX0FilterRun 
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubOptions} filter ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0FilterRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleX0ClusteringOutput
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleX0FilterLog sh -c "[ -f ${testdir}/output/logs/filter-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0FilterLog PROPERTIES DEPENDS TestJobsubExampleX0FilterRun)

    # we now expect an average of 22.9 clusters per event (after filtering)
    ADD_TEST( TestJobsubExampleX0FilterOutput sh -c "[ -f ${testdir}/output/lcio/run${PaddedRunNr}-clustering-filtered.slcio ] && lcio_check_col_elements --average --expelements 23 filtered_cluster_m26 ${testdir}/output/lcio/run${PaddedRunNr}-clustering-filtered.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0FilterOutput PROPERTIES DEPENDS TestJobsubExampleX0FilterRun)



#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 3: HITMAKER
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    ADD_TEST( NAME TestJobsubExampleX0HitmakerRun 
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubOptions} hitmaker ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0HitmakerRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleX0ClusteringOutput
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleX0HitmakerLog sh -c "[ -f ${testdir}/output/logs/hitmaker-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0HitmakerLog PROPERTIES DEPENDS TestJobsubExampleX0HitmakerRun)

    ADD_TEST( TestJobsubExampleX0HitmakerHisto sh -c "[ -f ${testdir}/output/histograms/run${PaddedRunNr}-hitmaker.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0HitmakerHisto PROPERTIES DEPENDS TestJobsubExampleX0HitmakerRun)

    ADD_TEST( TestJobsubExampleX0HitmakerPrealign sh -c "[ -f ${testdir}/output/database/run${PaddedRunNr}-prealignment.slcio ] && lcio_check_col_elements --expelements 6  alignment  ${testdir}/output/database/run${PaddedRunNr}-prealignment.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0HitmakerPrealign PROPERTIES DEPENDS TestJobsubExampleX0HitmakerRun)

    # we expect an average hit number of 23 for run 97 (wide geometry) using the example configuration
    ADD_TEST( TestJobsubExampleX0HitmakerOutput sh -c "[ -f ${testdir}/output/lcio/run${PaddedRunNr}-hitmaker.slcio ] && lcio_check_col_elements -a --expelements 23 hit ${testdir}/output/lcio/run${PaddedRunNr}-hitmaker.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0HitmakerOutput PROPERTIES DEPENDS TestJobsubExampleX0HitmakerRun)


#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 4: ALIGNMENT
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    # all this regular expressions must be matched for the test to pass
    SET( align_pass_regex_1 "Initialising Mille" )
    SET( align_pass_regex_2 "Pede successfully finished" )

    ADD_TEST( NAME TestJobsubExampleX0AlignRun
    	      WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubOptions} aligndaf ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0AlignRun PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${align_pass_regex_1}.*${marlin_pass_regex_1}.*${align_pass_regex_2}.*${jobsub_pass_regex_2}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleX0HitmakerOutput
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleX0AlignLog sh -c "[ -f ${testdir}/output/logs/align-${PaddedRunNr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0AlignLog PROPERTIES DEPENDS TestJobsubExampleX0AlignRun)

    ADD_TEST( TestJobsubExampleX0AlignHisto sh -c "[ -f ${testdir}/output/histograms/run${PaddedRunNr}-alignment.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0AlignHisto PROPERTIES DEPENDS TestJobsubExampleX0AlignRun)

    ADD_TEST( TestJobsubExampleX0AlignDB sh -c "[ -f ${testdir}/output/database/run${PaddedRunNr}-alignment.slcio ] && lcio_check_col_elements --expelements 6  alignment  ${testdir}/output/database/run${PaddedRunNr}-alignment.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0AlignDB PROPERTIES DEPENDS TestJobsubExampleX0AlignRun)

    ADD_TEST( TestJobsubExampleX0AlignOutput sh -c "[ -f ${testdir}/output/database/run${PaddedRunNr}-align-mille.bin -a -f ${testdir}/output/database/run${PaddedRunNr}-pede-steer.txt ] " )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0AlignOutput PROPERTIES DEPENDS TestJobsubExampleX0AlignRun)



#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  step 5: fitter
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    #set( fit_pass_regex_1 "total number of reconstructed tracks *[0-9][0-9][0-9][0-9][0-9]+" )

    ADD_TEST( NAME TestJobsubExampleX0FitterRun
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsuboptions} fitter ${Runnr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0FitterRun PROPERTIES 
        # test will pass if all of the following expressions are matched
	PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        # test will fail if any of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	# test DEPENDS on earlier steps
	DEPENDS TestJobsubExampleX0AlignRun
	# Fitter step sometimes takes a bit longer (in s)
	timeout 2500
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleX0Fitterlog sh -c "[ -f ${testdir}/output/logs/fitter-${paddedRunnr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0Fitterlog PROPERTIES DEPENDS TestJobsubExampleX0FitterRun)

    ADD_TEST( TestJobsubExampleX0Fitterhisto sh -c "[ -f ${testdir}/output/histograms/Run${paddedRunnr}-Fitter.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0Fitterhisto PROPERTIES DEPENDS TestJobsubExampleX0FitterRun)

    # we expect to see between 1 and 3 tracks in every event 
    # but tolerate if this is not the case in 40% of the events (empty events are counted)
    ADD_TEST( TestJobsubExampleX0Fitteroutput sh -c "[ -f ${testdir}/output/lcio/Run${paddedRunnr}-track.slcio ] && lcio_check_col_elements --pedantic --expelements 2 --abselementerror 1 --releventerror .40 track0 ${testdir}/output/lcio/Run${paddedRunnr}-track.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0Fitteroutput PROPERTIES DEPENDS TestJobsubExampleX0FitterRun)
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  step 6: x0 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    #set( fit_pass_regex_1 "total number of reconstructed tracks *[0-9][0-9][0-9][0-9][0-9]+" )

    ADD_TEST( NAME TestJobsubExampleX0X0Run 
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsuboptions} x0 ${Runnr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0X0Run PROPERTIES
        # test will pass if all of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        # test will fail if any of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	# test DEPENDS on earlier steps
	DEPENDS TestJobsubExampleX0FitterRun
	# X0 step sometimes takes a bit longer (in s)
	timeout 2500
	)
    # now check if the expected output files exist and look ok
    ADD_TEST( TestJobsubExampleX0X0log sh -c "[ -f ${testdir}/output/logs/X0-${paddedRunnr}.zip ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0X0log PROPERTIES DEPENDS TestJobsubExampleX0X0Run)

    ADD_TEST( TestJobsubExampleX0X0histo sh -c "[ -f ${testdir}/output/histograms/Run${paddedRunnr}-X0.root ]" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0X0histo PROPERTIES DEPENDS TestJobsubExampleX0X0Run)

    # we expect to see between 1 and 3 tracks in every event 
    # but tolerate if this is not the case in 40% of the events (empty events are counted)
    ADD_TEST( TestJobsubExampleX0X0output sh -c "[ -f ${testdir}/output/lcio/Run${paddedRunnr}-track.slcio ] && lcio_check_col_elements --pedantic --expelements 2 --abselementerror 1 --releventerror .40 track0 ${testdir}/output/lcio/Run${paddedRunnr}-track.slcio" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0X0output PROPERTIES DEPENDS TestJobsubExampleX0X0Run)

#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 7: StatTest
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
    SET( executable "python ${stattestdir}/runtests.py" )

    # all this regular expressions must be matched for the test to pass
    SET( fit_pass_regex_1 "SUCCESS" )
    SET( fit_fail_regex "FAILED" "NOT PASSED" "Error" "segmentation violation")

    # run stattest tool on output from previous step and test it against reference file; test are configured in specified config file (*.qa)

    ADD_TEST( TestJobsubExampleX0StatTestClustering sh -c "PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH ${executable} --cdash -g ${testdir}/output/stattest_report_clustering.pdf ${referencedatadir}/StatTestConf_X0Clustering.qa ${testdir}/output/histograms/run${PaddedRunNr}-clustering.root ${referencedatadir}/run${PaddedRunNr}-clustering.root" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0StatTestClustering PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${fit_pass_regex_1}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${fit_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleX0ClusteringRun
	)


    ADD_TEST( TestJobsubExampleX0StatTestAlign sh -c "PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH ${executable} --cdash -g ${testdir}/output/stattest_report_align.pdf ${referencedatadir}/StatTestConf_X0Align.qa ${testdir}/output/histograms/run${PaddedRunNr}-alignment.root ${referencedatadir}/run${PaddedRunNr}-alignment.root" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0StatTestAlign PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${fit_pass_regex_1}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${fit_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleX0AlignRun
	)


    ADD_TEST( TestJobsubExampleX0StatTestFitter sh -c "PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH ${executable} --cdash  -g${testdir}/output/stattest_report_fitter.pdf ${referencedatadir}/StatTestConf_X0Fitter.qa ${testdir}/output/histograms/run${PaddedRunNr}-fitter.root ${referencedatadir}/run${PaddedRunNr}-fitter.root" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0StatTestFitter PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${fit_pass_regex_1}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${fit_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleX0FitterRun
	)


    ADD_TEST( TestJobsubExampleX0StatTestx0 sh -c "PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH ${executable} --cdash  -g${testdir}/output/stattest_report_x0.pdf ${referencedatadir}/StatTestConf_X0x0.qa ${testdir}/output/histograms/run${PaddedRunNr}-x0.root ${referencedatadir}/run${PaddedRunNr}-x0.root" )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0StatTestx0 PROPERTIES
        # test will pass if ALL of the following expressions are matched
        PASS_REGULAR_EXPRESSION "${fit_pass_regex_1}"
        # test will fail if ANY of the following expressions is matched 
        FAIL_REGULAR_EXPRESSION "${fit_fail_regex}"
	# test depends on earlier steps
	DEPENDS TestJobsubExampleX0x0Run
	)

#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  STEP 7: MemChecks
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
  # STEP 1-6 VARIANTS USED FOR MEMCHECKS ONLY:
    SET( executable python -tt ${jobsubdir}/jobsub.py )
    # options for memcheck runs: reduced run range, plain output for valgrind parsing
    SET( jobsubMemCheckOptions --config=${exampledir}/config.cfg -csv ${exampledir}/runs.csv -o NativePath=${datadir} -o MaxRecordNumber=2000 --plain)

  # Converter run with reduced run range
    ADD_TEST( NAME TestJobsubExampleX0ConverterRunMemCheck
              WORKING_DIRECTORY "${testdir}"
	      COMMAND ${executable} ${jobsubMemCheckOptions} converter ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0ConverterRunMemCheck PROPERTIES
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
    )

    ADD_TEST( NAME TestJobsubExampleX0ClusteringRunMemCheck
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubMemCheckOptions} clustering ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0ClusteringRunMemCheck PROPERTIES
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	)

    ADD_TEST( NAME TestJobsubExampleX0FilterRunMemCheck
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubMemCheckOptions} filter ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0FilterRunMemCheck PROPERTIES
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	)

    ADD_TEST( NAME TestJobsubExampleX0HitmakerRunMemCheck
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubMemCheckOptions} hitmaker ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0HitmakerRunMemCheck PROPERTIES
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	)

    ADD_TEST( NAME TestJobsubExampleX0AlignRunMemCheck
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubMemCheckOptions} aligndaf ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0AlignRunMemCheck PROPERTIES
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	)

    ADD_TEST( NAME TestJobsubExampleX0FitterRunMemCheck
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubMemCheckOptions} fitter ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0FitterRunMemCheck PROPERTIES
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	)
    ADD_TEST( NAME TestJobsubExampleX0x0RunMemCheck
              WORKING_DIRECTORY ${testdir} 
	      COMMAND ${executable} ${jobsubMemCheckOptions} x0 ${RunNr} )
    SET_TESTS_PROPERTIES (TestJobsubExampleX0x0RunMemCheck PROPERTIES
        PASS_REGULAR_EXPRESSION "${jobsub_pass_regex_1}.*${marlin_pass_regex_1}.*${jobsub_pass_regex_2}"
        FAIL_REGULAR_EXPRESSION "${generic_fail_regex}"
	)
