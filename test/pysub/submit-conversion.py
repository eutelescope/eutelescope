#!/usr/bin/env python
## @package Pysub
# Submit a conversion job
#
# This Python script can be used by both normal and power users for
# conversion job submission. Converting files from the native DAQ format
# to the LCIO one is the very first step in the analysis and has to be done
# for every acquired run.
#
# This script can execute Marlin either local or remotely on the GRID, All
# the default parameters are saved in a configuration file (template/config.py).
# The user can select a different configuration file either using the --config-file
# option or using the enviromental variable SUBMIT_CONFIG.
#
#
# Call ./submit-conversion.py --help for all the options.
#
# @Version $Id: submit-conversion.py,v 1.9 2009-05-08 08:23:27 bulgheroni Exp $
# @Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>


import sys
import string
import re
import os
import popen2
from optparse import OptionParser
from optparse import OptionGroup



blue="\033[1;34m"
black="\033[0m"
red="\033[1;31m"
green="\033[1;32m"

## Submit a conversion job
#    This Python script can be used by both normal and power users for
#    conversion job submission. Converting files from the native DAQ format
#    to the LCIO one is the very first step in the analysis and has to be done
#    for every acquired run.
#
#    This script can execute Marlin either local or remotely on the GRID, All
#    the default parameters are saved in a configuration file (template/config.py).
#    The user can select a different configuration file either using the --config-file
#    option or using the enviromental variable SUBMIT_CONFIG.
#
#
#    Call ./submit-conversion.py --help for all the options.
#
#    @Version $Id: submit-conversion.py,v 1.9 2009-05-08 08:23:27 bulgheroni Exp $
#    @Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
def submit_conversion():
    # parse options and arguments
    usage = "usage: %prog [options] run-list"
    ver   = "%prog $Revision: 1.9 $"

    parser = OptionParser( usage=usage, version=ver )
    executionGroup = OptionGroup( parser, "Execution options",
                                  "Use these options to select where and how the jobs have to executed")


    executionHelp = """
    Select where to execute the job. all-local means: input and output files are stored locally and the job is executed on the local CPU.
    all-grid means: input and output files are taken from the storage element and the job will be submitted to the GRID.
    cpu-local means: input and output files are taken from the GRID SE, but the job will be executed on the local CPU
    """
    executionGroup.add_option( "-x", "--execution",
                               type="choice",
                               action="store",
                               choices=['all-local', 'all-grid', 'cpu-local','only-generate'],
                               dest="execution",
                               help=executionHelp,
                               metavar="EXECUTION")

    parser.set_defaults(execution="all-local")

    executionGroup.add_option( "-l", "--all-local",
                               action="store_const",
                               const="all-local",
                               dest="execution",
                               help="Same as -x all-local")

    executionGroup.add_option( "-r", "--all-grid",
                               action="store_const",
                               const="all-grid",
                               dest="execution",
                               help="Same as -x all-grid" )

    executionGroup.add_option( "-c", "--cpu-local",
                               action="store_const",
                               const="cpu-local",
                               dest="execution",
                               help="Same as -x cpu-local")

    executionGroup.add_option( "-s", "--only-generate",
                               action="store_const",
                               const="only-generate",
                               dest="execution",
                               help="Same as -x only-generate")

    parser.add_option_group( executionGroup )


    ioGroup = OptionGroup( parser, "Input and output files related options",
                           "Use these options to specify whether or not to keep the input and output files"
                           "These values are overwriting the default configurations depending on the execution mode:"
                           "  all-local = input and output files are kept "
                           "  cpu-local = input and output files are removed. ")

    ioGroup.add_option( "--keep-input",
                        action="store_true",
                        dest="force_keep_input",
                        help="Force not to delete the input files after finishing, independently of the execution mode")

    ioGroup.add_option( "--remove-input",
                        action="store_false",
                        dest="force_keep_input",
                        help="Force to delete the input files after finishing, independently of the execution mode")

    ioGroup.add_option( "--keep-output",
                        action="store_true",
                        dest="force_keep_output",
                        help="Force not to delete the output files after finishing, independently of the execution mode")

    ioGroup.add_option( "--remove-output",
                        action="store_false",
                        dest="force_keep_output",
                        help="Force to delete the output files after finishing, independently of the execution mode")


    parser.set_defaults(force_keep_input=True)
    parser.set_defaults(force_keep_output=True)

    parser.add_option_group( ioGroup )

    configurationGroup = OptionGroup( parser, "Configuration option",
                                      "Use these options to specify which configuration file you want to use,"
                                      "and other additional external files" )

    configurationGroup.add_option( "-g", "--gear-file",
                                   action="store",
                                   dest="gear_file",
                                   help="Specify the GEAR file to be used")

    configurationGroup.add_option( "--config-file",
                                   action="store",
                                   dest="config_file",
                                   help="Specify the configuration file to be used")

    parser.add_option_group( configurationGroup )
    option, args = parser.parse_args()

    # end of options and arguments parsing.
    # now checking if there are errors in the arguments

    if len(args) == 0:
        print red, "No input run provided", black
        sys.exit( 0 )

    # verify that all the numbers are integers numbers
    runList = []
    for i in args :
        try:
            runList.append( int ( i ) )

        except ValueError:
            print red, "Invalid run number %(i)s" % { "i": i }
            sys.exit( 1 )


    # this is just for now
    if option.execution == "all-grid" :
        print red, "Full GRID submission not yet implemented", black
        sys.exit( 1 )

    # load the configuration file
    config_file = ""
    if option.config_file == None :
        # this means that the user didn't provide any config file
        # check if the enviromental variable SUBMIT_CONFIG is defined
        try:
            envirConfigFile = os.environ['SUBMIT_CONFIG']
            if len(envirConfigFile ) == 0:
                config_file = "template/config.py"
            else:
                config_file = envirConfigFile
        except KeyError:
            # the enviromental variable is not even defined
            # use the default
            config_file = "template/config.py"
        else:
            config_file = option.config_file

    # verify that the config file does exists!
    if not os.path.exists( config_file ) :
        print red, "The configuration file", config_file, "doesn't exist!", black
        sys.exit ( 2 );
    else:
        execfile( config_file )

    # set the grid paths to the default one
    gridFolderNative  = defaultGRIDFolderNative
    gridFolderLcioRaw = defaultGRIDFolderLcioRaw
    gridFolderConvertJobout = defaultGRIDFolderConvertJoboutput


    # verify if the user provided a specific gear file
    gear_file = ""
    if option.gear_file == None :
        # no, use default then
        gear_file = defaultGEARFile
    else:
        gear_file = option.gear_file

    # in case the execution mode is really doing something, check that the gear file
    # exists before starting
    if option.execution != "only-generate" :
        if not os.path.exists( gear_file ) :
            print red, "GEAR file", gear_file, "doesn't exist!", black
            sys.exit( 3 )


    # let's start now!
    print red, "Conversion job submitter" , black

    # clear these variables
    copiedLCIOFile = []
    copiedTARFile = []


    # set the default io configuration depending on the execution mode
    optionKeepInput  = True
    optionKeepOutput = True

    if option.execution == "all-local" :
        optionKeepInput  = True
        optionKeepOuput  = True

    elif option.execution == "all-grid" :
        # doens't really matter but set all to true
        optionKeepInput  = True
        optionKeepOuput  = True

    elif option.execution == "cpu-local" :
        optionKeepInput  = False
        optionKeepOuput  = False

    # now overwrite the previous setting if the user wants
    if option.force_keep_input == True:
        optionKeepInput = True
    else:
        optionKeepInput = False

    if option.force_keep_output == True:
        optionKeepOutput = True
    else:
        optionKeepOutput = False

    # ask to the users for confirmation about file removal
    if optionKeepInput == False and option.execution != "only-generate" and option.execution != "all-grid" :
        goodAnswer = False
        while goodAnswer == False:
            print green, "This script is going to delete the input file(s) when finished.\n" ,\
                  " Are you sure to continue? [y/n]", black
            answer = raw_input( "--> ").lower()
            if answer != "y" and answer != "yes" and answer != "n" and answer != "no":
                print red, "Invalid answer, please type y or n", black
                answer = raw_input( "--> " ).lower()
            elif answer == "y" or answer == "yes":
                goodAnswer = True
            elif answer == "n" or answer == "no":
                goodAnswer = True
                sys.exit("Aborted by the user")

    if optionKeepOutput == False and option.execution != "only-generate" and option.execution != "all-grid":
        goodAnswer = False
        while goodAnswer == False:
            print green, "This script is going to delete the output file(s) when finished.\n" ,\
                  " Are you sure to continue? [y/n]", black
            answer = raw_input( "--> ").lower()
            if answer != "y" and answer != "yes" and answer != "n" and answer != "no":
                print red, "Invalid answer, please type y or n", black
                answer = raw_input( "--> " ).lower()
            elif answer == "y" or answer == "yes":
                goodAnswer = True
            elif answer == "n" or answer == "no":
                goodAnswer = True
                sys.exit("Aborted by the user")

    for run in runList[:]:
        returnvalue = 0
        runString = "%(#)06d" % { "#": run }

        print blue, "Generating the steering file... (universal-%(run)s.xml)" % { "run": runString }, black

        # open the template steering file for reading
        templateSteeringFile = open( "./template/universal-tmp.xml", "r")

        # read the whole content of the template and put it into a string
        templateSteeringString = templateSteeringFile.read()

        actualSteeringString = templateSteeringString.replace("@RunNumber@", runString )
        actualSteeringString = actualSteeringString.replace("@GearFile@", gear_file )

        actualSteeringFile = open( "universal-%(run)s.xml" % { "run": runString }, "w" )
        actualSteeringFile.write( actualSteeringString )

        actualSteeringFile.close()

        templateSteeringFile.close()

        # if only-generate... then we are done already...
        # continue to the next run
        if option.execution == "only-generate" :
            continue

        # this is the case for all-grid submission
        if option.execution == "all-grid" :
            # in principle we should never get to this point for the time being
            # continue anyway
            continue

        # this is the case for all-local and cpu-local
        # since they are quite similar, I'll put the two in the same state
        # of the state machine
        if option.execution == "all-local" or option.execution == "cpu-local" :


            # if cpu-local, we need the input file from the GRID
            if option.execution == "cpu-local":
                print blue, "Getting the native raw file from the GRID...", black
                command = "lcg-cp -v lfn:%(nativeFolder)s/run%(run)s.raw file:$PWD/native/run%(run)s.raw" % \
                      { "nativeFolder": gridFolderNative , "run": runString }
                returnvalue = os.system( command )
                if returnvalue != 0:
                    print red, "Problem getting the input file! (errno %(returnvalue)s)" % {"returnvalue":returnvalue }, black
                    continue



            # run Marlin in both all local and cpu local configuration
            print blue, "Running Marlin...",black

            # to properly get the return code I need to execute Marlin into a separete pipe
            logfile = open( "universal-%(run)s.log " % { "run": runString }, "w" )
            marlin = popen2.Popen4( "Marlin universal-%(run)s.xml" % { "run": runString } )
            while marlin.poll() == -1:
                line = marlin.fromchild.readline()
                print line.strip()
                logfile.write( line )

            logfile.close()
            returnvalue = marlin.poll()

            if returnvalue != 0:
                print red, "Problem executing Marlin! (errno %(returnvalue)s)" % {"returnvalue":returnvalue } , black
                continue

            # Copy to the GRID only in cpu-local mode
            if option.execution == "cpu-local" :
                print blue, "Copying register the output file to the GRID...", black
                command = "lcg-cr -v -l lfn:%(lciorawFolder)s/run%(run)s.slcio file:$PWD/lcio-raw/run%(run)s.slcio" % \
                          { "lciorawFolder": gridFolderLcioRaw , "run": runString }
                returnvalue = os.system( command )
                if returnvalue != 0:
                    print red, "Problem copying the lcio-raw file to the GRID! (errno %(returnvalue)s)" % {"returnvalue":returnvalue }, black
                else:
                    copiedLCIOFile.append( "%(lciorawFolder)s/run%(run)s.slcio" % { "lciorawFolder": gridFolderLcioRaw , "run": runString } )

                # Prepare the tarbal in any local configurations (only in cpu-local)
                print blue, "Preparing the joboutput tarbal...", black

                # create a temporary folder
                command = "mkdir /tmp/universal-%(run)s" % { "run": runString }
                returnvalue = os.system( command )
                if returnvalue != 1 and returnvalue != 0:
                    print red, "Problem creating the temporary folder! (errno %(returnvalue)s)" % {"returnvalue":returnvalue }, black
                else:

                    # copy there all the tarbal stuff
                    command = "cp universal-%(run)s* %(gear)s /tmp/universal-%(run)s" % { "gear": gear_file , "run" : runString }
                    returnvalue = os.system( command )
                    if returnvalue != 0:
                        print red, "Problem copying the joboutput files in the temporary folder! (errno %(returnvalue)s)" % {"returnvalue":returnvalue }, black
                    else:

                        # tar gzipped the temp folder
                        command = "tar czvf universal-%(run)s.tar.gz /tmp/universal-%(run)s"  % { "run" : runString }
                        returnvalue = os.system( command )
                        if returnvalue != 0:
                            print red, "Problem tar gzipping the temporary folder! (errno %(returnvalue)s)" % {"returnvalue":returnvalue }, black
                        else:

                            # copy the tarbal to the GRID
                            if optionCPULocal == 1:
                                print blue, "Copying the joboutput tarbal to the GRID...", black
                                command = "lcg-cr -v -l lfn:%(joboutputFolder)s/universal-%(run)s.tar.gz file:$PWD/universal-%(run)s.tar.gz" % \
                                          { "joboutputFolder": gridFolderConvertJobout , "run": runString }
                                returnvalue = os.system( command )
                                if returnvalue != 0:
                                    print red, "Problem copying to the GRID! (errno %(returnvalue)s)" % {"returnvalue":returnvalue }, black
                                else:
                                    # adding this file to the list of copied files
                                    copiedTARFile.append( "%(joboutputFolder)s/universal-%(run)s.tar.gz" \
                                                          %  { "joboutputFolder": gridFolderConvertJobout , "run": runString } )
                                    # remove all the xml and log files from the local machine including the tarbal just copied
                                    command = "rm -vrf  universal-%(run)s.tar.gz" % { "run" : runString }
                                    returnvalue = os.system( command )
                                    if returnvalue != 0:
                                        print red, "Problem removing the tarbal", black

            # clean up the execution enviroment (both all-local and cpu-local)
            print blue, "Cleaning up enviroment", black
            command = "rm -vrf /tmp/universal-%(run)s universal-%(run)s*" % { "run": runString }
            returnvalue = os.system( command )
            if returnvalue != 0:
                print red, "Problem removing temporary files! (errno %(returnvalue)s)" % {"returnvalue":returnvalue }, black

            # copying the tarbal in the log folder only in all-local mode
            if option.execution == "all-local":
                print blue, "Copying the tarbal in the local log folder", black
                command = "mv  universal-%(run)s.tar.gz logs/." % { "run" : runString }
                returnvalue = os.system( command )
                if returnvalue != 0:
                    print red, "Problem moving the tarbal to the local log folder (errno %(returnvalue)s)" % {"returnvalue":returnvalue }, black

            if optionKeepOutput == 0 and optionKeepInput == 0:
                command = "rm -v lcio-raw/run%(run)s.slcio native/run%(run)s.raw" % { "run": runString }
                returnvalue = os.system( command )
                if returnvalue != 0:
                    print red, "Problem removing the input and output files! (errno %(returnvalue)s)" % {"returnvalue":returnvalue }, black
                    continue

            elif optionKeepOutput == 1 and optionKeepInput == 0:
                command = "rm -v native/run%(run)s.raw" % { "run": runString }
                returnvalue = os.system( command )
                if returnvalue != 0:
                    print red, "Problem removing the input file! (errno %(returnvalue)s)" % {"returnvalue":returnvalue }, black
                    continue

            elif optionKeepOutput == 0 and optionKeepInput == 1:
                command = "rm -v lcio-raw/run%(run)s.slcio" % { "run": runString }
                returnvalue = os.system( command )
                if returnvalue != 0:
                    print red, "Problem removing the output file! (errno %(returnvalue)s)" % {"returnvalue":returnvalue }, black
                    continue

    if len(copiedLCIOFile) != 0 or len(copiedTARFile) != 0:
        print green, "Summary"
        print " The following lcio-raw file(s) were copied to the GRID:", blue
        for file in copiedLCIOFile[:] :
            print "\t ", file

        print "\n", green , "The following tarbal(s) were copied to the GRID:", blue
        for file in copiedTARFile[:] :
            print "\t ", file

    print black



if __name__ = "__main__" :
    submit_conversion()
