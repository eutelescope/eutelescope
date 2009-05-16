import os
import sys
import shutil
import sha
import glob
import tarfile
import popen2
import ConfigParser
import logging
import logging.handlers
import math
import datetime
import time
from submitbase import SubmitBase
from error import *

## Submit cluster filtering jobs
#
# This calss is responsible to submit jobs for pedestal and noise calculation.
# It is inheriting from SubmitBase and it is called by the submit-filter.py script
#
#
# @version $Id: submitfilter.py,v 1.11 2009-05-16 13:05:44 bulgheroni Exp $
# @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
#
class SubmitFilter( SubmitBase ):

    ## Version
    # This number is used when printing out the version of the package.
    #
    # Static member.
    #
    cvsVersion = "$Revision: 1.11 $"

    ## Name
    # This is the namer of the class. It is used in flagging all the log entries
    # and preparing all the file
    #
    # Static member.
    #
    name = "filter"

    ## General configure
    #
    # This method is called by the constructor itself and performs all the
    # basic configurations from the configuration file.
    # In particular is calling the configureLogger method to start the logging
    # system in its full glory
    #
    def configure( self ):

        # first of all we need to see if we have a user defined configuration
        # file different from the template/config.cfg
        #
        # The configuration file can be passed either as a command line option
        # or via an enviromental variable
        #
        SubmitBase.configure( self )

        # now we can properly set the logger.
        SubmitBase.configureLogger( self, self.name )

        # print the welcome message!
        self._logger.log( 15, "**********************************************************" )
        self._logger.log( 15, "Started submit-%(name)s" % { "name": self.name }  )
        self._logger.log( 15, "**********************************************************" )

        # now print the configuration to the log
        SubmitBase.logConfigurationFile( self )

        # now log the run list
        SubmitBase.logRunList( self )

        # now check the I/O options
        # default values depend on the execution mode
        self._keepInput  = True
        self._keepOutput = True

        if self._option.execution == "all-grid" :
            # it doesn't really matter since all the inputs and outputs will be on the GRID
            # and not locally. Leave both to yes
            pass
        elif self._option.execution == "all-local" :
            # usually you want to keep both!
            pass
        elif self._option.execution == "cpu-local" :
            # remove everything
            self._keepInput  = False
            self._keepOutput = False
        elif self._option.execution == "only-generate" :
            # doesn't matter but leave both on
            pass

        # do some checks on the command line options.
        if  self._option.force_keep_input and self._option.force_remove_input :
            self._logger.critical( "Keep and remove input file options are mutually exclusive" )
            self._optionParser.error( "Keep and remove input file options are mutually exclusive" )

        if  self._option.force_keep_output and self._option.force_remove_output :
            self._logger.critical( "Keep and remove output file options are mutually exclusive" )
            self._optionParser.error( "Keep and remove output file options are mutually exclusive" )

        if  self._option.verify_output and self._option.execution != "cpu-local" :
            self._logger.warning( "Selected option verify output in an execution mode different from cpu-local is meaningless" )
            self._option.verift_output = False

        # now check if the user overwrites this setting in the options
        if  self._option.force_keep_input and not self._keepInput:
            self._keepInput = True
            self._logger.info( "User forces to keep the input files" )

        if  self._option.force_remove_input and self._keepInput:
            self._keepInput = False
            self._logger.info( "User forces to remove the input files" )

        if  self._option.force_keep_output and not self._keepOutput:
            self._keepOutput = True
            self._logger.info( "User forces to keep the output files" )

        if  self._option.force_remove_output and self._keepOutput:
            self._keepOutput = False
            self._logger.info( "User forces to remove the output files" )

        # now in case the user wants to interact and the situation is dangerous
        # i.e. removing files, ask confirmation
        try :
            interactive = self._configParser.getboolean( "General", "Interactive" )
        except ConfigParser.NoOptionError :
            interactive = True

        if self._keepInput == False and interactive :
            self._logger.warning("This script is going to delete the input file(s) when finished." )
            self._logger.warning("Are you sure to continue? [y/n]")
            message = "--> "
            if self.askYesNo(prompt = message ) == True:
                self._logger.info("User decide to continue removing input files")
            else:
                self._logger.info("Aborted by user")
                sys.exit( 4 )

        if self._keepOutput == False and interactive :
            self._logger.warning("This script is going to delete the output file(s) when finished." )
            self._logger.warning("Are you sure to continue? [y/n]")
            message = "--> "
            if self.askYesNo(prompt = message ) == True:
                self._logger.info("User decide to continue removing output files")
            else:
                self._logger.info("Aborted by user")
                sys.exit( 4 )


    ## Execute method
    #
    # This is the real method, responsible for the job submission
    # Actually this is just a sort of big switch calling the real submitter
    # depending of the execution mode
    #
    def execute( self ) :

        # first of all verify that the pedestal run was provided, otherwise stop immediately
        self._pedeString = ""
        if self._option.pedestal == None:
            message = "Pedestal run not provided. Please use option -p to specify the pedestal run"
            self._logger.critical( message )
            raise StopExecutionError( message )

        else :
            self._pedeString = "%(pede)06d" % { "pede":  self._option.pedestal }

        # check the pedestal file
        try :
            self.checkPedestalFile()

        except MissingPedestalFileError, error:
            raise StopExecutionError("Missing pedestal file (locally)")

        except MissingPedestalFileOnGRIDError, error:
            raise StopExecutionError("Missing pedestal file (on GRID)")


        # now the execution is different in case we go for merging or not
        if self._option.merge:
            self.executeMerge( )
        else:
            self.executeSingleRun()


    ## Execute in case of merge
    def executeMerge( self ) :

        self._runStringList = []
        for i in self._args:
            try:
                self._runStringList.append( "%(run)06d" % { "run": int( i ) } )
                entry = i, "Unknown", "Unknown", "Unknown", "Unknown", "Unknown"

                # the current entry to the ntuple
                self._summaryNTuple.append( entry )
            except ValueError:
                message = "Invalid run number %(i)s. Skipping to the next." % { "i": i }
                self._logger.warning( message )

        entry = self._option.output, "N/A",  "Unknown", "Unknown", "Unknown", "Unknown"
        self._summaryNTuple.append( entry )

        self._logger.info( "Starting the merging run %(output)s" % { "output": self._option.output } )

        # now do something different depending on the execution option
        try:
            if self._option.execution == "all-grid" :
                self.executeMergeAllGRID( )

            elif self._option.execution == "all-local" :
                self.executeMergeAllLocal( )

            elif self._option.execution == "cpu-local":
                self.executeMergeCPULocal( )

            elif self._option.execution == "only-generate":
                self.executeMergeOnlyGenerate( )

        except MissingSteeringTemplateError, error:
            message = "Steering template %(file)s unavailble. Quitting!" % { "file": error._filename }
            self._logger.critical( message )
            raise StopExecutionError( message )


    ## Execute in case of single runs
    #
    # this is the real part
    # convert all the argumets into a run string compatible with the file name convention
    def executeSingleRun( self ):
        self._runList = []
        for i in self._args:
            try:
                self._runList.append( int( i  ) )

                # if it is a good run number than we can prepare an entry for the summary ntuple
                # the ntuple contains 6 variables:
                # runNumber ; inputFileStatus ; marlinStatus ; outputFileStatus ; histoFileStatus ; joboutputFileStatus
                entry = i , "Unknown", "Unknown", "Unknown", "Unknown", "Unknown"

                # the current entry to the ntuple
                self._summaryNTuple.append( entry )

                # do the same also for the GRID NTuple
                entry = i , "Unknown"
                self._gridJobNTuple.append( entry )

            except ValueError:
                message = "Invalid run number %(i)s" % { "i": i }
                self._logger.critical( message )
                raise StopExecutionError( message )

        for index, run in enumerate( self._runList ) :
            # prepare a string such 123456
            runString = "%(run)06d" % { "run" : run }

            message = "Now processing run %(run)s [ %(i)d / %(n)d ] " % {
                "run" : runString, "i": index + 1, "n": len(self._runList ) }
            self._logger.info( message )

            try:

                # now do something different depending on the execution option
                if self._option.execution == "all-grid" :
                    self.executeSingleRunAllGRID( index , runString )

                elif self._option.execution == "all-local" :
                    self.executeSingleRunAllLocal( index , runString )

                elif self._option.execution == "cpu-local":
                    self.executeSingleRunCPULocal( index , runString )

                elif self._option.execution == "only-generate":
                    self.executeSingleRunOnlyGenerate( index , runString )

            except MissingInputFileError, error:
                message = "Missing input file %(file)s" % { "file": error._filename }
                self._logger.error( message )
                self._logger.error("Skipping to the next run ")
                run, input, marlin, output, histo, tarball = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, "Missing", "Skipped", output, histo, tarball

            except MissingInputFileOnGRIDError, error:
                message = "Missing input file %(file)s" % { "file": error._filename }
                self._logger.error( message )
                self._logger.error("Skipping to the next run ")
                run, input, marlin, output, histo, tarball = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, "Missing", "Skipped", output, histo, tarball

            except OutputFileAlreadyOnGRIDError, error:
                message = "Output file %(file)s already on GRID" % { "file": error._filename }
                self._logger.error( message )
                self._logger.error("Skipping to the next run ")
                run, input, marlin, output, histo, tarball = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, input, "Skipped", "GRID", histo, tarball

            except HistogramFileAlreadyOnGRIDError, error:
                message = "Histogram file %(file)s already on GRID" % { "file": error._filename }
                self._logger.error( message )
                self._logger.error("Skipping to the next run ")
                run, input, marlin, output, histo, tarball = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, input, "Skipped", output, "GRID", tarball

            except JoboutputFileAlreadyOnGRIDError, error:
                message = "Joboutput file %(file)s already on GRID" % { "file": error._filename }
                self._logger.error( message )
                self._logger.error("Skipping to the next run ")
                run, input, marlin, output, histo, tarball = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, input, "Skipped", output, histo, "GRID"

            except MissingSteeringTemplateError, error:
                message = "Steering template %(file)s unavailble. Quitting!" % { "file": error._filename }
                self._logger.critical( message )
                raise StopExecutionError( message )

            except MissingGEARFileError, error:
                message = "Missing GEAR file %(file)s. Quitting! "  % { "file": error._filename }
                self._logger.critical( message )
                raise StopExecutionError( message )

            except MarlinError, error:
                message = "Error with Marlin execution (%(msg)s - errno = %(errno)s )" \
                    % { "msg": error._message, "errno": error._errno }
                self._logger.error( message )
                self._logger.error("Skipping to the next run ")
                run, input, marlin, output, histo, tarball = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, input, "Failed", "Missing", "Missing", "Missing"

            except MissingOutputFileError, error:
                message = "The output file (%(file)s) was not properly generated, possible failure" % { "file": error._filename }
                self._logger.error( message )
                run, b, c, d, e, f = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, b, c, "Missing", d, f

            except MissingHistogramFileError, error:
                message = "The output file (%(file)s) was not properly generated, possible failure" % { "file": error._filename }
                self._logger.error( message )
                run, b, c, d, e, f = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, b, c, d, "Missing", f

            except MissingJoboutputFileError, error:
                message = "The joboutput tarball (%(file)s) is missing, possible failure" % { "file": error._filename }
                self._logger.error( message )
                run, b, c, d, e, f = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, b, c, d, e, "Missing"

            except GRID_LCG_CPError, error:
                message = "Problem copying the input file (%(file)s)" % { "file": error._filename }
                self._logger.error( message )
                self._logger.error( "Skipping to the next run " )
                run, input, marlin, output, histo, tarball = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, "Missing", marlin, output, histo, tarball

    ## Generate only submitter with merging
    def executeMergeOnlyGenerate( self ) :

        # just need to generate the steerign file
        self.generateSteeringFileMerge()

    ## Generate only submitter
    #
    # This methods is responsibile of dry-run with only steering file
    # generation
    #
    def executeSingleRunOnlyGenerate( self, index , runString ):

        # just need to generate the steering file
        self.generateSteeringFileSingleRun( runString )


    ## All Local submitter
    #
    # This methods represents the sequence of commands that should be done while
    # submitting jobs on the local computer.
    #
    def executeSingleRunAllLocal( self, index , runString ):

        # before any futher, check we have the input file for this run
        self.checkInputFile( index, runString )

        # first generate the steering file
        self._steeringFileName = self.generateSteeringFileSingleRun( runString )

        # prepare a separate file for logging the output of Marlin
        self._logFileName = "%(name)s-%(run)s.log" % { "name": self.name, "run" : runString }

        # run marlin
        self.runMarlin( index, runString)

        # advice the user that Marlin is over
        self._logger.info( "Marlin finished successfully")

        # verify the presence of the output files
        self.checkOutputFile( index, runString )

        # verify the presence of the histogram files
        self.checkHistogramFile( index, runString )

        # prepare a tarbal for the records
        self.prepareTarball( runString )

        # check the presence of the joboutput tarball
        # this should be named something like
        # filter-123456.tar.gz
        self.checkJoboutputFile( index, runString )

        # clean up the local pc
        self.cleanup( runString )

    ## All local submitter with merge option
    def executeMergeAllLocal( self ):

        # before any futher, check all the input files
        for index, runString in enumerate( self._runStringList ) :
            try:
                self.checkInputFile( index, runString )
            except MissingInputFileError, error:
                self._logger.error( "The input file %(file)s is not available" % { "file": error._filename } )

                if self._configParser.get( "General", "Interactive" ):
                    if self.askYesNo( "Would you like to continue w/o this file? [y/n] " ):
                        self._logger.info( "User decided to continue w/o %(file)s" % { "file": error._filename } )
                        run, input, marlin, output, histo, tarball = self._summaryNTuple[ index ]
                        self._summaryNTuple[ index ] = run, "Skipped", "Skipped", "Skipped", "Skipped", "Skipped"
                        self._runStringList[ index ] = "DEADFACE"
                    else :
                        raise StopExecutionError( "Aborted by user" )
                else :
                    self._logger.warning( "Skipping %(file)s because missing" % { "file": error._filename } )
                    run, input, marlin, output, histo, tarball = self._summaryNTuple[ index ]
                    self._summaryNTuple[ index ] = run, "Skipped", "Skipped", "Skipped", "Skipped", "Skipped"


        # generate the steering file
        self._steeringFileName = self.generateSteeringFileMerge( )

        # prepare a separate file for logging the output of Marlin. This should be in the form:
        # %(name)s-%(output)s.log
        self._logFileName = "%(name)s-%(output)s.log" % { "name": self.name, "output": self._option.output }

        # run Marlin, this is exactly the same for the single run but these two 
        # small differences:
        #
        # 1. the index it is used only for filling the summaryNTuple, 
        #    let's use the very last index of the table.
        # 2. the runString is not used at all. Passing ""
        try :
            self.runMarlin( len( self._summaryNTuple ) - 1, None )

        except MarlinError, error:
            # this is critical, because for sure we can't continue 
            self._logger.critical( "Error with Marlin execution (%(msg)s - errno = %(errno)s )" \
                    % { "msg": error._message, "errno": error._errno } )

            # fix the summary table and stop
            for index, entry in enumerate ( self._summaryNTuple ):
                run, input, marlin, output, histo, tarball = entry
                if input != "Skipped" :
                    marlin  = "Failed"
                    output  = "Missing"
                    histo   = "Missing"
                    tarball = "Missing"
                    self._summaryNTuple[ index ] = run, input, marlin, output, histo, tarball

            # stop here
            raise StopExecutionError( "Error with Marlin. Terminating!" )

        # we need anyway to clean up the summaryNTuple
        for index, entry in enumerate ( self._summaryNTuple ):
            run, input, marlin, output, histo, tarball = entry
            if input != "Skipped" :
                marlin = "OK"
                if index < len( self._summaryNTuple ) -1 :
                    output  = "Merged"
                    histo   = "Merged"
                    tarball = "Merged"
                self._summaryNTuple[ index ] = run, input, marlin, output, histo, tarball

        # check the presence of the output file. This time it is one only
        try :
            self.checkOutputFile( len( self._summaryNTuple ) - 1, None )
        except MissingOutputFileError, error:
            self._logger.critical( "Missing output file %(file)s " % { "file": error._filename } )
            # fixing the summary table and stop the execution
            run, input, marlin, output, histo, tarball =  self._summaryNTuple[ len( self._summaryNTuple ) - 1]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, marlin, "Missing", histo, tarball
            raise StopExecutionError( "No output file was produced. Terminating" )

        # check the presence of the histogram file
        try :
            self.checkHistogramFile( len( self._summaryNTuple ) - 1, None )
        except MissingHistogramFileError, error:
            self._logger.critical( "Missing histogram file %(file)s " % { "file": error._filename } )
            # fixing the summary table and stop the execution
            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, marlin, output, "Missing", tarball
            raise StopExecutionError( "No histogram file was produced. Terminating" )


        # prepare the tarball
        self.prepareTarball( self._option.output )

        # check the presence of the joboutput file
        try:
            self.checkJoboutputFile( len( self._summaryNTuple ) - 1,  self._option.output  )
        except MissingJoboutputFileError, error:
            # fixing the summary table and stop the execution
            self._logger.critical( "Missing joboutput file %(file)s " % { "file": error._filename } )
            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, marlin, output, histo, "Missing"
            raise StopExecutionError( "No joboutput file was produced. Terminating" )

        # cleanup
        self.cleanup( self._option.output )

    ## CPU Local submitter
    #
    # This methods represents the sequence of commands that should be done while
    # submitting jobs on the local computer but using remote data
    #
    def executeSingleRunCPULocal( self, index , runString ):

        # do preliminary checks
        self.doPreliminaryTest( index, runString )

        # get the input file from the GRID
        self.getRunFromGRID( index, runString )

        # double check the presence of the input file
        self.checkInputFile( index, runString )

        #  generate the steering file
        self._steeringFileName = self.generateSteeringFileSingleRun( runString )

        # prepare a separate file for logging the output of Marlin
        self._logFileName = "%(name)s-%(run)s.log" % { "name": self.name,"run" : runString }

        # run marlin
        self.runMarlin( index, runString )

        # advice the user that Marlin is over
        self._logger.info( "Marlin finished successfully")

        # verify the presence of the output files
        self.checkOutputFile( index, runString )

        # prepare a tarbal for the records
        self.prepareTarball( runString )

        try :
            # copy the DB file to the GRID
            self.putOutputOnGRID( index, runString )

            # copy the histogram file to the GRID
            self.putHistogramOnGRID( index, runString )

            # copy the joboutput to the GRID
            self.putJoboutputOnGRID( index, runString )

            # clean up the local pc
            self.cleanup( runString )

        except GRID_LCG_CRError, error:
            message = "The file (%(file)s) couldn't be copied on the GRID"  % { "file": error._filename }
            self._logger.error( message )

    ## CPU Local submitter with merge
    # 
    def executeMergeCPULocal( self ) :

        # do preliminary test over all input runs
        self._doFullCheck = True
        try:
            for index, runString in enumerate( self._runStringList ) :
                try :
                    self.doPreliminaryTest( index, runString, self._doFullCheck )
                    self._doFullCheck = False
                except MissingInputFileOnGRIDError, error:
                    self._logger.error( "The input file %(file)s is not available" % { "file": error._filename } )

                    if self._configParser.get( "General", "Interactive" ):
                        if self.askYesNo( "Would you like to continue w/o this file? [y/n] " ):
                            self._logger.info( "User decided to continue w/o %(file)s" % { "file": error._filename } )
                            run, input, marlin, output, histo, tarball = self._summaryNTuple[ index ]
                            self._summaryNTuple[ index ] = run, "Skipped", "Skipped", "Skipped", "Skipped", "Skipped"
                            self._runStringList[ index ] = "DEADFACE"
                        else :
                            raise StopExecutionError( "Aborted by user" )
                    else :
                        self._logger.warning( "Skipping %(file)s because missing" % { "file": error._filename } )
                        run, input, marlin, output, histo, tarball = self._summaryNTuple[ index ]
                        self._summaryNTuple[ index ] = run, "Skipped", "Skipped", "Skipped", "Skipped", "Skipped"



        except OutputFileAlreadyOnGRIDError, error: 
            self._logger.critical( "Output file %(file)s already on GRID" % { "file": error._filename } )
            raise StopExecutionError( "Cannot continue" )

        except JoboutputFileAlreadyOnGRIDError, error: 
            self._logger.critical( "Joboutput file %(file)s already on GRID" % { "file": error._filename } )
            raise StopExecutionError( "Cannot continue" )

        except HistogramFileAlreadyOnGRIDError, error: 
            self._logger.critical( "Histogram file %(file)s already on GRID" % { "file": error._filename } )
            raise StopExecutionError( "Cannot continue" )
        # end of preliminary checks


    ## Get the input run from the GRID
    def getRunFromGRID(self, index, runString ) :

        self._logger.info(  "Getting the input file from the GRID" )

        try :
            inputPath = self._configParser.get( "GRID", "GRIDFolderClusearchResults" )
        except ConfigParser.NoOptionError :
            message = "GRIDFolderClusearchResults missing in the configuration file. Quitting."
            self._logger.critical( message )
            raise StopExecutionError( message )

        try :
            localPath = self._configParser.get( "LOCAL", "LocalFolderClusearchResults" )
        except ConfigParser.NoOptionError :
            localPath = "$PWD/results"

        baseCommand = "lcg-cp "
        if self._option.verbose :
            baseCommand = baseCommand + " -v "

        command = "%(base)s lfn:%(gridPath)s/run%(run)s-clu-p%(pede)s.slcio file:%(localPath)s/run%(run)s.slcio" %  \
            { "base": baseCommand, "gridPath" : inputPath, "pede": self._pedeString, "run": runString, "localPath": localPath }
        if os.system( command ) != 0:
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, "Missing", c, d, e, f
            raise GRID_LCG_CPError( "lfn:%(gridPath)s/run%(run)s-clu-p%(pede)s.slcio" % \
                                        { "gridPath" : inputPath, "run": runString, "pede": self._pedeString } )
        else:
            self._logger.info("Input file successfully copied from the GRID")
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, "OK", c, d, e, f


    ## Put the DB to the GRID
    def putOutputOnGRID( self, index, runString ):

        self._logger.info(  "Putting the output file to the GRID" )

        try :
            gridPath = self._configParser.get( "GRID", "GRIDFolderFilterResults")
        except ConfigParser.NoOptionError :
            message = "GRIDFolderFilterResults missing in the configuration file. Quitting."
            self._logger.critical( message )
            raise StopExecutionError( message )

        try :
            localPath = self._configParser.get( "LOCAL", "LocalFolderFilterResults" )
        except ConfigParser.NoOptionError :
            localPath = "$PWD/results"

        baseCommand = "lcg-cr "
        if self._option.verbose :
            baseCommand = baseCommand + " -v "

        command = "%(base)s -l lfn:%(gridFolder)s/run%(run)s-filter-p%(pede)s.slcio file:%(localFolder)s/run%(run)s-filter-p%(pede)s.slcio" % \
            { "base": baseCommand, "gridFolder": gridPath, "localFolder": localPath, "run" : runString , "pede": self._pedeString }
        if os.system( command ) != 0 :
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, c, "LOCAL", e, f
            raise GRID_LCG_CRError( "lfn:%(gridFolder)s/run%(run)s-filter-p%(pede)s.slcio" % \
                                        { "gridFolder": gridPath, "run" : runString , "pede": self._pedeString} )
        else:
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, c, "GRID", e, f
            self._logger.info( "Output file successfully copied to the GRID" )

    ## Put the histograms file to the GRID
    def putHistogramOnGRID( self, index, runString ):

        self._logger.info( "Putting the histogram file to the GRID" )

        try:
            gridPath = self._configParser.get( "GRID", "GRIDFolderFilterHisto")
        except ConfigParser.NoOptionError :
            message = "GRIDFolderFilterHisto missing in the configuration file. Quitting."
            self._logger.critical( message )
            raise StopExecutionError( message )

        try :
            localPath = self._configParser.get( "LOCAL", "LocalFolderFilterHisto" )
        except ConfigParser.NoOptionError :
            localPath = "$PWD/histo"

        baseCommand = "lcg-cr "
        if self._option.verbose :
            baseCommand = baseCommand + " -v "

        command = "%(base)s -l lfn:%(gridFolder)s/run%(run)s-filter-histo.root file:%(localFolder)s/run%(run)s-filter-histo.root" % \
            { "base": baseCommand, "gridFolder": gridPath, "localFolder": localPath, "run" : runString }
        if os.system( command ) != 0 :
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, c, d, "LOCAL", f
            raise GRID_LCG_CRError( "lfn:%(gridFolder)s/run%(run)s-filter-histo.root" % \
                                        { "gridFolder": gridPath, "run" : runString } )
        else:
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, c, d, "GRID", f
            self._logger.info( "Histogram file successfully copied to the GRID" )

    ## Put the joboutput file to the GRID
    def putJoboutputOnGRID( self, index, runString ):

        self._logger.info( "Putting the joboutput file to the GRID" )

        try:
            gridPath = self._configParser.get( "GRID", "GRIDFolderFilterJoboutput")
        except ConfigParser.NoOptionError :
            message = "GRIDFolderFilterJoboutout missing in the configuration file. Quitting."
            self._logger.critical( message )
            raise StopExecutionError( message )

        try :
            localPath = self._configParser.get( "LOCAL", "LocalFolderFilterJoboutput" )
        except ConfigParser.NoOptionError :
            localPath = "$PWD/log"

        baseCommand = "lcg-cr "
        if self._option.verbose :
            baseCommand = baseCommand + " -v "

        command = "%(base)s -l lfn:%(gridFolder)s/%(name)s-%(run)s.tar.gz file:%(localFolder)s/%(name)s-%(run)s.tar.gz" % \
            { "name": self.name, "base": baseCommand, "gridFolder": gridPath, "localFolder": localPath, "run" : runString }
        if os.system( command ) != 0 :
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, c, d, e, "LOCAL"
            raise GRID_LCG_CRError( "lfn:%(gridFolder)s/%(name)s-%(run)s.tar.gz" % \
                                        { "name": self.name, "gridFolder": gridPath, "run" : runString } )
        else:
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, c, d, e, "GRID"
            self._logger.info( "Joboutput file successfully copied to the GRID" )

    ## Check the input file
    #
    def checkInputFile( self, index, runString ) :
        # the input file should be something like:
        # results/run123456-clu-p654321.slcio

        self._logger.info( "Checking the input file" )

        try :
            inputFilePath = self._configParser.get( "LOCAL", "LocalFolderClusearchResults" )
        except ConfigParser.NoOptionError :
            inputFilePath = "results"
        inputFileName = "run%(run)s-clu-p%(pede)s.slcio" % { "run": runString , "pede": self._pedeString }
        if not os.access( os.path.join(inputFilePath, inputFileName) , os.R_OK ):
            message = "Problem accessing the input file (%(file)s), trying next run" % {"file": inputFileName }
            self._logger.error( message )
            raise MissingInputFileError( inputFileName )
        else :
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, "OK", c, d, e, f


    ## Generate the steering file
    def generateSteeringFileSingleRun( self, runString  ) :
        message = "Generating the steering file (%(name)s-%(run)s.xml) " % { "name": self.name, "run" : runString }
        self._logger.info( message )

        steeringFileTemplate =  ""
        try:
            steeringFileTemplate = self._configParser.get( "SteeringTemplate", "FilterSteeringFile" )
        except ConfigParser.NoOptionError :
            steeringFileTemplate = "template/%(name)s-tmp.xml" % { "name": self.name }

        message = "Using %(file)s as steering template" %{ "file": steeringFileTemplate }
        self._logger.debug( message )

        # check if the template exists
        if not os.path.exists( steeringFileTemplate ) :
            raise MissingSteeringTemplateError ( steeringFileTemplate )

        # open the template for reading
        templateSteeringFile = open( steeringFileTemplate , "r")

        # read the whole content
        templateSteeringString = templateSteeringFile.read()

        # make all the changes
        actualSteeringString = templateSteeringString

        # replace the file paths
        #
        # first the gear path
        if self._option.execution == "all-grid" :
            self._gearPath = "."
        else:
            try :
                self._gearPath = self._configParser.get( "LOCAL", "LocalFolderGear" )
            except ConfigParser.NoOptionError :
                self._gearPath = ""

        actualSteeringString = actualSteeringString.replace("@GearPath@", self._gearPath )

        # find the gear file, first check the configuration file, then the commad line options
        # and last use a default gear_telescope.xml
        self._gear_file = ""
        try:
            self._gear_file = self._configParser.get( "General", "GEARFile" )
        except ConfigParser.NoOptionError :
            self._logger.debug( "No GEAR file in the configuration file" )

        if self._option.gear_file != None :
            # this means that the user wants to override the configuration file
            self._gear_file = self._option.gear_file
            self._logger.debug( "Using command line GEAR file" )


        if self._gear_file == "" :
            # using default GEAR file
            defaultGEARFile = "gear_telescope.xml"
            self._gear_file = defaultGEARFile
            message = "Using default GEAR file %(gear)s" %{ "gear": defaultGEARFile }
            self._logger.warning( message )

        actualSteeringString = actualSteeringString.replace("@GearFile@", self._gear_file )

        # now repeat the game with the histo info xml file even if it is not compulsory to have it
        if self._option.execution == "all-grid" :
            self._histoinfoPath = "."
        else:
            try :
                self._histoinfoPath = self._configParser.get( "LOCAL", "LocalFolderHistoinfo" )
            except ConfigParser.NoOptionError :
                self._histoinfoPath = ""

        actualSteeringString = actualSteeringString.replace("@HistoInfoPath@", self._histoinfoPath )

        self._histoinfoFilename = ""
        try:
            self._histoinfoFilename = self._configParser.get( "General", "Histoinfo" )
        except ConfigParser.NoOptionError :
            self._logger.debug( "No histoinfo file in the configuration file, using default histoinfo.xml" )
            self._histoinfoFilename = "histoinfo.xml"

        actualSteeringString = actualSteeringString.replace("@HistoInfo@", self._histoinfoFilename )

        # now replace the input folder path
        if self._option.execution == "all-grid" :
            inputFolder = "results"
        else :
            try:
                inputFolder = self._configParser.get( "LOCAL", "LocalFolderClusearchResults" )
            except ConfigParser.NoOptionError :
                inputFolder = "results"
        actualSteeringString = actualSteeringString.replace("@InputPath@", inputFolder )

        # now replace the pedestal folder path
        if self._option.execution == "all-grid" :
            dbFolder = "db"
        else:
            try:
                dbFolder = self._configParser.get("LOCAL", "LocalFolderDBPede")
            except ConfigParser.NoOptionError :
                dbFolder = "db"
        actualSteeringString = actualSteeringString.replace("@DBPath@" ,dbFolder )

        # that's a good place to also replace the pedestal run number
        actualSteeringString = actualSteeringString.replace("@PedeRunNumber@", self._pedeString )

        # now replace the histo folder path
        if self._option.execution == "all-grid" :
            histoFolder = "histo"
        else:
            try:
                histoFolder = self._configParser.get("LOCAL", "LocalFolderFilterHisto")
            except ConfigParser.NoOptionError :
                histoFolder = "histo"
        actualSteeringString = actualSteeringString.replace("@HistoPath@" ,histoFolder )

        # now replace the output folder path
        if self._option.execution == "all-grid" :
            resultFolder = "results"
        else :
            try :
                resultFolder = self._configParser.get("LOCAL", "LocalFolderFilterResults")
            except ConfigParser.NoOptionError :
                resultFolder = "results"
        actualSteeringString = actualSteeringString.replace( "@ResultPath@", resultFolder )

        # now fix the issue with the name of the outputs
        actualSteeringString = actualSteeringString.replace( "@Base@", "run" )

        # finally replace the run string !
        actualSteeringString = actualSteeringString.replace("@RunNumber@", runString )

        # open the new steering file for writing
        steeringFileName = "%(name)s-%(run)s.xml" % { "name": self.name, "run" : runString }
        actualSteeringFile = open( steeringFileName, "w" )

        # write the new steering file
        actualSteeringFile.write( actualSteeringString )

        # close both files
        actualSteeringFile.close()
        templateSteeringFile.close()

        return steeringFileName

    ## Generate the steering file for merging
    def generateSteeringFileMerge( self ):
        self._logger.info( "Generating the steering file (%(name)s-%(output)s.xml) " %
                           { "name": self.name, "output": self._option.output } )

        steeringFileTemplate =  ""
        try:
            steeringFileTemplate = self._configParser.get( "SteeringTemplate", "FilterSteeringFile" )
        except ConfigParser.NoOptionError :
            steeringFileTemplate = "template/%(name)s-tmp.xml" % { "name": self.name }

        message = "Using %(file)s as steering template" %{ "file": steeringFileTemplate }
        self._logger.debug( message )

        # check if the template exists
        if not os.path.exists( steeringFileTemplate ) :
            raise MissingSteeringTemplateError ( steeringFileTemplate )

        # open the template for reading
        templateSteeringFile = open( steeringFileTemplate , "r")

        # read the whole content
        templateSteeringString = templateSteeringFile.read()

        # make all the changes
        actualSteeringString = templateSteeringString

        # replace the file paths
        #
        # first the gear path
        if self._option.execution == "all-grid" :
            self._gearPath = "."
        else:
            try :
                self._gearPath = self._configParser.get( "LOCAL", "LocalFolderGear" )
            except ConfigParser.NoOptionError :
                self._gearPath = ""

        actualSteeringString = actualSteeringString.replace("@GearPath@", self._gearPath )

        # find the gear file, first check the configuration file, then the commad line options
        # and last use a default gear_telescope.xml
        self._gear_file = ""
        try:
            self._gear_file = self._configParser.get( "General", "GEARFile" )
        except ConfigParser.NoOptionError :
            self._logger.debug( "No GEAR file in the configuration file" )

        if self._option.gear_file != None :
            # this means that the user wants to override the configuration file
            self._gear_file = self._option.gear_file
            self._logger.debug( "Using command line GEAR file" )


        if self._gear_file == "" :
            # using default GEAR file
            defaultGEARFile = "gear_telescope.xml"
            self._gear_file = defaultGEARFile
            message = "Using default GEAR file %(gear)s" %{ "gear": defaultGEARFile }
            self._logger.warning( message )

        actualSteeringString = actualSteeringString.replace("@GearFile@", self._gear_file )

        # now repeat the game with the histo info xml file even if it is not compulsory to have it
        if self._option.execution == "all-grid" :
            self._histoinfoPath = "."
        else:
            try :
                self._histoinfoPath = self._configParser.get( "LOCAL", "LocalFolderHistoinfo" )
            except ConfigParser.NoOptionError :
                self._histoinfoPath = ""

        actualSteeringString = actualSteeringString.replace("@HistoInfoPath@", self._histoinfoPath )

        self._histoinfoFilename = ""
        try:
            self._histoinfoFilename = self._configParser.get( "General", "Histoinfo" )
        except ConfigParser.NoOptionError :
            self._logger.debug( "No histoinfo file in the configuration file, using default histoinfo.xml" )
            self._histoinfoFilename = "histoinfo.xml"

        actualSteeringString = actualSteeringString.replace("@HistoInfo@", self._histoinfoFilename )

        # now replace the pedestal folder path
        if self._option.execution == "all-grid" :
            dbFolder = "db"
        else:
            try:
                dbFolder = self._configParser.get("LOCAL", "LocalFolderDBPede")
            except ConfigParser.NoOptionError :
                dbFolder = "db"
        actualSteeringString = actualSteeringString.replace("@DBPath@" ,dbFolder )

        # now replace the histo folder path
        if self._option.execution == "all-grid" :
            histoFolder = "histo"
        else:
            try:
                histoFolder = self._configParser.get("LOCAL", "LocalFolderFilterHisto")
            except ConfigParser.NoOptionError :
                histoFolder = "histo"
        actualSteeringString = actualSteeringString.replace("@HistoPath@" ,histoFolder )

        # now replace the output folder path
        if self._option.execution == "all-grid" :
            resultFolder = "results"
        else :
            try :
                resultFolder = self._configParser.get("LOCAL", "LocalFolderFilterResults")
            except ConfigParser.NoOptionError :
                resultFolder = "results"
        actualSteeringString = actualSteeringString.replace( "@ResultPath@", resultFolder )


        # now fix the issue with the name of the outputs
        actualSteeringString = actualSteeringString.replace( "@Base@@RunNumber@", self._option.output )

        # now we simply need to fix the input files
        goodMerging = False
        for index, runString in enumerate(self._runStringList):
            if runString != "DEADFACE":
                goodMerging = goodMerging or True
                if self._option.execution == "all-grid":
                    inputFolder = "results"
                else:
                    try:
                        inputFolder = self._configParser.get( "LOCAL", "LocalFolderClusearchResults" )
                    except ConfigParser.NoOptionError :
                        inputFolder = "results"
                actualSteeringString = actualSteeringString.replace( "@InputPath@", inputFolder )
                actualSteeringString = actualSteeringString.replace( "@RunNumber@", runString )
                actualSteeringString = actualSteeringString.replace( "@PedeRunNumber@", self._pedeString )
                file = "%(path)s/run%(run)s-clu-p%(pede)s.slcio" % { "path" : inputFolder, "run": runString, "pede": self._pedeString }
                actualSteeringString = actualSteeringString.replace( file, file + " @InputPath@/run@RunNumber@-clu-p@PedeRunNumber@.slcio" )

        if not goodMerging:
            message = "No available input files!"
            self._logger.critical( message )
            raise StopExecutionError( message )

        actualSteeringString = actualSteeringString.replace( " @InputPath@/run@RunNumber@-clu-p@PedeRunNumber@.slcio", "" )

        # in case there are still some un-replaced pederunnumber do it now!
        actualSteeringString = actualSteeringString.replace( "@PedeRunNumber@", self._pedeString )

        # open the new steering file for writing
        steeringFileName = "%(name)s-%(base)s.xml" % { "name": self.name, "base" : self._option.output }
        actualSteeringFile = open( steeringFileName, "w" )

        # write the new steering file
        actualSteeringFile.write( actualSteeringString )

        # close both files
        actualSteeringFile.close()
        templateSteeringFile.close()

        return steeringFileName

    ## Execute Marlin
    def runMarlin( self, index, runString  ) :

        self._logger.info( "Running Marlin" )

        # first check that the gear file exists
        if not os.path.exists( os.path.join( self._gearPath, self._gear_file )) :
            raise MissingGEARFileError(   os.path.join( self._gearPath, self._gear_file ) )

        # do some tricks for having the logfile
        logFile = open( self._logFileName, "w")
        marlin  = popen2.Popen4( "Marlin %(steer)s" % { "steer": self._steeringFileName } )
        while marlin.poll() == -1:
            line = marlin.fromchild.readline()
            print line.strip()
            logFile.write( line )

        logFile.close()
        returnValue = marlin.poll()
        if returnValue != 0:
            raise MarlinError( "", returnValue )
        else :
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, "OK", d, e, f
            return returnValue

    ## Check the output file
    #
    def checkOutputFile( self, index, runString ) :
        # this should be named something like
        # results/run123456-filter-p654321.slcio or 
        # results/%(output)s-filter-p654321.slcio

        try :
            outputFilePath = self._configParser.get( "LOCAL", "LocalFolderFilterResults" )
        except ConfigParser.NoOptionError :
            outputFilePath = "results"

        if self._option.merge :
            outputFileName = "%(output)s-filter-p%(pede)s.slcio"  % { "output": self._option.output, "pede": self._pedeString }
        else :
            outputFileName = "run%(run)s-filter-p%(pede)s.slcio" % { "run": runString, "pede": self._pedeString }
        self._logger.log(15, "Checking for file %(file)s " % { "file":  os.path.join( outputFilePath , outputFileName) })
        if not os.access( os.path.join( outputFilePath , outputFileName) , os.R_OK ):
            raise MissingOutputFileError( outputFileName )
        else :
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, c, "OK", e, f

    ## Check the histo file
    #
    def checkHistogramFile( self, index, runString ) :
        # this should be named something like
        # histo/run123456-filter-histo.root

        try :
            histoFilePath = self._configParser.get( "LOCAL", "LocalFolderFilterHisto" )
        except ConfigParser.NoOptionError :
            histoFilePath = "histo"
        if self._option.merge :
            histoFileName = "%(output)s-filter-histo.root" % { "output": self._option.output }
        else:
            histoFileName = "run%(run)s-filter-histo.root" % { "run": runString }

        if not os.access( os.path.join( histoFilePath , histoFileName ) , os.R_OK ):
            raise MissingHistogramFileError( histoFileName )
        else:
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, c, d, "OK", f

    ## Prepare the joboutput tarball
    def prepareTarball( self, runString ) :
        self._logger.info("Preparing the joboutput tarball" )

        # first prepare a folder to store them temporary
        destFolder = "%(name)s-%(run)s" %{ "name": self.name, "run": runString }
        shutil.rmtree( destFolder, True )
        os.mkdir( destFolder )

        # prepare the list of files we want to copy.
        listOfFiles = []

        # the gear file, the histo info and the config.cfg
        listOfFiles.append( os.path.join( self._gearPath, self._gear_file ) )
        try :
            histoinfoPath = self._configParser.get( "LOCAL", "LocalFolderHistoinfo" )
        except ConfigParser.NoOptionError :
            histoinfoPath = ""
        if os.access( os.path.join( histoinfoPath, self._histoinfoFilename ), os.R_OK ) :
            listOfFiles.append( os.path.join( histoinfoPath, self._histoinfoFilename ))
        listOfFiles.append( self._configFile )

        # all files starting with filter-123456 in the local folder
        for file in glob.glob( "%(name)s-*.*" % {"name": self.name } ):
            message = "Adding %(file)s to the joboutput tarball" % { "file": file }
            self._logger.debug( message )
            listOfFiles.append( file )

        # the histogram file
        if self._option.execution != "all-grid":
            try :
                histoFilePath = self._configParser.get( "LOCAL", "LocalFolderFilterHisto" )
            except ConfigParser.NoOptionError :
                histoFilePath = "histo"
                if self._option.merge:
                    histoFileName = "%(output)s-filter-histo.root" % { "output": self._option.output }
                else :
                    histoFileName =  "run%(run)s-filter-histo.root" % { "run": runString }
                listOfFiles.append( os.path.join( histoFilePath, histoFileName ) )

        # copy everything into a temporary folder
        for file in listOfFiles :
            shutil.copy( file, destFolder )

        # do the tarball
        self._tarballFileName = "%(name)s-%(run)s.tar.gz" % { "name": self.name, "run": runString }
        tarball = tarfile.open( self._tarballFileName, "w:gz" )
        tarball.add( destFolder )
        tarball.close()

        # remove the temporary folder
        shutil.rmtree( destFolder )

        # copy the tarball in the log folder
        try:
            localFolder = self._configParser.get( "LOCAL", "LocalFolderFilterJoboutput")
        except ConfigParser.NoOptionError :
            self._logger.debug( "LocalFolderFilterJoboutput not available, using $PWD/log ")
            localFolder = "log/"

        shutil.move( self._tarballFileName, localFolder )


    ## Check the joboutput file
    #
    def checkJoboutputFile( self, index, runString) :
        # this should be named something like
        # log/filter-123456.tar.gz

        try :
            outputFilePath = self._configParser.get( "LOCAL", "LocalFolderFilterJoboutput" )
        except ConfigParser.NoOptionError :
            outputFilePath = "log"
        tarballFileName = "%(name)s-%(run)s.tar.gz" % { "name": self.name, "run": runString }
        if not os.access( os.path.join( outputFilePath, tarballFileName), os.R_OK ):
            raise MissingJoboutputFileError( tarballFileName )
        else:
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, c, d, e, "OK"


    ## Cleanup after each run pedestal calculation
    def cleanup( self, runString ):

        self._logger.info( "Cleaning up the local pc" )

        # remove the log file and the steering file
        for file in glob.glob( "%(name)s-*" % {"name": self.name}):
            os.remove( file )


        # remove the input and output file
        if self._keepInput == False:
            try :
                inputFilePath = self._configParser.get( "LOCAL", "LocalFolderClusearchResults" )
            except ConfigParser.NoOptionError :
                inputFilePath = "results"

            if self._option.merge :
                for runString in self._runStringList:
                    if runString != "DEADFACE" :
                        inputFile  = "run%(run)s-clu-p%(pede)s.slcio" % { "run" : runString, "pede": self._pedeString }
                        os.remove( os.path.join( inputFilePath, inputFile ))
            else :
                inputFile  = "run%(run)s-clu-p%(pede)s.slcio" % { "run" : runString, "pede": self._pedeString }
                os.remove( os.path.join( inputFilePath, inputFile ))

        if self._keepOutput == False :
            try :
                outputFilePath = self._configParser.get( "LOCAL", "LocalFolderFilterResults" )
            except ConfigParser.NoOptionError :
                outputFilePath = "results"

            if self._option.merge:
                outputFile = "%(run)s-filter-p%(pede)s.slcio" % { "run" : runString, "pede": self._pedeString }
            else:
                outputFile = "run%(run)s-filter-p%(pede)s.slcio" % { "run" : runString, "pede": self._pedeString }
            for file in glob.glob( os.path.join( outputFilePath , outputFile ) ):
                os.remove( file )

            try :
                histoFilePath = self._configParser.get( "LOCAL", "LocalFolderFilterHisto" )
            except ConfigParser.NoOptionError :
                histoFilePath = "histo"

            if self._option.merge:
                histoFile = "%(run)s-filter-histo.root" % { "run": runString }
            else:
                histoFile = "run%(run)s-filter-histo.root" % { "run": runString }
            os.remove( os.path.join( histoFilePath, histoFile ) )


    def end( self ) :

        if self._keepInput == False :
            # remove also the pedestal run
            try:
                dbFilePath = self._configParser.get( "LOCAL", "LocalFodlerDB" )
            except ConfigParser.NoOptionError :
                dbFilePath = "db"

            os.remove( os.path.join( dbFilePath , self._pedeFilename ) )

        if self._option.execution == "all-grid" :
            self.prepareJIDFile()
            self.logGRIDJobs( )

        SubmitBase.end( self )


    def logGRIDJobs( self ):
        self._logger.info( "" )
        self._logger.info( "== GRID JOB ID ==============================================================" )
        for entry in self._gridJobNTuple :
            run, jid = entry
            message = "| %(run)6s | %(jid)64s |" % { "run" : run, "jid":jid.strip() }
            self._logger.info( message )

        self._logger.info("=============================================================================" )
        self._logger.info( "" )

    def prepareJIDFile( self ):
        unique = datetime.datetime.fromtimestamp( self._timeBegin ).strftime("%Y%m%d-%H%M%S")
        self._logger.info("Preparing the JID for this submission (%(name)s-%(date)s.jid)" % {
                "name": self.name, "date" : unique } )

        try :
            localPath = self._configParser.get( "LOCAL", "LocalFolderFilterJoboutput")
        except ConfigParser.NoOptionError :
            localPath = "log/"

        jidFile = open( os.path.join( localPath, "%(name)s-%(date)s.jid" % { "name": self.name, "date": unique } ), "w" )
        for run, jid in self._gridJobNTuple:
            if jid != "Unknown":
                jidFile.write( jid )

        jidFile.close()



    ## Preliminary checks
    #
    # This method performs some preliminary checks
    # before starting a full GRID or a CPU local submission
    # In particular, it checks that the proxy exists and it is still valid,
    # the input file is available on the GRID at the specified path and that
    # all the output files are not already on the Storage Element.
    #
    # In case the proxy is missing or expired a StopExecutionError will be raised,
    # while other errors will be thrown in case of non existing input files or
    # already existing output files.
    #
    # @throw StopExecutionError
    # @throw MissingFileOnGRIDError
    # @throw FileAlreadyOnGRIDError
    #
    def doPreliminaryTest( self, index, runString, fullCheck = True ) :

        if fullCheck :
            self._logger.log(15, "Performing a full check" )

        if fullCheck :
            # first log the voms-proxy-info
            self._logger.info( "Logging the voms-proxy-info" )
            info = popen2.Popen4("voms-proxy-info -all")
            while info.poll() == -1:
                line = info.fromchild.readline()
                self._logger.info( line.strip() )

            if info.poll() != 0:
                message = "Problem with the GRID_UI"
                self._logger.critical( message )
                raise StopExecutionError( message )

            # also check that the proxy is still valid
            if os.system( "voms-proxy-info -e" ) != 0:
                message = "Expired proxy"
                self._logger.critical( message )
                raise StopExecutionError( message )
            else:
                self._logger.info( "Valid proxy found" )

            # get all the needed path from the configuration file
            try :
                self._inputPathGRID     = self._configParser.get("GRID", "GRIDFolderClusearchResults")
                self._pedePathGRID      = self._configParser.get("GRID", "GRIDFolderDBPede" )
                self._outputPathGRID    = self._configParser.get("GRID", "GRIDFolderFilterResults" )
                self._joboutputPathGRID = self._configParser.get("GRID", "GRIDFolderFilterJoboutput")
                self._histogramPathGRID = self._configParser.get("GRID", "GRIDFolderFilterHisto")
                folderList =  [ self._outputPathGRID, self._joboutputPathGRID, self._histogramPathGRID ]
            except ConfigParser.NoOptionError:
                message = "Missing path from the configuration file"
                self._logger.critical( message )
                raise StopExecutionError( message )


        # check if the input file is on the GRID, otherwise go to next run
        # the existence of the pedestal run has been assured already.
        command = "lfc-ls %(inputPathGRID)s/run%(run)s-clu-p%(pede)s.slcio" % { "inputPathGRID" : self._inputPathGRID,  
                                                                                "pede": self._pedeString, "run": runString }

        self._logger.log(15, "Check input file %(path)s/run%(run)s-clu-p%(pede)s.slcio" % { "path" : self._inputPathGRID,  
                                                                                           "pede": self._pedeString, "run": runString } )
        lfc = popen2.Popen4( command )
        while lfc.poll() == -1:
            pass
        if lfc.poll() == 0:
            self._logger.info( "Input file found on the SE" )
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, "GRID", c, d, e, f
        else:
            self._logger.error( "Input file NOT found on the SE. Trying next run" )
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, "Missing", c, d, e, f
            raise MissingInputFileOnGRIDError( "%(inputPathGRID)s/run%(run)s-clu-p%(pede)s.slcio" % { "inputPathGRID" : self._inputPathGRID,  
                                                                                                      "pede": self._pedeString, "run": runString } )

        if fullCheck :
            # check the existence of the folders
            try :
                for folder in folderList:
                    self.checkGRIDFolder( folder )

            except MissingGRIDFolderError, error :
                message = "Folder %(folder)s is unavailable. Quitting" % { "folder": error._filename }
                self._logger.critical( message )
                raise StopExecutionError( message )


        # check if the output file already exists
        if self._option.merge:
            outputFilename = "%(run)s-filter-p%(pede)s.slcio" % { "pede": self._pedeString, "run": self._option.output }
        else:
            outputFilename = "run%(run)s-filter-p%(pede)s.slcio" % { "pede": self._pedeString, "run": runString }

        self._logger.log(15, "Check output file %(path)s/%(file)s" % {"path":  self._outputPathGRID, "file": outputFilename } )

        if fullCheck:
            command = "lfc-ls %(outputPathGRID)s/%(file)s" % { "outputPathGRID": self._outputPathGRID, "file": outputFilename }
            lfc = popen2.Popen4( command )
            while lfc.poll() == -1:
                pass
            if lfc.poll() == 0:
                self._logger.warning( "Output file %(outputPathGRID)s/%(file)s" % { "outputPathGRID": self._outputPathGRID, "file": outputFilename })
                if self._configParser.get("General","Interactive" ):
                    if self.askYesNo( "Would you like to remove it?  [y/n] " ):
                        self._logger.info( "User decided to remove %(outputPathGRID)s/%(file)s from the GRID"
                                           % { "outputPathGRID": self._outputPathGRID, "file": outputFilename } )
                        command = "lcg-del -a lfn:%(outputPathGRID)s/%(file)s" % { "outputPathGRID": self._outputPathGRID, "file": outputFilename }
                        os.system( command )
                    else :
                        raise OutputFileAlreadyOnGRIDError( "%(outputPathGRID)s/%(file)s on the GRID"
                                                  % { "outputPathGRID": self._outputPathGRID, "file": outputFilename } )
                else :
                    raise OutputFileAlreadyOnGRIDError( "%(outputPathGRID)s/%(file)s on the GRID"
                                                        % { "outputPathGRID": self._outputPathGRID, "file": outputFilename } )

            # check if the job output file already exists
            if self._option.merge:
                outputFilename =  "%(name)s-%(run)s.tar.gz" % { "name": self.name, "run": self._option.output }
            else:
                outputFilename =  "%(name)s-%(run)s.tar.gz" % { "name": self.name, "run": runString }

            self._logger.log(15, "Check joboutput file %(path)s/%(file)s" % {"path":  self._joboutputPathGRID, "file": outputFilename } )
            command = "lfc-ls %(outputPathGRID)s/%(file)s" % { "outputPathGRID": self._joboutputPathGRID, "file": outputFilename }
            lfc = popen2.Popen4( command )
            while lfc.poll() == -1:
                pass
            if lfc.poll() == 0:
                self._logger.warning( "Joboutput file %(outputPathGRID)s/%(file)s already exists"
                                      % {"outputPathGRID": self._joboutputPathGRID, "file": outputFilename } )
                if self._configParser.get("General","Interactive" ):
                    if self.askYesNo( "Would you like to remove it?  [y/n] " ):
                        self._logger.info( "User decided to remove %(outputPathGRID)s/%(file)s from the GRID"
                                           % { "outputPathGRID": self._joboutputPathGRID, "file": outputFilename } )
                        command = "lcg-del -a lfn:%(outputPathGRID)s/%(file)s" % { "outputPathGRID": self._joboutputPathGRID, "file": outputFilename }
                        os.system( command )
                    else :
                        raise JoboutputFileAlreadyOnGRIDError( "%(outputPathGRID)s/%(file)s on the GRID"
                                                               % { "outputPathGRID": self._joboutputPathGRID, "file": outputFilename } )
                else :
                    raise JoboutputFileAlreadyOnGRIDError( "%(outputPathGRID)s/%(file)s on the GRID"
                                                           % {"outputPathGRID": self._joboutputPathGRID,"file": outputFilename} )


            # check if the histogram file already exists
            if self._option.merge:
                outputFilename =  "%(run)s-filter-histo.root" % { "run": self._option.output }
            else:
                outputFilename =  "run%(run)s-filter-histo.root" % { "run": runString }

            self._logger.log(15, "Check histogram file %(path)s/%(file)s" % {"path":  self._histogramPathGRID, "file": outputFilename } )
            command = "lfc-ls %(outputPathGRID)s/%(file)s" % { "outputPathGRID": self._histogramPathGRID, "file": outputFilename }
            lfc = popen2.Popen4( command )
            while lfc.poll() == -1:
                pass
            if lfc.poll() == 0:
                self._logger.warning( "Histogram file %(outputPathGRID)s/%(file)s already exists"
                                      % { "outputPathGRID": self._histogramPathGRID, "file": outputFilename } )
                if self._configParser.get("General","Interactive" ):
                    if self.askYesNo( "Would you like to remove it?  [y/n] " ):
                        self._logger.info( "User decided to remove %(outputPathGRID)s/%(file)s from the GRID"
                                           % { "outputPathGRID": self._histogramPathGRID, "file": outputFilename } )
                        command = "lcg-del -a lfn:%(outputPathGRID)s/%(file)s" % { "outputPathGRID": self._histogramPathGRID, "file": outputFilename }
                        os.system( command )
                    else :
                        raise HistogramFileAlreadyOnGRIDError( "%(outputPathGRID)s/%(file)s on the GRID"
                                                               % { "outputPathGRID": self._histogramPathGRID, "file": outputFilename } )
                else :
                    raise HistogramFileAlreadyOnGRIDError( "%(outputPathGRID)s/%(file)s on the GRID"
                                                           % { "outputPathGRID": self._histogramPathGRID, "file": outputFilename } )


    ## Execute all GRID
    #
    def executeSingleRunAllGRID( self, index, runString ) :

        # do preliminary checks
        self.doPreliminaryTest( index, runString )

        # now prepare the jdl file from the template
        self.generateJDLFile( index, runString )

        # now generate the run job script
        self.generateRunjobFile( index, runString )

        # don't forget to generate the steering file!
        self._steeringFileName = self.generateSteeringFileSingleRun( runString )

        # finally ready to submit! Let's do it...!
        self.submitJDL( index, runString )

        # prepare a tarball with all the ancillaries
        self.prepareTarball( runString )

        # cleaning up the system
        self.cleanup( runString )




    ## Generate the run job
    #
    # This method is used to generate the run job script
    #
    def generateRunjobFile( self, index, runString ):
        message = "Generating the executable (%(name)s-%(run)s.sh)" % { "name": self.name, "run": runString }
        self._logger.info( message )
        try :
            runTemplate = self._configParser.get( "SteeringTemplate", "FilterGRIDScript" )
        except ConfigParser.NoOptionError:
            message = "Unable to find a valid executable template"
            self._logger.critical( message )
            raise StopExecutionError( message )

        runTemplateString = open( runTemplate, "r" ).read()
        runActualString = runTemplateString

        # replace the runString
        runActualString = runActualString.replace( "@RunString@", runString )

        # replace the pedeString as well
        runActualString = runActualString.replace( "@PedeString@", self._pedeString )

        # replace the job name
        runActualString = runActualString.replace( "@Name@", self.name )

        variableList = [ "GRIDCE", "GRIDSE", "GRIDStoreProtocol", "GRIDVO",
                         "GRIDFolderBase", "GRIDFolderClusearchResults", "GRIDFolderDBPede", "GRIDFolderFilterResults",
                         "GRIDFolderFilterJoboutput", "GRIDFolderFilterHisto", "GRIDLibraryTarball", "GRIDILCSoftVersion" ]
        for variable in variableList:
            try:
                value = self._configParser.get( "GRID", variable )
                if variable == "GRIDCE":
                    self._gridCE = value
                runActualString = runActualString.replace( "@%(value)s@" % {"value":variable} , value )
            except ConfigParser.NoOptionError:
                message = "Unable to find variable %(var)s in the config file" % { "var" : variable }
                self._logger.critical( message )
                raise StopExecutionError( message )

        # prepare the run list
        if not self._option.merge:
            runActualString = runActualString.replace( "@RunList@", runString)
            runActualString = runActualString.replace( "@Merge@", "no" )
            runActualString = runActualString.replace( "@Base@", runString )
        else :
            for run in self._runStringList :
                runActualString = runActualString.replace( "@RunList@", "%(run)s @RunList@" % { "run": run } )
            runActualString = runActualString.replace( "@RunList@", "" )
            runActualString = runActualString.replace( "@Merge@", "yes" )
            runActualString = runActualString.replace( "@Base@", self._option.output )

        self._runScriptFilename = "%(name)s-%(run)s.sh" % { "name": self.name, "run": runString }
        runActualFile = open( self._runScriptFilename, "w" )
        runActualFile.write( runActualString )
        runActualFile.close()

        # change the mode of the run script to 0777
        os.chmod(self._runScriptFilename, 0777)


    ## Do the real submission
    #
    # This is doing the job submission
    #
    def submitJDL( self, index, runString ) :
        self._logger.info("Submitting the job to the GRID")
        command = "glite-wms-job-submit -a -r %(GRIDCE)s -o %(name)s-%(run)s.jid %(name)s-%(run)s.jdl" % {
            "name": self.name, "run": runString , "GRIDCE":self._gridCE }
        glite = popen2.Popen4( command )
        while glite.poll() == -1:
            message = glite.fromchild.readline().strip()
            self._logger.log(15, message )

        if glite.poll() == 0:
            self._logger.info( "Job successfully submitted to the GRID" )
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, "Submit'd", d, e, f
        else :
            raise GRIDSubmissionError ( "%(name)s-%(run)s.jdl" % { "name": self.name,"run": runString } )

        # read back the the job id file
        # this is made by two lines only the first is comment, the second
        # is what we are interested in !
        jidFile = open( "%(name)s-%(run)s.jid" % { "name": self.name, "run": runString } )
        jidFile.readline()
        run, b = self._gridJobNTuple[ index ]
        self._gridJobNTuple[ index ] = run, jidFile.readline()
        jidFile.close()


    ## Check the existence of the pedestal file
    # 
    # Differently from the check input and output methods, this is done
    # once only at the very beginning before entering in the loop on runs.
    #
    def checkPedestalFile( self ) :
        self._logger.info( "Checking the pedestal file (run%(file)s-ped-db.slcio)" % { "file": self._pedeString } )

        # the pedestal file should be something like this
        # run123456-ped-db.slcio
        self._pedeFilename = "run%(run)s-ped-db.slcio" % { "run": self._pedeString }

        # where to check, depends from the execution mode! 

        if self._option.execution == "all-local":
            self.checkPedestalFileLocally()
        elif self._option.execution == "cpu-local" or self._option.execution == "all-grid" :
            self.checkPedestalFileGRID()


    ## Check locally for the pedestal file
    #
    def checkPedestalFileLocally( self ):

        # the pedestal file should be something like this
        # db/run123456-ped-db.slcio

        try :
            localPath = self._configParser.get( "LOCAL", "LocalFolderDBPede" )
        except ConfigParser.NoOptionError :
            localPath = "db"

        if not os.access( os.path.join( localPath, self._pedeFilename ), os.R_OK ):
            message = "Missing pedestal file %(pede)s" % { "pede": self._pedeFilename }
            self._logger.critical ( message )
            raise MissingPedestalFileError( self._pedeFilename )


    ## Check on GRID for the pedestal file
    #
    def checkPedestalFileGRID( self ):

        try :
            gridPath = self._configParser.get( "GRID" , "GRIDFolderDBPede" )
        except ConfigParser.NoOptionError :
            message = "Unable to find the GRIDFolderDBPede."
            self._logger.critical( message )
            raise StopExecutionError( message )

        try :
            localPath = self._configParser.get( "LOCAL", "LocalFolderDBPede" )
        except ConfigParser.NoOptionError :
            localPath = "db"

        command = "lfc-ls %(gridPath)s/%(file)s > /dev/null 2>&1 " % {
            "gridPath": gridPath, "file": self._pedeFilename }
        if os.system( command ) != 0:
            message = "Missing pedestal file %(file)s" %{ "file": self._pedeFilename }
            self._logger.critical ( message )
            raise MissingPedestalFileOnGRIDError( self._pedeFilename )

        if self._option.execution == "cpu-local":

            # get the file then
            baseCommand = "lcg-cp "
            if self._option.verbose :
                baseCommand = baseCommand + " -v "
                command = baseCommand + "  lfn:%(gridPath)s/%(file)s file:%(localPath)s/%(file)s " % {
                    "gridPath" : gridPath, "file": self._pedeFilename, "localPath": localPath }

            self._logger.info( "Getting the pedestal file %(file)s" % { "file": self._pedeFilename } )
            if os.system( command ) != 0:
                message = "Problem getting the pedestal file from the GRID"
                self._logger.critical( message  )
                raise StopExecutionError( message )
