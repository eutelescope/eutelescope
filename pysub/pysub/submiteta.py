import os
import sys
import shutil
import sha
import glob
import tarfile
import popen2
import commands
import ConfigParser
import logging
import logging.handlers
import math
import datetime
from submitbase import SubmitBase
from error import *

## Submit eta calculation jobs
#
# This calss is responsible to submit jobs for eta function calculation.
# It is inheriting from SubmitBase and it is called by the submit-eta.py script
#
#
# @version $Id: submiteta.py,v 1.14 2009-07-28 09:24:41 bulgheroni Exp $
# @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
#
class SubmitEta( SubmitBase ):

    ## Version
    # This number is used when printing out the version of the package.
    #
    # Static member.
    #
    cvsVersion = "$Revision: 1.14 $"

    ## Name
    # This is the namer of the class. It is used in flagging all the log entries
    # and preparing all the file
    #
    # Static member.
    #
    name = "eta"

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

        # issue a critical error in case the user didn't provide any output base name
        if self._option.output == None:
            self._logger.critical( "No output base name provided. Please use -o option to specify it" )
            sys.exit( 5 )

        # issua a critical error also in case there are more than one input files, because 
        # the current ETA calculation is able to handle one input file only.
        if len( self._args ) != 1 :
            self._logger.critical( "Current ETA calculation can handle one and one only input file!" )
            sys.exit( 6 )

    ## Execute method
    #
    # This is the real method, responsible for the job submission
    # Actually this is just a sort of big switch calling the real submitter
    # depending of the execution mode
    #
    def execute( self ) :

        # there is one and one only input runs.
        # so no reason to loop over _args.

        # the ntuple contains 6 variables:
        # runNumber ; inputFileStatus ; marlinStatus ; outputFileStatus ; histoFileStatus ; joboutputFileStatus
        entry = self._option.output, "Unknown", "Unknown", "Unknown", "Unknown", "Unknown"
        self._summaryNTuple.append( entry )

        # do the same for the GRID NTuple
        entry = self._option.output , "Unknown"
        self._gridJobNTuple.append( entry )

        try:

            # now do something different depending on the execution option
            if self._option.execution == "all-grid" :
                self.executeAllGRID( 0, self._option.output )

            elif self._option.execution == "all-local" :
                self.executeAllLocal( 0 , self._option.output )

            elif self._option.execution == "cpu-local":
                self.executeCPULocal( 0 , self._option.output )

            elif self._option.execution == "only-generate":
                self.executeOnlyGenerate( 0 , self._option.output )

        except MissingInputFileError, error:
            message = "Missing input file %(file)s" % { "file": error._filename }
            self._logger.critical( message )
            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, "Missing", "Skipped", output, histo, tarball
            raise StopExecutionError( message )

        except MissingInputFileOnGRIDError, error:
            message = "Missing input file %(file)s" % { "file": error._filename }
            self._logger.error( message )
            self._logger.error("Skipping to the next run ")
            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, "Missing", "Skipped", output, histo, tarball
            raise StopExecutionError( message )

        except OutputFileAlreadyOnGRIDError, error:
            message = "Output file %(file)s already on GRID" % { "file": error._filename }
            self._logger.error( message )
            self._logger.error("Skipping to the next run ")
            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, "Skipped", "GRID", histo, tarball
            raise StopExecutionError( message )

        except HistogramFileAlreadyOnGRIDError, error:
            message = "Histogram file %(file)s already on GRID" % { "file": error._filename }
            self._logger.error( message )
            self._logger.error("Skipping to the next run ")
            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, "Skipped", output, "GRID", tarball
            raise StopExecutionError( message )

        except JoboutputFileAlreadyOnGRIDError, error:
            message = "Joboutput file %(file)s already on GRID" % { "file": error._filename }
            self._logger.error( message )
            self._logger.error("Skipping to the next run ")
            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, "Skipped", output, histo, "GRID"
            raise StopExecutionError( message )

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
            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, "Failed", "Missing", "Missing", "Missing"
            raise StopExecutionError( message )

        except MissingOutputFileError, error:
            message = "The output file (%(file)s) was not properly generated, possible failure" % { "file": error._filename }
            self._logger.error( message )
            run, b, c, d, e, f = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, b, c, "Missing", d, f
            raise StopExecutionError( message )

        except MissingHistogramFileError, error:
            message = "The output file (%(file)s) was not properly generated, possible failure" % { "file": error._filename }
            self._logger.error( message )
            run, b, c, d, e, f = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, b, c, d, "Missing", f
            raise StopExecutionError( message )

        except MissingJoboutputFileError, error:
            message = "The joboutput tarball (%(file)s) is missing, possible failure" % { "file": error._filename }
            self._logger.error( message )
            run, b, c, d, e, f = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, b, c, d, e, "Missing"
            raise StopExecutionError( message )

        except GRID_LCG_CPError, error:
            message = "Problem copying the input file (%(file)s)" % { "file": error._filename }
            self._logger.error( message )
            self._logger.error( "Skipping to the next run " )
            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, "Missing", marlin, output, histo, tarball
            raise StopExecutionError( message )

        except GRIDSubmissionError, error:
            message = "Problem submitting the job to the GRID (%(value)s)" % { "value": error._message }
            self._logger.critical( message )
            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, "Failed", marlin, output, histo, tarball
            raise StopExecutionError( message )

    ## Generate only submitter
    #
    # This methods is responsibile of dry-run with only steering file
    # generation
    #
    def executeOnlyGenerate( self, index , runString ):

        # just need to generate the steering file
        self.generateSteeringFile( runString )


    ## All Local submitter
    #
    # This methods represents the sequence of commands that should be done while
    # submitting jobs on the local computer.
    #
    def executeAllLocal( self, index , runString ):

        # before any futher, check we have the input file for this run
        # runString is not used in this case, since the input file is provided 
        # directly by the user via the command line
        self.checkInputFile( index, None )

        # first generate the steering file
        self._steeringFileName = self.generateSteeringFile( runString )

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
        # clusearch-123456.tar.gz
        self.checkJoboutputFile( index, runString )

        # clean up the local pc
        self.cleanup( runString )

    ## CPU Local submitter
    #
    # This methods represents the sequence of commands that should be done while
    # submitting jobs on the local computer but using remote data
    #
    def executeCPULocal( self, index , runString ):

        # do preliminary checks
        self.doPreliminaryTest( index, runString )

        # get the input file from the GRID
        self.getRunFromGRID( index, runString )

        # double check the presence of the input file
        self.checkInputFile( index, runString )

        #  generate the steering file
        self._steeringFileName = self.generateSteeringFile( runString )

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

    ## Get the input run from the GRID
    def getRunFromGRID(self, index, runString ) :

        self._logger.info(  "Getting the input file from the GRID" )

        try :
            inputPath = self._configParser.get( "GRID", "GRIDFolderFilterResults" )
        except ConfigParser.NoOptionError :
            message = "GRIDFolderFilterResults missing in the configuration file. Quitting."
            self._logger.critical( message )
            raise StopExecutionError( message )

        try :
            localPath = self._configParser.get( "LOCAL", "LocalFolderFilterResults" )
        except ConfigParser.NoOptionError :
            localPath = "$PWD/results"

        baseCommand = "lcg-cp "
        if self._option.verbose :
            baseCommand = baseCommand + " -v "

        # the input file name is given by the user via the command line.
        # it could be that the user provided file with a piece of path
        # attached. So remove it before proceeding...
        trash, file = os.path.split( self._args[0] )

        command = "%(base)s lfn:%(gridPath)s/%(file)s file:%(localPath)s/%(file)s" %  \
            { "base": baseCommand, "gridPath" : inputPath, "file": file, "localPath": localPath }
        if os.system( command ) != 0:
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, "Missing", c, d, e, f
            raise GRID_LCG_CPError( "lfn:%(gridPath)s/%(file)s" % \
                                        { "gridPath" : inputPath, "run": runString } )
        else:
            self._logger.info("Input file successfully copied from the GRID")
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, "OK", c, d, e, f


    ## Put the DB to the GRID
    def putOutputOnGRID( self, index, runString ):

        self._logger.info(  "Putting the output file to the GRID" )

        try :
            gridPath = self._configParser.get( "GRID", "GRIDFolderDBEta")
        except ConfigParser.NoOptionError :
            message = "GRIDFolderDBEta missing in the configuration file. Quitting."
            self._logger.critical( message )
            raise StopExecutionError( message )

        try :
            localPath = self._configParser.get( "LOCAL", "LocalFolderDBEta" )
        except ConfigParser.NoOptionError :
            localPath = "$PWD/db"

        baseCommand = "lcg-cr "
        if self._option.verbose :
            baseCommand = baseCommand + " -v "

        # the input file name is given by the user via the command line.
        # it could be that the user provided file with a piece of path
        # attached. So remove it before proceeding....
        file = "%(output)s-eta-db.slcio" % { "output" : self._option.output }

        command = "%(base)s -l lfn:%(gridFolder)s/%(file)s file:%(localFolder)s/%(file)s" % \
            { "base": baseCommand, "gridFolder": gridPath, "localFolder": localPath, "file" : file  }
        if os.system( command ) != 0 :
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, c, "LOCAL", e, f
            raise GRID_LCG_CRError( "lfn:%(gridFolder)s/%(file)s" % { "gridFolder": gridPath, "file" : file } )
        else:
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, c, "GRID", e, f
            self._logger.info( "Output file successfully copied to the GRID" )

        if self._option.verify_output:
            self._logger.info( "Verifying the output file integrity on the GRID" )
            filename = file
            localCopy = open( os.path.join( localPath, filename ) ).read()
            localCopyHash = sha.new( localCopy ).hexdigest() 
            self._logger.log( 15, "Local copy hash is %(hash)s" % { "hash" : localCopyHash } )

            # now copying back the just copied file.
            baseCommand = "lcg-cp "
            if self._option.verbose :
                baseCommand = baseCommand + " -v "

            filenametest = "%(output)s-eta-db-test.slcio" % { "output" : self._option.output }
            command = "%(base)s lfn:%(gridFolder)s/%(file)s file:%(localFolder)s/%(filetest)s" % \
                { "base": baseCommand, "gridFolder": gridPath, "localFolder": localPath, "file" : filename, "filetest":filenametest }
            if os.system( command ) != 0 :
                run, input, marlin, output, histogram, tarball = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, input, marlin, "GRID - Fail", histogram, tarball
                self._logger.error( "Problem with the verification!" )
                raise GRID_LCG_CRError( "lfn:%(gridFolder)s/%(run)s" % { "gridFolder": gridPath, "run" : filename } )

            remoteCopy = open( os.path.join( localPath, filenametest ) ).read()
            remoteCopyHash = sha.new( remoteCopy ).hexdigest()
            self._logger.log( 15, "Remote copy hash is %(hash)s" % { "hash" : remoteCopyHash } )

            if remoteCopyHash == localCopyHash:
                run, input, marlin, output, histogram, tarball = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, input, marlin, "GRID - Ver", histogram, tarball
                self._logger.info( "Verification successful" )
            else :
                run, input, marlin, output, histogram, tarball = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, input, marlin, "GRID - Fail", histogram, tarball
                self._logger.error( "Problem with the verification!" )

            os.remove( os.path.join( localPath, filenametest ) )

    ## Put the histograms file to the GRID
    def putHistogramOnGRID( self, index, runString ):

        self._logger.info( "Putting the histogram file to the GRID" )

        try:
            gridPath = self._configParser.get( "GRID", "GRIDFolderEtaHisto")
        except ConfigParser.NoOptionError :
            message = "GRIDFolderEtaHisto missing in the configuration file. Quitting."
            self._logger.critical( message )
            raise StopExecutionError( message )

        try :
            localPath = self._configParser.get( "LOCAL", "LocalFolderEtaHisto" )
        except ConfigParser.NoOptionError :
            localPath = "$PWD/histo"

        baseCommand = "lcg-cr "
        if self._option.verbose :
            baseCommand = baseCommand + " -v "

        command = "%(base)s -l lfn:%(gridFolder)s/%(run)s-eta-histo.root file:%(localFolder)s/%(run)s-eta-histo.root" % \
            { "base": baseCommand, "gridFolder": gridPath, "localFolder": localPath, "run" : runString }
        if os.system( command ) != 0 :
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, c, d, "LOCAL", f
            raise GRID_LCG_CRError( "lfn:%(gridFolder)s/%(run)s-eta-histo.root" %   { "gridFolder": gridPath, "run" : runString } )
        else:
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, c, d, "GRID", f
            self._logger.info( "Histogram file successfully copied to the GRID" )

        if self._option.verify_output:
            self._logger.info( "Verifying the histogram integrity on the GRID" )
            filename = "%(run)s-eta-histo.root"  % { "run" : runString }

            localCopy = open( os.path.join( localPath, filename ) ).read()
            localCopyHash = sha.new( localCopy ).hexdigest() 
            self._logger.log( 15, "Local copy hash is %(hash)s" % { "hash" : localCopyHash } )

            # now copying back the just copied file.
            baseCommand = "lcg-cp "
            if self._option.verbose :
                baseCommand = baseCommand + " -v "

            filenametest = "%(run)s-eta-histo-test.root"  % { "run" : runString }

            command = "%(base)s lfn:%(gridFolder)s/%(file)s file:%(localFolder)s/%(filetest)s" % \
                { "base": baseCommand, "gridFolder": gridPath, "localFolder": localPath, "filetest": filenametest, "file" : filename }
            if os.system( command ) != 0 : 
                run, input, marlin, output, histogram, tarball = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, input, marlin, output, "GRID - Fail", tarball
                self._logger.error( "Problem with the verification!" )
                raise GRID_LCG_CRError( "lfn:%(gridFolder)s/%(run)s" % { "gridFolder": gridPath, "run" : filename } )

            remoteCopy = open( os.path.join( localPath, filenametest ) ).read()
            remoteCopyHash = sha.new( remoteCopy ).hexdigest()
            self._logger.log( 15, "Remote copy hash is %(hash)s" % { "hash" : remoteCopyHash } )

            if remoteCopyHash == localCopyHash:
                run, input, marlin, output, histogram, tarball = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, input, marlin, output, "GRID - Ver", tarball
                self._logger.info( "Verification successful" )
            else :
                run, input, marlin, output, histogram, tarball = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, input, marlin, output, "GRID - Fail", tarball
                self._logger.error( "Problem with the verification!" )

            os.remove( os.path.join( localPath, filenametest ) )

    ## Put the joboutput file to the GRID
    def putJoboutputOnGRID( self, index, runString ):

        self._logger.info( "Putting the joboutput file to the GRID" )

        try:
            gridPath = self._configParser.get( "GRID", "GRIDFolderEtaJoboutput")
        except ConfigParser.NoOptionError :
            message = "GRIDFolderEtaJoboutput missing in the configuration file. Quitting."
            self._logger.critical( message )
            raise StopExecutionError( message )

        try :
            localPath = self._configParser.get( "LOCAL", "LocalFolderEtaJoboutput" )
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

        if self._option.verify_output:
            self._logger.info( "Verifying the joboutput integrity on the GRID" )
            filename = "%(name)s-%(run)s.tar.gz" % { "name": self.name, "run" : runString }
            localCopy = open( os.path.join( localPath, filename ) ).read()
            localCopyHash = sha.new( localCopy ).hexdigest()
            self._logger.log( 15, "Local copy hash is %(hash)s" % { "hash" : localCopyHash } )

            # now copying back the just copied file.
            baseCommand = "lcg-cp "
            if self._option.verbose :
                baseCommand = baseCommand + " -v "

            filenametest = "%(name)s-%(run)s-test.tar.gz" % { "name": self.name,  "run" : runString }
            command = "%(base)s lfn:%(gridFolder)s/%(file)s file:%(localFolder)s/%(filetest)s" % \
                { "base": baseCommand, "gridFolder": gridPath, "localFolder": localPath, "filetest": filenametest,"file" : filename }
            if os.system( command ) != 0 : 
                run, input, marlin, output, histogram, tarball = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, input, marlin, output, histogram, "GRID - Fail"
                self._logger.error( "Problem with the verification!" )
                raise GRID_LCG_CRError( "lfn:%(gridFolder)s/%(run)s" % { "gridFolder": gridPath, "run" : filename } )

            remoteCopy = open( os.path.join( localPath, filenametest ) ).read()
            remoteCopyHash = sha.new( remoteCopy ).hexdigest()
            self._logger.log( 15, "Remote copy hash is %(hash)s" % { "hash" : remoteCopyHash } )

            if remoteCopyHash == localCopyHash:
                run, input, marlin, output, histogram, tarball = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, input, marlin, output, histogram, "GRID - Ver"
                self._logger.info( "Verification successful" )
            else :
                run, input, marlin, output, histogram, tarball = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, input, marlin, output, histogram, "GRID - Fail"
                self._logger.error( "Problem with the verification!" )

            os.remove( os.path.join( localPath, filenametest ) )

    ## Check the input file
    #
    def checkInputFile( self, index, runString ) :
        # the input file should be something like:
        # results/%(file)s

        self._logger.info( "Checking the input file" )

        try :
            inputFilePath = self._configParser.get( "LOCAL", "LocalFolderFilterResults" )
        except ConfigParser.NoOptionError :
            inputFilePath = "results"

        # the input file name is the one provided by the user in the command line
        inputFileName = self._args[0]
        if self._option.execution == "cpu-local" :
            path = os.path.join( inputFilePath, inputFileName )
        else:
            path = os.path.join(inputFilePath, os.path.abspath( inputFileName ) )

        if not os.access( path  , os.R_OK ):
            message = "Problem accessing the input file (%(file)s), trying next run" % {"file": inputFileName }
            self._logger.error( message )
            raise MissingInputFileError( inputFileName )
        else :
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, "OK", c, d, e, f


    ## Generate the steering file
    def generateSteeringFile( self, runString  ) :
        message = "Generating the steering file (%(name)s-%(run)s.xml) " % { "name": self.name, "run" : runString }
        self._logger.info( message )

        steeringFileTemplate =  ""
        try:
            steeringFileTemplate = self._configParser.get( "SteeringTemplate", "EtaSteeringFile" )
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


        # now replace the input file name. This is a bit peculiar compared with the previous.
        # the full path will be assembled here and replaced in one go
        if self._option.execution == "all-grid" :
            folder = "results/"

            # check how the input file was given
            trash, file = os.path.split( self._args[0] )

            # the fully qualified input file name will be:
            self._fqInputFile = folder + file
        else :
            try:
                folder = self._configParser.get( "LOCAL", "LocalFolderFilterResults" )
            except ConfigParser.NoOptionError :
                folder = "results"

            # compose the fully qualified input file
            if self._option.execution == "cpu-local":
                self._fqInputFile = os.path.join( folder,  self._args[0]  )
            else:
                self._fqInputFile = os.path.join( folder, os.path.abspath( self._args[0] ) )

        actualSteeringString = actualSteeringString.replace("@InputFile@", self._fqInputFile )

        # now replace the eta folder path
        if self._option.execution == "all-grid" :
            dbFolder = "db"
        else:
            try:
                dbFolder = self._configParser.get("LOCAL", "LocalFolderDBEta")
            except ConfigParser.NoOptionError :
                dbFolder = "db"
        actualSteeringString = actualSteeringString.replace("@DBPath@" ,dbFolder )

        # now replace the histo folder path
        if self._option.execution == "all-grid" :
            histoFolder = "histo"
        else:
            try:
                histoFolder = self._configParser.get("LOCAL", "LocalFolderEtaHisto")
            except ConfigParser.NoOptionError :
                histoFolder = "histo"
        actualSteeringString = actualSteeringString.replace("@HistoPath@" ,histoFolder )

        # finally replace the run string !
        actualSteeringString = actualSteeringString.replace("@RunNumber@", self._option.output )

        # open the new steering file for writing
        steeringFileName = "%(name)s-%(run)s.xml" % { "name": self.name, "run" : runString }
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
            l = line.strip()
            if(len(l) !=0):
              print l
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
        # db/output-eta-db.slcio

        try :
            outputFilePath = self._configParser.get( "LOCAL", "LocalFolderDBEta" )
        except ConfigParser.NoOptionError :
            outputFilePath = "db"

        outputFileName = "%(run)s-eta-db.slcio" % { "run": runString, }
        if not os.access( os.path.join( outputFilePath , outputFileName) , os.R_OK ):
            raise MissingOutputFileError( outputFileName )
        else :
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, c, "OK", e, f

    ## Check the histo file
    #
    def checkHistogramFile( self, index, runString ) :
        # this should be named something like
        # histo/output-eta-histo.root

        try :
            histoFilePath = self._configParser.get( "LOCAL", "LocalFolderEtaHisto" )
        except ConfigParser.NoOptionError :
            histoFilePath = "histo"
        histoFileName = "%(run)s-eta-histo.root" % { "run": runString }
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
        try :
            gearPath = self._configParser.get( "LOCAL", "LocalFolderGear" )
        except ConfigParser.NoOptionError :
            gearPath = ""
        listOfFiles.append( os.path.join( gearPath, self._gear_file ) )
       
        try :
            histoinfoPath = self._configParser.get( "LOCAL", "LocalFolderHistoinfo" )
        except ConfigParser.NoOptionError :
            histoinfoPath = ""
        listOfFiles.append( self._configFile )

        # all files starting with clusearch-123456 in the local folder
        for file in glob.glob( "%(name)s-*.*" % {"name": self.name } ):
            message = "Adding %(file)s to the joboutput tarball" % { "file": file }
            self._logger.debug( message )
            listOfFiles.append( file )

        # the histogram file
        if self._option.execution != "all-grid":
            try :
                histoFilePath = self._configParser.get( "LOCAL", "LocalFolderEtaHisto" )
            except ConfigParser.NoOptionError :
                histoFilePath = "histo"
                listOfFiles.append( os.path.join( histoFilePath, "%(run)s-eta-histo.root" % { "run": runString } ) )

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
            localFolder = self._configParser.get( "LOCAL", "LocalFolderEtaJoboutput")
        except ConfigParser.NoOptionError :
            self._logger.debug( "LocalFolderClusearchJoboutput not available, using $PWD/log ")
            localFolder = "log/"

        shutil.move( self._tarballFileName, localFolder )


    ## Check the joboutput file
    #
    def checkJoboutputFile( self, index, runString) :
        # this should be named something like
        # log/eta-output.tar.gz

        try :
            outputFilePath = self._configParser.get( "LOCAL", "LocalFolderEtaJoboutput" )
        except ConfigParser.NoOptionError :
            outputFilePath = "log"
        tarballFileName = "%(name)s-%(run)s.tar.gz" % { "name": self.name, "run": runString }
        if not os.access( os.path.join( outputFilePath, tarballFileName), os.R_OK ):
            raise MissingJoboutputFileError( tarballFileName )
        else:
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, c, d, e, "OK"


    ## Cleanup after each calculation
    def cleanup( self, runString ):

        self._logger.info( "Cleaning up the local pc" )

        # remove the log file and the steering file
        for file in glob.glob( "%(name)s-*" % {"name": self.name}):
            os.remove( file )


        # remove the input and output file
        if self._keepInput == False:
            try :
                inputFilePath = self._configParser.get( "LOCAL", "LocalFolderFilterResults" )
            except ConfigParser.NoOptionError :
                inputFilePath = "results"

            # the input file is taken from the command line
            trash, inputFile = os.path.split( self._args[0] )
            os.remove( os.path.join( inputFilePath, inputFile ))

        if self._keepOutput == False :
            try :
                outputFilePath = self._configParser.get( "LOCAL", "LocalFolderDBEta" )
            except ConfigParser.NoOptionError :
                outputFilePath = "db"

            outputFile = "%(run)s-eta-db.slcio" % { "run" : runString }
            for file in glob.glob( os.path.join( outputFilePath , outputFile ) ):
                os.remove( file )

            try :
                histoFilePath = self._configParser.get( "LOCAL", "LocalFolderEtaHisto" )
            except ConfigParser.NoOptionError :
                histoFilePath = "histo"

            histoFile = "%(run)s-eta-histo.root" % { "run": runString }
            os.remove( os.path.join( histoFilePath, histoFile ) )


    def end( self ) :

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
            localPath = self._configParser.get( "LOCAL", "LocalFolderClusearchJoboutput")
        except ConfigParser.NoOptionError :
            localPath = "log/"

        jidFile = open( os.path.join( localPath, "%(name)s-%(date)s.jid" % { "name": self.name, "date": unique } ), "w" )
        currentJIDFile = open( "current.jid" , "a" )
        for run, jid in self._gridJobNTuple:
            if jid != "Unknown" and jid != "See below":
                jidFile.write( jid )
                currentJIDFile.write( "# %(name)s %(output)s %(unique)s\n" % {"name":self.name, "output": self._option.output, "unique": unique } ) 
                currentJIDFile.write( jid )
                
        jidFile.close()
        currentJIDFile.close()


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
    def doPreliminaryTest( self, index, runString ) :
        # first log the voms-proxy-info
        self._logger.info( "Logging the voms-proxy-info" )
        command = "voms-proxy-info -all"
        status, output = commands.getstatusoutput( command )
        for line in output.splitlines():
            self._logger.info( line.strip() )

        if status != 0:
            message = "Problem with the GRID_UI"
            self._logger.critical( message )
            raise StopExecutionError( message )

        # also check that the proxy is still valid
        command = "voms-proxy-info -e"
        status, output = commands.getstatusoutput( command )
        if status != 0:
            message = "Expired proxy"
            self._logger.critical( message )
            raise StopExecutionError( message )
        else:
            self._logger.info( "Valid proxy found" )

        # check if we have already a delegation proxy
        self.checkDelegationProxy()

        # get all the needed path from the configuration file
        try :
            self._inputPathGRID     = self._configParser.get("GRID", "GRIDFolderFilterResults")
            self._outputPathGRID    = self._configParser.get("GRID", "GRIDFolderDBEta" )
            self._joboutputPathGRID = self._configParser.get("GRID", "GRIDFolderEtaJoboutput")
            self._histogramPathGRID = self._configParser.get("GRID", "GRIDFolderEtaHisto")
            folderList =  [ self._outputPathGRID, self._joboutputPathGRID, self._histogramPathGRID ]
        except ConfigParser.NoOptionError:
            message = "Missing path from the configuration file"
            self._logger.critical( message )
            raise StopExecutionError( message )

        # check if the input file is on the GRID,
        # the input file name is given by the user via the command line.
        # it could be that the user provided file with a piece of path
        # attached. So remove it before proceeding...
        trash, file = os.path.split( self._args[0] )
        command = "lfc-ls %(inputPathGRID)s/%(file)s" % { "inputPathGRID" : self._inputPathGRID,  "file": file }

        lfc = popen2.Popen4( command )
        while lfc.poll() == -1:
            pass
        if lfc.poll() == 0:
            self._logger.info( "Input file (%(file)s) found on the SE" % { "file": file } )
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, "GRID", c, d, e, f
        else:
            self._logger.error( "Input file (%(file)s) NOT found on the SE. Trying next run" % { "file": file } )
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, "Missing", c, d, e, f
            raise MissingInputFileOnGRIDError( "%(inputPathGRID)s/%(file)s" % { "inputPathGRID" : self._inputPathGRID,  "file": file } )

        # check the existence of the folders
        try :
            for folder in folderList:
                self.checkGRIDFolder( folder )

        except MissingGRIDFolderError, error :
            message = "Folder %(folder)s is unavailable. Quitting" % { "folder": error._filename }
            self._logger.critical( message )
            raise StopExecutionError( message )


        # check if the output file already exists
        command = "lfc-ls %(outputPathGRID)s/%(run)s-eta-db.slcio" % { "outputPathGRID": self._outputPathGRID, "run": runString }
        lfc = popen2.Popen4( command )
        while lfc.poll() == -1:
            pass
        if lfc.poll() == 0:
            self._logger.warning( "Output file %(outputPathGRID)s/%(run)s-eta-db.slcio already exists"
                                  % { "outputPathGRID": self._outputPathGRID, "run": runString } )
            if self._configParser.get("General","Interactive" ):
                if self.askYesNo( "Would you like to remove it?  [y/n] " ):
                    self._logger.info( "User decided to remove %(outputPathGRID)s/%(run)s-eta-db.slcio from the GRID"
                                       % { "outputPathGRID": self._outputPathGRID, "run": runString } )
                    command = "lcg-del -a lfn:%(outputPathGRID)s/%(run)s-eta-db.slcio" % { "outputPathGRID": self._outputPathGRID,"run": runString }
                    os.system( command )
                else :
                    raise OutputFileAlreadyOnGRIDError( "%(outputPathGRID)s/%(run)s-eta-db.slcio on the GRID"
                                                  % { "outputPathGRID": self._outputPathGRID, "run": runString } )
            else :
                raise OutputAlreadyOnGRIDError( "%(outputPathGRID)s/%(run)s-eta-db.slcio on the GRID"
                                              % { "outputPathGRID": self._outputPathGRID, "run": runString } )

        # check if the job output file already exists
        command = "lfc-ls %(outputPathGRID)s/%(name)s-%(run)s.tar.gz" % { 
            "name": self.name, "outputPathGRID": self._joboutputPathGRID, "run": runString }
        lfc = popen2.Popen4( command )
        while lfc.poll() == -1:
            pass
        if lfc.poll() == 0:
            self._logger.warning( "Joboutput file %(outputPathGRID)s/%(name)s-%(run)s.tar.gz already exists"
                                  % { "name": self.name, "outputPathGRID": self._joboutputPathGRID, "run": runString } )
            if self._configParser.get("General","Interactive" ):
                if self.askYesNo( "Would you like to remove it?  [y/n] " ):
                    self._logger.info( "User decided to remove %(outputPathGRID)s/%(name)s-%(run)s.tar.gz from the GRID"
                                       % { "name": self.name, "outputPathGRID": self._joboutputPathGRID, "run": runString } )
                    command = "lcg-del -a lfn:%(outputPathGRID)s/%(name)s-%(run)s.tar.gz" % { 
                        "name": self.name, "outputPathGRID": self._joboutputPathGRID, "run": runString }
                    os.system( command )
                else :
                    raise JoboutputFileAlreadyOnGRIDError( "%(outputPathGRID)s/%(name)s-%(run)s.tar.gz on the GRID"
                                                  % { "name": self.name, "outputPathGRID": self._joboutputPathGRID, "run": runString } )
            else :
                raise JoboutputFileAlreadyOnGRIDError( "%(outputPathGRID)s/%(name)s-%(run)s.tar.gz on the GRID"
                                              % { "name": self.name, "outputPathGRID": self._joboutputPathGRID, "run": runString } )


        # check if the histogram file already exists
        command = "lfc-ls %(outputPathGRID)s/%(run)s-eta-histo.root" % { "outputPathGRID": self._histogramPathGRID, "run": runString }
        lfc = popen2.Popen4( command )
        while lfc.poll() == -1:
            pass
        if lfc.poll() == 0:
            self._logger.warning( "Histogram file %(outputPathGRID)s/%(run)s-eta-histo.root already exists"
                                  % { "outputPathGRID": self._histogramPathGRID, "run": runString } )
            if self._configParser.get("General","Interactive" ):
                if self.askYesNo( "Would you like to remove it?  [y/n] " ):
                    self._logger.info( "User decided to remove %(outputPathGRID)s/%(run)s-eta-histo.root from the GRID"
                                       % { "outputPathGRID": self._histogramPathGRID, "run": runString } )
                    command = "lcg-del -a lfn:%(outputPathGRID)s/%(run)s-eta-histo.root" % { "outputPathGRID": self._histogramPathGRID, "run": runString }
                    os.system( command )
                else :
                    raise HistogramFileAlreadyOnGRIDError( "%(outputPathGRID)s/run%(run)s-clu-histo.root on the GRID"
                                                           % { "outputPathGRID": self._histogramPathGRID, "run": runString } )
            else :
                raise HistogramFileAlreadyOnGRIDError( "%(outputPathGRID)s/%(run)s-eta-histo.root on the GRID"
                                                       % { "outputPathGRID": self._histogramPathGRID, "run": runString } )
    ## Execute all GRID
    #
    def executeAllGRID( self, index, runString ) :

        # do preliminary checks
        self.doPreliminaryTest( index, runString )

        # now prepare the jdl file from the template
        self.generateJDLFile( index, runString )

        # now generate the run job script
        self.generateRunjobFile( index, runString )

        # don't forget to generate the steering file!
        self._steeringFileName = self.generateSteeringFile( runString )

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
            runTemplate = self._configParser.get( "SteeringTemplate", "EtaGRIDScript" )
        except ConfigParser.NoOptionError:
            message = "Unable to find a valid executable template"
            self._logger.critical( message )
            raise StopExecutionError( message )

        runTemplateString = open( runTemplate, "r" ).read()
        runActualString = runTemplateString

        # replace the input file name. Again strip away any possible path from this.
        trash, filename = os.path.split( self._args[0] )
        runActualString = runActualString.replace ("@Filename@", filename )

        # replace the output base
        runActualString = runActualString.replace( "@OutputBase@", self._option.output )

        # replace the job name
        runActualString = runActualString.replace( "@Name@", self.name )

        variableList = [ "GRIDCE", "GRIDSE", "GRIDStoreProtocol", "GRIDVO",
                         "GRIDFolderBase", "GRIDFolderFilterResults", "GRIDFolderDBEta", "GRIDFolderEtaJoboutput",
                         "GRIDFolderEtaJoboutput", "GRIDFolderEtaHisto", "GRIDLibraryTarball",
                         "GRIDLibraryTarballPath", "GRIDILCSoftVersion" ]
        for variable in variableList:
            try:
                value = self._configParser.get( "GRID", variable )
                if variable == "GRIDCE":
                    self._gridCE = value
                if variable == "GRIDLibraryTarballPath":
                    if  value.startswith( "lfn:" ) :
                        runActualString = runActualString.replace( "@HasLocalGRIDLibraryTarball@", "no" )
                        value = value.lstrip("lfn:")
                    else:
                        runActualString = runActualString.replace( "@HasLocalGRIDLibraryTarball@", "yes" )

                runActualString = runActualString.replace( "@%(value)s@" % {"value":variable} , value )

            except ConfigParser.NoOptionError:
                message = "Unable to find variable %(var)s in the config file" % { "var" : variable }
                self._logger.critical( message )
                raise StopExecutionError( message )

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
        command = "glite-wms-job-submit %(del)s -r %(GRIDCE)s -o %(name)s-%(run)s.jid %(name)s-%(run)s.jdl" % {
            "name": self.name, "run": runString , "GRIDCE":self._gridCE, "del": self._jobDelegation }
        status, output = commands.getstatusoutput( command )
        for line in output.splitlines():
            self._logger.log( 15, line.strip() )
        if status == 0:
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
