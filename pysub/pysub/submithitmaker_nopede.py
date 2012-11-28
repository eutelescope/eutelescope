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

## Submit hit making jobs
#
# This calss is responsible to submit jobs for hit making, i.e. going from the 
# detector frame of reference to the global telescope one.
# It is inheriting from SubmitBase and it is called by the submit-hitmaker.py script
#
#
# @version $Id: submithitmaker.py,v 1.14 2009-07-28 09:24:41 bulgheroni Exp $
# @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
#
class SubmitHitMaker( SubmitBase ):

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
    name = "hitmaker"

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
                self._logger.info("User decided to continue removing input files")
            else:
                self._logger.info("Aborted by user")
                sys.exit( 4 )

        if self._keepOutput == False and interactive :
            self._logger.warning("This script is going to delete the output file(s) when finished." )
            self._logger.warning("Are you sure to continue? [y/n]")
            message = "--> "
            if self.askYesNo(prompt = message ) == True:
                self._logger.info("User decided to continue removing output files")
            else:
                self._logger.info("Aborted by user")
                sys.exit( 4 )

        # issue a critical error in case the user didn't provide any output base name
        if self._option.output == None:
            self._logger.critical( "No output base name provided. Please use -o option to specify it" )
            sys.exit( 5 )

        # issue a warning in case the user didn't provide an eta file
        if self._option.eta == None:
            self._logger.warning("No ETA file provided. Please use -e option to specify it (if needed)" )
            self._option.eta = ""
#            sys.exit( 6 )


    ## Execute method
    #
    # This is the real method, responsible for the job submission
    # Actually this is just a sort of big switch calling the real submitter
    # depending of the execution mode
    #
    def execute( self ) :

        # we may have several input files provided as arguments through the command line
        # anyway it will be working in a sort of merging mode, it means that the steering
        # file will be generated having all the input files and one only output file
        # chaining all the hits
        self._inputFileList = []
        self._fqInputFileList = []
        self._justInputFileList = []
        for index, inputFile in enumerate( self._args ):

            # this is the list of basename
            self._justInputFileList.append(  os.path.basename( inputFile ) )

            # the ntuple contains 6 variables:
            # runNumber ; inputFileStatus ; marlinStatus ; outputFileStatus ; histoFileStatus ; joboutputFileStatus
            entry = os.path.basename( inputFile )[:6], "Unknown", "Unknown", "Unknown", "Unknown", "Unknown"
            self._summaryNTuple.append( entry )

            # do the same for the GRID NTuple
            entry = os.path.basename( inputFile )[:6] , "Unknown"
            self._gridJobNTuple.append( entry )

            self._inputFileList.append( inputFile ) 

            # now replace the input file name. This is a bit peculiar compared with the previous.
            # the full path will be assembled here and replaced in one go
            if self._option.execution == "all-grid" :
                folder = "results/"

                # check how the input file was given
                file = os.path.basename( inputFile )

                # the fully qualified input file name will be:
                self._fqInputFileList.append(folder + file)
            else :
                try:
                    folder = self._configParser.get( "LOCAL", "LocalFolderFilterResults" )
                except ConfigParser.NoOptionError :
                    folder = "results"

                if self._option.execution != "cpu-local":
                    self._fqInputFileList.append( os.path.join( folder,  inputFile  ) )
                else:
                    self._fqInputFileList.append( os.path.join( folder, os.path.abspath( inputFile ) ) )


        entry = self._option.output,  "Unknown", "Unknown", "Unknown", "Unknown", "Unknown"
        self._summaryNTuple.append( entry )

        entry = self._option.output,  "Unknown"
        self._gridJobNTuple.append( entry )

        try:

            # now do something different depending on the execution option
            if self._option.execution == "all-grid" :
                self.executeAllGRID(  )

            elif self._option.execution == "all-local" :
                self.executeAllLocal(  )

            elif self._option.execution == "cpu-local":
                self.executeCPULocal(  )

            elif self._option.execution == "only-generate":
                self.executeOnlyGenerate(  )

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
            self._logger.critical( message )
            raise StopExecutionError( message )

        except MissingOutputFileError, error:
            message = "The output file (%(file)s) was not properly generated, possible failure" % { "file": error._filename }
            self._logger.critical( message )
            raise StopExecutionError( message )

        except MissingHistogramFileError, error:
            message = "The histogram file (%(file)s) was not properly generated, possible failure" % { "file": error._filename }
            self._logger.error( message )
            raise StopExecutionError( message )

        except MissingJoboutputFileError, error:
            message = "The joboutput tarball (%(file)s) is missing, possible failure" % { "file": error._filename }
            self._logger.error( message )
            raise StopExecutionError( message )

        except OutputFileAlreadyOnGRIDError, error:
            message = "Output file %(file)s already on GRID" % { "file": error._filename }
            self._logger.critical( message )
            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len(self._summaryNTuple) - 1 ]
            self._summaryNTuple[ len(self._summaryNTuple) - 1 ] = run, input, "Skipped", "GRID", histo, tarball
            raise StopExecutionError( message )

        except HistogramFileAlreadyOnGRIDError, error:
            message = "Histogram file %(file)s already on GRID" % { "file": error._filename }
            self._logger.critical( message )
            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len(self._summaryNTuple) - 1 ]
            self._summaryNTuple[ len(self._summaryNTuple) - 1 ] = run, input, "Skipped", output, "GRID", tarball
            raise StopExecutionError( message )

        except JoboutputFileAlreadyOnGRIDError, error:
            message = "Joboutput file %(file)s already on GRID" % { "file": error._filename }
            self._logger.critical( message )
            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len(self._summaryNTuple) - 1 ]
            self._summaryNTuple[ len(self._summaryNTuple) - 1 ] = run, input, "Skipped", output, histo, "GRID"
            raise StopExecutionError( message )


        except GRIDSubmissionError, error:
            message ="Unable to perform the job submission (%(error)s) " % { "error": error._message }
            self._logger.critical( message )  
            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len(self._summaryNTuple) - 1 ]
            self._summaryNTuple[ len(self._summaryNTuple) - 1 ] = run, input, "Failed", output, histo, tarball
            raise StopExecutionError( message )

        except GRID_LCG_CPError, error:
            message = "Problem copying the input file (%(file)s)" % { "file": error._filename }
            self._logger.error( message )
            self._logger.error( "Skipping to the next run " )
            run, input, marlin, output, histo, tarball = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, "Missing", marlin, output, histo, tarball
            raise StopExecutionError( message )

    ## Generate only submitter
    #
    # This methods is responsibile of dry-run with only steering file
    # generation
    #
    def executeOnlyGenerate( self ):

        # just need to generate the steering file
        self.generateSteeringFile(  )


    ## All Local submitter
    #
    # This methods represents the sequence of commands that should be done while
    # submitting jobs on the local computer.
    #
    def executeAllLocal( self ):

        # before any futher, check we have the input file for this run
        # runString is not used in this case, since the input file is provided 
        # directly by the user via the command line
        self.checkInputFile(  )

        # first generate the steering file
        self._steeringFileName = self.generateSteeringFile( )

        # prepare a separate file for logging the output of Marlin
        self._logFileName = "%(name)s-%(run)s.log" % { "name": self.name, "run" : self._option.output }

        # run marlin
        self.runMarlin( )

        # advice the user that Marlin is over
        self._logger.info( "Marlin finished successfully")

        # verify the presence of the output files
        self.checkOutputFile(  )

        # verify the presence of the histogram files
        self.checkHistogramFile(  )

        # prepare a tarbal for the records
        self.prepareTarball(  )

        # check the presence of the joboutput tarball
        # this should be named something like
        # clusearch-123456.tar.gz
        self.checkJoboutputFile( )

        # clean up the local pc
        self.cleanup(  )

    ## CPU Local submitter
    #
    # This methods represents the sequence of commands that should be done while
    # submitting jobs on the local computer but using remote data
    #
    def executeCPULocal( self ):

        # do preliminary checks
        self.doPreliminaryTest(  )

        # get the input files from the GRID
        self.getRunFromGRID( )

        # double check the presence of the input files
        self.checkInputFile( )

        #  generate the steering file
        self._steeringFileName = self.generateSteeringFile(  )

        # prepare a separate file for logging the output of Marlin
        self._logFileName = "%(name)s-%(run)s.log" % { "name": self.name,"run" : self._option.output }

        # run marlin
        self.runMarlin( )

        # advice the user that Marlin is over
        self._logger.info( "Marlin finished successfully")

        # verify the presence of the output files
        self.checkOutputFile(  )

        # prepare a tarbal for the records
        self.prepareTarball(  )

        try :
            # copy the DB file to the GRID
            self.putOutputOnGRID( )

            # copy the histogram file to the GRID
            self.putHistogramOnGRID(  )

            # copy the joboutput to the GRID
            self.putJoboutputOnGRID(  )

            # clean up the local pc
            self.cleanup(  )

        except GRID_LCG_CRError, error:
            message = "The file (%(file)s) couldn't be copied on the GRID"  % { "file": error._filename }
            self._logger.critical( message )
            raise StopExecutionError( message )

    ## Get the input run from the GRID
    def getRunFromGRID( self ) :

        self._logger.info( "Getting the Eta file from the GRID" )
        try :
            inputPath = self._configParser.get( "GRID", "GRIDFolderDBEta" )
        except ConfigParser.NoOptionError :
            message = "GRIDFolderDBEta missing in the configuration file. Quitting."
            self._logger.critical( message )
            raise StopExecutionError( message )

        try :
            localPath = self._configParser.get( "LOCAL", "LocalFolderDBEta" )
        except ConfigParser.NoOptionError :
            localPath = "$PWD/db"

        baseCommand = "lcg-cp "
        if self._option.verbose :
            baseCommand = baseCommand + " -v "

        command = "%(base)s lfn:%(gridPath)s/%(file)s file:%(localPath)s/%(file)s" %  \
                      { "base": baseCommand, "gridPath" : inputPath, "file": os.path.basename( self._option.eta ), "localPath": localPath }

        if os.system( command ) != 0:
            message = "Unable to get the Eta file. Terminating!"
            self._logger.critical( message )
#            raise StopExecutionError( message )
        else:
            self._logger.info( "Eta file successfully copied from the GRID" )

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

        for index, inputFile in enumerate( self._inputFileList) :
            if inputFile != "DEADFACE" :
                command = "%(base)s lfn:%(gridPath)s/%(file)s file:%(localPath)s/%(file)s" %  \
                    { "base": baseCommand, "gridPath" : inputPath, "file": self._justInputFileList[ index ], "localPath": localPath }
                if os.system( command ) != 0:
                    run, b, c, d, e, f = self._summaryNTuple[ index ]
                    self._summaryNTuple[ index ] = run, "Missing", c, d, e, f
                    self._logger.error( "Error copying file %(file)s, skipping it" % { "file": self._justInputFileList[ index ] } )
                    self._inputFileList[ index ] = "DEADFACE" 
                else:
                    self._logger.info("Input file %(file)s successfully copied from the GRID" % { "file": self._justInputFileList[ index ] } )
                    run, b, c, d, e, f = self._summaryNTuple[ index ]
                    self._summaryNTuple[ index ] = run, "OK", c, d, e, f

        run, b, c, d, e, f = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
        self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, "See above", c, d, e, f


    ## Put the DB to the GRID
    def putOutputOnGRID( self ):

        self._logger.info(  "Putting the output file to the GRID" )

        try :
            gridPath = self._configParser.get( "GRID", "GRIDFolderHitmakerResults")
        except ConfigParser.NoOptionError :
            message = "GRIDFolderHitmakerResults missing in the configuration file. Quitting."
            self._logger.critical( message )
            raise StopExecutionError( message )

        try :
            localPath = self._configParser.get( "LOCAL", "LocalFolderHitmakerResults" )
        except ConfigParser.NoOptionError :
            localPath = "$PWD/results"

        baseCommand = "lcg-cr "
        if self._option.verbose :
            baseCommand = baseCommand + " -v "

        # the input file name is given by the user via the command line.
        # it could be that the user provided file with a piece of path
        # attached. So remove it before proceeding....
        file = "%(output)s-hit.slcio" % { "output" : self._option.output }

        command = "%(base)s -l lfn:%(gridFolder)s/%(file)s file:%(localFolder)s/%(file)s" % \
            { "base": baseCommand, "gridFolder": gridPath, "localFolder": localPath, "file" : file  }
        if os.system( command ) != 0 :
            self._logger.critical( "Problem copying the output file %(file)s to the GRID" % {"file" : file} )
            for index, entry in enumerate ( self._summaryNTuple ):
                run, input, marlin, output, histo, tarball = entry
                self._summaryNTuple[ index ] = run, input, marlin, "See below", histo, tarball

            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, marlin, "Local", histo, tarball
            raise GRID_LCG_CRError( "lfn:%(gridFolder)s/%(file)s" % { "gridFolder": gridPath, "file" : file } )
        else:
            self._logger.info( "Output file successfully copied to the GRID" )
            for index, entry in enumerate( self._summaryNTuple ):
                run, input, marlin, output, histo, tarball = entry
                self._summaryNTuple[ index ] = run, input, marlin, "See below", histo, tarball

            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, marlin, "GRID", histo, tarball

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

            filenametest = "%(output)s-hit-test.slcio" % { "output" : self._option.output }
            command = "%(base)s lfn:%(gridFolder)s/%(file)s file:%(localFolder)s/%(filetest)s" % \
                { "base": baseCommand, "gridFolder": gridPath, "localFolder": localPath, "file" : filename, "filetest":filenametest }
            if os.system( command ) != 0 :
                run, input, marlin, output, histogram, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
                self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, marlin, "GRID - Fail", histogram, tarball
                self._logger.error( "Problem with the verification!" )
                raise GRID_LCG_CRError( "lfn:%(gridFolder)s/%(run)s" % { "gridFolder": gridPath, "run" : filename } )

            remoteCopy = open( os.path.join( localPath, filenametest ) ).read()
            remoteCopyHash = sha.new( remoteCopy ).hexdigest()
            self._logger.log( 15, "Remote copy hash is %(hash)s" % { "hash" : remoteCopyHash } )

            if remoteCopyHash == localCopyHash:
                run, input, marlin, output, histogram, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
                self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, marlin, "GRID - Ver", histogram, tarball
                self._logger.info( "Verification successful" )
            else :
                run, input, marlin, output, histogram, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
                self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, marlin, "GRID - Fail", histogram, tarball
                self._logger.error( "Problem with the verification!" )

            os.remove( os.path.join( localPath, filenametest ) )

    ## Put the histograms file to the GRID
    def putHistogramOnGRID( self ):

        self._logger.info( "Putting the histogram file to the GRID" )

        try:
            gridPath = self._configParser.get( "GRID", "GRIDFolderHitmakerHisto")
        except ConfigParser.NoOptionError :
            message = "GRIDFolderHitmakerHisto missing in the configuration file. Quitting."
            self._logger.critical( message )
            raise StopExecutionError( message )

        try :
            localPath = self._configParser.get( "LOCAL", "LocalFolderHitmakerHisto" )
        except ConfigParser.NoOptionError :
            localPath = "$PWD/histo"

        baseCommand = "lcg-cr "
        if self._option.verbose :
            baseCommand = baseCommand + " -v "

        command = "%(base)s -l lfn:%(gridFolder)s/%(run)s-hit-histo.root file:%(localFolder)s/%(run)s-hit-histo.root" % \
            { "base": baseCommand, "gridFolder": gridPath, "localFolder": localPath, "run" : self._option.output }
        if os.system( command ) != 0 :
            self._logger.critical( "Problem copying the histogram file %(run)s-hit-histo.root to the GRID" % {"run" : self._option.output} )
            for index, entry in enumerate( self._summaryNTuple ):
                run, input, marlin, output, histo, tarball = entry
                self._summaryNTuple[ index ] = run, input, marlin, output, "See below", tarball

            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, marlin, output, "Local", tarball
            raise GRID_LCG_CRError( "lfn:%(gridFolder)s/%(run)s-hit-histo.root" %   { "gridFolder": gridPath, "run" : self._option.output } )
        else:
            self._logger.info("Histogram file sucessfully copied to the GRID" )
            for index, entry in enumerate( self._summaryNTuple ):
                run, input, marlin, output, histo, tarball = entry
                self._summaryNTuple[ index ] = run, input, marlin, output, "See below", tarball

            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, marlin, output, "GRID", tarball

        if self._option.verify_output:
            self._logger.info( "Verifying the histogram integrity on the GRID" )
            filename = "%(run)s-hit-histo.root"  % { "run" : self._option.output }

            localCopy = open( os.path.join( localPath, filename ) ).read()
            localCopyHash = sha.new( localCopy ).hexdigest() 
            self._logger.log( 15, "Local copy hash is %(hash)s" % { "hash" : localCopyHash } )

            # now copying back the just copied file.
            baseCommand = "lcg-cp "
            if self._option.verbose :
                baseCommand = baseCommand + " -v "

            filenametest = "%(run)s-hit-histo-test.root"  % { "run" : self._option.output }

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
                run, input, marlin, output, histogram, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
                self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, marlin, output, "GRID - Ver", tarball
                self._logger.info( "Verification successful" )
            else :
                run, input, marlin, output, histogram, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
                self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, marlin, output, "GRID - Fail", tarball
                self._logger.error( "Problem with the verification!" )

            os.remove( os.path.join( localPath, filenametest ) )

    ## Put the joboutput file to the GRID
    def putJoboutputOnGRID( self ):

        self._logger.info( "Putting the joboutput file to the GRID" )

        try:
            gridPath = self._configParser.get( "GRID", "GRIDFolderHitmakerJoboutput")
        except ConfigParser.NoOptionError :
            message = "GRIDFolderHitmakerJoboutput missing in the configuration file. Quitting."
            self._logger.critical( message )
            raise StopExecutionError( message )

        try :
            localPath = self._configParser.get( "LOCAL", "LocalFolderHitmakerJoboutput" )
        except ConfigParser.NoOptionError :
            localPath = "$PWD/log"

        baseCommand = "lcg-cr "
        if self._option.verbose :
            baseCommand = baseCommand + " -v "

        command = "%(base)s -l lfn:%(gridFolder)s/%(name)s-%(run)s.tar.gz file:%(localFolder)s/%(name)s-%(run)s.tar.gz" % \
            { "name": self.name, "base": baseCommand, "gridFolder": gridPath, "localFolder": localPath, "run" : self._option.output }
        if os.system( command ) != 0 :
            self._logger.critical( "Problem copying the joboutput file %(name)s-%(run)s.tar.gz to the GRID" % {
                     "name": self.name, "run" : self._option.output } )
            for index, entry in enumerate( self._summaryNTuple ):
                run, input, marlin, output, histo, tarball = entry
                self._summaryNTuple[ index ] = run, input, marlin, output, histo, "See below"

            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, marlin, output, histo, "Local"
            raise GRID_LCG_CRError( "lfn:%(gridFolder)s/%(run)s-hit-histo.root" %   { "gridFolder": gridPath, "run" : self._option.output } )

        else:
            self._logger.info("Joboutput file sucessfully copied to the GRID" )
            for index, entry in enumerate( self._summaryNTuple ):
                run, input, marlin, output, histo, tarball = entry
                self._summaryNTuple[ index ] = run, input, marlin, output, histo, "See below"

            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, marlin, output, histo, "GRID"

        if self._option.verify_output:
            self._logger.info( "Verifying the joboutput integrity on the GRID" )
            filename = "%(name)s-%(run)s.tar.gz" % { "name": self.name, "run" : self._option.output  }
            localCopy = open( os.path.join( localPath, filename ) ).read()
            localCopyHash = sha.new( localCopy ).hexdigest()
            self._logger.log( 15, "Local copy hash is %(hash)s" % { "hash" : localCopyHash } )

            # now copying back the just copied file.
            baseCommand = "lcg-cp "
            if self._option.verbose :
                baseCommand = baseCommand + " -v "

            filenametest = "%(name)s-%(run)s-test.tar.gz" % { "name": self.name,  "run" : self._option.output  }
            command = "%(base)s lfn:%(gridFolder)s/%(file)s file:%(localFolder)s/%(filetest)s" % \
                { "base": baseCommand, "gridFolder": gridPath, "localFolder": localPath, "filetest": filenametest,"file" : filename }
            if os.system( command ) != 0 : 
                run, input, marlin, output, histogram, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
                self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, marlin, output, histogram, "GRID - Fail"
                self._logger.error( "Problem with the verification!" )
                raise GRID_LCG_CRError( "lfn:%(gridFolder)s/%(run)s" % { "gridFolder": gridPath, "run" : filename } )

            remoteCopy = open( os.path.join( localPath, filenametest ) ).read()
            remoteCopyHash = sha.new( remoteCopy ).hexdigest()
            self._logger.log( 15, "Remote copy hash is %(hash)s" % { "hash" : remoteCopyHash } )

            if remoteCopyHash == localCopyHash:
                run, input, marlin, output, histogram, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
                self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, marlin, output, histogram, "GRID - Ver"
                self._logger.info( "Verification successful" )
            else :
                run, input, marlin, output, histogram, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
                self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, marlin, output, histogram, "GRID - Fail"
                self._logger.error( "Problem with the verification!" )

            os.remove( os.path.join( localPath, filenametest ) )

    ## Check the input file
    #
    def checkInputFile( self ) :
        # the input files are coming from the command line

        try :
            inputFilePath = self._configParser.get( "LOCAL", "LocalFolderFilterResults" )
        except ConfigParser.NoOptionError :
            inputFilePath = "results"

        # loop over all the input files
        for index, inputFile in enumerate( self._inputFileList ):
            if inputFile == "DEADFACE":
                continue

            self._logger.info( "Checking the input file %(file)s" % { "file": self._justInputFileList[ index ] } )

            if not os.access( self._fqInputFileList[ index ] , os.R_OK ):
                self._logger.error( "Problem accessing the input file %(file)s, trying nex one" % { "file": self._fqInputFileList[ index ] } )
                self._inputFileList[ index ] = "DEADFACE"
                run, input, marlin, output, histo, tarball = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, "Missing", marlin, output, histo, tarball

            else :
                run, input, marlin, output, histo, tarball = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, "OK", marlin, output, histo, tarball

        run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
        self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, "See above", marlin, output, histo, tarball

    ## Generate the steering file
    def generateSteeringFile( self ) :
        message = "Generating the steering file (%(name)s-%(output)s.xml) " % { "name": self.name, "output" : self._option.output }
        self._logger.info( message )

        steeringFileTemplate =  ""
        try:
            steeringFileTemplate = self._configParser.get( "SteeringTemplate", "HitmakerSteeringFile" )
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

        # now replace the offset folder path
        if self._option.execution == "all-grid" :
#            try:
#                dbFolder = self._configParser.get("GRID", "GRIDFolderDBPede")
#            except ConfigParser.NoOptionError :
                dbFolder = "db"
        else:
            try:
                dbFolder = self._configParser.get("LOCAL", "LocalFolderDBPede")
            except ConfigParser.NoOptionError :
                dbFolder = "db"
        actualSteeringString = actualSteeringString.replace("@DBPath@" ,dbFolder )

        # set the run number (offsets to be taken from)
        # can be on run for many in the run block

        if self._option.offsetRunNumber == None :
            # this means that the user does not want to use hot pixel db file
            print "offset"
            print self._option.offsetRunNumber
        
            actualSteeringString = actualSteeringString.replace("@UseOffsetDBCommentLeft@", "!--" )
            actualSteeringString = actualSteeringString.replace("@UseOffsetDBCommentRight@", "--" )
 
#            actualSteeringString = actualSteeringString.replace("@OffsetRunNumber@", "000000" )
           
            self._logger.debug( "NOT Using offset db file" )
        else:    
            # this means that the user DOES want to use hot pixel db file
            print "offset YES"
            
            actualSteeringString = actualSteeringString.replace("@UseOffsetDBCommentLeft@", "" )
            actualSteeringString = actualSteeringString.replace("@UseOffsetDBCommentRight@", "" )
 
            offsetRunNumber = "0%(run)s" % { "run": self._option.offsetRunNumber }
            actualSteeringString = actualSteeringString.replace("@OffsetRunNumber@", offsetRunNumber )
     
            self._logger.debug( "Using offset db file" )
        
         

#        print actualSteeringString



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

        for index, inputFile in enumerate( self._inputFileList ):
            if inputFile != "DEADFACE":
                actualSteeringString = actualSteeringString.replace("@InputFile@", "%(fqfn)s @InputFile@" % { "fqfn": self._fqInputFileList[ index] } )

        actualSteeringString = actualSteeringString.replace( "@InputFile@", "" )

        # the eta file is provided by the user through the command line option -e. 
        # be careful in what you do!

        # now replace the eta folder path
        if self._option.execution == "all-grid" :
            etaFolder = "db"
            trash, etaFile = os.path.split( self._option.eta )
            fqEtaFile = os.path.join( etaFolder, etaFile )
        else:
            try:
                etaFolder = self._configParser.get("LOCAL", "LocalFolderDBEta")
            except ConfigParser.NoOptionError :
                etaFolder = "db"

            if self._option.execution != "cpu-local" :
                trash, etaFile = os.path.split( self._option.eta )
                fqEtaFile = os.path.join( etaFolder, etaFile )
            else :
                fqEtaFile = os.path.join( etaFolder, os.path.abspath( self._option.eta ) )

        actualSteeringString = actualSteeringString.replace("@EtaFile@" , fqEtaFile )

        # now replace the output folder path
        if self._option.execution == "all-grid" :
            outputFolder = "results"
        else:
            try :
                outputFolder =  self._configParser.get("LOCAL", "LocalFolderHitmakerResults")
            except ConfigParser.NoOptionError :
                outputFolder = "results"
        actualSteeringString = actualSteeringString.replace ("@ResultsFolder@", outputFolder )


        # now replace the histo folder path
        if self._option.execution == "all-grid" :
            histoFolder = "histo"
        else:
            try:
                histoFolder = self._configParser.get("LOCAL", "LocalFolderHitmakerHisto")
            except ConfigParser.NoOptionError :
                histoFolder = "histo"
        actualSteeringString = actualSteeringString.replace("@HistoPath@" ,histoFolder )

        # finally replace the run string !
        actualSteeringString = actualSteeringString.replace("@Output@", self._option.output )

        # replace the event range
#        print self._eventRange_begin
#        print self._eventRange_end
#        actualSteeringString = actualSteeringString.replace( "@FirstEvent@", self._eventRange_begin  )
#        actualSteeringString = actualSteeringString.replace( "@LastEvent@",  self._eventRange_end  )
 
        # open the new steering file for writing
        steeringFileName = "%(name)s-%(run)s.xml" % { "name": self.name, "run" : self._option.output }
        actualSteeringFile = open( steeringFileName, "w" )

        # write the new steering file
        actualSteeringFile.write( actualSteeringString )

        # close both files
        actualSteeringFile.close()
        templateSteeringFile.close()

        return steeringFileName

    ## Execute Marlin
    def runMarlin( self ) :

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
            for index, entry in enumerate( self._summaryNTuple ) :
                run, input, marlin, output, histo, tarball = entry
                self._summaryNTuple[ index ] = run, input, "Error", output, histo, tarball

            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, "Error", output, histo, tarball
            raise MarlinError( "", returnValue )
        else :
            for index, entry in enumerate( self._summaryNTuple ) :
                run, input, marlin, output, histo, tarball = entry
                self._summaryNTuple[ index ] = run, input, "OK", output, histo, tarball

            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, "OK", output, histo, tarball


    ## Check the output file
    #
    def checkOutputFile( self ) :

        try :
            outputFilePath = self._configParser.get( "LOCAL", "LocalFolderHitmakerResults" )
        except ConfigParser.NoOptionError :
            outputFilePath = "results"

        outputFileName = "%(run)s-hit.slcio" % { "run": self._option.output  }
        if not os.access( os.path.join( outputFilePath , outputFileName) , os.R_OK ):
            for index, entry in enumerate( self._summaryNTuple ):
                run, input, marlin, output, histo, tarball = entry
                self._summaryNTuple[ index ] = run, input, marlin, "Missing", histo, tarball
            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, marlin, "Missing", histo, tarball
            raise MissingOutputFileError( outputFileName )

        else :
            for index, entry in enumerate( self._summaryNTuple ) :
                run, input, marlin, output, histo, tarball = entry
                self._summaryNTuple[ index ] = run, input, marlin, "OK", histo, tarball
            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, marlin, "OK", histo, tarball

    ## Check the histo file
    #
    def checkHistogramFile( self ) :


        try :
            histoFilePath = self._configParser.get( "LOCAL", "LocalFolderHitmakerHisto" )
        except ConfigParser.NoOptionError :
            histoFilePath = "histo"

        histoFileName = "%(run)s-hit-histo.root" % { "run": self._option.output }
        if not os.access( os.path.join( histoFilePath , histoFileName ) , os.R_OK ):
            for index, entry in enumerate( self._summaryNTuple ):
                run, input, marlin, output, histo, tarball = entry
                self._summaryNTuple[ index ] = run, input, marlin, output, "Missing", tarball
            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, marlin, output, "Missing", tarball
            raise MissingHistogramFileError( histoFileName )
        else:
            for index, entry in enumerate( self._summaryNTuple ):
                run, input, marlin, output, histo, tarball = entry
                self._summaryNTuple[ index ] = run, input, marlin, output, "OK", tarball
            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, marlin, output, "OK", tarball

    ## Prepare the joboutput tarball
    def prepareTarball( self ) :
        self._logger.info("Preparing the joboutput tarball" )

        # first prepare a folder to store them temporary
        destFolder = "%(name)s-%(run)s" %{ "name": self.name, "run": self._option.output }
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
                histoFilePath = self._configParser.get( "LOCAL", "LocalFolderHitmakerHisto" )
            except ConfigParser.NoOptionError :
                histoFilePath = "histo"
            listOfFiles.append( os.path.join( histoFilePath, "%(run)s-hit-histo.root" % { "run": self._option.output } ) )

        # copy everything into a temporary folder
        for file in listOfFiles :
            shutil.copy( file, destFolder )

        # do the tarball
        self._tarballFileName = "%(name)s-%(run)s.tar.gz" % { "name": self.name, "run": self._option.output }
        tarball = tarfile.open( self._tarballFileName, "w:gz" )
        tarball.add( destFolder )
        tarball.close()

        # remove the temporary folder
        shutil.rmtree( destFolder )

        # copy the tarball in the log folder
        try:
            localFolder = self._configParser.get( "LOCAL", "LocalFolderHitmakerJoboutput")
        except ConfigParser.NoOptionError :
            self._logger.debug( "LocalFolderHitmakerJoboutput not available, using $PWD/log ")
            localFolder = "log/"

        dest = os.path.join( localFolder, self._tarballFileName )
        if os.path.exists(dest):
          os.remove(dest)

        shutil.move( self._tarballFileName, localFolder )


    ## Check the joboutput file
    #
    def checkJoboutputFile( self ) :

        try :
            outputFilePath = self._configParser.get( "LOCAL", "LocalFolderHitmakerJoboutput" )
        except ConfigParser.NoOptionError :
            outputFilePath = "log"
        tarballFileName = "%(name)s-%(run)s.tar.gz" % { "name": self.name, "run": self._option.output }

        if not os.access( os.path.join( outputFilePath, tarballFileName), os.R_OK ):
            for index, entry in enumerate( self._summaryNTuple ):
                run, input, marlin, output, histo, tarball = entry
                self._summaryNTuple[ index ] = run, input, marlin, output, histo, "Missing"
            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, marlin, output, histo, "Missing"
            raise MissingJoboutputFileError( tarballFileName )
        else:
            for index, entry in enumerate( self._summaryNTuple ) :
                run, input, marlin, output, histo, tarball = entry
                self._summaryNTuple[ index ] = run, input, marlin, output, histo, "OK"
            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) -1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, marlin, output, histo, "OK"


    ## Cleanup after each calculation
    def cleanup( self ):

        self._logger.info( "Cleaning up the local pc" )

        # remove the log file and the steering file
        for file in glob.glob( "%(name)s-*" % {"name": self.name}):
            os.remove( file )


        # remove the input and output file
        if self._keepInput == False:

            for index, file in enumerate( self._inputFileList ):
                if file != "DEADFACE" :
                    os.remove( self._fqInputFileList[ index ] )

        if self._keepOutput == False :
            try :
                outputFilePath = self._configParser.get( "LOCAL", "LocalFolderHitmakerResults" )
            except ConfigParser.NoOptionError :
                outputFilePath = "results"

            outputFile = "%(run)s-hit.slcio" % { "run" : self._option.output }
            for file in glob.glob( os.path.join( outputFilePath , outputFile ) ):
                os.remove( file )

            try :
                histoFilePath = self._configParser.get( "LOCAL", "LocalFolderHitmakerHisto" )
            except ConfigParser.NoOptionError :
                histoFilePath = "histo"

            histoFile = "%(run)s-hit-histo.root" % { "run": self._option.output }
            os.remove( os.path.join( histoFilePath, histoFile ) )


    def end( self ) :

        if self._keepInput == False :
            # remove also the eta file
            try:
                dbFilePath = self._configParser.get( "LOCAL", "LocalFolderDBEta" )
            except ConfigParser.NoOptionError :
                dbFilePath = "db"

            if self._option.execution != "cpu-local":
                etaFile = os.path.basename( self._option.eta )
            else:
                etaFile = os.path.abspath( self._option.eta)

            if os.access( os.path.join( dbFilePath , etaFile ), os.W_OK ):
                os.remove( os.path.join( dbFilePath , etaFile ) )


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
    def doPreliminaryTest( self ) :
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
            self._outputPathGRID    = self._configParser.get("GRID", "GRIDFolderHitmakerResults" )
            self._joboutputPathGRID = self._configParser.get("GRID", "GRIDFolderHitmakerJoboutput")
            self._histogramPathGRID = self._configParser.get("GRID", "GRIDFolderHitmakerHisto")
            self._etaPathGRID       = self._configParser.get("GRID", "GRIDFolderDBEta")
            self._offsetPathGRID    = self._configParser.get("GRID", "GRIDFolderDBOffset" )
            self._dbHotPixelPathGRID    = self._configParser.get("GRID", "GRIDFolderDBHotPixel" )
            folderList =  [ self._outputPathGRID, self._joboutputPathGRID, self._histogramPathGRID, self._offsetPathGRID, self._dbHotPixelPathGRID ]
        except ConfigParser.NoOptionError:
            message = "Missing path from the configuration file"
            self._logger.critical( message )
            raise StopExecutionError( message )

        # check if the eta file is on the GRID
        etaFile = os.path.basename( self._option.eta )
        command = "lfc-ls %(etaPathGRID)s/%(etaFile)s" % {"etaPathGRID" : self._etaPathGRID , "etaFile": etaFile }
#        print command
        lfc = popen2.Popen4( command )
        if lfc.poll() == 0:
            self._logger.info( "Eta file found on SE" )
        else:
            message = "Eta file NOT found on SE. Terminating"
            self._logger.critical( message )
#            raise StopExecutionError( message )

        # check if the input files is on the GRID
        for index, inputFile in enumerate( self._inputFileList ):
            if inputFile != "DEADFACE" :
                justFile = self._justInputFileList[ index ] 
                command = "lfc-ls %(inputPathGRID)s/%(file)s" % { "inputPathGRID" : self._inputPathGRID,  "file": justFile  }
                print command

                lfc = popen2.Popen4( command )
                while lfc.poll() == -1:
                    pass

                if lfc.poll() == 0:
                    self._logger.info( "Input file %(justFile)s found on the SE" % {"justFile": justFile } )
                    run, b, c, d, e, f = self._summaryNTuple[ index ]
                    self._summaryNTuple[ index ] = run, "GRID", c, d, e, f
                else:
                    self._logger.error( "Input file %(justFile)s NOT found on the SE" % {"justFile": justFile } )
                    run, b, c, d, e, f = self._summaryNTuple[ index ]
                    self._summaryNTuple[ index ] = run, "Missing", c, d, e, f
                    self._inputFileList[ index ] = "DEADFACE"
                    if self._configParser.get("General","Interactive" ):
                        if not self.askYesNo( "Would you like to skip it and continue? [y/n] " ) :
                            message = "User decided to stop here"
                            self._logger.critical( message )
                            raise StopExecutionError( message )
                        else :
                            self._logger.info( "Skipping to the next run" )
                    else:
                        self._logger.info( "Skipping to the next run" )

        # check the existence of the folders
        try :
            for folder in folderList:
                self.checkGRIDFolder( folder )

        except MissingGRIDFolderError, error :
            message = "Folder %(folder)s is unavailable. Quitting" % { "folder": error._filename }
            self._logger.critical( message )
            raise StopExecutionError( message )


        # check if the output file already exists
        command = "lfc-ls %(outputPathGRID)s/%(run)s-hit.slcio" % { "outputPathGRID": self._outputPathGRID, "run": self._option.output }
        lfc = popen2.Popen4( command )
        while lfc.poll() == -1:
            pass
        if lfc.poll() == 0:
            self._logger.warning( "Output file %(outputPathGRID)s/%(run)s-hit.slcio already exists"
                                  % { "outputPathGRID": self._outputPathGRID, "run": self._option.output } )
            if self._configParser.get("General","Interactive" ):
                if self.askYesNo( "Would you like to remove it?  [y/n] " ):
                    self._logger.info( "User decided to remove %(outputPathGRID)s/%(run)s-hit.slcio from the GRID"
                                       % { "outputPathGRID": self._outputPathGRID, "run": self._option.output } )
                    command = "lcg-del -a lfn:%(outputPathGRID)s/%(run)s-hit.slcio" % { "outputPathGRID": self._outputPathGRID,"run": self._option.output }
                    os.system( command )
                else :
                    raise OutputFileAlreadyOnGRIDError( "%(outputPathGRID)s/%(run)s-hit.slcio on the GRID"
                                                  % { "outputPathGRID": self._outputPathGRID, "run": self._option.output } )
            else :
                raise OutputAlreadyOnGRIDError( "%(outputPathGRID)s/%(run)s-hit.slcio on the GRID"
                                              % { "outputPathGRID": self._outputPathGRID, "run": self._option.output } )

        # check if the job output file already exists
        command = "lfc-ls %(outputPathGRID)s/%(name)s-%(run)s.tar.gz" % { 
            "name": self.name, "outputPathGRID": self._joboutputPathGRID, "run": self._option.output }
        lfc = popen2.Popen4( command )
        while lfc.poll() == -1:
            pass
        if lfc.poll() == 0:
            self._logger.warning( "Joboutput file %(outputPathGRID)s/%(name)s-%(run)s.tar.gz already exists"
                                  % { "name": self.name, "outputPathGRID": self._joboutputPathGRID, "run": self._option.output } )
            if self._configParser.get("General","Interactive" ):
                if self.askYesNo( "Would you like to remove it?  [y/n] " ):
                    self._logger.info( "User decided to remove %(outputPathGRID)s/%(name)s-%(run)s.tar.gz from the GRID"
                                       % { "name": self.name, "outputPathGRID": self._joboutputPathGRID, "run": self._option.output } )
                    command = "lcg-del -a lfn:%(outputPathGRID)s/%(name)s-%(run)s.tar.gz" % { 
                        "name": self.name, "outputPathGRID": self._joboutputPathGRID, "run": self._option.output }
                    os.system( command )
                else :
                    raise JoboutputFileAlreadyOnGRIDError( "%(outputPathGRID)s/%(name)s-%(run)s.tar.gz on the GRID"
                                                  % { "name": self.name, "outputPathGRID": self._joboutputPathGRID, "run": self._option.output } )
            else :
                raise JoboutputFileAlreadyOnGRIDError( "%(outputPathGRID)s/%(name)s-%(run)s.tar.gz on the GRID"
                                              % { "name": self.name, "outputPathGRID": self._joboutputPathGRID, "run": self._option.output } )


        # check if the histogram file already exists
        command = "lfc-ls %(outputPathGRID)s/%(run)s-hit-histo.root" % { "outputPathGRID": self._histogramPathGRID, "run": self._option.output }
        lfc = popen2.Popen4( command )
        while lfc.poll() == -1:
            pass
        if lfc.poll() == 0:
            self._logger.warning( "Histogram file %(outputPathGRID)s/%(run)s-hit-histo.root already exists"
                                  % { "outputPathGRID": self._histogramPathGRID, "run": self._option.output } )
            if self._configParser.get("General","Interactive" ):
                if self.askYesNo( "Would you like to remove it?  [y/n] " ):
                    self._logger.info( "User decided to remove %(outputPathGRID)s/%(run)s-hit-histo.root from the GRID"
                                       % { "outputPathGRID": self._histogramPathGRID, "run": self._option.output } )
                    command = "lcg-del -a lfn:%(outputPathGRID)s/%(run)s-hit-histo.root" % { "outputPathGRID": self._histogramPathGRID, "run": self._option.output }
                    os.system( command )
                else :
                    raise HistogramFileAlreadyOnGRIDError( "%(outputPathGRID)s/run%(run)s-hit-histo.root on the GRID"
                                                           % { "outputPathGRID": self._histogramPathGRID, "run": self._option.output } )
            else :
                raise HistogramFileAlreadyOnGRIDError( "%(outputPathGRID)s/%(run)s-hit-histo.root on the GRID"
                                                       % { "outputPathGRID": self._histogramPathGRID, "run": self._option.output } )
    ## Execute all GRID
    #
    def executeAllGRID( self  ) :

        # do preliminary checks
        self.doPreliminaryTest( )

        # now prepare the jdl file from the template
        self.generateJDLFile( 0, self._option.output )

        # now generate the run job script
        self.generateRunjobFile(  )

        # don't forget to generate the steering file!
        self._steeringFileName = self.generateSteeringFile( )

        # finally ready to submit! Let's do it...!
        self.submitJDL(  )

        # prepare a tarball with all the ancillaries
        self.prepareTarball(  )

        # cleaning up the system
        self.cleanup(  )




    ## Generate the run job
    #
    # This method is used to generate the run job script
    #
    def generateRunjobFile( self ):
        message = "Generating the executable (%(name)s-%(run)s.sh)" % { "name": self.name, "run": self._option.output }
        self._logger.info( message )
        try :
            runTemplate = self._configParser.get( "SteeringTemplate", "HitmakerGRIDScript" )
        except ConfigParser.NoOptionError:
            message = "Unable to find a valid executable template"
            self._logger.critical( message )
            raise StopExecutionError( message )

        runTemplateString = open( runTemplate, "r" ).read()
        runActualString = runTemplateString

        # replace the name prefix
        runActualString = runActualString.replace( "@Name@", self.name )

        # replace the output prefix.
        runActualString = runActualString.replace( "@Output@", self._option.output )

        # replace the event range
        # runActualString = runActualString.replace( "@FirstEvent@",        self._eventRange_begin  )
        # runActualString = runActualString.replace( "@LastEvent@",        self._eventRange_end  )
        
        # replace the input file names
        for index, inputFile in enumerate( self._inputFileList ) :
            if inputFile != "DEADFACE" :
                runActualString = runActualString.replace( "@InputFileList@", "%(file)s @InputFileList@" % {"file": self._justInputFileList[index] } )
        runActualString = runActualString.replace( "@InputFileList@" , "" )

        # replace the eta file. Again strip away any possible path from this
        etaFile = os.path.basename( self._option.eta )
        runActualString = runActualString.replace( "@EtaFile@", etaFile )

        if self._option.offsetRunNumber == None :
           offsetRunNumber = 0
        else:
           offsetRunNumber = "%(run)s" % { "run": self._option.offsetRunNumber }
           
        fqOffsetFile = "run0%(run)s-offset-db.slcio" %{ "run": offsetRunNumber }
        runActualString = runActualString.replace("@OffsetFile@", fqOffsetFile )
 
        variableList = [ "GRIDCE", "GRIDSE", "GRIDStoreProtocol", "GRIDVO",
                         "GRIDFolderBase", "GRIDFolderFilterResults", "GRIDFolderDBEta", 
                         "GRIDFolderDBOffset",  "GRIDFolderDBHotPixel", "GRIDFolderHitmakerResults",
                         "GRIDFolderHitmakerJoboutput", "GRIDFolderHitmakerHisto", "GRIDLibraryTarball",
                         "GRIDLibraryTarballPath","GRIDILCSoftVersion" ]
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

        self._runScriptFilename = "%(name)s-%(run)s.sh" % { "name": self.name, "run": self._option.output }
        runActualFile = open( self._runScriptFilename, "w" )
        runActualFile.write( runActualString )
        runActualFile.close()

        # change the mode of the run script to 0777
        os.chmod(self._runScriptFilename, 0777)

    ## Do the real submission
    #
    # This is doing the job submission
    #
    def submitJDL( self ) :
        self._logger.info("Submitting the job to the GRID")
        command = "glite-wms-job-submit %(del)s -r %(GRIDCE)s -o %(name)s-%(run)s.jid %(name)s-%(run)s.jdl" % {
            "name": self.name, "run": self._option.output , "GRIDCE":self._gridCE, "del": self._jobDelegation }
        status, output = commands.getstatusoutput( command )
        for line in output.splitlines():
            self._logger.log( 15, line.strip() )
        if status == 0:
            self._logger.info( "Job successfully submitted to the GRID" )
            for index, inputFile in enumerate( self._inputFileList ):
                if inputFile != "DEADFACE":
                    run, b, c, d, e, f = self._summaryNTuple[ index ]
                    self._summaryNTuple[ index ] = run, b, "See below", d, e, f
            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, "Submit'd", output, histo, tarball

        else :
            for index, inputFile in enumerate( self._inputFileList ):
                if inputFile != "DEADFACE":
                    run, b, c, d, e, f = self._summaryNTuple[ index ]
                    self._summaryNTuple[ index ] = run, b, "See below", d, e, f
            run, input, marlin, output, histo, tarball = self._summaryNTuple[ len( self._summaryNTuple ) - 1 ]
            self._summaryNTuple[ len( self._summaryNTuple ) - 1 ] = run, input, "Failed", output, histo, tarball
            raise GRIDSubmissionError ( "%(name)s-%(run)s.jdl" % { "name": self.name,"run": self._option.output } )

        # read back the the job id file
        # this is made by two lines only the first is comment, the second
        # is what we are interested in !
        jidFile = open( "%(name)s-%(run)s.jid" % { "name": self.name, "run": self._option.output} )
        jidFile.readline()
        for index, inputFile in enumerate ( self._inputFileList ) :
            if inputFile != "DEADFACE":
                run, b = self._gridJobNTuple[ index ]
                self._gridJobNTuple[ index ] = run, "See below"
        run, b = self._gridJobNTuple[ len( self._gridJobNTuple ) - 1 ]
        self._gridJobNTuple[ len( self._gridJobNTuple ) - 1  ] = run, jidFile.readline()
        jidFile.close()





