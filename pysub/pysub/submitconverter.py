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

## Submit conversion jobs
#
#  This class is responsible to submit jobs for native to lcio-raw format conversion.
#  It is inheriting from SubmitBase and it is called by the submit-converter.py script
#
#
#
#  @version $Id: submitconverter.py,v 1.40 2009-07-28 09:24:41 bulgheroni Exp $
#  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
#
class SubmitConverter( SubmitBase ) :

    ## Version
    # This number is used when printing out the version of the package.
    #
    # Static member.
    #
    cvsVersion = "$Revision: 1.40 $"

    ## Name
    # This is the namer of the class. It is used in flagging all the log entries
    # and preparing all the file
    #
    # Static member.
    #
    name = "converter"

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
        self._logger.log( 15, "Started submit-%(name)s" % { "name": self.name } )
        self._logger.log( 15, "**********************************************************" )

        # now print the configuration to the log
        SubmitBase.logConfigurationFile( self )

        # now log the run list
        SubmitBase.logRunList( self )

        # now check the I/O options
        # default values depend on the execution mode
        self._keepInput  = True
        self._keepOutput = True

        # set default for writing out hot pixel db
        self._keepHotPixelOutput   = False
        self._keepHotPixelDBNumber = 0 


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

        if  self._option.force_hotpixel_output != 0 and not self._keepHotPixelOutput:
            self._keepHotPixelOutput = True
            self._logger.info( "User forces to save hot pixel db (even if it is not a off-beam run)" )
            self._keepHotPixelDBNumber = self._option.force_hotpixel_output 





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

        # this is the real part
        # convert all the argumets into a run string compatible with the file name convention
        self._runList = [];
        for i in self._args:
            try:
                if( self._runList.count(int(i)) == 0 ):
                  self._runList.append( int(i) )

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
                    self.executeAllGRID( index , runString )

                elif self._option.execution == "all-local" :
                    self.executeAllLocal( index , runString )

                elif self._option.execution == "cpu-local":
                    self.executeCPULocal( index , runString )

                elif self._option.execution == "only-generate":
                    self.executeOnlyGenerate( index , runString )

            except GRID_LCG_CPError, error:
                message = "Unable to copy %(file)s" % { "file": error._filename }
                self._logger.error( message )
                self._logger.error("Skipping to the next run ")
                run, b, c, d, e, f = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, "Missing", "Skipped", d,e,f

            except MissingInputFileError, error:
                message = "Missing input file %(file)s" % { "file": error._filename }
                self._logger.error( message )
                self._logger.error("Skipping to the next run ")
                run, b, c, d, e, f = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, "Missing", "Skipped", d,e,f

            except MissingInputFileOnGRIDError, error:
                message = "Missing input file %(file)s" % { "file": error._filename }
                self._logger.error( message )
                self._logger.error("Skipping to the next run ")
                run, input, marlin, output, histo, tarball = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, "Missing", "Skipped", output, histo, tarball

            except OutputFileAlreadyOnGRIDError, error:
                message ="File %(file)s already on SE" % { "file":error._filename }
                self._logger.error( message )
                self._logger.error( "Skipping to the next run" )
                run, b, c, d, e, f = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, b, "Skipped", "GRID",  "N/A", f

            except JoboutputFileAlreadyOnGRIDError, error:
                message ="File %(file)s already on SE" % { "file":error._filename }
                self._logger.error( message )
                self._logger.error( "Skipping to the next run" )
                run, b, c, d, e, f = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, b, "Skipped", "d",  "N/A", "GRID"

            except MissingSteeringTemplateError, error:
                message = "Steering template %(file)s unavailble. Quitting!" % { "file": error._filename }
                self._logger.critical( message )
                raise StopExecutionError( message )

            except MissingGEARFileError, error:
                message = "Missing GEAR file %(file)s. Quitting! "  % { "file": error._filename }
                self._logger.critical( message )
                raise StopExecutionError( message )

            except MarlinError, error:
                message = "Error with Marlin execution (%(msg)s - errno = %(errno)s )" % { "msg": error._message, "errno": error._errno }
                self._logger.error( message )
                self._logger.error("Skipping to the next run ")
                run, input, marlin, output, histo, tarball = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, input, "Failed", "Missing", "Missing", "Missing"

            except GRIDSubmissionError, error:
                message = error._message
                self._logger.error( message )
                self._logger.error( "Skipping to the next run ")
                run, b, c, d, e, f = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, b, "Failed", f, e, f

            except MissingOutputFileError, error:
                message = "The output file (%(file)s) was not properly generated, possible failure" % { "file": error._filename }
                self._logger.error( message )
                run, b, c, d, e, f = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, b, c, "Missing", "N/A", f

            except GRID_LCG_CRError, error:
                message = "The file (%(file)s) couldn't be copied on the GRID"  % { "file": error._filename }
                self._logger.error( message )

            except MissingJoboutputFileError, error:
                message = "The joboutput tarball (%(file)s) is missing, possible failure" % { "file": error._filename }
                self._logger.error( message )
                run, b, c, d, e, f = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, b, c, d, e, "Missing"

            except MissingLibraryFileError, error:
                self._logger.critical( "The GRID lib (%(file)s) is not avaialable. Can't continue" % { "file": error._filename } )
                raise StopExecutionError( "Missing library" )


    ## CPU Local submitter
    #
    # This methods represents the sequence of commands that should be done while
    # submitting jobs on the local computer but using remote data
    #
    def executeCPULocal( self, index , runString ):

        # perform preliminary checks
        self.doPreliminaryTest( index, runString )

        # get the input file from the GRID
        self.getRunFromGRID( index, runString )

        # double check the presence of the input file
        self.checkInputFile( index, runString )

        #  generate the steering file
        self._steeringFileName = self.generateSteeringFile( runString )

        # prepare a separate file for logging the output of Marlin
        self._logFileName = "%(name)s-%(run)s.log" % { "name": self.name, "run" : runString }

        # run marlin
        self.runMarlin( index, runString )

        # advice the user that Marlin is over
        self._logger.info( "Marlin finished successfully")

        # verify the presence of the output files
        self.checkOutputFile( index, runString )

        # prepare a tarbal for the records
        self.prepareTarball( runString )

        try :
            # copy the output LCIO file to the GRID
            self.putRunOnGRID( index, runString )

            # copy the joboutput to the GRID
            self.putJoboutputOnGRID( index, runString )

            # clean up the local pc
            self.cleanup( runString )

        except GRID_LCG_CRError, error:
            message = "The file (%(file)s) couldn't be copied on the GRID"  % { "file": error._filename }
            self._logger.error( message )



    ## Local submitter
    #
    # This methods represents the sequence of commands that should be done while
    # submitting jobs on the local computer.
    #
    def executeAllLocal( self, index , runString ):

        # before any futher, check we have the input file for this run
        self.checkInputFile( index, runString )

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

        # prepare a tarbal for the records
        self.prepareTarball( runString )

        # check the presence of the joboutput tarball
        # this should be named something like
        # converter-123456.tar.gz
        self.checkJoboutputFile( index, runString )

        # clean up the local pc
        self.cleanup( runString )

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


    ## Preliminary checks
    #
    # This method performs some preliminary checks
    # before starting a full GRID submission.
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
            self._inputPathGRID     = self._configParser.get("GRID", "GRIDFolderNative")
            self._outputPathGRID    = self._configParser.get("GRID", "GRIDFolderLcioRaw" )
            self._joboutputPathGRID = self._configParser.get("GRID", "GRIDFolderConvertJoboutput")
            self._dbHotPixelGRID    = self._configParser.get("GRID", "GRIDFolderDBHotPixel")
            folderList = [self._inputPathGRID , self._outputPathGRID, self._joboutputPathGRID, self._dbHotPixelGRID ]
        except ConfigParser.NoOptionError:
            message = "Missing path from the configuration file"
            self._logger.critical( message )
            raise StopExecutionError( message )

        # check if the input file is on the GRID, otherwise go to next run
        command = "lfc-ls %(inputPathGRID)s/run%(run)s.raw" % { "inputPathGRID" : self._inputPathGRID,  "run": runString }
        status, output = commands.getstatusoutput( command )
        if status == 0:
            self._logger.info( "Input file found on the SE" )
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, "GRID", c, d, "N/A", f
        else:
            self._logger.error( "Input file NOT found on the SE. Trying next run" )
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, "Missing", c, d, "N/A", f
            raise MissingInputFileOnGRIDError( "%(inputPathGRID)s/run%(run)s.raw" % { "inputPathGRID" : self._inputPathGRID,  "run": runString } )

        # verify that those folders really exists
        self._logger.info( "Checking the presence of the output folders on the GRID" )
        try :
            for folder in folderList:
                self.checkGRIDFolder( folder )

        except MissingGRIDFolderError, error :
            message = "Folder %(folder)s is unavailable. Quitting" % { "folder": error._filename }
            self._logger.critical( message )
            raise StopExecutionError( message )

        # check if the output file already exists
        command = "lfc-ls %(outputPathGRID)s/run%(run)s.slcio" % { "outputPathGRID": self._outputPathGRID, "run": runString }
        status, output = commands.getstatusoutput( command )
        if status == 0:
            self._logger.warning( "Output file %(outputPathGRID)s/run%(run)s.slcio already exists" % { "outputPathGRID": self._outputPathGRID, "run": runString } )
            if self._configParser.get("General","Interactive" ):
                if self.askYesNo( "Would you like to remove it?  [y/n] " ):
                    self._logger.info( "User decided to remove %(outputPathGRID)s/run%(run)s.slcio from the GRID"
                                       % { "outputPathGRID": self._outputPathGRID, "run": runString } )
                    command = "lcg-del -a lfn:%(outputPathGRID)s/run%(run)s.slcio" % { "outputPathGRID": self._outputPathGRID, "run": runString }
                    os.system( command )
                else :
                    raise OutputFileAlreadyOnGRIDError( "%(outputPathGRID)s/run%(run)s.slcio on the GRID"
                                                  % { "outputPathGRID": self._outputPathGRID, "run": runString } )
            else :
                raise OutputFileAlreadyOnGRIDError( "%(outputPathGRID)s/run%(run)s.slcio on the GRID"
                                              % { "outputPathGRID": self._outputPathGRID, "run": runString } )

        # check if the job output file already exists
        command = "lfc-ls %(outputPathGRID)s/%(name)s-%(run)s.tar.gz" \
            % { "name": self.name, "outputPathGRID": self._joboutputPathGRID, "run": runString }
        status, output = commands.getstatusoutput( command )
        if status == 0:
            self._logger.warning( "Output file %(outputPathGRID)s/%(name)s-%(run)s.tar.gz already exists"
                                  % { "name": self.name, "outputPathGRID": self._joboutputPathGRID, "run": runString } )
            if self._configParser.get("General","Interactive" ):
                if self.askYesNo( "Would you like to remove it?  [y/n] " ):
                    self._logger.info( "User decided to remove %(outputPathGRID)s/%(name)s-%(run)s.tar.gz from the GRID"
                                       % { "name": self.name, "outputPathGRID": self._joboutputPathGRID, "run": runString } )
                    command = "lcg-del -a lfn:%(outputPathGRID)s/%(name)s-%(run)s.tar.gz" % \
                        { "name": self.name, "outputPathGRID": self._joboutputPathGRID, "run": runString }
                    os.system( command )
                else :
                    raise JoboutputFileAlreadyOnGRIDError( "%(outputPathGRID)s/%(name)s-%(run)s.tar.gz on the GRID"
                                                  % { "name": self.name, "outputPathGRID": self._joboutputPathGRID, 
                                                      "run": runString } )
            else :
                raise JoboutputFileAlreadyOnGRIDError( "%(outputPathGRID)s/%(name)s-%(run)s.tar.gz on the GRID"
                                              % { "name": self.name, "outputPathGRID": self._joboutputPathGRID, "run": runString } )



    ## Generate the run job
    #
    # This method is used to generate the run job script
    #
    def generateRunjobFile( self, index, runString ):
        message = "Generating the executable (%(name)s-%(run)s.sh)" % { "name": self.name, "run": runString }
        self._logger.info( message )
        try :
            runTemplate = self._configParser.get( "SteeringTemplate", "ConverterGRIDScript" )
        except ConfigParser.NoOptionError:
            message = "Unable to find a valid executable template"
            self._logger.critical( message )
            raise StopExecutionError( message )

        runTemplateString = open( runTemplate, "r" ).read()
        runActualString = runTemplateString

        # replace the runString
        runActualString = runActualString.replace( "@RunString@", runString )

        # replace the hotpixelString as well
        runActualString = runActualString.replace( "@HotPixelRunNumber@", self._keepHotPixelDBNumber )
        message = "Replacing HotPixelRunNumber with %(var)s          " % { "var" : self._keepHotPixelDBNumber }
        self._logger.info( message )
#        self._logger.info( runActualString ) 

        # replace the job name
        runActualString = runActualString.replace( "@Name@", self.name )

        variableList = [ "GRIDCE", "GRIDSE", "GRIDStoreProtocol", "GRIDVO",
                         "GRIDFolderBase", "GRIDFolderNative", "GRIDFolderLcioRaw", "GRIDFolderDBHotPixel",
                         "GRIDFolderConvertJoboutput", "GRIDLibraryTarball", "GRIDLibraryTarballPath", "GRIDILCSoftVersion" ]
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
        self._logger.info( command )
        status, output = commands.getstatusoutput( command )
        for line in output.splitlines():
            self._logger.log( 15, line.strip() )

        if status == 0 :
            self._logger.info( "Job successfully submitted to the GRID" )
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, "Submit'd", d, "N/A", f
        else :
            raise GRIDSubmissionError ( "%(name)s-%(run)s.jdl" % { "name": self.name, "run": runString } )

        # read back the the job id file
        # this is made by two lines only the first is comment, the second
        # is what we are interested in !
        jidFile = open( "%(name)s-%(run)s.jid" % { "name": self.name, "run": runString } )
        jidFile.readline()
        run, b = self._gridJobNTuple[ index ]
        self._gridJobNTuple[ index ] = run, jidFile.readline()
        jidFile.close()

    ## Generate only submitter
    #
    # This methods is responsibile of dry-run with only steering file
    # generation
    #
    def executeOnlyGenerate( self, index , runString ):

        # just need to generate the steering file
        self.generateSteeringFile( runString )

    ## Get the input run from the GRID
    def getRunFromGRID(self, index, runString ) :

        self._logger.info(  "Getting the input file from the GRID" )

        try :
            gridNativePath = self._configParser.get( "GRID", "GRIDFolderNative" )
        except ConfigParser.NoOptionError :
            message = "GRIDFolderNative missing in the configuration file. Quitting."
            self._logger.critical( message )
            raise StopExecutionError( message )

        try :
            localPath = self._configParser.get( "LOCAL", "LocalFolderNative" )
        except ConfigParser.NoOptionError :
            localPath = "$PWD/native"

        baseCommand = "lcg-cp "
        if self._option.verbose :
            baseCommand = baseCommand + " -v "

        command = "%(base)s lfn:%(gridNativePath)s/run%(run)s.raw file:%(localPath)s/run%(run)s.raw" %  \
            { "base": baseCommand, "gridNativePath" : gridNativePath, "run": runString, "localPath": localPath }
        if os.system( command ) != 0:
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, "Missing", c, d, "N/A", f
            raise GRID_LCG_CPError( "lfn:%(gridNativePath)s/run%(run)s.raw" %   { "gridNativePath" : gridNativePath, "run": runString } )
        else:
            self._logger.info("Input file successfully copied from the GRID")
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, "OK", c, d, "N/A", f


    ## Put the output run to the GRID
    def putRunOnGRID( self, index, runString ):

        self._logger.info(  "Putting the LCIO file to the GRID" )

        try :
            gridLcioRawPath = self._configParser.get( "GRID", "GRIDFolderLcioRaw")
        except ConfigParser.NoOptionError :
            message = "GRIDFolderLcioRaw missing in the configuration file. Quitting."
            self._logger.critical( message )
            raise StopExecutionError( message )

        try :
            localPath = self._configParser.get( "LOCAL", "LocalFolderLcioRaw" )
        except ConfigParser.NoOptionError :
            localPath = "$PWD/lcio-raw"

        baseCommand = "lcg-cr "
        if self._option.verbose :
            baseCommand = baseCommand + " -v "

        command = "%(base)s -l lfn:%(gridFolder)s/run%(run)s.slcio file:%(localFolder)s/run%(run)s.slcio" % \
            { "base": baseCommand, "gridFolder": gridLcioRawPath, "localFolder": localPath, "run" : runString }
        if os.system( command ) != 0 :
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, c, "LOCAL", "N/A", f
            raise GRID_LCG_CRError( "lfn:%(gridFolder)s/run%(run)s.slcio" % \
                                        { "gridFolder": gridLcioRawPath, "run" : runString } )
        else:
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, c, "GRID", "N/A", f
            self._logger.info( "Output file successfully copied to the GRID" )

            # verify?
            # if so, first calculate the sha of the local file
            if self._option.verify_output :
                self._logger.info ("Verifying the LCIO file integrity on the GRID" )
                localCopy = open( "%(localFolder)s/run%(run)s.slcio" % { "localFolder": localPath, "run" : runString }, "r" ).read()
                localCopyHash = sha.new( localCopy ).hexdigest()
                message = "Local copy hash %(hash)s" % { "hash": localCopyHash }
                self._logger.log( 15, message )

                # now copy back the remote file
                # if so we need to copy back the file
                baseCommand = "lcg-cp "
                if self._option.verbose :
                    baseCommand = baseCommand + " -v "

                command = "%(base)s lfn:%(gridFolder)s/run%(run)s.slcio file:%(localFolder)s/run%(run)s-test.slcio" % \
                          { "base": baseCommand, "gridFolder": gridLcioRawPath, "localFolder": localPath, "run" : runString }
                if os.system( command ) != 0 :
                    run, b, c, d, e, f = self._summaryNTuple[ index ]
                    self._summaryNTuple[ index ] = run, b, c, "LOCAL", "N/A", f
                    raise GRID_LCG_CRError( "lfn:%(gridFolder)s/run%(run)s.slcio" % \
                                    { "gridFolder": gridLcioRawPath, "run" : runString } )

                remoteCopy = open( "%(localFolder)s/run%(run)s-test.slcio" % { "localFolder": localPath, "run" : runString }, "r" ).read()
                remoteCopyHash = sha.new( remoteCopy ).hexdigest()
                message = "Remote copy hash %(hash)s" % { "hash": remoteCopyHash }
                self._logger.log( 15, message )

                if remoteCopyHash == localCopyHash:
                    run, b, c, d, e, f = self._summaryNTuple[ index ]
                    self._summaryNTuple[ index ] = run, b, c, "GRID - Ver!", "N/A", f
                    self._logger.info( "Verification successful" )
                else:
                    run, b, c, d, e, f = self._summaryNTuple[ index ]
                    self._summaryNTuple[ index ] = run, b, c, "GRID - Fail!", "N/A", f
                    self._logger.error( "Problem with the verification!" )

                # remove the test file
                os.remove( "%(localFolder)s/run%(run)s-test.slcio" % { "localFolder": localPath, "run" : runString } )


    ## Put the joboutput to the GRID
    def putJoboutputOnGRID( self, index, runString ):

        self._logger.info(  "Putting the joboutput file to the GRID" )

        try :
            gridFolder = self._configParser.get( "GRID", "GRIDFolderConvertJoboutput")
        except ConfigParser.NoOptionError :
            message = "GRIDFolderConvertJoboutput missing in the configuration file. Quitting."
            self._logger.critical( message )
            raise StopExecutionError( message )

        try :
            localPath = self._configParser.get( "LOCAL", "LocalFolderConvertJoboutput")
        except ConfigParser.NoOptionError :
            localPath = "log/"

        baseCommand = "lcg-cp "
        if self._option.verbose :
            baseCommand = baseCommand + " -v "

        command = "lcg-cr -v -l lfn:%(gridFolder)s/%(name)s-%(run)s.tar.gz file:%(localFolder)s/%(name)s-%(run)s.tar.gz" % \
            { "name": self.name, "gridFolder": gridFolder, "localFolder": localPath, "run" : runString }

        if os.system( command ) != 0 :
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, c, d, e, "LOCAL"
            raise GRID_LCG_CRError( "lfn:%(gridFolder)s/%(name)s-%(run)s.tar.gz"% \
                                        { "name": self.name, "gridFolder": gridFolder, "run" : runString } )
        else:
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, c, d, e, "GRID"
            self._logger.info("Jobouput file successfully copied from the GRID")

            # verify?
            # if so, first calculate the sha of the local file
            if self._option.verify_output :
                self._logger.info ("Verifying the joboutput file integrity on the GRID" )
                localCopy = open( "%(localFolder)s/%(name)s-%(run)s.tar.gz" % { 
                        "name": self.name, "localFolder": localPath, "run" : runString }, "r" ).read()
                localCopyHash = sha.new( localCopy ).hexdigest()
                message = "Local copy hash %(hash)s" % { "hash": localCopyHash }
                self._logger.log( 15, message )

                # now copy back the remote file
                # if so we need to copy back the file
                baseCommand = "lcg-cp "
                if self._option.verbose :
                    baseCommand = baseCommand + " -v "

                command = "%(base)s lfn:%(gridFolder)s/%(name)s-%(run)s.tar.gz file:%(localFolder)s/%(name)s-%(run)s-test.tar.gz" % \
                          {"name": self.name, "base": base, "gridFolder": gridFolder, "localFolder": localPath, "run" : runString }
                if os.system( command ) != 0 :
                    run, b, c, d, e, f = self._summaryNTuple[ index ]
                    self._summaryNTuple[ index ] = run, b, c, d, "N/A", "GRID - Fail!"
                    raise GRID_LCG_CRError( "lfn:%(gridFolder)s/run%(run)s.slcio" % \
                                    { "gridFolder": gridLcioRawPath, "run" : runString } )

                remoteCopy = open( "%(localFolder)s/%(name)s-%(run)s-test.tar.gz" % { 
                        "name": self.name, "localFolder": localPath, "run" : runString }, "r" ).read()
                remoteCopyHash = sha.new( remoteCopy ).hexdigest()
                message = "Remote copy hash %(hash)s" % { "hash": remoteCopyHash }
                self._logger.log( 15, message )

                if remoteCopyHash == localCopyHash:
                    run, b, c, d, e, f = self._summaryNTuple[ index ]
                    self._summaryNTuple[ index ] = run, b, c, "GRID - Ver!", "N/A", f
                    self._logger.info( "Verification successful" )
                else:
                    run, b, c, d, e, f = self._summaryNTuple[ index ]
                    self._summaryNTuple[ index ] = run, b, c, "GRID - Fail!", "N/A", f
                    self._logger.error( "Problem with the verification!" )

                # anyway remove the test file
                os.remove( "%(localFolder)s/%(name)s-%(run)s-test.tar.gz" % { "name": self.name, "localFolder": localPath, "run" : runString } )

    ## Generate the steering file
    def generateSteeringFile( self, runString  ) :
        message = "Generating the steering file (%(name)s-%(run)s.xml) " % { "name": self.name, "run" : runString }
        self._logger.info( message )

        steeringFileTemplate =  ""
        try:
            steeringFileTemplate = self._configParser.get( "SteeringTemplate", "ConverterSteeringFile" )
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

        # decide if User wants to write out hot pixel db file
        # involves AutoPedestal and HotPixelKiller  processors
        if self._keepHotPixelOutput == False:
             self._logger.info( "_keepHotPixelOutput = FALSE" )
             actualSteeringString = actualSteeringString.replace("@RunHotPixelKillerCommentLeft@", "!--" )
             actualSteeringString = actualSteeringString.replace("@RunHotPixelKillerCommentRight@", "--" )
        else:    
             self._logger.info( "_keepHotPixelOutput = TRUE" )
             actualSteeringString = actualSteeringString.replace("@RunHotPixelKillerCommentLeft@", "" )
             actualSteeringString = actualSteeringString.replace("@RunHotPixelKillerCommentRight@", "" )
            
        
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

        # now replace the native folder path
        if self._option.execution == "all-grid" :
            nativeFolder = "native"
        else :
            try:
                nativeFolder = self._configParser.get( "LOCAL", "LocalFolderNative" )
            except ConfigParser.NoOptionError :
                nativeFolder = "native"
        actualSteeringString = actualSteeringString.replace("@NativeFolder@", nativeFolder )

        # now replace the lcio-raw folder path
        if self._option.execution == "all-grid" :
            lcioRawFolder = "lcio-raw"
        else:
            try:
                lcioRawFolder = self._configParser.get("LOCAL", "LocalFolderLcioRaw")
            except ConfigParser.NoOptionError :
                lcioRawFolder = "lcio-raw"
        actualSteeringString = actualSteeringString.replace("@LcioRawFolder@" ,lcioRawFolder )

        # finally replace the run string !
        actualSteeringString = actualSteeringString.replace("@RunNumber@", runString )

        # if one runs the HotPixelKiller at this stage (with converter) 
        # now replace the histo folder path
        if self._option.execution == "all-grid" :
            self._logger.debug( " all-grid is true" )      
            histoFolder = "histo"
        else:
            self._logger.debug( " all-grid is NOT true" )      
            try:
                histoFolder = self._configParser.get("LOCAL", "LocalFolderConvertHisto")
            except ConfigParser.NoOptionError :
                histoFolder = "histo"
        actualSteeringString = actualSteeringString.replace("@HistoPath@" ,histoFolder )

        # now replace the DB folder path
        if self._option.execution == "all-grid" :
            dbPath = "db"
        else:
            try:
                dbPath = self._configParser.get("LOCAL", "LocalFolderDBAlign" )
            except  ConfigParser.NoOptionError :
                dbPath = "db"
        actualSteeringString = actualSteeringString.replace("@DBPath@", dbPath )
 
        # replace the DB run number (output)   
        if( self._keepHotPixelDBNumber != 0 ):
           actualSteeringString = actualSteeringString.replace("@HotPixelRunNumber@", self._keepHotPixelDBNumber )
           
            


        # replace the run string !
        actualSteeringString = actualSteeringString.replace("@Output@", "%(output)s" %
                                                            { "output": self._option.output } )

        # replace the run string !
        actualSteeringString = actualSteeringString.replace("@RunNumber@", runString )


        #-----------------------------

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
        marlin  = popen2.Popen4( "Marlin %(steer)s; Marlin -d %(steer)s %(steer)s.dot; " % { "steer": self._steeringFileName } )
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

    ## Prepare the joboutput tarball
    def prepareTarball( self, runString ):

        self._logger.info("Preparing the joboutput tarball" )

        # first prepare a folder to store them temporary
        destFolder = "%(name)s-%(run)s" %{ "name": self.name, "run": runString}
        shutil.rmtree( destFolder, True )
        os.mkdir( destFolder )

        # prepare the list of files we want to copy.
        listOfFiles = []
        try :
            gearPath = self._configParser.get( "LOCAL", "LocalFolderGear" )
        except ConfigParser.NoOptionError :
            gearPath = ""
        listOfFiles.append( os.path.join( gearPath, self._gear_file ) )
        listOfFiles.append( self._configFile )
        for file in glob.glob( "%(name)s-*.*" % {"name": self.name} ):
            message = "Adding %(file)s to the joboutput tarball" % { "file": file } 
            self._logger.debug( message )
            listOfFiles.append( file )

        for file in listOfFiles :
            shutil.copy( file , destFolder )

        # do the tarball
        self._tarballFileName = "%(name)s-%(run)s.tar.gz" % {"name": self.name,"run": runString }
        tarball = tarfile.open( self._tarballFileName, "w:gz" )
        tarball.add( destFolder )
        tarball.close()

        # remove the temporary folder
        shutil.rmtree( destFolder )

        # copy the tarball in the log folder
        try:
            localFolder = self._configParser.get( "LOCAL", "LocalFolderConvertJoboutput")
        except ConfigParser.NoOptionError :
            localFolder = "log/"

        dest = os.path.join( localFolder, self._tarballFileName )
        if os.path.exists(dest):
          os.remove(dest)

        shutil.move( self._tarballFileName, localFolder )



    ## Cleanup after each run conversion
    def cleanup( self, runString ):

        self._logger.info( "Cleaning up the local pc" )

        # remove the log file and the steering file
        for file in glob.glob( "%(name)s-*"  % {"name": self.name}):
            os.remove( file )


        # remove the input and output file
        try :
            inputFilePath = self._configParser.get( "LOCAL", "LocalFolderNative" )
        except ConfigParser.NoOptionError :
            inputFilePath = "native"
        inputFile  = "run%(run)s.raw" % { "run" : runString }
        if self._keepInput == False:
            os.remove( os.path.join( inputFilePath, inputFile ))

        try :
            outputFilePath = self._configParser.get( "LOCAL", "LocalFolderLcioRaw" )
        except ConfigParser.NoOptionError :
            outputFilePath = "lcio-raw"
        outputFile = "run%(run)s.slcio" % { "run" : runString }
        if self._keepOutput == False :
            os.remove( os.path.join(outputFilePath, outputFile ))

    ## Check the input file
    #
    def checkInputFile( self, index, runString ) :
        # the input file should be something like:
        # native/run123456.raw

        try :
            inputFilePath = self._configParser.get( "LOCAL", "LocalFolderNative" )
        except ConfigParser.NoOptionError :
            inputFilePath = "native"
        inputFileName = "run%(run)s.raw" % { "run": runString }
        if not os.access( os.path.join(inputFilePath, inputFileName) , os.R_OK ):
            message = "Problem accessing the input file (%(file)s), trying next run" % {"file": inputFileName }
            self._logger.error( message )
            raise MissingInputFileError( inputFileName )
        else :
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, "OK", c, d, e, f

    ## Check the output file
    #
    def checkOutputFile( self, index, runString ) :
        # this should be named something like
        # lcio-raw/run123456.slcio

        try :
            outputFilePath = self._configParser.get( "LOCAL", "LocalFolderLcioRaw" )
        except ConfigParser.NoOptionError :
            outputFilePath = "lcio-raw"
        outputFileName = "run%(run)s.slcio" % { "run": runString }
        if not os.access( os.path.join( outputFilePath , outputFileName) , os.R_OK ):
            raise MissingOutputFileError( outputFileName )
        else :
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, c, "OK", "N/A", f

    ## Check the joboutput file
    #
    def checkJoboutputFile( self, index, runString) :
        # this should be named something like
        # log/converter-123456.tar.gz

        try :
            outputFilePath = self._configParser.get( "LOCAL", "LocalFolderConvertJoboutput" )
        except ConfigParser.NoOptionError :
            outputFilePath = "log"
        tarballFileName = "%(name)s-%(run)s.tar.gz" % { "name": self.name,"run": runString }
        if not os.access( os.path.join( outputFilePath, tarballFileName), os.R_OK ):
            raise MissingJoboutputFileError( tarballFileName )
        else:
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, c, d, e, "OK"

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
        self._logger.info("Preparing the JID for this submission (%(name)s-%(date)s.jid)" % { "name": self.name, "date" : unique } )

        try :
            localPath = self._configParser.get( "LOCAL", "LocalFolderConvertJoboutput")
        except ConfigParser.NoOptionError :
            localPath = "log/"

        jidFile = open( os.path.join( localPath, "%(name)s-%(date)s.jid" % { "name": self.name, "date": unique } ), "w" )
        currentJIDFile = open( "current.jid" , "a" )
        for run, jid in self._gridJobNTuple:
            if jid != "Unknown" and jid != "See below":
                jidFile.write( jid )
                currentJIDFile.write( "# %(name)s %(run)s %(unique)s\n" % {"name":self.name, "run":run , "unique": unique } )
                currentJIDFile.write( jid )
        jidFile.close()
        currentJIDFile.close()
                                                                        
