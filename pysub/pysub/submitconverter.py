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
from submitbase import SubmitBase
from error import *

## Submit conversion jobs
#
#  This class is responsible to submit jobs for native to lcio-raw format conversion.
#  It is inheriting from SubmitBase and it is called by the submit-converter.py script
#
#
#
#  @version $Id: submitconverter.py,v 1.16 2009-05-11 12:47:51 bulgheroni Exp $
#  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
#
class SubmitConverter( SubmitBase ) :


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
        self._configFile = ""
        if ( self._option.config_file == None ) :
            # not from the option, check from environ
            try:
                self._configFile = os.environ['SUBMIT_CONFIG']
            except KeyError:
                # no way, just use the default one
                self._configFile = "config/config.cfg"
        else:
            self._configFile = self._option.config_file

        # before proceeding check if the configuration file
        # really exists!
        if not os.path.exists( self._configFile ):
            logging.critical( "Configuration file %(cfg)s doesn't exist!" % { "cfg": self._configFile } )

        # if it exists, then I can read it!!!
        self._configParser = ConfigParser.SafeConfigParser()
        self._configParser.read( self._configFile )

        # now we can properly set the logger.
        self.configureLogger()

        # now print the configuration to the log
        self.logConfigurationFile()

        # now log the run list
        self.logRunList()

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

    ## Logger configurator
    #
    # This is another configure method that is called by the main configure
    # to properly set up the logging system.
    # Before calling this method only a very simple logging system is available
    #
    def configureLogger( self ):

        if not self._configParser.has_section( "Logger" ) :
            logging.warning( "No Logger section in the configuration file, using default logging" )
            return

        # rename the logger according to this module
        self._logger = logging.getLogger( "Converter" )
        # set the global logging level
        self._logger.setLevel( self._configParser.getint( "Logger", "GlobalLoggerLevel" ))

        # configuring the console handler
        self._logger.removeHandler( self._consoleHandler )
        if self._configParser.getboolean( "Logger", "ConsoleHandler" ):
            self._consoleHandler = logging.StreamHandler()
            self._consoleHandler.setLevel( self._configParser.getint( "Logger", "ConsoleHandlerLevel" ) )
            self._consoleHandler.setFormatter(logging.Formatter("%(asctime)s - %(name)-10s [%(levelname)-8s]: %(message)s","%a, %d %b %Y %H:%M:%S"
                                                                ))
            self._logger.addHandler( self._consoleHandler )

        # configuring the time rotating file handler
        if self._configParser.getboolean( "Logger", "RotatingFileHandler" ):
            try:
                rotatingHandler = logging.handlers.RotatingFileHandler(
                    filename = self._configParser.get( "Logger", "RotatingFileHandlerFileName" ),
                    maxBytes = self._configParser.getint( "Logger", "RotatingFileHandlerSize" ),
                    backupCount = 10)
                rotatingHandler.setLevel( self._configParser.getint( "Logger", "RotatingFileHandlerLevel" ) )
                rotatingHandler.setFormatter( logging.Formatter("%(asctime)s - %(name)-10s [%(levelname)-8s]: %(message)s","%a, %d %b %Y %H:%M:%S") )
                self._logger.addHandler( rotatingHandler )
            except IOError, detail:
                message = "IOError: %(detail)s" % { "detail":detail }
                self._logger.error( message )
                self._logger.error( "Impossible to open the rotating handler for logging" )
        self._logger.log( 15, "**********************************************************" )
        self._logger.log( 15, "Started submit-converter" )
        self._logger.log( 15, "**********************************************************" )

    ## Print the configuration file
    #
    # This method is used to print to the logger the options specified in
    # the configuration file.
    # The log level for this is ALL = 15 and it is usually not inclueded into the
    # stdout
    #
    def logConfigurationFile( self ):
        self._logger.log(15,"")
        self._logger.log(15, "Logging the configuration file " )
        for section in self._configParser.sections():
            self._logger.log(15,"")
            message = "[ Section: %(section)s ]" % {'section':section}
            self._logger.log(15, message)
            for option in self._configParser.options(section):
                message = "    %(option)-35s = %(value)s " % { "option":option, "value":self._configParser.get( section, option) }
                self._logger.log(15,  message )

        self._logger.log(15, "---------------------------------------")
        self._logger.log(15,"")


    ## Print the run list to the logger
    #
    # Simple method to print out on the log all the runs we are about to process
    #
    def logRunList( self ):
        self._logger.log(15, "" )
        self._logger.log(15, "Logging the run list ")
        message = ""
        i = 0
        for j, run in enumerate( self._args ) :
            message = ""
            for k in range(0,4):
                if i < len( self._args) :
                    message = message + "%(run)10s " % { "run": self._args[i] }
                    i = i + 1
            if len(message) != 0:
                self._logger.log(15, message )

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
                self._runList.append( int( i  ) )

                # if it is a good run number than we can prepare an entry for the summary ntuple
                # the ntuple contains 6 variables:
                # runNumber ; inputFileStatus ; marlinStatus ; outputFileStatus ; histoFileStatus ; joboutputFileStatus
                entry = i , "Unknown", "Unknown", "Unknown", "Unknown", "Unknown"

                # the current entry to the ntuple
                self._summaryNTuple.append( entry )

            except ValueError:
                message = "Invalid run number %(i)s" % { "i": i }
                self._logger.critical( message )
                sys.exit( 1 )


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

    ## This is the real all-grid submitter
    def allGridSubmission( self ) :
        pass

    ## CPU Local submitter
    #
    # This methods represents the sequence of commands that should be done while
    # submitting jobs on the local computer but using remote data
    #
    def executeCPULocal( self, index , runString ):

        # get the input file from the GRID
        self.getRunFromGRID( index, runString )

        # double check the presence of the input file
        self.checkInputFile( index, runString )

        #  generate the steering file
        self._steeringFileName = self.generateSteeringFile( runString )

        # prepare a separate file for logging the output of Marlin
        self._logFileName = "universal-%(run)s.log" % { "run" : runString }

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
        self._logFileName = "universal-%(run)s.log" % { "run" : runString }

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
        # universal-123456.tar.gz
        self.checkJoboutputFile( index, runString )

        # clean up the local pc
        self.cleanup( runString )

    ## Execute all GRID
    #
    def executeAllGRID( self, index, runString ) :
        pass

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

        command = "lcg-cp -v lfn:%(gridNativePath)s/run%(run)s.raw file:%(localPath)s/run%(run)s.raw" %  \
            { "gridNativePath" : gridNativePath, "run": runString, "localPath": localPath }
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

        command = "lcg-cr -v -l lfn:%(gridFolder)s/run%(run)s.slcio file:%(localFolder)s/run%(run)s.slcio" % \
            { "gridFolder": gridLcioRawPath, "localFolder": localPath, "run" : runString }
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
                command = "lcg-cp -v lfn:%(gridFolder)s/run%(run)s.slcio file:%(localFolder)s/run%(run)s-test.slcio" % \
                          { "gridFolder": gridLcioRawPath, "localFolder": localPath, "run" : runString }
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

                # anyway remove the test file
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

        command = "lcg-cr -v -l lfn:%(gridFolder)s/universal-%(run)s.tar.gz file:%(localFolder)s/universal-%(run)s.tar.gz" % \
            { "gridFolder": gridFolder, "localFolder": localPath, "run" : runString }

        if os.system( command ) != 0 :
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, c, d, e, "LOCAL"
            raise GRID_LCG_CRError( "lfn:%(gridFolder)s/universal-%(run)s.tar.gz"% \
                                        { "gridFolder": gridFolder, "run" : runString } )
        else:
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, c, d, e, "GRID"
            self._logger.info("Jobouput file successfully copied from the GRID")

            # verify?
            # if so, first calculate the sha of the local file
            if self._option.verify_output :
                self._logger.info ("Verifying the joboutput file integrity on the GRID" )
                localCopy = open( "%(localFolder)s/universal-%(run)s.tar.gz" % { "localFolder": localPath, "run" : runString }, "r" ).read()
                localCopyHash = sha.new( localCopy ).hexdigest()
                message = "Local copy hash %(hash)s" % { "hash": localCopyHash }
                self._logger.log( 15, message )

                # now copy back the remote file
                # if so we need to copy back the file
                command = "lcg-cp -v lfn:%(gridFolder)s/universal-%(run)s.tar.gz file:%(localFolder)s/universal-%(run)s-test.tar.gz" % \
                          { "gridFolder": gridFolder, "localFolder": localPath, "run" : runString }
                if os.system( command ) != 0 :
                    run, b, c, d, e, f = self._summaryNTuple[ index ]
                    self._summaryNTuple[ index ] = run, b, c, d, "N/A", "GRID - Fail!"
                    raise GRID_LCG_CRError( "lfn:%(gridFolder)s/run%(run)s.slcio" % \
                                    { "gridFolder": gridLcioRawPath, "run" : runString } )

                remoteCopy = open( "%(localFolder)s/universal-%(run)s-test.tar.gz" % { "localFolder": localPath, "run" : runString }, "r" ).read()
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
                os.remove( "%(localFolder)s/universal-%(run)s-test.tar.gz" % { "localFolder": localPath, "run" : runString } )

    ## Generate the steering file
    def generateSteeringFile( self, runString  ) :
        message = "Generating the steering file (universal-%(run)s.xml) " % { "run" : runString }
        self._logger.info( message )

        steeringFileTemplate =  ""
        try:
            steeringFileTemplate = self._configParser.get( "SteeringTemplate", "ConverterSteeringFile" )
        except ConfigParser.NoOptionError :
            steeringFileTemplate = "template/universal-tmp.xml"

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
        try:
            nativeFolder = self._configParser.get( "LOCAL", "LocalFolderNative" )
        except ConfigParser.NoOptionError :
            nativeFolder = "native"
        actualSteeringString = actualSteeringString.replace("@NativeFolder@", nativeFolder )

        # now replace the lcio-raw folder path
        try:
            lcioRawFolder = self._configParser.get("LOCAL", "LocalFolderLcioRaw")
        except ConfigParser.NoOptionError :
            lcioRawFolder = "lcio-raw"
        actualSteeringString = actualSteeringString.replace("@LcioRawFolder@" ,lcioRawFolder )

        # finally replace the run string !
        actualSteeringString = actualSteeringString.replace("@RunNumber@", runString )

        # open the new steering file for writing
        steeringFileName = "universal-%(run)s.xml" % { "run" : runString }
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

    ## Prepare the joboutput tarball
    def prepareTarball( self, runString ):

        self._logger.info("Preparing the joboutput tarball" )

        # first prepare a folder to store them temporary
        destFolder = "universal-%(run)s" %{ "run": runString}
        shutil.rmtree( destFolder, True )
        os.mkdir( destFolder )

        # prepare the list of files we want to copy.
        listOfFiles = []
        listOfFiles.append( os.path.join( self._gearPath, self._gear_file ) )
        listOfFiles.append( self._steeringFileName )
        listOfFiles.append( self._logFileName )
        listOfFiles.append( self._configFile )

        for file in listOfFiles :
            shutil.copy( file , destFolder )

        # do the tarball
        self._tarballFileName = "universal-%(run)s.tar.gz" % {"run": runString }
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

        shutil.move( self._tarballFileName, localFolder )

    ## Cleanup after each run conversion
    def cleanup( self, runString ):

        self._logger.info( "Cleaning up the local pc" )

        # remove the log file and the steering file
        os.remove( self._steeringFileName )
        os.remove( self._logFileName )

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
        # universal-123456.tar.gz

        try :
            outputFilePath = self._configParser.get( "LOCAL", "LocalFolderConvertJoboutput" )
        except ConfigParser.NoOptionError :
            outputFilePath = "log"
        tarballFileName = "universal-%(run)s.tar.gz" % { "run": runString }
        if not os.access( os.path.join( outputFilePath, tarballFileName), os.R_OK ):
            raise MissingJoboutputFileError( tarballFileName )
        else:
            run, b, c, d, e, f = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, b, c, d, e, "OK"
