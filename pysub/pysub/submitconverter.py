import os
import sys
import shutil
import glob
import tarfile
import popen2
import ConfigParser
import logging
import logging.handlers
from submitbase import SubmitBase



## Submit conversion jobs
#
# This is the class responsible of executing job conversions
class SubmitConverter( SubmitBase ) :
    """ This is a test submit class to understand how it works! """

    ## General configure
    #
    # This method is called by the constructor itself and performs all the
    # basic configurations from the configuration file.
    # In particular is calling the configureLogger method to start the logging
    # system in its full glory
    def configure( self ):

        # first of all we need to see if we have a user defined configuration 
        # file different from the template/config.cfg
        # 
        # The configuration file can be passed either as a command line option
        # or via an enviromental variable
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
        self._logger.log(1, "Configuration file " )
        for section in self._configParser.sections():
            message = "Section: %(section)s" % {'section':section}
            self._logger.log(1, message) 
            for option in self._configParser.options(section):
                message = " %(option)s = %(value)s " % { "option":option, "value":self._configParser.get( section, option) }
                self._logger.log(1,  message )


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

        # now check if the user overwrites this setting in the options
        if  self._option.force_keep_input == True:
            self._keepInput = True
            self._logger.info( "User forces to keep the input files" )
        else: 
            self._keepInput = False
            self._logger.info( "User forces to remove the input files" )

        if  self._option.force_keep_output == True:
            self._keepOutput = True
            self._logger.info( "User forces to keep the output files" )
        else:
            self._keepOutput = False
            self._logger.info( "User forces to remove the output files" )

        # now in case the user wants to interact and the situation is dangerous
        # i.e. removing files, ask confirmation
        try :
            interactive = self._configParser.getboolean( "General", "Interactive" )
        except ConfigParser.NoOptionError :
            interactive = True

        if self._keepInput == False and interactive :
            goodAnswer = False
            while goodAnswer == False:
                self._logger.warning("This script is going to delete the input file(s) when finished." )
                self._logger.warning("Are you sure to continue? [y/n]")
                answer = raw_input( "--> ").lower()
                if answer != "y" and answer != "yes" and answer != "n" and answer != "no":
                    self._logger.warning("Invalid answer, please type y or n")
                    answer = raw_input( "--> " ).lower()
                elif answer == "y" or answer == "yes":
                    goodAnswer = True
                    self._logger.info("User decide to continue removing input files")
                elif answer == "n" or answer == "no":
                    goodAnswer = True
                    self._logger.info("Aborted by user")
                    sys.exit( 4 )

        if self._keepOutput == False and interactive :
            goodAnswer = False
            while goodAnswer == False:
                self._logger.warning("This script is going to delete the output file(s) when finished." )
                self._logger.warning("Are you sure to continue? [y/n]")
                answer = raw_input( "--> ").lower()
                if answer != "y" and answer != "yes" and answer != "n" and answer != "no":
                    self._logger.warning("Invalid answer, please type y or n")
                    answer = raw_input( "--> " ).lower()
                elif answer == "y" or answer == "yes":
                    goodAnswer = True
                    self._logger.info("User decide to continue removing output files")
                elif answer == "n" or answer == "no":
                    goodAnswer = True
                    self._logger.info("Aborted by user")
                    sys.exit( 4 )

    ## Logger configurator
    #
    # This is another configure method that is called by the main configure
    # to properly set up the logging system.
    # Before calling this method only a very simple logging system is available
    def configureLogger( self ):

        if not self._configParser.has_section( "Logger" ) :
            logging.warning( "No Logger section in the configuration file, using default logging" )
            return

        # set the global logging level
        self._logger.setLevel( self._configParser.getint( "Logger", "GlobalLoggerLevel" ))

        # configuring the console handler
        self._logger.removeHandler( self._consoleHandler )
        if self._configParser.getboolean( "Logger", "ConsoleHandler" ):
            self._consoleHandler = logging.StreamHandler()
            self._consoleHandler.setLevel( self._configParser.getint( "Logger", "ConsoleHandlerLevel" ) )
            self._consoleHandler.setFormatter( logging.Formatter("%(asctime)s - %(levelname)s - %(message)s") )
            self._logger.addHandler( self._consoleHandler )

        # configuring the time rotating file handler
        if self._configParser.getboolean( "Logger", "RotatingFileHandler" ):
            timeRotatingHandler = logging.handlers.RotatingFileHandler(
                filename = self._configParser.get( "Logger", "RotatingFileHandlerFileName" ),
                maxBytes = self._configParser.getint( "Logger", "RotatingFileHandlerSize" ),
                backupCount = 10)
            timeRotatingHandler.setLevel( self._configParser.getint( "Logger", "RotatingFileHandlerLevel" ) )
            timeRotatingHandler.setFormatter( logging.Formatter("%(asctime)s - %(levelname)s - %(message)s") )
            self._logger.addHandler( timeRotatingHandler )

        self._logger.log( 1, "**********************************************************" )
        self._logger.log( 1, "Started submit-converter" )


    ## Execute method
    #
    # This is the real method, responsible for the job submission
    # Actually this is just a sort of big switch calling the real submitter
    # depending of the execution mode
    def execute( self ) :

        # this is the real part
        # convert all the argumets into a run string compatible with the file name convention
        self._runList = [];
        for i in self._args:
            try:
                self._runList.append( int( i  ) )
            except ValueError:
                message = "Invalid run number %(i)s" % { "i": i }
                self._logger.critical( message )
                sys.exit( 1 )


        # now do something different depending on the execution option
        if self._option.execution == "all-grid" :
            self.allGridSubmission()

        elif self._option.execution == "all-local" :
            self.allLocalSubmission()

        elif self._option.execution == "cpu-local":
            self.cpuLocalSubmission()

        elif self._option.execution == "only-generate":
            self.onlyGenerateSubmission()


    ## This is the real all-grid submitter
    def allGridSubmission( self ) :
        pass

    ## CPU Local submitter
    #
    # This methods represents the sequence of commands that should be done while
    # submitting jobs on the local computer but using remote data
    def cpuLocalSubmission( self ):

        for run in self._runList :
            runString = "%(#)06d" % { "#" : run }

            # get the input file from the GRID
            if self.getRunFromGRID( runString ) != 0:
                self._logger.error("Problem copying input file from the GRID, trying next!" )
                continue

            self._logger.info("Input file successfully copied from the GRID")

            #  generate the steering file
            self._steeringFileName = self.generateSteeringFile( runString )

            # prepare a separate file for logging the output of Marlin
            self._logFileName = "universal-%(run)s.log" % { "run" : runString }

            # run marlin takes as inputs the steering file and the log file
            if self.runMarlin( ) != 0 :
                self._logger.error( "Problem executing Marlin, skipping to the next run" )
                continue

            # advice the user that Marlin is over
            self._logger.info( "Marlin finished successfully") 

            # prepare a tarbal for the records
            self.prepareTarball( runString )

            # copy the output LCIO file to the GRID
            if self.putRunOnGRID( runString ) != 0:
                self._logger.error("Problem copying output LCIO file to the GRID, trying next!" )
                self._logger.error("Input and output files will be kept" )
                continue
            self._logger.info("Ouput file successfully copied from the GRID")

            # copy the joboutput file to the GRID
            if self.putJoboutputOnGRID( runString ) != 0:
                self._logger.error("Problem copying joboutput file to the GRID, trying next!" )
                self._logger.error("Input and output files will be kept" )
                continue
            self._logger.info("Jobouput file successfully copied from the GRID")

            # clean up the local pc
            self.cleanup( runString )


    ## Local submitter
    #
    # This methods represents the sequence of commands that should be done while
    # submitting jobs on the local computer.
    def allLocalSubmission( self ):

        for run in self._runList :
            runString = "%(#)06d" % { "#" : run }

            # first generate the steering file
            self._steeringFileName = self.generateSteeringFile( runString )

            # prepare a separate file for logging the output of Marlin
            self._logFileName = "universal-%(run)s.log" % { "run" : runString }

            # run marlin takes as inputs the steering file and the log file
            if self.runMarlin( ) != 0 :
                self._logger.error( "Problem executing Marlin, skipping to the next run" )
                continue

            # advice the user that Marlin is over
            self._logger.info( "Marlin finished successfully") 


            # prepare a tarbal for the records
            self.prepareTarball( runString )

            # clean up the local pc
            self.cleanup( runString )

    ## Generate only submitter
    #
    # This methods is responsibile of dry-run with only steering file
    # generation
    def onlyGenerateSubmission( self ):
        for run in self._runList :
            runString = "%(#)06d" % { "#" : run }

            # just need to generate the steering file
            self.generateSteeringFile( runString )

    ## Get the input run from the GRID
    def getRunFromGRID(self, runString ) :

        self._logger.info(  "Getting the input file from the GRID" )

        try :
            gridNativePath = self._configParser.get( "GRID", "GRIDFolderNative" )
        except ConfigParser.NoOptionError :
            self._logger.critical( "GRIDFolderNative missing in the configuration file. Quitting." )
            sys.exit( 5 )

        localPath = "$PWD/native"
        command = "lcg-cp -v lfn:%(gridNativePath)s/run%(run)s.raw file:%(localPath)s/run%(run)s.raw" %  \
                  { "gridNativePath" : gridNativePath, "run": runString, "localPath": localPath }
        return  os.system( command )

    ## Put the output run to the GRID
    def putRunOnGRID( self, runString ): 

        self._logger.info(  "Putting the LCIO file to the GRID" )

        try :
            gridLcioRawPath = self._configParser.get( "GRID", "GRIDFolderLcioRaw")
        except ConfigParser.NoOptionError :
            self._logger.critical( "GRIDFolderLcioRaw missing in the configuration file. Quitting." )
            sys.exit( 5 )

        localPath = "$PWD/lcio-raw"
        command = "lcg-cr -v -l lfn:%(gridFolder)s/run%(run)s.slcio file:%(localFolder)s/run%(run)s.slcio" % \
            { "gridFolder": gridLcioRawPath, "localFolder": localPath, "run" : runString }
        return os.system( command )

    ## Put the joboutput to the GRID
    def putJoboutputOnGRID( self, runString ):

        self._logger.info(  "Putting the joboutput file to the GRID" )

        try :
            gridFolder = self._configParser.get( "GRID", "GRIDFolderConvertJoboutput")
        except ConfigParser.NoOptionError :
            self._logger.critical( "GRIDFolderConvertJoboutput missing in the configuration file. Quitting." )
            sys.exit( 5 )

        localPath = "$PWD"
        command = "lcg-cr -v -l lfn:%(gridFolder)s/universal-%(run)s.tar.gz file:%(localFolder)s/universal-%(run)s.tar.gz" % \
            { "gridFolder": gridFolder, "localFolder": localPath, "run" : runString }

        return os.system( command )

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
            message ="The steering template %(tmp)s is missing. Quitting!" % { "tmp":steeringFileTemplate }
            self._logger.critical( message )
            sys.exit( 2 )

        # open the template for reading
        templateSteeringFile = open( steeringFileTemplate , "r")

        # read the whole content
        templateSteeringString = templateSteeringFile.read()

        # make all the changes
        actualSteeringString = templateSteeringString
        actualSteeringString = actualSteeringString.replace("@RunNumber@", runString )

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
    def runMarlin( self ) :

        self._logger.info( "Running Marlin" )

        # first check that the gear file exists
        if not os.path.exists( self._gear_file ) :
            message = "The GEAR file %(gear)s doesn't exist!" % {"gear" : self._gear_file } 
            self._logger.error( message )
            return 3

        # do some tricks for having the logfile
        logFile = open( self._logFileName, "w")
        marlin  = popen2.Popen4( "Marlin %(steer)s" % { "steer": self._steeringFileName } )
        while marlin.poll() == -1:
            line = marlin.fromchild.readline()
            print line.strip()
            logFile.write( line )

        logFile.close()
        return marlin.poll()

    ## Prepare the joboutput tarball
    def prepareTarball( self, runString ):

        self._logger.info("Preparing the joboutput tarball" )

        # first prepare a folder to store them temporary
        destFolder = "universal-%(run)s" %{ "run": runString}
        shutil.rmtree( destFolder, True )
        os.mkdir( destFolder )

        # prepare the list of files we want to copy.
        listOfFiles = []
        listOfFiles.append( self._gear_file )
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


    ## Cleanup after each run conversion
    def cleanup( self, runString ):

        self._logger.info( "Cleaning up the local pc" )

        # copy the tarball in the log folder
        shutil.move( self._tarballFileName, "log/")

        # remove the log file and the steering file
        os.remove( self._steeringFileName )
        os.remove( self._logFileName )

        # remove the input and output file
        inputFile  = "native/run%(run)s.raw" % { "run" : runString }
        outputFile = "lcio-raw/run%(run)s.slcio" % { "run" : runString }

        # but before removing check if this is needed
        if self._keepInput == False :
            os.remove( inputFile )

        if self._keepOutput == False :
            os.remove( outputFile )
