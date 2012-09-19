from optparse import OptionParser
from error    import *
import ConfigParser
import logging
import time
import os
import sys
import popen2
import commands

## SubmitBase
# This is the base class for all submitters
#
# It is just defining a sort of template for the submit class that are
# inheriting from this.
#
# @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
# @version $Id: submitbase.py,v 1.29 2009-08-01 15:14:15 bulgheroni Exp $
#
class SubmitBase :

    ## Version
    # This number is used when printing out the version of the package.
    #
    # Static member.
    #
    cvsVersion = "$Revision: 1.29 $"

    ## Name
    # This is the namer of the class. It is used in flagging all the log entries
    # and preparing all the file
    #
    # Static member.
    #
    name = "base"

    ## The base class constructor
    #
    # Unless the daughter class has its own constructor, this is the
    # only one is called whenever a new submitter is created.
    #
    # It creates a very simple console logger with medium verbosity
    # level and then call the configure method.
    def __init__( self, optionParser ):

        # assign the option from outside to the local optionparser
        self._optionParser = optionParser
        self._option, self._args = self._optionParser.parse_args()

        # create a logger, for the time being the logger is a simple
        # console handler, after the configuration is read, it will be
        # properly configured
        self._logger = logging.getLogger( self.name )
        logging.addLevelName(15, "ALL" )
        self._logger.setLevel( logging.INFO )
        self._consoleHandler = logging.StreamHandler()
        self._consoleHandler.setLevel( logging.INFO )
        self._logger.addHandler( self._consoleHandler )

        # defining the t_0 of the submitter
        self._timeBegin = time.time()

        # load the configuration
        try:
            self.configure()
        except MissingConfigurationFileError, error:
            logging.critical( "Configuration file %(cfg)s doesn't exist!" % { "cfg": error._filename } )
            sys.exit( 1 )

        # initialize a summary ntuple
        self._summaryNTuple = [];

        # initialize the grid job ntuple
        self._gridJobNTuple = [];


    ## The configure method
    #
    # This method, in particular its overload in the daughter class
    # should contain all the need instruction to read and parse the
    # configuration file.
    #
    def configure( self ) :

        # get the config file 
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

        # get the range to process 
        self._eventRange = ""
        if ( self._option.event_range == None ) :
            # not from the option, check from environ
            self._eventRange_begin = 0
            self._eventRange_end   = 10000000
            self._eventRange = 10000000
        else:
            self._eventRange_string = self._option.event_range.split("-")
            self._eventRange_begin  = self._eventRange_string[0]
            self._eventRange_end    = self._eventRange_string[1]

        # before proceeding check if the configuration file
        # really exists!
        if not os.path.exists( self._configFile ):
            raise MissingConfigurationFileError( self._configFile )

        # set up dictionary to be used for defaults in the config file
        # populate defaults e.g. with environment variables
        config_defaults = {
            'env_eutelescope_path' : os.environ['EUTELESCOPE'],
            'env_home_path' : os.environ['HOME']
            }

        # if it exists, then I can read it!!!
        self._configParser = ConfigParser.SafeConfigParser(config_defaults)
        self._configParser.read( self._configFile )

        # set a couple of usefull information from the configuration file
        self._isInteractive = False;
        self._isForceYes    = False;
        try :
            self._isInteractive = self._configParser.getboolean( "General", "Interactive" )
            self._isForceYes    = self._configParser.getboolean( "General", "ForceYes" )
        except ConfigParser.NoOptionError:
            self._logger.debug( "Unable to find interactive keys in the configuration file" )


    ## Logger configurator
    #
    # This is another configure method that is called by the main configure
    # to properly set up the logging system.
    # Before calling this method only a very simple logging system is available
    #
    def configureLogger( self, name ):

        if not self._configParser.has_section( "Logger" ) :
            logging.warning( "No Logger section in the configuration file, using default logging" )
            return

        # rename the logger according to this module
        self._logger = logging.getLogger( name )
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


    ## The execute method
    #
    # This is the real method in which the job is finally
    # submitted. It can be made by several other methods depending on
    # the complexity of the job and of the execution mode
    def execute( self ) :
        pass

    ## The end metod
    #
    # This is the last method called in the sequence. It is doing some
    # housekeeping leaving the system clean and print the last line in
    # the log file
    def end( self ) :

        self._logger.info( "Finished loop on input runs" )
        self.logSummary()

        self._timeEnd = time.time()
        message = "Submission completed in %(time)d seconds" % {
            "time": self._timeEnd - self._timeBegin }
        self._logger.info( message )

    ## Ask yes or no
    #
    def askYesNo( self, prompt, complaint= "Yes or no please", retries = 4 ):

        try :
            self._isForceYes

        except AttributeError:

            pass

        else:

            if self._isForceYes :
                return True

        while True:
            ok = raw_input( prompt )
            ok = ok.lower()
            if ok in ['y','yes','yea','yeap']  : return True
            if ok in ['n','no',"nope" ] : return False
            retries = retries - 1
            if retries < 0 : raise IOError, "Wrong answer"
            print complaint

    ## Log summary
    #
    # This method is used to log the summary of the submission
    # Remember that the ntuple is made in this way:
    # run; inputFileStatus; marlinStatus ; outputFileStatus ; histoFileStatus ; joboutputFileStatus
    #
    def logSummary( self ):

        if len( self._summaryNTuple) == 0 :
            pass

        else:
            self._logger.info( "" )
            self._logger.info( "== SUBMISSION SUMMARY =======================================================" )
            message = "| %(run)6s | %(inputFileStatus)10s | %(marlinStatus)10s | %(outputFileStatus)11s | %(histoFileStatus)11s | %(joboutputFileStatus)10s |" \
                % { "run": "Run", "inputFileStatus" : "Input File", "marlinStatus": "Marlin",
                    "outputFileStatus": "Output File", "histoFileStatus": "Histo File", "joboutputFileStatus": "Joboutput" }
            self._logger.info( message )
            for entry in self._summaryNTuple :
                run, inputFileStatus, marlinStatus, outputFileStatus, histoFileStatus, joboutputFileStatus = entry
                message = "| %(run)6s | %(inputFileStatus)10s | %(marlinStatus)10s | %(outputFileStatus)11s | %(histoFileStatus)11s | %(joboutputFileStatus)10s |" \
                    % { "run": run, "inputFileStatus" : inputFileStatus, "marlinStatus":marlinStatus,
                        "outputFileStatus": outputFileStatus, "histoFileStatus": histoFileStatus, "joboutputFileStatus":joboutputFileStatus }
                self._logger.info( message )
            self._logger.info("=============================================================================" )
            self._logger.info( "" )


    ## Verify the existence of a file on the GRID
    #
    #
    def checkGRIDFile( self, filename  ):
        self._logger.info( "Checking file %(file)s" % { "file": filename } )
        command = "lfc-ls %(file)s" % { "file": filename }
        lfc = popen2.Popen4( command )
        while lfc.poll() == -1:
            pass
        if lfc.poll() == 0:
            self._logger.warning( "%(file)s already on GRID" % { "file": filename } )
            if self._configParser.get("General","Interactive" ):
                if self.askYesNo( "Would you like to remove it?  [y/n] " ):
                    self._logger.info( "User decided to remove %(file)s " % { "file": filename } )
                    command = "lcg-del -a lfn:%(file)s" % { "file": filename }
                    self._logger.info( command )
                    os.system( command )
                    return True
                else :
                    raise OutputFileAlreadyOnGRIDError( "%(file)s on the GRID" % { "file": filename } )
            else :
                raise OutputFileAlreadyOnGRIDError( "%(file)s" % {"file": filename } )


    ## Verify the existence of a GRID folder
    #
    # @return True if the folder exists
    # @throw MissingGRIDFolder in case it doesn't
    #
    def checkGRIDFolder( self, folder ) :

        # first check if it is in interactive mode or it is just a ROBOT
        try :
            interactive = self._configParser.getboolean( "General", "Interactive" )
        except ConfigParser.NoOptionError :
            interactive = True

        command = "lfc-ls %(folder)s > /dev/null 2>&1" % { "folder": folder }
        if os.system( command ) == 0:
            self._logger.log(15, "Found folder %(folder)s" % { "folder": folder } )
            return True
        else :
            if interactive :
                # ask the user if he wants to create the folder or not
                self._logger.error( "Unable to find folder %(folder)s" % { "folder": folder })
                if self.askYesNo( "Do you want to create it? " ) :
                    self._logger.info( "User decided to create folder %(folder)s" % { "folder": folder } )
                    command = "lfc-mkdir -p %(folder)s " % { "folder" :folder }
                    if os.system( command ) == 0 :
                        return True
                    else :
                        self._logger.critical("Unable to create folder %(folder)s" % { "folder": folder } )
                        raise MissingGRIDFolderError( folder )
                else:
                    raise MissingGRIDFolderError( folder )
            else :
                raise MissingGRIDFolderError( folder )

    ## Generate JDL file
    #
    # This method is called to generate a JDL file
    #
    def generateJDLFile( self, index, runString, *args ):
        message = "Generating the JDL file (%(name)s-%(run)s.jdl)" % { "name": self.name,"run": runString }
        self._logger.info( message )
        try :
            jdlTemplate = self._configParser.get("GRID", "GRIDJDLTemplate" )
        except ConfigParser.NoOptionError:
            jdlTemplate = "grid/jdl-tmp.jdl"
            if os.path.exists ( jdlTemplate ) :
                message = "Using JDL template (%(file)s) " % { "file": jdlTemplate }
                self._logger.info( message )
            else:
                message = "Unable to find a valid JDL template"
                self._logger.critical( message )
                raise StopExecutionError( message )

        jdlTemplateString = open( jdlTemplate, "r" ).read()
        jdlActualString = jdlTemplateString

        # modify the executable name
        jdlActualString = jdlActualString.replace( "@Executable@", "%(name)s-%(run)s.sh" % { "name": self.name, "run": runString } )

        # the gear file!
        # try to read it from the config file, then from the command line option
        try :
            self._gearPath = self._configParser.get( "LOCAL", "LocalFolderGear" )
        except ConfigParser.NoOptionError :
            self._gearPath = ""
        jdlActualString = jdlActualString.replace("@GearPath@", self._gearPath )

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

        jdlActualString = jdlActualString.replace( "@GearFile@", "%(path)s/%(gear)s"
                                                   % { "path": self._gearPath, "gear": self._gear_file } )

        # replace the steering file
        jdlActualString = jdlActualString.replace( "@SteeringFile@", "%(name)s-%(run)s.xml" % { "name": self.name, "run" : runString } )

        # replace the GRIDLib
        try :
            gridLibraryTarball = self._configParser.get( "GRID", "GRIDLibraryTarball" )
            gridLibraryTarballPath = self._configParser.get( "GRID", "GRIDLibraryTarballPath" )
        except ConfigParser.NoOptionError :
            message = "GRID library tarball unavailable!"
            self._logger.critical( message )
            raise StopExecutionError( message )

        # guess if the library is locally on the computer, or if it is already on the GRID
        if gridLibraryTarballPath.startswith( "lfn:" ):
            # it's on the storage element, check if it is there:
            command = "lfc-ls  %(path)s/%(file)s" % { "path": gridLibraryTarballPath.lstrip("lfn:"), "file": gridLibraryTarball }
            status, output = commands.getstatusoutput( command )
            if self._option.verbose:
                for line in output.splitlines():
                    self._logger.info( line.strip() )
            if status != 0:
                raise MissingLibraryFileError ( "%(path)s/%(file)s" % { "path": gridLibraryTarballPath, "file": gridLibraryTarball } )

            jdlActualString = jdlActualString.replace( "@GRIDLibraryTarball@", "" )

        else :

            jdlActualString = jdlActualString.replace( "@GRIDLibraryTarball@", ", \"%(path)s/%(file)s\" " %
                                                       { "path": gridLibraryTarballPath, "file":gridLibraryTarball } )

        # replace the histoinfo file
        try:
            histoinfo = self._configParser.get("General","Histoinfo")
        except ConfigParser.NoOptionError :
            self._logger.debug( "No histoinfo file in the configuration file" )
            histoinfo = ""

        try:
            histoinfoPath = self._configParser.get("LOCAL", "LocalFolderHistoinfo")
        except ConfigParser.NoOptionError :
            self._logger.debug( "No histoinfo path file in the configuration file" )
            histoinfoPath = "./"

        # check if it exists:
        if os.access( os.path.join( histoinfoPath, histoinfo), os.R_OK ) :
            jdlActualString = jdlActualString.replace( "@HistoInfo@", os.path.join( histoinfoPath, histoinfo) )
        else:
            jdlActualString = jdlActualString.replace( "@HistoInfo@", "" )

        # replace the VO
        try:
            vo = self._configParser.get( "GRID" , "GRIDVO" )
        except ConfigParser.NoOptionError :
            self._logger.warning( "Unable to find the GRIDVO. Using ilc" )
            vo = "ilc"
        jdlActualString = jdlActualString.replace( "@GRIDVO@", vo )

        # replace the ILCSoftVestion
        try :
            ilcsoftVersion = self._configParser.get( "GRID" , "GRIDILCSoftVersion" )
        except ConfigParser.NoOptionError :
            self._logger.warning( "Unable to find the GRIDILCSoftVersion. Using v01-06" )
            ilcsoftVersion = "v01-06"
        jdlActualString = jdlActualString.replace( "@GRIDILCSoftVersion@", ilcsoftVersion )

        # replace any other additional arguments:
        for arg in args:
            jdlActualString = jdlActualString.replace( "@Others@", ", \"%(arg)s\" @Others@" % { "arg": arg } )
        # remove spare others
        jdlActualString = jdlActualString.replace( "@Others@", "" )

        self._jdlFilename = "%(name)s-%(run)s.jdl" % { "name": self.name, "run": runString }
        jdlActualFile = open( self._jdlFilename, "w" )
        jdlActualFile.write( jdlActualString )
        jdlActualFile.close()


    def checkDelegationProxy( self ) :
        self._logger.info( "Checking for an existing delegation proxy")
        delegation = popen2.Popen4( "glite-wms-job-info -d $USER" , -1 )
#        delegation = popen2.Popen4(  "glite-wms-job-delegate-proxy -d $USER", -1  )
        rtr = delegation.wait()
        output = delegation.fromchild.read()
        for line in output.splitlines():
            self._logger.log( 15, line.strip() )

        if rtr != 0:
            self._logger.info( "Generating a new delegate proxy" )
            # prepare a delegation proxy
            delegation = popen2.Popen4(  "glite-wms-job-delegate-proxy -d $USER", -1  )
            rtr = delegation.wait()
            output = delegation.fromchild.read()
            for line in output.splitlines():
                self._logger.log( 15, line.strip() )
            if rtr != 0:
                self._logger.warning( "Unable to do job proxy delegation. Using auto-delegation" )
                self._jobDelegation  = " -a "
            else :
                self._jobDelegation  = " -d $USER "
        else:
            self._jobDelegation = " -d $USER "

