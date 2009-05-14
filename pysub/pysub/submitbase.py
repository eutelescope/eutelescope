from optparse import OptionParser
from error    import *
import ConfigParser
import logging
import time
import os
import sys
import popen2


## SubmitBase
# This is the base class for all submitters
#
# It is just defining a sort of template for the submit class that are
# inheriting from this.
#
# @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
# @version $Id: submitbase.py,v 1.11 2009-05-14 09:29:51 bulgheroni Exp $
#
class SubmitBase :

    cvsVersion = "$Revision: 1.11 $"

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
        self._logger = logging.getLogger( 'Base' )
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
            raise MissingConfigurationFileError( self._configFile )

        # if it exists, then I can read it!!!
        self._configParser = ConfigParser.SafeConfigParser()
        self._configParser.read( self._configFile )


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

        command = "lfc-ls %(folder)s" % { "folder": folder }
        lfc = popen2.Popen4( command )
        while lfc.poll() == -1:
            pass
        if lfc.poll() == 0 :
            self._logger.info( "Found folder %(folder)s" % { "folder": folder } )
            return True
        else :
            if interactive : 
                # ask the user if he wants to create the folder or not
                self._logger.error( "Unable to find folder %(folder)s" % { "folder": folder })
                if self.askYesNo( "Do you want to create it?" ) :
                    self._logger.info( "User deficed to create folder %(folder)s" % { "folder": folder } 
                    command = "lfc-mkdir -p %(folder)s " % { "folder" :folder } )
                    if os.system( command ) == 0 :
                        return True
                    else :
                        self._logger.critical("Unable to create folder %(folder)s" % { "folder": folder } )
                        raise MissingGRIDFolderError( folder )
                else:
                    raise MissingGRIDFolderError( folder )
            else :
                raise MissingGRIDFolderError( folder ) 
