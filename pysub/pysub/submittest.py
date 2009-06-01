import os
from submitbase import SubmitBase
import ConfigParser
import logging
import logging.handlers
import time

## SubmitTest
# @class SubmitTest
#
# This is a test submit class to understand how it works!
class SubmitTest( SubmitBase ) :
    """ This is a test submit class to understand how it works! """

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
            self._consoleHandler.setFormatter( logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
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

        # do some verification on the command line options
        if self._option.plus and self._option.minus :
            self._optionParser.error( "Options --minus and --plus are mutually exclusive" )

        # before proceeding check if the configuration file 
        # really exists!
        if not os.path.exists( self._configFile ):
            logging.critical( "Configuration file %(cfg)s doesn't exist!" % { "cfg": self._configFile } )

        # if it exists, then I can read it!!!
        self._configParser = ConfigParser.SafeConfigParser()
        self._configParser.read( self._configFile )

        # now we can properly set the logger.
        self.configureLogger()
        
    def execute( self ) :
        print "Working very hard..."
        time.sleep( 1 )

