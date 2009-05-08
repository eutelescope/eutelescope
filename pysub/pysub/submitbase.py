from optparse import OptionParser
import ConfigParser
import logging

## SubmitBase
#
# This is the base class for all submitter
class SubmitBase :
    """ The base class for all submitters  """


    ## The base class constructor
    def __init__( self, optionParser ):
        """ The consustructor """

        # assign the option from outside to the local optionparser
        self._optionParser = optionParser
        self._option, self._args = self._optionParser.parse_args()

        # create a logger, for the time being the logger is a simple 
        # console handler, after the configuration is read, it will be
        # properly configured
        self._logger = logging.getLogger( 'Logger' )
        self._logger.setLevel( logging.INFO )
        self._consoleHandler = logging.StreamHandler()
        self._consoleHandler.setLevel( logging.INFO )
        self._logger.addHandler( self._consoleHandler )

        # load the configuration
        self.configure()


    ## The configure method
    def configure( self ) :
        """ This is the place where the configuration file is loaded """
        pass


    ## The execute method
    def execute( self ) :
        """ This is the place where the real job is done """
        pass

    ## The end metod
    def end( self ) :
        """ This is the place where the cleaning up is done """
        pass
