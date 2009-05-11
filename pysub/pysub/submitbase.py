## @package pysub
# Python Submission scripts
#
# In principle running the EUTelescope processors to get your data
# analyzed is not so difficult, but one has to remember which
# processors are best suited for a particular procedure, the correct
# order and reasonable parameter default values.
#
# Moreover if you think that Marlin and EUTelescope can also be
# executed on the GRID this is introducing an extra degree of
# complexity!
#
# To simplify all this, we are now providing together with EUTelescope
# itselft and the eutelHistogram package, also a set of very powerful
# submission script, written in python with the main goal to allow any
# user or even a sofware robot to submit his/her own job either
# locally or on the GRID without having to modify a single line of
# steering file!
#
# Of course this come at a price! Marlin and EUTelescope were born to
# be as flexibile as possible but all this flexibility is making an
# automatic job submission nearly impossible.
#
# pysub is a sort of compromise between automatization and
# flexibility. The idea is that all the main parameters of pysub
# itself, like the GRID configuration or the verbosity level can be
# adjusted changing their value in a simple ASCII configuration file
# (config/config.cfg). The execution parameters instead are given in
# some steering file used as templates by the script.
#
# Our advice is the following:
#  @li Everytime you start a new analysis
#  copy the config/config.cfg into a new file, say
#  config/MyConfig.cfg.
#  @li In the new configuration properly set all the GRID paramters in
#  the GRID section.
#  @li Again in the configuration file modify the section concerning the
#  steering template file name in order to point to your own modified
#  version, if this is required.
#  @li Tell pysub to use your own config file and not the default one
#  either specifing the configuration file as a command line option
#  (--config-file = CONFIGFILE ) or by setting the enviromental
#  variable SUBMIT_CONFIG.
#
# From a programmer point of view, the pysub package is designed like
# this:
# <h4> The base class </h4>
#  The base submitter class (SubmitBase) is a very simple class
# contaning just the overall scheme of all the operation that a
# submitter has to have. Those are <b>configure</b>, <b>execute</b>
# and <b>end</b>.
#
# <h4> The submitter class </h4>
#  This is a class inheriting from SubmitBase, like SubmitConverter or
# SubmitTest implementing the call back the base class.
#
# <h4> The submitter script </h4>
#  By it self the submitter class is not running, but just defining a
# functionality that has to be executed. This is the main job of the
# submitter script like submit-converter.py that is actually
# responsible to create a new instance of the submitter class and make
# it working.
#
# @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
# @version $Id: submitbase.py,v 1.5 2009-05-11 10:06:39 bulgheroni Exp $

from optparse import OptionParser
import ConfigParser
import logging
import time

## SubmitBase
# This is the base class for all submitters
#
# It is just defining a sort of template for the submit class that are
# inheriting from this.
#
# @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
# @version $Id: submitbase.py,v 1.5 2009-05-11 10:06:39 bulgheroni Exp $
#
class SubmitBase :

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
        self.configure()

        # initialize a summary ntuple
        self._summaryNTuple = [];

    ## The configure method
    #
    # This method, in particular its overload in the daughter class
    # should contain all the need instruction to read and parse the
    # configuration file.
    def configure( self ) :
        pass


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
