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
import datetime
from submitbase import SubmitBase
from error import *

## Submit native copy
#
# This class is responsible to copy native format file from the local PC
# to the SE on the GRID.
#
# The idea is that this class is called by the ./submit-native-copy.py script
# and it can be done in an almost automatic way by an external watchdog checking
# which are the newly added files
#
#
# @version $Id: submitnativecopy.py,v 1.2 2009-05-13 09:21:01 bulgheroni Exp $
# @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
#
class SubmitNativeCopy( SubmitBase ) :

    cvsVersion = "$Revision: 1.2 $"

    ## General configure
    #
    # This method is called by the constructor itself and it
    # mainly calls the base class methods to properly set up
    # all the bit and pieces
    def configure( self ):

        # do some basic configuration, like reading the config.cfg file
        SubmitBase.configure( self )

        # set up the logger according to the configuration file just read
        SubmitBase.configureLogger( self, "NativeCopy" )

        # print the welcome message!
        self._logger.log( 15, "**********************************************************" )
        self._logger.log( 15, "Started submit-native-copy" )
        self._logger.log( 15, "**********************************************************" )

        # now print the configuration to the log
        SubmitBase.logConfigurationFile( self )

        # now log the run list
        SubmitBase.logRunList( self )
