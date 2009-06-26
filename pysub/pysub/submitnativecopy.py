import os
import sys
import shutil
import sha
import glob
import tarfile
import commands
import ConfigParser
import logging
import logging.handlers
import datetime
import time
import tempfile
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
# @version $Id: submitnativecopy.py,v 1.9 2009-06-26 09:45:14 jbehr Exp $
# @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
#
class SubmitNativeCopy( SubmitBase ) :

    ## Version
    # This number is used when printing out the version of the package.
    #
    # Static member.
    #
    cvsVersion = "$Revision: 1.9 $"

        ## Name
    # This is the namer of the class. It is used in flagging all the log entries
    # and preparing all the file
    #
    # Static member.
    #
    name = "nativecopy"

    ## General configure
    #
    # This method is called by the constructor itself and it
    # mainly calls the base class methods to properly set up
    # all the bit and pieces
    def configure( self ):

        # do some basic configuration, like reading the config.cfg file
        SubmitBase.configure( self )

        # set up the logger according to the configuration file just read
        SubmitBase.configureLogger( self, self.name )

        # print the welcome message!
        self._logger.log( 15, "**********************************************************" )
        self._logger.log( 15, "Started submit-%(name)s" % { "name" : self.name } )
        self._logger.log( 15, "**********************************************************" )

        # now print the configuration to the log
        SubmitBase.logConfigurationFile( self )

        # now log the run list
        SubmitBase.logRunList( self )

    ## Execute method
    #
    # This is the real method, responsible for the job submission
    # Actually this is just a sort of big switch calling the real submitter
    # depending of the execution mode
    #
    def execute( self ) : 

        # this is the real part 
        # convert all the arguments into an integer number
        self._runList = [] ;
        for i in self._args:
            try:
                self._runList.append( int( i  ) )

                # prepare also the entry for the summaryNTuple 
                # this is formatted in the way
                # runNumber, localFile, GRIDFile, Verifyon
                entry = i, "Unknown",  "Unknown",  "Unknown"
                self._summaryNTuple.append( entry )

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
                # check the presence of the input file
                self.checkInputFile( index, runString )

                # check if the file already exists on the GRID
                try :
                    self.checkAlreadyOnGRID( index, runString )
                except  FileAlreadyOnGRIDError, error:
                    message = "File %(file)s already on SE" % { "file":error._filename }
                    self._logger.warning( message )
                    self._logger.warning( "Skipping the copy and do directly the verification" )
                    run, input, output, verification = self._summaryNTuple[ index ]
                    self._summaryNTuple[ index ] = run, input, "Already", verification

                else :

                    # do the copy
                    self.copyRunOnGRID( index, runString )


                # make the verification
                self.verify( index, runString )

            except MissingInputFileError, error:
                message = "Missing input file %(file)s" % { "file": error._filename }
                self._logger.error( message )
                self._logger.error("Skipping to the next run ")
                run, input, output, verification = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, "Missing", output, verification

            except GRID_LCG_CRError, error:
                message = "Unable to copy %(file)s" % { "file": error._filename }
                self._logger.error( message )
                self._logger.error("Skipping to the next run ")
                run, input, output, verification = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, input, "Failed", "Failed"

            except GRID_LCG_CPError, error:
                message = "Unable to copy %(file)s" % { "file": error._filename }
                self._logger.error( message )
                self._logger.error("Skipping to the next run ")
                run, input, output, verification = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, input, output, "Failed"

            except VerificationError, error:
                message = "The verification failed!"
                self._logger.error( message )
                self._logger.error("Skipping to the next run ")
                run, input, output, verification = self._summaryNTuple[ index ]
                self._summaryNTuple[ index ] = run, input, output, "Failed"

    ## Check the input file
    #
    def checkInputFile( self, index, runString ) :
        # the input file should be something like:
        # native/run123456.raw

        self._logger.info( "Checking if the input run exists" )
        try :
            self._inputFilePath = self._configParser.get( "LOCAL", "LocalFolderNative" )
        except ConfigParser.NoOptionError :
            self._inputFilePath = "native"
            self._logger.warning( "Unable to find the local native folder in the configuration file" )
            self._logger.warning( "Using default one (%(path)s)" % { "path": self._inputFilePath } )
        inputFileName = "run%(run)s.raw" % { "run": runString }
        if not os.access( os.path.join( self._inputFilePath, inputFileName) , os.R_OK ):
            message = "Problem accessing the input file (%(file)s), trying next run" % {"file": inputFileName }
            self._logger.error( message )
            raise MissingInputFileError( inputFileName )
        else :
            run, input, output, verification = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, "OK", output, verification


    ## Check if a file already is already on the GRID
    #
    def checkAlreadyOnGRID( self, index, runString ):

        self._logger.info( "Checking if the input run is already on the GRID" )
        try :
            self._outputPathGRID    = self._configParser.get("GRID", "GRIDFolderNative" )
        except ConfigParser.NoOptionError:
            message = "Missing path from the configuration file"
            self._logger.critical( message )
            raise StopExecutionError( message )

        command = "lfc-ls %(outputPathGRID)s/run%(run)s.raw" % { "outputPathGRID": self._outputPathGRID, "run": runString }
        status, output = commands.getstatusoutput( command )

        if status == 0:
            self._logger.warning( "File %(outputPathGRID)s/run%(run)s.raw already on GRID"
                                  % { "outputPathGRID": self._outputPathGRID, "run": runString } )

            if self._configParser.get("General","Interactive" ):
                if self.askYesNo( "Would you like to remove it?  [y/n] " ):
                    self._logger.info( "User decided to remove %(outputPathGRID)s/run%(run)s.raw from the GRID"
                                       % { "outputPathGRID": self._outputPathGRID, "run": runString } )
                    command = "lcg-del -a lfn:%(outputPathGRID)s/run%(run)s.raw" % { "outputPathGRID": self._outputPathGRID, "run": runString }
                    os.system( command )
                else :
                    raise FileAlreadyOnGRIDError( "%(outputPathGRID)s/run%(run)s.raw from the GRID"
                                                  % { "outputPathGRID": self._outputPathGRID, "run": runString } )
            else :
                raise FileAlreadyOnGRIDError( "%(outputPathGRID)s/run%(run)s.raw from the GRID"
                                              % { "outputPathGRID": self._outputPathGRID, "run": runString } )

    ## Copy on GRID
    #
    def copyRunOnGRID( self, index, runString ) :

        self._logger.info( "Started copy to the GRID" )
        try :
            self._gridse = self._configParser.get( "GRID", "GRIDSE" )
        except ConfigParser.NoOptionError :
            self._logger.warning("Unable to get the GRIDSE from the configuration file.")
            self._gridse = "dcache-se-desy.desy.de"
            self._logger.warning("Using the default one (%(gridse)s)." % { "gridse": self._gridse } )

        baseCommand = "lcg-cr "
        if self._option.verbose:
            baseCommand = baseCommand + "-v "

        # check if we have a specific destination
        if self._gridse.startswith( "srm://"):
             self._gridse = "-d " + self._gridse.rstrip("/") + "/run%(runString)s.raw" % { "runString": runString }
        else:
             self._gridse = ""
 
 
        command = "%(base)s %(gridse)s -l lfn:%(gridFolder)s/run%(runString)s.raw file:%(inputFolder)s/run%(runString)s.raw" \
            % {  "base": baseCommand, "gridse": self._gridse, "gridFolder" : self._outputPathGRID,   "runString": runString, 
                 "inputFolder" : self._inputFilePath }
        if os.system( command ) != 0:
            raise GRID_LCG_CRError( "lfn:%(gridFolder)s/run%(runString)s.raw" % { "gridFolder" : self._outputPathGRID, "runString": runString } )
        else :
            run, input, output, verification = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, input, "OK", verification

    ## Verify the copy
    #
    def verify( self, index, runString ):
        self._logger.info( "Doing verification..." )

        # calculate the sha1sum of the local copy
        localCopy = open( "%(localFolder)s/run%(run)s.raw" % { "localFolder" : self._inputFilePath, "run": runString }, "r" ).read()
        localCopyHash = sha.new( localCopy ).hexdigest()
        message = "Local copy hash %(hash)s" % { "hash": localCopyHash }
        self._logger.log( 15, message )

        # now copy back the remote file
        # if so we need to copy back the file
        baseCommand = "lcg-cp "
        if self._option.verbose:
            baseCommand = baseCommand + "-v "

        tempFolder = tempfile.mkdtemp( prefix =self.name )
        command = "%(base)s lfn:%(gridFolder)s/run%(run)s.raw file:%(localFolder)s/run%(run)s-test.raw" % \
            { "base": baseCommand, "gridFolder": self._outputPathGRID, "localFolder": tempFolder, "run" : runString }

        if os.system( command ) != 0 :
            raise GRID_LCG_CRError( "lfn:%(gridFolder)s/run%(run)s.raw" % \
                                        { "gridFolder": self._inputFilePath, "run" : runString } )

        # calculate the sha1sum of the remote file 
        remoteCopy = open( "%(localFolder)s/run%(run)s-test.raw" % { "localFolder": tempFolder, "run" : runString }, "r" ).read()
        remoteCopyHash = sha.new( remoteCopy ).hexdigest()
        message = "Remote copy hash %(hash)s" % { "hash": remoteCopyHash }
        self._logger.log( 15, message )

        if remoteCopyHash == localCopyHash:
            run, input, output, verification = self._summaryNTuple[ index ]
            self._summaryNTuple[ index ] = run, input, output, "OK"
        else :
            raise VerificationError

        # if successful remove the test file
        shutil.rmtree( tempFolder )

    ## Final method
    #
    def end( self ) :

        self._logger.info( "Finished loop on input runs" )

        # print summary
        self.logSummary()

        self._timeEnd = time.time()
        message = "Submission completed in %(time)d seconds" % {
            "time": self._timeEnd - self._timeBegin }
        self._logger.info( message )



    ## Log the summary
    #
    def logSummary( self ):

        if len( self._summaryNTuple) == 0 :
            pass

        else:
            self._logger.info( "" ) 
            self._logger.info( "== SUBMISSION SUMMARY =======================================================" )
            message = "| %(run)13s | %(inputFileStatus)17s | %(outputFileStatus)17s | %(verification)17s |" \
                % { "run": "Run", "inputFileStatus" : "Input File",
                    "outputFileStatus": "Output File", "verification": "Verification" }
            self._logger.info( message )
            for run, input, output, verification in self._summaryNTuple :
                message = "| %(run)13s | %(inputFileStatus)17s | %(outputFileStatus)17s | %(verification)17s |" \
                    % { "run": run, "inputFileStatus" : input, "outputFileStatus": output, "verification": verification }
                self._logger.info( message )
            self._logger.info("=============================================================================" )
            self._logger.info( "" )
