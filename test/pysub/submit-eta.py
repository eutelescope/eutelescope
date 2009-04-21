#!/usr/bin/env python

import sys
import string
import re
import os
import popen2
import glob

blue="\033[1;34m"
black="\033[0m"
red="\033[1;31m"
green="\033[1;32m"


def usage( commandname ):

    print green , "Usage:"
    print commandname, "[ OPTIONS ] -o | --output-basename basename  inputfile-name"
    print """
          List of OPTIONS

             inputfile-name            This is the path to a file LCIO file contaning the cluster collection to be
                                       used for ETA calculation

             -o or --output-basename   The output file name for the ETA database and for the corresponding histogram
                                       file will be built using this output base name according to the following
                                       naming convention:

                                       DB File     --> db/basename-eta-db.slcio
                                       Histo File  --> histo/basename-eta-histo.root

             -h or --help         Print this help message

             --only-generate      Simply genereate the steering file

    """
    print black


# here starts the main
print red, "Eta calculation job submitter" , black

# default options
optionGenerateOnly      = 0

# parse the argmuments
narg = 0
inputFileNameFullPath = ""
inputFileName         = ""
outputFileBaseName    = ""
dbFileNameFullPath    = ""
dbFileName            = ""
histoFileNameFullPath = ""
histoFileName         = ""
outputFileIndex       = -1 

goodArg = sys.argv[1:]

for i,arg in enumerate( goodArg ):
    narg = narg + 1

    # don't count twice the pedestal run
    if i == outputFileIndex:
        continue

    if arg == "-o" or arg == "--output-basename":
        outputFileBaseName = goodArg[ i + 1]
        outputFileIndex = i + 1

    elif arg == "--only-generate":
        optionGenerateOnly = 1

    elif arg == "-h" or arg == "--help":
        usage(sys.argv[0])
        sys.exit(0)

    else :

        # if I'm here it means that this is the input file
        inputFileNameFullPath = os.path.abspath( goodArg[ i ] )
        basepath, inputFileName = os.path.split( inputFileNameFullPath )
        if not os.path.exists( inputFileNameFullPath ) :
            print red, "The input file doesn't exist. " , black
            sys.exit( -1 )


if narg == 0:
    usage(sys.argv[0])
    sys.exit(0)

if outputFileIndex == -1:
    print red, "Output file basename not specified. Please use option -o to specify the output file name", black
    sys.exit(1)

# preparing the output file names
dbFileName         = outputFileBaseName + "-eta-db.slcio"
dbFileNameFullPath = "db/" + dbFileName

histoFileName          = outputFileBaseName + "-eta-histo"
histoFileNameFullPath  = "histo/" + histoFileName

print blue, "Generating the steering file...", black

#open the template steering file for reading
templateSteeringFile = open( "./template/eta-tmp.xml", "r")

# read the whole content of the template and put it into a string
actualSteeringString = templateSteeringFile.read()

actualSteeringString = actualSteeringString.replace( "@InputFile@", "%(input)s" % { "input": inputFileNameFullPath } )
actualSteeringString = actualSteeringString.replace( "@DBFile@", "%(db)s"  % { "db": dbFileNameFullPath } )
actualSteeringString = actualSteeringString.replace( "@HistoFile@", "%(histo)s"  % { "histo": histoFileNameFullPath } )

actualSteeringFile = open( "eta-%(output)s.xml" % {"output" : outputFileBaseName }, "w" )
actualSteeringFile.write( actualSteeringString )
actualSteeringFile.close()
templateSteeringFile.close()

if optionGenerateOnly == 0:

    # run Marlin in both all local and cpu local configuration
    print blue, "Running Marlin...",black

    # to properly get the return code I need to execute Marlin into a separete pipe
    logfile = open( "eta-%(output)s.log " %  {"output" : outputFileBaseName }, "w" )
    marlin = popen2.Popen4( "Marlin eta-%(output)s.xml" % {"output" : outputFileBaseName } )
    while marlin.poll() == -1:
        line = marlin.fromchild.readline()
        print line.strip()
        logfile.write( line )

    logfile.close()
    returnvalue = marlin.poll()

    if returnvalue != 0:
        print red, "Problem executing Marlin! (errno %(returnvalue)s)" % {"returnvalue":returnvalue }, black
        sys.exit( 1 )

    print blue, "Removing the steering and log files...", black
    command = "rm -r eta-%(output)s.*" % {"output" : outputFileBaseName }
    returnvalue = os.system( command )
    if returnvalue != 0:
        print red, "Problem removing the steering and log files...", black
