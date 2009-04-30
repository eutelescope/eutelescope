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

             inputfile-names           This is the path to one or more LCIO files contaning the cluster collection to be
                                       converted into space points. The results of all input files will be merged into a single
                                       output hit file

             -o or --output-basename   The output file name for the hit file and the corresponding histogram file
                                       according to the following naming convention

                                       hit File     --> restults/basename-hit.slcio
                                       Histo File   --> histo/basename-hit-histo.root
                                       
             -g or --gear-file         Specify a GEAR file to be used. In case it is not specified, gear_telescope.xml will be used
             -e or --eta-file          Specify a ETA DB file to be used for cluster center correction. In case it is not specified, the eta
                                       correction will not be applied

             -h or --help              Print this help message

             --only-generate           Simply genereate the steering file

             -l or --all-local         Execute the job entirely locally.
             -r or --all-grid          Execute the job on the GRID
    """
    print black


# here starts the main
print red, "Hit maker job submitter" , black

# default options
optionGenerateOnly      = 0
optionExecution         = "all-local"

# parse the argmuments
narg = 0
inputFileList         = []
outputFileBaseName    = ""
etaFileName           = ""
gearFileName          = "gear_telescope.xml"
histoFileName         = ""
outputFileName        = ""


goodArg = sys.argv[1:]
skipNext = False

for i,arg in enumerate( goodArg ):
    narg = narg + 1

    if skipNext == True:
        skipNext = False
        continue

    if arg == "-o" or arg == "--output-basename":
        outputFileBaseName = goodArg[ i + 1]
        skipNext = True

    elif arg == "-g" or arg == "--gear-file" :
        gearFileName = goodArg[ i + 1 ]
        skipNext = True

    elif arg == "-e" or arg == "--eta-file" :
        etaFileName = goodArg[ i + 1 ]
        skipNext = True

    elif arg == "--only-generate":
        optionGenerateOnly = 1

    elif arg == "-h" or arg == "--help":
        usage(sys.argv[0])
        sys.exit(0)

    elif arg == "-l" or arg == "--all-local":
        optionExecution = "all-local"

    elif arg == "-r" or arg == "-all-grid":
        optionExecution = "all-grid"

    else :

        # if I'm here it means that this is the input file
        tempList = glob.glob( goodArg[ i ] )
        for item in tempList:
            inputFileName = item
            if not os.path.exists( inputFileName ) :
                print red, "The input file", inputFileName, "doesn't exist!", black
                sys.exit( 1 );
            else :
                inputFileList.append( inputFileName )

if narg == 0:
    usage(sys.argv[0])
    sys.exit(0)

# verify that we have the output file basename
if outputFileBaseName == "" :
    print red, "Output file basename not specified. Please use option -o to specify the output file base name", black
    sys.exit( 1 )

# verify that we have at least one input file
if len( inputFileList ) == 0 :
    print red, "No input file found. No reason to continue!", black
    sys.exit( 0 )

# verify that the gear file exists
if not os.path.exists( gearFileName ) :
    print red, "The specified gear file (%(gear)s) doesn't exist!" % { "gear": gearFileName }, black
    sys.exit( 1 )

if gearFileName == "gear_telescope.xml":
    print red, "W A R N I N G: Using generic gear file! " , black

# verify if the eta file exists
if not os.path.exists(  etaFileName )  :
    if etaFileName == "" :
        print red, "W A R N I N G:  ETA file  not specified, the eta function will not be applied" , black
    else:
        print red, "W A R N I N G:  ETA file (%(eta)s) not found, the eta function will not be applied" % { "eta" : etaFileName } , black


# verify the optionExecution
if optionExecution == "all-grid" :
    print red, "Sorry only local execution available for the time being", black
    sys.exit( 1 )

# prepare the histoFileName and the outputFileName
histoFileName          = outputFileBaseName + "-hit-histo"
outputFileName         = outputFileBaseName + "-hit.slcio"

# steering file name and logfile
steeringFileName       = "hitmaker-" + outputFileBaseName + ".xml"
logFileName            = "hitmaker-" + outputFileBaseName + ".log"

print blue, "Generating the steering file...", black

#open the template steering file for reading
templateSteeringFile = open( "./template/hitmaker-tmp.xml", "r")

# read the whole content of the template and put it into a string
actualSteeringString = templateSteeringFile.read()

# add all the input files
for inputFile in inputFileList:
    actualSteeringString = actualSteeringString.replace( "@InputFile@", "%(input)s @InputFile@" % { "input": inputFile } )

# remove the last @InputFile@ not replaced
actualSteeringString = actualSteeringString.replace( "@InputFile@", "" )

# add the gear file
actualSteeringString = actualSteeringString.replace( "@GearFile@", "%(gear)s" % { "gear": gearFileName })

# add the eta file
actualSteeringString = actualSteeringString.replace( "@EtaFile@", "%(eta)s" % { "eta": etaFileName } )

# add the histogram file name
actualSteeringString = actualSteeringString.replace( "@HistoFile@", "%(histo)s" % {"histo": histoFileName } )

# add the output file name
actualSteeringString = actualSteeringString.replace( "@OutputFile@", "%(output)s" % {"output": outpuFileName } )

# write the new file on disk
actualSteeringFile = open( steeringFileName, "w" )
actualSteeringFile.write( actualSteeringString )
actualSteeringFile.close()
templateSteeringFile.close()

if optionGenerateOnly == 0:

    # run Marlin in both all local and cpu local configuration
    print blue, "Running Marlin...",black

    # to properly get the return code I need to execute Marlin into a separete pipe
    logfile = open( logFileName, "w" )
    marlin = popen2.Popen4( "Marlin %(steer)s" % {"steer" : steeringFileName } )
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
    command = "rm -r %(steer)s %(log)s" % {"steer" : steeringFileName , "log": logFileName }
    returnvalue = os.system( command )
    if returnvalue != 0:
        print red, "Problem removing the steering and log files...", black
