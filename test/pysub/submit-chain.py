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
    print commandname, "[ OPTIONS ] -o | --output-file outputfile-name list-of-runs-to-be-merged"
    print """
          List of OPTIONS

             -o or --outputfile   The name of the output file

             -h or --help         Print this help message

             --only-generate      Simply genereate the steering file

             --keep-input         Option not to delete the input files in case they are taken from
                                  a remote localtion. By default remote input files are removed when
                                  chaining is over. In case the files are local they are never deleted

             --keep-output        Option to keep a copy of the output in case it is copied somewhere
                                  remotely. By default remote output file is always deleted, while local
                                  output file is not.

             The list-of-runs-to-be-merged is a file list that can be either local or remote.
             Remote files are identified by a prefix defining the transfer protocol.
             Available protocols are with the corresponding prefix are:

                 Procol    Prefix      Example

                - GRID --> lfn:        lfn:$LFC_HOME/2008/tb-cern-summer/lcio-raw-depfet/run001000.slcio
                - WEB  --> http://     http://www.roma3.infn.it/~toto/run1000.slcio
                       --> https://
                - SSH  --> ssh:        ssh:bulgheroni@blade.roma.infn.it:~/my_data/runs/run001999.slcio

             Remote files are copied locally before processing. 
             Local files are either identified using the file: prefix or even no prefix at all.


             The output file name can be again either local or remote. In case of remote output file a temporary
             file will be generated on the local machine and then copied to the remote location afterwards.
             Available protocols for remote output are GRID (lfn:) and SSH (ssh:).
    """
    print black


# here starts the main
print red, "Chaining job submitter" , black

# default options
optionGenerateOnly = 0
optionKeepInput    = 0
optionKeepOutput   = 0
# defaul GRID paths

# parse the argmuments
narg = 0
fileList           = []
fileListLocal      = []
fileListRemoteGRID = []
fileListRemoteSSH  = []
fileListRemoteWEB  = []
fileMapRemote      = []

outputFileIndex    = -1

goodArg = sys.argv[1:]

remoteCounter = 0;

for i,arg in enumerate( goodArg ):
    narg = narg + 1

    # don't count twice the pedestal run
    if i == outputFileIndex:
        continue

    if arg == "-o" or arg == "--output-file":
        outputFileArg = goodArg[ i + 1]
        outputFileIndex = i + 1

    elif arg == "--only-generate":
        optionGenerateOnly = 1

    elif arg == "-h" or arg == "--help":
        usage(sys.argv[0])
        sys.exit(0)

    elif arg == "--keep-input":
        optionKeepInput = 1

    elif arg == "--keep-output":
        optionKeepOutput = 1

    else :

        # if I'm here this should be one of the file to be merged
        # first guess which protocol we need to use
        if arg[:len("lfn:") ] == "lfn:":
            # this is a grid file
            fileListRemoteGRID.append( arg )

        elif arg[:len("http://") ]   == "http://" or \
             arg[:len("https://") ]  == "https://" :
            # this is from the web
            fileListRemoteWEB.append( arg )

        elif arg[:len("ssh:")  ] == "ssh:" :
            # this is available thourg ssh
            fileListRemoteSSH.append( arg )

        elif arg[:len("file:")  ] == "file:" :
            # this is a local file with the prefix
            # it can contains * or ? or any other regular expression.
            # So explicit them before continuing
            tempList = glob.glob( arg[len("file:")  :] )
            for item in tempList:
                fileListLocal.append( item )

        else :
            # this is a local file w/o prefix
            tempList = glob.glob( arg )
            for item in tempList:
                fileListLocal.append( item )


if narg == 0:
    usage(sys.argv[0])
    sys.exit(0)

if outputFileIndex == -1:
    print red, "Output file name not specified. Please use option -o to specify the output file name", black
    sys.exit(1)

if len( fileListLocal ) + len( fileListRemoteSSH ) + len(fileListRemoteWEB) + len(fileListRemoteGRID) == 0:
    print red, "No input file specified. No need to chain", black
    sys.exit(0)

if  len( fileListRemoteSSH ) + len(fileListRemoteWEB) + len(fileListRemoteGRID) == 0:
    isInputFileAllLocal = 1
else:
    isInputFileAllLocal = 0


# getting the remote file first. We need this only if optionGenereateOnly == 0
if len( fileListRemoteGRID ) != 0 :
    print blue, "Getting all the files from the GRID", black

    for file in fileListRemoteGRID :

        if optionGenerateOnly == 1 :
            fileList.append( "/tmp/remote-%(num)06d.slcio" % { "num": remoteCounter } )
            remoteCounter = remoteCounter + 1
        else: 
            command = "lcg-cp -v %(file)s file:/tmp/remote-%(num)06d.slcio" % { "file": file, "num": remoteCounter }
            returnvalue = os.system( command )
            if returnvalue != 0:
                print red, "Problem getting %(file)s" % { "file": file }, black
            else:
                fileList.append( "/tmp/remote-%(num)06d.slcio" % { "num": remoteCounter } )
                fileMapRemote.append( ( "/tmp/remote-%(num)06d.slcio" % { "num": remoteCounter } , "%(file)s" % { "file":file } ) )
                remoteCounter = remoteCounter + 1

if len( fileListRemoteWEB ) != 0 :
    print blue, "Getting all the files from the WEB", black

    for file in fileListRemoteWEB :

        if optionGenerateOnly == 1 :
            fileList.append( "/tmp/remote-%(num)06d.slcio" % { "num": remoteCounter } )
            remoteCounter = remoteCounter + 1
        else: 
        
            command = "wget -O /tmp/remote-%(num)06d.slcio %(file)s"  % { "file": file, "num": remoteCounter }
            returnvalue = os.system( command )
            if returnvalue != 0:
                print red, "Problem getting %(file)s" % { "file": file }, black
            else:
                fileList.append( "/tmp/remote-%(num)06d.slcio" % { "num": remoteCounter } )
                fileMapRemote.append( ( "/tmp/remote-%(num)06d.slcio" % { "num": remoteCounter } , "%(file)s" % { "file":file } ) )
                remoteCounter = remoteCounter + 1

if len( fileListRemoteSSH ) != 0 :
    print blue, "Getting all the files from SSH (password maybe needed!)", black

    for file in fileListRemoteSSH :

        if optionGenerateOnly == 1 :
            fileList.append( "/tmp/remote-%(num)06d.slcio" % { "num": remoteCounter } )
            remoteCounter = remoteCounter + 1
        else: 
            
            command = "scp %(file)s /tmp/remote-%(num)06d.slcio" % { "file": file[len("ssh:") :], "num": remoteCounter }
            returnvalue = os.system( command )
            if returnvalue != 0:
                print red,  "Problem getting %(file)s" % { "file": file }, black
            else:
                fileList.append( "/tmp/remote-%(num)06d.slcio" % { "num": remoteCounter } )
                fileMapRemote.append( ( "/tmp/remote-%(num)06d.slcio" % { "num": remoteCounter } , "%(file)s" % { "file":file } ) )
                remoteCounter = remoteCounter + 1

if len( fileListLocal ) != 0 :
    for file in fileListLocal :
        fileList.append( file )

# handling the output file
# guess if it is local or remote
if outputFileArg[:len("lfn:") ] == "lfn:" or \
   outputFileArg[:len("ssh:") ] == "ssh:" :
    isOutputFileLocal = 0
    outputFileName = "/tmp/output.slcio"
    outputFileNameShort = "output.slcio"


elif outputFileArg[:len("file:") ] == "file:" :
    isOutputFileLocal = 1
    outputFileName = outputFileArg[len("file:") :]
    index = outputFileName.rfind( "/")
    if index != -1:
        outputFileNameShort = outputFileName[index + 1:]
    else:
        outputFileNameShort = outputFileName

else :
    isOutputFileLocal = 1
    outputFileName = outputFileArg
    index = outputFileName.rfind( "/")
    if index != -1:
        outputFileNameShort = outputFileName[index + 1:]
    else:
        outputFileNameShort = outputFileName

print blue, "Generating the steering file...", black

#open the template steering file for reading
templateSteeringFile = open( "./template/chain-tmp.xml", "r")

# read the whole content of the template and put it into a string
actualSteeringString = templateSteeringFile.read()

for file in fileList :

    actualSteeringString = actualSteeringString.replace( "@RunList@", "%(run)s @RunList@" % { "run": file } )

actualSteeringString = actualSteeringString.replace( "@RunList@", "" )
actualSteeringString = actualSteeringString.replace( "@OutputFile@", "%(output)s" % { "output": outputFileName } )

actualSteeringFile = open( "chain-%(output)s.xml" % {"output" : outputFileNameShort }, "w" )
actualSteeringFile.write( actualSteeringString )
actualSteeringFile.close()
templateSteeringFile.close()

if optionGenerateOnly == 0:

    # run Marlin in both all local and cpu local configuration
    print blue, "Running Marlin...",black

    # to properly get the return code I need to execute Marlin into a separete pipe
    logfile = open( "chain-%(run)s.log " % { "run": outputFileNameShort }, "w" )
    marlin = popen2.Popen4( "Marlin chain-%(run)s.xml" % { "run": outputFileNameShort } )
    while marlin.poll() == -1:
        line = marlin.fromchild.readline()
        print line.strip()
        logfile.write( line )

    logfile.close()
    returnvalue = marlin.poll()

    if returnvalue != 0:
        print red, "Problem executing Marlin! (errno %(returnvalue)s)" % {"returnvalue":returnvalue }, black
        sys.exit( 1 )

    if isOutputFileLocal == 0 :
        # we have to copy the output file somewhere
        # guess the protocol
        if outputFileArg[:len("lfn:") ] == "lfn:" :
            print blue, "Copying and registering the output file to the GRID...", black
            command =  "lcg-cr -v -l %(lfn)s file:%(file)s" % { "lfn" : outputFileArg, "file": outputFileName } 
            returnvalue = os.system( command )
            
            if returnvalue != 0:
                print red, "Problem copying the output file to the GRID! (errno %(returnvalue)s)" % {"returnvalue":returnvalue }, black
            else :
                if optionKeepOutput == 0:
                    # remove the output file
                    print blue, "Removing the local copy of the output file", black
                    command = "rm -v %(file)s" % { "file": outputFileName }
                    returnvalue = os.system( command )
                    if returnvalue != 0:
                        print red, "Problem removing the local copy of the output file", black
                else :
                    print blue, "The local copy of the output file is", outputFileName, black
                    
        elif outputFileArg[:len("ssh:")] == "ssh:" :
            print blue, "Copying the output file via SSH...", black
            command = "scp %(file)s %(ssh)s" % {"file":outputFileName , "ssh": outputFileArg[ len("ssh:") : ] }
            returnvalue = os.system( command )

            if returnvalue != 0:
                    print red, "Problem copying the output file via SSH! (errno %(returnvalue)s)" % {"returnvalue":returnvalue }, black
            else :
                if optionKeepOutput == 0:
                    # remove the output file
                    print blue, "Removing the local copy of the output file", black
                    command = "rm -v %(file)s" % { "file": outputFileName }
                    returnvalue = os.system( command )
                    if returnvalue != 0:
                        print red, "Problem removing the local copy of the output file", black
                else :
                    print blue, "The local copy of the output file is", outputFileName, black

    # removing the non local input files
    if isInputFileAllLocal == 0:
        if optionKeepInput == 0:
            print blue, "Removing the local copies of the remote input files...", black
            command = "rm -r /tmp/remote-*.slcio"
            returnvalue = os.system( command )
            if returnvalue != 0:
                print red, "Problem removing the local copies of the remote input files!  (errno %(returnvalue)s)" % {"returnvalue":returnvalue }, black
        else :
            print blue, "Local copy of the input files:"
            for local,remote in fileMapRemote :
                print "%(local)s ==> %(remote)s" % {"local":local, "remote":remote}
            
    print blue, "Removing the steering and log files...", black
    command = "rm -r chain-%(output)s.*" % {"output" : outputFileNameShort }
    returnvalue = os.system( command )
    if returnvalue != 0:
        print red, "Problem removing the steering and log files...", black

        
