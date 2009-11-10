#!/usr/bin/env python

import os
import os.path
from optparse import OptionParser
import ConfigParser
import shutil
import tarfile
import tempfile
import glob
import commands
import sys


## Ask yes or no
#
def askYesNo(  prompt, complaint= "Yes or no please", retries = 4 ):
    while True:
        ok = raw_input( prompt )
        ok = ok.lower()
        if ok in ['y','yes','yea','yeap']  : return True
        if ok in ['n','no',"nope" ] : return False
        retries = retries - 1
        if retries < 0 : raise IOError, "Wrong answer"
        print complaint

def main() :

    usage = "usage: %prog [options] additional-files"
    cvsVersion = "$Revision: 1.9 $"
    version = "%prog version" +  cvsVersion[10:len(cvsVersion)-1]
    parser = OptionParser( usage = usage, version = version )

    parser.add_option( "-o", "--output", type="string", action="store", dest="output",
                       help="The name of the output GRIDLib tarball. If this is starting with lfn:"
                       " then it will be copied to the default storage element and registered with this logical filename" )

    parser.add_option("-c", "--config-file", type="string", action="store", dest="config",
                      help="The configuration file")

    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="Make the script verbose")

    parser.set_defaults( output="gridlib.tar.gz" )

    options, args = parser.parse_args()

    if options.config == None:
        # check if the user set the configuration file via a envirom
        # variable
        try:
            configFile = os.environ['SUBMIT_CONFIG']
        except KeyError:
            configFile = os.path.join( os.getcwd(), "config/config.cfg" )
    else :
        configFile = os.path.abspath( options.config )

    if not os.access( configFile, os.R_OK ):
        print "Problem accessing the configuration file. Please use either the -c option or the SUBMIT_CONFIG variable"
        sys.exit( 1 )

    configParser = ConfigParser.SafeConfigParser()
    configParser.read( configFile )


    goodList = []
    linkList = []
    fileList = configParser.items( "GRIDLibraryContent" )
    for key, file in fileList:
        if len(glob.glob( os.path.expandvars( os.path.expanduser( file )))) == 0:
            print "Problem accessing file %(key)s = %(file)s. Skipping " % { "key":key, "file":file }
        # in principle each file in the list maybe a glob
        for file1 in glob.glob( os.path.expandvars( os.path.expanduser( file ) )):
            if os.path.islink( file1 ):
                linkList.append( file1 )
            else :
                goodList.append( file1 )


    for otherFile in args:
        if not os.access( otherFile, os.R_OK ):
            print "Problem accessing file %(file)s " %{ "file":otherFile }
            sys.exit( 3 )
        else :
            goodList.append( otherFile )

    tempDir = tempfile.mkdtemp()

    # first copy the real file
    for goodFile in goodList :
        shutil.copy( goodFile, tempDir )

    # now copy the links
    for link in linkList:
        dest   = os.path.join( tempDir , os.path.basename ( link )  )
        linkto = os.readlink( link )
        os.symlink(linkto, dest)

    # check if the ouptut is a local file or a lfn
    isTarballLocal = True
    if options.output.startswith( "lfn:" ):
        isTarballLocal = False
        filename = os.path.basename( options.output )
    else :
        isTarballLocal = True
        filename = options.output

    tarball = tarfile.open( filename, "w:gz" )
    tarball.dereference=True

    for i in glob.glob( "%(dir)s/*" % { "dir": tempDir } ):
        if options.verbose :
            print "Adding %(i)s" % { "i": os.path.basename( i ) }
        tarball.add( i , os.path.basename( i ) )

    tarball.close()
    shutil.rmtree( tempDir )

    if not isTarballLocal:
        # check if the file already exists on the GRID
        command = "lfc-ls %(file)s" % { "file": options.output.lstrip("lfn:") }
        status, output = commands.getstatusoutput( command )
        if options.verbose:
            for line in output.splitlines():
                print line.strip()

        if status == 0:
            prompt = "File %(file)s is already on the GRID, do you want to remove it? [y/n]" % { "file": options.output.lstrip("lfn:") }
            if askYesNo( prompt ):
                command = "lcg-del -a %(lfn)s" %{ "lfn": options.output }
                status, output = commands.getstatusoutput( command )
                if options.verbose:
                    for line in output.splitlines():
                        print line.strip()

            else:
                sys.exit( 1 )

        print "Copying to the GRID..."
        command = "lcg-cr -v -l %(lfn)s file:$PWD/%(local)s" % { "lfn": options.output, "local": filename }
        status, output = commands.getstatusoutput( command )

        if options.verbose :
            for line in output.splitlines():
                print line.strip()


if __name__ == "__main__" :
    main()




