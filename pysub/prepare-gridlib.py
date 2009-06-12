#!/usr/bin/env python

import os
import os.path
from optparse import OptionParser
import ConfigParser
import shutil
import tarfile
import tempfile
import glob
import sys

def main() :

    usage = "usage: %prog [options] additional-files"
    cvsVersion = "$Revision: 1.3 $"
    version = "%prog version" +  cvsVersion[10:len(cvsVersion)-1] 
    parser = OptionParser( usage = usage, version = version )

    parser.add_option( "-o", "--output", type="string", action="store", dest="output",
                       help="The name of the output GRIDLib tarball" )

    parser.add_option("-c", "--config-file", type="string", action="store", dest="config",
                      help="The configuration file")

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
    fileList = configParser.items( "GRIDLibraryContent" )
    for key, file in fileList:
        if not os.access( file, os.R_OK ):
            print "Problem accessing file %(key)s = %(file)s " % { "key":key, "file":file }
            sys.exit( 2 )
        else :
            goodList.append( file )


    for otherFile in args:
        if not os.access( otherFile, os.R_OK ):
            print "Problem accessing file %(file)s " %{ "file":otherFile }
            sys.exit( 3 )
        else :
            goodList.append( otherFile )

    tempDir = tempfile.mkdtemp()
    for goodFile in goodList :
        shutil.copy( goodFile, tempDir )
        
    tarball = tarfile.open( options.output, "w:gz" )
 
    for i in glob.glob( "%(dir)s/*" % { "dir": tempDir } ):
        tarball.add( i , os.path.basename( i ) )

    tarball.close()
    shutil.rmtree( tempDir )
            
    

if __name__ == "__main__" :
    main()
    


