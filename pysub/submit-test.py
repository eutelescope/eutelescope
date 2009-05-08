#!/usr/bin/python
import sys
import os.path
import commands
from optparse import OptionParser
from pysub import *

def main() :

    usage = "usage: %prog [options] run-list"

    parser = OptionParser( usage = usage )
    parser.add_option( "--config-file", action="store", dest="config_file",
                       help="Specify the configuration file to be used" )

    # create the new submitter, this will automatically configure it
    submitTest = SubmitTest( parser )

    # now execute it! 

if __name__ == "__main__" :
    main()
