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
    parser.add_option( "--plus",  action="store_true", dest="plus",  help="add a plus sign" )
    parser.add_option( "--minus", action="store_true", dest="minus", help="add a minus sign" )

    # create the new submitter, this will automatically configure it
    submitTest = SubmitTest( parser )

    # now execute it! 
    submitTest.execute()

    # finish
    submitTest.end()

if __name__ == "__main__" :
    main()
