#! /usr/bin/env python
import os
from pysub import SubmitNativeCopy
from pysub import SubmitBase
from pysub import StopExecutionError
from optparse import OptionParser
from optparse import OptionGroup

## The native copy submitter script.
#
# This script is simply defining all the input options for the
# SubmiNativeCopy and create an instance of this object
#
# @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
# @version $Id: submit-native-copy.py,v 1.4 2009-05-13 12:41:18 bulgheroni Exp $

def main() :

    usage   = """
%prog is part of pysub the Job Sumbitter environment of EUTelescope.

usage: %prog [options] [configuration-options] run-list
"""
    cvsVersion = "$Revision: 1.4 $"
    submitNativeCopyCVSVersion = SubmitNativeCopy.cvsVersion
    submitBaseCVSVersion = SubmitBase.cvsVersion
    version = "%prog version" + cvsVersion[10:len(cvsVersion)-1] + \
        "\nSubmitNativeCopy version " + submitNativeCopyCVSVersion[10:len(submitNativeCopyCVSVersion)-1] + \
        "\nSubmitBase class version " + submitBaseCVSVersion[10:len(submitBaseCVSVersion)-1] + \
        "\ncompiled on a " + os.name + " system"

    parser = OptionParser( version=version, usage = usage )

    # adding the verbose option
    parser.add_option( "-v", "--verbose",
                       action="store_true",
                       dest="verbose",
                       help="sake the output of GRID commands verbose" )


    # adding a group will all the configuration options
    configurationGroup = OptionGroup( parser, "Configuration option",
                                      "Use these options to specify which configuration file you want to use")

    configurationGroup.add_option( "--config-file",
                                   action="store",
                                   dest="config_file",
                                   help="specify the configuration file to be used")

    parser.add_option_group( configurationGroup )

    submitNativeCopy = SubmitNativeCopy( parser )

    try :
        submitNativeCopy.execute()

    except StopExecutionError, error:
        submitNativeCopy._logger.critical( "Cannot continue with execution" )

    # good bye!
    submitNativeCopy.end()


if __name__ == "__main__":
    main()
