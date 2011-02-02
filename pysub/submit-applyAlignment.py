#! /usr/bin/env python
import os
from pysub import SubmitApplyAlignment
from pysub import SubmitBase
from pysub import StopExecutionError
from optparse import OptionParser
from optparse import OptionGroup


## The fitting submitter script.
#
# This script is simply defining all the input options for
# the SubmitApplyAlignment and create an instance of this object.
#
# @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
# @version $Id: submit-applyAlignment.py,v 1.2 2009-06-06 12:38:40 bulgheroni Exp $

def main() :

    usage = """
%prog is part of pysub the Job Sumbitter environment of EUTelescope.

usage: %prog [execution-options] [io-options] [configuration-options] [applyAlignment-options] -o basename-output -a alignment-file input-files
"""
    cvsVersion = "$Revision: 1.2 $"
    submitApplyAlignmentCVSVersion = SubmitApplyAlignment.cvsVersion
    submitBaseCVSVersion = SubmitBase.cvsVersion
    version = "%prog version" + cvsVersion[10:len(cvsVersion)-1] + \
        "\nSubmitApplyAlignment class version " + \
        submitApplyAlignmentCVSVersion[10:len(submitApplyAlignmentCVSVersion)-1] + \
        "\nSubmitBase class version " + submitBaseCVSVersion[10:len(submitBaseCVSVersion)-1] + \
        "\ncompiled on a " + os.name + " system"


    parser = OptionParser( version=version, usage = usage )

    parser.add_option( "-o", "--output", type="string", action="store", dest="output", metavar="OUTPUT",
                       help = "This is the base name of the output file. This string will be used to generate"
                       " the name of the output files according the following naming convention:"
                       " Track LCIO file ==> OUTPUT-alignedHit.slcio, histogram file ==> OUTPUT-alignedHit-histo.root "
                       )

    parser.add_option( "-a", "--alignment-file", type="string", action="store", dest="alignment_file", 
                       help ="This is the name of the alignment file to be used." )


    # adding a group with alignment related options
    applyAlignmentGroup = OptionGroup( parser, "ApplyAlignment options", "Use these options to modify the behaviour of the applyAlignment processor")

    applyAlignmentGroup.add_option( "--split-job", action="store", type="int", dest="split_job", help="Use this option to define how many alignment job you want to submit."
                            " This option can be used only when working in all-grid mode and SPLIT_JOB jobs will be submitted each of them analyzing a SPLIT_SIZE number "
                            "of events.")

    applyAlignmentGroup.add_option( "--split-size", action="store", type="int", dest="split_size",
                            help="Use this option to define how many records should be analyzed by each fitting job" )

    applyAlignmentGroup.add_option( "--skip", action="store", type="int", dest="skip", help="Use this parameter to skip some records from the beginning of the input files" )


    parser.set_defaults( split_job=1 )
    parser.set_defaults( split_size=1000 )
    parser.set_defaults( skip=0 )

    parser.add_option_group( applyAlignmentGroup )

    # adding in a group with all the execution options.
    executionGroup = OptionGroup( parser, "Execution options",
                                  "Use these options to select where and how the jobs have to executed")

    executionHelp = """
Select where to execute the job.
all-local means: input and output files are stored locally
and the job is executed on the local CPU.
all-grid means: input and output files are taken
from the storage element and the job will be submitted to the GRID.
cpu-local means: input and output files are taken
from the GRID SE, but the job will be executed on the local CPU
    """
    executionGroup.add_option( "-x", "--execution",
                               type="choice",
                               action="store",
                               choices=['all-local', 'all-grid', 'cpu-local','only-generate'],
                               dest="execution",
                               help=executionHelp,
                               metavar="EXECUTION")

    parser.set_defaults(execution="all-local")

    executionGroup.add_option( "-l", "--all-local",
                               action="store_const",
                               const="all-local",
                               dest="execution",
                               help="Same as -x all-local")

    executionGroup.add_option( "-r", "--all-grid",
                               action="store_const",
                               const="all-grid",
                               dest="execution",
                               help="Same as -x all-grid" )

    executionGroup.add_option( "-c", "--cpu-local",
                               action="store_const",
                               const="cpu-local",
                               dest="execution",
                               help="Same as -x cpu-local")

    executionGroup.add_option( "-s", "--only-generate",
                               action="store_const",
                               const="only-generate",
                               dest="execution",
                               help="Same as -x only-generate")

    parser.add_option_group( executionGroup )

    # adding a group with all the I/O options
    ioGroup = OptionGroup( parser, "Input and output files related options",
                           "Use these options to specify whether or not to keep the input and output files"
                           "These values are overwriting the default configurations depending on the execution mode:"
                           "  all-local = input and output files are kept "
                           "  cpu-local = input and output files are removed. ")

    ioGroup.add_option( "--keep-input",
                        action="store_true",
                        dest="force_keep_input",
                        help="Force not to delete the input files after finishing, independently of the execution mode")

    ioGroup.add_option( "--remove-input",
                        action="store_true",
                        dest="force_remove_input",
                        help="Force to delete the input files after finishing, independently of the execution mode")

    ioGroup.add_option( "--keep-output",
                        action="store_true",
                        dest="force_keep_output",
                        help="Force not to delete the output files after finishing, independently of the execution mode")

    ioGroup.add_option( "--remove-output",
                        action="store_true",
                        dest="force_remove_output",
                        help="Force to delete the output files after finishing, independently of the execution mode")

    ioGroup.add_option( "--verify-output",
                        action="store_true",
                        dest="verify_output",
                        help="Verify that the output files copied to the GRID are equal to the local ones" +
                        "Option available on when working in cpu-local mode and requires more time" )

    ioGroup.add_option( "-v", "--verbose",
                        action="store_true",
                        dest="verbose",
                        help="sake the output of GRID commands verbose" )

    parser.set_defaults(force_keep_input=False)
    parser.set_defaults(force_keep_output=False)
    parser.set_defaults(force_remove_input=False)
    parser.set_defaults(force_remove_output=False)
    parser.set_defaults(verify_output=False)
    parser.set_defaults(verbose=False)

    parser.add_option_group( ioGroup )

    # adding a group will all the configuration options
    configurationGroup = OptionGroup( parser, "Configuration option",
                                      "Use these options to specify which configuration file you want to use,"
                                      "and other additional external files" )

    configurationGroup.add_option( "-g", "--gear-file",
                                   action="store",
                                   dest="gear_file",
                                   help="Specify the GEAR file to be used")

    configurationGroup.add_option( "--config-file",
                                   action="store",
                                   dest="config_file",
                                   help="Specify the configuration file to be used")

    configurationGroup.add_option( "--event-range",
                                   action="store",
                                   dest="event_range",
                                   help="Set the event range to process")


    parser.add_option_group( configurationGroup )

    # end of options

    # create the new submitter!
    submitApplyAlignment = SubmitApplyAlignment( parser )

    # execute it!
    try :
        submitApplyAlignment.execute()

    except StopExecutionError, error:
        submitApplyAlignment._logger.critical( "Cannot continue the execution" )

    # good bye!
    submitApplyAlignment.end()

if __name__ == "__main__" :
    main()
