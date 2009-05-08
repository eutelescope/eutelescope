#! /usr/bin/env python
from pysub import *
from optparse import OptionParser
from optparse import OptionGroup


def main() :

    usage   = "usage: %prog [execution-options] [io-options] [configuration-options] run-list"
    version = "$Revision: 1.1 $"
    parser = OptionParser( version=version, usage = usage )


    # adding in a group all the execution options.
    executionGroup = OptionGroup( parser, "Execution options",
                                  "Use these options to select where and how the jobs have to executed")

    executionHelp = """
    Select where to execute the job. all-local means: input and output files are stored locally and the job is executed on the local CPU.
    all-grid means: input and output files are taken from the storage element and the job will be submitted to the GRID.
    cpu-local means: input and output files are taken from the GRID SE, but the job will be executed on the local CPU
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
                        action="store_false",
                        dest="force_keep_input",
                        help="Force to delete the input files after finishing, independently of the execution mode")

    ioGroup.add_option( "--keep-output",
                        action="store_true",
                        dest="force_keep_output",
                        help="Force not to delete the output files after finishing, independently of the execution mode")

    ioGroup.add_option( "--remove-output",
                        action="store_false",
                        dest="force_keep_output",
                        help="Force to delete the output files after finishing, independently of the execution mode")


    parser.set_defaults(force_keep_input=True)
    parser.set_defaults(force_keep_output=True)

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

    parser.add_option_group( configurationGroup )

    parser.add_option_group( ioGroup )

    # end of options


    # create the new submitter! 
    submitConversion = SubmitConverter( parser )

    # execute it!
    submitConversion.execute()

if __name__ == "__main__" :
    main()

