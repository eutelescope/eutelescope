#!/bin/env python
#Main program

import sys
def usage():
    print 'Usage: runtests.py [-h] [-V] [--list] [--colors] [-g pdfname] conffile inputfile [referencefile]'
    print ' Performs statistical tests for distributions contained in inputfile'
    print ' tests are configured in conffile'
    print '    conffile : configuration file'
    print '    inputfile : input file containing distributions'
    print '    referencefile (optional) : reference file name that overwrites the one contained in the conf file'
    print '    -h : this help'
    print '    -v : verbose'
    print '    -V : very verbose'
    print '    -g : produce report.pdf with plots'
    print '    --cdash : if producing plots also produce .png files and necessary xml output for CDash integration'
    print '    --list : list all available tests'
    print '    --colors : use ANSI colors for output results table'
    print ''
    print ' Example: runtests.py testconf.qa test.root'
    print '          See content of example directory for an example of a configuration file'

def printList():
    from Tests import _testsMap
    print 'Known Tests:'
    for k,v in _testsMap.iteritems():
        print "Test: %s"%k
        print v.__doc__
        
def run( argv = sys.argv ):
    from Interface import runROOT
    colors = False #By default use ANSI colors for output
    cdash = False

    if len(argv)>=1:
        import getopt
        try:
            opts,args = getopt.getopt( argv ,"hvVg:",["list","colors","cdash"])
        except getopt.GetoptError,err:
            sys.stderr.write( str(err) )
            usage()
            sys.exit(2)

        sys.stderr.write( "Number of arguments: %d\n"%len(args))
            
        doplots = "report.pdf"
        for op,ar in opts:
            if op == '-h':
                usage()
                sys.exit(0)
            elif op == '-v':
                from Interface import _logger
                _logger.setLevel('INFO')
            elif op == '-V':
                from Interface import _logger
                _logger.setLevel('DEBUG')
            elif op == '--list':
                printList()
                sys.exit(0)
            elif op == "--colors":
                colors = True
            elif op == "--cdash":
                cdash = True
            elif op == '-g':
                doplots=ar
        if len(args)>=2:
            try:
                if len(args)>=3:
                    errcode = runROOT( args[0],args[1] , doPlots = doplots , defaultReferenceFile = args[2] , useColorsForOutput = colors, doCDashOutput = cdash)
                else:
                    errcode = runROOT( args[0] , args[1] ,doPlots = doplots , useColorsForOutput = colors, doCDashOutput = cdash)
            except Exception,e:
                sys.stderr.write("Error: %s\n"%e)
                sys.exit(1)
        else:
            usage()
            sys.exit(0)
            
    else:
        usage()
        sys.exit(0)
    if errcode > 0:
        sys.stderr.write( "Number of tests FAILED or NOT PASSED: %d\n"%errcode)
        sys.exit(1)
    sys.exit( 0 )

if __name__ == "__main__":
    run( sys.argv[1:])
