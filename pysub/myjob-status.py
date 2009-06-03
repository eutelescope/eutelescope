#! /usr/bin/env python
#
#  Python script to display the status of GRID jobs
#
#  This script should work with any file containing JIDs
#  but it works *better* when using the JID produced by
#  pysub
#
#  Author: Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
#  Version: $Id: myjob-status.py,v 1.2 2009-06-03 09:22:07 bulgheroni Exp $

from optparse import OptionParser
import os
import popen2
import time

def main() :
    usage = "%prog [options] JID-files"
    version = "$Revision: 1.2 $"
    parser = OptionParser( usage=usage, version=version[10:])


    parser.add_option( "-t", "--type" ,  type="choice",
                       action="store",
                       choices=['detailed', 'summary'],
                       dest="type",
                       help="Select how you want to see the status. Possible choices are \'detailed\' and \'summary\'",
                       metavar="TYPE")

    parser.add_option( "-d", "--detailed", action="store_const", const="detailed", dest="type",
                       help="Print the status of every job in the file with a description")

    parser.add_option( "-s", "--summary", action="store_const", const="summary", dest="type",
                       help="Print a summary of the status of all jobs in the list")

    parser.add_option( "-c", "--repeat-for-ever", action="store_true", dest="repeat", help="Loop indefinetely" )

    parser.add_option( "-f", "--filter", action="store", dest="filter", help="Apply a filter to the input files")


    parser.set_defaults( type="detailed" )
    parser.set_defaults( repeat=False )

    options, args = parser.parse_args()

    # check if the input files exist or not
    goodFile = []
    for file in args:
        if not os.access( file , os.R_OK ):
            print "Error: file %(file)s doens't exists" % { "file": file }
        else:
            goodFile.append( file )

    repeat = True
    first = True
    while repeat :

        if options.type == "detailed" :

            print "=" * 94
            print "| %(id)4s | %(desc)-60s | %(status)-20s |" % { "id": "ID", "desc": "Description", "status": "Status" }
            print "-" * 94

            id = 1
            for filename in goodFile:
                file = open( filename )
                while ( True ):
                    line = file.readline()
                    if line == '' :
                        break
                    if line.startswith( "#" ) :
                        description = line.strip()
                        jid = file.readline()
                    elif line.startswith("https") :
                        description  = "None"
                        jid = line
                    else :
                        print "Error parsing the input file %(file)s "% {"file": filename }
                        exit( 1 )

                    command = "glite-wms-job-status %(jid)s" % {"jid": jid }
                    statusCommand = popen2.Popen4( command, -1 ) 

                    output = ""
                    status = None
                    rtr = statusCommand.wait()
                    output = statusCommand.fromchild.read() 

                    if rtr != 0:
                        status = "Error"
                    else:
                        lineList = output.splitlines()
                        for currLine in lineList:
                            if  currLine.find("Current Status:") != -1 :
                                status = currLine.replace("Current Status:","").strip()
                                break

                    print "| %(id)4d | %(desc)-60s | %(status)-20s |" % { "id": id, "desc": description, "status": status }

                    id = id + 1

                file.close()
            print "=" * 94
            
        elif options.type == "summary" :

            if first:
                print "="*73
                print "| %(scheduled)-15s | %(running)-15s | %(done)-15s | %(aborted)-15s | " \
                      % { "scheduled": "Scheduled", "running": "Running", "done": "Done", "aborted": "Aborted" }
                print "-"*73
                first = False

            scheduled = 0
            done = 0
            running = 0
            aborted = 0
            for filename in goodFile:
                file = open( filename )
                command = "glite-wms-job-status --logfile /dev/null --noint -i %(filename)s > /tmp/status.log" % { "filename": filename }
                rtr = os.system( command )
                
                if rtr != 0:
                    print "Error polling for status in file %(file)s " % {"file": filename }
                else:
                    outputFile = open( "/tmp/status.log" )
                    output = outputFile.read()
                    scheduled = scheduled + output.count( "Scheduled" )
                    done = done + output.count( "Done" )
                    running = running + output.count("Running" )
                    aborted = aborted + output.count("Aborted" )
                    outputFile.close()
                    
                file.close()
            print "| %(scheduled)-15d | %(running)-15d | %(done)-15d | %(aborted)-15d | " % { "scheduled": scheduled, "running": running, "done": done, "aborted": aborted }
            print "-"*73

        repeat = options.repeat
        if repeat:
            time.sleep( 10 )
        
if __name__ == "__main__":
    main()

