#! /usr/bin/env python
#
#  Python script to display the status of GRID jobs
#
#  This script should work with any file containing JIDs
#  but it works *better* when using the JID produced by
#  pysub
#
#  Author: Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
#  Version: $Id: myjob-status.py,v 1.8 2009-08-01 10:40:23 bulgheroni Exp $

from optparse import OptionParser
import os
import popen2
import time
import datetime

def main() :
    usage = "%prog [options] JID-files"
    version = "$Revision: 1.8 $"
    version = version.replace("$Revision:", "")
    version = version.replace("$", "")
    parser = OptionParser( usage=usage, version=version.strip())


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

    parser.add_option( "-w", "--wait", action="store", type="int", dest="wait", help="Time to wait between one loop and the next when working in continuos mode")

    parser.add_option( "-i", action="store_true", help="Does nothing, only for backward compatibility" )


    parser.set_defaults( type="detailed" )
    parser.set_defaults( repeat=False )
    parser.set_defaults( wait=10 )


    options, args = parser.parse_args()

    if options.wait < 10 :
        options.wait = 10

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

            status = "Status @ %(time)s" % { "time":datetime.datetime.now().strftime("%H:%M:%S") }

            print "=" * 94
            print "| %(id)4s | %(desc)-60s | %(status)-20s |" % { "id": "ID", "desc": "Description", "status": status }
            print "-" * 94

            id = 1
            for filename in goodFile:
                file = open( filename )
                while ( True ):
                    line = file.readline()
                    if line == '' :
                        break
                    if line.startswith( "#" ) :
                        description = line.strip()[0:58]
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
            time.sleep( options.wait )
        
if __name__ == "__main__":
    main()

