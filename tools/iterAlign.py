#!/usr/bin/env python2
import os
import argparse
import sys
import os
import subprocess
"""
iterAlign: Will run the iterative alignment procedure from any GBL directory. The directory must have required structure. See GBL examples.
"""
def checkDollar(var):
    """
    Replace input variables of bash symbols with actual values. 
    If the symbol can not be found then we exit.
    Return the string with the actual values. 
    """
    newParts =[]
    print "Here is the variable ", var
    output = ""
    parts = var.split(" ")
    for part in parts:
      #  print "part: " , part
        if part.find("$") >= 0 :
       #     print "Found dollar"
        #    print part
            ind = part.replace('$',"")
            if os.environ.get(ind) == None:
                print "Can not find the variable to replace? Must end. Sorry!"
                sys.exit(-1)
            newParts.append(os.environ[ind])
         #   print newParts
        else:
            newParts.append(part)

  #  print "Here is newParts: " , newParts
    for newpart in newParts:
 #       print "Inside output", newpart
        output += newpart + " "

 #   print "Here is the dollar removed variable: " , output        
    return output

def findIterScripts():
    eutel = os.environ['EUTELESCOPE']
    output = eutel + "/jobsub/examples/GBL/iterativeAlignmentScripts"
    return output


def main():
    """
    This will take in the required arguments. Find the location of files needed and export all variables to be used in bash scripts
    The it will run the first bash script to start the process.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', help="The number of iterations", default=1)
    parser.add_argument('-i', help="This is the string to identify the output histograms and geometry files",default="DEFAULT" )
    parser.add_argument('-r', help="The run number to align",default="DEFAULT")
    parser.add_argument('-o', help="option: End after single run with only XYShifts. 1/0",default="0")
    args = parser.parse_args()
    if args.r == "DEFAULT":
        print "No run number given"
        print "Use the flag -h to see help"
        sys.exit(-1)
    if args.i == "DEFAULT":
        print "No identifier given"
        print "Use the flag -h to see help"
        sys.exit(-1)
    #Set the variables to export to bash scripts 
    os.environ["RUN"] = args.r
    os.environ["numberOfIterations"] = args.n
    os.environ["outputIdentifier"] = args.i
    os.environ["singleLoop"] = args.o
    #Get the scripts which will actually do the work and location of the example and add to the system path.
    itScriptLoc=findIterScripts()
    itScriptLocPython=itScriptLoc + "/pythonScripts"
    sys.path.append(itScriptLocPython)
    from configRead import readConfig
    from findVariable import findVar

    base=findVar("BasePath")
    string = str(base[1]).replace(" ","")
    string = string.replace('\n',"")
    os.environ["exampleLocation"] = string

    vars = readConfig(string)

    for var in vars:
        print "var: ",var
        value = str(checkDollar(var[1]))
        os.environ[var[0]] = value

    bashPath = itScriptLoc + "/initialiseAndRun.sh"
    os.environ["scriptsLocation"] = itScriptLoc

    print bashPath
    subprocess.call(bashPath, shell=True)

#RUN THE PROGRAM HERE!
main()
