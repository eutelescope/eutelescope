#!/usr/bin/python
#Input to this function is as follows: resolutions (All), excluded planes, planes(including duts, therefore all),the factor
import fileinput
import sys 
import logging
import os
from tempfile import mkstemp
from shutil import move
from os import remove, close
#Functions
def main():
	#Before we do anything we need to find the line number to insert the results into.
	start=0
	end=0
	start=findLineNumberFromString("Parameter",sys.argv[1])
	end=findLineNumberFromString("method inversion", sys.argv[1])
	print "Start %d"  %start
	print "End %d"  %end
	#Fill resultsLines before we loop through steering file
	with open(sys.argv[2]) as result:
		resultLines = result.readlines()
	#Create temp file. This will become the steer.txt file after.
	fh ,os.curdir = mkstemp()
	temp = open(os.curdir,'w')
	#Open old file to begin looping over
	steer = open(sys.argv[1])
	counter=1
	resultLineCounter=0
	for lineSteer in  steer:
		if(counter < start):
			temp.write(lineSteer)
		elif( counter >= start and counter < (end-1)):
			temp.write(resultLines[resultLineCounter])
			resultLineCounter += 1
		elif( counter >= (end-1)):
			temp.write(lineSteer)
		else:
			print "Some thing is wrong in resToSteer python script"
		counter += 1
	temp.close()
	close(fh)
	steer.close()

	remove(sys.argv[1])
	move(os.curdir, sys.argv[1])

def findLinesToReadBetween(start,end,filename):
	lookup1 = 'Parameter'
	start=findLineNumberFromString(lookup1,filename)
	lookup2 = 'method inversion'
	end=findLineNumberFromString(lookup2,filename)
	print end


def findLineNumberFromString(lookup,filename):
	with open(filename) as myFile:
		for num, line in enumerate(myFile, 1):
			if lookup in line:
				print "found line: %s at %d"  %(line, num)
				return num

#END OF FUNCTIONS////////////////////////////////////////////////////////////////////
print "Executing Results into steering file python script"
main()

