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
	#Fill resultsLines before we loop through steering file
	with open(sys.argv[2]) as result:
		resultLines = result.readlines()
	#Create temp file. This will become the steer.txt file after.
	fh ,os.curdir = mkstemp()
	temp = open(os.curdir,'w')
	#Open old file to begin looping over
	steer = open(sys.argv[1])
	counter=0
	resultLineCounter=0
	for lineSteer in  steer:
		if(counter < 3):
			temp.write(lineSteer)
		elif( counter >= 3 and counter < 39):
			temp.write(resultLines[resultLineCounter])
			resultLineCounter += 1
		elif( counter >= 39):
			temp.write(lineSteer)
		else:
			print "Some thing is wrong in resToSteer python script"
		counter += 1
	temp.close()
	close(fh)
	steer.close()

	remove(sys.argv[1])
	move(os.curdir, sys.argv[1])

#END OF FUNCTIONS////////////////////////////////////////////////////////////////////
print "Executing Results into steering file python script"
main()

