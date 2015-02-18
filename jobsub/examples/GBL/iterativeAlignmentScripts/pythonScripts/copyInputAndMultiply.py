#!/usr/bin/python

import sys
resolution=str(sys.argv[1])
numberOfDUTs=str(sys.argv[2])
output=""
for i in range(0,int(numberOfDUTs)):
	output+=str(resolution)
	output+=" "
print output
