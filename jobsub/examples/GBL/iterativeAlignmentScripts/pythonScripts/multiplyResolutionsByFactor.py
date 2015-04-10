#!/usr/bin/python
#Input to this function is as follows: resolutions (All), excluded planes, planes(including duts, therefore all),the factor
import sys
import logging
#Functions
def findResolutionExcludeAndFactorTerms(input,resolutions,exclude,allDuts):
	slashCount=0
	for i in range(1,int( len(input))):
		if input[i] == "/":
			slashCount+=1
		else:
			logging.debug('Entering else statement to fill vectors')
			if slashCount == 0:
				resolutions.append(float(input[i]))
			elif slashCount == 1:
				exclude.append(int(input[i]))
			elif slashCount == 2:
				allDuts.append(int(input[i]))
			elif slashCount == 3:
				return float(input[i])
			else:
				sys.exit("We have went to far. Problem in findResolutionExcludeAndFactorTerms()")

def newResolution(resolutions,exclude,allDuts,factor):
	if len(resolutions) != len(allDuts):
		print("Here is the resolution vector size and number of DUTs: %s and %s  " % (len(resolutions), len(allDuts)))
		sys.exit("The resolution is not the same as the number of duts? In newResolution()")
	newResolution=[]
	if exclude:
		for dut, res  in zip(allDuts, resolutions):
			for exc in exclude:
				if dut == exc:
					newResolution.append(res)
					break
				elif dut != exc and exc == exclude[-1]:
					newResolution.append(res*factor)
	elif not exclude: 
		for dut, res  in zip(allDuts, resolutions):
			newResolution.append(res*factor)

	return newResolution

def stringResolution(newResolution):
	stringResolution=str(newResolution)
	edit1 = stringResolution.replace("[","")
	edit2 = edit1.replace("]","")
	edit3 = edit2.replace(",","")

	return edit3
#End of functions.
input=sys.argv
resolutions=[]
exclude=[]
allDuts=[]
#Find the initial position of the different inputs from the vector.
#Must return integer since does not set with factor=input[i] in function
factor=findResolutionExcludeAndFactorTerms(input,resolutions,exclude,allDuts)
newResolutions=newResolution(resolutions,exclude,allDuts,factor)
stringResolutions=stringResolution(newResolutions)
print stringResolutions


