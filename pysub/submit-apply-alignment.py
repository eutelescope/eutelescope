#!/usr/bin/python

# -----------------------------------------------------------------------------------------------------------------#
# script to edit steering template file for apply-alignment                                                        #
#                                                                                                                  #
# usage: ./submit-apply-alignment.py [options]                                                                     #
#                                                                                                                  #
# input lcio file: LocalFolderHitmakerResults/@RunNumber@-@InputCollectionName@.slcio                              #
# alignment file: LocalFolderDBAlign/@AlignmentFile@                                                               #
# output lcio file: LocalFolderHitmakerResults/@RunNumber@-@OutputCollectionName@.slcio                            #
# output histo file: LocalFolderHitmakerHisto/@RunNumber@-@OutputCollectionName@-histo.root                        #
# LocalFolderHitmakerResults, LocalFolderDBAlign, LocalFolderHitmakerResults, LocalFolderHitmakerHisto             #
# are taken from the config file                                                                                   #
#                                                                                                                  #
# V. Libov (libov@mail.desy.de), 30 June 2010                                                                      #
# Last modified: 05 September 2010                                                                                 #
#------------------------------------------------------------------------------------------------------------------#

import commands
from optparse import OptionParser
from os import path
import sys
import os
import ConfigParser
# -----  options parsing ---------
parser = OptionParser()

parser.add_option("-r", "--RunNumber", dest="RunNumber", type="int", default="10307", help="RunNumber to be processed")
parser.add_option("-i", "--InputCollectionName", dest="InputCollectionName", default="hit", help="Input hit collection name")
parser.add_option("-o", "--OutputCollectionName", dest="OutputCollectionName", default="aligned_iter1-hit", help="Ouptut hit collection name")
parser.add_option("-a", "--AlignmentFile", dest="AlignmentFile", default="dummy", help="Alignment file name")
parser.add_option("--InputLCIOFile", dest="InputLCIOFile", default="dummy", help="Input lcio file name")
parser.add_option("-c", "--runCorrelator", dest="runCorrelator", action="store_true", default="false", help="Run Correlator processor or not (time consuming)")
parser.add_option("-f", "--config", dest="configFile", default="config.cfg", help="name of configuration file; assumed to be in config/ subdirectory")

(options, args) = parser.parse_args()

RunNumber = options.RunNumber
InputCollectionName = options.InputCollectionName
OutputCollectionName = options.OutputCollectionName
AlignmentFile = options.AlignmentFile
runCorrelator = options.runCorrelator
# -- config file
configFile = options.configFile
# check if config file is available, otherwise terminate
if not os.path.exists( configFile ):
	print "configuration file ", configFile, " not found. Terminating."
	sys.exit(-1)
# now we can get some stuff from it
configParser = ConfigParser.SafeConfigParser()
configParser.read( configFile )
# -- gear file
GearPath = configParser.get("LOCAL", "LocalFolderGear")
GearFile=configParser.get("General", "GEARFile" )
# directory where input/output is placed
HitmakerResults= configParser.get("LOCAL", "LocalFolderHitmakerResults")
# ok, now full name of input/output files can be set
OutputLCIOFile="%(hit)s/%(run)d-%(out)s.slcio" % { "hit" :HitmakerResults, "run" : RunNumber, "out" : OutputCollectionName }
InputLCIOFile=HitmakerResults+'/'+options.InputLCIOFile
# path to alignment database
AlignmentDBPath=configParser.get("LOCAL", "LocalFolderDBAlign")
# path to histograms
HistoPath = configParser.get("LOCAL", "LocalFolderHitmakerHisto")
# --------end option parsing -------

# ---- create actual xml steering file from the template
OutputFileName = "apply-alignment-%(run)s.xml" % { "run" : RunNumber }
RMFile = "rm -rf %(hit)s/%(run)d-%(out)s.000.slcio " % { "hit" : HitmakerResults, "run" : RunNumber, "out" : OutputCollectionName} 
commands.getoutput(RMFile)
#commands.getoutput("rm -rf "+HitmakerResults+"/"+RunNumber+"-"+OutputCollectionName+".000.slcio")
# commands.getoutput("rm -rf "+OutputLCIOFile)

fileContent = open('template_ibl/apply-alignment-tmp.xml', 'r').read()
sRunNumber = "%(run)d" % { "run" : RunNumber }
fileContent = fileContent.replace("@RunNumber@", sRunNumber)
fileContent = fileContent.replace("@GearPath@", GearPath)
fileContent = fileContent.replace("@GearFile@", GearFile)
fileContent = fileContent.replace("@InputCollectionName@", InputCollectionName)
fileContent = fileContent.replace("@OutputCollectionName@", OutputCollectionName)
fileContent = fileContent.replace("@OutputLCIOFile@", OutputLCIOFile)
fileContent = fileContent.replace("@InputLCIOFile@", InputLCIOFile)
fileContent = fileContent.replace("@AlignmentDBPath@", AlignmentDBPath)
fileContent = fileContent.replace("@AlignmentFile@", AlignmentFile)
fileContent = fileContent.replace("@Results@", HitmakerResults)
fileContent = fileContent.replace("@HistoPath@", HistoPath)
if runCorrelator==True:
	fileContent = fileContent.replace("@LeftComment@", "")
	fileContent = fileContent.replace("@RightComment@", "")
else:
	fileContent = fileContent.replace("@LeftComment@", "!--")
	fileContent = fileContent.replace("@RightComment@", "--")

newfile = open(OutputFileName, 'write')
newfile.write(fileContent)
newfile.close()

# run Marlin and write log file
log = commands.getoutput("Marlin "+OutputFileName)
#logs=os.getenv('logs')

#print log
#logpath = logs + "/apply-alignment." + sRunNumber  
#print logpath
#         "%(run)s-align-mille.bin" % { "run": self._option.output  }

logpath = "logs/apply-alignment.%(run)s-%(out)s.log" % { "run" : sRunNumber, "out" : OutputCollectionName }
print logpath
logFile = open(logpath+"-"+OutputCollectionName+".log", "write")
logFile.write(log)
logFile.close()

# cleaning up - remove actual steering file
commands.getoutput("rm -rf "+OutputFileName)


