#!/bin/sh
# A template of pedestal job
#
# @author Antonio Bulgheroni <mailto:antonio.bulgheroni@gmail.com>
# @version $Id: runjob-clusearch-tmp.sh,v 1.3 2009-07-18 17:24:06 bulgheroni Exp $
#
# errno  0: No error.
# errno  1: Unable to get the input file from the SE.
# errno  2: Unable to get the pedestal file from the SE.
# errno 20: Problem during Marlin execution.
# errno 30: Problem copying and registering the LCIO output to the SE.
# errno 31: Problem copying and registering the Joboutput to the SE.
# errno 32: Problem copying and registering the histogram to the SE
#


#############
# Defining a function to output a command line message
#
# It takes exactly one string
doCommand() {
    echo "> $1 "
    $1
    r=$?
    echo
    return $r
}


#############
# defining a function to get files from GRID
#
# It takes exactly two arguments:
#
# $1: the LFN
# $2: the local file name
#
getFromGRID(){

    doCommand "lcg-cp -v lfn:$1 file:$2"
    r=$?
    return $r


}

#############
# Defining a function to put and register files on GRID
#
# It takes exactly three arguments
#
# $1: the local file name
# $2: the LFN
# $3: the SE
#
putOnGRID() {
    doCommand "lcg-cr -v -d $3 -l lfn:$2 file:$1"
    return $?
}


## Run the job
#
# This is the real executable. It doesn't take any argument and
# it return 0 in case of successful execution or the following error
# codes in case of problems
#
# errno  0: No error.
# errno  1: Unable to get the input file from the SE.
# errno  2: Unable to get the pedestal file from the SE.
# errno 20: Problem during Marlin execution.
# errno 30: Problem copying and registering the LCIO output to the SE.
# errno 31: Problem copying and registering the Joboutput to the SE.
# errno 32: Problem copying and registering the histogram to the SE
#

# To be replaced with the runString in the format %(run)06d
RunString="@RunString@"

# To be replaced with the pedeString in the format %(pede)06d
PedeString="@PedeString@"

# To be replace with the job name used for the identification of
# all files. It should be something like converter
Name="@Name@"

# Define here all the variables modified by the submitter
GRIDCE="@GRIDCE@"
GRIDSE="@GRIDSE@"
GRIDStoreProtocol="@GRIDStoreProtocol@"
GRIDVO="@GRIDVO@"
LFC_HOME="/grid/ilc/eudet-jra1"
GRIDFolderBase="@GRIDFolderBase@"
GRIDFolderLcioRaw="@GRIDFolderLcioRaw@"
GRIDFolderDBPede="@GRIDFolderDBPede@"
GRIDFolderClusearchResults="@GRIDFolderClusearchResults@"
GRIDFolderClusearchJoboutput="@GRIDFolderClusearchJoboutput@"
GRIDFolderClusearchHisto="@GRIDFolderClusearchHisto@"
GRIDLibraryTarball="@GRIDLibraryTarball@"
GRIDILCSoftVersion="@GRIDILCSoftVersion@"

InputLcioRawLFN=$GRIDFolderLcioRaw/run$RunString.slcio
InputPedeLFN=$GRIDFolderDBPede/run$PedeString-ped-db.slcio
OutputLcioLFN=$GRIDFolderClusearchResults/run$RunString-clu-p$PedeString.slcio
OutputJoboutputLFN=$GRIDFolderClusearchJoboutput/$Name-$RunString.tar.gz
OutputHistoLFN=$GRIDFolderClusearchHisto/run$RunString-clu-histo.root

InputLcioRawLocal=$PWD/lcio-raw/run$RunString.slcio
InputPedeLocal=$PWD/db/run$PedeString-ped-db.slcio
OutputLcioLocal=$PWD/results/run$RunString-clu-p$PedeString.slcio
OutputJoboutputLocal=$PWD/log/$Name-$RunString.tar.gz
OutputHistoLocal=$PWD/histo/run$RunString-clu-histo.root
SteeringFile=$Name-$RunString.xml
LogFile=$Name-$RunString.log

echo
echo "########################################################################"
echo "# Starting $Name-$RunString at `date `"
echo "########################################################################"
echo

# prepare the directory structure as local
echo
echo "########################################################################"
echo "# Preparing the directory structure..."
echo "########################################################################"
echo

doCommand "mkdir native"
doCommand "mkdir lcio-raw"
doCommand "mkdir results"
doCommand "mkdir histo"
doCommand "mkdir pics"
doCommand "mkdir db"
doCommand "mkdir log"

# unpack the library
echo
echo "########################################################################"
echo "# Uncompressing the job tarball..."
echo "########################################################################"
echo
doCommand "tar xzvf $GRIDLibraryTarball"

# rename the simjob.slcio because otherwise it gets delete
doCommand "mv simjob.slcio simjob.slcio.keepme"

# from now on doing things to get access to ESA
doCommand "source ./ilc-grid-config.sh"
doCommand "$BASH ./ilc-grid-test-sys.sh || abort \"system tests failed!\" "
doCommand ". $VO_ILC_SW_DIR/initILCSOFT.sh $GRIDILCSoftVersion"
doCommand "$BASH ./ilc-grid-test-sw.sh"
r=$?
if [ $r -ne 0 ] ; then
    echo "******* ERROR: " ; cat ./ilc-grid-test-sw.log
    abort "software tests failed!"
fi
# ILCSOFT ready
echo
echo "########################################################################"
echo "# ILCSOFT ready to use"
echo "########################################################################"
echo

# now it's safe to rename the simjob to the original
doCommand "mv simjob.slcio.keepme simjob.slcio"

# set the list of Marlin plugins and the LD_LIBRARY_PATH
doCommand "export MARLIN_DLL=$PWD/libEutelescope.so"
doCommand "export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH"

echo
echo "########################################################################"
echo "# Getting the input file ${InputLcioRawLocal}"
echo "########################################################################"
echo
# get the input raw file
doCommand "getFromGRID ${InputLcioRawLFN} ${InputLcioRawLocal}"
r=$?
if [ $r -ne 0 ] ; then
    echo "Problem copying ${InputLcioRawLFN}. Exiting with error."
    exit 1
fi

echo
echo "########################################################################"
echo "# Getting the pedestal DB file ${InputPedeLocal}"
echo "########################################################################"
echo
doCommand "getFromGRID ${InputPedeLFN} ${InputPedeLocal}"
r=$?
if [ $r -ne 0 ] ; then
    echo "Problem copying ${InputPedeLFN}. Exiting with error."
    exit 2
fi

# list all the files available
doCommand "ls -al"

# ready to run marlin
echo
echo "########################################################################"
echo "# Starting Marlin `date`"
echo "########################################################################"
echo
c="Marlin $SteeringFile"
echo $c
$c
r=$?

if [ $r -ne 0 ] ; then
    echo "****** Problem running Marlin"
    exit 20
fi

echo
echo "########################################################################"
echo "# Marlin successfully finished `date `"
echo "########################################################################"
echo


# remove the input files
doCommand "rm ${InputLcioRawLocal}"
doCommand "rm ${InputPedeLocal}"

# fix the histogram file
doCommand "hadd -f temp.root empty.root ${OutputHistoLocal}"
doCommand "mv temp.root ${OutputHistoLocal}"

# put back the files to the GRID
echo
echo "########################################################################"
echo "# Copying and registering the output file to SE"
echo "########################################################################"
echo
doCommand "putOnGRID ${OutputLcioLocal} ${OutputLcioLFN} ${GRIDSE}"
r=$?
if [ $r -ne 0 ] ; then
    echo "****** Problem copying the ${OutputLcioLocal} to the GRID"
    exit 30
fi

echo
echo "########################################################################"
echo "# Preparing the joboutput tarball"
echo "########################################################################"
echo
doCommand "tar czvf ${OutputJoboutputLocal} *.log *.xml out err histo/*root"

echo
echo "########################################################################"
echo "# Copying and registering the tarball to SE"
echo "########################################################################"
echo
doCommand "putOnGRID  ${OutputJoboutputLocal} ${OutputJoboutputLFN} ${GRIDSE}"
r=$?
if [ $r -ne 0 ] ; then
    echo "****** Problem copying the ${OutputJoboutputLocal} to the GRID"
    exit 31
fi

echo
echo "########################################################################"
echo "# Copying and registering the histogram file to SE"
echo "########################################################################"
echo
doCommand "putOnGRID ${OutputHistoLocal} ${OutputHistoLFN} ${GRIDSE}"
r=$?
if [ $r -ne 0 ] ; then
    echo "****** Problem copying the ${OutputHistoLocal} to the GRID"
    exit 32
fi

# Job finished
echo
echo "########################################################################"
echo "# Job ($Name-$RunString) finished at `date`"
echo "########################################################################"
echo

