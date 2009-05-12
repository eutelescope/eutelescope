#!/bin/sh

# A template of converter job
#
# @author Antonio Bulgheroni <mailto:antonio.bulgheroni@gmail.com>
# @version $Id: runjob-converter-tmp.sh,v 1.2 2009-05-12 16:48:20 bulgheroni Exp $
#

#############
# Defining a function to output a command line message
# 
# It takes exactly one string
doCommand() {
    echo "> $1 "
    $1
    return $?
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
    echo "Copying $1 in $2 "

    doCommand "lcg-cp -v lfn:$1 file:$2"
    return $?

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
    echo "Copying and registering $1 in $2 on $3"

    doCommand "lcg-cr -v -d $3 -l lfn:$2 file:$1"
    return $?
}


# This is the begin!
echo "Starting universal-$RunString at `date `"

# To be replaced with the runString in the format %(run)06d
RunString="@RunString@"

# Define here all the variables modified by the submitter
GRIDCE="@GRIDCE@"
GRIDSE="@GRIDSE@"
GRIDStoreProtocol="@GRIDStoreProtocol@"
GRIDVO="@GRIDVO@"
LFC_HOME="/grid/ilc/eudet-jra1"
GRIDFolderBase="@GRIDFolderBase@"
GRIDFolderNative="@GRIDFolderNative@"
GRIDFolderLcioRaw="@GRIDFolderLcioRaw@"
GRIDFolderConvertJoboutput="@GRIDFolderConvertJoboutput@"
GRIDLibraryTarball="@GRIDLibraryTarball@"
GRIDILCSoftVersion="@GRIDILCSoftVersion@"

InputRawLFN=$GRIDFolderNative/run$RunString.raw
OutputLcioLFN=$GRIDFolderLcioRaw/run$RunString.slcio
OutputJoboutputLFN=$GRIDFolderConvertJoboutput/universal-$RunString.tar.gz


InputRawLocal=$PWD/native/run$RunString.raw
OutputLcioLocal=$PWD/lcio-raw/run$RunString.slcio
OutputJoboutputLocal=$PWD/log/universal-$RunString.tar.gz
SteeringFile=universal-$RunString.xml
LogFile=universal-$RunString.log


# prepare the directory structure as local
echo "Preparing the directory structure..."
doCommand "mkdir native"
doCommand "mkdir lcio-raw"
doCommand "mkdir results"
doCommand "mkdir histo"
doCommand "mkdir pics"
doCommand "mkdir db"
doCommand "mkdir log"

# unpack the library
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
    echo "ERROR: " ; cat ./ilc-grid-test-sw.log
    abort "software tests failed!"
fi
# ILCSOFT ready
echo " "
echo " ILCSOFT ready to use "
echo " ------------------------------------------"
echo " "

# now it's safe to rename the simjob to the original
doCommand "mv simjob.slcio.keepme simjob.slcio"

# set the list of Marlin plugins and the LD_LIBRARY_PATH
doCommand "export MARLIN_DLL=$PWD/libEutelescope.so.0.0.8"
doCommand "export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH"

# get the input raw file
doCommand "getFromGRID ${InputRawLFN} ${InputRawLocal}"
r=$?
if [ $r -ne 0 ] ; then
    echo "Problem copying ${InputRawLFN}. Exiting with error."
    exit 1
fi

# list all the files available
doCommand "ls -al"

# ready to run marlin
echo "Starting Marlin `date`"
c="Marlin $SteeringFile | tee $LogFile"
doCommand $c

# remove the input file
doCommand "rm ${InputRawLocal}"

# put back the files to the GRID
echo "Copying and registering the output file to SE"
doCommand "putOnGRID ${OutputLcioLocal} ${OutputLcioLFN} ${GRIDSE}"
r=$?
if [ $r -ne 0 ] ; then
    echo "Problem copying the ${OutputLcioLocal} to the GRID"
    exit 2
fi

echo "Preparing the joboutput tarball"
doCommand "tar czvf ${OutputJoboutputLocal} *.log *.xml"

echo "Copying and registering the tarball to SE"
doCommand "putOnGRID  ${OutputJoboutputLocal} ${OutputJoboutputLNF} ${GRIDSE}"
r=$?
if [ $r -ne 0 ] ; then
    echo "Problem copying the ${OutputJoboutputLocal} to the GRID"
    exit 3
fi

# Job finished
echo "Job finished at `date`"
