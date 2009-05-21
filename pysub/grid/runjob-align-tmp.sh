#!/bin/sh
# A template of alignment job
#
# @author Antonio Bulgheroni <mailto:antonio.bulgheroni@gmail.com>
# @version $Id: runjob-align-tmp.sh,v 1.2 2009-05-21 12:26:00 bulgheroni Exp $
#
# errno  0: No error.
# errno  1: Unable to get the input file from the SE.
# errno 20: Problem during Marlin execution.
# errno 30: Problem copying and registering the DB output to the SE.
# errno 31: Problem copying and registering the pede steering file
# errno 32: Problem copying and registering the mille binary file
# errno 33: Problem copying and registering the joboutput file
# errno 34: Problem copying and registering the histogram file
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
# errno 20: Problem during Marlin execution.
# errno 30: Problem copying and registering the DB output to the SE.
# errno 31: Problem copying and registering the pede steering file
# errno 32: Problem copying and registering the mille binary file
# errno 33: Problem copying and registering the joboutput file
# errno 34: Problem copying and registering the histogram file
#

# To be replaced with the output suffix
Output="@Output@"

# To be replace with the job name used for the identification of
# all files. It should be something like converter
Name="@Name@"

# This is the list of input files
InputFileList="@InputFileList@"

# To be replace with yes or no
RunPede="@RunPede@"

# Define here all the variables modified by the submitter
GRIDCE="@GRIDCE@"
GRIDSE="@GRIDSE@"
GRIDStoreProtocol="@GRIDStoreProtocol@"
GRIDVO="@GRIDVO@"
LFC_HOME="/grid/ilc/eudet-jra1"
GRIDFolderBase="@GRIDFolderBase@"
GRIDFolderHitmakerResults="@GRIDFolderHitmakerResults@"
GRIDFolderDBAlign="@GRIDFolderDBAlign@"
GRIDFolderAlignResults="@GRIDFolderAlignResults@"
GRIDFolderAlignJoboutput="@GRIDFolderAlignJoboutput@"
GRIDFolderAlignHisto="@GRIDFolderAlignHisto@"
GRIDLibraryTarball="@GRIDLibraryTarball@"
GRIDILCSoftVersion="@GRIDILCSoftVersion@"

# end of things to be replaced from the main script

OutputDBLFN=$GRIDFolderDBAlign/$Output-align-db.slcio
OutputSteerLFN=$GRIDFolderAlignResults/$Output-pede-steer.txt
OutputMilleLFN=$GRIDFolderAlignResults/$Output-align-mille.bin
OutputJoboutputLFN=$GRIDFolderAlignJoboutput/$Name-$Output.tar.gz
OutputHistoLFN=$GRIDFolderAlignHisto/$Output-align-histo.root

OutputDBLocal=$PWD/db/$Output-align-db.slcio
OutputSteerLocal=$PWD/results/$Output-pede-steer.txt
OutputMilleLocal=$PWD/results/$Output-align-mille.bin
OutputJoboutputLocal=$PWD/log/$Name-$Output.tar.gz
OutputHistoLocal=$PWD/histo/$Output-align-histo.root
SteeringFile=$Name-$Output.xml
LogFile=$Name-$Output.log

echo
echo "########################################################################"
echo "# Starting $Name-$Output at `date `"
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
doCommand "chmod 777 pede"

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
doCommand "export MARLIN_DLL=$PWD/libEutelescope.so.0.0.8"
doCommand "export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH"
doCommand "export PATH=$PWD:$PATH"

echo
echo "########################################################################"
echo "# Getting all the input files"
echo "########################################################################"
echo
echo
for file in $InputFileList; do 

    echo "--> $file"
    echo
    InputLFN=$GRIDFolderHitmakerResults/$file
    InputLocal=$PWD/results/$file

    doCommand "getFromGRID ${InputLFN} ${InputLocal}"
    r=$?
    if [ $r -ne 0 ] ; then
        echo "Problem copying ${InputLFN}. Exiting with error."
        exit 1
    fi

done 


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


# remove the input file
for file in $InputFileList; do
    InputLocal=$PWD/results/$file
    doCommand "rm ${InputLocal}"
done



# put back the files to the GRID
if [ $RunPede == "yes" ] ; then
    echo
    echo "########################################################################"
    echo "# Copying and registering the output DB file to SE"
    echo "########################################################################"
    echo
    doCommand "putOnGRID ${OutputDBLocal} ${OutputDBLFN} ${GRIDSE}"
    r=$?
    if [ $r -ne 0 ] ; then
        echo "****** Problem copying the ${OutputDBLocal} to the GRID"
        exit 30
    fi
fi

echo
echo "########################################################################"
echo "# Copying and registering the output pede steering to SE"
echo "########################################################################"
echo
doCommand "putOnGRID ${OutputSteerLocal} ${OutputSteerLFN} ${GRIDSE}"
r=$?
if [ $r -ne 0 ] ; then
    echo "****** Problem copying the ${OutputSteerLocal} to the GRID"
    exit 31
fi

echo
echo "########################################################################"
echo "# Copying and registering the output mille binary to SE"
echo "########################################################################"
echo
doCommand "putOnGRID ${OutputMilleLocal} ${OutputMilleLFN} ${GRIDSE}"
r=$?
if [ $r -ne 0 ] ; then
    echo "****** Problem copying the ${OutputMilleLocal} to the GRID"
    exit 32
fi

echo
echo "########################################################################"
echo "# Preparing the joboutput tarball"
echo "########################################################################"
echo
doCommand "tar czvf ${OutputJoboutputLocal} *.log *.xml out err mille* histo/*root"

echo
echo "########################################################################"
echo "# Copying and registering the tarball to SE"
echo "########################################################################"
echo
doCommand "putOnGRID  ${OutputJoboutputLocal} ${OutputJoboutputLFN} ${GRIDSE}"
r=$?
if [ $r -ne 0 ] ; then
    echo "****** Problem copying the ${OutputJoboutputLocal} to the GRID"
    exit 33
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
    exit 34
fi

# Job finished
echo
echo "########################################################################"
echo "# Job ($Name-$Output) finished at `date`"
echo "########################################################################"
echo

