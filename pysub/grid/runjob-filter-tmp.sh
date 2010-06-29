#!/bin/sh
# A template of filtering job
#
# @author Antonio Bulgheroni <mailto:antonio.bulgheroni@gmail.com>
# @version $Id: runjob-filter-tmp.sh,v 1.8 2009-08-01 10:44:30 bulgheroni Exp $
#
# errno  0: No error.
# errno  1: Unable to get the GRID library tarball from the SE
# errno  2: Unable to get the input file from the SE.
# errno  3: Unable to get the pedestal file from the SE.
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
# errno  1: Unable to get the GRID library tarball from the SE
# errno  2: Unable to get the input file from the SE.
# errno  3: Unable to get the pedestal file from the SE.
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

# Run list (actually used only when merging...)
RunList="@RunList@"

# Merging ?
Merge="@Merge@"

# Base name (actually used only when merging...)
Base="@Base@"

# Depending on the presence or not of the DUT, set the following
# variables
IsTelescopeOnly="@IsTelescopeOnly@"
IsTelescopeAndDUT="@IsTelescopeAndDUT@"
IsDUTOnly="@IsDUTOnly@"

# DUTSuffix
DUTSuffix="@DUTSuffix@"

# Define here all the variables modified by the submitter
GRIDCE="@GRIDCE@"
GRIDSE="@GRIDSE@"
GRIDStoreProtocol="@GRIDStoreProtocol@"
GRIDVO="@GRIDVO@"
LFC_HOME="/grid/ilc/eudet-jra1"
GRIDFolderBase="@GRIDFolderBase@"
GRIDFolderClusearchResults="@GRIDFolderClusearchResults@"
GRIDFolderDBPede="@GRIDFolderDBPede@"
GRIDFolderFilterResults="@GRIDFolderFilterResults@"
GRIDFolderFilterJoboutput="@GRIDFolderFilterJoboutput@"
GRIDFolderFilterHisto="@GRIDFolderFilterHisto@"
GRIDILCSoftVersion="@GRIDILCSoftVersion@"

# GRID Tarball
# LocalGRIDLibraryTarball --> "yes" means that the tarball is uploaded along with the JDL file
#                         --> "no" means that it has to be downloaded from a SE
HasLocalGRIDLibraryTarball="@HasLocalGRIDLibraryTarball@"
GRIDLibraryTarball="@GRIDLibraryTarball@"
GRIDLibraryTarballPath="@GRIDLibraryTarballPath@"
GRIDLibraryLocal=$PWD/$GRIDLibraryTarball
GRIDLibraryLFN=$GRIDLibraryTarballPath/$GRIDLibraryTarball

InputPedeLFN=$GRIDFolderDBPede/run$PedeString-ped-db.slcio
InputPedeLocal=$PWD/db/run$PedeString-ped-db.slcio

# This runjob needs to have a special definition of input files because they could be many!
# moreover the naming convention is rather different from the standard one when merging. 

if [ $Merge == "no" ] ; then

    if [ $IsDUTOnly == "yes" ] ; then
        OutputLcioLFN=$GRIDFolderFilterResults/run$RunString-filter-$DUTSuffix-p$PedeString.slcio
        OutputJoboutputLFN=$GRIDFolderFilterJoboutput/$Name-$RunString-$DUTSuffix.tar.gz
        OutputHistoLFN=$GRIDFolderFilterHisto/run$RunString-filter-$DUTSuffix-histo.root
#
        OutputLcioLocal=$PWD/results/run$RunString-filter-$DUTSuffix-p$PedeString.slcio
        OutputJoboutputLocal=$PWD/log/$Name-$RunString-$DUTSuffix.tar.gz
        OutputHistoLocal=$PWD/histo/run$RunString-filter-$DUTSuffix-histo.root
    else
        OutputLcioLFN=$GRIDFolderFilterResults/run$RunString-filter-p$PedeString.slcio
        OutputJoboutputLFN=$GRIDFolderFilterJoboutput/$Name-$RunString.tar.gz
        OutputHistoLFN=$GRIDFolderFilterHisto/run$RunString-filter-histo.root
#
        OutputLcioLocal=$PWD/results/run$RunString-filter-p$PedeString.slcio
        OutputJoboutputLocal=$PWD/log/$Name-$RunString.tar.gz
        OutputHistoLocal=$PWD/histo/run$RunString-filter-histo.root
    fi

    SteeringFile=$Name-$RunString.xml
    LogFile=$Name-$RunString.log
else

    if [ $IsDUTOnly == "yes" ] ; then
        OutputLcioLFN=$GRIDFolderFilterResults/$Base-filter-$DUTSuffix-p$PedeString.slcio
        OutputJoboutputLFN=$GRIDFolderFilterJoboutput/$Name-$Base-$DUTSuffix.tar.gz
        OutputHistoLFN=$GRIDFolderFilterHisto/$Base-filter-$DUTSuffix-histo.root
#
        OutputLcioLocal=$PWD/results/$Base-filter-$DUTSuffix-p$PedeString.slcio
        OutputJoboutputLocal=$PWD/log/$Name-$Base-$DUTSuffix.tar.gz
        OutputHistoLocal=$PWD/histo/$Base-filter-$DUTSuffix-histo.root
    else
        OutputLcioLFN=$GRIDFolderFilterResults/$Base-filter-p$PedeString.slcio
        OutputJoboutputLFN=$GRIDFolderFilterJoboutput/$Name-$Base.tar.gz
        OutputHistoLFN=$GRIDFolderFilterHisto/$Base-filter-histo.root
#
        OutputLcioLocal=$PWD/results/$Base-filter-p$PedeString.slcio
        OutputJoboutputLocal=$PWD/log/$Name-$Base.tar.gz
        OutputHistoLocal=$PWD/histo/$Base-filter-histo.root
    fi

    SteeringFile=$Name-$Base.xml
    LogFile=$Name-$Base.log
fi

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

# check if we need to get the tarbal or not
if [ $HasLocalGRIDLibraryTarball == "no" ] ; then

    echo
    echo "########################################################################"
    echo "# Getting the lib tarball..."
    echo "########################################################################"
    echo

    doCommand "getFromGRID ${GRIDLibraryLFN} ${GRIDLibraryLocal} "
    r=$?
    if [ $r -ne 0 ] ; then
        echo "Problem copying ${GRIDLibraryLFN}. Exiting with error"
        exit 1
    fi
fi

# unpack the library
echo
echo "########################################################################"
echo "# Uncompressing the job tarball..."
echo "########################################################################"
echo
doCommand "tar xzvf $GRIDLibraryTarball"


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


# set the list of Marlin plugins and the LD_LIBRARY_PATH
doCommand "export MARLIN_DLL=$PWD/libEutelescope.so"
doCommand "export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH"

echo
echo "########################################################################"
echo "# Getting all the input files"
echo "########################################################################"
echo
echo
for run in $RunList; do

    echo "--> $run"
    echo
    InputLcioLFN=$GRIDFolderClusearchResults/run$run-clu-p$PedeString.slcio
    InputLcioLocal=$PWD/results/run$run-clu-p$PedeString.slcio

    doCommand "getFromGRID ${InputLcioLFN} ${InputLcioLocal}"
    r=$?
    if [ $r -ne 0 ] ; then
        echo "Problem copying ${InputLcioLFN}. Exiting with error."
        exit 2
    fi

done

echo
echo "########################################################################"
echo "# Getting the pedestal DB file ${InputPedeLocal}"
echo "########################################################################"
echo
doCommand "getFromGRID ${InputPedeLFN} ${InputPedeLocal}"
r=$?
if [ $r -ne 0 ] ; then
   echo "Problem copying ${InputPedeLFN}. Warning! If the data under question is
   from an Analog source this Warning should be treated as an ERROR."
#    exit 3
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

for run in $RunList; do
    InputLcioLocal=$PWD/results/run$run-clu-p$PedeString.slcio
    doCommand "rm ${InputLcioLocal}"
done

doCommand "rm ${InputPedeLocal}"

# fixing the problem with the histogram file
# i.e. hadd the output file with an empty one.
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

