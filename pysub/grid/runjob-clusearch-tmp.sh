#!/bin/sh
# A template of pedestal job
#
# @author Antonio Bulgheroni <mailto:antonio.bulgheroni@gmail.com>
# @version $Id: runjob-clusearch-tmp.sh,v 1.5 2009-07-28 15:54:50 bulgheroni Exp $
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

# To be replaced with the hotpixelString in the format %(pede)s (does not have
# to be integer)
HotPixelRunString="@HotPixelRunNumber@"

# To be replace with the job name used for the identification of
# all files. It should be something like converter
Name="@Name@"

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
GRIDFolderLcioRaw="@GRIDFolderLcioRaw@"
GRIDFolderDBOffset="@GRIDFolderDBOffset@"
GRIDFolderDBPede="@GRIDFolderDBPede@"
GRIDFolderClusearchResults="@GRIDFolderClusearchResults@"
GRIDFolderClusearchJoboutput="@GRIDFolderClusearchJoboutput@"
GRIDFolderClusearchHisto="@GRIDFolderClusearchHisto@"
GRIDILCSoftVersion="@GRIDILCSoftVersion@"

# GRID Tarball
# LocalGRIDLibraryTarball --> "yes" means that the tarball is uploaded along with the JDL file
#                         --> "no" means that it has to be downloaded from a SE
HasLocalGRIDLibraryTarball="@HasLocalGRIDLibraryTarball@"
GRIDLibraryTarball="@GRIDLibraryTarball@"
GRIDLibraryTarballPath="@GRIDLibraryTarballPath@"
GRIDLibraryLocal=$PWD/$GRIDLibraryTarball
GRIDLibraryLFN=$GRIDLibraryTarballPath/$GRIDLibraryTarball

if [ $IsDUTOnly == "yes" ] ; then

    # lfn
    OutputLcioLFN=$GRIDFolderClusearchResults/run$RunString-clu-$DUTSuffix-p$PedeString.slcio
    OutputHistoLFN=$GRIDFolderClusearchHisto/run$RunString-clu-$DUTSuffix-histo.root
    OutputJoboutputLFN=$GRIDFolderClusearchJoboutput/$Name-$RunString-$DUTSuffix.tar.gz

    # local
    OutputLcioLocal=$PWD/results/run$RunString-clu-$DUTSuffix-p$PedeString.slcio
    OutputHistoLocal=$PWD/histo/run$RunString-clu-$DUTSuffix-histo.root
    OutputJoboutputLocal=$PWD/log/$Name-$RunString-$DUTSuffix.tar.gz

else

    # lfn
    OutputLcioLFN=$GRIDFolderClusearchResults/run$RunString-clu-p$PedeString.slcio
    OutputLcioLFNShort=$GRIDFolderClusearchResults/run$RunString-clu-p$PedeString
    OutputLcioLFNDIR=$GRIDFolderClusearchResults/
    OutputHistoLFN=$GRIDFolderClusearchHisto/run$RunString-clu-histo.root
    OutputJoboutputLFN=$GRIDFolderClusearchJoboutput/$Name-$RunString.tar.gz
    OutputDBLFN=$GRIDFolderDBOffset/run$RunString-offset-db.slcio


    # local
    OutputLcioLocal=$PWD/results/run$RunString-clu-p$PedeString.slcio
    OutputLcioLocalShort=$PWD/results/run$RunString-clu-p$PedeString
    OutputLcioLocalDIR=$PWD/results/
    OutputHistoLocal=$PWD/histo/run$RunString-clu-histo.root
    OutputJoboutputLocal=$PWD/log/$Name-$RunString.tar.gz
    OutputDBLocal=$PWD/db/run$RunString-offset-db.slcio


fi

# lfn
InputLcioRawLFN=$GRIDFolderLcioRaw/run$RunString.slcio
InputPedeLFN=$GRIDFolderDBPede/run$PedeString-ped-db.slcio
InputHotPixelLFN=$GRIDFolderDBPede/run$HotPixelRunString-hotpixel-db.slcio

# local
InputLcioRawLocal=$PWD/lcio-raw/run$RunString.slcio
InputPedeLocal=$PWD/db/run$PedeString-ped-db.slcio
InputHotPixelLocal=$PWD/db/run$HotPixelRunString-hotpixel-db.slcio

# Steering and log
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
#doCommand ". $VO_ILC_SW_DIR/initILCSOFT.sh $GRIDILCSoftVersion"
echo "where is VO_ILC_SW_DIR: $VO_ILC_SW_DIR"
doCommand ". $VO_ILC_SW_DIR/ilcsoft/x86_64_gcc41_sl5/init_ilcsoft.sh $GRIDILCSoftVersion"
#doCommand "$BASH ./ilc-grid-test-sw.sh"

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
echo "# Getting the input file ${InputLcioRawLocal}"
echo "########################################################################"
echo
# get the input raw file
doCommand "getFromGRID ${InputLcioRawLFN} ${InputLcioRawLocal}"
r=$?
if [ $r -ne 0 ] ; then
    echo "Problem copying ${InputLcioRawLFN}. Exiting with error."
    exit 2
fi

echo
echo "########################################################################"
echo "# Getting the hotpixel DB file ${InputHotPixelLocal}"
echo "########################################################################"
echo
# 
doCommand "getFromGRID ${InputHotPixelLFN} ${InputHotPixelLocal}"
r=$?
if [ $r -ne 0 ] ; then
   echo "Problem copying ${InputHotPixelLFN}.    Warning! No hotpixel db file found at given location."
#    echo "Please, check your input to make sure you all values are set properly. If you do not want to apply hotpixel masks you can ignore this warning."
###    exit 3
fi
#

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
touch temp.err

echo                                                                            | tee -a temp.err
echo "########################################################################" | tee -a temp.err
echo "# Starting Marlin `date`"                                                 | tee -a temp.err
echo "########################################################################" | tee -a temp.err
echo                                                                            | tee -a temp.err
find ${OutputLcioLocalDIR} | tee -a temp.err
find ./                    | tee -a temp.err
#mail -s "$Name; run=$RunString; starting Marlin " ${USER}@mail.desy.de        < temp.err

c="Marlin $SteeringFile"
echo $c
$c
r=$?

cat "" > temp.err
#mail -s "$Name; run=$RunString; Marlin $r " ${USER}@mail.desy.de        < temp.err

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
#doCommand "hadd -f temp.root empty.root ${OutputHistoLocal}"
#doCommand "mv temp.root ${OutputHistoLocal}"

# put back the files to the GRID
echo
echo "########################################################################"
echo "# Copying and registering the output file to SE"
echo "########################################################################"
echo
cat "" > temp.err
find ${OutputLcioLocalDIR} | tee -a temp.err
find ./                    | tee -a temp.err
#mail -s "$Name; run=$RunString; copy files to GRID" ${USER}@mail.desy.de        < temp.err

doCommand "putOnGRID ${OutputLcioLocal} ${OutputLcioLFN} ${GRIDSE}"
r=$?
if [ $r -ne 0 ] ; then
    echo "****** Problem copying the ${OutputLcioLocal} to the GRID"
    ls ${OutputLcioLocalShort}.*.slcio | tee slcio.txt
    runs=`grep -c slcio slcio.txt`
    echo "*** $runs found "

    for (( i=0; i<$runs; i++ ))
    do
      printf -v ii "%03i" $i
      echo "current run = $ii"
      OutputLocal=${OutputLcioLocalShort}.$ii.slcio
      OutputLFN=${OutputLcioLFNShort}.$ii.slcio
      echo "putOnGRID ${OutputLocal} ${OutputLFN} ${GRIDSE}"
      doCommand "putOnGRID ${OutputLocal} ${OutputLFN} ${GRIDSE}"
      r=$?
      if [ $r -ne 0 ] ; then
        echo "****** Problem copying $runs files from ${OutputLcioLocalShort}.xxx.slcio to the GRID"
        exit 30
      fi  
    done
#   exit 30
fi

# put back the files to the GRID
echo
echo "########################################################################"
echo "# Copying and registering the output DB offset file to SE"
echo "########################################################################"
echo "lfc-ls -l ${OutputDBLFN}"
lfc-ls -l ${OutputDBLFN}
echo "ls -ltr ${OutputDBLocal}"
ls -ltr ${OutputDBLocal}
echo doCommand "putOnGRID ${OutputDBLocal} ${OutputDBLFN} ${GRIDSE}"
 
doCommand "putOnGRID ${OutputDBLocal} ${OutputDBLFN} ${GRIDSE}"
r=$?
if [ $r -ne 0 ] ; then
    echo "****** Problem copying the ${OutputDBLocal} to the GRID"
#    exit 30
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
#    exit 31
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
#    exit 32
fi

# Job finished
echo
echo "########################################################################"
echo "# Job ($Name-$RunString) finished at `date`"
echo "########################################################################"
echo

