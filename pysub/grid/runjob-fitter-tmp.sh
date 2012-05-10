#!/bin/sh
# A template of fitter job
#
# @author Antonio Bulgheroni <mailto:antonio.bulgheroni@gmail.com>
# @version $Id: runjob-fitter-tmp.sh,v 1.4 2009-07-28 00:13:59 bulgheroni Exp $
#
# errno  0: No error.
# errno  1: Unable to get the GRID library tarball from the SE
# errno  2: Unable to get the input file from the SE.
# errno  3: Unable to get the alignment file from the SE.
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
    doCommand "lcg-del -a lfn:$2 "
    doCommand "lcg-cr -v -d $3 -l lfn:$2 file:$1"
    return $?
}


## Run the job
#
# This is the real executable. It doesn't take any argument and
# it return 0 in case of successful execution or the following error
# codes in case of problems
#
# errno  1: Unable to get the GRID library tarball from the SE
# errno  2: Unable to get the input file from the SE.
# errno  3: Unable to get the alignment file from the SE.
# errno 20: Problem during Marlin execution
# errno 30: Problem copying and registering the LCIO output to the SE
# errno 31: Problem copying and registering the Joboutput to the SE


# To be replaced with the output suffix
Output="@Output@"

# To be replace with the job name used for the identification of
# all files. It should be something like fitter
Name="@Name@"

HotPixelRunNumber="@Output@"


# This is the list of input files
InputFileList="@InputFileList@"

# This is the alignment file w/o path
AlignFile="@AlignFile@"

# Define here all the variables modified by the submitter
GRIDCE="@GRIDCE@"
GRIDSE="@GRIDSE@"
GRIDStoreProtocol="@GRIDStoreProtocol@"
GRIDVO="@GRIDVO@"
LFC_HOME="/grid/ilc/eudet-jra1"
GRIDFolderBase="@GRIDFolderBase@"
GRIDFolderHitmakerResults="@GRIDFolderHitmakerResults@"
GRIDFolderDBAlign="@GRIDFolderDBAlign@"
GRIDFolderFitterResults="@GRIDFolderFitterResults@"
GRIDFolderFitterJoboutput="@GRIDFolderFitterJoboutput@"
GRIDFolderFitterHisto="@GRIDFolderFitterHisto@"
GRIDILCSoftVersion="@GRIDILCSoftVersion@"

# GRID Tarball
# LocalGRIDLibraryTarball --> "yes" means that the tarball is uploaded along with the JDL file
#                         --> "no" means that it has to be downloaded from a SE
HasLocalGRIDLibraryTarball="@HasLocalGRIDLibraryTarball@"
GRIDLibraryTarball="@GRIDLibraryTarball@"
GRIDLibraryTarballPath="@GRIDLibraryTarballPath@"
GRIDLibraryLocal=$PWD/$GRIDLibraryTarball
GRIDLibraryLFN=$GRIDLibraryTarballPath/$GRIDLibraryTarball

# end of things to be replaced from the main script

InputAlignLFN=$GRIDFolderDBAlign/$AlignFile
OutputJoboutputLFN=$GRIDFolderFitterJoboutput/$Name-$Output.tar.gz
OutputHistoLFN=$GRIDFolderFitterHisto/$Output-track-histo.root

InputAlignLocal=$PWD/db/$AlignFile
OutputJoboutputLocal=$PWD/log/$Name-$Output.tar.gz
OutputHistoLocal=$PWD/histo/$Output-track-histo.root
SteeringFile=$Name-$Output.xml
LogFile=$Name-$Output.log

InputHotPixelLFN=$GRIDFolderDBHotPixel/run$HotPixelRunNumber-hotpixel-db.slcio
InputHotPixelLocal=$PWD/db/run$HotPixelRunNumber-hotpixel-db.slcio
LocalPWD=$PWD


OutputTBTrackLocal=$PWD/histo/tbtrack$Output.root
OutputTBTrackLFN=$GRIDFolderFitterHisto/tbtrack$Output.root

RefhitLFN=$GRIDFolderDBOffset/$Output-refhit-db.slcio
RefhitLocal=$PWD/db/$Output-refhit-db.slcio


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
echo doCommand ". $VO_ILC_SW_DIR/initILCSOFT.sh $GRIDILCSoftVersion"
echo "ILCSOFT : "$ILCSOFT
#


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
for file in $InputFileList; do 

    echo "--> $file"
    echo
    InputLFN=$GRIDFolderHitmakerResults/$file
    InputLocal=$PWD/results/$file

    doCommand "getFromGRID ${InputLFN} ${InputLocal}"
    r=$?
    if [ $r -ne 0 ] ; then
        echo "Problem copying ${InputLFN}. Exiting with error."
        exit 2
    fi

done 


RefhitLFN=$GRIDFolderDBAlign/$Output-refhit-db.slcio
RefhitLocal=$PWD/db/$Output-refhit-db.slcio

echo
echo "########################################################################"
echo "# Getting the refhit DB file ${PreAlignLocal}"
echo "########################################################################"
echo
doCommand "getFromGRID ${RefhitLFN} ${RefhitLocal}"
r=$?
if [ $r -ne 0 ] ; then
    echo "Problem copying ${RefhitLFN}. Exiting with error."
    exit 3
fi


PreAlignLFN=$GRIDFolderDBAlign/$Output-prealign-db.slcio
PreAlignLocal=$PWD/db/$Output-prealign-db.slcio

echo
echo "########################################################################"
echo "# Getting the prealign DB file ${PreAlignLocal}"
echo "########################################################################"
echo
doCommand "getFromGRID ${PreAlignLFN} ${PreAlignLocal}"
r=$?
if [ $r -ne 0 ] ; then
    echo "Problem copying ${PreAlignLFN}. Exiting with error."
    exit 3
fi

# list all the files available
doCommand "ls -al"

echo
echo "########################################################################"
echo "# Getting the align DB file ${InputAlignLocal}"
echo "########################################################################"
echo
doCommand "getFromGRID ${InputAlignLFN} ${InputAlignLocal}"
r=$?
if [ $r -ne 0 ] ; then
    echo "Problem copying ${InputAlignLFN}. Exiting with error."
    exit 3
fi


echo
echo "########################################################################"
echo "# Getting the hot pixel db file ${InputHotPixelLocal}"
echo "########################################################################"
echo
#  

doCommand "getFromGRID  ${InputHotPixelLFN} ${InputHotPixelLocal}"
r=$?
if [ $r -ne 0 ] ; then
   echo "Problem copying ${InputHotPixelLocal} into ${InputHotPixelLFN}.    Warning! No hotpixel db file found at given location."
   echo "Please, check your input to make sure you all values are set properly. If you do not want to apply hotpixel masks you can ignore this warning."
###    exit 3
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


# remove the input file
for file in $InputFileList; do
    InputLocal=$PWD/results/$file
    doCommand "rm ${InputLocal}"
done

# remove the eta file
doCommand "rm ${InputAlignLocal}"

# fix the histogram file
#doCommand "hadd -f temp.root empty.root ${OutputHistoLocal}"
#doCommand "mv temp.root ${OutputHistoLocal}"

# put back the files to the GRID
echo
echo "########################################################################"
echo "# Copying and registering the output file to SE"
echo "########################################################################"
echo
# the output LCIO can be more than one.
# just loop over all possible values!
OutputLocal=$PWD/results/$Output-track.slcio
OutputLFN=$GRIDFolderFitterResults/$Output-track.slcio

doCommand "putOnGRID ${OutputLocal} ${OutputLFN} ${GRIDSE}"
r=$?
if [ $r -eq 0 ] ; then
  echo "****** ${OutputLocal} copied to the GRID"

else
  echo "****** Problem copying the ${OutputLocal} to the GRID, try with .xxx.slcio "

  for file in `ls $PWD/results/$Output-track.[0-9][0-9][0-9].slcio` ; do
    
    basefile=`basename $file`
    echo "--> $basefile"
    echo

    OutputLocal=$file
    OutputLFN=$GRIDFolderFitterResults/$basefile

    doCommand "putOnGRID ${OutputLocal} ${OutputLFN} ${GRIDSE}"
    r=$?
    if [ $r -ne 0 ] ; then
        echo "****** Problem copying the ${OutputLocal} to the GRID"
        exit 30
    fi
  done
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
    ls -ltrah ${OutputHistoLocal} 
    find ./
    echo "****** Problem copying the ${OutputHistoLocal} to the GRID"
    exit 32
fi

doCommand "putOnGRID ${OutputTBTrackLocal} ${OutputTBTrackLFN} ${GRIDSE}"
r=$?
if [ $r -ne 0 ] ; then
    ls -ltrah ${OutputTBTrackLocal}
    find ./
    echo "****** Problem copying the ${OutputTBTrackLocal} to the GRID"
    exit 32
fi


# Job finished
echo
echo "########################################################################"
echo "# Job ($Name-$Output) finished at `date`"
echo "########################################################################"
echo

