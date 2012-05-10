#!/bin/sh
# A template of alignment job
#
# @author Antonio Bulgheroni <mailto:antonio.bulgheroni@gmail.com>
# @version $Id: runjob-align-tmp.sh,v 1.9 2009-07-28 00:13:59 bulgheroni Exp $
#
# errno  0: No error.
# errno  1: Unable to get the GRID library tarball from the SE
# errno  2: Unable to get the input file from the SE.
# errno 20: Problem during Marlin execution.
# errno 21: Problem during pede execution.
# errno 30: Problem copying and registering the DB output to the SE.
# errno 31: Problem copying and registering the pede steering file.
# errno 32: Problem copying and registering the mille binary file.
# errno 33: Problem copying and registering the joboutput file.
# errno 34: Problem copying and registering the histogram file.
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
# errno 20: Problem during Marlin execution.
# errno 21: Problem during pede execution.
# errno 30: Problem copying and registering the DB output to the SE.
# errno 31: Problem copying and registering the pede steering file.
# errno 32: Problem copying and registering the mille binary file.
# errno 33: Problem copying and registering the joboutput file.
# errno 34: Problem copying and registering the histogram file.
#

# To be replaced with the output suffix
Output="@Output@"
HotPixelRunNumber="@Output@"


# To be replace with the job name used for the identification of
# all files. It should be something like converter
Name="@Name@"

# This is the list of input files
InputFileList="@InputFileList@"

# To be replace with yes or no
RunPede="@RunPede@"
PedeSteerTemplate="@PedeSteerTemplate@"

# Define here all the variables modified by the submitter
GRIDCE="@GRIDCE@"
GRIDSE="@GRIDSE@"
GRIDStoreProtocol="@GRIDStoreProtocol@"
GRIDVO="@GRIDVO@"
LFC_HOME="/grid/ilc/eudet-jra1"
GRIDFolderBase="@GRIDFolderBase@"
GRIDFolderHitmakerResults="@GRIDFolderHitmakerResults@"
GRIDFolderDBAlign="@GRIDFolderDBAlign@"
GRIDFolderDBHotPixel="@GRIDFolderDBHotPixel@"
GRIDFolderAlignResults="@GRIDFolderAlignResults@"
GRIDFolderAlignJoboutput="@GRIDFolderAlignJoboutput@"
GRIDFolderAlignHisto="@GRIDFolderAlignHisto@"
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

OutputDBLFN=$GRIDFolderDBAlign/$Output-align-db.slcio
OutputSteerLFN=$GRIDFolderAlignResults/$Output-pede-steer.txt
OutputMilleLFN=$GRIDFolderAlignResults/$Output-align-mille.bin
OutputJoboutputLFN=$GRIDFolderAlignJoboutput/$Name-$Output.tar.gz
OutputHistoLFN=$GRIDFolderAlignHisto/$Output-align-histo.root

OutputDBLocal=$PWD/db/$Output-align-db.slcio
OutputSteerLocal=$PWD/results/$Output-pede-steer.txt
OutputMilleLocal=$PWD/results/$Output-align-mille.bin
OutputMilleShort=results/$Output-align-mille.bin
OutputJoboutputLocal=$PWD/log/$Name-$Output.tar.gz
OutputHistoLocal=$PWD/histo/$Output-align-histo.root
SteeringFile=$Name-$Output.xml
LogFile=$Name-$Output.log

InputHotPixelLFN=$GRIDFolderDBHotPixel/run$HotPixelRunNumber-hotpixel-db.slcio
InputHotPixelLocal=$PWD/db/run$HotPixelRunNumber-hotpixel-db.slcio
LocalPWD=$PWD


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
doCommand "chmod 777 pede pede2lcio"

# rename the simjob.slcio because otherwise it gets delete
doCommand "mv simjob.slcio simjob.slcio.keepme"

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

# now it's safe to rename the simjob to the original
doCommand "mv simjob.slcio.keepme simjob.slcio"

# set the list of Marlin plugins and the LD_LIBRARY_PATH
doCommand "export MARLIN_DLL=$PWD/libEutelescope.so"
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
#    echo "Please, check your input to make sure you all values are set properly. If you do not want to apply hotpixel masks you can ignore this warning."
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

#
# sometimes this iteration of pede (on GRID only!) misbehaves 
# cure: comment this condition OUT 
# This is a temp solution!
#
#if [ $RunPede == "yes" ] ; then
if [ 1 -eq 0 ] ; then

    echo
    echo "########################################################################"
    echo "# Rerunning pede `date`"
    echo "########################################################################"
    echo

    echo "Preparing the new steering file ${PedeSteerTemplate}"
    if  test -f $PedeSteerTemplate ; then
        echo "Pede steering template found"
     else
        echo "****** Missing pede steering template"
        exit 22
    fi

    echo "Replacing the mille binary file"
    sed -e "s|@MilleBinFile@|$OutputMilleShort|" < $PedeSteerTemplate > pede-steer-tmp.working

    echo "Replacing the starting points"
    if test -f millepede.res ; then 
        echo "Previous pede iteration results file found"
    else 
        echo "****** Missing millepede.res file"
        exit 23
    fi
    cat millepede.res | while read line ; do
        value=$line"\n@Parameters@"
        sed -e "s|@Parameters@|$value|" < pede-steer-tmp.working > pede-steer-tmp.working2
        mv pede-steer-tmp.working2 pede-steer-tmp.working
    done

    # now remove any other additional @Parameters@
    echo "Cleaning up the steering file"
    sed -e "s|@Parameters@||" < pede-steer-tmp.working > pede-steer-tmp.final

    echo "New steering file ready to be used"
    # now rename the initial pede steering file in something else, that will be archieved
    doCommand "mv ${OutputSteerLocal} $PWD/results/$Output-pede-steer-iter-0.txt"

    # rename the final steering in the current one
    doCommand "mv pede-steer-tmp.final ${OutputSteerLocal}"

    # remove working copy of the steering template
    doCommand "rm pede-steer-tmp.working*"

    # running pede
    # I'm not sure if pede is really using the return value
    # but anyway read it back
    c="pede ${OutputSteerLocal}"
    echo $c
    $c
    r=$?
    if [ $r -ne 0 ] ; then
        echo "****** Problem running pede"
        exit 21
    fi

    # copying the old db file in a backup file 
    doCommand "mv db/$Output-align-db.slcio db/$Output-align-db-iter-0.slcio"
    
    # run the pede2lcio utility
    doCommand "pede2lcio millepede.res db/$Output-align-db.slcio"
    r=$?
    if [ $r -ne 0 ] ; then
        echo "****** Problem running pede2lcio"
        exit 24
    fi

    echo
    echo "########################################################################"
    echo "# Second iteration of  pede finished at `date`"
    echo "########################################################################"
    echo
fi 

# remove the input file
for file in $InputFileList; do
    InputLocal=$PWD/results/$file
    doCommand "rm ${InputLocal}"
done

# fixing the problem with the histogram file
#doCommand "hadd -f temp.root empty.root ${OutputHistoLocal}"
#doCommand "mv temp.root ${OutputHistoLocal}"

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
doCommand "tar czvf ${OutputJoboutputLocal} *.log *.xml *.txt out err mille* db/*iter* histo/*root"

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

