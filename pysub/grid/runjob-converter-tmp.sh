#!/bin/sh
# A template of converter job
#
# @author Antonio Bulgheroni <mailto:antonio.bulgheroni@gmail.com>
# @version $Id: runjob-converter-tmp.sh,v 1.10 2009-07-28 00:13:59 bulgheroni Exp $
#
# errno  0: No error.
# errno  1: Unable to get the GRID library tarball from the SE
# errno  2: Unable to get the input file from the SE.
# errno 20: Problem during Marlin execution.
# errno 30: Problem copying and registering the LCIO output to the SE.
# errno 31: Problem copying and registering the Joboutput to the SE.
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
#    echo      "lcg-cp -v lfn:$1 file:$2"    
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
    echo      "lcg-cr -v -d $3 -l lfn:$2 file:$1"   
    if [ -e $1 ] ; then
      doCommand "lcg-del -a lfn:$2 "
      doCommand "lcg-cr -v -d $3 -l lfn:$2 file:$1"
    else
      for i in $(seq 1 999)
      do
         stringZfrom=$1
         stringZto=$2
         
         if [ $i -lt 10 ] ; then
           a="00"$i
         elif [ $i -lt 100 ] ; then
           a="0"$i
         elif [ $i -lt 1000 ] ; then
           a=$i
         fi  
          
         stringAfrom=${stringZfrom/slcio/$a".slcio"}
         stringAto=${stringZto/slcio/$a".slcio"}
         
         if [ -e $stringAfrom ] ;then
#           echo "doCommand lcg-cr -v -d $3 -l lfn:$stringAto file:$stringAfrom"            
           doCommand "lcg-del -a lfn:$stringAto "
           doCommand "lcg-cr -v -d $3 -l lfn:$stringAto file:$stringAfrom"            
         fi
#         echo $stringA         
      done      
    fi
    
    return $?
}


## Run the job
#
# This is the real executable. It doesn't take any argument and
# it return 0 in case of successful execution or the following error
# codes in case of problems
#
#
# errno  1: Unable to get the GRID library tarball from the SE
# errno  2: Unable to get the input file from the SE.
# errno 20: Problem during Marlin execution
# errno 30: Problem copying and registering the LCIO output to the SE
# errno 31: Problem copying and registering the Joboutput to the SE


# To be replaced with the runString in the format %(run)06d
RunString="@RunString@"

# To be replaced with the RunHotPixelNumber value in the format %(run)s can be a string
HotPixelRunNumber="@HotPixelRunNumber@"


# To be replace with the job name used for the identification of
# all files. It should be something like converter
Name="@Name@"

# Define here all the variables modified by the submitter
GRIDCE="@GRIDCE@"
GRIDSE="@GRIDSE@"
GRIDStoreProtocol="@GRIDStoreProtocol@"
GRIDVO="@GRIDVO@"
LFC_HOME="/grid/ilc/aida-wp9"
GRIDFolderBase="@GRIDFolderBase@"
GRIDFolderNative="@GRIDFolderNative@"
GRIDFolderLcioRaw="@GRIDFolderLcioRaw@"
GRIDFolderDBHotPixel="@GRIDFolderDBHotPixel@"
GRIDFolderConvertJoboutput="@GRIDFolderConvertJoboutput@"
GRIDILCSoftVersion="@GRIDILCSoftVersion@"

# GRID Tarball
# LocalGRIDLibraryTarball --> "yes" means that the tarball is uploaded along with the JDL file
#                         --> "no" means that it has to be downloaded from a SE
HasLocalGRIDLibraryTarball="@HasLocalGRIDLibraryTarball@"
GRIDLibraryTarball="@GRIDLibraryTarball@"
GRIDLibraryTarballPath="@GRIDLibraryTarballPath@"
GRIDLibraryLocal=$PWD/$GRIDLibraryTarball
GRIDLibraryLFN=$GRIDLibraryTarballPath/$GRIDLibraryTarball

# I / O files (Local and LFN)
InputRawLFN=$GRIDFolderNative/run$RunString.raw
OutputLcioLFN=$GRIDFolderLcioRaw/run$RunString.slcio
OutputJoboutputLFN=$GRIDFolderConvertJoboutput/$Name-$RunString.tar.gz
InputHotPixelLFN=$GRIDFolderDBHotPixel/run$HotPixelRunNumber-hotpixel-db.slcio

InputRawLocal=$PWD/native/run$RunString.raw
OutputLcioLocal=$PWD/lcio-raw/run$RunString.slcio
OutputJoboutputLocal=$PWD/log/$Name-$RunString.tar.gz
InputHotPixelLocal=$PWD/db/run$HotPixelRunNumber-hotpixel-db.slcio
LocalPWD=$PWD

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
doCommand "tar xzvf $GRIDLibraryLocal"

# from now on doing things to get access to ESA
doCommand "source ./ilc-grid-config.sh"
doCommand "$BASH ./ilc-grid-test-sys.sh || abort \"system tests failed!\" "
#doCommand ". $VO_ILC_SW_DIR/initILCSOFT.sh $GRIDILCSoftVersion"
#doCommand ". $VO_ILC_SW_DIR/ilcsoft/x86_64_gcc41_sl5/init_ilcsoft.sh v01-10"
echo "where is VO_ILC_SW_DIR: $VO_ILC_SW_DIR"
#doCommand ". $VO_ILC_SW_DIR/ilcsoft/x86_64_gcc41_sl5/init_ilcsoft.sh v01-12-01"
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

# get the input raw file
doCommand "getFromGRID ${InputRawLFN} ${InputRawLocal}"
r=$?
if [ $r -ne 0 ] ; then
    echo "Problem copying ${InputRawLFN}. Exiting with error."
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
c="ldd lib*"
echo $c
$c


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

doCommand "find ./"

echo
echo "########################################################################"
echo "# Marlin successfully finished `date `"
echo "########################################################################"
echo


# remove the input file
doCommand "rm ${InputRawLocal}"

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
echo "# Copying and registering file ${InputHotPixelLocal}" on SE
echo "########################################################################"
echo
#  
doCommand "find ${InputHotPixelLocal}"
doCommand "find ${LocalPWD}"

doCommand "putOnGRID ${InputHotPixelLocal} ${InputHotPixelLFN} ${GRIDSE} "
r=$?
if [ $r -ne 0 ] ; then
   echo "Problem copying ${InputHotPixelLocal} into ${InputHotPixelLFN}.    Warning! No hotpixel db file found at given location."
#    echo "Please, check your input to make sure you all values are set properly. If you do not want to apply hotpixel masks you can ignore this warning."
###    exit 3
fi
#

echo
echo "########################################################################"
echo "# Preparing the joboutput tarball"
echo "########################################################################"
echo
doCommand "tar czvf ${OutputJoboutputLocal} *.log *.xml out err"

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

# Job finished
echo
echo "########################################################################"
echo "# Job ($Name-$RunString) finished at `date`"
echo "########################################################################"
echo

