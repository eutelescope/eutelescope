#!/bin/sh
# A template of pedestal job
#
# @author Antonio Bulgheroni <mailto:antonio.bulgheroni@gmail.com>
# @version $Id: runjob-pedestal-tmp.sh,v 1.10 2009-07-28 00:13:59 bulgheroni Exp $
#
# errno  0: No error.
# errno  1: Unable to get the GRID library tarball from the SE
# errno  2: Unable to get the input file from the SE.
# errno 20: Problem during Marlin execution (telescope only)
# errno 21: Problem during Marlin execution (telescope part)
# errno 22: Problem during Marlin execution (dut part)
# errno 23: Problem during Marlin execution (dut only)
# errno 24: Problem running pedestalmerge
# errno 30: Problem copying and registering the DB (full file) output to the SE.
# errno 31: Problem copying and registering the DB (telescope file) output to the SE.
# errno 32: Problem copying and registering the DB (DUT file) output to the SE.
# errno 33: Problem copying and registering the Joboutput to the SE.
# errno 34: Problem copying and registering the histogram (full file) to the SE
# errno 35: Problem copying and registering the histogram (full file) to the SE
# errno 36: Problem copying and registering the histogram (full file) to the SE
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
# errno 20: Problem during Marlin execution (telescope only)
# errno 21: Problem during Marlin execution (telescope part)
# errno 22: Problem during Marlin execution (dut part)
# errno 23: Problem during Marlin execution (dut only)
# errno 24: Problem running pedestalmerge
# errno 30: Problem copying and registering the DB (full file) output to the SE.
# errno 31: Problem copying and registering the DB (telescope file) output to the SE.
# errno 32: Problem copying and registering the DB (DUT file) output to the SE.
# errno 33: Problem copying and registering the Joboutput to the SE.
# errno 34: Problem copying and registering the histogram (full file) to the SE
# errno 35: Problem copying and registering the histogram (full file) to the SE
# errno 36: Problem copying and registering the histogram (full file) to the SE


# To be replaced with the runString in the format %(run)06d
RunString="@RunString@"

# To be replace with the job name used for the identification of
# all files. It should be something like converter
Name="@Name@"

# Depending on the presence or not of the DUT, set the following
# variables
IsTelescopeOnly="@IsTelescopeOnly@"
IsTelescopeAndDUT="@IsTelescopeAndDUT@"
IsDUTOnly="@IsDUTOnly@"

# the DUT Suffix
DUTSuffix="@DUTSuffix@"

# Define here all the variables modified by the submitter
GRIDCE="@GRIDCE@"
GRIDSE="@GRIDSE@"
GRIDStoreProtocol="@GRIDStoreProtocol@"
GRIDVO="@GRIDVO@"
LFC_HOME="/grid/ilc/eudet-jra1"
GRIDFolderBase="@GRIDFolderBase@"
GRIDFolderLcioRaw="@GRIDFolderLcioRaw@"
GRIDFolderDBPede="@GRIDFolderDBPede@"
GRIDFolderPedestalJoboutput="@GRIDFolderPedestalJoboutput@"
GRIDFolderPedestalHisto="@GRIDFolderPedestalHisto@"
GRIDILCSoftVersion="@GRIDILCSoftVersion@"

# GRID Tarball
# LocalGRIDLibraryTarball --> "yes" means that the tarball is uploaded along with the JDL file
#                         --> "no" means that it has to be downloaded from a SE
HasLocalGRIDLibraryTarball="@HasLocalGRIDLibraryTarball@"
GRIDLibraryTarball="@GRIDLibraryTarball@"
GRIDLibraryTarballPath="@GRIDLibraryTarballPath@"
GRIDLibraryLocal=$PWD/$GRIDLibraryTarball
GRIDLibraryLFN=$GRIDLibraryTarballPath/$GRIDLibraryTarball

# Input file local and lfn
InputLcioRawLocal=$PWD/lcio-raw/run$RunString.slcio
InputLcioRawLFN=$GRIDFolderLcioRaw/run$RunString.slcio

# output db local and lfn
OutputDBLocal=""
OutputDB_TEL_Local=""
OutputDB_DUT_Local=""
OutputDBLFN=""
OutputDB_TEL_LFN=""
OutputDB_DUT_LFN=""

# output histo file local and lfn
OutputHistoLocal=""
OutputHisto_TEL_Local=""
OutputHisto_DUT_Local=""
OutputHistoLFN=""
OutputHisto_TEL_LFN=""
OutputHisto_DUT_LFN=""

# steering file
SteeringFile=$Name-$RunString.xml
SteeringFile_TEL=$Name-$RunString-telescope.xml
SteeringFile_DUT=$Name-$RunString-$DUTSuffix.xml

if [ $IsTelescopeOnly == "yes" ] ; then
    OutputDBLocal=$PWD/db/run$RunString-ped-db.slcio
    OutputDBLFN=$GRIDFolderDBPede/run$RunString-ped-db.slcio
#
    OutputHistoLocal=$PWD/histo/run$RunString-ped-histo.root
    OutputHistoLFN=$GRIDFolderPedestalHisto/run$RunString-ped-histo.root
fi

if [ $IsTelescopeAndDUT == "yes" ] ; then
    OutputDBLocal=$PWD/db/run$RunString-ped-db.slcio
    OutputDB_TEL_Local=$PWD/db/run$RunString-ped-telescope-db.slcio
    OutputDB_DUT_Local=$PWD/db/run$RunString-ped-$DUTSuffix-db.slcio
    OutputDBLFN=$GRIDFolderDBPede/run$RunString-ped-db.slcio
    OutputDB_TEL_LFN=$GRIDFolderDBPede/run$RunString-ped-telescope-db.slcio
    OutputDB_DUT_LFN=$GRIDFolderDBPede/run$RunString-ped-$DUTSuffix-db.slcio
#
    OutputHistoLocal=$PWD/histo/run$RunString-ped-histo.root
    OutputHisto_TEL_Local=$PWD/histo/run$RunString-ped-telescope-histo.root
    OutputHisto_DUT_Local=$PWD/histo/run$RunString-ped-$DUTSuffix-histo.root
    OutputHistoLFN=$GRIDFolderPedestalHisto/run$RunString-ped-histo.root
    OutputHisto_TEL_LFN=$GRIDFolderPedestalHisto/run$RunString-ped-telescope-histo.root
    OutputHisto_DUT_LFN=$GRIDFolderPedestalHisto/run$RunString-ped-$DUTSuffix-histo.root
fi

if [ $IsDUTOnly == "yes" ] ; then
    OutputDB_DUT_Local=$PWD/db/run$RunString-ped-$DUTSuffix-db.slcio
    OutputDB_DUT_LFN=$GRIDFolderDBPede/run$RunString-ped-$DUTSuffix-db.slcio
#
    OutputHisto_DUT_Local=$PWD/histo/run$RunString-ped-$DUTSuffix-histo.root
    OutputHisto_DUT_LFN=$GRIDFolderPedestalHisto/run$RunString-ped-$DUTSuffix-histo.root
fi

# joboutput file local and lfn
OutputJoboutputLocal=$PWD/log/$Name-$RunString.tar.gz
OutputJoboutputLFN=$GRIDFolderPedestalJoboutput/$Name-$RunString.tar.gz

# log file name
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
doCommand "export PATH=$PWD:$PATH"

# get the input raw file
doCommand "getFromGRID ${InputLcioRawLFN} ${InputLcioRawLocal}"
r=$?
if [ $r -ne 0 ] ; then
    echo "Problem copying ${InputLcioRawLFN}. Exiting with error."
    exit 0
fi

# list all the files available
doCommand "ls -al"

# ready to run marlin
echo
echo "########################################################################"
echo "# Starting Marlin `date`"
echo "########################################################################"
echo

if [ $IsTelescopeOnly == "yes" ] ; then
    c="Marlin $SteeringFile"
    echo $c
    $c
    r=$?


    if [ $r -ne 0 ] ; then
        echo "****** Problem running Marlin (telescope only)"
        exit 20
    fi
fi

if [ $IsTelescopeAndDUT == "yes" ] ; then 
    c="Marlin $SteeringFile_TEL"
    echo $c
    $c
    r=$?

    if [ $r -ne 0 ] ; then
        echo "****** Problem running Marlin for the telescope part"
        exit 21
    fi

    c="Marlin $SteeringFile_DUT"
    echo $c
    $c
    r=$?

    if [ $r -ne 0 ] ; then
        echo "****** Problem running Marlin for the DUT part"
        exit 22
    fi

    # merge the db
    doCommand "pedestalmerge -v -o ${OutputDBLocal} ${OutputDB_TEL_Local} ${OutputDB_DUT_Local}"
    r=$?
    if [ $r -ne 0 ] ; then
        echo "****** Problem running the pedestal merge"
        exit 24
    fi
fi

if [ $IsDUTOnly == "yes" ] ; then
    c="Marlin $SteeringFile_DUT"
    echo $c
    $c
    r=$?

    if [ $r -ne 0 ] ; then
        echo "****** Problem running Marlin for the DUT part"
        exit 23
    fi
fi

echo
echo "########################################################################"
echo "# Marlin successfully finished `date `"
echo "########################################################################"
echo


# remove the input file
doCommand "rm ${InputLcioRawLocal}"

# fixing the histograms
if [ $IsTelescopeOnly == "yes" ] ; then
    doCommand "hadd -f temp.root empty.root ${OutputHistoLocal}"
    doCommand "mv temp.root ${OutputHistoLocal}"
fi

if [ $IsTelescopeAndDUT == "yes" ] ; then
    doCommand "hadd -f ${OutputHistoLocal} ${OutputHisto_TEL_Local} ${OutputHisto_DUT_Local}"
fi

if [ $IsDUTOnly == "yes" ] ; then
    doCommand "hadd -f temp.root empty.root ${OutputHisto_DUT_Local}"
    doCommand "mv temp.root ${OutputHisto_DUT_Local}"
fi

# put back the files to the GRID
echo
echo "########################################################################"
echo "# Copying and registering the output file to SE"
echo "########################################################################"
echo

if [ $IsTelescopeOnly == "yes" ] ; then
    doCommand "putOnGRID ${OutputDBLocal} ${OutputDBLFN} ${GRIDSE}"
    r=$?
    if [ $r -ne 0 ] ; then
        echo "****** Problem copying the ${OutputDBLocal} to the GRID"
        exit 30
    fi
fi

if [ $IsTelescopeAndDUT == "yes" ] ; then
    doCommand "putOnGRID ${OutputDBLocal} ${OutputDBLFN} ${GRIDSE}"
    r=$?
    if [ $r -ne 0 ] ; then
        echo "****** Problem copying the ${OutputDBLocal} to the GRID"
        exit 30
    fi

    doCommand "putOnGRID ${OutputDB_TEL_Local} ${OutputDB_TEL_LFN} ${GRIDSE}"
    r=$?
    if [ $r -ne 0 ] ; then
        echo "****** Problem copying the ${OutputDB_TEL_Local} to the GRID"
        exit 31
    fi

    doCommand "putOnGRID ${OutputDB_DUT_Local} ${OutputDB_DUT_LFN} ${GRIDSE}"
    r=$?
    if [ $r -ne 0 ] ; then
        echo "****** Problem copying the ${OutputDB_DUT_Local} to the GRID"
        exit 32
    fi
fi

if [ $IsDUTOnly == "yes" ] ; then
    doCommand "putOnGRID ${OutputDB_DUT_Local} ${OutputDB_DUT_LFN} ${GRIDSE}"
    r=$?
    if [ $r -ne 0 ] ; then
        echo "****** Problem copying the ${OutputDB_DUT_Local} to the GRID"
        exit 32
    fi
fi


echo
echo "########################################################################"
echo "# Preparing the joboutput tarball"
echo "########################################################################"
echo
doCommand "tar czvf ${OutputJoboutputLocal} *.log *.xml out err db/*dat histo/*root"

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

if [ $IsTelescopeOnly == "yes" ] ; then
    doCommand "putOnGRID ${OutputHistoLocal} ${OutputHistoLFN} ${GRIDSE}"
    r=$?
    if [ $r -ne 0 ] ; then
        echo "****** Problem copying the ${OutputHistoLocal} to the GRID"
        exit 34
    fi
fi

if [ $IsTelescopeAndDUT == "yes" ] ; then
    doCommand "putOnGRID ${OutputHistoLocal} ${OutputHistoLFN} ${GRIDSE}"
    r=$?
    if [ $r -ne 0 ] ; then
        echo "****** Problem copying the ${OutputHistoLocal} to the GRID"
        exit 34
    fi

    doCommand "putOnGRID ${OutputHisto_TEL_Local} ${OutputHisto_TEL_LFN} ${GRIDSE}"
    r=$?
    if [ $r -ne 0 ] ; then
        echo "****** Problem copying the ${OutputHisto_TEL_Local} to the GRID"
        exit 35
    fi

    doCommand "putOnGRID ${OutputHisto_DUT_Local} ${OutputHisto_DUT_LFN} ${GRIDSE}"
    r=$?
    if [ $r -ne 0 ] ; then
        echo "****** Problem copying the ${OutputHisto_DUT_Local} to the GRID"
        exit 36
    fi
fi

if [ $IsDUTOnly == "yes" ] ; then
    doCommand "putOnGRID ${OutputHisto_DUT_Local} ${OutputHisto_DUT_LFN} ${GRIDSE}"
    r=$?
    if [ $r -ne 0 ] ; then
        echo "****** Problem copying the ${OutputHisto_DUT_Local} to the GRID"
        exit 36
    fi
fi

# Job finished
echo
echo "########################################################################"
echo "# Job ($Name-$RunString) finished at `date`"
echo "########################################################################"
echo

