#!/bin/sh

# A template of converter job
#
# @author Antonio Bulgheroni <mailto:antonio.bulgheroni@gmail.com>
# @version $Id: runjob-converter-tmp.sh,v 1.1 2009-05-11 17:16:17 bulgheroni Exp $
#

# To be replaced with the runString in the format %(run)06d
RunString="@RunString@"

# Define here all the variables modified by the submitter
GRIDCE="@GRIDCE@"
GRIDSE="@GRIDSE@"
GRIDStoreProtocol="@GRIDStoreProtocol@"
GRIDVO="@GRIDVO@"
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
mkdir native
mkdir lcio-raw
mkdir results
mkdir histo
mkdir pics
mkdir db
mkdir log


# unpack the library
echo "> tar xzvf $GRIDLibraryTarball"
tar xzvf $GRIDLibraryTarball

# rename the simjob.slcio because otherwise it gets delete
echo "> mv simjob.slcio simjob.slcio.keepme "
mv simjob.slcio simjob.slcio.keepme
echo " "

# from now on doing things to get access to ESA
echo "> source ./ilc-grid-config.sh"
source ./ilc-grid-config.sh
echo " " 

echo "> $BASH ./ilc-grid-test-sys.sh || abort \"system tests failed!\""
$BASH ./ilc-grid-test-sys.sh || abort "system tests failed!"
echo " " 

echo "> . new $VO_ILC_SW_DIR/initILCSOFT.sh $GRIDILCSoftVersion "
. $VO_ILC_SW_DIR/initILCSOFT.sh $GRIDILCSoftVersion
echo " "


echo "> $BASH ./ilc-grid-test-sw.sh ; r=$?"
$BASH ./ilc-grid-test-sw.sh ; r=$?
if [ $r -ne 0 ] ; then 
    echo "ERROR: " ; cat ./ilc-grid-test-sw.log
    abort "software tests failed!"
fi
echo " "


# ILCSOFT ready
echo " "
echo " ILCSOFT ready to use "
echo " ------------------------------------------"
echo " "

# now it's safe to rename the simjob to the original
echo "> mv simjob.slcio.keepme simjob.slcio"
mv simjob.slcio.keepme simjob.slcio
echo " "

# set the list of Marlin plugins
c="export MARLIN_DLL=$PWD/libEutelescope.so.0.0.8"
echo "> $c"
$c
echo " "

echo "> export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH
echo " "

# setting permissions for scripts.
echo "> chmod a+x ./get_file_by_lfn.sh"
chmod  a+x ./get_file_by_lfn.sh
echo " "

echo "> chmod a+x ./put_file.sh"
chmod  a+x ./put_file.sh
echo " "

return 0

echo "> ./get_file_by_lfn.sh ${InputRawLFN} ${InputRawLocal}"
./get_file_by_lfn.sh ${InputRawLFN} ${InputRawLocal}"
rc=$?
if [ $[rc!=0] = "1" ]; then
  echo "##### Copying '${InputRawLFN}' failed. Exiting with error."
  exit 0
else
  echo "+++++ Copying '${InputRawLFN}' succeeded. Continuing."
fi

# list all the files available
echo "> ls -l"
ls -al
echo " "

echo "> ls -l"
ls -al
echo " "

# ready to run marlin
echo "Starting Marlin `date`"
echo " "

# c="Marlin $SteeringFile | tee $LogFile "
c="Marlin $SteeringFile | tee $LogFile"
echo "> $c"
$c
echo " "

# remove the input file
echo "> rm ${InputRawLocal}"
rm ${InputRawLocal}

echo "Preparing the joboutput tarball"
tar czvf ${OutputJoboutputLocal} *.log *.xml 

echo "Copying and registering the tarball to SE"
./put_file.sh ${OutputJoboutputLNF} ${OutputJoboutputLocal} ${GRIDSE}
if [ $[rc!=0] = "1" ]; then
    echo "##### Copying '$LFN' failed. No output tar-ball written."
    exit 0
else
    echo "+++++ Copying '$LFN' succeeded. Output tar-ball written."
fi

echo "Copying and registering the output file to SE"
./put_file.sh ${OutputLcioLFN} ${OutputLcioLocal} ${GRIDSE}
if [ $[rc!=0] = "1" ]; then
    echo "##### Copying '$LFN' failed."
    exit 0
else
    echo "+++++ Copying '$LFN' succeeded."
fi


# Job finished
echo "Job finished at `date`"
