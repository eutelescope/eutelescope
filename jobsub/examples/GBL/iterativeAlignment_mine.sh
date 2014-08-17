#!/bin/bash

# parameters:
# 1 - runnumber         $RUN  
# 2 - runlist csv file: $RUNLIST

while getopts c:l:u:r: option
do
        case "${option}"
        in
                c) CONFIG=${OPTARG};;
                l) RUNLIST=${OPTARG};;
                u) Chi2Cut=${OPTARG};;
                r) RUN=${OPTARG};;
       esac
done

MaxRecordNumber="1000" 

echo "Input recieved"
echo "Config: " $CONFIG
echo "Runlist: " $RUNLIST
echo "Chi2Cut: " $Chi2Cut
echo "Run (MAKE SURE YOU HAVE THE CORRECT NUMBER OF ZEROS!): " $RUN #TO DO
echo "MaxRecordNumber: " $MaxRecordNumber

ExcludePlanes=""

#This is the alignment mode. It sets the size of the alignment jacobian dimensions.
amode="7";

pede="chiscut 50. 5. "

r="0.25000";
xres="$r $r $r $r $r $r";
yres="$r $r $r $r $r $r";

prev="$r";
echo "prev:$prev and r:$r"
fileTrack="/scratch/ilcsoft/v01-17-05/Eutelescope/master/jobsub/examples/GBL/output/logs/GBLTrackFit-0000${RUN}.zip"

if [ $# -ne 8 ]
then
 echo "$# parameters: $CONFIG $RUNLIST $Chi2Cut $RUN $file $gear10"
 return
fi
Fxr="0 1 2 3 4 5"
Fxs="0         "
Fyr="0 1 2 3 4 5"
Fys="0         "
Fzr="0 1 2 3 4 5"
Fzs="0 1 2 3 4 5"


$do jobsub.py  $DRY -c $CONFIG -csv $RUNLIST -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" GBLTrackFit  $RUN $DR
AverageChi2=`unzip  -p  $fileTrack |grep "This is the average chi2 -" |cut -d '-' -f2`; 
echo "The average chi2ndf is :  $AverageChi2"

