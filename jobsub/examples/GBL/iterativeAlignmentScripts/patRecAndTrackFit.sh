#!/bin/bash
while getopts i:h: option
do
        case "${option}"
        in
								i) inputGear=${OPTARG};;
								h) histoNameInput=${OPTARG};;
       esac
done

echo "These are the input parameters to patRec and TrackFit."
echo "Input Gear: $inputGear"
echo "Iteration: $histoNameInput"

if [  -z "$inputGear" ]; then
	echo "The input gear is empty to patRec and TrackFit!"
	exit
fi
if [  -z "$histoNameInput" ]; then
	echo "The histogram input string to patRec and TrackFit!"
	exit
fi

#THIS PART WILL RUN PATTERN RECOGNTION AND TRACKFITTING AGAIN WITH THE FINAL GEAR FILE.
$do jobsub.py -c $CONFIG -csv $RUNLIST  -o planeDimensions="${planeDimensions}" -o Verbosity="$Verbosity" -o MaxRecordNumber="$MaxRecordNumber" -o GearFile="$inputGear"  -o ExcludePlanes="$ExcludePlanes" $PatRec $RUN 
#TO DO: This be set here might be confusing should move to iterativeAlignment
r="0.005"; #Correct mimosa resolution
export dutX="0.250" #Correct DUT resolution.
export dutY="0.05"
export dutXs="$dutX $dutX" 
export dutYs="$dutY $dutY" 
xres="$r $r $r $dutXs $r $r $r";
yres="$r $r $r $dutYs $r $r $r";

#Then we create the first tracks using pattern recognition tracks.
$do jobsub.py  -c $CONFIG -csv $RUNLIST -o Verbosity="$Verbosity" -o GearFile="$inputGear" -o histoName="$histoNameInput" -o lcioInputName="trackcand"  -o inputCollectionName="track_candidates" -o lcioOutputName="GBLtracks" -o outputCollectionName="tracks"  -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" $TrackFit  $RUN
