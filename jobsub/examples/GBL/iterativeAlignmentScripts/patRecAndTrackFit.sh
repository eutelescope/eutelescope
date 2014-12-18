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
$do jobsub.py -c $CONFIG -csv $RUNLIST  -o GearFile="$inputGear" $PatRec $RUN 

#Then we create the first tracks using pattern recognition tracks.
$do jobsub.py  -c $CONFIG -csv $RUNLIST -o GearFile="$inputGear" -o histoName="$histoNameInput" $TrackFit  $RUN
