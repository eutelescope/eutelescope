#!/bin/bash
#This file does a quick pat->track using the gear specified in config.
while getopts r:h: option
do
        case "${option}"
        in
								r) RUN=${OPTARG};;
								h) histoNameInput=${OPTARG};;
       esac
done

echo "These are the input parameters to patRec and TrackFit."
echo "Run number $RUN"
echo "HistoNameInput $histoNameInput"

PatRec="patternRecognition"
DRY="--dry-run"
#THIS PART WILL RUN PATTERN RECOGNTION AND TRACKFITTING AGAIN WITH THE FINAL GEAR FILE.
$do jobsub.py -c config/config.cfg -csv runlist/runlist.csv $PatRec $RUN 
TrackFit="GBLTrackFit"
#TrackFit="GBLTrackFit-Multi"
#Then we create the first tracks using pattern recognition tracks.
$do jobsub.py  -c config/config.cfg -csv runlist/runlist.csv -o histoName="$histoNameInput" $TrackFit  $RUN
