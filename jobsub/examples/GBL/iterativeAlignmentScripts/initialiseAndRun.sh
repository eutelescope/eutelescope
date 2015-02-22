#!/bin/bash
#Check if python exists
if [ -z `which python`  ]; then 
	echo "Python not installed"; 
	exit;	
fi;

#FUNCTIONS.
#THIS FUNCTION WILL ADD THE CORRECT NUMBER OF ZEROS TO THE END OF THE NUMBER
function adding_zeros_to_RUN {
	RUN_Before=$1
	char_length=${#RUN_Before}
	if [ $char_length -eq 1 ]
	then 
		 five_zeros="00000"
		 RUN_After="$five_zeros$RUN_Before"
	 elif [ $char_length -eq 2 ]   
	 then
		 four_zeros="0000"
		 RUN_After="$four_zeros$RUN_Before"
	 elif [ $char_length -eq 3 ]
	 then 
		 three_zeros="000"
		 RUN_After="$three_zeros$RUN_Before"
	 elif [ $char_length -eq 4 ]
	 then
		 two_zeros="00"
		 RUN_After="$two_zeros$RUN_Before"
	 elif [ $char_length -eq 5 ]
	 then
		 one_zeros="0"
		 RUN_After="$one_zeros$RUN_Before"
	 else
		 echo "There is a problem with run number"
	 fi
 echo $RUN_After
}

#VARIABLES TO INITIALISE. THERE IS NO NEED TO CHANGE THESE.
export CONFIG="$exampleLocation/config/config.cfg"
export RUNLIST="$exampleLocation/runlist/runlist.csv"
export directory="$exampleLocation/output/logs"
export pythonLocation="$scriptsLocation/pythonScripts"
export outputGearFinal="gear-${outputIdentifier}-${RUN}.xml" #This is name of the gear after all iterations of alignment. 
export histoNameInputFinal="Alignment-Runs-${outputIdentifier}-${RUN}" #This is the name of the histograms which will use the final gear to produce the tracks.
export amode="7";
export patRecMultiplicationFactor=2 #This is the factor which we increase the window of acceptance by if too few tracks found.
export PatRec=patternRecognition #This is the name of the pattern recognition steering file
export TrackFit=GBLTrackFit #This is the name of the GBLTrack fitting steering file.
export Align=GBLAlign #This is the alignment steering file name
export inputCollectionName="track_candidates"
export lcioInputName="trackcand"
export RUN=$(adding_zeros_to_RUN $RUN)
#CHECK INPUT
if [ ! -d "$directory" ]; then
	echo "The directory we will look in to find logs files does not exist. Must end now"
	exit
fi
#PRINT VARIBLES TO SCREEN
echo "These are all the input parameters to iterative alignment."
echo "Run number: $RUN"
echo "Config file: $CONFIG"
echo "Runlist file: $RUNLIST"
echo "This is the resolutions X/Y:  $xres/$yres."

#THIS WILL RUN THE ALIGNMENT PROCESS AS MANY TIME AS YOU LIKE TO IMPROVE ALIGNMENT
$scriptsLocation/howManyIterationsDecider.sh -n "$numberOfIterations"
#$scriptsLocation/patRecAndTrackFit.sh -i "$outputGearFinal" -h "$histoNameInputFinal"  
