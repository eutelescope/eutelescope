#!/bin/bash
#VARIABLES. IMPORTANT CONSTANT VARIABLES THROUGH THE WHOLE ALIGNMENT PROCESS ARE SET HERE. IMPORTANT NOT ALL VARIABLES ARE HERE LOOK IN STEERING FILES/CONFIG
RUN="105" #This is the run number. Zeros are added later and then exported
export CONFIG="/afs/phas.gla.ac.uk/user/a/amorton/ilcsoft/v01-17-05/Eutelescope/trunk/jobsub/examples/GBL/config/config.cfg"
export RUNLIST="/afs/phas.gla.ac.uk/user/a/amorton/ilcsoft/v01-17-05/Eutelescope/trunk/jobsub/examples/GBL/runlist/runlist.csv"
export directory="/afs/phas.gla.ac.uk/user/a/amorton/ilcsoft/v01-17-05/Eutelescope/trunk/jobsub/examples/GBL/output/logs"
export PatRec=patternRecognition #This is the name of the pattern recognition steering file
export TrackFit=GBLTrackFit #This is the name of the GBLTrack fitting steering file.
export Align=GBLAlign #This is the alignment steering file name
export lcioPatternCollection="hit_filtered_m26" #This is the lcio name given to the hit collection. We may need to change this if we are give only lcio file with hits.
export planeDimensions="2 2 2 2 2 2 2 2" #This specified if the planes is pixel (2) or strip (1)
export MaxMissingHitsPerTrack="0" #This is the number of missing hits that a track can have on the planes
export AllowedSharedHitsOnTrackCandidate="0" #This is the number of hits two tracks can have in common. One track is discarded if over the limit.
export minTracksPerEventAcceptance=1 #This is the number of tracks that is needed per event before we stop pattern recognition 
export ExcludePlanes="20 21" #These planes are completely excluded from the analysis. The scattering from the plane however is still taken into account.
export ResidualsRMax="4" #This is the window size on the next plane that we will accept hits from. This will increase if less than 1 track per event is found.
export Verbosity="MESSAGE5"
export r="1"; #This is the resolution of the mimosa sensors taking into account the systematic uncertainty of the alignment. After alignment should be 0.005 mm
export dutX="1 1" #This is the resolution of the DUT in the x LOCAL direction taking into account the misalignment
export dutY="1 1" #This is the resolution of the DUT in the x LOCAL direction taking into account the misalignment
export MaxRecordNumber="10000" 
#export inputGearInitial="gear-1T.xml" #This is the initial input gear. 
#export inputGearInitial="gear_desy2012_150mm.xml"
#export inputGearInitial="gear-stripSensor-noDUT.xml"
export inputGearInitial="gear-quad.xml"
export outputGearFinal="gear-final-${RUN}.xml" #This is name of the gear after all iterations of alignment. 
export histoNameInputFinal="GBLtrack-final-${RUN}" #This is the name of the histograms which will use the final gear to produce the tracks.

#VARIABLES TO INITIALISE. THERE IS NO NEED TO CHANGE THESE.
export xres="$r $r $r $dutX $r $r $r";
export yres="$r $r $r $dutY $r $r $r";
export amode="7";
export patRecMultiplicationFactor=2 #This is the factor which we increase the window of acceptance by if too few tracks found.

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
./patRecToAlignmentMultiLoop.sh -i "$inputGearInitial" -o "gear-finished-iteration-1.xml" -n 1
./patRecToAlignmentMultiLoop.sh -i "gear-finished-iteration-1.xml" -o "gear-finished-iteration-2.xml" -n 2
./patRecToAlignmentMultiLoop.sh -i "gear-finished-iteration-2.xml" -o "$outputGearFinal" -n 3


#./patRecToAlignmentMultiLoop.sh -i "$inputGearInitial" -o "gear-finished-iteration-1.xml" -n 1
#./patRecToAlignmentMultiLoop.sh -i "gear-finished-iteration-1.xml" -o "$outputGearFinal" -n 2

#./patRecToAlignmentMultiLoop.sh -i "$inputGearInitial" -o "$outputGearFinal" -n 1

./patRecAndTrackFit.sh -i "$outputGearFinal" -h "$histoNameInputFinal"  
