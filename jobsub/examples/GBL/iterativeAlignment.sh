#!/bin/bash
#VARIABLES. IMPORTANT CONSTANT VARIABLES THROUGH THE WHOLE ALIGNMENT PROCESS IS SET HERE. IMPORTANT NOT ALL VARIABLES ARE HERE LOOK IN STEERING FILES/CONFIG
RUN="290" #This is the run number. Zeros are added later and then export
export CONFIG="/afs/phas.gla.ac.uk/user/a/amorton/ilcsoft/v01-17-05/Eutelescope/trunk/jobsub/examples/GBL/config/config.cfg"
export RUNLIST="/afs/phas.gla.ac.uk/user/a/amorton/ilcsoft/v01-17-05/Eutelescope/trunk/jobsub/examples/GBL/runlist/runlist.csv"
export directory="/afs/phas.gla.ac.uk/user/a/amorton/ilcsoft/v01-17-05/Eutelescope/trunk/jobsub/examples/GBL/output/logs"
export PatRec=patternRecognition #This is the name of the pattern recognition steering file
export TrackFit=GBLTrackFit #This is the name of the GBLTrack fitting steering file.
export Align=GBLAlign #This is the alignment steering file name
export lcioPatternCollection="hit_filtered_m26" #This is the lcio name given to the hit collection. We may need to change this if we are give only lcio with hits.
export planeDimensions="2 2 2 2 2 2"
export MaxMissingHitsPerTrack="0" #This is the number of missing hits that a track can have on the planes
export AllowedSharedHitsOnTrackCandidate="0" #This is the number of hits two tracks can have in common. One is discarded if over.
export minTracksPerEventAcceptance=1 #This is the number of tracks that is needed per event before we stop pattern recognition 
export ExcludePlanes="" #These planes are completely excluded from the analysis. The scattering from the plane however is still taken into account.
export minChi2AlignAcceptance=1
export ResidualsRMax="0.5" #This is the window size on the next plane that we will accept hits from.
export Verbosity="MESSAGE5"
export r="1"; #This is the resolution of the mimosa sensors. This can be changed to taken into account of the misalignment but in most cases leave alone.
export dutX="" #This is the resolution of the DUT in the x LOCAL direction. This should just a rough estimate of the expect resolution.
export dutY="" #This is the resolution of the DUT in the y LOCAL direction. This should just a rough estimate of the expect resolution.
export MaxRecordNumber="5000" 
export inputGearInitial="gear-1T.xml" #This is the initial input gear. 
#export inputGearInitial="gear_desy2012_150mm.xml"
export outputGearFinal="gear-final-${RUN}.xml" #This is name of the gear after all iterations of alignment. 
export histoNameInputFinal="GBLtrack-final-${RUN}" #This is the name of the histograms which will will use the final gear to produce the tracks.

#VARIABLES TO INITIALISE. THERE IS NO NEED TO CHANGE THESE.
export xres="$r $r $r $dutX $r $r $r";
export yres="$r $r $r $dutY $r $r $r";
export amode="7";
export patRecMultiplicationFactor=2 #This is the factor which we increase the window of acceptance by if too few tracks.


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
./iterativeAlignment2.sh -i $inputGearInitial -o "gear-finished-iteration-1" -n 1
./iterativeAlignment2.sh -i "gear-finished-iteration-1" -o "$outputGearFinal" -n 2

#THIS PART WILL RUN PATTERN RECOGNTION AND TRACKFITTING AGAIN WITH THE FINAL GEAR FILE.
$do jobsub.py -c $CONFIG -csv $RUNLIST -o HitInputCollectionName="$lcioPatternCollection" -o Verbosity="$Verbosity" -o MaxRecordNumber="$MaxRecordNumber" -o GearFile="$outputGearFinal"  -o ExcludePlanes="$ExcludePlanes" $PatRec $RUN 
r="0.005";
xres="$r $r $r $dut $r $r $r";
yres="$r $r $r $dut $r $r $r";

#Then we create the first tracks using pattern recognition tracks.
$do jobsub.py  -c $CONFIG -csv $RUNLIST -o Verbosity="$Verbosity" -o GearFile="$outputGearFinal" -o histoName="$histoNameInputFinal" -o lcioInputName="trackcand"  -o inputCollectionName="track_candidates" -o lcioOutputName="GBLtracks" -o outputCollectionName="tracks"  -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" $TrackFit  $RUN
