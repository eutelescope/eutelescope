#!/bin/bash
#VARIABLES. EVERYTHING THAT IS SET FOR ALL ALIGNMENT STEPS IS SET HERE.IMPORTANT NOT ALL VARIABLES ARE HERE LOOK IN STEERING FILES
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
export ResidualsRMax="0.5" #This is the window size on the next plane that we will accept hits from.
export minTracksPerEventAcceptance=1 #This is the number of tracks that is needed per event before we stop pattern recognition 
export Fxr="0 1 2 3 4 5" #This is the fixed planes for rotations round the x axis
export Fxs="0         5" #This is the fixed planes for shifts in the x axis
export Fyr="0 1 2 3 4 5" #This is the fixed planes for rotations round the y axis
export Fys="0         5" #This is the fixed planes for shifts in the y axis
export Fzr="0 1 2 3 4 5" #This is the fixed planes for rotations round the z axis
export Fzs="0 1 2 3 4 5" #This is the fixed planes for shifts in the z axis
export ExcludePlanes="" #These planes are completely excluded from the analysis. The scattering from the plane however is still taken into account.
export minChi2AlignAcceptance=1

export Verbosity="MESSAGE5"
export inputGear="gear-1T.xml"
#export inputGear="gear_desy2012_150mm.xml"
export outputGear="gear-final-XYshiftTest1-${RUN}.xml"
export histoNameInput="GBLtrack-XYshiftTest1-${RUN}"

export r="1"; #This is the resolution of the mimosa sensors. This can be changed to taken into account of the misalignment but in most cases leave alone.
export dutX="" #This is the resolution of the DUT in the x LOCAL direction. This should just a rough estimate of the expect resolution.
export dutY="" #This is the resolution of the DUT in the y LOCAL direction. This should just a rough estimate of the expect resolution.
export MaxRecordNumber="5000" 

#VARIABLES TO INITIALISE. THERE IS NO NEED TO CHANGE THESE.
export xres="$r $r $r $dutX $r $r $r";
export yres="$r $r $r $dutY $r $r $r";
export amode="7";
export patRecMultiplicationFactor=2 #This is the factor which we increase the window of acceptance by if too few tracks.

#FUNCTIONS. 
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
#This is the loop to produce the new gear file for x/y shifts. 
./singleIterationAlignment.sh

export Fxr="0 1 2 3 4 5" #This is the fixed planes for rotations round the x axis
export Fxs="0         5" #This is the fixed planes for shifts in the x axis
export Fyr="0 1 2 3 4 5" #This is the fixed planes for rotations round the y axis
export Fys="0         5" #This is the fixed planes for shifts in the y axis
export Fzr="0" #This is the fixed planes for rotations round the z axis
export Fzs="0 1 2 3 4 5" #This is the fixed planes for shifts in the z axis
export inputGear="$outputGear"
export outputGear="gear-final-ZRotations-${RUN}.xml"
export histoNameInput="GBLtrack-ZRotations-${RUN}"
./singleIterationAlignment.sh

$do jobsub.py -c $CONFIG -csv $RUNLIST -o HitInputCollectionName="$lcioPatternCollection" -o Verbosity="$Verbosity" -o MaxRecordNumber="$MaxRecordNumber" -o GearFile="$outputGear"  -o ExcludePlanes="$ExcludePlanes" $PatRec $RUN 
r="0.005";
xres="$r $r $r $dut $r $r $r";
yres="$r $r $r $dut $r $r $r";

#Then we create the first tracks using pattern recognition tracks.
$do jobsub.py  -c $CONFIG -csv $RUNLIST -o Verbosity="$Verbosity" -o GearFile="$outputGear" -o histoName="$histoNameInput" -o lcioInputName="trackcand"  -o inputCollectionName="track_candidates" -o lcioOutputName="GBLtracks" -o outputCollectionName="tracks"  -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" $TrackFit  $RUN
