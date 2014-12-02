#We have a bash script for each example to store the setting to align. Note these will be different from what is stored in configuration file. Since the configuration file will store the settings to fit tracks. 
#!/bin/bash
while getopts n: option
do
        case "${option}"
        in
								n) numberOfIterations=${OPTARG};;
       esac
done

if [  -z "$numberOfIterations" ]; then
	echo "The number of iterations is empty in iterativeAlignment bash script!"
	exit
fi

#VARIABLES. IMPORTANT CONSTANT VARIABLES THROUGH THE WHOLE ALIGNMENT PROCESS ARE SET HERE. IMPORTANT NOT ALL VARIABLES ARE HERE LOOK IN STEERING FILES/CONFIG
export RUN="613" 
export numberOfIterations="$numberOfIterations"
export exampleLocation="/afs/phas.gla.ac.uk/user/a/amorton/ilcsoft/v01-17-05/Eutelescope/trunk/jobsub/examples/GBL/stripSensorExample"
export scriptsLocation="/afs/phas.gla.ac.uk/user/a/amorton/ilcsoft/v01-17-05/Eutelescope/trunk/GBL/iterativeAlignmentScripts"
export planeDimensions="2 2 2 1 1 2 2 2" #This specified if the planes is pixel (2) or strip (1)
export MaxMissingHitsPerTrack="0" #This is the number of missing hits that a track can have on the planes
export AllowedSharedHitsOnTrackCandidate="0" #This is the number of hits two tracks can have in common. One track is discarded if over the limit.
export minTracksPerEventAcceptance=0.5 #This is the number of tracks that is needed per event before we stop pattern recognition. Note value should depend on other cuts. 
export dutPlanes="6 7" #Since the mimosa planes are always named the same but dut can be different, must note dut. This does NOT indicate that it will be used in analysis.
export ExcludePlanes="6" #These planes are completely excluded from the analysis. The scattering from the plane however is still taken into account.
export ResidualsRMax="1" #This is the window size on the next plane that we will accept hits from. This will increase if less than 1 track per event is found.
export Verbosity="MESSAGE5"
export r="1"; #Make resolution large so we begin with small chi2 and factor improvement to get to chi2/ndf=1 on next iteration. 
export dutX=0.5 #Need to add duts resolution like this since in alignment we times this by some factor during the process.
export dutY=1000000000000
export dutXs="$dutX $dutX" #This is the resolution of the DUT in the x LOCAL direction taking into account the misalignment
export dutYs="$dutY $dutY" #This is the resolution of the DUT in the y LOCAL direction taking into account the misalignment
export allMimosaPlanesFixed=" " #Planes 0 and 5 are always fixed. These are the additional planes you want to fix.
export MaxRecordNumber="30000" 
export inputGearInitial="gear-stripSensor-Aligned-With-noDUT.xml"
#export inputGearInitial="gear-final-noDUT-${RUN}.xml"
export outputIdentifier="final-DUT6-30k-Events-R=1" #Use this string to identify final gear/histogram and all iterations before.


$scriptsLocation/initialiseAndRun.sh
