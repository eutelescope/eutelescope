#We have a bash script for each example to store the setting to align. Note these will be different from what is stored in configuration file. Since the configuration file will store the settings to fit tracks. 
#!/bin/bash
while getopts n:i:r: option
do
        case "${option}"
        in
								n) numberOfIterations=${OPTARG};;
								i) identifier=${OPTARG};;
								r) RUN=${OPTARG};;
       esac
done

if [  -z "$numberOfIterations" ]; then
	echo "The number of iterations is empty in iterativeAlignment bash script!"
	exit
fi
#The location of the example and the scripts
export exampleLocation="/afs/phas.gla.ac.uk/user/a/amorton/ilcsoft/v01-17-05/Eutelescope/trunk/jobsub/examples/GBL/noDUTExample"
export scriptsLocation="/afs/phas.gla.ac.uk/user/a/amorton/ilcsoft/v01-17-05/Eutelescope/trunk/GBL/iterativeAlignmentScripts"

#VARIABLES WHICH VARY THROUGH ALIGNMENT ARE SET HERE. IMPORTANT NOT ALL VARIABLES ARE HERE LOOK IN CONFIG-ALIGNMENT.
export r="0.006"; #Make resolution large so we begin with small chi2 and factor improvement to get to chi2/ndf=1 on next iteration. 
export xres="0.006 $r $r $r $r 0.006";
export yres="0.006 $r $r $r $r 0.006";
export minTracksPerEventAcceptance=0.001 #This is the number of tracks that is needed per event before we stop pattern recognition. Note value should depend on other cuts. 
export ResidualsRMax="0.3" #This is the window size on the next plane that we will accept hits from. This will increase if less than 1 track per event is found.
export inputGearInitial="gear_desy2012_150mm.xml" #Note the gear file changes through the process so must be placed here.
export allPlanes="0 1 2 3 4 5"
export allPlanesFixed="0  5"  


#Export input variables to other bash files.  
export RUN="$RUN" 
export numberOfIterations="$numberOfIterations"
export outputIdentifier="$identifier" #Use this string to identify final gear/histogram and all iterations before.


$scriptsLocation/initialiseAndRun.sh
