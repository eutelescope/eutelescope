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

r="0.125";
xres="$r $r $r $r $r $r";
yres="$r $r $r $r $r $r";

prev="$r";
echo "prev:$prev and r:$r"
fileTrack="/scratch/ilcsoft/v01-17-05/Eutelescope/master/jobsub/examples/GBL/output/logs/GBLTrackFit-0000${RUN}.zip"

if [ $# -ne 8 ]
then
 echo "$# parameters: $CONFIG $RUNLIST $Chi2Cut $RUN $file $gear10"
 exit
fi

#We first run pattern recognition
$do jobsub.py  $DRY -c $CONFIG -csv $RUNLIST -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" tracksearchHelix $RUN 

#Then we create the first tracks using pattern recognition tracks.
$do jobsub.py  $DRY -c $CONFIG -csv $RUNLIST -o lcioInputName="trackcand"  -o inputCollectionName="track_candidates" -o lcioOutputName="GBLtracks" -o outputCollectionName="tracks"  -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" GBLTrackFit  $RUN 

averageChi2=`unzip  -p  $fileTrack |grep "This is the average chi2 -" |cut -d '-' -f2`; 
echo "The average chi2ndf is :  $averageChi2"
if [[ $averageChi2 == "" ]]; then
	echo "ERROR!!!!!!!!! string for chi2 GBL not found. Check output log file is in the correct place. Furthermore check the string is there. "
  exit
fi
if [[ $(echo "$averageChi2 < 0.5"|bc) -eq 1 ]]; then
	echo "The average chi2 is: " $averageChi2. "So decrease resolution."		
	r=$(echo "scale=4;$r*0.5"|bc);
	xres="$r $r $r $r $r $r";
	yres="$r $r $r $r $r $r";
	echo "New resolutions are for (X/Y):" $xres"/"$yres
elif [[ $(echo "$averageChi2 > 3"|bc) -eq 1 ]]; then
	echo "The average chi2 is: " $averageChi2. "So increase resolution."		
	r=$(echo "scale=4;$r*2"|bc);
	xres="$r $r $r $r $r $r";
	yres="$r $r $r $r $r $r";
	echo "New resolutions are for (X/Y):" $xres"/"$yres
else 
	echo "The average chi2 is: " $averageChi2. "So keep resolution the same."		
fi

#Now we enter the GBLTrack loop. Here we use the GBLTrack to fit new track with a chi2 tending to 1.
while : 
	do
	echo "Resolution inside GBLTrack loop beginning (X/Y):" $xres"/"$yres
	$do jobsub.py -c $CONFIG -csv $RUNLIST -o lcioInputName="GBLtracks" -o lcioOutputName="GBLtracks2"  -o inputCollectionName="tracks" -o  outputCollectionName="tracks2" -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" -o dropCollectionName="tracks tracks_states"  GBLTrackFit  $RUN 
	rm "/scratch/ilcsoft/v01-17-05/Eutelescope/master/jobsub/examples/GBL/output/lcio/run0000${RUN}-GBLtracks.slcio"
	averageChi2=`unzip  -p  $fileTrack |grep "This is the average chi2 -" |cut -d '-' -f2`; 
	echo "The average chi2ndf is :  $averageChi2"
	if [[ $averageChi2 == "" ]]; then
		echo "ERROR!!!!!!!!! string for chi2 GBL not found. Check output log file is in the correct place. Furthermore check the string is there. "
		exit
	fi
	if [[ $(echo "$averageChi2 < 0.5"|bc) -eq 1 ]]; then
		echo "The average chi2 is: " $averageChi2. "So decrease resolution."		
		r=$(echo "scale=4;$r*0.5"|bc);
		xres="$r $r $r $r $r $r";
		yres="$r $r $r $r $r $r";
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	elif [[ $(echo "$averageChi2 > 3"|bc) -eq 1 ]]; then
		echo "The average chi2 is: " $averageChi2. "So increase resolution."		
		r=$(echo "scale=4;$r*2"|bc);
		xres="$r $r $r $r $r $r";
		yres="$r $r $r $r $r $r";
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	else 
		echo "The average chi2 is: " $averageChi2. "So keep resolution the same."		
	fi
	echo "We are at the second part of the track loop"
	$do jobsub.py  $DRY -c $CONFIG -csv $RUNLIST -o dropCollectionName="tracks2 tracks2_states" -o lcioInputName="GBLtracks2" -o lcioOutputName="GBLtracks" -o inputCollectionName="tracks2" -o  outputCollectionName="tracks"  -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" GBLTrackFit  $RUN 
	rm "/scratch/ilcsoft/v01-17-05/Eutelescope/master/jobsub/examples/GBL/output/lcio/run0000${RUN}-GBLtracks2.slcio"

	averageChi2=`unzip  -p  $fileTrack |grep "This is the average chi2 -" |cut -d '-' -f2`; 
	echo "The average chi2ndf is :  $averageChi2"
	if [[ $averageChi2 == "" ]]; then
		echo "ERROR!!!!!!!!! string for chi2 GBL not found. Check output log file is in the correct place. Furthermore check the string is there. "
		exit
	fi
	if [[ $(echo "$averageChi2 < 0.5"|bc) -eq 1 ]]; then
		echo "The average chi2 is: " $averageChi2. "So decrease resolution."		
		r=$(echo "scale=4;$r*0.5"|bc);
		xres="$r $r $r $r $r $r";
		yres="$r $r $r $r $r $r";
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	elif [[ $(echo "$averageChi2 > 3"|bc) -eq 1 ]]; then
		echo "The average chi2 is: " $averageChi2. "So increase resolution."		
		r=$(echo "scale=4;$r*2"|bc);
		xres="$r $r $r $r $r $r";
		yres="$r $r $r $r $r $r";
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	else 
		echo "The average chi2 is: " $averageChi2. "So keep resolution the same and exit loop."		
		break
	fi
done
#Entering alignment steps########################################################################################
fileAlign="/scratch/ilcsoft/v01-17-05/Eutelescope/master/jobsub/examples/GBL/output/logs/GBLAlign-0000${RUN}.zip"

Fxr="0 1 2 3 4 5"
Fxs="0         5"
Fyr="0 1 2 3 4 5"
Fys="0         5"
Fzr="0 1 2 3 4 5"
Fzs="0 1 2 3 4 5"

#This is the alignment mode. It sets the size of the alignment jacobian dimensions.
amode="3";

pede="chiscut 50. 5. " #This is the input that tell millepede what tracks to discard.  
Chi2Cut="5." #This is the input tracks we allow to be used in alignment. Note that this relatively low value is possible since we loop over the track fit many times to get low chi2

#Do one iteration of alignment to create new loop gear
$do jobsub.py -c $CONFIG -csv $RUNLIST -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" -o GearAlignedFile="gear-final-0000${RUN}.xml" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o MilleMaxChi2Cut="$Chi2Cut" -o pede="$pede" GBLAlign  $RUN 


#Now we enter the alignment loop.
while :
	do
	echo "Resolution inside ALignment loop beginning (X/Y):" $xres"/"$yres
	$do jobsub.py -c $CONFIG -csv $RUNLIST -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" -o GearFile="gear-final-0000${RUN}.xml" -o GearAlignedFile="gear-loop-0000${RUN}.xml" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" -o AlignmentMode="$amode" -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}" -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o MilleMaxChi2Cut="$Chi2Cut" -o pede="$pede"  GBLAlign  $RUN

	rejected=`unzip  -p  $fileAlign |grep "Too many rejects" |cut -d '-' -f2`; 
	echo "Rejects word found: " $rejected
	if [[ $rejected != "" ]];then
		echo "Too many rejects. Resolution must increase."
		r=$(echo "scale=4;$r*2"|bc);
		xres="$r $r $r $r $r $r";
		yres="$r $r $r $r $r $r";
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	fi
	factor=`unzip  -p  $fileAlign |grep "multiply all input standard deviations by factor" |cut -d ' ' -f2`;
	echo "factor word found: " $factor
	if [[ $factor != "" ]];then
		echo "Too many rejects. Resolution must increase."
		r=$(echo "scale=4;$r*$factor"|bc);
		xres="$r $r $r $r $r $r";
		yres="$r $r $r $r $r $r";
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	fi
	$do jobsub.py -c $CONFIG -csv $RUNLIST -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" -o GearFile="gear-loop-0000${RUN}.xml " -o GearAlignedFile="gear-final-0000${RUN}.xml" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o MilleMaxChi2Cut="$Chi2Cut" -o pede="$pede"  GBLAlign  $RUN 
	rejected=`unzip  -p  $fileAlign |grep "Too many rejects" |cut -d '-' -f2`; 
	echo "Rejects word found: " $rejected
	if [[ $rejected != "" ]];then
		echo "Too many rejects. Resolution must increase."
		r=$(echo "scale=4;$r*2"|bc);
		xres="$r $r $r $r $r $r";
		yres="$r $r $r $r $r $r";
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	fi
	factor=`unzip  -p  $fileAlign |grep "multiply all input standard deviations by factor" |cut -d ' ' -f2`;
	echo "factor word found: " $factor
	if [[ $factor != "" ]];then
		echo "Too many rejects. Resolution must increase."
		r=$(echo "scale=4;$r*$factor"|bc);
		xres="$r $r $r $r $r $r";
		yres="$r $r $r $r $r $r";
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	fi
	if [[ $(echo "$averageChi2 < 0.5"|bc) -eq 1 ]]; then
		echo "The average chi2 is: " $averageChi2. "So decrease resolution."		
		echo "New resolutions are for (X/Y):" $xres"/"$yres
		echo "Continue"
	elif [[ $(echo "$averageChi2 > 3"|bc) -eq 1 ]]; then
		echo "The average chi2 is: " $averageChi2. "So increase resolution."		
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	else 
		echo "The average chi2 is: " $averageChi2. "So keep resolution the same and exit loop."		
		break
	fi
done 
