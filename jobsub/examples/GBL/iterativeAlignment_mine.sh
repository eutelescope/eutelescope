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

Fxr="0 1 2 3 4 5"
Fxs="0         5"
Fyr="0 1 2 3 4 5"
Fys="0         5"
Fzr="0 1 2 3 4 5"
Fzs="0 1 2 3 4 5"

#inputGear="gear_desy2012_150mm.xml"
inputGear="gear_lam_1T.xml"
outputGear="gear-final-XYshift-000${RUN}.xml"
histoNameInput="GBLtrack-XYshift-000${RUN}"

#inputGear="gear-final-XYshift-000${RUN}.xml"
#outputGear="gear-final-Zrotations2-000${RUN}.xml"
#histoNameInput="GBLtrack-zRotation2-000${RUN}"

#inputGear="gear-final-Zrotations-000${RUN}.xml"
#outputGear="gear-final-XRotation-000${RUN}.xml"
#histoNameInput="GBLtrack-xRotation"

#This is the alignment mode. It sets the size of the alignment jacobian dimensions.
amode="2";

pede="chiscut 15. 7. " #This is the input that tell millepede what tracks to discard.  

ExcludePlanes=""

r="0.0864 ";
xres="$r $r $r $r $r $r";
yres="$r $r $r $r $r $r";

prev="$r";
echo "prev:$prev and r:$r"
directoryTrack="/scratch/ilcsoft/v01-17-05/Eutelescope/master/jobsub/examples/GBL/output/logs"

if [ $# -ne 8 ]
then
 echo "$# parameters: $CONFIG $RUNLIST $Chi2Cut $RUN $file $gear10"
 exit
fi
#TO DO:For some reason when I drop collections this causes a segfault in alignment. I have no clue why. So I create many lcio files in this process.
#We first run pattern recognition
$do jobsub.py -c $CONFIG -csv $RUNLIST -o MaxRecordNumber="$MaxRecordNumber" -o GearFile="$inputGear"  -o ExcludePlanes="$ExcludePlanes" tracksearchHelix $RUN 

#Then we create the first tracks using pattern recognition tracks.
$do jobsub.py  -c $CONFIG -csv $RUNLIST -o GearFile="$inputGear" -o lcioInputName="trackcand"  -o inputCollectionName="track_candidates" -o lcioOutputName="GBLtracks-1" -o outputCollectionName="tracks1"  -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" GBLTrackFit  $RUN 
fileName="GBLTrackFit-000${RUN}.zip"
fullPath="$directoryTrack/$fileName"
echo "The full path to the log file is: $fullPath" 
averageChi2=`unzip  -p  $fullPath |grep "This is the average chi2 -" |cut -d '-' -f2`; 
echo "The average chi2ndf is :  $averageChi2"
if [[ $averageChi2 == "" ]]; then
	echo "ERROR!!!!!!!!! string for chi2 GBL not found. Check output log file is in the correct place. Furthermore check the string is there. "
  exit
fi
if [[ $(echo "$averageChi2 < 0.8"|bc) -eq 1 ]]; then
	echo "The average chi2 is: " $averageChi2. "So decrease resolution."		
	r=$(echo "scale=4;$r*0.8"|bc);
	xres="$r $r $r $r $r $r";
	yres="$r $r $r $r $r $r";
	echo "New resolutions are for (X/Y):" $xres"/"$yres
elif [[ $(echo "$averageChi2 > 1.2"|bc) -eq 1 ]]; then
	echo "The average chi2 is: " $averageChi2. "So increase resolution."		
	r=$(echo "scale=4;$r*1.2"|bc);
	xres="$r $r $r $r $r $r";
	yres="$r $r $r $r $r $r";
	echo "New resolutions are for (X/Y):" $xres"/"$yres
else 
	echo "The average chi2 is: " $averageChi2. "So keep resolution the same."		
fi

#Now we enter the GBLTrack loop. Here we use the GBLTrack to fit new track with a chi2 tending to 1.
lcioNumber=0
for x in {1..10}; do
	echo "Resolution inside GBLTrack loop beginning (X/Y):" $xres"/"$yres
	xnext=$(($x+1))
	echo "This is x and xnext at the start of the loop x: $x and xnext: $xnext "
	$do jobsub.py -c $CONFIG -csv $RUNLIST -o histoName="$histoNameInput" -o lcioInputName="GBLtracks-$x" -o lcioOutputName="GBLtracks-$xnext"  -o inputCollectionName="tracks$x" -o  outputCollectionName="tracks$xnext" -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" GBLTrackFit  $RUN 
	fileName="GBLTrackFit-000${RUN}.zip"
	fullPath="$directoryTrack/$fileName"
	echo "The full path to the log file is: $fullPath" 
#	rm "/scratch/ilcsoft/v01-17-05/Eutelescope/master/jobsub/examples/GBL/output/lcio/run000${RUN}-GBLtracks.slcio"
	averageChi2=`unzip  -p  $fullPath |grep "This is the average chi2 -" |cut -d '-' -f2`; 
	echo "The average chi2ndf is :  $averageChi2"
	if [[ $averageChi2 == "" ]]; then
		echo "ERROR!!!!!!!!! string for chi2 GBL not found. Check output log file is in the correct place. Furthermore check the string is there. "
		exit
	fi
	if [[ $(echo "$averageChi2 < 0.8"|bc) -eq 1 ]]; then
		echo "The average chi2 is: " $averageChi2. "So decrease resolution."		
		r=$(echo "scale=4;$r*0.8"|bc);
		xres="$r $r $r $r $r $r";
		yres="$r $r $r $r $r $r";
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	elif [[ $(echo "$averageChi2 > 1.2"|bc) -eq 1 ]]; then
		echo "The average chi2 is: " $averageChi2. "So increase resolution."		
		r=$(echo "scale=4;$r*1.2"|bc);
		xres="$r $r $r $r $r $r";
		yres="$r $r $r $r $r $r";
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	else 
		echo "The average chi2 is: " $averageChi2. "So keep resolution the same."		
		lcioNumber=$x
		break
	fi
#	echo "We are at the second part of the track loop"
#	$do jobsub.py  -c $CONFIG -csv $RUNLIST -o dropCollectionName="tracks2 tracks2_states" -o lcioInputName="GBLtracks2" -o lcioOutputName="GBLtracks" -o inputCollectionName="tracks2" -o  outputCollectionName="tracks"  -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" GBLTrackFit  $RUN 
##	rm "/scratch/ilcsoft/v01-17-05/Eutelescope/master/jobsub/examples/GBL/output/lcio/run000${RUN}-GBLtracks2.slcio"
#
#	averageChi2=`unzip  -p  $fileTrack |grep "This is the average chi2 -" |cut -d '-' -f2`; 
#	echo "The average chi2ndf is :  $averageChi2"
#	if [[ $averageChi2 == "" ]]; then
#		echo "ERROR!!!!!!!!! string for chi2 GBL not found. Check output log file is in the correct place. Furthermore check the string is there. "
#		exit
#	fi
#	if [[ $(echo "$averageChi2 < 0.5"|bc) -eq 1 ]]; then
#		echo "The average chi2 is: " $averageChi2. "So decrease resolution."		
#		r=$(echo "scale=4;$r*0.5"|bc);
#		xres="$r $r $r $r $r $r";
#		yres="$r $r $r $r $r $r";
#		echo "New resolutions are for (X/Y):" $xres"/"$yres
#	elif [[ $(echo "$averageChi2 > 3"|bc) -eq 1 ]]; then
#		echo "The average chi2 is: " $averageChi2. "So increase resolution."		
#		r=$(echo "scale=4;$r*2"|bc);
#		xres="$r $r $r $r $r $r";
#		yres="$r $r $r $r $r $r";
#		echo "New resolutions are for (X/Y):" $xres"/"$yres
#	else 
#		echo "The average chi2 is: " $averageChi2. "So keep resolution the same and exit loop."		
#		break
#	fi
done
if [[ $lcioNumber -eq 0 ]]; then
	echo "There is a problem with setting the lcio number for alignment"
	exit
fi
#Entering alignment steps########################################################################################
r="0.1 ";
xres="$r $r $r $r $r $r";
yres="$r $r $r $r $r $r";

fileAlign="/scratch/ilcsoft/v01-17-05/Eutelescope/master/jobsub/examples/GBL/output/logs/GBLAlign-000${RUN}.zip"

#Do one iteration of alignment to create new loop gear
$do jobsub.py -c $CONFIG -csv $RUNLIST -o lcioInputName="GBLtracks-$lcioNumber" -o inputCollectionName="tracks$lcioNumber" -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" -o GearFile="$inputGear" -o GearAlignedFile="$outputGear" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o MilleMaxChi2Cut="$Chi2Cut" -o pede="$pede" GBLAlign  $RUN 

error=`unzip  -p  $fileAlign |grep "Backtrace for this error:" | awk '{ print $NF }'`;
if [[ $error != "" ]];then
	echo "We have a segfault" 
	exit
fi

#Now we enter the alignment loop.
while :
	do
	echo "Resolution inside ALignment loop beginning (X/Y):" $xres"/"$yres
	$do jobsub.py -c $CONFIG -csv $RUNLIST -o lcioInputName="GBLtracks-$lcioNumber" -o inputCollectionName="tracks$lcioNumber" -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" -o GearFile="$inputGear" -o GearAlignedFile="gear-loop-000${RUN}.xml" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" -o AlignmentMode="$amode" -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}" -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o MilleMaxChi2Cut="$Chi2Cut" -o pede="$pede"  GBLAlign  $RUN

	rejected=`unzip  -p  $fileAlign |grep "Too many rejects" |cut -d '-' -f2`; 
	echo "Rejects word:  $rejected "
	if [[ $rejected != "" ]];then
		echo "Too many rejects. Resolution must increase by factor 2."
		r=$(echo "scale=4;$r*2"|bc);
		xres="$r $r $r $r $r $r";
		yres="$r $r $r $r $r $r";
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	fi
	factor=`unzip  -p  $fileAlign |grep "multiply all input standard deviations by factor" | awk '{ print $NF }'`;
	echo "factor word: " $factor
	if [[ $factor != "" ]];then
		echo "Factor word found! Resolution must increase by factor $factor "
		r=$(echo "scale=4;  $r*$factor" |bc);
		xres="$r $r $r $r $r $r";
		yres="$r $r $r $r $r $r";
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	fi
	$do jobsub.py -c $CONFIG -csv $RUNLIST -o lcioInputName="GBLtracks-$lcioNumber" -o inputCollectionName="tracks$lcioNumber" -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" -o GearFile="$inputGear" -o GearAlignedFile="$outputGear" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o MilleMaxChi2Cut="$Chi2Cut" -o pede="$pede"  GBLAlign  $RUN 
	rejected=`unzip  -p  $fileAlign |grep "Too many rejects" |cut -d '-' -f2`; 

	averageChi2Mille=`unzip -p $fileAlign |grep "Chi^2/Ndf" | awk '{ print $(NF-5) }'`;

	echo "Rejects word:  $rejected "
	if [[ $rejected != "" ]];then
		echo "Too many rejects. Resolution must increase by factor 2."
		r=$(echo "scale=4;$r*2"|bc);
		xres="$r $r $r $r $r $r";
		yres="$r $r $r $r $r $r";
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	fi
	factor=`unzip  -p  $fileAlign |grep "multiply all input standard deviations by factor" | awk '{ print $NF }'`;
	echo "factor word: " $factor
	if [[ $factor != "" ]];then
		echo "Factor word found! Resolution must increase by $factor."
		r=$(echo "scale=4;$r*$factor"|bc);
		xres="$r $r $r $r $r $r";
		yres="$r $r $r $r $r $r";
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	fi
	if [[ $averageChi2Mille == "" ]] && [[ $factor == "" ]] && [[ $rejected == "" ]];then
		echo "Mille chi2 is non existant. Here it is: $averageChi2Mille"
		echo "We can not find this or factor or rejects. Sit chi2 to 1 and exit"
		averageChi2Mille=1
	fi
#	averageChi2Mille=1
	if [[ $(echo "$averageChi2Mille < 0.8"|bc) -eq 1 ]]; then
		echo "The average chi2 is: " $averageChi2Mille. "So decrease resolution."		
		echo "New resolutions are for (X/Y):" $xres"/"$yres
		echo "Continue"
	elif [[ $(echo "$averageChi2Mille > 1.2"|bc) -eq 1 ]]; then
		echo "The average chi2 is: " $averageChi2Mille. "So increase resolution."		
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	else 
		echo "The average chi2 is: " $averageChi2Mille. "So keep resolution the same and exit loop."		
		break
	fi
done 


$do jobsub.py -c $CONFIG -csv $RUNLIST -o MaxRecordNumber="$MaxRecordNumber" -o GearFile="$outputGear"  -o ExcludePlanes="$ExcludePlanes" tracksearchHelix $RUN 
r="0.005";
xres="$r $r $r $r $r $r";
yres="$r $r $r $r $r $r";

#Then we create the first tracks using pattern recognition tracks.
$do jobsub.py  -c $CONFIG -csv $RUNLIST -o GearFile="$outputGear" -o histoName="$histoNameInput" -o lcioInputName="trackcand"  -o inputCollectionName="track_candidates" -o lcioOutputName="GBLtracks-1" -o outputCollectionName="tracks1"  -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" GBLTrackFit  $RUN 
