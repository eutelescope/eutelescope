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
###################################################################################Functions
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
########################################################################3
RUN=$(adding_zeros_to_RUN $RUN)

PatRec=tracksearchHelix-apix
TrackFit=GBLTrackFit-apix
Align=GBLAlign-apix
#PatRec=tracksearchHelix
#TrackFit=GBLTrackFit
#Align=GBLAlign

MaxRecordNumber="3" 

echo "Input recieved"
echo "Config: " $CONFIG
echo "Runlist: " $RUNLIST
echo "Chi2Cut: " $Chi2Cut
echo "Run: " $RUN 
echo "MaxRecordNumber: " $MaxRecordNumber

#Fxr="0 1 2 3 4 5"
#Fxs="0         5"
#Fyr="0 1 2 3 4 5"
#Fys="0         5"
#Fzr="0 "
#Fzs="0 1 2 3 4 5"

Fxr="0 1 2 3 4 5 20"
Fxs="0           "
Fyr="0 1 2 3 4 5 20"
Fys="0           "
Fzr="0           20"
Fzs="0 1 2 3 4 5 20"

#inputGear="gear_desy2012_150mm.xml"
#inputGear="gear_lam_1T.xml"
#inputGear="gear_150mm_1fei4_30.xml"
#outputGear="gear-final-XYshift-${RUN}.xml"
#histoNameInput="GBLtrack-XYshift-${RUN}"

#inputGear="gear-final-XYshift-${RUN}.xml"

#outputGear="gear-final-XYshiftS2-${RUN}.xml"
#histoNameInput="GBLtrack-XYshiftS2-${RUN}"

#inputGear="gear-final-Zrotations-${RUN}.xml"
inputGear="gear-final-ZRotation-${RUN}.xml"
#outputGear="gear-final-ZRotation-${RUN}.xml"
outputGear="gear-final-XYShift-DUT-${RUN}.xml"
histoNameInput="GBLtrack-XYshiftS-DUT-${RUN}"

#histoNameInput="GBLtrack-zRotation2-${RUN}"

#This is the alignment mode. It sets the size of the alignment jacobian dimensions.
amode="7";

pede="chiscut 10. 5. " #This is the input that tell millepede what tracks to discard.  

#ExcludePlanes="20"
ExcludePlanes=

r="0.0828";
xres="$r $r $r $r $r $r $r";
yres="$r $r $r $r $r $r $r";

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
$do jobsub.py -c $CONFIG -csv $RUNLIST -o MaxRecordNumber="$MaxRecordNumber" -o GearFile="$inputGear"  -o ExcludePlanes="$ExcludePlanes" $PatRec $RUN 


#Then we create the first tracks using pattern recognition tracks. We iterate until the chi2 is close to one. ############################
$do jobsub.py  -c $CONFIG -csv $RUNLIST -o GearFile="$inputGear" -o lcioInputName="trackcand"  -o inputCollectionName="track_candidates" -o lcioOutputName="GBLtracks-1" -o outputCollectionName="tracks1"  -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" $TrackFit  $RUN 
fileName="$TrackFit-${RUN}.zip"
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
	xres="$r $r $r $r $r $r $r";
	yres="$r $r $r $r $r $r $r";
	echo "New resolutions are for (X/Y):" $xres"/"$yres
elif [[ $(echo "$averageChi2 > 1.2"|bc) -eq 1 ]]; then
	echo "The average chi2 is: " $averageChi2. "So increase resolution."		
	r=$(echo "scale=4;$r*1.2"|bc);
	xres="$r $r $r $r $r $r $r";
	yres="$r $r $r $r $r $r $r";
	echo "New resolutions are for (X/Y):" $xres"/"$yres
else 
	echo "The average chi2 is: " $averageChi2. "So keep resolution the same."		
fi

#Now we enter the GBLTrack loop. Here we use the GBLTrack to fit new track with a chi2 tending to 1.
for x in {1..10}; do
	echo "Resolution inside GBLTrack loop beginning (X/Y):" $xres"/"$yres
	xnext=$(($x+1))
	echo "This is x and xnext at the start of the loop x: $x and xnext: $xnext "
	$do jobsub.py -c $CONFIG -csv $RUNLIST -o histoName="$histoNameInput" -o lcioInputName="GBLtracks-$x" -o lcioOutputName="GBLtracks-$xnext"  -o inputCollectionName="tracks$x" -o  outputCollectionName="tracks$xnext" -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" $TrackFit  $RUN 
	fileName="$TrackFit-${RUN}.zip"
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
		xres="$r $r $r $r $r $r $r";
		yres="$r $r $r $r $r $r $r";
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	elif [[ $(echo "$averageChi2 > 1.2"|bc) -eq 1 ]]; then
		echo "The average chi2 is: " $averageChi2. "So increase resolution."		
		r=$(echo "scale=4;$r*1.2"|bc);
		xres="$r $r $r $r $r $r $r";
		yres="$r $r $r $r $r $r $r";
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	else 
		echo "The average chi2 is: " $averageChi2. "So keep resolution the same."		
		break
	fi
done

############################################################################################################################################
######################################################################################################################Now we make sure the correction to the GBL tracks are close to 0 as possible
#echo "Now to reduce the corrections to near 1."
#while :
#	do
#	echo "Resolution inside GBLTrack corrections loop beginning (X/Y):" $xres"/"$yres
#	echo "We are on the $x iteration of lcio file" 
#	xnext=$(($x+1))
#	echo "This is x and xnext at the start of the loop x: $x and xnext: $xnext "
#	$do jobsub.py -c $CONFIG -csv $RUNLIST -o histoName="$histoNameInput" -o lcioInputName="GBLtracks-$x" -o lcioOutputName="GBLtracks-$xnext"  -o inputCollectionName="tracks$x" -o  outputCollectionName="tracks$xnext" -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" $TrackFit  $RUN 
#	fileName="$TrackFit-${RUN}.zip"
#	fullPath="$directoryTrack/$fileName"
#	echo "The full path to the log file is: $fullPath" 
#	correctionOmega=`unzip  -p  $fullPath |grep "This is the average correction for omega:" |cut -d '-' -f2`; 
#	correctionXZInclination=`unzip  -p  $fullPath |grep "This is the average correction for local xz inclination:" |cut -d ':' -f2`; 
#	correctionYZInclination=`unzip  -p  $fullPath |grep "This is the average correction for local yz inclination:" |cut -d ':' -f2`; 
#	correctionX=`unzip  -p  $fullPath |grep "This is the average correction for local x: " |cut -d ':' -f2`; 
#	correctionY=`unzip  -p  $fullPath |grep "This is the average correction for local y: " |cut -d ':' -f2`; 
#	echo "The average corrections:  $correctionOmega ,  $correctionXZInclination , $correctionYZInclination , $correctionX , $correctionY"   
#	if [[ $correctionOmega == "" ]] || [[ $correctionXZInclination == "" ]] || [[ $correctionYZInclination == "" ]] || [[ $correctionX == "" ]] || [[ $correctionY == "" ]];then
#		echo "One of the corrections is empty this must be wrong. Exit!!!!!11"
#		exit
#	fi
#	echo "Now to check if the corrections are small enough"
#	if [[ $(echo "$correctionOmega < 0.001"|bc) -eq 1 ]];then # && [[ $(echo "$correctionYZInclination < 0.001"|bc) -eq 1  ]] && [[ $(echo "$correctionX < 0.001"|bc) -eq 1  ]] && [[ $(echo "$correctionY < 0.001"|bc) -eq 1  ]];then
#		echo "The corrections are now small enough. Lets now align. "
#	break	
#	fi
#	x=$(($x+1))
#	echo "The corrections are not small enough. The new x is $x"
#done
#Entering alignment steps########################################################################################
echo "Now begin alignment"
fileAlign="/scratch/ilcsoft/v01-17-05/Eutelescope/master/jobsub/examples/GBL/output/logs/$Align-${RUN}.zip"


#Now we enter the alignment loop.
while :
	do
	echo "Resolution inside ALignment loop beginning (X/Y):" $xres"/"$yres
	$do jobsub.py -c $CONFIG -csv $RUNLIST -o lcioInputName="GBLtracks-$x" -o inputCollectionName="tracks$x" -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" -o GearFile="$inputGear" -o GearAlignedFile="$outputGear" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o MilleMaxChi2Cut="$Chi2Cut" -o pede="$pede"  $Align  $RUN 

	error=`unzip  -p  $fileAlign |grep "Backtrace for this error:" | awk '{ print $NF }'`;
	if [[ $error != "" ]];then
		echo "We have a segfault" 
		exit
	fi
	rejected=`unzip  -p  $fileAlign |grep "Too many rejects" |cut -d '-' -f2`; 
	averageChi2Mille=`unzip -p $fileAlign |grep "Chi^2/Ndf" | awk '{ print $(NF-5) }'`;
	echo "Rejects word:  $rejected "
	if [[ $rejected != "" ]];then
		echo "Too many rejects. Resolution must increase by factor 2."
		r=$(echo "scale=4;$r*2"|bc);
		xres="$r $r $r $r $r $r $r";
		yres="$r $r $r $r $r $r $r";
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	fi
	factor=`unzip  -p  $fileAlign |grep "multiply all input standard deviations by factor" | awk '{ print $NF }'`;
	echo "factor word: " $factor
	if [[ $factor != "" ]];then
		echo "Factor word found! Resolution must increase by $factor."
		r=$(echo "scale=4;$r*$factor"|bc);
		xres="$r $r $r $r $r $r $r";
		yres="$r $r $r $r $r $r $r";
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	fi
	if [[ $averageChi2Mille == "" ]] && [[ $factor == "" ]] && [[ $rejected == "" ]];then
		echo "Mille chi2 is non existant. Here it is: $averageChi2Mille"
		echo "We can not find this or factor or rejects. Sit chi2 to 1 and exit"
		averageChi2Mille=1
	fi
#	averageChi2Mille=1
	if [[ $(echo "$averageChi2Mille < 0.8"|bc) -eq 1 ]] && [[ $averageChi2Mille != "" ]]; then
		echo "The average chi2 is: " $averageChi2Mille. "So decrease resolution."		
		echo "New resolutions are for (X/Y):" $xres"/"$yres
		echo "Continue"
	elif [[ $(echo "$averageChi2Mille > 1.2"|bc) -eq 1 ]] && [[ $averageChi2Mille != "" ]]; then
		echo "The average chi2 is: " $averageChi2Mille. "So increase resolution."		
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	else 
		echo "The average chi2 is: " $averageChi2Mille. "So keep resolution the same and exit loop."		
		break
	fi
done 


$do jobsub.py -c $CONFIG -csv $RUNLIST -o MaxRecordNumber="$MaxRecordNumber" -o GearFile="$outputGear"  -o ExcludePlanes="$ExcludePlanes" $PatRec $RUN 
r="0.005";
xres="$r $r $r $r $r $r $r";
yres="$r $r $r $r $r $r $r";

#Then we create the first tracks using pattern recognition tracks.
x=$(($x+1))
$do jobsub.py  -c $CONFIG -csv $RUNLIST -o GearFile="$outputGear" -o histoName="$histoNameInput" -o lcioInputName="trackcand"  -o inputCollectionName="track_candidates" -o lcioOutputName="GBLtracks-$x" -o outputCollectionName="tracks$x"  -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" $TrackFit  $RUN
