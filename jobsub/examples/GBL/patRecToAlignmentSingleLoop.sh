#!/bin/bash
#This is the bash script that performs a single iteration of alignment. 
#It does: 
#1) It will run pattern recognition and check that the number of found track candidates is reasonable if not then it will widen its acceptance window and try again.
#2) Using these track candidates we will then produce tracks with the estimates resolution of the planes and DUT. 0.005 for mimosa planes. 
#3) With these tracks we then use this to align and perform a chi2 cut. We repeat with different chi2 cuts until the chi2 is approximately the same each iteration.
#		As we reduce the chi2 cut we reduce the tracks provided to millepede AND the threshold when millepede will reject tracks. 
#		Therefore we would expect the chi2 to reduce constantly BUT falling very quicky when we are including rubbish tracks.
#		We would expect this since the error on an actual track is due to the intrinsic resolution and misalignment which is roughly constant
#   However noise will be completely random and will most likely be much larger than the error due to resolution and misalignment.
#   Another way to know you have discarded the background is if the chi2 is below one. Assuming you residual is ONLY the intrisic resolution.
#4) The repeat if the chi2 of alignment is too high or has changed drastically from the iteration before. Otherwise end.

#THIS IS PART (1)
for x in {1..10}; do

	$do jobsub.py -c $CONFIG -csv $RUNLIST -o MaxMissingHitsPerTrack="$MaxMissingHitsPerTrack" -o AllowedSharedHitsOnTrackCandidate="$AllowedSharedHitsOnTrackCandidate" -o ResidualsRMax="$ResidualsRMax" -o HitInputCollectionName="$lcioPatternCollection" -o Verbosity="$Verbosity" -o Verbosity="$Verbosity" -o planeDimensions="${planeDimensions}" -o MaxRecordNumber="$MaxRecordNumber" -o GearFile="$inputGear"  -o ExcludePlanes="$ExcludePlanes" $PatRec $RUN 

	fileName="$PatRec-${RUN}.zip"
	fullPath="$directory/$fileName"
	echo "The full path to the log file is: $fullPath" 
	averageTracksPerEvent=`unzip  -p  $fullPath |grep "The average number of tracks per event: " |cut -d ':' -f2`; 
	echo "The average number of tracks per event from the log file is $averageTracksPerEvent"
	if [[ $averageTracksPerEvent == "" ]];then
		echo "averagetrackPerEvent variable is empty"		
	fi
	if [[ $(echo "$averageTracksPerEvent >$minTracksPerEventAcceptance "|bc) -eq 1 ]]; then
		break
	fi
	#If we reach this part of the code we need to increase the window size.
	export ResidualsRMax=$(echo "scale=4;$ResidualsRMax*$patRecMultiplicationFactor"|bc)    
done

#THIS IS PART (2)
$do jobsub.py  -c $CONFIG -csv $RUNLIST -o Verbosity="$Verbosity" -o GearFile="$inputGear" -o lcioInputName="trackcand"  -o inputCollectionName="track_candidates" -o lcioOutputName="GBLtracks" -o outputCollectionName="tracks"  -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" $TrackFit  $RUN 

fileName="$TrackFit-${RUN}.zip"
fullPath="$directory/$fileName"
echo "The full path to the log file is: $fullPath" 
averageChi2=`unzip  -p  $fullPath |grep "This is the average chi2 -" |cut -d '-' -f2`; 
echo "The average chi2ndf from the log is :  $averageChi2"
if [[ $averageChi2 == "" ]]; then
	echo "ERROR!!!!!!!!! string for chi2 GBL not found. Check output log file is in the correct place. Furthermore check the string is there. "
  exit
fi
Chi2Cut="$averageChi2"
#We must make sure this cut is close to the average so we do not cut too many
#tracks
fraction=$(echo "scale=4;$Chi2Cut*0.0333"|bc); #Divided by 30(0.0333) since
#this is close to the value of chi2 that is 3 standard deviations away. 
export pede="!chiscut $fraction  $fraction" #! denotes a comment in the steering file we remove this to activate this functionality. TO DO: Must comment below as well must fix
#this.
#This is the factor in iteration 1 and
#iteration 2 we times by the chi2 corresponding to max/min value that
#would be 3 standard derviation away from the mean of n normal gaussians. 
#Remember if our errors are purely statistically and not systematic like
#alignment uncertainties then our 6 telescople planes and DUT will act like
# 6 + 1 = 7 normal gaussian distributions (The 1 of course comes from the DUT).
# Note that the chi2cut effectively just cuts tracks that are deemed bad
# overall. Also note that these cuts are done after the track fit and do not
# change how determine that chi2. Outliers on the otherhand will fit individual
# tracks again to improve the chi2 by finding good tracks labelled as bad by a
# single data point of noise that produces a large chi2 which is above the 5 or 2.5 cut
# above.
export outlierdownweighting="!outlierdownweighting 0"


#THIS IS PART (3)
fileAlign="$directory/$Align-${RUN}.zip"
for x in {1..10}; do
	echo "We are on loop number $x"
	echo "The Chi2Cut variable is $Chi2Cut"
	echo "The pede input string is:  $pede"
	$do jobsub.py -c $CONFIG -csv $RUNLIST -o Verbosity="$Verbosity" -o lcioInputName="GBLtracks" -o inputCollectionName="tracks" -o MaxRecordNumber="$MaxRecordNumber" -o ExcludePlanes="$ExcludePlanes" -o GearFile="$inputGear" -o GearAlignedFile="$outputGear" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o chi2Cut="$Chi2Cut" -o pede="$pede" -o outlierdownweighting="$outlierdownweighting"  $Align  $RUN 

	error=`unzip  -p  $fileAlign |grep "Backtrace for this error:" | awk '{ print $NF }'`;
	if [[ $error != "" ]];then
		echo "We have a segfault" 
		exit
	fi
	factor=`unzip  -p  $fileAlign |grep "multiply all input standard deviations by factor" | awk '{ print $NF }'`;
	factor=`echo ${factor} | sed -e 's/[eE]+*/\\*10\\^/'`
	echo "factor word: " $factor
	if [[ $factor != "" ]];then
		echo "Factor word found! Resolution must increase by $factor."
		r=$(echo "scale=4;$r*$factor"|bc);
		xres="$r $r $r $dut $r $r $r";
		yres="$r $r $r $dut $r $r $r";
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	fi
	rejected=`unzip  -p  $fileAlign |grep "Too many rejects" |cut -d '-' -f2`; 
	averageChi2Mille=`unzip -p $fileAlign |grep "Chi^2/Ndf" | awk '{ print $(NF-5) }'`;
	echo "Rejects word:  $rejected "
	if [[ $rejected != "" ]];then #This makes sure that we do not cut too many tracks.
		echo "Too many rejects. Resolution must increase by factor 2."
		r=$(echo "scale=4;$r*2"|bc);
		xres="$r $r $r $dut $r $r $r";
		yres="$r $r $r $dut $r $r $r";
		echo "New resolutions are for (X/Y):" $xres"/"$yres
		echo "Too many rejects. Chi2 cut must increase."
		export Chi2Cut=$(echo "scale=4;$Chi2Cut*1.2"|bc);
		fraction=$(echo "scale=4;$Chi2Cut*0.0333"|bc);  
		export pede="!chiscut $fraction  $fraction" 
	fi
#	if [[ $(echo "$averageChi2Mille <$minChi2AlignAcceptance "|bc) -eq 1 ]] && [[ $averageChi2Mille != "" ]]; then
#		echo "The average chi2 is:  $averageChi2Mille. This is acceptable so finish."		
#		break
#	fi
	if [[ $averageChi2Mille == "" ]] && [[ $factor == "" ]] && [[ $rejected == "" ]];then
		echo "Mille chi2 is non existant. Here it is: $averageChi2Mille"
		echo "We can not find this or factor or rejects."
		break
	fi

done 
