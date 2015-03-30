#!/bin/bash
#This is the bash script that performs a single iteration of alignment. 
#The steps are: 
#1) It will run pattern recognition and check that the number of found track candidates is reasonable if not then it will widen its acceptance window and try again.
#2) Using these track candidates we will then produce tracks with the estimates resolution of the planes and DUT. This is approximately 1mm since it is misaligned. We also determine a average chi2 to use as our cut in alignment. 
#3) With these tracks we then use this to align. Alignment uses 2 different kinds of cuts. We cut tracks with high chi2 before alignment procedure and during.   
#4) If our chi2/ndf is not rougly equal to 1 or we have a high number of rejects then we run again with a different resolution.  
#TO DO: If we change the resolution chi2 cut should also change. However at the moment we keep this constant for many iteration of alignment.
#THIS IS PART (1)
echo "The input to single loop on iteration $number"
echo "Input gear: $inputGear" 
echo "Output gear: $outputGear"
echo "This is the resolutions X/Y:  $xres/$yres."
#for x in {1..10}; do 
#	echo "PATTERN RECOGNTION ATTEMPT $x ON ITERATION $number"
	$do jobsub.py -c $CONFIG -csv $RUNLIST  -o ResidualsRMax="$ResidualsRMax" -o GearFile="$inputGear" -o TrackCandHitOutputCollectionName="$PatOutGBLIn" $PatRec $RUN  

#	fileName="$PatRec-${RUN}.zip"
#	fullPath="$directory/$fileName"
#	echo "The full path to the log file is: $fullPath" 
#	averageTracksPerEvent=`unzip  -p  $fullPath |grep "The average number of tracks per event: " |cut -d ':' -f2`; 
#	echo "The average number of tracks per event from the log file is $averageTracksPerEvent"
#	if [[ $averageTracksPerEvent == "" ]];then
#		echo "averagetrackPerEvent variable is empty"		
#	fi
#	if [[ $(echo "$averageTracksPerEvent >$minTracksPerEventAcceptance "|bc) -eq 1 ]]; then
#		break
#	fi
	#If we reach this part of the code we need to increase the window size.
#	echo "$ResidualsRMax is the size of the radius of acceptance in this iteration"
#	export ResidualsRMax=$(echo "scale=4;$ResidualsRMax*$patRecMultiplicationFactor"|bc)    
#	echo "$ResidualsRMax is the size of the radius of acceptance in the next iteration"
	#Did we find the correct number of track within x=10 iterations.
#	if [[ $x -eq 10 ]];then
#		echo "We are are iteration 10 and still have not found enough tracks in pattern recogntion"
#		exit 1
#	fi
#done
#THIS IS PART (2)
#This part should analyse the output of pattern recogntion and remove tracks which are clearly poor quality. This is difficult with the systematic problems due to misalignment. 
#Should be a separate processor so we can possible compare pattern recognition techniques or chain. Future work. 
#echo "GBLTRACKS CREATED FOR ALIGNMENT ON ITERATION $number"
#$do jobsub.py  -c $CONFIG -csv $RUNLIST  -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" -o GearFile="$inputGear"  $TrackFit  $RUN  

#THIS IS PART (3)
upperLimit=1 #So we only accept decreases in resolution form millepede. 
#if there is a estimate of -90%, an increase in resolution, then use this result. 
fileAlign="$directory/$Align-${RUN}.zip"
numberRejectedAlignmentAttempts=0 #We set this since we do not want to fall in a loop were chi<1 the rejects then chi<1 .....
dry="--dry-run"
alignment=false
rejectFactor=2
echo "Enter while loop here"
while [ "$alignment" == "false" ]; do 
	echo "THE NUMBER OF FAILED ALIGNMENT ATTEMPTS $numberRejectedAlignmentAttempts"
	$do jobsub.py -c $CONFIG -csv $RUNLIST -o GearFile="$inputGear" -o GearAlignedFile="$outputGear" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres"  -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o lcioInputName="trackcand" -o inputCollectionName="$PatOutGBLIn" $Align  $RUN  
	#Fill variables.
	error1=`unzip  -p  $fileAlign |grep "Backtrace for this error:" | awk '{ print $NF }'`;
	error2=`unzip  -p  $fileAlign |grep "This is the entire stack" | awk '{ print $NF }'`;
	error3=`unzip  -p  $fileAlign |grep "problem with stream or EOF reading" | awk '{ print $NF }'`;
	error4=`unzip  -p  $fileAlign |grep "STOP STOP" | awk '{ print $NF }'`;

	rejected=`unzip  -p  $fileAlign |grep "Too many rejects" |cut -d '-' -f2`; 
	averageChi2Mille=`unzip -p $fileAlign |grep "Chi^2/Ndf" | awk '{ print $(NF-5) }'`;
	averageChi2Mille=$(echo $averageChi2Mille | cut -f1 -d' ')
	factor=`unzip  -p  $fileAlign |grep "multiply all input standard deviations by factor" | awk '{ print $NF }'`;
	factor=$(echo $factor | cut -f1 -d' ')
	factor=`echo ${factor} | sed -e 's/[eE]+*/\\*10\\^/'`;
	noCut=`echo " $averageChi2Mille < 10.0" | bc`; #Do not want to decrease the resolution too much
	notUnderEstimated=`echo " $averageChi2Mille > 1.0" | bc ` #Do not under estimate the errors.
	#We can either have a error, rejected, factor, or nothing.
	echo $error1 $error2 $error3
	if [ "$error1" != "" ] || [ "$error2" != "" ] || [ "$error3" != "" ] || [ "$error4" != "" ]  ;then #Must put quotes if we expect error1/2 to be empty
		echo "We have a segfault" 
		break
	elif [ "$rejected" != "" ];then
		export	numberRejectedAlignmentAttempts=$(($numberRejectedAlignmentAttempts+1))
		echo "Too many rejects. Resolution must increase by factor $rejectFactor."
		xInput=$xres;
		yInput=$yres;
		unset xres;
		unset yres;
		echo "Rejection factor is: $rejectFactor"
		xres=`python $pythonLocation/multiplyResolutionsByFactor.py $xInput /   / $allPlanes / $rejectFactor` 
		yres=`python $pythonLocation/multiplyResolutionsByFactor.py $yInput /   / $allPlanes / $rejectFactor`
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	elif [ "$noCut" == "1" ] && [ "$notUnderEstimated" == "1" ];then #Chi2 cut and must not understimate. Range [0,10]
		echo "Chi2: $averageChi2Mille and boolean: $noCut $notUnderEstimated"
		#If the factor word is out of range then we use this fit anyway.
		alignment=true
		echo "Set alignment to equal: $alignment"
	else
		echo "Chi2: $averageChi2Mille and boolean: $noCut $notUnderEstimated"
		if [[ "$factor" != "" ]];then
			echo "Here is the factor: $factor " 
			passUpper=`expr "$factor" '<' "$upperLimit"` #Do not want to decrease the resolution too much
			echo $passUpper
			if [ "$passUpper" -ne "0"  ];then
				echo "Factor word found! Resolution must increase by $factor."
				#for some reason the output of python will not overwrite the xres or yres? So must unset then set.
				xInput=$xres
				yInput=$yres
				unset xres;
				unset yres;
				xres=`python $pythonLocation/multiplyResolutionsByFactor.py $xInput / / $allPlanes / $factor`
				yres=`python $pythonLocation/multiplyResolutionsByFactor.py $yInput / / $allPlanes / $factor`
				echo "New resolutions are for (X/Y):" $xres"/"$yres
				#Also decrease the rejection factor since we must be close to a fit.
				rejectFactor=1.01
			else
				#If there is no factor to decrease resolution then just use this fit.
				alignment=true
			fi
		else
			alignment=true
		fi
	fi
	if [[ $numberRejectedAlignmentAttempts -eq 10 ]]
	then
		echo "We have already rejected $numberRejectedAlignmentAttempts times."
		break
	fi
done 
#We exit with an error if no alignment.
if [ "$alignment" == "false" ];then
	exit 1
fi
