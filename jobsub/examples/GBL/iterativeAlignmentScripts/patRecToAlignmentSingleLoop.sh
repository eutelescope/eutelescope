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
for x in {1..10}; do 
	echo "PATTERN RECOGNTION ATTEMPT $x ON ITERATION $number"
	$do jobsub.py -c $CONFIG -csv $RUNLIST  -o ResidualsRMax="$ResidualsRMax" -o GearFile="$inputGear"  $PatRec $RUN  

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
	echo "$ResidualsRMax is the size of the radius of acceptance in this iteration"
	export ResidualsRMax=$(echo "scale=4;$ResidualsRMax*$patRecMultiplicationFactor"|bc)    
	echo "$ResidualsRMax is the size of the radius of acceptance in the next iteration"
	#Did we find the correct number of track within x=10 iterations.
	if [[ $x -eq 10 ]];then
		echo "We are are iteration 10 and still have not found enough tracks in pattern recogntion"
		exit
	fi
done
#THIS IS PART (2)
#This part should analyse the output of pattern recogntion and remove tracks which are clearly poor quality. This is difficult with the systematic problems due to misalignment. 
#Should be a separate processor so we can possible compare pattern recognition techniques or chain. Future work. 
#At the moment we just run this through the GBL fitter to improve the tracks. 
#echo "GBLTRACKS CREATED FOR ALIGNMENT ON ITERATION $number"
$do jobsub.py  -c $CONFIG -csv $RUNLIST  -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" -o GearFile="$inputGear"  $TrackFit  $RUN  

#fileName="$TrackFit-${RUN}.zip"
#fullPath="$directory/$fileName"
#echo "The full path to the log file is: $fullPath" 
#averageChi2=`unzip  -p  $fullPath |grep "This is the average chi2 -" |cut -d '-' -f2`; 
#echo "The average chi2ndf from the log is :  $averageChi2"
#if [[ $averageChi2 == "" ]]; then
#	echo "ERROR!!!!!!!!! string for chi2 GBL not found. Check output log file is in the correct place. Furthermore check the string is there. "
#  exit
#fi


#THIS IS PART (3)
fileAlign="$directory/$Align-${RUN}.zip"
numberRejectedAlignmentAttempts=0 #We set this since we do not want to fall in a loop were chi<1 the rejects then chi<1 .....
tooManyRejectsExitLoopBool=false
dry="--dry-run"
#We use the last successful attempt when we have two rejected otherwise if we have one then we increase the resolution and continue.
# TO DO: Must change this to while loop since we exit down below as well. 
for x in {1..10}; do 
	echo "GBLALIGN ATTEMPT $x ON ITERATION $number"
	echo "THE NUMBER OF FAILED ALIGNMENT ATTEMPTS $numberRejectedAlignmentAttempts"
	$do jobsub.py -c $CONFIG -csv $RUNLIST -o GearFile="$inputGear" -o GearAlignedFile="$outputGear" -o xResolutionPlane="$xres" -o yResolutionPlane="$yres"  -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}"  $Align  $RUN  
	#Check that we have not seg fault within millepede.
	error=`unzip  -p  $fileAlign |grep "Backtrace for this error:" | awk '{ print $NF }'`;
	if [[ $error != "" ]];then
		echo "We have a segfault" 
		exit
	fi
	#What overall chi2 did we get from the fit. This is different from the chi2 of the individual tracks. 
	averageChi2Mille=`unzip -p $fileAlign |grep "Chi^2/Ndf" | awk '{ print $(NF-5) }'`;
	echo "The chi2/ndf of the fit is $averageChi2Mille" 
	if $tooManyRejectsExitLoopBool;then
		echo "WE HAVE TOO MANY REJECTED ATTEMPTS AND NOW HAVE USED THE LAST WORKING EXAMPLE AND EXITING"
		exit
	fi
	factor=`unzip  -p  $fileAlign |grep "multiply all input standard deviations by factor" | awk '{ print $NF }'`;
	factor=`echo ${factor} | sed -e 's/[eE]+*/\\*10\\^/'`
	echo "factor word: " $factor
	if [[ $factor != "" ]];then
		export xresWorking=$xres; #Must be set before the new resolution is set which may cause too many rejects
		export yresWorking=$yres;
		echo "Factor word found! Resolution must increase by $factor."
		#for some reason the output of python will not overwrite the xres or yres? So must unset then set.
		xInput=$xres
		yInput=$yres
		unset xres;
		unset yres;
		xres=`python $pythonLocation/multiplyResolutionsByFactor.py $xInput / / $allPlanes / $factor`
		yres=`python $pythonLocation/multiplyResolutionsByFactor.py $yInput / / $allPlanes / $factor`
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	fi
	rejected=`unzip  -p  $fileAlign |grep "Too many rejects" |cut -d '-' -f2`; 
	echo "Rejects word:  $rejected "
	if [[ $rejected != "" ]];then #This makes sure that we do not cut too many tracks.
		export	numberRejectedAlignmentAttempts=$(($numberRejectedAlignmentAttempts+1))
		echo "Too many rejects. Resolution must increase by factor 10."
		xInput=$xres;
		yInput=$yres;
		unset xres;
		unset yres;
		xres=`python $pythonLocation/multiplyResolutionsByFactor.py $xInput /   / $allPlanes / 2` #TO DO: This script breaks if you provide no fixed planes
		yres=`python $pythonLocation/multiplyResolutionsByFactor.py $yInput /   / $allPlanes / 2`
		echo "New resolutions are for (X/Y):" $xres"/"$yres
	fi
	if [[ $numberRejectedAlignmentAttempts -eq 10 ]]
	then
		echo "We have already rejected $numberRejectedAlignmentAttempts times. However have found the factor $factor to improve alignment to continue."
	fi
	if [ $numberRejectedAlignmentAttempts -eq 10 ]   
	then
		echo "Too many rejects found set x/y resolution to last working alignment fit, run, then exit."
		export xres=$xresWorking;
		export yres=$yresWorking;
		tooManyRejectsExitLoopBool=true
	fi
	#If the loop passes after 5 attempts with a small number of reject use this. 
	#Otherwise we use the millepede esitmation to guide our fit.
	if [[ $numberRejectedAlignmentAttempts -lt 5 ]] && [[ $rejected == "" ]];then
		echo "Enough have passed the chi2 cut. Use this iteration to update gear parameters"
		break
	elif [[ $factor == "" ]] && [[ $rejected == "" ]];then
		echo "Mille chi2 is non existant. Here it is: $averageChi2Mille"
		echo "We can not find this or factor or rejects. Break alignment loop."
		break
	fi
done 
