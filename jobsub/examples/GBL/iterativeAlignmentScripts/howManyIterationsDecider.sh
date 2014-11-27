while getopts n: option
do
        case "${option}"
        in
								n) numberOfIterations=${OPTARG};;
       esac
done

if [  -z "$numberOfIterations" ]; then
	echo "The number of iterations is empty to howManyIterationsDecider bash script!"
	exit
fi

if [ $numberOfIterations -eq 1 ]
then 
	$scriptsLocation/patRecToAlignmentMultiLoop.sh -i "$inputGearInitial" -o "$outputGearFinal" -n 1
elif [ $numberOfIterations -eq 2 ]   
then
	#THIS WILL RUN THE ALIGNMENT PROCESS AS MANY TIME AS YOU LIKE TO IMPROVE ALIGNMENT
	$scriptsLocation/patRecToAlignmentMultiLoop.sh -i "$inputGearInitial" -o "gear-finished-iteration-${outputIdentifier}-1.xml" -n 1
	$scriptsLocation/patRecToAlignmentMultiLoop.sh -i "gear-finished-iteration-${outputIdentifier}-1.xml" -o "$outputGearFinal" -n 2
elif [ $numberOfIterations -eq 3 ]
then 
	$scriptsLocation/patRecToAlignmentMultiLoop.sh -i "$inputGearInitial" -o "gear-finished-iteration-${outputIdentifier}-1.xml" -n 1
	$scriptsLocation/patRecToAlignmentMultiLoop.sh -i "gear-finished-iteration-${outputIdentifier}-1.xml" -o "gear-finished-iteration-${outputIdentifier}-2.xml" -n 2
	$scriptsLocation/patRecToAlignmentMultiLoop.sh -i "gear-finished-iteration-${outputIdentifier}-2.xml" -o "$outputGearFinal" -n 3
else
 echo "Too Many iterations chosen!"
fi



