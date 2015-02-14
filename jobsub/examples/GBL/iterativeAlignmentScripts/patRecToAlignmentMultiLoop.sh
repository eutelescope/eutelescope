#!/bin/bash
#This file will run all the different alignment modes. So we first run X/Y shifts and then ZRotations. If you have trouble with alignment you can add other step to this file. However you must make sure that the first gear input is the#intial gear and the last is the final gear. 
while getopts i:o:n: option
do
        case "${option}"
        in
								i) inputGear=${OPTARG};;
								o) outputGearLast=${OPTARG};;
								n) number=${OPTARG};; #Number is used to specify what iteration number we are on. 
       esac
done
echo "Entering MultiLoop on iteration $number"
echo "These are the input parameters in multi loop."
echo "Input Gear: $inputGear"
echo "Output Gear: $outputGearLast"
echo "Iteration: $number"

if [  -z "$inputGear" ]; then
	echo "The input gear is empty to multiLoop bash script!"
	exit
fi
if [  -z "$outputGearLast" ]; then
	echo "The output gear is empty to multiLoop bash script!"
	exit
fi
if [  -z "$number" ]; then
	echo "The number gear is empty to multiLoop bash script!"
	exit
fi
export number="$number" #We export this to echo in single loop bash script.
#THESE VARIABLES ARE ALSO CHANGED BELOW TO ALIGN WITH DIFFERENT SHIFTS AND ROTATIONS. DOES NOT HAVE TO BE IN Z-ORDERING.
export Fxr="$allPlanes" #This is the fixed planes for rotations round the x axis. Note we use dutPlanes since rotation round x axis always fixed parameter.
export Fxs="$allPlanesFixed" #This is the fixed planes for shifts in the x axis. Note we use exclude planes here so we can sometimes let dut alignment parameters free
export Fyr="$allPlanes" #This is the fixed planes for rotations round the y axis
export Fys="$allPlanesFixed" #This is the fixed planes for shifts in the y axis
export Fzr="$allPlanes" #This is the fixed planes for rotations round the z axis
export Fzs="$allPlanes" #This is the fixed planes for shifts in the z axis
#First input gear comes from the input of the bash script. IMPORTANT: If you edit this file then the first alignment must use the inital gear.
export inputGear="$inputGear"
export outputGear="gear-XYShift-${outputIdentifier}-Iteration-$number-${RUN}.xml"
export histoNameInput="Alignment-XYShift-${outputIdentifier}-Iteration-$number-${RUN}" 
#This is the loop to produce the new gear file for x/y shifts. 
$scriptsLocation/patRecToAlignmentSingleLoop.sh
$scriptsLocation/patRecAndTrackFit.sh -i "$outputGear" -h "$histoNameInput"  
echo "We have produced new histogram $histoNameInput"
export Fxr="$allPlanes" #This is the fixed planes for rotations round the x axis
export Fxs="$allPlanesFixed" #This is the fixed planes for shifts in the x axis
export Fyr="$allPlanes" #This is the fixed planes for rotations round the y axis
export Fys="$allPlanesFixed" #This is the fixed planes for shifts in the y axis
export Fzr="$allPlanesFixed" #This is the fixed planes for rotations round the z axis
export Fzs="$allPlanes" #This is the fixed planes for shifts in the z axis
export inputGear="$outputGear" #We now use the gear file produced from the last iteration.
export outputGear="gear-ZRotations-XYShift-${outputIdentifier}-Iteration-$number-${RUN}.xml"  
export histoNameInput="Alignment-ZRotations-XYShifts-${outputIdentifier}-Iteration-$number-${RUN}"
$scriptsLocation/patRecToAlignmentSingleLoop.sh
$scriptsLocation/patRecAndTrackFit.sh -i "$outputGear" -h "$histoNameInput"  

##We keep all fixed and free z-shifts. Otherwise we find that we can not find a sensible solution to the z shifts.
#echo "We have produced new histogram $histoNameInput"
#export Fxr="0 1 2 3 4 5 $dutPlanes" #This is the fixed planes for rotations round the x axis
#export Fxs="0 1 2 3 4 5 $dutPlanes" #This is the fixed planes for shifts in the x axis
#export Fyr="0 1 2 3 4 5 $dutPlanes" #This is the fixed planes for rotations round the y axis
#export Fys="0 1 2 3 4 5 $dutPlanes" #This is the fixed planes for shifts in the y axis
#export Fzr="0 1 2 3 4 5 $dutPlanes" #This is the fixed planes for rotations round the z axis
#export Fzs="$allPlanesFixed" #This is the fixed planes for shifts in the z axis
#export inputGear="$outputGear"
#export outputGear="gear-ZRotations-XYZShift-${outputIdentifier}-Iteration-$number-${RUN}.xml" 
#export histoNameInput="Alignment-ZRotations-XYZShifts-${outputIdentifier}-Iteration-$number-${RUN}"
#$scriptsLocation/patRecToAlignmentSingleLoop.sh
#$scriptsLocation/patRecAndTrackFit.sh -i "$outputGear" -h "$histoNameInput"  

##We keep all fixed and free x/y-rotations. 
#echo "We have produced new histogram $histoNameInput"
#export Fxr=" 0 1 2 3 4 5 $dutPlanes" #This is the fixed planes for rotations round the x axis
#export Fxs="0 1 2 3 4 5 $dutPlanes" #This is the fixed planes for shifts in the x axis
#export Fyr="$allPlanesFixed" #This is the fixed planes for rotations round the y axis
#export Fys="0 1 2 3 4 5 $dutPlanes" #This is the fixed planes for shifts in the y axis
#export Fzr="0 1 2 3 4 5 $dutPlanes" #This is the fixed planes for rotations round the z axis
#export Fzs="0 1 2 3 4 5 $dutPlanes" #This is the fixed planes for shifts in the z axis
#export inputGear="$outputGear"
#export outputGear="gear-YZRotations-XYZShift-${outputIdentifier}-Iteration-$number-${RUN}.xml" #IMPORTANT: The last alignment must output the new gear to the output gear name specified. 
#export histoNameInput="Alignment-YZRotations-XYZShifts-${outputIdentifier}-Iteration-$number-${RUN}"
#$scriptsLocation/patRecToAlignmentSingleLoop.sh
#$scriptsLocation/patRecAndTrackFit.sh -i "$outputGear" -h "$histoNameInput"  
#
##We keep all fixed and free x/y-rotations. 
#echo "We have produced new histogram $histoNameInput"
#export Fxr="$allPlanesFixed" #This is the fixed planes for rotations round the x axis
#export Fxs="0 1 2 3 4 5 $dutPlanes" #This is the fixed planes for shifts in the x axis
#export Fyr="0 1 2 3 4 5 $dutPlanes" #This is the fixed planes for rotations round the y axis
#export Fys="0 1 2 3 4 5 $dutPlanes" #This is the fixed planes for shifts in the y axis
#export Fzr="0 1 2 3 4 5 $dutPlanes" #This is the fixed planes for rotations round the z axis
#export Fzs="0 1 2 3 4 5 $dutPlanes" #This is the fixed planes for shifts in the z axis
#export inputGear="$outputGear"
#export outputGear="$outputGearLast" #IMPORTANT: The last alignment must output the new gear to the output gear name specified. 
#export histoNameInput="Alignment-XYZRotations-XYZShifts-${outputIdentifier}-Iteration-$number-${RUN}"
#$scriptsLocation/patRecToAlignmentSingleLoop.sh
#$scriptsLocation/patRecAndTrackFit.sh -i "$outputGear" -h "$histoNameInput"  
#
