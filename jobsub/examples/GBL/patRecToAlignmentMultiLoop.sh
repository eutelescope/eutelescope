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
#THESE VARIABLES ARE ALSO CHANGED BELOW TO ALIGN WITH DIFFERENT SHIFTS AND ROTATIONS.
export Fxr="0 1 2 3 4 5" #This is the fixed planes for rotations round the x axis
export Fxs="0         5" #This is the fixed planes for shifts in the x axis
export Fyr="0 1 2 3 4 5" #This is the fixed planes for rotations round the y axis
export Fys="0         5" #This is the fixed planes for shifts in the y axis
export Fzr="0 1 2 3 4 5" #This is the fixed planes for rotations round the z axis
export Fzs="0 1 2 3 4 5" #This is the fixed planes for shifts in the z axis
#First input gear comes from the input of the bash script. IMPORTANT: If you edit this file then the first alignment must use the inital gear.
export inputGear="$inputGear"
export outputGear="gear-XYShift-Iteration-$number-${RUN}.xml"
export histoNameInput="GBLtrack-XYShift-Iteration-$number-${RUN}" #Note we must add another bash script which will get determine for each step the tracks for r=0.005 mm

#This is the loop to produce the new gear file for x/y shifts. 
./patRecToAlignmentSingleLoop.sh
./patRecAndTrackFit.sh -i "$outputGear" -h "$histoNameInput"  
export Fxr="0 1 2 3 4 5" #This is the fixed planes for rotations round the x axis
export Fxs="0         5" #This is the fixed planes for shifts in the x axis
export Fyr="0 1 2 3 4 5" #This is the fixed planes for rotations round the y axis
export Fys="0         5" #This is the fixed planes for shifts in the y axis
export Fzr="0" #This is the fixed planes for rotations round the z axis
export Fzs="0 1 2 3 4 5" #This is the fixed planes for shifts in the z axis
export inputGear="$outputGear" #We now use the gear file produced from the last iteration.
#export outputGear="gear-ZRotations-XYShifts-Iteration-$number-${RUN}.xml"
export outputGear="$outputGearLast" #IMPORTANT: The last alignment must output the new gear to the output gear name specified. 
export histoNameInput="GBLtrack-ZRotations-XYShifts-Iteration-$number-${RUN}"
./patRecToAlignmentSingleLoop.sh
./patRecAndTrackFit.sh -i "$outputGear" -h "$histoNameInput"  
#export Fxr="0 1 2 3 4 5" #This is the fixed planes for rotations round the x axis
#export Fxs="0         5" #This is the fixed planes for shifts in the x axis
#export Fyr="0 1 2 3 4 5" #This is the fixed planes for rotations round the y axis
#export Fys="0         5" #This is the fixed planes for shifts in the y axis
#export Fzr="0" #This is the fixed planes for rotations round the z axis
#export Fzs="0  5" #This is the fixed planes for shifts in the z axis
#export inputGear="$outputGear"
#export outputGear="$outputGearLast" #IMPORTANT: The last alignment must output the new gear to the output gear name specified. 
#export histoNameInput="GBLtrack-ZRotations-XYZShifts-Iteration-$number-${RUN}"
#./patRecToAlignmentSingleLoop.sh
#./patRecAndTrackFit.sh -i "$outputGear" -h "$histoNameInput"  

