#!/bin/bash
while getopts i:o:n: option
do
        case "${option}"
        in
								i) inputGear=${OPTARG};;
								o) outputGearLast=${OPTARG};;
								n) number=${OPTARG};;
       esac
done

#THESE VARIABLES ARE ALSO CHANGED BELOW TO ALIGN WITH DIFFERENT SHIFTS AND ROTATIONS.
export Fxr="0 1 2 3 4 5" #This is the fixed planes for rotations round the x axis
export Fxs="0         5" #This is the fixed planes for shifts in the x axis
export Fyr="0 1 2 3 4 5" #This is the fixed planes for rotations round the y axis
export Fys="0         5" #This is the fixed planes for shifts in the y axis
export Fzr="0 1 2 3 4 5" #This is the fixed planes for rotations round the z axis
export Fzs="0 1 2 3 4 5" #This is the fixed planes for shifts in the z axis
#First input gear comes from input.
export outputGear="gear-XYShift-$number-${RUN}.xml"
export histoNameInput="GBLtrack-XYShift-$number-${RUN}"

#This is the loop to produce the new gear file for x/y shifts. 
./singleIterationAlignment.sh

export Fxr="0 1 2 3 4 5" #This is the fixed planes for rotations round the x axis
export Fxs="0         5" #This is the fixed planes for shifts in the x axis
export Fyr="0 1 2 3 4 5" #This is the fixed planes for rotations round the y axis
export Fys="0         5" #This is the fixed planes for shifts in the y axis
export Fzr="0" #This is the fixed planes for rotations round the z axis
export Fzs="0 1 2 3 4 5" #This is the fixed planes for shifts in the z axis
export inputGear="$outputGear"
export outputGear="gear-ZRotations-XYShifts-$number-${RUN}.xml"
export histoNameInput="GBLtrack-ZRotations-XYShifts-$number-${RUN}"
./singleIterationAlignment.sh
export Fxr="0 1 2 3 4 5" #This is the fixed planes for rotations round the x axis
export Fxs="0         5" #This is the fixed planes for shifts in the x axis
export Fyr="0 1 2 3 4 5" #This is the fixed planes for rotations round the y axis
export Fys="0         5" #This is the fixed planes for shifts in the y axis
export Fzr="0" #This is the fixed planes for rotations round the z axis
export Fzs="0  5" #This is the fixed planes for shifts in the z axis
export inputGear="$outputGear"
export outputGear="$outputGearLast"
export histoNameInput="GBLtrack-ZRotations-XYShifts-$number-${RUN}"
./singleIterationAlignment.sh

