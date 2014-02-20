#!/bin/bash

# parameters:
# 1 - runnumber         $RUN  
# 2 - runlist csv file: $RUNLIST

while getopts r:l: option
do
        case "${option}"
        in
                r) RUN=${OPTARG};;
                l) RUNLIST=${OPTARG};;
        esac
done

    for x in {1..40}; do

gear[x]="gear-$RUN-20-$x.xml"
echo ${gear[x]}
    done    

AlignPlaneIds="0 1 2 20 3 4 5"
Planes="0 1 2 20 3 4 5"
MaxRecordNumber="10000"

r="0.035"
rfei4="0.10"
ResolutionPlane="$r $r $r $rfei4 $r $r $r"
#
amode="7"; 


if [ $# -ne 4 ]
then
 echo "$# parameters: $RUN $RUNLIST $file $gear10"
 exit
fi

Chi2Cut="5000"

do="echo "

Fxr="0 1 2 20 3 4 5"
Fxs="0 1 2    3 4 5"
Fyr="0 1 2 20 3 4 5"
Fys="0 1 2    3 4 5"
Fzr="0 1 2    3 4 5"
Fzs="0 1 2 20 3 4 5"
DRY="--dry-run"

$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o MaxRecordNumber="$MaxRecordNumber" -o AlignPlaneIds="$AlignPlaneIds" -o Planes="$Planes" -o GearFile="gear-${RUN}-20.xml"   -o GearAlignedFile="${gear[1]}" -o ResolutionPlane="$ResolutionPlane" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o Chi2Cut="$Chi2Cut"  -o pede="$pede" aligngbl $RUN
# reduce Chi2Cut
Chi2Cut="30"
####


exit 0;

#  for x in `seq 10 30 2`; do
#echo $x
Fxr="0 1 2 20 3 4 5"
Fxs="0 1 2    3 4 5"
Fyr="0 1 2 20 3 4 5"
Fys="0 1 2    3 4 5"
Fzr="0 1 2    3 4 5"
Fzs="0 1 2 20 3 4 5"

echo "starting XY shifts/rotations"
#do=""
 for x in {1..10}; do
gear1=${gear[x]}
gear2=${gear[x+1]}
echo $gear1" to "$gear2
#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o MaxRecordNumber="$MaxRecordNumber" -o AlignPlaneIds="$AlignPlaneIds" -o Planes="$Planes" -o GearFile="${gear1}"  -o GearAlignedFile="${gear2}" -o ResolutionPlane="$ResolutionPlane" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o Chi2Cut="$Chi2Cut"  -o pede="$pede" aligngbl $RUN
#########################
   done


#
#do=""
#
Fxr="0 1 2    3 4 5"
Fxs="0 1 2 20 3 4 5"
Fyr="0 1 2    3 4 5"
Fys="0 1 2 20 3 4 5"
Fzr="0 1 2 20 3 4 5"
Fzs="0 1 2 20 3 4 5"

 for x in {11..20}; do
gear1=${gear[x]}
gear2=${gear[x+1]}
echo $gear1" to "$gear2
#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o MaxRecordNumber="$MaxRecordNumber" -o AlignPlaneIds="$AlignPlaneIds" -o Planes="$Planes" -o GearFile="${gear1}"  -o GearAlignedFile="${gear2}" -o ResolutionPlane="$ResolutionPlane" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o Chi2Cut="$Chi2Cut"  -o pede="$pede" aligngbl $RUN
#########################
   done



