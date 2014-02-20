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

gear[x]="gear-$RUN-$x.xml"
echo ${gear[x]}
    done    


MaxRecordNumber="10000"
AlignPlaneIds="0 1 2 20 3 4 5"
Planes="0 1 2 20 3 4 5"

#
amode="7";

r="0.0500";
rfei4="10.001";
res="$r $r $r $rfei4 $r $r $r";
prev="$r";
echo "prev:$prev and r:$r"
file="output/logs/aligngbl-00${RUN}.zip"

if [ $# -ne 4 ]
then
 echo "$# parameters: $RUN $RUNLIST $file $gear10"
 exit
fi


#  for x in `seq 10 30 2`; do
#echo $x
Fxr="0 1 2 20 3 4 5"
Fxs="0     20     5"
Fyr="0 1 2 20 3 4 5"
Fys="0     20     5"
Fzr="0     20      "
Fzs="0 1 2 20 3 4 5"



Chi2Cut="5000"

#do="echo "
#DRY="--dry-run"


$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o MaxRecordNumber="$MaxRecordNumber" -o AlignPlaneIds="$AlignPlaneIds" -o Planes="$Planes"                         -o GearAlignedFile="${gear[1]}" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o Chi2Cut="$Chi2Cut"  -o pede="$pede" aligngbl $RUN
# reduce Chi2Cut
Chi2Cut="30"
####

echo "starting XY shifts/rotations"
#do=""
 for x in {1..6}; do
gear1=${gear[x]}
gear2=${gear[x+1]}
echo $gear1" to "$gear2
#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o MaxRecordNumber="$MaxRecordNumber" -o AlignPlaneIds="$AlignPlaneIds" -o Planes="$Planes" -o GearFile="${gear1}"  -o GearAlignedFile="${gear2}" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o Chi2Cut="$Chi2Cut"  -o pede="$pede" aligngbl $RUN
####
echo "file: $file"
multi=`unzip  -p  $file |grep "multiply all input standard deviations" |cut -d 'r' -f4`; 
multi=${multi/[eE]+/*10^+};
multi=${multi/[eE]-/*10^-};

echo "multi:$multi  prev: $prev"; 
if [[ -n $multi && -n $prev && $(echo "$prev > 0.001"|bc) -eq 1 ]];then
r=$(echo "scale=4;$prev*$multi"|bc);
prev=$r; 
res="$r $r $r 10.000 $r $r $r"
else
echo "multi:$multi  prev: $prev";
fi
#echo "resolution $res"
#########################
   done

#
MaxRecordNumber="1000"
#
#do=" "
#
Fxr="0     20     5"
Fxs="0 1 2 20 3 4 5"
Fyr="0     20     5"
Fys="0 1 2 20 3 4 5"
Fzr="0 1 2 20 3 4 5"
Fzs="0 1 2 20 3 4 5"

for x in {7..10}; do

gear1=${gear[x]}
gear2=${gear[x+1]}
echo ${gear1}" to "$gear2

#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o MaxRecordNumber="$MaxRecordNumber"  -o AlignPlaneIds="$AlignPlaneIds" -o Planes="$Planes" -o GearFile="$gear1"  -o GearAlignedFile="$gear2" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o Chi2Cut="$Chi2Cut"  -o pede="$pede" aligngbl $RUN
####

done



