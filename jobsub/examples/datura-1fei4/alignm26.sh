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
                c) Chi2Cut=${OPTARG};;
       esac
done


    for x in {1..40}; do

gear[x]="gear-$RUN-$x.xml"
echo ${gear[x]}
    done    


MaxRecordNumber="1000"
AlignPlaneIds="0 1 2 20 3 4 5"
Planes="0 1 2 20 3 4 5"

#
amode="7";

#pede="chiscut 500. 50. "
Chi2Cut="500."


r="0.25000";
rfei4x="0.050"
rfei4y="0.250"
xres="$r $r $r $rfei4x $r $r $r";
yres="$r $r $r $rfei4y $r $r $r";

prev="$r";
echo "prev:$prev and r:$r"
file="output/logs/aligngbl-00${RUN}.zip"

if [ $# -ne 6 ]
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

ResidualsRMax="10.0"


#do="echo "
#DRY="--dry-run"

for x in {1..10}; do

$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o ResidualsRMax="$ResidualsRMax" -o MaxRecordNumber="$MaxRecordNumber" -o AlignPlaneIds="$AlignPlaneIds" -o Planes="$Planes"                         -o GearAlignedFile="${gear[1]}"  -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o Chi2Cut="$Chi2Cut"  -o pede="$pede" -o suffix="-XY-0" aligngbl $RUN

chi2string=`unzip  -p  $file |grep -i chi`
echo $chi2string

rejected=`unzip  -p  $file |grep "Too many rejects" |cut -d '-' -f2`; 
echo " rejected?$rejected;"

r=$(echo "scale=4;$prev*2.0"|bc);
prev=$r; 
xres="$r $r $r $rfei4x $r $r $r"
yres="$r $r $r $rfei4y $r $r $r"
echo "x=$x :: xres: $xres"
if [[ $rejected -eq "" ]];then
 echo " rejected?$rejected;;"
 echo "alignment is converging continue reducing sensor resolution "
 break;
fi

done

Fxr="0 1 2 20 3 4 5"
Fxs="0     20     5"
Fyr="0 1 2 20 3 4 5"
Fys="0     20     5"
Fzr="0     20      "
Fzs="0 1 2 20 3 4 5"


ResidualsRMax="0.5" 

echo "starting XY shifts/rotations"
#do=""
 for x in {1..10}; do
gear1=${gear[x]}
gear2=${gear[x+1]}
echo $gear1" to "$gear2
#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o ResidualsRMax="$ResidualsRMax" -o MaxRecordNumber="$MaxRecordNumber" -o AlignPlaneIds="$AlignPlaneIds" -o Planes="$Planes" -o GearFile="${gear1}"  -o GearAlignedFile="${gear2}"  -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o Chi2Cut="$Chi2Cut"  -o pede="$pede" -o suffix="-XY-$x" aligngbl $RUN

rejected=`unzip  -p  $file |grep "Too many rejects" |cut -d '-' -f2`; 
echo "at $x rejected?$rejected;;;"

STOP=`unzip  -p  $file |grep "STOP" |cut -d '-' -f6`; 
echo " STOP:1:$STOP:"

if [[ $STOP != "" ]]; then
 r=$(echo "scale=4;$prev*2.1"|bc);
 STOP="";
 prev=$r; 
 xres="$r $r $r $rfei4x $r $r $r"
 yres="$r $r $r $rfei4y $r $r $r"
 echo "STOP detected: taking new x=$x :: xres: $xres"
fi

if [[ $rejected = " stop"  ]]; then
 r=$(echo "scale=4;$prev*2.1"|bc);
 STOP="";
 prev=$r; 
 xres="$r $r $r $rfei4x $r $r $r"
 yres="$r $r $r $rfei4y $r $r $r"
 echo "rejected?$rejected?taking new x=$x :: xres: $xres"
else
 echo " rejected?$rejected: continue reducing sensor resolution"
 echo "alignment is converging continue reducing sensor resolution "
#
 if [[  $(echo "$prev > 0.010"|bc) -eq 1 ]];then
   r=$(echo "scale=4;$prev/2.0"|bc);
   prev=$r; 
   xres="$r $r $r $rfei4x $r $r $r"
   yres="$r $r $r $rfei4y $r $r $r"
   echo "not rejected: OK: x=$x :: xres : $xres  "
 fi
#
fi
#fi
#echo "resolution $res"
#########################
   done

exit 0;

#
Fxr="0     20     5"
Fxs="0 1 2 20 3 4 5"
Fyr="0     20     5"
Fys="0 1 2 20 3 4 5"
Fzr="0 1 2 20 3 4 5"
Fzs="0 1 2 20 3 4 5"


r="0.050";

#pede="chiscut 500. 50. "
Chi2Cut="500."


for x in {11..20}; do

gear1=${gear[x]}
gear2=${gear[x+1]}
echo ${gear1}" to "$gear2

#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o ResidualsRMax="$ResidualsRMax" -o MaxRecordNumber="$MaxRecordNumber"  -o AlignPlaneIds="$AlignPlaneIds" -o Planes="$Planes" -o GearFile="$gear1"  -o GearAlignedFile="$gear2"  -o xResolutionPlane="$xres" -o yResolutionPlane="$yres"  -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o Chi2Cut="$Chi2Cut"  -o pede="$pede" -o suffix="-tilts-$x" aligngbl $RUN
####
rejected=`unzip  -p  $file |grep "Too many rejects" |cut -d '-' -f2`; 
echo " rejected? $rejected"
STOP=`unzip  -p  $file |grep "STOP" |cut -d '-' -f6`; 
echo " STOP:$STOP:"

if [[ $STOP != "" ]]; then
 r=$(echo "scale=4;$prev*2.1"|bc);
 STOP="";
 prev=$r; 
 xres="$r $r $r $rfei4x $r $r $r"
 yres="$r $r $r $rfei4y $r $r $r"
 echo "STOP detected: taking new x=$x :: xres: $xres"
fi
if [[ $rejected = " stop" ]]; then
 r=$(echo "scale=4;$prev*2.1"|bc);
 STOP="";
 prev=$r; 
 xres="$r $r $r $rfei4x $r $r $r"
 yres="$r $r $r $rfei4y $r $r $r"
 echo "rejected:$rejected: taking new x=$x :: xres: $xres"
else
 echo " rejected? "$rejected
 echo "alignment is converging continue reducing sensor resolution "
#
 if [[  $(echo "$prev > 0.010"|bc) -eq 1 ]];then
   r=$(echo "scale=4;$prev/2.0"|bc);
   prev=$r; 
   xres="$r $r $r $rfei4x $r $r $r"
   yres="$r $r $r $rfei4y $r $r $r"
   echo "not rejected: OK: x=$x :: xres : $xres  "
 fi
 echo "keeping: not rejected: OK: x=$x :: xres : $xres  "

fi

done


exit 0;
# do not do z-alignment

#
#MaxRecordNumber="10000"
#
#do=" "
#
Fxr="0 1 2 20 3 4 5"
Fxs="0 1 2 20 3 4 5"
Fyr="0 1 2 20 3 4 5"
Fys="0 1 2 20 3 4 5"
Fzr="0 1 2 20 3 4 5"
Fzs="0     20      "

for x in {21..30}; do

gear1=${gear[x]}
gear2=${gear[x+1]}
echo ${gear1}" to "$gear2

#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o ResidualsRMax="$ResidualsRMax" -o MaxRecordNumber="$MaxRecordNumber"  -o AlignPlaneIds="$AlignPlaneIds" -o Planes="$Planes" -o GearFile="$gear1"  -o GearAlignedFile="$gear2"  -o xResolutionPlane="$xres" -o yResolutionPlane="$yres"  -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o Chi2Cut="$Chi2Cut"  -o pede="$pede" -o suffix="-z-$x" aligngbl $RUN
####

done



