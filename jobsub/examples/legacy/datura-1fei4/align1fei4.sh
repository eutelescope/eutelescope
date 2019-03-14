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
  gear[x]="gear-$RUN-20-$x.xml"
  echo ${gear[x]}
done    

AlignPlaneIds="0 1 2 20 3 4 5"
Planes="0 1 2 20 3 4 5"
MaxRecordNumber="1000"

r="0.0035"
xrfei4="1.000"
yrfei4="4.000"

xres="$r $r $r $xrfei4 $r $r $r"
yres="$r $r $r $yrfei4 $r $r $r"

xprev="$xrfei4"
yprev="$yrfei4"

MaxMissingHitsPerTrack="0"
ResidualsRMax="10.0"

#
amode="7"; 
file="output/logs/aligngbl-00${RUN}.zip"


if [ $# -ne 6 ]
then
 echo "$# parameters: $RUN $RUNLIST $file $gear10"
 exit
fi


#do="echo "

Fxr="0 1 2 20 3 4 5"
Fxs="0 1 2    3 4 5"
Fyr="0 1 2 20 3 4 5"
Fys="0 1 2    3 4 5"
Fzr="0 1 2    3 4 5"
Fzs="0 1 2 20 3 4 5"
#DRY="--dry-run"

for x in {1..10}; do

$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o ResidualsRMax="$ResidualsRMax" -o MaxRecordNumber="$MaxRecordNumber" -o AlignPlaneIds="$AlignPlaneIds" -o Planes="$Planes"  -o GearFile="gear-$RUN-10.xml"                        -o GearAlignedFile="${gear[1]}"  -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o Chi2Cut="$Chi2Cut"  -o pede="$pede" -o suffix="-1fei4-XY-0" aligngbl $RUN

chi2string=`unzip  -p  $file |grep -i chi`
echo $chi2string

rejected=`unzip  -p  $file |grep "Too many rejects" |cut -d '-' -f2`; 
echo "$x rejected?$rejected;"

xrfei4=$(echo "scale=4;$xprev*2.0"|bc);
echo "new xrfei4: $xrfei4";
yrfei4=$(echo "scale=4;$yprev*2.0"|bc);
echo "new yrfei4: $yrfei4";

xprev=$xrfei4
yprev=$yrfei4
xres="$r $r $r $xrfei4 $r $r $r"
yres="$r $r $r $yrfei4 $r $r $r"
echo "x=$x :: xres: $xres"
if [[ $rejected = "" ]]; 
then
 echo "$x rejected?$rejected;;"
 echo "alignment is converging continue reducing sensor resolution "
 break;
fi

done


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
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o MaxRecordNumber="$MaxRecordNumber" -o AlignPlaneIds="$AlignPlaneIds" -o Planes="$Planes" -o GearFile="${gear1}"  -o GearAlignedFile="${gear2}"  -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o Chi2Cut="$Chi2Cut"  -o pede="$pede" -o MaxMissingHitsPerTrack="$MaxMissingHitsPerTrack"  -o ResidualsRMax="$ResidualsRMax"  -o suffix="-1fei4-XY-$x" aligngbl $RUN
#########################
# echo "file: $file"
#multi=`unzip  -p  $file |grep "multiply all input standard deviations" |cut -d 'r' -f4`; 
#multi=${multi/[eE]+/*10^+};
#multi=${multi/[eE]-/*10^-};

#echo "multi:$multi  prev: $xprev"; 
#if [[ -n $multi && -n $xprev && $(echo "$xprev > 0.001"|bc) -eq 1 ]];then

if [[  $(echo "$xprev > 0.020"|bc) -eq 1 ]];then
xrfei4=$(echo "scale=4;$xprev/1.5"|bc);
xprev=$xrfei4; 

yrfei4=$(echo "scale=4;$yprev/1.5"|bc);
yprev=$yrfei4; 

xres="$r $r $r $xrfei4 $r $r $r"
yres="$r $r $r $yrfei4 $r $r $r"

echo "xres: "$xres
echo "yres: "$yres
fi

#else
#echo "multi:$multi   xprev: $xprev   yprev: $yprev";
#fi
#echo "resolution $res"
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
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o MaxRecordNumber="$MaxRecordNumber" -o AlignPlaneIds="$AlignPlaneIds" -o Planes="$Planes" -o GearFile="${gear1}"  -o GearAlignedFile="${gear2}"  -o xResolutionPlane="$xres" -o yResolutionPlane="$yres" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o Chi2Cut="$Chi2Cut"  -o pede="$pede" -o MaxMissingHitsPerTrack="$MaxMissingHitsPerTrack"  -o ResidualsRMax="$ResidualsRMax"  -o suffix="-1fei4-tilts-$x" aligngbl $RUN
#########################
   done



