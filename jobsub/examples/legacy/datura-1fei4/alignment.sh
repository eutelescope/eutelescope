#!/bin/sh

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

gear1="gear-$RUN-001.xml"
gear2="gear-$RUN-002.xml"
gear3="gear-$RUN-003.xml"
gear4="gear-$RUN-004.xml"
gear5="gear-$RUN-005.xml"
gear6="gear-$RUN-006.xml"
gear7="gear-$RUN-007.xml"
gear8="gear-$RUN-008.xml"
gear9="gear-$RUN-009.xml"
gear10="gear-$RUN-010.xml"

amode="7"; r="0.100";res="$r $r $r $r $r $r";
Fxr="0 1 2 3 4 5"
Fxs="0         5"
Fyr="0 1 2 3 4 5"
Fys="0         5"
Fzr="0          "
Fzs="0 1 2 3 4 5"

r="0.100";res="$r $r $r $r $r $r";prev="$r";
echo "prev:$prev and r:$r"
file="output/logs/aligngbl-0000${RUN}.zip"

if [ $# -ne 4 ]
then
 echo "$# parameters: $RUN $RUNLIST $file $gear10"
 exit
fi

#do="echo"
#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST                      -o GearAlignedFile="$gear1" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o Chi2Cut="5000"  -o pede="$pede" aligngbl $RUN
####
echo "file: $file"
multi=`unzip  -p  $file |grep "multiply all input standard deviations" |cut -d 'r' -f4`; 
multi=${multi/[eE]+/*10^+};
multi=${multi/[eE]-/*10^-};

echo "multi:$multi  prev: $prev"; 
if [[ -n $multi && -n $prev && $(echo "$prev > 0.004"|bc) -eq 1 ]];then
r=$(echo "scale=4;$prev*$multi"|bc);
prev=$r; 
res="$r $r $r  $r $r $r"
else
echo "multi:$multi  prev: $prev";
fi
echo "resolution $res"
#########################

#multi="0.33"; 

#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST                      -o GearAlignedFile="$gear1" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o Chi2Cut="5000"  -o pede="$pede" aligngbl $RUN
####
multi=`unzip  -p  $file |grep "multiply all input standard deviations" |cut -d 'r' -f4`; 
multi=${multi/[eE]-/*10^-};

echo "multi:$multi  prev: $prev"; 
if [[ -n $multi && -n $prev && $(echo "$prev > 0.004"|bc) -eq 1 ]];then
r=$(echo "scale=4;$prev*$multi"|bc);
prev=$r; 
res="$r $r $r $r $r $r"
else
echo "multi:$multi  prev: $prev";
fi
echo "resolution $res"
#########################

pede="chiscut  15. 5.";

#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o GearFile="$gear1" -o GearAlignedFile="$gear2" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o Chi2Cut="30"  -o pede="$pede" aligngbl $RUN
####
multi=`unzip  -p  $file |grep "multiply all input standard deviations" |cut -d 'r' -f4`; 
multi=${multi/[eE]-/*10^-};

echo "multi:$multi  prev: $prev"; 
if [[ -n $multi && -n $prev && $(echo "$prev > 0.004"|bc) -eq 1 ]];then
r=$(echo "scale=4;$prev*$multi"|bc);
prev=$r; 
res="$r $r $r $r $r $r"
else
echo "multi:$multi  prev: $prev";

fi
echo "resolution $res"
#########################



#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o GearFile="$gear2" -o GearAlignedFile="$gear3" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o Chi2Cut="30"  -o pede="$pede" aligngbl $RUN
####
multi=`unzip  -p  $file |grep "multiply all input standard deviations" |cut -d 'r' -f4`; 
multi=${multi/[eE]-/*10^-};

echo "multi:$multi  prev: $prev"; 
if [[ -n $multi && -n $prev && $(echo "$prev > 0.004"|bc) -eq 1 ]];then
r=$(echo "scale=4;$prev*$multi"|bc);
prev=$r; 
res="$r $r $r $r $r $r"
else
echo "multi:$multi  prev: $prev";
fi
echo "resolution $res"
#########################

# do rotations:

Fxr="0         5"
Fxs="0 1 2 3 4 5"
Fyr="0         5"
Fys="0 1 2 3 4 5"
Fzr="0 1 2 3 4 5"
Fzs="0 1 2 3 4 5"

#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o GearFile="$gear3" -o GearAlignedFile="$gear4" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}"  -o Chi2Cut="30" -o pede="$pede" 	aligngbl $RUN
####
multi=`unzip  -p  $file |grep "multiply all input standard deviations" |cut -d 'r' -f4`; 
multi=${multi/[eE]-/*10^-};

echo "multi:$multi  prev: $prev"; 
if [[ -n $multi && -n $prev && $(echo "$prev > 0.004"|bc) -eq 1 ]];then
r=$(echo "scale=4;$prev*$multi"|bc);
prev=$r; 
res="$r $r $r $r $r $r"
else
echo "multi:$multi  prev: $prev";
fi
echo "resolution $res"
#########################


#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o GearFile="$gear4" -o GearAlignedFile="$gear5" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}"  -o Chi2Cut="30" -o pede="$pede"	aligngbl $RUN
####
multi=`unzip  -p  $file |grep "multiply all input standard deviations" |cut -d 'r' -f4`; 
multi=${multi/[eE]-/*10^-};

echo "multi:$multi  prev: $prev"; 
if [[ -n $multi && -n $prev && $(echo "$prev > 0.004"|bc) -eq 1 ]];then
r=$(echo "scale=4;$prev*$multi"|bc);
prev=$r; 
res="$r $r $r $r $r $r"
else
echo "multi:$multi  prev: $prev";
fi
echo "resolution $res"
#########################

#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o GearFile="$gear5" -o GearAlignedFile="$gear6" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}"  -o Chi2Cut="30" -o pede="$pede"	aligngbl $RUN
####
multi=`unzip  -p  $file |grep "multiply all input standard deviations" |cut -d 'r' -f4`; 
multi=${multi/[eE]-/*10^-};

echo "multi:$multi  prev: $prev"; 
if [[ -n $multi && -n $prev && $(echo "$prev > 0.004"|bc) -eq 1 ]];then
r=$(echo "scale=4;$prev*$multi"|bc);
prev=$r; 
res="$r $r $r $r $r $r"
else
echo "multi:$multi  prev: $prev";

fi
echo "resolution $res"
#########################

# do z shifts
Fxr="0 1 2 3 4 5"
Fxs="0 1 2 3 4 5"
Fyr="0 1 2 3 4 5"
Fys="0 1 2 3 4 5"
Fzr="0 1 2 3 4 5"
Fzs="0          "

#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o GearFile="$gear6" -o GearAlignedFile="$gear7" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}"  -o Chi2Cut="30" -o pede="$pede"	aligngbl $RUN
####
multi=`unzip  -p  $file |grep "multiply all input standard deviations" |cut -d 'r' -f4`; 
multi=${multi/[eE]-/*10^-};

echo "multi:$multi  prev: $prev"; 
if [[ -n $multi && -n $prev && $(echo "$prev > 0.004"|bc) -eq 1 ]];then
r=$(echo "scale=4;$prev*$multi"|bc);
prev=$r; 
res="$r $r $r $r $r $r"
else
echo "multi:$multi  prev: $prev";
fi
echo "resolution $res"
#########################

#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o GearFile="$gear7" -o GearAlignedFile="$gear8" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}"  -o Chi2Cut="30" -o pede="$pede"	aligngbl $RUN
####
multi=`unzip  -p  $file |grep "multiply all input standard deviations" |cut -d 'r' -f4`; 
multi=${multi/[eE]-/*10^-};

echo "multi:$multi  prev: $prev"; 
if [[ -n $multi && -n $prev && $(echo "$prev > 0.004"|bc) -eq 1 ]];then
r=$(echo "scale=4;$prev*$multi"|bc);
prev=$r; 
res="$r $r $r $r $r $r"
else
echo "multi:$multi  prev: $prev";
fi
echo "resolution $res"
#########################

#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o GearFile="$gear8" -o GearAlignedFile="$gear9" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}"  -o Chi2Cut="30" -o pede="$pede"	aligngbl $RUN
####
multi=`unzip  -p  $file |grep "multiply all input standard deviations" |cut -d 'r' -f4`; 
multi=${multi/[eE]-/*10^-};

echo "multi:$multi  prev: $prev"; 
if [[ -n $multi && -n $prev && $(echo "$prev > 0.004"|bc) -eq 1 ]];then
r=$(echo "scale=4;$prev*$multi"|bc);
prev=$r; 
res="$r $r $r $r $r $r"
else
echo "multi:$multi  prev: $prev";
fi
echo "resolution $res"
#########################

#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o GearFile="$gear9" -o GearAlignedFile="$gear10" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}"  -o Chi2Cut="30" -o pede="$pede"	aligngbl $RUN
####
multi=`unzip  -p  $file |grep "multiply all input standard deviations" |cut -d 'r' -f4`; 
multi=${multi/[eE]-/*10^-};

echo "multi:$multi  prev: $prev"; 
if [[ -n $multi && -n $prev && $(echo "$prev > 0.004"|bc) -eq 1 ]];then
r=$(echo "scale=4;$prev*$multi"|bc);
prev=$r; 
res="$r $r $r $r $r $r"
else
echo "multi:$multi  prev: $prev";
fi
echo "resolution $res"
#########################

