#!/bin/sh

first="97"
last="97"
RUNLIST="runlist-150.csv"

#first="33"
#last="33"
#RUNLIST="runlist-20.csv"

#modus="straight"
#modus="daf"
modus="gbl"

#DRY=--dry-run

# GearAlignedFile = gear-@RunNumber@-aligned.xml 
# GearFile    	= @GearGeoFile@


#do="echo"

for i in `seq $first $last`; do

gear1="gear-$i-001.xml"
gear2="gear-$i-002.xml"
gear3="gear-$i-003.xml"
gear4="gear-$i-004.xml"
gear5="gear-$i-005.xml"
gear6="gear-$i-006.xml"
gear7="gear-$i-007.xml"
gear8="gear-$i-008.xml"
gear9="gear-$i-009.xml"
gear10="gear-$i-010.xml"
gear11="gear-$i-011.xml"

amode="7"; r="0.100";res="$r $r $r $r $r $r";
Fxr="0 1 2 3 4 5"
Fxs="0         5"
Fyr="0 1 2 3 4 5"
Fys="0         5"
Fzr="0          "
Fzs="0 1 2 3 4 5"

r="0.100";res="$r $r $r $r $r $r";

file="output/logs/aligngbl-0000${i}.zip"
#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST                      -o GearAlignedFile="$gear1" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o Chi2Cut="5000"  -o pede="$pede" aligngbl $i
####
multi=`unzip  -p  $file |grep "multiply all input standard deviations" |cut -d 'r' -f4`; 
echo "multi:" $multi; 
r=$(echo "$r*$multi"|bc); 
res="$r $r $r $r $r $r";
#########################

#multi="0.33"; 

#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST                      -o GearAlignedFile="$gear1" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o Chi2Cut="5000"  -o pede="$pede" aligngbl $i
####
multi=`unzip  -p  $file |grep "multiply all input standard deviations" |cut -d 'r' -f4`; 
echo "multi:" $multi; 
r=$(echo "$r*$multi"|bc); 
res="$r $r $r $r $r $r";
#########################

#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o GearFile="$gear1" -o GearAlignedFile="$gear2" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o Chi2Cut="5000"  -o pede="$pede" aligngbl $i
####
multi=`unzip  -p  $file |grep "multiply all input standard deviations" |cut -d 'r' -f4`; 
echo "multi:" $multi; 
r=$(echo "$r*$multi"|bc); 
res="$r $r $r $r $r $r";
#########################



#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o GearFile="$gear2" -o GearAlignedFile="$gear3" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}" -o Chi2Cut="5000"  -o pede="$pede" aligngbl $i
####
multi=`unzip  -p  $file |grep "multiply all input standard deviations" |cut -d 'r' -f4`; 
echo "multi:" $multi; 
r=$(echo "$r*$multi"|bc); 
res="$r $r $r $r $r $r";
#########################

# do rotations:
pede="chiscut  50. 10.";

Fxr="0         5"
Fxs="0 1 2 3 4 5"
Fyr="0         5"
Fys="0 1 2 3 4 5"
Fzr="0 1 2 3 4 5"
Fzs="0 1 2 3 4 5"

#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o GearFile="$gear3" -o GearAlignedFile="$gear4" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}"  -o Chi2Cut="5000" -o pede="$pede" 	aligngbl $i
####
multi=`unzip  -p  $file |grep "multiply all input standard deviations" |cut -d 'r' -f4`; 
echo "multi:" $multi; 
r=$(echo "$r*$multi"|bc); 
res="$r $r $r $r $r $r";
#########################


#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o GearFile="$gear4" -o GearAlignedFile="$gear5" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}"  -o Chi2Cut="5000" -o pede="$pede"	aligngbl $i
####
multi=`unzip  -p  $file |grep "multiply all input standard deviations" |cut -d 'r' -f4`; 
echo "multi:" $multi; 
r=$(echo "$r*$multi"|bc); 
res="$r $r $r $r $r $r";
#########################

#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o GearFile="$gear5" -o GearAlignedFile="$gear6" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}"  -o Chi2Cut="5000" -o pede="$pede"	aligngbl $i
####
multi=`unzip  -p  $file |grep "multiply all input standard deviations" |cut -d 'r' -f4`; 
echo "multi:" $multi; 
r=$(echo "$r*$multi"|bc); 
res="$r $r $r $r $r $r";
#########################

# do z shifts
Fxr="0 1 2 3 4 5"
Fxs="0 1 2 3 4 5"
Fyr="0 1 2 3 4 5"
Fys="0 1 2 3 4 5"
Fzr="0 1 2 3 4 5"
Fzs="0          "

#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o GearFile="$gear6" -o GearAlignedFile="$gear7" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}"  -o Chi2Cut="5000" -o pede="$pede"	aligngbl $i
####
multi=`unzip  -p  $file |grep "multiply all input standard deviations" |cut -d 'r' -f4`; 
echo "multi:" $multi; 
r=$(echo "$r*$multi"|bc); 
res="$r $r $r $r $r $r";
#########################

#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o GearFile="$gear7" -o GearAlignedFile="$gear8" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}"  -o Chi2Cut="5000" -o pede="$pede"	aligngbl $i
####
multi=`unzip  -p  $file |grep "multiply all input standard deviations" |cut -d 'r' -f4`; 
echo "multi:" $multi; 
r=$(echo "$r*$multi"|bc); 
res="$r $r $r $r $r $r";
#########################

#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o GearFile="$gear8" -o GearAlignedFile="$gear9" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}"  -o Chi2Cut="5000" -o pede="$pede"	aligngbl $i
####
multi=`unzip  -p  $file |grep "multiply all input standard deviations" |cut -d 'r' -f4`; 
echo "multi:" $multi; 
r=$(echo "$r*$multi"|bc); 
res="$r $r $r $r $r $r";
#########################

#########################
$do jobsub.py  $DRY -c config.cfg -csv $RUNLIST -o GearFile="$gear9" -o GearAlignedFile="$gear10" -o ResolutionPlane="$res" -o AlignmentMode="$amode"   -o FixXrot="${Fxr}" -o FixXshifts="${Fxs}"  -o FixYrot="${Fyr}" -o FixYshifts="${Fys}" -o FixZrot="${Fzr}" -o FixZshifts="${Fzs}"  -o Chi2Cut="5000" -o pede="$pede"	aligngbl $i
####
multi=`unzip  -p  $file |grep "multiply all input standard deviations" |cut -d 'r' -f4`; 
echo "multi:" $multi; 
r=$(echo "$r*$multi"|bc); 
res="$r $r $r $r $r $r";
#########################




done

#
