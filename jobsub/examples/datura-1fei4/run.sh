#!/bin/sh

first="3546"
#last="3546"

#first="3552"
first="3568"
last="3568"
#last="3569"
#last="3552"
RUNLIST="runlist-1fei4.csv"

#first="33"
#last="33"
#RUNLIST="runlist-20.csv"

#modus="straight"
#modus="daf"
modus="gbl"

#DRY=--dry-run

for i in `seq $first $last`; do

# jobsub.py $DRY -c config.cfg -csv $RUNLIST converter  $i
# jobsub.py  $DRY -c config.cfg -csv $RUNLIST clustering $i
##### jobsub.py  $DRY -c config.cfg -csv $RUNLIST filter $i

if [[ $modus == "straight" ]]; then
 jobsub.py $DRY -c config.cfg -csv $RUNLIST hitmaker   $i
# alignment using straight line assumption
 jobsub.py $DRY -c config.cfg -csv $RUNLIST align      $i
# fitter using broken line implementation by F.Zarnezki
 jobsub.py $DRY -c config.cfg -csv $RUNLIST trackTestFitter $i

elif [[ $modus == "daf" ]]; then
 jobsub.py $DRY -c config.cfg -csv $RUNLIST hitmaker   $i
 jobsub.py $DRY -c config.cfg -csv $RUNLIST aligndaf   $i
 jobsub.py $DRY -c config.cfg -csv $RUNLIST trackdaf   $i

elif [[ $modus == "gbl" ]]; then
# FOR NEW gbl ONE NEEDS TO GET HITS IN LOCAL COORDINATE SYSTEM:
# jobsub.py  $DRY -c config.cfg -csv $RUNLIST hitlocal $i
# Exhautsive and Helix track search results are not identical - to be investigated (perhaps trivially explained)

#### jobsub.py  $DRY -c config.cfg -csv $RUNLIST tracksearchExh $i
# jobsub.py  $DR -o MaxMissingHitsPerTrack="1"  -c config.cfg -csv $RUNLIST tracksearchHelix $i
# echo jobsub.py  $DRY -c config.cfg -csv $RUNLIST aligngbl $i
#
# bash alignm26.sh -r ${i} -l ${RUNLIST} -c 5000
# bash align1fei4.sh -r ${i} -l ${RUNLIST} 

#### jobsub.py  $DRY -c config.cfg -csv $RUNLIST aligngbl $i
#
# 10 iterations to aligne 6 planes (6D)
#   jobsub.py  $DRY -c config.cfg -o GearAlignedFile="gear-${i}-20-21.xml" -csv $RUNLIST trackgbl $i
 jobsub.py  $DRY -c config.cfg -o GearAlignedFile="gear-${i}-16.xml" -o Chi2Cut="5000" -csv $RUNLIST trackgbl $i

fi
done

#
