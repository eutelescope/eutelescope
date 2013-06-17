#!/bin/sh

first="233"
last="5080"

#modus="straight"
#modus="daf"
modus="gbl"

for i in `seq $first $last`; do

../../jobsub.py -c config.cfg -csv runlist.csv converter  $i
../../jobsub.py -c config.cfg -csv runlist.csv clustering $i
../../jobsub.py -c config.cfg -csv runlist.csv hitmaker   $i

if [[ $modus == "straight" ]]; then
  ../../jobsub.py -c config.cfg -csv runlist.csv align      $i
  ../../jobsub.py -c config.cfg -csv runlist.csv fitter     $i

elif [[ $modus == "daf" ]]; then
  ../../jobsub.py -c config.cfg -csv runlist.csv aligndaf   $i
  ../../jobsub.py -c config.cfg -csv runlist.csv fitter     $i

elif [[ $modus == "gbl" ]]; then
  ../../jobsub.py -c config.cfg -csv runlist.csv tracksearch $i
  ../../jobsub.py -c config.cfg -csv runlist.csv aligngbl $i
  sh ../../../tools/parsepede/parsemilleout.sh output/database/run00$i-pede-steer.txt millepede.res output/database/run00$i-align.slcio
  ../../jobsub.py -c config.cfg -csv runlist.csv trackfit $i

fi
done

#
