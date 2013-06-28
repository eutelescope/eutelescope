#!/bin/sh

first="97"
last="97"

modus="straight"
#modus="daf"
#modus="gbl"


for i in `seq $first $last`; do

jobsub.py -c config.cfg -csv runlist.csv converter  $i
jobsub.py -c config.cfg -csv runlist.csv clustering $i
jobsub.py -c config.cfg -csv runlist.csv hitmaker   $i

if [[ $modus == "straight" ]]; then
  jobsub.py -c config.cfg -csv runlist.csv align      $i
  jobsub.py -c config.cfg -csv runlist.csv fitter     $i

elif [[ $modus == "daf" ]]; then
  jobsub.py -c config.cfg -csv runlist.csv aligndaf   $i
  jobsub.py -c config.cfg -csv runlist.csv fitter     $i

elif [[ $modus == "gbl" ]]; then
  jobsub.py -c config.cfg -csv runlist.csv tracksearch $i
  jobsub.py -c config.cfg -csv runlist.csv aligngbl $i
  jobsub.py -c config.cfg -csv runlist.csv trackfit $i
  jobsub.py -c config.cfg -csv runlist.csv fitter $i

fi
done

#
