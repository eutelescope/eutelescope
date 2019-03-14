#!/bin/sh
#export RUNRANGE="5037-5080"
#export RUNRANGE="5058-5080"
export RUNRANGE="5079"


#
../../jobsub.py -c config.cfg -csv runlist.csv converter  $RUNRANGE
../../jobsub.py -c config.cfg -csv runlist.csv clustering $RUNRANGE
../../jobsub.py -c config.cfg -csv runlist.csv hitmaker   $RUNRANGE
../../jobsub.py -c config.cfg -csv runlist.csv align      $RUNRANGE
../../jobsub.py -c config.cfg -csv runlist.csv fitter     $RUNRANGE
#
