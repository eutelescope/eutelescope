#!/bin/sh

#################################
#                               #
# CMS MPA Analysis - EUDAQ Data #
#                               #
# GBL alignment                 #
#                               #
#################################

# usage: sh x_run.sh runnumber

jobsub.py -c config.cfg -csv runlist.csv converter $1
jobsub.py -c config.cfg -csv runlist.csv clustering $1
jobsub.py -c config.cfg -csv runlist.csv filter $1
jobsub.py -c config.cfg -csv runlist.csv hitmaker $1
jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-1 $1
jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-2 $1
jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-3 $1