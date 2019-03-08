#!/bin/sh

################################
#                              #
# ALiBaVa Analysis - Pedestal  #
#                              #
################################

# usage: sh x_pedestal.sh runnumber

jobsub.py -c config.cfg -csv runlist.csv converter $1
jobsub.py -c config.cfg -csv runlist.csv pedestal $1
jobsub.py -c config.cfg -csv runlist.csv commonmode $1
jobsub.py -c config.cfg -csv runlist.csv pedestal2 $1
jobsub.py -c config.cfg -csv runlist.csv pedestalhisto $1

#
