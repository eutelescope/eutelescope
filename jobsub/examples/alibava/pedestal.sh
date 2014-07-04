#!/bin/sh

################################
#                              #
# ALiBaVa Analysis - Pedestal  #
#                              #
################################

# usage: sh pedestal.sh runnumber

To calculate Pedestals
../../jobsub.py -c config.cfg -csv runlist.csv convert-ped $1
../../jobsub.py -c config.cfg pedestal $1
../../jobsub.py -c config.cfg -csv runlist.csv commonmode $1
../../jobsub.py -c config.cfg pedestal2 $1

