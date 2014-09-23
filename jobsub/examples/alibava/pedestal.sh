#!/bin/sh

################################
#                              #
# ALiBaVa Analysis - Pedestal  #
#                              #
################################

# usage: sh pedestal.sh runnumber

#To calculate Pedestals
jobsub -c config.cfg -csv runlist.csv alibava-convert-ped $1
jobsub -c config.cfg alibava-pedestal $1
jobsub -c config.cfg -csv runlist.csv alibava-commonmode $1
jobsub -c config.cfg alibava-pedestal2 $1

