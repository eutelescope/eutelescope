#!/bin/sh

##################################
#                                #
# ALiBaVa Analysis - Calibration #
#                                #
##################################

# usage: sh x_calibration.sh runnumber

jobsub.py -c config.cfg -csv runlist.csv converter $1
jobsub.py -c config.cfg -csv runlist.csv calibration $1

#
