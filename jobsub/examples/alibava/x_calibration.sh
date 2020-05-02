#!/bin/sh

##################################
#                                #
# ALiBaVa Analysis - Calibration #
#                                #
# Use for calibration data and   #
# header analysis                #
#                                #
# x_pedestal.sh should be called #
# before...                      #
#                                #
##################################

# usage: sh x_calibration.sh runnumber

jobsub.py -c config.cfg -csv runlist.csv converter $1
jobsub.py -c config.cfg -csv runlist.csv header $1
jobsub.py -c config.cfg -csv runlist.csv calibration $1

#
