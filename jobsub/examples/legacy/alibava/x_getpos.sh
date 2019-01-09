#!/bin/sh

########################################
#                                      #
# ALiBaVa Analysis - Get DUT Position  #
#                                      #
########################################

# usage: sh x_getpos.sh runnumber (4-digit)

jobsub.py -c config.cfg -csv runlist.csv createdummy $1
jobsub.py -c config.cfg -csv runlist.csv getdutposition $1
cat output/histograms/00$1-finalposition.txt

#
