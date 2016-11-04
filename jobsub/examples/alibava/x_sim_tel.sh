#!/bin/sh

################################
#                              #
# ALiBaVa Analysis - Telescope #
#                              #
# for single telescope runs    #
#                              #
################################

# usage: sh x_telescope-single.sh runnumber

jobsub.py -c config.cfg -csv runlist.csv telescope-clustering-sim $1
jobsub.py -c config.cfg -csv runlist.csv telescope-filter $1

#
