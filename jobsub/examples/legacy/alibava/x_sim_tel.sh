#!/bin/sh

####################################
#                                  #
# ALiBaVa Analysis - Telescope Sim #
#                                  #
####################################

# usage: sh x_sim_tel.sh runnumber

jobsub.py -c config.cfg -csv runlist.csv telescope-clustering-sim $1
jobsub.py -c config.cfg -csv runlist.csv telescope-filter $1

#
