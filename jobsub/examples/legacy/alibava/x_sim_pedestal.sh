#!/bin/sh

###################################
#                                 #
# ALiBaVa Analysis - Sim Pedestal #
#                                 #
###################################

# usage: sh x_sim_pedestal.sh runnumber

jobsub.py -c config.cfg -csv runlist.csv simconverter $1
jobsub.py -c config.cfg -csv runlist.csv pedestal $1
jobsub.py -c config.cfg -csv runlist.csv commonmode $1
jobsub.py -c config.cfg -csv runlist.csv pedestal2 $1
jobsub.py -c config.cfg -csv runlist.csv pedestalhisto $1

#
