#!/bin/sh

################################
#                              #
# ALiBaVa Analysis - Pedestal  #
#                              #
################################

# usage: sh reco.sh runnumber

# To get pedestal subtracted and common mode corrected signal values
jobsub -c config.cfg -csv runlist.csv alibava-converter $1
jobsub -c config.cfg -csv runlist.csv alibava-reco $1
jobsub -c config.cfg -csv runlist.csv alibava-datahisto $1
jobsub -c config.cfg alibava-commonmodecut $1
jobsub -c config.cfg -csv runlist.csv alibava-seedclustering $1
jobsub -c config.cfg -csv runlist.csv alibava-clusterhisto $1



