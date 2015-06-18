#!/bin/sh

################################
#                              #
# ALiBaVa Analysis - Pedestal  #
#                              #
################################

# usage: sh reco.sh runnumber

# To get pedestal subtracted and common mode corrected signal values
jobsub -c config/config.cfg -csv runlist/runlist.csv alibava-converter $1
jobsub -c config/config.cfg -csv runlist/runlist.csv alibava-reco $1
jobsub -c config/config.cfg -csv runlist/runlist.csv alibava-datahisto $1
jobsub -c config/config.cfg -csv runlist/runlist.csv alibava-commonmodecut $1
jobsub -c config/config.cfg -csv runlist/runlist.csv alibava-seedclustering $1
#jobsub -c config/config.cfg -csv runlist/runlist.csv alibava-clusterhisto $1



