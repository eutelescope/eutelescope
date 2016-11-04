#!/bin/sh

################################
#                              #
# ALiBaVa Analysis - Telescope #
#                              #
# for single telescope runs    #
# with reference plane         #
#                              #
################################

# usage: sh x_telescope-single-ref.sh runnumber

jobsub.py -c config.cfg -csv runlist.csv telescope-converter $1
jobsub.py -c config.cfg -csv runlist.csv telescope-clustering-ref $1
jobsub.py -c config.cfg -csv runlist.csv telescope-filter $1

#

