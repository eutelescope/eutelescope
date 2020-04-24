#!/bin/sh

################################
#                              #
# ALiBaVa Analysis - Data      #
#                              #
# GBL alignment - n-type       #
#                              #
################################

# usage: sh x_gbl_n.sh runnumber

################################
# Run the analysis chain.
# x_pedestal.sh and x_telescope-<multi/single>.sh should be called before...
################################

jobsub.py -c config.cfg -csv runlist.csv converter $1
jobsub.py -c config.cfg -csv runlist.csv reco $1
jobsub.py -c config.cfg -csv runlist.csv rghfilter $1
jobsub.py -c config.cfg -csv runlist.csv clustering-1-afterrgh $1
jobsub.py -c config.cfg -csv runlist.csv clustering-2-afterrgh $1
jobsub.py -c config.cfg -csv runlist.csv datahisto $1
jobsub.py -c config.cfg -csv runlist.csv merge $1
jobsub.py -c config.cfg -csv runlist.csv hitmaker $1
jobsub.py -c config.cfg -csv runlist.csv coordinator $1

jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-1 $1
jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-2 $1
jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-3 $1
jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-4 $1
jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-5 $1
jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-6 $1
jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-7 $1
jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-8 $1
jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-9 $1
jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-10 $1
jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-10 $1
jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-11 $1
jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-12 $1
jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-13 $1
jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-14 $1
jobsub.py -c config.cfg -csv runlist.csv tracking-gbl $1

################################
# event-viewer displays tracks
################################
################################
# Uncomment this for event display server:
################################

#	# check if glced is running, if not start it
#	SERVICE='glced'
#	if ps ax | grep -v grep | grep $SERVICE > /dev/null
#	then
#		echo "$SERVICE is running!"
#	else
#		echo "$SERVICE is not running, will start it!"
#	$SERVICE &
#	fi
#
#	jobsub.py -c config.cfg -csv runlist.csv event-viewer $1

#
