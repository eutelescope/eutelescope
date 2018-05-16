#!/bin/sh

################################
#                              #
# ALiBaVa Analysis - Data      #
#                              #
# DAF alignment                #
#                              #
################################

# usage: sh x_daf.sh runnumber

################################
# Run the analysis chain.
# x_pedestal.sh and x_telescope-<multi/single>.sh should be called before...
################################

jobsub.py -c config.cfg -csv runlist.csv converter $1
jobsub.py -c config.cfg -csv runlist.csv reco $1
jobsub.py -c config.cfg -csv runlist.csv clustering-1 $1
jobsub.py -c config.cfg -csv runlist.csv clustering-2 $1
jobsub.py -c config.cfg -csv runlist.csv datahisto $1
jobsub.py -c config.cfg -csv runlist.csv merge $1
jobsub.py -c config.cfg -csv runlist.csv hitmaker $1
jobsub.py -c config.cfg -csv runlist.csv coordinator $1

jobsub.py -c config.cfg -csv runlist.csv alignment-daf-1 $1
jobsub.py -c config.cfg -csv runlist.csv alignment-daf-2 $1
jobsub.py -c config.cfg -csv runlist.csv alignment-daf-3 $1
jobsub.py -c config.cfg -csv runlist.csv alignment-daf-4 $1
jobsub.py -c config.cfg -csv runlist.csv alignment-daf-5 $1
jobsub.py -c config.cfg -csv runlist.csv tracking-1 $1
jobsub.py -c config.cfg -csv runlist.csv alignment-daf-7 $1
jobsub.py -c config.cfg -csv runlist.csv tracking-2 $1
jobsub.py -c config.cfg -csv runlist.csv alignment-daf-9 $1
jobsub.py -c config.cfg -csv runlist.csv alignment-daf-10 $1
jobsub.py -c config.cfg -csv runlist.csv tracking-3 $1

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
