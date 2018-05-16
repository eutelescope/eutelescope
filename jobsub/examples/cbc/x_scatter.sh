#!/bin/sh

####################################
#                                  #
# CMS CBC Analysis - X0 Scattering #
#                                  #
# GBL alignment                    #
#                                  #
####################################

# usage: sh x_scatter.sh runnumber

jobsub.py -c config.cfg -csv runlist.csv telescope-converter $1
jobsub.py -c config.cfg -csv runlist.csv scatter-clustering $1
jobsub.py -c config.cfg -csv runlist.csv scatter-hitmaker $1
jobsub.py -c config.cfg -csv runlist.csv scatter-alignment-gbl-1 $1
jobsub.py -c config.cfg -csv runlist.csv scatter-alignment-gbl-2 $1
jobsub.py -c config.cfg -csv runlist.csv scatter-alignment-gbl-3 $1
jobsub.py -c config.cfg -csv runlist.csv scatter-alignment-gbl-4 $1
jobsub.py -c config.cfg -csv runlist.csv scatter-alignment-gbl-5 $1
jobsub.py -c config.cfg -csv runlist.csv scatter-alignment-gbl-6 $1
jobsub.py -c config.cfg -csv runlist.csv scatter-alignment-gbl-7 $1
jobsub.py -c config.cfg -csv runlist.csv scatter-tracking-gbl $1

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
