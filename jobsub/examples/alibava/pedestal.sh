#!/bin/sh

################################
#                              #
# ALiBaVa Analysis - Pedestal  #
#                              #
################################

# usage: sh pedestal.sh runnumber

../../jobsub.py -c config.cfg convert-ped $1
../../jobsub.py -c config.cfg pedestal $1
../../jobsub.py -c config.cfg commonmode $1
../../jobsub.py -c config.cfg pedestal2 $1

#