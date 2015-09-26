#!/bin/zsh

################################
#                              #
# ALiBaVa Analysis - Pedestal  #
#                              #
################################

# usage: source merge_and_createhits.sh runnumber(s)

# To merge alibava and telescope clusters and then create hits
# you will need to run these templates
TEMPLATELIST=('merger' 'hitmaker-local')

for TEMPLATE in $TEMPLATELIST; do
jobsub -c config/config.cfg -csv runlistfiles/runlist.csv $TEMPLATE $1
done


