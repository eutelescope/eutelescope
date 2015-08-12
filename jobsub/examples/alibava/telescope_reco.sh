#!/bin/zsh

################################
#                              #
# ALiBaVa Analysis - Pedestal  #
#                              #
################################

# usage: source telescope_reco.sh runnumber(s)

# To get telescope clusters to be formed and filtered
# you will need to run these templates
TEMPLATELIST=('telescope-converter' 'telescope-clustering' 'telescope-filter')

for TEMPLATE in $TEMPLATELIST; do
jobsub -c config.cfg -csv runlist.csv $TEMPLATE $@
done


