#!/bin/zsh

################################
#                              #
# ALiBaVa Analysis - Pedestal  #
#                              #
################################

# usage: source reco.sh runnumber(s)

# To get pedestal subtracted, common mode and cross talk corrected signal values
# you will need to run these templates
TEMPLATELIST=('alibava-converter' 'alibava-reco' 'alibava-datahisto' 'alibava-commonmodecut' 'alibava-seedclustering' 'alibava-crosstalk-it1' 'alibava-crosstalk-it2' 'alibava-crosstalk-it3')

for TEMPLATE in $TEMPLATELIST; do
jobsub -c config.cfg -csv runlist.csv $TEMPLATE $@
done


