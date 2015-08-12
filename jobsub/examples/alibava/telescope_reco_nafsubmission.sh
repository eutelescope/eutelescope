#!/bin/zsh

################################
#                              #
# ALiBaVa Analysis - Pedestal  #
#                              #
################################

# usage: source telescope_reco_nafsubmission.sh runnumber(s)

# To get telescope clusters to be formed and filtered
# you will need to run these templates
TEMPLATELIST=('telescope-converter' 'telescope-clustering' 'telescope-filter')

for TEMPLATE in $TEMPLATELIST; do
jobsub -c config.cfg -csv runlist.csv --subdir --naf qsubparameters.cfg $TEMPLATE $@

#check if there is any run submitted to the naf
#if there is wait for them to be finished 
while [[ -n `qstat` ]]; do
        echo Waiting for jobs to be finished
        sleep 30
done
if [[ -z `qstat` ]];then
        echo ready to start new job
fi

done

