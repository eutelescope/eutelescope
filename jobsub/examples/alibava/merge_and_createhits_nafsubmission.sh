#!/bin/zsh

################################
#                              #
# ALiBaVa Analysis - Pedestal  #
#                              #
################################

# usage: source merge_and_createhits_nafsubmission.sh runnumber(s)

# To merge alibava and telescope clusters and then create hits
# you will need to run these templates
TEMPLATELIST=('merger' 'hitmaker-local')

for TEMPLATE in $TEMPLATELIST; do
jobsub -c config.cfg -csv runlist.csv --subdir --naf qsubparameters.cfg $TEMPLATE $1

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

