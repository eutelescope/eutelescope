#!/bin/zsh

################################
#                              #
# ALiBaVa Analysis - Pedestal  #
#                              #
################################

# usage: source pedestal_nafsubmission.sh runnumber(s)

#To calculate Pedestals
# you will need to run these templates
TEMPLATELIST=('alibava-convert-ped' 'alibava-pedestal' 'alibava-commonmode' 'alibava-pedestal2')


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

