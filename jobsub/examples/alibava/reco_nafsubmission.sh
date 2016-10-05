#!/bin/zsh

###########################################################
#                                                         #
# ALiBaVa Analysis - Signal Reconstruction and clustering #
#                                                         #
###########################################################

# usage: source reco_nafsubmission.sh runnumber(s)

# To get pedestal subtracted, common mode and cross talk corrected signal values
# you will need to run these templates
TEMPLATELIST=('alibava-converter' 'alibava-reco' 'alibava-applyxtalk')

for TEMPLATE in $TEMPLATELIST; do
	jobsub -c config/config.cfg -csv runlistfiles/runlist.csv --subdir --naf qsubparameters.cfg $TEMPLATE $@

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

