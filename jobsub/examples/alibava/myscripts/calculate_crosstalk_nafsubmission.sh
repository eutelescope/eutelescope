#!/bin/zsh

#########################################################################################
#                                                                                       #
# ALiBaVa Analysis - Calculates CrossTalk and does Signal Reconstruction and clustering #
#                                                                                       #
#########################################################################################

# usage: source reco.sh runnumber(s)

# To get pedestal subtracted, common mode and cross talk corrected signal values
# you will need to run these templates
#TEMPLATELIST=('alibava-convert-ped' 'alibava-pedestal' 'alibava-commonmode' 'alibava-pedestal2')
#TEMPLATELIST=('alibava-pedestal' 'alibava-commonmode' 'alibava-pedestal2')

#TEMPLATELIST=('alibava-converter' 'alibava-reco' 'alibava-crosstalk-it1' 'alibava-crosstalk-it2' 'alibava-crosstalk-it3')
TEMPLATELIST=('alibava-reco' 'alibava-crosstalk-it1' 'alibava-crosstalk-it2' 'alibava-crosstalk-it3')

for TEMPLATE in $TEMPLATELIST; do
	jobsub -c config/dut4_cce_config.cfg -csv runlistfiles/runlist.csv --subdir --naf qsubparameters.cfg $TEMPLATE $@

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


