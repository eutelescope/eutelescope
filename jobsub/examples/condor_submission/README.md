# Job Submission to batch via HTCondor

To find out more about HTCondor submission, try:
- https://confluence.desy.de/pages/viewpage.action?pageId=67639562

To submit your EUTel reconstruction job to HTCondor, you add just the ```--condor``` flag plus the configuration file ```condorparameters.cfg``` to jobsub.

E.g. this can then look like

```
jobsub -c config.cfg -csv runlist.csv --condor condorparameters.cfg JOBTASK RUNNR
```

with the parameters JOBTASK as the wished reconstruction step and RUNNR as the desired run number.


This will send the job to HTCondor (using the script ```JOBTASK-RUNNR_jobsub.sh``` by executing ```condor_submit JOBTASK-RUNNR_jobsub.submit```).
The output logs are stored as usal in ```./output/logs```. In total there are three logs: ```*.log``` (normal execution log), ```*.error``` (error states), ```*.condor``` (log from submit node).
Furthermore, the submit scripts are stored there.

For trying out the HTCondor submission, please copy the file ```condorparameters.cfg``` to your reconstruction folder and use the flag ```--condor``` for the reconstruction step.

It is possible to make some adaptions to the submit procedure, for this please check the configuration file ```condorparameters.cfg```.

To check for the status of your job, you can use ```condor_q```