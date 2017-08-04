## Example for a EUTelescope analysis using the GBL processors for telescope data only.


Detailed information on how to use EUTELESCOPE software is available at

http://eutelescope.desy.de



Go to $EUTELESCOPE/jobsub/examples/gbl_local and create the directories needed for running the analysis:

`mkdir -p ./output/histograms && mkdir -p ./output/database && mkdir -p ./output/logs && mkdir -p ./output/lcio`



Source the build_env.sh from the EUTel main directory to set the correct environment:

`source $EUTELESCOPE/build_env.sh`



The analysis we want to perform is controlled by a config file (config.cfg), a csv-table (runlist.csv) and steering file templates (*.xml), all located in the working directory.

You can set the variables for your working directory here as well as input parameters for the different processors.



You should be able to run the analysis now. For this, execute the following commands one after another and check for each of the processors to successfully finish.

```
jobsub -c config.cfg -csv runlist.csv -g converter 117
jobsub -c config.cfg -csv runlist.csv -g clustering 117
jobsub -c config.cfg -csv runlist.csv -g hitmaker 117
jobsub -c config.cfg -csv runlist.csv -g aligngbl 117
jobsub -c config.cfg -csv runlist.csv -g aligngbl2 117
jobsub -c config.cfg -csv runlist.csv -g aligngbl3 117
jobsub -c config.cfg -csv runlist.csv -g trackgbltriplet 117
```


At every step a ROOT file is created, containing a set of histograms, which you can find at output/histograms/ after completion.
