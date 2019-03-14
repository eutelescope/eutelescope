# Example for a EUTelescope analysis using the GBL processors for telescope data only (empty telescope = no DUT).

Detailed information on how to use EUTELESCOPE software is available at http://eutelescope.desy.de

## Prerequisites

Source the build_env.sh from the EUTel main directory to set the correct environment:

`source $EUTELESCOPE/build_env.sh`

## Running

The analysis we want to perform is controlled by a config file (config.cfg), a csv-table (runlist.csv) and steering file templates (*.xml), all located in the working directory.

You can set the variables for your working directory here as well as input parameters for the different processors.


There a two different versions, one for raw-data taken with EUDAQ1 and one with already converted data files taken with EUDAQ2.

### EUDAQ2 version

You should be able to run the analysis now. For this, execute the following commands one after another and check for each of the processors to successfully finish.

```
jobsub -c config_DESY2018.cfg -csv runlist.csv -g noisypixel 701
jobsub -c config_DESY2018.cfg -csv runlist.csv -g clustering 701
jobsub -c config_DESY2018.cfg -csv runlist.csv -g hitmaker 701
jobsub -c config_DESY2018.cfg -csv runlist.csv -g alignGBL 701
jobsub -c config_DESY2018.cfg -csv runlist.csv -g alignGBL2 701
jobsub -c config_DESY2018.cfg -csv runlist.csv -g alignGBL3 701
jobsub -c config_DESY2018.cfg -csv runlist.csv -g fitGBL 701
```


### EUDAQ1 version

For this version, one has to convert the file in the first step. Furthermore, it is then needed to change the inputs for the following processes to the newly created LCIO file.
Please have a look in the steering-templates for noisypixel and clustering; here comment out the current LCIOInput and take the other provided one.
Then the analysis should work in the same way:

```
jobsub -c config_DESY2012.cfg -csv runlist.csv -g converter 117
jobsub -c config_DESY2012.cfg -csv runlist.csv -g noisypixel 117
jobsub -c config_DESY2012.cfg -csv runlist.csv -g clustering 117
jobsub -c config_DESY2012.cfg -csv runlist.csv -g hitmaker 117
jobsub -c config_DESY2012.cfg -csv runlist.csv -g alignGBL 117
jobsub -c config_DESY2012.cfg -csv runlist.csv -g alignGBL2 117
jobsub -c config_DESY2012.cfg -csv runlist.csv -g alignGBL3 117
jobsub -c config_DESY2012.cfg -csv runlist.csv -g fitGBL 117
```

## Output

At every step a ROOT file is created, containing a set of histograms, which you can be found at output/histograms/ after completion.