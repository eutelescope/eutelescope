# Example for a EUTelescope analysis using the DAF processors for Aconite telescope with 2 DUTs (FEI4Quad + FEI4Single).

Detailed information on how to use EUTELESCOPE software is available at http://eutelescope.desy.de

## Prerequisites

Source the build_env.sh from the EUTel main directory to set the correct environment:

`source $EUTELESCOPE/build_env.sh`

## Running

The analysis we want to perform is controlled by a config file (config.cfg), a csv-table (runlist.csv) and steering file templates (*.xml), all located in the working directory.

You can set the variables for your working directory here as well as input parameters for the different processors.

The data was taken with EUDAQ1, therefore as first step the conversion of the raw data has to be done.


You should be able to run the analysis now. For this, execute the following commands one after another and check for each of the processors to successfully finish.

```
jobsub -c config.cfg -csv runlist.csv -g converter 1085
jobsub -c config.cfg -csv runlist.csv -g clustering 1085
jobsub -c config.cfg -csv runlist.csv -g hitmaker 1085
jobsub -c config.cfg -csv runlist.csv -g alignDAF 1085
jobsub -c config.cfg -csv runlist.csv -g fitDAF 1085
```

## Output

At every step a ROOT file is created, containing a set of histograms, which you can be found at output/histograms/ after completion. Furthermore, the final track data is dumped to a ROOT
ntuple, called run001085_output.root .
