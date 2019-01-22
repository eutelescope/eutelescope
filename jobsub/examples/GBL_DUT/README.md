# Example for a EUTelescope analysis using the GBL processors and an active DUT (+ a pixel reference plane).

Detailed information on how to use EUTELESCOPE software is available at http://eutelescope.desy.de

This example consists in a 2018 data set including an ATLAS strip sensor as an active DUT and a FEI4 reference plane placed before the last telescope plane. The measurement was performed at DESY in TB 22 using the EUDET-type telescope DURANTA.

## Prerequisites

Source the build_env.sh from the EUTel main directory to set the correct environment:

`source $EUTELESCOPE/build_env.sh`

## Running

To run the example, execute the following commands one after another and check for each of the processors to successfully finish.

```
jobsub -c config.cfg -csv runlist.csv -g noisypixel 2365
jobsub -c config.cfg -csv runlist.csv -g clustering 2365
jobsub -c config.cfg -csv runlist.csv -g hitmaker 2365
jobsub -c config.cfg -csv runlist.csv -g alignGBL1 2365
jobsub -c config.cfg -csv runlist.csv -g alignGBL2 2365
jobsub -c config.cfg -csv runlist.csv -g alignGBL3 2365
jobsub -c config.cfg -csv runlist.csv -g fitGBL 2365
```

## Output

At every step a ROOT file is created, containing a set of histograms, which you can be found at output/histograms/ after completion. The final track and hit parameters are stored in a root NTupla called NTuple.root.
