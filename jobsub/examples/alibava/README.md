===============================================================================

examples/alibava

Thomas Eichhorn 2013-2020
thomas.eichhorn@desy.de

===============================================================================

This is the example full analysis and reconstruction of a strip sensor DUT
with an ALiBaVa readout.

===============================================================================

Folder contents:

/gearfiles:

  Contains gear files for various setup scenarios.

/output:

  The output created by the EUTelescope processors will be written to the 
subfolders herein.

/steering-templates:

  This folder contains the steering templates used to control the processors.

/config.cfg:

  The configuration file. Data paths should be changed, other settings should 
work for most scenarios.

/runlist.csv:

  A file containing the information of the runs to be processed. Run-specific 
parameters are read from this file.

/README:

  This file.

/x_<name>.sh:

  Scripts for performing a specific task on one run. They should be called 
with sh scriptname.sh <runnumber> . Further information is in the scripts.

===============================================================================

Example analysis steps:


AllPix-simulated data:

# Simulated telescope data:
$ sh x_sim_tel.sh 908
# Simulated pedestal data:
$ sh x_sim_ped.sh 900
# Simulated ALiBaVa data:
$ sh x_sim_dut.sh 908


Unirradiated p-type sensor:

# Multiple telescope runs:
$ sh x_telescope-multi.sh 12260
# Off-beam calibration:
$ sh x_calibration.sh 763
# Off-beam pedestal data:
$ sh x_pedestal.sh 19
# On-beam p-type ALiBaVa data:
$ sh x_p_dut.sh 22


Irradiated and rotated n-type sensor:

# Single telescope run:
$ sh x_telescope-single.sh 3692
# Off-beam calibration:
$ sh x_calibration.sh 453
# Off-beam pedestal data:
$ sh x_pedestal.sh 447
# On-beam n-type ALiBaVa data:
$ sh x_n_dut.sh 438

===============================================================================

