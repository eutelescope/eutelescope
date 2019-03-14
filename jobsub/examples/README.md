
This subdirectory contains example configurations and steering files
that can be used as a basis for your own analysis:

* GBL_noDUT
     The Datura telescope with six planes of Mimosa26 without DUT

* GBL_SUT
     The Datura telescope with six planes of Mimosa26 and a SUT for material budget imaging

* GBL_DUT
     The Datura telescope with six planes of Mimosa26 and a FEI4 detector as DUT plane

* mixed-mode
     Test data for EUDAQ2 and AIDA-TLU operations

* condor-submission
     Example file for the condor_submit parameters needed for NAF job submission

All examples are also being used for automated data-driven tests using
the CMake/CTest framework to verify the correct functionality of
EUTelescope; see 'datura-alone/testing.cmake' for an example setup for
such a test.

All examples are self-contained and only refer to other files within
their subfolders (with the exception of the raw data files).

Please consider to add your configuration as an example for others and
as basis for automated and constant verification: you would need to
provide the necessary steering-files together with a data sample and
documentation. Please contact the EUTelescope software coordinators
for further information.