#!/bin/bash
# -----------------------------------------------------------------------------
# Script for running full EUTelescope analysis chain
# IMPORTANT: Check this file carefully,
# in particular alignment steps you want to do.
# In this example, the following alignment steps are performed:
# 1. Align telescope planes with correlation bands + apply alignment constants
# 2. Improve it with standard residual cuts approach + apply alignment constants
# 3. Align each DUT separately with correlation bands
# 4. Apply alignment
# 5. Improve with std cuts
# 6. Apply alignment
# NOTE: Last output collection after alignment should be named alignedHit
# (default name of the input collection in the fitter)
#
# Author: Vladyslav Libov, libov@mail.desy.de
# Last modified: 30 August 2010
#
# Edited: by Igor Rubinskiy, 21-11-2010
#
# -----------------------------------------------------------------------------
# -- choose set of runs to be processed
#RUNSET=`seq 11061 11191`
RUNSET=$1

RUNHOT="10494-01"
echo set of runs to be processed: $RUNSET
# -- choose config file (should be placed in pysub/config/)
CONFIG=config-ibl-phase1-rot.cfg
# -----------------------------------------------------------------------------
# -- additional variables
# planes to be excluded during DUT alignment
NONE=""
# telescope planes
TELPLANES="0 1 2 3 4 5"
# fixed planes for telescope alignment
TELPLANESFIXED="0 5"
export logs="$EUTELESCOPE/pysub/logs"
#
# local folder name to be edited by each user, example:
RESULTS="/data/zenith229a/rubinsky/EUDET/TestBeam/Data/2010/tb-cern-autumn/ibl/results"
#
# -----------------------------------------------------------------------------
for RUN in $RUNSET;
do
	echo Running full analysis chain for run $RUN
   if [ -d $logs ]; then
		echo "Writing logfiles to $logs/"
	else
		echo "Directory $logs/ does not exist. Creating..."
		mkdir $logs
		echo "Writing logfiles to $logs/"
	fi
	rm -rf $logs/*.$RUN*.log

    #
    # standard convertion from RAW to LCIO data format
	#
    echo ./submit-converter.py $RUN --config=config/$CONFIG > $logs/converter.$RUN.log
    ./submit-converter.py $RUN --config=config/$CONFIG > $logs/converter.$RUN.log

    #
    # --hot parameter tells to use the htpixelDB file 
    #
    echo ./submit-clusearch.py $RUN --config=config/$CONFIG --hot $RUNHOT > $logs/clusearch.$RUN.log
    ./submit-clusearch.py $RUN --config=config/$CONFIG --hot $RUNHOT > $logs/clusearch.$RUN.log

    #
    # -n parameter tells to use the offsetDB values
    # offsetDB is equivalent to preAlignment
	#
    echo  ./submit-hitmaker.py -o $RUN run0$RUN-clu-p.slcio --config=config/$CONFIG -n    $RUN > $logs/hitmaker.$RUN.log
    ./submit-hitmaker.py -o $RUN run0$RUN-clu-p.slcio --config=config/$CONFIG -n    $RUN > $logs/hitmaker.$RUN.log

    #
    # align all Planes!
    # alignment with Millepede, fine alignment!
    #
    # it works only when it is possible to have a track goign through all sensor
    # planes at the same time, due to wrong sensor positioning during the beam
    # test setup it happens so that a pair of DUTs can be displaced.'
    #
    echo ./submit-align.py -f "$TELPLANESFIXED" -e "$NONE"  -o $RUN $RUN-hit.slcio  --config=config/$CONFIG > $logs/align.$RUN.log
    ./submit-align.py -f "$TELPLANESFIXED" -e "$NONE"  -o $RUN $RUN-hit.slcio  --config=config/$CONFIG > $logs/align.$RUN.log

	#
    # -- Finally! Fit the tracks! --
    #
    echo ./submit-fitter.py -o $RUN $RUN-hit.slcio -a $RUN-align-db.slcio --config=config/$CONFIG  > $logs/fitter.$RUN.log
    ./submit-fitter.py -o $RUN $RUN-hit.slcio -a $RUN-align-db.slcio --config=config/$CONFIG  > $logs/fitter.$RUN.log

    echo "done for run $RUN"
done
