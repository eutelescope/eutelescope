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
CONFIG=config-ibl-phase1.cfg
# -----------------------------------------------------------------------------
# -- additional variables
# planes to be excluded during DUT alignment
PLANE10="11 12 "
PLANE11="10 12 "
PLANE12="10 11 "
PLANES=( "$PLANE10" "$PLANE11" "$PLANE12" )
ALLDUTS="10 11 12 "
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

    # standard convertion
	./submit-converter.py $RUN --config=config/$CONFIG > $logs/converter.$RUN.log

    # --hot parameter tells to use the htpixelDB file 
    ./submit-clusearch.py $RUN --config=config/$CONFIG --hot $RUNHOT > $logs/clusearch.$RUN.log

	# 5-digit runnumber
    # -n parameter tells to use te offsetDB values
	./submit-hitmaker.py -o $RUN run0$RUN-clu-p.slcio --config=config/$CONFIG -n    $RUN > $logs/hitmaker.$RUN.log


    # 
    # Millepede alignment (start with M26 sensors)
    #
    echo ./submit-align.py -f \"$TELPLANESFIXED\" -e \"$ALLDUTS\" -i hit -o $RUN-iter1 $RUN-hit.slcio  --config=config/$CONFIG 
    ./submit-align.py -f "$TELPLANESFIXED" -e "$ALLDUTS" -i hit -o $RUN-iter1 $RUN-hit.slcio  --config=config/$CONFIG > $logs/align-iter1.$RUN.log

	# apply this alignment
	# --- REMINDER: parameters for ./submit-apply-alignment.py
	# -a - alignment file name, given by -o option of previous processor, i.e.  $RUN-iter1
	# -r RunNumber, -i InputHitCollection, -o OutputHitCollection
	# output lcio hit file name is constructed like RunNumber-OutputHitCollection.slcio
	# output root file name is constructed like RunNumber-OutputHitCollection.root ----
	echo ./submit-apply-alignment.py -a $RUN-iter1-align-db.slcio -r $RUN -i hit -o iter1_hit --InputLCIOFile $RUN-hit.slcio --config=config/$CONFIG
    ./submit-apply-alignment.py -a $RUN-iter1-align-db.slcio -r $RUN -i hit -o iter1_hit --InputLCIOFile $RUN-hit.slcio --config=config/$CONFIG



	rm $RESULTS/$RUN-iter1_hit.000.slcio


	# -----------
	# --- sequential DUT: alignment and apply-alignment ---#
	DUTindex=0
	for i in `seq  2  4`;
	do
		let j=i-1			# for input collections
		in=iter$j'_hit'		# input collections
		out=iter$i'_hit'
		outBase=$RUN-iter$i	# output base
		#OPTION=--use-residual-cuts
        echo  ./submit-align.py -f \"$TELPLANES\" -e \"${PLANES[$DUTindex]}\" -i $in -o $outBase $RUN-$in.000.slcio $OPTION --config=config/$CONFIG 
	    ./submit-align.py -f "$TELPLANES" -e "${PLANES[$DUTindex]}" -i $in -o $outBase $RUN-$in.000.slcio $OPTION --config=config/$CONFIG >$logs/align-iter$i.$RUN.log
		# apply alignment
		# the last applied alignment should have output collection name alignedHit - otherwise fitter steering file needs help ;)
		if [ $i -eq 4 ]; then out=alignedHit; fi
 		if [ $i -eq 4 ]; then out=alignedhit; fi
        echo ./submit-apply-alignment.py -a $outBase-align-db.slcio -r $RUN -i $in -o $out --InputLCIOFile $RUN-$in.000.slcio --config=config/$CONFIG
        ./submit-apply-alignment.py -a $outBase-align-db.slcio -r $RUN -i $in -o $out --InputLCIOFile $RUN-$in.000.slcio --config=config/$CONFIG
   		rm $RESULTS/$RUN-$in.000.slcio
		let DUTindex++
	done
	# ------------------------------------------------------------

	# -- Finally! Fit the tracks! --
	echo ./submit-fitter.py -o $RUN $RUN-alignedHit.000.slcio -a $RUN-iter4-align-db.slcio --config=config/$CONFIG  
    ./submit-fitter.py -o $RUN $RUN-alignedHit.000.slcio -a $RUN-iter4-align-db.slcio --config=config/$CONFIG  > $logs/fitter.$RUN.log
	echo "done for run $RUN"
done
