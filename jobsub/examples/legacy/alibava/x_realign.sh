#!/bin/bash

################################
#                              #
# ALiBaVa Realignment          #
#                              #
################################

# usage: sh x_realign.sh badrunnumber goodrunnumber

################################
# Applies the alignment of <goodrunnumber> to <badrunnumber>
# Alignment constants of <badrunnumber> are backuped beforehand
# Run numbers should be 4 digits!
################################

# select mode
mode=gbl
#mode=daf

if (($mode == gbl)) ; then

    cp output/database/run00$1-prealignment.slcio output/database/run00$1-original-prealignment.slcio
    cp output/database/run00$1-alignment-1-gbl.slcio output/database/run00$1-original-alignment-1-gbl.slcio
    cp output/database/run00$1-alignment-2-gbl.slcio output/database/run00$1-original-alignment-2-gbl.slcio
    cp output/database/run00$1-alignment-3-gbl.slcio output/database/run00$1-original-alignment-3-gbl.slcio
    cp output/database/run00$1-alignment-4-gbl.slcio output/database/run00$1-original-alignment-4-gbl.slcio
    cp output/database/run00$1-alignment-5-gbl.slcio output/database/run00$1-original-alignment-5-gbl.slcio
    cp output/database/run00$1-alignment-6-gbl.slcio output/database/run00$1-original-alignment-6-gbl.slcio
    cp output/database/run00$1-alignment-7-gbl.slcio output/database/run00$1-original-alignment-7-gbl.slcio
    cp output/database/run00$1-alignment-8-gbl.slcio output/database/run00$1-original-alignment-8-gbl.slcio
    cp output/database/run00$1-alignment-9-gbl.slcio output/database/run00$1-original-alignment-9-gbl.slcio
    cp output/database/run00$1-alignment-10-gbl.slcio output/database/run00$1-original-alignment-10-gbl.slcio
    cp output/database/run00$1-alignment-11-gbl.slcio output/database/run00$1-original-alignment-11-gbl.slcio
    cp output/database/run00$1-alignment-12-gbl.slcio output/database/run00$1-original-alignment-12-gbl.slcio
    cp output/database/run00$1-alignment-13-gbl.slcio output/database/run00$1-original-alignment-13-gbl.slcio
    cp output/database/run00$1-alignment-14-gbl.slcio output/database/run00$1-original-alignment-14-gbl.slcio

    echo "Backed up all old gbl alignment files"

    cp output/database/run00$2-prealignment.slcio output/database/run00$1-prealignment.slcio
    cp output/database/run00$2-alignment-1-gbl.slcio output/database/run00$1-alignment-1-gbl.slcio
    cp output/database/run00$2-alignment-2-gbl.slcio output/database/run00$1-alignment-2-gbl.slcio
    cp output/database/run00$2-alignment-3-gbl.slcio output/database/run00$1-alignment-3-gbl.slcio
    cp output/database/run00$2-alignment-4-gbl.slcio output/database/run00$1-alignment-4-gbl.slcio
    cp output/database/run00$2-alignment-5-gbl.slcio output/database/run00$1-alignment-5-gbl.slcio
    cp output/database/run00$2-alignment-6-gbl.slcio output/database/run00$1-alignment-6-gbl.slcio
    cp output/database/run00$2-alignment-7-gbl.slcio output/database/run00$1-alignment-7-gbl.slcio
    cp output/database/run00$2-alignment-8-gbl.slcio output/database/run00$1-alignment-8-gbl.slcio
    cp output/database/run00$2-alignment-9-gbl.slcio output/database/run00$1-alignment-9-gbl.slcio
    cp output/database/run00$2-alignment-10-gbl.slcio output/database/run00$1-alignment-10-gbl.slcio
    cp output/database/run00$2-alignment-11-gbl.slcio output/database/run00$1-alignment-11-gbl.slcio
    cp output/database/run00$2-alignment-12-gbl.slcio output/database/run00$1-alignment-12-gbl.slcio
    cp output/database/run00$2-alignment-13-gbl.slcio output/database/run00$1-alignment-13-gbl.slcio
    cp output/database/run00$2-alignment-14-gbl.slcio output/database/run00$1-alignment-14-gbl.slcio

    echo "Copied gbl reference alignment"

    echo "Running coordinator..."
    jobsub -c config.cfg -csv runlist.csv coordinator $1
    echo "Running tracking"
    jobsub -c config.cfg -csv runlist.csv tracking-gbl $1

elif (($mode == daf)) ; then

    cp output/database/run00$1-prealignment.slcio output/database/run00$1-original-prealignment.slcio
    cp output/database/run00$1-alignment-1.slcio output/database/run00$1-original-alignment-1.slcio
    cp output/database/run00$1-alignment-2.slcio output/database/run00$1-original-alignment-2.slcio
    cp output/database/run00$1-alignment-3.slcio output/database/run00$1-original-alignment-3.slcio
    cp output/database/run00$1-alignment-4.slcio output/database/run00$1-original-alignment-4.slcio
    cp output/database/run00$1-alignment-5.slcio output/database/run00$1-original-alignment-5.slcio
    cp output/database/run00$1-alignment-6.slcio output/database/run00$1-original-alignment-6.slcio
    cp output/database/run00$1-alignment-7.slcio output/database/run00$1-original-alignment-7.slcio
    cp output/database/run00$1-alignment-8.slcio output/database/run00$1-original-alignment-8.slcio
    cp output/database/run00$1-alignment-9.slcio output/database/run00$1-original-alignment-9.slcio
    cp output/database/run00$1-alignment-10.slcio output/database/run00$1-original-alignment-10.slcio

    echo "Backed up all old daf alignment files"

    cp output/database/run00$2-prealignment.slcio output/database/run00$1-prealignment.slcio
    cp output/database/run00$2-alignment-1.slcio output/database/run00$1-alignment-1.slcio
    cp output/database/run00$2-alignment-2.slcio output/database/run00$1-alignment-2.slcio
    cp output/database/run00$2-alignment-3.slcio output/database/run00$1-alignment-3.slcio
    cp output/database/run00$2-alignment-4.slcio output/database/run00$1-alignment-4.slcio
    cp output/database/run00$2-alignment-5.slcio output/database/run00$1-alignment-5.slcio
    cp output/database/run00$2-alignment-6.slcio output/database/run00$1-alignment-6.slcio
    cp output/database/run00$2-alignment-7.slcio output/database/run00$1-alignment-7.slcio
    cp output/database/run00$2-alignment-8.slcio output/database/run00$1-alignment-8.slcio
    cp output/database/run00$2-alignment-9.slcio output/database/run00$1-alignment-9.slcio
    cp output/database/run00$2-alignment-10.slcio output/database/run00$1-alignment-10.slcio

    echo "Copied daf reference alignment"

    echo "Running coordinator..."
    jobsub -c config.cfg -csv runlist.csv coordinator $1
    echo "Running tracking"
    jobsub -c config.cfg -csv runlist.csv tracking-3 $1

fi

#
