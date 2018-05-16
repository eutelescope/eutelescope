#!/bin/sh

################################
#                              #
# ALiBaVa Analysis - Telescope #
#                              #
# for concatenated telescope   #
# runs                         #
#                              #
################################

# usage: sh x_telescope-multi.sh runnumber

# Warning: if your lcio files do NOT lie in ./output/lcio, you should adjust the path!

jobsub.py -c config.cfg -csv runlist.csv telescope-converter $1
jobsub.py -c config.cfg -csv runlist.csv telescope-converter $(($1+1))
jobsub.py -c config.cfg -csv runlist.csv telescope-converter $(($1+2))
jobsub.py -c config.cfg -csv runlist.csv telescope-converter $(($1+3))
minussign='-'
end=$(($1+3))
range=$1$minussign$end
jobsub -c config.cfg --concatenate --option concatinput=output/lcio/run@RunRange@-converter.slcio -csv runlist.csv telescope-clustering-concat $range
jobsub.py -c config.cfg -csv runlist.csv telescope-filter $1

#
