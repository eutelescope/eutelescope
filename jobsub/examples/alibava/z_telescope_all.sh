#!/bin/bash
INPUT=runlist.csv
MAXJOBS=10
jobcount=0
T="$(date "+%d/%m/%y %T")"
touch logoutput_telescope.txt
echo "Starting script at $T !" > logoutput_telescope.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_telescope.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in ""|"c"|"d")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running converter on run $runnumber at $T!" >> logoutput_telescope.txt
    jobsub.py -c config.cfg -csv runlist.csv telescope-converter $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all converter runs at $T! Proceeding with clustering!" >> logoutput_telescope.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_telescope.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running clustering on run $runnumber at $T!" >> logoutput_telescope.txt
    jobsub.py -c config.cfg -csv runlist.csv telescope-clustering $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all clustering runs at $T! Proceeding with concat clustering!" >> logoutput_telescope.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_telescope.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "c")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running concat clustering on run $runnumber at $T!" >> logoutput_telescope.txt
    jobsub -c config.cfg --concatenate --option concatinput=output/lcio/run@RunRange@-converter.slcio -csv runlist.csv telescope-clustering-concat $runnumber-0`expr $runnumber + 3` &

  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all concat clustering runs at $T! Proceeding with telescope filtering!" >> logoutput_telescope.txt

LDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_telescope.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in ""|"c")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running filter on run $runnumber at $T!" >> logoutput_telescope.txt
    jobsub.py -c config.cfg -csv runlist.csv telescope-filter $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Finished all at $T!" >> logoutput_telescope.txt

#
