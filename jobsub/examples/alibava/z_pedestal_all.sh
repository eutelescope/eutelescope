#!/bin/bash
INPUT=runlist.csv
MAXJOBS=5
jobcount=0
T="$(date "+%d/%m/%y %T")"
touch logoutput_pedestal.txt
echo "Starting script at $T !" > logoutput_pedestal.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_pedestal.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "ped")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running pedestal converter on run $runnumber at $T!" >> logoutput_pedestal.txt
    jobsub.py -c config.cfg -csv runlist.csv converter $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all pedestal converter runs at $T! Proceeding with pedestal!" >> logoutput_pedestal.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_pedestal.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "ped")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running pedestal on run $runnumber at $T!" >> logoutput_pedestal.txt
    jobsub.py -c config.cfg -csv runlist.csv pedestal $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all pedestal runs at $T! Proceeding with commonmode!" >> logoutput_pedestal.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_pedestal.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "ped")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running commonmode on run $runnumber at $T!" >> logoutput_pedestal.txt
    jobsub.py -c config.cfg -csv runlist.csv commonmode $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all commonmode runs at $T! Proceeding with pedestal2!" >> logoutput_pedestal.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_pedestal.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "ped")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running pedestal2 on run $runnumber at $T!" >> logoutput_pedestal.txt
    jobsub.py -c config.cfg -csv runlist.csv pedestal2 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all pedestal2 runs at $T! Proceeding with pedestalhisto!" >> logoutput_pedestal.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_pedestal.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "ped")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running pedestalhisto on run $runnumber at $T!" >> logoutput_pedestal.txt
    jobsub.py -c config.cfg -csv runlist.csv pedestalhisto $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Finished all at $T!" >> logoutput_pedestal.txt

#
