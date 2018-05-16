#!/bin/bash
INPUT=runlists/n.csv
MAXJOBS=10
jobcount=0
T="$(date "+%d/%m/%y %T")"
touch logoutput_data_gbl_n.txt
echo "Starting script at $T !" > logoutput_data_gbl_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_gbl_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running converter on run $runnumber at $T!" >> logoutput_data_gbl_n.txt
    jobsub.py -c config.cfg -csv runlist.csv converter $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all converter runs at $T! Proceeding with reco!" >> logoutput_data_gbl_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_gbl_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running reco on run $runnumber at $T!" >> logoutput_data_gbl_n.txt
    jobsub.py -c config.cfg -csv runlist.csv reco $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all reco runs at $T! Proceeding with rghfilter!" >> logoutput_data_gbl_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_gbl_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running rghfilter on run $runnumber at $T!" >> logoutput_data_gbl_n.txt
    jobsub.py -c config.cfg -csv runlist.csv rghfilter $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all rghfilter runs at $T! Proceeding with clustering-1-afterrgh!" >> logoutput_data_gbl_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_gbl_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running clustering-1-afterrgh on run $runnumber at $T!" >> logoutput_data_gbl_n.txt
    jobsub.py -c config.cfg -csv runlist.csv clustering-1-afterrgh $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all clustering-1-afterrgh runs at $T! Proceeding with clustering-2-afterrgh!" >> logoutput_data_gbl_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_gbl_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running clustering-2-afterrgh on run $runnumber at $T!" >> logoutput_data_gbl_n.txt
    jobsub.py -c config.cfg -csv runlist.csv clustering-2-afterrgh $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all clustering-2-afterrgh runs at $T! Proceeding with merge!" >> logoutput_data_gbl_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_gbl_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running merge on run $runnumber at $T!" >> logoutput_data_gbl_n.txt
    jobsub.py -c config.cfg -csv runlist.csv merge $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all merge runs at $T! Proceeding with hitmaker!" >> logoutput_data_gbl_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_gbl_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running hitmaker on run $runnumber at $T!" >> logoutput_data_gbl_n.txt
    jobsub.py -c config.cfg -csv runlist.csv hitmaker $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all hitmaker runs at $T! Proceeding with coordinator!" >> logoutput_data_gbl_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_gbl_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running coordinator on run $runnumber at $T!" >> logoutput_data_gbl_n.txt
    jobsub.py -c config.cfg -csv runlist.csv coordinator $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all coordinator runs at $T! Proceeding with alignment-gbl-1!" >> logoutput_data_gbl_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_gbl_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running alignment-gbl-1 on run $runnumber at $T!" >> logoutput_data_gbl_n.txt
    jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-1 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all alignment-gbl-1 runs at $T! Proceeding with alignment-gbl-2!" >> logoutput_data_gbl_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_gbl_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running alignment-gbl-2 on run $runnumber at $T!" >> logoutput_data_gbl_n.txt
    jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-2 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all alignment-gbl-2 runs at $T! Proceeding with alignment-gbl-3!" >> logoutput_data_gbl_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_gbl_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running alignment-gbl-3 on run $runnumber at $T!" >> logoutput_data_gbl_n.txt
    jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-3 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all alignment-gbl-3 runs at $T! Proceeding with alignment-gbl-4!" >> logoutput_data_gbl_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_gbl_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running alignment-gbl-4 on run $runnumber at $T!" >> logoutput_data_gbl_n.txt
    jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-4 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all alignment-gbl-4 runs at $T! Proceeding with alignment-gbl-5!" >> logoutput_data_gbl_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_gbl_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running alignment-gbl-5 on run $runnumber at $T!" >> logoutput_data_gbl_n.txt
    jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-5 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all alignment-gbl-5 runs at $T! Proceeding with alignment-gbl-6!" >> logoutput_data_gbl_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_gbl_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running alignment-gbl-6 on run $runnumber at $T!" >> logoutput_data_gbl_n.txt
    jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-6 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all alignment-gbl-6 runs at $T! Proceeding with alignment-gbl-7!" >> logoutput_data_gbl_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_gbl_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running alignment-gbl-7 on run $runnumber at $T!" >> logoutput_data_gbl_n.txt
    jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-7 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all alignment-gbl-7 runs at $T! Proceeding with alignment-gbl-8!" >> logoutput_data_gbl_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_gbl_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running alignment-gbl-8 on run $runnumber at $T!" >> logoutput_data_gbl_n.txt
    jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-8 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all alignment-gbl-8 runs at $T! Proceeding with alignment-gbl-9!" >> logoutput_data_gbl_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_gbl_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running alignment-gbl-9 on run $runnumber at $T!" >> logoutput_data_gbl_n.txt
    jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-9 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all alignment-gbl-9 runs at $T! Proceeding with alignment-gbl-10!" >> logoutput_data_gbl_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_gbl_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running alignment-gbl-10 on run $runnumber at $T!" >> logoutput_data_gbl_n.txt
    jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-10 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all alignment-gbl-10 runs at $T! Proceeding with alignment-gbl-11!" >> logoutput_data_gbl_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_gbl_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running alignment-gbl-11 on run $runnumber at $T!" >> logoutput_data_gbl_n.txt
    jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-11 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all alignment-gbl-11 runs at $T! Proceeding with alignment-gbl-12!" >> logoutput_data_gbl_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_gbl_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running alignment-gbl-12 on run $runnumber at $T!" >> logoutput_data_gbl_n.txt
    jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-12 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all alignment-gbl-12 runs at $T! Proceeding with alignment-gbl-13!" >> logoutput_data_gbl_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_gbl_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running alignment-gbl-13 on run $runnumber at $T!" >> logoutput_data_gbl_n.txt
    jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-13 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all alignment-gbl-13 runs at $T! Proceeding with alignment-gbl-14!" >> logoutput_data_gbl_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_gbl_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running alignment-gbl-14 on run $runnumber at $T!" >> logoutput_data_gbl_n.txt
    jobsub.py -c config.cfg -csv runlist.csv alignment-gbl-14 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all alignment-gbl-14 runs at $T! Proceeding with tracking-gbl!" >> logoutput_data_gbl_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_gbl_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running tracking-gbl on run $runnumber at $T!" >> logoutput_data_gbl_n.txt
    jobsub.py -c config.cfg -csv runlist.csv tracking-gbl $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Finished all at $T!" >> logoutput_data_gbl_n.txt

#

