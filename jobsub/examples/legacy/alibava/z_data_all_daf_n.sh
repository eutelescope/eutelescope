#!/bin/bash
INPUT=runlists/n.csv
MAXJOBS=10
jobcount=0
T="$(date "+%d/%m/%y %T")"
touch logoutput_data_daf_n.txt
echo "Starting script at $T !" > logoutput_data_daf_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_daf_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running converter on run $runnumber at $T!" >> logoutput_data_daf_n.txt
    jobsub.py -c config.cfg -csv runlist.csv converter $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all converter runs at $T! Proceeding with reco!" >> logoutput_data_daf_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_daf_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running reco on run $runnumber at $T!" >> logoutput_data_daf_n.txt
    jobsub.py -c config.cfg -csv runlist.csv reco $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all reco runs at $T! Proceeding with rghfilter!" >> logoutput_data_daf_n.txt

wait

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_daf_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running rghfilter on run $runnumber at $T!" >> logoutput_data_daf_n.txt
    jobsub.py -c config.cfg -csv runlist.csv rghfilter $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all rghfilter runs at $T! Proceeding with clustering-1-afterrgh!" >> logoutput_data_daf_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_daf_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running clustering-1-afterrgh on run $runnumber at $T!" >> logoutput_data_daf_n.txt
    jobsub.py -c config.cfg -csv runlist.csv clustering-1-afterrgh $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all clustering-1-afterrgh runs at $T! Proceeding with clustering-2-afterrgh!" >> logoutput_data_daf_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_daf_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running clustering-2-afterrgh on run $runnumber at $T!" >> logoutput_data_daf_n.txt
    jobsub.py -c config.cfg -csv runlist.csv clustering-2-afterrgh $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all clustering-2-afterrgh runs at $T! Proceeding with merge!" >> logoutput_data_daf_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_daf_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running merge on run $runnumber at $T!" >> logoutput_data_daf_n.txt
    jobsub.py -c config.cfg -csv runlist.csv merge $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all merge runs at $T! Proceeding with hitmaker!" >> logoutput_data_daf_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_daf_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running hitmaker on run $runnumber at $T!" >> logoutput_data_daf_n.txt
    jobsub.py -c config.cfg -csv runlist.csv hitmaker $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all hitmaker runs at $T! Proceeding with alignment-daf-1!" >> logoutput_data_daf_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_daf_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running alignment-daf-1 on run $runnumber at $T!" >> logoutput_data_daf_n.txt
    jobsub.py -c config.cfg -csv runlist.csv alignment-daf-1 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all alignment-daf-1 runs at $T! Proceeding with alignment-daf-2!" >> logoutput_data_daf_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_daf_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running alignment-daf-2 on run $runnumber at $T!" >> logoutput_data_daf_n.txt
    jobsub.py -c config.cfg -csv runlist.csv alignment-daf-2 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all alignment-daf-2 runs at $T! Proceeding with alignment-daf-3!" >> logoutput_data_daf_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_daf_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running alignment-daf-3 on run $runnumber at $T!" >> logoutput_data_daf_n.txt
    jobsub.py -c config.cfg -csv runlist.csv alignment-daf-3 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all alignment-daf-3 runs at $T! Proceeding with alignment-daf-4!" >> logoutput_data_daf_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_daf_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running alignment-daf-4 on run $runnumber at $T!" >> logoutput_data_daf_n.txt
    jobsub.py -c config.cfg -csv runlist.csv alignment-daf-4 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all alignment-daf-4 runs at $T! Proceeding with alignment-daf-5!" >> logoutput_data_daf_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_daf_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running alignment-daf-5 on run $runnumber at $T!" >> logoutput_data_daf_n.txt
    jobsub.py -c config.cfg -csv runlist.csv alignment-daf-5 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all alignment-daf-5 runs at $T! Proceeding with tracking-1!" >> logoutput_data_daf_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_daf_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running tracking-1 on run $runnumber at $T!" >> logoutput_data_daf_n.txt
    jobsub.py -c config.cfg -csv runlist.csv tracking-1 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all tracking-1 runs at $T! Proceeding with alignment-daf-7!" >> logoutput_data_daf_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_daf_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running alignment-daf-7 on run $runnumber at $T!" >> logoutput_data_daf_n.txt
    jobsub.py -c config.cfg -csv runlist.csv alignment-daf-7 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all alignment-daf-7 runs at $T! Proceeding with tracking-2!" >> logoutput_data_daf_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_daf_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running tracking-2 on run $runnumber at $T!" >> logoutput_data_daf_n.txt
    jobsub.py -c config.cfg -csv runlist.csv tracking-2 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all tracking-2 runs at $T! Proceeding with alignment-daf-9!" >> logoutput_data_daf_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_daf_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running alignment-daf-9 on run $runnumber at $T!" >> logoutput_data_daf_n.txt
    jobsub.py -c config.cfg -csv runlist.csv alignment-daf-9 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all alignment-daf-9 runs at $T! Proceeding with alignment-daf-10!" >> logoutput_data_daf_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_daf_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running alignment-daf-10 on run $runnumber at $T!" >> logoutput_data_daf_n.txt
    jobsub.py -c config.cfg -csv runlist.csv alignment-daf-10 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Done all alignment-daf-10 runs at $T! Proceeding with tracking-3!" >> logoutput_data_daf_n.txt

OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read runnumber type fail
do
  if (( $jobcount == $MAXJOBS  )) ; then
    T="$(date "+%d/%m/%y %T")"
    echo "Waiting! Job limit reached at $T!" >> logoutput_data_daf_n.txt
    wait
    jobcount=0
  fi
  case "$runnumber" in \#*|"RunNumber"|"")
    continue
  ;; esac
  case "$type" in "dat")
    let jobcount+=1
    T="$(date "+%d/%m/%y %T")"
    echo "Running tracking-3 on run $runnumber at $T!" >> logoutput_data_daf_n.txt
    jobsub.py -c config.cfg -csv runlist.csv tracking-3 $runnumber &
  ;; esac
done < $INPUT
IFS=$OLDIFS

wait

T="$(date "+%d/%m/%y %T")"
echo "Finished all at $T!" >> logoutput_data_daf_n.txt

#

