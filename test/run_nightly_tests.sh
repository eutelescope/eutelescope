#!/bin/bash
# this script runs continuous test using CMake's testing mechanisms and submits the results
# to the dashboard.
# use GNU screen to start this script from ssh sessions and then detach the session.

WAKEUPAT="02:05" # time to wake up every day in HH:MM in UTC

if [ -z "$EUTELESCOPE" ]
then
    echo " Variable \$EUTELESCOPE not set, trying to parse 'build_env.sh'"
    if [ -f ./build_env.sh ]
	then
	source build_env.sh
    else 
	if [ -f ../build_env.sh ]
	then
	    source ../build_env.sh
	else
	    echo " Could not find build_env.sh - please run script from EUTelescope root directory!"
	    exit
	fi
    fi
fi

if [ -z "$DISPLAY" ]
then
    # if $DISPLAY is not set (e.g. on a VM or server), ROOT sents a warning to stderr e.g. when loading Marlin libraries; this might fail tests
    echo " Variable \$DISPLAY not set, set to 'localhost:0' (otherwise tests might fail due to ROOT error message)"
    export DISPLAY=localhost:0
fi

cd $EUTELESCOPE
if (( $? )); then
    {
        echo " Could not change into EUTelescope directory - please set up environment variables correctly!";
        exit;
    }
fi;
echo " Environment set up correctly "
cd build

# setup done!
echo " Waiting for my time to wake up ($WAKEUPAT UTC)... "

# infinite loop
while :; do
    now="$(date --utc +%H:%M)"
    if [[ "$now" = "$WAKEUPAT" ]]; then
	echo " it's $now, time to wake up!"
	echo " .. cleaning up .."
	make clean
	echo " .. running nightly checks .."
        ctest -D Nightly # also performs SVN update
	#ctest -D Experimental # does no SVN update before building
	echo " .. running memory usage checks .."
	ctest -D MemoryCheck
	echo " .. my job is done done for now, going to sleep on $(date) .. "
	sleep 59
    fi
    sleep 30
done
