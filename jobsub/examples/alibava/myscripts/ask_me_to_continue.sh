#!/bin/bash

input="n"
qstat

while [[ $(bash -c 'read -e -t 130 -p "Do you want to continue?[y/n]: " tmp; echo $tmp') != "y" ]]; do
echo " "
#echo "Sleeping for 30 sec"
#sleep 30 
qstat

done

