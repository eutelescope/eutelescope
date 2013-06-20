#!/bin/bash

###################################################################
#########
#########    Millipede output to lcio conversion 
#########
#########    Last change: 28.05.2013   Denys Lontkovskyi
#########		 Introduce basic functionality
###################################################################
#
#	The script takes milipede steering files produced by 
#	EUTelProcessorGBLTracking during alignment step and output
#	of pede (mille.res) and converts alignment constants 
#	into LCIO collection.
#
#	NOTE: Steeting file must have comments in special format
#	understood by pede2lcio executable. Be carefull with 
#	modifications.
#
###################################################################
# Program constants 
EXPECTED_ARGS=3
E_BADARGS=65
E_BADFILE=66

# 
function usage {
  echo "Usage: $0 mille-steer.txt mille.res alignment.slcio"
  exit $E_BADARGS
}


function is_file_bad {
	local file="$1"
	[[ -f "$file" ]] && return 1 || return 0
}

function is_file_millesteer {
	local file="$1"
	[[ $file == *.txt ]] && return 1 || return 0
}

function is_file_milleout {
	local file="$1"
	[[ $file == *.res ]] && return 1 || return 0
}
readonly -f usage

##################################################################
############# Script entry point
##################################################################

# Check number of arguments supplied
if [ $# -ne $EXPECTED_ARGS ]; then
  usage
fi

steering_file=$1
outres_file=$2
lcio_file=$3

# Chech if supplied files exist
if ( is_file_bad "$steering_file" )
then
 echo "File $steering_file not found. Terminate..."
 exit $E_BADFILE
fi

# Chech extension of the first supplied file (require mille steering)
if ( is_file_millesteer "$steering_file" )
then
 echo "File $steering_file in not a millepede steering file. Terminate..."
 exit $E_BADFILE
fi

# Chech extension of the first supplied file (require mille steering)
if ( is_file_milleout "$outres_file" )
then
 echo "File $outres_file in not a millepede output file. Terminate..."
 exit $E_BADFILE
fi

if ( is_file_bad "$outres_file" )
then
 echo "File $outres_file not found. Terminate..."
 exit $E_BADFILE
fi

# Merge MILLIPEDE steering and output files
sort -g $steering_file > file1
sort -g $outres_file   > file2
join file1 file2 | grep '^[0-9]' | awk '{ printf "%-10s%-4s%-20s%-15s%-15s%-15s\n",$1,$5,$6,$7,$8,$11}' | tee out.pede2lcio

# Convert to lcio collection file
pede2lcio out.pede2lcio $lcio_file
rm -f out.pede2lcio
rm -f file1 file2

# WARNING
echo ""
echo "ATTENTION!ATTENTION!ATTENTION!ATTENTION!ATTENTION!"
echo "Please, check the consitency of conversion step!"
echo "ATTENTION!ATTENTION!ATTENTION!ATTENTION!ATTENTION!"
echo ""
