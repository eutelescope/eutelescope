#!/bin/bash

###################################################################
#########
#########    Millipede output to lcio conversion 
#########
#########    change: 28.05.2013   Denys Lontkovskyi
#########		 Introduce basic functionality
#########    Last change: 11.10.2013   Denys Lontkovskyi
#########		 Intoduce possibility to generate gear file
###################################################################
#
#	The script takes milipede steering files produced by 
#	EUTelProcessorGBLTracking during alignment step and output
#	of pede (mille.res) and converts alignment constants 
#	into LCIO and/or GEAR collection.
#
#	NOTE: Steeting file must have comments in special format
#	understood by pede2lcio executable. Be carefull with 
#	modifications.
#
###################################################################
# Program constants 
E_BADARGS=65
E_BADFILE=66

# 
function usage {
  echo "Usage: $0 mille-steer.txt mille.res alignment.slcio oldGear.xml newGear.xml"
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

function is_file_gear {
        local file="$1"
        [[ $file == *.xml ]] && return 1 || return 0
}

readonly -f usage

##################################################################
############# Script entry point
##################################################################

# Check number of arguments supplied
if !( [ $# -eq 3 ] || [ $# -eq 5 ] ); then
  usage
fi

steering_file=$1
outres_file=$2
lcio_file=$3
oldgear_file=$4
newgear_file=$5

# Check if supplied files exist
if ( is_file_bad "$steering_file" )
then
 echo "File $steering_file not found. Terminate..."
 exit $E_BADFILE
fi

# Check extension of the first supplied file (require mille steering)
if ( is_file_millesteer "$steering_file" )
then
 echo "File $steering_file in not a millepede steering file. Terminate..."
 exit $E_BADFILE
fi

# Check extension of the first supplied file (require mille steering)
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
sort -n -u $steering_file | grep '^[[:blank:]]*[[:digit:]]' > tmp_parsepede_file1
sort -n -u $outres_file   | grep '^[[:blank:]]*[[:digit:]]'   > tmp_parsepede_file2
join tmp_parsepede_file1 tmp_parsepede_file2 | grep '^[0-9]' | awk '{ printf "%-15s%-15s%-25s%-15s%-15s%-15s\n",$1,$5,$6,$7,$8,$11}' | tee out.pede2lcio


# Convert to lcio collection file
if [ $# == 3 ]; then 
 echo  pede2lcio out.pede2lcio $lcio_file
 pede2lcio out.pede2lcio $lcio_file


 rm -f out.pede2lcio
 rm -f file1 file2

fi

# Convert to lcio collection file and new GEAR file
if [ $# == 5 ]; then 

 # Check extension of supplied old GEAR file
 if ( is_file_gear "$oldgear_file" )
 then
  echo "File $oldgear_file in not a gear file. Terminate..."
  exit $E_BADFILE
 fi

 echo pede2lcio -g out.pede2lcio $lcio_file $oldgear_file $newgear_file
 pede2lcio -g out.pede2lcio $lcio_file $oldgear_file $newgear_file

 rm -f out.pede2lcio
 rm -f file1 file2

fi
rm tmp_parsepede_file1 tmp_parsepede_file2

