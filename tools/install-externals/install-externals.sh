#!/bin/bash

# check number of arguments and display help text if not as expected
if [ ! $# -eq 1 ]
then
    echo "This script installs external dependencies for EUTelescope."
    echo "Please provide the directory path in which to install as argument e.g. '$0 $EUTELESCOPE/external'"
    exit
fi

if [ ! -d $1 ]
then
    echo "This script installs external dependencies for EUTelescope."
    echo "Please provide *an existing* directory path in which to install e.g. '$0 $EUTELESCOPE/external'"
    exit
fi

extpath=$1
pypath="$extpath/pymodules" # where to store external python modules
tempdir="$extpath/packages.tmp" # where to store the .tar.gz files and company

if [ ! -d $tempdir ]; then mkdir $tempdir; fi;
if [ ! -d $pypath ]; then mkdir $pypath; fi;

# argparse -- needed for python scripts on SL5 machines
argparsefile="argparse-1.2.1.tar.gz"
wget --no-check-certificate -P "$tempdir" 'http://argparse.googlecode.com/files/'$argparsefile
if [ ! -d "$pypath/argparse" ]; then mkdir "$pypath/argparse"; fi;
echo "Extracting tar archive..."
tar --strip-components 1 -C "$pypath/argparse" -xzf "$tempdir/$argparsefile" 
echo "... done with argparse python module"

# rootpy -- useful for interfacing python and ROOT
rootpyversion="f4600dfd30073a0702c39d4db2c7f307eff953f9" # this is after 0.7.1; should point to a (later) release when provided!
wget --no-check-certificate --output-document="$tempdir/${rootpyversion}.tar.gz" 'https://github.com/rootpy/rootpy/archive/'${rootpyversion}.tar.gz
if [ ! -d "$pypath/rootpy-dev" ]; then mkdir "$pypath/rootpy-dev"; fi;
echo "Extracting tar archive..."
tar --strip-components 1 -C "$pypath/rootpy-dev" -xzf "$tempdir/${rootpyversion}.tar.gz" 
echo "... done with rootpy python module"

# Eigen -- linear algebra package
eigenversion="3.2.2"
wget --no-check-certificate --output-document="$tempdir/${eigenversion}.tar.gz" 'http://bitbucket.org/eigen/eigen/get/'${eigenversion}.tar.gz
if [ ! -d "$extpath/Eigen" ]; then mkdir "$extpath/Eigen"; fi;
echo "Extracting tar archive..."
tar --strip-components 1 -C "$extpath/Eigen" -xzf "$tempdir/${eigenversion}.tar.gz" 
echo "... done with Eigen library"

# cmspixeldecoder -- needed for decoding of CMS pixel tracker data
cmxpixeldecoderversion="master" # this is HEAD; should point to a (later) release when provided!
wget --no-check-certificate --output-document="$tempdir/${cmxpixeldecoderversion}.tar.gz" 'https://github.com/simonspa/CMSPixelDecoder/archive/'${cmxpixeldecoderversion}.tar.gz
if [ ! -d "$extpath/CMSPixelDecoder" ]; then mkdir "$extpath/CMSPixelDecoder"; fi;
echo "Extracting tar archive..."
tar --strip-components 1 -C "$extpath/CMSPixelDecoder" -xzf "$tempdir/${cmxpixeldecoderversion}.tar.gz" 
echo "... done with cmxpixeldecoder library"

# clean up
rm $tempdir/*
rmdir $tempdir
