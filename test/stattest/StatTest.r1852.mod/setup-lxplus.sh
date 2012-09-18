#!/bin/bash
sys=x86_64-slc5-gcc43-opt
gccver=4.3.2
pyver=2.6.5p2
rootver=5.34.00

basedir=/afs/cern.ch/sw/lcg/external

#gcc
source ${basedir}/gcc/${gccver}/${sys}/setup.sh
#python
export PATH=${basedir}/Python/${pyver}/${sys}/bin:$PATH
export LD_LIBRARY_PATH=${basedir}/Python/${pyver}/${sys}/lib:$LD_LIBRARY_PATH
#ROOT
source /afs/cern.ch/sw/lcg/app/releases/ROOT/${rootver}/${sys}/root/bin/thisroot.sh