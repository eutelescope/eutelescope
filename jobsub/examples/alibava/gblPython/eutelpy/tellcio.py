'''
Read telescope hits from LCIO file with pyLCIO and write them to text file

Needs:
 - ROOT dictionary for lcio
 - LCIO shell variable
 - $ROOTSYS/lib:$LCIO/src/python in PYTHONPATH
 
Created on Jun 12, 2015

@author: kleinwrt
'''

## \file
# telescope lcio file 

from pyLCIO import IOIMPL  #@UnresolvedImport

## telescope hit.
class TelLcioFile(object):
  
  ## constructor.
  #
  #
  def __init__(self):
    ## reader 
    self.__reader = IOIMPL.LCFactory.getInstance().createLCReader()
    
  ## open file
  #
  # @param[in] fileName  file name
  #
  def open(self, fileName):
    # open LCIO file
    self.__reader.open(fileName)
    
  ## close file
  def close(self):
    # close KCIO file  
    self.__reader.close()
    
  ## read event
  def readEvent(self):
    # read next event
    return self.__reader.readNextEvent()
    
