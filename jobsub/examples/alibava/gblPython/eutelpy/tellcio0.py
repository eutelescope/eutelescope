'''
Created on Jun 12, 2015

@author: kleinwrt
'''

## \file
# DUMMY telescope lcio file 


## dummy telescope lcio file.
class TelLcioFile(object):
  
  ## constructor.
  #
  #
  def __init__(self):
    print "LCIO files not supported !!!"
    
  ## open file
  #
  # @param[in] fileName  file name
  #
  def open(self, fileName):
    print "   cannot read from ", fileName
    
  ## close file
  def close(self):
    pass
    
  ## read event
  def readEvent(self):
    return None
    
