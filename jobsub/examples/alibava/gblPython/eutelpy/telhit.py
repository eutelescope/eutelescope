'''
Created on Jun 8, 2015

@author: kleinwrt
'''

## \file
# telescope event 

import numpy as np

## telescope hit.
class TelHit(object):
  
  ## constructor.
  #
  # @param[in] det      detector
  # @param[in] index    index
  # @param[in] planeID: plane ID
  # @param[in] pos      (local) position 
  # @param[in] cluster  cluster information (optional)
  #
  def __init__(self, det, index, planeID, pos, cluster):
    ## detector plane ID
    self.__id = planeID
    # local position
    locPos = np.array(pos)
    ## global position
    self.__pos = det[self.__id].transformLocalToGlobal(locPos).tolist()
    ## used flag 
    self.__used = False 
    ## matches
    self.__matches = 0
    ## GBL label
    self.__label = 0
    ## index (in input hit collection)
    self.__index = index
    ## cluster information
    self.__cluster = cluster

  ## Set GBL label.
  #
  # @param[in] l  label
  #
  def setLabel(self, l):
    self.__label = l  

  ## Get GBL label.
  def getLabel(self):
    return self.__label  
 
  ## Get GBL label.
  def getIndex(self):
    return self.__index  
            
  ## Add match
  def addMatch(self):
    self.__matches += 1
    
  ## Get matches
  def getMatches(self):
    return self.__matches

  ## Set use flag.
  def setUsed(self):
    self.__used = True

  ## Get use flag.  
  def getUsed(self):
    return self.__used    
 
  ## Get X coordinate.
  def getX(self):
    return self.__pos[0]
  
  ## Get Y coordinate.
  def getY(self):
    return self.__pos[1]
    
  ## Get Z coordinate.
  def getZ(self):
    return self.__pos[2]

  ## Get coordinates.
  def getPos(self):
    return self.__pos
  
  ## Get plane.
  def getPlane(self):
    return self.__id
    
  ## Get cluster information
  def getCluster(self):
    return self.__cluster  
        
  ## Dump.
  def dump(self):
    print " id ", self.__id, " pos ", self.__pos
