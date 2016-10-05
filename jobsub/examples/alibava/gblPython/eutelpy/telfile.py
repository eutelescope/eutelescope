'''
Created on Jun 12, 2015

@author: kleinwrt
'''

## \file
# telescope file 

from eutelpy.tellcio import TelLcioFile

## telescope file (lcio or text).
class TelFile(object):
  
  ## constructor.
  #
  # @param[in] fileName    file name
  # @param[in] collection  hit collection for lcio, None for text file
  #
  def __init__(self, fileName, collection=None):
    ## file name
    self.__fileName = fileName
    ## collection
    self.__collection = collection 
    ## file
    self.__file = None
    
  ## open file
  def open(self):
    if self.__collection is None:
      print "Opening text file", self.__fileName
      self.__file = open(self.__fileName)
    else:
      print "Opening lcio file", self.__fileName
      self.__file = TelLcioFile()
      self.__file.open(self.__fileName)  
    print
    
  ## close file
  def close(self):
    self.__file.close()
    
  ## read file
  def read(self):
    # header
    header = None
    # list of hits
    hits = []
      
    if self.__collection is None:
      # read from text file
      index = 0
      while True:
        line = self.__file.readline()
        if line == "":
          break
        # line read
        fields = line.split()
        if len(fields) < 4:
          # last line with event information: run ,event number
          header = (int(fields[0]), int(fields[1]))      
          break
        else:
          # line with hit information: planeID, (local) x, y, z
          hits.append((index, int(fields[0]), float(fields[1]), float(fields[2]), float(fields[3]), []))
          index += 1 
    else:
      # read from lcio file
      event = self.__file.readEvent()
      if event:
        header = (event.getRunNumber(), event.getEventNumber())    
        # analyze hit collection
        for collectionName, collection in event:
          if collectionName == self.__collection:
            index = 0
            for h in collection:
              # get hit info
              pos = h.getPosition()
              id0 = int(h.getCellID0())
              cluster = []
              # for DUT get cluster
              if id0 > 5:                
                channelVec = h.getRawHits()[0].getChargeValues()
                cluster = [ c for c in channelVec]
              hits.append((index, id0, pos[0], pos[1], pos[2], cluster))
              index += 1
          
    return header, hits
