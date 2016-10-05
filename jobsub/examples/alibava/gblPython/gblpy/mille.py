'''
Input/output of MP-II binary records.

Created on Aug 1, 2011

@author: kleinwrt
'''

## \file
# Millepede-II (binary) record

import array, math

## Millepede-II (binary) record.
#  
#  Containing information for local (track) and global fit.
# 
#  The data blocks are collected in two arrays, a real array and
#  an integer array, of same length.  The content of the arrays::
#
#\verbatim
#         real array              integer array    
#     0   0.0                     error count (this record)  
#     1   RMEAS, measured value   0                            __iMeas   -+
#     2   local derivative        index of local derivative               |
#     3   local derivative        index of local derivative               |
#     4    ...                                                            | block
#         SIGMA, error (>0)       0                            __ iErr    |
#         global derivative       label of global derivative              |
#         global derivative       label of global derivative             -+
#         RMEAS, measured value   0                            __position
#         local derivative        index of local derivative
#         local derivative        index of local derivative
#         ...
#         SIGMA, error            0
#         global derivative       label of global derivative
#         global derivative       label of global derivative
#         ...
#         global derivative       label of global derivative   __recLen
#\endverbatim
#
class MilleRecord(object):

  ## Create MP-II binary record.
  #
  def __init__(self):
    ## position in record, usually start of next data block; int
    self.__position = 1 
    ## number of data blocks in record; int
    self.__numData = 0  
    ## record length; int
    self.__recLen = 0 
    ## position of value in current data block; int 
    self.__iMeas = 0    
    ## position of error in current data block; int
    self.__iErr = 0     
    ## array with markers (0) and labels; array(int32)
    self.__inder = array.array('i') 
    ## array with values, errors and derivatives; (float32)
    self.__glder = array.array('f') 
 
  ## Add data block to (end of) record.
  #  
  #  @param dataList list with measurement, error, labels and derivatives; list 
  #  
  def addData(self, dataList):
    if (self.__numData == 0): # first word is error counter
      self.__inder.append(0)
      self.__glder.append(0.) 
    self.__numData += 1     
 
    aMeas, aPrec, indLocal, derLocal, labGlobal, derGlobal = dataList
    self.__inder.append(0)
    self.__glder.append(aMeas)
    self.__inder.fromlist(indLocal)
    self.__glder.fromlist(derLocal)    
    self.__inder.append(0)
    self.__glder.append(1.0 / math.sqrt(aPrec)) # convert to error
    self.__inder.fromlist(labGlobal)
    self.__glder.fromlist(derGlobal)

  ## Get data block from current position in record.
  #  
  #  @return list with measurement, error, labels and derivatives; list
  #    
  def getData(self):
    aMeas = self.__glder[self.__iMeas]
    indLocal = []
    derLocal = []
    for i in range(self.__iMeas + 1, self.__iErr):
      indLocal.append(self.__inder[i])
      derLocal.append(self.__glder[i])
    aPrec = 1.0 / self.__glder[self.__iErr] ** 2 # convert to precision 
    indGlobal = []
    derGlobal = []
    for i in range(self.__iErr + 1, self.__position):
      indGlobal.append(self.__inder[i])
      derGlobal.append(self.__glder[i])    
    return aMeas, aPrec, indLocal, derLocal, indGlobal, derGlobal   
    
 
  ## Print record. 
  def printRecord(self):
    print " MilleRecord, len: ", len(self.__inder)
    print self.__inder
    print self.__glder
 
  ## Write record to file.
  #  
  #  @param aFile (binary) file
  #
  def writeRecord(self, aFile):
    header = array.array('i') # header with number of words
    header.append(len(self.__inder) * 2)
    header.tofile(aFile)
    self.__glder.tofile(aFile)
    self.__inder.tofile(aFile)
 
  ## Read record from file.
  #  
  #  @param aFile (binary) file
  #    
  def readRecord(self, aFile):
    header = array.array('i') # header with number of words
    header.fromfile(aFile, 1)
    self.__recLen = header[0] / 2
    self.__glder.fromfile(aFile, self.__recLen)
    self.__inder.fromfile(aFile, self.__recLen)
 
  ## Locate next data block.
  #  
  #  @return next block exists; bool 
  #   
  def moreData(self):
    if (self.__position < self.__recLen):
      while (self.__position < self.__recLen and self.__inder[self.__position] != 0):
        self.__position += 1
      self.__iMeas = self.__position 
      self.__position += 1
      while (self.__position < self.__recLen and self.__inder[self.__position] != 0):
        self.__position += 1
      self.__iErr = self.__position 
      self.__position += 1
      while (self.__position < self.__recLen and self.__inder[self.__position] != 0):
        self.__position += 1
      self.__numData += 1  
      return True
    else:
      return False

  ## Get special data tag from block.
  #  
  #  @return tag or -1 for ordinary data block; int    
  #
  def specialDataTag(self):
    aTag = -1
    if (self.__iMeas + 1 == self.__iErr and self.__glder[self.__iErr] < 0):
      aTag = int(-self.__glder[self.__iErr] * 10. + 0.5) % 10
    return aTag
