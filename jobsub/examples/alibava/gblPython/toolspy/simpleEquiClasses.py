'''
Created on Aug 22, 2013

@author: kleinwrt
'''
## \file
# Simple wquivalence classes.

# ######################################################################       

## Equivalence classes
#
# From indices of matching (list) objects.
#
# Adopted from (H1) Fortran version by V. Blobel
#
class SimpleEquivalenceClasses(object):

  ## Constructor.
  def __init__(self):
    ## indices
    self._index = {}
    ## matches
    self._matches = []
    
  ## Add index.
  def addIndex(self, aIndex):
    self._index[aIndex] = aIndex
    
  ## Add matching (index) pair.
  def addMatch(self, aIndexPair):
    self._matches.append(aIndexPair)
    
  ## Get classes.
  def getClasses(self, debug=False):
    #
    # determine equivalence classes - inline algorithm
    #
    if debug:
      print " Indices ", self._index
      print " matches ", self._matches
    for i1, i2 in self._matches:
      # trace first element up to is ancestor
      while self._index[i1] != i1:
        i1 = self._index[i1]
      # trace second element up to is ancestor
      while self._index[i2] != i2:
        i2 = self._index[i2]
      if i1 != i2:
        self._index[i1] = i2 
    # final sweep-up to highest ancestor
    for i1 in self._index.iterkeys():
      while self._index[i1] != self._index[self._index[i1]]:
        self._index[i1] = self._index[self._index[i1]]
    if debug:
      print " Ancestors ", self._index    
    classes = {}
    for i, a in self._index.iteritems():
      if a not in classes:
        classes[a] = []
      classes[a].append(i)
    if debug:
      print " Classes ", classes 
    return classes
