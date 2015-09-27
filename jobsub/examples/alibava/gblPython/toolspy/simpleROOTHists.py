'''
ROOT histograms

Needs:
 - $ROOTSYS/lib:

Created on Aug 31, 2015

@author: kleinwrt
'''

## \file
# Simple histograms with ROOT.
 
from ROOT import TFile, TH1F  #@UnresolvedImport
 
## Simple histograms.
#
#  Based on matplotlib and scipy.
#
class SimpleHists(object):

  ## Constructor.
  #
  # @param[in] labels  labels
  # 
  def __init__(self, labels):
    ## labels
    self.__labels = labels
    ## histograms
    self.__hist = {}
    for l in labels:
      self.__hist[l] = []
      
  ## Add entry
  #
  # @param[in] label  label
  # @param[in] x      value
  #
  def addEntry(self, l, x):
    if l in self.__hist:
      self.__hist[l].append(x)
 
  ## Dump histograms      
  def dump(self):
    for l in self.__labels:
      print " hist ", l, self.__hist[l]
 
  ## Plot histograms      
  def plot(self):
    print " Plotting of histograms with ROOT is not implemented, use 'save(fileName)' !"
      
  ## save histograms (to ROOT)
  #
  # @param[in] fileName  file name
  #
  def save(self, fileName):
    tFile = TFile(fileName, 'recreate')
    for i, l in enumerate(self.__labels):
      entries = len(self.__hist[l])
      if entries > 1:
        h = TH1F('h' + str(i), l, 100, min(self.__hist[l]), max(self.__hist[l]))
        for x in self.__hist[l]:
          h.Fill(x)
        h.Write()    
    tFile.Close()
