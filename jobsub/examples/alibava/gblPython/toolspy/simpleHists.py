'''
Created on May 18, 2015

@author: kleinwrt
'''

## \file
# Simple histograms.

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import norm

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
    for l in self.__labels:
      entries = len(self.__hist[l])
      if entries > 1:
        plt.hist(self.__hist[l], 100, histtype='step')
        mean, rms = norm.fit(self.__hist[l]) 
        plt.title(' %s, entries=%i mean=%.4g rms=%.4g' % (l, entries, mean, rms))
        plt.show()
      
  ## save histograms (to pdf)
  #
  # @param[in] fileName  file name
  #
  def save(self, fileName):
    pdf = PdfPages(fileName)
    for l in self.__labels:
      entries = len(self.__hist[l])
      if entries > 1:
        plt.hist(self.__hist[l], 100, histtype='step')
        mean, rms = norm.fit(self.__hist[l])
        plt.title(' %s, entries=%i mean=%.4g rms=%.4g' % (l, entries, mean, rms))
        pdf.savefig()
        plt.close()
    pdf.close()
        
