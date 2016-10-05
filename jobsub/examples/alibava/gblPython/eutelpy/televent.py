'''
Created on Jun 8, 2015

@author: kleinwrt
'''

## \file
# telescope event 
#
# read from text or lcio file

from eutelpy.telhit import TelHit
# import matplotlib.pyplot as plt          # remove to get rid of X11 error
import numpy as np

## telescope event
#
class TelEvent(object):
  
  ## constructor.
  #
  # @param[in] detector    detector
  # @param[in] beamEnergy  (test)beam energy
  # @param[in] bField  magnetic field
  #
  def __init__(self, detector, beamEnergy, bField):
    ## detector
    self.__det = detector
    ## run number
    self.__runNumber = 0
    ## event number
    self.__eventNumber = 0
    ## q/p
    self.__qbyp = 1. / beamEnergy
    ## magnetic field
    self.__bField = bField
    ## deflection strength
    self.__bfac = (-0.0003 * bField[1], 0.0003 * bField[0], 0.)  # -cBxZ
    ## hits
    self.__hits = []
    ## tracks
    self.__tracks = []
    ## trajectories
    self.__traj = []
    ## valid event
    self.__valid = False
    ## momentum
    self.__Pseg = []
    self.__Pgbl = []


  ## Read event from file.
  #
  # @param[in] telFile  TelFile (text or lcio)
  #  
  def read(self, telFile):
    # read event
    header, hits = telFile.read()
    if header is not None:
      self.__runNumber, self.__eventNumber = header
      self.__valid = True
      #
      for h in hits:
        self.__hits.append(TelHit(self.__det, h[0], h[1], (h[2], h[3], h[4]), h[5]))
   
  ## Event is valid?
  def isValid(self):
    return self.__valid
  
  ## Get run number.
  def getRunNumber(self):
    return self.__runNumber
  
  ## Get run number.
  def getEventNumber(self):
    return self.__eventNumber    

  ## Get run, event number.
  def getHeader(self):
    return self.__runNumber, self.__eventNumber   
 
  ## Get number of hits.
  def getNumHits(self):
    return len(self.__hits)  

  ## Get number of DUT hits.
  def getNumDUThits(self):
    count = 0
    for h in self.__hits:
      if h.getPlane() > 5:
        count += 1
    return count  
  
  ## Get hits.
  def getHits(self):
    return self.__hits  

  ## Get number of tracks.
  def getNumTracks(self):
    return len(self.__tracks)  
     
  ## Dump.
  def dump(self):
    print " run ", self.__runNumber, " event ", self.__eventNumber
    for h in self.__hits:
      h.dump()
 
  ## Plot hits.
  def plot(self, show=True):
    X1 = []; Y1 = []; Z1 = []
    X2 = []; Y2 = []; Z2 = []
    for h in self.__hits:
      l = h.getPlane()
      p = h.getPos()
      if l < 6:
        X1.append(p[0])
        Y1.append(p[1])
        Z1.append(p[2])
      else:    
        X2.append(p[0])
        Y2.append(p[1])
        Z2.append(p[2])
    fig = plt.figure(1)
    fig.suptitle("run %i event %i " % (self.__runNumber, self.__eventNumber))
    plt.subplot(211)      
    plt.plot(Z1, X1, 'r+')
    plt.plot(Z2, X2, 'b*')
    plt.plot([-50], [0.], '')
    plt.title('X vs Z')
    plt.subplot(212)                
    plt.plot(Z1, Y1, 'r+')  
    plt.plot(Z2, Y2, 'b*')  
    plt.plot([-50], [0.], '')
    plt.title('Y vs Z')
    if show: 
      plt.show()        
  
  ## Find tracks.    
  #
  # @param[in] finder   track finder
  # @param[in] zMag     start (in z) of magnetic field
  #
  def findTracks(self, finder, zMag):   
    finder.findTracks(self.__hits, self.__qbyp, self.__bfac, zMag)
    finder.matchDUT()
    self.__tracks = finder.getSegments()    
      
  ## Fit with GBL.
  #
  # @param[in] binaryFile   Millepede-II binary file
  # @param[in] consBL       beam line constraint
  #
  def fitGBL(self, binaryFile, consBL=None):
    # fit segments
    for track in self.__tracks:
      self.__traj.append(track.fitSegment(self.__det, self.__qbyp, self.__bField, binaryFile, consBL))

  ## Analyze event.
  #
  # Look at DUT hits
  #
  # @param[in] rootTree  Root TTree
  #
  def analyze(self, rootTree=None):
    for it, t in enumerate(self.__tracks):
      hits = t.getHits()
      for p, h in hits.iteritems():
	if p == 0:
	  globalPos_0 = np.array(h.getPos())
	if p == 1:
	  globalPos_1 = np.array(h.getPos())
	if p == 2:
	  globalPos_2 = np.array(h.getPos())
	if p == 3:
	  globalPos_3 = np.array(h.getPos())
	if p == 4:
	  globalPos_4 = np.array(h.getPos())
	if p == 5:
	  globalPos_5 = np.array(h.getPos())
        if p > 5:
          plane = self.__det[p]
          # from plane
          rot = plane.getRotation()
          # from hit
          z = h.getZ()
	  globalPosDUT = np.array(h.getPos())
          locMeas = plane.transformGlobalToLocal(np.array(h.getPos()))
          # from seed track
          position = t.getPosition(z)
          direction = t.getDirection(z)
          locDir = np.dot(rot.T, np.array(direction))
          locPos = plane.transformGlobalToLocal(np.array(position))
          # corrections from (fitted) GBL trajectory
          locCorr = self.__traj[it].getResults(h.getLabel())[0]
          locSlope = (locDir[0] / locDir[2] + locCorr[1], locDir[1] / locDir[2] + locCorr[2])
          locPos = (locPos[0] + locCorr[3], locPos[1] + locCorr[4])
          #print " analyze %i %i %i %i %.4f %.4f %.4f %.4f %.4g %.4g" % (self.__runNumber, self.__eventNumber, h.getIndex(), p, locMeas[0], locMeas[1], locPos[0], locPos[1], locSlope[0], locSlope[1])
          #print "   channelVec ", h.getCluster()
          if self.__bField[0] == 0 and self.__bField[1] == 0 and self.__bField[2] == 0: 
            pseg = 0 
            pgbl = 0
          else:
            pseg = self.__Pseg[it] 
            pgbl = self.__Pgbl[it] 
          if rootTree is not None:
            rootTree.fill(self.__runNumber, self.__eventNumber, pseg, pgbl, globalPos_0, globalPos_1, globalPos_2, globalPos_3, globalPos_4, globalPos_5, globalPosDUT, h.getIndex(), p, locPos, locSlope, h.getCluster())

  ## Analyze event.
  #
  # Look at DUT scattering
  #
  # @param[in] hists  histograms
  #
  def analyzeScat(self, hists=None):
    for t in self.__traj:
      # local corrections at first point          
      locCorr, locCov = t.getResults(1)
      if locCorr.shape[0] > 6:
        print " scatPar ", locCorr[5], locCorr[6], locCov[5, 5], locCov[6, 6]
        if hists is not None:
          hists.addEntry("xScat", locCorr[5])
          hists.addEntry("yScat", locCorr[6])

  ## Analyze event.
  #
  # Look at residuals
  #
  def analyzeRes(self):
    for t in self.__traj:
      # check residuals
      for i in range(t.getNumPoints()):
        numData, aResiduals, aMeasErr, aResErr, aDownWeight = t.getMeasResults(i + 1)
        for j in range(numData):
          print " measRes " , i, j, aResiduals[j], aMeasErr[j], aResErr[j], aDownWeight[j]
        numData, aResiduals, aMeasErr, aResErr, aDownWeight = t.getScatResults(i + 1)
        for j in range(numData):
          print " scatRes " , i, j, aResiduals[j], aMeasErr[j], aResErr[j], aDownWeight[j]

  ## Analyze event.
  #
  # Look at measured momentum
  #
  # @param[in] hists  histograms
  #
  def analyzeMom(self, hists):
    bfac2 = self.__bfac[0] * self.__bfac[0] + self.__bfac[1] * self.__bfac[1]
    if bfac2 <= 0.:
      return
    #  Bfield on
    for it, t in enumerate(self.__tracks):
      # q/p of segment (fron slope difference of triplets)
      curv = t.getCurvature()
      qbyp = (curv[0] * self.__bfac[0] + curv[1] * self.__bfac[1]) / bfac2
      # correction from GBL fit
      locCorr = self.__traj[it].getResults(1)[0]
      dqbyp = locCorr[0]
      #print " q/p ", it, qbyp, dqbyp
      self.__Pseg.append(1. / qbyp)
      self.__Pgbl.append(1. / (qbyp + dqbyp))

      if hists is not None:
          hists.addEntry("Pseg", self.__Pseg[it])
          hists.addEntry("Pgbl", self.__Pgbl[it])
  
