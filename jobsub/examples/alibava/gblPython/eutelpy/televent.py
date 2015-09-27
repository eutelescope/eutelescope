'''
Created on Jun 8, 2015

@author: kleinwrt
'''

## \file
# telescope event 
#
# read from text or lcio file

from eutelpy.telhit import TelHit
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
  #
  def fitGBL(self, binaryFile):
    # fit segments
    for track in self.__tracks:
      self.__traj.append(track.fitSegment(self.__det, self.__qbyp, self.__bField, binaryFile))

  ## Analyze event.
  #
  # Look at DUT hits
  #
  # @param[in] rootTree  Root TTree
  #
  def analyze(self, rootTree=None):
    for it, t in enumerate(self.__tracks):
      hits = t.getHits()
      for p in range(6, 9):
        if p in hits:
          plane = self.__det[p]
          # from plane
          rot = plane.getRotation()
          h = hits[p]
          # from hit
          z = h.getZ()
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
          if rootTree is not None:
            rootTree.fill(self.__runNumber, self.__eventNumber, h.getIndex(), h.getPlane(), locPos, locSlope, h.getCluster())
          
