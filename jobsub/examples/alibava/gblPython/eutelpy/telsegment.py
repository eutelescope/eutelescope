'''
Created on Jun 8, 2015

@author: kleinwrt
'''

## \file
# telescope segment

#import matplotlib.pyplot as plt          # remove to get rid of X11 error
import numpy as np
import math
from eutelpy.telsystem import globalZXYsystem
from gblpy.gblfit import GblPoint, GblTrajectory

## Segment (track candidate).
class TelSegment(object): 
    
  ## Constructor.
  #
  # Multiple hits per plane are not implemented.
  #
  # @param[in] header   event header
  # @param[in] hits     list of hits
  # @param[in] curv     curvatures in X(Z), Y(Z)
  #
  def __init__(self, header, hits, curv=(0., 0.)): 
    ## (event) header
    self.__hdr = header     
    ## hits
    self.__hits = {}
    for h in hits:
      self.__hits[h.getPlane()] = h
    # straight line with first and last hit
    ## offsets (position of first hit and z of last)
    self.__offset = (hits[0].getX(), hits[0].getY(), hits[0].getZ(), hits[-1].getZ())  
    ## (average) slope
    dz = self.__offset[3] - self.__offset[2]
    self.__slope = ((hits[-1].getX() - self.__offset[0]) / dz, (hits[-1].getY() - self.__offset[1]) / dz)
    ## curvature
    self.__curv = (curv[0], curv[1])
    ## matches
    self.__matches = 0
    
    #self.dump()                      
  
  ## Dump
  def dump(self):
    print " segment ", len(self.__hits), self.__offset[0], self.__offset[1], self.__offset[2], \
                       self.__offset[3], self.__slope[0], self.__slope[1], self.__curv[0], self.__curv[1]  
   
  ## Add match
  def addMatch(self):
    self.__matches += 1
    
  ## Get matches
  def getMatches(self):
    return self.__matches       
 
  ## Reset matches
  def resetMatches(self):
    self.__matches = 0       
    
  ## Add DUT.
  #
  # @param[in] hit  DUT hit
  #
  def addDUT(self, hit):
    self.__hits[hit.getPlane()] = hit
  
  ## Get hits.  
  def getHits(self):
    return self.__hits
    
  ## Get Z of track center
  def getZCenter(self):
    return 0.5 * (self.__offset[2] + self.__offset[3])  
    
  ## Get prediction.
  #
  # @param[in] z   z position
  #
  def getPrediction(self, z):
    # parabolic in z
    dz1 = z - self.__offset[2]  # distance (in Z) to first hit
    dz2 = z - self.__offset[3]  # distance (in Z) to last hit
    px = self.__offset[0] + dz1 * self.__slope[0] + 0.5 * dz1 * dz2 * self.__curv[0] 
    py = self.__offset[1] + dz1 * self.__slope[1] + 0.5 * dz1 * dz2 * self.__curv[1]
    return (px, py)
    
  ## Get position.
  #
  # @param[in] z   z position
  #
  def getPosition(self, z):
    px, py = self.getPrediction(z)
    return (px, py, z)
      
  ## Get slope.
  #
  # @param[in] z   z position
  #
  def getSlope(self, z):
    dz = z - 0.5 * (self.__offset[2] + self.__offset[3])
    return (self.__slope[0] + dz * self.__curv[0], self.__slope[1] + dz * self.__curv[1])      

  ## Get curvature.
  def getCurvature(self):
    return self.__curv      
      
  ## Get direction.
  #
  # @param[in] z   z position
  #
  def getDirection(self, z):
    slopeX, slopeY = self.getSlope(z)
    norm = math.sqrt(slopeX * slopeX + slopeY * slopeY + 1.)
    return (slopeX / norm, slopeY / norm, 1. / norm) 
  
  ## Check residuals (pre-alignment).
  def checkRes(self):
    # require complete segment
    if self.__hits[0].getLayer() != 0 or self.__hits[-1].getLayer() != 5:
      return

    for h in self.__hits[1:-1]:
      x = h.getX(); y = h.getY(); z = h.getZ()
      dz1 = z - self.__offset[2]  # distance (in Z) to first hit
      dz2 = z - self.__offset[3]  # distance (in Z) to last hit
      px = self.__offset[0] + dz1 * self.__slope[0] + 0.5 * dz1 * dz2 * self.__curv[0] 
      py = self.__offset[1] + dz1 * self.__slope[1] + 0.5 * dz1 * dz2 * self.__curv[1]
      print " res ", h.getLayer(), x, px, y, py, self.__slope[0], self.__slope[1]
      
  ## Fit segment.
  #
  #  Fit segment with GBL, parabolic in Z.
  #  As local (fit) systems the global ZXY system is used. 
  #
  # @param[in] det     detector (planes)
  # @param[in] qbyp    assumed Q/P (e.g. from beam energy)
  # @param[in] bField  magnetic field
  # @param[in] milleFile   Millepede-II binary file
  # @param[in] consBL    beam line constraint (for B=0)
  # @return GBL trajectory
  #
  def fitSegment(self, det, qbyp, bField, milleFile, consBL=None):
    # expand track at center
    refZ = self.getZCenter()
    refPos = self.getPosition(refZ)
    refDir = self.getDirection(refZ)
    # dZ/ds
    dZds = refDir[2] 
    # local system
    locsys = globalZXYsystem(qbyp, refDir, refPos)
    # create trajectory
    traj = GblTrajectory(bField[0] != 0. or bField[1] != 0.) 
    jacPointToPoint = np.eye(5)
    sOld = None
    sDUT = None
    #X0Air = None
    X0Air = 304200.  # [mm]
    #
    # loop over planes ordered by Z intersect (for overlapping planes prefer the one with a hit)
    overlap = 1.   
    for planeID in sorted(det.keys(), key=lambda p: det[p].getZintersect() - (overlap if p in self.__hits else 0.)):
      #print " plane ", planeID, planeID in self.__hits
      plane = det[planeID]
      # plane position
      position = plane.getPosition()
      # normal to plane
      normal = plane.getNormal()     
      # measurement direction
      measDir = plane.getMeasDir()
      # arc-length
      sArc = locsys.getArcLengthToPlane(position, normal)
      # skip scattering in air ?
      skipAir = False
      # dummy point in front of first point for beam line constraint
      # (slope with not change until first point: no bending, scattering in air) 
      if consBL is not None and sOld is None:
        skipAir = True
        sOld = sArc - overlap - 1.
        point = GblPoint(jacPointToPoint)
        # for beam line constraint create 4D measurement (u', v', u, v)
        meas4D = np.zeros(4)
        prec4D = np.zeros(4)
        if bField[0] == 0.:
          # no bending in YZ
          meas4D[1] = consBL[0][1] - self.__slope[1]
          prec4D[1] = 1. / consBL[1][1] ** 2
        if bField[1] == 0.:
          # no bending in XZ
          meas4D[0] = consBL[0][0] - self.__slope[0]
          prec4D[0] = 1. / consBL[1][0] ** 2
        #print " beam line constraint ", meas4D, prec4D          
        point.addMeasurement([None, meas4D, prec4D]) 
        # add point to trajectory      
        iLabel = traj.addPoint(point)
      # except for the first point
      if sOld is not None:
        # skip too close by plane
        if sArc < sOld + overlap:
          #print " skip ", planeID, sArc, sOld
          continue
        # scattering in air?
        if (X0Air is not None) and (not skipAir):
          thickness = 0.05
          sAir = sArc - sOld - thickness
          XbyX0 = sAir / X0Air
          sMean = sOld + thickness + 0.5 * sAir 
          sRMS = sAir / math.sqrt(12.)
          # two (identical) equivalent thin scatterers (in local system)
          for s, x in [(sMean - sRMS, XbyX0 * 0.5), (sMean + sRMS, XbyX0 * 0.5)]:
            # update jacobian
            jacPointToPoint = locsys.getSimpleJacobian(s - sOld, bField)
            #print " jac1 ", s
            #print jacPointToPoint
            sOld = s
            # point for scatterer
            point = GblPoint(jacPointToPoint)
            direction = self.getDirection(refZ + dZds * s)
            scatP = plane.getScatPrecision(qbyp, direction, x)   
            scat = np.array([ 0., 0.]) 
            # add scatterer to point
            point.addScatterer([scat, scatP])
            # add point to trajectory      
            traj.addPoint(point) 
        # update jacobian
        jacPointToPoint = locsys.getSimpleJacobian(sArc - sOld, bField)
        #print " jac2 ", s
        #print jacPointToPoint

      sOld = sArc 
      # Z position
      posZ = refZ + dZds * sArc
      # direction (on seed track)
      direction = self.getDirection(posZ)
      #print " sArc ", sArc, posZ
      # position (on seed track)
      position = self.getPosition(posZ)
      # point with (independent) measurements (in measurement system)
      point = GblPoint(jacPointToPoint)
      if planeID in self.__hits:
        # measurement (residuals)
        meas = plane.getResiduals(self.__hits[planeID], position)   
        # measurement resolution
        measPrec = plane.getPrecision()
        # projection from local (=ZXY) to measurement system
        proL2m = locsys.getProjection(measDir)  
        # add measurement to point
        point.addMeasurement([proL2m, meas, measPrec])
        # global derivatives (in local = global ZXY system)
        labGlobal, derGlobal = plane.getGlobalDerivativesLoc(position, direction, '111111')
        point.addGlobals(labGlobal, derGlobal)
        '''
        # to measure scattering in DUT with additional local parameters
        if sDUT is not None:
          addDer = (sArc - sDUT) * proL2m
          point.addLocals(addDer)

      if planeID > 5:
        sDUT = sArc
      '''
      # scattering error
      XbyX0 = plane.getXbyX0()
      if XbyX0 > 0.:
        scatP = plane.getScatPrecision(qbyp, direction, XbyX0)
        scat = np.array([ 0., 0.]) 
        # add scatterer to point
        point.addScatterer([scat, scatP])
      # add point to trajectory      
      iLabel = traj.addPoint(point)
      if planeID in self.__hits:
        self.__hits[planeID].setLabel(iLabel)
   
    # fit trajectory
    Chi2, Ndf = traj.fit()[:2]
    #print " GBL: Chi2, Ndf ", Chi2, Ndf

    if milleFile is not None:
      if Ndf > 0:
        traj.milleOut(milleFile)
    return traj
           
  ## Plot segment.
  def plot(self):
    X = []; Y = []; Z = []
    for h in sorted(self.__hits.itervalues(), key=lambda x: x.getZ()):
      X.append(h.getX())  
      Y.append(h.getY())  
      Z.append(h.getZ())  
    # ZX
    plt.subplot(211)  
    plt.plot(Z, X, 'm--')    
    # ZY 
    plt.subplot(212)     
    plt.plot(Z, Y, 'm--')  
