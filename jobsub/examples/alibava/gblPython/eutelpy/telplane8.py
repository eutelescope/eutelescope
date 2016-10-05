'''
Created on Feb 22, 2015

@author: kleinwrt
'''

## \file
# Telescope plane

import numpy as np
import math

## Telescope plane
class telescopePlane(object):

  ## Create new plane
  #
  # @param  idx    index
  # @param  pid    ID
  # @param  pos    position
  # @param  rot    rotations
  # @param  xbyx0  thickness / radiation length
  # @param  res    resolution
  # @param  npix   number of pixels (x,y)
  # @param  pitch  pitch 
  # @param  locrot local rotation
  # @param  aerr   alignment error (mimosa, DUT)
  #
  def __init__(self, idx, pid, pos, rot, xbyx0, res, npix, pitch, locrot, aerr):
    ## index
    self.__index = idx
    ## ID
    self.__id = pid
    ## position
    self.__pos = np.array(pos)
    # layer orientation (transformation from local to global system)
    alpha = rot[0] / 57.3; cos = math.cos(alpha); sin = math.sin(alpha) 
    rotx = np.array([[1., 0., 0.], [0., cos, -sin], [0., sin, cos]])  
    alpha = rot[1] / 57.3; cos = math.cos(alpha); sin = math.sin(alpha) 
    roty = np.array([[cos, 0., sin], [0., 1., 0.], [-sin, 0., cos]])  
    alpha = rot[2] / 57.3; cos = math.cos(alpha); sin = math.sin(alpha) 
    rotz = np.array([[cos, -sin, 0.], [sin, cos, 0.], [0., 0., 1.]]) 
    rotl = np.array([[locrot[0], locrot[1], 0.], [locrot[2], locrot[3], 0.], [0., 0., 1.]]) 
    ## rotation (inverse = transposed)
    self.__rot = np.dot(roty, np.dot(rotx, np.dot(rotz, rotl)))  # rot=roty*rotx*rotz*rotl 
    ## xbyx0
    self.__xbyx0 = xbyx0
    ## resolution
    ##self.__res = res
    if pid < 6:    
      #self.__res = [0.006, 0.006]  #  6 microns
      self.__res = [res, res]  #  from gear
      alignErr = aerr[0]
    else:
      self.__res = [pitch[0] / math.sqrt(12.), pitch[1] / math.sqrt(12.)]
      alignErr = aerr[1]
    # alignment error  
    if alignErr > 0:
      self.__res[0] = math.sqrt(self.__res[0] ** 2 + alignErr ** 2)   
      self.__res[1] = math.sqrt(self.__res[1] ** 2 + alignErr ** 2)   
    ## number of pixels
    self.__npix = npix

  ## Transform from local to global system.
  #
  # @param[in] locPos local position
  #
  def transformLocalToGlobal(self, locPos):
    return np.dot(self.__rot, locPos) + self.__pos
  
  ## Transform from global to local system.
  #
  # @param[in] gloPos    global position
  #
  def transformGlobalToLocal(self, gloPos):
    return np.dot(self.__rot.T, gloPos - self.__pos)
  
  ## get ID
  def getID(self):
    return self.__id
  
  ## get position
  def getPosition(self):
    return self.__pos      
    
  ## get intersection with Z-axis
  def getZintersect(self):
    return np.dot(self.__pos, self.__rot[:, 2])/self.__rot[2, 2]    
       
  ## get rotation
  def getRotation(self):
    return self.__rot
    
  ## get normal vector
  def getNormal(self):
    return self.__rot[:, 2]
   
  ## get measurement directions
  def getMeasDir(self):
    return self.__rot[:, :2]
  
  ## get (diagonal of) precision matrix
  #
  # Measurement directions with only 1 pixel are ignored (strip detector)
  #
  def getPrecision(self):
    prec = np.zeros(2)
    for i in range(2):
      if self.__npix[i] > 1:  # strips have only one sensitive direction (other: npix=1)
        prec[i] = 1. / self.__res[i] ** 2
    return prec
 
  ## get X/X0 
  def getXbyX0(self):
    return self.__xbyx0
    
  ## Get scattering precision
  #
  #  Multiple scattering precision matrix in local system for a thin scatterer.
  #
  #  @param  qbyp       q/p
  #  @param  direction  direction (of track) at hit
  #  @param  XbyX0      X/X0 (thickness / radiation lenght)
  #  @return (2*2) precision matrix
  #   
  def getScatPrecision(self, qbyp, direction, XbyX0):
    # scattering error
    cosIncident = np.dot(direction, self.__rot[:, 2])  # T*I
    #scatErr = 0.015 / energy * math.sqrt(XbyX0 / abs(cosIncident))
    scatErr = 0.0136 * qbyp * math.sqrt(XbyX0 / abs(cosIncident)) * 0.80  # with average log correction
    # Scattering precision matrix in local system  
    c1 = np.dot(direction, self.__rot[:, 0])  # T*J
    c2 = np.dot(direction, self.__rot[:, 1])  # T*K
    #print " c_i ", c1, c2, scatErr
    fac = (1 - c1 * c1 - c2 * c2) / (scatErr * scatErr)
    scatP = np.empty((2, 2))
    scatP[0, 0] = fac * (1 - c1 * c1)
    scatP[0, 1] = fac * (-c1 * c2)
    scatP[1, 0] = fac * (-c1 * c2)
    scatP[1, 1] = fac * (1 - c2 * c2)
    return scatP  
  
  ## get global (rigid body) derivatives (in measurement system)
  #
  # Global derivatives in the plane/measurement system for a rigid body.
  #
  # @param   hitPos     position of hit (measured or predicted)
  # @param   direction  direction (of track) at hit
  # @param   mask       mask for active parameters
  # @return  global labels and derivatives
  #
  def getGlobalDerivativesMeas(self, hitPos, direction, mask='111111'):
    # transform to (measurement) plane system
    posPlane = np.dot(self.__rot.T, hitPos - self.__pos)
    dirPlane = np.dot(self.__rot.T, direction)
    # dr/dm
    normal = np.array([0., 0., 1.])
    drdm = np.eye(3) - np.outer(dirPlane, normal) / np.dot(dirPlane, normal)
    # dm/dg
    posPlane[2] = 0.  # should be zero in local system
    dmdg = np.array([[1., 0., 0., 0., posPlane[2], -posPlane[1]], \
                     [0., 1., 0., -posPlane[2], 0., posPlane[0]], \
                     [0., 0., 1., posPlane[1], -posPlane[0], 0.]])             
    # dr/dg
    drdg = np.dot(drdm, dmdg)
    # select active parameters
    xder = []
    yder = []
    labels = []   
    for i in range(len(mask)):
      if mask[i] == '1' and self.__index >= 0:
        labels.append(self.__id * 10 + i + 1)
        xder.append(drdg[0, i])
        yder.append(drdg[1, i]) 
       
    return np.array([labels, labels], dtype=int), np.array([xder, yder])
    
  ## get global (rigid body) derivatives (in local system)
  #
  # Global derivatives in the plane/measurement system for a rigid body.
  #
  # @param   hitPos     position of hit (measured or predicted)
  # @param   direction  direction (of track) at hit
  # @param   mask       mask for active parameters
  # @return  global labels and derivatives
  #
  def getGlobalDerivativesLoc(self, hitPos, direction, mask='111111', planeID=-1):
    # rotate around (center) position
    if self.__id > 5:
      # DUT (rotate around intersection with Z axis)
      posRel = hitPos - np.array([0., 0., self.getZintersect()])
    else:
      # mimosa
      posRel = hitPos - self.__pos
    # drl/dm
    normal = self.__rot[:, 2]
    drldm = np.eye(3) - np.outer(direction, normal) / np.dot(direction, normal)
    # dm/dg
    dmdg = np.array([[1., 0., 0., 0., posRel[2], -posRel[1]], \
                     [0., 1., 0., -posRel[2], 0., posRel[0]], \
                     [0., 0., 1., posRel[1], -posRel[0], 0.]])
    # dr/dg
    drdg = np.dot(self.__rot.T, np.dot(drldm, dmdg))
    # select active parameters
    xder = []
    yder = []
    labels = [] 
    ioff = 0
    '''
    mask9=mask
    if planeID > 5:
      #ioff = 100*max(0,min(19,int(hitPos[0]+10.)))
      ioff = 100*max(0,min(9,int(hitPos[1]+5.)))
      mask9 = '101000'
    '''  
    for i in range(len(mask)):
      if mask[i] == '1' and self.__index >= 0:
        labels.append(self.__id * 10 + i + 1 + ioff)
        xder.append(drdg[0, i])
        yder.append(drdg[1, i]) 
       
    return np.array([labels, labels], dtype=int), np.array([xder, yder])
      
  ## get residuals
  #
  #  @param  hit         hit
  #  @param  prediction  prediction of hit position
  #  @return residuals vector (or zeros if no hit found in plane) 
  #
  def getResiduals(self, hit, prediction):
    return np.dot(self.__rot.T, np.array(hit.getPos()) - prediction)[:2] 
