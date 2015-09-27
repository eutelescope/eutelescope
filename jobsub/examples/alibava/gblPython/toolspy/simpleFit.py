'''
Created on Jun 10, 2013

@author: kleinwrt
'''
## \file
# Simple line and circle fits

import math
import numpy as np

# ######################################################################       

## Simple (XY) line fit: phi, dca
#
# Using residuals perpendicular to line.
#
class SimpleLineFitXY(object): 

  ## Constructor.
  #
  # @param[in] refPoint   reference point
  #
  def __init__(self, refPoint=(0., 0.)):
    ## reference point
    self._refPoint = refPoint   
    ## number of points
    self._numPoints = 0
    ## sum(X)        
    self._sx = 0. 
    ## sum(Y)        
    self._sy = 0. 
    ## sum(X*X)        
    self._sxx = 0. 
    ## sum(X*Y)        
    self._sxy = 0. 
    ## sum(Y*Y)        
    self._syy = 0. 
    ## sum(weight)        
    self._sw = 0.
    ## parameter  phi0, dca
    self._parameter = [0. for i in range(2)]
    ## covariance matrix
    self._covariance = np.zeros((2, 2))
        
  ## Add point.
  #
  # @param[in] position position of point
  # @param[in] weight   weight of point
  #
  def addPoint(self, position, weight):
    self._numPoints += 1
    x = position[0] - self._refPoint[0]
    y = position[1] - self._refPoint[1]
    # sum up
    self._sx += x * weight
    self._sxx += x * x * weight 
    self._sy += y * weight; 
    self._syy += y * y * weight 
    self._sxy += x * y * weight
    self._sw += weight
    
  ## Perform fit.
  def fit(self):
    # averages
    ax = self._sx / self._sw
    axx = self._sxx / self._sw 
    ay = self._sy / self._sw 
    ayy = self._syy / self._sw 
    axy = self._sxy / self._sw  
    # variances  
    cxx = axx - ax * ax 
    cyy = ayy - ay * ay 
    cxy = axy - ax * ay 
   
    q1 = cxy; q2 = cxx - cyy 
    phi = 0.5 * math.atan2(2.*q1, q2)
    sinphi = math.sin(phi)
    cosphi = math.cos(phi)
    # compare phi with initial track direction
    if cosphi * (ax + self._refPoint[0]) + sinphi * (ay + self._refPoint[1]) < 0.:
      # reverse direction
      phi -= cmp(phi, 0) * math.pi
      cosphi = -cosphi; sinphi = -sinphi
    # track parameters
    d = sinphi * ax - cosphi * ay
    self._parameter = [phi, d]
    # chi2
    chi2 = self._sw * \
      (sinphi * sinphi * cxx - 2.*sinphi * cosphi * cxy + cosphi * cosphi * cyy)
 
    # calculate covariance matrix
    aMatrix = [[0. for j in range(2)] for i in range(2)]
    aMatrix[0][0] = cosphi * cosphi * self._sxx + 2. * cosphi * sinphi * self._sxy + sinphi * sinphi * self._syy 
    aMatrix[0][1] = -(cosphi * self._sx + sinphi * self._sy)
    aMatrix[1][0] = aMatrix[0][1]
    aMatrix[1][1] = self._sw
    self._covariance = np.linalg.inv(aMatrix)
    return self._numPoints, chi2
 
  ## Get parameters.
  def getParameters(self):
    return self._parameter  
 
  ## Get covariance matrix.   
  def getCovarianceMatrix(self):
    return self._covariance
 
# ######################################################################       

## Simple (ZS) line fit: dz/ds, z0
#
class SimpleLineFitZS(object): 

  ## Constructor
  def __init__(self):
    
    ## number of points
    self._numPoints = 0
    ## sum(X)        
    self._sx = 0. 
    ## sum(Y)        
    self._sy = 0. 
    ## sum(X*X)        
    self._sxx = 0. 
    ## sum(X*Y)        
    self._sxy = 0. 
    ## sum(Y*Y)        
    self._syy = 0. 
    ## sum(weight)        
    self._sw = 0.
    ## parameter  dz/ds, z0
    self._parameter = [0. for i in range(2)]
    ## covariance matrix
    self._covariance = np.zeros((2, 2))
        
  ## Add poin.t
  #
  # @param[in] position position of point
  # @param[in] weight   weight of point
  #
  def addPoint(self, position, weight):
    self._numPoints += 1
    x = position[0] # s
    y = position[1] # z
    # sum up
    self._sx += x * weight
    self._sxx += x * x * weight 
    self._sy += y * weight; 
    self._syy += y * y * weight 
    self._sxy += x * y * weight
    self._sw += weight
    
  ## Perform fit.
  def fit(self):
    # calculate covariance matrix
    aMatrix = [[self._sxx, self._sx ], [self._sx, self._sw]]
    self._covariance = np.linalg.inv(aMatrix)
    # solution
    dzds = self._sxy * self._covariance[0][0] + self._sy * self._covariance[0][1]
    z0 = self._sxy * self._covariance[1][0] + self._sy * self._covariance[1][1]
    self._parameter = [dzds, z0]
    # chi2
    chi2 = dzds * dzds * self._sxx + z0 * z0 * self._sw + self._syy + \
           2.*dzds * z0 * self._sx - 2.*dzds * self._sxy - 2.*z0 * self._sy
 
    return self._numPoints, chi2
 
  ## Get parameters.   
  def getParameters(self):
    return self._parameter  

  ## Get covariance matrix.       
  def getCovarianceMatrix(self):
    return self._covariance  
 
# ######################################################################       

## Karimaki circle fit
#
class SimpleCircleFit(object): 

  ## Constructor.
  #
  # @param[in] refPoint   reference point
  #
  def __init__(self, refPoint=(0., 0.)):
    ## reference point
    self._refPoint = refPoint       
    ## number of points
    self._numPoints = 0
    ## sum(X)        
    self._sx = 0. 
    ## sum(Y)        
    self._sy = 0. 
    ## sum(X*X)        
    self._sxx = 0. 
    ## sum(X*Y)        
    self._sxy = 0. 
    ## sum(Y*Y)        
    self._syy = 0. 
    ## sum(X*R^2)        
    self._sxr = 0. 
    ## sum(Y*R^2)        
    self._syr = 0. 
    ## sum(R^2)        
    self._sr = 0. 
    ## sum(R^2*R^2)        
    self._srr = 0. 
    ## sum(weight)        
    self._sw = 0.
    ## parameter rho (=-1/R), phi0, dca
    self._parameter = [0. for i in range(3)]
    ## covariance matrix
    self._covariance = np.zeros((3, 3))
        
  ## Add point.
  #
  # @param[in] position position of point
  # @param[in] weight   weight of point
  #
  def addPoint(self, position, weight):
    self._numPoints += 1
    x = position[0] - self._refPoint[0]
    y = position[1] - self._refPoint[1]
    r2 = x * x + y * y
    # sum up
    self._sx += x * weight
    self._sxx += x * x * weight 
    self._sy += y * weight; 
    self._syy += y * y * weight 
    self._sxy += x * y * weight
    self._sxr += x * r2 * weight 
    self._syr += y * r2 * weight 
    self._sr += r2 * weight 
    self._srr += r2 * r2 * weight 
    self._sw += weight
    
  ## Perform fit.
  def fit(self):
    # averages
    ax = self._sx / self._sw
    axx = self._sxx / self._sw 
    ay = self._sy / self._sw 
    ayy = self._syy / self._sw 
    axy = self._sxy / self._sw  
    axr = self._sxr / self._sw 
    ayr = self._syr / self._sw 
    ar = self._sr / self._sw
    arr = self._srr / self._sw  
    # variances  
    cxx = axx - ax * ax 
    cyy = ayy - ay * ay 
    cxy = axy - ax * ay 
    cxr = axr - ax * ar
    cyr = ayr - ay * ar
    crr = arr - ar * ar
    #
    q1 = crr * cxy - cxr * cyr; q2 = crr * (cxx - cyy) - cxr * cxr + cyr * cyr 
    phi = 0.5 * math.atan2(2.*q1, q2)
    sinphi = math.sin(phi)
    cosphi = math.cos(phi)
    # compare phi with initial track direction
    if cosphi * (ax + self._refPoint[0]) + sinphi * (ay + self._refPoint[1]) < 0.:
      # reverse direction
      phi -= cmp(phi, 0) * math.pi
      cosphi = -cosphi; sinphi = -sinphi
    kappa = (sinphi * cxr - cosphi * cyr) / crr
    delta = -kappa * ar + sinphi * ax - cosphi * ay
    # track parameters
    rho = -2.*kappa / math.sqrt(1. - 4.*delta * kappa)
    d = 2.*delta / (1. + math.sqrt(1. - 4.*delta * kappa))
    self._parameter = [rho, phi, d]
    # chi2
    u = 1. - rho * d
    chi2 = self._sw * u * u * \
      (sinphi * sinphi * cxx - 2.*sinphi * cosphi * cxy + cosphi * cosphi * cyy - kappa * kappa * crr)
    #print " xyfit ", chi2, self._numPoints, ":", rho, phi, d
 
    # calculate covariance matrix
    sa = sinphi * self._sx - cosphi * self._sy
    sb = cosphi * self._sx + sinphi * self._sy
    sc = (sinphi * sinphi - cosphi * cosphi) * self._sxy + sinphi * cosphi * (self._sxx - self._syy)
    sd = sinphi * self._sxr - cosphi * self._syr
    saa = sinphi * sinphi * self._sxx - 2.*sinphi * cosphi * self._sxy + cosphi * cosphi * self._syy
    aMatrix = [[0. for j in range(3)] for i in range(3)]
    aMatrix[0][0] = 0.25 * self._srr - d * (sd - d * (saa + 0.5 * self._sr - d * (sa - 0.25 * d * self._sw)))
    aMatrix[0][1] = u * (0.5 * (cosphi * self._sxr + sinphi * self._syr) - d * (sc - 0.5 * d * sb)) 
    aMatrix[1][0] = aMatrix[0][1]
    aMatrix[1][1] = u * u * (cosphi * cosphi * self._sxx + 2. * cosphi * sinphi * self._sxy + sinphi * sinphi * self._syy) 
    aMatrix[0][2] = rho * (-0.5 * sd + d * saa) - 0.5 * u * self._sr + 0.5 * d * ((3 * u - 1.) * sa - u * d * self._sw)
    aMatrix[2][0] = aMatrix[0][2]
    aMatrix[1][2] = -u * (rho * sc + u * sb)
    aMatrix[2][1] = aMatrix[1][2]
    aMatrix[2][2] = rho * (rho * saa + 2 * u * sa) + u * u * self._sw
    self._covariance = np.linalg.inv(aMatrix)
    
    return self._numPoints, chi2

  ## Get parameters.        
  def getParameters(self):
    return self._parameter  
 
  ## Get covariance matrix.          
  def getCovarianceMatrix(self):
    return self._covariance
  
