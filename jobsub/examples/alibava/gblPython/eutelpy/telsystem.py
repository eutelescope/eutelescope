'''
Created on Feb 22, 2015

@author: kleinwrt
'''

## \file
# define local system(s)

import numpy as np
import math

## global ZXY system
#
# Define global (Z,X,Y) as local (I,J,K) system.
# Conventions according to A. Strandlie, W. Wittek, NIM A, 566 (2006) 687-698.
#
class globalZXYsystem(object):

  ## Constructor
  #
  #  Construct local system with track parameters
  #
  #  @param qbyp       q/p
  #  @param direction  direction
  #  @param position   position
  #
  def __init__(self, qbyp, direction, position):
    ## track position
    self.__tpos = np.array(position)
    ## track direction
    self.__tdir = np.array(direction)
    ## track parameters (q/p, dx/dz, dy/dz, x, y)
    self.__tpar = np.array([qbyp, direction[0] / direction[2], direction[1] / direction[2], \
                            position[0], position[1]])
    ## directions for x, y  (2*3)
    self.__xydir = np.array([[1., 0., -direction[0] / direction[2]], \
                             [0., 1., -direction[1] / direction[2]]])
        
  ## Get arc-length to plane
  #
  # Intersect straight line with plane.
  #
  # @param position  point on plane (any)
  # @param normal    normal vector of plane
  # @return arc-length (3D)
  #
  def getArcLengthToPlane(self, position, normal):
    return np.dot(position - self.__tpos, normal) / np.dot(self.__tdir, normal)
     
  ## Get simplified helix propagator (arbitrary magnetic field).
  #  
  #  Simple jacobian for track parameters (q/p, slopes, offsets) 
  #  in (Z,X,Y) system (offsets=(x,y), slopes=(dx/dz,dy/dz)),
  #  quadratic in arc length difference (first order in bending Q).
  #  Units are: T, GeV, mm 
  #
  #  Derived from full curvilinear propagator and transformations from/to
  #  local (I,J,K) system in the limit Q -> 0.
  #
  #  @param[in] ds     arc length difference; float
  #  @param[in] bfield magnetic field (vector) 
  #  @return (5*5) propagation matrix
  #
  def getSimpleJacobian(self, ds, bfield):
  
    sinLambda = self.__tdir[2]
    BxT = np.cross(bfield, self.__tdir)  # BxT
    bfac = -0.0002998 * np.dot(self.__xydir, BxT)
    jac = np.eye(5)
    jac[1, 0] = bfac[0] * ds / sinLambda
    jac[2, 0] = bfac[1] * ds / sinLambda
    jac[3, 0] = 0.5 * bfac[0] * ds * ds
    jac[4, 0] = 0.5 * bfac[1] * ds * ds
    jac[3, 1] = ds * sinLambda
    jac[4, 2] = ds * sinLambda 
    return jac        

  ## Get scattering precision
  #
  #  Multiple scattering precision matrix in local system for a thin scatterer.
  #
  #  @param  XbyX0   X/X0 (thickness / radiation lenght)
  #  @param  normal  normal vector to scattering plane
  #  @return (2*2) precision matrix
  #   
  def getScatPresicion(self, XbyX0, normal):
    # scattering error
    cosIncident = np.dot(self.__tdir, normal)  # T*I
    #scatErr = 0.015 * self.__tpar[0] * math.sqrt(XbyX0 / abs(cosIncident))
    scatErr = 0.0136 * self.__tpar[0] * math.sqrt(XbyX0 / abs(cosIncident)) * 0.80  # with avarage log corr
    # Scattering precision matrix in local system  
    c1 = self.__tdir[0]  # T*J
    c2 = self.__tdir[1]  # T*K
    fac = (1 - c1 * c1 - c2 * c2) / (scatErr * scatErr)
    scatP = np.empty((2, 2))
    scatP[0, 0] = fac * (1 - c1 * c1)
    scatP[0, 1] = fac * (-c1 * c2)
    scatP[1, 0] = fac * (-c1 * c2)
    scatP[1, 1] = fac * (1 - c2 * c2)
    return scatP

  ## get projection local to measurement system
  #
  #  Derived from projections from measurement and local to curvilinear system.
  #
  #  @param measDir  measurement directions (3*2 matrix)
  #  @return (2*2) projection matrix 
  #
  def getProjection(self, measDir):
    # projection measurement to local system
    proM2l = np.dot(self.__xydir, measDir)
    # projection local to measurement system
    proL2m = np.linalg.inv(proM2l)
    return proL2m
  
  ## get transformation measurement to local system
  #
  #  Derived from transformations from measurement and local to curvilinear system.
  #
  #  @param measDir  measurement directions (3*2 matrix)
  #  @param  normal  normal vector to measurement plane
  #  @return (5*5) transformation matrix for track parameters (q/p, slopes, offsets)
  #
  def getTransMeasToLocal(self, measDir, normal):
    # T*I
    cosIncident = np.dot(self.__tdir, normal)
    # projection measurement to local system
    proM2l = np.dot(self.__xydir, measDir)
    # transformation measurement to local system
    transM2l = np.eye(5)
    transM2l[1:3, 1:3] = proM2l * cosIncident / self.__tdir[2] 
    transM2l[3:5, 3:5] = proM2l 
    return transM2l
  
  ## get transformation local to measurement system
  #
  #  Derived from transformations from measurement and local to curvilinear system.
  #
  #  @param measDir  measurement directions (3*2 matrix)
  #  @param  normal  normal vector to measurement plane
  #  @return (5*5) transformation matrix for track parameters (q/p, slopes, offsets)
  #
  def getTransLocalToMeas(self, measDir, normal):
    # transformation local to measurement system
    return np.linalg.inv(self.getTransMeasToLocal(measDir, normal))
