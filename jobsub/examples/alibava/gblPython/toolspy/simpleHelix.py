'''
Created on Aug 22, 2013

@author: kleinwrt
'''
'''
Created on Jun 13, 2013

@author: kleinwrt
'''
## \file
# Simple helix.

import math
import numpy as np

# ######################################################################       
 
## Simple helix
#
# Assuming constant magnetic field in (positive) Z-direction
#
class SimpleHelix(object):  
  
  ## Constructor.
  #
  # @param[in] parameter helix parameter ((curv,) phi0, dca, dzds, z0)
  # @param[in] refPoint   reference point
  #
  # For comparison: Generalized circle equation: 
  #  n_0 + x*n_1 + y*n_2 + (x*x+y*y)*n_3 = 0, 
  #  n_0 ~= -dca, (n_1, n_2) = -(cos(phi_0), sin(phi_0)), n_3 = 0.5*rinv      
  #
  def __init__(self, parameter, refPoint=(0., 0.)):
    ## reference point
    self._refPoint = refPoint   
    ## curvature (in XY)
    self._rinv = parameter[-5] if len(parameter) > 4 else 0.
    ## flight direction at point of closest approach (in XY)
    self._phi0 = parameter[-4]
    ## direction vector at point of closest approach (in XY)
    self._dir0 = (math.cos(self._phi0), math.sin(self._phi0))
    ## distance of closest approach in (XY)
    self._dca = parameter[-3]
    ## dZ/ds
    self._dzds = parameter[-2]
    ## Z position at distance of closest approach
    self._z0 = parameter[-1]
    ## XY circle parameter: X position of center / R
    self._xRelCenter = -(1. - self._dca * self._rinv) * self._dir0[1]
    ## XY circle parameter: Y position of center / R
    self._yRelCenter = (1. - self._dca * self._rinv) * self._dir0[0]
    
  ## Dump helix. 
  def dump(self):
    print " helix ", self._rinv, self._phi0, self._dca, self._dzds, self._z0
    print "   rel. center ", self._xRelCenter, self._yRelCenter
    
  ## Get phi (of point on circle) for given radius (to ref. point)
  #
  # ( |dca| < radius < |rad-2*dca|, from H1/cjfphi, not restricted to -Pi .. +Pi )
  #
  # @param[in] aRadius radius
  # @return azimuth of point on circle
  #
  def getPhi(self, aRadius):
    arg = (0.5 * self._rinv * (aRadius * aRadius + self._dca * self._dca) - self._dca) \
        / (aRadius * (1.0 - self._rinv * self._dca))
    return math.asin(arg) + self._phi0    

  ## Get (2D) arc length for given radius (to ref. point)  
  #
  # ( |dca| < radius < |rad-2*dca|, from H1/cjfsxy )
  #
  # @param[in] aRadius radius
  # @return arc length from dca to point on circle
  #    
  def getArcLengthR(self, aRadius):
    arg = (0.5 * self._rinv * (aRadius * aRadius + self._dca * self._dca) - self._dca) \
        / (aRadius * (1.0 - self._rinv * self._dca))
    if (abs(arg) >= 1.):
      #print " bad arc ", aRadius, self._rinv, self._dca
      return None    
    # line
    if self._rinv == 0:     
      return math.sqrt(aRadius * aRadius - self._dca * self._dca) 
    # helix  
    sxy = math.asin(aRadius * self._rinv * math.sqrt(1.0 - arg * arg)) / self._rinv   
    if (0.5 * self._rinv * self._rinv * \
       (aRadius * aRadius - self._dca * self._dca) - 1. + self._rinv * self._dca > 0.):
      sxy = math.pi / abs(self._rinv) - sxy
    return sxy
   
  ## Get (2D) arc length for given point.
  #
  # Arc length from dca to point on circle on intersection with line
  # from circle center to given point
  #
  # @param[in] pos (XY) Position
  # @return arc length from dca to point on circle
  #  
  def getArcLengthXY(self, pos):
    xLoc = pos[0] - self._refPoint[0]; yLoc = pos[1] - self._refPoint[1]
    # line
    if self._rinv == 0:
      return self._dir0[0] * xLoc + self._dir0[1] * yLoc 
    # helix  
    dx = (xLoc * self._rinv - self._xRelCenter)
    dy = (yLoc * self._rinv - self._yRelCenter)
    dphi = math.atan2(dx, -dy) - self._phi0
    if (abs(dphi) > math.pi):
      dphi -= cmp(dphi, 0.) * 2.0 * math.pi 
    return dphi / self._rinv

  ## Get ZS direction (cosLambda, sinLambda).
  def getZSDirection(self):
    tanl = self._dzds
    cosl = 1. / math.sqrt(1. + tanl * tanl)
    return cosl, cosl * tanl
      
  ## Get expected position in plane. 
  #
  # Plane is defined by direction in XY and Z axis.
  # Intersection of circle and straight line in XY. 
  #
  # @param[in] xPos  X position on plane
  # @param[in] yPos  Y position on plane
  # @param[in] xDir  X direction of plane
  # @param[in] yDir  Y direction of plane
  # @param[in] zPos  Z position on plane
  # @param[in] zDir  Z direction (+/- 1)
  # @param[in] thick thickness of plane (for (Z of) double sided sensors)
  # @return position in plane, Z coordinate (and cos(beta), (2D) arc-length, phi)
  #
  def getExpectedPlanePos(self, xPos, yPos, xDir, yDir, zPos=0., zDir=1., thick=0.):
    xLoc = xPos - self._refPoint[0]; yLoc = yPos - self._refPoint[1]
    rDir = math.sqrt(xDir * xDir + yDir * yDir)
    ex = xDir / rDir
    ey = yDir / rDir
    # line
    if self._rinv == 0:
      cosb = self._dir0[0] * ey - self._dir0[1] * ex
      sinb = self._dir0[0] * ex + self._dir0[1] * ey
      sArc = (ey * xLoc - ex * yLoc - self._dca * sinb) / cosb
      xyPred = (self._dir0[1] * xLoc - self._dir0[0] * yLoc - self._dca) / cosb 
      zPred = (self._z0 + (sArc + thick / abs(cosb)) * self._dzds - zPos) / zDir
      return xyPred, zPred, cosb, sArc, self._phi0 
    # helix  
    dx = (self._xRelCenter - xLoc * self._rinv)
    dy = (self._yRelCenter - yLoc * self._rinv)
    A = ex * dx + ey * dy
    B = dx * dx + dy * dy - 1.
    #C = ex * dy - ey * dx 
    if A * A < B:
      return None
    cosb = A * math.sqrt(1. - B / (A * A))
    #sinb = C
    xyPred = (A - cosb) / self._rinv
    dx = (xLoc + ex * xyPred) * self._rinv - self._xRelCenter
    dy = (yLoc + ey * xyPred) * self._rinv - self._yRelCenter
    phi = math.atan2(dx, -dy)
    dphi = phi - self._phi0
    if (abs(dphi) > math.pi):
      dphi -= cmp(dphi, 0.) * 2.0 * math.pi 
    sArc = dphi / self._rinv 
    zPred = (self._z0 + (sArc + thick / abs(cosb)) * self._dzds - zPos) / zDir
    return xyPred, zPred, cosb, sArc, phi
    
  ## Get analytical helix propagator (in solenoidal magnetic field).
  # 
  # Adapted from TRPRFN.F (GEANT3), 
  # for curvilinear track parameters (q/p,lambda,phi,x_t,y_t).
  #
  # @param[in] phi1 azimutal direction at start point
  # @param[in] ds   (3D) arc length to end point
  # @param[in] bfac Magnetic field strength (*c)
  # @return (5*5) propagation matrix
  #
  def getPropagator(self, phi1, ds, bfac):
    
    ## Cross product (3*3)
    def cross(a, b):
      return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
    ## Dot product (3*3)
    def dot(a, b):
      return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

    curv = self._rinv
    # line
    if curv == 0.:
      ajac = np.eye(5)
      ajac[3][2] = ds / math.sqrt(1. + self._dzds * self._dzds)
      ajac[4][1] = ds    
      return ajac      
  
    # at start
    cosl1 = 1. / math.sqrt(1. + self._dzds * self._dzds)
    sinl1 = cosl1 * self._dzds
    # at end
    cosl2 = cosl1
    sinl2 = sinl1
    cosl2Inv = 1. / cosl2
    # direction of magnetic field 
    hn = [0., 0., 1.]
    # track direction vectors T0, T
    phi2 = phi1 + ds * cosl1 * curv
    t1 = [cosl1 * math.cos(phi1), cosl1 * math.sin(phi1), sinl1]
    t2 = [cosl2 * math.cos(phi2), cosl2 * math.sin(phi2), sinl2]
    #
    q = curv * cosl1
    theta = q * ds
    sint = math.sin(theta)
    cost = math.cos(theta)
    gamma = dot(hn, t2)  # H*T
    an1 = cross(hn, t1)  # HxT0
    an2 = cross(hn, t2)  # HxT
    # U0, V0
    au = 1. / math.sqrt(t1[0] * t1[0] + t1[1] * t1[1])
    u1 = [ -au * t1[1], au * t1[0], 0.]
    v1 = [ -t1[2] * u1[1], t1[2] * u1[0], t1[0] * u1[1] - t1[1] * u1[0] ]
    # U, V
    au = 1. / math.sqrt(t2[0] * t2[0] + t2[1] * t2[1])
    u2 = [ -au * t2[1], au * t2[0], 0.]
    v2 = [ -t2[2] * u2[1], t2[2] * u2[0], t2[0] * u2[1] - t2[1] * u2[0] ]
    #
    qp = -bfac
    pav = qp / cosl1 / curv
    anv = -dot(hn, u2)  # N*V=-H*U
    anu = dot(hn, v2)  # N*U= H*V
    omcost = 1. - cost
    tmsint = theta - sint
    # M0-M
    dx = [ -(gamma * tmsint * hn[0] + sint * t1[0] + omcost * an1[0]) / q,
          - (gamma * tmsint * hn[1] + sint * t1[1] + omcost * an1[1]) / q,
          - (gamma * tmsint * hn[2] + sint * t1[2] + omcost * an1[2]) / q ]
    # HxU0
    hu1 = cross(hn, u1)
    # HxV0
    hv1 = cross(hn, v1)
    # some dot products
    u1u2 = dot(u1, u2); u1v2 = dot(u1, v2); v1u2 = dot(v1, u2); v1v2 = dot(v1, v2) 
    hu1u2 = dot(hu1, u2); hu1v2 = dot(hu1, v2); hv1u2 = dot(hv1, u2); hv1v2 = dot(hv1, v2)
    hnu1 = dot(hn, u1); hnv1 = dot(hn, v1); hnu2 = dot(hn, u2); hnv2 = dot(hn, v2)
    t2u1 = dot(t2, u1); t2v1 = dot(t2, v1);
    t2dx = dot(t2, dx); u2dx = dot(u2, dx); v2dx = dot(v2, dx)
    an2u1 = dot(an2, u1); an2v1 = dot(an2, v1)
    # jacobian
    ajac = np.zeros((5, 5))
    # 1/P
    ajac[0][0] = 1.
    # Lambda
    ajac[1][0] = -qp * anv * t2dx
    ajac[1][1] = cost * v1v2 + sint * hv1v2 + omcost * hnv1 * hnv2 + \
                 anv * (-sint * t2v1 + omcost * an2v1 - gamma * tmsint * hnv1)
    ajac[1][2] = cosl1 * (cost * u1v2 + sint * hu1v2 + omcost * hnu1 * hnv2 + \
                 anv * (-sint * t2u1 + omcost * an2u1 - gamma * tmsint * hnu1))
    ajac[1][3] = -q * anv * t2u1
    ajac[1][4] = -q * anv * t2v1
    # Phi
    ajac[2][0] = -qp * anu * t2dx * cosl2Inv
    ajac[2][1] = cosl2Inv * (cost * v1u2 + sint * hv1u2 + omcost * hnv1 * hnu2 + \
                 anu * (-sint * t2v1 + omcost * an2v1 - gamma * tmsint * hnv1))
    ajac[2][2] = cosl2Inv * cosl1 * (cost * u1u2 + sint * hu1u2 + omcost * hnu1 * hnu2 + \
                 anu * (-sint * t2u1 + omcost * an2u1 - gamma * tmsint * hnu1))
    ajac[2][3] = -q * anu * t2u1 * cosl2Inv              
    ajac[2][4] = -q * anu * t2v1 * cosl2Inv
    # Xt
    ajac[3][0] = pav * u2dx
    ajac[3][1] = (sint * v1u2 + omcost * hv1u2 + tmsint * hnu2 * hnv1) / q             
    ajac[3][2] = (sint * u1u2 + omcost * hu1u2 + tmsint * hnu2 * hnu1) * cosl1 / q 
    ajac[3][3] = u1u2            
    ajac[3][4] = v1u2            
    # Yt
    ajac[4][0] = pav * v2dx
    ajac[4][1] = (sint * v1v2 + omcost * hv1v2 + tmsint * hnv2 * hnv1) / q             
    ajac[4][2] = (sint * u1v2 + omcost * hu1v2 + tmsint * hnv2 * hnu1) * cosl1 / q 
    ajac[4][3] = u1v2            
    ajac[4][4] = v1v2          
 
    return ajac
  
  ## Get simplified helix propagator (in solenoidal magnetic field).
  #
  # Parabola with curvilinear track parameters (q/p,lambda,phi,x_t,y_t).
  #
  # @param[in] ds   (3D) arc length from start to end point
  # @param[in] bfac Magnetic field strength (*c)
  # @return (5*5) propagation matrix
  #
  def getPropagatorSimple(self, ds, bfac):
    cosl = 1. / math.sqrt(1. + self._dzds * self._dzds)
    ajac = np.zeros((5, 5))

    ajac[0][0] = 1.0
    ajac[1][1] = 1.0
    ajac[2][0] = -bfac * ds
    ajac[2][2] = 1.0
    ajac[3][0] = -0.5 * bfac * ds * ds * cosl
    ajac[3][2] = ds * cosl
    ajac[3][3] = 1.0
    ajac[4][1] = ds
    ajac[4][4] = 1.0
    
    return ajac
 
  ## Transform from local to curvilinear system.
  #
  # @param[in] phi  azimutal direction at start point
  # @param[in] i    local system (normal vector)
  # @param[in] j    local system (plane vector 1)
  # @param[in] k    local system (plane vector 2)
  # @return (5*5) propagation matrix
  #
  def localToCurv(self, phi, i, j, k):
    
    ## Cross product (3*3)
    def cross(a, b):
      return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
    ## Dot product (3*3)
    def dot(a, b):
      return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

    curv = self._rinv

    # cos(lambda)
    cosl = 1. / math.sqrt(1. + self._dzds * self._dzds)
    sinl = cosl * self._dzds
    # direction of magnetic field 
    hn = [0., 0., 1.]
    # track direction vector T
    t = [cosl * math.cos(phi), cosl * math.sin(phi), sinl]
    #
    q = curv * cosl
    an = cross(hn, t)  # HxT
    # U, V
    au = 1. / math.sqrt(t[0] * t[0] + t[1] * t[1])
    u = [ -au * t[1], au * t[0], 0.]
    v = [ -t[2] * u[1], t[2] * u[0], t[0] * u[1] - t[1] * u[0] ]
    # some dot products
    ti = dot(t, i); tj = dot(t, j); tk = dot(t, k)
    uj = dot(u, j); uk = dot(u, k)
    vj = dot(v, j); vk = dot(v, k)
    anu = dot(an, u); anv = dot(an, v)
    # jacobian
    ajac = np.zeros((5, 5))
    # 1/P
    ajac[0][0] = 1.
    ajac[1][1] = ti * vj
    ajac[1][2] = ti * vk
    ajac[1][3] = -anv * q * tj
    ajac[1][4] = -anv * q * tk
    ajac[2][1] = ti * uj / cosl
    ajac[2][2] = ti * uk / cosl
    ajac[2][3] = -anu * q * tj / cosl
    ajac[2][4] = -anu * q * tk / cosl
    ajac[3][3] = uj
    ajac[3][4] = uk 
    ajac[4][3] = vj 
    ajac[4][4] = vk
    return ajac

  ## Change reference point
  #
  # @param[in] newRefPoint new reference point       
  # @return    new helix parameters      
  #
  def moveTo(self, newRefPoint):
    dx = self._refPoint[0] - newRefPoint[0]
    dy = self._refPoint[1] - newRefPoint[1]
    dz = 0.
    if len(newRefPoint) > 2: 
      if  len(self._refPoint) > 2:
        dz = self._refPoint[2] - newRefPoint[2]
      else:
        dz = -newRefPoint[2]       
    rho = self._rinv
    phi = self._phi0
    dca = self._dca
    dzds = self._dzds
    z0 = self._z0
    
    cosphi = math.cos(phi)
    sinphi = math.sin(phi)
    u = 1. - rho * dca
    dp = dx * sinphi - dy * cosphi + dca
    dl = dx * cosphi + dy * sinphi
    sa = 2. * dp - rho * (dp * dp + dl * dl)
    sb = -rho * dx + u * sinphi
    sc = rho * dy + u * cosphi
    sd = math.sqrt(1. - rho * sa)
    # transformed parameters
    if rho == 0.:      
      dca = dp
      sArc = -dl
      newPar = [phi, dca]
    else:
      phi = math.atan2(sb, sc)
      dca = sa / (1. + sd)
      dphi = phi - self._phi0 
      if abs(dphi) > math.pi: dphi -= cmp(dphi, 0.) * 2.0 * math.pi 
      sArc = dphi / rho
      newPar = [rho, phi, dca]  
    z0 += sArc * dzds - dz   
    newPar += [dzds, z0]
    
    return newPar
      
  ## Get state (parameters, covariance matrix) at point
  #
  # @param[in] aPoint   point       
  # @param[in] refCov   covariance matrix at reference point
  # @return             state (parameters, covariance matrix)      
  #
  def getStateAt(self, aPoint, refCov):
    dx = self._refPoint[0] - aPoint[0]
    dy = self._refPoint[1] - aPoint[1]
    rho = self._rinv
    phi = self._phi0
    dca = self._dca
    dzds = self._dzds
    z0 = self._z0
    
    cosphi = math.cos(phi)
    sinphi = math.sin(phi)
    u = 1. - rho * dca
    dp = dx * sinphi - dy * cosphi + dca
    dl = dx * cosphi + dy * sinphi
    sa = 2. * dp - rho * (dp * dp + dl * dl)
    sb = -rho * dx + u * sinphi
    sc = rho * dy + u * cosphi
    sd = math.sqrt(1. - rho * sa)
    # transformed parameters
    if rho == 0.:      
      dca = dp
      sArc = -dl
      newPar = [phi, dca]
    else:
      phi = math.atan2(sb, sc)
      dca = sa / (1. + sd)
      dphi = phi - self._phi0 
      if abs(dphi) > math.pi: dphi -= cmp(dphi, 0.) * 2.0 * math.pi 
      sArc = dphi / rho
      newPar = [rho, phi, dca]  
    z0 += sArc * dzds  #+ dz
    newPar += [dzds, z0]
    
    # define jacobian
    if rho == 0.:
      aJac = np.eye(4)
      aJac[1, 0] = -sArc
      aJac[3, 2] = sArc
    else:  
      fc = 1. / (sb * sb + sc * sc)
      fl = 0.5 * sa / (sd * (1. + sd) * (1. + sd))
      fm = -rho * fl + 1. / (sd * (1. + sd))
      aJac = np.zeros((5, 5))
      aJac[0][0] = 1.
      aJac[1][0] = -fc * dl
      aJac[1][1] = fc * u * (1. - rho * dp)
      aJac[1][2] = -fc * rho * rho * dl
      aJac[2][0] = -fm * (dp * dp + dl * dl) + fl * sa
      aJac[2][1] = 2. * fm * u * dl
      aJac[2][2] = 2. * fm * (1. - rho * dp)
      aJac[3][3] = 1.
      aJac[4][0] = dzds * (aJac[1][0] - sArc) / rho
      aJac[4][1] = dzds * (aJac[1][1] - 1.) / rho
      aJac[4][2] = dzds * aJac[1][2] / rho
      aJac[4][3] = sArc
      aJac[4][4] = 1.      
    #print " jac ", aJac
    # transformed covariance
    newCov = np.dot(aJac, np.dot(refCov, aJac.T))    
    
    return np.array(newPar), newCov
            
