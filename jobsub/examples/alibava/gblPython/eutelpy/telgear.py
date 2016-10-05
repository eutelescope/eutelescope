'''
Created on Feb 22, 2015

@author: kleinwrt
'''

## \file
# extract telescope setup from GEAR file 

from xml.dom.minidom import parse
from numpy import array
from telplane import telescopePlane

## GEAR setup.
class GearSetup(object):
  
  ## Parse GEAR file (with DOM).
  def __init__(self, gearFile, alignErr=(0., 0.)):
    ## DOM
    self.__dom = parse(gearFile)
    ## detector
    self.__det = {}
    layers = self.__dom.getElementsByTagName("layer")
    for i, layer in enumerate(layers):
      ladder = layer.getElementsByTagName("ladder")[0]
      lid = int(ladder.getAttribute("ID")) 
      pos = [float(ladder.getAttribute("positionX")), \
           float(ladder.getAttribute("positionY")), \
           float(ladder.getAttribute("positionZ"))]
      print " pos ", i, lid, pos     
      rot = [float(ladder.getAttribute("rotationZY")), \
           float(ladder.getAttribute("rotationZX")), \
           float(ladder.getAttribute("rotationXY"))]
      sens = layer.getElementsByTagName("sensitive")[0]
      xbyx0 = float(sens.getAttribute("thickness")) / float(sens.getAttribute("radLength"))
      res = float(sens.getAttribute("resolution"))
      npix = [int(sens.getAttribute("npixelX")), int(sens.getAttribute("npixelY"))]
      pitch = [float(sens.getAttribute("pitchX")), float(sens.getAttribute("pitchY"))]
      locrot = [ float(sens.getAttribute("rotation1")), float(sens.getAttribute("rotation2")), \
                 float(sens.getAttribute("rotation3")), float(sens.getAttribute("rotation4")) ]
      self.__det[lid] = telescopePlane(i, lid, pos, rot, xbyx0, res, npix, pitch, locrot, alignErr)
      
  ## get magnetic field from DOM.
  def getBField(self):
    bf = self.__dom.getElementsByTagName("BField")[0]
    field = [float(bf.getAttribute("x")), float(bf.getAttribute("y")), float(bf.getAttribute("z"))]
    print  bf.getAttribute("type"), field
    return array(field)  

  ## get detector from DOM.
  #
  # @return dict of \ref telplane::telescopePlane "telescope planes"
  #
  def getDetector(self):
    return self.__det 

    
