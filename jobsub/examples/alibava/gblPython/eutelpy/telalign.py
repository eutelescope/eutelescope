'''
Created on Jun 16, 2015

@author: kleinwrt
'''

## \file
# alignment with Millepede-II


## MP2 alignment
class TelMP2alignment(object):

  ## Constructor
  #
  #  @param[in] fileName   name of MP-II binary file
  #
  def __init__(self, fileName):
    ## file name
    self.__fileName = fileName
    ## binary file
    self.__binaryFile = None
    
  ## open binary file  
  def open(self):  
    self.__binaryFile = open(self.__fileName, "wb")
    
  ## get binary file
  def getBinaryFile(self):
    return self.__binaryFile
    
  ## create steering files
  #
  # @param[in] fileName    name of steering file
  # @param[in] det         detector
  # @param[in] parMimosa   alignment parameters for the Mimosa planes
  # @param[in] parDUT      alignment parameters for DUT planes
  # @param[in] DUTpair     pair of DUT (plane IDs) to combine into one plane
  #
  def createSteering(self, fileName, det, parMimosa, parDUT='111111', DUTpair=None):
    #
    f = open(fileName, 'w') 
    #
    f.write("*** Millepede-II steering file from eutelpy ***\n\n") 
    f.write("Cfiles ! C-type binary files\n%s\n\n" % (self.__fileName))
    f.write('* Global labels are plane-ID * 10 + 1..6 for\n* x-shift, y-shift, z-shift, x-rot, y-rot, z-rot\n\n')
    # check for 1D DUTs
    DUTpar = list(parDUT)
    for planeID, plane in det.iteritems():
      if planeID < 6:
        continue
      # DUT
      prec = plane.getPrecision()
      measDir = plane.getMeasDir()
      for i in range(2):
        if prec[i] == 0.:
          jDir = None
          for j in range(3):
            if abs(measDir[j][i]) == 1.:
              jDir = j
          # missing measurement is in some global direction?    
          if jDir is not None:
            if DUTpar[jDir] == '1':
              f.write("Parameter ! Fix missing measurement direction for plane %i\n" % (planeID))
              f.write(" %6i  0.  -1. ! Fixed\n" % (planeID * 10 + jDir + 1))   
              f.write('\n')
          else:   
            f.write("Constraint  0. ! Fix missing measurement direction for plane %i\n" % (planeID))
            for j in range(3):
              if measDir[j][i] != 0. and DUTpar[j] != '1':
                DUTpar[j] = '1'  #  ! forced on
              f.write(" %6i %f\n" % (planeID * 10 + j + 1, measDir[j][i]))
            f.write('\n')
    # combine 2 DUTs into a single plane
    if DUTpair is not None:
      # position (of planes)
      pos0 = det[DUTpair[0]].getPosition()
      pos1 = det[DUTpair[1]].getPosition()
      f.write("! Force DUT %i and %i to be in the same plane\n" % (DUTpair[0], DUTpair[1]))
      # constraints for Z-shift, x-rot, y-rot
      for iCon, comment in enumerate(["z-shift", "x-rot", "y-rot"]):
        if DUTpar[iCon + 2] == '1':
          f.write("Constraint  0. ! %s\n" % (comment))
          f.write(" %6i %f\n %6i %f\n" % (DUTpair[0] * 10 + iCon + 3, 1., DUTpair[1] * 10 + iCon + 3, -1.))
          if iCon == 0:
            # same dZ in common system
            f.write(" %6i %f\n %6i %f\n" % (DUTpair[0] * 10 + 4, -pos0[1], DUTpair[1] * 10 + 4, pos1[1]))
            f.write(" %6i %f\n %6i %f\n" % (DUTpair[0] * 10 + 5, pos0[0], DUTpair[1] * 10 + 5, -pos1[0]))
          f.write('\n')
    # fix parameters in reference planes and undetermined degrees of freedoms
    f.write("Parameter ! fixed parameters\n") 
    for planeID in det.iterkeys():
      if planeID < 6:
      # Mimosa
        for i in range(6):  
          if parMimosa[planeID][i] in ('r', 'R'):
            f.write(" %6i  0.  -1. ! Reference\n" % (planeID * 10 + i + 1))
          elif parMimosa[planeID][i] != '1':
            f.write(" %6i  0.  -1. ! Unused\n" % (planeID * 10 + i + 1))              
      else:
      # DUT  
        for i in range(6):  
          if DUTpar[i] != '1':
            f.write(" %6i  0.  -1. ! Unused\n" % (planeID * 10 + i + 1))   
    f.write('\n')  
    # pede commands
    f.write("method inversion 5 0.1 ! Solution method\n")
    f.write("outlierdownweighting 4 ! M-estimators to downweight outliers")
    f.write("\nchiscut 30. 6. ! Chi2-cuts for local (track) fits\n")
    f.write("! Other options\nprintcounts\nskipemptycons\nend")
    f.close()
