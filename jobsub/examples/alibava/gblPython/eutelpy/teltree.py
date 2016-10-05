'''
ROOT tree with DUT information

Needs:
 - $ROOTSYS/lib:
 
Created on Jul 23, 2015

@author: kleinwrt
'''

## \file
# telescope ROOT tree for DUT information

from ROOT import TFile, TTree  #@UnresolvedImport
from array import array

## telescope hit.
class TelTree(object):
  
  ## constructor.
  #
  # @param[in] fileName        ROOT file name
  # @param[in] maxClusterSize  maximal cluster size
  #
  def __init__(self, fileName, maxClusterSize=20):
    ## max cluster Size
    self.__maxSize = maxClusterSize
    ## TFile
    self.__tFile = TFile(fileName, 'recreate')
    ## TTree
    self.__tTree = TTree('treeDUT', 'tree with DUT hit and cluster information')
    ## array for run number
    self.__run = array('i', [ 0 ])
    ## array for event number
    self.__evt = array('i', [ 0 ])
    ## array for Pseg 
    self.__Pseg = array('f', [ 0 ])
    ## array for Pgbl
    self.__Pgbl = array('f', [ 0 ])
    ## array for tel plane global pos
    self.__tel_globalPos_0 = array('f', 3 * [ 0 ])
    ## array for tel plane global pos
    self.__tel_globalPos_1 = array('f', 3 * [ 0 ])
    ## array for tel plane global pos
    self.__tel_globalPos_2 = array('f', 3 * [ 0 ])
    ## array for tel plane global pos
    self.__tel_globalPos_3 = array('f', 3 * [ 0 ])
    ## array for tel plane global pos
    self.__tel_globalPos_4 = array('f', 3 * [ 0 ])
    ## array for tel plane global pos
    self.__tel_globalPos_5 = array('f', 3 * [ 0 ])
    ## array for dut plane global pos
    self.__dut_globalPos = array('f', 3 * [ 0 ])
    ## array for hit number
    self.__hit = array('i', [ 0 ])
    ## array for dut id number
    self.__dutID = array('i', [ 0 ])
    ## array for local X position
    self.__locXpos = array('f', [ 0 ])
    ## array for local Y position
    self.__locYpos = array('f', [ 0 ])
    ## array for local X slope
    self.__locXslope = array('f', [ 0 ])
    ## array for local Y slope
    self.__locYslope = array('f', [ 0 ])
    ## array for cluster size
    self.__clusterSize = array('i', [ 0 ])
    ## array for cluster X channel
    self.__clusterXchan = array('i', maxClusterSize * [ 0 ])
    ## array for cluster Y channel
    self.__clusterYchan = array('i', maxClusterSize * [ 0 ])
    ## array for cluster charge
    self.__clusterCharge = array('f', maxClusterSize * [ 0 ])
    # add branches
    self.__tTree.Branch('runnum', self.__run, 'runum/I')
    self.__tTree.Branch('evtnum', self.__evt, 'evtnum/I')
    self.__tTree.Branch('Pseg', self.__Pseg, 'Pseg/F')
    self.__tTree.Branch('Pgbl', self.__Pgbl, 'Pgbl/F')
    self.__tTree.Branch('tel_globalPos_0', self.__tel_globalPos_0, 'tel_globalPos_0[3]/F')
    self.__tTree.Branch('tel_globalPos_1', self.__tel_globalPos_1, 'tel_globalPos_1[3]/F')
    self.__tTree.Branch('tel_globalPos_2', self.__tel_globalPos_2, 'tel_globalPos_2[3]/F')
    self.__tTree.Branch('tel_globalPos_3', self.__tel_globalPos_3, 'tel_globalPos_3[3]/F')
    self.__tTree.Branch('tel_globalPos_4', self.__tel_globalPos_4, 'tel_globalPos_4[3]/F')
    self.__tTree.Branch('tel_globalPos_5', self.__tel_globalPos_5, 'tel_globalPos_5[3]/F')
    self.__tTree.Branch('dut_globalPos', self.__dut_globalPos, 'dut_globalPos[3]/F')
    self.__tTree.Branch('hitnum', self.__hit, 'hitnum/I')
    self.__tTree.Branch('dutID', self.__dutID, 'dutID/I')
    self.__tTree.Branch('locxpos', self.__locXpos, 'locxpos/F')
    self.__tTree.Branch('locypos', self.__locYpos, 'locypos/F')
    self.__tTree.Branch('locxslope', self.__locXslope, 'locxpos/F')
    self.__tTree.Branch('locyslope', self.__locYslope, 'locypos/F')
    self.__tTree.Branch('clustersize', self.__clusterSize, 'clustersize/I')
    self.__tTree.Branch('clusterxchan', self.__clusterXchan, 'clusterxchan[clustersize]/I')
    self.__tTree.Branch('clusterychan', self.__clusterYchan, 'clusterychan[clustersize]/I')
    self.__tTree.Branch('clustercharge', self.__clusterCharge, 'clustercharge[clustersize]/F')
       
  ## fill tree
  #
  # @param[in] run      run number
  # @param[in] evt      event number
  # @param[in] hit      hit number
  # @param[in] locPos   local position
  # @param[in] locSlope local slope
  # @param[in] cluster  cluster information
  #
  def fill(self, run , evt, pseg, pgbl, globalPos_0, globalPos_1, globalPos_2, globalPos_3, globalPos_4, globalPos_5, globalPosDUT, hit, dutID, locPos, locSlope, cluster):
    # fill arrays
    self.__run[0] = run
    self.__evt[0] = evt
    self.__Pseg[0] = pseg
    self.__Pgbl[0] = pgbl
    for i in range(0,3):
      self.__tel_globalPos_0[i] = globalPos_0[i]
      self.__tel_globalPos_1[i] = globalPos_1[i]
      self.__tel_globalPos_2[i] = globalPos_2[i]
      self.__tel_globalPos_3[i] = globalPos_3[i]
      self.__tel_globalPos_4[i] = globalPos_4[i]
      self.__tel_globalPos_5[i] = globalPos_5[i]
      self.__dut_globalPos[i] = globalPosDUT[i] 
    self.__hit[0] = hit
    self.__dutID[0] = dutID
    self.__locXpos[0] = locPos[0]
    self.__locYpos[0] = locPos[1]
    self.__locXslope[0] = locSlope[0]
    self.__locYslope[0] = locSlope[1]
    clusterSize = min(len(cluster) / 4, self.__maxSize)
    self.__clusterSize[0] = clusterSize
    for i in range(clusterSize):
      self.__clusterXchan[i] = int(cluster[4 * i])
      self.__clusterYchan[i] = int(cluster[4 * i + 1])
      self.__clusterCharge[i] = cluster[4 * i + 2]
    # fill tree
    self.__tTree.Fill()
     
  ## write and close file  
  def write(self):
    self.__tFile.Write()
    self.__tFile.Close()
    
