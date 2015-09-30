'''
Created on Jun 8, 2015

@author: kleinwrt
'''
## \file
# main program
#
# telescope track finding and fitting

from eutelpy.telgear import GearSetup
from eutelpy.televent import TelEvent
from eutelpy.teltriplets import TelTriplets
from eutelpy.telfile import TelFile
from eutelpy.telalign import TelMP2alignment
from eutelpy.teltree import TelTree
from toolspy.simpleROOTHists import SimpleHists
import time
import sys, getopt
import os

def main(argv):
  start = time.clock()
  runnumber = 0
  gearfile = ''

  try:
    opts, args = getopt.getopt(argv,"hg:r:",["gearfile=","runnumber="])
  except getopt.GetoptError:
    print 'telmain.py -g <gearfile> -r <runnumber>'
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print 'telmain.py -g <gearfile> -r <runnumber>'
      sys.exit()
    elif opt in ("-g", "--gearfile"):
      gearfile = arg
    elif opt in ("-r", "--runnumber"):
      runnumber = arg.zfill(6)
  print 'Gear file is "', gearfile
  print 'Run number is "', runnumber

  # create output directory if it doesn't exist
  outputdir = '/nfs/dust/atlas/user/yeda/ilcsoft/v01-17-05/Eutelescope/trunk/jobsub/examples/alibava/gblPython/output/run' + runnumber
  if not os.path.exists(outputdir):
    os.makedirs(outputdir)

  # maximum number of events
  #maxEvt = 5000
  maxEvt = 10000000
  # beam energy (*q)
  beamEnergy = -4.4  # 613
  # lcio data file
  lciofile = '/nfs/dust/atlas/user/yeda/EUTelescopeOutput/output/lcio/runXXXXXX-hitmaker-local.slcio'
  dataFile = TelFile(lciofile.replace("XXXXXX",runnumber), 'merged_hits')
  dataFile.open() 
  
  # GEAR file
  #gearfile = "gear-start.xml"

  # Millepede-II alignment
  millepedefile = outputdir + '/milleBinary.dat'
  mp2 = TelMP2alignment(millepedefile)
  mp2.open()  # open file for alignment

  # estimated alignment error (mimosa, DUT)
  #alignmentError = (0.1, 0.05)
  #alignmentError = (0.005, 0.01)
  #alignmentError = (0.005, 0.1)
  alignmentError = (0.005, 0.05)
  #alignmentError = (0.05, 1)

  # number of events to show event display
  displayEvents = 0

  # ROOT tree file
  #tree = None
  treefile = outputdir + '/dutTree.root'
  tree = TelTree(treefile)
  
  # parse gear file
  g = GearSetup(gearfile, alignmentError)
  bField = g.getBField()
  detector = g.getDetector()
  # create MP2 steering files
  # for mimosa planes determine positions and XY rotations, fix first plane and position of last as Reference
  #parMimosa = ['RRRRRR', '111001', '111001', '111001', '111001', 'RRR001'] 
  parMimosa = ['RRRRRR', '000000', '000000', '000000', '000000', 'RRRRRR'] 
  if bField[1] != 0.:
    # for B field in Y fix X position of plane 3 in addition
    parMimosa[3] = 'R11001'
  # for DUT determine positions and (all) rotations
  #parDUT = '111111'
  parDUT = '101011'
  #parDUT = '101001'
  # combined alignment of DUTs (in common plane, need to have same intersection point with Z-axis in gear file)
  combDUT = (6, 7)
  # don't combine DUTs
  #combDUT = None    
  steerfile = outputdir + '/steer.txt'
  mp2.createSteering(steerfile, detector, parMimosa, parDUT, combDUT)
  #
  # Cuts in (X,Y) for triplet finder, DUT matching 
  #   (peak from true, background from random matches; tail in bending plane for B<>0, electrons)
  #   Doublet definition: X/Y-distances (mean depends on beam direction and Z-distance)
  #   Triplet definition: X/Y-distances of third hit to doublet (RMS = triplet resolution)
  #   Track definition: X/Y-slope differences of pair of triplets (RMS depends on Z-distances)
  #   Track definition: X/Y-position differences of pair of triplets (RMS depends on Z-distances)
  #   DUT matching: X/Y-distances of DUT hit to track
  #cuts = ((2., 1.), (1., 1.), (0.01, 0.01), (5., 5.))  # iteration 1, no alignment, B on
  #cuts = ((1., 1.), (1., 1.), (0.01, 0.01), (5., 5.))  # iteration 1, no alignment, B off
  #cuts = ((2., 1.), (0.25, 0.25), (0.01, 0.0025), (1., 1.))  # iteration 1, some alignment
  #cuts = ((2., 1.), (0.1, 0.025), (0.01, 0.001), (0.1, 0.1))  # 286, B on
  #cuts = ((1., 1.), (0.025, 0.025), (0.01, 0.01), (0.1, 0.1), (0.2, 10.0))  # 613, B off

  # my cuts to align telescope only
  #cuts = ((1., 1.), (0.5, 0.5), (0.02, 0.02), (5.0, 5.0), (0, 0)) 
  #cuts = ((1., 1.), (0.1 ,0.1), (0.01, 0.01), (0.5, 0.5), (0, 0)) 
  #cuts = ((0.5, 1.), (0.02,0.02), (0.01, 0.01), (0.1, 0.1), (0, 0))  
  #cuts = ((0.5, 0.5), (0.02,0.02), (0.01, 0.01), (0.1, 0.1), (0, 0))  
  #cuts = ((0.5, 0.5), (0.02,0.02), (0.01, 0.01), (0.1, 0.1), (10.2, 10))  
  cuts = ((0.2, 0.2), (0.01,0.01), (0.002, 0.002), (0.1, 0.1), (10.2, 10))  
  #cuts = ((0.5, 0.5), (0.02,0.02), (0.01, 0.01), (0.1, 0.1), (0.2, 10))  


  # histograms for cut values
  #hists = None
  hists = SimpleHists(("doubletDx", "doubletDy", "tripletDx", "tripletDy", "dslopeX", "dslopeY", \
                      "dposX", "dposY", "DUT-dx", "DUT-dy", "nTriplets", "nMatches", "match/triplet", "nSegments"))
  #
  numEvt = 0
  numTrack = 0
  while numEvt < maxEvt:
    event = TelEvent(detector, beamEnergy, bField)
    event.read(dataFile)
    if not event.isValid():
      break
    # event read
    #if event.getNumHits() <= 0:
    if event.getNumDUThits() <= 0:
      continue
    numEvt += 1 
    #event.dump()
    # triplet finder
    finder = TelTriplets(event.getHeader(), detector.keys(), cuts, hists)
    # find tracks and match DUT hits
    event.findTracks(finder, -40.)
    numTrack += event.getNumTracks()
    # fit segments
    event.fitGBL(mp2.getBinaryFile())
    # display event and finder results
    if numEvt < displayEvents:      
      event.plot(False)
      finder.plotTriplets(False) 
      finder.plotSegments() 

    # analyze event
    event.analyze(tree)
    
  dataFile.close() 
  # write and close tree file
  if tree is not None:
    tree.write()
  # plot/save histograms   
  if hists is not None:
    #hists.dump()
    #hists.plot()
    histfile = outputdir + '/cuts.root'
    hists.save(histfile)
  print " Events ", numEvt
  print " Tracks ", numTrack
  end = time.clock()
  print " Time [s] ", end - start  

if __name__ == "__main__":
  main(sys.argv[1:]) 
