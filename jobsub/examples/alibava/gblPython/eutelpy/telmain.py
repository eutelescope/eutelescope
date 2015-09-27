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
from toolspy.simpleHists import SimpleHists
import time

if __name__ == '__main__':
  start = time.clock()
  # maximum number of events
  maxEvt = 1000
  # beam energy (*q)
  #beamEnergy = -3.0 # 286
  beamEnergy = -4.4  # 613
  # lcio data file
  #dataFile = TelFile('run000613-hitmaker.slcio', 'merged_hits')
  # text data file
  dataFile = TelFile('allhits-613.dat')
  dataFile.open() 
  # GEAR file
  gearFile = "gear-613.xml"
  # Millepede-II alignment
  mp2 = TelMP2alignment("milleBinary.dat")
  mp2.open()  # open file for alignment
  # estimated alignment error (mimosa, DUT)
  alignmentError = (0.005, 0.01)
  # number of events to show event display
  displayEvents = 0
  
  # parse gear file
  g = GearSetup(gearFile, alignmentError)
  bField = g.getBField()
  detector = g.getDetector()
  # create MP2 steering files
  mp2.createSteering("steer.txt", detector, bField, '111001', '111111', (6, 7))
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
  cuts = ((1., 1.), (0.025, 0.025), (0.01, 0.01), (0.1, 0.1), (0.2, 10.0))  # 613, B off
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
    if event.getNumHits() <= 0:
    #if event.getNumDUThits() <= 0:
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
    event.analyze()
    
  dataFile.close()  
  if hists is not None:
    #hists.dump()
    #hists.plot()
    hists.save("cuts.pdf")
  print " Events ", numEvt
  print " Tracks ", numTrack
  end = time.clock()
  print " Time [s] ", end - start   
