'''
Created on Apr 23, 2015

@author: kleinwrt
'''

## \file
# telescope triplets

from telsegment import TelSegment

## TPC Doublet
#
class TelDoublet(object):

  ## Constructor.
  #
  # @param[in] hit1    first hit
  # @param[in] hit2    second hit
  # @param[in] curv    assumed curvature
  # @param[in] zMag   start (in z) of magnetic field
  # 
  def __init__(self, hit1, hit2, curv, zMag):
    # position and direction (corrected for curvature)
    x1, y1, z1 = hit1.getPos()
    x1 -= 0.5 * curv[0] * (z1 - zMag) ** 2
    y1 -= 0.5 * curv[1] * (z1 - zMag) ** 2
    x2, y2, z2 = hit2.getPos()
    x2 -= 0.5 * curv[0] * (z2 - zMag) ** 2
    y2 -= 0.5 * curv[1] * (z2 - zMag) ** 2
    ## curvature (in XZ, YZ)
    self.__curv = curv
    ## start of magnetic field
    self.__zMag = zMag
    # mean position
    xm = (x1 + x2) / 2.; ym = (y1 + y2) / 2.; zm = (z1 + z2) / 2.
    ## position
    self.__position = (xm, ym, zm)
    #
    dx = x2 - x1; dy = y2 - y1; dz = z2 - z1
    ## slope
    self.__slope = (dx / dz, dy / dz)
    ## dx
    self.__dx = dx
    ## dy
    self.__dy = dy

  ## get distances     
  def getDists(self):
    return self.__dx, self.__dy
    
  ## get position
  def getPos(self):
    return self.__position

  ## get slope
  def getSlope(self):
    return self.__slope
  
  # match with hit 
  #
  # @param[in]  hit  hit
  # @param[in]  distCut  distance cuts (in x, y)
  #   
  def match(self, hit, distCut):
    # triplet residuals
    x, y, z = hit.getPos()
    x -= 0.5 * self.__curv[0] * (z - self.__zMag) ** 2
    y -= 0.5 * self.__curv[1] * (z - self.__zMag) ** 2    
    dz = z - self.__position[2]
    dx = x - self.__position[0] - dz * self.__slope[0]
    dy = y - self.__position[1] - dz * self.__slope[1]
    #
    if abs(dx) > distCut[0] or abs(dy) > distCut[1]:
      return None
    #
    return (dx, dy)

## telescope triplet
#
class TelTriplet(object):

  ## Constructor.
  #
  # @param[in] hits    hits
  # @param[in] dbl     seeding doublet
  # @param[in] trp     pair of 2D triplets in XY, ZS 
  # 
  def __init__(self, hits, dbl, trp):
    ## hits
    self.__hits = hits
    ## pair of 2D triplets in XY, ZS
    self.__triplet = trp 
    ## position
    self.__position = dbl.getPos()
    ## slope
    self.__slope = dbl.getSlope()
    ## matches
    self.__matches = 0
  
  ## Add match
  def addMatch(self):
    self.__matches += 1
    
  ## Get matches
  def getMatches(self):
    return self.__matches
       
  ## Get hits.
  def getHits(self):
    return self.__hits

  ## Get triplet.
  def getTriplet(self):
    return self.__triplet

  ## Get position.
  def getPos(self):
    return self.__position

  ## Get position at different z value.
  #
  # @param[in] z   z position
  #
  def getPosAt(self, z):
    dz = z - self.__position[2] 
    x = self.__position[0] + self.__slope[0] * dz
    y = self.__position[1] + self.__slope[1] * dz
    return (x, y, z)

  ## get slope
  def getSlope(self):
    return self.__slope
    
  ## dump
  def dump(self):
    print " triplet ", self.__hits[1].getLayer(), self.__matches,
    for p in self.__position:
      print p,
    for s in self.__slope:
      print s,
    for t in self.__triplet:
      print t,
    print  
             
## Triplets for track finding.
class TelTriplets(object):    
  
  ## Constructor.
  #
  # @param[in] header   event header
  # @param[in] planeIDs IDs of detector planes
  # @param[in] cuts     cuts
  # @param[in] hists    histograms
  #
  def __init__(self, header, planeIDs, cuts, hists):
    ## (event) header
    self.__hdr = header 
    ## hits (in layers)
    self.__hits = {}
    for p in planeIDs:
      self.__hits[p] = []
    ## triplets
    self.__triplets = [[], []]   
    ## segments
    self.__seg = []
    ## cuts for triplets
    self.__cuts = cuts
    ## histograms
    self.__hists = hists
  
  ## Get segments.
  def getSegments(self):
    return self.__seg
        
  ## Plot segments.
  #
  # @param[in] show  flag for showing plot
  #
  def plotSegments(self, show=True):
    plt.figure(1)
    for segment in self.__seg:
      segment.plot()
             
    if show: 
      plt.show()  
      
  ## Find tracks.    
  #
  #  Find triplets in first and last three mimosa planes and combine unique matches 
  #  of those to tracks. Requires a single hit in all six mimosa planes on each track. 
  #
  # @param[in] hits   list of hits
  # @param[in] qbyp   assumed Q/P (e.g. from beam energy)
  # @param[in] bfac   deflection strength (c*BxZ)
  # @param[in] zMag   start (in z) of magnetic field
  #
  def findTracks(self, hits, qbyp=0., bfac=(0., 0.), zMag=0.):   
    # cuts
    doubletCut, tripletCut, slopeCut, positionCut = self.__cuts[:4]
    # curvature
    curv = [bfac[0] * qbyp, bfac[1] * qbyp]
    #
    numFound = 0
    #print " hits ", len(hits)
    for h in hits:
      self.__hits[h.getPlane()].append(h)
    #for i in range(6):
    #  print " layer ", i, len(self.__hits[i])
    # look for triplets in layers (0,1,2) and (3,4,5)
    nTrp = 0
    for lcenter in [1, 4]:
      for hit1 in self.__hits[lcenter - 1]:
        for hit3 in self.__hits[lcenter + 1]:
          aDoublet = TelDoublet(hit1, hit3, curv, zMag)
          dx, dy = aDoublet.getDists()
          if abs(dx) > doubletCut[0] or abs(dy) > doubletCut[1]:
            continue
          if self.__hists is not None:
            self.__hists.addEntry("doubletDx", dx)
            self.__hists.addEntry("doubletDy", dy)
          for hit2 in self.__hits[lcenter]:
            trp = aDoublet.match(hit2, tripletCut)
            if trp is not None:
              nTrp += 1
              #print " trp ", lcenter, dx, dy, trp[0], trp[1]
              self.__triplets[lcenter / 3].append(TelTriplet([hit1, hit2, hit3], aDoublet, trp))
              if self.__hists is not None:
                self.__hists.addEntry("tripletDx", trp[0])
                self.__hists.addEntry("tripletDy", trp[1])

    #print " nTrp ", nTrp
    if self.__hists is not None:
      self.__hists.addEntry("nTriplets", nTrp)        
    # combine triplets
    matches = []
    for i1, trp1 in enumerate(self.__triplets[0]):
      slope1 = trp1.getSlope()
      z1 = trp1.getPos()[2]
      for i2, trp2 in enumerate(self.__triplets[1]):
        slope2 = trp2.getSlope()
        dslopeX = slope1[0] - slope2[0]; dslopeY = slope1[1] - slope2[1] 
        #print " slopes ", slope1, slope2, dslope
        if abs(dslopeX) > slopeCut[0] or abs(dslopeY) > slopeCut[1]:
          continue
        if self.__hists is not None:
          self.__hists.addEntry("dslopeX", dslopeX)
          self.__hists.addEntry("dslopeY", dslopeY)
        z2 = trp2.getPos()[2]
        zmean = 0.5 * (z1 + z2)
        #print " zmean ", zmean
        pos1 = trp1.getPosAt(zmean)
        pos2 = trp2.getPosAt(zmean)
        #print " match ", dslopeX, dslopeY, pos1[0], pos1[1], pos2[0], pos2[1]
        dposX = pos1[0] - pos2[0]; dposY = pos1[1] - pos2[1]
        if abs(dposX) < positionCut[0] and abs(dposY) < positionCut[1]:
          matches.append((i1, i2))
          trp1.addMatch()
          trp2.addMatch()
          if self.__hists is not None:
            self.__hists.addEntry("dposX", dposX)
            self.__hists.addEntry("dposY", dposY)
          
    #print " matches ", matches
    if self.__hists is not None:
      self.__hists.addEntry("nMatches", len(matches))
      for triplets in self.__triplets:
        for trp in triplets:
          self.__hists.addEntry("match/triplet", trp.getMatches()) 
    # look for unique matches
    for i1, i2 in matches:
      trp1 = self.__triplets[0][i1] 
      trp2 = self.__triplets[1][i2] 
      #print " i ", i1, i2, trp1.getMatches(), trp2.getMatches() 
      if  trp1.getMatches() == 1 and trp2.getMatches() == 1:
        hits = trp1.getHits() + trp2.getHits()
        # improve curvature from slope difference
        if bfac[0] != 0. or bfac[1] != 0.:
          slope1 = trp1.getSlope(); z1 = trp1.getPos()[2]  
          slope2 = trp2.getSlope(); z2 = trp2.getPos()[2]
          # derivatives (d(slope1-slope2)/d(qbyp))
          derx = bfac[0] * (z1 - z2); dery = bfac[1] * (z1 - z2)
          # qbyp correction
          dqbyp = (derx * (slope1[0] - slope2[0]) + dery * (slope1[1] - slope2[1])) \
                / (derx ** 2 + dery ** 2)
          #print " dqbyp(slope) ", dqbyp, dqbyp * bfac[0], dqbyp * bfac[1] 
          curv[0] += dqbyp * bfac[0]; curv[1] += dqbyp * bfac[1]
        self.__seg.append(TelSegment(self.__hdr, hits, curv)) 
        numFound += 1
        
    #print " nSeg ", len(self.__seg)
    if self.__hists is not None:
      self.__hists.addEntry("nSegments", len(self.__seg))      
        
  ## match DUT
  #
  # Unique matching of DUT hits (plane ID > 5) and segments.
  #
  def matchDUT(self):
    # DUT matching cuts defined ?
    if len(self.__cuts) < 5:
      return
    
    distCut = self.__cuts[4]
    for l, hits in self.__hits.iteritems():
      if l < 6:
        continue
      # look for unique matches
      matches = []
      for seg in self.__seg:
        seg.resetMatches()
      for i, h in enumerate(hits):  
        x, y, z = h.getPos()
        for j, seg in enumerate(self.__seg):
          px, py = seg.getPrediction(z)
          dx = abs(x - px)
          dy = abs(y - py)
          if dx < distCut[0] and dy < distCut[1]:
            matches.append((i, j))
            h.addMatch()
            seg.addMatch()
            if self.__hists is not None:
              if h.getPlane() == 6:
                self.__hists.addEntry("DUT6-dx", x - px)
                self.__hists.addEntry("DUT6-dy", y - py)
              elif h.getPlane() == 7:
                self.__hists.addEntry("DUT7-dx", x - px)
                self.__hists.addEntry("DUT7-dy", y - py)
              else:
                print "DUT not identified"
      for i, j in matches:
        h = hits[i]
        seg = self.__seg[j]
        if h.getMatches() == 1 and seg.getMatches() == 1:
          seg.addDUT(h)
        
  ## plot triplets
  #
  # @param[in] show  flag for showing plot
  #
  def plotTriplets(self, show=True):
    plt.figure(1)
    for superLayer in self.__triplets:
      for triplet in superLayer:
        X = []
        Y = []
        Z = []
        for h in triplet.getHits():
          X.append(h.getX())      
          Y.append(h.getY())      
          Z.append(h.getZ())      
        plt.subplot(211)     
        plt.plot(Z, X, 'g-')
        plt.subplot(212)     
        plt.plot(Z, Y, 'g-')
    if show: 
      plt.show()  
      
  ## dump triplets
  def dumpTriplets(self):
    for superLayer in self.__triplets:
      for triplet in superLayer: 
        triplet.dump()
                     
  ## write hits on segments to text file
  #
  # @param[in]  fout  output file
  # @param[in]  be    beam energy
  #
  def writeHitsOnSegments(self, fout, be):
    event = self.__hdr.getRunEvent()[1]
    for seg in self.__seg:
      fout.write("%i %.2f \n" % (event, be))
      for h in seg.getHits():
        l = h.getLayer()
        x, y, z = h.getPos()
        fout.write("%i %.4f %.4f %.4f\n" % (l, x, y, z))
        
