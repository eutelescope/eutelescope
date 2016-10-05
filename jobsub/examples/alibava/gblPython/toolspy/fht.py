'''
Created on Dec 1, 2013

@author: kleinwrt
'''

# # \file
# Fast Hough transformation

import math
import matplotlib.pyplot as plt

## cube counter
nCube = 0

## Fast Hough transformation
#
# Fast Hough transformation to construct tracks from space points.
#
# The average track direction is in (local) X direction. 
# The measured coordinates have to be shifted and rescaled to: X->u in [-1, +1],
# Y,Z -> v,w in [-0.5, 0.5]. The resolution in v and w should be comparable.
# The track projections v(u), w(u) 
# are parametrized by series of Legendre polynomials L_i, in the (XY) bending plane up to 
# order 2 and in the (XZ) plane with order 1.:
#   - v(u) = a_0 * L_0(u) + a_1 * L_1(u) (+ a_2 *L_2(u))
#   - w(u) = b_0 * L_0(u) + b_1 * L_1(u)
#   - with L_0(u) = 1., L_1(u) = u, L_2(u) = 1.5*u*u - 0.5
#
# The track parameters a_i, b_i are in [-0.5, +0.5] (the initial hyper cube). 
# Each space point defines 2 hyper planes (v(u)-v_meas=0, w(u)-w_meas=0) in the 5 (or 4)
# dimensional parameter space (a_i, b_i) with normal vectors (L_i(u), 0.) and (0., L_i(u)).
# The unit normal vector is used to calculate the distance 
# of a hyper plane to the center of a hyper cube. They intersect if the 
# absolute value of the distance times the largest component of the unit normal vector is
# smaller than 0.5. The root cube is recursivly split for each dimension in two child cubes.
# Only childs containing a minimum number of rows are considered further until the hit
# density in a cube is compatible with a single track. The search stops with the first found 
# track candidate.
# 
class FastHoughTrans(object):

  ## Constructor.
  #
  # Construct hyper planes and calculate entries in root cube.
  #  
  # @param[in] setup   number of parameters in (XY, XZ) 
  # @param[in] hitInfo hit information (index, row, u, v, w)
  #
  def __init__(self, setup, hitInfo):
    ## setup
    self.__setup = setup
    ## hyper planes
    self.__planes = []
    ## entries (distances and hyper plane index)
    self.__entries = []
        
    entries = []
    nDim = self.__setup[0] + self.__setup[1]
    numChilds = 1 << nDim
    #print " dim ", nDim, numChilds
    for i, info in enumerate(hitInfo):
        # map from feature to parameter space
        inside = True
        rad = []
        directions = []
        for p in range(2):
          if self.__setup[p] > 0:
            u = info[2]; v = info[3 + p] 
            a = [ 1.0 ]  # (a1)
            if  self.__setup[p] > 1:
              a.append((u))  # (a1, a2)
              if  self.__setup[p] > 2:
                a.append((1.5 * u * u - 0.5))  # (a1, a2, a3)
            q2 = 0.
            for x in a:
              q2 += x * x  
            q = math.sqrt(q2)
            r = -v / q  # a0
            for j in range(len(a)):
              a[j] /= q
            rad.append(r)
            directions.append(a)
            inside = inside and (abs(r * a[0]) < 0.5)  # better cut ?
            #print " p ", p, u, v, a, r, inside

        hp = HyperPlane(info[0], info[1], info[2:5], directions)
        #hp.dump()
        self.__planes.append(hp)
        if inside:
          entries.append((rad[:], i))

    # order by row
    self.__entries = sorted(entries, key=lambda e: self.__planes[e[1]].getRow())    
    #print " ordered entries " 
    #for r, i in self.__entries:
    #  print r, i, self.__planes[i].getRow()  

  ## Run fast Hough transformation.
  #
  # @param[in] rowScale row scale (for row combination to 'meta' rows)  
  # @param[in] minRows  minimum number of rows for a track candidate  
  # @param[in] effCut   effciency cut for a track candidate  
  # @param[in] purCut   purity cut for a track candidate
  # @param[in] minLevel minimum number of subdivision of root cube 
  # @param[in] maxLevel maxinum number of subdivision of root cube 
  # @param[in] maxCube  maxinum number child cubes to check 
  #
  def run(self, rowScale, minRows, effCut=0.9, purCut=1.1, minLevel=5, maxLevel=8, maxCube=250):
          
    root = HyperCube(self.__setup, self.__planes, self.__entries)
    #root.dump()    
    global nCube
    nCube = 0
    steering = (maxCube, minRows, minLevel, maxLevel, effCut, purCut, rowScale)
    cand = root.divide(steering) 
    print " cubes ", nCube

    return cand  
      
## Hyper plane
#
# Build from the unit normals of the (1 or 2) hyper planes belonging to a space point.
# Contains the steps for updating the distance to a cube center for splitting in 
# any direction and the cut for the intersection.
#    
class HyperPlane(object):
  
  ## constructor
  #
  # @param[in] ihit   hit index
  # @param[in] row    row
  # @param[in] pos    pos
  # @param[in] dirs   unit normals
  # @param[in] dCut   distance cut (allowing for some overlap with neighbor cubes)
  #
  def __init__(self, ihit, row, pos, dirs, dCut=0.75):
    ## hit (index)
    self.__ihit = ihit
    ## row (number)
    self.__row = row
    ## position
    self.__pos = pos
    ## direction
    self.__dir = dirs[:]
    ## steps
    self.__steps = []
    ## distance cut
    self.__cut = []
    for d in dirs:
      # cut
      self.__cut.append(dCut / abs(d[0]))  # d[0] is max component
      # steps  
      nSteps = 2 ** len(d)
      steps = [ 0. for i in range(nSteps) ]
      for i in range(nSteps):
        ii = i
        for a in d:
          steps[i] += 0.5 * a if ii % 2 == 1 else -0.5 * a
          ii /= 2
      self.__steps.append(steps)
 
  ## dump 
  def dump(self, verbose=False):
    print " Hyperplane ", self.__ihit, self.__row, self.__cut 
    print "   directions ", self.__dir
    if verbose: 
      print "   steps      ", self.__steps
 
  ## get position
  def getPos(self):
    return self.__pos
    
  ## get direction 
  #
  # @param[in] idim dimension (0 for v, 1 for w)
  # @return unit normal vector
  #  
  def getDir(self, idim):
    return self.__dir[idim]
     
  ## get steps
  #
  # @return unit normal vector
  #  
  def getSteps(self):
    return self.__steps

  ## get cut
  #
  # @return cut (for maximal component of distance to cube center)
  #  
  def getCut(self):
    return self.__cut

  ## get row   
  #
  # @return row number
  # 
  def getRow(self):
    return self.__row

## Hypercube
#
# Defined by list of space points inside cube (with distances, hyper plane index).
#    
class HyperCube(object):
  
  ## constructor
  #
  # @param[in] dim     number of parameters in (XY, XZ) 
  # @param[in] planes  list of hyper planes
  # @param[in] entries list of space points inside cube (with distances, hyper plane index)
  # @param[in,out] descent cube splitting path
  #
  def __init__(self, dim, planes, entries, descent=[]):
    ## dimension
    self.__dimension = dim
    ## number of childs
    self.__numChilds = 1 << (dim[0] + dim[1])
    ## hyper planes 
    self.__planes = planes
    ## entries (distances and index)
    self.__entries = entries
    ## descent
    self.__descent = descent
    #
    global nCube
    nCube += 1

  ## dump    
  def dump(self, flag=False):
    print "  HyperCube ", len(self.__entries), self.__descent, self.getCenter()
    if flag:
      print self.__entries  

  ## get number of entries
  #
  # @return number of entries
  #
  def getNumEntries(self):
    return len(self.__entries)
  
  ## divide cube
  #
  # Divide cube into two child cubes in each dimension. Order childs by number of rows 
  # and spread of distances. Check only childs containing some minimum number of rows.
  # Stop with first accepted track candidate (based on hit density and level).
  #
  # @param[in] steering steering parameters
  # @return None of cube with track candidate
  #
  def divide(self, steering):
    level = len(self.__descent)
    rowScale = steering[6]
    #print " level ", level, " childs ", self.__descent, self.center()
    childs = [ [] for i in range(self.__numChilds) ]
    rows = [ [0, -1] for i in range(self.__numChilds) ]
    sums = [ [0.001, 0., 0., 0., 0.] for i in range(self.__numChilds) ]
    # patterns
    for rOld, ip in self.__entries:
      hp = self.__planes[ip]
      row = hp.getRow() / rowScale
      cut = hp.getCut()
      steps = hp.getSteps()
      # factorized subspaces?
      if len(rOld) > 1:
        # two factorized subspaces
        rNew2 = []
        for b1, step1 in enumerate(steps[1]):
          rChild1 = rOld[1] * 2. + step1
          if abs(rChild1) > cut[1]:
            continue
          rNew2.append((b1, rChild1))        
        for b0, step0 in enumerate(steps[0]):
          rChild0 = rOld[0] * 2. + step0
          if abs(rChild0) > cut[0]:
            continue
          for b1, rChild1 in rNew2:
            b = b0 + (b1 << self.__dimension[0])
            childs[b].append(([ rChild0, rChild1 ], ip))  
            if row > rows[b][1]:
              rows[b][0] += 1; rows[b][1] = row
            sums[b][0] += 1.; sums[b][1] += rChild0; sums[b][2] += rChild0 * rChild0
            sums[b][3] += rChild1; sums[b][4] += rChild1 * rChild1
      else:
        # single subspace
        for b0, step0 in enumerate(steps[0]):
          rChild0 = rOld[0] * 2. + step0
          #if level<2:
          #  print " r ", level, row, rOld[0], b0, step0, rChild0, cut[0]
          if abs(rChild0) > cut[0]:
            continue
          b = b0
          childs[b].append(([ rChild0 ], ip))
          if row > rows[b][1]:
            rows[b][0] += 1; rows[b][1] = row              
          sums[b][0] += 1.; sums[b][1] += rChild0; sums[b][2] += rChild0 * rChild0 

    #for b in range(self.__numChilds):
    #  print " b ", level, b, rows[b][0]
    #  self.plot(self.__entries, 'b+')
    #  self.plotCenter([b])
    #  self.plot(childs[b], 'r*', True)
    #for b in sorted(range(self.__numChilds), key=lambda x: rows[x][0] - sums[x][2] / sums[x][0] + (sums[x][1] / sums[x][0]) ** 2, reverse=True):
    #  print " b ", b,  rows[b][0],   sums[b][2] / sums[b][0] - (sums[b][1] / sums[b][0]) ** 2  
    #for n in counts:
    #  maxCount = max(maxCount,n)
    #for b in sorted(range(self.__numChilds), key=lambda x: rows[x][0] - sums[x][2] / sums[x][0] + (sums[x][1] / sums[x][0]) ** 2, reverse=True):
    for b in sorted(range(self.__numChilds), key=lambda x: rows[x][0] - sums[x][2] / sums[x][0] - sums[x][4] / sums[x][0], reverse=True):
      nRow = rows[b][0]
      if nRow < steering[1]:
        return None
      newPlanes = childs[b]
      newCube = HyperCube(self.__dimension, self.__planes, newPlanes, self.__descent + [b])
      if nCube > steering[0]:
        return None  
      rowLength = self.__planes[newPlanes[-1][1]].getRow() - self.__planes[newPlanes[0][1]].getRow() + 1
      ratio = len(newPlanes) / float(rowLength)  # hit density (number of hits / row length) 
      #print "                  "[:level],
      #print " sum ", level, b, rows[b][0], len(newPlanes), rowLength, ratio  #, \
      #  sums[b][1] / sums[b][0], sums[b][2] / sums[b][0] - (sums[b][1] / sums[b][0]) ** 2, sums[b][2] / sums[b][0], \
      #  sums[b][3] / sums[b][0], sums[b][4] / sums[b][0] - (sums[b][3] / sums[b][0]) ** 2, sums[b][4] / sums[b][0]
      if level >= steering[2]:
        if steering[4] < ratio < steering[5]:
          return newCube
          #continue   
      if level < steering[3]:            
        cand = newCube.divide(steering)
        if cand is not None:
          return cand    
 
  ## plot
  def plot(self, entries, mode, show=False):
    X = []; Y = []; Z = []
    for r, ip in entries:
      hp = self.__planes[ip]
      u, v, w = hp.getPos()
      X.append(u)
      Y.append(v)
      Z.append(w)
      
    plt.figure(1)
    plt.subplot(211)      
    plt.plot(X, Y, mode)
    #plt.title('V vs U')
    plt.subplot(212)            
    plt.plot(X, Z, mode)    
    #plt.title('W vs U')
    if show: 
      plt.show()   
 
  ## plot center
  def plotCenter(self, b=[]):
    center = self.getCenter(b)
    #print " center ", center
    X = []; Y = []; Z = []
    for i in range(21):
      u = i * 0.1 - 1.
      lp = [1., u, 1.5 * u * u - 0.5]
      v = 0.
      for d in range(self.__dimension[0]):
        v += center[d] * lp[d]
      w = 0.
      for d in range(self.__dimension[1]):
        w += center[d + self.__dimension[0]] * lp[d]        
      X.append(u)
      Y.append(v)  
      Z.append(w)  
    plt.subplot(211)
    if self.__dimension[0] > 0 :plt.plot(X, Y, 'm--')  
    string = "cube "
    for i in self.__descent + b: string += str(i) + " "
    plt.title(string) 
    plt.subplot(212)
    if self.__dimension[1] > 0 :plt.plot(X, Z, 'm--')  
                        
  ## get center of cube 
  #
  # @return center of cube
  #    
  def getCenter(self, b=[]):
    nDim = self.__dimension[0] + self.__dimension[1]
    best = [0. for i in range(nDim)]
    for child in reversed(self.__descent + b):
      c = child
      for i in range(nDim):
        best[i] = 0.5 * best[i] + 0.25 * (2 * (c % 2) - 1)
        c /= 2
    return best 

  ## get indices
  #
  # @return hyper plane indices
  #
  def getIndices(self):
    return [ p[1] for p in self.__entries]
 
  ## get distances 
  #
  # @return distances
  #  
  def getDistances(self):
    return [ p[0] for p in self.__entries]
 
  ## get entries
  #
  # @return entries (distances, hyper planes index)
  #  
  def getEntries(self):
    return self.__entries

  ## get level  
  #
  # @return cube level (number of splittings)
  #   
  def getLevel(self):
    return len(self.__descent)
  
