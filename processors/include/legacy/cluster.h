#ifndef cluster_h
#define cluster_h 1

#include <algorithm>
#include <iostream>
#include <map>
#include <vector>

class Cluster {
public:
  Cluster();
  void set_values(size_t, std::vector<int>, std::vector<int>);
  Cluster mirrorX();
  Cluster mirrorY();
  Cluster rotate90();
  void NeighbourPixels(int x, int y, std::vector<int> xOriginal,
                       std::vector<int> yOriginal, std::vector<int> &xNeighbour,
                       std::vector<int> &yNeighbour);
  void FindReferenceClusters(std::vector<Cluster> &clusterVec, size_t sizeMax);
  std::map<int, int> SymmetryPairs(std::vector<Cluster> clusterVec,
                                   const char *type);
  std::vector<std::vector<int>> sameShape(std::vector<Cluster> clusterVec);
  int WhichClusterShape(Cluster cluster, std::vector<Cluster> clusterVec);
  std::vector<int> getX() { return x; }
  std::vector<int> getY() { return y; }
  size_t Size() { return size; }
  bool operator==(Cluster c2);
  void getCenterOfGravity(float &xCenter, float &yCenter);

protected:
  size_t size;
  std::vector<int> x;
  std::vector<int> y;
};

#endif
