#ifndef cluster_h
#define cluster_h 1

#include <iostream>
#include <vector>
#include <algorithm>
#include <map>

using namespace std;

class Cluster {
  public:
    Cluster();
    void set_values(int,vector<int>,vector<int>);
    Cluster mirrorX();
    Cluster mirrorY();
    Cluster rotate90();
    void NeighbourPixels(int x, int y, vector<int> xOriginal, vector<int> yOriginal, vector<int> &xNeighbour, vector<int> &yNeighbour);
    void FindReferenceClusters(vector<Cluster> &clusterVec, int sizeMax);
    std::map<int,int> SymmetryPairs(vector<Cluster> clusterVec, const char* type);
    vector< vector<int> > sameShape(vector<Cluster> clusterVec);
    int WhichClusterShape(Cluster cluster, vector<Cluster> clusterVec);
    vector<int> getX() {return x;}
    vector<int> getY() {return y;}
    int Size() {return size;}
    bool operator==(Cluster c2); 
    void getCenterOfGravity(float &xCenter, float &yCenter);
 protected:
    int size;
    vector<int> x;
    vector<int> y;
};

#endif
