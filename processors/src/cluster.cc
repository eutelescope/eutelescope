#ifdef USE_GEAR

#include "cluster.h"

using namespace std;

Cluster aCluster;

Cluster::Cluster()
    :size(0),
    x(0),
    y(0){
} 

bool Cluster::operator==(Cluster c2) {
  int x1Min = *min_element(x.begin(), x.end());                       
  int x2Min = *min_element(c2.x.begin(), c2.x.end());                 
  int y1Min = *min_element(y.begin(), y.end());                       
  int y2Min = *min_element(c2.y.begin(), c2.y.end());                 
  if (size != c2.Size()) return false;                                
  bool isEqual = false;                                               
  for (int i=0; i<size; i++)                                          
  {                                                                   
    for (int j=0; j<size; j++)                                        
    {                                                                 
      if (x[i]-x1Min == c2.x[j]-x2Min && y[i]-y1Min == c2.y[j]-y2Min) 
      {                                                               
        isEqual = true;                                               
        break;                                                        
      }                                                               
    }                                                                 
    if (isEqual)                                                      
    {                                                                 
      isEqual = false;                                                
      continue;                                                       
    }                                                                 
    else return false;                                                
  }                                                                   
  return true;                                                        
}                                                                     



Cluster Cluster::mirrorX() {
  vector<int> xNew(size);
  vector<int> yNew(size);
  int yMax = *max_element(y.begin(), y.end());
  for (int i=0; i<size; i++)
  {
    xNew[i]=x[i];
    yNew[i]=yMax-y[i];
  }
  Cluster clusterNew;
  clusterNew.set_values(size,xNew,yNew);
  return clusterNew;
}

Cluster Cluster::mirrorY() {
  vector<int> xNew(size);
  vector<int> yNew(size);
  int xMax = *max_element(x.begin(), x.end());
  for (int i=0; i<size; i++)
  {
    xNew[i]=xMax-x[i];
    yNew[i]=y[i];
  }
  Cluster clusterNew;
  clusterNew.set_values(size,xNew,yNew);
  return clusterNew;
}

Cluster Cluster::rotate90(){
  vector<int> xNew(size);
  vector<int> yNew(size);
  for (int i=0; i<size; i++)
  {
    xNew[i]=-1*y[i];
    yNew[i]=x[i];
  }
  int xMin = *min_element(xNew.begin(), xNew.end());
  int yMin = *min_element(yNew.begin(), yNew.end());
  for (int i=0; i<size; i++)
  {
    xNew[i]-=xMin;
    yNew[i]-=yMin;
  }
  Cluster clusterNew;
  clusterNew.set_values(size,xNew,yNew);
  return clusterNew;

}


void Cluster::set_values (int s, vector<int> a, vector<int> b) {
  size = s;
  if (x.size() == 0)
    for (int i=0; i<size; i++)
    {
      x.push_back(a[i]);
      y.push_back(b[i]);
    }
  else 
    for (int i=0; i<size; i++)
      if (x.size() <= (unsigned int)size)
      {
        x[i] = a[i];
        y[i] = b[i];
      }
      else
      {
        x.push_back(a[i]);
        y.push_back(b[i]);
      }
}

void Cluster::NeighbourPixels(int x, int y, vector<int> xOriginal, vector<int> yOriginal, vector<int> &xNeighbour, vector<int> &yNeighbour)
{
  int yTmp = 0;
  for (int xTmp=x-1; xTmp<=x+1; xTmp++)
  {
    bool pixelExists = false;
    if (xTmp==x-1 || xTmp==x+1) 
    {
      yTmp = y;
      for (unsigned int i=0; i<xOriginal.size(); i++)
      {
        if (xTmp == xOriginal[i] && yTmp == yOriginal[i]) 
        {
          pixelExists = true;
          break;
        }
      }
      if (pixelExists) continue;
      xNeighbour.push_back(xTmp);
      yNeighbour.push_back(yTmp);
      pixelExists = false;
    }
    else
    {
      for (int j=0; j<2; j++)
      {
        yTmp = (j==0?y+1:y-1);
        for (unsigned int i=0; i<xOriginal.size(); i++)
        {
          if (xTmp == xOriginal[i] && yTmp == yOriginal[i]) 
          {
            pixelExists = true;
            break;
          }
        }
        if (pixelExists) 
        {
          pixelExists = false; 
          continue;
        }
        xNeighbour.push_back(xTmp);
        yNeighbour.push_back(yTmp);
        pixelExists = false;
      }
    }
  }
}

void Cluster::FindReferenceClusters(vector<Cluster> &clusterVec, int sizeMax)
{
  vector<int> xTmp(1);
  vector<int> yTmp(1);
  xTmp[0]=0;
  yTmp[0]=0;
  Cluster cTmp;
  cTmp.set_values(1,xTmp,yTmp);
  clusterVec.push_back(cTmp);
  for (int size=2; size<=sizeMax; size++)
  {
    cout << "Looking for clusters with size: " << size << endl;
    int clusterVecSize = clusterVec.size();
    for (int iCluster=0; iCluster<clusterVecSize; iCluster++)
    {
      if (clusterVec[iCluster].Size() < size-1) continue;
      Cluster cluster; 
      vector<int> x = clusterVec[iCluster].getX();
      vector<int> y = clusterVec[iCluster].getY();
      for (int iPixel=0; iPixel<clusterVec[iCluster].Size(); iPixel++)
      {
        vector<int> xNeighbour;
        vector<int> yNeighbour;
        NeighbourPixels(x[iPixel],y[iPixel],x,y,xNeighbour,yNeighbour);
        for (unsigned int i=0; i<xNeighbour.size();i++)
        {
          bool areadyExists = false;
          if (x.size() < (unsigned int)size) 
          {
            x.push_back(xNeighbour[i]);
            y.push_back(yNeighbour[i]);
          }
          else
          {
            x[size-1] = xNeighbour[i];
            y[size-1] = yNeighbour[i];
          }
          cluster.set_values(size,x,y);
          for (unsigned int k=0; k<clusterVec.size(); k++)
            if (cluster == clusterVec[k]) 
            {
              areadyExists = true; 
              break;
            }
          if (!areadyExists)
            clusterVec.push_back(cluster);
        }
      }
    }
  }
  for (unsigned int iCluster = 0; iCluster<clusterVec.size();iCluster++)
  {
    Cluster cluster = clusterVec[iCluster];
    vector<int> X = cluster.getX();
    int minX = *min_element(X.begin(),X.end());
    if (minX < 0 )
      for (unsigned int iPixel=0; iPixel<X.size(); iPixel++)
        X[iPixel] -= minX;
    vector<int> Y = cluster.getY();
    int minY = *min_element(Y.begin(),Y.end());
    if (minY < 0 )
      for (unsigned int iPixel=0; iPixel<Y.size(); iPixel++)
        Y[iPixel] -= minY;
    cluster.set_values(X.size(),X,Y);
    clusterVec[iCluster] = cluster;
    cout << "Shapes: ID " << iCluster << " X size " << X.size() << " Y size " << Y.size() << endl; 
    cout << "(X,Y) :" ;
    for (unsigned int iPixel=0; iPixel<X.size(); iPixel++) cout << "(" << X[iPixel] << "," << Y[iPixel] << ") " ;
    cout << endl;
  }
  cout << "All shapes found!" << endl;
}

std::map<int,int> Cluster::SymmetryPairs(vector<Cluster> clusterVec, const char* type){
  std::map<int,int> pair;
  const char* typeX = "x";
  const char* typeY = "y";
  if (type != typeX && type != typeY)   
  {
    cerr << "Type has to be y or x, assuming x" << endl;
    type = "x";
  }
  for (unsigned int i=0; i<clusterVec.size(); i++)
  {
    Cluster cluster;
    if (type == typeX) cluster = clusterVec[i].mirrorX();
    else if (type == typeY) cluster = clusterVec[i].mirrorY();
    for (unsigned int j=0; j<clusterVec.size(); j++)
    {
      bool alreadyAdded = false;
      for(map<int,int>::iterator it = pair.begin(); it != pair.end(); ++it)
        if (it->first == (int)j)
        {
          alreadyAdded = true;
          break;
        }
      if (alreadyAdded) continue;
      if (cluster == clusterVec[j]) 
      {
        if (i == j) break;
        pair.insert(make_pair(i,j));
      }
    }
  }
  return pair;
}

vector< vector<int> > Cluster::sameShape(vector<Cluster> clusterVec){
  vector< vector<int> > symmetryGroups;
  for (unsigned int i=0; i<clusterVec.size(); i++)
  {
    vector<int> group;
    Cluster clusterOriginal = clusterVec[i];
    bool alreadyAdded = false;
    for (unsigned int k=0; k<symmetryGroups.size(); k++)
    {
      for (unsigned int l=0; l<symmetryGroups[k].size(); l++)
        if (symmetryGroups[k][l] == (int)i)
        {
          alreadyAdded = true;
          break;
        }
      if (alreadyAdded) break;
    }
    if (alreadyAdded) continue;
    else group.push_back(i);
    Cluster cluster2;
    for (int mir=0; mir<2; mir++)
    {
      if (mir == 0)  cluster2 = clusterOriginal;
      else  cluster2 = clusterOriginal.mirrorX();
      for (int rot=0; rot<4; rot++)
      {
        cluster2 = cluster2.rotate90();
        if (cluster2 == clusterOriginal) break;
        for (unsigned int j=0; j<clusterVec.size(); j++)
        {
          if (clusterVec[j] == cluster2)
          {
            alreadyAdded = false;
            for (unsigned int k=0; k<group.size(); k++)
              if (group[k] == (int)j)
              {
                alreadyAdded = true;
                break;
              }
            for (unsigned int k=0; k<symmetryGroups.size() && !alreadyAdded; k++)
            {
              for (unsigned int l=0; l<symmetryGroups[k].size(); l++)
                if (symmetryGroups[k][l] == (int)j)
                {
                  alreadyAdded = true;
                  break;
                }
              if (alreadyAdded) break;
            }
            if (!alreadyAdded) group.push_back(j);
          }
        }
      }
    }
    sort(group.begin(), group.end());
    symmetryGroups.push_back(group);
  }
  return symmetryGroups;
}


int Cluster::WhichClusterShape(Cluster cluster, vector<Cluster> clusterVec)
{
  for (unsigned int iCluster=0; iCluster<clusterVec.size(); iCluster++)
    if (cluster == clusterVec[iCluster]) return iCluster;
  return -1;
}

void Cluster::getCenterOfGravity(float &xCenter, float &yCenter)
{
  xCenter = 0;
  yCenter = 0;
  for (unsigned int iX=0; iX<x.size(); iX++)
    xCenter += x[iX];
  xCenter = xCenter/x.size();
  for (unsigned int iY=0; iY<y.size(); iY++)
    yCenter += y[iY];
  yCenter = yCenter/y.size();
}
#endif
