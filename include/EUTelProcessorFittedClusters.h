// Version: $Id: EUTelDafBase.h 2925 2013-09-02 11:02:00Z hamnett $
//! Author Havard Gjersdal <haavagj@fys.uio.no>
/*
 *   This source code is part of the Eutelescope package of Marlin.
 */

#ifndef EUTELDAFBASE_H
#define EUTELDAFBASE_H

// built only if GEAR is available
#ifdef USE_GEAR

// eutelescope includes
#include "EUTelAlignmentConstant.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>

// system includes <>
#include <string>
#include <fstream>
#include <vector>
#include <map>

namespace eutelescope {
  class ClusterInfo;
  class EUTelFittedClusters : public marlin::Processor {
  public:
    EUTelFittedClusters ();
    EUTelFittedClusters (std::string);
    virtual void init ();
    virtual void processRunHeader (LCRunHeader * run);
    virtual void processEvent (LCEvent * evt);
    virtual void end();
    virtual int guessSensorID( double* hit);

  protected:
    LCCollectionVec *_trackCollection;
    LCCollectionVec *_clusterCollection;
    std::string _trackCollectionName;
    std::string _clusterCollectionName;
    double _minBinX;
    double _minBinY;
    double _maxBinX;
    double _maxBinY;
    double _binSizeX;
    double _binSizeY;

    //! Silicon planes parameters as described in GEAR
    gear::SiPlanesParameters * _siPlanesParameters;
    gear::SiPlanesLayerLayout * _siPlanesLayerLayout;

    //The following map represents each plane
    std::map< int, ClusterInfo* > _clusteringInfomationForEachPlane;

    TH2D *_clusterSizeTotalVsPositionXAllPlanesCombined;
    TH2D *_clusterSizeTotalVsPositionYAllPlanesCombined;
    TH2D *_clusterSizeXVsPositionXAllPlanesCombined;
    TH2D *_clusterSizeXVsPositionYAllPlanesCombined;
    TH2D *_clusterSizeYVsPositionXAllPlanesCombined;
    TH2D *_clusterSizeYVsPositionYAllPlanesCombined;
    std::vector< TH2D* > _clusterSizeTotalVsPositionXIndividualPlanes;
    std::vector< TH2D* > _clusterSizeTotalVsPositionYIndividualPlanes;
    std::vector< TH2D* > _clusterSizeXVsPositionXIndividualPlanes;
    std::vector< TH2D* > _clusterSizeXVsPositionYIndividualPlanes;
    std::vector< TH2D* > _clusterSizeYVsPositionXIndividualPlanes;
    std::vector< TH2D* > _clusterSizeYVsPositionYIndividualPlanes;
    TH2D *_clusterSizeTotalVsResidualXAllPlanesCombined;
    TH2D *_clusterSizeTotalVsResidualYAllPlanesCombined;
    TH2D *_clusterSizeXVsResidualXAllPlanesCombined;
    TH2D *_clusterSizeXVsResidualYAllPlanesCombined;
    TH2D *_clusterSizeYVsResidualXAllPlanesCombined;
    TH2D *_clusterSizeYVsResidualYAllPlanesCombined;
    std::vector< TH2D* > _clusterSizeTotalVsResidualXIndividualPlanes;
    std::vector< TH2D* > _clusterSizeTotalVsResidualYIndividualPlanes;
    std::vector< TH2D* > _clusterSizeXVsResidualXIndividualPlanes;
    std::vector< TH2D* > _clusterSizeXVsResidualYIndividualPlanes;
    std::vector< TH2D* > _clusterSizeYVsResidualXIndividualPlanes;
    std::vector< TH2D* > _clusterSizeYVsResidualYIndividualPlanes;
    TH2D *_clusterSizeTotalVsChi2;
    TH2D *_clusterSizeXVsChi2;
    TH2D *_clusterSizeYVsChi2;
   
    TH1D *_averageClusterSizeTotalVsPositionXAllPlanesCombined;
    TH1D *_averageClusterSizeTotalVsPositionYAllPlanesCombined;
    TH1D *_averageClusterSizeXVsPositionXAllPlanesCombined;
    TH1D *_averageClusterSizeXVsPositionYAllPlanesCombined;
    TH1D *_averageClusterSizeYVsPositionXAllPlanesCombined;
    TH1D *_averageClusterSizeYVsPositionYAllPlanesCombined;
    std::vector< TH1D* > _averageClusterSizeTotalVsPositionXIndividualPlanes;
    std::vector< TH1D* > _averageClusterSizeTotalVsPositionYIndividualPlanes;
    std::vector< TH1D* > _averageClusterSizeXVsPositionXIndividualPlanes;
    std::vector< TH1D* > _averageClusterSizeXVsPositionYIndividualPlanes;
    std::vector< TH1D* > _averageClusterSizeYVsPositionXIndividualPlanes;
    std::vector< TH1D* > _averageClusterSizeYVsPositionYIndividualPlanes;
    TH1D *_averageClusterSizeTotalVsResolutionXAllPlanesCombined;
    TH1D *_averageClusterSizeTotalVsResolutionYAllPlanesCombined;
    TH1D *_averageClusterSizeXVsResolutionXAllPlanesCombined;
    TH1D *_averageClusterSizeXVsResolutionYAllPlanesCombined;
    TH1D *_averageClusterSizeYVsResolutionXAllPlanesCombined;
    TH1D *_averageClusterSizeYVsResolutionYAllPlanesCombined;
    std::vector< TH1D* > _averageClusterSizeTotalVsResolutionXIndividualPlanes;
    std::vector< TH1D* > _averageClusterSizeTotalVsResolutionYIndividualPlanes;
    std::vector< TH1D* > _averageClusterSizeXVsResolutionXIndividualPlanes;
    std::vector< TH1D* > _averageClusterSizeXVsResolutionYIndividualPlanes;
    std::vector< TH1D* > _averageClusterSizeYVsResolutionXIndividualPlanes;
    std::vector< TH1D* > _averageClusterSizeYVsResolutionYIndividualPlanes;
    TH1D *_averageClusterSizeTotalVsChi2;
    TH1D *_averageClusterSizeXVsChi2;
    TH1D *_averageClusterSizeYVsChi2;

    TH2D *_clusterSizeTotalVsAngleX;
    TH2D *_clusterSizeXVsAngleX;
    TH2D *_clusterSizeYVsAngleX;
    TH2D *_clusterSizeTotalVsAngleY;
    TH2D *_clusterSizeXVsAngleY;
    TH2D *_clusterSizeYVsAngleY;
    TH1D *_averageClusterSizeTotalVsAngleX;
    TH1D *_averageClusterSizeXVsAngleX;
    TH1D *_averageClusterSizeYVsAngleX;
    TH1D *_averageClusterSizeTotalVsAngleY;
    TH1D *_averageClusterSizeXVsAngleY;
    TH1D *_averageClusterSizeYVsAngleY;
  };

class ClusterInfo{
  public:
    //The following maps have the integer as the cluster size and the vector as the name of the variable
    std::map< int, std::vector < double > > _xPositionVsSizeTotal;
    std::map< int, std::vector < double > > _yPositionVsSizeTotal;
    std::map< int, std::vector < double > > _chi2VsSizeTotal;
    std::map< int, std::vector < double > > _angleXVsSizeTotal;
    std::map< int, std::vector < double > > _angleYVsSizeTotal;
    std::map< int, std::vector < double > > _residualXVsSizeTotal;
    std::map< int, std::vector < double > > _residualYVsSizeTotal;
    std::map< int, std::vector < double > > _xPositionVsSizeX;
    std::map< int, std::vector < double > > _yPositionVsSizeX;
    std::map< int, std::vector < double > > _chi2VsSizeX;
    std::map< int, std::vector < double > > _angleXVsSizeX;
    std::map< int, std::vector < double > > _angleYVsSizeX;
    std::map< int, std::vector < double > > _residualXVsSizeX;
    std::map< int, std::vector < double > > _residualYVsSizeX;
    std::map< int, std::vector < double > > _xPositionVsSizeY;
    std::map< int, std::vector < double > > _yPositionVsSizeY;
    std::map< int, std::vector < double > > _chi2VsSizeY;
    std::map< int, std::vector < double > > _angleXVsSizeY;
    std::map< int, std::vector < double > > _angleYVsSizeY;
    std::map< int, std::vector < double > > _residualYVsSizeY;
    std::map< int, std::vector < double > > _residualYVsSizeY;
};

}
#endif
#endif

