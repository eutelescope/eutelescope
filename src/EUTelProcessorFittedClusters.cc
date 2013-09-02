// Version: $Id: EUTelFittedClusters.cc 2925 2013-09-02 11:02:00Z hamnett $

// Author Havard Gjersdal, UiO(haavagj@fys.uio.no)
/*!
 * This is a track fitting processor for the Eutelescope package. 
 *
 * It preforms track finding and fitting on a supplied hit collection.
 *
 * The track finder works by propagating all hits to plane 0, currently assuming straight
 * line fits, then running a cluster finder. Hit clusters above some set value are considered
 * track candidates.
 *
 * This track candidate is then fitted using a implementation of a Deterministic Annealing
 * Filter (DAF), that in short is a Kalman Filter running iteratively over a set of weighted
 * measurements, reweighing the measurements after each fit based on the residuals and a
 * supplied chi2 cut off.
 *
 * This package uses the Eigen library for linear algebra. This package is very quick when
 * compiled properly, but very slow when compiled for debugging. Make sure to compile
 * properly before running productions.
 *
 * Running 'cmake -i' inside the build folder, and then when it asks
 * Variable Name: CMAKE_CXX_FLAGS_RELEASE
 * Description: Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files).                          
 *
 * enter:
 * New Value (Enter to keep current value): -O3 -msse2 -ftree-vectorize -DNDEBUG
 *
 * When it asks
 * Variable Name: CMAKE_BUILD_TYPE
 * enter:
 * New Value (Enter to keep current value): Release
 *
 * If youc cpu supports it, you could try -msse4 or -msse3 aswell.
 */

// built only if GEAR and MARLINUTIL are used
#if defined(USE_GEAR)
// eutelescope includes ".h"
#include "EUTelFittedClusters.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelExceptions.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelReferenceHit.h"


// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"
#include "marlin/Exceptions.h"
#include "marlin/AIDAProcessor.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// lcio includes <.h>
#include <IO/LCWriter.h>
#include <UTIL/LCTime.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCObject.h>
#include <EVENT/TrackerPulse.h>
#include <EVENT/TrackerData.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/SimTrackerHitImpl.h>


#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>

// system includes <>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>

#include <TVector3.h>
#include <TH1D.h>
#include <TF1.h>
#include <Eigen/Geometry> 

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;


EUTelFittedClusters::EUTelFittedClusters(std::string name) : marlin::Processor(name)
: _trackCollection(NULL),
  _clusterCollection(NULL),
  _trackCollectionName(""),
  _clusterCollectionName(""),
  _minBinX(0.0),
  _minBinY(0.0),
  _maxBinX(0.0),
  _maxBinY(0.0),
  _binSizeX(0.0),
  _binSizeY(0.0){
  registerProcessorParameter( "TrackCollectionName", "Names of input track collection", _trackCollectionName,"default");
  registerProcessorParameter ("ClusterCollectionName", "Name of the cluster collection to use", _clusterCollectionName, "default" );
  registerOptionalParameter("MinimumBinX","The smallest size in the x direction", _minBinX, -10.75);
  registerOptionalParameter("MaximumBinX","The largest size in the x direction", _maxBinX, 10.75);
  registerOptionalParameter("MinimumBinY","The smallest size in the y direction", _minBinY, -6.85);
  registerOptionalParameter("MaximumBinY","The largest size in the y direction", _maxBinY, 6.85);
  registerOptionalParameter("BinSizeX","The bin size in the x direction", _binSizeX, 0.1);
  registerOptionalParameter("BinSizeY","The bin size in the y direction", _binSizeY, 0.05);
}

void EUTelFittedClusters::init() {
  
}

void EUTelFittedClusters::processRunHeader (LCRunHeader * rdr) {
  auto_ptr<EUTelRunHeaderImpl> header ( new EUTelRunHeaderImpl (rdr) );
  header->addProcessor( type() ) ;
  ++_iRun;
}

int EUTelFittedClusters::guessSensorID( double * hit ) 
{

  int sensorID = -1;
  double minDistance =  numeric_limits< double >::max() ;
  if( ReferenceHitVecIsSet() )
  {
    streamlog_out( MESSAGE5 ) << "_referenceHitVec is empty" << endl;
    return 0;
  }

      for(size_t ii = 0 ; ii <  static_cast< unsigned int >(_referenceHitVec->getNumberOfElements()); ii++)
      {
        EUTelReferenceHit* refhit = static_cast< EUTelReferenceHit*> ( _referenceHitVec->getElementAt(ii) ) ;
        TVector3 hit3d( hit[0], hit[1], hit[2] );
        TVector3 hitInPlane( refhit->getXOffset(), refhit->getYOffset(), refhit->getZOffset());
        TVector3 norm2Plane( refhit->getAlpha(), refhit->getBeta(), refhit->getGamma() );
 
        double distance = abs( norm2Plane.Dot(hit3d-hitInPlane) );
        if ( distance < minDistance ) 
        {
           minDistance = distance;
           sensorID = ii;                    // number in the GEAR file z ordered
        }    

      }

  return sensorID;
}


void EUTelFittedClusters::processEvent(LCEvent * event){

  try{
    _clusterCollection = dynamic_cast < LCCollectionVec * > (event->getCollection( _clusterCollectionName));
    _trackCollection = dynamic_cast < LCCollectionVec * > (event->getCollection( _trackCollectionName));
    size_t numberoftracks = _trackCollection->getNumberOfElements();
    for(size_t i = 0; i < numberoftracks; ++i){
      fillPlots(track);
    }
  } catch(...){
    streamlog_out(MESSAGE2) << "No cluster in this event" << endl;
  }

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);
  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG2 ) << "EORE found: nothing else to do." << endl;
    return;
  }

}

void EUTelFittedClusters::fillPlots(EVENT::Track *track){
  vector< TrackerHit* > hits = track->getTrackerHits();
  double chi2 = track->getChi2();
  double ndof = track->getNdf();
  size_t tracksize = hits.size()/2;
  double clustersizetotal(0);
  double clustersizex(0);
  double clustersizey(0);
  for(size_t z = 0; z < tracksize; ++z){
    for(int i = 0; i < _clusterCollection->getNumberOfElements(); ++i){
      TrackerPulse * cluster = dynamic_cast< TrackerPulse* >(_clusterCollection->getElementAt(i));
      streamlog_out(DEBUG0) << "Element " << i << " has CellID1: " << cluster->getCellID1() << endl;
      std::vector< float > covmatrix = cluster->getCovMatrix();
      for(size_t j = 0; j < covmatrix.size(); ++j){
        streamlog_out(DEBUG0) << covmatrix[j] << ", ";
      }
      streamlog_out(DEBUG0) << endl;
      TrackerData *clusterdata = dynamic_cast< TrackerData* >(cluster->getTrackerData());
      int   rhs = 0;
      lcio::long64 mask  = 0x1F;
      lcio::long64 cell0 = static_cast<lcio::long64> (clusterdata->getCellID0());
      streamlog_out(DEBUG0) << "The sensor ID is: " << static_cast<int> ( ( cell0 & mask ) >> rhs ) << endl;
      streamlog_out(DEBUG0) << "clusterdata has CellID0: " << clusterdata->getCellID0() << endl;
      streamlog_out(DEBUG0) << "clusterdata has CellID1: " << clusterdata->getCellID1() << endl;
      streamlog_out(DEBUG0) << "clusterdata has time: " << clusterdata->getTime() << endl;
      std::vector< float > chargevalues = clusterdata->getChargeValues();

      double clusterseedx(0);
      double clusterseedy(0);
      double ClusterHitCut = 0.0184*20;
      streamlog_out(DEBUG0) << "chargevalues has size: " << chargevalues.size() << endl;
      clustersizetotal = static_cast<double>(chargevalues.size())/3.0;
      for(size_t j = 0; j < chargevalues.size(); ++j){
        if((j+1)%3 != 0){
          if((j+1) % 3 == 2){
            if(clusterseedy == 0){
              clusterseedy = 576*0.0184/2.0 - chargevalues[j]*0.0184;
            }
            streamlog_out(DEBUG0) << "y = " << chargevalues[j]*0.0184 - 576*0.0184/2.0 << ", ";//This is where the hit positions are actually stored
          } else{
            if(clusterseedx == 0){
              clusterseedx = 1152*0.0184/2.0 - chargevalues[j]*0.0184;
            }
            streamlog_out(DEBUG0) << "x = " << chargevalues[j]*0.0184 - 1152*0.0184/2.0 << ", ";//This is where the hit positions are actually stored
          }
        } else{
          streamlog_out(DEBUG0) << endl;
        }
      }
    }
    _clusteringInfomationForEachPlane[z]->_xPositionVsSizeTotal[clustersizetotal].push_back(xpos);
  }
}

double getAverage(std::vector< double > z){
  double mean(0);  
  for(vector< double >::iterator it = z.begin(); it != z.end(); ++it){
    mean += *it;
  }
  mean /= static_cast<double>(z.size());
  return mean;
}

double getError(std::vector< double > z, double mean){
  double error(0);
  for(vector<double>::iterator it = z.begin(); it != z.end(); ++it){
    error += (mean - *it)*(mean - *it);
  }
  error = sqrt(error);
  error /= sqrt(z.size());
  return error;
}

double GetAverageResolution(std::vector< double > averageresidual){
  TH1D *histogram = new TH1D("histogram","histogram",100,-0.04,0.04);
  for(vector< double >::iterator i = averageresidual.begin(); i != averageresidual.end(); ++i){
    histogram->Fill(*i);
  }
  TF1 *fit = new TF1("fit","gaus");
  histogram->Fit(fit,"Q");
  return fit->GetParameter(2);
}

void EUTelFittedClusters::end() {
}
#endif // USE_GEAR
