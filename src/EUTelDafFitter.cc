// Version: $Id$
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
#include "EUTelDafFitter.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelExceptions.h"
#include "EUTelReferenceHit.h"


// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"
#include "marlin/Exceptions.h"
#include "marlin/AIDAProcessor.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <IO/LCWriter.h>
#include <UTIL/LCTime.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <EVENT/TrackerPulse.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
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

// ROOT includes
#include "TVector3.h"

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;


EUTelDafFitter::EUTelDafFitter () : EUTelDafBase("EUTelDafFitter"){
  //cout << "DafFitter " << "51" << endl;
    //Child spesific params and description
  dafParams();
  //cout << "DafFitter " << "52" << endl;
}

void EUTelDafFitter::dafParams(){
  //cout << "DafFitter " << "53" << endl;
  _description = "This processor preforms track reconstruction. The tracks are as final track fit for analysis.";

  //Tracker system options
  registerOptionalParameter("AddToLCIO", "Should plots be made and filled?", _addToLCIO, static_cast<bool>(true));
  registerOptionalParameter("FitDuts","Set this to true if you want DUTs to be included in the track fit", _fitDuts, static_cast<bool>(false)); 
  //Track fitter options
  registerOutputCollection(LCIO::TRACK,"TrackCollectionName", "Collection name for fitted tracks", _trackCollectionName, string ("fittracks"));
  //cout << "DafFitter " << "54" << endl;
}

void EUTelDafFitter::dafInit() {
  //cout << "DafFitter " << "55" << endl;
//printf("EUTelDafFitter::dafInit()\n"); 
  if(_fitDuts){
    for( size_t ii = 0; ii< _system.planes.size(); ii++){
      if( find(_dutPlanes.begin(), _dutPlanes.end(), _system.planes.at(ii).getSensorID()) != _dutPlanes.end()){ 
	_system.planes.at(ii).include(); 
      }
    }
  }

//  _isFirstEvent = true; 

  //cout << "DafFitter " << "56" << endl;
}


void EUTelDafFitter::dafEvent (LCEvent * event) {
  //cout << "DafFitter " << "57" << endl;
// printf("EUTelDafFitter::dafEvent()\n"); 
/*
  if( _isFirstEvent )
  {
    try
    {
     _referenceHitVec     = dynamic_cast < LCCollectionVec * > (event->getCollection( _referenceHitCollectionName));  
    }
    catch(...)
    {
      streamlog_out( ERROR5 ) << "Critical error; the referennce hit collection was not found, pls check your steering file." << endl;
    }
 
  }
*/

  //Prepare track collection
  if(_addToLCIO){
    // Define output track and hit collections
    _fittrackvec = new LCCollectionVec(LCIO::TRACK);
    _fitpointvec = new LCCollectionVec(LCIO::TRACKERHIT);
    // Set flag for storing track hits in track collection
    LCFlagImpl flag(_fittrackvec->getFlag());
    flag.setBit( LCIO::TRBIT_HITS );
    _fittrackvec->setFlag(flag.getFlag());
  }
  
  //Check found tracks
  for(size_t ii = 0; ii < _system.getNtracks(); ii++ ){
//printf("EUTelDafFitter::dafEvent track %3d \n", ii);
    //run track fitte
    _nClusters++;
    _system.fitPlanesInfoDaf(_system.tracks.at(ii));
//printf("EUTelDafFitter::dafEvent track %3d info is OK \n", ii);
    //Check resids, intime, angles
    if(not checkTrack( _system.tracks.at(ii))) { continue;};
    int inTimeHits = checkInTime();
    if( inTimeHits < _nDutHits) { continue;}
//printf("EUTelDafFitter::dafEvent track %3d is OK \n", ii);
 
    dumpToAscii();
    //Fill plots
    if(_histogramSwitch){ 
      fillPlots( _system.tracks.at(ii) ); 
      fillDetailPlots( _system.tracks.at(ii) ); 
    }
    //Dump to LCIO
    if( _addToLCIO) { addToLCIO(_system.tracks.at(ii)); }
    _nTracks++;
  }

 //Add track collection
  if(_addToLCIO)
  { 
    event->addCollection(_fittrackvec,_trackCollectionName); 
    std::string  sfitpoints = "" ;  

    
    for(int i = 0; i<2000; i++) //TODO (Phillip Hamnett) Why is this hard coded to 1000? Ric changed to 2000.
    {
      sfitpoints = "fitpoints" + i;
      try
      {
       dynamic_cast < LCCollectionVec * > ( event->getCollection( sfitpoints ) )  ;
      }
      catch(...)
      {
        break;
      } 
    }
    event->addCollection(_fitpointvec, sfitpoints );
  }

//   if( _isFirstEvent )
//   {
//     _isFirstEvent = false;
//   }

  //cout << "DafFitter " << "58" << endl;
}

void EUTelDafFitter::addToLCIO(daffitter::TrackCandidate* track){
  //cout << "DafFitter " << "59" << endl;
  TrackImpl * fittrack = new TrackImpl();
  // Impact parameters are useless and set to 0
  fittrack->setD0(0.);        // impact paramter of the track in (r-phi)
  fittrack->setZ0(0.);        // impact paramter of the track in (r-z)
  fittrack->setTanLambda(0.); // dip angle of the track at reference point

  daffitter::TrackEstimate* est = track->estimates.at(0);

  //No good way of storing the track angles, so
  fittrack->setOmega( est->getXdz()); // Storing dxdz as Omega
  fittrack->setPhi( est->getYdz() );   // Storing dx/dy as phi

  fittrack->setChi2(track->chi2);
  fittrack->setNdf(int(track->ndof + 0.2f) );
  
//  fittrack->setIsReferencePointPCA(false);
  float refpoint[3];
  
  for(size_t plane = 0; plane < _system.planes.size(); plane++){
    daffitter::FitPlane& pl = _system.planes.at(plane);
    daffitter::TrackEstimate* estim = track->estimates.at( plane );
    TrackerHitImpl * fitpoint = new TrackerHitImpl();
    // Hit type 32 means a fitted hit
    fitpoint->setType(32);
    double pos[3];
    pos[0]= estim->getX() / 1000.0;
    pos[1]= estim->getY() / 1000.0;
    pos[2]= pl.getMeasZ() / 1000.0;
// overload z coordinate calculation -> important for proper sensor Identification by the hit coordinates based onthe refhit collection
    pos[2] = getZfromRefHit(plane, pos);    
//
    fitpoint->setPosition(pos);
    // Covariance matrix of the fitted position
    // (stored as lower triangle matrix, i.e.  cov(xx),cov(y,x),cov(y,y) ).
    float cov[TRKHITNCOVMATRIX];
    cov[0]= estim->cov(0,0);
    cov[1]= estim->cov(0,1);
    cov[2]= estim->cov(1,1);
    //Error of z position of fit is unclear to me, this would be a systematic alignment
    //error. Set to 0 along with all covariances.
    cov[3]=cov[4]=cov[5]=0.;
    fitpoint->setCovMatrix(cov);
    _fitpointvec->push_back(fitpoint);
    fittrack->addHit(fitpoint);

/* this code causes seg fault ....
    //At this point, we have a fitted hit position called fitpoint, here we will add the clustering information
    streamlog_out(DEBUG0) << "Track hit is at position " << fitpoint->getPosition()[0] << ", " << fitpoint->getPosition()[1] << ", " << fitpoint->getPosition()[2] << endl;
/*    streamlog_out(DEBUG0) << "_clusterVec contains: " << _clusterVec->getNumberOfElements() << endl;
    for(int i = 0; i < _clusterVec->getNumberOfElements(); ++i){
      TrackerPulse * cluster = dynamic_cast< TrackerPulse* >(_clusterVec->getElementAt(i));
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
      streamlog_out(DEBUG0) << "OMG WE FOUND IT, THE SENSOR ID IS..." << static_cast<int> ( ( cell0 & mask ) >> rhs ) << endl;
      streamlog_out(DEBUG0) << "clusterdata has CellID0: " << clusterdata->getCellID0() << endl;
      streamlog_out(DEBUG0) << "clusterdata has CellID1: " << clusterdata->getCellID1() << endl;
      streamlog_out(DEBUG0) << "clusterdata has time: " << clusterdata->getTime() << endl;
      std::vector< float > chargevalues = clusterdata->getChargeValues();

      double clusterseedx(0);
      double clusterseedy(0);
      double ClusterHitCut = 0.0184*20;
      streamlog_out(DEBUG0) << "chargevalues has size: " << chargevalues.size() << endl;
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
      if(static_cast<unsigned int> ( ( cell0 & mask ) >> rhs ) == plane){
	double trackpositionx = fitpoint->getPosition()[0];
	double trackpositiony = fitpoint->getPosition()[1];
        if(abs(clusterseedy - trackpositiony) < ClusterHitCut &&
	   abs(clusterseedx - trackpositionx) < ClusterHitCut){
          int ClusterSize = static_cast< int >(chargevalues.size()/3.0);
          if(ClusterSize > MAXCLUSTERSIZE){
            MAXCLUSTERSIZE = ClusterSize;
          }
          double Chi2OverNdof = fittrack->getChi2()/fittrack->getNdf();
          streamlog_out(DEBUG0) << "Beginning to fill the first clustering histograms" << endl;
          _aidaHistoMap["ClusterSize"]->fill(static_cast< double >(ClusterSize)); 
          _aidaHistoMap2D["ClusterSizeVsXPosition"]->fill(trackpositionx,static_cast< double >(ClusterSize)); 
          _aidaHistoMap2D["ClusterSizeVsYPosition"]->fill(trackpositiony,static_cast< double >(ClusterSize)); 
          streamlog_out(DEBUG0) << "Filled all histograms except Chi2" << endl;
          streamlog_out(DEBUG0) << "ClusterSize is " << ClusterSize << endl;
          streamlog_out(DEBUG0) << "Chi2OverNdof is " << Chi2OverNdof << endl;
          _aidaHistoMap2D["ClusterSizeVsChi2"]->fill(Chi2OverNdof,static_cast< double >(ClusterSize)); 
          streamlog_out(DEBUG0) << "Filled the first clustering histograms" << endl;



//  fittrack->setChi2(track->chi2);
//  fittrack->setNdf(int(track->ndof + 0.2f) );
 

          double newminbinx(minx);
          double newminbiny(miny);
          int actualbinx(0);
          int actualbiny(0);
          bool notfoundx = true;
          bool notfoundy = true;

          //int acutalbin(0);

          std::pair< int, int > position;
          streamlog_out(DEBUG0) << "Starting the crazy while loops" << endl;
          while(notfoundx){
            if(trackpositionx > maxx){
              string xistoobig = "The x coordinate of the hit is larger than the size of the sensor. Skipping this track";
              throw xistoobig;
            } //end of: if(trackpositionx > maxx)
            else if(trackpositionx < newminbinx){
              string xistoosmall = "The x coordinate of the hit is either smaller than the size of the sensor, or should have been caught in the previous bin. Skipping this track";
              throw xistoosmall;
            } //end of: else if(x < newminbinx)
            else if(trackpositionx < newminbinx + binsizex){
              notfoundx = false;
        
              while(notfoundy){
                if(trackpositiony > maxy){
                  string yistoobig = "The y coordinate of the hit is larger than the size of the sensor. Skipping this track";
                  throw yistoobig;
                } //end of: if(y > mayy)
                else if(trackpositiony < newminbiny){
                  string yistoosmall = "The y coordinate of the hit is either smaller than the size of the sensor, or should have been caught in the previous bin. Skipping this track";
                  throw yistoosmall;
                } //end of: else if(y < newminbiny)
                else if(trackpositiony < newminbiny + binsizey){
                  notfoundy = false;
                  _xPositionForClustering[actualbinx].push_back(ClusterSize);
                  _yPositionForClustering[actualbiny].push_back(ClusterSize);

                } //end of: else if(y < newminbiny + binsizey)
                else{
                  newminbiny += binsizey;
                  actualbiny++;
                } //end of: else [else if(y < newminbiny + binsizey)]
              } // end of: while(notfoundy)
            } //end of: else if(x < newminbinx + binsizex)
            else{
              newminbinx += binsizex;
              actualbinx++;
            } //end of: else [else if(x < newminbinx + binsizex)]
          } // end of: while(notfoundx)
          streamlog_out(DEBUG0) << "Ending the crazy while loops" << endl;
        
          _Chi2sForAverage[ClusterSize].push_back(Chi2OverNdof);

//Plotting the residuals for fitted tracks 
          for( size_t ii = 0; ii < _system.planes.size() ; ii++){
            daffitter::FitPlane& plane = _system.planes.at(ii);
            char iden[4];
            sprintf(iden, "%d", plane.getSensorID());
            string bname = static_cast< string >("pl") + iden + "_";
            daffitter::TrackEstimate* estim = track->estimates.at(ii);
            for(size_t w = 0; w < plane.meas.size(); w++){
              if( plane.weights(w) < 0.5f ) {  continue; }
              daffitter::Measurement& meas = plane.meas.at(w);
              _aidaHistoMap2D["ResidualXVsClusterSize"]->fill(ClusterSize,(estim->getX() - meas.getX())*1e-3);
              _aidaHistoMap2D["ResidualYVsClusterSize"]->fill(ClusterSize,(estim->getY() - meas.getY())*1e-3);
              _resolutionXForClustering[ClusterSize].push_back((estim->getX() - meas.getX())*1e-3);
              _resolutionYForClustering[ClusterSize].push_back((estim->getY() - meas.getY())*1e-3);
            }
          }

	}
      }
      streamlog_out(DEBUG0) << endl;
    }
*/
    if(plane == 0){
      refpoint[0] = pos[0];
      refpoint[1] = pos[1];
      refpoint[2] = pos[2];
    }
    //Also add measurement point
    for(size_t mm = 0; mm < pl.meas.size(); mm++){
      if( pl.weights(mm) < 0.5f){ continue; }
      if( pl.isExcluded()) { continue; }
      TrackerHitImpl* meashit = static_cast<TrackerHitImpl*> ( _hitCollection->getElementAt( pl.meas.at(mm).getIden()) );
      fittrack->addHit(meashit);
    }
  }
  fittrack->setReferencePoint(refpoint);
  _fittrackvec->addElement(fittrack);
  //cout << "DafFitter " << "60" << endl;
}

double EUTelDafFitter::getZfromRefHit(int plane, double *pos){
  //cout << "DafFitter " << "61" << endl;
//rintf("plane %5d estim X:%5.2f Y:%5.2f Z:%5.2f \n", plane,  pos[0], pos[1], pos[2] );      

 try
 {
   double doesNothing = pos[2];
   doesNothing++;  //Stops a warning about unused variable
 }
 catch(...)
 {
   printf(" input array pos{3} access error, perhaps it has less then 3 elements => check yoru input !! \n"); 
   return 0.;
 }
   
  if( ReferenceHitVecIsSet() )
  {
    streamlog_out( MESSAGE5 ) << "_referenceHitVec is empty" << endl;
    return 0;
  }
  
  EUTelReferenceHit* refhit = static_cast< EUTelReferenceHit*> ( _referenceHitVec->getElementAt(plane) ) ;
//        printf(" _referenceHitVec %p refhit %p \n", _referenceHitVec, refhit);
        
  TVector3 lpoint( pos[0], pos[1], pos[2] );
  TVector3 lvector( 0., 0., 1. );
  TVector3 hitInPlane( refhit->getXOffset(), refhit->getYOffset(), refhit->getZOffset());
  TVector3 norm2Plane( refhit->getAlpha(), refhit->getBeta(), refhit->getGamma() );

  TVector3 point( 1.,1.,1. );
          
  double linecoord_numenator   = norm2Plane.Dot(hitInPlane-lpoint);
  double linecoord_denumenator = norm2Plane.Dot(lvector);
  point = (linecoord_numenator/linecoord_denumenator)*lvector + lpoint;

//  printf("plane %5d estim X:%5.2f Y:%5.2f Z:%5.2f -> (%5.2f, %5.2f, %5.2f)\n", plane,  pos[0], pos[1], pos[2], point(0), point(1), point(2) );      

  //cout << "DafFitter " << "62" << endl;
return point(2);
}

void EUTelDafFitter::dafEnd() {
  //cout << "DafFitter " << "63" << endl;


}
#endif // USE_GEAR
