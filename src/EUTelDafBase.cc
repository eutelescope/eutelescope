
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
#include "EUTelDafBase.h"
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

// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <IO/LCWriter.h>
#include <UTIL/LCTime.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
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

#include <TVector3.h>
#include <Eigen/Geometry> 

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;


EUTelDafBase::EUTelDafBase(std::string name) : marlin::Processor(name) {
  //Universal DAF params
  
  // input collection
  std::vector<std::string > HitCollectionNameVecExample;
  HitCollectionNameVecExample.push_back("hit");
  registerInputCollections(LCIO::TRACKERHIT,"HitCollectionName", "Names of input hit collections", _hitCollectionName,HitCollectionNameVecExample);
  
  //Tracker system options
  registerOptionalParameter("MakePlots", "Should plots be made and filled?", _histogramSwitch, static_cast<bool>(false));
  registerProcessorParameter("TelescopePlanes","List of sensor IDs for the telescope planes. These planes are used for the track finder, and track fitter.", _telPlanes ,std::vector<int>());
  registerOptionalParameter("DutPlanes","List of sensor IDs for the DUT planes. Used to make the decision on whether ro accept the track or not. These planes are not used in track finder, and not in the track fitter unless option 'useDutsInFit' is set.", _dutPlanes ,std::vector<int>());
  registerOptionalParameter("Ebeam", "Beam energy [GeV], used to calculate amount of scatter", _eBeam,  static_cast < float > (120.0));
  registerOptionalParameter("TelResolutionX", "Sigma of telescope resolution in the global X plane,", _telResX,  static_cast < float > (5.3));
  registerOptionalParameter("TelResolutionY", "Sigma of telescope resolution in the global Y plane,", _telResY,  static_cast < float > (5.3));
  registerOptionalParameter("DutResolutionX", "Sigma of telescope resolution in the global X plane,", _dutResX,  static_cast < float > (115.4));
  registerOptionalParameter("DutResolutionY", "Sigma of telescope resolution in the global Y plane,", _dutResY,  static_cast < float > (14.4));
  registerOptionalParameter("ScaleScattering","Scale thickness of DUT planes", _scaleScatter, static_cast<float>(1.0f)); 
  registerOptionalParameter("RadiationLengths","Radiation lengths of planes.", _radLength, std::vector<float>());
  
  //Track finder options
  registerOptionalParameter("FinderRadius","Track finding: The maximum allowed distance between to hits in the xy plane for inclusion in track candidate", _clusterRadius, static_cast<float>(300.0));
  registerOptionalParameter("Chi2Cutoff","DAF fitter: The cutoff value for a measurement to be included in the fit.", _chi2cutoff, static_cast<float>(300.0f));
  registerOptionalParameter("RequireNTelPlanes","How many telescope planes do we require to be included in the fit?",_nSkipMax ,static_cast <float> (0.0f));
  registerOptionalParameter("NominalDxdz", "dx/dz assumed by track finder", _nXdz, static_cast<float>(0.0f));
  registerOptionalParameter("NominalDydz", "dy/dz assumed by track finder", _nYdz, static_cast<float>(0.0f));
  
  // 
  registerOptionalParameter("ReferenceCollection","reference hit collection name ", _referenceHitCollectionName, static_cast <string> ("reference") );
  registerOptionalParameter("ApplyToReferenceCollection","Do you want the reference hit collection to be corrected by the shifts and tilts from the alignment collection? (default - false )",  _applyToReferenceHitCollection, static_cast< bool   > ( false ));
 

  //Track quality parameters
  registerOptionalParameter("MaxChi2OverNdof", "Maximum allowed global chi2/ndof", _maxChi2, static_cast<float> ( 9999.0));
  registerOptionalParameter("TrackAsciiName", "Filename for fitted tracks", _asciiName, string ("tracks.txt"));
  registerOptionalParameter("NDutHits", "How many DUT hits do we need in order to accept track?", _nDutHits, static_cast <int>(0));
  registerOptionalParameter("AlignmentCollectionNames", "Names of alignment collections, should be in same order as application", _alignColNames, std::vector<std::string>());
}

bool EUTelDafBase::defineSystemFromData(){
  bool gotIt = true;
  bool gotPlane = false;
  for(size_t plane = 0; plane < _system.planes.size(); plane++){
    daffitter::FitPlane& pl = _system.planes.at(plane);
    gotPlane = false;
    for(size_t meas = 0; meas < pl.meas.size(); meas++){
//printf("plane %5d meas %5d ref %5d \n", plane, meas, _nRef.at(plane) );
      if( _nRef.at(plane) > 2){ gotPlane = true; continue; }
      if( _nRef.at(plane) == 0 ){
//printf("pl.getRef0: %8.3f %8.3f %8.3f \n", pl.getRef0()(0), pl.getRef0()(1), pl.getRef0()(2) );
	pl.setRef0( Vector3f(pl.meas.at(meas).getX(), pl.meas.at(meas).getY(), pl.meas.at(meas).getZ()));
//printf("pl.getRef0: %8.3f %8.3f %8.3f \n", pl.getRef0()(0), pl.getRef0()(1), pl.getRef0()(2) );
	_nRef.at(plane)++;
	gotPlane = false;
//	printf("FALSE plane %5d meas %5d ref %5d \n", plane, meas, _nRef.at(plane) );
 continue;
      }
      if( fabs(pl.meas.at(meas).getX() - pl.getRef0()(0) ) < 500) { continue; }
      if( fabs(pl.meas.at(meas).getY() - pl.getRef0()(1) ) < 500) { continue; }
      if( _nRef.at(plane) == 1 ){
	pl.setRef1( Vector3f(pl.meas.at(meas).getX(), pl.meas.at(meas).getY(), pl.meas.at(meas).getZ()));
	_nRef.at(plane)++;
	gotPlane = false;
//	printf("FALSE plane %5d meas %5d ref %5d \n", plane, meas, _nRef.at(plane) );
 continue;
      }
      if( fabs(pl.meas.at(meas).getX() - pl.getRef1()(0) ) < 500) {  continue; }
      if( fabs(pl.meas.at(meas).getY() - pl.getRef1()(1) ) < 500) {  continue; }
      if( _nRef.at(plane) == 2 ){
	pl.setRef2( Vector3f(pl.meas.at(meas).getX(), pl.meas.at(meas).getY(), pl.meas.at(meas).getZ()));
	_nRef.at(plane)++;
	getPlaneNorm(pl);
	pl.print();
	gotPlane = true;
//	printf("TRUE plane %5d meas %5d ref %5d \n", plane, meas, _nRef.at(plane) );
 continue;
      }
    }
    if(not gotPlane) { gotIt = false;}
//   printf("gotPlane  %5d  gotIt %5d \n",gotPlane, gotIt);

}
  return(gotIt);
}

void EUTelDafBase::gearRotate(size_t index, size_t gearIndex){
  daffitter::FitPlane& pl = _system.planes.at(index);

  double gRotation[3] = { 0., 0., 0.};
  gRotation[0] = _siPlanesLayerLayout->getLayerRotationXY(gearIndex); // Euler alpha ;
  gRotation[1] = _siPlanesLayerLayout->getLayerRotationZX(gearIndex); // Euler alpha ;
  gRotation[2] = _siPlanesLayerLayout->getLayerRotationZY(gearIndex); // Euler alpha ;
  // transform into radians
  gRotation[0] =  gRotation[0]*3.1415926/180.; 
  gRotation[1] =  gRotation[1]*3.1415926/180.; 
  gRotation[2] =  gRotation[2]*3.1415926/180.; 
  
  //Reference points define plane
  //Transform ref, cloning hitmaker logic
  TVector3 ref0 = TVector3(0.0, 0.0, 0.0);
  TVector3 ref1 = TVector3(10.0, 0.0, 0.0);
  TVector3 ref2 = TVector3(0.0, 10.0, 0.0);

  double zZero        = _siPlanesLayerLayout->getSensitivePositionZ(gearIndex); // mm
  double zThickness   = _siPlanesLayerLayout->getLayerThickness(gearIndex); // mm

  double nomZ = zZero + 0.5 * zThickness;

  if( TMath::Abs( gRotation[0] ) > 1e-6 ){
    ref0.RotateZ( gRotation[0] );
    ref1.RotateZ( gRotation[0] );
    ref2.RotateZ( gRotation[0] );
  }
  if( TMath::Abs( gRotation[1] )> 1e-6 ) {
    ref0.RotateY( gRotation[1] );
    ref1.RotateY( gRotation[1] );
    ref2.RotateY( gRotation[1] );
  }
  if( TMath::Abs( gRotation[2] ) > 1e-6 ){
    ref0.RotateX( gRotation[2] );
    ref1.RotateX( gRotation[2] );
    ref2.RotateX( gRotation[2] );
  }


  pl.setRef0( Vector3f( ref0.X() * 1000.0f, ref0.Y() * 1000.0f, (ref0.Z() + nomZ) * 1000.0f ));
  pl.setRef1( Vector3f( ref1.X() * 1000.0f, ref1.Y() * 1000.0f, (ref1.Z() + nomZ) * 1000.0f ));
  pl.setRef2( Vector3f( ref2.X() * 1000.0f, ref2.Y() * 1000.0f, (ref2.Z() + nomZ) * 1000.0f ));
//printf("gearRotate[%2d] %8.3f %8.3f %8.3f \n", index, ref0.X(), ref0.Y(), ref0.Z()+nomZ);
//printf("gearRotate[%2d] %8.3f %8.3f %8.3f \n", index, ref1.X(), ref1.Y(), ref1.Z()+nomZ);
//printf("gearRotate[%2d] %8.3f %8.3f %8.3f \n", index, ref2.X(), ref2.Y(), ref2.Z()+nomZ);

  //Tracks are propagated to glob xy plane => Errors are in glob xy plane. scales like cosine
  //Errors not corrected for xy rotation
  pl.scaleErrors( std::fabs(std::cos(gRotation[1])), std::fabs(std::cos(gRotation[0])));
  getPlaneNorm(pl);
} 
/*
Vector3f EUTelDafBase::applyAlignment(EUTelAlignmentConstant* alignment, Vector3f point){
  Vector3f outpoint;
  outpoint(0) = point(0) * (1.0 + alignment->getAlpha())  + alignment->getGamma() * point(1);
  outpoint(1) = point(1) * (1.0 + alignment->getBeta())   - alignment->getGamma() * point(0);
  outpoint(2) = point(2);
  // Shifts
  outpoint(0) -= 1000.0f * alignment->getXOffset();
  outpoint(1) -= 1000.0f * alignment->getYOffset();
  outpoint(2) -= 1000.0f * alignment->getZOffset();

//printf("point %5.3f %5.3f %5.3f \n", point(0), point(1), point(2) );
//printf("outut %5.3f %5.3f %5.3f \n", outpoint(0), outpoint(1), outpoint(2) );

  return(outpoint);
}*/

Vector3f EUTelDafBase::applyAlignment(EUTelAlignmentConstant* alignment, Vector3f point){
  Vector3f outpoint;
  double alpha = alignment->getAlpha();
  double beta  = alignment->getBeta();  
  double gamma = alignment->getGamma();
  double z_sensor = point(2);

// sync with updated sign convention of the rotation angles:
//
  outpoint(0) =                point(0) + (-1)*gamma * point(1) +      beta  * (point(2) - z_sensor) ;
  outpoint(1) =        gamma * point(0) +              point(1) + (-1)*alpha * (point(2) - z_sensor) ;
  outpoint(2) = (-1) * beta  * point(0) +      alpha * point(1) +              (point(2) - z_sensor) ;

  // second the shift
  outpoint(0) -= 1000.0f * alignment->getXOffset();
  outpoint(1) -= 1000.0f * alignment->getYOffset();
  outpoint(2) -= 1000.0f * alignment->getZOffset();

            
  outpoint(2) += z_sensor ;
//printf("input %5.3f %5.3f %5.3f   %5.3f %5.3f %5.3f \n", point(0), point(1), point(2),  alpha, beta, gamma );
//printf("outut %5.3f %5.3f %5.3f \n", outpoint(0), outpoint(1), outpoint(2) );

  return(outpoint);
}

void EUTelDafBase::alignRotate(std::string collectionName, LCEvent* event) {
  cout << "Reading in alignment collection " << collectionName << endl;
  LCCollectionVec * alignmentCollectionVec;
  try {
    alignmentCollectionVec     = dynamic_cast < LCCollectionVec * > (event->getCollection(collectionName));
  } catch (DataNotAvailableException& e) {
    throw runtime_error("Unable to open alignment collection " + collectionName);
  }
  for( size_t plane = 0; plane < _system.planes.size() ; plane++){
// printf("plane %5d \n", plane);
    daffitter::FitPlane& pl = _system.planes.at(plane);
    int iden = pl.getSensorID();
    for ( size_t ii = 0; ii < alignmentCollectionVec->size(); ++ii ) {
      try{
      EUTelAlignmentConstant * alignment = static_cast< EUTelAlignmentConstant * >
	( alignmentCollectionVec->getElementAt( ii ) );
      if( alignment->getSensorID() != iden) { continue; }
      pl.setRef0( applyAlignment(alignment, pl.getRef0()) );
      pl.setRef1( applyAlignment(alignment, pl.getRef1()) );
      pl.setRef2( applyAlignment(alignment, pl.getRef2()) );
      //Errors not corrected for xy rotation
      getPlaneNorm(pl);
      pl.scaleErrors(alignment->getAlpha() + 1.0f, alignment->getBeta() + 1.0f);
      pl.print();
      }
      catch(...)
      {
        printf("could not find sensor in %s at %5d \n", collectionName.c_str(), ii);
      }
    }
  }
}
void EUTelDafBase::getPlaneNorm(daffitter::FitPlane& pl){
  Vector3f l1 = pl.getRef1() - pl.getRef0();
  Vector3f l2 = pl.getRef2() - pl.getRef0();
  //Calculate plane normal vector from ref points
  pl.setPlaneNorm( l2.cross(l1));
  // cout << "Ref0 " << endl << pl.getRef0() << endl;
  // cout << "Norm vec " << endl << pl.getPlaneNorm() << endl;
  // cout << "Norm dot r0 - r1 " << pl.getPlaneNorm().dot(pl .getRef1() - pl.getRef0()) << endl;
  // cout << "Norm dot r2 - r0 " << pl.getPlaneNorm().dot(pl.getRef2() - pl.getRef0()) << endl;
}

void EUTelDafBase::init() {
  //cout << "Trying to open " << _asciiName << endl;
  trackstream.open(_asciiName.c_str());
  
  printParameters ();

  _iRun = 0; _iEvt = 0; _nTracks = 0; _nClusters =0;
  n_passedNdof =0; n_passedChi2OverNdof = 0; n_passedIsnan = 0;
  _initializedSystem = false;

  //Geometry description
  _siPlanesParameters  = const_cast<gear::SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

  //Use map to sort planes by z
  _zSort.clear();
  for(int plane = 0; plane < _siPlanesLayerLayout->getNLayers(); plane++){
    _zSort[ _siPlanesLayerLayout->getLayerPositionZ(plane) * 1000.0 ] = plane;
  }

  //Add Planes to tracker system,
  map<float, int>::iterator zit = _zSort.begin();
  size_t index(0), nActive(0);
  for( ; zit != _zSort.end(); index++, zit++){
//printf("index= %5d \n",index);
    _nRef.push_back(3);
    int sensorID = _siPlanesLayerLayout->getID( (*zit).second );
    //float zPos  = (*zit).first + _siPlanesLayerLayout->getSensitiveThickness( (*zit).second );
    //Read sensitive as 0, in case the two are different
    float zPos  = _siPlanesLayerLayout->getSensitivePositionZ( (*zit).second )* 1000.0
      + 0.5 * 1000.0 *  _siPlanesLayerLayout->getSensitiveThickness( (*zit).second) ; // um
    //Figure out what kind of plane we are dealing with
    float errX(0.0f), errY(0.0f);
    bool excluded = true;

    // Get scatter using x / x0
    float radLength = _siPlanesLayerLayout->getLayerThickness( (*zit).second ) /  _siPlanesLayerLayout->getLayerRadLength( (*zit).second );
    if( _radLength.size() > index){
      radLength = _radLength.at(index);
    }
    float scatter = getScatterThetaVar( radLength );
    
    //Is current plane a telescope plane?
    if( find(_telPlanes.begin(), _telPlanes.end(), sensorID) != _telPlanes.end()){
      _nRef.at(index) = 0;
      errX = _telResX; errY = _telResY;
      excluded = false;
    }
    //Is current plane a DUT plane?
    if( find(_dutPlanes.begin(), _dutPlanes.end(), sensorID) != _dutPlanes.end()){
      _nRef.at(index) = 0;
      errX = _dutResX; errY = _dutResY;
      radLength = _scaleScatter * radLength;
      scatter = getScatterThetaVar( radLength );
      //scatter *= _scaleScatter;
    }
    //If plane is neither Tel nor Dut, all we need is the zpos and scatter amount.

    //Add plane to tracker system
    if(not excluded){ nActive++;}
// printf("sensorID %5d zpos %8.3f %5d \n", sensorID, zPos, _nRef.at(index) );
   _system.addPlane(sensorID, zPos , errX, errY, scatter, excluded);
    gearRotate(index, (*zit).second);
  }
 

  //Prepare track finder
  _system.setClusterRadius(_clusterRadius);
  _system.setNominalXdz(_nXdz);
  _system.setNominalYdz(_nYdz);
  //Prepare and preallocate memory for track fitter
  _system.setChi2OverNdofCut(_maxChi2);
  _system.setDAFChi2Cut(_chi2cutoff);
  _system.init();

  //Fuzzy assignment by DAF might make a plane only partially included, This means ndof is
  //not a integer. Everything above (ndof - 0.5) is assumed to include at least ndof degrees of
  //freedom. 
  _ndofMin = -4 + _nSkipMax * 2 - 0.5;
  streamlog_out ( MESSAGE ) << "NDOF min is " << _ndofMin << endl;
  if( _ndofMin < 0.5) {
    streamlog_out ( ERROR ) << "Too few active planes(" << nActive << ") when " << _nSkipMax << " planes can be skipped." 
			    << "Please check your configuration." << endl;
    exit(1);
  }
  dafInit();
    
  if(_histogramSwitch) {
    bookHistos(); 
    bookDetailedHistos();
  }

  //Define region for edge masking
  for(size_t ii = 0; ii < _dutPlanes.size(); ii++){
    int iden = _dutPlanes.at(ii);
    int xMin = _colMin.size() > ii ? _colMin.at(ii) : -9999999;
    int xMax = _colMax.size() > ii ? _colMax.at(ii) : 9999999;
    int yMin = _rowMin.size() > ii ? _rowMin.at(ii) : -9999999;
    int yMax = _rowMax.size() > ii ? _rowMax.at(ii) : 9999999;
    _colMinMax[iden] = make_pair(xMin, xMax);
    _rowMinMax[iden] = make_pair(yMin, yMax);
  }
}

void EUTelDafBase::processRunHeader (LCRunHeader * rdr) {
  auto_ptr<EUTelRunHeaderImpl> header ( new EUTelRunHeaderImpl (rdr) );
  header->addProcessor( type() ) ;
  ++_iRun;
}


float EUTelDafBase::getScatterThetaVar( float radLength ){
  //From pdg live
  float scatterTheta = 0.0136f/_eBeam * sqrt( radLength ) *  (1.0f + 0.038f * std::log(radLength) );
  return(scatterTheta * scatterTheta);
}

size_t EUTelDafBase::getPlaneIndex(float zPos){
  //Get plane index from z-position of hit
  map<float,int>::iterator it = _zSort.begin();
  size_t index(0);
  bool foundIt(false);
  for(;it != _zSort.end(); index++, it++){
    if( fabs((*it).first - zPos) < 30000.0){ 
      foundIt = true; break;
    }
  }
  if(not foundIt){ 
    streamlog_out (ERROR ) << "Found hit at z=" << zPos << " , not able to assign to any plane!" << endl; 
    return(-1);
  }else{
//    streamlog_out (MESSAGE ) << "Found hit at z=" << zPos << " , assign index = " << index << endl;     
  }

  return(index);
}

int EUTelDafBase::guessSensorID(const double * hit ) 
{

  int sensorID = -1;
  double minDistance =  numeric_limits< double >::max() ;
//  double * hitPosition = const_cast<double * > (hit->getPosition());

//  LCCollectionVec * referenceHitVec     = dynamic_cast < LCCollectionVec * > (evt->getCollection( _referenceHitCollectionName));
  if( ReferenceHitVecIsSet() )
  {
    streamlog_out(MESSAGE) << "_referenceHitVec is empty" << endl;
    return 0;
  }

      for(size_t ii = 0 ; ii <  _referenceHitVec->getNumberOfElements(); ii++)
      {
        EUTelReferenceHit* refhit = static_cast< EUTelReferenceHit*> ( _referenceHitVec->getElementAt(ii) ) ;
//        printf(" _referenceHitVec %p refhit %p \n", _referenceHitVec, refhit);
        
        TVector3 hit3d( hit[0], hit[1], hit[2] );
        TVector3 hitInPlane( refhit->getXOffset(), refhit->getYOffset(), refhit->getZOffset());
        TVector3 norm2Plane( refhit->getAlpha(), refhit->getBeta(), refhit->getGamma() );
 
        double distance = abs( norm2Plane.Dot(hit3d-hitInPlane) );
//        printf("iPlane %5d   hitPos:  [%8.3f;%8.3f%8.3f]  distance: %8.3f \n", refhit->getSensorID(), hitPosition[0],hitPosition[1],hitPosition[2], distance  );
        if ( distance < minDistance ) 
        {
           minDistance = distance;
//           sensorID = refhit->getSensorID(); // proper ID
           sensorID = ii;                    // number in the GEAR file z ordered
//           printf("sensorID: %5d \n", sensorID );
        }    

      }

  return sensorID;
}


void EUTelDafBase::readHitCollection(LCEvent* event){
  //Dump LCIO hit collection to tracker system
  //Extract hits from collection, add to tracker system
//printf("EUTelDafBase::readHitCollection \n");
  for(size_t i =0;i < _hitCollectionName.size();i++){
    try {
      _hitCollection = event->getCollection(_hitCollectionName[i]);
    } catch (DataNotAvailableException& e) {
      streamlog_out ( WARNING2 ) << "No input collection " << _hitCollectionName[i] << " found for event " << event->getEventNumber()
				 << " in run " << event->getRunNumber() << endl;
      throw SkipEventException( this );
    }
// printf("EUTelDafBase::readHitCollection add all hits from collection %5d \n",i);
    //Add all hits in collection to corresponding plane
    for ( int iHit = 0; iHit < _hitCollection->getNumberOfElements(); iHit++ ) {
      TrackerHitImpl* hit = static_cast<TrackerHitImpl*> ( _hitCollection->getElementAt(iHit) );
      const double * pos = hit->getPosition();
//      int planeIndex = getPlaneIndex( pos[2]  * 1000.0f);
      int planeIndex = guessSensorID( pos );
//      printf("planeIndex %5d of %5d \n", planeIndex, _system.planes.size() ); 
      bool region = checkClusterRegion( hit, _system.planes.at(planeIndex).getSensorID() );
      //if( not region){ continue; }
      if(planeIndex >=0 ) { 
	_system.addMeasurement( planeIndex, (float) pos[0] * 1000.0f, (float) pos[1] * 1000.0f, (float) pos[2] * 1000.0f,  region, iHit);
      }
    }
// printf("EUTelDafBase::readHitCollection add all hits from collection %5d DONE \n",i);
 
  }
// printf("EUTelDafBase::readHitCollection -- DONE \n");
 
}

bool EUTelDafBase::checkClusterRegion(lcio::TrackerHitImpl* hit, int iden){
  bool goodRegion(true);
  if( hit->getType() == kEUTelAPIXClusterImpl ){
    auto_ptr<EUTelVirtualCluster> cluster( new EUTelSparseClusterImpl< EUTelAPIXSparsePixel >
  					   ( static_cast<TrackerDataImpl *> ( hit->getRawHits()[0] )));
    // float xPos(0), yPos(0);
    // cluster->getCenterOfGravity(xPos, yPos);
    int xSeed(0), ySeed(0);
    cluster->getCenterCoord(xSeed, ySeed);
    int xSize(0), ySize(0);
    cluster->getClusterSize(xSize, ySize);
    //if( iden == 10) std::cout << iden << " c: " << xSeed << " r: " << ySeed << std::endl;
    std::pair<int, int> &colMinMax = _colMinMax[iden]; 
    if( (xSeed - xSize / 2 ) < colMinMax.first ) { goodRegion = false;}
    if( (xSeed + xSeed / 2) > colMinMax.second ) { goodRegion = false;}
    std::pair<int, int> &rowMinMax = _rowMinMax[iden]; 
    if( (ySeed - ySize / 2) < rowMinMax.first ) { goodRegion = false;}
    if( (ySeed + ySize / 2) > rowMinMax.second ) { goodRegion = false;}
  }
  return(goodRegion);
}

int EUTelDafBase::checkInTime(daffitter::TrackCandidate * track){
  size_t nMatches(0);
  for( size_t ii = 0; ii < _system.planes.size() ; ii++){
    daffitter::FitPlane& plane = _system.planes.at(ii);
    int sensorID = plane.getSensorID();
//    printf("plane %5d with ID %5d \n", ii, sensorID);
    //Check if any DUT plane is "In time" with the 
    if( find(_dutPlanes.begin(), _dutPlanes.end(), sensorID) == _dutPlanes.end()){ continue; }
    //In timeness can be checked by seeing if the plane has assigned weight
    for(size_t w = 0; w < plane.meas.size(); w++){
//      printf("plane %5d with ID %5d --- weight = %5.2f\n", ii, sensorID, plane.weights(w) );
      if( plane.weights(w) < 0.5f ) {  continue; }
      nMatches++; break;
    }
  }
  return(nMatches);
}

void EUTelDafBase::processEvent(LCEvent * event){
  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);
  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG2 ) << "EORE found: nothing else to do." << endl;
    return;
  }
  if( isFirstEvent() ){
    for(size_t ii = 0; ii < _alignColNames.size(); ii++){
//printf("aligncollection %5d \n", ii);
      alignRotate(_alignColNames.at(ii), event);
    }
    if ( _applyToReferenceHitCollection ) 
    {
      _referenceHitVec = dynamic_cast < LCCollectionVec * > (event->getCollection( _referenceHitCollectionName));
    }
  }
  //Prepare tracker system for new data, new tracks
  _system.clear();

  //Dump hit collection to collection sorted by plane
  readHitCollection(event);

  if( not _initializedSystem ){ 
//    _initializedSystem = defineSystemFromData();
//    if(not _initializedSystem) { return; }
//    cout << "Initialized system at event " << event->getEventNumber() << endl;
  }
//printf("EUTelDafBase::processEvent \n");
  //Run track finder
  _system.clusterTracker();
  
  //Child specific actions
  dafEvent(event);

  if(event->getEventNumber() % 1000 == 0){
    streamlog_out ( MESSAGE ) << "Accepted " << _nTracks <<" tracks at event " << event->getEventNumber() << endl;
  }
}

bool EUTelDafBase::checkTrack(daffitter::TrackCandidate * track){

//printf("EUTelDafBase::checkTrack %7.3f %8.3f \n", track->ndof, track->chi2 );
//track->print();
  if( track->ndof < _ndofMin) { return(false); }
  n_passedNdof++;
  if( (track->chi2 / track->ndof) > _maxChi2 ) { return(false); }
  n_passedChi2OverNdof++;
  if( isnan(track->ndof)) { return(false); }
  n_passedIsnan++;

//printf("EUTelDafBase::checkTrack %7.3f %8.3f ;; return true;\n", track->ndof, track->chi2 );
  return(true);
}

void EUTelDafBase::dumpToAscii(daffitter::TrackCandidate *track){
  trackstream << _nTracks << endl;
  for( size_t ii = 0; ii < _system.planes.size() ; ii++){
    daffitter::FitPlane& plane = _system.planes.at(ii);
    for(size_t w = 0; w < plane.meas.size(); w++){
      if( plane.weights(w) < 0.5f ) {  continue; }
      daffitter::Measurement& meas = plane.meas.at(w);
      trackstream << meas.getX() << "\t" << meas.getY() << "\t"
		  << plane.getMeasZ() << endl;
    }
  }
  trackstream << endl;
}

void EUTelDafBase::fillPlots(daffitter::TrackCandidate *track){
  _aidaHistoMap["chi2"]->fill( track->chi2);
  _aidaHistoMap["logchi2"]->fill( std::log10(track->chi2));
  _aidaHistoMap["ndof"]->fill( track->ndof);
  _aidaHistoMap["chi2overndof"]->fill( track->chi2 / track->ndof);
  //Fill plots per plane
  for( size_t ii = 0; ii < _system.planes.size() ; ii++){
    daffitter::FitPlane& plane = _system.planes.at(ii);
    char iden[4];
    sprintf(iden, "%d", plane.getSensorID());
    string bname = (string)"pl" + iden + "_";
    //Plot resids, angles for all hits with > 50% includion in track.
    //This should be one measurement per track
    daffitter::TrackEstimate* estim = track->estimates.at(ii);
    for(size_t w = 0; w < plane.meas.size(); w++){
      if( plane.weights(w) < 0.5f ) {  continue; }
      daffitter::Measurement& meas = plane.meas.at(w);
      //Resids 
      _aidaHistoMap[bname + "residualX"]->fill( estim->getX() - meas.getX() );
      _aidaHistoMap[bname + "residualY"]->fill( estim->getY() - meas.getY() );
      //Resids 
      _aidaHistoMapProf1D[bname + "residualdXvsX"]->fill(estim->getX(), estim->getX() - meas.getX() );
      _aidaHistoMapProf1D[bname + "residualdYvsX"]->fill(estim->getX(), estim->getY() - meas.getY() );
      _aidaHistoMapProf1D[bname + "residualdXvsY"]->fill(estim->getY(), estim->getX() - meas.getX() );
      _aidaHistoMapProf1D[bname + "residualdYvsY"]->fill(estim->getY(), estim->getY() - meas.getY() );
      _aidaHistoMapProf1D[bname + "residualdZvsX"]->fill(estim->getX(), plane.getMeasZ() - meas.getZ()  );
      _aidaHistoMapProf1D[bname + "residualdZvsY"]->fill(estim->getY(), plane.getMeasZ() - meas.getZ()  );

/*
      _aidaHistoMap2D[bname + "residualXdX"]->fill(estim->getX(), estim->getX() - meas.getX() );
      _aidaHistoMap2D[bname + "residualYdX"]->fill(estim->getX(), estim->getY() - meas.getY() );
      _aidaHistoMap2D[bname + "residualXdY"]->fill(estim->getY(), estim->getX() - meas.getX() );
      _aidaHistoMap2D[bname + "residualYdY"]->fill(estim->getY(), estim->getY() - meas.getY() );
      _aidaHistoMap2D[bname + "residualdZvsY"]->fill(estim->getY(), plane.getMeasZ() - meas.getZ()  );
*/
      _aidaHistoMap2D[bname + "residualmeasZvsmeasX"]->fill(  meas.getZ()/1000., meas.getX()  );
      _aidaHistoMap2D[bname + "residualmeasZvsmeasY"]->fill(  meas.getZ()/1000., meas.getY()  );
      _aidaHistoMap2D[bname + "residualfitZvsmeasX"]->fill( plane.getMeasZ()/1000., meas.getX() );
      _aidaHistoMap2D[bname + "residualfitZvsmeasY"]->fill( plane.getMeasZ()/1000., meas.getY() );
 
      _aidaHistoMap2D[ "AllResidmeasZvsmeasX"]->fill(  meas.getZ()/1000., meas.getX()  );
      _aidaHistoMap2D[ "AllResidmeasZvsmeasY"]->fill(  meas.getZ()/1000., meas.getY()  );
      _aidaHistoMap2D[ "AllResidfitZvsmeasX"]->fill( plane.getMeasZ()/1000., meas.getX() );
      _aidaHistoMap2D[ "AllResidfitZvsmeasY"]->fill( plane.getMeasZ()/1000., meas.getY() );
    //Angles
      _aidaHistoMap[bname + "dxdz"]->fill( estim->getXdz() );
      _aidaHistoMap[bname + "dydz"]->fill( estim->getYdz() );
      if( ii != 4) { continue; }
     //plane.setMeasZ( plane.getMeasZ() - 175.0);
      _aidaZvHitX->fill(estim->getX(), meas.getZ() - plane.getZpos());
      _aidaZvFitX->fill(estim->getX(), (plane.getMeasZ() - plane.getZpos()) - (meas.getZ() - plane.getZpos()));
      _aidaZvHitY->fill(estim->getY(), meas.getZ() - plane.getZpos());
//      _aidaZvHitY->fill(estim->getY(), plane.getMeasZ() - plane.getZpos());
      _aidaZvFitY->fill(estim->getY(), (plane.getMeasZ() - plane.getZpos()) - (meas.getZ() - plane.getZpos()));
    }
  }
}

void EUTelDafBase::fillDetailPlots(daffitter::TrackCandidate *track){
  for( size_t ii = 0; ii < _system.planes.size() ; ii++){
    daffitter::FitPlane& plane = _system.planes.at(ii);

    daffitter::TrackEstimate* estim = track->estimates.at(ii);

    char iden[4];
    sprintf(iden, "%d", plane.getSensorID());
    string bname = (string)"pl" + iden + "_";

    //Plot resids, angles for all hits with > 50% includion in track.
    //This should be one measurement per track
    for(size_t w = 0; w < plane.meas.size(); w++){
      daffitter::Measurement& meas = plane.meas.at(w);
      if( plane.weights(w) < 0.5f) { continue; }
      //Resids 
      float resX = ( estim->getX() - meas.getX() );
      resX *= resX;
      resX /= plane.getSigmaX() *  plane.getSigmaX() + estim->cov(0,0);
      float resY = ( estim->getY() - meas.getY() );
      resY *= resY;
      resY /= plane.getSigmaY() *  plane.getSigmaY() + estim->cov(1,1);
      _aidaHistoMap[bname + "hitChi2"]->fill( resX + resY );
      
      _aidaHistoMap[bname + "sigmaX"]->fill( sqrt(estim->cov(0,0)) );
      _aidaHistoMap[bname + "sigmaY"]->fill( sqrt(estim->cov(1,1)) );
      
      float pullX =  ( estim->getX() - meas.getX() ) / sqrt(plane.getSigmaX() * plane.getSigmaX() + estim->cov(0,0));
      float pullY =  ( estim->getY() - meas.getY() ) / sqrt(plane.getSigmaY() * plane.getSigmaY() + estim->cov(1,1));
      _aidaHistoMap[bname + "pullX"]->fill( pullX );
      _aidaHistoMap[bname + "pullY"]->fill( pullY );
    }
  }
}

void EUTelDafBase::bookHistos(){
  int maxNdof = -4 + _system.planes.size() * 2 + 1;
  _aidaHistoMap["chi2"] = AIDAProcessor::histogramFactory(this)->createHistogram1D("chi2", 1000, 0, maxNdof * _maxChi2);
  _aidaHistoMap["logchi2"] = AIDAProcessor::histogramFactory(this)->createHistogram1D("logchi2", 1000, 0, std::log10(maxNdof * _maxChi2));
  if( _aidaHistoMap["chi2"] == NULL){
    streamlog_out ( ERROR2 ) << "Problem with histo booking. Check paths!" << std::endl;
    _histogramSwitch = false;
    return;
  }
  _aidaHistoMap["ndof"] = AIDAProcessor::histogramFactory(this)->createHistogram1D("ndof", maxNdof * 10, 0, maxNdof);
  _aidaHistoMap["chi2overndof"] = AIDAProcessor::histogramFactory(this)->createHistogram1D("Chi2OverNdof", maxNdof * 10, 0, _maxChi2);

  _aidaZvFitX = AIDAProcessor::histogramFactory(this)->createHistogram2D("ZvHitX", 1000,  -5000.0,  5000.0, 1000,   -100.0,   100.0);
  _aidaZvHitX = AIDAProcessor::histogramFactory(this)->createHistogram2D("ZvFitX", 1000, -10000.0, 10000.0, 1000, -10000.0, 10000.0);
  _aidaZvFitY = AIDAProcessor::histogramFactory(this)->createHistogram2D("ZvHitY", 1000, -10000.0, 10000.0, 1000,   -100.0,   100.0);
  _aidaZvHitY = AIDAProcessor::histogramFactory(this)->createHistogram2D("ZvFitY", 1000, -10000.0, 10000.0, 1000, -10000.0, 10000.0);

  _aidaHistoMap2D["AllResidmeasZvsmeasX"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( "AllResidmeasZvsmeasX",1000 ,0., 1000., 1000 ,-10000., 10000.);
  _aidaHistoMap2D["AllResidmeasZvsmeasY"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( "AllResidmeasZvsmeasY",1000 ,0., 1000., 1000 ,-10000., 10000.);
  _aidaHistoMap2D["AllResidfitZvsmeasX"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( "AllResidfitZvsmeasX",1000 ,0., 1000., 1000 ,-10000., 10000.);
  _aidaHistoMap2D["AllResidfitZvsmeasY"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( "AllResidfitZvsmeasY",1000 ,0., 1000., 1000 ,-10000., 10000.);


  for( size_t ii = 0; ii < _system.planes.size() ; ii++){
    daffitter::FitPlane& plane = _system.planes.at(ii);
    char iden[4];
    sprintf(iden, "%d", plane.getSensorID());
    string bname = (string)"pl" + iden + "_";
    //Resids 
    _aidaHistoMap[bname + "residualX"] =  AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "residualX", 1000 , -1000., 1000.);
    _aidaHistoMap[bname + "residualY"] =  AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "residualY", 1000 , -1000., 1000.);
    //Resids 2D // profiles
    _aidaHistoMapProf1D[bname+"residualdXvsX"]= AIDAProcessor::histogramFactory(this)->createProfile1D(bname+"dXvsX", 1000, -10000., 10000.,  -1000., 1000. );
    _aidaHistoMapProf1D[bname+"residualdYvsX"]= AIDAProcessor::histogramFactory(this)->createProfile1D(bname+"dXvsY", 1000, -10000., 10000.,  -1000., 1000. );
    _aidaHistoMapProf1D[bname+"residualdXvsY"]= AIDAProcessor::histogramFactory(this)->createProfile1D(bname+"dYvsX", 1000, -10000., 10000.,  -1000., 1000. );
    _aidaHistoMapProf1D[bname+"residualdYvsY"]= AIDAProcessor::histogramFactory(this)->createProfile1D(bname+"dYvsY", 1000, -10000., 10000.,  -1000., 1000. );
    _aidaHistoMapProf1D[bname+"residualdZvsX"]= AIDAProcessor::histogramFactory(this)->createProfile1D(bname+"dZvsX", 1000, -10000., 10000.,  -100., 100. );
    _aidaHistoMapProf1D[bname+"residualdZvsY"]= AIDAProcessor::histogramFactory(this)->createProfile1D(bname+"dZvsY", 1000, -10000., 10000.,  -100., 100. );
 
//    _aidaHistoMap2D[bname + "residualXdX"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( bname + "residualXdX", 1000 , -10000., 10000., 1000 , -1000., 1000.);
//    _aidaHistoMap2D[bname + "residualYdX"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( bname + "residualYdX", 1000 , -10000., 10000., 1000 , -1000., 1000.);
//    _aidaHistoMap2D[bname + "residualXdY"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( bname + "residualXdY", 1000 , -10000., 10000., 1000 , -1000., 1000.);
//    _aidaHistoMap2D[bname + "residualYdY"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( bname + "residualYdY", 1000 , -10000., 10000., 1000 , -1000., 1000.);
//  _aidaHistoMap2D[bname + "residualdZvsX"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( bname + "residualdZvsX",1000 ,-10000., 10000., 1000 ,-50., 50.);
//  _aidaHistoMap2D[bname + "residualdZvsY"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( bname + "residualdZvsY",1000 ,-10000., 10000., 1000 ,-50., 50.);

    // residuals
    _aidaHistoMap2D[bname + "residualmeasZvsmeasX"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( bname + "residualmeasZvsmeasX",1000 ,0., 1000., 1000 ,-10000., 10000.);
    _aidaHistoMap2D[bname + "residualmeasZvsmeasY"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( bname + "residualmeasZvsmeasY",1000 ,0., 1000., 1000 ,-10000., 10000.);
    _aidaHistoMap2D[bname + "residualfitZvsmeasX"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( bname + "residualfitZvsmeasX",1000 ,0., 1000., 1000 ,-10000., 10000.);
    _aidaHistoMap2D[bname + "residualfitZvsmeasY"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( bname + "residualfitZvsmeasY",1000 ,0., 1000., 1000 ,-10000., 10000.);

    //Angles
    _aidaHistoMap[bname + "dxdz"] = AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "dxdz", 10000, -0.1, 0.1);
    _aidaHistoMap[bname + "dydz"] = AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "dydz", 10000, -0.1, 0.1);
    //Zpositions
  }
}

void EUTelDafBase::bookDetailedHistos(){
  for( size_t ii = 0; ii < _system.planes.size() ; ii++){
    daffitter::FitPlane& plane = _system.planes.at(ii);
    char iden[4];
    sprintf(iden, "%d", plane.getSensorID());
    string bname = (string)"pl" + iden + "_";
    _aidaHistoMap[bname + "sigmaX"] =  AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "sigmaX", 10000, 0.0f, 100);
    _aidaHistoMap[bname + "sigmaY"] =  AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "sigmaY", 10000, 0.0f, 100);
    _aidaHistoMap[bname + "hitChi2"] =  AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "hitChi2", 300, 0, 300);
    _aidaHistoMap[bname + "pullX"] =  AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "pullX", 10000, -100, 100);
    _aidaHistoMap[bname + "pullY"] =  AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "pullY", 10000, -100, 100);
  }
}

void EUTelDafBase::end() {
  cout << "Trying to close " << _asciiName << endl;
  trackstream.close();

  dafEnd();
  
  streamlog_out ( MESSAGE ) << endl;
  streamlog_out ( MESSAGE ) << "Number of found hit clusters: " << _nClusters << endl;
  streamlog_out ( MESSAGE ) << "Tracks with ok ndof: " << n_passedNdof << endl;
  streamlog_out ( MESSAGE ) << "Tracks with ok chi2/ndof: " << n_passedChi2OverNdof << endl;
  streamlog_out ( MESSAGE ) << "Tracks with no NaNs: " << n_passedIsnan<< endl;
  streamlog_out ( MESSAGE ) << "Number of fitted tracks: " << _nTracks << endl;
  streamlog_out ( MESSAGE ) << "Successfully finished" << endl;
  for( size_t ii = 0; ii < _system.planes.size() ; ii++){
    daffitter::FitPlane& plane = _system.planes.at(ii);
    char iden[4];
    sprintf(iden, "%d", plane.getSensorID());
    string bname = (string)"pl" + iden + "_";
    streamlog_out ( MESSAGE ) << "plane:" << ii <<
                               "  x-stat :" <<  _aidaHistoMap[bname + "residualX"]->allEntries() <<
                               "  x-mean:"  <<  _aidaHistoMap[bname + "residualX"]->mean() << 
                               "  x-rms :"  <<  _aidaHistoMap[bname + "residualX"]->rms() << 
                               "  y-stat :" <<  _aidaHistoMap[bname + "residualY"]->allEntries() <<
                               "  y-mean:"  <<  _aidaHistoMap[bname + "residualY"]->mean() << 
                               "  y-rms :"  <<  _aidaHistoMap[bname + "residualY"]->rms() << endl;
  }


}
#endif // USE_GEAR
