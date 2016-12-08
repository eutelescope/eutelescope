// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#if defined USE_GEAR

// EUTelescope includes:
#include "EUTelTripletGBL.h"
#include "EUTelTripletGBLInstance.h"

#include "EUTELESCOPE.h"
//#include "EUTelSparseDataImpl.h"
#include "EUTelBaseSparsePixel.h"
#include "EUTelGenericSparsePixel.h"
#include "EUTelTrackerDataInterfacerImpl.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelPStream.h" // process streams redi::ipstream

// for clustersize
#include "EUTelSimpleVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelGenericSparseClusterImpl.h"
#include "EUTelGeometricClusterImpl.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

// AIDA histogram package (on top of ROOT):

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IAxis.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/IProfile2D.h>
#include <AIDA/ITree.h>
#endif

// GBL:
#include "include/GblTrajectory.h"
#include "include/MilleBinary.h"


// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/ProcessorMgr.h"
#include "marlin/Global.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// LCIO includes <.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>

// system includes <>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <memory>
#include <string.h>
#include <map>
#include <cstdlib>
#include <limits>

// ROOT includes ".h"
#include <TMath.h>
#include <TVectorD.h>
#include <TF1.h>
#include <TMatrixD.h>
#include <TVector3.h>
#include <TRotation.h>
#include "TH1D.h"

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;


EUTelTripletGBL::EUTelTripletGBL() : Processor("EUTelTripletGBL"), _siPlanesParameters(), _siPlanesLayerLayout(), _inputCollectionTelescope(""), _isFirstEvent(0), _eBeam(0), _nEvt(0), _nTelPlanes(0), _dut_plane(3), _cutx(0.15), _cuty(0.15), _planeSort(), _planeID(), _planePosition(), _planeThickness(), _planeX0(), _planeResolution() {

  // modify processor description
  _description = "Analysis for DATURA reference analysis ";

  // processor parameters
  registerInputCollection( LCIO::TRACKERHIT,
      "InputCollectionTelescope" ,
      "Name of the input TrackerHit collection of the telescope",
      _inputCollectionTelescope,
      std::string("") ); // no collection defaulted forcing the user to pass something meaningful

  registerProcessorParameter( "Ebeam",
      "Beam energy [GeV]",
      _eBeam, static_cast < double >( 0.0));

  registerOptionalParameter( "triCut", "Upstream/Downstream triplet residual cut [mm]", _triCut, 0.1 );

  registerProcessorParameter( "matching_cut_x",
      "cut for matching in x coordinate in mm",
      _cutx, static_cast < double >(0.15));
  registerProcessorParameter( "matching_cut_y",
      "cut for matching in y coordinate in mm",
      _cuty, static_cast < double >(0.15));

  registerProcessorParameter( "slope_cut_x",
      "cut for track slopes in x coordinate in rad",
      _slope_cut_x, static_cast < double >(0.002));
  registerProcessorParameter( "slope_cut_y",
      "cut for track slopes in y coordinate in rad",
      _slope_cut_y, static_cast < double >(0.002));

  registerProcessorParameter( "dut_plane",
      "plane to be considered the DUT and excluded from the track fit",
      _dut_plane, static_cast <int>(3));

  registerProcessorParameter( "eff_radius",
      "radius on DUT plane to accept match with triplet",
      _eff_radius, static_cast <double>(0.1)); // in mm. 0.1 is for 6 GeV, 20 mm. This is scaled with E and dz

  registerProcessorParameter( "kappa",
      "global factor to Highland formula",
      _kappa, static_cast <double>(1.0)); // 1.0 means HL as is, 1.2 means 20% additional scattering

  registerProcessorParameter( "aluthickum",
      "thickness of alu target, if present",
      _aluthickum, static_cast <double>(0.0));

  registerProcessorParameter( "probchi2_cut",
      "Cut on Prob(chi2,ndf) rejecting bad tracks with prob < cut",
      _probchi2_cut, static_cast <double>(.01)); 

  registerOptionalParameter("Resolution",
      "resolution parameter for each Cluster size, same for all planes. first value is average of all CSes. Up to CS6 plus larger than 6, hence in total 8 numbers. Disable with -1, e.g. (3.5e-3, -1, -1, ...., -1)",
      _resolution,  FloatVec (static_cast <double> (8), 3.5*1e-3));

  registerOptionalParameter("Thickness","thickness parameter for each plane. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_thickness,  FloatVec (static_cast <double> (6), 50*1e-3));


}


void EUTelTripletGBL::init() {

  // usually a good idea to
  printParameters();

  _nEvt = 0;

  _isFirstEvent = true;

  streamlog_out(MESSAGE0) << "Beam energy " << _eBeam << " GeV" <<  std::endl;

  // Read geometry information from GEAR

  streamlog_out( MESSAGE0 ) << "Reading telescope geometry description from GEAR " << std::endl;

  _siPlanesParameters  = const_cast<gear::SiPlanesParameters* >( &(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*>(  &(_siPlanesParameters->getSiPlanesLayerLayout() ));

  // Take all layers defined in GEAR geometry
  _nTelPlanes = _siPlanesLayerLayout->getNLayers();

  _planeSort = new int[_nTelPlanes];
  _planePosition = new double[_nTelPlanes]; // z pos
  _planeID         = new int[_nTelPlanes];
  _planeThickness  = new double[_nTelPlanes];
  _planeX0         = new double[_nTelPlanes];
  _planeResolution = new double[_nTelPlanes];

  for(int i = 0; i < _nTelPlanes; i++) {

    _planeID[i]=_siPlanesLayerLayout->getID(i);
    _planePosition[i]   =_siPlanesLayerLayout->getLayerPositionZ(i);
    //_planeThickness[i]  =_siPlanesLayerLayout->getLayerThickness(i);
    _planeThickness[i]  = _thickness[i];
    streamlog_out( MESSAGE6 ) << "  Thickness Plane " << i << " = " << _thickness[i] << std::endl;
    _planeX0[i]         =_siPlanesLayerLayout->getLayerRadLength(i);
    //_planeResolution[i] = _siPlanesLayerLayout->getSensitiveResolution(i); // Get it from input
    _planeResolution[i] = _resolution[0]; // pass the average here, which is the first entry of the vector
    streamlog_out( MESSAGE6 ) << "  Avg. reso Plane " << i << " = " << _resolution[0] << std::endl;
  }

  streamlog_out( MESSAGE2 ) <<  "Telescope configuration with " << _nTelPlanes << " planes" << std::endl;

  for( int ipl = 0; ipl < _nTelPlanes; ipl++) {
    std::stringstream ss;
    ss << "  ID = " << _planeID[ipl]
      << "  at Z [mm] = " << _planePosition[ipl]
      << " dZ [um] = " << _planeThickness[ipl]*1000.;

    ss << "  average Res [um] = " << _planeResolution[ipl]*1000.;

    streamlog_out( MESSAGE2 ) <<  ss.str() << std::endl;

  }


  if( _resolution.size() != 8 )
  {
    throw InvalidParameterException("WARNING, length of resolution vector is != 8. You need to pass 8 resolutions (-1 if not known). \n");
  }
  if ( _resolution.at(0) < 0) {
    streamlog_out( ERROR2 ) << " Found avg resolution < 0, namely = " << _resolution.size() << 
      "\n   will have to termiante." << std::endl;
    return;
  }
  for (int i = 1; i < 8; i++) if(_resolution[i] < 0) _resolution.at(i) = _resolution.at(0);

  for (int i = 1; i < 8; i++) // not from 0 , 0 is average
    streamlog_out( MESSAGE6 ) <<  "CS resolutions: CS" << i << " = " << _resolution.at(i) << std::endl;

  //_triCut = _triCut *6. / _eBeam * (_planePosition[1] - _planePosition[0]) / 20.;

  // Book histograms:
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  bookHistos();
#endif
}//init

//------------------------------------------------------------------------------
void EUTelTripletGBL::processRunHeader( LCRunHeader* runHeader) {

  std::auto_ptr<EUTelRunHeaderImpl> eutelHeader( new EUTelRunHeaderImpl( runHeader ) );
  eutelHeader->addProcessor( type() );

  // Decode and print out Run Header information - just a check

  _nRun = runHeader->getRunNumber();

  streamlog_out( MESSAGE2 )  << "Processing run header"
    << ", run nr " << runHeader->getRunNumber() << std::endl;

  const std::string detectorName = runHeader->getDetectorName();
  const std::string detectorDescription = runHeader->getDescription();

  streamlog_out( MESSAGE0 ) << detectorName << " : " << detectorDescription << std::endl;

} // processRunHeader

//----------------------------------------------------------------------------
void EUTelTripletGBL::processEvent( LCEvent * event ) {

  if( _nEvt % 1000 == 0 ) {
    streamlog_out( MESSAGE2 ) << "Processing event "
      << setw(6) << setiosflags(ios::right)
      << event->getEventNumber() << " in run "
      << setw(6) << setiosflags(ios::right)
      << event->getRunNumber()
      << endl;
  }

  EUTelEventImpl * euEvent = static_cast<EUTelEventImpl*>( event );

  if( euEvent->getEventType() == kEORE ) {
    streamlog_out( DEBUG5 ) <<  "EORE found: nothing else to do." << std::endl;
    return;
  }

  int runNumber = event->getRunNumber();

  if( _isFirstEvent ) {

    // apply all GEAR/alignment offsets to get corrected X,Y,Z
    // position of the sensor center

    for( int iplane = 0; iplane < _siPlanesLayerLayout->getNLayers(); iplane++ ) {

      std::map< unsigned int , double > _planeCenter;
      std::map< unsigned int , double > _planeNormal;

      // start filling the map with Gear values:
      _planeCenter[ 0 ] =  _siPlanesLayerLayout->getLayerPositionX(iplane); // X
      _planeCenter[ 1 ] =  _siPlanesLayerLayout->getLayerPositionY(iplane); // Y
      _planeCenter[ 2 ] =  _siPlanesLayerLayout->getLayerPositionZ(iplane); // Z
      _planeNormal[ 0 ] =  0.; // X
      _planeNormal[ 1 ] =  0.; // Y
      _planeNormal[ 2 ] =  1.; // Z

      TVector3  _normalTVec( _planeNormal[0], _planeNormal[1], _planeNormal[2]); 

      // do initial rotation from GEAR
      try{
	double gRotation[3] = { 0., 0., 0.}; // not rotated

	gRotation[0] = _siPlanesLayerLayout->getLayerRotationXY(iplane); // Euler alpha
	gRotation[1] = _siPlanesLayerLayout->getLayerRotationZX(iplane); // Euler alpha
	gRotation[2] = _siPlanesLayerLayout->getLayerRotationZY(iplane); // Euler alpha

	// input angles are in degree, translate into radian
	gRotation[0] =  gRotation[0]*3.1415926/180.;
	gRotation[1] =  gRotation[1]*3.1415926/180.;
	gRotation[2] =  gRotation[2]*3.1415926/180.;

	TRotation r;
	r.RotateX( gRotation[2] );
	r.RotateY( gRotation[1] );
	r.RotateZ( gRotation[0] );
	_normalTVec.Transform( r );

	_planeNormal[0] = _normalTVec[0];
	_planeNormal[1] = _normalTVec[1];
	_planeNormal[2] = _normalTVec[2];
      }
      catch(...) {
	streamlog_out(DEBUG5) << "no sensor rotation is given in the GEAR steering file, assume NONE" << std::endl;
      }
    }
    if( isFirstEvent() ) _isFirstEvent = false;

  } // FirstEvent

  // Increase event count
  _nEvt++;

  //----------------------------------------------------------------------------
  // check input collection (aligned hits):

  LCCollection *collection;
  try {
    collection = event->getCollection( _inputCollectionTelescope );
  }
  catch( lcio::DataNotAvailableException& e) {
    streamlog_out( DEBUG1 ) << "Not able to get collections "
      << _inputCollectionTelescope << " "
      << "\nfrom event " << event->getEventNumber()
      << " in run " << runNumber  << std::endl;

    return;
  }


  //----------------------------------------------------------------------------
  // Copy hits to local table
  // Assign hits to sensor planes

  if(_nEvt < 10) streamlog_out( WARNING2 )  << "Total of " << collection->getNumberOfElements() << " tracker hits in input collection " << std::endl;
  nAllHitHisto->fill(collection->getNumberOfElements());

  //----------------------------------------------------------------------------

  std::vector<hit> hits;

  // Extract all hits from the LCIO collection and calculate the
  // position uncertainty if necessary:
  for( int ihit = 0; ihit < collection->getNumberOfElements(); ihit++ ) {

    hit newhit;
    TrackerHit * meshit = dynamic_cast<TrackerHit*>( collection->getElementAt(ihit) );
    const double * pos = meshit->getPosition();
    const EVENT::FloatVec cov = meshit->getCovMatrix();
    float charge = 0.;

    TrackerDataImpl* clusterVector = static_cast<TrackerDataImpl*>( meshit->getRawHits()[0]);
    EUTelSimpleVirtualCluster * cluster=0;

    float locx = -1;
    float locy = -1;
    if ( meshit->getType() == kEUTelSparseClusterImpl ) 
    {
      cluster = new EUTelSparseClusterImpl< EUTelGenericSparsePixel > ( clusterVector );
      charge = cluster->getTotalCharge();
      cluster->getCenterOfGravity(locx, locy);
      streamlog_out(DEBUG4) << " charge is " << charge << std::endl;
    }


    // Write the position:
    newhit.x = pos[0];
    newhit.y = pos[1];
    newhit.z = pos[2];
    newhit.locx = locx;
    newhit.locy = locy;

    //if(_nEvt < 10) streamlog_out( WARNING2 ) << "hit x = " << meshit->getPosition()[0] << endl;

    // Write clustersize
    newhit.clustersize = charge;
    delete cluster;

    // Find Plane ID to which the hit belongs by minimizing its
    // distance in Z
    //FIXME to be fixed with more general geometry description!
    double distMin = 99; // [mm]

    bool foundplane = false;
    for( int ipl = 0; ipl < _nTelPlanes; ipl++ ) {
      if( abs(newhit.z - _planePosition[ipl]) < distMin ) {
	newhit.plane = ipl;
	foundplane = true;
	distMin = abs(newhit.z - _planePosition[ipl]);
      }
    }

    if(!foundplane) {
      streamlog_out(ERROR5) << "Could not associate hit with a telescope plane. Skipping event." << std::endl;
      return;
    }

    // ...and uncertainty:
    if( cov.at(0) > 0. ) newhit.ex = sqrt(cov.at(0));
    else newhit.ex = _planeResolution[newhit.plane];

    if( cov.at(2) > 0. ) newhit.ey = sqrt(cov.at(2));
    else newhit.ey = _planeResolution[newhit.plane];

    streamlog_out(DEBUG3) << "Hit " << ihit << ": " << newhit << std::endl;

    // Use only hits with clustersize < 5 (>4 means there was a noise hit in the cluster or something else was fishy)
    //if(newhit.clustersize > 4) continue;

    // Add it to the vector of telescope hits:
    hits.push_back(newhit);

    //delete clusterVector;
    //delete meshit;
    //delete pos;

  } // loop over hits

  streamlog_out(DEBUG4) << "Event " << event->getEventNumber()
    << " contains " << hits.size() << " hits" << std::endl;


  // Fill the telescope plane correlation plots:
  TelescopeCorrelationPlots(hits);


  // -----------------------------------------------------------------
  // Downstream Telescope Triplets ("driplets")

  // Generate new triplet set for the Telescope Downstream Arm:
  std::vector<triplet> downstream_triplets;
  FindTriplets(hits, 3, 4, 5, downstream_triplets);

  // Iterate over all found downstream triplets to fill histograms and match them to the REF and DUT:
  for( std::vector<triplet>::iterator drip = downstream_triplets.begin(); drip != downstream_triplets.end(); drip++ ){

    // Fill some histograms for downstream triplets:
    dridxHisto->fill( (*drip).getdx(4)*1E3 ); 
    dridxvsx->fill( (*drip).base().x, (*drip).getdx(4)*1E3 ); // check for rot
    dridxvsy->fill( (*drip).base().y, (*drip).getdx(4)*1E3 );
    dridxvstx->fill( (*drip).slope().x*1E3, (*drip).getdx(4)*1E3 ); // check for z shift
    dridxvsty->fill( (*drip).slope().y*1E3, (*drip).getdx(4)*1E3 );

    dridyHisto->fill( (*drip).getdy(4)*1E3 );
    dridyvsx->fill( (*drip).base().x, (*drip).getdy(4)*1E3 );
    dridyvsy->fill( (*drip).base().y, (*drip).getdy(4)*1E3 );
    dridyvstx->fill( (*drip).slope().x*1E3, (*drip).getdy(4)*1E3 );
    dridyvsty->fill( (*drip).slope().y*1E3, (*drip).getdy(4)*1E3 );

    drixHisto->fill( -(*drip).gethit(4).x ); // -x = x_DP = out
    driyHisto->fill( -(*drip).gethit(4).y ); // -y = y_DP =  up
    drixyHisto->fill( -(*drip).gethit(4).x, -(*drip).gethit(4).y );
    dritxHisto->fill( (*drip).slope().x*1E3 );
    drityHisto->fill( (*drip).slope().y*1E3 );

  }//loop over driplets

  ndriHisto->fill( downstream_triplets.size() );


  //##########################  0-1-2 upstream triplets #######################
  // here again a telescope hit triplet is formed, from planes 0 and 2
  // and the correlation to plane 1. the found triplets are
  // extrapolated and matched to the DUT plane and exhaustive
  // histograms are written.
  // Furthermore the hits of a triplet matched to the DUT are
  // collected and fed into GBL for a track fit.

  // Generate new triplet set for the Telescope Upstream Arm:
  std::vector<triplet> upstream_triplets;
  FindTriplets(hits, 0, 1, 2, upstream_triplets);

  // Iterate over all found upstream triplets to fill histograms and match them to the REF and DUT:
  for( std::vector<triplet>::iterator trip = upstream_triplets.begin(); trip != upstream_triplets.end(); trip++ ) {

    // Fill some histograms for the upstream triplets:
    tridxHisto->fill( (*trip).getdx(1)*1E3 );
    tridyHisto->fill( (*trip).getdy(1)*1E3 );
    tridx1Histo->fill( (*trip).getdx(1)*1E0 );
    tridy1Histo->fill( (*trip).getdy(1)*1E0 );
    tridxvsx->fill( (*trip).base().x, (*trip).getdx(1)*1E3 ); // check for rot
    tridxvsy->fill( (*trip).base().y, (*trip).getdx(1)*1E3 );
    tridxvstx->fill( (*trip).slope().x*1E3, (*trip).getdx(1)*1E3 ); // check for z shift
    tridxvsty->fill( (*trip).slope().y*1E3, (*trip).getdx(1)*1E3 );
    tridyvsx->fill( (*trip).base().x, (*trip).getdy(1)*1E3 );
    tridyvsy->fill( (*trip).base().y, (*trip).getdy(1)*1E3 );
    tridyvstx->fill( (*trip).slope().x*1E3, (*trip).getdy(1)*1E3 );
    tridyvsty->fill( (*trip).slope().y*1E3, (*trip).getdy(1)*1E3 );
    trixHisto->fill( -(*trip).gethit(1).x );
    triyHisto->fill( -(*trip).gethit(1).y );
    trixyHisto->fill( -(*trip).gethit(1).x, -(*trip).gethit(1).y );
    tritxHisto->fill( (*trip).slope().x*1E3 );
    trityHisto->fill( (*trip).slope().y*1E3 );


    // Extrapolate Upstream triplet to Downstream planes 3,4,5: Resolution studies
    for( std::vector<hit>::iterator lhit = hits.begin(); lhit != hits.end(); lhit++ ){

      if( (*lhit).plane <= 2 ) continue; // want 3,4, or 5

      // Fill residuals of triplet and hit in the selected plane:
      if( (*lhit).plane == 3 ) {
	tridx3Histo->fill( (*trip).getdx((*lhit))*1E0 );
	tridy3Histo->fill( (*trip).getdy((*lhit))*1E0 ); // 65 um at 4.7 GeV with CMS
	tridx3bHisto->fill( (*trip).getdx((*lhit))*1E3 ); // finer binning
	tridy3bHisto->fill( (*trip).getdy((*lhit))*1E3 ); // 
      }
      else if( (*lhit).plane == 4 ) {
	tridx4Histo->fill( (*trip).getdx((*lhit))*1E0 );
	tridy4Histo->fill( (*trip).getdy((*lhit))*1E0 ); //174 um at 4.7 GeV
	tridx4bHisto->fill( (*trip).getdx((*lhit))*1E3 ); // finer binning
	tridy4bHisto->fill( (*trip).getdy((*lhit))*1E3 ); // 
      }
      else if( (*lhit).plane == 5 ) {
	tridx5Histo->fill( (*trip).getdx((*lhit))*1E0 );
	tridy5Histo->fill( (*trip).getdy((*lhit))*1E0 ); //273 um at 4.7 GeV
	tridx5bHisto->fill( (*trip).getdx((*lhit))*1E3 ); // finer binning
	tridy5bHisto->fill( (*trip).getdy((*lhit))*1E3 ); // 
      }
    }// Resolution studies



  }// iterate over upstream triplets
  ntriHisto->fill( upstream_triplets.size() );



  //----------------------------------------------------------------------------
  // calculate efficiency of plane 3 by forming a triplet from planes 0, 1, 2; 2, 4, 5.
  // Then try to find a match on DUT (plane 3)
  //
  //
  double DUTz = _planePosition[3];
  // This scales the radius with beam energy spacing (efficiencty shouldnt depend on amount of scattering!). Define radius at 6 GeV, 20 mm
  double eff_radius = _eff_radius * 6. / _eBeam * (_planePosition[1] - _planePosition[0]) / 20.; 


  // Generate new triplet set with planes 0, 1, 2; 2,4,5:
  std::vector<triplet> eff_triplets_UP = upstream_triplets;
  //FindTriplets(hits, 0, 1, 2, eff_triplets_UP);

  std::vector<triplet> eff_triplets_DOWN;
  FindTriplets(hits, 2, 4, 5, eff_triplets_DOWN);

  std::vector<triplet> eff_triplets;

  // Iterate over all found eff-triplets to match them to the DUT (plane3):
  //int n_matched_trips = 0;
  //int n_unmatched_trips = 0;

  //std::cout << " n eff triplets UP   = " << eff_triplets_UP->size() << std::endl;
  //std::cout << " n eff triplets DOWN = " << eff_triplets_DOWN->size() << std::endl;

  for( std::vector<triplet>::iterator trip = eff_triplets_UP.begin(); trip != eff_triplets_UP.end(); trip++ ) {

    // Track impact position at Matching Point from Upstream:
    double xA = (*trip).getx_at(DUTz); // triplet impact point at matching position
    double yA = (*trip).gety_at(DUTz);

    // check if trip is isolated
    bool IsolatedTrip = IsTripletIsolated(trip, eff_triplets_UP, DUTz, _triCut*2.01);

    for( std::vector<triplet>::iterator drip = eff_triplets_DOWN.begin(); drip != eff_triplets_DOWN.end(); drip++ ){

      // Track impact position at Matching Point from Downstream:
      double xB = (*drip).getx_at(DUTz); // triplet impact point at matching position
      double yB = (*drip).gety_at(DUTz);

      // check if drip is isolated
      bool IsolatedDrip = IsTripletIsolated(drip, eff_triplets_DOWN, DUTz, _triCut*2.01);

      // driplet - triplet
      double dx = xB - xA; 
      double dy = yB - yA;

      // match driplet and triplet:
      if( abs(dx) > _cutx) continue;
      if( abs(dy) > _cuty) continue;

      //std::cout << " intersec ";

      // check isolation
      if( !IsolatedTrip || !IsolatedDrip ) continue;
      //std::cout << " , isolated ";

      // apply fiducial cut
      if ( fabs(xA) >  9.0) continue;
      if (     -yA  < -4.0) continue;

      //std::cout << " , fiducial " << std::endl;

      eff_triplets.push_back(*trip);


    } // Downstream
  } // Upstream

  //std::cout << " n eff triplets = " << eff_triplets.size() << std::endl;

  for( std::vector<triplet>::iterator trip = eff_triplets.begin(); trip != eff_triplets.end(); trip++ ) {

    double ddAMin = -1000.0;
    double xTrip = (*trip).getx_at(DUTz);
    double yTrip = (*trip).gety_at(DUTz);

    for( std::vector<hit>::iterator lhit = hits.begin(); lhit != hits.end(); lhit++ ){

      if( (*lhit).plane <= 2 || (*lhit).plane  > 3) continue; // want 3

      // Fill residuals of triplet and hit in the selected plane:
      if( (*lhit).plane == 3 ) {
	double xHit = (*lhit).x;
	double yHit = (*lhit).y;

	double ddA = sqrt( fabs(xHit - xTrip)*fabs(xHit - xTrip) 
	    + fabs(yHit - yTrip)*fabs(yHit - yTrip) );
	if(ddAMin < 0 || ddA < ddAMin) ddAMin = ddA;
      }
    } // end loop over hits

    // if distance is smaller then limit, accept this as matched Hit
    if(fabs(ddAMin) < eff_radius) {
      //n_matched_trips++;
      effix3->fill(-xTrip, 1.);
      effiy3->fill(-yTrip, 1.);
    } else {
      effix3->fill(-xTrip, 0.);
      effiy3->fill(-yTrip, 0.);
    }

  }
  
  eff_triplets_UP.clear();
  eff_triplets_DOWN.clear();
  eff_triplets.clear();

  //----------------------------------------------------------------------------
  // calculate efficiency of plane 2 by forming a triplet from planes 0, 1, 3; 3, 4, 5.
  // Then try to find a match on DUT (plane 2)
  //
  //
  DUTz = _planePosition[3];
  // This scales the radius with beam energy spacing (efficiencty shouldnt depend on amount of scattering!). Define radius at 6 GeV, 20 mm


  // Generate new triplet set with planes 0, 1, 3; 3,4,5:
  FindTriplets(hits, 0, 1, 3, eff_triplets_UP);

  eff_triplets_DOWN = downstream_triplets;

  //std::vector<triplet> eff_triplets;

  // Iterate over all found eff-triplets to match them to the DUT (plane3):
  //int n_matched_trips = 0;
  //int n_unmatched_trips = 0;

  //std::cout << " n eff triplets UP   = " << eff_triplets_UP->size() << std::endl;
  //std::cout << " n eff triplets DOWN = " << eff_triplets_DOWN->size() << std::endl;

  for( std::vector<triplet>::iterator trip = eff_triplets_UP.begin(); trip != eff_triplets_UP.end(); trip++ ) {

    // Track impact position at Matching Point from Upstream:
    double xA = (*trip).getx_at(DUTz); // triplet impact point at matching position
    double yA = (*trip).gety_at(DUTz);

    // check if trip is isolated
    bool IsolatedTrip = IsTripletIsolated(trip, eff_triplets_UP, DUTz, _triCut*2.01);

    for( std::vector<triplet>::iterator drip = eff_triplets_DOWN.begin(); drip != eff_triplets_DOWN.end(); drip++ ){

      // Track impact position at Matching Point from Downstream:
      double xB = (*drip).getx_at(DUTz); // triplet impact point at matching position
      double yB = (*drip).gety_at(DUTz);

      // check if drip is isolated
      bool IsolatedDrip = IsTripletIsolated(drip, eff_triplets_DOWN, DUTz, _triCut*2.01);

      // driplet - triplet
      double dx = xB - xA; 
      double dy = yB - yA;

      // match driplet and triplet:
      if( abs(dx) > _cutx) continue;
      if( abs(dy) > _cuty) continue;

      //std::cout << " intersec ";

      // check isolation
      if( !IsolatedTrip || !IsolatedDrip ) continue;
      //std::cout << " , isolated ";

      // apply fiducial cut
      if ( fabs(xA) >  9.0) continue;
      if (     -yA  < -4.0) continue;

      //std::cout << " , fiducial " << std::endl;

      eff_triplets.push_back(*trip);


    } // Downstream
  } // Upstream

  //std::cout << " n eff triplets = " << eff_triplets.size() << std::endl;

  for( std::vector<triplet>::iterator trip = eff_triplets.begin(); trip != eff_triplets.end(); trip++ ) {

    double ddAMin = -1000.0;
    double xTrip = (*trip).getx_at(DUTz);
    double yTrip = (*trip).gety_at(DUTz);

    for( std::vector<hit>::iterator lhit = hits.begin(); lhit != hits.end(); lhit++ ){

      if( (*lhit).plane <= 1 || (*lhit).plane  > 2) continue; // want 2

      // Fill residuals of triplet and hit in the selected plane:
      if( (*lhit).plane == 2 ) {
	double xHit = (*lhit).x;
	double yHit = (*lhit).y;

	double ddA = sqrt( fabs(xHit - xTrip)*fabs(xHit - xTrip) 
	    + fabs(yHit - yTrip)*fabs(yHit - yTrip) );
	if(ddAMin < 0 || ddA < ddAMin) ddAMin = ddA;
      }
    } // end loop over hits

    // if distance is smaller then limit, accept this as matched Hit
    if(fabs(ddAMin) < eff_radius) {
      //n_matched_trips++;
      effix2->fill(-xTrip, 1.);
      effiy2->fill(-yTrip, 1.);
    } else {
      effix2->fill(-xTrip, 0.);
      effiy2->fill(-yTrip, 0.);
    }

  }

  eff_triplets_UP.clear();
  eff_triplets_DOWN.clear();
  eff_triplets.clear();

  //----------------------------------------------------------------------------
  // calculate efficiency of plane 1 by forming a triplet from planes 0, 2, 3; 3, 4, 5.
  // Then try to find a match on DUT (plane 1)
  //
  //
  DUTz = _planePosition[3];
  // This scales the radius with beam energy spacing (efficiencty shouldnt depend on amount of scattering!). Define radius at 6 GeV, 20 mm


  // Generate new triplet set with planes 0, 2, 3; 3,4,5:
  FindTriplets(hits, 0, 2, 3, eff_triplets_UP);

  eff_triplets_DOWN = downstream_triplets;

  //std::vector<triplet> eff_triplets;

  // Iterate over all found eff-triplets to match them to the DUT (plane3):
  //int n_matched_trips = 0;
  //int n_unmatched_trips = 0;

  //std::cout << " n eff triplets UP   = " << eff_triplets_UP->size() << std::endl;
  //std::cout << " n eff triplets DOWN = " << eff_triplets_DOWN->size() << std::endl;

  for( std::vector<triplet>::iterator trip = eff_triplets_UP.begin(); trip != eff_triplets_UP.end(); trip++ ) {

    // Track impact position at Matching Point from Upstream:
    double xA = (*trip).getx_at(DUTz); // triplet impact point at matching position
    double yA = (*trip).gety_at(DUTz);

    // check if trip is isolated
    bool IsolatedTrip = IsTripletIsolated(trip, eff_triplets_UP, DUTz, _triCut*2.01);

    for( std::vector<triplet>::iterator drip = eff_triplets_DOWN.begin(); drip != eff_triplets_DOWN.end(); drip++ ){

      // Track impact position at Matching Point from Downstream:
      double xB = (*drip).getx_at(DUTz); // triplet impact point at matching position
      double yB = (*drip).gety_at(DUTz);

      // check if drip is isolated
      bool IsolatedDrip = IsTripletIsolated(drip, eff_triplets_DOWN, DUTz, _triCut*2.01);

      // driplet - triplet
      double dx = xB - xA; 
      double dy = yB - yA;

      // match driplet and triplet:
      if( abs(dx) > _cutx) continue;
      if( abs(dy) > _cuty) continue;

      //std::cout << " intersec ";

      // check isolation
      if( !IsolatedTrip || !IsolatedDrip ) continue;
      //std::cout << " , isolated ";

      // apply fiducial cut
      if ( fabs(xA) >  9.0) continue;
      if (     -yA  < -4.0) continue;

      //std::cout << " , fiducial " << std::endl;

      eff_triplets.push_back(*trip);


    } // Downstream
  } // Upstream

  //std::cout << " n eff triplets = " << eff_triplets.size() << std::endl;

  for( std::vector<triplet>::iterator trip = eff_triplets.begin(); trip != eff_triplets.end(); trip++ ) {

    double ddAMin = -1000.0;
    double xTrip = (*trip).getx_at(DUTz);
    double yTrip = (*trip).gety_at(DUTz);

    for( std::vector<hit>::iterator lhit = hits.begin(); lhit != hits.end(); lhit++ ){

      if( (*lhit).plane <= 0 || (*lhit).plane  > 1) continue; // want 1

      // Fill residuals of triplet and hit in the selected plane:
      if( (*lhit).plane == 1 ) {
	double xHit = (*lhit).x;
	double yHit = (*lhit).y;

	double ddA = sqrt( fabs(xHit - xTrip)*fabs(xHit - xTrip) 
	    + fabs(yHit - yTrip)*fabs(yHit - yTrip) );
	if(ddAMin < 0 || ddA < ddAMin) ddAMin = ddA;
      }
    } // end loop over hits

    // if distance is smaller then limit, accept this as matched Hit
    if(fabs(ddAMin) < eff_radius) {
      //n_matched_trips++;
      effix1->fill(-xTrip, 1.);
      effiy1->fill(-yTrip, 1.);
    } else {
      effix1->fill(-xTrip, 0.);
      effiy1->fill(-yTrip, 0.);
    }

  }

  eff_triplets_UP.clear();
  eff_triplets_DOWN.clear();
  eff_triplets.clear();

  //----------------------------------------------------------------------------
  // calculate efficiency of plane 0 by forming a triplet from planes 1, 2, 3; 3, 4, 5.
  // Then try to find a match on DUT (plane 0)
  //
  //
  DUTz = _planePosition[3];
  // This scales the radius with beam energy spacing (efficiencty shouldnt depend on amount of scattering!). Define radius at 6 GeV, 20 mm


  // Generate new triplet set with planes 1, 2, 3; 3,4,5:
  FindTriplets(hits, 1, 2, 3, eff_triplets_UP);

  eff_triplets_DOWN = downstream_triplets;

  //std::vector<triplet> eff_triplets;

  // Iterate over all found eff-triplets to match them to the DUT (plane3):
  //int n_matched_trips = 0;
  //int n_unmatched_trips = 0;

  //std::cout << " n eff triplets UP   = " << eff_triplets_UP->size() << std::endl;
  //std::cout << " n eff triplets DOWN = " << eff_triplets_DOWN->size() << std::endl;

  for( std::vector<triplet>::iterator trip = eff_triplets_UP.begin(); trip != eff_triplets_UP.end(); trip++ ) {

    // Track impact position at Matching Point from Upstream:
    double xA = (*trip).getx_at(DUTz); // triplet impact point at matching position
    double yA = (*trip).gety_at(DUTz);

    // check if trip is isolated
    bool IsolatedTrip = IsTripletIsolated(trip, eff_triplets_UP, DUTz, _triCut*2.01);

    for( std::vector<triplet>::iterator drip = eff_triplets_DOWN.begin(); drip != eff_triplets_DOWN.end(); drip++ ){

      // Track impact position at Matching Point from Downstream:
      double xB = (*drip).getx_at(DUTz); // triplet impact point at matching position
      double yB = (*drip).gety_at(DUTz);

      // check if drip is isolated
      bool IsolatedDrip = IsTripletIsolated(drip, eff_triplets_DOWN, DUTz, _triCut*2.01);

      // driplet - triplet
      double dx = xB - xA; 
      double dy = yB - yA;

      // match driplet and triplet:
      if( abs(dx) > _cutx) continue;
      if( abs(dy) > _cuty) continue;

      //std::cout << " intersec ";

      // check isolation
      if( !IsolatedTrip || !IsolatedDrip ) continue;
      //std::cout << " , isolated ";

      // apply fiducial cut
      if ( fabs(xA) >  9.0) continue;
      if (     -yA  < -4.0) continue;

      //std::cout << " , fiducial " << std::endl;

      eff_triplets.push_back(*trip);


    } // Downstream
  } // Upstream

  //std::cout << " n eff triplets = " << eff_triplets.size() << std::endl;

  for( std::vector<triplet>::iterator trip = eff_triplets.begin(); trip != eff_triplets.end(); trip++ ) {

    double ddAMin = -1000.0;
    double xTrip = (*trip).getx_at(DUTz);
    double yTrip = (*trip).gety_at(DUTz);

    for( std::vector<hit>::iterator lhit = hits.begin(); lhit != hits.end(); lhit++ ){

      if( (*lhit).plane  > 0) continue; // want 0

      // Fill residuals of triplet and hit in the selected plane:
      if( (*lhit).plane == 0 ) {
	double xHit = (*lhit).x;
	double yHit = (*lhit).y;

	double ddA = sqrt( fabs(xHit - xTrip)*fabs(xHit - xTrip) 
	    + fabs(yHit - yTrip)*fabs(yHit - yTrip) );
	if(ddAMin < 0 || ddA < ddAMin) ddAMin = ddA;
      }
    } // end loop over hits

    // if distance is smaller then limit, accept this as matched Hit
    if(fabs(ddAMin) < eff_radius) {
      //n_matched_trips++;
      effix0->fill(-xTrip, 1.);
      effiy0->fill(-yTrip, 1.);
    } else {
      effix0->fill(-xTrip, 0.);
      effiy0->fill(-yTrip, 0.);
    }

  }

  eff_triplets_UP.clear();
  eff_triplets_DOWN.clear();
  eff_triplets.clear();

  //----------------------------------------------------------------------------
  // calculate efficiency of plane 4 by forming a triplet from planes 0, 1, 2; 2, 3, 5.
  // Then try to find a match on DUT (plane 3)
  //
  //
  DUTz = _planePosition[3];


  // Generate new triplet set with planes 0, 1, 2; 2,4,5:
  eff_triplets_UP = upstream_triplets;
  FindTriplets(hits, 2, 3, 5, eff_triplets_DOWN);

  for( std::vector<triplet>::iterator trip = eff_triplets_UP.begin(); trip != eff_triplets_UP.end(); trip++ ) {

    // Track impact position at Matching Point from Upstream:
    double xA = (*trip).getx_at(DUTz); // triplet impact point at matching position
    double yA = (*trip).gety_at(DUTz);

    // check if trip is isolated
    bool IsolatedTrip = IsTripletIsolated(trip, eff_triplets_UP, DUTz, _triCut*2.01);

    for( std::vector<triplet>::iterator drip = eff_triplets_DOWN.begin(); drip != eff_triplets_DOWN.end(); drip++ ){

      // Track impact position at Matching Point from Downstream:
      double xB = (*drip).getx_at(DUTz); // triplet impact point at matching position
      double yB = (*drip).gety_at(DUTz);

      // check if drip is isolated
      bool IsolatedDrip = IsTripletIsolated(drip, eff_triplets_DOWN, DUTz, _triCut*2.01);

      // driplet - triplet
      double dx = xB - xA; 
      double dy = yB - yA;

      // match driplet and triplet:
      if( abs(dx) > _cutx) continue;
      if( abs(dy) > _cuty) continue;

      //std::cout << " intersec ";

      // check isolation
      if( !IsolatedTrip || !IsolatedDrip ) continue;
      //std::cout << " , isolated ";

      // apply fiducial cut
      if ( fabs(xA) >  9.0) continue;
      if (     -yA  < -4.0) continue;

      //std::cout << " , fiducial " << std::endl;

      eff_triplets.push_back(*trip);


    } // Downstream
  } // Upstream

  //std::cout << " n eff triplets = " << eff_triplets.size() << std::endl;

  for( std::vector<triplet>::iterator trip = eff_triplets.begin(); trip != eff_triplets.end(); trip++ ) {

    double ddAMin = -1000.0;
    double xTrip = (*trip).getx_at(DUTz);
    double yTrip = (*trip).gety_at(DUTz);

    for( std::vector<hit>::iterator lhit = hits.begin(); lhit != hits.end(); lhit++ ){

      if( (*lhit).plane <= 3 || (*lhit).plane  > 4) continue; // want 4

      // Fill residuals of triplet and hit in the selected plane:
      if( (*lhit).plane == 4 ) {
	double xHit = (*lhit).x;
	double yHit = (*lhit).y;

	double ddA = sqrt( fabs(xHit - xTrip)*fabs(xHit - xTrip) 
	    + fabs(yHit - yTrip)*fabs(yHit - yTrip) );
	if(ddAMin < 0 || ddA < ddAMin) ddAMin = ddA;
      }
    } // end loop over hits

    // if distance is smaller then limit, accept this as matched Hit
    if(fabs(ddAMin) < eff_radius) {
      //n_matched_trips++;
      effix4->fill(-xTrip, 1.);
      effiy4->fill(-yTrip, 1.);
    } else {
      effix4->fill(-xTrip, 0.);
      effiy4->fill(-yTrip, 0.);
    }

  }
  
  eff_triplets_UP.clear();
  eff_triplets_DOWN.clear();
  eff_triplets.clear();

  //----------------------------------------------------------------------------
  // calculate efficiency of plane 5 by forming a triplet from planes 0, 1, 2; 2, 3, 4.
  // Then try to find a match on DUT (plane 3)
  //
  //
  DUTz = _planePosition[3];


  // Generate new triplet set with planes 0, 1, 2; 2,3,4:
  eff_triplets_UP = upstream_triplets;
  FindTriplets(hits, 2, 3, 4, eff_triplets_DOWN);

  for( std::vector<triplet>::iterator trip = eff_triplets_UP.begin(); trip != eff_triplets_UP.end(); trip++ ) {

    // Track impact position at Matching Point from Upstream:
    double xA = (*trip).getx_at(DUTz); // triplet impact point at matching position
    double yA = (*trip).gety_at(DUTz);

    // check if trip is isolated
    bool IsolatedTrip = IsTripletIsolated(trip, eff_triplets_UP, DUTz, _triCut*2.01);

    for( std::vector<triplet>::iterator drip = eff_triplets_DOWN.begin(); drip != eff_triplets_DOWN.end(); drip++ ){

      // Track impact position at Matching Point from Downstream:
      double xB = (*drip).getx_at(DUTz); // triplet impact point at matching position
      double yB = (*drip).gety_at(DUTz);

      // check if drip is isolated
      bool IsolatedDrip = IsTripletIsolated(drip, eff_triplets_DOWN, DUTz, _triCut*2.01);

      // driplet - triplet
      double dx = xB - xA; 
      double dy = yB - yA;

      // match driplet and triplet:
      if( abs(dx) > _cutx) continue;
      if( abs(dy) > _cuty) continue;

      //std::cout << " intersec ";

      // check isolation
      if( !IsolatedTrip || !IsolatedDrip ) continue;
      //std::cout << " , isolated ";

      // apply fiducial cut
      if ( fabs(xA) >  9.0) continue;
      if (     -yA  < -4.0) continue;

      //std::cout << " , fiducial " << std::endl;

      eff_triplets.push_back(*trip);


    } // Downstream
  } // Upstream

  //std::cout << " n eff triplets = " << eff_triplets.size() << std::endl;

  for( std::vector<triplet>::iterator trip = eff_triplets.begin(); trip != eff_triplets.end(); trip++ ) {

    double ddAMin = -1000.0;
    double xTrip = (*trip).getx_at(DUTz);
    double yTrip = (*trip).gety_at(DUTz);

    for( std::vector<hit>::iterator lhit = hits.begin(); lhit != hits.end(); lhit++ ){

      if( (*lhit).plane <= 4 || (*lhit).plane  > 5) continue; // want 5

      // Fill residuals of triplet and hit in the selected plane:
      if( (*lhit).plane == 5 ) {
	double xHit = (*lhit).x;
	double yHit = (*lhit).y;

	double ddA = sqrt( fabs(xHit - xTrip)*fabs(xHit - xTrip) 
	    + fabs(yHit - yTrip)*fabs(yHit - yTrip) );
	if(ddAMin < 0 || ddA < ddAMin) ddAMin = ddA;
      }
    } // end loop over hits

    // if distance is smaller then limit, accept this as matched Hit
    if(fabs(ddAMin) < eff_radius) {
      //n_matched_trips++;
      effix5->fill(-xTrip, 1.);
      effiy5->fill(-yTrip, 1.);
    } else {
      effix5->fill(-xTrip, 0.);
      effiy5->fill(-yTrip, 0.);
    }

  }
  
  eff_triplets_UP.clear();
  eff_triplets_DOWN.clear();
  eff_triplets.clear();


  //----------------------------------------------------------------------------
  // six: triplets A and driplets B
  // matching and GBL fit
  // kinks: triplets A vs driplets B
  // scattering point = DUT:

  //FIXME: DUTz is the place where they should be matched!
  // now hardcoded to center of telescope...
  DUTz = _planePosition[2] + (_planePosition[3] - _planePosition[2])/2;

  // Match the Telescope Upstream and Downstream Arm triplets to get tracks:
  std::vector<track> telescope_tracks;
  MatchTriplets(upstream_triplets,downstream_triplets,DUTz, telescope_tracks);

  //delete downstream_triplets;
  //delete upstream_triplets;

  int ngbl = 0;
  for( std::vector<track>::iterator tr = telescope_tracks.begin(); tr != telescope_tracks.end(); tr++ ){

    triplet trip = (*tr).get_upstream();
    triplet drip = (*tr).get_downstream();
    triplet srip((*tr).gethit(0), (*tr).gethit(2), (*tr).gethit(5)); // seed triplet -> srip
    if(_aluthickum > 1.) srip = trip;

    std::vector<double> xAplanes(_nTelPlanes);
    std::vector<double> yAplanes(_nTelPlanes);
    for (int i = 0; i < _nTelPlanes; i++){
      xAplanes[i] = srip.getx_at(_planePosition[i]);
      yAplanes[i] = srip.gety_at(_planePosition[i]);
    }

    // Track kinks as difference in triplet slopes:
    //      double kx = (*drip).slope().x - (*trip).slope().x; //kink
    //      double ky = (*drip).slope().y - (*trip).slope().y;
    double kx = (*tr).kink_x();
    double ky = (*tr).kink_y();

    // Track impact position at DUT from Downstream:
    double xB = drip.getx_at(DUTz);
    double yB = drip.gety_at(DUTz);

    // Track impact position at DUT from Upstream:
    double xA = srip.getx_at(DUTz);
    double yA = srip.gety_at(DUTz);

    double dx = xB - xA; // driplet - triplet
    double dy = yB - yA;

    // GBL with triplet A as seed:

    std::vector<gbl::GblPoint> traj_points;

    // build up trajectory:
    std::vector<double> sPoint;

    // plane 0:
    double s = 0;

    TMatrixD proL2m(2,2);
    proL2m.UnitMatrix();

    TVectorD meas(2);

    //double res = 3.42E-3; // [mm] Anemone telescope intrinsic resolution
    //res = 4.5E-3; // EUDET

    // scatter:
    TVectorD scat(2);
    scat.Zero(); //mean is zero

    double p = _eBeam; // beam momentum
    double epsSi = -1;
    double epsAir = -1.; // define later when dz is known
    double epsAlu = _aluthickum/1000./88.97; // Alu target

    double sumeps = 0.0;
    double tetSi = -1;
    double tetAir = -1.;


    TVectorD wscatSi(2);
    TVectorD wscatAir(2);

    TMatrixD alDer( 2, 3 ); // alignment derivatives
    alDer[0][0] = 1.0; // dx/dx
    alDer[1][0] = 0.0; // dy/dx
    alDer[0][1] = 0.0; // dx/dy
    alDer[1][1] = 1.0; // dy/dy

    std::vector<unsigned int> ilab; // 0-5 = telescope

    // plane 0-5:

    double rx[6];
    double ry[6];
    double trackhitx[6];
    double trackhity[6];
    double trackhitxloc[6];
    double trackhityloc[6];
    //double zprev = _planePosition[0];
    double step = 0.;
    s = -10.;

    // loop over all scatterers first to calculate sumeps
    for( int ipl = 0; ipl < 6; ++ipl ){
      epsSi = _thickness[ipl]/ 93.66 + 0.050 / 285.6; // Si + Kapton (Kapton ist 2 * 25  = 50 um thick, but has 1/3 rad length = 17 um)
      sumeps += epsSi;
      if( ipl < 5) {
	double distplane = _planePosition[ipl+1] - _planePosition[ipl];
	epsAir =   distplane  / 304200.; 
	sumeps += epsAir;
      }

    }
    sumeps += epsAlu;
   // done with calculating sum eps


    int DUT_label;
    int centre_label = -100;
    for( int ipl = 0; ipl < 6; ++ipl ){
      if(ipl == 0){
	// one dummy point:

	gbl::GblPoint * point = new gbl::GblPoint( JacobianPointToPoint( step ) );
	traj_points.push_back(*point);
	s += step;
	sPoint.push_back(s);
	step = 10; // [mm]
	delete point;
      }

      // Get the corresponding hit from the track:
      hit trackhit = (*tr).gethit(ipl);

      if (ipl == 0) clustersize0->fill(trackhit.clustersize);
      if (ipl == 1) clustersize1->fill(trackhit.clustersize);
      if (ipl == 2) clustersize2->fill(trackhit.clustersize);
      if (ipl == 3) clustersize3->fill(trackhit.clustersize);
      if (ipl == 4) clustersize4->fill(trackhit.clustersize);
      if (ipl == 5) clustersize5->fill(trackhit.clustersize);

      double dz = trackhit.z - srip.base().z;
      double xs = srip.base().x + srip.slope().x * dz; // Ax at plane
      double ys = srip.base().y + srip.slope().y * dz; // Ay at plane

      trackhitx[ipl] = trackhit.x;
      trackhity[ipl] = trackhit.y;
      trackhitxloc[ipl] = trackhit.locx;
      trackhityloc[ipl] = trackhit.locy;

      rx[ipl] = trackhit.x - xs;
      ry[ipl] = trackhit.y - ys;

      if( ipl == 0 ) {
	sixx0Histo->fill( -trackhit.x );
	sixy0Histo->fill( -trackhit.y );
      }
      if( ipl == 1 ) {
	sixx1Histo->fill( -trackhit.x );
	sixy1Histo->fill( -trackhit.y );
      }
      if( ipl == 2 ) {
	sixx2Histo->fill( -trackhit.x );
	sixy2Histo->fill( -trackhit.y );
      }
      if( ipl == 3 ) {
	sixx3Histo->fill( -trackhit.x );
	sixy3Histo->fill( -trackhit.y );
      }
      if( ipl == 4 ) {
	sixx4Histo->fill( -trackhit.x );
	sixy4Histo->fill( -trackhit.y );
      }
      if( ipl == 5 ) {
	sixx5Histo->fill( -trackhit.x );
	sixy5Histo->fill( -trackhit.y );
      }

      gbl::GblPoint * point = new gbl::GblPoint( JacobianPointToPoint( step ) );

      meas[0] = rx[ipl];
      meas[1] = ry[ipl];
      //meas[0] = trackhit.x;
      //meas[1] = trackhit.y;

      epsSi = _thickness[ipl]/93.66 + 0.050 / 286.6; // Si + Kapton (Kapton ist 2 * 25  = 50 um thick, but has 1/3 rad length = 17 um)
      tetSi = _kappa * 0.0136 * sqrt(epsSi) / p * ( 1 + 0.038*std::log(sumeps) );

      wscatSi[0] = 1.0 / ( tetSi * tetSi ); //weight
      wscatSi[1] = 1.0 / ( tetSi * tetSi );

      TVectorD measPrec(2); // precision = 1/resolution^2
      //measPrec[0] = 1.0 / _resolution[ipl] / _resolution[ipl];
      //measPrec[1] = 1.0 / _resolution[ipl] / _resolution[ipl];
      double _resolution_tmp = -1.0;
      if(trackhit.clustersize == 1) _resolution_tmp = _resolution[1];
      if(trackhit.clustersize == 2) _resolution_tmp = _resolution[2];
      if(trackhit.clustersize == 3) _resolution_tmp = _resolution[3];
      if(trackhit.clustersize == 4) _resolution_tmp = _resolution[4];
      if(trackhit.clustersize == 5) _resolution_tmp = _resolution[5];
      if(trackhit.clustersize == 6) _resolution_tmp = _resolution[6];
      if(trackhit.clustersize >  6) _resolution_tmp = _resolution[7];
      //_resolution_tmp = 3.24*1e-3;

      measPrec[0] = 1.0 / _resolution_tmp / _resolution_tmp;
      measPrec[1] = 1.0 / _resolution_tmp / _resolution_tmp;

      if(ipl != _dut_plane) { point->addMeasurement( proL2m, meas, measPrec ); }

      point->addScatterer( scat, wscatSi );

      std::vector<int> globalLabels(3);
      globalLabels[0] = 10 + ipl; // dx
      globalLabels[1] = 20 + ipl; // dy
      globalLabels[2] = 40 + ipl; // rot
      alDer[0][2] = -ys; // dx/rot
      alDer[1][2] =  xs; // dy/rot

      traj_points.push_back(*point);
      s += step;
      sPoint.push_back( s );
      DUT_label = sPoint.size();
      ilab.push_back(DUT_label);
      delete point;




      if( ipl < 5) {
	double distplane = _planePosition[ipl+1] - _planePosition[ipl];

	step = 0.21*distplane;
	epsAir =   0.5*distplane  / 304200.; 
	tetAir = _kappa * 0.0136 * sqrt(epsAir) / p * ( 1 + 0.038*std::log(sumeps) ); // add 10% for trial

	wscatAir[0] = 1.0 / ( tetAir * tetAir ); // weight
	wscatAir[1] = 1.0 / ( tetAir * tetAir ); 

	gbl::GblPoint * point = new gbl::GblPoint( JacobianPointToPoint( step ) );
	point->addScatterer( scat, wscatAir );

	s += step;
	traj_points.push_back(*point);
	sPoint.push_back( s );
	delete point;

	step = 0.58*distplane;
	if(ipl ==2){ // insert point at centre
	  step = step / 2.;
	  gbl::GblPoint * pointcentre = new gbl::GblPoint( JacobianPointToPoint( step ) );

	  if(_aluthickum > 1.) { // at least 1 um target
	    double tetAlu = _kappa*0.0136 * sqrt(epsAlu) / p * ( 1 + 0.038*std::log(sumeps) );

	    TVectorD wscatAlu(2);
	    wscatAlu[0] = 1.0 / ( tetAlu * tetAlu ); // weight
	    wscatAlu[1] = 1.0 / ( tetAlu * tetAlu ); 

	    pointcentre->addScatterer( scat, wscatAlu );
	  }
	  s += step;
	  traj_points.push_back(*pointcentre);
	  sPoint.push_back( s );
          centre_label = sPoint.size();
	  delete pointcentre;
	  // now again step/2
        }

	gbl::GblPoint * point1 = new gbl::GblPoint( JacobianPointToPoint( step ) );
	point1->addScatterer( scat, wscatAir );

	s += step;
	traj_points.push_back(*point1);
	sPoint.push_back( s );
	delete point1;

	step = 0.21*distplane; // remaing distance to next plane
      }
      else{

	// one more dummy point:

	step = 10; // [mm]
	gbl::GblPoint * point = new gbl::GblPoint( JacobianPointToPoint( step ) );
	traj_points.push_back(*point);
	s += step;
	sPoint.push_back(s);
	delete point;
      }

    } // loop over planes

    // monitor what we put into GBL:
    selxHisto->fill( -xA ); // triplet at DUT
    selyHisto->fill( -yA );
    selaxHisto->fill( srip.slope().x*1E3 );
    selayHisto->fill( srip.slope().y*1E3 );
    seldxHisto->fill( dx*1E3 ); // triplet-driplet match
    seldyHisto->fill( dy*1E3 );
    selkxHisto->fill( kx*1E3 ); // triplet-driplet kink
    selkyHisto->fill( ky*1E3 );

    seldx1Histo->fill( rx[1]*1E3 ); // triplet interpol
    seldy1Histo->fill( ry[1]*1E3 );
    seldx3Histo->fill( rx[3]*1E3 ); // triplet extrapol
    seldy3Histo->fill( ry[3]*1E3 );
    seldx4Histo->fill( rx[4]*1E3 );
    seldy4Histo->fill( ry[4]*1E3 );
    seldx5Histo->fill( rx[5]*1E3 );
    seldy5Histo->fill( ry[5]*1E3 );

    double Chi2;
    int Ndf;
    double lostWeight;

    gbl::GblTrajectory traj(traj_points, false ); // curvature = false
    std::string fit_optionList = "";
    traj.fit( Chi2, Ndf, lostWeight, fit_optionList );
    //traj.getLabels(ilab); // instead pushback sPoint.size() when adding plane

    // debug:

    
    if(_nEvt < 100){
      cout << "traj with " << traj.getNumPoints() << " points:" << endl;
      for( int ipl = 0; ipl < 6; ++ipl ){
        cout << "  plane " << ipl << ", lab " << ilab[ipl];
        cout << "  z " << sPoint[ilab[ipl]-1];
        cout << " dx " << rx[ipl];
        cout << " dy " << ry[ipl];
        cout << " chi2 " << Chi2;
        cout << " ndf " << Ndf;
        cout << endl;
      }

       std::cout << " Is traj valid? " << traj.isValid() << std::endl;
       traj.printPoints();
      //traj.printTrajectory();
      //traj.printData();
    }

    ngbl++;

    gblndfHisto->fill( Ndf );
    if( Ndf == 8 ) 
      gblchi2aHisto->fill( Chi2 );
    else
      gblchi2bHisto->fill( Chi2 );

    double probchi = 0;
    probchi = TMath::Prob( Chi2, Ndf );
    gblprbHisto->fill( probchi );

    gblprbxHisto->fill(-xA, probchi);
    gblprbyHisto->fill(-yA, probchi);

    // bad fits:

    if( probchi < 0.01 ) {

      badxHisto->fill( -xA ); // triplet at DUT
      badyHisto->fill( -yA );
      badaxHisto->fill( srip.slope().x*1E3 );
      badayHisto->fill( srip.slope().y*1E3 );
      baddxHisto->fill( dx*1E3 ); // triplet-driplet match
      baddyHisto->fill( dy*1E3 );
      badkxHisto->fill( kx*1E3 ); // triplet-driplet kink
      badkyHisto->fill( ky*1E3 );

      baddx1Histo->fill( rx[1]*1E3 ); // triplet interpol
      baddy1Histo->fill( ry[1]*1E3 );
      baddx3Histo->fill( rx[3]*1E3 ); // triplet extrapol
      baddy3Histo->fill( ry[3]*1E3 );
      baddx4Histo->fill( rx[4]*1E3 );
      baddy4Histo->fill( ry[4]*1E3 );
      baddx5Histo->fill( rx[5]*1E3 );
      baddy5Histo->fill( ry[5]*1E3 );
    }// bad fit
    else {

      goodxHisto->fill( -xA ); // triplet at DUT
      goodyHisto->fill( -yA );
      goodx1Histo->fill( rx[1]*1E3 ); // triplet interpol
      goody1Histo->fill( ry[1]*1E3 );
    } // OK fit



    // look at good fit:
    //
    //double chi2_cut = 0.1;

    if( probchi > _probchi2_cut){

      TVectorD aCorrection(5);
      TMatrixDSym aCovariance(5);

      double ax[8];
      // double ay[8];
      unsigned int k = 0;


      unsigned int ndim = 2;
      TVectorD aResiduals(ndim);
      TVectorD aMeasErrors(ndim);
      TVectorD aResErrors(ndim);
      TVectorD aDownWeights(ndim);


      TVectorD aKinks(ndim);
      TVectorD aKinkErrors(ndim);
      TVectorD kResErrors(ndim);
      TVectorD kDownWeights(ndim);

      double pixel_size = 18.4e-3;

      unsigned int ndata = 2;
      //track = q/p, x', y', x, y
      //        0,   1,  2,  3, 4
      // at plane DUT:
      unsigned int ipos = ilab[0];
      traj.getResults( ipos, aCorrection, aCovariance );
      traj.getMeasResults( ipos, ndata, aResiduals, aMeasErrors, aResErrors, aDownWeights );
      traj.getScatResults( ipos, ndata, aKinks, aKinkErrors, kResErrors, kDownWeights );
      aResiduals[0] = rx[0] - aCorrection[3];
      aResiduals[1] = ry[0] - aCorrection[4];
      gblax0Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
      gbldx0Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
      gbldx01Histo->fill( aCorrection[3] ); // shift x [mm]
      gblrx0Histo->fill( ( rx[0] - aCorrection[3] ) * 1E3 ); // residual x [um]
      gblry0Histo->fill( ( ry[0] - aCorrection[4] ) * 1E3 ); // residual y [um]
      gblpx0Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
      gblpy0Histo->fill( aResiduals[1] / aResErrors[1] ); // pull
      if(_dut_plane == 0) gblpx0_unbHisto->fill( (rx[0] - aCorrection[3]) / sqrt(_resolution[0]*_resolution[0] + aCovariance(3,3)) ); // unbiased pull
      if(_dut_plane == 0) gblpy0_unbHisto->fill( (ry[0] - aCorrection[4]) / sqrt(_resolution[0]*_resolution[0] + aCovariance(4,4)) ); // unbiased pull
      gblqx0Histo->fill( aKinks[0]*1E3 ); // kink RESIDUAL (measured kink - fit kink)
      ax[k] = aCorrection[1]; // angle correction at plane, for kinks
      // TProfile for res_x a.f.o. x
      gblrxvsx0->fill( xAplanes.at(k), sqrt(TMath::Pi()/2.)*fabs(aResiduals[0]));
      gblryvsy0->fill( yAplanes.at(k), sqrt(TMath::Pi()/2.)*fabs(aResiduals[1]));
      gblrxvsx01->fill((trackhitx[0] - aResiduals[0]), sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); // seed corrected
      gblryvsy01->fill((trackhity[0] - aResiduals[1]), sqrt(TMath::Pi()/2.)*fabs(aResiduals[1]));

      double corrPos[2];
      corrPos[0] = trackhitx[0] - aResiduals[0]; //+ .4e-3;  should be zero ! :/
      corrPos[1] = trackhity[0] - aResiduals[1]; // + .15e-3; should be zero (no alignment on plane 0, not even pre-align)      
      int nx = (corrPos[0]) / pixel_size;
      int ny = (corrPos[1]) / pixel_size;
      int invsignx = -(corrPos[0]) / fabs((corrPos[0]));
      int invsigny = -(corrPos[1]) / fabs((corrPos[1]));

      // do again for extrapolated srip position, userd ONLY for gblrxvsxpix0 and gblryvsypix0
      int nx1 = (xAplanes[0]) / pixel_size;
      int ny1 = (yAplanes[0]) / pixel_size;
      int invsignx1 = -(xAplanes[0]) / fabs((xAplanes[0]));
      int invsigny1 = -(yAplanes[0]) / fabs((yAplanes[0]));

      gblrxvsxpix0->fill(                    (xAplanes[0] +invsignx1*(abs(nx1) +1.)*pixel_size)*1e3, sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); 
      gblryvsypix0->fill(                    (yAplanes[0] +invsigny1*(abs(ny1) +1.)*pixel_size)*1e3, sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])); 
      gblrxvsxpix01->fill( ((trackhitx[0] - aResiduals[0])+invsignx*(abs(nx) +1.)*pixel_size)*1e3, sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); 
      gblryvsypix01->fill( ((trackhity[0] - aResiduals[1])+invsigny*(abs(ny) +1.)*pixel_size)*1e3, sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])); 
      // clustersize-specific plots
      int CS0 = (*tr).gethit(0).clustersize;

      double modPos[2];
      modPos[0] = ((corrPos[0])+(invsignx*(abs(nx) +.5) + 0.5) *pixel_size)*1e3;
      modPos[1] = ((corrPos[1])+(invsigny*(abs(ny) +.5) + 0.5) *pixel_size)*1e3;
      if (CS0 == 1){
        gblrxvsxpix01_CS1->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); 
        gblryvsypix01_CS1->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])); 
      } 
      if (CS0 == 2){
        gblrxvsxpix01_CS2->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); 
        gblryvsypix01_CS2->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])); 
      } 
      if (CS0 == 3){
        gblrxvsxpix01_CS3->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); 
        gblryvsypix01_CS3->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])); 
      } 
      if (CS0 == 4){
        gblrxvsxpix01_CS4->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); 
        gblryvsypix01_CS4->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])); 
      } 


      k++;

      ipos = ilab[1];
      traj.getResults( ipos, aCorrection, aCovariance );
      traj.getMeasResults(static_cast<unsigned int>(ipos), ndata, aResiduals, aMeasErrors, aResErrors, aDownWeights );
      traj.getScatResults(static_cast<unsigned int>(ipos), ndata, aKinks, aKinkErrors, kResErrors, kDownWeights );
      aResiduals[0] = rx[1] - aCorrection[3];
      aResiduals[1] = ry[1] - aCorrection[4];
      gblax1Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
      gbldx1Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
      gbldx11Histo->fill( aCorrection[3] ); // shift x [mm]
      gblrx1Histo->fill( ( rx[1] - aCorrection[3] ) * 1E3 ); // residual x [um]
      gblry1Histo->fill( ( ry[1] - aCorrection[4] ) * 1E3 ); // residual y [um]
      gblpx1Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
      gblpy1Histo->fill( aResiduals[1] / aResErrors[1] ); // pull
      if(_dut_plane == 1) gblpx1_unbHisto->fill( (rx[1] - aCorrection[3]) / sqrt(_resolution[0]*_resolution[0] + aCovariance(3,3)) ); // unbiased pull
      if(_dut_plane == 1) gblpy1_unbHisto->fill( (ry[1] - aCorrection[4]) / sqrt(_resolution[0]*_resolution[0] + aCovariance(4,4)) ); // unbiased pull
      gblqx1Histo->fill( aKinks[0]*1E3 ); // kink-residual (NOT KINK itself!)
      gblsx1Histo->fill( aKinks[0]/aKinkErrors[0] ); // x kink pull
      gbltx1Histo->fill( aKinks[0]/kResErrors[0] ); // x kink pull
      ax[k] = aCorrection[1]; // angle correction at plane, for kinks
      //ay[k] = aCorrection[2]; // angle correction at plane, for kinks
      //std::cout << " aResiduals[0] = " << aResiduals[0] << " aResErrors[0] = " << aResErrors[0] << std::endl;  
      k++;

      ipos = ilab[2];
      traj.getResults( ipos, aCorrection, aCovariance );
      traj.getMeasResults(static_cast<unsigned int>(ipos), ndata, aResiduals, aMeasErrors, aResErrors, aDownWeights );
      traj.getScatResults(static_cast<unsigned int>(ipos), ndata, aKinks, aKinkErrors, kResErrors, kDownWeights );
      aResiduals[0] = rx[2] - aCorrection[3];
      aResiduals[1] = ry[2] - aCorrection[4];
      gblax2Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
      gbldx2Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
      gbldx21Histo->fill( aCorrection[3] ); // shift x [mm]
      gblrx2Histo->fill( ( rx[2] - aCorrection[3] ) * 1E3 ); // residual x [um]
      gblry2Histo->fill( ( ry[2] - aCorrection[4] ) * 1E3 ); // residual y [um]
      gblpx2Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
      gblpy2Histo->fill( aResiduals[1] / aResErrors[1] ); // pull
      if(_dut_plane == 2) gblpx2_unbHisto->fill( (rx[2] - aCorrection[3]) / sqrt(_resolution[0]*_resolution[0] + aCovariance(3,3)) ); // unbiased pull
      if(_dut_plane == 2) gblpy2_unbHisto->fill( (ry[2] - aCorrection[4]) / sqrt(_resolution[0]*_resolution[0] + aCovariance(4,4)) ); // unbiased pull
      gblqx2Histo->fill( aKinks[0]*1E3 ); // kink-residual
      gblsx2Histo->fill( aKinks[0]/aKinkErrors[0] ); // x kink res over kinkError
      gbltx2Histo->fill( aKinks[0]/kResErrors[0] ); // x kink pull
      ax[k] = aCorrection[1]; // angle correction at plane, for kinks
      //ay[k] = aCorrection[2]; // angle correction at plane, for kinks
      k++;

      ipos = ilab[3];
      traj.getResults( ipos, aCorrection, aCovariance );
      traj.getMeasResults(static_cast<unsigned int>(ipos), ndata, aResiduals, aMeasErrors, aResErrors, aDownWeights );
      traj.getScatResults(static_cast<unsigned int>(ipos), ndata, aKinks, aKinkErrors, kResErrors, kDownWeights );
      aResiduals[0] = rx[3] - aCorrection[3];
      aResiduals[1] = ry[3] - aCorrection[4];
      gblax3Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
      gbldx3Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
      gbldx31Histo->fill( aCorrection[3] ); // shift x [mm]
      gblrx3Histo->fill( ( rx[3] - aCorrection[3] ) * 1E3 ); // residual x [um]
      gblry3Histo->fill( ( ry[3] - aCorrection[4] ) * 1E3 ); // residual y [um]
      gblpx3Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
      gblpy3Histo->fill( aResiduals[1] / aResErrors[1] ); // pull
      gblqx3Histo->fill( aKinks[0]*1E3 ); // kink
      gblsx3Histo->fill( aKinks[0]/aKinkErrors[0] ); // x kink pull
      gbltx3Histo->fill( aKinks[0]/kResErrors[0] ); // x kink pull
      ax[k] = aCorrection[1]; // angle correction at plane, for kinks
      //ay[k] = aCorrection[2]; // angle correction at plane, for kinks
      //
      if(_dut_plane == 3) gblpx3_unbHisto->fill( (rx[3] - aCorrection[3]) / sqrt(_resolution[0]*_resolution[0] + aCovariance(3,3)) ); // unbiased pull
      if(_dut_plane == 3) gblpy3_unbHisto->fill( (ry[3] - aCorrection[4]) / sqrt(_resolution[0]*_resolution[0] + aCovariance(4,4)) ); // unbiased pull

      // clustersize-specific plots
      int CS3 = (*tr).gethit(3).clustersize;
      //corrPos[0] = trackhitx[3] - aResiduals[0] +  7.077e-3 + 0.35e-3; // run000117, get automatically 
      //corrPos[1] = trackhity[3] - aResiduals[1] + 18.128e-3 + 0.1e-3;  // run000117, get automatically      
      corrPos[0] = trackhitx[3] - aResiduals[0] +  6.645e-3; // run002140
      corrPos[1] = trackhity[3] - aResiduals[1] + 10.324e-3;  // run002140

      // correct for rotation
      double corrPos2[2];
      //double gamma3 = 5.3e-3; // run000117, get automatically
      double gamma3 = 1.94969e-03; // run002140
      corrPos2[0] = cos(gamma3)*corrPos[0] - sin(gamma3)*corrPos[1];
      corrPos2[1] = sin(gamma3)*corrPos[0] + cos(gamma3)*corrPos[1];
      
      corrPos[0] = corrPos2[0];
      corrPos[1] = corrPos2[1];

      nx = fabs((corrPos[0])) / pixel_size;
      ny = fabs((corrPos[1])) / pixel_size;
      invsignx = -(corrPos[0]) / fabs((corrPos[0]));
      invsigny = -(corrPos[1]) / fabs((corrPos[1]));

      modPos[0] = ((corrPos[0])+(invsignx*(abs(nx) +.5) + 0.5)*pixel_size)*1e3;
      modPos[1] = ((corrPos[1])+(invsigny*(abs(ny) +.5) + 0.5)*pixel_size)*1e3;

      //gblrxvsxpix31->fill( (trackhitx[3] - aResiduals[0])+invsignx*(abs(nx) +0.5)*pixel_size, sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); 
      //gblryvsypix31->fill( (trackhity[3] - aResiduals[1])+invsigny*(abs(ny) +0.5)*pixel_size, sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])); 
      gblrxvsxpix31->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); 
      gblryvsypix31->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])); 

      // overlay of all CSs
      gblnxy->fill(modPos[0], modPos[1], CS3 );

      if( CS3 == 1 ) {
	gblrx3_cs1Histo->fill( aResiduals[0]* 1E3 ); 
	gblry3_cs1Histo->fill( aResiduals[1]* 1E3 ); 
	gblpx3_cs1Histo->fill( aResiduals[0] / aResErrors[0] ); 
	gblpy3_cs1Histo->fill( aResiduals[1] / aResErrors[1] ); 

	gblnCS1xy->fill(modPos[0], modPos[1] );
	gblnCS1xy1->fill((trackhitx[3] - aResiduals[0]),
	    (trackhity[3] - aResiduals[1]));

	gblrxvsxpix3cs1->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); 
	gblryvsypix3cs1->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])); 
      }
      if( CS3 == 2 ) {
	gblrx3_cs2Histo->fill( aResiduals[0]* 1E3 ); 
	gblry3_cs2Histo->fill( aResiduals[1]* 1E3 ); 
	gblpx3_cs2Histo->fill( aResiduals[0] / aResErrors[0] ); 
	gblpy3_cs2Histo->fill( aResiduals[1] / aResErrors[1] ); 

	gblnCS2xy->fill(modPos[0], modPos[1] );
	gblnCS2xy1->fill((trackhitx[3] - aResiduals[0]),
	    (trackhity[3] - aResiduals[1]));

	gblrxvsxpix3cs2->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); 
	gblryvsypix3cs2->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])); 
      }
      if( CS3 == 3 ) {
	gblrx3_cs3Histo->fill( aResiduals[0]* 1E3 ); 
	gblry3_cs3Histo->fill( aResiduals[1]* 1E3 ); 
	gblpx3_cs3Histo->fill( aResiduals[0] / aResErrors[0] ); 
	gblpy3_cs3Histo->fill( aResiduals[1] / aResErrors[1] ); 

	gblnCS3xy->fill(modPos[0], modPos[1] );
	gblnCS3xy1->fill((trackhitx[3] - aResiduals[0]),
	    (trackhity[3] - aResiduals[1]));

	gblrxvsxpix3cs3->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); 
	gblryvsypix3cs3->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])); 
      }
      if( CS3 == 4 ) {
	gblrx3_cs4Histo->fill( aResiduals[0]* 1E3 ); 
	gblry3_cs4Histo->fill( aResiduals[1]* 1E3 ); 
	gblpx3_cs4Histo->fill( aResiduals[0] / aResErrors[0] ); 
	gblpy3_cs4Histo->fill( aResiduals[1] / aResErrors[1] ); 

	gblnCS4xy->fill(modPos[0], modPos[1] );
	gblnCS4xy1->fill((trackhitx[3] - aResiduals[0]),
	    (trackhity[3] - aResiduals[1]));

	gblrxvsxpix3cs4->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); 
	gblryvsypix3cs4->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])); 
      }
      if( CS3 == 5 ) {
	gblrx3_cs5Histo->fill( aResiduals[0]* 1E3 ); 
	gblry3_cs5Histo->fill( aResiduals[1]* 1E3 ); 
	gblpx3_cs5Histo->fill( aResiduals[0] / aResErrors[0] ); 
	gblpy3_cs5Histo->fill( aResiduals[1] / aResErrors[1] ); 

	gblnCS5xy->fill(modPos[0], modPos[1] );
	gblnCS5xy1->fill((trackhitx[3] - aResiduals[0]),
	    (trackhity[3] - aResiduals[1]));

	gblrxvsxpix3cs5->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); 
	gblryvsypix3cs5->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])); 
      }
      if( CS3 == 6 ) {
	gblrx3_cs6Histo->fill( aResiduals[0]* 1E3 ); 
	gblry3_cs6Histo->fill( aResiduals[1]* 1E3 ); 
	gblpx3_cs6Histo->fill( aResiduals[0] / aResErrors[0] ); 
	gblpy3_cs6Histo->fill( aResiduals[1] / aResErrors[1] ); 

	gblnCS6xy->fill(modPos[0], modPos[1] );
	gblnCS6xy1->fill((trackhitx[3] - aResiduals[0]),
	    (trackhity[3] - aResiduals[1]));

	gblrxvsxpix3cs6->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); 
	gblryvsypix3cs6->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])); 
      }
      if( CS3 > 6 ) {
	gblpx3_cs7Histo->fill( aResiduals[0] / aResErrors[0] ); 
	gblpy3_cs7Histo->fill( aResiduals[1] / aResErrors[1] ); 

	//gblnCS6xy->fill(modPos[0], modPos[1] );

	//gblrxvsxpix3cs6->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); 
	//gblryvsypix3cs6->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])); 
      }

      kinkpixvsxy->fill( modPos[0], modPos[1] , sqrt((aCorrection[1]*aCorrection[1] + aCorrection[2]*aCorrection[2]))*1E3 ); //sqrt(<kink^2>) [mrad]


      k++;

      ipos = ilab[4];
      traj.getResults( ipos, aCorrection, aCovariance );
      traj.getMeasResults(static_cast<unsigned int>(ipos), ndata, aResiduals, aMeasErrors, aResErrors, aDownWeights );
      traj.getScatResults(static_cast<unsigned int>(ipos), ndata, aKinks, aKinkErrors, kResErrors, kDownWeights );
      aResiduals[0] = rx[4] - aCorrection[3];
      aResiduals[1] = ry[4] - aCorrection[4];
      gblax4Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
      gbldx4Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
      gbldx41Histo->fill( aCorrection[3] ); // shift x [mm]
      gblrx4Histo->fill( ( rx[4] - aCorrection[3] ) * 1E3 ); // residual x [um]
      gblry4Histo->fill( ( ry[4] - aCorrection[4] ) * 1E3 ); // residual y [um]
      gblpx4Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
      gblpy4Histo->fill( aResiduals[1] / aResErrors[1] ); // pull
      if(_dut_plane == 4) gblpx4_unbHisto->fill( (rx[4] - aCorrection[3]) / sqrt(_resolution[0]*_resolution[0] + aCovariance(3,3)) ); // unbiased pull
      if(_dut_plane == 4) gblpy4_unbHisto->fill( (ry[4] - aCorrection[4]) / sqrt(_resolution[0]*_resolution[0] + aCovariance(4,4)) ); // unbiased pull
      gblqx4Histo->fill( aKinks[0]*1E3 ); // kink
      gblsx4Histo->fill( aKinks[0]/aKinkErrors[0] ); // x kink pull
      gbltx4Histo->fill( aKinks[0]/kResErrors[0] ); // x kink pull
      ax[k] = aCorrection[1]; // angle correction at plane, for kinks
      //ay[k] = aCorrection[2]; // angle correction at plane, for kinks
      k++;

      ipos = ilab[5];
      traj.getResults( ipos, aCorrection, aCovariance );
      traj.getMeasResults(static_cast<unsigned int>(ipos), ndata, aResiduals, aMeasErrors, aResErrors, aDownWeights );
      traj.getScatResults(static_cast<unsigned int>(ipos), ndata, aKinks, aKinkErrors, kResErrors, kDownWeights );
      aResiduals[0] = rx[5] - aCorrection[3];
      aResiduals[1] = ry[5] - aCorrection[4];
      gblax5Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
      gbldx5Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
      gbldx51Histo->fill( aCorrection[3] ); // shift x [mm]
      gblrx5Histo->fill( ( rx[5] - aCorrection[3] ) * 1E3 ); // residual x [um]
      gblry5Histo->fill( ( ry[5] - aCorrection[4] ) * 1E3 ); // residual y [um]
      gblpx5Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
      gblpy5Histo->fill( aResiduals[1] / aResErrors[1] ); // pull
      if(_dut_plane == 5) gblpx5_unbHisto->fill( (rx[5] - aCorrection[3]) / sqrt(_resolution[0]*_resolution[0] + aCovariance(3,3)) ); // unbiased pull
      if(_dut_plane == 5) gblpy5_unbHisto->fill( (ry[5] - aCorrection[4]) / sqrt(_resolution[0]*_resolution[0] + aCovariance(4,4)) ); // unbiased pull
      gblqx5Histo->fill( aKinks[0]*1E3 ); // kink
      ax[k] = aCorrection[1]; // angle correction at plane, for kinks
      //ay[k] = aCorrection[2]; // angle correction at plane, for kinks
      //

      corrPos[0] = trackhitx[5] - aResiduals[0] + 2.028e-3 ;  // alignment const (from pre-align only for plane 5) modulo 18.4
      corrPos[1] = trackhity[5] - aResiduals[1] + 1.930e-3;       
      
      nx = (corrPos[0]) / pixel_size;
      ny = (corrPos[1]) / pixel_size;
      invsignx = -(corrPos[0]) / fabs((corrPos[0]));
      invsigny = -(corrPos[1]) / fabs((corrPos[1]));

      modPos[0] = ((corrPos[0])+(invsignx*(abs(nx) +.5) + 0.5) *pixel_size)*1e3;
      modPos[1] = ((corrPos[1])+(invsigny*(abs(ny) +.5) + 0.5) *pixel_size)*1e3;

      // clustersize-specific plots
      int CS5 = (*tr).gethit(5).clustersize;

      if (CS5 == 1){
        gblrxvsxpix51_CS1->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); 
        gblryvsypix51_CS1->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])); 
      } 
      if (CS5 == 2){
        gblrxvsxpix51_CS2->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); 
        gblryvsypix51_CS2->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])); 
      } 
      if (CS5 == 3){
        gblrxvsxpix51_CS3->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); 
        gblryvsypix51_CS3->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])); 
      } 
      if (CS5 == 4){
        gblrxvsxpix51_CS4->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); 
        gblryvsypix51_CS4->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])); 
      } 



      k++;
      // do some more analysis with kinks. 
      //   Calculate the two angles as corr_n-1 - corr_n-2, corr_n - corr_n-1, corr_n+1 - corr_n -> sum of all =  corr_n+1 - corr_n-2

      ipos = centre_label-2;
      traj.getResults( ipos, aCorrection, aCovariance );
      double kink_upstream = aCorrection[1];

      ipos = centre_label+1; // 1 downstream of centre
      traj.getResults( ipos, aCorrection, aCovariance );
      double kink_downstream = aCorrection[1];

      gblkxCentreHisto->fill( (kink_downstream - kink_upstream)*1E3 ); // kink at air/alu (sum of neighbours) [mrad]
      gblkxCentre1Histo->fill((kink_downstream - kink_upstream) ); // kink at air/alu (sum of neighbours) [rad]
      gblkx1Histo->fill( (ax[1] - ax[0])*1E3 ); // kink at 1 [mrad]
      gblkx2Histo->fill( (ax[2] - ax[1])*1E3 ); // kink at 2 [mrad]
      gblkx3Histo->fill( (ax[3] - ax[2])*1E3 ); // kink at 3 [mrad]
      gblkx4Histo->fill( (ax[4] - ax[3])*1E3 ); // kink at 4 [mrad]
      gblkx5Histo->fill( (ax[5] - ax[4])*1E3 ); // kink at 5 [mrad] // THIS IS NULL as ax[5] is not updated
      //gblkx6Histo->fill( (ax[6] - ax[5])*1E3 ); // kink at 6 [mrad]

    } // end if good fit 


    //------------------------------------------------------------------------
    // intersect point in z:
    hit intersect = (*tr).intersect();
    double zx = intersect.x;
    double zy = intersect.y;

    if( abs(dy) < 0.1 ) { // no cut on dx
      if( abs( kx ) > 0.003 ) { // 
	sixzx3Histo->fill( zx - _planePosition[2] );
      }
      if( abs( kx ) > 0.002 ) { // 
	sixzx2Histo->fill( zx - _planePosition[2] );
      }
    }

    if( abs(dx) < 0.1 ) { // no cut on dy
      if( abs( ky ) > 0.003 ) { // 
	sixzy3Histo->fill( zy - _planePosition[2] );
      }
      if( abs( ky ) > 0.002 ) { // 
	sixzy2Histo->fill( zy - _planePosition[2] );
      }
    }

    //------------------------------------------------------------------------
    // z intersect:

    if( abs(dx) < 0.2 && abs(dy) < 0.2 ) { // looser cut allows more z range

      // measure scattering angle x after cuts in y:
      // cut on ky creates bias in kx

      if( abs( kx ) > 0.001 ) {
	sixzx1Histo->fill( zx - _planePosition[2] );
	if( abs( zx - DUTz ) < 30 ) {
	  sixkyzxHisto->fill( ky*1E3 );
	  sixkxzxHisto->fill( kx*1E3 ); // plot with gap, fittp0g.C("sixkxzx")
	}
      }

      if( abs( ky ) > 0.001 ) {
	sixzy1Histo->fill( zy - _planePosition[2] );
	if( abs( zy - DUTz ) < 30 ) {
	  sixkxzyHisto->fill( kx*1E3 );
	  sixkyzyHisto->fill( ky*1E3 ); // plot with gap
	}
      }

    }//match

    //------------------------------------------------------------------------

  } // Loop over found tracks

  nsixHisto->fill( telescope_tracks.size() );

  //prevdutrefddt = dutrefddt;
  // Clear memory
  //delete telescope_tracks;
}


//------------------------------------------------------------------------------
void EUTelTripletGBL::check( LCEvent * /* evt */  ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


//------------------------------------------------------------------------------
void EUTelTripletGBL::end(){

  // Print the summary:
  streamlog_out(MESSAGE5)
    << "---------------------------------------------------------------------------------------------------------" << std::endl
    << std::endl
    << "Processed events:    "
    << std::setw(10) << std::setiosflags(std::ios::right)
    << _nEvt << std::resetiosflags(std::ios::right) << std::endl;

} // end end


TMatrixD EUTelTripletGBL::JacobianPointToPoint( double ds ) {
  /* for GBL:
     Jacobian for straight line track
     track = q/p, x', y', x, y
     0,   1,  2,  3, 4
     */
  TMatrixD jac( 5, 5 );
  jac.UnitMatrix();
  jac[3][1] = ds; // x = x0 + xp * ds
  jac[4][2] = ds; // y = y0 + yp * ds
  return jac;
}


void EUTelTripletGBL::TelescopeCorrelationPlots(std::vector<hit> &telescopehits) {
  for( std::vector<hit>::iterator ihit = telescopehits.begin(); ihit != telescopehits.end(); ihit++ ){

    int ipl = (*ihit).plane;

    for( std::vector<hit>::iterator jhit = telescopehits.begin(); jhit != telescopehits.end(); jhit++ ){

      int jpl = (*jhit).plane;
      double dx = (*jhit).x - (*ihit).x;
      double dy = (*jhit).y - (*ihit).y;

      if( ipl == 0 ) {
	if( jpl == 1 ) {
	  dx01Histo->fill( dx );
	  dy01Histo->fill( dy );
	  if( abs(dy) < 1 ) du01Histo->fill( dx );
	}
	if( jpl == 2 ) {
	  dx02Histo->fill( dx );
	}
	if( jpl == 3 ) {
	  dx03Histo->fill( dx );
	}
	if( jpl == 4 ) {
	  dx04Histo->fill( dx );
	}
	if( jpl == 5 ) {
	  dx05Histo->fill( dx );
	}

      }//ipl==0

      if( ipl == 1 ) {
	if( jpl == 2 ) {
	  dx12Histo->fill( dx );
	  dy12Histo->fill( dy );
	  if( abs(dy) < 1 ) du12Histo->fill( dx );
	}
      }

      if( ipl == 2 ) {
	if( jpl == 3 ) {
	  dx23Histo->fill( dx );
	  dy23Histo->fill( dy );
	  if( abs(dy) < 1 ) du23Histo->fill( dx );
	}
      }

      if( ipl == 3 ) {
	if( jpl == 4 ) {
	  dx34Histo->fill( dx );
	  dy34Histo->fill( dy );
	  if( abs(dy) < 1 ) du34Histo->fill( dx );
	}
      }

      if( ipl == 4 ) {
	if( jpl == 5 ) {
	  dx45Histo->fill( dx );
	  dy45Histo->fill( dy );
	  if( abs(dy) < 1 ) du45Histo->fill( dx );
	}
      }
    }
  }
}


void EUTelTripletGBL::FindTriplets(std::vector<hit> &hits, unsigned int plane0, unsigned int plane1, unsigned int plane2, std::vector<EUTelTripletGBL::triplet> &triplets) {

  //std::vector<triplet> * triplets = new std::vector<triplet>;

  // Cut on the triplet track angle: +- 10 mrad
  //double triplet_angle_cut = 0.010;
  // Cut on the triplet residual on the middle plane:
  //double triplet_residual_cut = 0.1; // [mm]

  for( std::vector<hit>::iterator ihit = hits.begin(); ihit != hits.end(); ihit++ ){
    if( (*ihit).plane != plane0 ) continue; // First plane

    for( std::vector<hit>::iterator jhit = hits.begin(); jhit != hits.end(); jhit++ ){
      if( (*jhit).plane != plane2 ) continue; // Last plane

      //float sum_res = 999.;
      //float sum_res_old = 999.;
      //bool IsFirst = true;

      for( std::vector<hit>::iterator khit = hits.begin(); khit != hits.end(); khit++ ){
	if( (*khit).plane != plane1 ) continue; // Middle plane

	// Create new preliminary triplet from the three hits:
	triplet new_triplet((*ihit),(*khit),(*jhit));

	//if(_nEvt < 10) if( plane0 == 0 && plane1 == 1 && plane2 == 2) streamlog_out( WARNING2 ) << "dx triplet = " << new_triplet.getdx(plane1)*1e3 << endl;
	//if(_nEvt < 10) if( plane0 == 3 && plane1 == 4 && plane2 == 5) streamlog_out( WARNING2 ) << "dx driplet = " << new_triplet.getdx(plane1)*1e3 << endl;

	// Setting cuts on the triplet track angle:
	if( abs(new_triplet.getdx()) > _slope_cut_x * new_triplet.getdz()) continue;
	if( abs(new_triplet.getdy()) > _slope_cut_y * new_triplet.getdz()) continue;

	// Setting cuts on the triplet residual on the middle plane
	if( abs(new_triplet.getdx(plane1)) > _triCut) continue;
	if( abs(new_triplet.getdy(plane1)) > _triCut) continue;

	/*
	// For low threshold (high noise) and/or high occupancy, use only the triplet with the smallest sum of residuals on plane1
	sum_res = abs(new_triplet.getdx(plane1)) + abs(new_triplet.getdy(plane1));
	if(sum_res < sum_res_old || IsFirst){

	// Remove the last one since it fits worse, not if its the first
	if(!IsFirst) triplets->pop_back();
	// The triplet is accepted, push it back:
	triplets->push_back(new_triplet);
	IsFirst = false;
	streamlog_out(DEBUG2) << new_triplet;
	sum_res_old = sum_res;
	}
	*/

	// The triplet is accepted, push it back:
	triplets.push_back(new_triplet);
	streamlog_out(DEBUG2) << new_triplet;

      }//loop over hits
    }//loop over hits
  }// loop over hits

  //return triplets;
}

bool EUTelTripletGBL::IsTripletIsolated(std::vector<triplet>::iterator it, std::vector<triplet> &trip, double z_match, double isolation) { // isolation is defaulted to 0.3 mm
  bool IsolatedTrip = true;

  // check first if trip is isolated
  double xA = (*it).getx_at(z_match); // triplet impact point at matching position
  double yA = (*it).gety_at(z_match);


  double ddAMin = -1.0;
  for( std::vector<triplet>::iterator tripIsoCheck = trip.begin(); tripIsoCheck != trip.end(); tripIsoCheck++ ) {
    if(it != tripIsoCheck){
      double xAIsoCheck = (*tripIsoCheck).getx_at(z_match);
      double yAIsoCheck = (*tripIsoCheck).gety_at(z_match);
      double ddA = sqrt( fabs(xAIsoCheck - xA)*fabs(xAIsoCheck - xA) 
	  + fabs(yAIsoCheck - yA)*fabs(yAIsoCheck - yA) );
      if(ddAMin < 0 || ddA < ddAMin) ddAMin = ddA;
    }
  }

  //triddaMindutHisto->fill(ddAMin);
  if(ddAMin < isolation) IsolatedTrip = false;

  return IsolatedTrip;
}

void EUTelTripletGBL::MatchTriplets(std::vector<triplet> &up, std::vector<triplet> &down, double z_match, std::vector<EUTelTripletGBL::track> &tracks) {

  //std::vector<track> * tracks = new std::vector<track>;
  // Cut on the matching of two triplets [mm]
  //double intersect_residual_cut = 0.1;

  for( std::vector<triplet>::iterator trip = up.begin(); trip != up.end(); trip++ ){

    // Track impact position at Matching Point from Upstream:
    double xA = (*trip).getx_at(z_match); // triplet impact point at matching position
    double yA = (*trip).gety_at(z_match);

    // check if trip is isolated
    bool IsolatedTrip = IsTripletIsolated(trip, up, z_match, _triCut*2.01);

    for( std::vector<triplet>::iterator drip = down.begin(); drip != down.end(); drip++ ){

      // Track impact position at Matching Point from Downstream:
      double xB = (*drip).getx_at(z_match); // triplet impact point at matching position
      double yB = (*drip).gety_at(z_match);

      // check if drip is isolated
      bool IsolatedDrip = IsTripletIsolated(drip, down, z_match, _triCut*2.01);


      // Build a track candidate from one upstream and one downstream triplet:
      track newtrack((*trip),(*drip));

      // Track kinks as difference in triplet slopes:
      double kx = newtrack.kink_x();
      double ky = newtrack.kink_y();

      // driplet - triplet
      double dx = xB - xA; 
      double dy = yB - yA;

      sixkxHisto->fill( kx*1E3 );
      sixkyHisto->fill( ky*1E3 );
      sixdxHisto->fill( dx );
      sixdyHisto->fill( dy );


      if( abs(dy) < 0.4 ) sixdxcHisto->fill( dx*1E3 ); // sig = 17 um at 5 GeV
      if( abs(dx) < 0.4 ) sixdycHisto->fill( dy*1E3 );

      // match driplet and triplet:
      if( abs(dx) > _cutx) continue;
      if( abs(dy) > _cuty) continue;

      // check isolation
      if( !IsolatedTrip || !IsolatedDrip ) {
	hIso->fill(0);
	//continue;
      }
      else hIso->fill(1);

      sixkxcHisto->fill( kx*1E3 );
      sixkycHisto->fill( ky*1E3 );
      sixxHisto->fill( -xA ); // -xA = x_DP = out
      sixyHisto->fill( -yA ); // -yA = y_DP = up
      sixxyHisto->fill( -xA, -yA ); // DP: x_out, y_up

      // Fill kink map histogram:
      if( abs( kx ) > 0.002 || abs( ky ) > 0.002 ) sixxycHisto->fill( -xA, -yA );

      kinkvsxy->fill( -xA, -yA, (kx*kx + ky*ky)*1E6 ); //<kink^2> [mrad^2]

      // apply fiducial cut
      if ( fabs(xA) >  9.0) continue;
      if (     -yA  < -4.0) continue;

      // Add the track to the vector if trip/drip are isolated, the triplets are matched, and all other cuts are passed
      tracks.push_back(newtrack);

    } // Downstream
  } // Upstream

  streamlog_out(DEBUG2) << "Found " << tracks.size() << " tracks from matched triplets." << std::endl;
  //return tracks;
}



EUTelTripletGBL::track::track(triplet up, triplet down) : upstream(up), downstream(down) {}

double EUTelTripletGBL::track::kink_x() {
  return (downstream.slope().x - upstream.slope().x);
}

double EUTelTripletGBL::track::kink_y() {
  return (downstream.slope().y - upstream.slope().y);
}

EUTelTripletGBL::hit EUTelTripletGBL::track::intersect() {
  hit inter;
  // Re-check what this actually is...
  // and simplifie using triplet class members...
  inter.x = ( upstream.base().x - upstream.slope().x * upstream.base().z - downstream.base().x + downstream.slope().x * downstream.base().z ) / kink_x();
  inter.y = ( upstream.base().y - upstream.slope().y * upstream.base().z - downstream.base().y + downstream.slope().y * downstream.base().z ) / kink_y();
  return inter;
}

EUTelTripletGBL::triplet EUTelTripletGBL::track::get_upstream() {
  return upstream;
}

EUTelTripletGBL::triplet EUTelTripletGBL::track::get_downstream() {
  return downstream;
}

EUTelTripletGBL::hit EUTelTripletGBL::track::gethit(int plane) {
  if(plane < 3) return upstream.gethit(plane);
  else return downstream.gethit(plane);
}



EUTelTripletGBL::triplet::triplet() : linked_dut(false), hits() {
  // Empty default constructor
}

EUTelTripletGBL::triplet::triplet(hit hit0, hit hit1, hit hit2) : linked_dut(false), hits() {
  triplet();
  filltriplet(hit0, hit1, hit2);
}

EUTelTripletGBL::hit EUTelTripletGBL::triplet::getpoint_at(double z) {
  hit impact;
  impact.z = z - base().z;
  impact.x = base().x + slope().x * impact.z;
  impact.y = base().y + slope().y * impact.z;
  return impact;
}

double EUTelTripletGBL::triplet::getx_at(double z) {
  return this->base().x + this->slope().x * (z - this->base().z);
}

double EUTelTripletGBL::triplet::getdx() {
  return this->hits.rbegin()->second.x - this->hits.begin()->second.x;
}

double EUTelTripletGBL::triplet::getdx(int ipl) {
  return this->hits[ipl].x - this->base().x - this->slope().x * (this->hits[ipl].z - this->base().z);
}

double EUTelTripletGBL::triplet::getdx(hit point) {
  return point.x - this->base().x - this->slope().x * (point.z - this->base().z);
}

double EUTelTripletGBL::triplet::gety_at(double z) {
  return this->base().y + this->slope().y * (z - this->base().z);
}

double EUTelTripletGBL::triplet::getdy() {
  return this->hits.rbegin()->second.y - this->hits.begin()->second.y;
}

double EUTelTripletGBL::triplet::getdy(int ipl) {
  return this->hits[ipl].y - this->base().y - this->slope().y * (this->hits[ipl].z - this->base().z);
}

double EUTelTripletGBL::triplet::getdy(hit point) {
  return point.y - this->base().y - this->slope().y * (point.z - this->base().z);
}

double EUTelTripletGBL::triplet::getdz() {
  return this->hits.rbegin()->second.z - this->hits.begin()->second.z;
}

EUTelTripletGBL::hit EUTelTripletGBL::triplet::gethit(int plane) {
  return this->hits[plane];
}

EUTelTripletGBL::hit EUTelTripletGBL::triplet::base() {
  hit center;
  center.x = 0.5*( this->hits.begin()->second.x + this->hits.rbegin()->second.x );
  center.y = 0.5*( this->hits.begin()->second.y + this->hits.rbegin()->second.y );
  center.z = 0.5*( this->hits.begin()->second.z + this->hits.rbegin()->second.z );
  return center;
}

EUTelTripletGBL::hit EUTelTripletGBL::triplet::slope() {
  hit sl;
  double dz = (this->hits.rbegin()->second.z - this->hits.begin()->second.z);
  sl.x = (this->hits.rbegin()->second.x - this->hits.begin()->second.x) / dz;
  sl.y = (this->hits.rbegin()->second.y - this->hits.begin()->second.y) / dz;
  return sl;
}

#endif
