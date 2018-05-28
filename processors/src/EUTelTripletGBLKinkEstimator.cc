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
#include "EUTelTripletGBLKinkEstimator.h"
#include "EUTelTripletGBLKinkEstimatorInstance.h"
#include "EUTelTripletGBLUtility.h"

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
#include <AIDA/IHistogram3D.h>
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

// marlin util includes
#include "Mille.h"

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
#include <sstream>
#include <memory>
#include <map>

// ROOT includes ".h"
#include <TMath.h>
#include <TRotation.h>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;

// MillePede-II binary file
gbl::MilleBinary * milleGBL2; 
//gbl::MilleBinary * milleGBLgood; 

EUTelTripletGBLKinkEstimator::EUTelTripletGBLKinkEstimator() : Processor("EUTelTripletGBLKinkEstimator"), _siPlanesParameters(), _siPlanesLayerLayout(), _inputCollectionTelescope(""), _isFirstEvent(0), _eBeam(0), _nEvt(0), _nTelPlanes(0), _dut_plane(-1), _track_match_cut(0.15),  _planeSort(), _planeID(), _planePosition(), _planeThickness(), _planeX0(), _planeResolution() {

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
      _eBeam, 0.0);

  registerOptionalParameter( "triResCut", "Upstream/Downstream triplet residual cut [mm]", _triplet_res_cut, 0.1 );

  registerProcessorParameter( "matchingCut",
      "cut for matching in x coordinate in mm",
      _track_match_cut, 0.15);

  registerProcessorParameter( "slopeCut",
      "cut for track slopes in x coordinate in rad",
      _slope_cut, 0.002);

  registerProcessorParameter( "dut_plane",
      "aluminium plane to be considered as passive DUT",
      _dut_plane, -1);

  registerProcessorParameter( "eff_radius",
      "radius on DUT plane to accept match with triplet",
      _eff_radius, 0.1); // in mm. 0.1 is for 6 GeV, 20 mm. This is scaled with E and dz

  registerProcessorParameter( "kappa",
      "global factor to Highland formula",
      _kappa, 1.0); // 1.0 means HL as is, 1.2 means 20% additional scattering

  registerProcessorParameter( "targetthick",
      "thickness of alu target, if present",
      _targetthick, 0.0);

  registerProcessorParameter( "probchi2Cut",
      "Cut on Prob(chi2,ndf) rejecting bad tracks with prob < cut",
      _probchi2_cut, .01); 

  registerOptionalParameter("Resolution",
      "resolution parameter for each Cluster size, same for all planes. first value is average of all CSes. Up to CS6 plus larger than 6, hence in total 8 numbers. Disable with -1, e.g. (3.5e-3, -1, -1, ...., -1)",
      _resolution,  FloatVec( 8, 3.5*1e-3));

  registerOptionalParameter("Thickness","thickness parameter for each plane. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_thickness,  FloatVec( 7, 50*1e-3));


}


void EUTelTripletGBLKinkEstimator::init() {

  // usually a good idea to
  printParameters();

  _nEvt = 0;
  _ngbl = 0;

  _isFirstEvent = true;

  streamlog_out(MESSAGE0) << "Beam energy " << _eBeam << " GeV" <<  std::endl;

  // Read geometry information from GEAR

  streamlog_out(MESSAGE0) << "Reading telescope geometry description from GEAR " << std::endl;

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
  _planePosActive  = new int[_nTelPlanes];
  _planePosActive[0] = 0;
  _planePosActive[1] = 1;
  _planePosActive[2] = 2;
  _planePosActive[3] = -1;
  _planePosActive[4] = 3;
  _planePosActive[5] = 4;
  _planePosActive[6] = 5;

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

  _triplet_res_cut = _triplet_res_cut *6. / _eBeam * (_planePosition[1] - _planePosition[0]) / 20.;

  // Book histograms:
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  bookHistos();
#endif
  
  // for mille binary
  std::string _binaryFilename = "milleKINK.bin";
  std::string _binaryFilenamegood = "milleKINKgood.bin";
  unsigned int reserveSize = 8000;
  milleGBL2 = new gbl::MilleBinary( _binaryFilename, reserveSize );
  //milleGBLgood = new gbl::MilleBinary( _binaryFilenamegood, reserveSize );
}//init

//------------------------------------------------------------------------------
void EUTelTripletGBLKinkEstimator::processRunHeader( LCRunHeader* runHeader) {

  auto eutelHeader = std::make_unique<EUTelRunHeaderImpl>( runHeader );
  eutelHeader->addProcessor( type() );

  // Decode and print out Run Header information - just a check

  _nRun = runHeader->getRunNumber();

  streamlog_out( MESSAGE2 )  << "Processing run header"
    << ", run nr " << runHeader->getRunNumber() << std::endl;

  const std::string detectorName = runHeader->getDetectorName();
  const std::string detectorDescription = runHeader->getDescription();

  streamlog_out( MESSAGE0 ) << detectorName << " : " << detectorDescription << std::endl;

  // we want the gblutil class (which is NOT a marlin processor) to be able to write its own histograms, but into the same root file of the processor class, which uses the util class. Therefore, we let gblutil know its parent, which knows about the AIDA histogram handle
  gblutil.setParent(this);
  // and we book some histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  gblutil.bookHistos();
#endif

} // processRunHeader

//----------------------------------------------------------------------------
void EUTelTripletGBLKinkEstimator::processEvent( LCEvent * event ) {

  if( _nEvt % 1000 == 0 ) {
    streamlog_out( MESSAGE2 ) << "Processing event "
      << setw(6) << setiosflags(ios::right)
      << event->getEventNumber() << " in run "
      << setw(6) << setiosflags(ios::right)
      << event->getRunNumber()
      << ", currently having "
      << _ngbl << " good GBL tracks "
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

  if(_nEvt < 10) streamlog_out( MESSAGE6 )  << "Total of " << collection->getNumberOfElements() << " tracker hits in input collection " << std::endl;
  nAllHitHisto->fill(collection->getNumberOfElements());

  //----------------------------------------------------------------------------

  std::vector<EUTelTripletGBLUtility::hit> hits;

  // Extract all hits from the LCIO collection and calculate the
  // position uncertainty if necessary:
  for( int ihit = 0; ihit < collection->getNumberOfElements(); ihit++ ) {

    EUTelTripletGBLUtility::hit newhit;
    TrackerHit * meshit = dynamic_cast<TrackerHit*>( collection->getElementAt(ihit) );
    const double * pos = meshit->getPosition();
    const EVENT::FloatVec cov = meshit->getCovMatrix();
    float charge = 0.;

    TrackerDataImpl* clusterVector = static_cast<TrackerDataImpl*>( meshit->getRawHits()[0]);
    EUTelSimpleVirtualCluster* cluster = nullptr;

    float locx = -1;
    float locy = -1;
    if ( meshit->getType() == kEUTelSparseClusterImpl ) 
    {
      cluster = new EUTelSparseClusterImpl< EUTelGenericSparsePixel > ( clusterVector );
      charge = cluster->getTotalCharge();
      cluster->getCenterOfGravity(locx, locy);
      streamlog_out(DEBUG2) << " charge is " << charge << std::endl;
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
    double distMin = 2; // [mm]

    bool foundplane = false;
    for( int ipl = 0; ipl < _nTelPlanes; ipl++ ) {
      streamlog_out(DEBUG3) << "hit z coord = " << newhit.z << "  plane Pos = " << _planePosition[ipl] << "  and abs dist = " << abs(newhit.z - _planePosition[ipl]) << std::endl;
      if( abs(newhit.z - _planePosition[ipl]) < distMin ) {
	newhit.plane = _planePosActive[ipl];
        newhit.id = _planeID[ipl];
	foundplane = true;
	distMin = abs(newhit.z - _planePosition[ipl]);
      }
    }
    streamlog_out(DEBUG3)  << " Hit is at plane " << newhit.plane << " with ID " << newhit.id << std::endl;

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
  //
  // a word on gblutil:
  // for future use, the (contructor of?) gblutil could know about the AIDAProcessor, so hitso dont need to be passed, but live in the Util class.
  // This opens the possiblity to create histograms within methods of the Util class, rather then uglily hacking them in the TripletGBL class.
  // Cut values could be passed and stored as private member and used in the methods, rather than passed during fct call.
  // tbs...



  // -----------------------------------------------------------------
  // Downstream Telescope Triplets ("driplets")

  // Generate new triplet set for the Telescope Downstream Arm:
  std::vector<EUTelTripletGBLUtility::triplet> downstream_triplets;
  gblutil.FindTriplets(hits, std::array<size_t,3>{3, 4, 5}, _triplet_res_cut, 5*_slope_cut, downstream_triplets);
  streamlog_out(DEBUG4) << "number of found driplets = " << downstream_triplets.size() << std::endl;


  // Iterate over all found downstream triplets to fill histograms and match them to the REF and DUT:
  for( std::vector<EUTelTripletGBLUtility::triplet>::iterator drip = downstream_triplets.begin(); drip != downstream_triplets.end(); drip++ ){

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
  std::vector<EUTelTripletGBLUtility::triplet> upstream_triplets;
  gblutil.FindTriplets(hits, std::array<size_t,3>{0, 1, 2}, _triplet_res_cut, _slope_cut, upstream_triplets);
  streamlog_out(DEBUG4) << "number of found triplets = " << upstream_triplets.size() << std::endl;

  // Iterate over all found upstream triplets to fill histograms and match them to the REF and DUT:
  for( std::vector<EUTelTripletGBLUtility::triplet>::iterator trip = upstream_triplets.begin(); trip != upstream_triplets.end(); trip++ ) {

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
    for( std::vector<EUTelTripletGBLUtility::hit>::iterator lhit = hits.begin(); lhit != hits.end(); lhit++ ){

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
  // six: triplets A and driplets B
  // matching and GBL fit
  // kinks: triplets A vs driplets B
  // scattering point = DUT:

  // now hardcoded to target
  streamlog_out(DEBUG4) << " alu plane (_dut_plane) is number " << _dut_plane << endl;
  double DUTz = _planePosition[_dut_plane];
  streamlog_out(DEBUG4) << " DUTz position is " << DUTz << endl;

  // Match the Telescope Upstream and Downstream Arm triplets to get tracks:
  std::vector<EUTelTripletGBLUtility::track> telescope_tracks;
  gblutil.MatchTriplets(upstream_triplets,downstream_triplets, DUTz, _track_match_cut, telescope_tracks);

  //delete downstream_triplets;
  //delete upstream_triplets;

  for( std::vector<EUTelTripletGBLUtility::track>::iterator tr = telescope_tracks.begin(); tr != telescope_tracks.end(); tr++ ){

    EUTelTripletGBLUtility::triplet trip = (*tr).get_upstream();
    EUTelTripletGBLUtility::triplet drip = (*tr).get_downstream();
    EUTelTripletGBLUtility::triplet srip((*tr).gethit(0), (*tr).gethit(2), (*tr).gethit(5)); // seed triplet is called 'srip'

    std::vector<double> xAplanes(_nTelPlanes);
    std::vector<double> yAplanes(_nTelPlanes);
    for (int i = 0; i < _nTelPlanes; i++){
      xAplanes[i] = trip.getx_at(_planePosition[i]);
      yAplanes[i] = trip.gety_at(_planePosition[i]);
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
    double xA = trip.getx_at(DUTz);
    double yA = trip.gety_at(DUTz);

    double dx = xB - xA; // driplet - triplet
    double dy = yB - yA;

    // GBL with triplet A as seed:

    std::vector<gbl::GblPoint> traj_points;

    // build up trajectory:
    std::vector<double> sPoint;

    // plane 0:
    double s = 0.0;

    Eigen::Matrix2d proL2m = Eigen::Matrix2d::Identity();
    Eigen::Matrix<double,2,4> addDer = Eigen::Matrix<double,2,4>::Zero();

    //double res = 3.42E-3; // [mm] Anemone telescope intrinsic resolution
    //res = 4.5E-3; // EUDET

    // scatter:
    Eigen::Vector2d scat = Eigen::Vector2d::Zero(); //mean is zero

    double p = _eBeam; // beam momentum
    double epsSi = -1;
    double epsAir = -1.; // define later when dz is known
    double epsAlu = (_targetthick)/88.97; // Alu target // FIXME this is used only in the log-correction, but needs assumption of X0 :/
    // another way would be so solve this iteratively: start wihout target correcttion in the log-correction, in second iteration use result from first, ...

    double sumeps = 0.0;
    double tetSi = -1;
    double tetAir = -1.;

    std::vector<unsigned int> ilab; // 0-5 = telescope

    // plane 0-5:

    double rx[6];
    double ry[6];
    double trackhitx[6];
    double trackhity[6];
    //double trackhitxloc[6];
    //double trackhityloc[6];
    //double zprev = _planePosition[0];
    double step = 0.;

    // loop over all scatterers first to calculate sumeps
    for( int ipl = 0; ipl < 7; ++ipl ){ // FIXME include estimate for unknown scatterer into log corr ?! for now, set thick = 0 in gear
      if(_thickness[ipl] > 0.0){
        epsSi = _thickness[ipl]/ 93.66 + 0.050 / 285.6; // Si + Kapton (Kapton ist 2 * 25  = 50 um thick, but has 1/3 rad length = 17 um)
        sumeps += epsSi;
      }
      if( ipl < 6) {
	double distplane = _planePosition[ipl+1] - _planePosition[ipl];
	epsAir =   distplane  / 304200.; 
	sumeps += epsAir;
      }

    }
    sumeps += epsAlu;
    streamlog_out(DEBUG4) << "sumeps after alu = " << sumeps << endl;
    // done with calculating sum eps

    
    int DUT_label = -1;
    int siplane_label = -1;
    double sDUT = -1;
    int ipl =0;
    for( int iplprime = 0; iplprime < _nTelPlanes; ++iplprime ){


      // passive DUTs are listed with IDs starting 100, 101, ..., at least in my convention
      if(_planeID[iplprime] <100 ){ 
	
	// Get the corresponding hit from the track:
	EUTelTripletGBLUtility::hit trackhit = (*tr).gethit(ipl);

	if (ipl == 0) clustersize0->fill(trackhit.clustersize);
	if (ipl == 1) clustersize1->fill(trackhit.clustersize);
	if (ipl == 2) clustersize2->fill(trackhit.clustersize);
	if (ipl == 3) clustersize3->fill(trackhit.clustersize);
	if (ipl == 4) clustersize4->fill(trackhit.clustersize);
	if (ipl == 5) clustersize5->fill(trackhit.clustersize);

	double dz = trackhit.z - srip.base().z; // here we use the srip rather then the trip (but srip can be = trip)
	double xs = srip.base().x + srip.slope().x * dz; // Ax at plane
	double ys = srip.base().y + srip.slope().y * dz; // Ay at plane

	trackhitx[ipl] = trackhit.x;
	trackhity[ipl] = trackhit.y;
	//trackhitxloc[ipl] = trackhit.locx;
	//trackhityloc[ipl] = trackhit.locy;

	rx[ipl] = trackhit.x - xs; // distance to SEED triplet srip (which can either be equal to trip, or constructed from e.g. planes (0,2,5)
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

	gbl::GblPoint * point = new gbl::GblPoint( gblutil.JacobianPointToPoint( step ) );

    Eigen::Vector2d meas(2);
	meas[0] = rx[ipl];
	meas[1] = ry[ipl];
	//meas[0] = trackhit.x;
	//meas[1] = trackhit.y;

	epsSi = _thickness[iplprime]/93.66 + 0.050 / 286.6; // Si + Kapton (Kapton ist 2 * 25  = 50 um thick, but has 1/3 rad length = 17 um)
	tetSi = _kappa * 0.0136 * sqrt(epsSi) / p * ( 1 + 0.038*std::log(sumeps) );

    Eigen::Vector2d wscatSi(2);
	wscatSi[0] = 1.0 / ( tetSi * tetSi ); //weight
	wscatSi[1] = 1.0 / ( tetSi * tetSi );

	double _resolution_tmp = -1.0;
	if(trackhit.clustersize == 1) _resolution_tmp = _resolution[1];
	if(trackhit.clustersize == 2) _resolution_tmp = _resolution[2];
	if(trackhit.clustersize == 3) _resolution_tmp = _resolution[3];
	if(trackhit.clustersize == 4) _resolution_tmp = _resolution[4];
	if(trackhit.clustersize == 5) _resolution_tmp = _resolution[5];
	if(trackhit.clustersize == 6) _resolution_tmp = _resolution[6];
	if(trackhit.clustersize >  6) _resolution_tmp = _resolution[7];
	//_resolution_tmp = 3.24*1e-3;
 
	Eigen::Vector2d measPrec; // precision = 1/resolution^2
	measPrec[0] = 1.0 / _resolution_tmp / _resolution_tmp;
	measPrec[1] = 1.0 / _resolution_tmp / _resolution_tmp;

	point->addMeasurement( proL2m, meas, measPrec );

	point->addScatterer( scat, wscatSi );

	s += step;
	streamlog_out(DEBUG3) << " s = " << s  << std::endl;
	sPoint.push_back( s );
	siplane_label = sPoint.size();
	ilab.push_back(siplane_label);
        // add local derrivative if after DUT
	// Matrix
	//     dkink_11     dkink_21    dkink12    dkink_22
	//  x    s-s1           0         s-s2       0
	//  y     0            s-s1        0        s-s2
	if(sDUT > 0){
	  addDer(0,0) = (s - (sDUT + _targetthick/sqrt(12))); // First scatterer in target
	  addDer(1,1) = (s - (sDUT + _targetthick/sqrt(12))); //
	  addDer(0,2) = (s - (sDUT - _targetthick/sqrt(12))); // second scatterer in target
	  addDer(1,3) = (s - (sDUT - _targetthick/sqrt(12))); //
	  streamlog_out(DEBUG3) << " lever arm left DUT-point = " << (s - (sDUT + _targetthick/sqrt(12))) << " and right DUT-point = " << (s - (sDUT - _targetthick/sqrt(12))) << std::endl;

	  point->addLocals(addDer);
	}

	traj_points.push_back(*point);
	delete point;

        ipl++;
      } else {

        // now the plane should be the passive one
	gbl::GblPoint * point = new gbl::GblPoint( gblutil.JacobianPointToPoint( step ) );

        // we dont add a scatterer here, as we want it to be unbiased (addition of local derivative)
	/*double tetAlu = _kappa*0.0136 * sqrt(epsAlu) / p * ( 1 + 0.038*std::log(sumeps) );
	if(_targetthick < 1.) tetAlu = 1.e-10;
        std::cout << " tetAlu = " << tetAlu << std::endl;

	TVectorD wscatAlu(2);
	wscatAlu[0] = 1.0 / ( tetAlu * tetAlu ); // weight
	wscatAlu[1] = 1.0 / ( tetAlu * tetAlu ); 

	point->addScatterer( scat, wscatAlu );
	*/

	traj_points.push_back(*point);
	s += step;
	sPoint.push_back( s );
	DUT_label = sPoint.size();
	ilab.push_back(DUT_label);

	sDUT = s;
	streamlog_out(DEBUG3)  << " s_DUT = " << sDUT  << std::endl;
	delete point;
      }






      if( iplprime < 6) {
	double distplane = _planePosition[iplprime+1] - _planePosition[iplprime];

	streamlog_out(DEBUG3) << " dist plane = " << distplane  << std::endl;
	step = 0.21*distplane;
	epsAir =   0.5*distplane  / 304200.; 
	tetAir = _kappa * 0.0136 * sqrt(epsAir) / p * ( 1 + 0.038*std::log(sumeps) ); 

    Eigen::Vector2d wscatAir(2);
	wscatAir[0] = 1.0 / ( tetAir * tetAir ); // weight
	wscatAir[1] = 1.0 / ( tetAir * tetAir ); 

	gbl::GblPoint * point = new gbl::GblPoint( gblutil.JacobianPointToPoint( step ) );
	point->addScatterer( scat, wscatAir );

	s += step;
	traj_points.push_back(*point);
	sPoint.push_back( s );
	delete point;

	step = 0.58*distplane;
	gbl::GblPoint * point1 = new gbl::GblPoint( gblutil.JacobianPointToPoint( step ) );
	point1->addScatterer( scat, wscatAir );

	s += step;
	traj_points.push_back(*point1);
	sPoint.push_back( s );
	delete point1;

	step = 0.21*distplane; // remaing distance to next plane
      }


    } // loop over planes
    // monitor what we put into GBL:
    selxHisto->fill( -xA ); // triplet at DUT
    selyHisto->fill( -yA );
    selaxHisto->fill( trip.slope().x*1E3 );
    selayHisto->fill( trip.slope().y*1E3 );
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

    
    if(_nEvt < 10){
      streamlog_out(MESSAGE4) << "traj with " << traj.getNumPoints() << " points:" << endl;
      for( int iplb = 0; iplb < 7; iplb++ ){
        if(_planeID[iplb] < 100){
	  streamlog_out(MESSAGE4) << "  plane " << _planePosActive[iplb] << ", label " << ilab[iplb];
	  streamlog_out(MESSAGE4) << "  z " << sPoint[ilab[iplb]-1];
	  streamlog_out(MESSAGE4) << " dx " << rx[_planePosActive[iplb]];
	  streamlog_out(MESSAGE4) << " dy " << ry[_planePosActive[iplb]];
	  streamlog_out(MESSAGE4) << " chi2 " << Chi2;
	  streamlog_out(MESSAGE4) << " ndf " << Ndf;
	  streamlog_out(MESSAGE4) << endl;
	}
      }

      streamlog_out(DEBUG4)  << " Is traj valid? " << traj.isValid() << std::endl;
      traj.printPoints();
      //traj.printTrajectory();
      //traj.printData();
    }


    gblndfHisto->fill( Ndf );
    gblchi2Histo->fill( Chi2 );

    double probchi = 0;
    probchi = TMath::Prob( Chi2, Ndf );
    gblprbHisto->fill( probchi );

    gblprbxHisto->fill(-xA, probchi);
    gblprbyHisto->fill(-yA, probchi);

    // bad fits:

    if( probchi < 0.01 ) {

      badxHisto->fill( -xA ); // triplet at DUT
      badyHisto->fill( -yA );
      badaxHisto->fill( trip.slope().x*1E3 );
      badayHisto->fill( trip.slope().y*1E3 );
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

    traj.milleOut( *milleGBL2 );
    if( probchi > _probchi2_cut){

      _ngbl++;

      Eigen::VectorXd aCorrection(9);
      Eigen::MatrixXd aCovariance(9,9);
      
      double ax[8];
      // double ay[8];
      unsigned int k = 0;


      Eigen::VectorXd aResiduals;
      Eigen::VectorXd aMeasErrors;
      Eigen::VectorXd aResErrors;
      Eigen::VectorXd aDownWeights;


      Eigen::VectorXd aKinks;
      Eigen::VectorXd aKinkErrors;
      Eigen::VectorXd kResErrors;
      Eigen::VectorXd kDownWeights;

      unsigned int ndata = 2;
      //track = q/p, x', y', x, y
      //        0,   1,  2,  3, 4
      // at plane DUT:
      //
      //std::cout << " size of ilab = " << ilab.size() << std::endl;

      unsigned int ipos = ilab[0];
      traj.getResults( ipos, aCorrection, aCovariance );
      traj.getMeasResults( ipos, ndata, aResiduals, aMeasErrors, aResErrors, aDownWeights );
      traj.getScatResults( ipos, ndata, aKinks, aKinkErrors, kResErrors, kDownWeights );
      gblax0Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
      gbldx0Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
      gbldx01Histo->fill( aCorrection[3] ); // shift x [mm]
      gblrx0Histo->fill( ( aResiduals[0] ) * 1E3 ); // residual x [um]
      gblry0Histo->fill( ( aResiduals[1] ) * 1E3 ); // residual y [um]
      gblpx0Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
      gblpy0Histo->fill( aResiduals[1] / aResErrors[1] ); // pull
      gblqx0Histo->fill( aKinks[0]*1E3 ); // kink RESIDUAL (measured kink - fit kink)
      ax[k] = aCorrection[1]; // angle correction at plane, for kinks
      // TProfile for res_x a.f.o. x
      gblrxvsx0->fill( xAplanes.at(k), sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); // triplet extrapolation
      gblryvsy0->fill( yAplanes.at(k), sqrt(TMath::Pi()/2.)*fabs(aResiduals[1]));
      gblrxvsx01->fill((trackhitx[0] - aResiduals[0]), sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); // track fit
      gblryvsy01->fill((trackhity[0] - aResiduals[1]), sqrt(TMath::Pi()/2.)*fabs(aResiduals[1]));

      //std::cout << "Plane 0: \n aCorr: [0] = " << aCorrection[0]<< "  aCoor: [1] = " << aCorrection[1]<< "  aCoor: [2] = " << aCorrection[2]<< "  aCoor: [3] = " << aCorrection[3]<< "  aCoor: [4] = " << aCorrection[4]<< "  aCoor: [5] = " << aCorrection[5]<< "  aCoor: [6] = " << aCorrection[6] << std::endl;
      k++;

      ipos = ilab[1];
      traj.getResults( ipos, aCorrection, aCovariance );
      traj.getMeasResults(ipos, ndata, aResiduals, aMeasErrors, aResErrors, aDownWeights );
      traj.getScatResults(ipos, ndata, aKinks, aKinkErrors, kResErrors, kDownWeights );
      gblax1Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
      gbldx1Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
      gbldx11Histo->fill( aCorrection[3] ); // shift x [mm]
      gblrx1Histo->fill( ( aResiduals[0] ) * 1E3 ); // residual x [um]
      gblry1Histo->fill( ( aResiduals[1] ) * 1E3 ); // residual y [um]
      gblpx1Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
      gblpy1Histo->fill( aResiduals[1] / aResErrors[1] ); // pull
      gblqx1Histo->fill( aKinks[0]*1E3 ); // kink-residual (NOT KINK itself!)
      gblsx1Histo->fill( aKinks[0]/aKinkErrors[0] ); // x kink pull
      gbltx1Histo->fill( aKinks[0]/kResErrors[0] ); // x kink pull
      ax[k] = aCorrection[1]; // angle correction at plane, for kinks
      //ay[k] = aCorrection[2]; // angle correction at plane, for kinks
      //std::cout << " aResiduals[0] = " << aResiduals[0] << " aResErrors[0] = " << aResErrors[0] << std::endl;  
      k++;


      ipos = ilab[2];
      traj.getResults( ipos, aCorrection, aCovariance );
      traj.getMeasResults(ipos, ndata, aResiduals, aMeasErrors, aResErrors, aDownWeights );
      traj.getScatResults(ipos, ndata, aKinks, aKinkErrors, kResErrors, kDownWeights );
      gblax2Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
      gbldx2Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
      gbldx21Histo->fill( aCorrection[3] ); // shift x [mm]
      gblrx2Histo->fill( ( aResiduals[0] ) * 1E3 ); // residual x [um]
      gblry2Histo->fill( ( aResiduals[1] ) * 1E3 ); // residual y [um]
      gblpx2Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
      gblpy2Histo->fill( aResiduals[1] / aResErrors[1] ); // pull
      gblqx2Histo->fill( aKinks[0]*1E3 ); // kink-residual
      gblsx2Histo->fill( aKinks[0]/aKinkErrors[0] ); // x kink res over kinkError
      gbltx2Histo->fill( aKinks[0]/kResErrors[0] ); // x kink pull
      ax[k] = aCorrection[1]; // angle correction at plane, for kinks
      //ay[k] = aCorrection[2]; // angle correction at plane, for kinks
      k++;

      ipos = ilab[3];
      traj.getResults( ipos, aCorrection, aCovariance );
      gblax6Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
      gblay6Histo->fill( aCorrection[2]*1E3 ); // angle y [mrad]
      gblax6primeHisto->fill( aCorrection[5]*1E3 ); // angle x [mrad]
      gblay6primeHisto->fill( aCorrection[6]*1E3 ); // angle y [mrad]
      gblax6prime2Histo->fill( aCorrection[7]*1E3 ); // angle x [mrad]
      gblay6prime2Histo->fill( aCorrection[8]*1E3 ); // angle y [mrad]
      gblax6primeDiffHisto->fill( (aCorrection[7] + aCorrection[5])*1E3 ); // angle x [mrad]
      gblaxprimeayprime->fill(aCorrection[5]*1e3, aCorrection[6]*1e3);
      gblaxy6Histo->fill( (fabs(aCorrection[5]) + fabs(aCorrection[6]))/2.*1E3 ); // angle x [mrad]
      // get the cov 5 (and 6) here for variance of additional local derivative
      gblDUTkinkuncertHisto->fill(  sqrt(aCovariance(5,5))*1e3 ); // angle x [mrad]
      gblDUTkinkprimeuncertHisto->fill(  sqrt(aCovariance(7,7))*1e3 ); // angle x [mrad]
      gblDUTdecorrkinkuncertHisto->fill(  sqrt( aCovariance(5,5) + aCovariance(7,7) +2*aCovariance(5,7) )*1e3 ); // angle x [mrad] // CovMatrix is symmetric

      //std::cout << " sqrt aCov5 *1e3 = " << sqrt(aCovariance(5,5))*1e3 << std::endl;
      //std::cout << " sigma kink_DUT = " << sqrt(aCovariance(5,5) + aCovariance(6,6)) << std::endl;
      gblaxvsxy->fill( -xA, -yA, fabs(aCorrection[5])*1E3 ); //sqrt(<kink^2>) [mrad]
      gblayvsxy->fill( -xA, -yA, fabs(aCorrection[6])*1E3 ); //sqrt(<kink^2>) [mrad]
      gblaxyvsxy->fill( -xA, -yA, ( fabs(aCorrection[5]) + fabs(aCorrection[6]) )/2.*1E3 ); // [mrad]
      gblaxprime3D->fill(-xA, -yA, (aCorrection[7] + aCorrection[5])*1e3,1.);
      gblayprime3D->fill(-xA, -yA, (aCorrection[8] + aCorrection[6])*1e3,1.);
      gblaxprime6vsx->fill( -xA, sqrt(TMath::Pi()/2.)*fabs((aCorrection[7] + aCorrection[5])*1e3)); // trip extrapolation, kink angle
      ax[k] = aCorrection[1]; // angle correction at plane, for kinks

      //std::cout << "Plane DUT: \n aCorr: [0] = " << aCorrection[0]<< "  aCoor: [1] = " << aCorrection[1]<< "  aCoor: [2] = " << aCorrection[2]<< "  aCoor: [3] = " << aCorrection[3]<< "  aCoor: [4] = " << aCorrection[4]<< "  aCoor: [5] = " << aCorrection[5]<< "  aCoor: [6] = " << aCorrection[6] << std::endl;
      k++;

      ipos = ilab[4];
      traj.getResults( ipos, aCorrection, aCovariance );
      traj.getMeasResults(ipos, ndata, aResiduals, aMeasErrors, aResErrors, aDownWeights );
      traj.getScatResults(ipos, ndata, aKinks, aKinkErrors, kResErrors, kDownWeights );
      gblax3Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
      gbldx3Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
      gbldx31Histo->fill( aCorrection[3] ); // shift x [mm]
      gblrx3Histo->fill( ( aResiduals[0] ) * 1E3 ); // residual x [um]
      gblry3Histo->fill( ( aResiduals[1] ) * 1E3 ); // residual y [um]
      gblpx3Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
      gblpy3Histo->fill( aResiduals[1] / aResErrors[1] ); // pull
      gblqx3Histo->fill( aKinks[0]*1E3 ); // kink
      gblsx3Histo->fill( aKinks[0]/aKinkErrors[0] ); // x kink pull
      gbltx3Histo->fill( aKinks[0]/kResErrors[0] ); // x kink pull
      ax[k] = aCorrection[1]; // angle correction at plane, for kinks
      //
      //std::cout << "Plane 3: \n aCorr: [0] = " << aCorrection[0]<< "  aCoor: [1] = " << aCorrection[1]<< "  aCoor: [2] = " << aCorrection[2]<< "  aCoor: [3] = " << aCorrection[3]<< "  aCoor: [4] = " << aCorrection[4]<< "  aCoor: [5] = " << aCorrection[5]<< "  aCoor: [6] = " << aCorrection[6] << std::endl;
      k++;

      ipos = ilab[5];
      traj.getResults( ipos, aCorrection, aCovariance );
      traj.getMeasResults(ipos, ndata, aResiduals, aMeasErrors, aResErrors, aDownWeights );
      traj.getScatResults(ipos, ndata, aKinks, aKinkErrors, kResErrors, kDownWeights );
      gblax4Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
      gbldx4Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
      gbldx41Histo->fill( aCorrection[3] ); // shift x [mm]
      gblrx4Histo->fill( ( aResiduals[0] ) * 1E3 ); // residual x [um]
      gblry4Histo->fill( ( aResiduals[1] ) * 1E3 ); // residual y [um]
      gblpx4Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
      gblpy4Histo->fill( aResiduals[1] / aResErrors[1] ); // pull
      gblqx4Histo->fill( aKinks[0]*1E3 ); // kink
      gblsx4Histo->fill( aKinks[0]/aKinkErrors[0] ); // x kink pull
      gbltx4Histo->fill( aKinks[0]/kResErrors[0] ); // x kink pull
      ax[k] = aCorrection[1]; // angle correction at plane, for kinks
      k++;

      ipos = ilab[6];
      traj.getResults( ipos, aCorrection, aCovariance );
      traj.getMeasResults(ipos, ndata, aResiduals, aMeasErrors, aResErrors, aDownWeights );
      traj.getScatResults(ipos, ndata, aKinks, aKinkErrors, kResErrors, kDownWeights );
      gblax5Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
      gbldx5Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
      gbldx51Histo->fill( aCorrection[3] ); // shift x [mm]
      gblrx5Histo->fill( ( aResiduals[0] ) * 1E3 ); // residual x [um]
      gblry5Histo->fill( ( aResiduals[1] ) * 1E3 ); // residual y [um]
      gblpx5Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
      gblpy5Histo->fill( aResiduals[1] / aResErrors[1] ); // pull
      gblqx5Histo->fill( aKinks[0]*1E3 ); // kink
      ax[k] = aCorrection[1]; // angle correction at plane, for kinks
      //
      k++;

      // do some more analysis with kinks. 
      //   Calculate the two angles as corr_n-1 - corr_n-2, corr_n - corr_n-1, corr_n+1 - corr_n -> sum of all =  corr_n+1 - corr_n-2

      //ipos = DUT_label-2;
      //traj.getResults( ipos, aCorrection, aCovariance );
      //double kink_upstream = aCorrection[1];

      //ipos = DUT_label+1; // 1 downstream of centre
      //traj.getResults( ipos, aCorrection, aCovariance );
      //double kink_downstream = aCorrection[1];

      //gblkxCentreHisto->fill( (kink_downstream - kink_upstream) ); // kink at air/alu (sum of neighbours) [rad]
      //gblkxCentre1Histo->fill((kink_downstream - kink_upstream)*1e3 ); // kink at air/alu (sum of neighbours) [mrad]
      //gblkxCentreHisto->fill( (drip.slope().x - trip.slope().x) ); // kink at air/alu (sum of neighbours) [rad]
      //gblkxCentre1Histo->fill((drip.slope().x - trip.slope().x)*1e3 ); // kink at air/alu (sum of neighbours) [mrad]

      gblkx1Histo->fill( (ax[1] - ax[0])*1E3 ); // kink at 1 [mrad]
      gblkx2Histo->fill( (ax[2] - ax[1])*1E3 ); // kink at 2 [mrad]
      gblkx3Histo->fill( (ax[3] - ax[2])*1E3 ); // kink at 3 [mrad]
      gblkx4Histo->fill( (ax[4] - ax[3])*1E3 ); // kink at 4 [mrad]
      gblkx5Histo->fill( (ax[5] - ax[4])*1E3 ); // kink at 5 [mrad] 
      gblkx6Histo->fill( (ax[6] - ax[5])*1E3 ); // kink at 6 [mrad]
      
      
      //traj.milleOut( *milleGBLgood );
    } // end if good fit 

    

    //------------------------------------------------------------------------
    // intersect point in z:
    EUTelTripletGBLUtility::hit intersect = (*tr).intersect();
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
void EUTelTripletGBLKinkEstimator::check( LCEvent * /* evt */  ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


//------------------------------------------------------------------------------
void EUTelTripletGBLKinkEstimator::end(){

  // Print the summary:
  streamlog_out(MESSAGE5)
    << "---------------------------------------------------------------------------------------------------------" << std::endl
    << std::endl
    << "Processed events:    "
    << std::setw(10) << std::setiosflags(std::ios::right)
    << _nEvt << std::resetiosflags(std::ios::right) << std::endl;

  // close the output file:
  delete milleGBL2;
} // end end


void EUTelTripletGBLKinkEstimator::TelescopeCorrelationPlots(std::vector<EUTelTripletGBLUtility::hit> &telescopehits) {
  for( std::vector<EUTelTripletGBLUtility::hit>::iterator ihit = telescopehits.begin(); ihit != telescopehits.end(); ihit++ ){

    int ipl = (*ihit).plane;

    for( std::vector<EUTelTripletGBLUtility::hit>::iterator jhit = telescopehits.begin(); jhit != telescopehits.end(); jhit++ ){

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

#endif
