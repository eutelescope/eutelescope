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
#include "EUTelGBLFitter.h"
#include "EUTelTripletGBL.h"
#include "EUTelTripletGBLUtility.h"

#include "EUTELESCOPE.h"
#include "EUTelTrackerDataInterfacerImpl.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"

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

// LCIO includes <.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>

// system includes <>
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


EUTelGBLFitter::EUTelGBLFitter() : Processor("EUTelGBLFitter"), _inputCollectionTelescope(""), _isFirstEvent(0), _eBeam(0), _nEvt(0), _nPlanes(0), _dut_plane(3), _track_match_cut(0.15),  _planePosition() {
  // modify processor description
  _description = "Analysis for DATURA reference analysis ";

  // processor parameters
  registerInputCollection( LCIO::TRACKERHIT, "InputCollection", "Name of the input TrackerHit collection of the telescope",
      _inputCollectionTelescope, std::string("") ); // no collection defaulted forcing the user to pass something meaningful

  registerProcessorParameter( "Ebeam", "Beam energy [GeV]", _eBeam, static_cast<double>(0.0));

  registerOptionalParameter( "triResCut", "Upstream/Downstream triplet residual cut [mm]", _triplet_res_cut, 0.1);

  registerProcessorParameter( "matchingCut", "cut for matching in x coordinate in mm", _track_match_cut, static_cast<double>(0.15));

  registerProcessorParameter( "slopeCut", "cut for track slopes in x coordinate in rad", _slope_cut, static_cast<double>(0.002));

  registerProcessorParameter( "dut_plane", "plane to be considered the DUT and excluded from the track fit",
      _dut_plane, static_cast <int>(3));

  //in mm. 0.1 is for 6 GeV, 20 mm. This is scaled with E and dz
  registerProcessorParameter( "eff_radius", "radius on DUT plane to accept match with triplet", _eff_radius, static_cast <double>(0.1)); 

  //1.0 means HL as is, 1.2 means 20% additional scattering
  registerProcessorParameter( "kappa", "global factor to Highland formula", _kappa, static_cast<double>(1.0));

  registerProcessorParameter( "probchi2Cut", "Cut on Prob(chi2,ndf) rejecting bad tracks with prob < cut", _probchi2_cut, static_cast<double>(.01)); 

  registerOptionalParameter( "Resolution", "resolution parameter for each Cluster size, same for all planes. first value is average of all CSes. Up to CS6 plus larger than 6, hence in total 8 numbers. Disable with -1, e.g. (3.5e-3, -1, -1, ...., -1)", _resolution,  
  FloatVec(static_cast<double>(8), 3.5*1e-3));

  registerOptionalParameter( "Thickness","thickness parameter for each plane. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.", _thickness, FloatVec(static_cast<double>(6), 50*1e-3));

  registerProcessorParameter( "ClusterSizeSwitch", "boolian to switch b/w cluster charge (0, false) and cluster dimension (1, true) used as cluster size",
  _CSswitch, static_cast<bool>(true));
}

void EUTelGBLFitter::init() {
  // usually a good idea to
  printParameters();
  _isFirstEvent = true;
  _nEvt = 0;
  _printEventCounter= 0;
  _ngbl = 0;

  //This is the vector of sensorIDs ordered alogn the gloabl z-axis, 
  //this is guranteed by the framework 
  _sensorIDVec = geo::gGeometry().sensorIDsVec();
  _nPlanes = _sensorIDVec.size();

  for(auto& sensorID: _sensorIDVec) {
    _planePosition.emplace_back( geo::gGeometry().siPlaneZPosition(sensorID) );
    auto z = geo::gGeometry().siPlaneZSize(sensorID);
    auto rad = geo::gGeometry().siPlaneRadLength(sensorID);
    if(sensorID < 6) {
        _planeRadLength.emplace_back(55e-3 / 93.66 + 0.050 / 286.6); // Si + Kapton
    } else {
        _planeRadLength.emplace_back(z/rad);
    }
  }

  //to compute the total radiation length we will loop over all planes and add radiation
  //length for air from the first to the last plane
  double totalRadLength = (_planePosition.back()-_planePosition.front())/304200; //Radiation length for air is 30420 cm
  for(auto& radLen: _planeRadLength) {
    totalRadLength += radLen;
  }

  for(auto& radLen: _planeRadLength) {
      // Paper showed HL predicts too high angle, at least for biased measurement. 
      double tetSi = _kappa*0.0136 * sqrt(radLen) / _eBeam * ( 1 + 0.038*std::log(totalRadLength) );
      _planeWscatSi.emplace_back( 1.0/(tetSi*tetSi), 1.0/(tetSi*tetSi) );
  }

  for(size_t ipl = 0; ipl < _nPlanes-1; ipl++) {
    auto distplane = _planePosition[ipl+1] - _planePosition[ipl];
    double epsAir =   0.5*distplane  / 304200.;  // in [mm]
    double tetAir = _kappa*0.0136 * sqrt(epsAir) / _eBeam * ( 1 + 0.038*std::log(totalRadLength) );
    _planeWscatAir.emplace_back( 1.0/(tetAir*tetAir), 1.0/(tetAir*tetAir) );
  }

  streamlog_out(MESSAGE0) << "Beam energy " << _eBeam << " GeV" <<  std::endl;

  for(auto& sensorID: _sensorIDVec) {
    //streamlog_out( MESSAGE6 ) << "  Thickness Plane " << i << " = " << _thickness[i] << std::endl;
    //_planeResolution[i] = _siPlanesLayerLayout->getSensitiveResolution(i); // Get it from input
    //_planeResolution[i] = _resolution[0]; // pass the average here, which is the first entry of the vector
    streamlog_out( MESSAGE6 ) << "  Avg. reso Plane " << sensorID << " = " << _resolution[0] << std::endl;
    _planeResolutionX[sensorID] = _resolution;
    _planeResolutionY[sensorID] = _resolution;
  }

  streamlog_out( MESSAGE2 ) <<  "Telescope configuration with " << _nPlanes << " planes" << std::endl;

  for( size_t ipl = 0; ipl < _nPlanes; ipl++) {
    std::stringstream ss;
    ss << "  ID = " << "<ID>"
      << "  at Z [mm] = " << _planePosition[ipl]
      << " dZ [um] = " << "<planeThickness[ipl]*1000.>";
    //ss << "  average Res [um] = " << _planeResolution[ipl]*1000.;
    streamlog_out( MESSAGE2 ) <<  ss.str() << std::endl;
  }

/*
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
*/
  _triplet_res_cut = _triplet_res_cut *6. / _eBeam * (_planePosition[1] - _planePosition[0]) / 20.;

  // Book histograms:
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  bookHistos();
#endif
}//init

//------------------------------------------------------------------------------
void EUTelGBLFitter::processRunHeader( LCRunHeader* runHeader) {
  // we want the gblutil class (which is NOT a marlin processor) to be able to write its own histograms, but into the same root file of the processor class, which uses the util class. Therefore, we let gblutil know its parent, which knows about the AIDA histogram handle
  gblutil.setParent(this);
  gblutil.bookHistos();
  auto header = std::make_unique<EUTelRunHeaderImpl>(runHeader);
  header->addProcessor( type() ) ;
  // increment the run counter
  //++_iRun;
  // Decode and print out Run Header information - just a check
  _nRun = runHeader->getRunNumber();
  streamlog_out( MESSAGE2 )  << "Processing run header, run nr: " << runHeader->getRunNumber() << std::endl;
} // processRunHeader

//----------------------------------------------------------------------------
void EUTelGBLFitter::processEvent( LCEvent * event ) {

  if( _nEvt % 1000 == 0 ) {
    streamlog_out( MESSAGE2 ) << "Processing event "
      << setw(6) << setiosflags(ios::right)
      << event->getEventNumber() << " in run "
      << setw(6) << setiosflags(ios::right)
      << event->getRunNumber()
      << ", currently having "
      << _ngbl << " good tracks "
      << endl;
  }

  auto euEvent = static_cast<EUTelEventImpl*>( event );

  if( euEvent->getEventType() == kEORE ) {
    streamlog_out( DEBUG5 ) <<  "EORE found: nothing else to do." << std::endl;
    return;
  }

  int runNumber = event->getRunNumber();

  // Increase event count
  _nEvt++;

  //----------------------------------------------------------------------------
  // check input collection (aligned hits):

  LCCollection* collection = nullptr;
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

  CellIDDecoder<TrackerHit> hitCellDecoder(EUTELESCOPE::HITENCODING);
  std::vector<EUTelTripletGBLUtility::hit> hits;
  std::vector<EUTelTripletGBLUtility::hit> DUThits;

  for( int iHit = 0; iHit < collection->getNumberOfElements(); iHit++ ) {
    auto hit = static_cast<TrackerHitImpl*>( collection->getElementAt(iHit) );
    auto sensorID = hitCellDecoder(hit)["sensorID"];
    auto hitPosition = hit->getPosition();

    int cluX = 0, cluY = 0;
    float locX = 0, locY = 0;

    EUTelTripletGBLUtility::hit newhit(hitPosition, sensorID);

    auto rawData = static_cast<TrackerDataImpl*>(hit->getRawHits()[0]);
    if( hit->getType() == kEUTelSparseClusterImpl ){
      auto const & cluster = EUTelSparseClusterImpl<EUTelGenericSparsePixel>(rawData);
      cluster.getClusterSize(cluX, cluY);
      cluster.getCenterOfGravity(locX, locY);      
  
      newhit.locx = locX;
      newhit.locy = locY;
      newhit.clustersize = cluster.size();
      newhit.clustersizex = cluX;
      newhit.clustersizey = cluY;

      auto& resVecX = _planeResolutionX[sensorID];
      auto& resVecY = _planeResolutionY[sensorID];
      //The resolution vectors hold data like: <average> <size 1> <size 2> ... <size n> <greater than n> 
      //newhit.ex = (cluX >= resVecX.size()) ? resVecX.back() : resVecX[cluX];
      newhit.ex = resVecX.front();
      //newhit.ey = (cluY >= resVecY.size()) ? resVecY.back() : resVecY[cluY];
      newhit.ey = resVecY.front(); 
    } else {
      //throw cluster error
    }

    if(sensorID <= 6) {
      hits.emplace_back(newhit);
    } else {
      DUThits.emplace_back(newhit);
    }
  } // end loop over all hits in given collection

  streamlog_out(DEBUG4) << "Event " << event->getEventNumber() << " contains " << hits.size() << " hits" << std::endl;

  // Fill the telescope plane correlation plots:
  TelescopeCorrelationPlots(hits);
  // -----------------------------------------------------------------
  //
  // a word on gblutil:
  // for future use, the (contructor of?) gblutil could know about the AIDAProcessor, so histos dont need to be passed, but live in the Util class.
  // This opens the possiblity to create histograms within methods of the Util class, rather then uglily hacking them in the TripletGBL class.
  // DONE
  // Cut values could be passed and stored as private member and used in the methods, rather than passed during fct call.
  // tbd ...

  // -----------------------------------------------------------------
  // Downstream Telescope Triplets ("driplets")
  // Generate new triplet set for the Telescope Downstream Arm:
  std::vector<EUTelTripletGBLUtility::triplet> downstream_triplets;
  gblutil.FindTriplets(hits, 3, 4, 5, _triplet_res_cut, _slope_cut, downstream_triplets);
  streamlog_out(DEBUG4) << "Found " << downstream_triplets.size() << " driplets." << endl;

  // Iterate over all found downstream triplets to fill histograms and match them to the REF and DUT:
  for( auto& drip : downstream_triplets ){
    // Fill some histograms for downstream triplets:
    dridxHisto->fill( drip.getdx(4)*1E3 ); 
    dridxvsx->fill(  drip.base().x, drip.getdx(4)*1E3 ); // check for rot
    dridxvsy->fill( drip.base().y, drip.getdx(4)*1E3 );
    dridxvstx->fill( drip.slope().x*1E3, drip.getdx(4)*1E3 ); // check for z shift
    dridxvsty->fill( drip.slope().y*1E3, drip.getdx(4)*1E3 );

    dridyHisto->fill( drip.getdy(4)*1E3 );
    dridyvsx->fill( drip.base().x, drip.getdy(4)*1E3 );
    dridyvsy->fill( drip.base().y, drip.getdy(4)*1E3 );
    dridyvstx->fill( drip.slope().x*1E3, drip.getdy(4)*1E3 );
    dridyvsty->fill( drip.slope().y*1E3, drip.getdy(4)*1E3 );

    drixHisto->fill( -drip.gethit(4).x ); // -x = x_DP = out
    driyHisto->fill( -drip.gethit(4).y ); // -y = y_DP =  up
    drixyHisto->fill( -drip.gethit(4).x, -drip.gethit(4).y );
    dritxHisto->fill( drip.slope().x*1E3 );
    drityHisto->fill( drip.slope().y*1E3 );

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
  gblutil.FindTriplets(hits, 0, 1, 2, _triplet_res_cut, _slope_cut, upstream_triplets);
  streamlog_out(DEBUG4) << "Found " << upstream_triplets.size() << " triplets." << endl;

  // Iterate over all found upstream triplets to fill histograms and match them to the REF and DUT:
  for( auto& trip : upstream_triplets ) {
    // Fill some histograms for the upstream triplets:
    tridxHisto->fill( trip.getdx(1)*1E3 );
    tridyHisto->fill( trip.getdy(1)*1E3 );
    tridx1Histo->fill( trip.getdx(1)*1E0 );
    tridy1Histo->fill( trip.getdy(1)*1E0 );
    tridxvsx->fill( trip.base().x, trip.getdx(1)*1E3 ); // check for rot
    tridxvsy->fill( trip.base().y, trip.getdx(1)*1E3 );
    tridxvstx->fill( trip.slope().x*1E3, trip.getdx(1)*1E3 ); // check for z shift
    tridxvsty->fill( trip.slope().y*1E3, trip.getdx(1)*1E3 );
    tridyvsx->fill( trip.base().x, trip.getdy(1)*1E3 );
    tridyvsy->fill( trip.base().y, trip.getdy(1)*1E3 );
    tridyvstx->fill( trip.slope().x*1E3, trip.getdy(1)*1E3 );
    tridyvsty->fill( trip.slope().y*1E3, trip.getdy(1)*1E3 );
    trixHisto->fill( -trip.gethit(1).x );
    triyHisto->fill( -trip.gethit(1).y );
    trixyHisto->fill( -trip.gethit(1).x, -trip.gethit(1).y );
    tritxHisto->fill( trip.slope().x*1E3 );
    trityHisto->fill( trip.slope().y*1E3 );

    // Extrapolate Upstream triplet to Downstream planes 3,4,5: Resolution studies
    for( auto& lhit : hits ){
      if( lhit.plane <= 2 ) continue; // want 3,4, or 5

      // Fill residuals of triplet and hit in the selected plane:
      if( lhit.plane == 3 ) {
      	tridx3Histo->fill( trip.getdx(lhit)*1E0 );
      	tridy3Histo->fill( trip.getdy(lhit)*1E0 ); // 65 um at 4.7 GeV with CMS
      	tridx3bHisto->fill( trip.getdx(lhit)*1E3 ); // finer binning
      	tridy3bHisto->fill( trip.getdy(lhit)*1E3 ); // 
      } else if( lhit.plane == 4 ) {
      	tridx4Histo->fill( trip.getdx(lhit)*1E0 );
       	tridy4Histo->fill( trip.getdy(lhit)*1E0 ); //174 um at 4.7 GeV
      	tridx4bHisto->fill( trip.getdx(lhit)*1E3 ); // finer binning
      	tridy4bHisto->fill( trip.getdy(lhit)*1E3 ); // 
      } else if( lhit.plane == 5 ) {
      	tridx5Histo->fill( trip.getdx(lhit)*1E0 );
      	tridy5Histo->fill( trip.getdy(lhit)*1E0 ); //273 um at 4.7 GeV
      	tridx5bHisto->fill( trip.getdx(lhit)*1E3 ); // finer binning
      	tridy5bHisto->fill( trip.getdy(lhit)*1E3 ); // 
      }
    }// Resolution studies
  }// iterate over upstream triplets
  ntriHisto->fill( upstream_triplets.size() );

  //----------------------------------------------------------------------------
  // calculate efficiency of plane 3 by forming a triplet from planes 0, 1, 2; 2, 4, 5.
  // Then try to find a match on DUT (plane 3)

  // This scales the radius with beam energy spacing (efficiencty shouldnt depend on amount of scattering!). Define radius at 6 GeV, 20 mm
  double eff_radius = _eff_radius * 6. / _eBeam * (_planePosition[1] - _planePosition[0]) / 20.; 
  double track_match_z = _planePosition[3];
  double DUTz = _planePosition[3];

  // Generate new triplet set with planes 0, 1, 2; 2,4,5:
  std::vector<EUTelTripletGBLUtility::triplet> eff_triplets_UP = upstream_triplets;
  //gblutil.FindTriplets(hits, 0, 1, 2, _triplet_res_cut, _slope_cut, eff_triplets_UP);

  std::vector<EUTelTripletGBLUtility::triplet> eff_triplets_DOWN;
  gblutil.FindTriplets(hits, 2, 4, 5, _triplet_res_cut, _slope_cut, eff_triplets_DOWN);

  std::vector<AIDA::IProfile1D*> profiles;
  profiles.push_back(effix3);
  profiles.push_back(effiy3);
  gblutil.PlaneEfficiency(eff_triplets_UP, eff_triplets_DOWN, hits, 3, track_match_z, DUTz, _track_match_cut, eff_radius, profiles);
  
  //----------------------------------------------------------------------------
  // calculate efficiency of plane 2 by forming a triplet from planes 0, 1, 3; 3, 4, 5.
  // Then try to find a match on DUT (plane 2)
  track_match_z = _planePosition[3];
  DUTz = _planePosition[2];
  // This scales the radius with beam energy spacing (efficiencty shouldnt depend on amount of scattering!). Define radius at 6 GeV, 20 mm

  // Generate new triplet set with planes 0, 1, 3; 3,4,5:
  gblutil.FindTriplets(hits, 0, 1, 3, _triplet_res_cut, _slope_cut, eff_triplets_UP);
  // use existing one for down stream
  eff_triplets_DOWN = downstream_triplets;

  profiles.clear();
  profiles.push_back(effix2);
  profiles.push_back(effiy2);
  gblutil.PlaneEfficiency(eff_triplets_UP, eff_triplets_DOWN, hits, 2, track_match_z, DUTz, _track_match_cut, eff_radius, profiles);

  //----------------------------------------------------------------------------
  // calculate efficiency of plane 1 by forming a triplet from planes 0, 2, 3; 3, 4, 5.
  // Then try to find a match on DUT (plane 1)
  track_match_z = _planePosition[3];
  DUTz = _planePosition[1];
  // This scales the radius with beam energy spacing (efficiencty shouldnt depend on amount of scattering!). Define radius at 6 GeV, 20 mm

  // Generate new triplet set with planes 0, 2, 3; 3,4,5:
  gblutil.FindTriplets(hits, 0, 2, 3, _triplet_res_cut, _slope_cut, eff_triplets_UP);
  // use existing one for down stream
  eff_triplets_DOWN = downstream_triplets;

  profiles.clear();
  profiles.push_back(effix1);
  profiles.push_back(effiy1);
  gblutil.PlaneEfficiency(eff_triplets_UP, eff_triplets_DOWN, hits, 1, track_match_z, DUTz, _track_match_cut, eff_radius, profiles);

  //----------------------------------------------------------------------------
  // calculate efficiency of plane 0 by forming a triplet from planes 1, 2, 3; 3, 4, 5.
  // Then try to find a match on DUT (plane 0)
  track_match_z = _planePosition[3];
  DUTz = _planePosition[0];
  // This scales the radius with beam energy spacing (efficiencty shouldnt depend on amount of scattering!). Define radius at 6 GeV, 20 mm

  // Generate new triplet set with planes 1, 2, 3; 3,4,5:
  gblutil.FindTriplets(hits, 1, 2, 3, _triplet_res_cut, _slope_cut, eff_triplets_UP);
  eff_triplets_DOWN = downstream_triplets;

  profiles.clear();
  profiles.push_back(effix0);
  profiles.push_back(effiy0);
  gblutil.PlaneEfficiency(eff_triplets_UP, eff_triplets_DOWN, hits, 0, track_match_z, DUTz, _track_match_cut, eff_radius, profiles);

  //----------------------------------------------------------------------------
  // calculate efficiency of plane 4 by forming a triplet from planes 0, 1, 2; 2, 3, 5.
  // Then try to find a match on DUT (plane 3)
  track_match_z = _planePosition[3];
  DUTz = _planePosition[4];

  // Generate new triplet set with planes 0, 1, 2; 2,4,5:
  eff_triplets_UP = upstream_triplets;
  gblutil.FindTriplets(hits, 2, 3, 5, _triplet_res_cut, _slope_cut, eff_triplets_DOWN);

  profiles.clear();
  profiles.push_back(effix4);
  profiles.push_back(effiy4);
  gblutil.PlaneEfficiency(eff_triplets_UP, eff_triplets_DOWN, hits, 4, track_match_z, DUTz, _track_match_cut, eff_radius, profiles);

  //----------------------------------------------------------------------------
  // calculate efficiency of plane 5 by forming a triplet from planes 0, 1, 2; 2, 3, 4.
  // Then try to find a match on DUT (plane 3)
  track_match_z = _planePosition[3];
  DUTz = _planePosition[5];

  // Generate new triplet set with planes 0, 1, 2; 2,3,4:
  eff_triplets_UP = upstream_triplets;
  gblutil.FindTriplets(hits, 2, 3, 4, _triplet_res_cut, _slope_cut, eff_triplets_DOWN);

  profiles.clear();
  profiles.push_back(effix5);
  profiles.push_back(effiy5);
  gblutil.PlaneEfficiency(eff_triplets_UP, eff_triplets_DOWN, hits, 5, track_match_z, DUTz, _track_match_cut, eff_radius, profiles);

  //----------------------------------------------------------------------------
  // six: triplets A and driplets B
  // matching and GBL fit
  // kinks: triplets A vs driplets B
  // scattering point = DUT:

  //FIXME: DUTz is the place where they should be matched!
  // now hardcoded to center of telescope...
  DUTz = _planePosition[2] + (_planePosition[3] - _planePosition[2])/2;

  // Match the Telescope Upstream and Downstream Arm triplets to get tracks:
  std::vector<EUTelTripletGBLUtility::track> telescope_tracks;
  gblutil.MatchTriplets(upstream_triplets,downstream_triplets, DUTz, _track_match_cut, telescope_tracks);

  //delete downstream_triplets;
  //delete upstream_triplets;

  streamlog_out(DEBUG4) << "Found " << telescope_tracks.size() << " tracks from matching t/driplets." << endl;
  
  for( auto& tr: telescope_tracks ){
    auto const & trip = tr.get_upstream();
    auto const & drip = tr.get_downstream();
    EUTelTripletGBLUtility::triplet srip(tr.gethit(0), tr.gethit(2), tr.gethit(5)); // seed triplet is called srip

    std::vector<double> xAplanes(_nPlanes);
    std::vector<double> yAplanes(_nPlanes);

    for (int i = 0; i < _nPlanes; i++){
      xAplanes[i] = srip.getx_at(_planePosition[i]);
      yAplanes[i] = srip.gety_at(_planePosition[i]);
    }

    // Track kinks as difference in triplet slopes:
    //      double kx = drip.slope().x - trip.slope().x; //kink
    //      double ky = drip.slope().y - trip.slope().y;
    double kx = tr.kink_x();
    double ky = tr.kink_y();

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

    Eigen::Matrix2d proL2m = Eigen::Matrix2d::Identity();

    //double res = 3.42E-3; // [mm] Anemone telescope intrinsic resolution
    //res = 4.5E-3; // EUDET

    // scatter:
    Eigen::Vector2d scat = Eigen::Vector2d::Zero(); //mean is zero
    
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
    s = 0.;

    int DUT_label;
  
    for( int ipl = 0; ipl < 6; ++ipl ){

      // Get the corresponding hit from the track:
      EUTelTripletGBLUtility::hit trackhit = tr.gethit(ipl);

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
      } else if( ipl == 1 ) {
        sixx1Histo->fill( -trackhit.x );
        sixy1Histo->fill( -trackhit.y );
      } else if( ipl == 2 ) {
        sixx2Histo->fill( -trackhit.x );
        sixy2Histo->fill( -trackhit.y );
      } else if( ipl == 3 ) {
        sixx3Histo->fill( -trackhit.x );
        sixy3Histo->fill( -trackhit.y );
      } else if( ipl == 4 ) {
        sixx4Histo->fill( -trackhit.x );
        sixy4Histo->fill( -trackhit.y );
      } else if( ipl == 5 ) {
        sixx5Histo->fill( -trackhit.x );
        sixy5Histo->fill( -trackhit.y );
      }

      auto point = gbl::GblPoint( gblutil.JacobianPointToPoint( step ) );

      Eigen::Vector2d meas;
      meas[0] = rx[ipl];
      meas[1] = ry[ipl];
      //meas[0] = trackhit.x;
      //meas[1] = trackhit.y;

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
      
   	  Eigen::Vector2d measPrec; // precision = 1/resolution^2
      measPrec[0] = 1.0 / _resolution_tmp / _resolution_tmp;
      measPrec[1] = 1.0 / _resolution_tmp / _resolution_tmp;

      if(ipl != _dut_plane) { point.addMeasurement( proL2m, meas, measPrec ); }

      point.addScatterer( scat, _planeWscatSi[ipl] );

      std::vector<int> globalLabels(3);
      globalLabels[0] = 10 + ipl; // dx
      globalLabels[1] = 20 + ipl; // dy
      globalLabels[2] = 40 + ipl; // rot

      traj_points.emplace_back(point);
      s += step;
      sPoint.push_back( s );
      DUT_label = sPoint.size();
      ilab.push_back(DUT_label);

      if( ipl < 5) {
        double distplane = _planePosition[ipl+1] - _planePosition[ipl];
        step = 0.21*distplane;

        auto point = gbl::GblPoint( gblutil.JacobianPointToPoint( step ) );
        point.addScatterer( scat, _planeWscatAir[ipl] );
        s += step;
        traj_points.push_back(point);
        sPoint.push_back( s );

        step = 0.58*distplane;
        auto point1 = gbl::GblPoint( gblutil.JacobianPointToPoint( step ) );
        point1.addScatterer( scat, _planeWscatAir[ipl] );
        s += step;
        traj_points.push_back(point1);
        sPoint.push_back( s );
	      
        step = 0.21*distplane; // remaing distance to next plane
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
    if(_printEventCounter < 10){
      streamlog_out(MESSAGE4) << "traj with " << traj.getNumPoints() << " points:" << endl;
      for( unsigned int ipl = 0; ipl < sPoint.size(); ++ipl ){
        streamlog_out(DEBUG4) << "  GBL point " << ipl;
        streamlog_out(DEBUG4) << "  z " << sPoint[ipl]; 
        streamlog_out(DEBUG4) << endl;
      }
      for( unsigned int ipl = 0; ipl < 6; ++ipl ){
        streamlog_out(DEBUG4) << " plane " << ipl << ", lab " << ilab[ipl];
        streamlog_out(DEBUG4) << " z " << sPoint[ilab[ipl]-1];
        streamlog_out(DEBUG4) << "  dx " << rx[ipl];
        streamlog_out(DEBUG4) << "  dy " << ry[ipl];
        streamlog_out(DEBUG4) << "  chi2 " << Chi2;
        streamlog_out(DEBUG4) << "  ndf " << Ndf;
        streamlog_out(DEBUG4) << endl;
      }

      streamlog_out(DEBUG4)  << " Is traj valid? " << traj.isValid() << std::endl;
      //traj.printPoints();
      //traj.printTrajectory();
      //traj.printData();
      _printEventCounter++;
    }

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
      
      _ngbl++;

      Eigen::VectorXd aCorrection(2);
      Eigen::MatrixXd aCovariance(5,5);

      double ax[8];
      // double ay[8];
      unsigned int k = 0;

      unsigned int ndim = 2;
      Eigen::VectorXd aResiduals(ndim);
      Eigen::VectorXd aMeasErrors(ndim);
      Eigen::VectorXd aResErrors(ndim);
      Eigen::VectorXd aDownWeights(ndim);

      Eigen::VectorXd aKinks(ndim);
      Eigen::VectorXd aKinkErrors(ndim);
      Eigen::VectorXd kResErrors(ndim);
      Eigen::VectorXd kDownWeights(ndim);

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
      int CS0 = tr.gethit(0).clustersize;

      double modPos[2];
      modPos[0] = ((corrPos[0])+(invsignx*(abs(nx) +.5) + 0.5) *pixel_size)*1e3;
      modPos[1] = ((corrPos[1])+(invsigny*(abs(ny) +.5) + 0.5) *pixel_size)*1e3;
      
      // overlay of all CSs
      gblnxy_plane0->fill(modPos[0], modPos[1], CS0 );
      gblnxy1_plane0->fill(modPos[0], modPos[1], 1 );

      if (CS0 == 1){
	gblnCS1xy_plane0->fill(modPos[0], modPos[1] );
	//gblnCS1xy1_plane0->fill((trackhitx[3] - aResiduals[0]), (trackhity[3] - aResiduals[1]));

        gblrxvsxpix01_CS1->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])*1e3); 
        gblryvsypix01_CS1->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])*1e3); 
      } 
      if (CS0 == 2){
	gblnCS2xy_plane0->fill(modPos[0], modPos[1] );
	// gblnCS2xy1_plane0->fill((trackhitx[3] - aResiduals[0]),   (trackhity[3] - aResiduals[1]));
	
        gblrxvsxpix01_CS2->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])*1e3); 
        gblryvsypix01_CS2->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])*1e3); 
      } 
      if (CS0 == 3){
	gblnCS3xy_plane0->fill(modPos[0], modPos[1] );
	//gblnCS3xy1_plane0->fill((trackhitx[3] - aResiduals[0]),(trackhity[3] - aResiduals[1]));

        gblrxvsxpix01_CS3->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])*1e3); 
        gblryvsypix01_CS3->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])*1e3); 
      } 
      if (CS0 == 4){
	gblnCS4xy_plane0->fill(modPos[0], modPos[1] );
	//gblnCS4xy1_plane0->fill((trackhitx[3] - aResiduals[0]),   (trackhity[3] - aResiduals[1]));

        gblrxvsxpix01_CS4->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])*1e3); 
        gblryvsypix01_CS4->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])*1e3); 
      } 
      if( CS0 == 5 ) {
	gblnCS5xy_plane0->fill(modPos[0], modPos[1] );
	//gblnCS5xy1_plane0->fill((trackhitx[3] - aResiduals[0]), (trackhity[3] - aResiduals[1]));
      }
      if( CS0 == 6 ) {
	gblnCS6xy_plane0->fill(modPos[0], modPos[1] );
	//gblnCS6xy1_plane0->fill((trackhitx[3] - aResiduals[0]),  (trackhity[3] - aResiduals[1]));
      }
      if( CS0 > 6 ) {
	gblnCS7xy_plane0->fill(modPos[0], modPos[1] );
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
      //
      // FIXME get alignment constants from GEAR file (rather then opening all 4 alignment files ....)
      int CS3 = tr.gethit(3).clustersize;
      corrPos[0] = trackhitx[3] - aResiduals[0] +  7.077e-3 + 0.35e-3; // run000117, better get automatically, with new GEAR file
      corrPos[1] = trackhity[3] - aResiduals[1] + 18.128e-3 + 0.1e-3;  // run000117
      //corrPos[0] = trackhitx[3] - aResiduals[0] -  6.645e-3; // run002140
      //corrPos[1] = trackhity[3] - aResiduals[1] - 10.324e-3;  // run002140

      // correct for rotation
      double corrPos2[2];
      double gamma3 = 5.3e-3; // run000117, get automatically
      //double gamma3 = 1.94969e-03; // run002140
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
      gblrxvsxpix31->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])*1e3); 
      gblryvsypix31->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])*1e3); 

      // overlay of all CSs
      gblnxy_plane3->fill(modPos[0], modPos[1], CS3 );
      gblnxy1_plane3->fill(modPos[0], modPos[1], 1 );

      if( CS3 == 1 ) {
	gblrx3_cs1Histo->fill( aResiduals[0]* 1E3 ); 
	gblry3_cs1Histo->fill( aResiduals[1]* 1E3 ); 
	gblpx3_cs1Histo->fill( aResiduals[0] / aResErrors[0] ); 
	gblpy3_cs1Histo->fill( aResiduals[1] / aResErrors[1] ); 

	gblnCS1xy_plane3->fill(modPos[0], modPos[1] );
	gblnCS1xy1_plane3->fill((trackhitx[3] - aResiduals[0]),
	    (trackhity[3] - aResiduals[1]));

	gblrxvsxpix3cs1->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])*1e3); 
	gblryvsypix3cs1->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])*1e3); 
      }
      if( CS3 == 2 ) {
	gblrx3_cs2Histo->fill( aResiduals[0]* 1E3 ); 
	gblry3_cs2Histo->fill( aResiduals[1]* 1E3 ); 
	gblpx3_cs2Histo->fill( aResiduals[0] / aResErrors[0] ); 
	gblpy3_cs2Histo->fill( aResiduals[1] / aResErrors[1] ); 

	gblnCS2xy_plane3->fill(modPos[0], modPos[1] );
	gblnCS2xy1_plane3->fill((trackhitx[3] - aResiduals[0]),
	    (trackhity[3] - aResiduals[1]));

	gblrxvsxpix3cs2->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])*1e3); 
	gblryvsypix3cs2->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])*1e3); 
      }
      if( CS3 == 3 ) {
	gblrx3_cs3Histo->fill( aResiduals[0]* 1E3 ); 
	gblry3_cs3Histo->fill( aResiduals[1]* 1E3 ); 
	gblpx3_cs3Histo->fill( aResiduals[0] / aResErrors[0] ); 
	gblpy3_cs3Histo->fill( aResiduals[1] / aResErrors[1] ); 

	gblnCS3xy_plane3->fill(modPos[0], modPos[1] );
	gblnCS3xy1_plane3->fill((trackhitx[3] - aResiduals[0]),
	    (trackhity[3] - aResiduals[1]));

	gblrxvsxpix3cs3->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])*1e3); 
	gblryvsypix3cs3->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])*1e3); 
      }
      if( CS3 == 4 ) {
	gblrx3_cs4Histo->fill( aResiduals[0]* 1E3 ); 
	gblry3_cs4Histo->fill( aResiduals[1]* 1E3 ); 
	gblpx3_cs4Histo->fill( aResiduals[0] / aResErrors[0] ); 
	gblpy3_cs4Histo->fill( aResiduals[1] / aResErrors[1] ); 

	gblnCS4xy_plane3->fill(modPos[0], modPos[1] );
	gblnCS4xy1_plane3->fill((trackhitx[3] - aResiduals[0]),
	    (trackhity[3] - aResiduals[1]));

	gblrxvsxpix3cs4->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])*1e3); 
	gblryvsypix3cs4->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])*1e3); 
      }
      if( CS3 == 5 ) {
	gblrx3_cs5Histo->fill( aResiduals[0]* 1E3 ); 
	gblry3_cs5Histo->fill( aResiduals[1]* 1E3 ); 
	gblpx3_cs5Histo->fill( aResiduals[0] / aResErrors[0] ); 
	gblpy3_cs5Histo->fill( aResiduals[1] / aResErrors[1] ); 

	gblnCS5xy_plane3->fill(modPos[0], modPos[1] );
	gblnCS5xy1_plane3->fill((trackhitx[3] - aResiduals[0]),
	    (trackhity[3] - aResiduals[1]));

	gblrxvsxpix3cs5->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])*1e3); 
	gblryvsypix3cs5->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])*1e3); 
      }
      if( CS3 == 6 ) {
	gblrx3_cs6Histo->fill( aResiduals[0]* 1E3 ); 
	gblry3_cs6Histo->fill( aResiduals[1]* 1E3 ); 
	gblpx3_cs6Histo->fill( aResiduals[0] / aResErrors[0] ); 
	gblpy3_cs6Histo->fill( aResiduals[1] / aResErrors[1] ); 

	gblnCS6xy_plane3->fill(modPos[0], modPos[1] );
	gblnCS6xy1_plane3->fill((trackhitx[3] - aResiduals[0]),
	    (trackhity[3] - aResiduals[1]));

	gblrxvsxpix3cs6->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])*1e3); 
	gblryvsypix3cs6->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])*1e3); 
      }
      if( CS3 > 6 ) {
	gblpx3_cs7Histo->fill( aResiduals[0] / aResErrors[0] ); 
	gblpy3_cs7Histo->fill( aResiduals[1] / aResErrors[1] ); 

	gblnCS7xy_plane3->fill(modPos[0], modPos[1] );

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
      int CS5 = tr.gethit(5).clustersize;

      if (CS5 == 1){
        gblrxvsxpix51_CS1->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])*1e3); 
        gblryvsypix51_CS1->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])*1e3); 
      } 
      if (CS5 == 2){
        gblrxvsxpix51_CS2->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])*1e3); 
        gblryvsypix51_CS2->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])*1e3); 
      } 
      if (CS5 == 3){
        gblrxvsxpix51_CS3->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])*1e3); 
        gblryvsypix51_CS3->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])*1e3); 
      } 
      if (CS5 == 4){
        gblrxvsxpix51_CS4->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])*1e3); 
        gblryvsypix51_CS4->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])*1e3); 
      } 


      k++;
      // do some more analysis with kinks. 
      //   Calculate the two angles as corr_n-1 - corr_n-2, corr_n - corr_n-1, corr_n+1 - corr_n -> sum of all =  corr_n+1 - corr_n-2

      ipos = ilab[2]+1; // 1 point upstream of centre = 1 downstream of plane 2
      traj.getResults( ipos, aCorrection, aCovariance );
      double kink_upstream = aCorrection[1];

      ipos = ilab[3]-1; // 1 point downstream of centre = 1 upstream of plane 3
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
    EUTelTripletGBLUtility::hit intersect = tr.intersect();
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
void EUTelGBLFitter::check( LCEvent * /* evt */  ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


//------------------------------------------------------------------------------
void EUTelGBLFitter::end(){

  // Print the summary:
  streamlog_out(MESSAGE5)
    << "---------------------------------------------------------------------------------------------------------" << std::endl
    << std::endl
    << "Processed events:    "
    << std::setw(10) << std::setiosflags(std::ios::right)
    << _nEvt << std::resetiosflags(std::ios::right) << std::endl;

} // end end


void EUTelGBLFitter::TelescopeCorrelationPlots(std::vector<EUTelTripletGBLUtility::hit> const & telescopehits) {
  for( auto& ihit: telescopehits ){
    int ipl = ihit.plane;
    for( auto& jhit: telescopehits ){
      int jpl = jhit.plane;

      double dx = jhit.x - ihit.x;
      double dy = jhit.y - ihit.y;

      if( ipl == 0 ) {
        if( jpl == 1 ) {
	        dx01Histo->fill( dx );
	        dy01Histo->fill( dy );
	        if( abs(dy) < 1 ) 
            du01Histo->fill( dx );
	      } else if( jpl == 2 ) {
	        dx02Histo->fill( dx );
	      } else if( jpl == 3 ) {
          dx03Histo->fill( dx );
        } else if( jpl == 4 ) {
          dx04Histo->fill( dx );
        } else if( jpl == 5 ) {
          dx05Histo->fill( dx );
        }
      } else if( ipl == 1 ) {
        if( jpl == 2 ) {
	        dx12Histo->fill( dx );
	        dy12Histo->fill( dy );
	        if( abs(dy) < 1 ) 
            du12Histo->fill( dx );
	      }
      } else if( ipl == 2 ) {
        if( jpl == 3 ) {
	        dx23Histo->fill( dx );
	        dy23Histo->fill( dy );
	        if( abs(dy) < 1 )
            du23Histo->fill( dx );
	      }
      } else if( ipl == 3 ) {
	      if( jpl == 4 ) {
	        dx34Histo->fill( dx );
	        dy34Histo->fill( dy );
	        if( abs(dy) < 1 )
            du34Histo->fill( dx );
    	  }
      } else if( ipl == 4 ) {
        if( jpl == 5 ) {
          dx45Histo->fill( dx );
          dy45Histo->fill( dy );
          if( abs(dy) < 1 ) 
            du45Histo->fill( dx );
	      }
      }
    }
  }
}

//------------------------------------------------------------------------------
void EUTelGBLFitter::bookHistos()
{

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  // telescope hits per plane:
  AIDAProcessor::tree(this)->mkdir("Telescope");

  nAllHitHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/nallhit", 201, -0.5, 200.5 );
  nAllHitHisto->setTitle( "Telescope hits/event;telescope hits;events" );

  // telescope dx:

  dx01Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx01", 100, -1, 1 );
  dx01Histo->setTitle( "x1-x0;x_{1}-x_{0} [mm];hit pairs" );

  dy01Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dy01", 100, -1, 1 );
  dy01Histo->setTitle( "y1-y0;y_{1}-y_{0} [mm];hit pairs" );

  du01Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/du01", 100, -1, 1 );
  du01Histo->setTitle( "x1-x0, |dy| < 1;x_{1}-x_{0} [mm];hit pairs" );

  dx02Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx02", 100, -1, 1 );
  dx02Histo->setTitle( "x2-x0;x_{2}-x_{0} [mm];hit pairs" );

  dx03Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx03", 100, -1, 1 );
  dx03Histo->setTitle( "x3-x0;x_{3}-x_{0} [mm];hit pairs" );

  dx04Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx04", 100, -1, 1 );
  dx04Histo->setTitle( "x4-x0;x_{4}-x_{0} [mm];hit pairs" );

  dx05Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx05", 100, -1, 1 );
  dx05Histo->setTitle( "x5-x0;x_{5}-x_{0} [mm];hit pairs" );

  dx12Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx12", 100, -1, 1 );
  dx12Histo->setTitle( "x2-x1;x_{2}-x_{1} [mm];hit pairs" );

  dy12Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dy12", 100, -1, 1 );
  dy12Histo->setTitle( "y2-y1;y_{2}-y_{1} [mm];hit pairs" );

  du12Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/du12", 100, -1, 1 );
  du12Histo->setTitle( "x2-x1, |dy| < 1;x_{2}-x_{1} [mm];hit pairs" );

  dx23Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx23", 100, -1, 1 );
  dx23Histo->setTitle( "x3-x2;x_{3}-x_{2} [mm];hit pairs" );

  dy23Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dy23", 100, -1, 1 );
  dy23Histo->setTitle( "y3-y2;y_{3}-y_{2} [mm];hit pairs" );

  du23Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/du23", 100, -1, 1 );
  du23Histo->setTitle( "x3-x2, |dy| < 1;x_{3}-x_{2} [mm];hit pairs" );

  dx34Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx34", 100, -1, 1 );
  dx34Histo->setTitle( "x4-x3;x_{4}-x_{3} [mm];hit pairs" );

  dy34Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dy34", 100, -1, 1 );
  dy34Histo->setTitle( "y4-y3;y_{4}-y_{3} [mm];hit pairs" );

  du34Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/du34", 100, -1, 1 );
  du34Histo->setTitle( "x4-x3, |dy| < 1;x_{4}-x_{3} [mm];hit pairs" );

  dx45Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx45", 100, -1, 1 );
  dx45Histo->setTitle( "x5-x4;x_{5}-x_{4} [mm];hit pairs" );

  dy45Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dy45", 100, -1, 1 );
  dy45Histo->setTitle( "y5-y4;y_{5}-y_{4} [mm];hit pairs" );

  du45Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/du45", 100, -1, 1 );
  du45Histo->setTitle( "x5-x4, |dy| < 1;x_{5}-x_{4} [mm];hit pairs" );

  // triplets:
  AIDAProcessor::tree(this)->mkdir("Upstream");

  dzcvsxy = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "Upstream/dzcvsxy", 120, -12, 12, 60, -6, 6, -999, 999 );
  dzcvsxy->setTitle( "DUT plane;telescope track x_{DUT} [mm];telescope track y_{DUT} [mm];z_{DUT} [mm]" );

  z3vsxy = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "Upstream/z3vsxy", 120, -12, 12, 60, -6, 6, -999, 999 );
  z3vsxy->setTitle( "DUT plane;telescope track x_{DUT} [mm];telescope track y_{DUT} [mm];z_{DUT} [mm]" );

  tridxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx", 100, -100, 100 );
  tridxHisto->setTitle( "triplet dx;x_{1}-x_{m} [#mum];telescope triplets" );

  tridyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy", 100, -100, 100 );
  tridyHisto->setTitle( "triplet dy;y_{1}-y_{m} [#mum];telescope triplets" );

  tridx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx1", 100, -1, 1 );
  tridx1Histo->setTitle( "triplet dx;x_{1}-x_{t} [mm];telescope triplets" );

  tridy1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy1", 100, -1, 1 );
  tridy1Histo->setTitle( "triplet dy;y_{1}-y_{t} [mm];telescope triplets" );

  tridxvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridxvsx", 100, -10, 10, -100, 100 );
  tridxvsx->setTitle( "triplet x resid vs x;x [mm];triplet <#Deltax> [#mum]" );

  tridxvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridxvsy", 50, -5, 5, -100, 100 );
  tridxvsy->setTitle( "triplet x resid vs y;y [mm];triplet <#Deltax> [#mum]" );

  tridxvstx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridxvstx", 80, -2, 2, -100, 100 );
  tridxvstx->setTitle( "triplet x resid vs tx;t_{x} [mrad];triplet <#Deltax> [#mum]" );

  tridxvsty = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridxvsty", 80, -2, 2, -100, 100 );
  tridxvsty->setTitle( "triplet x resid vs ty;t_{y} [mrad];triplet <#Deltax> [#mum]" );

  tridyvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridyvsx", 100, -10, 10, -100, 100 );
  tridyvsx->setTitle( "triplet y resid vs x;x [mm];triplet <#Deltay> [#mum]" );

  tridyvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridyvsy", 50, -5, 5, -100, 100 );
  tridyvsy->setTitle( "triplet y resid vs y;y [mm];triplet <#Deltay> [#mum]" );

  tridyvstx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridyvstx", 80, -2, 2, -100, 100 );
  tridyvstx->setTitle( "triplet y resid vs tx;t_{x} [mrad];triplet <#Deltay> [#mum]" );

  tridyvsty = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridyvsty", 80, -2, 2, -100, 100 );
  tridyvsty->setTitle( "triplet y resid vs ty;t_{y} [mrad];triplet <#Deltay> [#mum]" );

  tridx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx3", 100, -1, 1 );
  tridx3Histo->setTitle( "triplet dx;x_{3}-x_{t} [mm];telescope triplets" );

  tridy3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy3", 100, -1, 1 );
  tridy3Histo->setTitle( "triplet dy;y_{3}-y_{t} [mm];telescope triplets" );

  tridx3bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx3b", 100, -250, 250 );
  tridx3bHisto->setTitle( "triplet dx;x_{3}-x_{t} [um];telescope triplets" );

  tridy3bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy3b", 100, -250, 250 );
  tridy3bHisto->setTitle( "triplet dy;y_{3}-y_{t} [um];telescope triplets" );

  tridx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx4", 100, -2, 2 );
  tridx4Histo->setTitle( "triplet dx;x_{4}-x_{t} [mm];telescope triplets" );

  tridy4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy4", 100, -2, 2 );
  tridy4Histo->setTitle( "triplet dy;y_{4}-y_{t} [mm];telescope triplets" );

  tridx4bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx4b", 100, -400, 400 );
  tridx4bHisto->setTitle( "triplet dx;x_{4}-x_{t} [um];telescope triplets" );

  tridy4bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy4b", 100, -400, 400 );
  tridy4bHisto->setTitle( "triplet dy;y_{4}-y_{t} [um];telescope triplets" );

  tridx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx5", 100, -3, 3 );
  tridx5Histo->setTitle( "triplet dx;x_{5}-x_{t} [mm];telescope triplets" );

  tridy5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy5", 100, -3, 3 );
  tridy5Histo->setTitle( "triplet dy;y_{5}-y_{t} [mm];telescope triplets" );

  tridx5bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx5b", 100, -1000, 1000 );
  tridx5bHisto->setTitle( "triplet dx;x_{5}-x_{t} [um];telescope triplets" );

  tridy5bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy5b", 100, -1000, 1000 );
  tridy5bHisto->setTitle( "triplet dy;y_{5}-y_{t} [um];telescope triplets" );

  trixHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/trix", 240, -12, 12 );
  trixHisto->setTitle( "triplet x1;x1_{out} [mm];telescope triplets" );

  triyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/triy", 120, -6, 6 );
  triyHisto->setTitle( "triplet y1;y1_{up} [mm];telescope triplets" );

  trixyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "Upstream/trixy", 240, -12, 12, 120, -6, 6 );
  trixyHisto->setTitle( "triplet y1 vs x1;x1_{out} [mm];y1_{up} [mm];telescope triplets" );

  tritxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tritx", 100, -10, 10 );
  tritxHisto->setTitle( "triplet slope x;#theta_{x} [mrad];telescope triplets" );

  trityHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/trity", 100, -10, 10 );
  trityHisto->setTitle( "triplet slope y;#theta_{y} [mrad];telescope triplets" );

  trixdutHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/trixdut", 240, -12, 12 );
  trixdutHisto->setTitle( "triplet at DUT;triplet x_{out} at DUT [mm];telescope triplets" );

  triydutHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/triydut", 120, -6, 6 );
  triydutHisto->setTitle( "triplet at DUT;triplet y_{up} at DUT [mm];telescope triplets" );

  trixydutHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "Upstream/trixydut", 240, -12, 12, 120, -6, 6 );
  trixydutHisto->setTitle( "triplet at DUT;triplet x_{out} at DUT [mm];triplet y_{up} at DUT [mm];telescope triplets" );

  ntriHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/ntri", 31, -0.5, 30.5 );
  ntriHisto->setTitle( "telescope triplets;0-1-2 triplets;events" );


  // driplets:
  AIDAProcessor::tree(this)->mkdir("Downstream");

  dridxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/dridx", 100, -100, 100 );
  dridxHisto->setTitle( "driplet dx;x_{4}-x_{m} [#mum];driplets" );

  dridyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/dridy", 100, -100, 100 );
  dridyHisto->setTitle( "driplet dy;y_{4}-y_{m} [#mum];driplets" );

  dridxvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridxvsx", 100, -10, 10, -100, 100 );
  dridxvsx->setTitle( "driplet x resid vs x;x [mm];driplet <#Deltax> [#mum]" );

  dridxvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridxvsy", 50, -5, 5, -100, 100 );
  dridxvsy->setTitle( "driplet x resid vs y;y [mm];driplet <#Deltax> [#mum]" );

  dridxvstx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridxvstx", 80, -2, 2, -100, 100 );
  dridxvstx->setTitle( "driplet x resid vs tx;t_{x} [mrad];driplet <#Deltax> [#mum]" );

  dridxvsty = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridxvsty", 80, -2, 2, -100, 100 );
  dridxvsty->setTitle( "driplet x resid vs ty;t_{y} [mrad];driplet <#Deltax> [#mum]" );

  dridyvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridyvsx", 100, -10, 10, -100, 100 );
  dridyvsx->setTitle( "driplet y resid vs x;x [mm];driplet <#Deltay> [#mum]" );

  dridyvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridyvsy", 50, -5, 5, -100, 100 );
  dridyvsy->setTitle( "driplet y resid vs y;y [mm];driplet <#Deltay> [#mum]" );

  dridyvstx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridyvstx", 80, -2, 2, -100, 100 );
  dridyvstx->setTitle( "driplet y resid vs tx;t_{x} [mrad];driplet <#Deltay> [#mum]" );

  dridyvsty = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridyvsty", 80, -2, 2, -100, 100 );
  dridyvsty->setTitle( "driplet y resid vs ty;t_{y} [mrad];driplet <#Deltay> [#mum]" );

  drixHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/drix", 240, -12, 12 );
  drixHisto->setTitle( "driplet x4;x4_{out} [mm];driplets" );

  driyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/driy", 120, -6, 6 );
  driyHisto->setTitle( "driplet y4;y4_{up} [mm];driplets" );

  drixyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "Downstream/drixy", 240, -12, 12, 120, -6, 6 );
  drixyHisto->setTitle( "driplet y4 vs x4;x4_{out} [mm];y4_{up} [mm];driplets" );

  dritxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/dritx", 100, -10, 10 );
  dritxHisto->setTitle( "driplet slope x;#theta_{x} [mrad];driplets" );

  drityHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/drity", 100, -10, 10 );
  drityHisto->setTitle( "driplet slope y;#theta_{y} [mrad];driplets" );



  bacsxaHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "bacsxa", 100, -10, 10 );
  bacsxaHisto->setTitle( "DUT + driplet x;DUT cluster + driplet #Sigmax [mm];DUT clusters" );

  bacdyaHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "bacdya", 100, -10, 10 );
  bacdyaHisto->setTitle( "DUT - driplet y;DUT cluster - driplet #Deltay [mm];DUT clusters" );

  bacsxcHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "bacsxc", 200, -500, 500 );
  bacsxcHisto->setTitle( "DUT + driplet x;DUT cluster + driplet #Sigmax [mm];DUT clusters" );

  bacdycHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "bacdyc", 200, -500, 500 );
  bacdycHisto->setTitle( "DUT - driplet y;DUT cluster - driplet #Deltay [mm];DUT clusters" );

  bacsxcqHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "bacsxcq", 200, -500, 500 );
  bacsxcqHisto->setTitle( "DUT + driplet x;DUT cluster + driplet #Sigmax [mm];DUT clusters" );

  bacdycqHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "bacdycq", 100, -500, 500 );
  bacdycqHisto->setTitle( "DUT - driplet y;DUT cluster - driplet #Deltay [mm];DUT clusters" );


  ndriHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/ndri", 31, -0.5, 30.5 );
  ndriHisto->setTitle( "telescope driplets;3-4-5 driplets;events" );

  // efficiency plane 3 with triplet 024
  AIDAProcessor::tree(this)->mkdir("Effi");

 effix0 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Effi/effix0", 120, -12, 12, -0.1, 1.1 );
  effix0->setTitle( "trip-effi vs x at plane 0;x [mm]; efficiency" );

  effiy0 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Effi/effiy0", 60, -6, 6, -0.1, 1.1 );
  effiy0->setTitle( "trip-effi vs y at plane 0;y [mm]; efficiency" );

 effix1 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Effi/effix1", 120, -12, 12, -0.1, 1.1 );
  effix1->setTitle( "trip-effi vs x at plane 1;x [mm]; efficiency" );

  effiy1 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Effi/effiy1", 60, -6, 6, -0.1, 1.1 );
  effiy1->setTitle( "trip-effi vs y at plane 1;y [mm]; efficiency" );


  effix2 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Effi/effix2", 120, -12, 12, -0.1, 1.1 );
  effix2->setTitle( "trip-effi vs x at plane 2;x [mm]; efficiency" );

  effiy2 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Effi/effiy2", 60, -6, 6, -0.1, 1.1 );
  effiy2->setTitle( "trip-effi vs y at plane 2;y [mm]; efficiency" );

  effix3 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Effi/effix3", 120, -12, 12, -0.1, 1.1 );
  effix3->setTitle( "trip-effi vs x at plane 3;x [mm]; efficiency" );

  effiy3 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Effi/effiy3", 60, -6, 6, -0.1, 1.1 );
  effiy3->setTitle( "trip-effi vs y at plane 3;y [mm]; efficiency" );

  effix4 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Effi/effix4", 120, -12, 12, -0.1, 1.1 );
  effix4->setTitle( "trip-effi vs x at plane 4;x [mm]; efficiency" );

  effiy4 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Effi/effiy4", 60, -6, 6, -0.1, 1.1 );
  effiy4->setTitle( "trip-effi vs y at plane 4;y [mm]; efficiency" );

  effix5 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Effi/effix5", 120, -12, 12, -0.1, 1.1 );
  effix5->setTitle( "trip-effi vs x at plane 5;x [mm]; efficiency" );

  effiy5 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Effi/effiy5", 60, -6, 6, -0.1, 1.1 );
  effiy5->setTitle( "trip-effi vs y at plane 5;y [mm]; efficiency" );





  //driplets-triplets
  // Tracks
  AIDAProcessor::tree(this)->mkdir("Tracks");

  nsixHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/nsix", 21, -0.5, 20.5 );
  nsixHisto->setTitle( "telescope six-plane-tracks;six-plane-tracks;events" );


  clustersize0 = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/clustersize0", 100, 0, 100 );
  clustersize0->setTitle( "clustersize at 0; clustersize;six-plane tracks" );

  sixx0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx0", 240, -12, 12 );
  sixx0Histo->setTitle( "six x at 0;six x_{out} at 0 [mm];six-plane tracks" );

  sixy0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy0", 120, -6, 6 );
  sixy0Histo->setTitle( "six y at 0;six y_{up} at 0 [mm];six-plane tracks" );

  clustersize1 = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/clustersize1", 100, 0, 100 );
  clustersize1->setTitle( "clustersize at 1; clustersize;six-plane tracks" );

  sixx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx1", 240, -12, 12 );
  sixx1Histo->setTitle( "six x at 1;six x_{out} at 1 [mm];six-plane tracks" );

  sixy1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy1", 120, -6, 6 );
  sixy1Histo->setTitle( "six y at 1;six y_{up} at 1 [mm];six-plane tracks" );

  clustersize2 = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/clustersize2", 100, 0, 100 );
  clustersize2->setTitle( "clustersize at 2; clustersize;six-plane tracks" );

  sixx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx2", 240, -12, 12 );
  sixx2Histo->setTitle( "six x at 2;six x_{out} at 2 [mm];six-plane tracks" );

  sixy2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy2", 120, -6, 6 );
  sixy2Histo->setTitle( "six y at 2;six y_{up} at 2 [mm];six-plane tracks" );

  clustersize3 = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/clustersize3", 100, 0, 100 );
  clustersize3->setTitle( "clustersize at 3; clustersize;six-plane tracks" );

  sixx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx3", 240, -12, 12 );
  sixx3Histo->setTitle( "six x at 3;six x_{out} at 3 [mm];six-plane tracks" );

  sixy3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy3", 120, -6, 6 );
  sixy3Histo->setTitle( "six y at 3;six y_{up} at 3 [mm];six-plane tracks" );

  clustersize4 = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/clustersize4", 100, 0, 100 );
  clustersize4->setTitle( "clustersize at 4; clustersize;six-plane tracks" );

  sixx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx4", 240, -12, 12 );
  sixx4Histo->setTitle( "six x at 4;six x_{out} at 4 [mm];six-plane tracks" );

  sixy4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy4", 120, -6, 6 );
  sixy4Histo->setTitle( "six y at 4;six y_{up} at 4 [mm];six-plane tracks" );

  clustersize5 = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/clustersize5", 100, 0, 100 );
  clustersize5->setTitle( "clustersize at 5; clustersize;six-plane tracks" );

  sixx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx5", 240, -12, 12 );
  sixx5Histo->setTitle( "six x at 5;six x_{out} at 5 [mm];six-plane tracks" );

  sixy5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy5", 120, -6, 6 );
  sixy5Histo->setTitle( "six y at 5;six y_{up} at 5 [mm];six-plane tracks" );

  // GBL:
  AIDAProcessor::tree(this)->mkdir("GBL");

  derxtiltHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/derxtilt", 100, -0.1, 0.1 );
  derxtiltHisto->setTitle( "ddx/dtilt;ddx/dtilt [mm/rad];align hits" );

  derytiltHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/derytilt", 100, -10, 10 );
  derytiltHisto->setTitle( "ddy/dtilt;ddy/dtilt [mm/rad];align hits" );

  derxturnHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/derxturn", 100, -10, 10 );
  derxturnHisto->setTitle( "ddx/dturn;ddx/dturn [mm/rad];align hits" );

  deryturnHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/deryturn", 100, -1, 1 );
  deryturnHisto->setTitle( "ddy/dturn;ddy/dturn [mm/rad];align hits" );

  selxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/selx", 240, -12, 12 );
  selxHisto->setTitle( "x at DUT, sel GBL;six x_{out} at DUT [mm];selected tracks" );

  selyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/sely", 120, -6, 6 );
  selyHisto->setTitle( "y at DUT, sel GBL;six y_{up} at DUT [mm];selected tracks" );

  selaxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/selax", 100, -5, 5 );
  selaxHisto->setTitle( "track angle x, sel GBL;x angle [mrad];tracks" );

  selayHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/selay", 100, -5, 5 );
  selayHisto->setTitle( "track angle y, sel GBL;y angle [mrad];tracks" );

  seldxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldx", 100, -150, 150 );
  seldxHisto->setTitle( "track match x, sel GBL;#Deltax [#mum];tracks" );

  seldyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldy", 100, -150, 150 );
  seldyHisto->setTitle( "track match y, sel GBL;#Deltay [#mum];tracks" );

  selkxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/selkx", 100, -10, 10 );
  selkxHisto->setTitle( "kink x, sel GBL;kink x [mrad];tracks" );

  selkyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/selky", 100, -10, 10 );
  selkyHisto->setTitle( "kink y, sel GBL;kink y [mrad];tracks" );

  seldx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldx1", 100, -100, 100 );
  seldx1Histo->setTitle( "triplet resid x at 1, sel GBL;#Deltax [#mum];tracks" );

  seldy1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldy1", 100, -100, 100 );
  seldy1Histo->setTitle( "triplet resid y at 1, sel GBL;#Deltay [#mum];tracks" );

  seldx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldx3", 100, -1000, 1000 );
  seldx3Histo->setTitle( "triplet resid x at 3, sel GBL;#Deltax [#mum];tracks" );

  seldy3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldy3", 100, -1000, 1000 );
  seldy3Histo->setTitle( "triplet resid y at 3, sel GBL;#Deltay [#mum];tracks" );

  seldx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldx4", 100, -1500, 1500 );
  seldx4Histo->setTitle( "triplet resid x at 4, sel GBL;#Deltax [#mum];tracks" );

  seldy4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldy4", 100, -1500, 1500 );
  seldy4Histo->setTitle( "triplet resid y at 4, sel GBL;#Deltay [#mum];tracks" );

  seldx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldx5", 100, -3000, 3000 );
  seldx5Histo->setTitle( "triplet resid x at 5, sel GBL;#Deltax [#mum];tracks" );

  seldy5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldy5", 100, -3000, 3000 );
  seldy5Histo->setTitle( "triplet resid y at 5, sel GBL;#Deltay [#mum];tracks" );

  /*seldx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldx6", 100, -500, 500 );
    seldx6Histo->setTitle( "triplet resid x at DUT, sel GBL;#Deltax [#mum];tracks" );

    seldy6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldy6", 100, -500, 500 );
    seldy6Histo->setTitle( "triplet resid y at DUT, sel GBL;#Deltay [#mum];tracks" );*/

  gblndfHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblndf", 16, -0.5, 15.5 );
  gblndfHisto->setTitle( "GBL fit NDF;GBL NDF;tracks" );

  gblchi2aHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblchi2a", 100, 0, 100 );
  gblchi2aHisto->setTitle( "GBL fit chi2, DoF 8;GBL chi2;tracks" );

  gblchi2bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblchi2b", 100, 0, 100 );
  gblchi2bHisto->setTitle( "GBL fit chi2;GBL chi2;tracks" );

  gblprbHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblprb", 1000, 0, 1 );
  gblprbHisto->setTitle( "GBL fit probability;GBL fit probability;tracks" );

  gblprbxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "GBL/gblprbx", 240, -12, 12, 1000, 0, 1 );
  gblprbxHisto->setTitle( "GBL fit probability vs. x at DUT;x at DUT [mm]; GBL fit probability;tracks" );

  gblprbyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "GBL/gblprby", 120, -6, 6, 1000, 0, 1 );
  gblprbyHisto->setTitle( "GBL fit probability vs. y at DUT;y at DUT [mm]; GBL fit probability;tracks" );

  gblnxy_plane0 = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "GBL/gblnxy_plane0", 35, 0., 18.4, 35, 0., 18.4,0,5 );
  gblnxy_plane0->setTitle( "GBL intra-pixel weighted occurrence of all CSs;GBL track x at plane0 [#mum];GBL track y at plane0 [mm];events" );

  gblnxy1_plane0 = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "GBL/gblnxy1_plane0", 35, 0., 18.4, 35, 0., 18.4);
  gblnxy1_plane0->setTitle( "GBL intra-pixel occurrence of all CSs;GBL track x at plane0 [#mum];GBL track y at plane0 [mm];events" );

  gblnxy_plane3 = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "GBL/gblnxy_plane3", 35, 0., 18.4, 35, 0., 18.4,0,5 );
  gblnxy_plane3->setTitle( "GBL intra-pixel weighted occurrence of all CSs;GBL track x at plane3 [#mum];GBL track y at plane3 [mm];events" );

  gblnxy1_plane3 = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "GBL/gblnxy1_plane3", 35, 0., 18.4, 35, 0., 18.4);
  gblnxy1_plane3->setTitle( "GBL intra-pixel occurrence of all CSs;GBL track x at plane3 [#mum];GBL track y at plane3 [mm];events" );

  gblnCS1xy_plane0 = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "GBL/gblnCS1xy_plane0", 35, 0., 18.4, 35, 0., 18.4 );
  gblnCS1xy_plane0->setTitle( "GBL in-pixel occurrence of CS1;GBL track x at plane0 [#mum];GBL track y at plane0 [mm];events" );

  gblnCS2xy_plane0 = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "GBL/gblnCS2xy_plane0", 35, 0., 18.4, 35, 0., 18.4 );
  gblnCS2xy_plane0->setTitle( "GBL in-pixel occurrence of CS2;GBL track x at plane0 [#mum];GBL track y at plane0 [mm];events" );

  gblnCS3xy_plane0 = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "GBL/gblnCS3xy_plane0", 35, 0., 18.4, 35, 0., 18.4 );
  gblnCS3xy_plane0->setTitle( "GBL in-pixel occurrence of CS3;GBL track x at plane0 [#mum];GBL track y at plane0 [mm];events" );

  gblnCS4xy_plane0 = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "GBL/gblnCS4xy_plane0", 35, 0., 18.4, 35, 0., 18.4 );
  gblnCS4xy_plane0->setTitle( "GBL in-pixel occurrence of CS4;GBL track x at plane0 [#mum];GBL track y at plane0 [mm];events" );

  gblnCS5xy_plane0 = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "GBL/gblnCS5xy_plane0", 35, 0., 18.4, 35, 0., 18.4 );
  gblnCS5xy_plane0->setTitle( "GBL in-pixel occurrence of CS5;GBL track x at plane0 [#mum];GBL track y at plane0 [mm];events" );

  gblnCS6xy_plane0 = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "GBL/gblnCS6xy_plane0", 35, 0., 18.4, 35, 0., 18.4 );
  gblnCS6xy_plane0->setTitle( "GBL in-pixel occurrence of CS6;GBL track x at plane0 [#mum];GBL track y at plane0 [mm];events" );

  gblnCS7xy_plane0 = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "GBL/gblnCS7xy_plane0", 35, 0., 18.4, 35, 0., 18.4 );
  gblnCS7xy_plane0->setTitle( "GBL in-pixel occurrence of CS>6;GBL track x at plane0 [#mum];GBL track y at plane0 [mm];events" );

  gblnCS1xy_plane3 = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "GBL/gblnCS1xy_plane3", 35, 0., 18.4, 35, 0., 18.4 );
  gblnCS1xy_plane3->setTitle( "GBL in-pixel occurrence of CS1;GBL track x at plane3 [#mum];GBL track y at plane3 [mm];events" );

  gblnCS2xy_plane3 = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "GBL/gblnCS2xy_plane3", 35, 0., 18.4, 35, 0., 18.4 );
  gblnCS2xy_plane3->setTitle( "GBL in-pixel occurrence of CS2;GBL track x at plane3 [#mum];GBL track y at plane3 [mm];events" );

  gblnCS3xy_plane3 = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "GBL/gblnCS3xy_plane3", 35, 0., 18.4, 35, 0., 18.4 );
  gblnCS3xy_plane3->setTitle( "GBL in-pixel occurrence of CS3;GBL track x at plane3 [#mum];GBL track y at plane3 [mm];events" );

  gblnCS4xy_plane3 = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "GBL/gblnCS4xy_plane3", 35, 0., 18.4, 35, 0., 18.4 );
  gblnCS4xy_plane3->setTitle( "GBL in-pixel occurrence of CS4;GBL track x at plane3 [#mum];GBL track y at plane3 [mm];events" );

  gblnCS5xy_plane3 = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "GBL/gblnCS5xy_plane3", 35, 0., 18.4, 35, 0., 18.4 );
  gblnCS5xy_plane3->setTitle( "GBL in-pixel occurrence of CS5;GBL track x at plane3 [#mum];GBL track y at plane3 [mm];events" );

  gblnCS6xy_plane3 = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "GBL/gblnCS6xy_plane3", 35, 0., 18.4, 35, 0., 18.4 );
  gblnCS6xy_plane3->setTitle( "GBL in-pixel occurrence of CS6;GBL track x at plane3 [#mum];GBL track y at plane3 [mm];events" );

  gblnCS7xy_plane3 = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "GBL/gblnCS7xy_plane3", 35, 0., 18.4, 35, 0., 18.4 );
  gblnCS7xy_plane3->setTitle( "GBL in-pixel occurrence of CS>6;GBL track x at plane3 [#mum];GBL track y at plane3 [mm];events" );

 gblnCS1xy1_plane3 = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "GBL/gblnCS1xy1_plane3", 110, -11, 11, 60, -6, 6 );
  gblnCS1xy1_plane3->setTitle( "GBL occurrence of CS1;GBL track x at plane3 [mm];GBL track y at plane3 [mm];events" );

  gblnCS2xy1_plane3 = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "GBL/gblnCS2xy1_plane3", 110, -11, 11, 60, -6, 6 );
  gblnCS2xy1_plane3->setTitle( "GBL occurrence of CS2;GBL track x at plane3 [mm];GBL track y at plane3 [mm];events" );

  gblnCS3xy1_plane3 = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "GBL/gblnCS3xy1_plane3", 110, -11, 11, 60, -6, 6 );
  gblnCS3xy1_plane3->setTitle( "GBL occurrence of CS3;GBL track x at plane3 [mm];GBL track y at plane3 [mm];events" );

  gblnCS4xy1_plane3 = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "GBL/gblnCS4xy1_plane3", 110, -11, 11, 60, -6, 6 );
  gblnCS4xy1_plane3->setTitle( "GBL occurrence of CS4;GBL track x at plane3 [mm];GBL track y at plane3 [mm];events" );

  gblnCS5xy1_plane3 = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "GBL/gblnCS5xy1_plane3", 110, -11, 11, 60, -6, 6 );
  gblnCS5xy1_plane3->setTitle( "GBL occurrence of CS5;GBL track x at plane3 [mm];GBL track y at plane3 [mm];events" );

  gblnCS6xy1_plane3 = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "GBL/gblnCS6xy1_plane3", 110, -11, 11, 60, -6, 6 );
  gblnCS6xy1_plane3->setTitle( "GBL occurrence of CS>5;GBL track x at plane3 [mm];GBL track y at plane3 [mm];events" );


  // bad fits:

  badxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/badx", 240, -12, 12 );
  badxHisto->setTitle( "x at DUT, bad GBL;six x_{out} at DUT [mm];bad tracks" );

  badyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/bady", 120, -6, 6 );
  badyHisto->setTitle( "y at DUT, bad GBL;six y_{up} at DUT [mm];bad tracks" );

  badaxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/badax", 100, -5, 5 );
  badaxHisto->setTitle( "track angle x, bad GBL;x angle [mrad];tracks" );

  badayHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baday", 100, -5, 5 );
  badayHisto->setTitle( "track angle y, bad GBL;y angle [mrad];tracks" );

  baddxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddx", 100, -150, 150 );
  baddxHisto->setTitle( "track match x, bad GBL;#Deltax [#mum];tracks" );

  baddyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddy", 100, -150, 150 );
  baddyHisto->setTitle( "track match y, bad GBL;#Deltay [#mum];tracks" );

  badkxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/badkx", 100, -10, 10 );
  badkxHisto->setTitle( "kink x, bad GBL;kink x [mrad];tracks" );

  badkyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/badky", 100, -10, 10 );
  badkyHisto->setTitle( "kink y, bad GBL;kink y [mrad];tracks" );

  baddx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddx1", 100, -100, 100 );
  baddx1Histo->setTitle( "triplet resid x at 1, bad GBL;#Deltax [#mum];tracks" );

  baddy1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddy1", 100, -100, 100 );
  baddy1Histo->setTitle( "triplet resid y at 1, bad GBL;#Deltay [#mum];tracks" );

  baddx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddx3", 100, -1000, 1000 );
  baddx3Histo->setTitle( "triplet resid x at 3, bad GBL;#Deltax [#mum];tracks" );

  baddy3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddy3", 100, -1000, 1000 );
  baddy3Histo->setTitle( "triplet resid y at 3, bad GBL;#Deltay [#mum];tracks" );

  baddx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddx4", 100, -1500, 1500 );
  baddx4Histo->setTitle( "triplet resid x at 4, bad GBL;#Deltax [#mum];tracks" );

  baddy4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddy4", 100, -1500, 1500 );
  baddy4Histo->setTitle( "triplet resid y at 4, bad GBL;#Deltay [#mum];tracks" );

  baddx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddx5", 100, -3000, 3000 );
  baddx5Histo->setTitle( "triplet resid x at 5, bad GBL;#Deltax [#mum];tracks" );

  baddy5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddy5", 100, -3000, 3000 );
  baddy5Histo->setTitle( "triplet resid y at 5, bad GBL;#Deltay [#mum];tracks" );

  /*baddx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddx6", 100, -250, 250 );
    baddx6Histo->setTitle( "triplet resid x at DUT, bad GBL;#Deltax [#mum];tracks" );

    baddy6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddy6", 100, -250, 250 );
    baddy6Histo->setTitle( "triplet resid y at DUT, bad GBL;#Deltay [#mum];tracks" );*/

  // good fits:

  goodxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/goodx", 240, -12, 12 );
  goodxHisto->setTitle( "x at DUT, good GBL;six x_{out} at DUT [mm];good tracks" );

  goodyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/goody", 120, -6, 6 );
  goodyHisto->setTitle( "y at DUT, good GBL;six y_{up} at DUT [mm];good tracks" );

  goodx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/goodx1", 100, -100, 100 );
  goodx1Histo->setTitle( "triplet resid x at 1, good GBL;#Deltax [#mum];tracks" );

  goody1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/goody1", 100, -100, 100 );
  goody1Histo->setTitle( "triplet resid y at 1, good GBL;#Deltay [#mum];tracks" );

  /*goodx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/goodx6", 100, -250, 250 );
    goodx6Histo->setTitle( "triplet resid x at 6, good GBL;#Deltax [#mum];tracks" );

    goody6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/goody6", 100, -250, 250 );
    goody6Histo->setTitle( "triplet resid y at 6, good GBL;#Deltay [#mum];tracks" );*/

  // look at fit:

  gblax0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax0", 100, -1, 1 );
  gblax0Histo->setTitle( "GBL x angle at plane 0;x angle at plane 0 [mrad];tracks" );

  gbldx0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx0", 100, -10, 10 );
  gbldx0Histo->setTitle( "GBL x shift at plane 0;x shift at plane 0 [#mum];tracks" );

  gbldx01Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx01", 100, -10, 10 );
  gbldx01Histo->setTitle( "GBL x shift at plane 0;x shift at plane 0 [mm];tracks" );

  gblrx0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx0", 250, -25, 25 );
  gblrx0Histo->setTitle( "GBL x resid at plane 0;x resid at plane 0 [#mum];tracks" );

  gblry0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblry0", 250, -25, 25 );
  gblry0Histo->setTitle( "GBL y resid at plane 0;y resid at plane 0 [#mum];tracks" );

  gblpx0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx0", 100, -10, 10 );
  gblpx0Histo->setTitle( "GBL x pull at plane 0;x pull at plane 0 [#sigma];tracks" );

  gblpx0_unbHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx0_unb", 100, -10, 10 );
  gblpx0_unbHisto->setTitle( "GBL x unbiased pull at plane 0;x unbiased pull at plane 0 [#sigma];tracks" );

  gblpy0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpy0", 100, -10, 10 );
  gblpy0Histo->setTitle( "GBL y pull at plane 0;y pull at plane 0 [#sigma];tracks" );

  gblpy0_unbHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpy0_unb", 100, -10, 10 );
  gblpy0_unbHisto->setTitle( "GBL y unbiased pull at plane 0;y unbiased pull at plane 0 [#sigma];tracks" );

  gblqx0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx0", 200, -1, 1 );
  gblqx0Histo->setTitle( "GBL x kink at plane 0;x kink at plane 0 [mrad];tracks" );

  gblrxvsx0 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblrxvsx0", 110, -11, 11, -100, 100 );
  gblrxvsx0->setTitle( "gbl x resid vs x;x [mm];<#sigma_{x}> [#mum]" );

  gblryvsy0 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblryvsy0", 110, -11, 11, -100, 100 );
  gblryvsy0->setTitle( "gbl y resid vs y;y [mm];<#sigma_{y}> [#mum]" );

  gblrxvsx01 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblrxvsx01", 110, -11, 11, -100, 100 );
  gblrxvsx01->setTitle( "gbl x resid vs x;x [mm];<#sigma_{x}> [#mum]" );

  gblryvsy01 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblryvsy01", 110, -11, 11, -100, 100 );
  gblryvsy01->setTitle( "gbl y resid vs y;y [mm];<#sigma_{y}> [#mum]" );

  gblrxvsxpix0 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblrxvsxpix0", 200, 0., 18.4, -100, 100 );
  gblrxvsxpix0->setTitle( "gbl x resid vs x;x [mm];<#sigma_{x}> [#mum]" );

  gblryvsypix0 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblryvsypix0", 200, 0., 18.4, -100, 100 );
  gblryvsypix0->setTitle( "gbl y resid vs y;y [mm];<#sigma_{y}> [#mum]" );

  gblrxvsxpix01 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblrxvsxpix01", 200, 0., 18.4, -100, 100 );
  gblrxvsxpix01->setTitle( "gbl x resid vs x at plane 0;x [#mum];<#sigma_{x}> [#mum]" );

  gblryvsypix01 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblryvsypix01", 200, 0., 18.4, -100, 100 );
  gblryvsypix01->setTitle( "gbl y resid vs y at plane 0;y [#mum];<#sigma_{y}> [#mum]" );

  gblrxvsxpix01_CS1 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblrxvsxpix01_CS1", 200, 0., 18.4, -100, 100 );
  gblrxvsxpix01_CS1->setTitle( "gbl x resid vs x at plane 0 for CS1;x [#mum];<#sigma_{x}> [#mum]" );

  gblryvsypix01_CS1 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblryvsypix01_CS1", 200, 0., 18.4, -100, 100 );
  gblryvsypix01_CS1->setTitle( "gbl y resid vs y at plane 0 for CS1;y [#mum];<#sigma_{y}> [#mum]" );

  gblrxvsxpix01_CS2 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblrxvsxpix01_CS2", 200, 0., 18.4, -100, 100 );
  gblrxvsxpix01_CS2->setTitle( "gbl x resid vs x at plane 0 for CS2;x [#mum];<#sigma_{x}> [#mum]" );

  gblryvsypix01_CS2 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblryvsypix01_CS2", 200, 0., 18.4, -100, 100 );
  gblryvsypix01_CS2->setTitle( "gbl y resid vs y at plane 0 for CS2;y [#mum];<#sigma_{y}> [#mum]" );

  gblrxvsxpix01_CS3 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblrxvsxpix01_CS3", 200, 0., 18.4, -100, 100 );
  gblrxvsxpix01_CS3->setTitle( "gbl x resid vs x at plane 0 for CS3;x [#mum];<#sigma_{x}> [#mum]" );

  gblryvsypix01_CS3 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblryvsypix01_CS3", 200, 0.,18.4 , -100, 100 );
  gblryvsypix01_CS3->setTitle( "gbl y resid vs y at plane 0 for CS3;y [#mum];<#sigma_{y}> [#mum]" );

  gblrxvsxpix01_CS4 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblrxvsxpix01_CS4", 200, 0.,18.4 , -100, 100 );
  gblrxvsxpix01_CS4->setTitle( "gbl x resid vs x at plane 0 for CS4;x [#mum];<#sigma_{x}> [#mum]" );

  gblryvsypix01_CS4 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblryvsypix01_CS4", 200, 0.,18.4 , -100, 100 );
  gblryvsypix01_CS4->setTitle( "gbl y resid vs y at plane 0 for CS4;y [#mum];<#sigma_{y}> [#mum]" );

  gblrxvsxpix51_CS1 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblrxvsxpix51_CS1", 200, 0.,18.4 , -100, 100 );
  gblrxvsxpix51_CS1->setTitle( "gbl x resid vs x at plane 5 for CS1;x [#mum];<#sigma_{x}> [#mum]" );

  gblryvsypix51_CS1 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblryvsypix51_CS1", 200, 0.,18.4 , -100, 100 );
  gblryvsypix51_CS1->setTitle( "gbl y resid vs y at plane 5 for CS1;y [#mum];<#sigma_{y}> [#mum]" );

  gblrxvsxpix51_CS2 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblrxvsxpix51_CS2", 200, 0.,18.4 , -100, 100 );
  gblrxvsxpix51_CS2->setTitle( "gbl x resid vs x at plane 5 for CS2;x [#mum];<#sigma_{x}> [#mum]" );

  gblryvsypix51_CS2 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblryvsypix51_CS2", 200, 0.,18.4 , -100, 100 );
  gblryvsypix51_CS2->setTitle( "gbl y resid vs y at plane 5 for CS2;y [#mum];<#sigma_{y}> [#mum]" );

  gblrxvsxpix51_CS3 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblrxvsxpix51_CS3", 200, 0.,18.4 , -100, 100 );
  gblrxvsxpix51_CS3->setTitle( "gbl x resid vs x at plane 5 for CS3;x [#mum];<#sigma_{x}> [#mum]" );

  gblryvsypix51_CS3 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblryvsypix51_CS3", 200, 0.,18.4 , -100, 100 );
  gblryvsypix51_CS3->setTitle( "gbl y resid vs y at plane 5 for CS3;y [#mum];<#sigma_{y}> [#mum]" );

  gblrxvsxpix51_CS4 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblrxvsxpix51_CS4", 200, 0.,18.4 , -100, 100 );
  gblrxvsxpix51_CS4->setTitle( "gbl x resid vs x at plane 5 for CS4;x [#mum];<#sigma_{x}> [#mum]" );

  gblryvsypix51_CS4 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblryvsypix51_CS4", 200, 0.,18.4 , -100, 100 );
  gblryvsypix51_CS4->setTitle( "gbl y resid vs y at plane 5 for CS4;y [#mum];<#sigma_{y}> [#mum]" );


  gblrxvsxpix31 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblrxvsxpix31", 200, 0.,18.4 , -100, 100 );
  gblrxvsxpix31->setTitle( "gbl x resid vs x at plane 3;x [#mum];<#sigma_{x}> [#mum]" );

  gblryvsypix31 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblryvsypix31", 200, 0.,18.4 , -100, 100 );
  gblryvsypix31->setTitle( "gbl y resid vs y at plane 3;y [#mum];<#sigma_{y}> [#mum]" );

  gblrxvsxpix3cs1 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblrxvsxpix3cs1", 200, 0.,18.4 , -100, 100 );
  gblrxvsxpix3cs1->setTitle( "gbl x resid vs x at plane 3 for CS 1;x [#mum];<#sigma_{x}> [#mum]" );

  gblryvsypix3cs1 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblryvsypix3cs1", 200, 0.,18.4 , -100, 100 );
  gblryvsypix3cs1->setTitle( "gbl y resid vs y at plane 3 for CS 1;y [#mum];<#sigma_{y}> [#mum]" );

  gblrxvsxpix3cs2 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblrxvsxpix3cs2", 200, 0.,18.4 , -100, 100 );
  gblrxvsxpix3cs2->setTitle( "gbl x resid vs x at plane 3 for CS 2;x [#mum];<#sigma_{x}> [#mum]" );

  gblryvsypix3cs2 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblryvsypix3cs2", 200, 0.,18.4 , -100, 100 );
  gblryvsypix3cs2->setTitle( "gbl y resid vs y at plane 3 for CS 2;y [#mum];<#sigma_{y}> [#mum]" );

  gblrxvsxpix3cs3 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblrxvsxpix3cs3", 200, 0.,18.4 , -100, 100 );
  gblrxvsxpix3cs3->setTitle( "gbl x resid vs x at plane 3 for CS 3;x [#mum];<#sigma_{x}> [#mum]" );

  gblryvsypix3cs3 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblryvsypix3cs3", 200, 0.,18.4 , -100, 100 );
  gblryvsypix3cs3->setTitle( "gbl y resid vs y at plane 3 for CS 3;y [#mum];<#sigma_{y}> [#mum]" );

  gblrxvsxpix3cs4 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblrxvsxpix3cs4", 200, 0.,18.4 , -100, 100 );
  gblrxvsxpix3cs4->setTitle( "gbl x resid vs x at plane 3 for CS 4;x [#mum];<#sigma_{x}> [#mum]" );

  gblryvsypix3cs4 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblryvsypix3cs4", 200, 0.,18.4 , -100, 100 );
  gblryvsypix3cs4->setTitle( "gbl y resid vs y at plane 3 for CS 4;y [#mum];<#sigma_{y}> [#mum]" );

  gblrxvsxpix3cs5 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblrxvsxpix3cs5", 200, 0.,18.4 , -100, 100 );
  gblrxvsxpix3cs5->setTitle( "gbl x resid vs x at plane 3 for CS 5;x [#mum];<#sigma_{x}> [#mum]" );

  gblryvsypix3cs5 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblryvsypix3cs5", 200, 0.,18.4 , -100, 100 );
  gblryvsypix3cs5->setTitle( "gbl y resid vs y at plane 3 for CS 5;y [#mum];<#sigma_{y}> [#mum]" );

  gblrxvsxpix3cs6 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblrxvsxpix3cs6", 200, 0.,18.4 , -100, 100 );
  gblrxvsxpix3cs6->setTitle( "gbl x resid vs x at plane 3 for CS 6;x [#mum];<#sigma_{x}> [#mum]" );

  gblryvsypix3cs6 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "GBL/gblryvsypix3cs6", 200, 0.,18.4 , -100, 100 );
  gblryvsypix3cs6->setTitle( "gbl y resid vs y at plane 3 for CS 6;y [#mum];<#sigma_{y}> [#mum]" );

  gblax1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax1", 100, -1, 1 );
  gblax1Histo->setTitle( "GBL x angle at plane 1;x angle at plane 1 [mrad];tracks" );

  gbldx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx1", 10, -10, 10 );
  gbldx1Histo->setTitle( "GBL x shift at plane 1;x shift at plane 1 [#mum];tracks" );

  gbldx11Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx11", 100, -10, 10 );
  gbldx11Histo->setTitle( "GBL x shift at plane 1;x shift at plane 1 [mm];tracks" );

  gblrx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx1", 250, -25, 25 );
  gblrx1Histo->setTitle( "GBL x resid at plane 1;x resid at plane 1 [#mum];tracks" );

  gblry1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblry1", 250, -25, 25 );
  gblry1Histo->setTitle( "GBL y resid at plane 1;y resid at plane 1 [#mum];tracks" );

  gblpx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx1", 100, -10, 10 );
  gblpx1Histo->setTitle( "GBL x pull at plane 1;x pull at plane 1 [#sigma];tracks" );

  gblpx1_unbHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx1_unb", 100, -10, 10 );
  gblpx1_unbHisto->setTitle( "GBL x unbiased pull at plane 1;x unbiased pull at plane 1 [#sigma];tracks" );

  gblpy1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpy1", 100, -10, 10 );
  gblpy1Histo->setTitle( "GBL y pull at plane 1;y pull at plane 1 [#sigma];tracks" );

  gblpy1_unbHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpy1_unb", 100, -10, 10 );
  gblpy1_unbHisto->setTitle( "GBL y unbiased pull at plane 1;y unbiased pull at plane 1 [#sigma];tracks" );

  gblqx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx1", 200, -1, 1 );
  gblqx1Histo->setTitle( "GBL x kink at plane 1;x kink at plane 1 [mrad];tracks" );

  gblsx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblsx1", 200, -10, 10 );
  gblsx1Histo->setTitle( "GBL x kink at plane 1/error;x kink at plane 1/error;tracks" );

  gbltx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbltx1", 200, -10, 10 );
  gbltx1Histo->setTitle( "GBL x kink pull at plane 1;x kink pull at plane 1;tracks" );


  gblax2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax2", 100, -1, 1 );
  gblax2Histo->setTitle( "GBL x angle at plane 2;x angle at plane 2 [mrad];tracks" );

  gbldx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx2", 100, -10, 10 );
  gbldx2Histo->setTitle( "GBL x shift at plane 2;x shift at plane 2 [#mum];tracks" );

  gbldx21Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx21", 100, -10, 10 );
  gbldx21Histo->setTitle( "GBL x shift at plane 2;x shift at plane 2 [mm];tracks" );

  gblrx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx2", 250, -25, 25 );
  gblrx2Histo->setTitle( "GBL x resid at plane 2;x resid at plane 2 [#mum];tracks" );

  gblry2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblry2", 250, -25, 25 );
  gblry2Histo->setTitle( "GBL y resid at plane 2;y resid at plane 2 [#mum];tracks" );

  gblpx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx2", 100, -10, 10 );
  gblpx2Histo->setTitle( "GBL x pull at plane 2;x pull at plane 2 [#sigma];tracks" );

  gblpx2_unbHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx2_unb", 100, -10, 10 );
  gblpx2_unbHisto->setTitle( "GBL x unbiased pull at plane 2;x unbiased pull at plane 2 [#sigma];tracks" );

  gblpy2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpy2", 100, -10, 10 );
  gblpy2Histo->setTitle( "GBL y pull at plane 2;y pull at plane 2 [#sigma];tracks" );

  gblpy2_unbHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpy2_unb", 100, -10, 10 );
  gblpy2_unbHisto->setTitle( "GBL y unbiased pull at plane 2;y unbiased pull at plane 2 [#sigma];tracks" );

  gblqx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx2", 200, -1, 1 );
  gblqx2Histo->setTitle( "GBL x kink at plane 2;x kink at plane 2 [mrad];tracks" );

  gblsx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblsx2", 200, -10, 10 );
  gblsx2Histo->setTitle( "GBL x kink at plane 2/error;x kink at plane 2/error;tracks" );

  gbltx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbltx2", 200, -10, 10 );
  gbltx2Histo->setTitle( "GBL x kink pull at plane 2;x kink pull at plane 2;tracks" );


  gblax3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax3", 100, -1, 1 );
  gblax3Histo->setTitle( "GBL x angle at plane 3;x angle at plane 3 [mrad];tracks" );

  gbldx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx3", 100, -10, 10 );
  gbldx3Histo->setTitle( "GBL x shift at plane 3;x shift at plane 3 [#mum];tracks" );

  gbldx31Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx31", 100, -10, 10 );
  gbldx31Histo->setTitle( "GBL x shift at plane 3;x shift at plane 3 [mm];tracks" );

  gblrx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx3", 250, -25, 25 );
  gblrx3Histo->setTitle( "GBL x resid at plane 3;x resid at plane 3 [#mum];tracks" );

  gblry3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblry3", 250, -25, 25 );
  gblry3Histo->setTitle( "GBL y resid at plane 3;y resid at plane 3 [#mum];tracks" );

  gblrx3_cs1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx3_cs1", 250, -25, 25 );
  gblrx3_cs1Histo->setTitle( "GBL x resid at plane 3_cs1;x resid at plane 3 [#mum];tracks" );

  gblry3_cs1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblry3_cs1", 250, -25, 25 );
  gblry3_cs1Histo->setTitle( "GBL y resid at plane 3_cs1;y resid at plane 3 [#mum];tracks" );

  gblrx3_cs2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx3_cs2", 250, -25, 25 );
  gblrx3_cs2Histo->setTitle( "GBL x resid at plane 3_cs2;x resid at plane 3 [#mum];tracks" );

  gblry3_cs2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblry3_cs2", 250, -25, 25 );
  gblry3_cs2Histo->setTitle( "GBL y resid at plane 3_cs2;y resid at plane 3 [#mum];tracks" );

  gblrx3_cs3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx3_cs3", 250, -25, 25 );
  gblrx3_cs3Histo->setTitle( "GBL x resid at plane 3_cs3;x resid at plane 3 [#mum];tracks" );

  gblry3_cs3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblry3_cs3", 250, -25, 25 );
  gblry3_cs3Histo->setTitle( "GBL x resid at plane 3_cs3;x resid at plane 3 [#mum];tracks" );

  gblrx3_cs4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx3_cs4", 250, -25, 25 );
  gblrx3_cs4Histo->setTitle( "GBL y resid at plane 3_cs4;y resid at plane 3 [#mum];tracks" );

  gblry3_cs4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblry3_cs4", 250, -25, 25 );
  gblry3_cs4Histo->setTitle( "GBL y resid at plane 3_cs4;y resid at plane 3 [#mum];tracks" );

  gblrx3_cs5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx3_cs5", 250, -25, 25 );
  gblrx3_cs5Histo->setTitle( "GBL y resid at plane 3_cs5y resid at plane 3 [#mum];tracks" );

  gblry3_cs5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblry3_cs5", 250, -25, 25 );
  gblry3_cs5Histo->setTitle( "GBL y resid at plane 3_cs5;y resid at plane 3 [#mum];tracks" );

  gblrx3_cs6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx3_cs6", 250, -25, 25 );
  gblrx3_cs6Histo->setTitle( "GBL y resid at plane 3_cs>5;y resid at plane 3 [#mum];tracks" );

  gblry3_cs6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblry3_cs6", 250, -25, 25 );
  gblry3_cs6Histo->setTitle( "GBL y resid at plane 3_cs>5;y resid at plane 3 [#mum];tracks" );

  gblpx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx3", 100, -10, 10 );
  gblpx3Histo->setTitle( "GBL x pull at plane 3;x pull at plane 3 [#sigma];tracks" );

  gblpy3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpy3", 100, -10, 10 );
  gblpy3Histo->setTitle( "GBL y pull at plane 3;y pull at plane 3 [#sigma];tracks" );

  gblpx3_cs1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx3_cs1", 100, -10, 10 );
  gblpx3_cs1Histo->setTitle( "GBL x pull at plane 3_cs1;x pull at plane 3 [#sigma];tracks" );

  gblpy3_cs1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpy3_cs1", 100, -10, 10 );
  gblpy3_cs1Histo->setTitle( "GBL y pull at plane 3_cs1;y pull at plane 3 [#sigma];tracks" );

  gblpx3_cs2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx3_cs2", 100, -10, 10 );
  gblpx3_cs2Histo->setTitle( "GBL x pull at plane 3_cs2;x pull at plane 3 [#sigma];tracks" );

  gblpy3_cs2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpy3_cs2", 100, -10, 10 );
  gblpy3_cs2Histo->setTitle( "GBL y pull at plane 3_cs2;y pull at plane 3 [#sigma];tracks" );

  gblpx3_cs3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx3_cs3", 100, -10, 10 );
  gblpx3_cs3Histo->setTitle( "GBL x pull at plane 3_cs3;x pull at plane 3 [#sigma];tracks" );

  gblpy3_cs3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpy3_cs3", 100, -10, 10 );
  gblpy3_cs3Histo->setTitle( "GBL y pull at plane 3_cs3;y pull at plane 3 [#sigma];tracks" );

  gblpx3_cs4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx3_cs4", 100, -10, 10 );
  gblpx3_cs4Histo->setTitle( "GBL x pull at plane 3_cs4;x pull at plane 3 [#sigma];tracks" );

  gblpy3_cs4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpy3_cs4", 100, -10, 10 );
  gblpy3_cs4Histo->setTitle( "GBL y pull at plane 3_cs4;y pull at plane 3 [#sigma];tracks" );

  gblpx3_cs5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx3_cs5", 100, -10, 10 );
  gblpx3_cs5Histo->setTitle( "GBL x pull at plane 3_cs5;x pull at plane 3 [#sigma];tracks" );

  gblpy3_cs5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpy3_cs5", 100, -10, 10 );
  gblpy3_cs5Histo->setTitle( "GBL y pull at plane 3_cs5;y pull at plane 3 [#sigma];tracks" );

  gblpx3_cs6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx3_cs6", 100, -10, 10 );
  gblpx3_cs6Histo->setTitle( "GBL x pull at plane 3_cs>5;x pull at plane 3 [#sigma];tracks" );

  gblpy3_cs6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpy3_cs6", 100, -10, 10 );
  gblpy3_cs6Histo->setTitle( "GBL y pull at plane 3_cs>5;y pull at plane 3 [#sigma];tracks" );

  gblpx3_cs7Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx3_cs7", 100, -10, 10 );
  gblpx3_cs7Histo->setTitle( "GBL x pull at plane 3_cs>6;x pull at plane 3 [#sigma];tracks" );

  gblpy3_cs7Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpy3_cs7", 100, -10, 10 );
  gblpy3_cs7Histo->setTitle( "GBL y pull at plane 3_cs>6;y pull at plane 3 [#sigma];tracks" );

  gblpx3_unbHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx3_unb", 100, -10, 10 );
  gblpx3_unbHisto->setTitle( "GBL x unbiased pull at plane 3;x unbiased pull at plane 3 [#sigma];tracks" );

  gblpy3_unbHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpy3_unb", 100, -10, 10 );
  gblpy3_unbHisto->setTitle( "GBL y unbiased pull at plane 3;y unbiased pull at plane 3 [#sigma];tracks" );

  gblqx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx3", 200, -1, 1 );
  gblqx3Histo->setTitle( "GBL x kink at plane 3;x kink at plane 3 [mrad];tracks" );

  gblsx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblsx3", 200, -10, 10 );
  gblsx3Histo->setTitle( "GBL x kink at plane 3/error;x kink at plane 3/error;tracks" );

  gbltx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbltx3", 200, -10, 10 );
  gbltx3Histo->setTitle( "GBL x kink pull at plane 3;x kink pull at plane 3;tracks" );


  gblax4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax4", 100, -1, 1 );
  gblax4Histo->setTitle( "GBL x angle at plane 4;x angle at plane 4 [mrad];tracks" );

  gbldx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx4", 100, -10, 10 );
  gbldx4Histo->setTitle( "GBL x shift at plane 4;x shift at plane 4 [#mum];tracks" );

  gbldx41Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx41", 100, -10, 10 );
  gbldx41Histo->setTitle( "GBL x shift at plane 4;x shift at plane 4 [nm];tracks" );

  gblrx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx4", 250, -25, 25 );
  gblrx4Histo->setTitle( "GBL x resid at plane 4;x resid at plane 4 [#mum];tracks" );

  gblry4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblry4", 250, -25, 25 );
  gblry4Histo->setTitle( "GBL y resid at plane 4;y resid at plane 4 [#mum];tracks" );

  gblpx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx4", 100, -10, 10 );
  gblpx4Histo->setTitle( "GBL x pull at plane 4;x pull at plane 4 [#sigma];tracks" );

  gblpx4_unbHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx4_unb", 100, -10, 10 );
  gblpx4_unbHisto->setTitle( "GBL x unbiased pull at plane 4;x unbiased pull at plane 4 [#sigma];tracks" );

  gblpy4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpy4", 100, -10, 10 );
  gblpy4Histo->setTitle( "GBL y pull at plane 4;y pull at plane 4 [#sigma];tracks" );

  gblpy4_unbHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpy4_unb", 100, -10, 10 );
  gblpy4_unbHisto->setTitle( "GBL y unbiased pull at plane 4;y unbiased pull at plane 4 [#sigma];tracks" );

  gblqx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx4", 200, -1, 1 );
  gblqx4Histo->setTitle( "GBL x kink at plane 4;x kink at plane 4 [mrad];tracks" );

  gblsx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblsx4", 200, -10, 10 );
  gblsx4Histo->setTitle( "GBL x kink at plane 4/error;x kink at plane 4/error;tracks" );

  gbltx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbltx4", 200, -10, 10 );
  gbltx4Histo->setTitle( "GBL x kink pull at plane 4;x kink pull at plane 4;tracks" );


  gblax5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax5", 100, -1, 1 );
  gblax5Histo->setTitle( "GBL x angle at plane 5;x angle at plane 5 [mrad];tracks" );

  gbldx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx5", 100, -10, 10 );
  gbldx5Histo->setTitle( "GBL x shift at plane 5;x shift at plane 5 [#mum];tracks" );

  gbldx51Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx51", 100, -10, 10 );
  gbldx51Histo->setTitle( "GBL x shift at plane 5;x shift at plane 5 [mm];tracks" );

  gblrx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx5", 250, -25, 25 );
  gblrx5Histo->setTitle( "GBL x resid at plane 5;x resid at plane 5 [#mum];tracks" );

  gblry5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblry5", 250, -25, 25 );
  gblry5Histo->setTitle( "GBL y resid at plane 5;y resid at plane 5 [#mum];tracks" );

  gblpx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx5", 100, -10, 10 );
  gblpx5Histo->setTitle( "GBL x pull at plane 5;x pull at plane 5 [#sigma];tracks" );

  gblpx5_unbHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx5_unb", 100, -10, 10 );
  gblpx5_unbHisto->setTitle( "GBL x unbiased pull at plane 5;x unbiased pull at plane 5 [#sigma];tracks" );

  gblpy5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpy5", 100, -10, 10 );
  gblpy5Histo->setTitle( "GBL y pull at plane 5;y pull at plane 5 [#sigma];tracks" );

  gblpy5_unbHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpy5_unb", 100, -10, 10 );
  gblpy5_unbHisto->setTitle( "GBL y unbiased pull at plane 5;y unbiased pull at plane 5 [#sigma];tracks" );

  gblqx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx5", 200, -1, 1 );
  gblqx5Histo->setTitle( "GBL x kink at plane 5;x kink at plane 5 [mrad];tracks" );


  /*gblax6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax6", 100, -1, 1 );
    gblax6Histo->setTitle( "GBL x angle at DUT;x angle at DUT [mrad];tracks" );

    gbldx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx6", 100, -1000, 1000 );
    gbldx6Histo->setTitle( "GBL x shift at DUT;x shift at DUT [#mum];tracks" );

    gbldy6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldy6", 100, -1000, 1000 );
    gbldy6Histo->setTitle( "GBL y shift at DUT;y shift at DUT [#mum];tracks" );

    gblrx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx6", 100, -250, 250 );
    gblrx6Histo->setTitle( "GBL x resid at DUT;x resid at DUT [#mum];tracks" );

    gblry6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblry6", 100, -100, 100 );
    gblry6Histo->setTitle( "GBL y resid at DUT;y resid at DUT [#mum];tracks" );

    gblpx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx6", 100, -10, 10 );
    gblpx6Histo->setTitle( "GBL x pull at DUT;x pull at DUT [#sigma];tracks" );

    gblpy6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpy6", 100, -10, 10 );
    gblpy6Histo->setTitle( "GBL y pull at DUT;y pull at DUT [#sigma];tracks" );

    gblqx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx6", 200, -10, 10 );
    gblqx6Histo->setTitle( "GBL x kink at DUT;x kink at DUT [mrad];tracks" );

    gblsx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblsx6", 100, -10, 10 );
    gblsx6Histo->setTitle( "GBL x kink at DUT/error;x kink at DUT/error;tracks" );

    gbltx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbltx6", 100, -10, 10 );
    gbltx6Histo->setTitle( "GBL x kink pull at DUT;x kink pull at DUT;tracks" );*/


  gblkxCentreHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkxCentre", 200, -0.2, 0.2 );
  gblkxCentreHisto->setTitle( "GBL kink angle at centre;centre kink [mrad];tracks" );

  gblkxCentre1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkxCentre1", 200, -0.2, 0.2 );
  gblkxCentre1Histo->setTitle( "GBL kink angle at centre in rad;centre kink [rad];tracks" );

  gblkx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkx1", 200, -0.2, 0.2 );
  gblkx1Histo->setTitle( "GBL kink angle at plane 1;plane 1 kink [mrad];tracks" );

  gblkx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkx2", 200, -0.2, 0.2 );
  gblkx2Histo->setTitle( "GBL kink angle at plane 2;plane 2 kink [mrad];tracks" );

  gblkx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkx3", 200, -0.2, 0.2 );
  gblkx3Histo->setTitle( "GBL kink angle at plane 3;plane 3 kink [mrad];tracks" );

  gblkx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkx4", 200, -0.2, 0.2 );
  gblkx4Histo->setTitle( "GBL kink angle at plane 4;plane 4 kink [mrad];tracks" );

  gblkx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkx5", 200, -0.2, 0.2 );
  gblkx5Histo->setTitle( "GBL kink angle at plane 5;plane 5 kink [mrad];tracks" );

  /*gblkx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkx6", 100, -1, 1 );
    gblkx6Histo->setTitle( "GBL kink angle at plane 6;plane 6 kink [mrad];tracks" );*/

  kinkpixvsxy = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "GBL/kinkpixvsxy", 15, 0., 18.4, 15, 0., 18.4, 0, 100);
  kinkpixvsxy->setTitle( "GBL intra-pixel kink;GBL track x at plane3 [#mum];GBL track y at plane3 [#mum];sqrt(<kink^{2}>) [mrad]" );

  // z intersect:

  sixzx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixzx3", 100, -50, 250 );
  sixzx3Histo->setTitle( "intersect z-x, kink > 3 mrad;intersect z(x) - z_{2} [mm];tracks" );

  sixzy3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixzy3", 100, -50, 250 );
  sixzy3Histo->setTitle( "intersect z-y, kink > 3 mrad;intersect z(y) - z_{2} [mm];tracks" );

  sixzx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixzx2", 100, -50, 250 );
  sixzx2Histo->setTitle( "intersect z-x, kink > 2 mrad;intersect z(x) - z_{2} [mm];tracks" );

  sixzy2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixzy2", 100, -50, 250 );
  sixzy2Histo->setTitle( "intersect z-y, kink > 2 mrad;intersect z(y) - z_{2} [mm];tracks" );

  sixzx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixzx1", 100, -50, 250 );
  sixzx1Histo->setTitle( "intersect z-x, kink > 1 mrad;intersect z(x) - z_{2} [mm];tracks" );

  sixzy1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixzy1", 100, -50, 250 );
  sixzy1Histo->setTitle( "intersect z-y, kink > 1 mrad;intersect z(y) - z_{2} [mm];tracks" );

  sixkxzyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkxzy", 200, -20, 20 );
  sixkxzyHisto->setTitle( "kink x at DUT;kink x [mrad];tracks" );

  sixkyzxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkyzx", 200, -20, 20 );
  sixkyzxHisto->setTitle( "kink y at DUT;kink y [mrad];tracks" );

  sixkxzxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkxzx", 200, -20, 20 );
  sixkxzxHisto->setTitle( "kink x at DUT;kink x [mrad];tracks" );

  sixkyzyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkyzy", 200, -20, 20 );
  sixkyzyHisto->setTitle( "kink y at DUT;kink y [mrad];tracks" );

  hIso = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/nIso", 2, -0.5, 1.5 );
  sixkyzyHisto->setTitle( "n Isolated t/driplets;isolation;tracks" );


  streamlog_out( DEBUG2 ) << "Histogram booking completed \n\n" << std::endl;

#else

  streamlog_out( MESSAGE4 ) << "No histogram produced because Marlin doesn't use AIDA" << std::endl;

#endif

  return;
}






#endif
