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


EUTelGBLFitter::EUTelGBLFitter() : Processor("EUTelGBLFitter"), _inputCollectionTelescope(""), _isFirstEvent(0), _eBeam(0), _nEvt(0), _nPlanes(0), _track_match_cut(0.15),  _planePosition() {
  // modify processor description
  _description = "Analysis for DATURA reference analysis ";

  // processor parameters
  registerInputCollection( LCIO::TRACKERHIT, "InputCollection", "Name of the input TrackerHit collection of the telescope", _inputCollectionTelescope, std::string("") );

  registerProcessorParameter( "Ebeam", "Beam energy [GeV]", _eBeam, 0.0);

  registerOptionalParameter( "triResCut", "Upstream/Downstream triplet residual cut [mm]", _triplet_res_cut, 0.1);

  registerProcessorParameter( "matchingCut", "cut for matching in x coordinate in mm", _track_match_cut, 0.15);

  registerProcessorParameter( "slopeCut", "cut for track slopes in x coordinate in rad", _slope_cut, 0.002);

  registerProcessorParameter( "excludedPlanes", "Planes to be excluded from the track fit", _excluded_planes, EVENT::IntVec() );

  //in mm. 0.1 is for 6 GeV, 20 mm. This is scaled with E and dz
  registerProcessorParameter( "eff_radius", "radius on DUT plane to accept match with triplet", _eff_radius, 0.1); 

  //1.0 means HL as is, 1.2 means 20% additional scattering
  registerProcessorParameter( "kappa", "global factor to Highland formula", _kappa, 1.0);

  registerProcessorParameter( "probchi2Cut", "Cut on Prob(chi2,ndf) rejecting bad tracks with prob < cut", _probchi2_cut, .01); 

  registerProcessorParameter( "TelescopeResolution", "Resolution parameter for each cluster size (CS) for all telescope planes. First value is average of all CSes, subsequently for CS=1 them CS=2 and so on. The last value is for all CSes larger than the previous ones. I.e. if you provide five values: <avg> <CS 1> <CS 2> <CS 3> <CS greater 3>", _telResolution, FloatVec(8, 3.5*1e-3));

  registerOptionalParameter( "DUTXResolutions", "Same as TelescopeResolution, but now only for y-direction. Also, there needs to be an additional leading NEGATIVE number for the sensorID. E.g. -20 0.5 0.7 0.4 0.3 would correspond to <sensorID (20)> <avg> <CS 1> <CS 2> <CS greater 2>, this could be followed by a further section which again starts with a negative number for the next sensorID", _dutResolutionX, FloatVec(8, 3.5*1e-3));

  registerOptionalParameter( "DUTYResolutions", "Same as DUTXResolutions but in y-direction.", _dutResolutionY, FloatVec(8., 3.5*1e-3));
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
    auto const & pos = geo::gGeometry().getPlaneZPosition(sensorID);
    _planePosition.emplace_back( pos );
    auto const & z = geo::gGeometry().getPlaneZSize(sensorID);
    auto const & rad = geo::gGeometry().getPlaneRadiationLength(sensorID);

    if(sensorID < 6) {
        _planeRadLength.emplace_back(z/rad + 0.050 / 286.6); // Plane from GEAR + Kapton
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
    streamlog_out( MESSAGE6 ) << "  Avg. reso Plane " << sensorID << " = " << _telResolution[0] << std::endl;
    _planeResolutionX[sensorID] = _telResolution;
    _planeResolutionY[sensorID] = _telResolution;
  }

  std::vector<float>* currentResVector = nullptr;
  for(auto it = _dutResolutionX.begin(); it !=_dutResolutionX.end(); it++) {
    if(*it < 0) {
      auto sensorID = static_cast<int>(-(*it));
      if(std::find(_sensorIDVec.begin(), _sensorIDVec.end(), sensorID) == _sensorIDVec.end()) {
        streamlog_out(ERROR5) << "Something went horribly wrong - probably the sensor ID is not present in GEAR file" << std::endl;
      } else {
        currentResVector = &_planeResolutionX[sensorID];
        currentResVector->clear();
      }
    } else {
      if(currentResVector) {
        currentResVector->emplace_back(*it);
      } else {
        streamlog_out(ERROR5) << "Something went horribly wrong - probably no sensor ID was given" << std::endl;
      }
    }
  }
  currentResVector = nullptr;
  for(auto it = _dutResolutionY.begin(); it !=_dutResolutionY.end(); it++) {
    if(*it < 0) {
      auto sensorID = static_cast<int>(-(*it));
      if(std::find(_sensorIDVec.begin(), _sensorIDVec.end(), sensorID) == _sensorIDVec.end()) {
        streamlog_out(ERROR5) << "Something went horribly wrong - probably the sensor ID is not present in GEAR file" << std::endl;
      } else {
        currentResVector = &_planeResolutionY[sensorID];
        currentResVector->clear();
      }
    } else {
      if(currentResVector) {
        currentResVector->emplace_back(*it);
      } else {
        streamlog_out(ERROR5) << "Something went horribly wrong - probably no sensor ID was given" << std::endl;
      }
    }
  }

  //Resolution Info
  for(auto sensorID: _sensorIDVec) {
    auto& resX = _planeResolutionX[sensorID];
    auto& resY = _planeResolutionY[sensorID];
    std::stringstream ss;
    int const width = 8;
    ss  << "---- Plane: " << sensorID << " ----" 
        << std::endl << std::left << setw(20) << setfill(' ') << "Resolution X [um]: ";
    for(size_t i = 0; i < resX.size(); ++i) {
      if(i == 0) ss << std::left << setw(width) << setfill(' ') << "<avg>";
      else if( i == resX.size()-1) ss << setw(width) << setfill(' ') << "<CS >"+std::to_string(i-1)+">\n";
      else ss << setw(width) << setfill(' ') << "<CS "+std::to_string(i)+">";
    }
    ss << std::left << setw(20) << setfill(' ') << ' ';
    std::for_each(resX.begin(), resX.end(), [&](double const & n){ ss << std::left << setw(width) << setfill(' ') << n*1000;});
    ss  << std::endl <<  std::left << setw(20) << setfill(' ') <<  "Resolution Y [um]: ";
    for(size_t i = 0; i < resY.size(); ++i) {
      if(i == 0) ss << std::left << setw(width) << setfill(' ') << "<avg>";
      else if( i == resY.size()-1) ss << setw(width) << setfill(' ') << "<CS >"+std::to_string(i-1)+">\n";
      else ss << setw(width) << setfill(' ') << "<CS "+std::to_string(i)+">";
    }
    ss << std::left << setw(20) << setfill(' ') << ' ';
    std::for_each(resY.begin(), resY.end(), [&](double const & n){ ss << std::left << setw(width) << setfill(' ') << n*1000;});
    streamlog_out( MESSAGE2 ) <<  ss.str() << std::endl;
  }

  //Material Info
  for(size_t ix = 0; ix < _sensorIDVec.size(); ++ix) {
    auto sensorID = _sensorIDVec[ix];
    auto const & pos = _planePosition[ix];
    auto const & rad = _planeRadLength[ix];
    streamlog_out( MESSAGE2 ) << "Plane " << sensorID << " is located at (z-position) " << pos << " mm.\nIt has a radiation length of " << rad << " [mm/mm]\n"; 
  }

  for(auto sensorID: _sensorIDVec) {
    if( std::find( std::begin(_excluded_planes), std::end(_excluded_planes), sensorID ) == _excluded_planes.end() ) {
      _excludedSensorMap.insert( std::make_pair(sensorID, false) ); 
    } else {
      _excludedSensorMap.insert( std::make_pair(sensorID, true) ); 
    }  
  }

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

  nAllTelHitHisto->fill( hits.size() );
  nAllDUTHitHisto->fill( DUThits.size() );

  streamlog_out(DEBUG4) << "Event " << event->getEventNumber() << " contains " << hits.size() << " telscope and " << DUThits.size() << " DUT hits" << std::endl;

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
  gblutil.FindTriplets(hits, std::array<size_t,3>{3, 4, 5}, _triplet_res_cut, _slope_cut, downstream_triplets);
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
  gblutil.FindTriplets(hits, std::array<size_t,3>{0, 1, 2}, _triplet_res_cut, _slope_cut, upstream_triplets);

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
  //double eff_radius = _eff_radius * 6. / _eBeam * (_planePosition[1] - _planePosition[0]) / 20.; 
/*  
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

*/

  //----------------------------------------------------------------------------
  // six: triplets A and driplets B
  // matching and GBL fit
  // kinks: triplets A vs driplets B

  //Match half way between the upstream and downstream arms
  double zMid = 0.5*(_planePosition[_nPlanes-3] + _planePosition[2]);
  std::vector<EUTelTripletGBLUtility::track> telescope_tracks;
  gblutil.MatchTriplets(upstream_triplets,downstream_triplets, zMid, _track_match_cut, telescope_tracks);

  streamlog_out(DEBUG4) << "Found " << telescope_tracks.size() << " tracks from matching t/driplets." << endl;
  
  for( auto& tr: telescope_tracks ){
    //unused: auto const & trip = tr.get_upstream();
    auto const & drip = tr.get_downstream();
    EUTelTripletGBLUtility::triplet srip(tr.gethit(0), tr.gethit(2), tr.gethit(5)); // seed triplet is called srip

    std::vector<double> xAplanes(_nPlanes);
    std::vector<double> yAplanes(_nPlanes);

    for (size_t i = 0; i < _nPlanes; i++){
      xAplanes[i] = srip.getx_at(_planePosition[i]);
      yAplanes[i] = srip.gety_at(_planePosition[i]);
    }

    // Track kinks as difference in triplet slopes:
    //      double kx = drip.slope().x - trip.slope().x; //kink
    //      double ky = drip.slope().y - trip.slope().y;
    double kx = tr.kink_x();
    double ky = tr.kink_y();

    // Track impact position at DUT from Downstream:
    double xB = drip.getx_at(zMid);
    double yB = drip.gety_at(zMid);

    // Track impact position at DUT from Upstream:
    double xA = srip.getx_at(zMid);
    double yA = srip.gety_at(zMid);

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
    
    std::vector<unsigned int> ilab;
    std::map<size_t, size_t> ilabToSensorID;

    size_t DUTCount = _nPlanes-6;
    std::vector<double> rx (_nPlanes, -1.0);
    std::vector<double> ry (_nPlanes, -1.0);
    std::vector<double> trackhitx (_nPlanes, -1.0);
    std::vector<double> trackhity (_nPlanes, -1.0);
    std::vector<double> trackhitxloc (_nPlanes, -1.0);
    std::vector<double> trackhityloc (_nPlanes, -1.0);
    std::vector<bool> hasHit (_nPlanes, false);

    double step = 0.;
    s = 0.;
    for( size_t ipl = 0; ipl < _nPlanes; ++ipl ){

      //We have to add all the planes, the up and downstream arm of the telescope will definitely have
      //hits, the DUTs might not though! The first and last three planes are the telescope.
      EUTelTripletGBLUtility::hit const * trackhit = nullptr;
      if(ipl < 3) {
        trackhit = &tr.gethit(ipl);
      } else if( ipl < 3+DUTCount) {
//        auto sensorID = _sensorIDVec[ipl];
//        if(triplet.has_DUT(sensorID)) trackhit = &triplet.get_DUT_Hit(sensorID);
//        else if(driplet.has_DUT(sensorID)) trackhit = &driplet.get_DUT_Hit(sensorID);
      } else {
        trackhit = &tr.gethit(ipl-DUTCount);
      }

      //if there is no hit we take the plane position from the geo description
      //double zz = trackhit ? trackhit->z : _planePosition[ipl];// [mm]

      auto point = gbl::GblPoint( gblutil.JacobianPointToPoint( step ) );

      if(trackhit){
        //fill the trackhit relevant histograms
        fillTrackhitHisto(*trackhit, ipl);

        double dz = trackhit->z - srip.base().z;
        double xs = srip.base().x + srip.slope().x * dz; // Ax at plane
        double ys = srip.base().y + srip.slope().y * dz; // Ay at plane

        trackhitx[ipl] = trackhit->x;
        trackhity[ipl] = trackhit->y;
        trackhitxloc[ipl] = trackhit->locx;
        trackhityloc[ipl] = trackhit->locy;

        rx[ipl] = trackhit->x - xs;
        ry[ipl] = trackhit->y - ys;

        auto meas = Eigen::Vector2d(rx[ipl], ry[ipl]);

        auto const & resVecX = _planeResolutionX[trackhit->plane];
        auto const & resVecY = _planeResolutionY[trackhit->plane];
        auto cluX = static_cast<size_t>(trackhit->clustersizex);
        //auto cluX = static_cast<size_t>(trackhit->clustersize);
        auto cluY = static_cast<size_t>(trackhit->clustersizey);
        //auto cluY = static_cast<size_t>(trackhit->clustersize);
        double _x_resolution_tmp = (cluX >= resVecX.size()) ? resVecX.back() : resVecX[cluX];
        double _y_resolution_tmp = (cluY >= resVecY.size()) ? resVecY.back() : resVecY[cluY];
      
   	    auto measPrec = Eigen::Vector2d(1.0/_x_resolution_tmp/_x_resolution_tmp, 1.0/_y_resolution_tmp/_y_resolution_tmp);

        //The measurement is only included if it is a non excluded plane
        auto currentSensorID = _sensorIDVec[ipl];
        if( !_excludedSensorMap[currentSensorID]  ) {
          point.addMeasurement( proL2m, meas, measPrec );
        }

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
      }  
      point.addScatterer( scat, _planeWscatSi[ipl] );

      // streamlog_out(DEBUG4) << "Added Scatterer:\n" << _planeWscatSi[ipl] << std::endl; 
      // streamlog_out(DEBUG4) << "Meas Precision:\n" << measPrec << std::endl; 
      traj_points.emplace_back(point);
      s += step;
      sPoint.push_back( s );
      ilab.push_back( sPoint.size() );
      ilabToSensorID.insert( std::make_pair(ilab.back(), _sensorIDVec[ipl]) );

      if( ipl < 5) {
        double distplane = _planePosition[ipl+1] - _planePosition[ipl];
        step = 0.21*distplane;

        auto point = gbl::GblPoint( gblutil.JacobianPointToPoint( step ) );
        point.addScatterer( scat, _planeWscatAir[ipl] );
        s += step;
        traj_points.push_back(point);
        sPoint.push_back( s );

        // streamlog_out(DEBUG4) << "Added Air Scatterer:\n" << _planeWscatAir[ipl] << "\nat: " << s << std::endl; 
        step = 0.58*distplane;
        auto point1 = gbl::GblPoint( gblutil.JacobianPointToPoint( step ) );
        point1.addScatterer( scat, _planeWscatAir[ipl] );
        s += step;
        traj_points.push_back(point1);
        sPoint.push_back( s );
        // streamlog_out(DEBUG4) << "Added Air Scatterer:\n" << _planeWscatAir[ipl] << "\nat :" << s << std::endl; 
	      
        step = 0.21*distplane; // remaing distance to next plane
      }
    } // loop over planes

    double Chi2;
    int Ndf;
    double lostWeight;

    gbl::GblTrajectory traj(traj_points, false ); // curvature = false
    std::string fit_optionList = "";
    traj.fit( Chi2, Ndf, lostWeight, fit_optionList );

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
    if( probchi > _probchi2_cut){
      _ngbl++;

      Eigen::VectorXd aCorrection(2);
      Eigen::MatrixXd aCovariance(5,5);

      auto ax = std::vector<double>(_nPlanes, -2);
      // double ay[8];

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
      size_t ipos = ilab[0];
      size_t currentSensor = ilabToSensorID[ipos];

      for(size_t ipl = 0; ipl < ilab.size(); ++ipl) {
        ipos = ilab[ipl];
        currentSensor = ilabToSensorID[ipos];

        traj.getResults( ipos, aCorrection, aCovariance );
        traj.getMeasResults( ipos, ndata, aResiduals, aMeasErrors, aResErrors, aDownWeights );
        traj.getScatResults( ipos, ndata, aKinks, aKinkErrors, kResErrors, kDownWeights );

        aResiduals[0] = rx[ipl] - aCorrection[3];
        aResiduals[1] = ry[ipl] - aCorrection[4];

        gblaxHistos[ipl]->fill( aCorrection[1]*1E3 ); // angle x [mrad]
        gbldxHistos[ipl]->fill( aCorrection[3]*1E3 ); // shift x [um]

//      gbldx01Histo->fill( aCorrection[3] ); // shift x [mm]

        gblrxHistos[ipl]->fill( ( rx[ipl] - aCorrection[3] ) * 1E3 ); // residual x [um]
        gblryHistos[ipl]->fill( ( ry[ipl] - aCorrection[4] ) * 1E3 ); // residual y [um]
        gblpxHistos[ipl]->fill( aResiduals[0] / aResErrors[0] ); // pull
        gblpyHistos[ipl]->fill( aResiduals[1] / aResErrors[1] ); // pull
 //     if(_dut_plane == 0) gblpx0_unbHisto->fill( (rx[0] - aCorrection[3]) / sqrt(_telResolution[0]*_telResolution[0] + aCovariance(3,3)) ); // unbiased pull
 //     if(_dut_plane == 0) gblpy0_unbHisto->fill( (ry[0] - aCorrection[4]) / sqrt(_telResolution[0]*_telResolution[0] + aCovariance(4,4)) ); // unbiased pull

        gblqxHistos[ipl]->fill( aKinks[0]*1E3 ); // kink RESIDUAL (measured kink - fit kink)
      
        ax[ipl] = aCorrection[1]; // angle correction at plane, for kinks

        // TProfile for res_x a.f.o. x
        gblrxvsx[ipl]->fill( xAplanes.at(ipl), sqrt(TMath::Pi()/2.)*fabs(aResiduals[0]));
        gblryvsy[ipl]->fill( yAplanes.at(ipl), sqrt(TMath::Pi()/2.)*fabs(aResiduals[1]));
        gblrxvsx1[ipl]->fill((trackhitx[ipl] - aResiduals[0]), sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); // seed corrected
        gblryvsy1[ipl]->fill((trackhity[ipl] - aResiduals[1]), sqrt(TMath::Pi()/2.)*fabs(aResiduals[1]));
        
        std::array<double,2> corrPos;
        corrPos[0] = trackhitx[ipl] - aResiduals[0]; //+ .4e-3;  should be zero ! :/
        corrPos[1] = trackhity[ipl] - aResiduals[1]; // + .15e-3; should be zero (no alignment on plane 0, not even pre-align)      

        std::array<double, 3> globalCorrPos {corrPos[0], corrPos[1], _planePosition[ipl]};
        std::array<double, 3> localCorrPos;
 
        geo::gGeometry().master2Local(currentSensor, globalCorrPos, localCorrPos);
        corrPos[0] = localCorrPos[0];
        corrPos[1] = localCorrPos[1];

        int nx = (corrPos[0]) / pixel_size;
        int ny = (corrPos[1]) / pixel_size;
        int invsignx = -(corrPos[0]) / fabs((corrPos[0]));
        int invsigny = -(corrPos[1]) / fabs((corrPos[1])); 

        // do again for extrapolated srip position, userd ONLY for gblrxvsxpix0 and gblryvsypix0
        int nx1 = (xAplanes[ipl]) / pixel_size;
        int ny1 = (yAplanes[ipl]) / pixel_size;
        int invsignx1 = -(xAplanes[ipl]) / fabs((xAplanes[ipl]));
        int invsigny1 = -(yAplanes[ipl]) / fabs((yAplanes[ipl]));

        gblrxvsxpix[ipl]->fill(                    (xAplanes[ipl] +invsignx1*(abs(nx1) +1.)*pixel_size)*1e3, sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); 
        gblryvsypix[ipl]->fill(                    (yAplanes[ipl] +invsigny1*(abs(ny1) +1.)*pixel_size)*1e3, sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])); 
        gblrxvsxpix1[ipl]->fill( ((trackhitx[ipl] - aResiduals[0])+invsignx*(abs(nx) +1.)*pixel_size)*1e3, sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])); 
        gblryvsypix1[ipl]->fill( ((trackhity[ipl] - aResiduals[1])+invsigny*(abs(ny) +1.)*pixel_size)*1e3, sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])); 

        // clustersize-specific plots
        auto CStot = tr.gethit(currentSensor).clustersize;
        auto CSx = tr.gethit(currentSensor).clustersizex;
        auto CSy = tr.gethit(currentSensor).clustersizey;

        std::array<double,2> modPos;
        modPos[0] = ((corrPos[0])+(invsignx*(abs(nx) +.5) + 0.5) *pixel_size)*1e3;
        modPos[1] = ((corrPos[1])+(invsigny*(abs(ny) +.5) + 0.5) *pixel_size)*1e3;
  
        // overlay of all CSs
        gblnxy[ipl]->fill(modPos[0], modPos[1], CStot);
        gblnxy1[ipl]->fill(modPos[0], modPos[1], 1);
        gblcluxvscluy[ipl]->fill(CSx, CSy, 1);

        size_t CSTotIndex = (CStot > 6) ? 6 : CStot-1;
        size_t CSXIndex = (CSx > 6) ? 6 : CSx-1;
        size_t CSYIndex = (CSy > 6) ? 6 : CSy-1;

        gblnCSxy_tot[ipl][CSTotIndex]->fill(modPos[0], modPos[1] );    
        gblnCSxy_x[ipl][CSXIndex]->fill(modPos[0], modPos[1] );    
        gblnCSxy_y[ipl][CSYIndex]->fill(modPos[0], modPos[1] );    

        if(CSTotIndex < 4){
          gblrxvsxpix1CS[ipl][CSTotIndex]->fill( modPos[0], sqrt(TMath::Pi()/2.)*fabs(aResiduals[0])*1e3);
          gblryvsypix1CS[ipl][CSTotIndex]->fill( modPos[1], sqrt(TMath::Pi()/2.)*fabs(aResiduals[1])*1e3);
        }  
      } 
/*
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
*/
      for(size_t ix = 1; ix < _nPlanes; ++ix) {
        gblkxHistos[ix-1]->fill( (ax[ix] - ax[ix-1])*1E3 );
      }
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
	if( abs( zx - zMid ) < 30 ) {
	  sixkyzxHisto->fill( ky*1E3 );
	  sixkxzxHisto->fill( kx*1E3 ); // plot with gap, fittp0g.C("sixkxzx")
	}
      }

      if( abs( ky ) > 0.001 ) {
	sixzy1Histo->fill( zy - _planePosition[2] );
	if( abs( zy - zMid ) < 30 ) {
	  sixkxzyHisto->fill( kx*1E3 );
	  sixkyzyHisto->fill( ky*1E3 ); // plot with gap
	}
      }

    }//match

    //------------------------------------------------------------------------

  } // Loop over found tracks

  nsixHisto->fill( telescope_tracks.size() );
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
}

void EUTelGBLFitter::fillTrackhitHisto(EUTelTripletGBLUtility::hit const & hit, int ipl){
  clustersizeTotal[ipl]->fill(hit.clustersize);
  clustersizeX[ipl]->fill(hit.clustersizex);
  clustersizeY[ipl]->fill(hit.clustersizey);

  sixXHistos[ipl]->fill( -hit.x );
  sixYHistos[ipl]->fill( -hit.y );
}

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
void EUTelGBLFitter::bookHistos(){
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  // telescope and DUT hits per plane:
  AIDAProcessor::tree(this)->mkdir("Telescope");

  nAllTelHitHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/nalltelhit", 201, -0.5, 200.5 );
  nAllTelHitHisto->setTitle( "Telescope hits/event;telescope hits;events" );

  nAllDUTHitHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/nallduthit", 51, -0.5, 50.5 );
  nAllDUTHitHisto->setTitle( "DUT hits/event;DUT hits;events" );

  // telescope correlation plots:
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
  AIDAProcessor::tree(this)->mkdir("Tracks");

  nsixHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/nsix", 21, -0.5, 20.5 );
  nsixHisto->setTitle( "telescope six-plane-tracks;six-plane-tracks;events" );

  AIDAProcessor::tree(this)->mkdir("Tracks/ClusterSize");
  for(size_t ix = 0; ix < _sensorIDVec.size(); ++ix) {
    auto sensorIdString = std::to_string(_sensorIDVec[ix]);
    std::string histNameTotal = "Tracks/ClusterSize/ClusterSizeTotal_"+sensorIdString;
    std::string histNameX = "Tracks/ClusterSize/ClusterSizeX_"+sensorIdString;
    std::string histNameY = "Tracks/ClusterSize/ClusterSizeY_"+sensorIdString;

    clustersizeTotal.push_back(AIDAProcessor::histogramFactory(this)->
     createHistogram1D( histNameTotal, 101, -0.5, 100 ));
    clustersizeTotal.back()->setTitle( "total cluster size on plane "+sensorIdString+";#hit pixels in cluster;count" );

    clustersizeX.push_back(AIDAProcessor::histogramFactory(this)->
     createHistogram1D( histNameX, 101, -0.5, 100 ));
    clustersizeX.back()->setTitle( "total cluster size in x on plane "+sensorIdString+";#hit pixels in cluster in x;count" );

    clustersizeY.push_back(AIDAProcessor::histogramFactory(this)->
     createHistogram1D( histNameY, 101, -0.5, 100 ));
    clustersizeY.back()->setTitle( "total cluster size in y on plane "+sensorIdString+";#hit pixels in cluster in y;count" );
  }

  AIDAProcessor::tree(this)->mkdir("Tracks/HitPosition");
  for(size_t ix = 0; ix < _sensorIDVec.size(); ++ix) {
    auto sensorIdString = std::to_string(_sensorIDVec[ix]);
    std::string histNameX = "Tracks/HitPosition/X_"+sensorIdString;
    std::string histNameY = "Tracks/HitPosition/Y_"+sensorIdString;
  
    sixXHistos.push_back(AIDAProcessor::histogramFactory(this)->
     createHistogram1D( histNameX, 240, -12, 12 ));
    sixXHistos.back()->setTitle( "measured hit on track x on plane "+sensorIdString+";x position [mm]; tracks" );

    sixYHistos.push_back(AIDAProcessor::histogramFactory(this)->
     createHistogram1D( histNameY, 240, -12, 12 ));
    sixYHistos.back()->setTitle( "measured hit on track y on plane "+sensorIdString+";y position [mm]; tracks" );
  }

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

  // look at fit:
  AIDAProcessor::tree(this)->mkdir("GBL/Angles");
  AIDAProcessor::tree(this)->mkdir("GBL/Shifts");
  AIDAProcessor::tree(this)->mkdir("GBL/Residuals");
  AIDAProcessor::tree(this)->mkdir("GBL/Pulls");
  AIDAProcessor::tree(this)->mkdir("GBL/Kinks");
  AIDAProcessor::tree(this)->mkdir("GBL/ResVsPos");
  AIDAProcessor::tree(this)->mkdir("GBL/ClusterSize");
  AIDAProcessor::tree(this)->mkdir("GBL/ResVsPosVsCS");

  for(size_t ix = 0; ix < _sensorIDVec.size(); ++ix) {
    auto sensorIdString = std::to_string(_sensorIDVec[ix]);
    std::string histNameAngleX = "GBL/Angles/gblax_"+sensorIdString;
    std::string histNameShiftX = "GBL/Shifts/gbldx_"+sensorIdString;
    std::string histNameResX = "GBL/Residuals/gblrx_"+sensorIdString;
    std::string histNameResY = "GBL/Residuals/gblry_"+sensorIdString;
    std::string histNamePullX = "GBL/Pulls/gblpx_"+sensorIdString;
    std::string histNamePullY = "GBL/Pulls/gblpy_"+sensorIdString;
    std::string histNameKinkX = "GBL/Kinks/gblqx_"+sensorIdString;
    std::string histNameXResVsX = "GBL/ResVsPos/gblrxvsx_"+sensorIdString;  
    std::string histNameYResVsY = "GBL/ResVsPos/gblryvsy_"+sensorIdString;  
    std::string histNameXResVsX1 = "GBL/ResVsPos/gblrxvsx1_"+sensorIdString;  
    std::string histNameYResVsY1 = "GBL/ResVsPos/gblryvsy1_"+sensorIdString;  
    std::string histNameXResVsXPix = "GBL/ResVsPos/gblrxvsxpix_"+sensorIdString;
    std::string histNameYResVsYPix = "GBL/ResVsPos/gblryvsypix_"+sensorIdString;
    std::string histNameXResVsXPix1 = "GBL/ResVsPos/gblrxvsxpix1_"+sensorIdString;
    std::string histNameYResVsYPix1 = "GBL/ResVsPos/gblryvsypix1_"+sensorIdString;
    std::string histNameClusterSizes = "GBL/ClusterSize/gblnxy_"+sensorIdString;
    std::string histNameWeightClusterSizes = "GBL/ClusterSize/gblnxy1_"+sensorIdString;
    std::string histNameClusterXYCorrelation = "GBL/ClusterSize/gblcluxvscluy_"+sensorIdString;

    gblaxHistos.push_back(AIDAProcessor::histogramFactory(this)->
     createHistogram1D( histNameAngleX, 100, -1, 1 ));
    gblaxHistos.back()->setTitle( "GBL x angle at plane  "+sensorIdString+";x angle at plane [mrad];tracks" );

    gbldxHistos.push_back(AIDAProcessor::histogramFactory(this)->
     createHistogram1D( histNameShiftX, 100, -10, 10 ));
    gbldxHistos.back()->setTitle( "GBL x shift at plane "+sensorIdString+";x shift at plane [#mum];tracks" );

    gblrxHistos.push_back(AIDAProcessor::histogramFactory(this)->
     createHistogram1D( histNameResX, 250, -25, 25 ));
    gblrxHistos.back()->setTitle( "GBL x resid at plane "+sensorIdString+";x resid at plane [#mum];tracks" );

    gblryHistos.push_back(AIDAProcessor::histogramFactory(this)->
     createHistogram1D( histNameResY, 250, -25, 25 ));
    gblryHistos.back()->setTitle( "GBL y resid at plane "+sensorIdString+";y resid at plane [#mum];tracks" );

    gblpxHistos.push_back(AIDAProcessor::histogramFactory(this)->
     createHistogram1D( histNamePullX, 100, -10, 10 ));
    gblpxHistos.back()->setTitle( "GBL x pull at plane "+sensorIdString+";x pull at plane [#sigma];tracks" );

    gblpyHistos.push_back(AIDAProcessor::histogramFactory(this)->
     createHistogram1D( histNamePullY, 100, -10, 10 ));
    gblpyHistos.back()->setTitle( "GBL y pull at plane "+sensorIdString+";y pull at plane [#sigma];tracks" );

    gblqxHistos.push_back(AIDAProcessor::histogramFactory(this)->
     createHistogram1D( histNameKinkX, 200, -1, 1 ));
    gblqxHistos.back()->setTitle( "GBL x kink resid at plane "+sensorIdString+";x kink resid at plane [mrad];tracks" );

    gblrxvsx.push_back(AIDAProcessor::histogramFactory(this)->
     createProfile1D( histNameXResVsX, 110, -11, 11, -100, 100 ));
    gblrxvsx.back()->setTitle( "gbl x resid vs x plane"+sensorIdString+";x [mm];<#sigma_{x}> [#mum]" );

    gblryvsy.push_back(AIDAProcessor::histogramFactory(this)->
     createProfile1D( histNameYResVsY, 110, -11, 11, -100, 100 ));
    gblryvsy.back()->setTitle( "gbl y resid vs y plane"+sensorIdString+";y [mm];<#sigma_{y}> [#mum]" );

    gblrxvsx1.push_back(AIDAProcessor::histogramFactory(this)->
     createProfile1D( histNameXResVsX1, 110, -11, 11, -100, 100 ));
    gblrxvsx1.back()->setTitle( "gbl x resid vs x plane"+sensorIdString+";x [mm];<#sigma_{x}> [#mum]" );

    gblryvsy1.push_back(AIDAProcessor::histogramFactory(this)->
     createProfile1D( histNameYResVsY1, 110, -11, 11, -100, 100 ));
    gblryvsy1.back()->setTitle( "gbl y resid vs y plane"+sensorIdString+";y [mm];<#sigma_{y}> [#mum]" );

    gblrxvsxpix.push_back(AIDAProcessor::histogramFactory(this)->
     createProfile1D( histNameXResVsXPix, 200, 0., 18.4, -100, 100 ));
    gblrxvsxpix.back()->setTitle( "gbl x resid vs x plane"+sensorIdString+";x [mm];<#sigma_{x}> [#mum]" );

    gblryvsypix.push_back(AIDAProcessor::histogramFactory(this)->
     createProfile1D( histNameYResVsYPix, 200, 0., 18.4, -100, 100 ));
    gblryvsypix.back()->setTitle( "gbl y resid vs y plane"+sensorIdString+";y [mm];<#sigma_{y}> [#mum]" );

    gblrxvsxpix1.push_back(AIDAProcessor::histogramFactory(this)->
     createProfile1D( histNameXResVsXPix1, 200, 0., 18.4, -100, 100 ));
    gblrxvsxpix1.back()->setTitle( "gbl x resid vs x at plane"+sensorIdString+";x [mm];<#sigma_{x}> [#mum]" );

    gblryvsypix1.push_back(AIDAProcessor::histogramFactory(this)->
     createProfile1D( histNameYResVsYPix1, 200, 0., 18.4, -100, 100 ));
    gblryvsypix1.back()->setTitle( "gbl y resid vs y at plane"+sensorIdString+";y [mm];<#sigma_{y}> [#mum]" );

    gblnxy.push_back(AIDAProcessor::histogramFactory(this)->
     createProfile2D( histNameClusterSizes, 35, 0., 18.4, 35, 0., 18.4,0,5 ));
    gblnxy.back()->setTitle( "GBL intra-pixel weighted occurrence of all CSs;GBL track x at plane "+sensorIdString+" [#mum];GBL track y at plane [mm];events" );

    gblnxy1.push_back(AIDAProcessor::histogramFactory(this)->
     createHistogram2D( histNameWeightClusterSizes, 35, 0., 18.4, 35, 0., 18.4));
    gblnxy1.back()->setTitle( "GBL intra-pixel occurrence of all CSs;GBL track x at plane "+sensorIdString+" [#mum];GBL track y at plane [mm];events" );

    gblcluxvscluy.push_back(AIDAProcessor::histogramFactory(this)->
     createHistogram2D( histNameClusterXYCorrelation, 7, 0.5, 7.5, 7, 0.5, 7.5));
    gblcluxvscluy.back()->setTitle( "CS X vs. CS Y on plane "+sensorIdString+";Clustersize X;Clustersize Y;clusters" );

    std::string histNameCSRaw = "GBL/ClusterSize/CS";
    gblnCSxy_tot.emplace_back(std::vector<AIDA::IHistogram2D*>(7, nullptr));
    gblnCSxy_x.emplace_back(std::vector<AIDA::IHistogram2D*>(7, nullptr));
    gblnCSxy_y.emplace_back(std::vector<AIDA::IHistogram2D*>(7, nullptr));

    for(size_t iy = 0; iy < gblnCSxy_tot.back().size(); ++iy){
      std::string CSString = std::to_string(iy+1);
      std::string histNameTot = histNameCSRaw+"tot"+CSString+"_plane"+sensorIdString;
      std::string histNameX = histNameCSRaw+"X"+CSString+"_plane"+sensorIdString;
      std::string histNameY = histNameCSRaw+"Y"+CSString+"_plane"+sensorIdString;

      gblnCSxy_tot.back()[iy] = AIDAProcessor::histogramFactory(this)->
       createHistogram2D( histNameTot, 35, 0., 18.4, 35, 0., 18.4 );
      gblnCSxy_tot.back()[iy]->setTitle( "GBL in-pixel occurrence of CS_{tot} "+CSString+";GBL track x at plane "+sensorIdString+" [#mum];GBL track y at plane [mm];events" ); 

      gblnCSxy_x.back()[iy] = AIDAProcessor::histogramFactory(this)->
       createHistogram2D( histNameX, 35, 0., 18.4, 35, 0., 18.4 );
      gblnCSxy_x.back()[iy]->setTitle( "GBL in-pixel occurrence of CS_{x} "+CSString+";GBL track x at plane "+sensorIdString+" [#mum];GBL track y at plane [mm];events" ); 

      gblnCSxy_y.back()[iy] = AIDAProcessor::histogramFactory(this)->
       createHistogram2D( histNameY, 35, 0., 18.4, 35, 0., 18.4 );
      gblnCSxy_y.back()[iy]->setTitle( "GBL in-pixel occurrence of CS_{y} "+CSString+";GBL track x at plane "+sensorIdString+" [#mum];GBL track y at plane [mm];events" ); 
    }

    std::string histNameResVsPosByCSRawX = "GBL/ResVsPosVsCS/gblrxvsxpix_plane";
    std::string histNameResVsPosByCSRawY = "GBL/ResVsPosVsCS/gblryvsypix_plane";
    gblrxvsxpix1CS.emplace_back(std::vector<AIDA::IProfile1D*>(4, nullptr));
    gblryvsypix1CS.emplace_back(std::vector<AIDA::IProfile1D*>(4, nullptr));
    for(size_t iy = 0; iy < gblrxvsxpix1CS.back().size(); ++iy){
      std::string CSString = std::to_string(iy+1);
      std::string histNameX = histNameResVsPosByCSRawX+sensorIdString+"_CS"+CSString;
      std::string histNameY = histNameResVsPosByCSRawY+sensorIdString+"_CS"+CSString;

      gblrxvsxpix1CS.back()[iy] = AIDAProcessor::histogramFactory(this)->
       createProfile1D( histNameX, 200, 0., 18.4, -100, 100 );
      gblrxvsxpix1CS.back()[iy]->setTitle( "gbl x resid vs x at plane "+sensorIdString+" for CS"+CSString+";x [#mum];<#sigma_{x}> [#mum]" );

      gblryvsypix1CS.back()[iy] = AIDAProcessor::histogramFactory(this)->
       createProfile1D( histNameY, 200, 0., 18.4, -100, 100 );
      gblryvsypix1CS.back()[iy]->setTitle( "gbl y resid vs y at plane "+sensorIdString+" for CS"+CSString+";y [#mum];<#sigma_{y}> [#mum]" );
    }
  }
/*
  gblpx0_unbHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx0_unb", 100, -10, 10 );
  gblpx0_unbHisto->setTitle( "GBL x unbiased pull at plane 0;x unbiased pull at plane 0 [#sigma];tracks" );

  gblpy0_unbHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpy0_unb", 100, -10, 10 );
  gblpy0_unbHisto->setTitle( "GBL y unbiased pull at plane 0;y unbiased pull at plane 0 [#sigma];tracks" );
*/
  gblkxCentreHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkxCentre", 200, -0.2, 0.2 );
  gblkxCentreHisto->setTitle( "GBL kink angle at centre;centre kink [mrad];tracks" );

  gblkxCentre1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkxCentre1", 200, -0.2, 0.2 );
  gblkxCentre1Histo->setTitle( "GBL kink angle at centre in rad;centre kink [rad];tracks" );

  for(size_t ix = 1; ix < _nPlanes; ++ix){
    auto sensorIdString = std::to_string(_sensorIDVec[ix]);
    std::string histNameKinkX = "GBL/Kinks/gblkx"+sensorIdString;

    gblkxHistos.emplace_back(AIDAProcessor::histogramFactory(this)->
      createHistogram1D( histNameKinkX, 200, -0.2, 0.2 ));
    gblkxHistos.back()->setTitle( "GBL kink angle at plane "+sensorIdString+";plane "+sensorIdString+" kink [mrad];tracks" );
  }

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
