// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Authors:
// Philipp Roloff, DESY <mailto:philipp.roloff@desy.de>
// Joerg Behr, Hamburg Uni/DESY <joerg.behr@desy.de>
// Slava Libov, DESY <mailto:vladyslav.libov@desy.de>
// Igor Rubinskiy, DESY <mailto:igorrubinsky@gmail.com>
// Daniel Pitzl, DESY <mailto:daniel.pitzl@desy.de>
//
// Version: $Id: EUTelAlignGBL.cc,v 1.48 2009-08-01 10:49:46 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
// built only if GEAR and MARLINUTIL are used
#if defined(USE_GEAR) && defined(USE_MARLINUTIL)

// eutelescope includes ".h"
#include "EUTelAlignGBL.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"
#include "EUTelPStream.h" // process streams redi::ipstream
#include "EUTelGeometryTelescopeGeoDescription.h"

// GBL:
#include "include/GblTrajectory.h"
#include "include/MilleBinary.h"

#include "EUTelTripletGBLUtility.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"
#include "marlin/Exceptions.h"
#include "marlin/AIDAProcessor.h"

// marlin util includes
#include "Mille.h"

// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
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
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/LCCollectionVec.h>

// system includes <>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

AIDA::IHistogram1D * nAllHitHistoGBLAlign;

AIDA::IHistogram1D * tridxHistGBLAlign;
AIDA::IHistogram1D * tridyHistGBLAlign;
AIDA::IHistogram1D * tridxcHistGBLAlign;
AIDA::IHistogram1D * tridycHistGBLAlign;
AIDA::IHistogram1D * ntriHistGBLAlign;

AIDA::IHistogram1D * dridxHistGBLAlign;
AIDA::IHistogram1D * dridyHistGBLAlign;
AIDA::IHistogram1D * dridxcHistGBLAlign;
AIDA::IHistogram1D * dridycHistGBLAlign;
AIDA::IHistogram1D * ndriHistGBLAlign;

AIDA::IHistogram1D * sixdxHistGBLAlign;
AIDA::IHistogram1D * sixdyHistGBLAlign;
AIDA::IHistogram1D * sixkxHistGBLAlign;
AIDA::IHistogram1D * sixkyHistGBLAlign;

AIDA::IHistogram1D * selxHistGBLAlign;
AIDA::IHistogram1D * selyHistGBLAlign;
AIDA::IHistogram1D * selaxHistGBLAlign;
AIDA::IHistogram1D * selayHistGBLAlign;
AIDA::IHistogram1D * seldxHistGBLAlign;
AIDA::IHistogram1D * seldyHistGBLAlign;
AIDA::IHistogram1D * selkxHistGBLAlign;
AIDA::IHistogram1D * selkyHistGBLAlign;

AIDA::IHistogram1D * seldx1HistGBLAlign;
AIDA::IHistogram1D * seldy1HistGBLAlign;
AIDA::IHistogram1D * seldx3HistGBLAlign;
AIDA::IHistogram1D * seldy3HistGBLAlign;
AIDA::IHistogram1D * seldx4HistGBLAlign;
AIDA::IHistogram1D * seldy4HistGBLAlign;
AIDA::IHistogram1D * seldx5HistGBLAlign;
AIDA::IHistogram1D * seldy5HistGBLAlign;
AIDA::IHistogram1D * seldx6HistGBLAlign;
AIDA::IHistogram1D * seldy6HistGBLAlign;

AIDA::IHistogram1D * gblndfHistGBLAlign;
AIDA::IHistogram1D * gblchi2HistGBLAlign;
AIDA::IHistogram1D * gblprbHistGBLAlign;

AIDA::IHistogram1D * badxHistGBLAlign;
AIDA::IHistogram1D * badyHistGBLAlign;
AIDA::IHistogram1D * badaxHistGBLAlign;
AIDA::IHistogram1D * badayHistGBLAlign;
AIDA::IHistogram1D * baddxHistGBLAlign;
AIDA::IHistogram1D * baddyHistGBLAlign;
AIDA::IHistogram1D * badkxHistGBLAlign;
AIDA::IHistogram1D * badkyHistGBLAlign;

AIDA::IHistogram1D * baddx1HistGBLAlign;
AIDA::IHistogram1D * baddy1HistGBLAlign;
AIDA::IHistogram1D * baddx3HistGBLAlign;
AIDA::IHistogram1D * baddy3HistGBLAlign;
AIDA::IHistogram1D * baddx4HistGBLAlign;
AIDA::IHistogram1D * baddy4HistGBLAlign;
AIDA::IHistogram1D * baddx5HistGBLAlign;
AIDA::IHistogram1D * baddy5HistGBLAlign;
AIDA::IHistogram1D * baddx6HistGBLAlign;
AIDA::IHistogram1D * baddy6HistGBLAlign;

AIDA::IHistogram1D * goodx1HistGBLAlign;
AIDA::IHistogram1D * goody1HistGBLAlign;
AIDA::IHistogram1D * goodx6HistGBLAlign;
AIDA::IHistogram1D * goody6HistGBLAlign;

AIDA::IHistogram1D * gblax0HistGBLAlign;
AIDA::IHistogram1D * gbldx0HistGBLAlign;
AIDA::IHistogram1D * gblrx0HistGBLAlign;
AIDA::IHistogram1D * gblry0HistGBLAlign;
AIDA::IHistogram1D * gblpx0HistGBLAlign;
AIDA::IHistogram1D * gblpy0HistGBLAlign;

AIDA::IHistogram1D * gblax1HistGBLAlign;
AIDA::IHistogram1D * gbldx1HistGBLAlign;
AIDA::IHistogram1D * gblrx1HistGBLAlign;
AIDA::IHistogram1D * gblry1HistGBLAlign;

AIDA::IHistogram1D * gblax2HistGBLAlign;
AIDA::IHistogram1D * gbldx2HistGBLAlign;
AIDA::IHistogram1D * gblrx2HistGBLAlign;
AIDA::IHistogram1D * gblry2HistGBLAlign;

AIDA::IHistogram1D * gblax3HistGBLAlign;
AIDA::IHistogram1D * gbldx3HistGBLAlign;
AIDA::IHistogram1D * gblrx3HistGBLAlign;
AIDA::IHistogram1D * gblry3HistGBLAlign;

AIDA::IHistogram1D * gblax4HistGBLAlign;
AIDA::IHistogram1D * gbldx4HistGBLAlign;
AIDA::IHistogram1D * gblrx4HistGBLAlign;
AIDA::IHistogram1D * gblry4HistGBLAlign;

AIDA::IHistogram1D * gblax5HistGBLAlign;
AIDA::IHistogram1D * gbldx5HistGBLAlign;
AIDA::IHistogram1D * gblrx5HistGBLAlign;
AIDA::IHistogram1D * gblry5HistGBLAlign;

AIDA::IHistogram1D * gblax6HistGBLAlign;
AIDA::IHistogram1D * gbldx6HistGBLAlign;
AIDA::IHistogram1D * gbldy6HistGBLAlign;
AIDA::IHistogram1D * gblrx6HistGBLAlign;
AIDA::IHistogram1D * gblry6HistGBLAlign;

AIDA::IHistogram1D * gblkx1HistGBLAlign;
AIDA::IHistogram1D * gblkx2HistGBLAlign;
AIDA::IHistogram1D * gblkx3HistGBLAlign;
AIDA::IHistogram1D * gblkx4HistGBLAlign;
AIDA::IHistogram1D * gblkx5HistGBLAlign;
AIDA::IHistogram1D * gblkx6HistGBLAlign;

AIDA::IHistogram1D * nmHistGBLAlign;

#endif


//------------------------------------------------------------------------------
EUTelAlignGBL::EUTelAlignGBL(): Processor("EUTelAlignGBL") {

  // modify processor description
  _description = "EUTelAlignGBL uses the MILLE program to write data files for MILLEPEDE II.";

  // input collections
  std::vector<std::string > HitCollectionNameVecExample;
  HitCollectionNameVecExample.push_back("corrhits");

  registerInputCollections(LCIO::TRACKERHIT,"HitCollectionName","Hit collections name",_hitCollectionName,HitCollectionNameVecExample);
  registerProcessorParameter("Ebeam","Beam energy [GeV]",_eBeam, static_cast<double>(4.0));
  registerProcessorParameter("IsFirstAlignStep", "Bool: 1/0 (yes/no)",_IsFirstAlignStep, static_cast<int>(0));
  registerOptionalParameter("ExcludePlanes","Exclude planes from fit according to their sensor ids.",_excludePlanes_sensorIDs ,std::vector<int>());
  registerOptionalParameter("FixedPlanes","Fix sensor planes in the fit according to their sensor ids.",_FixedPlanes_sensorIDs ,std::vector<int>());
  registerOptionalParameter("MaxTrackCandidatesTotal","Maximal number of track candidates (Total).",_maxTrackCandidatesTotal, static_cast <int> (10000000));
  registerOptionalParameter("MaxTrackCandidates","Maximal number of track candidates.",_maxTrackCandidates, static_cast <int> (2000));
  registerOptionalParameter("BinaryFilename","Name of the Millepede binary file.",_binaryFilename, string ("mille.bin"));
  registerOptionalParameter("TelescopeResolution","Resolution of the telescope for Millepede (sigma_x=sigma_y.",_telescopeResolution, static_cast <float> (0.010));
  registerOptionalParameter("AlignMode","Number of alignment constants used. Available mode are: "
                              "\nXYZShifts - shifts in X and Y"
                              "\nXYShiftsRotZ - shifts in X and Y and rotation around the Z axis,"
                              "\nXYZShiftsRotZ - shifts in X,Y and Z and rotation around the Z axis",
                              _alignModeString, std::string("XYShiftsRotZ"));
  registerOptionalParameter("triCut", "Upstream triplet residual cut [um]", _triCut, 0.30);
  registerOptionalParameter("driCut", "Downstream triplet residual cut [um]", _driCut, 0.40);
  registerOptionalParameter("sixCut", "Upstream-Downstream Track matching cut [um]", _sixCut, 0.60);
  registerOptionalParameter("slopeCut", "t(d)riplet slope cut [radian]", _slopeCut, 0.01);
  registerOptionalParameter("GeneratePedeSteerfile","Generate a steering file for the pede program.",_generatePedeSteerfile, static_cast <int> (0));
  registerOptionalParameter("PedeSteerfileName","Name of the steering file for the pede program.",_pedeSteerfileName, string("steer_mille.txt"));
  registerProcessorParameter("kappa","Global factor to Highland formula", _kappa, static_cast <double>(1.0)); // 1.0 means HL as is, 1.2 means 20% additional scattering
}

//------------------------------------------------------------------------------
void EUTelAlignGBL::init() {
  geo::gGeometry().initializeTGeoDescription(EUTELESCOPE::GEOFILENAME,
                                             EUTELESCOPE::DUMPGEOROOT);

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

  streamlog_out(MESSAGE4) << "assumed beam energy " << _eBeam << " GeV" <<  endl;
  streamlog_out(MESSAGE4)     << "Summary:\n\t_planeRadLength has size: \t" << _planeRadLength.size() 
                            << "\n\t_planeWscatSi has size: \t" << _planeWscatSi.size() 
                            << "\n\t_planeWscatAir has size: \t" << _planeWscatAir.size() << '\n';


   // The measurement resolution depends on threshold, energy, dz, and number of iterations done ...
   // We make an accurate, heuristic  guess
   double distplane = _planePosition[1] - _planePosition[0];
   double res = 0;
   if( distplane > 100. ) {
    res = 100. - _eBeam*8;
    if (_eBeam > 11) {
        res = 8;
    }
   }
   else if( distplane < 30. ) {
    res = 20. - _eBeam*1.5; // if only tracks with prob(chi2,ndf) > 0.001 are passed to Mille
   }
   else if( distplane > 30. && distplane < 100. ) {
    res = 50 - _eBeam*3.8; 
    if (_eBeam > 10) {
        res = 6.;
    }
   }

   if(!_IsFirstAlignStep){
    // 2nd iteration 20
    if( distplane < 30. ) {
        res = 4.4 - _eBeam*0.1;
    }
    // 2nd iteration 150
    else if( distplane > 100. ) {
        res = 18 - _eBeam*1.5;
    }
   }

   res = res/1000.; // finally convert to [mm]
      
   streamlog_out( MESSAGE2 ) << "res x = " << res << endl;

   //For now we have the same x/y resolution and the same for telescope & DUT
   for(size_t ipl = 0; ipl < _nPlanes; ipl++) {
    _planeMeasPrec.emplace_back(1.0/res/res, 1.0/res/res); // precision = 1/resolution^2
   }

  // the user is giving sensor ids for the planes to be fixed.
  // These sensor ids have to be converted to a local index
  // according to the planes positions along the z axis, i.e. we
  // need to fill _FixedPlanes with the z-position which corresponds
  // to the postion in _sensorIDVec
  cout << "FixedPlanes " << _FixedPlanes_sensorIDs.size() << endl;
  for(auto& fixedPlaneID: _FixedPlanes_sensorIDs) {
    auto planeIt = std::find(_sensorIDVec.begin(), _sensorIDVec.end(), fixedPlaneID); 
    _FixedPlanes.emplace_back( planeIt - _sensorIDVec.begin() );
  }

  // same for excluded planes
  cout << "ExcludedPlanes " << _excludePlanes_sensorIDs.size() << endl;
  for(auto& excludedPlaneID: _excludePlanes_sensorIDs) {
    auto planeIt = std::find(_sensorIDVec.begin(), _sensorIDVec.end(), excludedPlaneID);
    _excludePlanes.emplace_back( planeIt - _sensorIDVec.begin() );
  }

  int counter = 0;
  for(auto& sensorID: _sensorIDVec) { 
    if(std::find(_excludePlanes_sensorIDs.begin(), _excludePlanes_sensorIDs.end(), sensorID) 
                 == _excludePlanes_sensorIDs.end()) {
        indexconverter.emplace_back(counter++);
    } else {
        indexconverter.emplace_back(-1);
    }
  }

#endif

  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;
  _printEventCounter = 0;

  // Initialize number of excluded planes
  _nExcludePlanes = _excludePlanes.size();

  streamlog_out( MESSAGE2 ) << "Number of planes excluded from the alignment fit: " << _nExcludePlanes << endl;

  // Initialise Mille statistics:
  _nMilleTracks = 0;

  // booking histograms
  bookHistos();

  streamlog_out( MESSAGE2 ) << "Initialising Mille..." << endl;

  unsigned int reserveSize = 8000;
  milleAlignGBL = std::make_unique<gbl::MilleBinary>( _binaryFilename, reserveSize );

  streamlog_out( MESSAGE2 ) << "The filename for the binary file is: " << _binaryFilename.c_str() << endl;

  // apply correction to cut if is not first alignment step
  if(!_IsFirstAlignStep){
    _triCut = _triCut*6./_eBeam;
    _driCut = _driCut*6./_eBeam;
    _sixCut = _sixCut*6./_eBeam;
  }

  if(_alignModeString.compare("XYShiftsRotZ") == 0 ) {
    _alignMode = Utility::alignMode::XYShiftsRotZ;
  } else if( _alignModeString.compare("XYShifts") == 0 ) {
    _alignMode = Utility::alignMode::XYShifts;
  } else if( _alignModeString.compare("XYZShiftsRotZ") == 0 ) {
    _alignMode = Utility::alignMode::XYZShiftsRotZ;
  } else {
    streamlog_out(ERROR) << "The chosen AlignMode: '" << _alignModeString << "' is invalid. Please correct your steering template and retry!" << std::endl;
    throw InvalidParameterException("AlignMode");
  }
  streamlog_out( MESSAGE2 ) << "end of init" << endl;
}

//------------------------------------------------------------------------------
void EUTelAlignGBL::processRunHeader( LCRunHeader* rdr ) {
  gblutil.setParent(this);
  gblutil.bookHistos();
  auto header = std::make_unique<EUTelRunHeaderImpl>(rdr);
  header->addProcessor( type() ) ;
  // increment the run counter
  ++_iRun;
}

//------------------------------------------------------------------------------
Eigen::Matrix<double,5,5> Jac55new( double ds ) {
  /* for GBL:
     Jacobian for straight line track
     track = q/p, x', y', x, y
     0,   1,  2,  3, 4
     */
  Eigen::Matrix<double,5,5> jac = Eigen::Matrix<double,5,5>::Identity();
  //jac.UnitMatrix();
  jac(3,1) = ds; // x = x0 + xp * ds
  jac(4,2) = ds; // y = y0 + yp * ds
  return jac;
}

//------------------------------------------------------------------------------
void EUTelAlignGBL::processEvent( LCEvent * event ) {

  if( _iEvt % 1000 == 0 ) {
    streamlog_out( MESSAGE2 ) << "Processing event "
      << setw(6) << setiosflags(ios::right)
      << event->getEventNumber() << " in run "
      << setw(6) << setiosflags(ios::right)
      << event->getRunNumber()
      << ", currently having "
      << _nMilleTracks << " tracks "
      << endl;
  }

  if( _nMilleTracks > (size_t)_maxTrackCandidatesTotal ) {
    throw StopProcessingException(this);
  }
  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event) ;

  if( evt->getEventType() == kEORE ) {
    streamlog_out( DEBUG2 ) << "EORE found: nothing else to do." << endl;
    return;
  }


  CellIDDecoder<TrackerHit> hitCellDecoder(EUTELESCOPE::HITENCODING);
  std::vector<EUTelTripletGBLUtility::hit> _hitsVec;

  for( size_t i = 0; i < _hitCollectionName.size(); i++ ) {
    LCCollection* collection = nullptr;
    try {
      collection = event->getCollection(_hitCollectionName[i]);
    } catch (DataNotAvailableException& e) {
      throw SkipEventException(this);
    }
    // loop over all hits in collection:

    if(_printEventCounter < NO_PRINT_EVENT_COUNTER) {
      streamlog_out(DEBUG2) << "Event " << event->getEventNumber() << " contains " 
      << collection->getNumberOfElements() << " hits" << endl;
    }
    
    nAllHitHistoGBLAlign->fill(collection->getNumberOfElements());

    for( int iHit = 0; iHit < collection->getNumberOfElements(); iHit++ ) {
      auto hit = static_cast<TrackerHitImpl*>( collection->getElementAt(iHit) );
      auto sensorID = hitCellDecoder(hit)["sensorID"];
      auto hitPosition = hit->getPosition();
      _hitsVec.emplace_back(hitPosition, sensorID);
      if(_printEventCounter < NO_PRINT_EVENT_COUNTER) std::cout << "Hit on plane " << sensorID << " at " << hitPosition[0] << "|" << hitPosition[1]  << '\n';
    } // end loop over all hits in given collection
  }//loop over all input hit collection

  // triplets 0-1-2:
  // driplets 3-4-5:
  int i0 = indexconverter[0]; // plane 0
  int i1 = indexconverter[1]; // plane 1
  int i2 = indexconverter[2]; // plane 2
  int i3 = indexconverter[3]; // plane 3
  int i4 = indexconverter[4]; // plane 4
  int i5 = indexconverter[5]; // plane 5
 
  if( i0*i1*i2*i3*i4*i5 >= 0 ) { // not excluded

    int nm = 0;
    auto tripletVec = std::vector<EUTelTripletGBLUtility::triplet>();
    auto dripletVec = std::vector<EUTelTripletGBLUtility::triplet>();
    gblutil.FindTriplets(_hitsVec, 0, 1, 2, _triCut, _slopeCut, tripletVec, false);
    gblutil.FindTriplets(_hitsVec, 3, 4, 5, _driCut, _slopeCut+0.012, dripletVec, false);

    if(_printEventCounter < NO_PRINT_EVENT_COUNTER){
      std::cout << "Triplets:\n";
      for(auto& triplet: tripletVec) {
        std::cout << triplet << '\n';
      }
      std::cout << "Driplets:\n";
      for(auto& driplet: dripletVec) {
        std::cout << driplet << '\n';
      }
    }

    ntriHistGBLAlign->fill( tripletVec.size()  );
    ndriHistGBLAlign->fill( dripletVec.size() );

    double zMid = 0.5*(_planePosition[3] + _planePosition[2]);
    auto matchedTripletVec = std::vector<EUTelTripletGBLUtility::track>();
    gblutil.MatchTriplets(tripletVec, dripletVec, zMid, _sixCut, matchedTripletVec);

    if(_printEventCounter < NO_PRINT_EVENT_COUNTER) std::cout << "Matched to:\n"; 
    for(auto& track: matchedTripletVec) {    
      //std::cout << "Dist 22:\n";
      auto& triplet = track.get_upstream();
      auto has22 = gblutil.AttachDUT(triplet, _hitsVec, 22, 3 );    
      //std::cout << "Dist 21:\n";
      auto& driplet = track.get_downstream();
      auto has21 = gblutil.AttachDUT(driplet, _hitsVec, 21, 3 );    
      if(_printEventCounter < NO_PRINT_EVENT_COUNTER){
         std::cout << "---pair--\n" << triplet << driplet << '\n';
        std::cout << "Expects hit on 22 at:\n";
        auto zPos = geo::gGeometry().siPlaneZPosition(22);
        auto trX = triplet.getx_at(zPos);
        auto trY = triplet.gety_at(zPos);
        std::cout << trX << "|" << trY << '\n';
        std::cout << "Expects hit on 21 at:\n";
        zPos = geo::gGeometry().siPlaneZPosition(21);
        trX = driplet.getx_at(zPos);
        trY = driplet.gety_at(zPos);
        std::cout << trX << "|" << trY << '\n';        
        for(auto& hit: _hitsVec) {
          if(hit.plane == 21 || hit.plane == 22) {
            std::cout << "Hit on " << hit.plane << " at: " << hit.x << "|" << hit.y << '\n';
          }
        }
      }
      if(has22 && has21) _nMilleTracks++;// std::cout << "Found triplet with both DUTs matched" << std::endl;
    }

    for(auto& track: matchedTripletVec) {
      // GBL point vector for the trajectory, all in [mm] !!
      // GBL with triplet A as seed
      std::vector<gbl::GblPoint> traj_points;
      // build up trajectory:
      std::vector<unsigned int> ilab; // 0-5 = telescope, 6 = DUT, 7 = REF
      vector<double> sPoint;

      // the arc length at the first measurement plane is 0.
      double s = 0;

      Eigen::Matrix2d proL2m = Eigen::Matrix2d::Identity();

      // scatter
      Eigen::Vector2d scat = Eigen::Vector2d::Zero(); //mean is zero

      auto triplet = track.get_upstream();
      auto driplet = track.get_downstream();
      //We need the _T_riplet slope to compute residual
      auto triSlope = triplet.slope();

      std::vector<EUTelTripletGBLUtility::hit> DUThits;
      for(auto it = triplet.DUT_begin(); it != triplet.DUT_end(); ++it){
        DUThits.emplace_back(it->second);
      }
      for(auto it = driplet.DUT_begin(); it != driplet.DUT_end(); ++it){
        DUThits.emplace_back(it->second);
      }

      if(_printEventCounter < NO_PRINT_EVENT_COUNTER) std::cout << "Track has " << DUThits.size() << " DUT hits\n";

      Eigen::Matrix2d alDer2; // alignment derivatives
      alDer2(0,0) = 1.0; // dx/dx GBL sign convetion
      alDer2(1,0) = 0.0; // dy/dx
      alDer2(0,1) = 0.0; // dx/dy
      alDer2(1,1) = 1.0; // dy/dy

      Eigen::Matrix<double,2,3> alDer3; // alignment derivatives
      alDer3(0,0) = 1.0; // dx/dx
      alDer3(1,0) = 0.0; // dy/dx
      alDer3(0,1) = 0.0; // dx/dy
      alDer3(1,1) = 1.0; // dy/dy

      Eigen::Matrix<double, 2,4> alDer4; // alignment derivatives
      alDer4(0,0) = 1.0; // dx/dx
      alDer4(1,0) = 0.0; // dy/dx
      alDer4(0,1) = 0.0; // dx/dy
      alDer4(1,1) = 1.0; // dy/dy
      alDer4(0,3) = triSlope.x; // dx/dz
      alDer4(1,3) = triSlope.y; // dx/dz

      size_t DUTCount = _nPlanes-6;
      std::vector<double> rx (_nPlanes, 0.0);
      std::vector<double> ry (_nPlanes, 0.0);
    
      double step = .0;
      unsigned int iLabel;

      for( size_t ipl = 0; ipl < _nPlanes; ++ipl ) {
        //We have to add all the planes, the up and downstream arm of the telescope will definitely have
        //hits, the DUTs might not though! The first and last three planes are the telescope.
        EUTelTripletGBLUtility::hit const * hit = nullptr;
        if(ipl < 3) {
          hit = &triplet.gethit(ipl);
        } else if( ipl < 3+DUTCount) {
          auto sensorID = _sensorIDVec[ipl];
          if(triplet.has_DUT(sensorID)) hit = &triplet.get_DUT_Hit(sensorID);
          else if(driplet.has_DUT(sensorID)) hit = &driplet.get_DUT_Hit(sensorID);
        } else {
          hit = &driplet.gethit(ipl-DUTCount);
        }

        //if there is no hit we take the plane position from the geo description
        double zz = hit ? hit->z : _planePosition[ipl];// [mm]

        //transport matrix in (q/p, x', y', x, y) space
        auto jacPointToPoint = Jac55new( step );
        auto point = gbl::GblPoint( jacPointToPoint );
        s += step;

        //of there is a hit we will add a measurement to the point
        if(hit){
          double xs = triplet.getx_at(zz);
          double ys = triplet.gety_at(zz);

          if(_printEventCounter < NO_PRINT_EVENT_COUNTER) std::cout << "xs = " << xs << "   ys = " << ys << std::endl;

          rx[ipl] = (hit->x - xs); // resid hit-triplet, in micrometer ...
          ry[ipl] = (hit->y - ys); // resid
      
          Eigen::Vector2d meas;
          meas[0] = rx[ipl]; // fill meas vector for GBL
          meas[1] = ry[ipl];
          point.addMeasurement( proL2m, meas, _planeMeasPrec[ipl] );

          if( _alignMode == Utility::alignMode::XYShifts ) { // only x and y shifts
            // global labels for MP:
            std::vector<int> globalLabels(2);
            globalLabels[0] = 1 + 2*ipl;
            globalLabels[1] = 2 + 2*ipl;
            point.addGlobals( globalLabels, alDer2 ); // for MillePede alignment
          }
          else if( _alignMode == Utility::alignMode::XYShiftsRotZ ) { // with rot
            std::vector<int> globalLabels(3);
            globalLabels[0] = 1 + 3*ipl; // x
            globalLabels[1] = 2 + 3*ipl; // y
            globalLabels[2] = 3 + 3*ipl; // rot
            alDer3(0,2) = -ys; // dx/dphi
            alDer3(1,2) =  xs; // dy/dphi
            point.addGlobals( globalLabels, alDer3 ); // for MillePede alignment
          }
          else if( _alignMode == Utility::alignMode::XYZShiftsRotZ ) { // with rot and z shift
            std::vector<int> globalLabels(4);
            globalLabels[0] = 1 + 4*ipl;
            globalLabels[1] = 2 + 4*ipl;
            globalLabels[2] = 3 + 4*ipl;
            globalLabels[3] = 4 + 4*ipl; // z
            alDer4(0,2) = -ys; // dx/dphi
            alDer4(1,2) =  xs; // dy/dphi
            point.addGlobals( globalLabels, alDer4 ); // for MillePede alignment
          }
          // monitor what we put into GBL:
          double xA = triplet.getx_at(zMid);
          double yA = triplet.gety_at(zMid);
          double xB = driplet.getx_at(zMid);
          double yB = driplet.gety_at(zMid);

          selxHistGBLAlign->fill( xA ); // triplet at mid
          selyHistGBLAlign->fill( yA );
          selaxHistGBLAlign->fill( triSlope.x*1E3 );//track slope
          selayHistGBLAlign->fill( triSlope.y*1E3 );
          seldxHistGBLAlign->fill( (xB-xA)*1E3 ); // triplet-driplet match
          seldyHistGBLAlign->fill( (yB-yA)*1E3 );
          selkxHistGBLAlign->fill( track.kink_x()*1E3 ); // triplet-driplet kink
          selkyHistGBLAlign->fill( track.kink_y()*1E3 );

          seldx1HistGBLAlign->fill( rx[1]*1E3 ); // triplet interpol
          seldy1HistGBLAlign->fill( ry[1]*1E3 );
          seldx3HistGBLAlign->fill( rx[3]*1E3 ); // triplet extrapol
          seldy3HistGBLAlign->fill( ry[3]*1E3 );
          seldx4HistGBLAlign->fill( rx[4]*1E3 );
          seldy4HistGBLAlign->fill( ry[4]*1E3 );
          seldx5HistGBLAlign->fill( rx[5]*1E3 );
          seldy5HistGBLAlign->fill( ry[5]*1E3 );
        }

        point.addScatterer( scat, _planeWscatSi[ipl] );
        sPoint.push_back( s );
        iLabel = sPoint.size();
        ilab.push_back(iLabel);
        traj_points.push_back(point);

        if( ipl < _nPlanes-1 ) {
          double distplane = _planePosition[ipl+1] - _planePosition[ipl];
          step = 0.21*distplane; // in [mm]
          auto point = gbl::GblPoint( Jac55new( step ) );
          point.addScatterer( scat, _planeWscatAir[ipl] );
          s += step;
          traj_points.push_back(point);
          sPoint.push_back( s );
          step = 0.58*distplane; // in [mm]
          auto point1 = gbl::GblPoint( Jac55new( step ) );
          point1.addScatterer( scat, _planeWscatAir[ipl] );
          s += step;
          traj_points.push_back(point1);
          sPoint.push_back( s );
          step = 0.21*distplane; // remaing distance to next plane, in [mm]
        }
      } // loop over planes

      double Chi2;
      int Ndf;
      double lostWeight;

      gbl::GblTrajectory traj(traj_points, false); // curvature = false
      traj.fit( Chi2, Ndf, lostWeight );
      //traj.getLabels(ilab); // instead pushback sPoint.size() when adding plane

      if(_printEventCounter < NO_PRINT_EVENT_COUNTER){
        streamlog_out(DEBUG4) << "traj with " << traj.getNumPoints() << " points:" << endl;
        for( size_t ipl = 0; ipl < sPoint.size(); ++ipl ){
          streamlog_out(DEBUG4) << "  GBL point " << ipl;
          streamlog_out(DEBUG4) << "  z " << sPoint[ipl]; 
          streamlog_out(DEBUG4) << endl;
        }
        for( size_t ipl = 0; ipl < 6; ++ipl ){
          streamlog_out(DEBUG4) << " plane " << ipl << ", lab " << ilab[ipl];
          streamlog_out(DEBUG4) << " z " << sPoint[ilab[ipl]-1];
          streamlog_out(DEBUG4) << "  dx " << rx[ipl];
          streamlog_out(DEBUG4) << "  dy " << ry[ipl];
          streamlog_out(DEBUG4) << "  chi2 " << Chi2;
          streamlog_out(DEBUG4) << "  ndf " << Ndf;
          streamlog_out(DEBUG4) << endl;
        }

        streamlog_out(DEBUG2)  << " Is traj valid? " << traj.isValid() << std::endl;
        traj.printPoints();
        //traj.printTrajectory();
        //traj.printData();
        _printEventCounter++;
      }

      //cout << " chi2 " << Chi2 << ", ndf " << Ndf << endl;

      gblndfHistGBLAlign->fill( Ndf );
      gblchi2HistGBLAlign->fill( Chi2 );
      double probchi = TMath::Prob( Chi2, Ndf );
      gblprbHistGBLAlign->fill( probchi );

      // bad fits:
      if( probchi < 0.01 ) {
        //badxHistGBLAlign->fill( xA ); // triplet at DUT
        //badyHistGBLAlign->fill( yA );
        badaxHistGBLAlign->fill( triSlope.x*1E3 );
        badayHistGBLAlign->fill( triSlope.y*1E3 );
        //baddxHistGBLAlign->fill( (xB-xA)*1E3 ); // triplet-driplet match
        //baddyHistGBLAlign->fill( (yB-yA)*1E3 );
        badkxHistGBLAlign->fill( track.kink_x()*1E3 ); // triplet-driplet kink
        badkyHistGBLAlign->fill( track.kink_y()*1E3 );

        baddx1HistGBLAlign->fill( rx[1]*1E3 ); // triplet interpol
        baddy1HistGBLAlign->fill( ry[1]*1E3 );
        baddx3HistGBLAlign->fill( rx[3]*1E3 ); // triplet extrapol
        baddy3HistGBLAlign->fill( ry[3]*1E3 );
        baddx4HistGBLAlign->fill( rx[4]*1E3 );
        baddy4HistGBLAlign->fill( ry[4]*1E3 );
        baddx5HistGBLAlign->fill( rx[5]*1E3 );
        baddy5HistGBLAlign->fill( ry[5]*1E3 );
      }// bad fit
      else {
        goodx1HistGBLAlign->fill( rx[1]*1E3 ); // triplet interpol
        goody1HistGBLAlign->fill( ry[1]*1E3 );
      } // OK fit

      // look at fit:
      Eigen::VectorXd aCorrection;
      Eigen::MatrixXd aCovariance;

      std::array<double,8> ax;
      std::array<double,8> ay;
      unsigned int k = 0;

      // at plane 0:
      unsigned int ndata = 2;
      unsigned int ndim = 2;
      Eigen::VectorXd aResiduals(ndim);
      Eigen::VectorXd aMeasErrors(ndim);
      Eigen::VectorXd aResErrors(ndim);
      Eigen::VectorXd aDownWeights(ndim);

      int ipos = ilab[0];
      traj.getResults( ipos, aCorrection, aCovariance );
      traj.getMeasResults( ipos, ndata, aResiduals, aMeasErrors, aResErrors, aDownWeights );

      //track = q/p, x', y', x, y
      //        0,   1,  2,  3, 4
      gblax0HistGBLAlign->fill( aCorrection[1]*1E3 ); // angle x [mrad]
      gbldx0HistGBLAlign->fill( aCorrection[3]*1e3 ); // shift x [um]
      gblrx0HistGBLAlign->fill( (aResiduals[0])*1e3 ); // residual x [um]
      gblry0HistGBLAlign->fill( (aResiduals[1])*1e3 ); // residual y [um]
      gblpx0HistGBLAlign->fill( (aResiduals[0] / aResErrors[0]) );
      gblpy0HistGBLAlign->fill( (aResiduals[1] / aResErrors[1]) ); 
      ax[k] = aCorrection[1]; // angle correction at plane, for kinks
      ay[k] = aCorrection[2]; // angle correction at plane, for kinks
      k++;

      ipos = ilab[1];
      traj.getResults( ipos, aCorrection, aCovariance );
      gblax1HistGBLAlign->fill( aCorrection[1]*1E3 ); // angle x [mrad]
      gbldx1HistGBLAlign->fill( aCorrection[3]*1e3 ); // shift x [um]
      gblrx1HistGBLAlign->fill( (rx[1] - aCorrection[3])*1e3 ); // residual x [um]
      gblry1HistGBLAlign->fill( (ry[1] - aCorrection[4])*1e3 ); // residual y [um]
      ax[k] = aCorrection[1]; // angle correction at plane, for kinks
      ay[k] = aCorrection[2]; // angle correction at plane, for kinks
      k++;

      ipos = ilab[2];
      traj.getResults( ipos, aCorrection, aCovariance );
      gblax2HistGBLAlign->fill( aCorrection[1]*1E3 ); // angle x [mrad]
      gbldx2HistGBLAlign->fill( aCorrection[3]*1e3 ); // shift x [um]
      gblrx2HistGBLAlign->fill( (rx[2] - aCorrection[3])*1e3 ); // residual x [um]
      gblry2HistGBLAlign->fill( (ry[2] - aCorrection[4])*1e3 ); // residual y [um]
      ax[k] = aCorrection[1]; // angle correction at plane, for kinks
      ay[k] = aCorrection[2]; // angle correction at plane, for kinks
      k++;

      ipos = ilab[3];
      traj.getResults( ipos, aCorrection, aCovariance );
      gblax3HistGBLAlign->fill( aCorrection[1]*1E3 ); // angle x [mrad]
      gbldx3HistGBLAlign->fill( aCorrection[3]*1e3 ); // shift x [um]
      gblrx3HistGBLAlign->fill( (rx[3] - aCorrection[3]*1e3) ); // residual x [um]
      gblry3HistGBLAlign->fill( (ry[3] - aCorrection[4])*1e3 ); // residual y [um]
      ax[k] = aCorrection[1]; // angle correction at plane, for kinks
      ay[k] = aCorrection[2]; // angle correction at plane, for kinks
      k++;

      ipos = ilab[4];
      traj.getResults( ipos, aCorrection, aCovariance );
      gblax4HistGBLAlign->fill( aCorrection[1]*1E3 ); // angle x [mrad]
      gbldx4HistGBLAlign->fill( aCorrection[3]*1e3 ); // shift x [um]
      gblrx4HistGBLAlign->fill( (rx[4] - aCorrection[3])*1e3 ); // residual x [um]
      gblry4HistGBLAlign->fill( (ry[4] - aCorrection[4])*1e3 ); // residual y [um]
      ax[k] = aCorrection[1]; // angle correction at plane, for kinks
      ay[k] = aCorrection[2]; // angle correction at plane, for kinks
      k++;

      ipos = ilab[5];
      traj.getResults( ipos, aCorrection, aCovariance );
      gblax5HistGBLAlign->fill( aCorrection[1]*1E3 ); // angle x [mrad]
      gbldx5HistGBLAlign->fill( aCorrection[3]*1e3 ); // shift x [um]
      gblrx5HistGBLAlign->fill( (rx[5] - aCorrection[3])*1e3 ); // residual x [um]
      gblry5HistGBLAlign->fill( (ry[5] - aCorrection[4])*1e3 ); // residual y [um]
      ax[k] = aCorrection[1]; // angle correction at plane, for kinks
      ay[k] = aCorrection[2]; // angle correction at plane, for kinks
      k++;

      // kinks: 1,2 = tele, 3 = DUT, 4,5 = tele
      gblkx1HistGBLAlign->fill( (ax[1] - ax[0])*1E3 ); // kink at 1 [mrad]
      gblkx2HistGBLAlign->fill( (ax[2] - ax[1])*1E3 ); // kink at 2 [mrad]
      gblkx3HistGBLAlign->fill( (ax[3] - ax[2])*1E3 ); // kink at 3 [mrad]
      gblkx4HistGBLAlign->fill( (ax[4] - ax[3])*1E3 ); // kink at 4 [mrad]
      gblkx5HistGBLAlign->fill( (ax[5] - ax[4])*1E3 ); // kink at 5 [mrad]

      // do not pass very bad tracks to mille
      if(probchi > 0.001) {
          traj.milleOut( *milleAlignGBL );
             nm++;
      }
    }
    nmHistGBLAlign->fill( nm );
    _nMilleTracks += nm;
  }//i3*i4*i5 valid
  // count events:
  ++_iEvt;
  if( isFirstEvent() ) _isFirstEvent = false;
}

//------------------------------------------------------------------------------
void EUTelAlignGBL::end() {
  milleAlignGBL.reset(nullptr);

  // if write the pede steering file
  if( _generatePedeSteerfile ) {

    streamlog_out( MESSAGE4 ) << endl << "Generating the steering file for the pede program..." << endl;

    ofstream steerFile;
    steerFile.open(_pedeSteerfileName.c_str());

    if( steerFile.is_open()) {

      // find first and last excluded plane
      auto firstnotexcl = _nPlanes;
      size_t lastnotexcl = 0;

      // loop over all planes:

      for( size_t ipl = 0; ipl < _nPlanes; ipl++) {

    size_t excluded = 0;

    // loop over excluded planes:

    for( size_t jpl = 0; jpl < _nExcludePlanes; jpl++ ) {
      if( ipl == _excludePlanes[jpl] ) excluded = 1;
    }

    if( excluded == 0 && firstnotexcl > ipl ) firstnotexcl = ipl;

    if( excluded == 0 && lastnotexcl < ipl ) lastnotexcl = ipl;

      } // end loop over all planes

      steerFile << "Cfiles" << endl;
      steerFile << _binaryFilename << endl;
      steerFile << endl;

      steerFile << "Parameter" << endl;

      int counter = 0;
      int nfix = 0;

      // loop over all planes:

      for( size_t ipl = 0; ipl < _nPlanes; ipl++) {

    int excluded = 0; // flag for excluded planes

    // check in list of excluded planes:

    for( size_t iex = 0; iex < _nExcludePlanes; iex++) {
      if( ipl == _excludePlanes[iex] )
        excluded = 1;
    }

    cout << "Plane " << ipl << " exclude = " << excluded << endl;

    if( excluded == 0 ) {

      bool fixed = false;
      for( size_t i = 0;i< _FixedPlanes.size(); i++ ) {
        if( _FixedPlanes[i] == ipl )
          fixed = true;
      }

      // if fixed planes

      if( fixed || (_FixedPlanes.size() == 0 && (ipl == firstnotexcl || ipl == lastnotexcl) ) ) {
        nfix++;
        if( _alignMode == Utility::alignMode::XYShifts ) {
          steerFile << (counter * 2 + 1) << "  0.0 -1.0" << endl;
          steerFile << (counter * 2 + 2) << "  0.0 -1.0" << endl;
        }
        if( _alignMode == Utility::alignMode::XYShiftsRotZ ) {
          steerFile << (counter * 3 + 1) << "  0.0 -1.0" << endl; // fix x
          steerFile << (counter * 3 + 2) << "  0.0 -1.0" << endl; // fix y
          steerFile << (counter * 3 + 3) << "  0.0 -1.0" << endl; // fix rot
        }
        if( _alignMode == Utility::alignMode::XYZShiftsRotZ ) {
          steerFile << (counter * 4 + 1) << "  0.0 -1.0" << endl;
          steerFile << (counter * 4 + 2) << "  0.0 -1.0" << endl;
          steerFile << (counter * 4 + 3) << "  0.0 -1.0" << endl;
        }
      }

      else {

        if( _alignMode == Utility::alignMode::XYShifts ) {
          steerFile << (counter * 2 + 1) << "  0.0  0.0" << endl;
          steerFile << (counter * 2 + 2) << "  0.0  0.0" << endl;
        }

        if( _alignMode == Utility::alignMode::XYShiftsRotZ ) {
          steerFile << (counter * 3 + 1) << "  0.0  0.0" << endl;
          steerFile << (counter * 3 + 2) << "  0.0  0.0" << endl;
          steerFile << (counter * 3 + 3) << "  0.0  0.0" << endl;
        }

        if( _alignMode == Utility::alignMode::XYZShiftsRotZ ) {
          steerFile << (counter * 4 + 1) << "  0.0  0.0" << endl;
          steerFile << (counter * 4 + 2) << "  0.0  0.0" << endl;
          steerFile << (counter * 4 + 3) << "  0.0  0.0" << endl;
        }

      }// not fixed

      // special for z shift:

      if( _alignMode == Utility::alignMode::XYZShiftsRotZ ) {
        if( ipl == 1 )
          steerFile << (counter * 4 + 4) << "  0.0 -1.0" << endl;
        else if( ipl == 4 )
          steerFile << (counter * 4 + 4) << "  0.0 -1.0" << endl;
        else
          steerFile << (counter * 4 + 4) << "  0.0  0.0" << endl;
      }

      counter++;

    } // end if plane not excluded

      } // end loop over all planes

      if( nfix < 2 ) {

    if( _alignMode == Utility::alignMode::XYShifts ) {

      steerFile << "Constraint 0 ! sum dx = 0" << endl;
      steerFile << " 1  1.0" << endl;
      steerFile << " 3  1.0" << endl;
      steerFile << " 5  1.0" << endl;
      steerFile << " 7  1.0" << endl;
      steerFile << " 9  1.0" << endl;
      steerFile << "11  1.0" << endl;

      steerFile << "Constraint 0 ! sum dy = 0" << endl;
      steerFile << " 2  1.0" << endl;
      steerFile << " 4  1.0" << endl;
      steerFile << " 6  1.0" << endl;
      steerFile << " 8  1.0" << endl;
      steerFile << "10  1.0" << endl;
      steerFile << "12  1.0" << endl;
    }

    if( _alignMode == Utility::alignMode::XYShiftsRotZ ) {

      steerFile << "Constraint 0 ! sum dx = 0" << endl;
      steerFile << " 1  1.0" << endl;
      steerFile << " 4  1.0" << endl;
      steerFile << " 7  1.0" << endl;
      steerFile << "10  1.0" << endl;
      steerFile << "13  1.0" << endl;
      steerFile << "16  1.0" << endl;

      steerFile << "Constraint 0 ! sum dy = 0" << endl;
      steerFile << " 2  1.0" << endl;
      steerFile << " 5  1.0" << endl;
      steerFile << " 8  1.0" << endl;
      steerFile << "11  1.0" << endl;
      steerFile << "14  1.0" << endl;
      steerFile << "17  1.0" << endl;

      steerFile << "Constraint 0 ! sum dphi = 0" << endl;
      steerFile << " 3  1.0" << endl;
      steerFile << " 6  1.0" << endl;
      steerFile << " 9  1.0" << endl;
      steerFile << "12  1.0" << endl;
      steerFile << "15  1.0" << endl;
      steerFile << "18  1.0" << endl;
    }

    if( _alignMode == Utility::alignMode::XYZShiftsRotZ ) {

      steerFile << "Constraint 0 ! sum dx = 0" << endl;
      steerFile << " 1  1.0" << endl;
      steerFile << " 5  1.0" << endl;
      steerFile << " 9  1.0" << endl;
      steerFile << "13  1.0" << endl;
      steerFile << "17  1.0" << endl;
      steerFile << "21  1.0" << endl;

      steerFile << "Constraint 0 ! sum dy = 0" << endl;
      steerFile << " 2  1.0" << endl;
      steerFile << " 6  1.0" << endl;
      steerFile << "10  1.0" << endl;
      steerFile << "14  1.0" << endl;
      steerFile << "18  1.0" << endl;
      steerFile << "22  1.0" << endl;

      steerFile << "Constraint 0 ! sum dphi = 0" << endl;
      steerFile << " 3  1.0" << endl;
      steerFile << " 7  1.0" << endl;
      steerFile << "11  1.0" << endl;
      steerFile << "15  1.0" << endl;
      steerFile << "19  1.0" << endl;
      steerFile << "23  1.0" << endl;

      steerFile << "Constraint 0 ! sum dz = 0" << endl;
      steerFile << " 4  1.0" << endl;
      steerFile << " 8  1.0" << endl;
      steerFile << "12  1.0" << endl;
      steerFile << "16  1.0" << endl;
      steerFile << "20  1.0" << endl;
      steerFile << "24  1.0" << endl;
    }

      }//nfix < 2

      steerFile << endl;
      steerFile << "! chiscut 5.0 2.5" << endl;
      steerFile << "outlierdownweighting 4" << endl;
      steerFile << "dwfractioncut 0.1" << endl;
      steerFile << endl;
      steerFile << "method inversion 10  0.1" << endl;
      // Use 10 OpenMP threads to process the data:
      steerFile << "threads 10 1" << endl;
      steerFile << endl;
      steerFile << "histprint" << endl;
      steerFile << endl;
      steerFile << "end" << endl;

      steerFile.close();

      streamlog_out( MESSAGE2 ) << "File " << _pedeSteerfileName << " written." << endl;

    }
    else {
      streamlog_out( ERROR2 ) << "Could not open steering file." << endl;
    }

  } // end if write the pede steering file

  streamlog_out( MESSAGE2 ) << endl;
  streamlog_out( MESSAGE2 ) << "Number of tracks used: " << _nMilleTracks << endl;

  streamlog_out( MESSAGE2 ) << endl;
  streamlog_out( MESSAGE2 ) << "Successfully finished" << endl;

}//end

//------------------------------------------------------------------------------
void EUTelAlignGBL::bookHistos() {

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  try {
    streamlog_out( MESSAGE2 ) << "Booking histograms..." << endl;

    //DP:
    
    nAllHitHistoGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "nallhit", 201, -0.5, 200.5 );
    nAllHitHistoGBLAlign->setTitle( "Telescope hits/event;telescope hits;events" );

    tridxHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "tridx", 100, -1000, 1000 );
    tridxHistGBLAlign->setTitle( "tridx;x_{1}-x_{02} [um];triplet" );

    tridyHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "tridy", 100, -1000, 1000 );
    tridyHistGBLAlign->setTitle( "tridy;y_{1}-y_{02} [um];triplet" );

    tridxcHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "tridxc", 100, -100, 100 );
    tridxcHistGBLAlign->setTitle( "tridxc;x_{1}-x_{02} [um];triplet" );

    tridycHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "tridyc", 100, -100, 100 );
    tridycHistGBLAlign->setTitle( "tridyc;y_{1}-y_{02} [um];triplet" );

    ntriHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "ntri", 21, -0.5, 20.5 );
    ntriHistGBLAlign->setTitle( "ntri;triplets;events" );

    dridxHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "dridx", 100, -1000, 1000 );
    dridxHistGBLAlign->setTitle( "dridx;x_{4}-x_{35} [um];triplet" );

    dridyHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "dridy", 100, -1000, 1000 );
    dridyHistGBLAlign->setTitle( "dridy;y_{4}-y_{35} [um];triplet" );

    dridxcHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "dridxc", 100, -100, 100 );
    dridxcHistGBLAlign->setTitle( "dridxc;x_{4}-x_{35} [um];triplet" );

    dridycHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "dridyc", 100, -100, 100 );
    dridycHistGBLAlign->setTitle( "dridyc;y_{4}-y_{35} [um];triplet" );

    ndriHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "ndri", 21, -0.5, 20.5 );
    ndriHistGBLAlign->setTitle( "ndri;driplets;events" );

    sixdxHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "sixdx", 250, -1000, 1000 );
    sixdxHistGBLAlign->setTitle( "sixdx;x_{A}-x_{B} [um];triplet pairs" );

    sixdyHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "sixdy", 250, -1000, 1000 );
    sixdyHistGBLAlign->setTitle( "sixdy;y_{A}-y_{B} [um];triplet pairs" );

    sixkxHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "sixkx", 100, -25, 25 );
    sixkxHistGBLAlign->setTitle( "sixkx;x kink angle [mrad];triplet pairs" );

    sixkyHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "sixky", 100, -25, 25 );
    sixkyHistGBLAlign->setTitle( "sixky;y kink angle [mrad];triplet pairs" );

    // GBL:

    selxHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "selx", 150, -15, 15 );
    selxHistGBLAlign->setTitle( "x at DUT, sel GBL;x [mm];tracks" );

    selyHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "sely", 100, -10, 10 );
    selyHistGBLAlign->setTitle( "y at DUT, sel GBL;y [mm];tracks" );

    selaxHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "selax", 100, -25, 25 );
    selaxHistGBLAlign->setTitle( "track angle x, sel GBL;x angle [mrad];tracks" );

    selayHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "selay", 100, -25, 25 );
    selayHistGBLAlign->setTitle( "track angle y, sel GBL;y angle [mrad];tracks" );

    seldxHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "seldx", 100, -5000, 5000 );
    seldxHistGBLAlign->setTitle( "track match x, sel GBL;#Deltax [#mum];tracks" );

    seldyHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "seldy", 100, -5000, 5000 );
    seldyHistGBLAlign->setTitle( "track match y, sel GBL;#Deltay [#mum];tracks" );

    selkxHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "selkx", 100, -25, 25 );
    selkxHistGBLAlign->setTitle( "kink x, sel GBL;kink x [mrad];tracks" );

    selkyHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "selky", 100, -25, 25 );
    selkyHistGBLAlign->setTitle( "kink y, sel GBL;kink y [mrad];tracks" );

    seldx1HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "seldx1", 100, -1000, 1000 );
    seldx1HistGBLAlign->setTitle( "triplet resid x at 1, sel GBL;#Deltax [#mum];tracks" );

    seldy1HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "seldy1", 100, -1000, 1000 );
    seldy1HistGBLAlign->setTitle( "triplet resid y at 1, sel GBL;#Deltay [#mum];tracks" );

    seldx3HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "seldx3", 100, -1000, 1000 );
    seldx3HistGBLAlign->setTitle( "triplet resid x at 3, sel GBL;#Deltax [#mum];tracks" );

    seldy3HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "seldy3", 100, -1000, 1000 );
    seldy3HistGBLAlign->setTitle( "triplet resid y at 3, sel GBL;#Deltay [#mum];tracks" );

    seldx4HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "seldx4", 100, -1000, 1000 );
    seldx4HistGBLAlign->setTitle( "triplet resid x at 4, sel GBL;#Deltax [#mum];tracks" );

    seldy4HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "seldy4", 100, -1000, 1000 );
    seldy4HistGBLAlign->setTitle( "triplet resid y at 4, sel GBL;#Deltay [#mum];tracks" );

    seldx5HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "seldx5", 100, -5000, 5000 );
    seldx5HistGBLAlign->setTitle( "triplet resid x at 5, sel GBL;#Deltax [#mum];tracks" );

    seldy5HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "seldy5", 100, -5000, 5000 );
    seldy5HistGBLAlign->setTitle( "triplet resid y at 5, sel GBL;#Deltay [#mum];tracks" );

    seldx6HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "seldx6", 100, -5000, 5000 );
    seldx6HistGBLAlign->setTitle( "triplet resid x at DUT, sel GBL;#Deltax [#mum];tracks" );

    seldy6HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "seldy6", 100, -5000, 5000 );
    seldy6HistGBLAlign->setTitle( "triplet resid y at DUT, sel GBL;#Deltay [#mum];tracks" );

    gblndfHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblndf", 16, -0.5, 15.5 );
    gblndfHistGBLAlign->setTitle( "GBL fit NDF;GBL NDF;tracks" );

    gblchi2HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblchi2", 100, 0, 100 );
    gblchi2HistGBLAlign->setTitle( "GBL fit chi2;GBL chi2;tracks" );

    gblprbHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblprb", 100, 0, 1 );
    gblprbHistGBLAlign->setTitle( "GBL fit probability;GBL fit probability;tracks" );

    // bad fits:

    badxHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "badx", 150, -15, 15 );
    badxHistGBLAlign->setTitle( "x at DUT, bad GBL;x [mm];tracks" );

    badyHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "bady", 100, -10, 10 );
    badyHistGBLAlign->setTitle( "y at DUT, bad GBL;y [mm];tracks" );

    badaxHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "badax", 100, -25, 25 );
    badaxHistGBLAlign->setTitle( "track angle x, bad GBL;x angle [mrad];tracks" );

    badayHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baday", 100, -25, 25 );
    badayHistGBLAlign->setTitle( "track angle y, bad GBL;y angle [mrad];tracks" );

    baddxHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baddx", 100, -5000, 5000 );
    baddxHistGBLAlign->setTitle( "track match x, bad GBL;#Deltax [#mum];tracks" );

    baddyHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baddy", 100, -5000, 5000 );
    baddyHistGBLAlign->setTitle( "track match y, bad GBL;#Deltay [#mum];tracks" );

    badkxHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "badkx", 100, -25, 25 );
    badkxHistGBLAlign->setTitle( "kink x, bad GBL;kink x [mrad];tracks" );

    badkyHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "badky", 100, -25, 25 );
    badkyHistGBLAlign->setTitle( "kink y, bad GBL;kink y [mrad];tracks" );

    baddx1HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baddx1", 100, -1000, 1000 );
    baddx1HistGBLAlign->setTitle( "triplet resid x at 1, bad GBL;#Deltax [#mum];tracks" );

    baddy1HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baddy1", 100, -1000, 1000 );
    baddy1HistGBLAlign->setTitle( "triplet resid y at 1, bad GBL;#Deltay [#mum];tracks" );

    baddx3HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baddx3", 100, -1000, 1000 );
    baddx3HistGBLAlign->setTitle( "triplet resid x at 3, bad GBL;#Deltax [#mum];tracks" );

    baddy3HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baddy3", 100, -1000, 1000 );
    baddy3HistGBLAlign->setTitle( "triplet resid y at 3, bad GBL;#Deltay [#mum];tracks" );

    baddx4HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baddx4", 100, -1500, 1500 );
    baddx4HistGBLAlign->setTitle( "triplet resid x at 4, bad GBL;#Deltax [#mum];tracks" );

    baddy4HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baddy4", 100, -1500, 1500 );
    baddy4HistGBLAlign->setTitle( "triplet resid y at 4, bad GBL;#Deltay [#mum];tracks" );

    baddx5HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baddx5", 100, -3000, 3000 );
    baddx5HistGBLAlign->setTitle( "triplet resid x at 5, bad GBL;#Deltax [#mum];tracks" );

    baddy5HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baddy5", 100, -3000, 3000 );
    baddy5HistGBLAlign->setTitle( "triplet resid y at 5, bad GBL;#Deltay [#mum];tracks" );

    baddx6HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baddx6", 100, -250, 250 );
    baddx6HistGBLAlign->setTitle( "triplet resid x at DUT, bad GBL;#Deltax [#mum];tracks" );

    baddy6HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baddy6", 100, -250, 250 );
    baddy6HistGBLAlign->setTitle( "triplet resid y at DUT, bad GBL;#Deltay [#mum];tracks" );

    // good fits:

    goodx1HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "goodx1", 100, -1000, 1000 );
    goodx1HistGBLAlign->setTitle( "triplet resid x at 1, good GBL;#Deltax [#mum];tracks" );

    goody1HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "goody1", 100, -1000, 1000 );
    goody1HistGBLAlign->setTitle( "triplet resid y at 1, good GBL;#Deltay [#mum];tracks" );

    goodx6HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "goodx6", 100, -250, 250 );
    goodx6HistGBLAlign->setTitle( "triplet resid x at 6, good GBL;#Deltax [#mum];tracks" );

    goody6HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "goody6", 100, -250, 250 );
    goody6HistGBLAlign->setTitle( "triplet resid y at 6, good GBL;#Deltay [#mum];tracks" );

    // look at fit:

    gblax0HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblax0", 100, -5, 5 );
    gblax0HistGBLAlign->setTitle( "GBL angle at plane 0;x angle at plane 0 [mrad];tracks" );

    gbldx0HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gbldx0", 500, -250, 250 );
    gbldx0HistGBLAlign->setTitle( "GBL shift at plane 0;x shift at plane 0 [#mum];tracks" );

    gblrx0HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblrx0", 500, -250, 250 );
    gblrx0HistGBLAlign->setTitle( "GBL resid at plane 0;x resid at plane 0 [#mum];tracks" );

    gblry0HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblry0", 500, -250, 250 );
    gblry0HistGBLAlign->setTitle( "GBL resid at plane 0;y resid at plane 0 [#mum];tracks" );

    gblpx0HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblpx0", 100, -5, 5 );
    gblpx0HistGBLAlign->setTitle( "GBL pull at plane 0;x pull at plane 0;tracks" );

    gblpy0HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblpy0", 100, -5, 5 );
    gblpy0HistGBLAlign->setTitle( "GBL pull at plane 0;y pull at plane 0;tracks" );


    gblax1HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblax1", 100, -5, 5 );
    gblax1HistGBLAlign->setTitle( "GBL angle at plane 1;x angle at plane 1 [mrad];tracks" );

    gbldx1HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gbldx1", 500, -250, 250 );
    gbldx1HistGBLAlign->setTitle( "GBL shift at plane 1;x shift at plane 1 [#mum];tracks" );

    gblrx1HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblrx1", 500, -250, 250 );
    gblrx1HistGBLAlign->setTitle( "GBL resid at plane 1;x resid at plane 1 [#mum];tracks" );

    gblry1HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblry1", 500, -250, 250 );
    gblry1HistGBLAlign->setTitle( "GBL resid at plane 1;y resid at plane 1 [#mum];tracks" );


    gblax2HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblax2", 100, -5, 5 );
    gblax2HistGBLAlign->setTitle( "GBL angle at plane 2;x angle at plane 2 [mrad];tracks" );

    gbldx2HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gbldx2", 500, -250, 250 );
    gbldx2HistGBLAlign->setTitle( "GBL shift at plane 2;x shift at plane 2 [#mum];tracks" );

    gblrx2HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblrx2", 500, -250, 250 );
    gblrx2HistGBLAlign->setTitle( "GBL resid at plane 2;x resid at plane 2 [#mum];tracks" );

    gblry2HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblry2", 500, -250, 250 );
    gblry2HistGBLAlign->setTitle( "GBL resid at plane 2;y resid at plane 2 [#mum];tracks" );


    gblax3HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblax3", 100, -5, 5 );
    gblax3HistGBLAlign->setTitle( "GBL angle at plane 3;x angle at plane 3 [mrad];tracks" );

    gbldx3HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gbldx3", 500, -250, 250 );
    gbldx3HistGBLAlign->setTitle( "GBL shift at plane 3;x shift at plane 3 [#mum];tracks" );

    gblrx3HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblrx3", 500, -250, 250 );
    gblrx3HistGBLAlign->setTitle( "GBL resid at plane 3;x resid at plane 3 [#mum];tracks" );

    gblry3HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblry3", 500, -250, 250 );
    gblry3HistGBLAlign->setTitle( "GBL resid at plane 3;y resid at plane 3 [#mum];tracks" );


    gblax4HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblax4", 100, -5, 5 );
    gblax4HistGBLAlign->setTitle( "GBL angle at plane 4;x angle at plane 4 [mrad];tracks" );

    gbldx4HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gbldx4", 500, -250, 250 );
    gbldx4HistGBLAlign->setTitle( "GBL shift at plane 4;x shift at plane 4 [#mum];tracks" );

    gblrx4HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblrx4", 500, -250, 250 );
    gblrx4HistGBLAlign->setTitle( "GBL resid at plane 4;x resid at plane 4 [#mum];tracks" );

    gblry4HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblry4", 500, -250, 250 );
    gblry4HistGBLAlign->setTitle( "GBL resid at plane 4;y resid at plane 4 [#mum];tracks" );


    gblax5HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblax5", 100, -5, 5 );
    gblax5HistGBLAlign->setTitle( "GBL angle at plane 5;x angle at plane 5 [mrad];tracks" );

    gbldx5HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gbldx5", 500, -250, 250 );
    gbldx5HistGBLAlign->setTitle( "GBL shift at plane 5;x shift at plane 5 [#mum];tracks" );

    gblrx5HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblrx5", 500, -250, 250 );
    gblrx5HistGBLAlign->setTitle( "GBL resid at plane 5;x resid at plane 5 [#mum];tracks" );

    gblry5HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblry5", 500, -250, 250 );
    gblry5HistGBLAlign->setTitle( "GBL resid at plane 5;y resid at plane 5 [#mum];tracks" );


    gblax6HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblax6", 100, -5, 5 );
    gblax6HistGBLAlign->setTitle( "GBL angle at DUT;x angle at DUT [mrad];tracks" );

    gbldx6HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gbldx6", 100, -250, 250 );
    gbldx6HistGBLAlign->setTitle( "GBL shift at DUT;x shift at DUT [#mum];tracks" );

    gbldy6HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gbldy6", 100, -1000, 1000 );
    gbldy6HistGBLAlign->setTitle( "GBL shift at DUT;y shift at DUT [#mum];tracks" );

    gblrx6HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblrx6", 100, -250, 250 );
    gblrx6HistGBLAlign->setTitle( "GBL resid at DUT;x resid at DUT [#mum];tracks" );

    gblry6HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblry6", 100, -250, 250 );
    gblry6HistGBLAlign->setTitle( "GBL resid at DUT;y resid at DUT [#mum];tracks" );


    gblkx1HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblkx1", 100, -5, 5 );
    gblkx1HistGBLAlign->setTitle( "GBL kink angle at plane 1;plane 1 kink [mrad];tracks" );

    gblkx2HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblkx2", 100, -5, 5 );
    gblkx2HistGBLAlign->setTitle( "GBL kink angle at plane 2;plane 2 kink [mrad];tracks" );

    gblkx3HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblkx3", 100, -5, 5 );
    gblkx3HistGBLAlign->setTitle( "GBL kink angle at plane 3;plane 3 kink [mrad];tracks" );

    gblkx4HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblkx4", 100, -5, 5 );
    gblkx4HistGBLAlign->setTitle( "GBL kink angle at plane 4;plane 4 kink [mrad];tracks" );

    gblkx5HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblkx5", 100, -5, 5 );
    gblkx5HistGBLAlign->setTitle( "GBL kink angle at plane 5;plane 5 kink [mrad];tracks" );

    gblkx6HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblkx6", 100, -5, 5 );
    gblkx6HistGBLAlign->setTitle( "GBL kink angle at plane 6;plane 6 kink [mrad];tracks" );

    nmHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "nm", 21, -0.5, 20.5 );
    nmHistGBLAlign->setTitle( "track matches;track matches;events" );

  }//try
  catch( lcio::Exception& e ) {

  }
#endif

}

