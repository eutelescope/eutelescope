// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-

// Author: A.F.Zarnecki, University of Warsaw <mailto:zarnecki@fuw.edu.pl>
// Version: $Id: EUTelTestFitter.cc,v 1.32 2008-10-04 21:08:14 bulgheroni Exp $
// Date 2007.06.04

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifdef USE_GEAR

// eutelescope inlcudes
#include "EUTelTestFitter.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelHistogramManager.h"
#include "EUTelExceptions.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/ITree.h>
#endif

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>


#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>

#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>


using namespace std;
using namespace lcio ;
using namespace marlin ;
using namespace eutelescope;

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelTestFitter::_linChi2HistoName   = "linChi2";
std::string EUTelTestFitter::_logChi2HistoName   = "logChi2";
std::string EUTelTestFitter::_firstChi2HistoName = "firstChi2";
std::string EUTelTestFitter::_bestChi2HistoName  = "bestChi2";
std::string EUTelTestFitter::_fullChi2HistoName  = "fullChi2";
std::string EUTelTestFitter::_nTrackHistoName    = "nTrack";
std::string EUTelTestFitter::_nAllHitHistoName   = "nAllHit";
std::string EUTelTestFitter::_nAccHitHistoName   = "nAccHit";
std::string EUTelTestFitter::_nHitHistoName      = "nHit";
std::string EUTelTestFitter::_nBestHistoName     = "nBest";
std::string EUTelTestFitter::_hitAmbiguityHistoName  = "nAmbig";
#endif


EUTelTestFitter::EUTelTestFitter() : Processor("EUTelTestFitter") {

  // modify processor description
  _description = "Analytical track fitting processor for EUDET telescope" ;


  // register steering parameters:
  //       name, description, class-variable, default value

  // input collection first:

  registerInputCollection( LCIO::TRACKERHIT,
                           "InputCollectionName" ,
                           "Name of the input TrackerHit collection"  ,
                           _inputColName ,
                           std::string("meshit") ) ;

  // output collection

  registerOutputCollection(LCIO::TRACK,"OutputTrackCollectionName",
                           "Collection name for fitted tracks",
                           _outputTrackColName, string ("testfittracks"));

  registerOutputCollection(LCIO::TRACKERHIT,"CorrectedHitCollectionName",
                           "Collection name for corrected particle positions",
                           _correctedHitColName, string ("corrfithits"));

  registerOutputCollection(LCIO::TRACKERHIT,"OutputHitCollectionName",
                           "Collection name for fitted particle hits (positions)",
                           _outputHitColName, string ("testfithits"));

  // compulsory parameters:

  registerProcessorParameter ("InputHitsInTrack",
                              "Flag for storing input (measured) hits in track",
                              _InputHitsInTrack,  static_cast < bool > (false));

  registerProcessorParameter ("OutputHitsInTrack",
                              "Flag for storing output (fitted) hits in track",
                              _OutputHitsInTrack,  static_cast < bool > (true));

  registerProcessorParameter ("DebugEventCount",
                              "Print out every DebugEnevtCount event",
                              _debugCount,  static_cast < int > (10));

  registerProcessorParameter ("AllowMissingHits",
                              "Allowed number of missing hits in the track",
                              _allowMissingHits,  static_cast < int > (0));

  registerProcessorParameter ("AllowSkipHits",
                              "Allowed number of hits removed from the track",
                              _allowSkipHits,  static_cast < int > (0));

  registerProcessorParameter ("MaxPlaneHits",
                              "Maximum number of considered hits per plane",
                              _maxPlaneHits,  static_cast < int > (5));

  registerProcessorParameter ("MissingHitPenalty",
                              "Chi2 penalty for missing hit in the track",
                              _missingHitPenalty,  static_cast < double > (0.));

  registerProcessorParameter ("SkipHitPenalty",
                              "Chi2 penalty for removing hit from the track",
                              _skipHitPenalty,  static_cast < double > (100.));

  registerProcessorParameter ("Chi2Max",
                              "Maximum Chi2 for accepted track fit",
                              _chi2Max,  static_cast < double > (1000.));

  registerProcessorParameter ("UseNominalResolution",
                              "Flag for using nominal resolution instead of position errors",
                              _useNominalResolution,  static_cast < bool > (false));

  registerProcessorParameter ("UseDUT",
                              "Flag for including DUT measurement in the fit",
                              _useDUT,  static_cast < bool > (false));

  registerProcessorParameter ("Ebeam",
                              "Beam energy [GeV]",
                              _eBeam,  static_cast < double > (6.0));

  registerProcessorParameter("HistoInfoFileName",
                             "Name of the histogram information file",
                             _histoInfoFileName, string( "histoinfo.xml" ) );

  // optional parameters

  std::vector<int > initLayerIDs;
  std::vector<float > initLayerShift;

  registerOptionalParameter ("SkipLayerIDs",
                             "Ids of layers which should NOT be included in the fit",
                             _SkipLayerIDs, initLayerIDs);

  registerOptionalParameter ("PassiveLayerIDs",
                             "Ids of layers which should be treated as passive in the fit",
                             _PassiveLayerIDs, initLayerIDs);

  registerOptionalParameter ("AlignLayerIDs",
                             "Ids of layers for which alignment corrections are given",
                             _AlignLayerIDs, initLayerIDs);

  registerOptionalParameter ("AlignLayerShiftX",
                             "Alignment corrections in X for these layers",
                             _AlignLayerShiftX, initLayerShift);

  registerOptionalParameter ("AlignLayerShiftY",
                             "Alignment corrections in Y for these layers",
                             _AlignLayerShiftY, initLayerShift);

  registerOptionalParameter ("AlignLayerRotZ",
                             "Rotation around Z for layer alignment",
                             _AlignLayerRotZ, initLayerShift);

  registerOptionalParameter ("WindowLayerIDs",
                             "Ids of layers for which position window cut are applied",
                             _WindowLayerIDs, initLayerIDs);

  registerOptionalParameter ("WindowMinX",
                             "Lower window edge in X",
                             _WindowMinX, initLayerShift);

  registerOptionalParameter ("WindowMaxX",
                             "Upper window edge in X",
                             _WindowMaxX, initLayerShift);

  registerOptionalParameter ("WindowMinY",
                             "Lower window edge in Y",
                             _WindowMinY, initLayerShift);

  registerOptionalParameter ("WindowMaxY",
                             "Upper window edge in Y",
                             _WindowMaxY, initLayerShift);


  registerOptionalParameter ("MaskLayerIDs",
                             "Ids of layers for which position masks are applied",
                             _MaskLayerIDs, initLayerIDs);

  registerOptionalParameter ("MaskMinX",
                             "Lower mask edge in X",
                             _MaskMinX, initLayerShift);

  registerOptionalParameter ("MaskMaxX",
                             "Upper mask edge in X",
                             _MaskMaxX, initLayerShift);

  registerOptionalParameter ("MaskMinY",
                             "Lower mask edge in Y",
                             _MaskMinY, initLayerShift);

  registerOptionalParameter ("MaskMaxY",
                             "Upper mask edge in Y",
                             _MaskMaxY, initLayerShift);


  registerOptionalParameter ("UseBeamConstraint",
                             "Flag for using beam direction constraint in the fit",
                             _useBeamConstraint,  static_cast < bool > (false));

  registerOptionalParameter ("BeamSpread",
                             "Assumed angular spread of the beam [rad]",
                             _beamSpread,  static_cast < double > (0.1));

  registerOptionalParameter ("SearchMultipleTracks",
                             "Flag for searching multiple tracks in events with multiple hits",
                             _searchMultipleTracks,  static_cast < bool > (false));

  registerOptionalParameter ("AllowAmbiguousHits",
                             "Allow same hit to be used in more than one track",
                             _allowAmbiguousHits, static_cast < bool > (false));

  // initialize all the counters
  _noOfEventWOInputHit   = 0;
  _noOfEventWOTrack = 0;
  _noOfTracks        = 0;

}


void EUTelTestFitter::init() {

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  // check if Marlin was built with GEAR support or not
#ifndef USE_GEAR

  streamlog_out ( ERROR2 ) << "Marlin was not built with GEAR support.\n"
    "You need to install GEAR and recompile Marlin with -DUSE_GEAR before continue." << endl;

  // I'm thinking if this is the case of throwing an exception or
  // not. This is a really error and not something that can
  // exceptionally happens. Still not sure what to do
  exit(-1);

#else

  // check if the GEAR manager pointer is not null!
  if ( Global::GEAR == 0x0 ) {
    streamlog_out ( ERROR2 ) <<  "The GearMgr is not available, for an unknown reason."  << endl;
    exit(-1);
  }

  // Read geometry information from GEAR

  streamlog_out ( MESSAGE0 ) << "Reading telescope geometry description from GEAR " << endl;

  _siPlanesParameters  = const_cast<gear::SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

#endif

// Test output

  if( _SkipLayerIDs.size() ) {
    streamlog_out ( MESSAGE0 )  <<  _SkipLayerIDs.size() << " layers should be skipped" << endl;
  }

  if( _PassiveLayerIDs.size() ) {
    streamlog_out ( MESSAGE0 ) <<  _PassiveLayerIDs.size() << " layers should be considered passive" << endl;
  }

// Active planes only:
//  _nTelPlanes = _siPlanesParameters->getSiPlanesNumber();

// Take all layers defined in GEAR geometry
  _nTelPlanes = _siPlanesLayerLayout->getNLayers();

// Check for DUT

  if( _siPlanesParameters->getSiPlanesType()==_siPlanesParameters->TelescopeWithDUT )    {
    _iDUT = _nTelPlanes ;
    _nTelPlanes++;
  }  else {
    _iDUT = -1 ;
  }

// Read position in Z (for sorting), skip layers if requested

  _planeSort = new int[_nTelPlanes];
  _planePosition   = new double[_nTelPlanes];

  int nSkip=0;

  for(int ipl=0; ipl <  _siPlanesLayerLayout->getNLayers(); ipl++)    {
    _planePosition[ipl-nSkip]=_siPlanesLayerLayout->getLayerPositionZ(ipl);
    _planeSort[ipl-nSkip]=ipl;

// Check if not on "skip list"

    int _plID = _siPlanesLayerLayout->getID(ipl);

    for(int spl=0; spl< (int)_SkipLayerIDs.size() ; spl++) {
      if ( _SkipLayerIDs.at(spl) == _plID) {
        streamlog_out ( MESSAGE0 ) <<  "Skipping layer ID " << _plID
                                   << " at Z = " << _siPlanesLayerLayout->getLayerPositionZ(ipl) << endl;
        nSkip++;
        break;
      }
    }
  }


  _nTelPlanes-=nSkip;

  if(_iDUT>0)  {
    _planePosition[_iDUT-nSkip]=_siPlanesLayerLayout->getDUTPositionZ();
    _planeSort[_iDUT-nSkip]=_iDUT;
  }

  // Binary sorting

  bool sorted;
  do {
    sorted=false;
    for(int iz=0; iz<_nTelPlanes-1 ; iz++) {
      if(_planePosition[iz]>_planePosition[iz+1])   {
        double _posZ = _planePosition[iz];
        _planePosition[iz] = _planePosition[iz+1];
        _planePosition[iz+1] = _posZ;

        int _idZ = _planeSort[iz];
        _planeSort[iz] = _planeSort[iz+1];
        _planeSort[iz+1] = _idZ;

        sorted=true;
      }
    }
  } while(sorted);

// Book local geometry arrays

  _planeID         = new int[_nTelPlanes];
  _planeShiftX     = new double[_nTelPlanes];
  _planeShiftY     = new double[_nTelPlanes];
  _planeRotZ       = new double[_nTelPlanes];
  _planeThickness  = new double[_nTelPlanes];
  _planeX0         = new double[_nTelPlanes];
  _planeResolution = new double[_nTelPlanes];
  _isActive        = new bool[_nTelPlanes];
  _planeWindowIDs  = new std::vector<int>[_nTelPlanes];
  _planeMaskIDs    = new std::vector<int>[_nTelPlanes];

  _nActivePlanes = 0 ;

// Fill arrays with parameters of layer, sorted in Z

  for(int iz=0; iz < _nTelPlanes ; iz++) {
    int ipl=_planeSort[iz];

    int iActive;
    double resolution;

// All dimensions are assumed to be in mm !!!

    if(ipl != _iDUT )    {
      _planeID[iz]=_siPlanesLayerLayout->getID(ipl);
      _planeThickness[iz]=_siPlanesLayerLayout->getLayerThickness(ipl);
      _planeX0[iz]=_siPlanesLayerLayout->getLayerRadLength(ipl);
      resolution = _siPlanesLayerLayout->getSensitiveResolution(ipl);
    }  else {
      _planeID[iz]=_siPlanesLayerLayout->getDUTID();
      _planeThickness[iz]=_siPlanesLayerLayout->getDUTThickness();
      _planeX0[iz]=_siPlanesLayerLayout->getDUTRadLength();
      resolution = _siPlanesLayerLayout->getDUTSensitiveResolution();
    }

    iActive = (resolution > 0);

// Check passive layer list

    for(int ppl=0; ppl< (int)_PassiveLayerIDs.size() ; ppl++) {
      if ( _PassiveLayerIDs.at(ppl) == _planeID[iz])    {
        streamlog_out ( MESSAGE0 ) <<  "Force passive layer ID " << _planeID[iz]
                                   << " at Z = " << _planePosition[iz] << endl ;
        iActive = false;
        break;
      }
    }

    if(iActive && (ipl != _iDUT || _useDUT ))   {
      _isActive[iz] = true ;
      _planeResolution[iz]=resolution;
      _nActivePlanes++ ;
    } else  {
      _isActive[iz] = false ;
      _planeResolution[iz]=0.;
    }

// No alignment corrections in GEAR file
// Look in input options

    _planeShiftX[iz]=0.;
    _planeShiftY[iz]=0.;
    _planeRotZ[iz]=0.;

    for(int apl=0; apl< (int)_AlignLayerIDs.size() ; apl++) {
      if ( _AlignLayerIDs.at(apl) == _planeID[iz]) {
        _planeShiftX[iz]=_AlignLayerShiftX.at(apl);
        _planeShiftY[iz]=_AlignLayerShiftY.at(apl);
        // Rotation can be skipped: check size
        if(apl < (int)_AlignLayerRotZ.size()) {
          _planeRotZ[iz]=_AlignLayerRotZ.at(apl);
        }
        break;
      }
    }

// Check, if there are additional cuts defined for this plane


    for(int wpl=0; wpl< (int)_WindowLayerIDs.size() ; wpl++) {
      if ( _WindowLayerIDs.at(wpl) == _planeID[iz]) {
        _planeWindowIDs[iz].push_back(wpl);
      }
    }
    for(int mpl=0; mpl< (int)_MaskLayerIDs.size() ; mpl++) {
      if ( _MaskLayerIDs.at(mpl) == _planeID[iz]) {
        _planeMaskIDs[iz].push_back(mpl);
      }
    }
    // End of plane loop
  }

  // Get new DUT position (after sorting)

  for(int iz=0;iz< _nTelPlanes ; iz++) {
    if(_planeSort[iz]==_iDUT)  {
      _iDUT=iz;
      break;
    }
  }


  // Print out geometry information
  streamlog_out ( MESSAGE2 ) <<  "Telescope configuration with " << _nTelPlanes << " planes" << endl;


  for(int ipl=0; ipl < _nTelPlanes; ipl++) {
    stringstream ss ;
    if(ipl == _iDUT) {
      ss << "D.U.T.  plane" ;
    } else {
      if(_isActive[ipl]) {
        ss << "Active  plane" ;
      } else {
        ss << "Passive plane" ;
      }
      ss << "  ID = " << _planeID[ipl]
         << "  at Z [mm] = " << _planePosition[ipl]
         << " dZ [um] = " << _planeThickness[ipl]*1000. ;

      if(_isActive[ipl])  {
        ss << "  Res [um] = " << _planeResolution[ipl]*1000. ;

        if(_planeShiftX[ipl] !=0. || _planeShiftY[ipl] !=0. ) {
          ss << "\n  alignment corrections:"
             <<  " dX [mm] = " << _planeShiftX[ipl]
             << " dY [mm] = " << _planeShiftY[ipl] ;
        }
        if(_planeRotZ[ipl] !=0.) {
          ss << " RotZ [rad] = " << _planeRotZ[ipl] ;
        }

        for(int wpl=0; wpl< (int)_planeWindowIDs[ipl].size() ; wpl++)  {
          int iwin= _planeWindowIDs[ipl].at(wpl);
          ss << "\n accepted window: X = " <<  _WindowMinX.at(iwin)
             << " to " <<  _WindowMaxX.at(iwin)
             << " Y = " <<  _WindowMinY.at(iwin)
             << " to " <<  _WindowMaxY.at(iwin);
        }

        for(int mpl=0; mpl< (int)_planeMaskIDs[ipl].size() ; mpl++)  {
          int imsk= _planeMaskIDs[ipl].at(mpl);
          ss << "\n imposed mask: X = " <<  _MaskMinX.at(imsk)
             << " to " <<  _MaskMaxX.at(imsk)
             << " Y = " <<  _MaskMinY.at(imsk)
             << " to " <<  _MaskMaxY.at(imsk);
        }
      }

      streamlog_out( MESSAGE2 ) <<  ss.str() << endl;
    }
  }
  streamlog_out( MESSAGE2 ) << "Total of " << _nActivePlanes << " active sensor planes " << endl;

  // Allocate arrays for track fitting

  _planeHits   = new int[_nTelPlanes];
  _planeChoice = new int[_nTelPlanes];
  _planeMod    = new int[_nTelPlanes];

  _planeX  = new double[_nTelPlanes];
  _planeEx = new double[_nTelPlanes];
  _planeY  = new double[_nTelPlanes];
  _planeEy = new double[_nTelPlanes];

  _planeDist = new double[_nTelPlanes];
  _planeScat = new double[_nTelPlanes];

  _fitX  = new double[_nTelPlanes];
  _fitEx = new double[_nTelPlanes];
  _fitY  = new double[_nTelPlanes];
  _fitEy = new double[_nTelPlanes];

  int arrayDim = _nTelPlanes * _nTelPlanes;

  _fitArray = new double[arrayDim];
  _nominalFitArray = new double[arrayDim];

  _nominalError = new double[_nTelPlanes];

  // Fill nominal fit matrices and
  // calculate expected precision of track fitting

  // Planes are ordered in position along the beam line !

  for(int ipl=0; ipl<_nTelPlanes ; ipl++) {
    if(ipl>0) {
      _planeDist[ipl-1]=1./(_planePosition[ipl] - _planePosition[ipl-1]) ;
    }
    _planeScat[ipl]= 0.0136/_eBeam * sqrt(_planeThickness[ipl]/_planeX0[ipl])
      * (1.+0.038*std::log(_planeThickness[ipl]/_planeX0[ipl])) ;

    if(ipl==0 && _useBeamConstraint) {
      _planeScat[ipl]= 1./(_planeScat[ipl]*_planeScat[ipl]+ _beamSpread*_beamSpread) ;
    } else {
      _planeScat[ipl]= 1./(_planeScat[ipl] * _planeScat[ipl]) ;
    }
    _fitX[ipl] =_fitY[ipl] = 0. ;
    _nominalError[ipl]= _planeResolution[ipl];
  }

  // Fit with nominal parameters

  int status = DoAnalFit(_fitX,_nominalError);

  if(status) {
    streamlog_out( ERROR2 ) << "\n Fit with nominal geometry failed !?!" << endl;
  }
  // Store fit matrix

  for(int imx=0; imx<arrayDim; imx++) {
    _nominalFitArray[imx] = _fitArray[imx];
  }

  stringstream ss;
  ss << "Expected position resolutions [um]: ";
  for(int ipl=0; ipl<_nTelPlanes ; ipl++) {
    ss << _nominalError[ipl]*1000. << "  " ;
  }
  streamlog_out ( MESSAGE2 ) << ss.str() << endl;


// Book histograms

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  bookHistos();
#endif

}

void EUTelTestFitter::processRunHeader( LCRunHeader* runHeader) {

  auto_ptr<EUTelRunHeaderImpl> eutelHeader( new EUTelRunHeaderImpl ( runHeader ) );
  eutelHeader->addProcessor( type() );

  _nRun++ ;

  // Decode and print out Run Header information - just a check

  int runNr = runHeader->getRunNumber();

  streamlog_out( MESSAGE2 )  << "Processing run header " << _nRun
                             << ", run nr " << runNr << endl;

  const std::string detectorName = runHeader->getDetectorName();
  const std::string detectorDescription = runHeader->getDescription();
  const std::vector<std::string> * subDets = runHeader->getActiveSubdetectors();

  streamlog_out( MESSAGE0 ) << detectorName << " : " << detectorDescription << endl;

  int nDet = subDets->size();

  streamlog_out( MESSAGE0 ) << nDet << " subdetectors defined :" << endl;
  stringstream ss;
  for(int idet=0;idet<nDet;idet++) {
    streamlog_out( MESSAGE0 )  << idet+1 << " : " << subDets->at(idet) << endl;
  }

}

void EUTelTestFitter::processEvent( LCEvent * event ) {

  bool debug = ( _debugCount>0 && _nEvt%_debugCount == 0);

  if ( _nEvt % 10 == 0 ) {
    streamlog_out( MESSAGE2 ) << "Processing event "
                              << setw(6) << setiosflags(ios::right) << event->getEventNumber() << " in run "
                              << setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber() << setfill(' ')
                              << " (Total = " << setw(10) << _nEvt << ")" << resetiosflags(ios::left) << endl;
  }
  _nEvt ++ ;

  EUTelEventImpl * euEvent = static_cast<EUTelEventImpl*> ( event );
  if ( euEvent->getEventType() == kEORE ) {
    streamlog_out ( DEBUG ) <<  "EORE found: nothing else to do." << endl;
    return;
  }


  LCCollection* col;
  try {
    col = event->getCollection( _inputColName ) ;
  } catch (lcio::DataNotAvailableException& e) {
    streamlog_out ( DEBUG ) << "Not able to get collection "
                            << _inputColName
                            << "\nfrom event " << event->getEventNumber()
                            << " in run " << event->getRunNumber()  << endl;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nTrackHistoName]))->fill(0);
#endif
    ++_noOfEventWOInputHit;
    return;
  }


  // Copy hits to local table
  // Assign hits to sensor planes
  // =============================


  int nHit = col->getNumberOfElements()  ;

  if ( debug ) {
    streamlog_out( TESTFITTERMESSAGE )  << "Total of " << nHit << " tracker hits in input collection " << endl;
  }

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nAllHitHistoName]))->fill(nHit);
#endif

  if(nHit + _allowMissingHits < _nActivePlanes) {
    streamlog_out ( TESTFITTERMESSAGE ) << "Not enough hits to perform the fit, exiting... " << endl;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nTrackHistoName]))->fill(0);
#endif
    ++_noOfEventWOTrack;
    return;
  }



  double * hitX  = new double[nHit];
  double * hitEx = new double[nHit];
  double * hitY  = new double[nHit];
  double * hitEy = new double[nHit];
  double * hitZ  = new double[nHit];
  int    * hitPlane = new int[nHit];
  int    * hitFits = new int[nHit];

  IntVec * planeHitID   = new IntVec[_nTelPlanes];

  double * bestX  = new double[_nTelPlanes];
  double * bestEx = new double[_nTelPlanes];
  double * bestY  = new double[_nTelPlanes];
  double * bestEy = new double[_nTelPlanes];

  // Loop over hits

  int nGoodHit = 0;

  for(int ihit=0; ihit< nHit ; ihit++) {
    TrackerHit * meshit = dynamic_cast<TrackerHit*>( col->getElementAt(ihit) ) ;

    // Hit position

    const double * pos = meshit->getPosition();

    hitZ[ihit] = pos[2];

    hitFits[ihit]=-1;

    // We have to find Plane ID of the hit
    // by looking at the Z position

    double distMin = 1.;
    hitPlane[ihit] = -1 ;

    for(int ipl=0;ipl<_nTelPlanes;ipl++)  {
      double dist =  hitZ[ihit] - _planePosition[ipl] ;

      if(dist*dist < distMin*distMin)  {
        hitPlane[ihit]=ipl;
        distMin=dist;
      }
    }

    // Ignore hits not matched to any plane

    if(hitPlane[ihit]<0) {

      streamlog_out( WARNING0 )  << "Reconstructed hit outside sensor planes z [mm] = "  << hitZ[ihit] << endl;
      continue;
    }

    // Ignore hit, if plane not declared as active (i.e. not used in the fit)

    if(! _isActive[hitPlane[ihit]]) {
      continue ;
    }

    // Ignore hit also, if maximum number of hits already matched to this plane

    if(_maxPlaneHits>0 && _maxPlaneHits-planeHitID[hitPlane[ihit]].size()<=0) {
      continue ;
    }

    // Hit will be used: correct X and Y position for plane alignment

    //      hitX[ihit] = pos[0];
    //      hitY[ihit] = pos[1];

    // Shift sign changed in 1.22 for consistency with alignment
    // processor

    hitX[ihit] = pos[0]*cos(_planeRotZ[hitPlane[ihit]])
      + pos[1]*sin(_planeRotZ[hitPlane[ihit]])
      - _planeShiftX[hitPlane[ihit]];

    hitY[ihit] = pos[1]*cos(_planeRotZ[hitPlane[ihit]])
      - pos[0]*sin(_planeRotZ[hitPlane[ihit]])
      - _planeShiftY[hitPlane[ihit]];

    // Check Window and Mask cuts, if defined

    bool hitcut = false;

    for(int wpl=0; wpl< (int)_planeWindowIDs[hitPlane[ihit]].size() ; wpl++)  {
      int iwin= _planeWindowIDs[hitPlane[ihit]].at(wpl);

      if(hitX[ihit] <  _WindowMinX.at(iwin)) hitcut = true;
      if(hitX[ihit] >  _WindowMaxX.at(iwin)) hitcut = true;
      if(hitY[ihit] <  _WindowMinY.at(iwin)) hitcut = true;
      if(hitY[ihit] >  _WindowMaxY.at(iwin)) hitcut = true;
    }

    if(hitcut) continue;


    for(int mpl=0; mpl< (int)_planeMaskIDs[hitPlane[ihit]].size() ; mpl++) {
      int imsk= _planeMaskIDs[hitPlane[ihit]].at(mpl);

      if(hitX[ihit] >  _MaskMinX.at(imsk) &&  hitX[ihit] <  _MaskMaxX.at(imsk) &&
         hitY[ihit] >  _MaskMinY.at(imsk) &&  hitY[ihit] <  _MaskMaxY.at(imsk)) {
        hitcut = true;
      }
    }

    if(hitcut)continue;


    // Add hit to hit list for given plane - to be used in track selection

    planeHitID[hitPlane[ihit]].push_back(ihit);
    nGoodHit++;

    // Position uncertainty. Use nominal resolution if not properly defined

    const EVENT::FloatVec cov = meshit->getCovMatrix();

    if(cov.at(0)>0.) {
      hitEx[ihit]=sqrt(cov.at(0));
    } else {
      hitEx[ihit]=_planeResolution[hitPlane[ihit]];
    }
    if(cov.at(2)>0.) {
      hitEy[ihit]=sqrt(cov.at(2));
    } else {
      hitEy[ihit]=_planeResolution[hitPlane[ihit]];
    }

    hitFits[ihit]=0;

    if ( debug ) {
      streamlog_out ( TESTFITTERMESSAGE ) << "Hit " << ihit
                                          << "   X = " << hitX[ihit] << " +/- " << hitEx[ihit]
                                          << "   Y = " << hitY[ihit] << " +/- " << hitEy[ihit]
                                          << "   Z = " << hitZ[ihit] << " (plane " << hitPlane[ihit] << ")" << endl;
    }
  }

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nAccHitHistoName]))->fill(nGoodHit);
#endif

  // Main analysis loop: finding multiple tracks (if allowed)

  // Define output track and hit collections
  LCCollectionVec     * fittrackvec = new LCCollectionVec(LCIO::TRACK);
  LCCollectionVec     * fitpointvec = new LCCollectionVec(LCIO::TRACKERHIT);
  LCCollectionVec     * corrpointvec = new LCCollectionVec(LCIO::TRACKERHIT);

  // Set flag for storing track hits in track collection

  LCFlagImpl flag(fittrackvec->getFlag());
  flag.setBit( LCIO::TRBIT_HITS );
  fittrackvec->setFlag(flag.getFlag());


  int nFittedTracks = 0 ;
  bool firstTrack = true ;
  int ibest;

  // In the current implementation ambiguity mode works only for full
  // tracks, i.e. when _allowMissingHits == 0

  bool ambiguityMode =  _allowAmbiguousHits && ( _allowMissingHits == 0 );

  do {

    // Count planes active in this event and number of fit possibilities

    int nFiredPlanes = 0;
    int nChoice = 1;
    ibest=-1;

    // Count from the last plane, to allow for "smart" track finding

    for(int ipl=_nTelPlanes-1; ipl>=0 ;ipl--)   {
      _planeHits[ipl] = planeHitID[ipl].size() ;

      if(_planeHits[ipl]>0)  {
        nFiredPlanes++;

        _planeChoice[ipl]=_planeHits[ipl]+1;
      }  else {
        _planeChoice[ipl]=1;
      }

      _planeMod[ipl]=nChoice;
      nChoice*=_planeChoice[ipl];

    }

    // Debug output

    if(firstTrack && debug ) {
      for(int ipl=0;ipl<_nTelPlanes;ipl++) {
        if( _isActive[ipl] )  {
          stringstream ss;
          ss << "Plane " << ipl << "  " << _planeHits[ipl] << " hit(s), hit IDs :";

          for( int ihit=0; ihit < (int) planeHitID[ipl].size() ; ihit ++) {
            ss << planeHitID[ipl].at(ihit) << " " ;
          }
          streamlog_out ( DEBUG )  << ss.str() << endl;
        }
      }
    }


    // Check if fit can be done

    if(nFiredPlanes + _allowMissingHits < _nActivePlanes) {
      if( firstTrack ) {
        if(debug) {
          streamlog_out ( TESTFITTERMESSAGE ) <<  "Not enough planes hit to perform the fit " << endl;
        }

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nTrackHistoName]))->fill(0);
#endif

        // before returning clean up the memory
        delete [] bestEy;
        delete [] bestY;
        delete [] bestEx;
        delete [] bestX;
        delete [] planeHitID;
        delete [] hitPlane;
        delete [] hitFits;
        delete [] hitZ;
        delete [] hitEy;
        delete [] hitY;
        delete [] hitEx;
        delete [] hitX;
        
        // increment the counter
        ++_noOfEventWOTrack;
        return ;
      }
    }


    if( firstTrack && debug ) {
      streamlog_out ( TESTFITTERMESSAGE ) << nFiredPlanes << " active sensor planes hit, checking "
                                          << nChoice << " fit possibilities "  << endl;
    }

    // Check all track possibilities

    double chi2min  = numeric_limits<double >::max();
    double chi2best = _chi2Max;
    double bestPenalty = 0;
    int nBestFired = 0;

    // Loop over fit possibilities
    // Start from one-hit track to allow for "smart" skipping of wrong matches

    int istart=0;
    int nmiss=_allowMissingHits;

    while(nmiss>0 || !_isActive[istart])  {
      if(_isActive[istart]) {
        nmiss--;
      }
      istart++;
    }

    for(int ichoice=nChoice-_planeMod[istart]-1; ichoice >=0 ; ichoice--)  {
      int nChoiceFired=0;
      double choiceChi2=-1.;
      double trackChi2=-1.;
      int ifirst=-1;
      int ilast=0;
      int nleft=0;

      // Fill position and error arrays for this hit configuration

      for(int ipl=0;ipl<_nTelPlanes;ipl++)  {
        _planeX[ipl]=_planeY[ipl]=_planeEx[ipl]=_planeEy[ipl]=0.;

        if(_isActive[ipl])  {
          int ihit=(ichoice/_planeMod[ipl])%_planeChoice[ipl];

          if(ihit<_planeHits[ipl])    {
            int jhit = planeHitID[ipl].at(ihit);
            _planeX[ipl]=hitX[jhit];
            _planeY[ipl]=hitY[jhit];
            _planeEx[ipl]=(_useNominalResolution)?_planeResolution[ipl]:hitEx[jhit];
            _planeEy[ipl]=(_useNominalResolution)?_planeResolution[ipl]:hitEy[jhit];
            if(ifirst<0) {
              ifirst=ipl;
            }
            ilast=ipl;
            nleft=0;
            nChoiceFired++;
          }  else {
            nleft++;
          }
        }
      }

      // Check number of selected hits
      // =============================

      // No fit to 1 hit :-)

      if(nChoiceFired < 2) continue;

      // Fit with 2 hits make sense only with beam constraint, or
      // when 2 point fit is allowed

      if(nChoiceFired==2 && !_useBeamConstraint
         && nChoiceFired + _allowMissingHits < _nActivePlanes) continue;

      // Skip also if the fit can not be extended to proper number
      // of planes; no need to check remaining planes !!!
      if(nChoiceFired + nleft < _nActivePlanes - _allowMissingHits ) {
        ichoice-=_planeMod[ilast]-1;
        continue;
      }


      // Select fit method
      // "Nominal" fit only if all active planes used

      if(_useNominalResolution && (nChoiceFired == _nActivePlanes)) {
        choiceChi2 = NominalFit();
      } else {
        if(_useNominalResolution)    choiceChi2 = SingleFit();
        else choiceChi2 = MatrixFit();
      }

      // Fit failed ?

      if(choiceChi2 < 0.)   {
        streamlog_out ( WARNING2 ) << "Fit to " << nChoiceFired
                                   << " planes failed for event " << event->getEventNumber()
                                   << " in run " << event->getRunNumber()  << endl;

        continue ;
      }

      // Penalty for missing or skiped hits

      double penalty = (_nActivePlanes-nFiredPlanes)*_missingHitPenalty
        +   (nFiredPlanes-nChoiceFired)*_skipHitPenalty ;


      trackChi2 = choiceChi2+penalty;

      if(nChoiceFired + _allowMissingHits >= _nActivePlanes &&
         nChoiceFired + _allowSkipHits    >= nFiredPlanes  &&
         trackChi2 < chi2min) {
        chi2min=trackChi2;
      }

      // Check if better than best fit (or chi2Max if hit ambiguity allowed)
      // If not: skip also all track possibilities which include
      // this hit selection !!!
      if((choiceChi2 >= _chi2Max) ||
         (choiceChi2 >=  chi2best  && !ambiguityMode) ) {
        ichoice-=_planeMod[ilast]-1;
        continue;
      }

      //
      // Skip fit if could not be accepted (too few planes fired)
      //
      if(nChoiceFired + _allowMissingHits < _nActivePlanes ||
         nChoiceFired + _allowSkipHits    < nFiredPlanes ) {
        continue;
      }

      // Best fit ?

      if(trackChi2<chi2best)   {
        chi2best=trackChi2;
        bestPenalty=penalty;
        ibest=ichoice;
        nBestFired=nChoiceFired;
        for(int ipl=0;ipl<_nTelPlanes;ipl++)  {
          bestX[ipl]=_fitX[ipl];
          bestY[ipl]=_fitY[ipl];
          bestEx[ipl]=_fitEx[ipl];
          bestEy[ipl]=_fitEy[ipl];
        }
      }

      // If hit ambiguity allowed (same hit can be used in many
      // tracks) - fill all tracks passing chi2 cut

      if(trackChi2<_chi2Max   && ambiguityMode ) {

        if(debug)  {
          streamlog_out ( TESTFITTERMESSAGE) <<  "Track reconstructed from " << nChoiceFired << " hits: " << endl;


          // print out hits contributing to the fit


          for(int ipl=0;ipl<_nTelPlanes;ipl++)  {
            if(_isActive[ipl]) {
              int ihit=(ichoice/_planeMod[ipl])%_planeChoice[ipl];

              if(ihit<_planeHits[ipl])  {
                int jhit = planeHitID[ipl].at(ihit);
                streamlog_out ( TESTFITTERMESSAGE) << "Hit " << jhit
                                                   << "   X = " << hitX[jhit]
                                                   << "   Y = " << hitY[jhit]
                                                   << "   Z = " << hitZ[jhit] << " (plane" << hitPlane[jhit] << ")"  << endl;
              }
            }
          }

          streamlog_out ( TESTFITTERMESSAGE) << " Fitted positions in telescope planes:" << endl;

          for(int ipl=0;ipl<_nTelPlanes;ipl++) {
            streamlog_out ( TESTFITTERMESSAGE) << "  X = " << _fitX[ipl] << " +/- " << _fitEx[ipl]
                                               << "  Y = " << _fitY[ipl] << " +/- " << _fitEy[ipl]
                                               << "  at Z = " << _planePosition[ipl]  << endl;
          }


          streamlog_out ( TESTFITTERMESSAGE) << " Fit chi2 = " << trackChi2 << " including penalties of " << penalty << endl;

        }

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
        // Fill Chi2 histograms

        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_linChi2HistoName]))->fill(trackChi2);

        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_logChi2HistoName]))->fill(log10(trackChi2));

        if(_allowMissingHits && nBestFired==_nActivePlanes) {
          (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_fullChi2HistoName]))->fill(log10(trackChi2));
        }

        // Fill hit histograms

        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nHitHistoName]))->fill(nChoiceFired);
#endif

        // Write fit result out

        TrackImpl * fittrack = new TrackImpl();

        // Following parameters are not used for Telescope
        // and are set to zero (just in case)
        fittrack->setOmega(0.);     // curvature of the track
        fittrack->setD0(0.);        // impact paramter of the track in (r-phi)
        fittrack->setZ0(0.);        // impact paramter of the track in (r-z)
        fittrack->setPhi(0.);       // phi of the track at reference point
        fittrack->setTanLambda(0.); // dip angle of the track at reference point

        // Used class members

        fittrack->setChi2(trackChi2);  // Chi2 of the fit (including penalties)
        fittrack->setNdf(nChoiceFired); // Number of planes fired (!)

        // Fitted position at DUT is stored as Reference Point
        // (see below in loop over telescope planes)
        // If no DUT present: use position in the first plane !

        fittrack->setIsReferencePointPCA(false);
        float refpoint[3];

        // Track points fitted in each plane are stored as track hits

        for(int ipl=0;ipl<_nTelPlanes;ipl++)  {
          TrackerHitImpl * fitpoint = new TrackerHitImpl;

          // Hit type is set to 32, to distinguish from measured
          // hits (hit type = cluster type = 0 ... 31)

          fitpoint->setType(32);

          // fitted position in a plane

          double pos[3];

          pos[0]=_fitX[ipl];
          pos[1]=_fitY[ipl];
          pos[2]=_planePosition[ipl];

          fitpoint->setPosition(pos);

          // Covariance matrix of the position
          // (stored as lower triangle matrix, i.e.  cov(xx),cov(y,x),cov(y,y) ).

          float cov[TRKHITNCOVMATRIX];

          cov[0]=_fitEx[ipl]*_fitEx[ipl];
          cov[1]=0.;
          cov[2]=_fitEy[ipl]*_fitEy[ipl];
          cov[3]=cov[4]=0.;
          cov[5]=_planeThickness[ipl]*_planeThickness[ipl]/12.;

          fitpoint->setCovMatrix(cov);

          // store fit point

          fitpointvec->push_back(fitpoint);

          //   add fitted point to track

          if(_OutputHitsInTrack)   fittrack->addHit(fitpoint);


          //                 << "   X = " << hitX[ihit] << " +/- " << hitEx[ihit]
          //       << "   Y = " << hitY[ihit] << " +/- " <<
          //       hitEy[ihit]

          // add measured point to track (if found)

          if(_InputHitsInTrack && _isActive[ipl])  {
            int ihit=(ichoice/_planeMod[ipl])%_planeChoice[ipl];

            if(ihit<_planeHits[ipl])   {
              int jhit = planeHitID[ipl].at(ihit);
              TrackerHitImpl * meshit = dynamic_cast<TrackerHitImpl*>( col->getElementAt(jhit) ) ;
              TrackerHitImpl * corrhit = new TrackerHitImpl;
              //
              // Copy input hit data
              //
              corrhit->setType(meshit->getType());
              corrhit->setTime(meshit->getTime());
              corrhit->setdEdx(meshit->getdEdx());
              corrhit->rawHits()=meshit->getRawHits();
              //
              // Use corrected position
              //
              pos[0]=hitX[jhit];
              pos[1]=hitY[jhit];
              pos[2]=_planePosition[ipl];

              corrhit->setPosition(pos);
              //
              // Include errors, as used in the fit
              //
              cov[0]=hitEx[jhit]*hitEx[jhit];
              cov[1]=0.;
              cov[2]=hitEy[jhit]*hitEy[jhit];
              cov[3]=cov[4]=0.;
              cov[5]=_planeThickness[ipl]*_planeThickness[ipl]/12.;

              corrhit->setCovMatrix(cov);

              corrpointvec->push_back(corrhit);

              fittrack->addHit(corrhit);
            }
          }

          // Use position at DUT as a track reference point.
          // If no DUT present: use position in the first plane !

          if(ipl==_iDUT || (_iDUT<0 && ipl==0)) {
            for(int iref=0;iref<3;iref++) {
              refpoint[iref]=pos[iref];
            }
          }
        }

        // Store track reference point.
        fittrack->setReferencePoint(refpoint);

        fittrackvec->addElement(fittrack);

        // Count number of tracks for each hit

        for(int ipl=0;ipl<_nTelPlanes;ipl++) {
          if(_isActive[ipl]) {
            int ihit=(ichoice/_planeMod[ipl])%_planeChoice[ipl];

            if(ihit<_planeHits[ipl]) {
              hitFits[planeHitID[ipl].at(ihit)]++;
            }
          }
        }

        nFittedTracks++;

      }  // end of track filling - if(choiceChi2+penalty<_chi2Max   && ambiguityMode )



    }
    // End of loop over track possibilities

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    if(firstTrack)
      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_firstChi2HistoName]))->fill(log10(chi2min));
#endif

    if(ibest<0 && firstTrack) {
      if(debug) {
        streamlog_out ( TESTFITTERMESSAGE ) << "No track fulfilling search criteria found ! " << endl;
      }

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nTrackHistoName]))->fill(0);
#endif

      // before throwing the exception I should clean up the
      // memory...
      delete [] bestEy;
      delete [] bestY;
      delete [] bestEx;
      delete [] bestX;
      delete [] planeHitID;
      delete [] hitPlane;
      delete [] hitFits;
      delete [] hitZ;
      delete [] hitEy;
      delete [] hitY;
      delete [] hitEx;
      delete [] hitX;
      return;
    }


    if(ibest>=0 && firstTrack && ambiguityMode )  {

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
      if(_searchMultipleTracks) {
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_bestChi2HistoName]))->fill(log10(chi2best));
      }
      if(_searchMultipleTracks) {
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nBestHistoName]))->fill(nBestFired);
      }
#endif

    }


    if(ibest>=0 && !ambiguityMode ) {

      if(debug)  {
        streamlog_out ( TESTFITTERMESSAGE)  << "Track reconstructed from " << nBestFired << " hits: " << endl;


        // print out hits contributing to the fit


        for(int ipl=0;ipl<_nTelPlanes;ipl++)  {
          if(_isActive[ipl])   {
            int ihit=(ibest/_planeMod[ipl])%_planeChoice[ipl];

            if(ihit<_planeHits[ipl])  {
              int jhit = planeHitID[ipl].at(ihit);
              streamlog_out ( TESTFITTERMESSAGE) << "Hit " << jhit
                                                 << "   X = " << hitX[jhit]
                                                 << "   Y = " << hitY[jhit]
                                                 << "   Z = " << hitZ[jhit] << " (plane" << hitPlane[jhit] << ")" << endl;
            }
          }
        }


        streamlog_out ( TESTFITTERMESSAGE) << " Fitted positions in telescope planes:" << endl;

        for(int ipl=0;ipl<_nTelPlanes;ipl++) {
          streamlog_out ( TESTFITTERMESSAGE) << "  X = " << bestX[ipl] << " +/- " << bestEx[ipl]
                                             << "  Y = " << bestY[ipl] << " +/- " << bestEy[ipl]
                                             << "  at Z = " << _planePosition[ipl] << endl;
        }

        streamlog_out ( TESTFITTERMESSAGE) << " Fit chi2 = " << chi2best << " including penalties of " << bestPenalty << endl;

      }

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
      // Fill Chi2 histograms

      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_linChi2HistoName]))->fill(chi2best);

      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_logChi2HistoName]))->fill(log10(chi2best));

      if(_searchMultipleTracks && firstTrack) {
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_bestChi2HistoName]))->fill(log10(chi2best));
      }

      if(_allowMissingHits && nBestFired==_nActivePlanes) {
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_fullChi2HistoName]))->fill(log10(chi2best));
      }

      // Fill hit histograms

      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nHitHistoName]))->fill(nBestFired);

      if(_searchMultipleTracks && firstTrack) {
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nBestHistoName]))->fill(nBestFired);
      }
#endif


      // Write fit result out

      TrackImpl * fittrack = new TrackImpl();

      // Following parameters are not used for Telescope
      // and are set to zero (just in case)
      fittrack->setOmega(0.);     // curvature of the track
      fittrack->setD0(0.);        // impact paramter of the track in (r-phi)
      fittrack->setZ0(0.);        // impact paramter of the track in (r-z)
      fittrack->setPhi(0.);       // phi of the track at reference point
      fittrack->setTanLambda(0.); // dip angle of the track at reference point

      // Used class members

      fittrack->setChi2(chi2best);  // Chi2 of the fit (including penalties)
      fittrack->setNdf(nBestFired); // Number of planes fired (!)

      // Fitted position at DUT is stored as Reference Point
      // (see below in loop over telescope planes)
      // If no DUT present: use position in the first plane !

      fittrack->setIsReferencePointPCA(false);
      float refpoint[3];

      // Track points fitted in each plane are stored as track hits

      for(int ipl=0;ipl<_nTelPlanes;ipl++)  {
        TrackerHitImpl * fitpoint = new TrackerHitImpl;

        // Hit type is set to 32, to distinguish from measured
        // hits (hit type = cluster type = 0 ... 31)

        fitpoint->setType(32);

        // fitted position in a plane

        double pos[3];

        pos[0]=bestX[ipl];
        pos[1]=bestY[ipl];
        pos[2]=_planePosition[ipl];

        fitpoint->setPosition(pos);

        // Covariance matrix of the position
        // (stored as lower triangle matrix, i.e.  cov(xx),cov(y,x),cov(y,y) ).

        float cov[TRKHITNCOVMATRIX];

        cov[0]=bestEx[ipl]*bestEx[ipl];
        cov[1]=0.;
        cov[2]=bestEy[ipl]*bestEy[ipl];
        cov[3]=cov[4]=0.;
        cov[5]=_planeThickness[ipl]*_planeThickness[ipl]/12.;

        fitpoint->setCovMatrix(cov);

        // store fit point

        fitpointvec->push_back(fitpoint);

        //   add fitted point to track

        if(_OutputHitsInTrack)  fittrack->addHit(fitpoint);


        //                 << "   X = " << hitX[ihit] << " +/- " << hitEx[ihit]
        //       << "   Y = " << hitY[ihit] << " +/- " <<
        //       hitEy[ihit]

        // add measured point to track (if found)

        if(_InputHitsInTrack && _isActive[ipl])  {
          int ihit=(ibest/_planeMod[ipl])%_planeChoice[ipl];

          if(ihit<_planeHits[ipl])  {
            int jhit = planeHitID[ipl].at(ihit);
            TrackerHitImpl * meshit = dynamic_cast<TrackerHitImpl*>( col->getElementAt(jhit) ) ;
            TrackerHitImpl * corrhit = new TrackerHitImpl;
            //
            // Copy input hit data
            //
            corrhit->setType(meshit->getType());
            corrhit->setTime(meshit->getTime());
            corrhit->setdEdx(meshit->getdEdx());
            corrhit->rawHits()=meshit->getRawHits();
            //
            // Use corrected position
            //
            pos[0]=hitX[jhit];
            pos[1]=hitY[jhit];
            pos[2]=_planePosition[ipl];

            corrhit->setPosition(pos);
            //
            // Include errors, as used in the fit
            //
            cov[0]=hitEx[jhit]*hitEx[jhit];
            cov[1]=0.;
            cov[2]=hitEy[jhit]*hitEy[jhit];
            cov[3]=cov[4]=0.;
            cov[5]=_planeThickness[ipl]*_planeThickness[ipl]/12.;

            corrhit->setCovMatrix(cov);

            corrpointvec->push_back(corrhit);

            fittrack->addHit(corrhit);
          }
        }

        // Use position at DUT as a track reference point.
        // If no DUT present: use position in the first plane !

        if(ipl==_iDUT || (_iDUT<0 && ipl==0)) {
          for(int iref=0;iref<3;iref++) {
            refpoint[iref]=pos[iref];
          }
        }
      }
      // Store track reference point.
      fittrack->setReferencePoint(refpoint);

      fittrackvec->addElement(fittrack);
      // increment the total track counter
      ++_noOfTracks;
      nFittedTracks++;
    }

    // If multiple tracks allowed: remove hits from fitted track from the list

    if(ibest>=0 && _searchMultipleTracks && !ambiguityMode )  {
      for(int ipl=0;ipl<_nTelPlanes;ipl++) {
        if(_isActive[ipl]) {
          int ihit=(ibest/_planeMod[ipl])%_planeChoice[ipl];
          if(ihit<_planeHits[ipl]){
            planeHitID[ipl].erase(planeHitID[ipl].begin()+ihit);
            nGoodHit--;
          }
        }
      }
    }
    firstTrack=false;
  } while(_searchMultipleTracks && ibest>=0 &&
          nGoodHit + _allowMissingHits >= _nActivePlanes &&  !ambiguityMode );

  // End of track loop

  if(nFittedTracks > 0 ) {
    event->addCollection(fittrackvec,_outputTrackColName);
    event->addCollection(fitpointvec,_outputHitColName);
    event->addCollection(corrpointvec,_correctedHitColName);
  } else {
    delete fittrackvec;
    delete fitpointvec;
    delete corrpointvec;
  }

  if ( fittrackvec->size() == 0 ) {
    ++_noOfEventWOTrack;
  }

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  // Number of reconstructed tracks

  (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nTrackHistoName]))->fill(nFittedTracks);
#endif

  // Hit ambiguity

  if ( ambiguityMode ) {
    for(int ihit=0; ihit< nHit ; ihit++) {
      if(_isActive[hitPlane[ihit]]) {

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_hitAmbiguityHistoName]))->fill(hitFits[ihit]);
#endif
      }
    }

  }


  // Clear all working arrays

  delete [] bestEy;
  delete [] bestY;
  delete [] bestEx;
  delete [] bestX;
  delete [] planeHitID;
  delete [] hitPlane;
  delete [] hitFits;
  delete [] hitZ;
  delete [] hitEy;
  delete [] hitY;
  delete [] hitEx;
  delete [] hitX;


  return;

}



void EUTelTestFitter::check( LCEvent * evt ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void EUTelTestFitter::end(){

  //   std::cout << "EUTelTestFitter::end()  " << name()
  //        << " processed " << _nEvt << " events in " << _nRun << " runs "
  //        << std::endl ;


  // Print the summer
  streamlog_out( MESSAGE ) << "Total number of processed events:    " << setw(10) << setiosflags(ios::right) << _nEvt << resetiosflags(ios::right) << endl
                           << "Total number of event w/o input hit: " << setw(10) << setiosflags(ios::right) << _noOfEventWOInputHit 
                           << resetiosflags(ios::right) << endl
                           << "Total number of event w/o track:     " << setw(10) << setiosflags(ios::right) << _noOfEventWOTrack
                           << resetiosflags(ios::right) << endl
                           << "Total number of reconstructed tracks " << setw(10) << setiosflags(ios::right) << _noOfTracks << resetiosflags(ios::right)
                           << endl;


  // Clean memory

  delete [] _planeSort ;
  delete [] _planeID ;
  delete [] _planeShiftX ;
  delete [] _planeShiftY ;
  delete [] _planeRotZ ;
  delete [] _planePosition ;
  delete [] _planeThickness  ;
  delete [] _planeX0  ;
  delete [] _planeResolution ;
  delete [] _planeDist ;
  delete [] _planeScat ;
  delete [] _isActive ;

  delete [] _planeHits ;
  delete [] _planeChoice ;
  delete [] _planeMod ;

  delete [] _planeX ;
  delete [] _planeEx ;
  delete [] _planeY ;
  delete [] _planeEy ;

  delete [] _fitX  ;
  delete [] _fitEx ;
  delete [] _fitY ;
  delete [] _fitEy ;
  delete [] _fitArray ;

  delete [] _nominalFitArray ;
  delete [] _nominalError ;
}



void EUTelTestFitter::bookHistos()
{

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  streamlog_out ( MESSAGE2 ) <<  "Booking histograms " << endl;


  streamlog_out ( MESSAGE2 ) << "Histogram information searched in " << _histoInfoFileName << endl;

  auto_ptr<EUTelHistogramManager> histoMgr( new EUTelHistogramManager( _histoInfoFileName ));
  EUTelHistogramInfo    * histoInfo;
  bool                    isHistoManagerAvailable;

  try {
    isHistoManagerAvailable = histoMgr->init();
  } catch ( ios::failure& e) {
    streamlog_out( ERROR ) << "I/O problem with " << _histoInfoFileName << "\n"
                           << "Continuing without histogram manager"    << endl;
    isHistoManagerAvailable = false;
  } catch ( ParseException& e ) {
    streamlog_out( ERROR ) << e.what() << "\n"
                           << "Continuing without histogram manager" << endl;
    isHistoManagerAvailable = false;
  }

  // Chi2 distribution for all accepted tracks (linear scale)

  int    chi2NBin  = 1000;
  double chi2Min   = 0.;
  double chi2Max   = 200.;
  string chi2Title = "Chi2 distribution";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_linChi2HistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG )  << (* histoInfo ) << endl;
          chi2NBin = histoInfo->_xBin;
          chi2Min  = histoInfo->_xMin;
          chi2Max  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) chi2Title = histoInfo->_title;
        }
    }


  AIDA::IHistogram1D * linChi2Histo = AIDAProcessor::histogramFactory(this)->createHistogram1D( _linChi2HistoName.c_str(),chi2NBin,chi2Min,chi2Max);

  linChi2Histo->setTitle(chi2Title.c_str());

  _aidaHistoMap.insert(make_pair(_linChi2HistoName, linChi2Histo));


  // log(Chi2) distribution for all accepted tracks

  chi2NBin  = 100;
  chi2Min   = -2.;
  chi2Max   = 8.;
  chi2Title = "log(Chi2) distribution";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_logChi2HistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG )   << (* histoInfo ) << endl;
          chi2NBin = histoInfo->_xBin;
          chi2Min  = histoInfo->_xMin;
          chi2Max  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) chi2Title = histoInfo->_title;
        }
    }


  AIDA::IHistogram1D * logChi2Histo = AIDAProcessor::histogramFactory(this)->createHistogram1D( _logChi2HistoName.c_str(),chi2NBin,chi2Min,chi2Max);

  logChi2Histo->setTitle(chi2Title.c_str());

  _aidaHistoMap.insert(make_pair(_logChi2HistoName, logChi2Histo));


// Additional Chi2 histogram for first track candidate (without chi2 cut)

  string firstchi2Title = chi2Title + ", first candidate (before cut)";

  AIDA::IHistogram1D * firstChi2Histo = AIDAProcessor::histogramFactory(this)->createHistogram1D( _firstChi2HistoName.c_str(),chi2NBin,chi2Min,chi2Max);

  firstChi2Histo->setTitle(firstchi2Title.c_str());


  _aidaHistoMap.insert(make_pair(_firstChi2HistoName, firstChi2Histo));


// Chi2 histogram for best tracks in an event - use same binning

  if(_searchMultipleTracks)
    {
      string bestchi2Title = chi2Title + ", best fit";

      AIDA::IHistogram1D * bestChi2Histo = AIDAProcessor::histogramFactory(this)->createHistogram1D( _bestChi2HistoName.c_str(),chi2NBin,chi2Min,chi2Max);

      bestChi2Histo->setTitle(bestchi2Title.c_str());


      _aidaHistoMap.insert(make_pair(_bestChi2HistoName, bestChi2Histo));
    }


// Another Chi2 histogram for all full tracks in an event - use same binning

  if(_allowMissingHits)
    {
      string fullchi2Title = chi2Title + ", full tracks";

      AIDA::IHistogram1D * fullChi2Histo = AIDAProcessor::histogramFactory(this)->createHistogram1D( _fullChi2HistoName.c_str(),chi2NBin,chi2Min,chi2Max);

      fullChi2Histo->setTitle(fullchi2Title.c_str());


      _aidaHistoMap.insert(make_pair(_fullChi2HistoName, fullChi2Histo));
    }


  // Distribution of number of tracks

  int    trkNBin  = 11;
  double trkMin   = -0.5;
  double trkMax   = 10.5;
  string trkTitle = "Number of reconstructed tracks in an event";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_nTrackHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG ) << (* histoInfo ) << endl;
          trkNBin = histoInfo->_xBin;
          trkMin  = histoInfo->_xMin;
          trkMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) trkTitle = histoInfo->_title;
        }
    }


  AIDA::IHistogram1D * nTrackHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( _nTrackHistoName.c_str(),trkNBin,trkMin,trkMax);

  nTrackHisto->setTitle(trkTitle.c_str());

  _aidaHistoMap.insert(make_pair(_nTrackHistoName, nTrackHisto));



  // Number of hits in input collection

  int    hitNBin  = 1000;
  double hitMin   = 0.;
  double hitMax   = 1000.;
  string hitTitle = "Number of hits in input collection";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_nAllHitHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG ) << (* histoInfo ) << endl;
          hitNBin = histoInfo->_xBin;
          hitMin  = histoInfo->_xMin;
          hitMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) hitTitle = histoInfo->_title;
        }
    }

  AIDA::IHistogram1D * nAllHitHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( _nAllHitHistoName.c_str(),hitNBin,hitMin,hitMax);

  nAllHitHisto->setTitle(hitTitle.c_str());

  _aidaHistoMap.insert(make_pair(_nAllHitHistoName, nAllHitHisto));



  // Number of accepted hits

  hitNBin  = 1000;
  hitMin   = 0.;
  hitMax   = 1000.;
  hitTitle = "Number of hits accepted from input collection";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_nAccHitHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG ) << (* histoInfo ) << endl;
          hitNBin = histoInfo->_xBin;
          hitMin  = histoInfo->_xMin;
          hitMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) hitTitle = histoInfo->_title;
        }
    }

  AIDA::IHistogram1D * nAccHitHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( _nAccHitHistoName.c_str(),hitNBin,hitMin,hitMax);

  nAccHitHisto->setTitle(hitTitle.c_str());

  _aidaHistoMap.insert(make_pair(_nAccHitHistoName, nAccHitHisto));



  // Number of hits per track

  hitNBin  = 11;
  hitMin   = -0.5;
  hitMax   = 10.5;
  hitTitle = "Number of hits used to reconstruct the track";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_nHitHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG ) << (* histoInfo ) << endl;
          hitNBin = histoInfo->_xBin;
          hitMin  = histoInfo->_xMin;
          hitMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) hitTitle = histoInfo->_title;
        }
    }

  AIDA::IHistogram1D * nHitHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( _nHitHistoName.c_str(),hitNBin,hitMin,hitMax);

  nHitHisto->setTitle(hitTitle.c_str());

  _aidaHistoMap.insert(make_pair(_nHitHistoName, nHitHisto));



// Additional hit number histogram for best tracks in an event - use same binning

  if(_searchMultipleTracks)
    {
      string bestTitle = hitTitle + ", best fit";

      AIDA::IHistogram1D * BestHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( _nBestHistoName.c_str(),hitNBin,hitMin,hitMax);

      BestHisto->setTitle(bestTitle.c_str());


      _aidaHistoMap.insert(make_pair(_nBestHistoName, BestHisto));
    }


// Additional histogram for number of tracks fitted to given hit - use same binning

  if(_allowAmbiguousHits && ( _allowMissingHits == 0 ))
    {
      string ambigTitle = "Number of tracks containing given hit (ambiguity)";

      AIDA::IHistogram1D * AmbigHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( _hitAmbiguityHistoName.c_str(),hitNBin,hitMin,hitMax);

      AmbigHisto->setTitle(ambigTitle.c_str());


      _aidaHistoMap.insert(make_pair(_hitAmbiguityHistoName, AmbigHisto));
    }


// List all booked histogram - check of histogram map filling

  streamlog_out ( MESSAGE0 ) <<  _aidaHistoMap.size() << " histograms booked" << endl;


  map<string, AIDA::IBaseHistogram *>::iterator mapIter;
  for(mapIter = _aidaHistoMap.begin(); mapIter != _aidaHistoMap.end() ; mapIter++ )
    streamlog_out ( DEBUG ) <<  mapIter->first << " : " <<  (mapIter->second)->title()  << endl;

  streamlog_out ( DEBUG ) << "Histogram booking completed \n\n" << endl;
#else
  streamlog_out ( MESSAGE4 ) << "No histogram produced because Marlin doesn't use AIDA" << endl;
#endif

  return;
}
//
// ===============================================================================
//
//  Private function members
//


double EUTelTestFitter::MatrixFit()
{
  for(int ipl=0; ipl<_nTelPlanes;ipl++)
    {
      _fitX[ipl]=_planeX[ipl];
      _fitEx[ipl]=_planeEx[ipl];
      _fitY[ipl]=_planeY[ipl];
      _fitEy[ipl]=_planeEy[ipl];
    }

  int status = DoAnalFit(_fitX,_fitEx);

  if(status)return -1. ;

  status = DoAnalFit(_fitY,_fitEy);

  if(status)return -1. ;

  double chi2=GetFitChi2();

  return chi2 ;
}

double EUTelTestFitter::SingleFit()
{

  for(int ipl=0; ipl<_nTelPlanes;ipl++)
    {
      _fitX[ipl]=_planeX[ipl];
      _fitEx[ipl]=_planeEx[ipl];
    }

  int status = DoAnalFit(_fitX,_fitEx);

  if(status)return -1. ;

  // Use same matrix to solve equation in Y

  for(int ipl=0; ipl<_nTelPlanes;ipl++)
    {
      _fitEy[ipl]=_fitEx[ipl];

      _fitY[ipl]=0. ;
      for(int jpl=0; jpl<_nTelPlanes;jpl++)
        if(_planeEy[jpl]>0.)
          _fitY[ipl]+=_fitArray[ipl+jpl*_nTelPlanes]*_planeY[jpl]/_planeEy[jpl]/_planeEy[jpl];
    }

  double chi2=GetFitChi2();

  return chi2 ;
}

double EUTelTestFitter::NominalFit()
{
  for(int ipl=0; ipl<_nTelPlanes;ipl++)
    {
      _fitEy[ipl]=_fitEx[ipl]=_nominalError[ipl];

      _fitX[ipl]=0. ;
      _fitY[ipl]=0. ;
      for(int jpl=0; jpl<_nTelPlanes;jpl++)
        {
          if(_planeEx[jpl]>0.)
            _fitX[ipl]+=_nominalFitArray[ipl+jpl*_nTelPlanes]*_planeX[jpl]/_planeEx[jpl]/_planeEx[jpl];
          if(_planeEy[jpl]>0.)
            _fitY[ipl]+=_nominalFitArray[ipl+jpl*_nTelPlanes]*_planeY[jpl]/_planeEy[jpl]/_planeEy[jpl];
        }
    }

  double chi2=GetFitChi2();

  return chi2 ;
}


int EUTelTestFitter::DoAnalFit(double * pos, double *err)
{
  for(int ipl=0; ipl<_nTelPlanes;ipl++)
    {
      if(_isActive[ipl] && err[ipl]>0)
        err[ipl]=1./err[ipl]/err[ipl] ;
      else
        err[ipl] = 0. ;

      pos[ipl]*=err[ipl];
    }


  for(int ipl=0; ipl<_nTelPlanes;ipl++)
    for(int jpl=0; jpl<_nTelPlanes;jpl++)
      {
        int imx=ipl+jpl*_nTelPlanes;

        _fitArray[imx] = 0.;

        if(jpl==ipl-2)
          _fitArray[imx] += _planeDist[ipl-2]*_planeDist[ipl-1]*_planeScat[ipl-1] ;

        if(jpl==ipl+2)
          _fitArray[imx] += _planeDist[ipl]*_planeDist[ipl+1]*_planeScat[ipl+1] ;

        if(jpl==ipl-1)
          {
            if(ipl>0 &&  ipl < _nTelPlanes-1)
              _fitArray[imx] -= _planeDist[ipl-1]*(_planeDist[ipl]+_planeDist[ipl-1])*_planeScat[ipl] ;
            if(ipl>1)
              _fitArray[imx] -= _planeDist[ipl-1]*(_planeDist[ipl-1]+_planeDist[ipl-2])*_planeScat[ipl-1] ;
          }

        if(jpl==ipl+1)
          {
            if(ipl>0 && ipl < _nTelPlanes-1)
              _fitArray[imx] -= _planeDist[ipl]*(_planeDist[ipl]+_planeDist[ipl-1])*_planeScat[ipl] ;
            if(ipl < _nTelPlanes-2)
              _fitArray[imx] -= _planeDist[ipl]*(_planeDist[ipl+1]+_planeDist[ipl])*_planeScat[ipl+1] ;
          }

        if(jpl==ipl)
          {
            _fitArray[imx] += err[ipl] ;

            if(ipl>0 && ipl<_nTelPlanes-1)
              _fitArray[imx] += _planeScat[ipl]*(_planeDist[ipl]+_planeDist[ipl-1])*(_planeDist[ipl]+_planeDist[ipl-1]) ;

            if(ipl > 1 )
              _fitArray[imx] += _planeScat[ipl-1]*_planeDist[ipl-1]*_planeDist[ipl-1] ;

            if(ipl < _nTelPlanes-2)
              _fitArray[imx] += _planeScat[ipl+1]*_planeDist[ipl]*_planeDist[ipl] ;
          }

        // For beam constraint

        if(ipl==jpl && ipl<2 && _useBeamConstraint)
          _fitArray[imx] += _planeScat[0]*_planeDist[0]*_planeDist[0] ;

        if(ipl+jpl==1 && _useBeamConstraint)
          _fitArray[imx] -= _planeScat[0]*_planeDist[0]*_planeDist[0] ;


      }

  int status=GaussjSolve(_fitArray,pos,_nTelPlanes) ;

  if(status)
    {
      cerr << "Singular matrix in track fitting algorithm ! " << endl;
      for(int ipl=0;ipl<_nTelPlanes;ipl++)
        err[ipl]=0. ;
    }
  else
    for(int ipl=0;ipl<_nTelPlanes;ipl++)
      err[ipl]=sqrt(_fitArray[ipl+ipl*_nTelPlanes]);

  return status ;
}



double EUTelTestFitter::GetFitChi2()
{
  double chi2=0. ;

  // Measurements

  for(int ipl=0; ipl<_nTelPlanes;ipl++)
    if(_isActive[ipl])
      {
        if(_planeEx[ipl]>0.)
          chi2+=(_fitX[ipl]-_planeX[ipl])*(_fitX[ipl]-_planeX[ipl])/_planeEx[ipl]/_planeEx[ipl] ;

        if(_planeEy[ipl]>0.)
          chi2+=(_fitY[ipl]-_planeY[ipl])*(_fitY[ipl]-_planeY[ipl])/_planeEy[ipl]/_planeEy[ipl] ;

      }

  // Scattering angles
  // Use approximate formulas, corresponding to the approximation
  // used in fitting algorithm

  for(int ipl=1; ipl<_nTelPlanes-1;ipl++)
    {
      double th1,th2,dth;

      th2=(_fitX[ipl+1]-_fitX[ipl])*_planeDist[ipl] ;
      th1=(_fitX[ipl]-_fitX[ipl-1])*_planeDist[ipl-1] ;
      //    dth=atan(th2)-atan(th1) ;
      dth= th2 - th1 ;
      chi2 += _planeScat[ipl] * dth * dth;

      th2=(_fitY[ipl+1]-_fitY[ipl])*_planeDist[ipl] ;
      th1=(_fitY[ipl]-_fitY[ipl-1])*_planeDist[ipl-1] ;
      //    dth=atan(th2)-atan(th1) ;
      dth= th2 - th1 ;
      chi2 += _planeScat[ipl] * dth * dth;
    }

  // Beam constraint

  if(_useBeamConstraint)
    {
      double dth;

      //    dth=atan((_fitX[1]-_fitX[0])*_planeDist[0]) ;
      dth=(_fitX[1]-_fitX[0])*_planeDist[0] ;
      chi2 += _planeScat[0] * dth * dth;

      //    dth=atan((_fitY[1]-_fitY[0])*_planeDist[0]) ;
      dth=(_fitY[1]-_fitY[0])*_planeDist[0] ;
      chi2 += _planeScat[0] * dth * dth;
    }


  return chi2 ;
}



int EUTelTestFitter::GaussjSolve(double *alfa,double *beta,int n)
{
  int *ipiv;
  int *indxr;
  int *indxc;
  int i,j,k;
  int irow=0;
  int icol=0;
  double abs,big,help,pivinv;

  ipiv = new int[n];
  indxr = new int[n];
  indxc = new int[n];

  for(i=0;i<n;ipiv[i++]=0);

  for(i=0;i<n;i++)
    {
      big=0.;
      for(j=0;j<n;j++)
        {
          if(ipiv[j]==1)continue;
          for(k=0;k<n;k++)
            {
              if(ipiv[k]!=0)continue;
              abs=fabs(alfa[n*j+k]);
              if(abs>big)
                {
                  big=abs;
                  irow=j;
                  icol=k;
                }
            }
        }
      ipiv[icol]++;

      if(ipiv[icol]>1)
        return 1;

      if(irow!=icol)
        {
          help=beta[irow];
          beta[irow]=beta[icol];
          beta[icol]=help;
          for(j=0;j<n;j++)
            {
              help=alfa[n*irow+j];
              alfa[n*irow+j]=alfa[n*icol+j];
              alfa[n*icol+j]=help;
            }
        }
      indxr[i]=irow;
      indxc[i]=icol;

      if(alfa[n*icol+icol]==0.)
        return 1;

      help=alfa[n*icol+icol];
      pivinv=1./help;
      alfa[n*icol+icol]=1.;
      for(j=0;j<n;alfa[n*icol+(j++)]*=pivinv);
      beta[icol]*=pivinv;

      for(j=0;j<n;j++)
        {
          if(j==icol)continue;
          help=alfa[n*j+icol];
          alfa[n*j+icol]=0.;
          for(k=0;k<n;k++)
            alfa[n*j+k]-=alfa[n*icol+k]*help;
          beta[j]-=beta[icol]*help;
        }
    }

  for(i=n-1;i>=0;i--)
    {
      if(indxr[i]==indxc[i])continue;
      for(j=0;j<n;j++)
        {
          help=alfa[n*j+indxr[i]];
          alfa[n*j+indxr[i]]=alfa[n*j+indxc[i]];
          alfa[n*j+indxc[i]]=help;
        }
    }

  delete [] ipiv;
  delete [] indxr;
  delete [] indxc;

  return 0;
}

#endif
