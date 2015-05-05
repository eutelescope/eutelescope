// Version: $Id$
        
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifdef USE_GEAR

// eutelescope includes
#include "EUTelTestFitter.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelHistogramManager.h"
#include "EUTelExceptions.h"
#include "EUTelReferenceHit.h"


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

// LCIO includes <.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/SimTrackerHitImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>

// system includes <>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <limits>

// ROOT includes ".h"
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TVector3.h>
#include <TRotation.h>
#else
#error *** You need ROOT to compile this code.  *** 
#endif

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


EUTelTestFitter::EUTelTestFitter()
: Processor("EUTelTestFitter"),
  _isFirstEvent(false),
  _siPlanesParameters(NULL),
  _siPlanesLayerLayout(NULL),
  _histoInfoFileName(""),
  _inputColName(""),
  _outputTrackColName(""),
  _correctedHitColName(""),
  _outputHitColName(""),
  _alignmentCollectionNames(),
  _InputHitsInTrack(false),
  _OutputHitsInTrack(false),
  _SkipLayerIDs(),
  _PassiveLayerIDs(),
  _AlignLayerIDs(),
  _AlignLayerShiftX(),
  _AlignLayerShiftY(),
  _AlignLayerRotZ(),
  _WindowLayerIDs(),
  _WindowMinX(),
  _WindowMaxX(),
  _WindowMinY(),
  _WindowMaxY(),
  _MaskLayerIDs(),
  _MaskMinX(),
  _MaskMaxX(),
  _MaskMinY(),
  _MaskMaxY(),
  _resolutionX(),
  _resolutionY(),
  _resolutionZ(),
  _referenceHitCollectionName(""),
  _useReferenceHitCollection(false),
  _referenceHitVec(NULL),
  _allowMissingHits(0),
  _allowSkipHits(0),
  _maxPlaneHits(0),
  _searchMultipleTracks(false),
  _allowAmbiguousHits(false),
  _maximumAmbiguousHits(false),
  _missingHitPenalty(0.0),
  _skipHitPenalty(0.0),
  _chi2Max(0.0),
  _chi2Min(0.0),
  _useNominalResolution(false),
  _useDUT(false),
  _useBeamConstraint(false),
  _beamSpread(0.0),
  _beamSlopeX(0.0),
  _beamSlopeY(0.0),
  _eBeam(0.0),
  _nTelPlanes(0),
  _nActivePlanes(0),
  _iDUT(0),
  _planeSort(NULL),
  _planeID(NULL),
  _planeShiftX(NULL),
  _planeShiftY(NULL),
  _planeRotZ(NULL),
  _planePosition(NULL),
  _planeThickness(NULL),
  _planeX0(NULL),
  _planeResolution(NULL),
  _isActive(NULL),
  _planeWindowIDs(NULL),
  _planeMaskIDs(NULL),
  _nRun(0),
  _nEvt(0),
  _planeHits(NULL),
  _planeChoice(NULL),
  _planeMod(NULL),
  _planeX(NULL),
  _planeEx(NULL),
  _planeY(NULL),
  _planeEy(NULL),
  _planeScatAngle(NULL),
  _planeDist(NULL),
  _planeScat(NULL),
  _fitX(NULL),
  _fitEx(NULL),
  _fitY(NULL),
  _fitEy(NULL),
  _fitArray(NULL),
  _nominalFitArrayX(NULL),
  _nominalErrorX(NULL),
  _nominalFitArrayY(NULL),
  _nominalErrorY(NULL),
  _noOfEventWOInputHit(0),
  _noOfEventWOTrack(0),
  _noOfTracks(0),
  _aidaHistoMap(),
  _aidaHistoMap1D(),
  _aidaHistoMap2D(),
  _UseSlope(false),
  _SlopeXLimit(0.0),
  _SlopeYLimit(0.0),
  _SlopeDistanceMax(0.0),
  _fittedXcorr(),
  _fittedYcorr(),
  _fittedZcorr(),
  _indexDUTneighbour(0),
  _zDUTneighbour(0.0),
  _siPlaneCenter(),
  _siPlaneNormal() {

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

  // alignment collections (might be important) 
  EVENT::StringVec	_alignmentCollectionExample;
  _alignmentCollectionExample.push_back("alignment");
  
  registerProcessorParameter ("alignmentCollectionNames",
                            "List of alignment collections which are neede to get track position on a Sensor surface ",
                            _alignmentCollectionNames, _alignmentCollectionExample );

  // compulsory parameters:

  registerProcessorParameter ("InputHitsInTrack",
                              "Flag for storing input (measured) hits in track",
                              _InputHitsInTrack,  static_cast < bool > (false));

  registerProcessorParameter ("OutputHitsInTrack",
                              "Flag for storing output (fitted) hits in track",
                              _OutputHitsInTrack,  static_cast < bool > (true));

  registerProcessorParameter ("AllowMissingHits",
                              "Allowed number of missing hits in the track",
                              _allowMissingHits,  static_cast < int > (0));

  registerProcessorParameter ("AllowSkipHits",
                              "Allowed number of hits removed from the track",
                              _allowSkipHits,  static_cast < int > (0));

  registerProcessorParameter ("MaxPlaneHits",
                              "Maximum number of considered hits per plane",
                              _maxPlaneHits,  static_cast < int > (100));

  registerProcessorParameter ("MissingHitPenalty",
                              "Chi2 penalty for missing hit in the track",
                              _missingHitPenalty,  static_cast < double > (0.));

  registerProcessorParameter ("SkipHitPenalty",
                              "Chi2 penalty for removing hit from the track",
                              _skipHitPenalty,  static_cast < double > (100.));

  registerProcessorParameter ("Chi2Min",
                              "Minimum Chi2 for accepted track fit",
                              _chi2Min,  static_cast < double > (0.0));

  registerProcessorParameter ("Chi2Max",
                              "Maximum Chi2 for accepted track fit",
                              _chi2Max,  static_cast < double > (100.));

  registerProcessorParameter ("UseNominalResolution",
                              "Flag for using nominal resolution instead of position errors",
                              _useNominalResolution,  static_cast < bool > (true));

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

  registerOptionalParameter("ReferenceCollection","reference hit collection name ", _referenceHitCollectionName, static_cast <string> ("referenceHit") );

  registerOptionalParameter("UseReferenceCollection","Do you want the reference hit collection to be used for coordinate transformations?",  _useReferenceHitCollection, static_cast< bool   > ( true ));
 
 
  // ------- Parameters added to allow correlation band info 02 August 2010 libov@mail.desy.de -------
  registerOptionalParameter("UseSlope","Use expected track direction to constraint number of considered hit combinations (track preselection).", _UseSlope, true );
  registerOptionalParameter("SlopeXLimit","Limit on track slope change when passing sensor layer (in X direction)", _SlopeXLimit, static_cast <float> (0.001));
  registerOptionalParameter("SlopeYLimit","Limit on track slope change when passing sensor layer (in Y direction)", _SlopeYLimit, static_cast <float> (0.001));
  registerOptionalParameter("SlopeDistanceMax","Maximum hit distance from the expected position, used for hit preselection in [mm]", _SlopeDistanceMax, static_cast <float> (1.));
  // -------------------------------------------------------------------------------------------------

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
                             "Assumed angular spread of the beam [rad] (for beam constraint)",
                             _beamSpread, static_cast < double > (0.0) );

  registerOptionalParameter ("BeamSlopeX",
                             "Beam direction tilt in X-Z plane [rad] (for beam constraint)",
                             _beamSlopeX,  static_cast < double > (0.));

  registerOptionalParameter ("BeamSlopeY",
                             "Beam direction tilt in Y-Z plane [rad] (for beam constraint)",
                             _beamSlopeY,  static_cast < double > (0.));

  registerOptionalParameter ("SearchMultipleTracks",
                             "Flag for searching multiple tracks in events with multiple hits",
                             _searchMultipleTracks,  static_cast < bool > (true));

  registerOptionalParameter ("AllowAmbiguousHits",
                             "Allow same hit to be used in more than one track",
                             _allowAmbiguousHits, static_cast < bool > (false));

  registerOptionalParameter ("MaximumAmbiguousHits",
                             "Maximum number of hits to be shared by more than one track",
                             _maximumAmbiguousHits, static_cast < int > (2));

  registerOptionalParameter("ResolutionX","X resolution parameter for each plane. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_resolutionX,  std::vector<float> (static_cast <int> (6), 10.));

  registerOptionalParameter("ResolutionY","Y resolution parameter for each plane. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_resolutionY,std::vector<float> (static_cast <int> (6), 10.));

  registerOptionalParameter("ResolutionZ","Z resolution parameter for each plane. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_resolutionZ,std::vector<float> (static_cast <int> (6), 10.));


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

  _isFirstEvent = true;

  _referenceHitVec = 0;


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

  if( !_SkipLayerIDs.empty() ) {
    streamlog_out ( MESSAGE0 )  <<  _SkipLayerIDs.size() << " layers should be skipped" << endl;
  }

  if( !_PassiveLayerIDs.empty() ) {
    streamlog_out ( MESSAGE0 ) <<  _PassiveLayerIDs.size() << " layers should be considered passive" << endl;
  }

  // Take all layers defined in GEAR geometry
  _nTelPlanes = _siPlanesLayerLayout->getNLayers();

  // Check for DUT
  if( _siPlanesParameters->getSiPlanesType()==_siPlanesParameters->TelescopeWithDUT )    {
    _iDUT = _nTelPlanes ;
    _nTelPlanes++;
  }  else {
    _iDUT = -1 ;
  }


  // _nTelPlanes is the total number of sensitive layers in the setup
  // summing both the telescopes and the DUT 
  // Read position in Z (for sorting), skip layers if requested

  _planeSort = new int[_nTelPlanes];
  _planePosition   = new double[_nTelPlanes];

  int nSkip=0;

  for(int ipl=0; ipl <  _siPlanesLayerLayout->getNLayers(); ipl++)    {
    _planePosition[ipl-nSkip]=_siPlanesLayerLayout->getLayerPositionZ(ipl);
    _planeSort[ipl-nSkip]=ipl;

       
    // Check if not on "skip list"

    int _plID = _siPlanesLayerLayout->getID(ipl);

    for(int spl=0; spl< static_cast<int>(_SkipLayerIDs.size()) ; spl++) {
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
  //
  // Note: all this code can in principle be replaced by a simple sort 
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
  //
  // N.B. the plane ID array is containing the sensorID sorted
  // according to the position along Z. We will be using this array
  // for histogram booking and filling.
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
      _planeThickness[iz]=_siPlanesLayerLayout->getLayerThickness(ipl) + _siPlanesLayerLayout->getSensitiveThickness(ipl) ;
      _planeX0[iz]=_siPlanesLayerLayout->getLayerRadLength(ipl) ;
      resolution = _siPlanesLayerLayout->getSensitiveResolution(ipl);
    }  else {
      _planeID[iz]=_siPlanesLayerLayout->getDUTID();
      _planeThickness[iz]=_siPlanesLayerLayout->getLayerThickness(ipl) + _siPlanesLayerLayout->getSensitiveThickness(ipl) ;
      _planeX0[iz]=_siPlanesLayerLayout->getDUTRadLength();
      resolution = _siPlanesLayerLayout->getDUTSensitiveResolution();
    }

    iActive = (resolution > 0);

    // Check passive layer list
    for(int ppl=0; ppl< static_cast<int>(_PassiveLayerIDs.size()) ; ppl++) {
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

    // Alignment corrections should already be added in HitMaker
    // But look in input options for possible adjustments

    _planeShiftX[iz]=0.;
    _planeShiftY[iz]=0.;
    _planeRotZ[iz]=0.;

    for(int apl=0; apl< static_cast<int>(_AlignLayerIDs.size()) ; apl++) {
      if ( _AlignLayerIDs.at(apl) == _planeID[iz]) {
        _planeShiftX[iz]=_AlignLayerShiftX.at(apl);
        _planeShiftY[iz]=_AlignLayerShiftY.at(apl);
        // Rotation can be skipped: check size
        if(apl < static_cast<int>(_AlignLayerRotZ.size())) {
          _planeRotZ[iz]=_AlignLayerRotZ.at(apl);
        }
        break;
      }
    }

    // Check, if there are additional cuts defined for this plane
    for(int wpl=0; wpl< static_cast<int>(_WindowLayerIDs.size()) ; wpl++) {
      if ( _WindowLayerIDs.at(wpl) == _planeID[iz]) {
        _planeWindowIDs[iz].push_back(wpl);
      }
    }
    for(int mpl=0; mpl< static_cast<int>(_MaskLayerIDs.size()) ; mpl++) {
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

        for(int wpl=0; wpl< static_cast<int>(_planeWindowIDs[ipl].size()) ; wpl++)  {
          int iwin= _planeWindowIDs[ipl].at(wpl);
          ss << "\n accepted window: X = " <<  _WindowMinX.at(iwin)
             << " to " <<  _WindowMaxX.at(iwin)
             << " Y = " <<  _WindowMinY.at(iwin)
             << " to " <<  _WindowMaxY.at(iwin);
        }

        for(int mpl=0; mpl< static_cast<int>(_planeMaskIDs[ipl].size()) ; mpl++)  {
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
  _planeMod    = new type_fitcount[_nTelPlanes];

  _planeX  = new double[_nTelPlanes];
  _planeEx = new double[_nTelPlanes];
  _planeY  = new double[_nTelPlanes];
  _planeEy = new double[_nTelPlanes];

  _planeScatAngle  = new double[_nTelPlanes];

  _planeDist = new double[_nTelPlanes];
  _planeScat = new double[_nTelPlanes];

  _fitX  = new double[_nTelPlanes];
  _fitEx = new double[_nTelPlanes];
  _fitY  = new double[_nTelPlanes];
  _fitEy = new double[_nTelPlanes];

  int arrayDim = _nTelPlanes * _nTelPlanes;

  _fitArray = new double[arrayDim];
  _nominalFitArrayX = new double[arrayDim];
  _nominalErrorX = new double[_nTelPlanes];

  _nominalFitArrayY = new double[arrayDim];
  _nominalErrorY = new double[_nTelPlanes];

  // Fill nominal fit matrices and
  // calculate expected precision of track fitting

  // Planes are ordered in position along the beam line !

  double totalScatAngle = 0.;

  for(int ipl=0; ipl<_nTelPlanes ; ipl++) 
  {
    if(ipl>0) 
    {
      _planeDist[ipl-1]=1./(_planePosition[ipl] - _planePosition[ipl-1]) ;
    }
    
    _planeScatAngle[ipl]= 0.0136/_eBeam * sqrt(_planeThickness[ipl]/_planeX0[ipl])
      * (1.+0.038*std::log(_planeThickness[ipl]/_planeX0[ipl])) ;
  
 
    if(ipl==0 && _useBeamConstraint) {
      _planeScat[ipl]= 1./(_planeScatAngle[ipl]*_planeScatAngle[ipl]+ _beamSpread*_beamSpread) ;
    } else {
      _planeScat[ipl]= 1./(_planeScatAngle[ipl] * _planeScatAngle[ipl]) ;
    }

    totalScatAngle+= _planeScatAngle[ipl] * _planeScatAngle[ipl];

    if(streamlog_level(DEBUG5)){
      streamlog_out( DEBUG5 ) << "Scattering angle in plane " << ipl << ": " << _planeScatAngle[ipl] << endl;
    }
    _fitX[ipl] =_fitY[ipl] = 0. ;
    if(static_cast<int>(_resolutionX.size()) < ipl+1 )
    {
      _nominalErrorX[ipl]= _planeResolution[ipl];
    }
    else
    {
      _nominalErrorX[ipl]= _resolutionX[ipl];
    }
    if(static_cast<int>(_resolutionY.size()) < ipl+1 )
    {
      _nominalErrorY[ipl]= _planeResolution[ipl];
    }
    else
    {
      _nominalErrorY[ipl]= _resolutionY[ipl];
    }
 
  }

  totalScatAngle = sqrt(totalScatAngle);

  // Fit with nominal parameters for X direction

  int status = DoAnalFit(_fitX,_nominalErrorX,_beamSlopeX);

  if(status) {
    streamlog_out( ERROR2 ) << "\n Fit in X with nominal geometry failed !?!" << endl;
  }
  // Store fit matrix

  for(int imx=0; imx<arrayDim; imx++) {
    _nominalFitArrayX[imx] = _fitArray[imx];
  }

  stringstream ss;
  ss << "Expected position resolutions in X [um]: ";
  for(int ipl=0; ipl<_nTelPlanes ; ipl++) {
    ss << _nominalErrorX[ipl]*1000. << "  " ;
  }

  ss << endl << "Expected scattering angle [mrad]: ";
  for(int ipl=0; ipl<_nTelPlanes ; ipl++) {
    ss << _planeScatAngle[ipl]*1000. << "  " ;
  }

  ss << endl << "Expected total scattering angle [mrad]: " << totalScatAngle*1000. ;

  streamlog_out ( MESSAGE2 ) << ss.str() << endl;


  // Fit with nominal parameters for Y direction

  status = DoAnalFit(_fitY,_nominalErrorY,_beamSlopeY);

  if(status) {
    streamlog_out( ERROR2 ) << "\n Fit in Y with nominal geometry failed !?!" << endl;
  }
  // Store fit matrix

  for(int imx=0; imx<arrayDim; imx++) {
    _nominalFitArrayY[imx] = _fitArray[imx];
  }

// Check if slope-based preselection parameter values are not too small

  if( _UseSlope && 
        _SlopeXLimit < 5.*totalScatAngle 
     )
    streamlog_out( ERROR2 ) << "SlopeXLimit cut probably too tight! Check parameters!" << endl; 
  

  if( _UseSlope && 
        _SlopeYLimit < 5.*totalScatAngle 
     )
    streamlog_out( ERROR2 ) << "SlopeXLimit cut probably too tight! Check parameters!" << endl; 
  
  // Take into account beam spread

  totalScatAngle = sqrt(totalScatAngle*totalScatAngle + _beamSpread*_beamSpread);

  if( _UseSlope && 
  _SlopeDistanceMax < 5.*totalScatAngle * (_planePosition[_nTelPlanes-1] - _planePosition[0])
     )
    streamlog_out( ERROR2 ) << "SlopeDistanceMax cut probably too tight! Check parameters!" << endl; 
  

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

  streamlog_out( MESSAGE4 )  << "Processing run header " << _nRun
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

  _nEvt ++ ;

  EUTelEventImpl * euEvent = static_cast<EUTelEventImpl*> ( event );
  if ( euEvent->getEventType() == kEORE ) 
  {
    streamlog_out ( DEBUG5 ) <<  "EORE found: nothing else to do." << endl;
    return;
  }


  LCCollection* col;
  try 
  {
    col = event->getCollection( _inputColName ) ;
  }
  catch (lcio::DataNotAvailableException& e) 
  {
    streamlog_out ( DEBUG5 ) << "Not able to get collection "
                            << _inputColName
                            << "\nfrom event " << event->getEventNumber()
                            << " in run " << event->getRunNumber()  << endl;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nTrackHistoName]))->fill(0);
#endif
    ++_noOfEventWOInputHit;
    return;
  }

  // setup cellIdDecoder to decode the sensor ID from the hits
  CellIDDecoder<TrackerHit>  hitCellDecoder(EUTELESCOPE::HITENCODING);


  if(_isFirstEvent)
  {
       if ( _useReferenceHitCollection ) 
       {
         _referenceHitVec = dynamic_cast < LCCollectionVec * > (event->getCollection( _referenceHitCollectionName));
       }
 
      // apply all GEAR/alignment offsets to get corrected X,Y,Z position of the
      // sensor center
      //
      //

      for(int iplane=0; iplane <  _siPlanesLayerLayout->getNLayers(); iplane++)    
      {

          map< unsigned int , double > _planeCenter;
          map< unsigned int , double > _planeNormal;
 
          int sensorID = _siPlanesLayerLayout->getID(iplane)                   ;

          // start filling the map with Gear values:

          _planeCenter[ 0 ] =  _siPlanesLayerLayout->getLayerPositionX(iplane) ; // X
          _planeCenter[ 1 ] =  _siPlanesLayerLayout->getLayerPositionY(iplane) ; // Y
          _planeCenter[ 2 ] =  _siPlanesLayerLayout->getLayerPositionZ(iplane) ; // Z

          _planeNormal[ 0 ] =  0.                                              ; // X
          _planeNormal[ 1 ] =  0.                                              ; // Y
          _planeNormal[ 2 ] =  1.                                              ; // Z

          TVector3  _normalTVec( _planeNormal[0], _planeNormal[1], _planeNormal[2]); 

          // do initial rotation (from GEAR)
          try
                      {
 
                          double gRotation[3] = { 0., 0., 0.}; // not rotated
                         
                          gRotation[0] = _siPlanesLayerLayout->getLayerRotationXY(iplane); // Euler alpha ;
                          gRotation[1] = _siPlanesLayerLayout->getLayerRotationZX(iplane); // Euler alpha ;
                          gRotation[2] = _siPlanesLayerLayout->getLayerRotationZY(iplane); // Euler alpha ;
                          
                          // input angles are in DEGREEs !!!
                          // translate into radians
                          gRotation[0] =  gRotation[0]*3.1415926/180.; // 
                          gRotation[1] =  gRotation[1]*3.1415926/180.; //
                          gRotation[2] =  gRotation[2]*3.1415926/180.; //

                          TRotation	r;
                          r.RotateZ( gRotation[0] );                          
                          r.RotateY( gRotation[1] );                            
                          r.RotateX( gRotation[2] );                          

                          _normalTVec.Transform( r );
                                  
                          _planeNormal[0] = _normalTVec[0];
                          _planeNormal[1] = _normalTVec[1];
                          _planeNormal[2] = _normalTVec[2];
                      }
                      catch(...)
                      {
			streamlog_out( MESSAGE5 ) << " no sensor rotation is given in the GEAR steering file, assume NONE" << endl;
                      }

          if( _alignmentCollectionNames.size() > 0 )
          {
              for( size_t i=0; i < _alignmentCollectionNames.size(); i++)
              {

                  //
                  // if offsets are given via the alignment collections apply:
                  //    

                  try
                  {
                      LCCollectionVec * alignmentCollectionVec     = dynamic_cast < LCCollectionVec * > (event->getCollection(_alignmentCollectionNames[i]));
                      // next, find the alignment constant corresponding to the DUT

                      EUTelAlignmentConstant * c=NULL;
                      for ( size_t iPos = 0; iPos < alignmentCollectionVec->size(); ++iPos ) 
                      {
                          c = static_cast< EUTelAlignmentConstant * > ( alignmentCollectionVec->getElementAt( iPos ) );
                          if (c -> getSensorID() == sensorID  )
                          {                
                              break;	         
                          }
                      }

		     if(c->getSensorID() == sensorID)
		     {
                       _planeCenter[ 0 ] += c->getXOffset() ;
                       _planeCenter[ 1 ] += c->getYOffset() ;
                       _planeCenter[ 2 ] += c->getZOffset() ;
                

                       double gRotation[3] = { 0., 0., 0.}; // not rotated

                       try
                       {
                           gRotation[0] = c->getGamma();
                           gRotation[1] = c->getBeta();
                           gRotation[2] = c->getAlpha();

                           TRotation	r;
                           r.RotateZ( gRotation[0] );                          
                           r.RotateY( gRotation[1] );                          
                           r.RotateX( gRotation[2] );                          
 
                           _normalTVec.Transform( r );
                                   
                        }
                        catch(...)
                        {
 			streamlog_out( MESSAGE5 ) << " no sensor rotation is given in the GEAR steering file, assume NONE" << endl;
 	                }

                       _planeNormal[0] = _normalTVec[0];
                       _planeNormal[1] = _normalTVec[1];
                       _planeNormal[2] = _normalTVec[2];

		    }
		    else {
			streamlog_out( DEBUG5 ) << "Wrong telescope plane, not applying offsets." << endl;
		    }                      
                  }
                  catch(...)
                  {
		    streamlog_out( WARNING2 ) << "Collection " << _alignmentCollectionNames[i].c_str() << " not found "  << endl;
                  }
              }            

                         
          }
 
          _siPlaneCenter[ sensorID]=  _planeCenter ;
          _siPlaneNormal[ sensorID]=  _planeNormal ;
     }
      
      if ( isFirstEvent() ) _isFirstEvent = false;
  }

  // Copy hits to local table
  // Assign hits to sensor planes
  // =============================


  int nHit = col->getNumberOfElements()  ;

  if(streamlog_level(DEBUG5)){
    streamlog_out( DEBUG5 )  << "Total of " << nHit << " tracker hits in input collection " << endl;
  }

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nAllHitHistoName]))->fill(nHit);
#endif

  if(nHit + _allowMissingHits < _nActivePlanes) 
  {
    streamlog_out ( DEBUG5 ) << "Not enough hits to perform the fit, exiting... " << endl;

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
  int    * hitFits  = new int[nHit];

  IntVec * planeHitID   = new IntVec[_nTelPlanes];

  // Loop over hits

  int nGoodHit = 0;

  for(int ihit=0; ihit< nHit ; ihit++) 
  {
    TrackerHit        * meshit = dynamic_cast<TrackerHit*>( col->getElementAt(ihit) ) ;
    SimTrackerHitImpl * simhit = dynamic_cast<SimTrackerHitImpl*>( col->getElementAt(ihit) ) ;

    // Hit position
    //
    double pos[3]={0.,0.,0.};
    if( meshit != 0 )
    {
      const double *pos0 = meshit->getPosition();
      pos[0] = pos0[0];
      pos[1] = pos0[1];
      pos[2] = pos0[2];
    } 
    else if( meshit == 0 && simhit != 0 )
    {
      const double *pos0 = simhit->getPosition(); 
      pos[0] = pos0[0];
      pos[1] = pos0[1];
      pos[2] = pos0[2];
    }
   
    hitZ[ihit] = pos[2];

    hitFits[ihit]=-1;

    // We have to find Plane ID of the hit
    // by looking at the Z position
    //
    double distMin =  1.;
    hitPlane[ihit] = -1 ;

    for(int ipl=0;ipl<_nTelPlanes;ipl++)  
    {
      double dist =  hitZ[ihit] - _planePosition[ipl] ;
      if( dist < 0 ) dist = -dist; // always positive defined !
       
      if(     dist  <   distMin    )  
      {
        hitPlane[ihit] = ipl;
        distMin        = dist ;
      }
    }

    // Ignore hits not matched to any plane
    //
    if(hitPlane[ihit]<0) 
    {
      if ( _isActive[hitPlane[ihit]] ) 
      {
          streamlog_out( WARNING0 )  << "Reconstructed hit outside sensor planes z [mm] = "  << hitZ[ihit] << endl;
      }
      continue;
    }

    // Ignore hit, if plane not declared as active (i.e. not used in the fit)
    // 
    if(! _isActive[hitPlane[ihit]]) 
    {
      continue ;
    }

    // Ignore hit also, if maximum number of hits already matched to this plane
    // 
    if(_maxPlaneHits>0 && _maxPlaneHits-planeHitID[hitPlane[ihit]].size()<=0) 
    {
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

    for(int wpl=0; wpl< static_cast<int>(_planeWindowIDs[hitPlane[ihit]].size()) ; wpl++)  
    {
      int iwin= _planeWindowIDs[hitPlane[ihit]].at(wpl);

      if(hitX[ihit] <  _WindowMinX.at(iwin)) hitcut = true;
      if(hitX[ihit] >  _WindowMaxX.at(iwin)) hitcut = true;
      if(hitY[ihit] <  _WindowMinY.at(iwin)) hitcut = true;
      if(hitY[ihit] >  _WindowMaxY.at(iwin)) hitcut = true;
    }

    if(hitcut) continue;


    for(int mpl=0; mpl< static_cast<int>(_planeMaskIDs[hitPlane[ihit]].size()) ; mpl++) 
    {
      int imsk= _planeMaskIDs[hitPlane[ihit]].at(mpl);

      if(
              hitX[ihit] >  _MaskMinX.at(imsk) &&  hitX[ihit] <  _MaskMaxX.at(imsk) &&
              hitY[ihit] >  _MaskMinY.at(imsk) &&  hitY[ihit] <  _MaskMaxY.at(imsk)
              ) 
      {
        hitcut = true;
      }
    }

    if(hitcut)continue;


    // Add hit to hit list for given plane - to be used in track selection
    //
    planeHitID[hitPlane[ihit]].push_back(ihit);
    nGoodHit++;

    // Position uncertainty. Use nominal resolution if not properly defined
    //
    EVENT::FloatVec cov(3);
    cov[0] = 0.;
    cov[1] = 0.;
    cov[2] = 0.;
    if( meshit != 0 )
    {
      const EVENT::FloatVec cov0 = meshit->getCovMatrix();
      cov[0]=cov0[0];
      cov[1]=cov0[1];
      cov[2]=cov0[2];
    } 
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

    if(streamlog_level(DEBUG5)){
      streamlog_out ( DEBUG5 ) << "Hit " << ihit
			       << "   X = " << hitX[ihit] << " +/- " << hitEx[ihit]
			       << "   Y = " << hitY[ihit] << " +/- " << hitEy[ihit]
			       << "   Z = " << hitZ[ihit] << " (plane " << hitPlane[ihit] << ")" << endl;
    }
  }

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nAccHitHistoName]))->fill(nGoodHit);
#endif

  // Main analysis loop: finding multiple tracks (if allowed)
  //=========================================================

  // Define output track and hit collections

  LCCollectionVec     *fittrackvec = new LCCollectionVec(LCIO::TRACK);
  LCCollectionVec     *fitpointvec = new LCCollectionVec(LCIO::TRACKERHIT);
  LCCollectionVec     *corrpointvec = new LCCollectionVec(LCIO::TRACKERHIT);

  // prepare an encoder for the hit collection
  CellIDEncoder<TrackerHitImpl> fitHitEncoder(EUTELESCOPE::HITENCODING, fitpointvec);
  CellIDEncoder<TrackerHitImpl> corrHitEncoder(EUTELESCOPE::HITENCODING, corrpointvec);


  // Set flag for storing track hits in track collection

  LCFlagImpl flag(fittrackvec->getFlag());
  flag.setBit( LCIO::TRBIT_HITS );
  fittrackvec->setFlag(flag.getFlag());


  //
  // New implementation: only one loop over all fit possibilities! 
  // =============================================================
  //
  // Method works in all cases, also when missing hits are allowed 
  // Duplicated tracks (if ambiguity is not allowed) rejected later


  // Store all fits passing cuts
  // Chi2 map will sort all possibilities according to Chi2 value

  std::multimap<double,int> fittedChi2;
  std::multimap<double,int>::iterator fitIterator,fitIterator2;


  // Vectors storing fit results (one number per fit)

  std::vector<double> fittedPenalty;
  std::vector<int> fittedFired;


  // Vectors storing fit results (_nTelPlanes numbers per fit)

  std::vector<double> fittedX;
  std::vector<double> fittedEx;
  std::vector<double> fittedY;
  std::vector<double> fittedEy;
  std::vector<int> fittedHits;

  // Total number of fitted tracks stored in vectors

  int nFittedTracks = 0 ;


    // Count planes active in this event and number of fit possibilities
    //
    int nFiredPlanes = 0;
    type_fitcount nChoice = 1;

    
    // Count from the last plane, to allow for "smart" track finding
    //
    for(int ipl=_nTelPlanes-1; ipl>=0 ;ipl--)   
    {
      _planeHits[ipl] = planeHitID[ipl].size() ;

      if(_planeHits[ipl]>0)  
      {
        nFiredPlanes++;

        _planeChoice[ipl]=_planeHits[ipl]+1;
      }
      else 
      {
        _planeChoice[ipl]=1;
      }

      _planeMod[ipl]=nChoice;
      nChoice*=_planeChoice[ipl];
    }
    

    // Debug output
    if(streamlog_level(DEBUG5)){
      for(int ipl=0;ipl<_nTelPlanes;ipl++) 
      {
        if( _isActive[ipl] )  
        {
          stringstream ss;
          ss << "Plane " << ipl << "  " << _planeHits[ipl] << " hit(s), hit IDs :";

          for( int ihit=0; ihit < static_cast<int>( planeHitID[ipl].size()) ; ihit ++) 
          {
            ss << planeHitID[ipl].at(ihit) << " " ;
          }
          streamlog_out ( DEBUG5 )  << ss.str() << endl;
        }
      }
    }


    // Check if fit can be done

    if(nFiredPlanes + _allowMissingHits < _nActivePlanes) 
    {
      if(streamlog_level(DEBUG5)){
	streamlog_out ( DEBUG5 ) <<  "Not enough planes hit to perform the fit " << endl;
      }

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nTrackHistoName]))->fill(0);
#endif

      // before returning clean up the memory
      delete fittrackvec;
      delete fitpointvec;
      delete corrpointvec;
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


    if(streamlog_level(DEBUG5))
    {
      streamlog_out ( DEBUG5 ) << nFiredPlanes << " active sensor planes hit, checking "
                                          << nChoice << " fit possibilities "  << endl;
    }

    // Check all track possibilities

    double chi2min  = numeric_limits<double >::max();

    // Loop over fit possibilities
    // Start from one-hit track to allow for "smart" skipping of wrong matches

    int istart=0;
    int nmiss=_allowMissingHits;

    while(nmiss>0 || !_isActive[istart])  
    {
      if(_isActive[istart]) 
      {
        nmiss--;
      }
      istart++;
    }

    
    for(type_fitcount ichoice = nChoice-_planeMod[istart]-1; ichoice >= 0; ichoice--)  
    {        
      int    nChoiceFired =  0 ;
      double choiceChi2   = -1.;
      double trackChi2    = -1.;
      int    ifirst       = -1 ;
      int    ilast        =  0 ;
      int    nleft        =  0 ;
     
      // New variables for preselection based on slope
      // will be set to plane number if
      //   - hit too far from the expected position (based on first
      //             plane + beam slope): hit missed
      //   - angle between track segments (slope change) too large:
      //                  track slope
      //
      // Value >0 gives first layer which failed the cut
      // 0 value means that preselection cuts were passed by all hits

      int firstHitMissed = 0;
      int firstTrackSlope = 0;
     
      // If beam constraint used: assume the track should go along
      // beam direction, otherwise beam is assumed to be perpendicular
      // to the sensor plane

      double expTrackSlopeX=0.;
      double expTrackSlopeY=0.;
 
      if(_useBeamConstraint)
	{
	  expTrackSlopeX=_beamSlopeX;
	  expTrackSlopeY=_beamSlopeY;
	}

      double lastSlopeX=0.;
      double lastSlopeY=0.;
 
      // Fill position and error arrays for this hit configuration

      for(int ipl=0;ipl<_nTelPlanes;ipl++)  
      {
         _planeX[ipl] = _planeY[ipl] = _planeEx[ipl] = _planeEy[ipl] = 0.;

        if(_isActive[ipl])  
        {
          int ihit   = (ichoice/_planeMod[ipl])%_planeChoice[ipl];
        

          if(ihit<_planeHits[ipl])    
          {
            int jhit      = planeHitID[ipl].at(ihit);
            
            _planeX[ipl]  = hitX[jhit];
            _planeY[ipl]  = hitY[jhit];
            _planeEx[ipl] = (_useNominalResolution)?_planeResolution[ipl]:hitEx[jhit];
            _planeEy[ipl] = (_useNominalResolution)?_planeResolution[ipl]:hitEy[jhit];
 
 
	    // Calculate distance from expected position
	    // starting from the second hit (when ifirst already set)

            if(_UseSlope && ifirst>=0 && firstHitMissed == 0)
	      {
              double expX = _planeX[ifirst] + expTrackSlopeX *(_planePosition[ipl]-_planePosition[ifirst]);
              double expY = _planeY[ifirst] + expTrackSlopeY *(_planePosition[ipl]-_planePosition[ifirst]);
              if(   abs( _planeX[ipl] - expX ) >  _SlopeDistanceMax/1000. 
                 || abs( _planeY[ipl] - expY ) >  _SlopeDistanceMax/1000.  
                   )  firstHitMissed = ipl;

	      }

	    // Calculate slope and check slope change w.r.t. previous slope

            if(_UseSlope && ifirst>=0 && firstTrackSlope==0  )
	      {
              double slopeX = (_planeX[ipl]-_planeX[ifirst])/
                               (_planePosition[ipl]-_planePosition[ifirst]);

              double slopeY = (_planeY[ipl]-_planeY[ifirst])/
                               (_planePosition[ipl]-_planePosition[ifirst]);

              if(ilast>ifirst && 
		 ( abs(slopeX - lastSlopeX) > _SlopeXLimit ||
                   abs(slopeY - lastSlopeY) > _SlopeYLimit )
		 ) firstTrackSlope=ipl;

              lastSlopeX=slopeX;
              lastSlopeY=slopeY;
	      }

         
            if(ifirst<0) 
            {            
              ifirst    = ipl;
            }
              
            ilast = ipl;
            nleft = 0;
            nChoiceFired++;
          } 
          else 
          {
            nleft++;        // Counts number of planes with missing
			    // hits after the last hit
          }
        }
      }
      // End of plane loop (decoding fit hypothesis)


      // Check number of selected hits
      // =============================

      // No fit to 1 hit :-)

      if(nChoiceFired < 2) 
        {
           continue;
        }
      // Fit with 2 hits make sense only with beam constraint, or
      // when 2 point fit is allowed

      if(       nChoiceFired==2 
                && !_useBeamConstraint
                && nChoiceFired + _allowMissingHits < _nActivePlanes     ) 
        {
           continue;
        }
      
      // Skip also if the fit can not be extended to proper number
      // of planes; no need to check remaining planes !!!

        if(nChoiceFired + nleft < _nActivePlanes - _allowMissingHits ) {
          ichoice-=_planeMod[ilast]-1;
         continue;
        }
     

      // Preselection added before full Chi2 calculation
      //
      // Cut on distance from expected position

	if(firstHitMissed>0){
          ichoice-=_planeMod[firstHitMissed]-1;
         continue;
        }
     
      // Cut on track slope changes

	if(firstTrackSlope>0){
          ichoice-=_planeMod[firstTrackSlope]-1;
         continue;
        }
     

 
      // Select fit method
      // "Nominal" fit only if all active planes used

      
      if(_useNominalResolution && (nChoiceFired == _nActivePlanes)) 
      {
        choiceChi2 = NominalFit();
      } else {
        if(_useNominalResolution && _beamSlopeX==_beamSlopeY) choiceChi2 = SingleFit();
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
      double penalty = 
          (_nActivePlanes-nFiredPlanes)*_missingHitPenalty
          +   
          (nFiredPlanes-nChoiceFired)*_skipHitPenalty ;


      trackChi2 = choiceChi2+penalty;

      if(
              nChoiceFired + _allowMissingHits >= _nActivePlanes 
              &&
              nChoiceFired + _allowSkipHits    >= nFiredPlanes  
              &&
              trackChi2 < chi2min
              ) 
      {
        chi2min=trackChi2;
      }

      // Check if better than chi2Max
      // If not: skip also all track possibilities which include
      // this hit selection !!!

      if( choiceChi2 >= _chi2Max  || choiceChi2 < _chi2Min ) 
      {        
        ichoice-=_planeMod[ilast]-1;
        continue;
      } 

      //
      // Skip fit if could not be accepted (too few planes fired)
      //

      if(
              nChoiceFired + _allowMissingHits < _nActivePlanes 
              ||
              nChoiceFired + _allowSkipHits    < nFiredPlanes 
              ) 
      {
        continue;
      }


      // Fill all tracks passing chi2 cut

      if( trackChi2 < _chi2Max && trackChi2 > _chi2Min ) 
      {

        fittedChi2.insert( make_pair( trackChi2, nFittedTracks ));

        fittedPenalty.push_back(penalty); 
        fittedFired.push_back(nChoiceFired);

        for(int ipl=0;ipl<_nTelPlanes;ipl++)  
        {
            int jhit=-1;

            if(_isActive[ipl])  
            {
                int ihit = (ichoice/_planeMod[ipl])%_planeChoice[ipl];
                
                if(ihit<_planeHits[ipl])
                {
                    jhit = planeHitID[ipl].at(ihit);
                }
            }

            fittedHits.push_back(jhit);

            fittedX.push_back(_fitX[ipl]);
            fittedY.push_back(_fitY[ipl]);
            fittedEx.push_back(_fitEx[ipl]);
            fittedEy.push_back(_fitEy[ipl]);
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    stringstream iden;
    iden << "pl" << _planeID[ipl] << "_";
    string bname = iden.str();
if(jhit>=0){
      _aidaHistoMap1D[bname + "fitX"]->fill( _fitX[ipl]  );
      _aidaHistoMap1D[bname + "fitY"]->fill( _fitY[ipl]  );
      _aidaHistoMap1D[bname + "hitX"]->fill(  hitX[jhit] );
      _aidaHistoMap1D[bname + "hitY"]->fill(  hitY[jhit] );
      _aidaHistoMap1D[bname + "residualX"]->fill( _fitX[ipl] - hitX[jhit] );
      _aidaHistoMap1D[bname + "residualY"]->fill( _fitY[ipl] - hitY[jhit] );
      //Resids 
      _aidaHistoMap2D[bname + "residualXdX"]->fill( _fitX[ipl]  , _fitX[ipl]    - hitX[jhit]  );
      _aidaHistoMap2D[bname + "residualYdX"]->fill( _fitX[ipl]  , _fitY[ipl]    - hitY[jhit]  );
      _aidaHistoMap2D[bname + "residualXdY"]->fill( _fitY[ipl]  , _fitX[ipl]    - hitX[jhit]  );
      _aidaHistoMap2D[bname + "residualYdY"]->fill( _fitY[ipl]  , _fitY[ipl]    - hitY[jhit]  );
 }
 
#endif

  
        }

        nFittedTracks++;

      }  // end of track filling 



    }
    // End of loop over track possibilities

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_firstChi2HistoName]))->fill(log10(chi2min));
#endif

    if(nFittedTracks==0) {
      if(streamlog_level(DEBUG5)){
	streamlog_out ( DEBUG5 ) << "No track fulfilling search criteria found ! " << endl;
      }

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nTrackHistoName]))->fill(0);
#endif

      // before throwing the exception I should clean up the
      // memory...
      delete fittrackvec;
      delete fitpointvec;
      delete corrpointvec;
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



    // Fitted track analysis
    // =====================


    // Take only first track (best Chi2), if multiple track finding disabled

    if(!_searchMultipleTracks) nFittedTracks=1;

    int itrk=0;
  
    int nStoredTracks=0;

    for(fitIterator=fittedChi2.begin(); fitIterator!=fittedChi2.end(); fitIterator++)
       {
       // Chi2 map is sorted by definition, we start from the best fit

       double choiceChi2 = fitIterator->first;
       int ifit = fitIterator->second;

       int nChoiceFired = fittedFired[ifit];
       double penalty = fittedPenalty[ifit];


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
       if(itrk==0 && _searchMultipleTracks) 
         {
     (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_bestChi2HistoName]))->fill(log10(choiceChi2));
     (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nBestHistoName]))->fill(nChoiceFired);
	 }
#endif


       // Skip tracks having too many ambiguous hits
       // Check better fits only!

       bool isAmbiguous = false;

       for(fitIterator2=fittedChi2.begin(); fitIterator2!=fitIterator; fitIterator2++)
	 {
           int ifit2 = fitIterator2->second;
           int nAmbiguous=0;

           for(int ipl=0;ipl<_nTelPlanes;ipl++)  
	       if(fittedHits[_nTelPlanes*ifit+ipl]==fittedHits[_nTelPlanes*ifit2+ipl]
                  && 
                  fittedHits[_nTelPlanes*ifit+ipl]>=0) 
                                                  nAmbiguous++;

           if( 
	      (  _allowAmbiguousHits && nAmbiguous > _maximumAmbiguousHits )
                  ||
	      ( !_allowAmbiguousHits && nAmbiguous > 0 )

              ) isAmbiguous=true;
	 }


       if(isAmbiguous)continue;

       // Count accepted track for contributing hits

       for(int ipl=0;ipl<_nTelPlanes;ipl++) 
         {
	   int jhit = fittedHits[_nTelPlanes*ifit+ipl];

	   if(jhit >= 0)
                  hitFits[jhit]++;
	 }

       if(streamlog_level(DEBUG5)){
         streamlog_out ( DEBUG5 )  << "Track reconstructed from " << nChoiceFired << " hits: " << endl;


        // print out hits contributing to the fit

        for(int ipl=0;ipl<_nTelPlanes;ipl++)  
          {
          int jhit = fittedHits[_nTelPlanes*ifit+ipl];
          if(jhit >= 0)
              streamlog_out ( DEBUG5 ) << "Hit " << jhit
                                                 << "   X = " << hitX[jhit]
                                                 << "   Y = " << hitY[jhit]
                                                 << "   Z = " << hitZ[jhit] 
                                                 << " (plane" << hitPlane[jhit] << ")" << endl;
          }


        streamlog_out ( DEBUG5) << " Fitted positions in telescope planes:" << endl;

        for(int ipl=0;ipl<_nTelPlanes;ipl++) {
          streamlog_out ( DEBUG5) << "  X = " << fittedX[_nTelPlanes*ifit+ipl] 
				  << " +/- " << fittedEx[_nTelPlanes*ifit+ipl]
				  << "  Y = " << fittedY[_nTelPlanes*ifit+ipl] 
				  << " +/- " << fittedEy[_nTelPlanes*ifit+ipl]
				  << "  at Z = " << _planePosition[ipl] << endl;
        }

        streamlog_out ( DEBUG5) << " Fit chi2 = " << choiceChi2 << " including penalties of " << penalty << endl;

       }

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
      // Fill Chi2 histograms

      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_linChi2HistoName]))->fill(choiceChi2);

      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_logChi2HistoName]))->fill(log10(choiceChi2));

      if(_allowMissingHits && nChoiceFired==_nActivePlanes) {
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_fullChi2HistoName]))->fill(log10(choiceChi2));
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

      fittrack->setChi2(choiceChi2);  // Chi2 of the fit (including penalties)
      fittrack->setNdf(nChoiceFired); // Number of planes fired (!)

      // Fitted position at DUT is stored as Reference Point
      // (see below in loop over telescope planes)
      // If no DUT present: use position in the first plane !

//      fittrack->setIsReferencePointPCA(false);
      float refpoint[3];

      // Track points fitted in each plane are stored as track hits

      for(int ipl=0;ipl<_nTelPlanes;ipl++)  
      {
        TrackerHitImpl * fitpoint = new TrackerHitImpl;

	// set sensorID
	fitHitEncoder["sensorID"] =  _planeID[ipl];

	// set the local/global and "fittedhit" bit flag properties for the hit
	fitHitEncoder["properties"] = kHitInGlobalCoord+kFittedHit;

	// store values
	fitHitEncoder.setCellID( fitpoint );

        // fitted position in a plane
        double pos[3];

        pos[0]=fittedX[_nTelPlanes*ifit+ipl];
        pos[1]=fittedY[_nTelPlanes*ifit+ipl];
        pos[2]=_planePosition[ipl];

        getFastTrackImpactPoint(pos[0], pos[1], pos[2], fittrack, event);
        fitpoint->setPosition(pos);

        // Covariance matrix of the position
        // (stored as lower triangle matrix, i.e.  cov(xx),cov(y,x),cov(y,y) ).

        float cov[TRKHITNCOVMATRIX];

        cov[0]=fittedEx[_nTelPlanes*ifit+ipl]*fittedEx[_nTelPlanes*ifit+ipl];
        cov[1]=0.;
        cov[2]=fittedEy[_nTelPlanes*ifit+ipl]*fittedEy[_nTelPlanes*ifit+ipl];
        cov[3]=cov[4]=0.;
        cov[5]=_planeThickness[ipl]*_planeThickness[ipl]/12.;

        fitpoint->setCovMatrix(cov);

        // store fit point

        fitpointvec->push_back(fitpoint);

        //   add fitted point to track
        //   
        if(_OutputHitsInTrack)  fittrack->addHit(fitpoint);

        // add measured point to track (if found)
        // 
        if(_InputHitsInTrack && _isActive[ipl])  
        {
          int jhit = fittedHits[_nTelPlanes*ifit+ipl];

          if(jhit >= 0)  
          {
            TrackerHitImpl    * meshit  = dynamic_cast<TrackerHitImpl*>( col->getElementAt(jhit) ) ;
            SimTrackerHitImpl * simhit  = dynamic_cast<SimTrackerHitImpl*>( col->getElementAt(jhit) ) ;
            TrackerHitImpl    * corrhit = new TrackerHitImpl;

            //
            // Copy input hit data
            //
            if( meshit != 0 ) 
            {
              fitpoint->setType(meshit->getTime());
              corrhit->setType(meshit->getType());
              corrhit->setTime(meshit->getTime());
              corrhit->setEDep(meshit->getEDep());
              corrhit->rawHits()=meshit->getRawHits();
	      corrhit->setCellID0(meshit->getCellID0());
	      corrhit->setCellID1(meshit->getCellID1());
            }
            else if( simhit != 0 )
            {
              fitpoint->setType(simhit->getTime());
              corrhit->setType(0);
              corrhit->setTime(simhit->getTime());
              corrhit->setEDep(simhit->getEDep());
	      corrhit->setCellID0(simhit->getCellID0());
	      corrhit->setCellID1(simhit->getCellID1());
            } else {
	      // set sensorID
	      corrHitEncoder["sensorID"] =  _planeID[ipl];
	    
	      // set the local/global and "fittedhit" bit flag properties for the hit
	      corrHitEncoder["properties"] = 0; // init
	      corrHitEncoder["properties"] = kHitInGlobalCoord;
	    
	      // store values
	      corrHitEncoder.setCellID( corrhit );
	    }
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

        if(ipl==_iDUT || (_iDUT<0 && ipl==0)) 
        {
          for(int iref=0;iref<3;iref++) 
          {
            refpoint[iref]=pos[iref];
          }
        }
      }
      // Store track reference point.
      fittrack->setReferencePoint(refpoint);
      fittrackvec->addElement(fittrack);
      // increment the total track counter

      _noOfTracks++;

      nStoredTracks++;

      itrk++;  // loop count
      }

  // End of track finding procedure

  if(nStoredTracks > 0 ) 
  {
    event->addCollection(fittrackvec,_outputTrackColName);
    event->addCollection(fitpointvec,_outputHitColName);
    event->addCollection(corrpointvec,_correctedHitColName);
  } 
  else 
  {
    delete fittrackvec;
    delete fitpointvec;
    delete corrpointvec;
  }


  if ( fittrackvec->size() == 0 ) 
  {
    ++_noOfEventWOTrack;
  }

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  // Number of reconstructed tracks
  // 
  (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nTrackHistoName]))->fill(nStoredTracks);
#endif

  // Hit ambiguity
  // 
  if (_allowAmbiguousHits) 
  {
    for(int ihit=0; ihit< nHit ; ihit++) 
    {
      if(_isActive[hitPlane[ihit]]) 
      {

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_hitAmbiguityHistoName]))->fill(hitFits[ihit]);
#endif
      }
    }
  }


  // Clear all working arrays

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



void EUTelTestFitter::check( LCEvent * /* evt */  ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void EUTelTestFitter::end(){

  if(streamlog_level(DEBUG5))
  {
    for(int ipl=0;ipl<_nTelPlanes;ipl++)  
      {
	stringstream iden;
	iden << "pl" << _planeID[ipl] << "_";
	string bname = iden.str();
	
	streamlog_out( DEBUG5 ) << "X: ["<< ipl << ":" << _planeID[ipl] <<"]" << 
	  _aidaHistoMap1D[bname + "residualX"]->allEntries()<< " " <<
	  _aidaHistoMap1D[bname + "residualX"]->mean()*1000. << " " <<
	  _aidaHistoMap1D[bname + "residualX"]->rms()*1000. << " " <<
	  _aidaHistoMap1D[bname + "residualY"]->allEntries()<< " " <<
	  _aidaHistoMap1D[bname + "residualY"]->mean()*1000. << " " <<
	  _aidaHistoMap1D[bname + "residualY"]->rms()*1000. << " " << endl;
      }
  }

  // Print the summary
  streamlog_out( MESSAGE5 ) << "Total number of processed events:     " << setw(10) << setiosflags(ios::right) << _nEvt << resetiosflags(ios::right) << endl
                           << "Total number of events w/o input hit: " << setw(10) << setiosflags(ios::right) << _noOfEventWOInputHit 
                           << resetiosflags(ios::right) << endl
                           << "Total number of events w/o track:     " << setw(10) << setiosflags(ios::right) << _noOfEventWOTrack
                           << resetiosflags(ios::right) << endl
                           << "Total number of reconstructed tracks  " << setw(10) << setiosflags(ios::right) << _noOfTracks << resetiosflags(ios::right)
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
  delete [] _planeScatAngle ;
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

  delete [] _nominalFitArrayX ;
  delete [] _nominalErrorX ;

  delete [] _nominalFitArrayY ;
  delete [] _nominalErrorY ;
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
    streamlog_out( ERROR5 ) << "I/O problem with " << _histoInfoFileName << "\n"
                           << "Continuing without histogram manager"    << endl;
    isHistoManagerAvailable = false;
  } catch ( ParseException& e ) {
    streamlog_out( ERROR5 ) << e.what() << "\n"
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
          streamlog_out ( DEBUG5 )  << (* histoInfo ) << endl;
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
          streamlog_out ( DEBUG5 )   << (* histoInfo ) << endl;
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

// plot plane by plane:
   for(int iz=0; iz < _nTelPlanes ; iz++) {
//plane id by      _planeID[iz]  
    stringstream iden;
    iden << "pl" << _planeID[iz] << "_";
    string bname = iden.str();
 

    int   limitXN  = 100;
    int   limitYN  = 100;
    int   limitZN  = 100;
    float limitY   = 10.0;
    float limitX   = 10.0; 
    float limitYr  = 0.1;
    float limitXr  = 0.1; 
    //float limitZ   = 50.0; 
    float limitZr  = 50.0; 
   //Resids 
    _aidaHistoMap1D[bname + "fitX"] =  AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "fitX", limitXN , -limitX, limitX);
    _aidaHistoMap1D[bname + "fitY"] =  AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "fitY", limitYN , -limitY, limitY);
    _aidaHistoMap1D[bname + "hitX"] =  AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "hitX", limitXN , -limitX, limitX); 
    _aidaHistoMap1D[bname + "hitY"] =  AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "hitY", limitYN , -limitY, limitY);
    _aidaHistoMap1D[bname + "residualX"] =  AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "residualX", limitXN , -limitXr, limitXr);
    _aidaHistoMap1D[bname + "residualY"] =  AIDAProcessor::histogramFactory(this)->createHistogram1D( bname + "residualY", limitYN , -limitYr, limitYr);
    //Resids 2D
    _aidaHistoMap2D[bname + "residualXdX"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( bname + "residualXdX", limitXN, -limitX, limitX, limitXN, -limitXr, limitXr);
    _aidaHistoMap2D[bname + "residualYdX"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( bname + "residualYdX", limitXN, -limitX, limitX, limitYN, -limitYr, limitYr);
    _aidaHistoMap2D[bname + "residualXdY"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( bname + "residualXdY", limitYN, -limitY, limitY, limitXN, -limitXr, limitXr);
    _aidaHistoMap2D[bname + "residualYdY"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( bname + "residualYdY", limitYN, -limitY, limitY, limitYN, -limitYr, limitYr);

    _aidaHistoMap2D[bname + "residualdZvsX"]        =  AIDAProcessor::histogramFactory(this)->createHistogram2D( bname + "residualdZvsX",limitXN, -limitX, limitX, limitZN ,-limitZr, limitZr);
    _aidaHistoMap2D[bname + "residualdZvsY"]        =  AIDAProcessor::histogramFactory(this)->createHistogram2D( bname + "residualdZvsY",limitYN, -limitY, limitY, limitZN ,-limitZr, limitZr);
    _aidaHistoMap2D[bname + "residualmeasZvsmeasX"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( bname + "residualmeasZvsmeasX",limitXN , -limitX, limitX, limitZN ,-limitZr, limitZr);
    _aidaHistoMap2D[bname + "residualmeasZvsmeasY"] =  AIDAProcessor::histogramFactory(this)->createHistogram2D( bname + "residualmeasZvsmeasY",limitYN , -limitY, limitY, limitZN ,-limitZr, limitZr);
    _aidaHistoMap2D[bname + "residualfitZvsmeasX"]  =  AIDAProcessor::histogramFactory(this)->createHistogram2D( bname + "residualfitZvsmeasX",limitXN , -limitX, limitX, limitZN ,-limitZr, limitZr);
    _aidaHistoMap2D[bname + "residualfitZvsmeasY"]  =  AIDAProcessor::histogramFactory(this)->createHistogram2D( bname + "residualfitZvsmeasY",limitYN , -limitY, limitY, limitZN ,-limitZr, limitZr);
  }

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
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
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
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
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
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
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
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
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
  if(_allowAmbiguousHits )
    {
      string ambigTitle = "Number of tracks containing given hit (ambiguity)";
      AIDA::IHistogram1D * AmbigHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( _hitAmbiguityHistoName.c_str(),hitNBin,hitMin,hitMax);
      AmbigHisto->setTitle(ambigTitle.c_str());
      _aidaHistoMap.insert(make_pair(_hitAmbiguityHistoName, AmbigHisto));
    }


  // List all booked histogram - check of histogram map filling
  streamlog_out ( DEBUG5 ) <<  _aidaHistoMap.size() << " histograms booked" << endl;

  map<string, AIDA::IBaseHistogram *>::iterator mapIter;
  for(mapIter = _aidaHistoMap.begin(); mapIter != _aidaHistoMap.end() ; mapIter++ )
    streamlog_out ( DEBUG5 ) <<  mapIter->first << " : " <<  (mapIter->second)->title()  << endl;

  streamlog_out ( DEBUG5 ) << "Histogram booking completed \n\n" << endl;

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

  int status = DoAnalFit(_fitX,_fitEx,_beamSlopeX);

  if(status)return -1. ;

  status = DoAnalFit(_fitY,_fitEy,_beamSlopeY);

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

  int status = DoAnalFit(_fitX,_fitEx,_beamSlopeX);

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
      _fitEx[ipl]=_nominalErrorX[ipl];
      _fitEy[ipl]=_nominalErrorY[ipl];

      _fitX[ipl]=0. ;
      _fitY[ipl]=0. ;

      for(int jpl=0; jpl<_nTelPlanes;jpl++)
        {
          if(_planeEx[jpl]>0.)
            _fitX[ipl]+=_nominalFitArrayX[ipl+jpl*_nTelPlanes]*_planeX[jpl]/_planeEx[jpl]/_planeEx[jpl];
          if(_planeEy[jpl]>0.)
            _fitY[ipl]+=_nominalFitArrayY[ipl+jpl*_nTelPlanes]*_planeY[jpl]/_planeEy[jpl]/_planeEy[jpl];
        }

      // Correction for beam slope

    if(_useBeamConstraint && _beamSlopeX!=0.)
      {
	_fitX[ipl]-=_nominalFitArrayX[ipl]*_beamSlopeX*_planeDist[0]*_planeScat[0];
	_fitX[ipl]+=_nominalFitArrayX[ipl+_nTelPlanes]*_beamSlopeX*_planeDist[0]*_planeScat[0];
      }

    if(_useBeamConstraint && _beamSlopeY!=0.)
      {
	_fitY[ipl]-=_nominalFitArrayY[ipl]*_beamSlopeY*_planeDist[0]*_planeScat[0];
	_fitY[ipl]+=_nominalFitArrayY[ipl+_nTelPlanes]*_beamSlopeY*_planeDist[0]*_planeScat[0];
      }

    }

  double chi2=GetFitChi2();

  return chi2 ;
}


int EUTelTestFitter::DoAnalFit(double * pos, double *err, double slope)
{
  for(int ipl=0; ipl<_nTelPlanes;ipl++)
    {
      if(_isActive[ipl] && err[ipl]>0)
        err[ipl]=1./err[ipl]/err[ipl] ;
      else
        err[ipl] = 0. ;

      pos[ipl]*=err[ipl];
    }

  // To take into account beam tilt

  if(_useBeamConstraint && slope!=0.)
    {
      pos[0] -= slope*_planeDist[0]*_planeScat[0];
      pos[1] += slope*_planeDist[0]*_planeScat[0];
    }


  for(int ipl=0; ipl<_nTelPlanes;ipl++)
  {
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
  // Beam slope taken into account.

  if(_useBeamConstraint)
    {
      double dth;

      // Use small angle approximation: atan(x) = x
      // Should be:
      //    dth=atan((_fitX[1]-_fitX[0])*_planeDist[0]) ;
      dth=(_fitX[1]-_fitX[0])*_planeDist[0]  -  _beamSlopeX;
      chi2 += _planeScat[0] * dth * dth;

      dth=(_fitY[1]-_fitY[0])*_planeDist[0]  -  _beamSlopeY ;
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

  for(i=0;i<n;i++)ipiv[i]=0;

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

      if(ipiv[icol]>1){
	// first clean up then bail out
	delete[] ipiv;
	delete[] indxr;
	delete[] indxc;
        return 1;
      }

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

      if(alfa[n*icol+icol]==0.){
	// first clean up then bail out
	delete[] ipiv;
	delete[] indxr;
	delete[] indxc;
        return 1;}

      help=alfa[n*icol+icol];
      pivinv=1./help;
      alfa[n*icol+icol]=1.;
      for(j=0;j<n;j++) alfa[n*icol+j]*=pivinv;

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

void EUTelTestFitter::getFastTrackImpactPoint(double & x, double & y, double & z, Track * /* tr */, LCEvent * /* ev */) {

    // given maps:
    //            _siPlaneCenter.insert( make_pair( sensorID,  _planeCenter ));
    //            _siPlaneNormal.insert( make_pair( sensorID,  _planeNormal ));


    // this is a simplified version of getTrackImpactPoint
    // 
    double pos[3]={x,y,z};
    int sensorID = guessSensorID( pos );

    //
    // assume no slope for a track
    //
    double slopeX = 0.;
    double slopeY = 0.;

    double offsetX = x;
    double offsetY = y;   

    // this is a normal vector to the plane
	TVectorD NormalVector(3);
	NormalVector[0] = _siPlaneNormal[sensorID][0] ;
	NormalVector[1] = _siPlaneNormal[sensorID][1] ;
	NormalVector[2] = _siPlaneNormal[sensorID][2] ;

	// this is a vector to the point in the plane
	TVectorD    r0Vector(3);
	r0Vector[0] = _siPlaneCenter[sensorID][0] ;
	r0Vector[1] = _siPlaneCenter[sensorID][1] ;
	r0Vector[2] = _siPlaneCenter[sensorID][2] ;


	// now have to solve the equation
	TVectorD	trackImpact(3);
	TMatrixD	equationMatrix(3,3);
	TVectorD	b(3);

        equationMatrix(0, 0) = NormalVector (0);
	equationMatrix(0, 1) = NormalVector (1);
	equationMatrix(0, 2) = NormalVector (2);
	equationMatrix(1, 0) = 1;
	equationMatrix(1, 1) = 0;
	equationMatrix(1, 2) = (-1)*slopeX;
	equationMatrix(2, 0) = 0;
	equationMatrix(2, 1) = 1;
	equationMatrix(2, 2) = (-1)*slopeY;

	b[0] = r0Vector(0) * NormalVector (0) + r0Vector(1) * NormalVector (1) + r0Vector(2) * NormalVector (2);
	b[1] = offsetX;
	b[2] = offsetY;

	trackImpact = equationMatrix.Invert() * b;

	x = trackImpact(0);
	y = trackImpact(1);
	z = trackImpact(2);

}


void EUTelTestFitter::getTrackImpactPoint(double & x, double & y, double & z, Track * tr, LCEvent * ev) {

  // what should be here:
  // get the center of sensor (X,Y,Z) and normal vector (a,b,c)
  // based on Gear information and all available alignment collections 
  // at the moment Millepede "alignment" and offset "preAlignment" collections
  //
  // solve linear algebra
  //  intersection of sensor Plane and Track
  //
  //  input: 
  //        Track   *tr
  //        LCEvent *ev
  //  output:
  //        double x,y,z
  //
  //
    
  int iplane = guessSensorID( tr );

  // angles only!
  TVector3 _GEAR_(0.,0.,0.);
  TVector3 _PreAlignment_(0.,0.,0.);
  TVector3 _Alignment_(0.,0.,0.);

  _GEAR_.SetX(   _siPlanesLayerLayout->getLayerRotationXY(iplane) ); // alfa
  _GEAR_.SetY(   _siPlanesLayerLayout->getLayerRotationZX(iplane) ); // beta
  _GEAR_.SetZ(   _siPlanesLayerLayout->getLayerRotationZY(iplane) ); // gamma

  // code taken from EUTelAPIXHistograms
  // small fixes to make it compile,
  // definitely not the way to write code
  // 
  double _zDUT           = 0.;
  double _zDUTneighbour  = 0.;
  int    _indexDUT       = 0 ;
  int    _manualDUTid    = 0 ;

  std::vector<EVENT::TrackerHit*>  trackhits = tr->getTrackerHits();
  int nHit =   trackhits.size();


  // first, get the track impact points at the DUT and the layer closest to it

  // coordinates of the track impact point at the closest to DUT layer
  double x1=0, y1=0, z1=0;
  // coordinates of the track impact point at the DUT layer
  double x2=0, y2=0, z2=0;

  // for sanity check
  bool	foundHitDUT = false;
  bool	foundHitNeighbour = false;

  // setup cellIdDecoder to decode the hit properties
  CellIDDecoder<TrackerHit>  hitCellDecoder(EUTELESCOPE::HITENCODING);
   
  double	dist;
  for(int ihit=0; ihit< nHit ; ihit++)
    {
      TrackerHit * meshit = trackhits.at(ihit);
 
      // Look at fitted hits only!
      if (meshit != 0){
	if ( (hitCellDecoder(meshit)["properties"] & kFittedHit) == 0 ){
	  continue;
	}
      }

      // Hit position
      //
      double pos[3]={0.,0.,0.};
      if( meshit != 0 )
	{
	  const double *pos0 = meshit->getPosition();
	  pos[0] = pos0[0];
	  pos[1] = pos0[1];
	  pos[2] = pos0[2];
	} 
      // look for a hit at DUT
      dist =  pos[2] - _zDUT ;
      if (dist * dist < 1) {
	if (foundHitDUT) {
	  cout << "hit at DUT layer already found! Terminating" << endl;
	  abort();
	}
	x2 = pos[0];
	y2 = pos[1];
	z2 = pos[2];
	foundHitDUT = true;
      }
      // look for a hit at DUT's neighbour
      dist =  pos[2] - _zDUTneighbour ;
      if (dist * dist < 1) {
	if (foundHitNeighbour) {
	  cout << "hit at DUT neigbour layer already found! Terminating" << endl;
	  abort();
	}
	x1 = pos[0];
	y1 = pos[1];
	z1 = pos[2];
      }
    
    }

  if (! (foundHitDUT && foundHitNeighbour) ) {
    cout << "Was not possible to find hits at dut and next-to-dut layers. Terminating." << endl;
    abort();
  }

  // now proceed to calculation of the intersection point
  // of this track and the rotated DUT layer

  // for the formulas see logbook 20/01/2011
  // track parametrization: X = offsetX + slopeX * Z, Y = offsetY + slopeY * Z,

  double	slopeX = (x2 - x1) / (z2-z1);
  double	offsetX = x1 - z1 * slopeX;

  double	slopeY = (y2 - y1) / (z2-z1);
  double	offsetY = y1 - z1 * slopeY;


  (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap["TrackSlopeX"]))->fill(slopeX);
  (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap["TrackOffsetX"]))->fill(offsetX);
  (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap["TrackSlopeY"]))->fill(slopeY);
  (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap["TrackOffsetY"]))->fill(offsetY);


  // this is a normal vector to the plane
  TVectorD	NormalVector(3);
  NormalVector(0) = 0.;
  NormalVector(1) = 0.;
  NormalVector(2) = 1.;

  double	z_sensor = _siPlanesLayerLayout->getSensitivePositionZ(_indexDUT) +
    0.5 * ( _siPlanesLayerLayout->getLayerThickness( _indexDUT ) + _siPlanesLayerLayout->getSensitiveThickness( _indexDUT ) );

  // this is a vector to the point in the plane
  TVectorD		r0Vector(3);
  r0Vector(0)=0;
  r0Vector(1)=0;
  r0Vector(2)=z_sensor;

  // this is a vector, pointing to the (0, 0, z_sensor )
  // it's needed to translate the coordinate system to
  // and to perform the rotation around that point, not the (0, 0, 0)
  // as it's done in the alignment/applyAlignment steps.
  TVectorD		auxVector(3);
  auxVector(0) = 0;
  auxVector(1) = 0;
  auxVector(2) = z_sensor;

  for ( unsigned i = 0; i< _alignmentCollectionNames.size(); i++) {

    LCCollectionVec * alignmentCollectionVec     = dynamic_cast < LCCollectionVec * > (ev->getCollection(_alignmentCollectionNames[i]));
    // next, find the alignment constant corresponding to the DUT
    EUTelAlignmentConstant * c=NULL;
    for ( size_t iPos = 0; iPos < alignmentCollectionVec->size(); ++iPos ) {

      c = static_cast< EUTelAlignmentConstant * > ( alignmentCollectionVec->getElementAt( iPos ) );
      if (c -> getSensorID() == _manualDUTid ) break;	// this means we found the alignment constant corresponding
      // to the DUT; the pointer to it is now stored in c and can
      // be furhter used
    }
    if ( c == NULL ) {
      cout << "Was not possible to found alignment constant, terminating" << endl;
      abort();
    }

    //------------------------------------------------------------------------------------
    TMatrixD	RotationMatrix(3,3);
    RotationMatrix(0,0) = 1;
    RotationMatrix(0,1) = c -> getGamma();
    RotationMatrix(0,2) = c -> getBeta();
    RotationMatrix(1,0) = (-1) * c -> getGamma();
    RotationMatrix(1,1) = 1;
    RotationMatrix(1,2) = c -> getAlpha();
    RotationMatrix(2,0) = (-1) * c -> getBeta();
    RotationMatrix(2,1) = (-1) * c -> getAlpha();
    RotationMatrix(2,2) = 1;

    // transform the normal vector (only the rotation)
    NormalVector = RotationMatrix * NormalVector;

    // transform the vector to the plane point (rotation+translations)
    r0Vector = RotationMatrix *  (r0Vector - auxVector) + auxVector;
    r0Vector(0) -= c -> getXOffset();
    r0Vector(1) -= c -> getYOffset();
    r0Vector(2) -= c -> getZOffset();

    //------------------------------------------------------------------------------------

    //deltaX -= c->getXOffset();
  }

  // now have to solve the equation
  TVectorD	trackImpact(3);
  TMatrixD	equationMatrix(3,3);
  TVectorD	b(3);

  equationMatrix(0, 0) = NormalVector (0);
  equationMatrix(0, 1) = NormalVector (1);
  equationMatrix(0, 2) = NormalVector (2);
  equationMatrix(1, 0) = 1;
  equationMatrix(1, 1) = 0;
  equationMatrix(1, 2) = (-1)*slopeX;
  equationMatrix(2, 0) = 0;
  equationMatrix(2, 1) = 1;
  equationMatrix(2, 2) = (-1)*slopeY;

  b(0) = r0Vector(0) * NormalVector (0) + r0Vector(1) * NormalVector (1) + r0Vector(2) * NormalVector (2);
  b(1) = offsetX;
  b(2) = offsetY;

  trackImpact = equationMatrix.Invert() * b;

  /*
  // very very naive approach
  // finally calculate track impact point
  z = ( z2 /tan((-1)*_beta) + deltaX - offsetX) / (slopeX + 1./tan((-1)*_beta));
  x = slopeX * z + offsetX;
  // in this simple approximation, no correction for y
  y = y2;*/

  x = trackImpact(0);
  y = trackImpact(1);
  z = trackImpact(2);

}



int EUTelTestFitter::guessSensorID( Track * /*track*/ ) {

  streamlog_out( ERROR5 ) << " guessSensorID(Track *) called but no longer implemented! Please use guessSensorID( double & x, double & y, double & z) instead. " << endl;

  return -1;
}


int EUTelTestFitter::guessSensorID( double & /*x*/, double & /*y*/, double & z) 
{

  int sensorID = -1;
  double minDistance =  numeric_limits< double >::max() ;

  for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); ++iPlane ) 
  {
      int _id = _siPlanesLayerLayout->getID(iPlane);
      double distance = std::abs( z - _siPlaneCenter[ _id  ][2] );
      if ( distance < minDistance ) 
      {
          minDistance = distance;
          sensorID = _id;
      }
  }

  if ( minDistance > 10 /* mm */ ) 
  {
    // advice the user that the guessing wasn't successful 
    streamlog_out( WARNING3 ) << "A hit was found " << minDistance << " mm far from the nearest plane\n"
      "Please check the consistency of the data with the GEAR file " << endl;
        
  }

  return sensorID;
}


int EUTelTestFitter::guessSensorID( double * hit ) 
{

  int sensorID = -1;
  double minDistance =  numeric_limits< double >::max() ;

  if( _referenceHitVec == 0 || _useReferenceHitCollection == false) {
      // use z information of planes instead of reference vector
      for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); ++iPlane ) {
	double distance = std::abs( hit[2] - _siPlanesLayerLayout->getLayerPositionZ(iPlane) );
	if ( distance < minDistance ) {
	  minDistance = distance;
	  sensorID = _siPlanesLayerLayout->getID( iPlane );
	}
      }
      if ( minDistance > 30  ) {
	// advice the user that the guessing wasn't successful 
	streamlog_out( WARNING3 ) << "A hit was found " << minDistance << " mm far from the nearest plane\n"
	  "Please check the consistency of the data with the GEAR file: hitPosition[2]=" << hit[2] <<       endl;
      }
    
      return sensorID;
    }

  for(int ii = 0 ; ii <  _referenceHitVec->getNumberOfElements(); ii++) {
        EUTelReferenceHit* refhit = static_cast< EUTelReferenceHit*> ( _referenceHitVec->getElementAt(ii) ) ;
        
        TVector3 hit3d( hit[0], hit[1], hit[2] );
        TVector3 hitInPlane( refhit->getXOffset(), refhit->getYOffset(), refhit->getZOffset());
        TVector3 norm2Plane( refhit->getAlpha(), refhit->getBeta(), refhit->getGamma() );
 
        double distance = abs( norm2Plane.Dot(hit3d-hitInPlane) );
        if ( distance < minDistance ) 
        {
           minDistance = distance;
           sensorID = refhit->getSensorID();
        }    

      }

  return sensorID;
}


#endif


