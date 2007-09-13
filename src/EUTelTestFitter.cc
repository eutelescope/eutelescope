// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-

// Author: A.F.Zarnecki, University of Warsaw <mailto:zarnecki@fuw.edu.pl>
// Version: $Id: EUTelTestFitter.cc,v 1.14 2007-09-13 17:09:43 zarnecki Exp $
// Date 2007.06.04

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope inlcudes
#include "EUTelTestFitter.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelHistogramManager.h"
#include "EUTelExceptions.h"

#ifdef MARLIN_USE_AIDA
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


using namespace std;
using namespace lcio ;
using namespace marlin ;
using namespace eutelescope;

// definition of static members mainly used to name histograms
#ifdef MARLIN_USE_AIDA
std::string EUTelTestFitter::_logChi2HistoName   = "logChi2";
std::string EUTelTestFitter::_firstChi2HistoName  = "firstChi2";
std::string EUTelTestFitter::_bestChi2HistoName  = "bestChi2";
std::string EUTelTestFitter::_fullChi2HistoName  = "fullChi2";
std::string EUTelTestFitter::_nTrackHistoName    = "nTrack";
std::string EUTelTestFitter::_nHitHistoName      = "nHit";
std::string EUTelTestFitter::_nBestHistoName     = "nBest";
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
			      _debugCount,  static_cast < int > (100));

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

  registerOptionalParameter ("UseBeamConstraint",
			     "Flag for using beam direction constraint in the fit",
                             _useBeamConstraint,  static_cast < bool > (false));

  registerOptionalParameter ("BeamSpread",
                             "Assumed angular spread of the beam [rad]",
                             _beamSpread,  static_cast < double > (0.1));

  registerOptionalParameter ("SearchMultipleTracks",
			     "Flag for searching multiple tracks in events with multiple hits",
                             _searchMultipleTracks,  static_cast < bool > (false));

}


void EUTelTestFitter::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  // check if Marlin was built with GEAR support or not
#ifndef USE_GEAR

  message<ERROR> ( "Marlin was not built with GEAR support." );
  message<ERROR> ( "You need to install GEAR and recompile Marlin with -DUSE_GEAR before continue.");
  
  // I'm thinking if this is the case of throwing an exception or
  // not. This is a really error and not something that can
  // exceptionally happens. Still not sure what to do
  exit(-1);

#else

  // check if the GEAR manager pointer is not null!
  if ( Global::GEAR == 0x0 ) {
    message<ERROR> ( "The GearMgr is not available, for an unknown reason." );
    exit(-1);
  }

  // Read geometry information from GEAR

  message<MESSAGE> ( log() << "Reading telescope geometry description from GEAR ") ;

  _siPlanesParameters  = const_cast<gear::SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

#endif

// Test output

 if( _SkipLayerIDs.size() )
      message<MESSAGE> ( log() <<  _SkipLayerIDs.size() << " layers should be skipped") ;

// Active planes only:
//  _nTelPlanes = _siPlanesParameters->getSiPlanesNumber();

// Take all layers defined in GEAR geometry
  _nTelPlanes = _siPlanesLayerLayout->getNLayers();

// Check for DUT

   if( _siPlanesParameters->getSiPlanesType()==_siPlanesParameters->TelescopeWithDUT )
        {
      _iDUT = _nTelPlanes ;
      _nTelPlanes++;
        }
   else
       _iDUT = -1 ;

// Read position in Z (for sorting), skip layers if requested

  _planeSort = new int[_nTelPlanes];
  _planePosition   = new double[_nTelPlanes];

  int nSkip=0;

  for(int ipl=0; ipl <  _siPlanesLayerLayout->getNLayers(); ipl++)
    {
    _planePosition[ipl-nSkip]=_siPlanesLayerLayout->getLayerPositionZ(ipl);
    _planeSort[ipl-nSkip]=ipl;

// Check if not on "skip list"

    int _plID = _siPlanesLayerLayout->getID(ipl);

    for(int spl=0; spl< (int)_SkipLayerIDs.size() ; spl++)
      if ( _SkipLayerIDs.at(spl) == _plID)
       {
       message<MESSAGE> ( log() <<  "Skipping layer ID " << _plID 
       << " at Z = " << _siPlanesLayerLayout->getLayerPositionZ(ipl) ) ;
       nSkip++;
       break;
       }

    }


  _nTelPlanes-=nSkip;
  
  if(_iDUT>0)
      {
      _iDUT-=nSkip;
      _planePosition[_iDUT]=_siPlanesLayerLayout->getDUTPositionZ();
      _planeSort[_iDUT]=_iDUT;
      }

 // Binary sorting

   bool sorted;
   do{
   sorted=false;
   for(int iz=0; iz<_nTelPlanes-1 ; iz++)
      if(_planePosition[iz]>_planePosition[iz+1])
         {
         double _posZ = _planePosition[iz];
         _planePosition[iz] = _planePosition[iz+1];
         _planePosition[iz+1] = _posZ;

         int _idZ = _planeSort[iz];
         _planeSort[iz] = _planeSort[iz+1];
         _planeSort[iz+1] = _idZ;

         sorted=true;
         }

     }while(sorted);

// Book local geometry arrays

  _planeID         = new int[_nTelPlanes];
  _planeShiftX     = new double[_nTelPlanes];
  _planeShiftY     = new double[_nTelPlanes];
  _planeRotZ       = new double[_nTelPlanes];
  _planeThickness  = new double[_nTelPlanes];
  _planeX0         = new double[_nTelPlanes];
  _planeResolution = new double[_nTelPlanes];
  _isActive        = new bool[_nTelPlanes];
  _nActivePlanes = 0 ;

// Fill arrays with parameters of layer, sorted in Z

  for(int iz=0; iz < _nTelPlanes ; iz++)
    {
      int ipl=_planeSort[iz];

      int iActive;
      double resolution;

// All dimensions are assumed to be in mm !!!

      if(ipl != _iDUT )
         {
      _planeID[iz]=_siPlanesLayerLayout->getID(ipl);
      _planeThickness[iz]=_siPlanesLayerLayout->getLayerThickness(ipl);
      _planeX0[iz]=_siPlanesLayerLayout->getLayerRadLength(ipl);
      resolution = _siPlanesLayerLayout->getSensitiveResolution(ipl);
        }
      else
         {
      _planeID[iz]=_siPlanesLayerLayout->getDUTID();
      _planeThickness[iz]=_siPlanesLayerLayout->getDUTThickness();
      _planeX0[iz]=_siPlanesLayerLayout->getDUTRadLength();
       resolution = _siPlanesLayerLayout->getDUTSensitiveResolution();
        }

      iActive = (resolution > 0);
 
      if(iActive && (ipl != _iDUT || _useDUT ))
	{
	  _isActive[iz] = true ;
	  _planeResolution[iz]=resolution;
	  _nActivePlanes++ ;
	}
      else
	{
          _isActive[iz] = false ;
          _planeResolution[iz]=0.;
	}

// No alignment corrections in GEAR file
// Look in input options

    _planeShiftX[iz]=0.;
    _planeShiftY[iz]=0.;
    _planeRotZ[iz]=0.;

     for(int apl=0; apl< (int)_AlignLayerIDs.size() ; apl++)
      if ( _AlignLayerIDs.at(apl) == _planeID[iz])
       {
       _planeShiftX[iz]=_AlignLayerShiftX.at(apl);
       _planeShiftY[iz]=_AlignLayerShiftY.at(apl);
       // Rotation can be skipped: check size
       if(apl < (int)_AlignLayerRotZ.size())
              _planeRotZ[iz]=_AlignLayerRotZ.at(apl);
       break;
       }

   }

  // Get new DUT position (after sorting)

  for(int iz=0;iz< _nTelPlanes ; iz++)
    if(_planeSort[iz]==_iDUT)
       {
        _iDUT=iz;
        break;
       }


  // Print out geometry information

  message<MESSAGE> ( log() << "Telescope configuration with " << _nTelPlanes << " planes" );


  for(int ipl=0; ipl < _nTelPlanes; ipl++)
    {
      stringstream ss ; 
      if(ipl == _iDUT)
	ss << "D.U.T.  plane" ;
      else
	if(_isActive[ipl])
	  ss << "Active  plane" ;
	else
	  ss << "Passive plane" ; 
      
      ss << "  ID = " << _planeID[ipl]
         << "  at Z [mm] = " << _planePosition[ipl] 
	 << " dZ [um] = " << _planeThickness[ipl]*1000. ;
      
      if(_isActive[ipl])
        {
	ss << "  Res [um] = " << _planeResolution[ipl]*1000. ;
      
        if(_planeShiftX[ipl] !=0. || _planeShiftY[ipl] !=0. ) 
	  ss << "\n  alignment corrections:" 
             <<  " dX [mm] = " << _planeShiftX[ipl] 
              << " dY [mm] = " << _planeShiftY[ipl] ;
        if(_planeRotZ[ipl] !=0.) 
           ss << " RotZ [rad] = " << _planeRotZ[ipl] ;
        }

      message<MESSAGE> ( log() << ss.str() );
    }

  message<MESSAGE> ( log() << "Total of " << _nActivePlanes << " active sensor planes " );
  
    // Allocate arrays for track fitting

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

  // Planes have to be ordered in position along the beam line !
  // This is not checked !!!

  for(int ipl=0; ipl<_nTelPlanes ; ipl++)
    {
      if(ipl>0)
	_planeDist[ipl-1]=1./(_planePosition[ipl] - _planePosition[ipl-1]) ;

      _planeScat[ipl]= 0.0136/_eBeam * sqrt(_planeThickness[ipl]/_planeX0[ipl])
	* (1.+0.038*std::log(_planeThickness[ipl]/_planeX0[ipl])) ;

      if(ipl==0 && _useBeamConstraint)
	_planeScat[ipl]= 1./(_planeScat[ipl]*_planeScat[ipl]+ _beamSpread*_beamSpread) ; 
      else
	_planeScat[ipl]= 1./(_planeScat[ipl] * _planeScat[ipl]) ; 

      _fitX[ipl] =_fitY[ipl] = 0. ;
      _nominalError[ipl]= _planeResolution[ipl];
    }

  // Fit with nominal parameters

  int status = DoAnalFit(_fitX,_nominalError);

  if(status)
    cerr << "\n Fit with nominal geometry failed !?!" << endl ;

  // Store fit matrix

  for(int imx=0; imx<arrayDim; imx++)
    _nominalFitArray[imx] = _fitArray[imx];

  message<MESSAGE> ( log() << "Expected position resolutions [um]: " ) ;
  stringstream ss;
  for(int ipl=0; ipl<_nTelPlanes ; ipl++) 
    ss << _nominalError[ipl]*1000. << "  " ;
  message<MESSAGE> ( log() << ss.str() );


// Book histograms

#ifdef MARLIN_USE_AIDA
    bookHistos();
#endif

}

void EUTelTestFitter::processRunHeader( LCRunHeader* runHeader) { 

  auto_ptr<EUTelRunHeaderImpl> eutelHeader( new EUTelRunHeaderImpl ( runHeader ) );
  eutelHeader->addProcessor( type() );
  
  _nRun++ ;

  // Decode and print out Run Header information - just a check

  int runNr = runHeader->getRunNumber();
  
  message<MESSAGE> ( log() << "Processing run header " << _nRun 
		     << ", run nr " << runNr );

  const std::string detectorName = runHeader->getDetectorName();
  const std::string detectorDescription = runHeader->getDescription();
  const std::vector<std::string> * subDets = runHeader->getActiveSubdetectors();

  message<MESSAGE> ( log() << detectorName << " : " << detectorDescription ) ;

  int nDet = subDets->size();

  message<MESSAGE> ( log() << nDet << " subdetectors defined :" );
  stringstream ss;
  for(int idet=0;idet<nDet;idet++)  message<MESSAGE> (log()  << idet+1 << " : " << subDets->at(idet) );


} 

void EUTelTestFitter::processEvent( LCEvent * event ) { 
  
  EUTelEventImpl * euEvent = static_cast<EUTelEventImpl*> ( event );
  if ( euEvent->getEventType() == kEORE ) {
    message<DEBUG> ( "EORE found: nothing else to do." );
    return;
  }

  bool debug = ( _debugCount>0 && _nEvt%_debugCount == 0);

  _nEvt ++ ;
  int evtNr = event->getEventNumber();


  if(debug)message<DEBUG> ( log() << "Processing record " << _nEvt << " == event " << evtNr );

  LCCollection* col;
  try {
    col = event->getCollection( _inputColName ) ;
  } catch (lcio::DataNotAvailableException& e) {
    message<ERROR> ( log() << "Not able to get collection " 
		     << _inputColName 
		     << "\nfrom event " << event->getEventNumber()
		     << " in run " << event->getRunNumber()  );
   (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nTrackHistoName]))->fill(0);
   throw SkipEventException(this);
  }
    

  // Copy hits to local table
  // Assign hits to sensor planes
  // =============================


  int nHit = col->getNumberOfElements()  ;

  if(debug)message<DEBUG> ( log() << "Total of " << nHit << " tracker hits in input collection " );

  if(nHit + _allowMissingHits < _nActivePlanes) {
    if(debug)message<DEBUG> ( log() << "Not enough hits to perform the fit, exiting... " );
    (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nTrackHistoName]))->fill(0);
    throw SkipEventException(this);
  }

  double * hitX  = new double[nHit];
  double * hitEx = new double[nHit];
  double * hitY  = new double[nHit];
  double * hitEy = new double[nHit];
  double * hitZ  = new double[nHit];
  int    * hitPlane = new int[nHit];

  int    * nPlaneHits   = new int[_nTelPlanes];
  int    * nPlaneChoice = new int[_nTelPlanes];
  IntVec * planeHitID   = new IntVec[_nTelPlanes];

  double * bestX  = new double[_nTelPlanes];
  double * bestEx = new double[_nTelPlanes];
  double * bestY  = new double[_nTelPlanes];
  double * bestEy = new double[_nTelPlanes];

  // Loop over hits

  int nGoodHit = 0;

  for(int ihit=0; ihit< nHit ; ihit++)
    {
      TrackerHit * meshit = dynamic_cast<TrackerHit*>( col->getElementAt(ihit) ) ;

      // Hit position

      const double * pos = meshit->getPosition();

      hitZ[ihit] = pos[2];


      // We have to find Plane ID of the hit
      // by looking at the Z position

      double distMin = 1.;
      hitPlane[ihit] = -1 ;

      for(int ipl=0;ipl<_nTelPlanes;ipl++)
	{
	  double dist =  hitZ[ihit] - _planePosition[ipl] ;

	  if(dist*dist < distMin*distMin)
	    {
	      hitPlane[ihit]=ipl;
	      distMin=dist;
	    }
	}

      // Ignore hits not matched to any plane

      if(hitPlane[ihit]<0) {

	message<ERROR> ( log() << "Reconstructed hit outside sensor planes z [mm] = "  << hitZ[ihit] );

      }

      // Ignore hit, if plane not declared as active (i.e. not used in the fit)

      if(! _isActive[hitPlane[ihit]])
	continue ;

      // Ignore hit also, if maximum number of hits already matched to this plane

      if(_maxPlaneHits>0 && 
         _maxPlaneHits-planeHitID[hitPlane[ihit]].size()<=0) 
	continue ;

      // Hit will be used: correct X and Y position for plane alignment

      //      hitX[ihit] = pos[0];
      //      hitY[ihit] = pos[1];

      hitX[ihit] = pos[0]*cos(_planeRotZ[hitPlane[ihit]])
                 + pos[1]*sin(_planeRotZ[hitPlane[ihit]])
                 + _planeShiftX[hitPlane[ihit]];

      hitY[ihit] = pos[1]*cos(_planeRotZ[hitPlane[ihit]])
                 - pos[0]*sin(_planeRotZ[hitPlane[ihit]])
                 + _planeShiftY[hitPlane[ihit]];

      // Add hit to hit list for given plane - to be used in track selection

      planeHitID[hitPlane[ihit]].push_back(ihit);
      nGoodHit++;

      // Position uncertainty. Use nominal resolution if not properly defined

      const EVENT::FloatVec cov = meshit->getCovMatrix();

      if(cov.at(0)>0.)
	hitEx[ihit]=sqrt(cov.at(0));
      else
	hitEx[ihit]=_planeResolution[hitPlane[ihit]];

      if(cov.at(2)>0.)
	hitEy[ihit]=sqrt(cov.at(2));
      else
	hitEy[ihit]=_planeResolution[hitPlane[ihit]];


      if(debug)message<DEBUG> ( log() << "Hit " << ihit
		       << "   X = " << hitX[ihit] << " +/- " << hitEx[ihit]  
		       << "   Y = " << hitY[ihit] << " +/- " << hitEy[ihit]  
		       << "   Z = " << hitZ[ihit] << " (plane " << hitPlane[ihit] << ")" );
    }


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

  do{

    // Count planes active in this event and number of fit possibilities

    int nFiredPlanes = 0;
    int nChoice = 1;
    ibest=-1;

    for(int ipl=0;ipl<_nTelPlanes;ipl++)
      {
	nPlaneHits[ipl] = planeHitID[ipl].size() ;

	if(nPlaneHits[ipl]>0)
	  {
	    nFiredPlanes++;

	    if(_allowSkipHits>0)
	      nPlaneChoice[ipl]=nPlaneHits[ipl]+1;
	    else
	      nPlaneChoice[ipl]=nPlaneHits[ipl];
	  }
	else
	  nPlaneChoice[ipl]=1;

	nChoice*=nPlaneChoice[ipl];

	if( _isActive[ipl] && firstTrack && debug )
	  {
	    stringstream ss;
	    ss << "Plane " << ipl << "  " << nPlaneHits[ipl] << " hit(s), hit IDs :";
	    for( int ihit=0; ihit < (int) planeHitID[ipl].size() ; ihit ++)
	      ss << planeHitID[ipl].at(ihit) << " " ;
	    message<DEBUG> ( log() << ss.str() );
	  }
      }
 
    // Check if fit can be done

    if(nFiredPlanes + _allowMissingHits < _nActivePlanes)
      {
        if( firstTrack ) {
	  if(debug)message<DEBUG> ( log() << "Not enough planes hit to perform the fit " );

          (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nTrackHistoName]))->fill(0);

	  // before throwing the exception I should clean up the
	  // memory...
	  delete [] bestEy;
	  delete [] bestY;
	  delete [] bestEx;
	  delete [] bestX;
	  delete [] planeHitID;
	  delete [] nPlaneChoice;
	  delete [] nPlaneHits;
	  delete [] hitPlane;
	  delete [] hitZ;
	  delete [] hitEy;
	  delete [] hitY;
	  delete [] hitEx;
	  delete [] hitX;
	  throw SkipEventException(this);
	}
      }


    if( firstTrack && debug ) 
           message<DEBUG> ( log() << nFiredPlanes << " active sensor planes hit, checking "
				       << nChoice << " fit possibilities " );

    // Check all track possibilities

    double chi2min  = 1E+100;
    double chi2best = _chi2Max;
    double bestPenalty = 0;
    int nBestFired = 0;

    // Loop over fit possibilities

    for(int ichoice=0; ichoice<nChoice; ichoice++)
      {
	int modchoice=ichoice;

	int nChoiceFired=0;
	double choiceChi2=-1.;

	// Fill position and error arrays for this hit configuration

	for(int ipl=0;ipl<_nTelPlanes;ipl++)
	  {
	    _planeX[ipl]=_planeY[ipl]=_planeEx[ipl]=_planeEy[ipl]=0.;

	    if(_isActive[ipl])
	      {
		int ihit=modchoice%nPlaneChoice[ipl];
		modchoice/=nPlaneChoice[ipl];

		if(ihit<nPlaneHits[ipl])
		  {
		    int jhit = planeHitID[ipl].at(ihit);
		    _planeX[ipl]=hitX[jhit];
		    _planeY[ipl]=hitY[jhit];
		    _planeEx[ipl]=(_useNominalResolution)?_planeResolution[ipl]:hitEx[jhit];
		    _planeEy[ipl]=(_useNominalResolution)?_planeResolution[ipl]:hitEy[jhit];
		    nChoiceFired++;
		  }
	      }
	  }

	// Check number of selected hits

	if(nChoiceFired + _allowMissingHits < _nActivePlanes ||
	   nChoiceFired + _allowSkipHits    < nFiredPlanes )
	  continue;

	// Select fit method
	// "Nominal" fit only if all active planes used

	if(_useNominalResolution && (nChoiceFired == _nActivePlanes))
          choiceChi2 = NominalFit();
	else
	  if(_useNominalResolution)
	    choiceChi2 = SingleFit();
	  else
	    choiceChi2 = MatrixFit();

	// Fit failed ?

	if(choiceChi2 < 0.)continue ;

	// Penalty for missing or skiped hits

	double penalty = (_nActivePlanes-nFiredPlanes)*_missingHitPenalty
	  +   (nFiredPlanes-nChoiceFired)*_skipHitPenalty ;


	if(choiceChi2+penalty<chi2min)
            chi2min=choiceChi2+penalty;   

	// Best fit ?

	if(choiceChi2+penalty<chi2best)
          {
	    chi2best=choiceChi2+penalty;
	    bestPenalty=penalty;
	    ibest=ichoice;
	    nBestFired=nChoiceFired;
	    for(int ipl=0;ipl<_nTelPlanes;ipl++)
	      {
		bestX[ipl]=_fitX[ipl];
		bestY[ipl]=_fitY[ipl];
		bestEx[ipl]=_fitEx[ipl];
		bestEy[ipl]=_fitEy[ipl];
	      }
	  }
      }
    // End of loop over track possibilities

    if(firstTrack) 
      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_firstChi2HistoName]))->fill(log10(chi2min));


    if(ibest<0 && firstTrack) {
      if(debug)message<DEBUG> ( log() << "No track fulfilling search criteria found ! " );


      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nTrackHistoName]))->fill(0);

      // before throwing the exception I should clean up the
      // memory...
      delete [] bestEy;
      delete [] bestY;
      delete [] bestEx;
      delete [] bestX;
      delete [] planeHitID;
      delete [] nPlaneChoice;
      delete [] nPlaneHits;
      delete [] hitPlane;
      delete [] hitZ;
      delete [] hitEy;
      delete [] hitY;
      delete [] hitEx;
      delete [] hitX;	
      throw SkipEventException(this);
    }
    
    
    if(ibest>=0)
      {

      if(debug)
        {
        message<DEBUG> ( log() << "Track reconstructed from " << nBestFired << " hits: " );


    // print out hits contributing to the fit

      
	int modchoice=ibest;

	for(int ipl=0;ipl<_nTelPlanes;ipl++)
	  {
	    if(_isActive[ipl])
	      {
		  int ihit=modchoice%nPlaneChoice[ipl];
		  modchoice/=nPlaneChoice[ipl];

		  if(ihit<nPlaneHits[ipl])
		    {
		    int jhit = planeHitID[ipl].at(ihit);
                    message<DEBUG> ( log() << "Hit " << jhit
		       << "   X = " << hitX[jhit] 
		       << "   Y = " << hitY[jhit] 
		       << "   Z = " << hitZ[jhit] << " (plane" << hitPlane[jhit] << ")" );
		    }
	      }
	  }


	message<DEBUG> (log() << " Fitted positions in telescope planes:");

	for(int ipl=0;ipl<_nTelPlanes;ipl++)
	  message<DEBUG> ( log() << "  X = " << bestX[ipl] << " +/- " << bestEx[ipl] 
			   << "  Y = " << bestY[ipl] << " +/- " << bestEy[ipl] 
			   << "  at Z = " << _planePosition[ipl] ) ;
	

	message<DEBUG> ( log() << " Fit chi2 = " << chi2best << " including penalties of " << bestPenalty );

        }

     // Fill Chi2 histograms

	(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_logChi2HistoName]))->fill(log10(chi2best));

	if(_searchMultipleTracks && firstTrack)
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_bestChi2HistoName]))->fill(log10(chi2best));

        if(_allowMissingHits && nBestFired==_nActivePlanes)
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_fullChi2HistoName]))->fill(log10(chi2best));

     // Fill hit histograms

	(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nHitHistoName]))->fill(nBestFired);

	if(_searchMultipleTracks && firstTrack)
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nBestHistoName]))->fill(nBestFired);



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

	int modchoice=ibest;

	for(int ipl=0;ipl<_nTelPlanes;ipl++)
	  {
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

            if(_OutputHitsInTrack)
              fittrack->addHit(fitpoint);


	    //		       << "   X = " << hitX[ihit] << " +/- " << hitEx[ihit]  
	    //       << "   Y = " << hitY[ihit] << " +/- " <<
	    //       hitEy[ihit]  

            // add measured point to track (if found)

 	    if(_InputHitsInTrack && _isActive[ipl])
	      {
		 int ihit=modchoice%nPlaneChoice[ipl];
		 modchoice/=nPlaneChoice[ipl];

		 if(ihit<nPlaneHits[ipl])
		    {
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

	    if(ipl==_iDUT || (_iDUT<0 && ipl==0))
	      for(int iref=0;iref<3;iref++)
		refpoint[iref]=pos[iref];
	  }
	
	// Store track reference point.
	fittrack->setReferencePoint(refpoint);
	
	fittrackvec->addElement(fittrack);
	nFittedTracks++;
      }
    
    // If multiple tracks allowed: remove hits from fitted track from the list
    
    if(ibest>=0 && _searchMultipleTracks)
      {
	int modchoice=ibest;
	
	for(int ipl=0;ipl<_nTelPlanes;ipl++)
	  if(_isActive[ipl])
	    {
	      int ihit=modchoice%nPlaneChoice[ipl];
	      modchoice/=nPlaneChoice[ipl];
	      
	      if(ihit<nPlaneHits[ipl])
		{
                  planeHitID[ipl].erase(planeHitID[ipl].begin()+ihit);
                  nGoodHit--;
		}
	    }
      }
    
    firstTrack=false;
  }
  while(_searchMultipleTracks && ibest>=0 && nGoodHit + _allowMissingHits >= _nActivePlanes);
    
  // End of track loop

  if(nFittedTracks > 0 )
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

  // Number of reconstructed tracks

  (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_nTrackHistoName]))->fill(nFittedTracks);

  // Clear all working arrays
  
  delete [] bestEy;
  delete [] bestY;
  delete [] bestEx;
  delete [] bestX;
  delete [] planeHitID;
  delete [] nPlaneChoice;
  delete [] nPlaneHits;
  delete [] hitPlane;
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
  // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
  // 	    << std::endl ;


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

#ifdef MARLIN_USE_AIDA

  message<MESSAGE> ( log() << "Booking histograms " );


  message<MESSAGE> ( log() << "Histogram information searched in " << _histoInfoFileName);

  auto_ptr<EUTelHistogramManager> histoMgr( new EUTelHistogramManager( _histoInfoFileName ));
  EUTelHistogramInfo    * histoInfo;
  bool                    isHistoManagerAvailable;

  try {
    isHistoManagerAvailable = histoMgr->init();
  } catch ( ios::failure& e) {
    message<ERROR> ( log() << "I/O problem with " << _histoInfoFileName << "\n"
		     << "Continuing without histogram manager"    );
    isHistoManagerAvailable = false;
  } catch ( ParseException& e ) {
    message<ERROR> ( log() << e.what() << "\n"
		     << "Continuing without histogram manager" );
    isHistoManagerAvailable = false;
  }

 // log(Chi2) distribution for all accepted tracks

    int    chi2NBin  = 100;
    double chi2Min   = -2.;
    double chi2Max   = 8.;
    string chi2Title = "log(Chi2) distribution";


    if ( isHistoManagerAvailable ) 
      {
      histoInfo = histoMgr->getHistogramInfo(_logChi2HistoName);
      if ( histoInfo ) 
         {
	 message<DEBUG> ( log() << (* histoInfo ) );
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
	 message<DEBUG> ( log() << (* histoInfo ) );
	 trkNBin = histoInfo->_xBin;
	 trkMin  = histoInfo->_xMin;
	 trkMax  = histoInfo->_xMax;
	 if ( histoInfo->_title != "" ) trkTitle = histoInfo->_title;
         }
      }


   AIDA::IHistogram1D * nTrackHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( _nTrackHistoName.c_str(),trkNBin,trkMin,trkMax);

     nTrackHisto->setTitle(trkTitle.c_str());

    _aidaHistoMap.insert(make_pair(_nTrackHistoName, nTrackHisto));



 // Number of hits per track

    int    hitNBin  = 11;
    double hitMin   = -0.5;
    double hitMax   = 10.5;
    string hitTitle = "Number of hits used to reconstruct the track";


    if ( isHistoManagerAvailable ) 
      {
      histoInfo = histoMgr->getHistogramInfo(_nHitHistoName);
      if ( histoInfo ) 
         {
	 message<DEBUG> ( log() << (* histoInfo ) );
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


// List all booked histogram - check of histogram map filling

  message<MESSAGE> ( log() <<  _aidaHistoMap.size() << " histograms booked");  


  map<string, AIDA::IBaseHistogram *>::iterator mapIter;
  for(mapIter = _aidaHistoMap.begin(); mapIter != _aidaHistoMap.end() ; mapIter++ )
     message<DEBUG> ( log() <<  mapIter->first << " : " <<  (mapIter->second)->title() ) ;

  message<DEBUG> ( log() << "Histogram booking completed \n\n");  
#else
  message<MESSAGE> ( log() << "No histogram produced because Marlin doesn't use AIDA" );
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
	  _fitArray[imx] += _planeScat[1]*_planeDist[1]*_planeDist[1] ;

	if(ipl+jpl==1 && _useBeamConstraint) 
	  _fitArray[imx] -= _planeScat[1]*_planeDist[1]*_planeDist[1] ;


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
