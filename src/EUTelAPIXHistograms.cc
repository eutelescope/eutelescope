
// Author: A.F.Zarnecki, University of Warsaw <mailto:zarnecki@fuw.edu.pl>
// @version: $Id$

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifdef USE_GEAR

// eutelescope includes ".h"
#include "EUTelAPIXHistograms.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelHistogramManager.h"
#include "EUTelExceptions.h"
#include "EUTelAPIXSparsePixel.h"
#include "EUTelSparseDataImpl.h"
#include "EUTelAPIXSparseClusterImpl.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes ".h"
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/IProfile2D.h>
#include <AIDA/ITree.h>
#endif

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// lcio includes ".h"
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>

// ROOT includes ".h"
#include <TVectorD.h>
#include <TMatrixD.h>

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

using namespace std;
using namespace lcio ;
using namespace marlin ;
using namespace eutelescope;

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelAPIXHistograms::_MeasuredXHistoName  = "measuredX";
std::string EUTelAPIXHistograms::_MeasuredYHistoName  = "measuredY";
std::string EUTelAPIXHistograms::_MeasuredXYHistoName = "measuredXY";

std::string EUTelAPIXHistograms::_MatchedXHistoName  = "matchedX";
std::string EUTelAPIXHistograms::_MatchedYHistoName  = "matchedY";
std::string EUTelAPIXHistograms::_MatchedXYHistoName = "matchedXY";

std::string EUTelAPIXHistograms::_UnMatchedXHistoName  = "unmatchedX";
std::string EUTelAPIXHistograms::_UnMatchedYHistoName  = "unmatchedY";
std::string EUTelAPIXHistograms::_UnMatchedXYHistoName = "unmatchedXY";

std::string EUTelAPIXHistograms::_FittedXHistoName  = "fittedX";
std::string EUTelAPIXHistograms::_FittedYHistoName  = "fittedY";
std::string EUTelAPIXHistograms::_FittedXYHistoName = "fittedXY";

std::string EUTelAPIXHistograms::_EfficiencyXHistoName  = "DUTeffiX";
std::string EUTelAPIXHistograms::_EfficiencyYHistoName  = "DUTeffiY";
std::string EUTelAPIXHistograms::_EfficiencyXYHistoName = "DUTeffiXY";

std::string EUTelAPIXHistograms::_BgEfficiencyXHistoName  = "BGeffiX";
std::string EUTelAPIXHistograms::_BgEfficiencyYHistoName  = "BGeffiY";
std::string EUTelAPIXHistograms::_BgEfficiencyXYHistoName = "BGeffiXY";

std::string EUTelAPIXHistograms::_NoiseXHistoName  = "DUTnoiseX";
std::string EUTelAPIXHistograms::_NoiseYHistoName  = "DUTnoiseY";
std::string EUTelAPIXHistograms::_NoiseXYHistoName = "DUTnoiseXY";

std::string EUTelAPIXHistograms::_ShiftXHistoName       = "DUTshiftX";
std::string EUTelAPIXHistograms::_ShiftYHistoName       = "DUTshiftY";
std::string EUTelAPIXHistograms::_ShiftXYHistoName      = "DUTshiftXY";

std::string EUTelAPIXHistograms::_BgShiftXHistoName       = "BGshiftX";
std::string EUTelAPIXHistograms::_BgShiftYHistoName       = "BGshiftY";
std::string EUTelAPIXHistograms::_BgShiftXYHistoName      = "BGshiftXY";

std::string EUTelAPIXHistograms::_ShiftXvsYHistoName    = "DUTshiftXvsY";
std::string EUTelAPIXHistograms::_ShiftYvsXHistoName    = "DUTshiftYvsX";
std::string EUTelAPIXHistograms::_ShiftXvsY2DHistoName    = "DUTshiftXvsY2D";
std::string EUTelAPIXHistograms::_ShiftYvsX2DHistoName    = "DUTshiftYvsX2D";

std::string EUTelAPIXHistograms::_ShiftYvsYHistoName    = "DUTshiftYvsY";
std::string EUTelAPIXHistograms::_ShiftXvsXHistoName    = "DUTshiftXvsX";
std::string EUTelAPIXHistograms::_ShiftXvsX2DHistoName    = "DUTshiftXvsX2D";
std::string EUTelAPIXHistograms::_ShiftYvsY2DHistoName    = "DUTshiftYvsY2D";

std::string EUTelAPIXHistograms::_EtaXHistoName      = "EtaX";
std::string EUTelAPIXHistograms::_EtaYHistoName      = "EtaY";
std::string EUTelAPIXHistograms::_EtaX2DHistoName    = "EtaX2D";
std::string EUTelAPIXHistograms::_EtaY2DHistoName    = "EtaY2D";
std::string EUTelAPIXHistograms::_EtaX3DHistoName    = "EtaX3D";
std::string EUTelAPIXHistograms::_EtaY3DHistoName    = "EtaY3D";


// some APIX histograms
// libov@mail.desy.de 05 August 2010
std::string	EUTelAPIXHistograms::_totAPIXmatchedHistoName				=	"APIX_ToT_matched";
std::string	EUTelAPIXHistograms::_totAPIXunmatchedHistoName			=	"APIX_ToT_unmatched";

std::string	EUTelAPIXHistograms::_lv1APIXmatchedHistoName				=	"APIX_lv1_matched";
std::string	EUTelAPIXHistograms::_lv1APIXunmatchedHistoName			=	"APIX_lv1_unmatched";

std::string	EUTelAPIXHistograms::_hitsInMatchedClusterHistoName		=	"APIX_hitsInCluster_matched";
std::string	EUTelAPIXHistograms::_hitsInUnmatchedClusterHistoName	=	"APIX_hitsInCluster_unmatched";

std::string	EUTelAPIXHistograms::_maxDifflv1MatchedHistoName			=	"APIX_maxlv1Diff_matched";
std::string	EUTelAPIXHistograms::_maxDifflv1UnmatchedHistoName		=	"APIX_maxlv1Diff_unmatched";

// Implementing DUT plots in its local FoR, 18 August
std::string EUTelAPIXHistograms::_EfficiencyXLOCALHistoName 			=	"APIX_DUTeffiXLOCAL";
std::string EUTelAPIXHistograms::_EfficiencyYLOCALHistoName 			=	"APIX_DUTeffiYLOCAL";
std::string EUTelAPIXHistograms::_EfficiencyXYLOCALHistoName			=	"APIX_DUTeffiXYLOCAL";

// 24 August 2010 libov@mail.desy.de
std::string EUTelAPIXHistograms::_MeasuredXLOCALHistoName				=	"APIX_measuredXLOCAL";
std::string EUTelAPIXHistograms::_MeasuredYLOCALHistoName				=	"APIX_measuredYLOCAL";
std::string EUTelAPIXHistograms::_MeasuredXYLOCALHistoName			=	"APIX_measuredXYLOCAL";

std::string EUTelAPIXHistograms::_FittedXLOCALHistoName				=	"APIX_fittedXLOCAL";
std::string EUTelAPIXHistograms::_FittedYLOCALHistoName				=	"APIX_fittedYLOCAL";
std::string EUTelAPIXHistograms::_FittedXYLOCALHistoName				=	"APIX_fittedXYLOCAL";

std::string	EUTelAPIXHistograms::_MatchedClusterSizeXHistoName		=	"APIX_ClusterSizeX_matched";
std::string	EUTelAPIXHistograms::_MatchedClusterSizeYHistoName		=	"APIX_ClusterSizeY_matched";

std::string	EUTelAPIXHistograms::_UnmatchedClusterSizeXHistoName	=	"APIX_ClusterSizeX_unmatched";
std::string	EUTelAPIXHistograms::_UnmatchedClusterSizeYHistoName	=	"APIX_ClusterSizeY_unmatched";

std::string	EUTelAPIXHistograms::_ChargeSharingProbXHistoName		=	"APIX_chargeSharingProbX";
std::string	EUTelAPIXHistograms::_ChargeSharingProbYHistoName		=	"APIX_chargeSharingProbY";
std::string	EUTelAPIXHistograms::_ChargeSharingProbXYHistoName		=	"APIX_chargeSharingProbXY";

// 06 September
std::string	EUTelAPIXHistograms::_totSinglePixelClustersXHistoName	=	"APIX_ToT_SinglePixelClustersX";
std::string	EUTelAPIXHistograms::_totSinglePixelClustersYHistoName	=	"APIX_ToT_SinglePixelClustersY";
std::string	EUTelAPIXHistograms::_totSinglePixelClustersXYHistoName	=	"APIX_ToT_SinglePixelClustersXY";

std::string	EUTelAPIXHistograms::_totAllClustersXYHistoName	=	"APIX_ToT_AllClustersXY";

std::string	EUTelAPIXHistograms::_MatchedHitsHistoName	=	"APIX_MatchedPerEvent";
std::string	EUTelAPIXHistograms::_MatchedHitsVSEventHistoName	=	"APIX_MatchedPerEventVSEvent";

std::string	EUTelAPIXHistograms::_NumberOfFittedTracksHistoName	=	"APIX_NumberOfFittedTracks";

std::string	EUTelAPIXHistograms::_NumberOfMatchedReferencesHistoName = "APIX_NumberOfMatchedReferences";
std::string	EUTelAPIXHistograms::_Lv1OfFirstMatchedReferenceHistoName = "APIX_Lv1OfFirstMatchedReference";



#endif


EUTelAPIXHistograms::EUTelAPIXHistograms() : Processor("EUTelAPIXHistograms") {

  // modify processor description
  _description = "Analysis of DUT performance based on the analytic track fit results" ;


  // register steering parameters:
  //       name, description, class-variable, default value

  // input collection first:

  registerInputCollection( LCIO::TRACK,
                           "InputTrackCollectionName" ,
                           "Name of the input Track collection"  ,
                           _inputTrackColName ,
                           std::string("testfittracks") ) ;

  registerInputCollection( LCIO::TRACKERHIT,
                           "InputHitCollectionName" ,
                           "Name of the input DUT hit collection"  ,
                           _inputHitColName ,
                           std::string("hit") ) ;

  // other processor parameters:

  registerProcessorParameter ("UseManualDUT",
                              "Flag for manual DUT selection",
                              _useManualDUT,  static_cast < bool > (false));

  registerProcessorParameter ("ManualDUTid",
                              "Id of telescope layer which should be used as DUT",
                              _manualDUTid,  static_cast < int > (0));

  registerProcessorParameter ("DistMax",
                              "Maximum allowed distance between fit and matched DUT hit",
                              _distMax,  static_cast < double > (0.1));


  registerProcessorParameter ("DUTpitchX",
                              "DUT sensor pitch in X",
                              _pitchX,  static_cast < double > (0.03));


  registerProcessorParameter ("DUTpitchY",
                              "DUT sensor pitch in Y",
                              _pitchY,  static_cast < double > (0.03));


  std::vector<float > initAlign;
  initAlign.push_back(0.);
  initAlign.push_back(0.);
  initAlign.push_back(0.);

  registerProcessorParameter ("DUTalignment",
                              "Alignment corrections for DUT: shift (in mm) in X, Y and rotation around Z",
                              _DUTalign, initAlign);

  registerProcessorParameter("HistoInfoFileName",
                             "Name of the histogram information file",
                             _histoInfoFileName, string( "histoinfo.xml" ) );

  registerProcessorParameter ("DebugEventCount",
                              "Print out every DebugEnevtCount event",
                              _debugCount,  static_cast < int > (100));

// 17 August 2010 libov@mail.desy.de

registerProcessorParameter ("onlyIntimeTracks",
                              "For the efficiency determination, use only tracks with have matched hit in reference plane",
                           _onlyIntimeTracks,  static_cast < bool > (false));

registerProcessorParameter ("referencePlaneID",
                              "Reference plane ID,",
                           _referencePlaneID,  static_cast < int > (13));

registerProcessorParameter ("DistMaxReference",
                              "Max distance from the track position to hit in a Reference Plane ID to be considered as matched",
                           _distMaxReference,  static_cast < double > (0.2));


StringVec	_alignmentCollectionSuffixExamples;
_alignmentCollectionSuffixExamples.push_back("iter3");

registerProcessorParameter ("alignmentCollectionNames",
                              "List of alignment collections that were applied to the DUT",
                           _alignmentCollectionSuffixes, _alignmentCollectionSuffixExamples);

registerProcessorParameter("Beta","Rotation Angle around Y axis",
                            _beta, static_cast< double > ( 0.00 ) );


}


void EUTelAPIXHistograms::init() {

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  // check if Marlin was built with GEAR support or not
#ifndef USE_GEAR

  message<ERROR5> ( "Marlin was not built with GEAR support." );
  message<ERROR5> ( "You need to install GEAR and recompile Marlin with -DUSE_GEAR before continue.");

  // I'm thinking if this is the case of throwing an exception or
  // not. This is a really error and not something that can
  // exceptionally happens. Still not sure what to do
  exit(-1);

#else

  // check if the GEAR manager pointer is not null!
  if ( Global::GEAR == 0x0 ) {
    message<ERROR5> ( "The GearMgr is not available, for an unknown reason." );
    exit(-1);
  }

  // Read geometry information from GEAR

  message<MESSAGE5> ( log() << "Reading telescope geometry description from GEAR ") ;

  _siPlanesParameters  = const_cast<gear::SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

#endif

  if(_useManualDUT)
    {
      bool _manualOK=false;

      for(int ipl=0; ipl <  _siPlanesLayerLayout->getNLayers(); ipl++)
        if(_siPlanesLayerLayout->getID(ipl)==_manualDUTid)
          {
            _iDUT=_manualDUTid;
            _zDUT=_siPlanesLayerLayout->getLayerPositionZ(ipl);
				_indexDUT = ipl;
            _manualOK=true;
          }
		cout << "_indexDUT = " << _indexDUT << endl;
      if(!_manualOK)
        {
          message<ERROR5> ( log() << "Manual DUT flag set, layer not found ID = "
                           << _manualDUTid
                           << "\n Program will terminate! Correct geometry description!");
          exit(-1);
        }
    }
  else
    if( _siPlanesParameters->getSiPlanesType()==_siPlanesParameters->TelescopeWithDUT )
      {
        _iDUT=_siPlanesLayerLayout->getDUTID();
        _zDUT=_siPlanesLayerLayout->getDUTPositionZ();
      }
    else
      {
        message<ERROR5> ( log() << "DUT analysis initialized, but no DUT found in GEAR \n"
                         << "Program will terminate! Correct geometry description!");
        exit(-1);
      }


// Print out geometry information

  message<MESSAGE5> ( log() << "D.U.T. plane  ID = " << _iDUT
                     << "  at Z [mm] = " << _zDUT );

	// 18 August
	_xPitch       = _siPlanesLayerLayout->getSensitivePitchX(_indexDUT);    // mm
	_yPitch       = _siPlanesLayerLayout->getSensitivePitchY(_indexDUT);    // mm

	a = _siPlanesLayerLayout->getSensitiveRotation1(_indexDUT); // was -1 ;
	b = _siPlanesLayerLayout->getSensitiveRotation2(_indexDUT); // was  0 ;
	c = _siPlanesLayerLayout->getSensitiveRotation3(_indexDUT); // was  0 ;
	d = _siPlanesLayerLayout->getSensitiveRotation4(_indexDUT); // was -1 ;
	// now rotate back -- add more description here or put a link -- see paper logbook 18 August 2010 libov@mail.desy.de
	double const_factor = static_cast< double >(1./((b*c)-(a*d)));
	_rot00 = const_factor * d * (-1);
	_rot01 = const_factor * b;
	_rot10 = const_factor * c;
	_rot11 = const_factor * a* (-1);
	
    // 21 January 2011
	// get an index of a layer, closest to the DUT (upstream)
	double	zMin = 1000000.;
	_indexDUTneighbour=-1;
	_zDUTneighbour=-1;
	for(int ipl=0; ipl <  _siPlanesLayerLayout->getNLayers(); ipl++) {

		double	z = _siPlanesLayerLayout->getLayerPositionZ(ipl);
		if  ( ( z < _zDUT ) && ( (_zDUT - z) < zMin) ) {
			_indexDUTneighbour = ipl;
			zMin = _zDUT - z;
		}
	}
	_zDUTneighbour = _siPlanesLayerLayout->getLayerPositionZ(_indexDUTneighbour);

	if ( (_indexDUTneighbour == -1) || (_zDUTneighbour == -1) ) {
		cout << "was not able to determine layer next to the 	DUT. Terminating!" << endl;
		abort();
	}

	cout << "The layer closest to the DUT: index= "<<_indexDUTneighbour <<  ", z= " << _zDUTneighbour << endl;


// Book histograms

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  bookHistos();
#endif

	_noHitYet = true;
	_eventsWithNoHit = 0;

	_transShiftX = 0;
	_transShiftY = 0;

	// transform the angle to radians
  _beta = _beta * (3.14159 / 180.);

}

void EUTelAPIXHistograms::processRunHeader( LCRunHeader* runHeader) {

  auto_ptr<EUTelRunHeaderImpl> eutelHeader( new EUTelRunHeaderImpl ( runHeader ) );
  eutelHeader->addProcessor( type() );

  _nRun++ ;

  // Decode and print out Run Header information - just a check

  int runNr = runHeader->getRunNumber();
  // convert to string
  char buf[256];
  sprintf(buf, "%i", runNr);
  std::string runNr_str(buf);

  message<MESSAGE5> ( log() << "Processing run header " << _nRun
                     << ", run nr " << runNr );

  const std::string detectorName = runHeader->getDetectorName();
  const std::string detectorDescription = runHeader->getDescription();

  message<MESSAGE5> ( log() << detectorName << " : " << detectorDescription ) ;


	// pick up correct alignment collection
	_alignmentCollectionNames.clear();
	for (unsigned i = 0; i < _alignmentCollectionSuffixes.size(); i++) {
		std::string	temp = "run"+runNr_str +_alignmentCollectionSuffixes[i];
		_alignmentCollectionNames.push_back(temp);
		cout << _alignmentCollectionNames[i] << endl;
	}

}

void EUTelAPIXHistograms::processEvent( LCEvent * event ) {


  EUTelEventImpl * euEvent = static_cast<EUTelEventImpl*> ( event );
  if ( euEvent->getEventType() == kEORE ) {
    message<DEBUG5> ( "EORE found: nothing else to do." );
    return;
  }

  bool debug = ( _debugCount>0 && _nEvt%_debugCount == 0);

  if ( _nEvt % 1000 == 0 ) {
    streamlog_out( MESSAGE2 ) << "Processing event "
                              << setw(6) << setiosflags(ios::right) << event->getEventNumber() << " in run "
                              << setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber() << setfill(' ')
                              << " (Total = " << setw(10) << _nEvt << ")" << resetiosflags(ios::left) << endl;
  }

  _nEvt ++ ;
  int evtNr = event->getEventNumber();


  if(debug)message<MESSAGE5> ( log() << "Processing record " << _nEvt << " == event " << evtNr );

  //
  // Get input collections
  //

  LCCollection* trackcol;
  try {
    trackcol = event->getCollection( _inputTrackColName ) ;
  } catch (lcio::DataNotAvailableException& e) {
    throw SkipEventException(this);
  }


  LCCollection* hitcol = NULL;
  bool _DUTok=true;
  try {
    hitcol = event->getCollection( _inputHitColName ) ;
  } catch (lcio::DataNotAvailableException& e) {
    message<ERROR5> ( log() << "Not able to get collection "
                     << _inputHitColName
                     << "\nfrom event " << event->getEventNumber()
                     << " in run " << event->getRunNumber()  );
    _DUTok=false;
    //
    // Do not skip event if DUT hits missing - efficiency and
    //   background calculations still have to be done!
    //
    //   throw SkipEventException(this);
  }

	LCCollection* original_zsdata = NULL;
	try {
		original_zsdata = event->getCollection( "cluster" ) ;
	} catch (lcio::DataNotAvailableException& e) {
		cout << "original_zsdata collection not available" << endl;
		cout << "Terminating... " << endl;
		abort();
    }



  // Clear local fit storage tables
  // 'fitted' table used for comparison with measured hits
  // 'bgfitted' used for comparison with hits from previous event

  _fittedX.clear();
  _fittedY.clear();

  _fittedXcorr.clear();
  _fittedYcorr.clear();
  _fittedZcorr.clear();



  _bgfittedX.clear();
  _bgfittedY.clear();


  // Loop over tracks in input track collection
  // Read fitted positions at DUT

  int nTrack = trackcol->getNumberOfElements()  ;


  if(debug)message<MESSAGE5> ( log() << "Total of " << nTrack << " tracks in input collection " );

  if (_manualDUTid >=10 ) 
  {	
      fillAPIXhits(hitcol, original_zsdata);
  }


  for(int itrack=0; itrack< nTrack ; itrack++)
    {
  
      Track * fittrack = dynamic_cast<Track*>( trackcol->getElementAt(itrack) ) ;
      std::vector<EVENT::TrackerHit*>  trackhits = fittrack->getTrackerHits();
      int nHit = 0;
      //----------------------------------------------------------------------------------------------
      // 17 August 2010 libov@mail.desy.de (1st implementation 14 July 2010)
      // The idea is to select only tracks that are matched to one of the APIX DUTs, so that its "intime"
      // check whether there is a hit in another plane near this track
      if (_onlyIntimeTracks)
      {
  		if ( ! ( hasMatchedHit ( fittrack ) ) ) continue;
      }
      //----------------------------------------------------------------------------------------------

      
      const float * por = fittrack->getReferencePoint();

      // Does track PoR match DUT position?

      double dist = por[2] - _zDUT ;

		double	fitX_corr = 0;
		double	fitY_corr = 0;
		double	fitZ_corr = 0;
		getTrackImpactPoint(fitX_corr, fitY_corr, fitZ_corr, fittrack, event);


      if(dist*dist < 1)
        {
          _fittedX.push_back(por[0]);
          _fittedY.push_back(por[1]);
          _fittedXcorr.push_back(fitX_corr);
          _fittedYcorr.push_back(fitY_corr);
          _fittedZcorr.push_back(fitZ_corr);
          _bgfittedX.push_back(por[0]);
          _bgfittedY.push_back(por[1]);
        }
      else
        {

          // Look at hits asigned to track
          // put this definition in the beginning of the loop - need to use there already
          //std::vector<EVENT::TrackerHit*>  trackhits = fittrack->getTrackerHits();
          // put nHit declaration above
          nHit =   trackhits.size();

          for(int ihit=0; ihit< nHit ; ihit++)
            {
              TrackerHit * meshit = trackhits.at(ihit);

              // Hit position

              const double * pos = meshit->getPosition();

              dist =  pos[2] - _zDUT ;

              // Look at fitted hits only!

              if( meshit->getType() >= 32 &&  dist*dist < 1 )
                {
                  _fittedX.push_back(pos[0]);
                  _fittedY.push_back(pos[1]);
                  _fittedXcorr.push_back(fitX_corr);
                  _fittedYcorr.push_back(fitY_corr);
                  _fittedZcorr.push_back(fitZ_corr);
                  _bgfittedX.push_back(pos[0]);
                  _bgfittedY.push_back(pos[1]);
                  break;
                }
            }
        }

      // End of loop over fitted tracks
    }

    
  if(debug)
  {
      message<MESSAGE5> ( log() << _fittedX.size()  << " fitted positions at DUT " );
  }

  if( _manualDUTid >= 10 )
  {
      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_NumberOfFittedTracksHistoName]))->fill(_fittedX.size());
  }
  
  if ( _fittedX.size() != 1 ) return;

  // Histograms of fitted positions

  for(int ifit=0; ifit < static_cast< int >(_fittedX.size()); ifit++)
    {
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA) 
      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_FittedXHistoName]))->fill(_fittedX[ifit]);        
        
      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_FittedYHistoName]))->fill(_fittedY[ifit]);

      (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_FittedXYHistoName]))->fill(_fittedX[ifit],_fittedY[ifit]);

      if(debug)message<MESSAGE5> ( log() << "Fit " << ifit
                                << "   X = " << _fittedX[ifit]
                                << "   Y = " << _fittedY[ifit]) ;
#endif
    }

  // Match fitted positions with hits from previous event (background estimate)
  // Skip first event

  if(_nEvt > 1)
    {
      int nMatch=0;
      double distmin;


      do{
        int bestfit=-1;
        int besthit=-1;

        distmin=_distMax*_distMax+1.;

        for(int ifit=0; ifit < static_cast< int >(_bgfittedX.size()); ifit++)
          for(int ihit=0; ihit< static_cast< int >(_bgmeasuredX.size()) ; ihit++)
            {
              double dist=
                (_bgmeasuredX[ihit]-_bgfittedX[ifit])*(_bgmeasuredX[ihit]-_bgfittedX[ifit])
                + (_bgmeasuredY[ihit]-_bgfittedY[ifit])*(_bgmeasuredY[ihit]-_bgfittedY[ifit]);

              if(dist<distmin)
                {
                  distmin=dist;
                  besthit=ihit;
                  bestfit=ifit;
                }
            }

        // Match found:

        if(distmin < _distMax*_distMax)
          {

            nMatch++;

            // Histograms of measured-fitted shifts
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

            (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_BgShiftXHistoName]))->fill(_bgmeasuredX[besthit]-_bgfittedX[bestfit]);

            (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_BgShiftYHistoName]))->fill(_bgmeasuredY[besthit]-_bgfittedY[bestfit]);

            (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_BgShiftXYHistoName]))->fill(_bgmeasuredX[besthit]-_bgfittedX[bestfit],_bgmeasuredY[besthit]-_bgfittedY[bestfit]);

            // Efficiency plots

            (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_BgEfficiencyXHistoName]))->fill(_bgfittedX[bestfit],1.);

            (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_BgEfficiencyYHistoName]))->fill(_bgfittedY[bestfit],1.);

            (dynamic_cast<AIDA::IProfile2D*> ( _aidaHistoMap[_BgEfficiencyXYHistoName]))->fill(_bgfittedX[bestfit],_bgfittedY[bestfit],1.);

#endif

            // Remove matched entries from the list (so the next matching pair
            // can be looked for)

            _bgfittedX.erase(_bgfittedX.begin()+bestfit);
            _bgfittedY.erase(_bgfittedY.begin()+bestfit);

            _bgmeasuredX.erase(_bgmeasuredX.begin()+besthit);
            _bgmeasuredY.erase(_bgmeasuredY.begin()+besthit);

          }

        // End of loop of matching background hits to fitted positions
      }
      while(distmin < _distMax*_distMax);


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

      // Background efficiency plots - tracks unmached to bg hits

      for(int ifit=0;ifit<static_cast< int >(_bgfittedX.size()); ifit++)
        {
          (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_BgEfficiencyXHistoName]))->fill(_bgfittedX[ifit],0.);

          (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_BgEfficiencyYHistoName]))->fill(_bgfittedY[ifit],0.);

          (dynamic_cast<AIDA::IProfile2D*> ( _aidaHistoMap[_BgEfficiencyXYHistoName]))->fill(_bgfittedX[ifit],_bgfittedY[ifit],0.);
        }

#endif

    }  // End of if(_nEvt > 1)  - end of background calculations



  // Clear local tables with measured position

  _measuredX.clear();
  _measuredY.clear();
  _measuredZ.clear();


  _localX.clear();
  _localY.clear();

  _bgmeasuredX.clear();
  _bgmeasuredY.clear();

  _lv1APIX.clear();
  _totAPIX.clear();
  _maxDifflv1.clear();
  _hitsInCluster.clear();

  _measuredXLOCAL.clear();
  _measuredYLOCAL.clear();

	// 24 august
	_clusterSizeX.clear();
	_clusterSizeY.clear();


  // Loop over hits in input collection
  // Read measured positions at DUT

  int nHit = 0;

  if(_DUTok) nHit = hitcol->getNumberOfElements();

  if(debug)message<MESSAGE5> ( log() << "Total of " << nHit << " tracker hits in input collection " );


  for(int ihit=0; ihit< nHit ; ihit++)
    {
      TrackerHit * meshit = dynamic_cast<TrackerHit*>( hitcol->getElementAt(ihit) ) ;

      // Hit position

      const double * pos = meshit->getPosition();

      double dist = pos[2] - _zDUT;

      if(dist*dist < 1)
        {
          // Apply alignment corrections

          double corrX  =   pos[0]*cos(_DUTalign.at(2))
            +pos[1]*sin(_DUTalign.at(2))
            +_DUTalign.at(0);

          double corrY  =  pos[1]*cos(_DUTalign.at(2))
            -pos[0]*sin(_DUTalign.at(2))
            +_DUTalign.at(1);

          _measuredX.push_back(corrX);
          _measuredY.push_back(corrY);
          _measuredZ.push_back(pos[2]);
          _bgmeasuredX.push_back(corrX);
          _bgmeasuredY.push_back(corrY);

			// ----  05 August 2010 libov@mail.desy.de --------------
			// Only for APIX at the moment
			if (_manualDUTid >= 10) {
				LCObjectVec clusterVec = (meshit->getRawHits());
				TrackerDataImpl * trackerData = dynamic_cast < TrackerDataImpl * > ( clusterVec[0] );
				EUTelSparseDataImpl <EUTelAPIXSparsePixel>  * sparseData = new EUTelSparseDataImpl <EUTelAPIXSparsePixel> (trackerData);
				EUTelAPIXSparsePixel * sparsePixel= new EUTelAPIXSparsePixel;
				float	tot=0;
				std::vector	<float> lv1Cluster;
				for (unsigned int i = 0; i < sparseData->size(); i++ ) {
					sparseData -> getSparsePixelAt (i, sparsePixel );
					tot += sparsePixel->getSignal();
					lv1Cluster.push_back ( sparsePixel->getTime() );
				}
				_totAPIX.push_back(tot);

				float	maxDifflv1 = 0;
				for (unsigned int i =1; i<lv1Cluster.size(); i++ ) {
					if (abs	(lv1Cluster[0] - lv1Cluster[i]) > maxDifflv1) {
						maxDifflv1 = abs(lv1Cluster[0] - lv1Cluster[i]);
					}
				}
				_lv1APIX.push_back(lv1Cluster[0]);
				_maxDifflv1.push_back(maxDifflv1);

				_hitsInCluster.push_back(sparseData->size());
				// this assumes SparseClustering algorithm!
				//EUTelSparseClusterImpl<EUTelAPIXSparsePixel> * sparseCluster = new EUTelSparseClusterImpl<EUTelAPIXSparsePixel> (trackerData);
				EUTelAPIXSparseClusterImpl * sparseCluster = new EUTelAPIXSparseClusterImpl (trackerData);
				if (sparseCluster -> getTotalCharge() != tot) {
					cout << "total charge got from cluster is not equal to calculated from its pixels. Check what you're doing. Terminating."<<endl;
					abort();
				}
				// get seed coordinates of a cluster
				float xCluSeed, yCluSeed;
				sparseCluster->getSeedCoord2(xCluSeed, yCluSeed);
				// get Center Of Gravity Shift - the procedure here should be the same as in the hitmaker!
				float	xShift, yShift;
				sparseCluster->getCenterOfGravityShift2( xShift, yShift);
				// Now transform cluster indices to coordinates in the local FoR
				//double measuredXLOCAL = ( static_cast<double> (xCluSeed) + static_cast<double>(xShift) + 0.5 ) * _xPitch ;
				//double measuredYLOCAL = ( static_cast<double> (yCluSeed) + static_cast<double>(yShift) + 0.5 ) * _yPitch ;
				double measuredXLOCAL = xCluSeed + xShift;
				double measuredYLOCAL = yCluSeed + yShift;

				_measuredXLOCAL.push_back(measuredXLOCAL);
				_measuredYLOCAL.push_back(measuredYLOCAL);

				// 07 September
				//float		xCluCenter, yCluCenter;
				//sparseCluster->getCenterCoord( xShift, yShift);		// no, because in hitmaker 3x3 center of gravity algorithm is used


				// 24 August
				int sizeX, sizeY;
				sparseCluster -> getClusterSize(sizeX, sizeY);
				_clusterSizeX.push_back(sizeX);
				_clusterSizeY.push_back(sizeY);

				delete sparseData;
				delete sparsePixel;
				delete sparseCluster;
			}
			// -----------------------------------------------------------
          // Local position should be taken from the cluster.
          // This is a temporary solution:

          double locX = pos[0];
          double locY = pos[1];

          // Subtract position of the central pixel

          int picX = static_cast< int >(locX/_pitchX);

          if(locX<0)picX--;

          locX-=(picX+0.5)*_pitchX;

          int picY = static_cast< int >(locY/_pitchY);

          if(locY<0)picY--;

          locY-=(picY+0.5)*_pitchY;

          _localX.push_back(locX);
          _localY.push_back(locY);
        }
    }

    if(debug)message<MESSAGE5> ( log() << _measuredX.size() << " hits at DUT " );

	if (_noHitYet && ( _measuredX.empty() ) ) {
		_eventsWithNoHit++;
		// anything to be done before returning?
		return;
	}
	if ((_noHitYet && (!_measuredX.empty()) ) || (_measuredX.size()>=1)){
		if (_noHitYet) cout << _eventsWithNoHit << "events before first hit" << endl;
		_noHitYet = false;
		// determine transformation ...
		getTransformationShifts();
	}

    // Histograms of measured positions

	if ( _manualDUTid >= 10 && _measuredX.size() != _measuredXLOCAL.size()) {
		cout<<"_measuredX.size() != _measuredXLOCAL.size()"<<endl;
		cout<< "Sizes of measured hits arrays in Telescope and Local FoR are not the same! Terminating! " << endl;
		abort();
	}


    for(int ihit=0; ihit < static_cast< int >(_measuredX.size()); ihit++)
    {
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_MeasuredXHistoName]))->fill(_measuredX[ihit]);


      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_MeasuredYHistoName]))->fill(_measuredY[ihit]);

      (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_MeasuredXYHistoName]))->fill(_measuredX[ihit],_measuredY[ihit]);

      if (_manualDUTid >= 10) {
          (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_MeasuredXLOCALHistoName]))->fill(_measuredXLOCAL[ihit]);
          (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_MeasuredYLOCALHistoName]))->fill(_measuredYLOCAL[ihit]);
          (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_MeasuredXYLOCALHistoName]))->fill(_measuredXLOCAL[ihit],_measuredYLOCAL[ihit]);
 
				double	z_sensor = _siPlanesLayerLayout->getSensitivePositionZ(_indexDUT)+ 0.5 * _siPlanesLayerLayout->getSensitiveThickness( _indexDUT );

				double	x = _measuredX[ihit];
				double	y = _measuredY[ihit];
				double	z = _measuredZ[ihit] - z_sensor;
				TransformToLocalFrame(x, y, z, event);
				(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap["diffX_correct"]))->fill(x - _measuredXLOCAL[ihit]);
				(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap["diffY_correct"]))->fill(y - _measuredYLOCAL[ihit]);
      }

#endif
      if(debug)message<MESSAGE5> ( log() << "Hit " << ihit
                                << "   X = " << _measuredX[ihit]
                                << "   Y = " << _measuredY[ihit]) ;
    }

	



  // Match measured and fitted positions

  int nMatch=0;
  double distmin;

  do{
    int bestfit=-1;
    int besthit=-1;

    distmin=_distMax*_distMax+1.;


    for(int ifit=0; ifit < static_cast< int >(_fittedX.size()); ifit++)
      for(int ihit=0; ihit< static_cast< int >(_measuredX.size()); ihit++)
        {
          double dist=
            (_measuredX[ihit]-_fittedX[ifit])*(_measuredX[ihit]-_fittedX[ifit])
            + (_measuredY[ihit]-_fittedY[ifit])*(_measuredY[ihit]-_fittedY[ifit]);


          if(dist<distmin)
            {
              distmin=dist;
              besthit=ihit;
              bestfit=ifit;
            }
        }

    // Match found:
    if(distmin < _distMax*_distMax)
      {

        nMatch++;

        // Matched hits positions

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_MatchedXHistoName]))->fill(_measuredX[besthit]);


        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_MatchedYHistoName]))->fill(_measuredY[besthit]);

        (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_MatchedXYHistoName]))->fill(_measuredX[besthit],_measuredY[besthit]);

        // Histograms of measured-fitted shifts

        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_ShiftXHistoName]))->fill(_measuredX[besthit]-_fittedX[bestfit]);

        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_ShiftYHistoName]))->fill(_measuredY[besthit]-_fittedY[bestfit]);

        (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_ShiftXYHistoName]))->fill(_measuredX[besthit]-_fittedX[bestfit],_measuredY[besthit]-_fittedY[bestfit]);

        (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_ShiftXvsYHistoName]))->fill(_fittedY[bestfit],_measuredX[besthit]-_fittedX[bestfit]);

        (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_ShiftYvsXHistoName]))->fill(_fittedX[bestfit],_measuredY[besthit]-_fittedY[bestfit]);

        (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_ShiftXvsX2DHistoName]))->fill(_fittedX[bestfit], _measuredX[besthit]-_fittedX[bestfit]);

        (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_ShiftXvsXHistoName]))->fill(_fittedX[bestfit], _measuredX[besthit]-_fittedX[bestfit]);

        (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_ShiftYvsY2DHistoName]))->fill(_fittedY[bestfit], _measuredY[besthit]-_fittedY[bestfit]);

        (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_ShiftYvsYHistoName]))->fill(_fittedY[bestfit], _measuredY[besthit]-_fittedY[bestfit]);

        (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_ShiftXvsY2DHistoName]))->fill(_fittedY[bestfit],_measuredX[besthit]-_fittedX[bestfit]);

        (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_ShiftYvsX2DHistoName]))->fill(_fittedX[bestfit],_measuredY[besthit]-_fittedY[bestfit]);


        // Eta function check plots


        (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_EtaXHistoName]))->fill(_localX[besthit],_measuredX[besthit]-_fittedX[bestfit]);

        (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_EtaYHistoName]))->fill(_localY[besthit],_measuredY[besthit]-_fittedY[bestfit]);

        (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_EtaX2DHistoName]))->fill(_localX[besthit],_measuredX[besthit]-_fittedX[bestfit]);

        (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_EtaY2DHistoName]))->fill(_localY[besthit],_measuredY[besthit]-_fittedY[bestfit]);

        (dynamic_cast<AIDA::IProfile2D*> ( _aidaHistoMap[_EtaX3DHistoName]))->fill(_localX[besthit],_localY[besthit],_measuredX[besthit]-_fittedX[bestfit]);

        (dynamic_cast<AIDA::IProfile2D*> ( _aidaHistoMap[_EtaY3DHistoName]))->fill(_localX[besthit],_localY[besthit],_measuredY[besthit]-_fittedY[bestfit]);

        // extend Eta histograms to 2 pitch range

#endif

        if(_localX[besthit]<0)
          _localX[besthit]+=_pitchX;
        else
          _localX[besthit]-=_pitchX;

        if(_localY[besthit]<0)
          _localY[besthit]+=_pitchY;
        else
          _localY[besthit]-=_pitchY;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

        (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_EtaXHistoName]))->fill(_localX[besthit],_measuredX[besthit]-_fittedX[bestfit]);

        (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_EtaYHistoName]))->fill(_localY[besthit],_measuredY[besthit]-_fittedY[bestfit]);

        (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_EtaX2DHistoName]))->fill(_localX[besthit],_measuredX[besthit]-_fittedX[bestfit]);

        (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_EtaY2DHistoName]))->fill(_localY[besthit],_measuredY[besthit]-_fittedY[bestfit]);


        // Efficiency plots
        (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_EfficiencyXHistoName]))->fill(_fittedX[bestfit],1.);

        (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_EfficiencyYHistoName]))->fill(_fittedY[bestfit],1.);

        (dynamic_cast<AIDA::IProfile2D*> ( _aidaHistoMap[_EfficiencyXYHistoName]))->fill(_fittedX[bestfit],_fittedY[bestfit],1.);


        // Noise plots

        (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_NoiseXHistoName]))->fill(_measuredX[besthit],0.);

        (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_NoiseYHistoName]))->fill(_measuredY[besthit],0.);

        (dynamic_cast<AIDA::IProfile2D*> ( _aidaHistoMap[_NoiseXYHistoName]))->fill(_measuredX[besthit],_measuredY[besthit],0.);

#endif

			// 05 August 2010 libov@mail.desy.de
			if (_manualDUTid >=10 ) {
				(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_totAPIXmatchedHistoName]))->fill(_totAPIX[besthit]);
				(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_hitsInMatchedClusterHistoName]))->fill(_hitsInCluster[besthit]);
				(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_lv1APIXmatchedHistoName]))->fill(_lv1APIX[besthit]);
				(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_maxDifflv1MatchedHistoName]))->fill(_maxDifflv1[besthit]);

				// --- 24 august ---
				(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_MatchedClusterSizeXHistoName]))->fill(_clusterSizeX[besthit]);
				(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_MatchedClusterSizeYHistoName]))->fill(_clusterSizeY[besthit]);
				// -----------------



				// -- 18 August
				// -- telescope frame of reference
				double	residualXTEL = _measuredX[besthit]-_fittedX[bestfit];
				double	residualYTEL = _measuredY[besthit]-_fittedY[bestfit];

				double	residualXLOCAL = _rot00 * residualXTEL + _rot01 * residualYTEL;
				double	residualYLOCAL = _rot10 * residualXTEL + _rot11 * residualYTEL;
				// -- consistency check
				//cout << (residualXTEL*residualXTEL+residualYTEL*residualYTEL) <<" <--> "<< (residualXLOCAL*residualXLOCAL + residualYLOCAL*residualYLOCAL)<<endl;
				if (((residualXTEL*residualXTEL+residualYTEL*residualYTEL) - (residualXLOCAL*residualXLOCAL + residualYLOCAL*residualYLOCAL)) > 0.0000001) {
					cout << "Bad matrix. Terminating "<< endl;
					abort();
				}
				double	fittedXLOCAL = _measuredXLOCAL[besthit] - residualXLOCAL;
				double	fittedYLOCAL = _measuredYLOCAL[besthit] - residualYLOCAL;


				// consistency check
				if (!_noHitYet) {


					double		fitX = _fittedXcorr[bestfit];
					double		fitY = _fittedYcorr[bestfit];


					double	z_sensor = _siPlanesLayerLayout->getSensitivePositionZ(_indexDUT)+ 0.5 * _siPlanesLayerLayout->getSensitiveThickness( _indexDUT );
					double	fitZ = (_fittedZcorr[bestfit] - z_sensor);

					(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap["Zdiff"]))->fill(1000 * (_measuredZ[besthit] - _fittedZcorr[bestfit]));

					TransformToLocalFrame(fitX, fitY, fitZ, event);

					/*cout << "measuredZ: "<<_measuredZ[besthit] <<" , fittedZ: "<<_fittedZcorr[bestfit]<<endl;
					cout << fitX << " " << fittedXLOCALMethod2 << " " << fittedXLOCAL << endl;
					cout << fitY << " " << fittedYLOCALMethod2 << " " << fittedYLOCAL << endl;*/

				(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap["diffX_crap"]))->fill(1000*(fitX - fittedXLOCAL));
				(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap["diffY_crap"]))->fill(1000*(fitY - fittedYLOCAL));

				fittedXLOCAL = fitX;
				fittedYLOCAL = fitY;


					/*if ((diffX>0.001) || (diffY>0.001)) {
						cout << "C r a p ! ! !" << endl;
						cout<< "fittedXLOCALMethod2 = "<<fittedXLOCALMethod2<<" fittedXLOCAL= "<<fittedXLOCAL<<" diff= "<<fittedXLOCALMethod2-fittedXLOCAL<<endl;
						cout<< "fittedYLOCALMethod2 = "<<fittedYLOCALMethod2<<" fittedYLOCAL= "<<fittedYLOCAL<<" diff= "<<fittedYLOCALMethod2-fittedYLOCAL<<endl;
						//abort();
					}*/
				}

				(dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_EfficiencyXLOCALHistoName]))->fill(fittedXLOCAL, 1.);
				(dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_EfficiencyYLOCALHistoName]))->fill(fittedYLOCAL, 1.);
				(dynamic_cast<AIDA::IProfile2D*> ( _aidaHistoMap[_EfficiencyXYLOCALHistoName]))->fill(fittedXLOCAL, fittedYLOCAL, 1.);

				// now stupidly done - to do: fill array before.  24 august
      		if (_manualDUTid >= 10) {
      			(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_FittedXLOCALHistoName]))->fill(fittedXLOCAL);
      			(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_FittedYLOCALHistoName]))->fill(fittedYLOCAL);
      			(dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_FittedXYLOCALHistoName]))->fill(fittedXLOCAL, fittedYLOCAL);

					if (_hitsInCluster[besthit] == 1) {
						if ( (_clusterSizeX[besthit] > 1) || (_clusterSizeY[besthit] > 1) ) {
							cout << "Bug! # of pixels = 1 but size X or size Y not! Terminating!"<<endl;
							abort();
						}
						// position of a fitted track inside a hit pixel
						float		X = (400./2) - (residualXLOCAL * 1000.);
						float		Y = (50./2) - (residualYLOCAL * 1000.);
						(dynamic_cast<AIDA::IProfile2D*> ( _aidaHistoMap[_totSinglePixelClustersXYHistoName]))->fill(X, Y, _totAPIX[besthit]);
					}
						// fittedXLOCAL was defined on line 1039


						bool		isMirrored = false;
						bool		longPixel = false;
						if ( (fittedXLOCAL < 0.6 ) || (fittedXLOCAL > 7.) ) longPixel = true;	// also the case if fitted posiiton OUTSIDE of the sensor in X!
						if ( ! longPixel ) {

							float		fittedIndexX = 1 + (fittedXLOCAL - 0.6) / _pitchX;

							float		fittedIndexY = ( _measuredYLOCAL[besthit] - residualYLOCAL ) / _pitchY;

							if (((static_cast< int >(floor(fittedIndexX + 1))) % 2) == 0) isMirrored = true;



							float		X = 1000 * (fittedIndexX - floor(fittedIndexX)) *_pitchX ;
							float		Y = 1000 * (fittedIndexY - floor(fittedIndexY)) *_pitchY ;

							if (isMirrored) X = 400. - X;

							(dynamic_cast<AIDA::IProfile2D*> ( _aidaHistoMap[_totAllClustersXYHistoName]))->fill(X, Y, _totAPIX[besthit]);


	        				if ( _hitsInCluster[besthit] > 1) {
								(dynamic_cast<AIDA::IProfile2D*> ( _aidaHistoMap[_ChargeSharingProbXYHistoName]))->fill(X, Y, 1.);
							} else if (_hitsInCluster[besthit] == 1) {
								(dynamic_cast<AIDA::IProfile2D*> ( _aidaHistoMap[_ChargeSharingProbXYHistoName]))->fill(X, Y, 0);
							} else {
								cout << "Bug, # of pixels in clusters < 0. Terminating! "<<endl;
								abort();
							}

							}
      		}

			}

        // Remove matched entries from the list (so the next matching pair
        // can be looked for)

        _fittedX.erase(_fittedX.begin()+bestfit);
        _fittedY.erase(_fittedY.begin()+bestfit);

        _fittedXcorr.erase(_fittedXcorr.begin()+bestfit);
        _fittedYcorr.erase(_fittedYcorr.begin()+bestfit);
        _fittedZcorr.erase(_fittedZcorr.begin()+bestfit);


		_measuredX.erase(_measuredX.begin()+besthit);
        _measuredY.erase(_measuredY.begin()+besthit);    
		_measuredZ.erase(_measuredZ.begin()+besthit);


        _localX.erase(_localX.begin()+besthit);
        _localY.erase(_localY.begin()+besthit);

			if (_manualDUTid >=10 ) {
				_totAPIX.erase(_totAPIX.begin()+besthit);
				_lv1APIX.erase(_lv1APIX.begin()+besthit);
				_hitsInCluster.erase(_hitsInCluster.begin()+besthit);
				_maxDifflv1.erase(_maxDifflv1.begin()+besthit);
				_measuredXLOCAL.erase(_measuredXLOCAL.begin()+besthit);
				_measuredYLOCAL.erase(_measuredYLOCAL.begin()+besthit);
				_clusterSizeX.erase(_clusterSizeX.begin()+besthit);
				_clusterSizeY.erase(_clusterSizeY.begin()+besthit);
			}
      }
    // End of loop of matching DUT hits to fitted positions
  }
  while(distmin < _distMax*_distMax);


  if(debug)
    {
      message<MESSAGE5> ( log() << nMatch << " DUT hits matched to fitted tracks ");
      message<MESSAGE5> ( log() << _measuredX.size() << " DUT hits not matched to any track ");
      message<MESSAGE5> ( log() << _fittedX.size() << " fitted tracks not matched to any DUT hit ");
    }

  if( _manualDUTid >= 10 )
  {
	(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_MatchedHitsHistoName])) -> fill (nMatch);
	(dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_MatchedHitsVSEventHistoName])) -> fill (evtNr, nMatch);
  // Efficiency plots - unmatched tracks
  }

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  for(int ifit=0; ifit < static_cast< int >(_fittedX.size()); ifit++)
    {
      (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_EfficiencyXHistoName]))->fill(_fittedX[ifit],0.);

      (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_EfficiencyYHistoName]))->fill(_fittedY[ifit],0.);

      (dynamic_cast<AIDA::IProfile2D*> ( _aidaHistoMap[_EfficiencyXYHistoName]))->fill(_fittedX[ifit],_fittedY[ifit],0.);


		// 19 August 2010 - plots in the local frame
	  if (_manualDUTid >= 10) {
        double	fittedXLOCALMethod2 = _rot00*(_fittedX[ifit] + _transShiftX) + _rot01 *  (_fittedY[ifit] + _transShiftY);
		double	fittedYLOCALMethod2 = _rot10*(_fittedX[ifit] + _transShiftX) + _rot11 *  (_fittedY[ifit] + _transShiftY);

		(dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_EfficiencyXLOCALHistoName]))->fill(fittedXLOCALMethod2, 0.);
		(dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_EfficiencyYLOCALHistoName]))->fill(fittedYLOCALMethod2, 0.);
		(dynamic_cast<AIDA::IProfile2D*> ( _aidaHistoMap[_EfficiencyXYLOCALHistoName]))->fill(fittedXLOCALMethod2, fittedYLOCALMethod2, 0.);

     	(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_FittedXLOCALHistoName]))->fill(fittedXLOCALMethod2);
      	(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_FittedYLOCALHistoName]))->fill(fittedYLOCALMethod2);
      	(dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_FittedXYLOCALHistoName]))->fill(fittedXLOCALMethod2, fittedYLOCALMethod2);
      }


    }


  // Noise plots - unmatched hits

  for(int ihit=0; ihit < static_cast< int >(_measuredX.size()); ihit++)
    {
      (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_NoiseXHistoName]))->fill(_measuredX[ihit],1.);

      (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_NoiseYHistoName]))->fill(_measuredY[ihit],1.);

      (dynamic_cast<AIDA::IProfile2D*> ( _aidaHistoMap[_NoiseXYHistoName]))->fill(_measuredX[ihit],_measuredY[ihit],1.);


      // Unmatched hit positions

      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_UnMatchedXHistoName]))->fill(_measuredX[ihit]);


      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_UnMatchedYHistoName]))->fill(_measuredY[ihit]);

      (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_UnMatchedXYHistoName]))->fill(_measuredX[ihit],_measuredY[ihit]);

		// 06 August 2010 libov@mail.desy.de
		if (_manualDUTid >=10 ) {
			(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_totAPIXunmatchedHistoName]))->fill(_totAPIX[ihit]);
			(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_hitsInUnmatchedClusterHistoName]))->fill(_hitsInCluster[ihit]);
			(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_lv1APIXunmatchedHistoName]))->fill(_lv1APIX[ihit]);
			(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_maxDifflv1UnmatchedHistoName]))->fill(_maxDifflv1[ihit]);

			(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_UnmatchedClusterSizeXHistoName]))->fill(_clusterSizeX[ihit]);
			(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_UnmatchedClusterSizeYHistoName]))->fill(_clusterSizeY[ihit]);

		}

    }

#endif

  return;
}



void EUTelAPIXHistograms::check( LCEvent * /* evt */ ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void EUTelAPIXHistograms::end(){

  // Nothing to do here
}



void EUTelAPIXHistograms::bookHistos()
{

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  message<MESSAGE5> ( log() << "Booking histograms " );


  message<MESSAGE5> ( log() << "Histogram information searched in " << _histoInfoFileName);

  auto_ptr<EUTelHistogramManager> histoMgr( new EUTelHistogramManager( _histoInfoFileName ));
  EUTelHistogramInfo    * histoInfo;
  bool                    isHistoManagerAvailable;

  try {
    isHistoManagerAvailable = histoMgr->init();
  } catch ( ios::failure& e) {
    message<ERROR5> ( log() << "I/O problem with " << _histoInfoFileName << "\n"
                     << "Continuing without histogram manager"    );
    isHistoManagerAvailable = false;
  } catch ( ParseException& e ) {
    message<ERROR5> ( log() << e.what() << "\n"
                     << "Continuing without histogram manager" );
    isHistoManagerAvailable = false;
  }



  // Measured position in X

  int    measXNBin  = 400;
  double measXMin   = -10.;
  double measXMax   =  10.;
  string measXTitle = "Measured particle position in X";
  string matchXTitle = "Matched hit position in X";
  string unmatchXTitle = "Unmatched hit position in X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_MeasuredXHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          measXNBin = histoInfo->_xBin;
          measXMin  = histoInfo->_xMin;
          measXMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) measXTitle = histoInfo->_title;
        }
    }


  AIDA::IHistogram1D * measXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_MeasuredXHistoName.c_str(),measXNBin,measXMin,measXMax);

  measXHisto->setTitle(measXTitle.c_str());

  _aidaHistoMap.insert(make_pair(_MeasuredXHistoName, measXHisto));


  AIDA::IHistogram1D * matchXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_MatchedXHistoName.c_str(),measXNBin,measXMin,measXMax);

  matchXHisto->setTitle(matchXTitle.c_str());

  _aidaHistoMap.insert(make_pair(_MatchedXHistoName, matchXHisto));


  AIDA::IHistogram1D * unmatchXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_UnMatchedXHistoName.c_str(),measXNBin,measXMin,measXMax);

  unmatchXHisto->setTitle(unmatchXTitle.c_str());

  _aidaHistoMap.insert(make_pair(_UnMatchedXHistoName, unmatchXHisto));


  // Measured position in Y

  int    measYNBin  = 400;
  double measYMin   = -10.;
  double measYMax   =  10.;
  string measYTitle = "Measured particle position in Y";
  string matchYTitle = "Matched hit position in Y";
  string unmatchYTitle = "Unmatched hit position in Y";



  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_MeasuredYHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          measYNBin = histoInfo->_xBin;
          measYMin  = histoInfo->_xMin;
          measYMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) measYTitle = histoInfo->_title;
        }
    }


  AIDA::IHistogram1D * measYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_MeasuredYHistoName.c_str(),measYNBin,measYMin,measYMax);

  measYHisto->setTitle(measYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_MeasuredYHistoName, measYHisto));



  AIDA::IHistogram1D * matchYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_MatchedYHistoName.c_str(),measYNBin,measYMin,measYMax);

  matchYHisto->setTitle(matchYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_MatchedYHistoName, matchYHisto));


  AIDA::IHistogram1D * unmatchYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_UnMatchedYHistoName.c_str(),measYNBin,measYMin,measYMax);

  unmatchYHisto->setTitle(unmatchYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_UnMatchedYHistoName, unmatchYHisto));


  // Measured position in X-Y

  measXNBin  = 100;
  measXMin   = -10.;
  measXMax   =  10.;
  measYNBin  = 100;
  measYMin   = -10.;
  measYMax   =  10.;
  string measXYTitle = "Measured particle position in XY";
  string matchXYTitle = "Matched hit position in XY";
  string unmatchXYTitle = "Unmatched hit position in XY";



  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_MeasuredXYHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          measXNBin = histoInfo->_xBin;
          measXMin  = histoInfo->_xMin;
          measXMax  = histoInfo->_xMax;
          measYNBin = histoInfo->_yBin;
          measYMin  = histoInfo->_yMin;
          measYMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) measXYTitle = histoInfo->_title;
        }
    }

  AIDA::IHistogram2D * measXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _MeasuredXYHistoName.c_str(),measXNBin,measXMin,measXMax,measYNBin,measYMin,measYMax);

  measXYHisto->setTitle(measXYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_MeasuredXYHistoName, measXYHisto));


  AIDA::IHistogram2D * matchXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _MatchedXYHistoName.c_str(),measXNBin,measXMin,measXMax,measYNBin,measYMin,measYMax);

  matchXYHisto->setTitle(matchXYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_MatchedXYHistoName, matchXYHisto));



  AIDA::IHistogram2D * unmatchXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _UnMatchedXYHistoName.c_str(),measXNBin,measXMin,measXMax,measYNBin,measYMin,measYMax);

  unmatchXYHisto->setTitle(unmatchXYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_UnMatchedXYHistoName, unmatchXYHisto));


  // Fitted position in X

  int    fitXNBin  = 400;
  double fitXMin   = -10.;
  double fitXMax   =  10.;
  string fitXTitle = "Fitted particle position in X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_FittedXHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          fitXNBin = histoInfo->_xBin;
          fitXMin  = histoInfo->_xMin;
          fitXMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) fitXTitle = histoInfo->_title;
        }
    }



  AIDA::IHistogram1D * fitXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_FittedXHistoName.c_str(),fitXNBin,fitXMin,fitXMax);

  fitXHisto->setTitle(fitXTitle.c_str());

  _aidaHistoMap.insert(make_pair(_FittedXHistoName, fitXHisto));



  // Fitted position in Y

  int    fitYNBin  = 400;
  double fitYMin   = -10.;
  double fitYMax   =  10.;
  string fitYTitle = "Fitted particle position in Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_FittedYHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          fitYNBin = histoInfo->_xBin;
          fitYMin  = histoInfo->_xMin;
          fitYMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) fitYTitle = histoInfo->_title;
        }
    }




  AIDA::IHistogram1D * fitYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_FittedYHistoName.c_str(),fitYNBin,fitYMin,fitYMax);

  fitYHisto->setTitle(fitYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_FittedYHistoName, fitYHisto));




  // Fitted position in X-Y

  fitXNBin  = 100;
  fitXMin   = -10.;
  fitXMax   =  10.;
  fitYNBin  = 100;
  fitYMin   = -10.;
  fitYMax   =  10.;
  string fitXYTitle = "Fitted particle position in XY";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_FittedXYHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          fitXNBin = histoInfo->_xBin;
          fitXMin  = histoInfo->_xMin;
          fitXMax  = histoInfo->_xMax;
          fitYNBin = histoInfo->_yBin;
          fitYMin  = histoInfo->_yMin;
          fitYMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) fitXYTitle = histoInfo->_title;
        }
    }



  AIDA::IHistogram2D * fitXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _FittedXYHistoName.c_str(),fitXNBin,fitXMin,fitXMax,fitYNBin,fitYMin,fitYMax);

  fitXYHisto->setTitle(fitXYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_FittedXYHistoName, fitXYHisto));




  // Measured - fitted position in X

  int    shiftXNBin  = 400;
  double shiftXMin   = -0.5;
  double shiftXMax   = 0.5;
  string shiftXTitle = "Measured - fitted X position";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ShiftXHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          shiftXNBin = histoInfo->_xBin;
          shiftXMin  = histoInfo->_xMin;
          shiftXMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) shiftXTitle = histoInfo->_title;
        }
    }


  AIDA::IHistogram1D * shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);

  shiftXHisto->setTitle(shiftXTitle.c_str());

  _aidaHistoMap.insert(make_pair(_ShiftXHistoName, shiftXHisto));


  // Corresponding background histogram

  shiftXTitle = "Measured - fitted X position (background)";

  AIDA::IHistogram1D * bgshiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_BgShiftXHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);

  bgshiftXHisto->setTitle(shiftXTitle.c_str());

  _aidaHistoMap.insert(make_pair(_BgShiftXHistoName, bgshiftXHisto));



  // Measured - fitted position in Y

  int    shiftYNBin  = 400;
  double shiftYMin   = -0.5;
  double shiftYMax   = 0.5;
  string shiftYTitle = "Measured - fitted Y position";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ShiftYHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          shiftYNBin = histoInfo->_xBin;
          shiftYMin  = histoInfo->_xMin;
          shiftYMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) shiftYTitle = histoInfo->_title;
        }
    }

  AIDA::IHistogram1D * shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYHistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);

  shiftYHisto->setTitle(shiftYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_ShiftYHistoName, shiftYHisto));


  // Corresponding background histogram

  shiftYTitle = "Measured - fitted Y position (background)";

  AIDA::IHistogram1D * bgshiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_BgShiftYHistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);

  bgshiftYHisto->setTitle(shiftYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_BgShiftYHistoName, bgshiftYHisto));




  // Measured - fitted position in X  vs Y

  shiftXNBin  = 200;
  shiftXMin   = -10.;
  shiftXMax   =  10.;
  double shiftVMin   = -0.5;
  double shiftVMax   = 0.5;
  shiftXTitle = "Measured - fitted X position vs Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ShiftXvsYHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          shiftXNBin = histoInfo->_xBin;
          shiftXMin  = histoInfo->_xMin;
          shiftXMax  = histoInfo->_xMax;
          shiftVMin  = histoInfo->_yMin;
          shiftVMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) shiftXTitle = histoInfo->_title;
        }
    }


  AIDA::IProfile1D * shiftXvsYHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_ShiftXvsYHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftVMin,shiftVMax);

  shiftXvsYHisto->setTitle(shiftXTitle.c_str());

  _aidaHistoMap.insert(make_pair(_ShiftXvsYHistoName, shiftXvsYHisto));


  // Measured - fitted position in X vs X
  shiftXTitle = "Measured - fitted X position vs X"; 
  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ShiftXvsXHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          shiftXNBin = histoInfo->_xBin;
          shiftXMin  = histoInfo->_xMin;
          shiftXMax  = histoInfo->_xMax;
          shiftVMin  = histoInfo->_yMin;
          shiftVMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) shiftXTitle = histoInfo->_title;
        }
    }
  
  AIDA::IProfile1D * shiftXvsXHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_ShiftXvsXHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftVMin,shiftVMax);

  shiftXvsXHisto->setTitle(shiftXTitle.c_str());

  _aidaHistoMap.insert(make_pair(_ShiftXvsXHistoName, shiftXvsXHisto));



  // Measured - fitted position in Y vs X

  shiftYNBin  = 200;
  shiftYMin   = -10.;
  shiftYMax   =  10.;
  shiftVMin   = -0.5;
  shiftVMax   = 0.5;
  shiftYTitle = "Measured - fitted Y position vs X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ShiftYvsXHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          shiftYNBin = histoInfo->_xBin;
          shiftYMin  = histoInfo->_xMin;
          shiftYMax  = histoInfo->_xMax;
          shiftVMin  = histoInfo->_yMin;
          shiftVMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) shiftYTitle = histoInfo->_title;
        }
    }

  AIDA::IProfile1D * shiftYvsXHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_ShiftYvsXHistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax,shiftVMin,shiftVMax);

  shiftYvsXHisto->setTitle(shiftYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_ShiftYvsXHistoName, shiftYvsXHisto));



  // Measured - fitted position in Y vs Y

  shiftYTitle = "Measured - fitted Y position vs Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ShiftYvsYHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          shiftYNBin = histoInfo->_xBin;
          shiftYMin  = histoInfo->_xMin;
          shiftYMax  = histoInfo->_xMax;
          shiftVMin  = histoInfo->_yMin;
          shiftVMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) shiftYTitle = histoInfo->_title;
        }
    }

  AIDA::IProfile1D * shiftYvsYHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_ShiftYvsYHistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax,shiftVMin,shiftVMax);

  shiftYvsYHisto->setTitle(shiftYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_ShiftYvsYHistoName, shiftYvsYHisto));


  // Measured - fitted position in X  vs Y (2D plot)

  shiftXNBin  = 200;
  shiftXMin   = -10.;
  shiftXMax   =  10.;
  int shiftVNBin  = 200;
  shiftVMin   = -0.5;
  shiftVMax   = 0.5;
  shiftXTitle = "Measured - fitted X position vs Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ShiftXvsY2DHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          shiftXNBin = histoInfo->_xBin;
          shiftXMin  = histoInfo->_xMin;
          shiftXMax  = histoInfo->_xMax;
          shiftVNBin = histoInfo->_yBin;
          shiftVMin  = histoInfo->_yMin;
          shiftVMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) shiftXTitle = histoInfo->_title;
        }
    }


  AIDA::IHistogram2D * shiftXvsY2DHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D(_ShiftXvsY2DHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftVNBin,shiftVMin,shiftVMax);

  shiftXvsY2DHisto->setTitle(shiftXTitle.c_str());

  _aidaHistoMap.insert(make_pair(_ShiftXvsY2DHistoName, shiftXvsY2DHisto));



  // Measured - fitted position in X vs X (2D plot) 

  shiftXTitle = "Measured - fitted X position vs X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ShiftXvsX2DHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          shiftXNBin = histoInfo->_xBin;
          shiftXMin  = histoInfo->_xMin;
          shiftXMax  = histoInfo->_xMax;
          shiftVNBin = histoInfo->_yBin;
          shiftVMin  = histoInfo->_yMin;
          shiftVMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) shiftXTitle = histoInfo->_title;
        }
    }


  AIDA::IHistogram2D * shiftXvsX2DHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D(_ShiftXvsX2DHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftVNBin,shiftVMin,shiftVMax);

  shiftXvsX2DHisto->setTitle(shiftXTitle.c_str());

  _aidaHistoMap.insert(make_pair(_ShiftXvsX2DHistoName, shiftXvsX2DHisto));



  // Measured - fitted position in Y vs X  (2D plot)

  shiftYNBin  = 200;
  shiftYMin   = -10.;
  shiftYMax   =  10.;
  shiftVNBin  = 500;
  shiftVMin   = -0.5;
  shiftVMax   = 0.5;
  shiftYTitle = "Measured - fitted Y position vs X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ShiftYvsX2DHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          shiftYNBin = histoInfo->_xBin;
          shiftYMin  = histoInfo->_xMin;
          shiftYMax  = histoInfo->_xMax;
          shiftVNBin = histoInfo->_yBin;
          shiftVMin  = histoInfo->_yMin;
          shiftVMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) shiftYTitle = histoInfo->_title;
        }
    }

  AIDA::IHistogram2D * shiftYvsX2DHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D(_ShiftYvsX2DHistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax,shiftVNBin,shiftVMin,shiftVMax);

  shiftYvsX2DHisto->setTitle(shiftYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_ShiftYvsX2DHistoName, shiftYvsX2DHisto));


  // Measured - fitted position in Y vs Y (2D plot)
  shiftYTitle = "Measured - fitted Y position vs Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ShiftYvsY2DHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          shiftYNBin = histoInfo->_xBin;
          shiftYMin  = histoInfo->_xMin;
          shiftYMax  = histoInfo->_xMax;
          shiftVNBin = histoInfo->_yBin;
          shiftVMin  = histoInfo->_yMin;
          shiftVMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) shiftYTitle = histoInfo->_title;
        }
    }

  AIDA::IHistogram2D * shiftYvsY2DHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D(_ShiftYvsY2DHistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax,shiftVNBin,shiftVMin,shiftVMax);

  shiftYvsY2DHisto->setTitle(shiftYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_ShiftYvsY2DHistoName, shiftYvsY2DHisto));


  // Measured - fitted position in X-Y

  shiftXNBin  = 500;
  shiftXMin   = -0.5;
  shiftXMax   = 0.5;
  shiftYNBin  = 500;
  shiftYMin   = -0.5;
  shiftYMax   = 0.5;
  string shiftXYTitle = "Measured - fitted position in XY";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ShiftXYHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          shiftXNBin = histoInfo->_xBin;
          shiftXMin  = histoInfo->_xMin;
          shiftXMax  = histoInfo->_xMax;
          shiftYNBin = histoInfo->_yBin;
          shiftYMin  = histoInfo->_yMin;
          shiftYMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) shiftXYTitle = histoInfo->_title;
        }
    }



  AIDA::IHistogram2D * shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);

  shiftXYHisto->setTitle(shiftXYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_ShiftXYHistoName, shiftXYHisto));



  // Corresponding background histogram

  shiftXYTitle = "Measured - fitted position in X-Y (background)";

  AIDA::IHistogram2D * bgshiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D(_BgShiftXYHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);

  bgshiftXYHisto->setTitle(shiftXYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_BgShiftXYHistoName, bgshiftXYHisto));



  // Efficiency as a function of the fitted position in X

  int    effiXNBin  = 200;
  double effiXMin   = -10.;
  double effiXMax   =  10.;
  string effiXTitle = "Efficiency vs particle position in X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_EfficiencyXHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          effiXNBin = histoInfo->_xBin;
          effiXMin  = histoInfo->_xMin;
          effiXMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) effiXTitle = histoInfo->_title;
        }
    }



  AIDA::IProfile1D * effiXHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_EfficiencyXHistoName.c_str(),effiXNBin,effiXMin,effiXMax);

  effiXHisto->setTitle(effiXTitle.c_str());

  _aidaHistoMap.insert(make_pair(_EfficiencyXHistoName, effiXHisto));


  // Corresponding background histogram


  effiXTitle = "Background match probability vs particle position in X";

  AIDA::IProfile1D * bgeffiXHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_BgEfficiencyXHistoName.c_str(),effiXNBin,effiXMin,effiXMax);

  bgeffiXHisto->setTitle(effiXTitle.c_str());

  _aidaHistoMap.insert(make_pair(_BgEfficiencyXHistoName, bgeffiXHisto));



  // Efficiency as a function of the fitted position in Y

  int    effiYNBin  = 200;
  double effiYMin   = -10.;
  double effiYMax   =  10.;
  string effiYTitle = "Efficiency vs particle position in Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_EfficiencyYHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          effiYNBin = histoInfo->_xBin;
          effiYMin  = histoInfo->_xMin;
          effiYMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) effiYTitle = histoInfo->_title;
        }
    }


  AIDA::IProfile1D * effiYHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_EfficiencyYHistoName.c_str(),effiYNBin,effiYMin,effiYMax);

  effiYHisto->setTitle(effiYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_EfficiencyYHistoName, effiYHisto));


  // Corresponding background histogram

  effiYTitle =  "Background match probability vs particle position in Y";

  AIDA::IProfile1D * bgeffiYHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_BgEfficiencyYHistoName.c_str(),effiYNBin,effiYMin,effiYMax);

  bgeffiYHisto->setTitle(effiYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_BgEfficiencyYHistoName, bgeffiYHisto));



  // Efficiency as a function of the fitted position in X-Y

  effiXNBin  = 50;
  effiXMin   = -10.;
  effiXMax   =  10.;
  effiYNBin  = 50;
  effiYMin   =  -5.;
  effiYMax   =   5.;
  string effiXYTitle = "Efficiency vs particle position in XY";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_EfficiencyXYHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          effiXNBin = histoInfo->_xBin;
          effiXMin  = histoInfo->_xMin;
          effiXMax  = histoInfo->_xMax;
          effiYNBin = histoInfo->_yBin;
          effiYMin  = histoInfo->_yMin;
          effiYMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) effiXYTitle = histoInfo->_title;
        }
    }



  AIDA::IProfile2D * effiXYHisto = AIDAProcessor::histogramFactory(this)->createProfile2D( _EfficiencyXYHistoName.c_str(),effiXNBin,effiXMin,effiXMax,effiYNBin,effiYMin,effiYMax);

  effiXYHisto->setTitle(effiXYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_EfficiencyXYHistoName, effiXYHisto));

  // Corresponding background histogram

  effiXYTitle =  "Background match probability vs particle position in XY";

  AIDA::IProfile2D * bgeffiXYHisto = AIDAProcessor::histogramFactory(this)->createProfile2D( _BgEfficiencyXYHistoName.c_str(),effiXNBin,effiXMin,effiXMax,effiYNBin,effiYMin,effiYMax);

  bgeffiXYHisto->setTitle(effiXYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_BgEfficiencyXYHistoName, bgeffiXYHisto));



  // Noise as a function of the measured position in X

  int    noiseXNBin  = 100;
  double noiseXMin   = -10.;
  double noiseXMax   =  10.;
  string noiseXTitle = "Noise fraction vs particle position in X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_NoiseXHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          noiseXNBin = histoInfo->_xBin;
          noiseXMin  = histoInfo->_xMin;
          noiseXMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) noiseXTitle = histoInfo->_title;
        }
    }



  AIDA::IProfile1D * noiseXHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_NoiseXHistoName.c_str(),noiseXNBin,noiseXMin,noiseXMax);

  noiseXHisto->setTitle(noiseXTitle.c_str());

  _aidaHistoMap.insert(make_pair(_NoiseXHistoName, noiseXHisto));



  // Noise as a function of the measured position in Y

  int    noiseYNBin  = 100;
  double noiseYMin   = -10.;
  double noiseYMax   =  10.;
  string noiseYTitle = "Noise fraction vs particle position in Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_NoiseYHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          noiseYNBin = histoInfo->_xBin;
          noiseYMin  = histoInfo->_xMin;
          noiseYMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) noiseYTitle = histoInfo->_title;
        }
    }




  AIDA::IProfile1D * noiseYHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_NoiseYHistoName.c_str(),noiseYNBin,noiseYMin,noiseYMax);

  noiseYHisto->setTitle(noiseYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_NoiseYHistoName, noiseYHisto));




  // Noise as a function of the measured position in X-Y

  noiseXNBin  = 50;
  noiseXMin   = -10.;
  noiseXMax   =  10.;
  noiseYNBin  = 50;
  noiseYMin   = -10.;
  noiseYMax   =  10.;
  string noiseXYTitle = "Noise fraction vs particle position in XY";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_NoiseXYHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          noiseXNBin = histoInfo->_xBin;
          noiseXMin  = histoInfo->_xMin;
          noiseXMax  = histoInfo->_xMax;
          noiseYNBin = histoInfo->_yBin;
          noiseYMin  = histoInfo->_yMin;
          noiseYMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) noiseXYTitle = histoInfo->_title;
        }
    }



  AIDA::IProfile2D * noiseXYHisto = AIDAProcessor::histogramFactory(this)->createProfile2D( _NoiseXYHistoName.c_str(),noiseXNBin,noiseXMin,noiseXMax,noiseYNBin,noiseYMin,noiseYMax);

  noiseXYHisto->setTitle(noiseXYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_NoiseXYHistoName, noiseXYHisto));




  // Eta function check: measured - fitted position in X  vs  local X

  int etaXNBin  = 200;
  double etaXMin   = -0.03;
  double etaXMax   = 0.03;
  double etaVMin   = -0.03;
  double etaVMax   = 0.03;
  string etaXTitle = "Measured - fitted X position vs local X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_EtaXHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          etaXNBin = histoInfo->_xBin;
          etaXMin  = histoInfo->_xMin;
          etaXMax  = histoInfo->_xMax;
          etaVMin  = histoInfo->_yMin;
          etaVMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) etaXTitle = histoInfo->_title;
        }
    }


  AIDA::IProfile1D * etaXHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_EtaXHistoName.c_str(),etaXNBin,etaXMin,etaXMax,etaVMin,etaVMax);

  etaXHisto->setTitle(etaXTitle.c_str());

  _aidaHistoMap.insert(make_pair(_EtaXHistoName, etaXHisto));



  // Eta function check: measured - fitted position in Y vs local Y

  int etaYNBin  = 200;
  double etaYMin   = -0.03;
  double etaYMax   = 0.03;
  etaVMin   = -0.03;
  etaVMax   = 0.03;
  string etaYTitle = "Measured - fitted Y position vs local Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_EtaYHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          etaYNBin = histoInfo->_xBin;
          etaYMin  = histoInfo->_xMin;
          etaYMax  = histoInfo->_xMax;
          etaVMin  = histoInfo->_yMin;
          etaVMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) etaYTitle = histoInfo->_title;
        }
    }

  AIDA::IProfile1D * etaYHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_EtaYHistoName.c_str(),etaYNBin,etaYMin,etaYMax,etaVMin,etaVMax);

  etaYHisto->setTitle(etaYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_EtaYHistoName, etaYHisto));




  // Eta function check: measured - fitted position in X  vs local X (2D plot)

  etaXNBin  = 200;
  etaXMin   = -0.03;
  etaXMax   = 0.03;
  int etaVNBin  = 200;
  etaVMin   = -0.03;
  etaVMax   = 0.03;
  etaXTitle = "Measured - fitted X position vs local X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_EtaX2DHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          etaXNBin = histoInfo->_xBin;
          etaXMin  = histoInfo->_xMin;
          etaXMax  = histoInfo->_xMax;
          etaVNBin = histoInfo->_yBin;
          etaVMin  = histoInfo->_yMin;
          etaVMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) etaXTitle = histoInfo->_title;
        }
    }


  AIDA::IHistogram2D * etaX2DHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D(_EtaX2DHistoName.c_str(),etaXNBin,etaXMin,etaXMax,etaVNBin,etaVMin,etaVMax);

  etaX2DHisto->setTitle(etaXTitle.c_str());

  _aidaHistoMap.insert(make_pair(_EtaX2DHistoName, etaX2DHisto));



  // Measured - fitted position in Y vs Y  (2D plot)

  etaYNBin  = 200;
  etaYMin   = -0.03;
  etaYMax   = 0.03;
  etaVNBin  = 200;
  etaVMin   = -0.03;
  etaVMax   = 0.03;
  etaYTitle = "Measured - fitted Y position vs local Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_EtaY2DHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          etaYNBin = histoInfo->_xBin;
          etaYMin  = histoInfo->_xMin;
          etaYMax  = histoInfo->_xMax;
          etaVNBin = histoInfo->_yBin;
          etaVMin  = histoInfo->_yMin;
          etaVMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) etaYTitle = histoInfo->_title;
        }
    }

  AIDA::IHistogram2D * etaY2DHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D(_EtaY2DHistoName.c_str(),etaYNBin,etaYMin,etaYMax,etaVNBin,etaVMin,etaVMax);

  etaY2DHisto->setTitle(etaYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_EtaY2DHistoName, etaY2DHisto));



  // Eta function check: measured - fitted position in X  vs  local X-Y ("3D")

  etaXNBin  = 200;
  etaXMin   = -0.03;
  etaXMax   = 0.03;
  etaYNBin  = 200;
  etaYMin   = -0.03;
  etaYMax   = 0.03;
  etaVMin   = -0.03;
  etaVMax   = 0.03;
  etaXTitle = "Measured - fitted X position vs local X-Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_EtaX3DHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          etaXNBin = histoInfo->_xBin;
          etaXMin  = histoInfo->_xMin;
          etaXMax  = histoInfo->_xMax;
          etaYNBin = histoInfo->_yBin;
          etaYMin  = histoInfo->_yMin;
          etaYMax  = histoInfo->_yMax;
          etaVMin  = histoInfo->_zMin;
          etaVMax  = histoInfo->_zMax;
          if ( histoInfo->_title != "" ) etaXTitle = histoInfo->_title;
        }
    }


  AIDA::IProfile2D * etaX3DHisto = AIDAProcessor::histogramFactory(this)->createProfile2D(_EtaX3DHistoName.c_str(),etaXNBin,etaXMin,etaXMax,etaYNBin,etaYMin,etaYMax,etaVMin,etaVMax);

  etaX3DHisto->setTitle(etaXTitle.c_str());

  _aidaHistoMap.insert(make_pair(_EtaX3DHistoName, etaX3DHisto));



  // Eta function check: measured - fitted position in Y vs local X-Y ("3D")

  etaXNBin  = 200;
  etaXMin   = -0.03;
  etaXMax   = 0.03;
  etaYNBin  = 200;
  etaYMin   = -0.03;
  etaYMax   = 0.03;
  etaVMin   = -0.03;
  etaVMax   = 0.03;
  etaYTitle = "Measured - fitted Y position vs local X-Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_EtaY3DHistoName);
      if ( histoInfo )
        {
          message<DEBUG5> ( log() << (* histoInfo ) );
          etaXNBin = histoInfo->_xBin;
          etaXMin  = histoInfo->_xMin;
          etaXMax  = histoInfo->_xMax;
          etaYNBin = histoInfo->_yBin;
          etaYMin  = histoInfo->_yMin;
          etaYMax  = histoInfo->_yMax;
          etaVMin  = histoInfo->_zMin;
          etaVMax  = histoInfo->_zMax;
          if ( histoInfo->_title != "" ) etaYTitle = histoInfo->_title;
        }
    }

  AIDA::IProfile2D * etaY3DHisto = AIDAProcessor::histogramFactory(this)->createProfile2D(_EtaY3DHistoName.c_str(),etaXNBin,etaXMin,etaXMax,etaYNBin,etaYMin,etaYMax,etaVMin,etaVMax);

  etaY3DHisto->setTitle(etaYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_EtaY3DHistoName, etaY3DHisto));

	// 05 August libov@mail.desy.de
	if (_manualDUTid >=10 ) {

		AIDA::IHistogram1D * totAPIXmatched = AIDAProcessor::histogramFactory(this)->createHistogram1D(_totAPIXmatchedHistoName.c_str(), 400, 0, 400);
		AIDA::IHistogram1D * totAPIXunmatched = AIDAProcessor::histogramFactory(this)->createHistogram1D(_totAPIXunmatchedHistoName.c_str(), 400, 0, 400);

		AIDA::IHistogram1D * lv1APIXmatched = AIDAProcessor::histogramFactory(this)->createHistogram1D(_lv1APIXmatchedHistoName.c_str(), 20, 0, 20);
		AIDA::IHistogram1D * lv1APIXunmatched = AIDAProcessor::histogramFactory(this)->createHistogram1D(_lv1APIXunmatchedHistoName.c_str(), 20, 0, 20);

		AIDA::IHistogram1D * hitsInMatchedCluster = AIDAProcessor::histogramFactory(this)->createHistogram1D(_hitsInMatchedClusterHistoName.c_str(), 20, 0, 20);
		AIDA::IHistogram1D * hitsInUnmatchedCluster = AIDAProcessor::histogramFactory(this)->createHistogram1D(_hitsInUnmatchedClusterHistoName.c_str(), 20, 0, 20);

		AIDA::IHistogram1D * maxDifflv1Matched = AIDAProcessor::histogramFactory(this)->createHistogram1D(_maxDifflv1MatchedHistoName.c_str(), 20, 0, 20);
		AIDA::IHistogram1D * maxDifflv1Unmatched = AIDAProcessor::histogramFactory(this)->createHistogram1D(_maxDifflv1UnmatchedHistoName.c_str(), 20, 0, 20);


		_aidaHistoMap.insert(make_pair(_totAPIXmatchedHistoName, totAPIXmatched));
		_aidaHistoMap.insert(make_pair(_totAPIXunmatchedHistoName, totAPIXunmatched));
		_aidaHistoMap.insert(make_pair(_lv1APIXmatchedHistoName, lv1APIXmatched));
		_aidaHistoMap.insert(make_pair(_lv1APIXunmatchedHistoName, lv1APIXunmatched));
		_aidaHistoMap.insert(make_pair(_hitsInMatchedClusterHistoName, hitsInMatchedCluster));
		_aidaHistoMap.insert(make_pair(_hitsInUnmatchedClusterHistoName, hitsInUnmatchedCluster));
		_aidaHistoMap.insert(make_pair(_maxDifflv1MatchedHistoName, maxDifflv1Matched));
		_aidaHistoMap.insert(make_pair(_maxDifflv1UnmatchedHistoName, maxDifflv1Unmatched));

		 AIDA::IProfile1D * effiXLOCALHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_EfficiencyXLOCALHistoName.c_str(),18, 0, 7.2);
		 AIDA::IProfile1D * effiYLOCALHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_EfficiencyYLOCALHistoName.c_str(),160, 0, 8);
		 AIDA::IProfile2D * effiXYLOCALHisto = AIDAProcessor::histogramFactory(this)->createProfile2D(_EfficiencyXYLOCALHistoName.c_str(), 18, 0, 7.2, 160, 0, 8);

		_aidaHistoMap.insert(make_pair(_EfficiencyXLOCALHistoName, effiXLOCALHisto));
		_aidaHistoMap.insert(make_pair(_EfficiencyYLOCALHistoName, effiYLOCALHisto));
		_aidaHistoMap.insert(make_pair(_EfficiencyXYLOCALHistoName, effiXYLOCALHisto));

		// -- more DUT plots, 24 August 2010, libov@mail.desy.de
		AIDA::IHistogram1D * MeasuredXLOCALHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_MeasuredXLOCALHistoName.c_str(), 18, 0, 7.2);
		AIDA::IHistogram1D * MeasuredYLOCALHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_MeasuredYLOCALHistoName.c_str(), 160, 0, 8);
		AIDA::IHistogram2D * MeasuredXYLOCALHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D(_MeasuredXYLOCALHistoName.c_str(), 18, 0, 7.2, 160, 0, 8);

		AIDA::IHistogram1D * FittedXLOCALHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_FittedXLOCALHistoName.c_str(), 18, 0, 7.2);
		AIDA::IHistogram1D * FittedYLOCALHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_FittedYLOCALHistoName.c_str(), 160, 0, 8);
		AIDA::IHistogram2D * FittedXYLOCALHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D(_FittedXYLOCALHistoName.c_str(), 18, 0, 7.2, 160, 0, 8);

		AIDA::IHistogram1D * _MatchedClusterSizeXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_MatchedClusterSizeXHistoName.c_str(), 10, 0, 10);
		AIDA::IHistogram1D * _MatchedClusterSizeYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_MatchedClusterSizeYHistoName.c_str(), 10, 0, 10);
		AIDA::IHistogram1D * _UnmatchedClusterSizeXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_UnmatchedClusterSizeXHistoName.c_str(), 10, 0, 10);
		AIDA::IHistogram1D * _UnmatchedClusterSizeYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_UnmatchedClusterSizeYHistoName.c_str(), 10, 0, 10);

		AIDA::IProfile1D * _ChargeSharingProbXHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_ChargeSharingProbXHistoName.c_str(),100,0,400);
		AIDA::IProfile1D * _ChargeSharingProbYHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_ChargeSharingProbYHistoName.c_str(),16,0,50);
		AIDA::IProfile2D * _ChargeSharingProbXYHisto = AIDAProcessor::histogramFactory(this)->createProfile2D(_ChargeSharingProbXYHistoName.c_str(),100,0, 400,16,0,50);

		_aidaHistoMap.insert(make_pair(_MeasuredXLOCALHistoName, MeasuredXLOCALHisto));
		_aidaHistoMap.insert(make_pair(_MeasuredYLOCALHistoName, MeasuredYLOCALHisto));
		_aidaHistoMap.insert(make_pair(_MeasuredXYLOCALHistoName, MeasuredXYLOCALHisto));

		_aidaHistoMap.insert(make_pair(_FittedXLOCALHistoName, FittedXLOCALHisto));
		_aidaHistoMap.insert(make_pair(_FittedYLOCALHistoName, FittedYLOCALHisto));
		_aidaHistoMap.insert(make_pair(_FittedXYLOCALHistoName, FittedXYLOCALHisto));
		
		_aidaHistoMap.insert(make_pair(_MatchedClusterSizeXHistoName, _MatchedClusterSizeXHisto));
		_aidaHistoMap.insert(make_pair(_MatchedClusterSizeYHistoName, _MatchedClusterSizeYHisto));
		_aidaHistoMap.insert(make_pair(_UnmatchedClusterSizeXHistoName, _UnmatchedClusterSizeXHisto));
		_aidaHistoMap.insert(make_pair(_UnmatchedClusterSizeYHistoName, _UnmatchedClusterSizeYHisto));

		_aidaHistoMap.insert(make_pair(_ChargeSharingProbXHistoName, _ChargeSharingProbXHisto));
		_aidaHistoMap.insert(make_pair(_ChargeSharingProbYHistoName, _ChargeSharingProbYHisto));
		_aidaHistoMap.insert(make_pair(_ChargeSharingProbXYHistoName, _ChargeSharingProbXYHisto));

		// 06 september
		AIDA::IProfile1D * _totSinglePixelClustersXHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_totSinglePixelClustersXHistoName.c_str(),100,0,400);

		AIDA::IProfile1D * _totSinglePixelClustersYHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_totSinglePixelClustersYHistoName.c_str(),16,0,50);

		AIDA::IProfile2D * _totSinglePixelClustersXYHisto = AIDAProcessor::histogramFactory(this)->createProfile2D(_totSinglePixelClustersXYHistoName.c_str(),100, 0, 400, 16, 0, 50);

		_aidaHistoMap.insert(make_pair(_totSinglePixelClustersXHistoName, _totSinglePixelClustersXHisto));
		_aidaHistoMap.insert(make_pair(_totSinglePixelClustersYHistoName, _totSinglePixelClustersYHisto));
		_aidaHistoMap.insert(make_pair(_totSinglePixelClustersXYHistoName, _totSinglePixelClustersXYHisto));

AIDA::IProfile2D * _totAllClustersXYHisto = AIDAProcessor::histogramFactory(this)->createProfile2D(_totAllClustersXYHistoName.c_str(),100, 0, 400, 16, 0, 50);

_aidaHistoMap.insert(make_pair(_totAllClustersXYHistoName, _totAllClustersXYHisto));

	// 06 october
	AIDA::IHistogram1D * _MatchedHitsHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_MatchedHitsHistoName.c_str(), 20, 0, 20 );
	AIDA::IProfile1D * _MatchedHitsVSEventHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_MatchedHitsVSEventHistoName.c_str(), 100000, 0, 100000 );

	_aidaHistoMap.insert(make_pair(_MatchedHitsHistoName, _MatchedHitsHisto));
	_aidaHistoMap.insert(make_pair(_MatchedHitsVSEventHistoName, _MatchedHitsVSEventHisto));

AIDA::IHistogram1D * _NumberOfFittedTracksHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_NumberOfFittedTracksHistoName.c_str(), 20, 0, 20 );
	_aidaHistoMap.insert(make_pair(_NumberOfFittedTracksHistoName, _NumberOfFittedTracksHisto));

	AIDA::IHistogram1D * _Lv1OfFirstMatchedReferenceHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_Lv1OfFirstMatchedReferenceHistoName.c_str(), 20, 0, 20);
	_aidaHistoMap.insert(make_pair(_Lv1OfFirstMatchedReferenceHistoName, _Lv1OfFirstMatchedReferenceHisto));

	AIDA::IHistogram1D * _NumberOfMatchedReferencesHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_NumberOfMatchedReferencesHistoName.c_str(), 20, 0, 20);
	_aidaHistoMap.insert(make_pair(_NumberOfMatchedReferencesHistoName, _NumberOfMatchedReferencesHisto));

	}

// List all booked histogram - check of histogram map filling



  message<MESSAGE5> ( log() <<  _aidaHistoMap.size() << " histograms booked");


  map<string, AIDA::IBaseHistogram *>::iterator mapIter;
  for(mapIter = _aidaHistoMap.begin(); mapIter != _aidaHistoMap.end() ; mapIter++ )
    message<DEBUG5> ( log() <<  mapIter->first << " : " <<  (mapIter->second)->title() ) ;

  message<DEBUG5> ( log() << "Histogram booking completed \n\n");
#else
  message<MESSAGE5> ( log() << "No histogram produced because Marlin doesn't use AIDA" );
#endif

  return;
}

#endif

void EUTelAPIXHistograms::getTransformationShifts() {

    if( _manualDUTid < 10 ) return;

    _transShiftX = (a * _measuredXLOCAL[0] + b * _measuredYLOCAL[0]) - _measuredX[0];
	_transShiftY = (c * _measuredXLOCAL[0] + d * _measuredYLOCAL[0]) - _measuredY[0];
}

bool EUTelAPIXHistograms::hasMatchedHit( Track* fittrack ) {
		std::vector<EVENT::TrackerHit*>  trackhits = fittrack->getTrackerHits();

		// loop over the hits, find fitted ones in the reference plane
		std::map<int, float> allDUtsfittedX, allDUtsfittedY;

		int nHit =   trackhits.size();
		for(int ihit=0; ihit < nHit ; ihit++) {
			TrackerHit * meshit = trackhits.at(ihit);
			if (meshit->getType() < 32) continue;
			// Hit position
			const double * pos = meshit -> getPosition();

			for(int ipl=0; ipl <  _siPlanesLayerLayout->getNLayers(); ipl++) {
				double dist =  pos[2] - ( _siPlanesLayerLayout -> getLayerPositionZ(ipl) ) ;
				if( dist*dist < 1 ) {
					allDUtsfittedX [ _siPlanesLayerLayout -> getID(ipl) ] = pos[0];
					allDUtsfittedY [ _siPlanesLayerLayout -> getID(ipl) ] = pos[1];
				}
			}
		}
		bool	matchedTrackToReferencePlane = false;
		int	nMatches=0;
		float	Lv1ofFirstMatch = -1;
		for (size_t i = 0; i< _allDUTsmeasuredX.size(); i++) {
			double deltaX = _allDUTsmeasuredX[i] - allDUtsfittedX[ _allDUTssensorID[i] ];
			double deltaY = _allDUTsmeasuredY[i] - allDUtsfittedY[ _allDUTssensorID[i] ];

			double	dist = deltaX * deltaX + deltaY * deltaY;
			if ( (dist < _distMaxReference*_distMaxReference) ) {
				nMatches++;
				if ( Lv1ofFirstMatch < 0 ) Lv1ofFirstMatch = _allDUTslv1[ i ];
			}
		}

		(dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_NumberOfMatchedReferencesHistoName] ) )->fill(nMatches);
		if ( nMatches > 0 ) (dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[_Lv1OfFirstMatchedReferenceHistoName] ))->fill(Lv1ofFirstMatch);

		if ( (nMatches > 0) && (Lv1ofFirstMatch >= 3 ) && (Lv1ofFirstMatch <= 7 ) ) matchedTrackToReferencePlane = true;

		return matchedTrackToReferencePlane;
}

void EUTelAPIXHistograms::fillAPIXhits(LCCollection* hitcol, LCCollection* original_zsdata) {
	// -- 02 October 2010
	CellIDDecoder<TrackerDataImpl> cellDecoder( original_zsdata );
	// fill measured hits for DUTs except DUT under consideration - for intime requirement
	_allDUTsmeasuredX.clear();
	_allDUTsmeasuredY.clear();
	_allDUTssensorID.clear();
	_allDUTslv1.clear();

	for(int ihit=0; ihit< hitcol->getNumberOfElements(); ihit++) {
		TrackerHit * meshit = dynamic_cast<TrackerHit*>( hitcol->getElementAt(ihit) ) ;
		LCObjectVec clusterVec = (meshit->getRawHits());
		TrackerDataImpl * trackerData = dynamic_cast < TrackerDataImpl * > ( clusterVec[0] );

		int sensorID             = static_cast<int > ( cellDecoder( trackerData )["sensorID"] );

		if ( (sensorID >=10) && (_manualDUTid >= 10) && (_manualDUTid != sensorID) ) { // only when DUT=apix, current sensorID is apix but not one
																												// under eff. measurement
				// Hit position
				const double * pos = meshit->getPosition();
				_allDUTsmeasuredX.push_back(pos[0]);
				_allDUTsmeasuredY.push_back(pos[1]);
				_allDUTssensorID.push_back(sensorID);

				EUTelSparseDataImpl <EUTelAPIXSparsePixel>  * sparseData = new EUTelSparseDataImpl <EUTelAPIXSparsePixel> (trackerData);
				EUTelAPIXSparsePixel * sparsePixel= new EUTelAPIXSparsePixel;
				sparseData -> getSparsePixelAt (0, sparsePixel );
				_allDUTslv1.push_back( sparsePixel -> getTime() );
				delete	sparseData;
				delete	sparsePixel;
		}
	}
	//--

}

void EUTelAPIXHistograms::getTrackImpactPoint(double & x, double & y, double & z, Track * tr, LCEvent * ev) {

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

	double	dist;
	for(int ihit=0; ihit< nHit ; ihit++)
	{
		TrackerHit * meshit = trackhits.at(ihit);

		// Look at fitted hits only!
		if (meshit->getType() < 32) continue;
		// Hit position
		const double * pos = meshit->getPosition();

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
	NormalVector(0) = sin( (-1)*_beta );
	NormalVector(1) = 0;
	NormalVector(2) = cos( (-1)*_beta );

	double	z_sensor = _siPlanesLayerLayout->getSensitivePositionZ(_indexDUT)+ 0.5 * _siPlanesLayerLayout->getSensitiveThickness( _indexDUT );

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

	x = trackImpact(0);
	y = trackImpact(1);
	z = trackImpact(2);

}




void EUTelAPIXHistograms::TransformToLocalFrame(double & x, double & y, double & z, LCEvent * ev) {

				// revert alignment, in an inverse order...
				for ( int i = _alignmentCollectionNames.size() - 1; i >= 0; i--) {
					revertAlignment (x, y, z, _alignmentCollectionNames[i], ev );
				}
				// revert beta rotations implememted in the hitmaker
				x = x / cos( _beta );
				z = z - (-1) * (-1) * x  * sin ( _beta );
				//cout << "z-z_sensor= " << z << endl;

				// revert setting (x,y) = (0,0) at the center of the sensor to the (row,col) = (0,0)
				double	xSize = _siPlanesLayerLayout->getSensitiveSizeX(_indexDUT);  // mm
				double	ySize = _siPlanesLayerLayout->getSensitiveSizeY(_indexDUT);  // mm
				// as in the hitmaker... -------
				double xPointing[2], yPointing[2];
				xPointing[0] = _siPlanesLayerLayout->getSensitiveRotation1(_indexDUT); // was -1 ;
				xPointing[1] = _siPlanesLayerLayout->getSensitiveRotation2(_indexDUT); // was  0 ;
				yPointing[0] = _siPlanesLayerLayout->getSensitiveRotation3(_indexDUT); // was  0 ;
				yPointing[1] = _siPlanesLayerLayout->getSensitiveRotation4(_indexDUT); // was -1 ;

				double sign = 0;
				if      ( xPointing[0] < -0.7 )       sign = -1 ;
				else if ( xPointing[0] > 0.7 )       sign =  1 ;
				else {
				if       ( xPointing[1] < -0.7 )    sign = -1 ;
				else if  ( xPointing[1] > 0.7 )    sign =  1 ;
				}
				x += sign * xSize/2;

				if      ( yPointing[0] < -0.7 )       sign = -1 ;
				else if ( yPointing[0] > 0.7 )       sign =  1 ;
				else {
				if       ( yPointing[1] < -0.7 )    sign = -1 ;
				else if  ( yPointing[1] > 0.7 )    sign =  1 ;
				}
				y += sign * ySize/2;
				//--------------

				// revert gear rotations
				double	x_temp = x;
				double	y_temp = y;

				x = _rot00 * x_temp + _rot01 * y_temp;
				y = _rot10 * x_temp + _rot11 * y_temp;
}

void EUTelAPIXHistograms::revertAlignment(double & x, double & y, double & z, std::string	collectionName, LCEvent * ev) {

	// in this function, some parts of the EUTelApplyAlignmentProcessor are used

	// get the alignment constant object
	// first, get the alignment collection
	LCCollectionVec * alignmentCollectionVec     = dynamic_cast < LCCollectionVec * > (ev->getCollection(collectionName));
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

	// now apply the constants

	// the way we apply constants is correctionMethod1 (first the rotations, second the shifts)
	// to revert the alignment transformation properly, the constants have to be applied in
	// a reverted way, i.e. first the shifts, second the rotations

	// note that the sign is different for all the constants with respect to what is done in EUTelApplyAlignmentProcessor -
	// the transformation is reverted

	// first the shifts
	x += c->getXOffset();
	y += c->getYOffset();
	z += c->getZOffset();

	double	x_temp = x;
	double	y_temp = y;
	double	z_temp = z;

	double	alpha = c -> getAlpha();
	double	beta = c -> getBeta();
	double	gamma = c -> getGamma();

	// second the rotation
	// for the inverse matrix derivation see paper log book 19/01/2011
	// libov@mail.desy.de

	x = x_temp * (1 + alpha * alpha ) + ( (-1) * gamma - alpha * beta) * y_temp + ( (-1) * beta + alpha * gamma) * z_temp;
	y = x_temp * (gamma - alpha * beta ) + (1 + beta * beta) * y_temp + ((-1) * alpha - beta * gamma) * z_temp;
	z = x_temp * (beta + alpha * gamma ) + (alpha - gamma * beta) * y_temp + ( 1 + gamma * gamma) * z_temp;

	double det = 1 + alpha * alpha + beta * beta + gamma * gamma;

	x = x / det;
	y = y / det;
	z = z / det;
}



