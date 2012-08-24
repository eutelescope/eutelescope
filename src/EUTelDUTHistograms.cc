// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-

// Author: A.F.Zarnecki, University of Warsaw <mailto:zarnecki@fuw.edu.pl>
// @version: $Id: EUTelDUTHistograms.cc,v 1.15 2009-07-29 11:05:02 bulgheroni Exp $

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifdef USE_GEAR

// ROOT includes:
#include "TVector3.h"

// eutelescope inlcudes
#include "EUTelDUTHistograms.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelHistogramManager.h"
#include "EUTelExceptions.h"
#include "EUTelReferenceHit.h"
// for cluster operations:
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelSparseCluster2Impl.h"
#
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
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

std::string EUTelDUTHistograms::_ClusterSizeXHistoName  = "clusterSizeX";
std::string EUTelDUTHistograms::_ClusterSizeYHistoName  = "clusterSizeY";
std::string EUTelDUTHistograms::_ClusterSizeXYHistoName = "clusterSuzeXY";

std::string EUTelDUTHistograms::_ClusterSizeXAHistoName  = "clusterSizeX submatrix A";
std::string EUTelDUTHistograms::_ClusterSizeYAHistoName  = "clusterSizeY submatrix A";
std::string EUTelDUTHistograms::_ClusterSizeXYAHistoName = "clusterSuzeXY submatrix A";

std::string EUTelDUTHistograms::_ClusterSizeXBHistoName  = "clusterSizeX submatrix B";
std::string EUTelDUTHistograms::_ClusterSizeYBHistoName  = "clusterSizeY submatrix B";
std::string EUTelDUTHistograms::_ClusterSizeXYBHistoName = "clusterSuzeXY submatrix B";

std::string EUTelDUTHistograms::_ClusterSizeXCHistoName  = "clusterSizeX submatrix C";
std::string EUTelDUTHistograms::_ClusterSizeYCHistoName  = "clusterSizeY submatrix C";
std::string EUTelDUTHistograms::_ClusterSizeXYCHistoName = "clusterSuzeXY submatrix C";

std::string EUTelDUTHistograms::_ClusterSizeXDHistoName  = "clusterSizeX submatrix D";
std::string EUTelDUTHistograms::_ClusterSizeYDHistoName  = "clusterSizeY submatrix D";
std::string EUTelDUTHistograms::_ClusterSizeXYDHistoName = "clusterSuzeXY submatrix D";


std::string EUTelDUTHistograms::_MeasuredXHistoName  = "measuredX";
std::string EUTelDUTHistograms::_MeasuredYHistoName  = "measuredY";
std::string EUTelDUTHistograms::_MeasuredXYHistoName = "measuredXY";

std::string EUTelDUTHistograms::_MatchedXHistoName  = "matchedX";
std::string EUTelDUTHistograms::_MatchedYHistoName  = "matchedY";
std::string EUTelDUTHistograms::_MatchedXYHistoName = "matchedXY";

std::string EUTelDUTHistograms::_UnMatchedXHistoName  = "unmatchedX";
std::string EUTelDUTHistograms::_UnMatchedYHistoName  = "unmatchedY";
std::string EUTelDUTHistograms::_UnMatchedXYHistoName = "unmatchedXY";

std::string EUTelDUTHistograms::_FittedXHistoName  = "fittedX";
std::string EUTelDUTHistograms::_FittedYHistoName  = "fittedY";
std::string EUTelDUTHistograms::_FittedXYHistoName = "fittedXY";

std::string EUTelDUTHistograms::_EfficiencyXHistoName  = "DUTeffiX";
std::string EUTelDUTHistograms::_EfficiencyYHistoName  = "DUTeffiY";
std::string EUTelDUTHistograms::_EfficiencyXYHistoName = "DUTeffiXY";

std::string EUTelDUTHistograms::_BgEfficiencyXHistoName  = "BGeffiX";
std::string EUTelDUTHistograms::_BgEfficiencyYHistoName  = "BGeffiY";
std::string EUTelDUTHistograms::_BgEfficiencyXYHistoName = "BGeffiXY";

std::string EUTelDUTHistograms::_NoiseXHistoName  = "DUTnoiseX";
std::string EUTelDUTHistograms::_NoiseYHistoName  = "DUTnoiseY";
std::string EUTelDUTHistograms::_NoiseXYHistoName = "DUTnoiseXY";

std::string EUTelDUTHistograms::_ShiftXHistoName       = "DUTshiftX";
std::string EUTelDUTHistograms::_ShiftYHistoName       = "DUTshiftY";
std::string EUTelDUTHistograms::_ShiftXYHistoName      = "DUTshiftXY";
// submatrix A
std::string EUTelDUTHistograms::_ShiftXAHistoName       = "DUTshiftXA";
std::string EUTelDUTHistograms::_ShiftYAHistoName       = "DUTshiftYA";
std::string EUTelDUTHistograms::_ShiftXYAHistoName      = "DUTshiftXYA";
// cluster size 1
std::string EUTelDUTHistograms::_ShiftXA1HistoName       = "DUTshiftXA1";
std::string EUTelDUTHistograms::_ShiftYA1HistoName       = "DUTshiftYA1";
std::string EUTelDUTHistograms::_ShiftXYA1HistoName      = "DUTshiftXYA1";
// cluster size 2
std::string EUTelDUTHistograms::_ShiftXA2HistoName       = "DUTshiftXA2";
std::string EUTelDUTHistograms::_ShiftYA2HistoName       = "DUTshiftYA2";
std::string EUTelDUTHistograms::_ShiftXYA2HistoName      = "DUTshiftXYA2";
// cluster size 3
std::string EUTelDUTHistograms::_ShiftXA3HistoName       = "DUTshiftXA3";
std::string EUTelDUTHistograms::_ShiftYA3HistoName       = "DUTshiftYA3";
std::string EUTelDUTHistograms::_ShiftXYA3HistoName      = "DUTshiftXYA3";
// cluster size 4
std::string EUTelDUTHistograms::_ShiftXA4HistoName       = "DUTshiftXA4";
std::string EUTelDUTHistograms::_ShiftYA4HistoName       = "DUTshiftYA4";
std::string EUTelDUTHistograms::_ShiftXYA4HistoName      = "DUTshiftXYA4";
// cluster size 5
std::string EUTelDUTHistograms::_ShiftXA5HistoName       = "DUTshiftXA5";
std::string EUTelDUTHistograms::_ShiftYA5HistoName       = "DUTshiftYA5";
std::string EUTelDUTHistograms::_ShiftXYA5HistoName      = "DUTshiftXYA5";
// cluster size 6
std::string EUTelDUTHistograms::_ShiftXA6HistoName       = "DUTshiftXA6";
std::string EUTelDUTHistograms::_ShiftYA6HistoName       = "DUTshiftYA6";
std::string EUTelDUTHistograms::_ShiftXYA6HistoName      = "DUTshiftXYA6";
// cluster size 7
std::string EUTelDUTHistograms::_ShiftXA7HistoName       = "DUTshiftXA7";
std::string EUTelDUTHistograms::_ShiftYA7HistoName       = "DUTshiftYA7";
std::string EUTelDUTHistograms::_ShiftXYA7HistoName      = "DUTshiftXYA7";

// submatrix B
std::string EUTelDUTHistograms::_ShiftXBHistoName       = "DUTshiftXB";
std::string EUTelDUTHistograms::_ShiftYBHistoName       = "DUTshiftYB";
std::string EUTelDUTHistograms::_ShiftXYBHistoName      = "DUTshiftXYB";
// cluster size 1
std::string EUTelDUTHistograms::_ShiftXB1HistoName       = "DUTshiftXB1";
std::string EUTelDUTHistograms::_ShiftYB1HistoName       = "DUTshiftYB1";
std::string EUTelDUTHistograms::_ShiftXYB1HistoName      = "DUTshiftXYB1";
// cluster size 2
std::string EUTelDUTHistograms::_ShiftXB2HistoName       = "DUTshiftXB2";
std::string EUTelDUTHistograms::_ShiftYB2HistoName       = "DUTshiftYB2";
std::string EUTelDUTHistograms::_ShiftXYB2HistoName      = "DUTshiftXYB2";
// cluster size 3
std::string EUTelDUTHistograms::_ShiftXB3HistoName       = "DUTshiftXB3";
std::string EUTelDUTHistograms::_ShiftYB3HistoName       = "DUTshiftYB3";
std::string EUTelDUTHistograms::_ShiftXYB3HistoName      = "DUTshiftXYB3";
// cluster size 4
std::string EUTelDUTHistograms::_ShiftXB4HistoName       = "DUTshiftXB4";
std::string EUTelDUTHistograms::_ShiftYB4HistoName       = "DUTshiftYB4";
std::string EUTelDUTHistograms::_ShiftXYB4HistoName      = "DUTshiftXYB4";
// cluster size 5
std::string EUTelDUTHistograms::_ShiftXB5HistoName       = "DUTshiftXB5";
std::string EUTelDUTHistograms::_ShiftYB5HistoName       = "DUTshiftYB5";
std::string EUTelDUTHistograms::_ShiftXYB5HistoName      = "DUTshiftXYB5";
// cluster size 6
std::string EUTelDUTHistograms::_ShiftXB6HistoName       = "DUTshiftXB6";
std::string EUTelDUTHistograms::_ShiftYB6HistoName       = "DUTshiftYB6";
std::string EUTelDUTHistograms::_ShiftXYB6HistoName      = "DUTshiftXYB6";
// cluster size 7
std::string EUTelDUTHistograms::_ShiftXB7HistoName       = "DUTshiftXB7";
std::string EUTelDUTHistograms::_ShiftYB7HistoName       = "DUTshiftYB7";
std::string EUTelDUTHistograms::_ShiftXYB7HistoName      = "DUTshiftXYB7";


// submatrix C
std::string EUTelDUTHistograms::_ShiftXCHistoName       = "DUTshiftXC";
std::string EUTelDUTHistograms::_ShiftYCHistoName       = "DUTshiftYC";
std::string EUTelDUTHistograms::_ShiftXYCHistoName      = "DUTshiftXYC";
// cluster size 1
std::string EUTelDUTHistograms::_ShiftXC1HistoName       = "DUTshiftXC1";
std::string EUTelDUTHistograms::_ShiftYC1HistoName       = "DUTshiftYC1";
std::string EUTelDUTHistograms::_ShiftXYC1HistoName      = "DUTshiftXYC1";
// cluster size 2
std::string EUTelDUTHistograms::_ShiftXC2HistoName       = "DUTshiftXC2";
std::string EUTelDUTHistograms::_ShiftYC2HistoName       = "DUTshiftYC2";
std::string EUTelDUTHistograms::_ShiftXYC2HistoName      = "DUTshiftXYC2";
// cluster size 3
std::string EUTelDUTHistograms::_ShiftXC3HistoName       = "DUTshiftXC3";
std::string EUTelDUTHistograms::_ShiftYC3HistoName       = "DUTshiftYC3";
std::string EUTelDUTHistograms::_ShiftXYC3HistoName      = "DUTshiftXYC3";
// cluster size 4
std::string EUTelDUTHistograms::_ShiftXC4HistoName       = "DUTshiftXC4";
std::string EUTelDUTHistograms::_ShiftYC4HistoName       = "DUTshiftYC4";
std::string EUTelDUTHistograms::_ShiftXYC4HistoName      = "DUTshiftXYC4";
// cluster size 5
std::string EUTelDUTHistograms::_ShiftXC5HistoName       = "DUTshiftXC5";
std::string EUTelDUTHistograms::_ShiftYC5HistoName       = "DUTshiftYC5";
std::string EUTelDUTHistograms::_ShiftXYC5HistoName      = "DUTshiftXYC5";
// cluster size 6
std::string EUTelDUTHistograms::_ShiftXC6HistoName       = "DUTshiftXC6";
std::string EUTelDUTHistograms::_ShiftYC6HistoName       = "DUTshiftYC6";
std::string EUTelDUTHistograms::_ShiftXYC6HistoName      = "DUTshiftXYC6";
// cluster size 7
std::string EUTelDUTHistograms::_ShiftXC7HistoName       = "DUTshiftXC7";
std::string EUTelDUTHistograms::_ShiftYC7HistoName       = "DUTshiftYC7";
std::string EUTelDUTHistograms::_ShiftXYC7HistoName      = "DUTshiftXYC7";


// submatrix D
std::string EUTelDUTHistograms::_ShiftXDHistoName       = "DUTshiftXD";
std::string EUTelDUTHistograms::_ShiftYDHistoName       = "DUTshiftYD";
std::string EUTelDUTHistograms::_ShiftXYDHistoName      = "DUTshiftXYD";
// cluster size 1
std::string EUTelDUTHistograms::_ShiftXD1HistoName       = "DUTshiftXD1";
std::string EUTelDUTHistograms::_ShiftYD1HistoName       = "DUTshiftYD1";
std::string EUTelDUTHistograms::_ShiftXYD1HistoName      = "DUTshiftXYD1";
// cluster size 2
std::string EUTelDUTHistograms::_ShiftXD2HistoName       = "DUTshiftXD2";
std::string EUTelDUTHistograms::_ShiftYD2HistoName       = "DUTshiftYD2";
std::string EUTelDUTHistograms::_ShiftXYD2HistoName      = "DUTshiftXYD2";
// cluster size 3
std::string EUTelDUTHistograms::_ShiftXD3HistoName       = "DUTshiftXD3";
std::string EUTelDUTHistograms::_ShiftYD3HistoName       = "DUTshiftYD3";
std::string EUTelDUTHistograms::_ShiftXYD3HistoName      = "DUTshiftXYD3";
// cluster size 4
std::string EUTelDUTHistograms::_ShiftXD4HistoName       = "DUTshiftXD4";
std::string EUTelDUTHistograms::_ShiftYD4HistoName       = "DUTshiftYD4";
std::string EUTelDUTHistograms::_ShiftXYD4HistoName      = "DUTshiftXYD4";
// cluster size 5
std::string EUTelDUTHistograms::_ShiftXD5HistoName       = "DUTshiftXD5";
std::string EUTelDUTHistograms::_ShiftYD5HistoName       = "DUTshiftYD5";
std::string EUTelDUTHistograms::_ShiftXYD5HistoName      = "DUTshiftXYD5";
// cluster size 6
std::string EUTelDUTHistograms::_ShiftXD6HistoName       = "DUTshiftXD6";
std::string EUTelDUTHistograms::_ShiftYD6HistoName       = "DUTshiftYD6";
std::string EUTelDUTHistograms::_ShiftXYD6HistoName      = "DUTshiftXYD6";
// cluster size 7
std::string EUTelDUTHistograms::_ShiftXD7HistoName       = "DUTshiftXD7";
std::string EUTelDUTHistograms::_ShiftYD7HistoName       = "DUTshiftYD7";
std::string EUTelDUTHistograms::_ShiftXYD7HistoName      = "DUTshiftXYD7";







std::string EUTelDUTHistograms::_BgShiftXHistoName       = "BGshiftX";
std::string EUTelDUTHistograms::_BgShiftYHistoName       = "BGshiftY";
std::string EUTelDUTHistograms::_BgShiftXYHistoName      = "BGshiftXY";

std::string EUTelDUTHistograms::_ShiftXvsYHistoName    = "DUTshiftXvsY";
std::string EUTelDUTHistograms::_ShiftYvsXHistoName    = "DUTshiftYvsX";
std::string EUTelDUTHistograms::_ShiftXvsY2DHistoName    = "DUTshiftXvsY2D";
std::string EUTelDUTHistograms::_ShiftYvsX2DHistoName    = "DUTshiftYvsX2D";

std::string EUTelDUTHistograms::_ShiftYvsYHistoName    = "DUTshiftYvsY";
std::string EUTelDUTHistograms::_ShiftXvsXHistoName    = "DUTshiftXvsX";
std::string EUTelDUTHistograms::_ShiftXvsX2DHistoName    = "DUTshiftXvsX2D";
std::string EUTelDUTHistograms::_ShiftYvsY2DHistoName    = "DUTshiftYvsY2D";

std::string EUTelDUTHistograms::_EtaXHistoName      = "EtaX";
std::string EUTelDUTHistograms::_EtaYHistoName      = "EtaY";
std::string EUTelDUTHistograms::_EtaX2DHistoName    = "EtaX2D";
std::string EUTelDUTHistograms::_EtaY2DHistoName    = "EtaY2D";
std::string EUTelDUTHistograms::_EtaX3DHistoName    = "EtaX3D";
std::string EUTelDUTHistograms::_EtaY3DHistoName    = "EtaY3D";

#endif


EUTelDUTHistograms::EUTelDUTHistograms() : Processor("EUTelDUTHistograms") {

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


  registerOptionalParameter("ReferenceCollection","reference hit collection name ", _referenceHitCollectionName, static_cast <string> ("reference") );
 
  registerOptionalParameter("ApplyToReferenceCollection","Do you want the reference hit collection to be corrected by the shifts and tilts from the alignment collection? (default - false )",  _applyToReferenceHitCollection, static_cast< bool   > ( false ));


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


}


void EUTelDUTHistograms::init() {

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  _referenceHitVec = 0;

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

  if(_useManualDUT)
    {
      bool _manualOK=false;

      for(int ipl=0; ipl <  _siPlanesLayerLayout->getNLayers(); ipl++)
        if(_siPlanesLayerLayout->getID(ipl)==_manualDUTid)
          {
            _iDUT=_manualDUTid;
            _zDUT=_siPlanesLayerLayout->getLayerPositionZ(ipl);
            _manualOK=true;
          }

      if(!_manualOK)
        {
          message<ERROR> ( log() << "Manual DUT flag set, layer not found ID = "
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
        message<ERROR> ( log() << "DUT analysis initialized, but no DUT found in GEAR \n"
                         << "Program will terminate! Correct geometry description!");
        exit(-1);
      }


// Print out geometry information

  message<MESSAGE> ( log() << "D.U.T. plane  ID = " << _iDUT
                     << "  at Z [mm] = " << _zDUT );



// Book histograms

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  bookHistos();
#endif

}

void EUTelDUTHistograms::processRunHeader( LCRunHeader* runHeader) {

  auto_ptr<EUTelRunHeaderImpl> eutelHeader( new EUTelRunHeaderImpl ( runHeader ) );
  eutelHeader->addProcessor( type() );

  _nRun++ ;

  // Decode and print out Run Header information - just a check

  int runNr = runHeader->getRunNumber();

  message<MESSAGE> ( log() << "Processing run header " << _nRun
                     << ", run nr " << runNr );

  const std::string detectorName = runHeader->getDetectorName();
  const std::string detectorDescription = runHeader->getDescription();
  //  const std::vector<std::string> * subDets = runHeader->getActiveSubdetectors();

  message<MESSAGE> ( log() << detectorName << " : " << detectorDescription ) ;


}

void EUTelDUTHistograms::processEvent( LCEvent * event ) {

  if ( isFirstEvent() )
  {
    if ( _applyToReferenceHitCollection ) 
    {
       _referenceHitVec = dynamic_cast < LCCollectionVec * > (event->getCollection( _referenceHitCollectionName));
    }
  }

  EUTelEventImpl * euEvent = static_cast<EUTelEventImpl*> ( event );
  if ( euEvent->getEventType() == kEORE ) {
    message<DEBUG> ( "EORE found: nothing else to do." );
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


  if(debug)message<MESSAGE> ( log() << "Processing record " << _nEvt << " == event " << evtNr );

  //
  // Get input collections
  //

  LCCollection* trackcol;
  try {
    trackcol = event->getCollection( _inputTrackColName ) ;
  } catch (lcio::DataNotAvailableException& e) {
    message<ERROR> ( log() << "Not able to get collection "
                     << _inputTrackColName
                     << "\nfrom event " << event->getEventNumber()
                     << " in run " << event->getRunNumber()  );
    return;
//    throw SkipEventException(this);
  }


  LCCollection* hitcol = NULL;
  bool _DUTok=true;
  try {
    hitcol = event->getCollection( _inputHitColName ) ;
  } catch (lcio::DataNotAvailableException& e) {
    message<ERROR> ( log() << "Not able to get collection "
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

  // Clear local fit storage tables
  // 'fitted' table used for comparison with measured hits
  // 'bgfitted' used for comparison with hits from previous event

  _fittedX.clear();
  _fittedY.clear();

  _bgfittedX.clear();
  _bgfittedY.clear();


  // Loop over tracks in input track collection
  // Read fitted positions at DUT

  int nTrack = trackcol->getNumberOfElements()  ;


  if(debug)message<MESSAGE> ( log() << "Total of " << nTrack << " tracks in input collection " );


  for(int itrack=0; itrack< nTrack ; itrack++)
    {
      Track * fittrack = dynamic_cast<Track*>( trackcol->getElementAt(itrack) ) ;


      const float * por = fittrack->getReferencePoint();

      // Does track PoR match DUT position?

      double dist = por[2] - _zDUT ;

      if(dist*dist < 1)
        {
          _fittedX.push_back(por[0]);
          _fittedY.push_back(por[1]);
          _bgfittedX.push_back(por[0]);
          _bgfittedY.push_back(por[1]);
        }
      else
        {

          // Look at hits asigned to track

          std::vector<EVENT::TrackerHit*>  trackhits = fittrack->getTrackerHits();

          int nHit =   trackhits.size();

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
                  _bgfittedX.push_back(pos[0]);
                  _bgfittedY.push_back(pos[1]);
                  break;
                }
            }
        }

      // End of loop over fitted tracks
    }

  if(debug)message<MESSAGE> ( log() << _fittedX.size()  << " fitted positions at DUT " );

  // Histograms of fitted positions



  for(int ifit=0;ifit<(int)_fittedX.size(); ifit++)
    {
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_FittedXHistoName]))->fill(_fittedX[ifit]);


      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_FittedYHistoName]))->fill(_fittedY[ifit]);

      (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_FittedXYHistoName]))->fill(_fittedX[ifit],_fittedY[ifit]);


      if(debug)message<MESSAGE> ( log() << "Fit " << ifit
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

        for(int ifit=0;ifit<(int)_bgfittedX.size(); ifit++)
          for(int ihit=0; ihit< (int)_bgmeasuredX.size() ; ihit++)
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

      for(int ifit=0;ifit<(int)_bgfittedX.size(); ifit++)
        {
          (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_BgEfficiencyXHistoName]))->fill(_bgfittedX[ifit],0.);

          (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_BgEfficiencyYHistoName]))->fill(_bgfittedY[ifit],0.);

          (dynamic_cast<AIDA::IProfile2D*> ( _aidaHistoMap[_BgEfficiencyXYHistoName]))->fill(_bgfittedX[ifit],_bgfittedY[ifit],0.);
        }

#endif

    }  // End of if(_nEvt > 1)  - end of background calculations


  // Clear local tables with measured position
  _clusterSizeX.clear();
  _clusterSizeY.clear();
  _subMatrix.clear();


  _measuredX.clear();
  _measuredY.clear();

  _localX.clear();
  _localY.clear();

  _bgmeasuredX.clear();
  _bgmeasuredY.clear();

  // Loop over hits in input collection
  // Read measured positions at DUT

  int nHit = 0;

  if(_DUTok) nHit = hitcol->getNumberOfElements();

  if(debug)message<MESSAGE> ( log() << "Total of " << nHit << " tracker hits in input collection " );


  for(int ihit=0; ihit< nHit ; ihit++)
    {
      TrackerHit * meshit = dynamic_cast<TrackerHit*>( hitcol->getElementAt(ihit) ) ;

      // Hit position

      const double * pos = meshit->getPosition();

      int   sensorID = guessSensorID( pos );

      double dist = pos[2] - _zDUT;

      if(dist*dist < 1)
        {
          int sizeX = -1;
          int sizeY = -1;
          int subMatrix = -1;
          getClusterSize( sensorID, static_cast<TrackerHitImpl*>(meshit), sizeX, sizeY, subMatrix);
          //printf("%d %d \n", sizeX,sizeY);  
          _clusterSizeX.push_back(sizeX);
          _clusterSizeY.push_back(sizeY);
          _subMatrix.push_back( subMatrix );


          // Apply alignment corrections

          double corrX  =   pos[0]*cos(_DUTalign.at(2))
            +pos[1]*sin(_DUTalign.at(2))
            +_DUTalign.at(0);

          double corrY  =  pos[1]*cos(_DUTalign.at(2))
            -pos[0]*sin(_DUTalign.at(2))
            +_DUTalign.at(1);

          _measuredX.push_back(corrX);
          _measuredY.push_back(corrY);
          _bgmeasuredX.push_back(corrX);
          _bgmeasuredY.push_back(corrY);

          // Local position should be taken from the cluster.
          // This is a temporary solution:

          double locX = pos[0];
          double locY = pos[1];

          // Subtract position of the central pixel

          int picX = (int)(locX/_pitchX);

          if(locX<0)picX--;

          locX-=(picX+0.5)*_pitchX;

          int picY = (int)(locY/_pitchY);

          if(locY<0)picY--;

          locY-=(picY+0.5)*_pitchY;

          _localX.push_back(locX);
          _localY.push_back(locY);
        }
    }


  if(debug)message<MESSAGE> ( log() << _measuredX.size() << " hits at DUT " );


  // Histograms of measured positions


  for(int ihit=0;ihit<(int)_measuredX.size(); ihit++)
    {
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA) 
      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_MeasuredXHistoName]))->fill(_measuredX[ihit]);


      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_MeasuredYHistoName]))->fill(_measuredY[ihit]);

      (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_MeasuredXYHistoName]))->fill(_measuredX[ihit],_measuredY[ihit]);
#endif
      if(debug)message<MESSAGE> ( log() << "Hit " << ihit
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

    for(int ifit=0;ifit<(int)_fittedX.size(); ifit++)
      for(int ihit=0; ihit< (int)_measuredX.size() ; ihit++)
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

        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_ClusterSizeXHistoName]))->fill(_clusterSizeX[besthit]+0.0);
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_ClusterSizeYHistoName]))->fill(_clusterSizeY[besthit]+0.0);
        (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_ClusterSizeXYHistoName]))->fill(_clusterSizeX[besthit]+0.0,_clusterSizeY[besthit]+0.0);

      if( _subMatrix[besthit] == 0)
       { 
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_ClusterSizeXAHistoName]))->fill(_clusterSizeX[besthit]+0.0);
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_ClusterSizeYAHistoName]))->fill(_clusterSizeY[besthit]+0.0);
        (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_ClusterSizeXYAHistoName]))->fill(_clusterSizeX[besthit]+0.0,_clusterSizeY[besthit]+0.0);
       }

      if( _subMatrix[besthit] == 1)
       { 
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_ClusterSizeXBHistoName]))->fill(_clusterSizeX[besthit]+0.0);
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_ClusterSizeYBHistoName]))->fill(_clusterSizeY[besthit]+0.0);
        (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_ClusterSizeXYBHistoName]))->fill(_clusterSizeX[besthit]+0.0,_clusterSizeY[besthit]+0.0);
       }

      if( _subMatrix[besthit] == 2)
       { 
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_ClusterSizeXCHistoName]))->fill(_clusterSizeX[besthit]+0.0);
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_ClusterSizeYCHistoName]))->fill(_clusterSizeY[besthit]+0.0);
        (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_ClusterSizeXYCHistoName]))->fill(_clusterSizeX[besthit]+0.0,_clusterSizeY[besthit]+0.0);
       }

      if( _subMatrix[besthit] == 3)
       { 
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_ClusterSizeXDHistoName]))->fill(_clusterSizeX[besthit]+0.0);
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_ClusterSizeYDHistoName]))->fill(_clusterSizeY[besthit]+0.0);
        (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_ClusterSizeXYDHistoName]))->fill(_clusterSizeX[besthit]+0.0,_clusterSizeY[besthit]+0.0);
       }



        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_MatchedXHistoName]))->fill(_measuredX[besthit]);


        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_MatchedYHistoName]))->fill(_measuredY[besthit]);

        (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_MatchedXYHistoName]))->fill(_measuredX[besthit],_measuredY[besthit]);

        // Histograms of measured-fitted shifts

        double shiftX =  _measuredX[besthit]-_fittedX[bestfit];
        double shiftY =  _measuredY[besthit]-_fittedY[bestfit];

        std::string histoX  = _ShiftXHistoName;
        std::string histoY  = _ShiftYHistoName;
        std::string histoXY = _ShiftXYHistoName;

// fill anyway the global one:
          (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[ histoX ]))->fill(shiftX);
          (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[ histoY ]))->fill(shiftY);
          (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[ histoXY ]))->fill(shiftX, shiftY);
        
        if( _subMatrix[besthit] == 0 )  
        { 
          histoX = _ShiftXAHistoName;
          histoY = _ShiftYAHistoName;
          histoXY= _ShiftXYAHistoName;
        }
        if( _subMatrix[besthit] == 1 )  
        { 
          histoX = _ShiftXBHistoName;
          histoY = _ShiftYBHistoName;
          histoXY= _ShiftXYBHistoName;
        }
        if( _subMatrix[besthit] == 2 )  
        { 
          histoX = _ShiftXCHistoName;
          histoY = _ShiftYCHistoName;
          histoXY= _ShiftXYCHistoName;
        }
        if( _subMatrix[besthit] == 3 )  
        { 
          histoX = _ShiftXDHistoName;
          histoY = _ShiftYDHistoName;
          histoXY= _ShiftXYDHistoName;
        }
         (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[ histoX ]))->fill(shiftX);
         (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[ histoY ]))->fill(shiftY);
         (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[ histoXY ]))->fill(shiftX, shiftY);
        
        if( _subMatrix[besthit] == 0 )  
        {
          if( _clusterSizeX[besthit] == 1 )  histoX = _ShiftXA1HistoName;
          if( _clusterSizeY[besthit] == 1 )  histoY = _ShiftYA1HistoName;
          if( _clusterSizeX[besthit] == 1 && _clusterSizeY[besthit] == 1 )  histoXY= _ShiftXYA1HistoName;

          if( _clusterSizeX[besthit] == 2 )  histoX = _ShiftXA2HistoName;
          if( _clusterSizeY[besthit] == 2 )  histoY = _ShiftYA2HistoName;
          if( _clusterSizeX[besthit] == 2 && _clusterSizeY[besthit] == 2 )  histoXY= _ShiftXYA2HistoName;
 
          if( _clusterSizeX[besthit] == 3 )  histoX = _ShiftXA3HistoName;
          if( _clusterSizeY[besthit] == 3 )  histoY = _ShiftYA3HistoName;
          if( _clusterSizeX[besthit] == 3 && _clusterSizeY[besthit] == 3 )  histoXY= _ShiftXYA3HistoName;
 
          if( _clusterSizeX[besthit] == 4 )  histoX = _ShiftXA4HistoName;
          if( _clusterSizeY[besthit] == 4 )  histoY = _ShiftYA4HistoName;
          if( _clusterSizeX[besthit] == 4 && _clusterSizeY[besthit] == 4 )  histoXY= _ShiftXYA4HistoName;
 
          if( _clusterSizeX[besthit] == 5 )  histoX = _ShiftXA5HistoName;
          if( _clusterSizeY[besthit] == 5 )  histoY = _ShiftYA5HistoName;
          if( _clusterSizeX[besthit] == 5 && _clusterSizeY[besthit] == 5 )  histoXY= _ShiftXYA5HistoName;
 
          if( _clusterSizeX[besthit] == 6 )  histoX = _ShiftXA6HistoName;
          if( _clusterSizeY[besthit] == 6 )  histoY = _ShiftYA6HistoName;
          if( _clusterSizeX[besthit] == 6 && _clusterSizeY[besthit] == 6 )  histoXY= _ShiftXYA6HistoName;
 
          if( _clusterSizeX[besthit] == 7 )  histoX = _ShiftXA7HistoName;
          if( _clusterSizeY[besthit] == 7 )  histoY = _ShiftYA7HistoName;
          if( _clusterSizeX[besthit] == 7 && _clusterSizeY[besthit] == 7 )  histoXY= _ShiftXYA7HistoName;
        }
        if( _subMatrix[besthit] == 1 )  
        { 
           if( _clusterSizeX[besthit] == 1 )  histoX = _ShiftXB1HistoName;
          if( _clusterSizeY[besthit] == 1 )  histoY = _ShiftYB1HistoName;
          if( _clusterSizeX[besthit] == 1 && _clusterSizeY[besthit] == 1 )  histoXY= _ShiftXYB1HistoName;

          if( _clusterSizeX[besthit] == 2 )  histoX = _ShiftXB2HistoName;
          if( _clusterSizeY[besthit] == 2 )  histoY = _ShiftYB2HistoName;
          if( _clusterSizeX[besthit] == 2 && _clusterSizeY[besthit] == 2 )  histoXY= _ShiftXYB2HistoName;
 
          if( _clusterSizeX[besthit] == 3 )  histoX = _ShiftXB3HistoName;
          if( _clusterSizeY[besthit] == 3 )  histoY = _ShiftYB3HistoName;
          if( _clusterSizeX[besthit] == 3 && _clusterSizeY[besthit] == 3 )  histoXY= _ShiftXYB3HistoName;
 
          if( _clusterSizeX[besthit] == 4 )  histoX = _ShiftXB4HistoName;
          if( _clusterSizeY[besthit] == 4 )  histoY = _ShiftYB4HistoName;
          if( _clusterSizeX[besthit] == 4 && _clusterSizeY[besthit] == 4 )  histoXY= _ShiftXYB4HistoName;
 
          if( _clusterSizeX[besthit] == 5 )  histoX = _ShiftXB5HistoName;
          if( _clusterSizeY[besthit] == 5 )  histoY = _ShiftYB5HistoName;
          if( _clusterSizeX[besthit] == 5 && _clusterSizeY[besthit] == 5 )  histoXY= _ShiftXYB5HistoName;
 
          if( _clusterSizeX[besthit] == 6 )  histoX = _ShiftXB6HistoName;
          if( _clusterSizeY[besthit] == 6 )  histoY = _ShiftYB6HistoName;
          if( _clusterSizeX[besthit] == 6 && _clusterSizeY[besthit] == 6 )  histoXY= _ShiftXYB6HistoName;
 
          if( _clusterSizeX[besthit] == 7 )  histoX = _ShiftXB7HistoName;
          if( _clusterSizeY[besthit] == 7 )  histoY = _ShiftYB7HistoName;
          if( _clusterSizeX[besthit] == 7 && _clusterSizeY[besthit] == 7 )  histoXY= _ShiftXYB7HistoName;
       }
        if( _subMatrix[besthit] == 2 )  
        { 
           if( _clusterSizeX[besthit] == 1 )  histoX = _ShiftXC1HistoName;
          if( _clusterSizeY[besthit] == 1 )  histoY = _ShiftYC1HistoName;
          if( _clusterSizeX[besthit] == 1 && _clusterSizeY[besthit] == 1 )  histoXY= _ShiftXYC1HistoName;

          if( _clusterSizeX[besthit] == 2 )  histoX = _ShiftXC2HistoName;
          if( _clusterSizeY[besthit] == 2 )  histoY = _ShiftYC2HistoName;
          if( _clusterSizeX[besthit] == 2 && _clusterSizeY[besthit] == 2 )  histoXY= _ShiftXYC2HistoName;
 
          if( _clusterSizeX[besthit] == 3 )  histoX = _ShiftXC3HistoName;
          if( _clusterSizeY[besthit] == 3 )  histoY = _ShiftYC3HistoName;
          if( _clusterSizeX[besthit] == 3 && _clusterSizeY[besthit] == 3 )  histoXY= _ShiftXYC3HistoName;
 
          if( _clusterSizeX[besthit] == 4 )  histoX = _ShiftXC4HistoName;
          if( _clusterSizeY[besthit] == 4 )  histoY = _ShiftYC4HistoName;
          if( _clusterSizeX[besthit] == 4 && _clusterSizeY[besthit] == 4 )  histoXY= _ShiftXYC4HistoName;
 
          if( _clusterSizeX[besthit] == 5 )  histoX = _ShiftXC5HistoName;
          if( _clusterSizeY[besthit] == 5 )  histoY = _ShiftYC5HistoName;
          if( _clusterSizeX[besthit] == 5 && _clusterSizeY[besthit] == 5 )  histoXY= _ShiftXYC5HistoName;
 
          if( _clusterSizeX[besthit] == 6 )  histoX = _ShiftXC6HistoName;
          if( _clusterSizeY[besthit] == 6 )  histoY = _ShiftYC6HistoName;
          if( _clusterSizeX[besthit] == 6 && _clusterSizeY[besthit] == 6 )  histoXY= _ShiftXYC6HistoName;
 
          if( _clusterSizeX[besthit] == 7 )  histoX = _ShiftXC7HistoName;
          if( _clusterSizeY[besthit] == 7 )  histoY = _ShiftYC7HistoName;
          if( _clusterSizeX[besthit] == 7 && _clusterSizeY[besthit] == 7 )  histoXY= _ShiftXYC7HistoName;
       }
        if( _subMatrix[besthit] == 3 )  
        { 
           if( _clusterSizeX[besthit] == 1 )  histoX = _ShiftXD1HistoName;
          if( _clusterSizeY[besthit] == 1 )  histoY = _ShiftYD1HistoName;
          if( _clusterSizeX[besthit] == 1 && _clusterSizeY[besthit] == 1 )  histoXY= _ShiftXYD1HistoName;

          if( _clusterSizeX[besthit] == 2 )  histoX = _ShiftXD2HistoName;
          if( _clusterSizeY[besthit] == 2 )  histoY = _ShiftYD2HistoName;
          if( _clusterSizeX[besthit] == 2 && _clusterSizeY[besthit] == 2 )  histoXY= _ShiftXYD2HistoName;
 
          if( _clusterSizeX[besthit] == 3 )  histoX = _ShiftXD3HistoName;
          if( _clusterSizeY[besthit] == 3 )  histoY = _ShiftYD3HistoName;
          if( _clusterSizeX[besthit] == 3 && _clusterSizeY[besthit] == 3 )  histoXY= _ShiftXYD3HistoName;
 
          if( _clusterSizeX[besthit] == 4 )  histoX = _ShiftXD4HistoName;
          if( _clusterSizeY[besthit] == 4 )  histoY = _ShiftYD4HistoName;
          if( _clusterSizeX[besthit] == 4 && _clusterSizeY[besthit] == 4 )  histoXY= _ShiftXYD4HistoName;
 
          if( _clusterSizeX[besthit] == 5 )  histoX = _ShiftXD5HistoName;
          if( _clusterSizeY[besthit] == 5 )  histoY = _ShiftYD5HistoName;
          if( _clusterSizeX[besthit] == 5 && _clusterSizeY[besthit] == 5 )  histoXY= _ShiftXYD5HistoName;
 
          if( _clusterSizeX[besthit] == 6 )  histoX = _ShiftXD6HistoName;
          if( _clusterSizeY[besthit] == 6 )  histoY = _ShiftYD6HistoName;
          if( _clusterSizeX[besthit] == 6 && _clusterSizeY[besthit] == 6 )  histoXY= _ShiftXYD6HistoName;
 
          if( _clusterSizeX[besthit] == 7 )  histoX = _ShiftXD7HistoName;
          if( _clusterSizeY[besthit] == 7 )  histoY = _ShiftYD7HistoName;
          if( _clusterSizeX[besthit] == 7 && _clusterSizeY[besthit] == 7 )  histoXY= _ShiftXYD7HistoName;
       }
          (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[ histoX ]))->fill(shiftX);
          (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[ histoY ]))->fill(shiftY);
          (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[ histoXY ]))->fill(shiftX, shiftY);
        

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

        // Remove matched entries from the list (so the next matching pair
        // can be looked for)

        _fittedX.erase(_fittedX.begin()+bestfit);
        _fittedY.erase(_fittedY.begin()+bestfit);

        _measuredX.erase(_measuredX.begin()+besthit);
        _measuredY.erase(_measuredY.begin()+besthit);

        _localX.erase(_localX.begin()+besthit);
        _localY.erase(_localY.begin()+besthit);

      }


    // End of loop of matching DUT hits to fitted positions
  }
  while(distmin < _distMax*_distMax);


  if(debug)
    {
      message<MESSAGE> ( log() << nMatch << " DUT hits matched to fitted tracks ");
      message<MESSAGE> ( log() << _measuredX.size() << " DUT hits not matched to any track ");
      message<MESSAGE> ( log() << _fittedX.size() << " fitted tracks not matched to any DUT hit ");
    }


  // Efficiency plots - unmatched tracks

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  for(int ifit=0;ifit<(int)_fittedX.size(); ifit++)
    {
      (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_EfficiencyXHistoName]))->fill(_fittedX[ifit],0.);

      (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_EfficiencyYHistoName]))->fill(_fittedY[ifit],0.);

      (dynamic_cast<AIDA::IProfile2D*> ( _aidaHistoMap[_EfficiencyXYHistoName]))->fill(_fittedX[ifit],_fittedY[ifit],0.);
    }


  // Noise plots - unmatched hits

  for(int ihit=0;ihit<(int)_measuredX.size(); ihit++)
    {
      (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_NoiseXHistoName]))->fill(_measuredX[ihit],1.);

      (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_NoiseYHistoName]))->fill(_measuredY[ihit],1.);

      (dynamic_cast<AIDA::IProfile2D*> ( _aidaHistoMap[_NoiseXYHistoName]))->fill(_measuredX[ihit],_measuredY[ihit],1.);


      // Unmatched hit positions

      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_UnMatchedXHistoName]))->fill(_measuredX[ihit]);


      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_UnMatchedYHistoName]))->fill(_measuredY[ihit]);

      (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_UnMatchedXYHistoName]))->fill(_measuredX[ihit],_measuredY[ihit]);

    }

#endif

  if ( isFirstEvent() ) _isFirstEvent = false;

  return;
}



void EUTelDUTHistograms::check( LCEvent * /* evt */ ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void EUTelDUTHistograms::end(){

  // Nothing to do here
}



void EUTelDUTHistograms::bookHistos()
{

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

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


  // cluster size in X and Y
  int    clusXNBin  =  12 ;
  double clusXMin   = -1.5;
  double clusXMax   =  10.5;
 
  int    clusYNBin  =  12 ;
  double clusYMin   = -1.5;
  double clusYMax   =  10.5; 
  string clusXTitle = "Fitted cluster size in X"       ;
  string clusYTitle = "Fitted cluster size in Y"       ;
  string clusXYTitle = "Fitted cluster size in X and Y" ;

  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ClusterSizeXHistoName);
      if ( histoInfo )
        {
          message<DEBUG> ( log() << (* histoInfo ) );
          clusXNBin = histoInfo->_xBin;
          clusXMin  = histoInfo->_xMin;
          clusXMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) clusXTitle = histoInfo->_title;
        }
    }
  AIDA::IHistogram1D * clusXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ClusterSizeXHistoName.c_str(),clusXNBin,clusXMin,clusXMax);
  clusXHisto->setTitle(clusXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ClusterSizeXHistoName, clusXHisto));

  clusXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ClusterSizeXAHistoName.c_str(),clusXNBin,clusXMin,clusXMax);
  clusXHisto->setTitle(clusXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ClusterSizeXAHistoName, clusXHisto));

  clusXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ClusterSizeXBHistoName.c_str(),clusXNBin,clusXMin,clusXMax);
  clusXHisto->setTitle(clusXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ClusterSizeXBHistoName, clusXHisto));

  clusXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ClusterSizeXCHistoName.c_str(),clusXNBin,clusXMin,clusXMax);
  clusXHisto->setTitle(clusXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ClusterSizeXCHistoName, clusXHisto));

  clusXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ClusterSizeXDHistoName.c_str(),clusXNBin,clusXMin,clusXMax);
  clusXHisto->setTitle(clusXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ClusterSizeXDHistoName, clusXHisto));


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ClusterSizeYHistoName);
      if ( histoInfo )
        {
          message<DEBUG> ( log() << (* histoInfo ) );
          clusYNBin = histoInfo->_xBin;
          clusYMin  = histoInfo->_xMin;
          clusYMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) clusYTitle = histoInfo->_title;
        }
    }
  AIDA::IHistogram1D * clusYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ClusterSizeYHistoName.c_str(),clusYNBin,clusYMin,clusYMax);
  clusYHisto->setTitle(clusYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ClusterSizeYHistoName, clusYHisto));

  clusYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ClusterSizeYAHistoName.c_str(),clusYNBin,clusYMin,clusYMax);
  clusYHisto->setTitle(clusYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ClusterSizeYAHistoName, clusYHisto));

  clusYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ClusterSizeYBHistoName.c_str(),clusYNBin,clusYMin,clusYMax);
  clusYHisto->setTitle(clusYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ClusterSizeYBHistoName, clusYHisto));

  clusYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ClusterSizeYCHistoName.c_str(),clusYNBin,clusYMin,clusYMax);
  clusYHisto->setTitle(clusYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ClusterSizeYCHistoName, clusYHisto));

  clusYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ClusterSizeYDHistoName.c_str(),clusYNBin,clusYMin,clusYMax);
  clusYHisto->setTitle(clusYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ClusterSizeYDHistoName, clusYHisto));

  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ClusterSizeXYHistoName);
      if ( histoInfo )
        {
          message<DEBUG> ( log() << (* histoInfo ) );
          clusYNBin = histoInfo->_yBin;
          clusYMin  = histoInfo->_yMin;
          clusYMax  = histoInfo->_yMax;
          clusXNBin = histoInfo->_xBin;
          clusXMin  = histoInfo->_xMin;
          clusXMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) clusXYTitle = histoInfo->_title;
        }
    }

  AIDA::IHistogram2D * clusXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D(_ClusterSizeXYHistoName.c_str(),clusXNBin,clusXMin,clusXMax,clusYNBin,clusYMin,clusYMax);
  clusXYHisto->setTitle(clusXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ClusterSizeXYHistoName, clusXYHisto));

  clusXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D(_ClusterSizeXYAHistoName.c_str(),clusXNBin,clusXMin,clusXMax,clusYNBin,clusYMin,clusYMax);
  clusXYHisto->setTitle(clusXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ClusterSizeXYAHistoName, clusXYHisto));

  clusXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D(_ClusterSizeXYBHistoName.c_str(),clusXNBin,clusXMin,clusXMax,clusYNBin,clusYMin,clusYMax);
  clusXYHisto->setTitle(clusXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ClusterSizeXYBHistoName, clusXYHisto));

  clusXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D(_ClusterSizeXYCHistoName.c_str(),clusXNBin,clusXMin,clusXMax,clusYNBin,clusYMin,clusYMax);
  clusXYHisto->setTitle(clusXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ClusterSizeXYCHistoName, clusXYHisto));

  clusXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D(_ClusterSizeXYDHistoName.c_str(),clusXNBin,clusXMin,clusXMax,clusYNBin,clusYMin,clusYMax);
  clusXYHisto->setTitle(clusXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ClusterSizeXYDHistoName, clusXYHisto));




  // Measured position in X

  int    measXNBin  =  300;
  double measXMin   = -15.;
  double measXMax   =  15.; 
  string measXTitle = "Measured particle position in X";
  string matchXTitle = "Matched hit position in X";
  string unmatchXTitle = "Unmatched hit position in X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_MeasuredXHistoName);
      if ( histoInfo )
        {
          message<DEBUG> ( log() << (* histoInfo ) );
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

  int    measYNBin  = 200;
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
          message<DEBUG> ( log() << (* histoInfo ) );
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

  measXNBin  = 150;
  measXMin   = -15.;
  measXMax   =  15.;
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
          message<DEBUG> ( log() << (* histoInfo ) );
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

  int    fitXNBin  = 150;
  double fitXMin   = -15.;
  double fitXMax   =  15.;
  string fitXTitle = "Fitted particle position in X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_FittedXHistoName);
      if ( histoInfo )
        {
          message<DEBUG> ( log() << (* histoInfo ) );
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

  int    fitYNBin  = 100;
  double fitYMin   = -10.;
  double fitYMax   =  10.;
  string fitYTitle = "Fitted particle position in Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_FittedYHistoName);
      if ( histoInfo )
        {
          message<DEBUG> ( log() << (* histoInfo ) );
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

  fitXNBin  = 150;
  fitXMin   = -15.;
  fitXMax   =  15.;
  fitYNBin  = 100;
  fitYMin   = -10.;
  fitYMax   =  10.;
  string fitXYTitle = "Fitted particle position in XY";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_FittedXYHistoName);
      if ( histoInfo )
        {
          message<DEBUG> ( log() << (* histoInfo ) );
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

  int    shiftXNBin  =  500;
  double shiftXMin   = -0.5;
  double shiftXMax   =  0.5;
  string shiftXTitle = "Measured - fitted X position";

  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ShiftXHistoName);
      if ( histoInfo )
        {
          message<DEBUG> ( log() << (* histoInfo ) );
          shiftXNBin = histoInfo->_xBin;
          shiftXMin  = histoInfo->_xMin;
          shiftXMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) shiftXTitle = histoInfo->_title;
        }
    }

  AIDA::IHistogram1D * shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXHistoName, shiftXHisto));

  // A
  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXAHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXAHistoName, shiftXHisto));
  
  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXA1HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXA1HistoName, shiftXHisto));
  
  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXA2HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXA2HistoName, shiftXHisto));
  
  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXA3HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXA3HistoName, shiftXHisto));
  
  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXA4HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXA4HistoName, shiftXHisto));
  
  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXA5HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXA5HistoName, shiftXHisto));
  
  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXA6HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXA6HistoName, shiftXHisto));

  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXA7HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXA7HistoName, shiftXHisto));
  // B
  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXBHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXBHistoName, shiftXHisto));

  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXB1HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXB1HistoName, shiftXHisto));
 
  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXB2HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXB2HistoName, shiftXHisto));
 
  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXB3HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXB3HistoName, shiftXHisto));
 
  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXB4HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXB4HistoName, shiftXHisto));
 
  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXB5HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXB5HistoName, shiftXHisto));
 
  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXB6HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXB6HistoName, shiftXHisto));
 
  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXB7HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXB7HistoName, shiftXHisto));
  // C
  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXCHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXCHistoName, shiftXHisto));

  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXC1HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXC1HistoName, shiftXHisto));
 
  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXC2HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXC2HistoName, shiftXHisto));
 
  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXC3HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXC3HistoName, shiftXHisto));
 
  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXC4HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXC4HistoName, shiftXHisto));
 
  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXC5HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXC5HistoName, shiftXHisto));
 
  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXC6HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXC6HistoName, shiftXHisto));
 
  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXC7HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXC7HistoName, shiftXHisto));
  // D
  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXDHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXDHistoName, shiftXHisto));

  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXD1HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXD1HistoName, shiftXHisto));

  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXD2HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXD2HistoName, shiftXHisto));

  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXD3HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXD3HistoName, shiftXHisto));

  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXD4HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXD4HistoName, shiftXHisto));

  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXD5HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXD5HistoName, shiftXHisto));

  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXD6HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXD6HistoName, shiftXHisto));

  shiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftXD7HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);
  shiftXHisto->setTitle(shiftXTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXD7HistoName, shiftXHisto));





  // Corresponding background histogram

  shiftXTitle = "Measured - fitted X position (background)";

  AIDA::IHistogram1D * bgshiftXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_BgShiftXHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);

  bgshiftXHisto->setTitle(shiftXTitle.c_str());

  _aidaHistoMap.insert(make_pair(_BgShiftXHistoName, bgshiftXHisto));



  // Measured - fitted position in Y

  int    shiftYNBin  =  500;
  double shiftYMin   = -0.5;
  double shiftYMax   =  0.5;
  string shiftYTitle = "Measured - fitted Y position";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ShiftYHistoName);
      if ( histoInfo )
        {
          message<DEBUG> ( log() << (* histoInfo ) );
          shiftYNBin = histoInfo->_xBin;
          shiftYMin  = histoInfo->_xMin;
          shiftYMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) shiftYTitle = histoInfo->_title;
        }
    }

  AIDA::IHistogram1D * shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYHistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYHistoName, shiftYHisto));


  // A
  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYAHistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYAHistoName, shiftYHisto));
  
  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYA1HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYA1HistoName, shiftYHisto));
  
  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYA2HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYA2HistoName, shiftYHisto));
  
  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYA3HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYA3HistoName, shiftYHisto));
  
  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYA4HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYA4HistoName, shiftYHisto));
  
  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYA5HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYA5HistoName, shiftYHisto));
  
  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYA6HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYA6HistoName, shiftYHisto));

  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYA7HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYA7HistoName, shiftYHisto));
  // B
  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYBHistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYBHistoName, shiftYHisto));

  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYB1HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYB1HistoName, shiftYHisto));
 
  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYB2HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYB2HistoName, shiftYHisto));
 
  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYB3HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYB3HistoName, shiftYHisto));
 
  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYB4HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYB4HistoName, shiftYHisto));
 
  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYB5HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYB5HistoName, shiftYHisto));
 
  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYB6HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYB6HistoName, shiftYHisto));
 
  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYB7HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYB7HistoName, shiftYHisto));
  // C
  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYCHistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYCHistoName, shiftYHisto));

  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYC1HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYC1HistoName, shiftYHisto));
 
  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYC2HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYC2HistoName, shiftYHisto));
 
  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYC3HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYC3HistoName, shiftYHisto));
 
  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYC4HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYC4HistoName, shiftYHisto));
 
  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYC5HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYC5HistoName, shiftYHisto));
 
  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYC6HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYC6HistoName, shiftYHisto));
 
  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYC7HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYC7HistoName, shiftYHisto));
  // D
  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYDHistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYDHistoName, shiftYHisto));

  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYD1HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYD1HistoName, shiftYHisto));

  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYD2HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYD2HistoName, shiftYHisto));

  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYD3HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYD3HistoName, shiftYHisto));

  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYD4HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYD4HistoName, shiftYHisto));

  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYD5HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYD5HistoName, shiftYHisto));

  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYD6HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYD6HistoName, shiftYHisto));

  shiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_ShiftYD7HistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);
  shiftYHisto->setTitle(shiftYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftYD7HistoName, shiftYHisto));



  // Corresponding background histogram

  shiftYTitle = "Measured - fitted Y position (background)";

  AIDA::IHistogram1D * bgshiftYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_BgShiftYHistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);

  bgshiftYHisto->setTitle(shiftYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_BgShiftYHistoName, bgshiftYHisto));




  // Measured - fitted position in X  vs Y

  shiftXNBin  = 150;
  shiftXMin   = -15.;
  shiftXMax   =  15.;
  double shiftVMin   = -1.0;
  double shiftVMax   =  1.0;
  shiftXTitle = "Measured - fitted X position vs Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ShiftXvsYHistoName);
      if ( histoInfo )
        {
          message<DEBUG> ( log() << (* histoInfo ) );
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
          message<DEBUG> ( log() << (* histoInfo ) );
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

  shiftYNBin  = 100;
  shiftYMin   = -10.;
  shiftYMax   =  10.;
  shiftVMin   = -1.0;
  shiftVMax   =  1.0;
  shiftYTitle = "Measured - fitted Y position vs X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ShiftYvsXHistoName);
      if ( histoInfo )
        {
          message<DEBUG> ( log() << (* histoInfo ) );
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
          message<DEBUG> ( log() << (* histoInfo ) );
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

  shiftXNBin  = 150;
  shiftXMin   = -15.;
  shiftXMax   =  15.; 
  int shiftVNBin  = 200;
  shiftVMin   = -1.0;
  shiftVMax   =  1.0;
  shiftXTitle = "Measured - fitted X position vs Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ShiftXvsY2DHistoName);
      if ( histoInfo )
        {
          message<DEBUG> ( log() << (* histoInfo ) );
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
          message<DEBUG> ( log() << (* histoInfo ) );
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

  shiftYNBin  = 150;
  shiftYMin   = -15.;
  shiftYMax   =  15.;
  shiftVNBin  =  200;
  shiftVMin   = -1.0;
  shiftVMax   =  1.0;
  shiftYTitle = "Measured - fitted Y position vs X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ShiftYvsX2DHistoName);
      if ( histoInfo )
        {
          message<DEBUG> ( log() << (* histoInfo ) );
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
          message<DEBUG> ( log() << (* histoInfo ) );
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

  shiftXNBin  =  500;
  shiftXMin   = -0.5;
  shiftXMax   =  0.5;
  shiftYNBin  =  500;
  shiftYMin   = -0.5;
  shiftYMax   =  0.5;
  string shiftXYTitle = "Measured - fitted position in XY";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ShiftXYHistoName);
      if ( histoInfo )
        {
          message<DEBUG> ( log() << (* histoInfo ) );
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

// A
  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYAHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYAHistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYA1HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYA1HistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYA2HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYA2HistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYA3HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYA3HistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYA4HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYA4HistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYA5HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYA5HistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYA6HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYA6HistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYA7HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYA7HistoName, shiftXYHisto));


// B
  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYBHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYBHistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYB1HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYB1HistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYB2HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYB2HistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYB3HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYB3HistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYB4HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYB4HistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYB5HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYB5HistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYB6HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYB6HistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYB7HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYB7HistoName, shiftXYHisto));


// C
  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYCHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYCHistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYC1HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYC1HistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYC2HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYC2HistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYC3HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYC3HistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYC4HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYC4HistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYC5HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYC5HistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYC6HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYC6HistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYC7HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYC7HistoName, shiftXYHisto));


// D
  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYDHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYDHistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYD1HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYD1HistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYD2HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYD2HistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYD3HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYD3HistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYD4HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYD4HistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYD5HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYD5HistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYD6HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYD6HistoName, shiftXYHisto));

  shiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _ShiftXYD7HistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);
  shiftXYHisto->setTitle(shiftXYTitle.c_str());
  _aidaHistoMap.insert(make_pair(_ShiftXYD7HistoName, shiftXYHisto));




  // Corresponding background histogram

  shiftXYTitle = "Measured - fitted position in X-Y (background)";

  AIDA::IHistogram2D * bgshiftXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D(_BgShiftXYHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);

  bgshiftXYHisto->setTitle(shiftXYTitle.c_str());

  _aidaHistoMap.insert(make_pair(_BgShiftXYHistoName, bgshiftXYHisto));



  // Efficiency as a function of the fitted position in X

  int    effiXNBin  = 150;
  double effiXMin   = -15.;
  double effiXMax   =  15.;
  string effiXTitle = "Efficiency vs particle position in X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_EfficiencyXHistoName);
      if ( histoInfo )
        {
          message<DEBUG> ( log() << (* histoInfo ) );
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

  int    effiYNBin  = 150;
  double effiYMin   = -15.;
  double effiYMax   =  15.; 
  string effiYTitle = "Efficiency vs particle position in Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_EfficiencyYHistoName);
      if ( histoInfo )
        {
          message<DEBUG> ( log() << (* histoInfo ) );
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

  effiXNBin  = 30;
  effiXMin   = -15.;
  effiXMax   =  15.;
  effiYNBin  = 20;
  effiYMin   = -10.;
  effiYMax   =  10.;
  string effiXYTitle = "Efficiency vs particle position in XY";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_EfficiencyXYHistoName);
      if ( histoInfo )
        {
          message<DEBUG> ( log() << (* histoInfo ) );
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

  int    noiseXNBin  = 150;
  double noiseXMin   = -15.;
  double noiseXMax   =  15.;
  string noiseXTitle = "Noise fraction vs particle position in X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_NoiseXHistoName);
      if ( histoInfo )
        {
          message<DEBUG> ( log() << (* histoInfo ) );
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

  int    noiseYNBin  = 150;
  double noiseYMin   = -15.;
  double noiseYMax   =  15.;
  string noiseYTitle = "Noise fraction vs particle position in Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_NoiseYHistoName);
      if ( histoInfo )
        {
          message<DEBUG> ( log() << (* histoInfo ) );
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

  noiseXNBin  = 30;
  noiseXMin   = -15.;
  noiseXMax   =  15.;
  noiseYNBin  = 20;
  noiseYMin   = -10.;
  noiseYMax   =  10.;
  string noiseXYTitle = "Noise fraction vs particle position in XY";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_NoiseXYHistoName);
      if ( histoInfo )
        {
          message<DEBUG> ( log() << (* histoInfo ) );
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

  int etaXNBin  =  60;
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
          message<DEBUG> ( log() << (* histoInfo ) );
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

  int etaYNBin  =  60;
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
          message<DEBUG> ( log() << (* histoInfo ) );
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

  etaXNBin  =  60;
  etaXMin   = -0.03;
  etaXMax   = 0.03;
  int etaVNBin  =  60;
  etaVMin   = -0.03;
  etaVMax   = 0.03;
  etaXTitle = "Measured - fitted X position vs local X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_EtaX2DHistoName);
      if ( histoInfo )
        {
          message<DEBUG> ( log() << (* histoInfo ) );
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

  etaYNBin  =  60;
  etaYMin   = -0.03;
  etaYMax   = 0.03;
  etaVNBin  =  60;
  etaVMin   = -0.03;
  etaVMax   = 0.03;
  etaYTitle = "Measured - fitted Y position vs local Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_EtaY2DHistoName);
      if ( histoInfo )
        {
          message<DEBUG> ( log() << (* histoInfo ) );
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

  etaXNBin  =  60;
  etaXMin   = -0.03;
  etaXMax   = 0.03;
  etaYNBin  =  60;
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
          message<DEBUG> ( log() << (* histoInfo ) );
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

  etaXNBin  =  60;
  etaXMin   = -0.03;
  etaXMax   = 0.03;
  etaYNBin  =  60;
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
          message<DEBUG> ( log() << (* histoInfo ) );
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

int EUTelDUTHistograms::guessSensorID(const double * hit ) 
{

  int sensorID = -1;
  double minDistance =  numeric_limits< double >::max() ;
//  double * hitPosition = const_cast<double * > (hit->getPosition());

//  LCCollectionVec * referenceHitVec     = dynamic_cast < LCCollectionVec * > (evt->getCollection( _referenceHitCollectionName));
  if( _referenceHitVec == 0)
  {
    streamlog_out(MESSAGE) << "_referenceHitVec is empty" << endl;
    return 0;
  }

      for(size_t ii = 0 ; ii <  _referenceHitVec->getNumberOfElements(); ii++)
      {
        EUTelReferenceHit* refhit = static_cast< EUTelReferenceHit*> ( _referenceHitVec->getElementAt(ii) ) ;
//        printf(" _referenceHitVec %p refhit %p at %5.2f %5.2f %5.2f wih %5.2f %5.2f %5.2f \n", _referenceHitVec, refhit,
//                refhit->getXOffset(), refhit->getYOffset(), refhit->getZOffset(),
//                refhit->getAlpha(), refhit->getBeta(), refhit->getGamma()  
//               );
        
        TVector3 hit3d( hit[0], hit[1], hit[2] );
        TVector3 hitInPlane( refhit->getXOffset(), refhit->getYOffset(), refhit->getZOffset());
        TVector3 norm2Plane( refhit->getAlpha(), refhit->getBeta(), refhit->getGamma() );
 
        double distance = abs( norm2Plane.Dot(hit3d-hitInPlane) );
//          printf("iPlane %5d   hitPos:  [%8.3f;%8.3f%8.3f]  distance: %8.3f \n", refhit->getSensorID(), hit[0],hit[1],hit[2], distance  );
        if ( distance < minDistance ) 
        {
           minDistance = distance;
           sensorID = refhit->getSensorID();
//           printf("sensorID: %5d \n", sensorID );
        }    

      }

  return sensorID;
}

int EUTelDUTHistograms::getClusterSize(int sensorID, TrackerHitImpl * hit, int& sizeX, int& sizeY, int& subMatrix ) {
  if(hit==0)
{
    streamlog_out( ERROR ) << "An invalid hit pointer supplied! will exit now\n" << endl;
    return -1;
}

        try
        {
            LCObjectVec clusterVector = hit->getRawHits();

            EUTelVirtualCluster * cluster=0;

            if ( hit->getType() == kEUTelBrickedClusterImpl ) {

               // fixed cluster implementation. Remember it
               //  can come from
               //  both RAW and ZS data
   
                cluster = new EUTelBrickedClusterImpl(static_cast<TrackerDataImpl *> ( clusterVector[0] ) );
                
            } else if ( hit->getType() == kEUTelDFFClusterImpl ) {
              
              // fixed cluster implementation. Remember it can come from
              // both RAW and ZS data
              cluster = new EUTelDFFClusterImpl( static_cast<TrackerDataImpl *> ( clusterVector[0] ) );
            } else if ( hit->getType() == kEUTelFFClusterImpl ) {
              
              // fixed cluster implementation. Remember it can come from
              // both RAW and ZS data
              cluster = new EUTelFFClusterImpl( static_cast<TrackerDataImpl *> ( clusterVector[0] ) );
            } 
            else if ( hit->getType() == kEUTelAPIXClusterImpl ) 
            {
              
//              cluster = new EUTelSparseClusterImpl< EUTelAPIXSparsePixel >
//                 ( static_cast<TrackerDataImpl *> ( clusterVector[ 0 ]  ) );

                // streamlog_out(MESSAGE4) << "Type is kEUTelAPIXClusterImpl" << endl;
                TrackerDataImpl * clusterFrame = static_cast<TrackerDataImpl*> ( clusterVector[0] );
                // streamlog_out(MESSAGE4) << "Vector size is: " << clusterVector.size() << endl;

                cluster = new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelAPIXSparsePixel >(clusterFrame);
	      
	        // CellIDDecoder<TrackerDataImpl> cellDecoder(clusterFrame);
                eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelAPIXSparsePixel > *apixCluster = new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelAPIXSparsePixel >(clusterFrame);
                
            }
            else if ( hit->getType() == kEUTelSparseClusterImpl ) 
            {
               cluster = new EUTelSparseClusterImpl< EUTelSimpleSparsePixel > ( static_cast<TrackerDataImpl *> ( clusterVector[0] ) );
            }

            if(cluster != 0)
            {
              float xlocal=-1.;
              float ylocal=-1.;
              cluster->getClusterSize(sizeX,sizeY);
              cluster->getCenterOfGravity(xlocal, ylocal);
              subMatrix = getSubMatrix(sensorID, xlocal);  
            //  printf("cluster x:%d y:%d    \n", sizeX, sizeY );
              return 0;         
            }
          }
          catch(...)
          {
            printf("guess SensorID crashed \n");
          }

return -1;
}

int EUTelDUTHistograms::getSubMatrix(int detectorID, float xlocal)
{
   // quarters : 0-287, 288-575, 576-863, 864-1151
   int fourlocal = static_cast<int>(xlocal*4.); 
   for ( int iLayer = 0; iLayer < _siPlanesLayerLayout->getNLayers(); iLayer++ ) {
      if ( _siPlanesLayerLayout->getID(iLayer) == detectorID ) {
         int subquarter =  fourlocal/static_cast<int>(_siPlanesLayerLayout->getSensitiveNpixelX(iLayer)); 
//         printf("detector %5d xlocal %8.3f %d %d\n", detectorID, xlocal, fourlocal, subquarter) ;        
         return subquarter;       
      }
   }

  return -1; 
}  

#endif

