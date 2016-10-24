// @version: $Id$

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// ROOT includes:
#include "TVector3.h"

// eutelescope includes
#include "EUTelDUTHistograms.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelHistogramManager.h"
#include "EUTelExceptions.h"
#include "EUTelReferenceHit.h"
// for cluster operations:
#include "EUTelSimpleVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelGenericSparseClusterImpl.h"
#include "EUTelGeometricClusterImpl.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/IProfile2D.h>
#include <AIDA/ITree.h>
#endif

// include ROOT
#include "TMath.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

#ifdef USE_GEAR

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

#endif

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>
#include <IMPL/SimTrackerHitImpl.h>
#include <UTIL/CellIDDecoder.h>


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

EUTelDUTHistograms::EUTelDUTHistograms() 
: Processor("EUTelDUTHistograms"),
  _referenceHitCollectionName(""),
  _useReferenceHitCollection(false),
  _referenceHitVec(NULL),
  _histoInfoFileName(""),
  _inputTrackColName(""),
  _inputHitColName(""),
  _inputRecHitColName(""),
  _inputFitHitColName(""),
  _iDUT(0),
  _nRun(0),
  _zDUT(0.0),
  _distMax(0.0),
  _pitchX(0.0),
  _pitchY(0.0),
  _clusterSizeX(),
  _clusterSizeY(),
  _subMatrix(),
  _maptrackid(0),
  _trackhitposX(),
  _trackhitposY(),
  _trackhitsizeX(),
  _trackhitsizeY(),
  _trackhitsubM(),
  _trackhitsensorID(),
  _cluSizeXCut(0),
  _cluSizeYCut(0),
  _trackNCluXCut(0),
  _trackNCluYCut(0),
  _measuredX(),
  _measuredY(),
  _bgmeasuredX(),
  _bgmeasuredY(),
  _localX(),
  _localY(),
  _fittedX(),
  _fittedY(),
  _bgfittedX(),
  _bgfittedY(),
  _DUTalign(),
_ClusterSizeHistos(),
_ShiftHistos(),
_MeasuredHistos(),
_MatchedHistos(),
_UnMatchedHistos(),
_FittedHistos(),
_EfficiencyHistos(),
_BgEfficiencyHistos(),
_NoiseHistos(),
_BgShiftHistos(),
_ShiftXvsYHisto(),
_ShiftYvsXHisto(),
_ShiftXvsXHisto(),
_ShiftYvsYHisto(),
_ShiftXvsY2DHisto(),
_ShiftYvsX2DHisto(),
_ShiftXvsX2DHisto(),
_ShiftYvsY2DHisto(),
_EtaXHisto(),
_EtaYHisto(),
_EtaX2DHisto(),
_EtaY2DHisto(),
_EtaX3DHisto(),
_EtaY3DHisto(),
_PixelEfficiencyHisto    (),
_PixelResolutionXHisto   (),
_PixelResolutionYHisto   (),
_PixelChargeSharingHisto ()

{

  // modify processor description
  _description = "Analysis of DUT performance based on the analytic track fit results";


  // register steering parameters:
  //       name, description, class-variable, default value

  // input collection first:

  registerInputCollection( LCIO::TRACK,
                           "InputTrackCollectionName" ,
                           "Name of the input Track collection"  ,
                           _inputTrackColName,
                           std::string("testfittracks") ) ;

  registerInputCollection( LCIO::TRACKERHIT,
                           "InputHitCollectionName" ,
                           "Name of the input DUT hit collection"  ,
                           _inputHitColName,
                           std::string("hit") ) ;

  registerInputCollection( LCIO::TRACKERHIT,
                           "InputRecHitCollectionName" ,
                           "Name of the input DUT hit collection"  ,
                           _inputRecHitColName,
                           std::string("localHit") ) ;

  registerInputCollection( LCIO::TRACKERHIT,
                           "InputFitHitCollectionName" ,
                           "Name of the input DUT hit collection"  ,
                           _inputFitHitColName,
                           std::string("dummy") ) ;

  // other processor parameters:

  registerProcessorParameter ("ManualDUTid",
                              "Id of telescope layer which should be used as DUT",
                              _iDUT,  static_cast < int > (0));

  registerProcessorParameter ("DistMax",
                              "Maximum allowed distance between fit and matched DUT hit in [mm]",
                              _distMax,  static_cast < double > (0.1));


  registerProcessorParameter ("DUTpitchX",
                              "DUT sensor pitch in X",
                              _pitchX,  static_cast < double > (0.0184));


  registerProcessorParameter ("DUTpitchY",
                              "DUT sensor pitch in Y",
                              _pitchY,  static_cast < double > (0.0184));


  registerOptionalParameter("ReferenceCollection","reference hit collection name ", _referenceHitCollectionName, static_cast <string> ("referenceHit") );
  registerOptionalParameter("UseReferenceCollection","Do you want the reference hit collection to be used for coordinate transformations?",  _useReferenceHitCollection, static_cast< bool   > ( true ));


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

  registerOptionalParameter("cluSizeXCut","cluster size X cut ", _cluSizeXCut, static_cast <int> (-1) );
 
  registerOptionalParameter("cluSizeYCut","cluster size Y cut ", _cluSizeYCut, static_cast <int> (-1) );
  
  registerOptionalParameter("trackNCluXCut","number of hit on a track with _cluSizeX cluster size ", _trackNCluXCut, static_cast <int> (0) );
 
  registerOptionalParameter("trackNCluYCut","number of hit on a track with _cluSizeY cluster size ", _trackNCluYCut, static_cast <int> (0) );

}

void EUTelDUTHistograms::init() {

  // usually a good idea to
  printParameters() ;
  _nRun = 0 ;
  _referenceHitVec = 0;
  _maptrackid = 0;

  _zDUT = geo::gGeometry().siPlaneZPosition(_iDUT);

// Print out geometry information

  message<MESSAGE5> ( log() << "D.U.T. plane  ID = " << _iDUT
                     << "  at Z [mm] = " << _zDUT );

// Book histograms

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  bookHistos();
#endif

}


void EUTelDUTHistograms::processRunHeader( LCRunHeader* runHeader) {

  auto eutelHeader = std::make_unique<EUTelRunHeaderImpl>(runHeader);
  eutelHeader->addProcessor(type());

  _nRun++ ;

  // Decode and print out Run Header information - just a check

  int runNr = runHeader->getRunNumber();

  message<DEBUG5> ( log() << "Processing run header " << _nRun
                     << ", run nr " << runNr );

  const std::string detectorName = runHeader->getDetectorName();
  const std::string detectorDescription = runHeader->getDescription();
  message<DEBUG5> ( log() << detectorName << " : " << detectorDescription ) ;


}

void EUTelDUTHistograms::processEvent( LCEvent * event ) {

  streamlog_out( DEBUG5 ) << "EUTelDUTHistograms::processEvent " << endl;

  {
    if ( _useReferenceHitCollection ) 
    {
       _referenceHitVec = dynamic_cast < LCCollectionVec * > (event->getCollection( _referenceHitCollectionName));
       
       if( streamlog_level( DEBUG5) ){
	 for(size_t ii = 0 ; ii < static_cast<size_t>(_referenceHitVec->getNumberOfElements()); ii++)
	   {
	     streamlog_out( DEBUG5 ) << " check output_refhit at : " << _referenceHitVec << " ";
	     EUTelReferenceHit* output_refhit = static_cast< EUTelReferenceHit*> ( _referenceHitVec->getElementAt(ii) ) ;
	     streamlog_out( DEBUG5 ) << " at : " <<  output_refhit << endl;     
	     streamlog_out( DEBUG5 ) << "CHK sensorID: " <<  output_refhit->getSensorID(   )     
				     << " x    :" <<        output_refhit->getXOffset(    )    
				     << " y    :" <<        output_refhit->getYOffset(    )    
				     << " z    :" <<        output_refhit->getZOffset(    )    
				     << " alfa :" <<        output_refhit->getAlpha()          
				     << " beta :" <<        output_refhit->getBeta()           
				     << " gamma:" <<        output_refhit->getGamma()        << endl ;
	   }
       }
    }
  }
  EUTelEventImpl * euEvent = static_cast<EUTelEventImpl*> ( event );
  if ( euEvent->getEventType() == kEORE ) {
    message<DEBUG5> ( "EORE found: nothing else to do." );
    return;
  }
 

// fill tracking info:
   if( _inputFitHitColName != "dummy" )
   {
     message<DEBUG5> ( log() << "inputFitHitColName = " << _inputFitHitColName << " (not dummy)" << endl); 
     if( read_track_from_collections( event ) > 0 ) 
     {
       return;
     }
   } 
   else
   {
     message<DEBUG5> ( log() << "inputFitHitColName = " << _inputFitHitColName << " (should be called dummy)" << endl); 
     if( read_track( event ) > 0 ) 
     {
       return;
     }
   }
  
  if(streamlog_level(DEBUG5))
  {
    message<DEBUG5> ( log() << _maptrackid  << " fitted tracks " ); 
    for( int itrack=0; itrack<_maptrackid; itrack++)
    {
      message<DEBUG5> ( log() <<" track " << itrack << " has " << _fittedX[itrack].size()  << " fitted positions at DUT " );
    }
  } 

  message<DEBUG5> ( log() << _measuredX.size() << " hits at DUT " );


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  // Histograms of fitted positions
  for( int itrack=0; itrack<_maptrackid; itrack++)
  {
    for( int ifit=0;ifit<static_cast<int>(_fittedX[itrack].size()); ifit++)
    {
      (dynamic_cast<AIDA::IHistogram1D*> ( _FittedHistos.at(projX)))->fill(_fittedX[itrack][ifit]);

      (dynamic_cast<AIDA::IHistogram1D*> ( _FittedHistos.at(projY)))->fill(_fittedY[itrack][ifit]);
      (dynamic_cast<AIDA::IHistogram2D*> ( _FittedHistos.at(projXY)))->fill(_fittedX[itrack][ifit],_fittedY[itrack][ifit]);
      if(streamlog_level(DEBUG5)){
	message<DEBUG5> ( log() << "Fit " << ifit << " [track:"<< itrack << "] "
			  << "   X = " << _fittedX[itrack][ifit]
			  << "   Y = " << _fittedY[itrack][ifit]) ;
      }
    }
  }
#endif

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA) 

  // Histograms of measured positions
  for(int ihit=0;ihit<static_cast<int>(_measuredX.size()); ihit++)
    {
      (dynamic_cast<AIDA::IHistogram1D*> ( _MeasuredHistos.at(projX)))->fill(_measuredX[ihit]);
      (dynamic_cast<AIDA::IHistogram1D*> ( _MeasuredHistos.at(projY)))->fill(_measuredY[ihit]);
      (dynamic_cast<AIDA::IHistogram2D*> ( _MeasuredHistos.at(projXY)))->fill(_measuredX[ihit],_measuredY[ihit]);
      if(streamlog_level(DEBUG5)){
	message<DEBUG5> ( log() << "Hit " << ihit
			  << "   X = " << _measuredX[ihit]
			  << "   Y = " << _measuredY[ihit]) ;
      }
    }
#endif


  // Match measured and fitted positions

  int nMatch=0;
  double distmin;

  for(int itrack=0; itrack< _maptrackid; itrack++)
  {
    int bestfit=-1;
    int besthit=-1;

    distmin = _distMax*_distMax + 10. ;
 
    if( static_cast<int>(_fittedX[itrack].size()) < 1 ) continue;
 
    for(int ifit=0;ifit<static_cast<int>(_fittedX[itrack].size()); ifit++)
    {
      if( _measuredX.empty() ) continue;

      for(int ihit=0; ihit< static_cast<int>(_measuredX.size()) ; ihit++)
        {

          double dist2rd=
            (_measuredX[ihit]-_fittedX[itrack][ifit])*(_measuredX[ihit]-_fittedX[itrack][ifit])
            + (_measuredY[ihit]-_fittedY[itrack][ifit])*(_measuredY[ihit]-_fittedY[itrack][ifit]);

	  if(streamlog_level(DEBUG5)){
	    message<DEBUG5> ( log() << "Fit ["<< itrack << ":" << _maptrackid <<"], ifit= " << ifit << " ["<< _fittedX[itrack][ifit] << ":" << _fittedY[itrack][ifit] << "]" << endl) ;
	    message<DEBUG5> ( log() << "rec " << ihit << " ["<< _measuredX[ihit] << ":" << _measuredY[ihit] << "]" << endl) ;
	    message<DEBUG5> ( log() << "distance : " << TMath::Sqrt( dist2rd )  << endl) ;
	  }
          if(dist2rd<distmin)
            {
              distmin = dist2rd;
              besthit = ihit;
              bestfit = ifit;
            }
        }
 
    }
 
    // Match found:

    if( distmin < _distMax*_distMax  )
      {

        nMatch++;

        // Matched hits positions

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

	// fill once for any matrix ("full detector")
        (dynamic_cast<AIDA::IHistogram1D*> ( _ClusterSizeHistos.at(projX).at(FullDetector)))->fill(_clusterSizeX[besthit]+0.0);
        (dynamic_cast<AIDA::IHistogram1D*> ( _ClusterSizeHistos.at(projY).at(FullDetector)))->fill(_clusterSizeY[besthit]+0.0);
        (dynamic_cast<AIDA::IHistogram2D*> ( _ClusterSizeHistos.at(projXY).at(FullDetector)))->fill(_clusterSizeX[besthit]+0.0,_clusterSizeY[besthit]+0.0);

	// .. and once for the submatrix (identified by the index)
        (dynamic_cast<AIDA::IHistogram1D*> ( _ClusterSizeHistos.at(projX).at(static_cast<detMatrix>(_subMatrix[besthit]))))
	  ->fill(_clusterSizeX[besthit]+0.0);
        (dynamic_cast<AIDA::IHistogram1D*> ( _ClusterSizeHistos.at(projY).at(static_cast<detMatrix>(_subMatrix[besthit]))))
	  ->fill(_clusterSizeY[besthit]+0.0);
        (dynamic_cast<AIDA::IHistogram2D*> ( _ClusterSizeHistos.at(projXY).at(static_cast<detMatrix>(_subMatrix[besthit]))))
	  ->fill(_clusterSizeX[besthit]+0.0,_clusterSizeY[besthit]+0.0);


        (dynamic_cast<AIDA::IHistogram1D*> ( _MatchedHistos.at(projX)))->fill(_measuredX[besthit]);
        (dynamic_cast<AIDA::IHistogram1D*> ( _MatchedHistos.at(projY)))->fill(_measuredY[besthit]);
        (dynamic_cast<AIDA::IHistogram2D*> ( _MatchedHistos.at(projXY)))->fill(_measuredX[besthit],_measuredY[besthit]);

        // Histograms of measured-fitted shifts
        double shiftX =  _measuredX[besthit]-_fittedX[itrack][bestfit];
        double shiftY =  _measuredY[besthit]-_fittedY[itrack][bestfit];

	// fill global: any matrix, any cluster size (cluster size 0 -> any cluster size)
	(dynamic_cast<AIDA::IHistogram1D*> (_ShiftHistos.at(projX).at(FullDetector).at(0)))->fill(shiftX);
	(dynamic_cast<AIDA::IHistogram1D*> (_ShiftHistos.at(projY).at(FullDetector).at(0)))->fill(shiftY);
	(dynamic_cast<AIDA::IHistogram2D*> (_ShiftHistos.at(projXY).at(FullDetector).at(0)))->fill(shiftX, shiftY);
        
	// fill for submatrix and any cluster size
	(dynamic_cast<AIDA::IHistogram1D*> (_ShiftHistos.at(projX).at(static_cast<detMatrix>(_subMatrix[besthit])).at(0)))->fill(shiftX);
	(dynamic_cast<AIDA::IHistogram1D*> (_ShiftHistos.at(projY).at(static_cast<detMatrix>(_subMatrix[besthit])).at(0)))->fill(shiftY);
	(dynamic_cast<AIDA::IHistogram2D*> (_ShiftHistos.at(projXY).at(static_cast<detMatrix>(_subMatrix[besthit])).at(0)))->fill(shiftX, shiftY);
	
	// check that the cluster size is within the limits of our multi diff. binning
	if (_clusterSizeX[besthit] <= HistoMaxClusterSize && _clusterSizeY[besthit] <= HistoMaxClusterSize){
	  // fill for any matrix
	  (dynamic_cast<AIDA::IHistogram1D*> (_ShiftHistos.at(projX).at(FullDetector).at(_clusterSizeX[besthit])))->fill(shiftX);
	  (dynamic_cast<AIDA::IHistogram1D*> (_ShiftHistos.at(projY).at(FullDetector).at(_clusterSizeY[besthit])))->fill(shiftY);
	  // for XY: only if cluster size identical in both x and y
	  if (_clusterSizeX[besthit]==_clusterSizeY[besthit]){
	    (dynamic_cast<AIDA::IHistogram2D*> (_ShiftHistos.at(projXY).at(FullDetector).at(_clusterSizeX[besthit])))->fill(shiftX, shiftY);}

	  // fill for submatrix
	  (dynamic_cast<AIDA::IHistogram1D*> (_ShiftHistos.at(projX).at(static_cast<detMatrix>(_subMatrix[besthit])).at(_clusterSizeX[besthit])))
	    ->fill(shiftX);
	  (dynamic_cast<AIDA::IHistogram1D*> (_ShiftHistos.at(projY).at(static_cast<detMatrix>(_subMatrix[besthit])).at(_clusterSizeY[besthit])))
	    ->fill(shiftY);
	  // for XY: only if cluster size identical in both x and y
	  if (_clusterSizeX[besthit]==_clusterSizeY[besthit]){
	    (dynamic_cast<AIDA::IHistogram2D*> (_ShiftHistos.at(projXY).at(static_cast<detMatrix>(_subMatrix[besthit])).at(_clusterSizeX[besthit])))->fill(shiftX, shiftY);}
	}


       if(  _clusterSizeX[besthit] == 1 &&  _clusterSizeY[besthit] == 1 )
       {
        _PixelEfficiencyHisto->fill( _localX[itrack][bestfit]*1000., _localY[itrack][bestfit]*1000., 1.);
        _PixelResolutionXHisto->fill( _localX[itrack][bestfit]*1000., _localY[itrack][bestfit]*1000., _measuredX[besthit]-_fittedX[itrack][bestfit] );
        _PixelResolutionYHisto->fill( _localX[itrack][bestfit]*1000., _localY[itrack][bestfit]*1000., _measuredY[besthit]-_fittedY[itrack][bestfit] );
       }


       _ShiftXvsYHisto->fill(_fittedY[itrack][bestfit],_measuredX[besthit]-_fittedX[itrack][bestfit]);
       _ShiftYvsXHisto->fill(_fittedX[itrack][bestfit],_measuredY[besthit]-_fittedY[itrack][bestfit]);
       _ShiftXvsX2DHisto->fill(_fittedX[itrack][bestfit], _measuredX[besthit]-_fittedX[itrack][bestfit]);
       _ShiftXvsXHisto->fill(_fittedX[itrack][bestfit], _measuredX[besthit]-_fittedX[itrack][bestfit]);
       
       _ShiftYvsY2DHisto->fill(_fittedY[itrack][bestfit], _measuredY[besthit]-_fittedY[itrack][bestfit]);
       
       _ShiftYvsYHisto->fill(_fittedY[itrack][bestfit], _measuredY[besthit]-_fittedY[itrack][bestfit]);
       
       _ShiftXvsY2DHisto->fill(_fittedY[itrack][bestfit],_measuredX[besthit]-_fittedX[itrack][bestfit]);
       
       _ShiftYvsX2DHisto->fill(_fittedX[itrack][bestfit],_measuredY[besthit]-_fittedY[itrack][bestfit]);


        // Eta function check plots
       if(  _clusterSizeX[besthit] == 1 &&  _clusterSizeY[besthit] == 1 ){
        _EtaXHisto->fill(_localX[itrack][bestfit],_measuredX[besthit]-_fittedX[itrack][bestfit]);
        _EtaYHisto->fill(_localY[itrack][bestfit],_measuredY[besthit]-_fittedY[itrack][bestfit]);
        _EtaX2DHisto->fill(_localX[itrack][bestfit],_measuredX[besthit]-_fittedX[itrack][bestfit]);
        _EtaY2DHisto->fill(_localY[itrack][bestfit],_measuredY[besthit]-_fittedY[itrack][bestfit]);
        _EtaX3DHisto->fill(_localX[itrack][bestfit],_localY[itrack][bestfit],_measuredX[besthit]-_fittedX[itrack][bestfit]);
        _EtaY3DHisto->fill(_localX[itrack][bestfit],_localY[itrack][bestfit],_measuredY[besthit]-_fittedY[itrack][bestfit]);
       } 
        // extend Eta histograms to 2 pitch range

#endif

        if(_localX[itrack][bestfit]<0)
          _localX[itrack][bestfit]+=_pitchX;
        else
          _localX[itrack][bestfit]-=_pitchX;

        if(_localY[itrack][bestfit]<0)
          _localY[itrack][bestfit]+=_pitchY;
        else
          _localY[itrack][bestfit]-=_pitchY;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

        _EtaXHisto->fill(_localX[itrack][bestfit],_measuredX[besthit]-_fittedX[itrack][bestfit]);
        _EtaYHisto->fill(_localY[itrack][bestfit],_measuredY[besthit]-_fittedY[itrack][bestfit]);
        _EtaX2DHisto->fill(_localX[itrack][bestfit],_measuredX[besthit]-_fittedX[itrack][bestfit]);
        _EtaY2DHisto->fill(_localY[itrack][bestfit],_measuredY[besthit]-_fittedY[itrack][bestfit]);

        // Efficiency plots
        (dynamic_cast<AIDA::IProfile1D*> ( _EfficiencyHistos.at(projX)))->fill(_fittedX[itrack][bestfit],1.);
        (dynamic_cast<AIDA::IProfile1D*> ( _EfficiencyHistos.at(projY)))->fill(_fittedY[itrack][bestfit],1.);
        (dynamic_cast<AIDA::IProfile2D*> ( _EfficiencyHistos.at(projXY)))->fill(_fittedX[itrack][bestfit],_fittedY[itrack][bestfit],1.);


        // Noise plots
        (dynamic_cast<AIDA::IProfile1D*> ( _NoiseHistos.at(projX)))->fill(_measuredX[besthit],0.);
        (dynamic_cast<AIDA::IProfile1D*> ( _NoiseHistos.at(projY)))->fill(_measuredY[besthit],0.);
        (dynamic_cast<AIDA::IProfile2D*> ( _NoiseHistos.at(projXY)))->fill(_measuredX[besthit],_measuredY[besthit],0.);

#endif

        // Remove matched entries from the list (so the next matching pair
        // can be looked for)

        _fittedX[itrack].erase(_fittedX[itrack].begin()+bestfit);
        _fittedY[itrack].erase(_fittedY[itrack].begin()+bestfit);

        _measuredX.erase(_measuredX.begin()+besthit);
        _measuredY.erase(_measuredY.begin()+besthit);

        _localX[itrack].erase(_localX[itrack].begin()+bestfit);
        _localY[itrack].erase(_localY[itrack].begin()+bestfit);


     }

    // End of loop of matching DUT hits to fitted positions

    if(streamlog_level(DEBUG5)){
      message<DEBUG5> ( log() << nMatch << " DUT hits matched to fitted tracks ");
      message<DEBUG5> ( log() << _measuredX.size() << " DUT hits not matched to any track ");
      message<DEBUG5> ( log() << "track "<<itrack<<" has " << _fittedX[itrack].size() << " _fittedX[itrack].size() not matched to any DUT hit ");
    }

  // Efficiency plots - unmatched tracks

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  for(int ifit=0;ifit<static_cast<int>(_localX[itrack].size()); ifit++)
    {
      _PixelEfficiencyHisto->fill(_localX[itrack][ ifit ]*1000.,_localY[itrack][ ifit ]*1000.,0.);
    }

  for(int ifit=0;ifit<static_cast<int>(_fittedX[itrack].size()); ifit++)
    {
      (dynamic_cast<AIDA::IProfile1D*> ( _EfficiencyHistos.at(projX)))->fill(_fittedX[itrack][ifit],0.);
      (dynamic_cast<AIDA::IProfile1D*> ( _EfficiencyHistos.at(projY)))->fill(_fittedY[itrack][ifit],0.);
      (dynamic_cast<AIDA::IProfile2D*> ( _EfficiencyHistos.at(projXY)))->fill(_fittedX[itrack][ifit],_fittedY[itrack][ifit],0.);
    }
  #endif
}


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  // Noise plots - unmatched hits

  for(int ihit=0;ihit<static_cast<int>(_measuredX.size()); ihit++){
      (dynamic_cast<AIDA::IProfile1D*> ( _NoiseHistos.at(projX)))->fill(_measuredX[ihit],1.);
      (dynamic_cast<AIDA::IProfile1D*> ( _NoiseHistos.at(projY)))->fill(_measuredY[ihit],1.);
      (dynamic_cast<AIDA::IProfile2D*> ( _NoiseHistos.at(projXY)))->fill(_measuredX[ihit],_measuredY[ihit],1.);

      // Unmatched hit positions
      (dynamic_cast<AIDA::IHistogram1D*> ( _UnMatchedHistos.at(projX)))->fill(_measuredX[ihit]);
      (dynamic_cast<AIDA::IHistogram1D*> ( _UnMatchedHistos.at(projY)))->fill(_measuredY[ihit]);
      (dynamic_cast<AIDA::IHistogram2D*> ( _UnMatchedHistos.at(projXY)))->fill(_measuredX[ihit],_measuredY[ihit]);

    }

#endif

  if ( isFirstEvent() ) _isFirstEvent = false;

  return;
}



void EUTelDUTHistograms::check( LCEvent * /* evt */ ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void EUTelDUTHistograms::end(){

	// fill global: any matrix, any cluster size (cluster size 0 -> any cluster size)
	streamlog_out( MESSAGE4 ) << "DUT " << 
        (dynamic_cast<AIDA::IHistogram1D*> (_ShiftHistos.at(projX).at(FullDetector).at(0)))->allEntries() << " " <<
	(dynamic_cast<AIDA::IHistogram1D*> (_ShiftHistos.at(projX).at(FullDetector).at(0)))->mean()*1000. << " " <<
	(dynamic_cast<AIDA::IHistogram1D*> (_ShiftHistos.at(projX).at(FullDetector).at(0)))->rms()*1000.  << " " <<
        (dynamic_cast<AIDA::IHistogram1D*> (_ShiftHistos.at(projY).at(FullDetector).at(0)))->allEntries() << " " <<
	(dynamic_cast<AIDA::IHistogram1D*> (_ShiftHistos.at(projY).at(FullDetector).at(0)))->mean()*1000. << " " <<
	(dynamic_cast<AIDA::IHistogram1D*> (_ShiftHistos.at(projY).at(FullDetector).at(0)))->rms()*1000.  << " " << endl;
      
}



void EUTelDUTHistograms::bookHistos()
{

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  std::string _ClusterSizeHistoBaseName  = "clusterSize"; // append [X, Y, XY][submatrix][A-D]
  std::string _MeasuredHistoBaseName  = "measured"; // append [X, Y, XY]
  std::string _MatchedHistoBaseName  = "matched"; // append [X, Y, XY]
  std::string _UnMatchedHistoBaseName  = "unmatched"; // append [X, Y, XY]
  std::string _FittedHistoBaseName  = "fitted"; // append [X, Y, XY]
  std::string _EfficiencyHistoBaseName  = "DUTeffi"; // append [X, Y, XY]
  std::string _BgEfficiencyBaseHistoName  = "BGeffi"; // append [X, Y, XY]
  std::string _NoiseHistoBaseName  = "DUTnoise"; // append [X, Y, XY]
  std::string _ShiftHistoBaseName       = "DUTshift"; // append [X, Y, XY][A-D][1-HistoMaxClusterSize]
  std::string _BgShiftHistoBaseName       = "BGshift"; // append [X, Y, XY]
  
  std::string _ShiftXvsYHistoName      = "DUTshiftXvsY";
  std::string _ShiftYvsXHistoName      = "DUTshiftYvsX";
  std::string _ShiftXvsY2DHistoName    = "DUTshiftXvsY2D";
  std::string _ShiftYvsX2DHistoName    = "DUTshiftYvsX2D";

  std::string _ShiftYvsYHistoName      = "DUTshiftYvsY";
  std::string _ShiftXvsXHistoName      = "DUTshiftXvsX";
  std::string _ShiftXvsX2DHistoName    = "DUTshiftXvsX2D";
  std::string _ShiftYvsY2DHistoName    = "DUTshiftYvsY2D";

  std::string _EtaXHistoName           = "EtaX";
  std::string _EtaYHistoName           = "EtaY";
  std::string _EtaX2DHistoName         = "EtaX2D";
  std::string _EtaY2DHistoName         = "EtaY2D";
  std::string _EtaX3DHistoName         = "EtaX3D";
  std::string _EtaY3DHistoName         = "EtaY3D";

  std::string _PixelEfficiencyHistoName       = "PixelEfficiency";
  std::string _PixelResolutionXHistoName      = "PixelResolutionX";
  std::string _PixelResolutionYHistoName      = "PixelResolutionY";
  std::string _PixelChargeSharingHistoName    = "PixelChargeSharing";


  message<MESSAGE5> ( log() << "Booking histograms " );


  message<MESSAGE5> ( log() << "Histogram information searched in " << _histoInfoFileName);

  auto histoMgr = std::make_unique<EUTelHistogramManager>( _histoInfoFileName);
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

  // helper map to construct histo names and titles
  std::map<projAxis,std::string>projectionName;
  projectionName.insert(make_pair(projX,"X"));
  projectionName.insert(make_pair(projY,"Y"));
  projectionName.insert(make_pair(projXY,"XY"));
  
  for (int thisProjection = 0; thisProjection<=2; thisProjection++){

    // helper maps: used to construct histo titles and histo names
    std::map< detMatrix, std::string> detMatrixDescription; // maps enumerator value to description e.g. for histo title
    detMatrixDescription.insert(make_pair(FullDetector,"(all matrices)"));
    detMatrixDescription.insert(make_pair(SubMatrixA,"(sub matrix A)"));
    detMatrixDescription.insert(make_pair(SubMatrixB,"(sub matrix B)"));
    detMatrixDescription.insert(make_pair(SubMatrixC,"(sub matrix C)"));
    detMatrixDescription.insert(make_pair(SubMatrixD,"(sub matrix D)"));
    
    std::map< detMatrix, std::string> detMatrixName; // maps enumerator value to name element e.g. for histo names
    detMatrixName.insert(make_pair(FullDetector,""));
    detMatrixName.insert(make_pair(SubMatrixA,"A"));
    detMatrixName.insert(make_pair(SubMatrixB,"B"));
    detMatrixName.insert(make_pair(SubMatrixC,"C"));
    detMatrixName.insert(make_pair(SubMatrixD,"D"));

    // set up histos for sub matrix studies
    for (int thisMatrix = 0; thisMatrix<5; thisMatrix++){
    
      {
	// cluster size
	int    NBinX  =  12 ; // defaults; 
	double MinX   =  -1.5;//  defaults; 
	double MaxX   =  10.5;//  defaults; 
	// for XY only:
	int    NBinY  =  12 ; // defaults; 
	double MinY   =  -1.5;//  defaults; 
	double MaxY   =  10.5;//  defaults; 
	
	// construct the name of the histogram we are looking for:
	std::string thisHistoName = "";
	if (static_cast<detMatrix>(thisMatrix) != FullDetector)
	  thisHistoName = _ClusterSizeHistoBaseName + projectionName.at(static_cast<projAxis>(thisProjection)) 
	    + "submatrix"+detMatrixName.at(static_cast<detMatrix>(thisMatrix));
	else
	  thisHistoName = _ClusterSizeHistoBaseName + projectionName.at(static_cast<projAxis>(thisProjection));
	
	// set default title
	std::string Title = "Fitted cluster size in " + projectionName.at(static_cast<projAxis>(thisProjection)) + " " + detMatrixDescription.at(static_cast<detMatrix>(thisMatrix));

	if (thisProjection!=projXY)
	  Title += ";cluster size [pixel];# of clusters";
	else
	  Title +=";cluster size [pixel];cluster size [pixel]";
	
	// replace info with user settings from histomanager
	if ( isHistoManagerAvailable ){
	  histoInfo = histoMgr->getHistogramInfo(thisHistoName);
	  if ( histoInfo ) {
	    message<DEBUG5> ( log() << (* histoInfo ) );
	    NBinX = histoInfo->_xBin;
	    MinX  = histoInfo->_xMin;
	    MaxX  = histoInfo->_xMax;
	    NBinY = histoInfo->_yBin;
	    MinY  = histoInfo->_yMin;
	    MaxY  = histoInfo->_yMax;
	    if ( histoInfo->_title != "" ) Title = histoInfo->_title;
	  }
	}
	
	// create 'outer' map if not existing
	_ClusterSizeHistos.insert( make_pair( static_cast<projAxis>(thisProjection), std::map<detMatrix,AIDA::IBaseHistogram*>()));
	AIDA::IBaseHistogram * thisHisto;
	if (static_cast<projAxis>(thisProjection)== projXY) 
	  thisHisto = dynamic_cast<AIDA::IBaseHistogram*>(AIDAProcessor::histogramFactory(this)->createHistogram2D(thisHistoName.c_str(),NBinX,MinX,MaxX,NBinY,MinY,MaxY));
	else 
	  thisHisto = dynamic_cast<AIDA::IBaseHistogram*>(AIDAProcessor::histogramFactory(this)->createHistogram1D(thisHistoName.c_str(),NBinX,MinX,MaxX));
	
	thisHisto->setTitle(Title);
	_ClusterSizeHistos.at(static_cast<projAxis>(thisProjection)).insert(make_pair(static_cast<detMatrix>(thisMatrix), thisHisto) );
      } // cluster size histos
  
      { 
	// Measured - fitted position
	
	// studied as function of cluster size
	for (int thisClusterSize=0;thisClusterSize<=HistoMaxClusterSize;thisClusterSize++){

	  int    NBinX  =  500;
	  double MinX   = -0.5;
	  double MaxX   =  0.5;
	  int    NBinY  =  500;
	  double MinY   = -0.5;
	  double MaxY   =  0.5;

	  // construct the name of the histogram we are looking for:
	  // append [X, Y, XY][A-D][1-HistoMaxClusterSize]
	  std::string thisHistoName = "";
	  if (thisClusterSize==0){
	      thisHistoName = _ShiftHistoBaseName + projectionName.at(static_cast<projAxis>(thisProjection)) 
		+detMatrixName.at(static_cast<detMatrix>(thisMatrix));
	    } else{
	      thisHistoName = _ShiftHistoBaseName + projectionName.at(static_cast<projAxis>(thisProjection)) 
		+detMatrixName.at(static_cast<detMatrix>(thisMatrix))+to_string(thisClusterSize);
	    }
	
	  // set default title
	  std::string Title = "";
	  if (thisClusterSize==0){
	     Title = " Measured - fitted " + projectionName.at(static_cast<projAxis>(thisProjection)) + " position " + detMatrixDescription.at(static_cast<detMatrix>(thisMatrix));
	  } else{
	    Title = " Measured - fitted " + projectionName.at(static_cast<projAxis>(thisProjection)) + " position " + detMatrixDescription.at(static_cast<detMatrix>(thisMatrix)) + "cl.sz. " + to_string(thisClusterSize);
	  }

	  if (thisProjection==projXY)
	    Title += "; #Delta X [mm]; #Delta Y [mm];";
	  else
	    Title +="; #Delta"+projectionName.at(static_cast<projAxis>(thisProjection))+" [mm];# of cluster";

	  // replace info with user settings from histomanager
	  if ( isHistoManagerAvailable ){
	    histoInfo = histoMgr->getHistogramInfo(thisHistoName);
	    if ( histoInfo ) {
	      message<DEBUG5> ( log() << (* histoInfo ) );
	      NBinX = histoInfo->_xBin;
	      MinX  = histoInfo->_xMin;
	      MaxX  = histoInfo->_xMax;
	      NBinY = histoInfo->_yBin;
	      MinY  = histoInfo->_yMin;
	      MaxY  = histoInfo->_yMax;
	      if ( histoInfo->_title != "" ) Title = histoInfo->_title;
	    }
	  }

	  // create 'outer' maps if not existing
	  _ShiftHistos.insert( make_pair( static_cast<projAxis>(thisProjection), std::map< detMatrix, std::map< int, AIDA::IBaseHistogram*> >()));
	  _ShiftHistos.at( static_cast<projAxis>(thisProjection)).insert(make_pair (static_cast<detMatrix>(thisMatrix) , std::map< int, AIDA::IBaseHistogram*>()));

	  AIDA::IBaseHistogram * thisHisto;
	  if (static_cast<projAxis>(thisProjection)== projXY) 
	    thisHisto = dynamic_cast<AIDA::IBaseHistogram*>(AIDAProcessor::histogramFactory(this)->createHistogram2D(thisHistoName.c_str(),NBinX,MinX,MaxX,NBinY,MinY,MaxY));
	  else 
	    thisHisto = dynamic_cast<AIDA::IBaseHistogram*>(AIDAProcessor::histogramFactory(this)->createHistogram1D(thisHistoName.c_str(),NBinX,MinX,MaxX));
	  
	  thisHisto->setTitle(Title);
	  _ShiftHistos.at(static_cast<projAxis>(thisProjection)).at(static_cast<detMatrix>(thisMatrix)).insert(make_pair(thisClusterSize, thisHisto) );
	} // cluster size
      } // shift histos

    } // submatrix loop

    {
      // Measured - fitted position background histos

      int    NBinX  =  500;
      double MinX   = -0.5;
      double MaxX   =  0.5;
      int    NBinY  =  500;
      double MinY   = -0.5;
      double MaxY   =  0.5;

      // construct the name of the histogram we are looking for:
      // append [X, Y, XY]
      std::string thisHistoName = "";
      thisHistoName = _BgShiftHistoBaseName + projectionName.at(static_cast<projAxis>(thisProjection));
	
      // set default title
      std::string Title = " Measured - fitted " + projectionName.at(static_cast<projAxis>(thisProjection)) + " position (background) ";

	if (thisProjection!=projXY)
	  Title += "; #Delta"+projectionName.at(static_cast<projAxis>(thisProjection))+" [mm];# of hits";
	else
	  Title +="#Delta X [mm]; #Delta Y [mm]";

	  
      AIDA::IBaseHistogram * thisHisto;
      if (static_cast<projAxis>(thisProjection)== projXY) 
	thisHisto = dynamic_cast<AIDA::IBaseHistogram*>(AIDAProcessor::histogramFactory(this)->createHistogram2D(thisHistoName.c_str(),NBinX,MinX,MaxX,NBinY,MinY,MaxY));
      else 
	thisHisto = dynamic_cast<AIDA::IBaseHistogram*>(AIDAProcessor::histogramFactory(this)->createHistogram1D(thisHistoName.c_str(),NBinX,MinX,MaxX));
	  
      thisHisto->setTitle(Title);
      _BgShiftHistos.insert(make_pair(static_cast<projAxis>(thisProjection), thisHisto) );

    } // shift histos: Corresponding background histogram
  

    { // Measured position
      int    NBinX  =  300 ; // defaults; 
      double MinX   =  -15.;//  defaults; 
      double MaxX   =  15.;//  defaults; 
      // for XY only:
      int    NBinY  =  200 ; // defaults; 
      double MinY   =  -10.;//  defaults; 
      double MaxY   =  10.;//  defaults; 
	
      if (static_cast<projAxis>(thisProjection) == projY){
	NBinX  =  200 ;
	MinX   =  -10.;
	MaxX   =  10.;
      }

      // construct the name of the histogram we are looking for:
      std::string thisHistoName = "";
      thisHistoName = _MeasuredHistoBaseName + projectionName.at(static_cast<projAxis>(thisProjection));
	
      // set default title
      std::string Title = "Measured particle position in " + projectionName.at(static_cast<projAxis>(thisProjection))
	+"; "+projectionName.at(static_cast<projAxis>(thisProjection))+" [mm];# of hits";

	if (thisProjection!=projXY)
	  Title += ";cluster size [pixel];# of clusters";
	else
	  Title +=";cluster size [pixel];cluster size [pixel]";
	
      // replace info with user settings from histomanager
      if ( isHistoManagerAvailable ){
	histoInfo = histoMgr->getHistogramInfo(thisHistoName);
	if ( histoInfo ) {
	  message<DEBUG5> ( log() << (* histoInfo ) );
	  NBinX = histoInfo->_xBin;
	  MinX  = histoInfo->_xMin;
	  MaxX  = histoInfo->_xMax;
	  NBinY = histoInfo->_yBin;
	  MinY  = histoInfo->_yMin;
	  MaxY  = histoInfo->_yMax;
	  if ( histoInfo->_title != "" ) Title = histoInfo->_title;
	}
      }
	
      AIDA::IBaseHistogram * thisHisto;
      if (static_cast<projAxis>(thisProjection)== projXY) 
	thisHisto = dynamic_cast<AIDA::IBaseHistogram*>(AIDAProcessor::histogramFactory(this)->createHistogram2D(thisHistoName.c_str(),NBinX,MinX,MaxX,NBinY,MinY,MaxY));
      else 
	thisHisto = dynamic_cast<AIDA::IBaseHistogram*>(AIDAProcessor::histogramFactory(this)->createHistogram1D(thisHistoName.c_str(),NBinX,MinX,MaxX));
	
      thisHisto->setTitle(Title);
      _MeasuredHistos.insert(make_pair(static_cast<projAxis>(thisProjection), thisHisto) );

    } // measured position



    { // matched position
      int    NBinX  =  300 ; // defaults; 
      double MinX   =  -15.;//  defaults; 
      double MaxX   =  15.;//  defaults; 
      // for XY only:
      int    NBinY  =  200 ; // defaults; 
      double MinY   =  -10.;//  defaults; 
      double MaxY   =  10.;//  defaults; 
	
      if (static_cast<projAxis>(thisProjection)== projY){
	NBinX  =  200 ;
	MinX   =  -10.;
	MaxX   =  10.;
      }

      // construct the name of the histogram we are looking for:
      std::string thisHistoName = "";
      thisHistoName = _MatchedHistoBaseName + projectionName.at(static_cast<projAxis>(thisProjection));
	
      // set default title
      std::string Title = "Matched particle position in " + projectionName.at(static_cast<projAxis>(thisProjection));

      if (thisProjection==projXY)
	Title +="X [mm]; Y [mm]";
      else
	Title +="; "+projectionName.at(static_cast<projAxis>(thisProjection))+" [mm];# of hits";

      // replace info with user settings from histomanager
      if ( isHistoManagerAvailable ){
	histoInfo = histoMgr->getHistogramInfo(thisHistoName);
	if ( histoInfo ) {
	  message<DEBUG5> ( log() << (* histoInfo ) );
	  NBinX = histoInfo->_xBin;
	  MinX  = histoInfo->_xMin;
	  MaxX  = histoInfo->_xMax;
	  NBinY = histoInfo->_yBin;
	  MinY  = histoInfo->_yMin;
	  MaxY  = histoInfo->_yMax;
	  if ( histoInfo->_title != "" ) Title = histoInfo->_title;
	}
      }
	
      AIDA::IBaseHistogram * thisHisto;
      if (static_cast<projAxis>(thisProjection)== projXY) 
	thisHisto = dynamic_cast<AIDA::IBaseHistogram*>(AIDAProcessor::histogramFactory(this)->createHistogram2D(thisHistoName.c_str(),NBinX,MinX,MaxX,NBinY,MinY,MaxY));
      else 
	thisHisto = dynamic_cast<AIDA::IBaseHistogram*>(AIDAProcessor::histogramFactory(this)->createHistogram1D(thisHistoName.c_str(),NBinX,MinX,MaxX));
	
      thisHisto->setTitle(Title);
      _MatchedHistos.insert(make_pair(static_cast<projAxis>(thisProjection), thisHisto) );

    } // matched position

    { // unmatched position
      int    NBinX  =  300 ; // defaults; 
      double MinX   =  -15.;//  defaults; 
      double MaxX   =  15.;//  defaults; 
      // for XY only:
      int    NBinY  =  200 ; // defaults; 
      double MinY   =  -10.;//  defaults; 
      double MaxY   =  10.;//  defaults; 
	
      if (static_cast<projAxis>(thisProjection)== projY){
	NBinX  =  200 ;
	MinX   =  -10.;
	MaxX   =  10.;
      }

      // construct the name of the histogram we are looking for:
      std::string thisHistoName = "";
      thisHistoName = _UnMatchedHistoBaseName + projectionName.at(static_cast<projAxis>(thisProjection));
	
      // set default title
      std::string Title = "Unmatched particle position in " + projectionName.at(static_cast<projAxis>(thisProjection));

      if (thisProjection==projXY)
	Title +="X [mm]; Y [mm]";
      else
	Title +="; "+projectionName.at(static_cast<projAxis>(thisProjection))+" [mm];# of hits";
	
      // replace info with user settings from histomanager
      if ( isHistoManagerAvailable ){
	histoInfo = histoMgr->getHistogramInfo(thisHistoName);
	if ( histoInfo ) {
	  message<DEBUG5> ( log() << (* histoInfo ) );
	  NBinX = histoInfo->_xBin;
	  MinX  = histoInfo->_xMin;
	  MaxX  = histoInfo->_xMax;
	  NBinY = histoInfo->_yBin;
	  MinY  = histoInfo->_yMin;
	  MaxY  = histoInfo->_yMax;
	  if ( histoInfo->_title != "" ) Title = histoInfo->_title;
	}
      }
	
      AIDA::IBaseHistogram * thisHisto;
      if (static_cast<projAxis>(thisProjection)== projXY) 
	thisHisto = dynamic_cast<AIDA::IBaseHistogram*>(AIDAProcessor::histogramFactory(this)->createHistogram2D(thisHistoName.c_str(),NBinX,MinX,MaxX,NBinY,MinY,MaxY));
      else 
	thisHisto = dynamic_cast<AIDA::IBaseHistogram*>(AIDAProcessor::histogramFactory(this)->createHistogram1D(thisHistoName.c_str(),NBinX,MinX,MaxX));
	
      thisHisto->setTitle(Title);
      _UnMatchedHistos.insert(make_pair(static_cast<projAxis>(thisProjection), thisHisto) );

    } // unmatched position

    { // fitted position
      int    NBinX  =  150 ; // defaults; 
      double MinX   =  -15.;//  defaults; 
      double MaxX   =  15.;//  defaults; 
      // for XY only:
      int    NBinY  =  100 ; // defaults; 
      double MinY   =  -10.;//  defaults; 
      double MaxY   =  10.;//  defaults; 
	
      if (static_cast<projAxis>(thisProjection)== projY){
	NBinX  =  200 ;
	MinX   =  -10.;
	MaxX   =  10.;
      }

      // construct the name of the histogram we are looking for:
      std::string thisHistoName = "";
      thisHistoName = _FittedHistoBaseName + projectionName.at(static_cast<projAxis>(thisProjection));
	
      // set default title
      std::string Title = "Fitted particle position in " + projectionName.at(static_cast<projAxis>(thisProjection));

      if (thisProjection==projXY)
	Title +="X [mm]; Y [mm]";
      else
	Title +="; "+projectionName.at(static_cast<projAxis>(thisProjection))+" [mm];# of hits";
	
      // replace info with user settings from histomanager
      if ( isHistoManagerAvailable ){
	histoInfo = histoMgr->getHistogramInfo(thisHistoName);
	if ( histoInfo ) {
	  message<DEBUG5> ( log() << (* histoInfo ) );
	  NBinX = histoInfo->_xBin;
	  MinX  = histoInfo->_xMin;
	  MaxX  = histoInfo->_xMax;
	  NBinY = histoInfo->_yBin;
	  MinY  = histoInfo->_yMin;
	  MaxY  = histoInfo->_yMax;
	  if ( histoInfo->_title != "" ) Title = histoInfo->_title;
	}
      }
	
      AIDA::IBaseHistogram * thisHisto;
      if (static_cast<projAxis>(thisProjection)== projXY) 
	thisHisto = dynamic_cast<AIDA::IBaseHistogram*>(AIDAProcessor::histogramFactory(this)->createHistogram2D(thisHistoName.c_str(),NBinX,MinX,MaxX,NBinY,MinY,MaxY));
      else 
	thisHisto = dynamic_cast<AIDA::IBaseHistogram*>(AIDAProcessor::histogramFactory(this)->createHistogram1D(thisHistoName.c_str(),NBinX,MinX,MaxX));
	
      thisHisto->setTitle(Title);
      _FittedHistos.insert(make_pair(static_cast<projAxis>(thisProjection), thisHisto) );

    } // fitted position



    { // efficiency as fcn of fitted position
      int    NBinX  =  150 ; // defaults; 
      double MinX   =  -15.;//  defaults; 
      double MaxX   =  15.;//  defaults; 
      // for XY only:
      int    NBinY  =  100 ; // defaults; 
      double MinY   =  -10.;//  defaults; 
      double MaxY   =  10.;//  defaults; 
	
      if (static_cast<projAxis>(thisProjection)== projY){
	NBinX  =  100 ;
	MinX   =  -10.;
	MaxX   =  10.;
      }

      // construct the name of the histogram we are looking for:
      std::string thisHistoName = "";
      thisHistoName = _EfficiencyHistoBaseName + projectionName.at(static_cast<projAxis>(thisProjection));
	
       // set default title
      std::string Title = "Efficiency vs particle position in " + projectionName.at(static_cast<projAxis>(thisProjection));

      if (thisProjection==projXY)
	Title +="X [mm]; Y [mm]";
      else
	Title +="; "+projectionName.at(static_cast<projAxis>(thisProjection))+" [mm];efficiency";
	
      // replace info with user settings from histomanager
      if ( isHistoManagerAvailable ){
	histoInfo = histoMgr->getHistogramInfo(thisHistoName);
	if ( histoInfo ) {
	  message<DEBUG5> ( log() << (* histoInfo ) );
	  NBinX = histoInfo->_xBin;
	  MinX  = histoInfo->_xMin;
	  MaxX  = histoInfo->_xMax;
	  NBinY = histoInfo->_yBin;
	  MinY  = histoInfo->_yMin;
	  MaxY  = histoInfo->_yMax;
	  if ( histoInfo->_title != "" ) Title = histoInfo->_title;
	}
      }
	
      AIDA::IBaseHistogram * thisHisto;
      if (static_cast<projAxis>(thisProjection)== projXY) 
	thisHisto = dynamic_cast<AIDA::IBaseHistogram*>(AIDAProcessor::histogramFactory(this)->createProfile2D(thisHistoName.c_str(),NBinX,MinX,MaxX,NBinY,MinY,MaxY));
      else 
	thisHisto = dynamic_cast<AIDA::IBaseHistogram*>(AIDAProcessor::histogramFactory(this)->createProfile1D(thisHistoName.c_str(),NBinX,MinX,MaxX));
	
      thisHisto->setTitle(Title);
      _EfficiencyHistos.insert(make_pair(static_cast<projAxis>(thisProjection), thisHisto) );

    } // efficiency as fcn of fitted position


    { // background match prob. as fcn of fitted position
      int    NBinX  =  150 ; // defaults; 
      double MinX   =  -15.;//  defaults; 
      double MaxX   =  15.;//  defaults; 
      // for XY only:
      int    NBinY  =  100 ; // defaults; 
      double MinY   =  -10.;//  defaults; 
      double MaxY   =  10.;//  defaults; 
	
      if (static_cast<projAxis>(thisProjection)== projY){
	NBinX  =  100 ;
	MinX   =  -10.;
	MaxX   =  10.;
      }

      // construct the name of the histogram we are looking for:
      std::string thisHistoName = "";
      thisHistoName = _BgEfficiencyBaseHistoName + projectionName.at(static_cast<projAxis>(thisProjection));
	
      // set default title
      std::string Title = "Bckgrd match prob vs part pos in " + projectionName.at(static_cast<projAxis>(thisProjection));

      if (thisProjection==projXY)
	Title +="X [mm]; Y [mm]";
      else
	Title +="; "+projectionName.at(static_cast<projAxis>(thisProjection))+" [mm];p_{bckgrd match}";
	

      // replace info with user settings from histomanager
      if ( isHistoManagerAvailable ){
	histoInfo = histoMgr->getHistogramInfo(thisHistoName);
	if ( histoInfo ) {
	  message<DEBUG5> ( log() << (* histoInfo ) );
	  NBinX = histoInfo->_xBin;
	  MinX  = histoInfo->_xMin;
	  MaxX  = histoInfo->_xMax;
	  NBinY = histoInfo->_yBin;
	  MinY  = histoInfo->_yMin;
	  MaxY  = histoInfo->_yMax;
	  if ( histoInfo->_title != "" ) Title = histoInfo->_title;
	}
      }
	
      AIDA::IBaseHistogram * thisHisto;
      if (static_cast<projAxis>(thisProjection)== projXY) 
	thisHisto = dynamic_cast<AIDA::IBaseHistogram*>(AIDAProcessor::histogramFactory(this)->createProfile2D(thisHistoName.c_str(),NBinX,MinX,MaxX,NBinY,MinY,MaxY));
      else 
	thisHisto = dynamic_cast<AIDA::IBaseHistogram*>(AIDAProcessor::histogramFactory(this)->createProfile1D(thisHistoName.c_str(),NBinX,MinX,MaxX));
	
      thisHisto->setTitle(Title);
      _BgEfficiencyHistos.insert(make_pair(static_cast<projAxis>(thisProjection), thisHisto) );

    } // background match prob. as fcn of fitted position

    { //  Noise as a function of the measured position
      int    NBinX  =  150 ; // defaults; 
      double MinX   =  -15.;//  defaults; 
      double MaxX   =  15.;//  defaults; 
      // for XY only:
      int    NBinY  =  100 ; // defaults; 
      double MinY   =  -10.;//  defaults; 
      double MaxY   =  10.;//  defaults; 

      if (static_cast<projAxis>(thisProjection)== projY){
	NBinX  =  100 ;
	MinX   =  -10.;
	MaxX   =  10.;
      }

      // construct the name of the histogram we are looking for:
      std::string thisHistoName = "";
      thisHistoName = _NoiseHistoBaseName + projectionName.at(static_cast<projAxis>(thisProjection));
	
      // set default title
      std::string Title = "Noise vs part pos in " + projectionName.at(static_cast<projAxis>(thisProjection));

      if (thisProjection==projXY)
	Title +="X [mm]; Y [mm]";
      else
	Title +="; "+projectionName.at(static_cast<projAxis>(thisProjection))+" [mm];# of noise hits";

	
      // replace info with user settings from histomanager
      if ( isHistoManagerAvailable ){
	histoInfo = histoMgr->getHistogramInfo(thisHistoName);
	if ( histoInfo ) {
	  message<DEBUG5> ( log() << (* histoInfo ) );
	  NBinX = histoInfo->_xBin;
	  MinX  = histoInfo->_xMin;
	  MaxX  = histoInfo->_xMax;
	  NBinY = histoInfo->_yBin;
	  MinY  = histoInfo->_yMin;
	  MaxY  = histoInfo->_yMax;
	  if ( histoInfo->_title != "" ) Title = histoInfo->_title;
	}
      }
	
      AIDA::IBaseHistogram * thisHisto = 0;
      if (static_cast<projAxis>(thisProjection)== projXY) 
	thisHisto = dynamic_cast<AIDA::IBaseHistogram*>(AIDAProcessor::histogramFactory(this)->createProfile2D(thisHistoName.c_str(),NBinX,MinX,MaxX,NBinY,MinY,MaxY));
      else 
	thisHisto = dynamic_cast<AIDA::IBaseHistogram*>(AIDAProcessor::histogramFactory(this)->createProfile1D(thisHistoName.c_str(),NBinX,MinX,MaxX));
	
      thisHisto->setTitle(Title);
      _NoiseHistos.insert(make_pair(static_cast<projAxis>(thisProjection), thisHisto) );

    } //  Noise as a function of the measured position

  } // projection axis loop

  // Measured - fitted position in X  vs Y
  int    shiftXNBin  = 150;
  double shiftXMin   = -15.;
  double shiftXMax   =  15.;
  double shiftVMin   = -1.0;
  double shiftVMax   =  1.0;
  std::string shiftXTitle = "Measured - fitted X position vs Y;Y; #Delta X [mm]";

  if ( isHistoManagerAvailable ){
    histoInfo = histoMgr->getHistogramInfo(_ShiftXvsYHistoName);
    if ( histoInfo ){
      message<DEBUG5> ( log() << (* histoInfo ) );
      shiftXNBin = histoInfo->_xBin;
      shiftXMin  = histoInfo->_xMin;
      shiftXMax  = histoInfo->_xMax;
      shiftVMin  = histoInfo->_yMin;
      shiftVMax  = histoInfo->_yMax;
      if ( histoInfo->_title != "" ) shiftXTitle = histoInfo->_title;
    }
  }

  _ShiftXvsYHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_ShiftXvsYHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftVMin,shiftVMax);
  _ShiftXvsYHisto->setTitle(shiftXTitle.c_str());

  // Measured - fitted position in X vs X
  shiftXTitle = "Measured - fitted X position vs X; X [mm];#Delta X [mm]"; 
  if ( isHistoManagerAvailable ){
    histoInfo = histoMgr->getHistogramInfo(_ShiftXvsXHistoName);
    if ( histoInfo ){
      message<DEBUG5> ( log() << (* histoInfo ) );
      shiftXNBin = histoInfo->_xBin;
      shiftXMin  = histoInfo->_xMin;
      shiftXMax  = histoInfo->_xMax;
      shiftVMin  = histoInfo->_yMin;
      shiftVMax  = histoInfo->_yMax;
      if ( histoInfo->_title != "" ) shiftXTitle = histoInfo->_title;
    }
  }
  
  _ShiftXvsXHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_ShiftXvsXHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftVMin,shiftVMax);
  _ShiftXvsXHisto->setTitle(shiftXTitle.c_str());


  // Measured - fitted position in Y vs X
  int shiftYNBin  = 100;
  double shiftYMin   = -10.;
  double shiftYMax   =  10.;
  shiftVMin   = -1.0;
  shiftVMax   =  1.0;
  std::string shiftYTitle = "Measured - fitted Y position vs X; X [mm]#Delta Y [mm];";

  if ( isHistoManagerAvailable ){
    histoInfo = histoMgr->getHistogramInfo(_ShiftYvsXHistoName);
    if ( histoInfo ){
      message<DEBUG5> ( log() << (* histoInfo ) );
      shiftYNBin = histoInfo->_xBin;
      shiftYMin  = histoInfo->_xMin;
      shiftYMax  = histoInfo->_xMax;
      shiftVMin  = histoInfo->_yMin;
      shiftVMax  = histoInfo->_yMax;
      if ( histoInfo->_title != "" ) shiftYTitle = histoInfo->_title;
    }
  }

  _ShiftYvsXHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_ShiftYvsXHistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax,shiftVMin,shiftVMax);

  _ShiftYvsXHisto->setTitle(shiftYTitle.c_str());


  // Measured - fitted position in Y vs Y
  shiftYTitle = "Measured - fitted Y position vs Y;Y [mm]; #Delta Y [mm]";
  if ( isHistoManagerAvailable ){
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

  _ShiftYvsYHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_ShiftYvsYHistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax,shiftVMin,shiftVMax);
  _ShiftYvsYHisto->setTitle(shiftYTitle.c_str());

  // Measured - fitted position in X  vs Y (2D plot)
  shiftXNBin  = 150;
  shiftXMin   = -15.;
  shiftXMax   =  15.; 
  int shiftVNBin  = 200;
  shiftVMin   = -1.0;
  shiftVMax   =  1.0;
  shiftXTitle = "Measured - fitted X position vs Y; Y [mm];#Delta X [mm]";

  if ( isHistoManagerAvailable ){
    histoInfo = histoMgr->getHistogramInfo(_ShiftXvsY2DHistoName);
    if ( histoInfo ){
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


  _ShiftXvsY2DHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D(_ShiftXvsY2DHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftVNBin,shiftVMin,shiftVMax);
  _ShiftXvsY2DHisto->setTitle(shiftXTitle.c_str());




  // Measured - fitted position in X vs X (2D plot) 
  shiftXTitle = "Measured - fitted X position vs X; X [mm]; #Delta X [mm]";

  if ( isHistoManagerAvailable ){
    histoInfo = histoMgr->getHistogramInfo(_ShiftXvsX2DHistoName);
    if ( histoInfo ){
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


  _ShiftXvsX2DHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D(_ShiftXvsX2DHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftVNBin,shiftVMin,shiftVMax);
  _ShiftXvsX2DHisto->setTitle(shiftXTitle.c_str());

  // Measured - fitted position in Y vs X  (2D plot)
  shiftYNBin  = 150;
  shiftYMin   = -15.;
  shiftYMax   =  15.;
  shiftVNBin  =  200;
  shiftVMin   = -1.0;
  shiftVMax   =  1.0;
  shiftYTitle = "Measured - fitted Y position vs X;  X [mm];#Delta Y [mm]";

  if ( isHistoManagerAvailable ){
    histoInfo = histoMgr->getHistogramInfo(_ShiftYvsX2DHistoName);
    if ( histoInfo ){
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

  _ShiftYvsX2DHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D(_ShiftYvsX2DHistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax,shiftVNBin,shiftVMin,shiftVMax);
  _ShiftYvsX2DHisto->setTitle(shiftYTitle.c_str());


  // Measured - fitted position in Y vs Y (2D plot)
  shiftYTitle = "Measured - fitted Y position vs Y;  Y [mm];#Delta Y [mm]";

  if ( isHistoManagerAvailable ){
    histoInfo = histoMgr->getHistogramInfo(_ShiftYvsY2DHistoName);
    if ( histoInfo ){
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

  _ShiftYvsY2DHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D(_ShiftYvsY2DHistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax,shiftVNBin,shiftVMin,shiftVMax);
  _ShiftYvsY2DHisto->setTitle(shiftYTitle.c_str());


  // Eta function check: measured - fitted position in X  vs  local X
  int etaXNBin  =  60;
  double etaXMin   = -0.03;
  double etaXMax   = 0.03;
  double etaVMin   = -0.03;
  double etaVMax   = 0.03;
  string etaXTitle = "Measured - fitted X position vs local X; local X [mm]; #Delta X";

  if ( isHistoManagerAvailable ){
    histoInfo = histoMgr->getHistogramInfo(_EtaXHistoName);
    if ( histoInfo ){
      message<DEBUG5> ( log() << (* histoInfo ) );
      etaXNBin = histoInfo->_xBin;
      etaXMin  = histoInfo->_xMin;
      etaXMax  = histoInfo->_xMax;
      etaVMin  = histoInfo->_yMin;
      etaVMax  = histoInfo->_yMax;
      if ( histoInfo->_title != "" ) etaXTitle = histoInfo->_title;
    }
  }


  _EtaXHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_EtaXHistoName.c_str(),etaXNBin,etaXMin,etaXMax,etaVMin,etaVMax);
  _EtaXHisto->setTitle(etaXTitle.c_str());

  // Eta function check: measured - fitted position in Y vs local Y
  int etaYNBin  =  60;
  double etaYMin   = -0.03;
  double etaYMax   = 0.03;
  etaVMin   = -0.03;
  etaVMax   = 0.03;
  string etaYTitle = "Measured - fitted Y position vs local Y;  local Y [mm];#Delta Y [mm]";

  if ( isHistoManagerAvailable ){
    histoInfo = histoMgr->getHistogramInfo(_EtaYHistoName);
    if ( histoInfo ){
      message<DEBUG5> ( log() << (* histoInfo ) );
      etaYNBin = histoInfo->_xBin;
      etaYMin  = histoInfo->_xMin;
      etaYMax  = histoInfo->_xMax;
      etaVMin  = histoInfo->_yMin;
      etaVMax  = histoInfo->_yMax;
      if ( histoInfo->_title != "" ) etaYTitle = histoInfo->_title;
    }
  }

  _EtaYHisto = AIDAProcessor::histogramFactory(this)->createProfile1D(_EtaYHistoName.c_str(),etaYNBin,etaYMin,etaYMax,etaVMin,etaVMax);
  _EtaYHisto->setTitle(etaYTitle.c_str());

  // Eta function check: measured - fitted position in X  vs local X (2D plot)
  etaXNBin  =  60;
  etaXMin   = -0.03;
  etaXMax   = 0.03;
  int etaVNBin  =  60;
  etaVMin   = -0.03;
  etaVMax   = 0.03;
  etaXTitle = "Measured - fitted X position vs local X; local x [mm]; #Delta X [mm]";

  if ( isHistoManagerAvailable ){
    histoInfo = histoMgr->getHistogramInfo(_EtaX2DHistoName);
    if ( histoInfo ){
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


  _EtaX2DHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D(_EtaX2DHistoName.c_str(),etaXNBin,etaXMin,etaXMax,etaVNBin,etaVMin,etaVMax);
  _EtaX2DHisto->setTitle(etaXTitle.c_str());

  // Measured - fitted position in Y vs Y  (2D plot)
  etaYNBin  =  60;
  etaYMin   = -0.03;
  etaYMax   = 0.03;
  etaVNBin  =  60;
  etaVMin   = -0.03;
  etaVMax   = 0.03;
  etaYTitle = "Measured - fitted Y position vs local Y; local Y [mm]; #Delta Y [mm]";

  if ( isHistoManagerAvailable ){
    histoInfo = histoMgr->getHistogramInfo(_EtaY2DHistoName);
    if ( histoInfo ){
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

  _EtaY2DHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D(_EtaY2DHistoName.c_str(),etaYNBin,etaYMin,etaYMax,etaVNBin,etaVMin,etaVMax);
  _EtaY2DHisto->setTitle(etaYTitle.c_str());




  // Eta function check: measured - fitted position in X  vs  local X-Y ("3D")

  etaXNBin  =  60;
  etaXMin   = -0.03;
  etaXMax   = 0.03;
  etaYNBin  =  60;
  etaYMin   = -0.03;
  etaYMax   = 0.03;
  etaVMin   = -0.03;
  etaVMax   = 0.03;
  etaXTitle = "Measured - fitted X position vs local X-Y; XY [mm];#Delta X [mm]";


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


  _EtaX3DHisto = AIDAProcessor::histogramFactory(this)->createProfile2D(_EtaX3DHistoName.c_str(),etaXNBin,etaXMin,etaXMax,etaYNBin,etaYMin,etaYMax,etaVMin,etaVMax);
  _EtaX3DHisto->setTitle(etaXTitle.c_str());


  // Eta function check: measured - fitted position in Y vs local X-Y ("3D")
  etaXNBin  =  60;
  etaXMin   = -0.03;
  etaXMax   = 0.03;
  etaYNBin  =  60;
  etaYMin   = -0.03;
  etaYMax   = 0.03;
  etaVMin   = -0.03;
  etaVMax   = 0.03;
  etaYTitle = "Measured - fitted Y position vs local X-Y; XY [mm];#Delta Y [mm]";

  if ( isHistoManagerAvailable ){
    histoInfo = histoMgr->getHistogramInfo(_EtaY3DHistoName);
    if ( histoInfo ){
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

  _EtaY3DHisto = AIDAProcessor::histogramFactory(this)->createProfile2D(_EtaY3DHistoName.c_str(),etaXNBin,etaXMin,etaXMax,etaYNBin,etaYMin,etaYMax,etaVMin,etaVMax);
  _EtaY3DHisto->setTitle(etaYTitle.c_str());


  // Pixel plots
  Int_t     pixYNBin  =  128;
  Double_t  pixYMin   = -18.4*2.;
  Double_t  pixYMax   =  18.4*2.;
  Int_t     pixXNBin  =  128;
  Double_t  pixXMin   = -18.4*2.;
  Double_t  pixXMax   =  18.4*2.;
  Double_t  pixVMin   = -1000000.;
  Double_t  pixVMax   =  1000000.;

  string pixTitle     = "In pixel efficiency; efficiency; pixel";
  // ---- // Efficiency
  if ( isHistoManagerAvailable ){
    histoInfo = histoMgr->getHistogramInfo( _PixelEfficiencyHistoName );
    if ( histoInfo ){
      message<DEBUG5> ( log() << (* histoInfo ) );
      pixXNBin = histoInfo->_xBin;
      pixXMin  = histoInfo->_xMin;
      pixXMax  = histoInfo->_xMax;
      pixYNBin = histoInfo->_yBin;
      pixYMin  = histoInfo->_yMin;
      pixYMax  = histoInfo->_yMax;
      pixVMin  = histoInfo->_zMin;
      pixVMax  = histoInfo->_zMax;
      if ( histoInfo->_title != "" ) pixTitle = histoInfo->_title;
    }
  }

  _PixelEfficiencyHisto = AIDAProcessor::histogramFactory(this)->createProfile2D(_PixelEfficiencyHistoName.c_str(), pixYNBin, pixYMin, pixYMax, pixXNBin, pixXMin, pixXMax, pixVMin, pixVMax );
  _PixelEfficiencyHisto->setTitle( pixTitle.c_str());


  // ---- // Resolution X
  if ( isHistoManagerAvailable ){
    histoInfo = histoMgr->getHistogramInfo( _PixelResolutionXHistoName );
    if ( histoInfo ){
      message<DEBUG5> ( log() << (* histoInfo ) );
      pixXNBin = histoInfo->_xBin;
      pixXMin  = histoInfo->_xMin;
      pixXMax  = histoInfo->_xMax;
      pixYNBin = histoInfo->_yBin;
      pixYMin  = histoInfo->_yMin;
      pixYMax  = histoInfo->_yMax;
      pixVMin  = histoInfo->_zMin;
      pixVMax  = histoInfo->_zMax;
      if ( histoInfo->_title != "" ) pixTitle = histoInfo->_title;
    }
  }

  _PixelResolutionXHisto = AIDAProcessor::histogramFactory(this)->createProfile2D(_PixelResolutionXHistoName.c_str(), pixYNBin, pixYMin, pixYMax, pixXNBin, pixXMin, pixXMax, pixVMin, pixVMax );
  _PixelResolutionXHisto->setTitle( pixTitle.c_str());


  // ---- // Resolution Y
  if ( isHistoManagerAvailable ){
    histoInfo = histoMgr->getHistogramInfo( _PixelResolutionYHistoName );
    if ( histoInfo ){
      message<DEBUG5> ( log() << (* histoInfo ) );
      pixXNBin = histoInfo->_xBin;
      pixXMin  = histoInfo->_xMin;
      pixXMax  = histoInfo->_xMax;
      pixYNBin = histoInfo->_yBin;
      pixYMin  = histoInfo->_yMin;
      pixYMax  = histoInfo->_yMax;
      pixVMin  = histoInfo->_zMin;
      pixVMax  = histoInfo->_zMax;
      if ( histoInfo->_title != "" ) pixTitle = histoInfo->_title;
    }
  }

  _PixelResolutionYHisto = AIDAProcessor::histogramFactory(this)->createProfile2D(_PixelResolutionYHistoName.c_str(), pixYNBin, pixYMin, pixYMax, pixXNBin, pixXMin, pixXMax, pixVMin, pixVMax );
  _PixelResolutionYHisto->setTitle( pixTitle.c_str());

  // ---- // Charge Sharing 
  if ( isHistoManagerAvailable ){
    histoInfo = histoMgr->getHistogramInfo( _PixelChargeSharingHistoName );
    if ( histoInfo ){
      message<DEBUG5> ( log() << (* histoInfo ) );
      pixXNBin = histoInfo->_xBin;
      pixXMin  = histoInfo->_xMin;
      pixXMax  = histoInfo->_xMax;
      pixYNBin = histoInfo->_yBin;
      pixYMin  = histoInfo->_yMin;
      pixYMax  = histoInfo->_yMax;
      pixVMin  = histoInfo->_zMin;
      pixVMax  = histoInfo->_zMax;
      if ( histoInfo->_title != "" ) pixTitle = histoInfo->_title;
    }
  }

  _PixelChargeSharingHisto = AIDAProcessor::histogramFactory(this)->createProfile2D(_PixelChargeSharingHistoName.c_str(), pixYNBin, pixYMin, pixYMax, pixXNBin, pixXMin, pixXMax, pixVMin, pixVMax );
  _PixelChargeSharingHisto->setTitle( pixTitle.c_str());


  message<DEBUG5> ( log() << "Histogram booking completed \n\n");
#else
  message<MESSAGE5> ( log() << "No histogram produced because Marlin doesn't use AIDA" );
#endif

  return;
}


int EUTelDUTHistograms::getClusterSize(int sensorID, TrackerHit * hit, int& sizeX, int& sizeY, int& subMatrix ) {

  if(hit==0)
  {
    streamlog_out( ERROR5 ) << "An invalid hit pointer supplied! will exit now\n" << endl;
    return -1;
  }

        try
        {
            TrackerDataImpl* clusterVector = static_cast<TrackerDataImpl*>( hit->getRawHits()[0]);

            EUTelSimpleVirtualCluster * cluster=0;

            if ( hit->getType() == kEUTelBrickedClusterImpl ) {

               // fixed cluster implementation. Remember it
               //  can come from
               //  both RAW and ZS data
   
                cluster = new EUTelBrickedClusterImpl( clusterVector );
                
            } else if ( hit->getType() == kEUTelDFFClusterImpl ) {
              
              // fixed cluster implementation. Remember it can come from
              // both RAW and ZS data
              cluster = new EUTelDFFClusterImpl( clusterVector );
            } else if ( hit->getType() == kEUTelFFClusterImpl ) {
              
              // fixed cluster implementation. Remember it can come from
              // both RAW and ZS data
              cluster = new EUTelFFClusterImpl( clusterVector );
            }
            else if ( hit->getType() == kEUTelSparseClusterImpl ) 
            {
               cluster = new EUTelSparseClusterImpl< EUTelGenericSparsePixel > ( clusterVector );
            }
	    else if ( hit->getType() ==	kEUTelGenericSparseClusterImpl )
	{

		CellIDDecoder<TrackerDataImpl> cellDecoder (EUTELESCOPE::ZSDATADEFAULTENCODING);
		SparsePixelType   type   = static_cast<SparsePixelType> ( static_cast<int> (cellDecoder( clusterVector )["sparsePixelType"]) );
		
		if(type == kEUTelGenericSparsePixel )
		{
			cluster = new EUTelGenericSparseClusterImpl<EUTelGenericSparsePixel>( clusterVector );
		}
		else if(type == kEUTelGeometricPixel )
		{
			cluster = new EUTelGenericSparseClusterImpl<EUTelGeometricPixel>( clusterVector );
		}
	}

            if(cluster != 0)
            {
              float xlocal=-1.;
              float ylocal=-1.;
              cluster->getClusterSize(sizeX,sizeY);
              cluster->getCenterOfGravity(xlocal, ylocal);
              subMatrix = getSubMatrix(sensorID, xlocal);  
              return 0;         
            }
          }
          catch(...)
          {
            streamlog_out(ERROR5)<< "guess SensorID crashed" << std::endl;
          }

return -1;
}

int EUTelDUTHistograms::getSubMatrix(int detectorID, float xlocal)
{
   // quarters : 0-287, 288-575, 576-863, 864-1151
   int fourlocal = static_cast<int>(xlocal*4.); 
   int subquarter =  fourlocal/geo::gGeometry().siPlaneXNpixels(detectorID);
   return subquarter;       
}  

// -------------------------------------------------------------------------------------------
int EUTelDUTHistograms::read_track_from_collections(LCEvent *event)
{
  // Clear local fit storage tables
  // 'fitted' table used for comparison with measured hits
  // 'bgfitted' used for comparison with hits from previous event

 //
  // Get input collections
  //

  // first check if there are tracks:
  LCCollection* trackcol;
  try {
    trackcol = event->getCollection( _inputTrackColName ) ;
  } catch (lcio::DataNotAvailableException& e) {
    streamlog_out( DEBUG5 ) << "Not able to get collection " << _inputTrackColName << " from event " << event->getEventNumber() << " in run " << event->getRunNumber() <<  endl;
    return 1;
  }
  
  // setup cellIdDecoder to decode the sensor ID from the hits
  CellIDDecoder<TrackerHit>  hitCellDecoder(EUTELESCOPE::HITENCODING);
  CellIDDecoder<SimTrackerHitImpl>  simhitCellDecoder(EUTELESCOPE::HITENCODING);


  // Loop over tracks in input track collection
  // Read fitted positions at DUT

  int nTrack = trackcol->getNumberOfElements()  ;

  message<DEBUG5> ( log() << "Total of " << nTrack << " tracks in input collection " );

  if(nTrack == 0 ) return 1;

  _fittedX.clear();
  _fittedY.clear();

  _bgfittedX.clear();
  _bgfittedY.clear();

  _localX.clear();
  _localY.clear();

  _trackhitposX.clear();   
  _trackhitposY.clear();
  _trackhitsizeX.clear();
  _trackhitsizeY.clear();
  _trackhitsubM.clear();
  _trackhitsensorID.clear();

  LCCollection* fit__col;
  try {
    fit__col = event->getCollection( _inputFitHitColName ) ;
  } catch (lcio::DataNotAvailableException& e) {
    streamlog_out( DEBUG5 ) << "Not able to get collection " << _inputFitHitColName << " from event " << event->getEventNumber() << " in run " << event->getRunNumber() <<  endl;
    return 1;
  }

  LCCollection* rec__col;
  try {
    rec__col = event->getCollection( _inputRecHitColName ) ;
  } catch (lcio::DataNotAvailableException& e) {
    streamlog_out( DEBUG5 ) << "Not able to get collection " << _inputRecHitColName << " from event " << event->getEventNumber() << " in run " << event->getRunNumber() <<  endl;
    return 2;
  }


  // Loop over tracks in input track collection
  // Read fitted positions at DUT

  int nTracks = fit__col->getNumberOfElements()  ;
  int nRecHits = rec__col->getNumberOfElements()  ;

  if(streamlog_level(DEBUG5)){
    message<DEBUG5> ( log() << "\n tracks " << nTracks << " \n hits " <<  nRecHits << endl );
    message<DEBUG5> ( log() << "Total of " << nTracks  <<" (fit) hits in input collection " );
    message<DEBUG5> ( log() << "Total of " << nRecHits <<" (rec) hits  in input collection " );
  }

// looking through a track info:
// initialise a class-internal track counter:
  _maptrackid = 0;
  for(int itrack=0; itrack< nTracks ; itrack++)
    {

      const double * pos = 0;
      int hsensorID      = 0;

     
      TrackerHit* fithit = dynamic_cast<TrackerHit*>( fit__col->getElementAt(itrack) ) ;
      SimTrackerHitImpl *fithit0 = dynamic_cast<SimTrackerHitImpl*>( fit__col->getElementAt(itrack) ) ;
 
      if(fithit != 0 ) 
      { 
        pos       = fithit->getPosition();
        hsensorID = hitCellDecoder(fithit)["sensorID"];
      } 
      else 
      if(fithit == 0 )
      {
        if(fithit0 != 0 ) 
        { 
          pos       = fithit0->getPosition();
          hsensorID = simhitCellDecoder(fithit0)["sensorID"];
        } 
      }
  
      // skip if for some reason the track collection is at NULL address
      if( fithit == 0 && fithit0 == 0 ) continue;

            {

              if( hsensorID == _iDUT  )  // get all fitted hits on board
                {

                  _fittedX[_maptrackid].push_back(pos[0]);
                  _fittedY[_maptrackid].push_back(pos[1]);
                  _bgfittedX[_maptrackid].push_back(pos[0]);
                  _bgfittedY[_maptrackid].push_back(pos[1]);

//
// using fitted position to calculate in pixel coordinates

          double locX = pos[0];
          double locY = pos[1];

          // Subtract position of the central pixel

          int picX = static_cast<int>(locX/_pitchX);

          if(locX<0)picX--;

          locX-=(picX+0.5)*_pitchX;

          int picY = static_cast<int>(locY/_pitchY);

          if(locY<0)picY--;

          locY-=(picY+0.5)*_pitchY;

          _localX[_maptrackid].push_back(locX);
          _localY[_maptrackid].push_back(locY);

	  if(streamlog_level(DEBUG5)){
	    message<DEBUG5> ( log() << "_fittedX element [" << _fittedX[_maptrackid].size()-1 <<  "]" << _fittedX[_maptrackid][ _fittedX.size()-1] << " " << _fittedY[_maptrackid][ _fittedX.size()-1] << " for DUT " << hsensorID << endl);
	  }

                }
            }

      // End of loop over fitted tracks
      _maptrackid++; 
   }

  if(streamlog_level(DEBUG5))
  {
    for(int ii=0;ii<_maptrackid;ii++)
    {
      message<MESSAGE5> ( log() << "for _maptrackid=" << ii << " found fithits " << _fittedX[ii].size() <<  endl);      
      for(unsigned int jj=0; jj < _fittedX[ii].size(); jj++)
      { 
         message<MESSAGE5> ( log() << "fit hits [" << jj << " of " << _fittedX[ii].size() << "] " << _fittedX[ii][jj] << " " <<  _fittedY[ii][jj] <<  endl);             
      }   
    }
  }

  message<DEBUG5> ( log() << "rechits " << endl );

   // Clear local tables with measured position
  _clusterSizeX.clear();
  _clusterSizeY.clear();
  _subMatrix.clear();


  _measuredX.clear();
  _measuredY.clear();

  _bgmeasuredX.clear();
  _bgmeasuredY.clear();
         // look at reconstructed hits only for all planes (excluding DUT !)       
          // hits that belong to track "_maptrackid"
          for(int ihit=0; ihit< nRecHits ; ihit++)
            {
              const double * pos = 0;
              int hsensorID      = 0;

              TrackerHit * meshit = dynamic_cast<TrackerHit*>( rec__col->getElementAt(ihit) ) ;
              SimTrackerHitImpl * meshit0 = dynamic_cast<SimTrackerHitImpl*>( rec__col->getElementAt(ihit) ) ;
 
      if(meshit != 0 ) 
      { 
        pos       = meshit->getPosition();
        hsensorID = hitCellDecoder(meshit)["sensorID"];
      } 
      else 
      if(meshit == 0 )
      {
        if(meshit0 != 0 ) 
        { 
          pos       = meshit0->getPosition();
          hsensorID = simhitCellDecoder(meshit0)["sensorID"];
        } 
      }
  
      // skip if for some reason the track collection is at NULL address
      if( meshit == 0 && meshit0 == 0 ) continue;


              if(hsensorID == _iDUT   ) 
                {
                  int sizeX = -1;
                  int sizeY = -1;
                  int subMatrix = -1; 
                  if(meshit != 0 )
                  {
                     getClusterSize( hsensorID, static_cast<TrackerHit*>(meshit), sizeX, sizeY, subMatrix);
                  } 
                  else if(meshit0 != 0 ) 
                  {
                     sizeX = 1;     
                     sizeY = 1;     
                     subMatrix = 0;     
                  }
		  _clusterSizeX.push_back(sizeX);
		  _clusterSizeY.push_back(sizeY);
		  _subMatrix.push_back( subMatrix );
		  
		  _measuredX.push_back(pos[0]);
		  _measuredY.push_back(pos[1]);
		  _bgmeasuredX.push_back(pos[0]);
		  _bgmeasuredY.push_back(pos[1]);
		  if(streamlog_level(DEBUG5)){
		    
		    message<DEBUG5> ( log() << "_measured element [" << _measuredX.size()-1 <<  "]"  << _measuredX[ _measuredX.size()-1] << " " << _measuredY[ _measuredX.size()-1] << " for DUT " << hsensorID << endl);
		  }
		}
            }

    message<DEBUG5> ( log() << "Total of " << _measuredX.size() << " hits at DUT " << _iDUT << endl);
      

 return 0;
}


// -----------------------------------------------------------------------------------------------------------
int EUTelDUTHistograms::read_track(LCEvent *event)
{
  // Clear local fit storage tables
  // 'fitted' table used for comparison with measured hits
  // 'bgfitted' used for comparison with hits from previous event

  _fittedX.clear();
  _fittedY.clear();

  _bgfittedX.clear();
  _bgfittedY.clear();

  _localX.clear();
  _localY.clear();

  _trackhitposX.clear();   
  _trackhitposY.clear();
  _trackhitsizeX.clear();
  _trackhitsizeY.clear();
  _trackhitsubM.clear();
  _trackhitsensorID.clear();

 //
  // Get input collections
  //

  LCCollection* trackcol;
  try {
    trackcol = event->getCollection( _inputTrackColName ) ;
  } catch (lcio::DataNotAvailableException& e) {
    streamlog_out( DEBUG5 ) << "Not able to get collection " << _inputTrackColName << " from event " << event->getEventNumber() << " in run " << event->getRunNumber() <<  endl;
    return 1;
  }

  // setup cellIdDecoder to decode the sensor ID from the hits
  CellIDDecoder<TrackerHit>  hitCellDecoder(EUTELESCOPE::HITENCODING);


  // Loop over tracks in input track collection
  // Read fitted positions at DUT

  int nTrack = trackcol->getNumberOfElements()  ;


  message<DEBUG5> ( log() << "Total of " << nTrack << " tracks in input collection " );

// looking through a track info:
// initialise a class-internal track counter:
  _maptrackid = 0;
  for(int itrack=0; itrack< nTrack ; itrack++)
    {
      Track * fittrack = dynamic_cast<Track*>( trackcol->getElementAt(itrack) ) ;
 
      // skip if for some reason the track collection is at NULL address
      if( fittrack == 0 ) continue;

        {
          // Look at hits assigned to track
          std::vector<EVENT::TrackerHit*>  trackhits = fittrack->getTrackerHits();

          int nHit =   trackhits.size();

          for(int ihit=0; ihit< nHit ; ihit++)
            {
              TrackerHit * meshit = trackhits.at(ihit);

              // Hit position

              const double * pos = meshit->getPosition();
              int hsensorID = hitCellDecoder(meshit)["sensorID"];

              // Look at fitted hits only!

              if( ((hitCellDecoder(meshit)["properties"] & kFittedHit) > 0)  && hsensorID == _iDUT  )  // get all fitted hits on board
                {

                  _fittedX[_maptrackid].push_back(pos[0]);
                  _fittedY[_maptrackid].push_back(pos[1]);
                  _bgfittedX[_maptrackid].push_back(pos[0]);
                  _bgfittedY[_maptrackid].push_back(pos[1]);

//
// using fitted position to calculate in pixel coordinates

          double locX = pos[0];
          double locY = pos[1];

          // Subtract position of the central pixel

          int picX = static_cast<int>(locX/_pitchX);

          if(locX<0)picX--;

          locX-=(picX+0.5)*_pitchX;

          int picY = static_cast<int>(locY/_pitchY);

          if(locY<0)picY--;

          locY-=(picY+0.5)*_pitchY;

          _localX[_maptrackid].push_back(locX);
          _localY[_maptrackid].push_back(locY);

                  break;
                }
            }
  
          // look at reconstructed hits only for all planes (excluding DUT !)       
          // hits that belong to track "_maptrackid"
          for(int ihit=0; ihit< nHit ; ihit++)
            {
              TrackerHit * meshit = trackhits.at(ihit);

              // Hit position
              const double * pos = meshit->getPosition();
              int hsensorID = hitCellDecoder(meshit)["sensorID"];

              if( (hitCellDecoder(meshit)["properties"] & kFittedHit) == 0  )
                {
                  int sizeX = -1;
                  int sizeY = -1;
                  int subMatrix = -1; 
                  getClusterSize( hsensorID, static_cast<TrackerHitImpl*>(meshit), sizeX, sizeY, subMatrix);
                 _trackhitposX[_maptrackid].push_back(pos[0]);
                  _trackhitposY[_maptrackid].push_back(pos[1]);
                  _trackhitsizeX[_maptrackid].push_back( sizeX );
                  _trackhitsizeY[_maptrackid].push_back( sizeY );
                  _trackhitsubM[_maptrackid].push_back(subMatrix);
                  _trackhitsensorID[_maptrackid].push_back( hsensorID );
                }
            }
        }

      // End of loop over fitted tracks
      _maptrackid++; 
   }

   // Clear local tables with measured position
  _clusterSizeX.clear();
  _clusterSizeY.clear();
  _subMatrix.clear();


  _measuredX.clear();
  _measuredY.clear();

  _bgmeasuredX.clear();
  _bgmeasuredY.clear();
 
  LCCollection* hitcol = NULL;
  try {
    hitcol = event->getCollection( _inputHitColName ) ;
  } catch (lcio::DataNotAvailableException& e) {
    streamlog_out( DEBUG5 ) << "Not able to get collection " << _inputHitColName << " from event " << event->getEventNumber() << " in run " << event->getRunNumber() <<  endl;
  }

  int nHit = 0;

  nHit = hitcol->getNumberOfElements();

  message<DEBUG5> ( log() << "Total of " << nHit << " tracker hits in input collection " );


  for(int ihit=0; ihit< nHit ; ihit++)
    {
      TrackerHit * meshit = dynamic_cast<TrackerHit*>( hitcol->getElementAt(ihit) ) ;

      // Hit position

      const double * pos = meshit->getPosition();

      int   sensorID = hitCellDecoder(meshit)["sensorID"];

     if( sensorID ==_iDUT ) // measured info only for DUT plane
        {
          int sizeX = -1;
          int sizeY = -1;
          int subMatrix = -1;
          getClusterSize( sensorID, static_cast<TrackerHitImpl*>(meshit), sizeX, sizeY, subMatrix);
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

       }
    }

    message<DEBUG5> ( log() << "Total of " << _measuredX.size() << " hits at DUT " << _iDUT << endl);

 return 0;
}


