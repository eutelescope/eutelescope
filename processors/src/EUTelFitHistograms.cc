
// Author: A.F.Zarnecki, University of Warsaw <mailto:zarnecki@fuw.edu.pl>
// Version: $Id$
// Date 2007.09.10

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// this processor is built only if USE_GEAR and USE_AIDA
#if defined(USE_GEAR) && ( defined(USE_AIDA) || defined(MARLIN_USE_AIDA) )

// eutelescope includes
#include "EUTelFitHistograms.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelHistogramManager.h"
#include "EUTelExceptions.h"

// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/IProfile2D.h>
#include <AIDA/ITree.h>


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
std::string EUTelFitHistograms::_ShiftXvsYHistoName    = "deltaXvsY";
std::string EUTelFitHistograms::_ShiftYvsXHistoName    = "deltaYvsX";

std::string EUTelFitHistograms::_MeasuredXHistoName  = "measuredX";
std::string EUTelFitHistograms::_MeasuredYHistoName  = "measuredY";
std::string EUTelFitHistograms::_MeasuredXYHistoName = "measuredXY";

std::string EUTelFitHistograms::_FittedXHistoName  = "fittedX";
std::string EUTelFitHistograms::_FittedYHistoName  = "fittedY";
std::string EUTelFitHistograms::_FittedXYHistoName = "fittedXY";

std::string EUTelFitHistograms::_ResidualXHistoName = "residualX";
std::string EUTelFitHistograms::_ResidualYHistoName = "residualY";
std::string EUTelFitHistograms::_ResidualXYHistoName = "residualXY";

std::string EUTelFitHistograms::_ScatXHistoName  = "scatteringX";
std::string EUTelFitHistograms::_ScatYHistoName  = "scatteringY";
std::string EUTelFitHistograms::_ScatXYHistoName = "scatteringXY";

std::string EUTelFitHistograms::_AngleXHistoName = "incangleX";
std::string EUTelFitHistograms::_AngleYHistoName = "incangleY";
std::string EUTelFitHistograms::_AngleXYHistoName = "incangleXY";

std::string EUTelFitHistograms::_clusterSignalHistoName  = "clusterSignal";
std::string EUTelFitHistograms::_meanSignalXHistoName  = "meanSignalX";
std::string EUTelFitHistograms::_meanSignalYHistoName  = "meanSignalY";
std::string EUTelFitHistograms::_meanSignalXYHistoName  = "meanSignalXY";

std::string EUTelFitHistograms::_beamShiftXHistoName  = "beamShiftX";
std::string EUTelFitHistograms::_beamShiftYHistoName  = "beamShiftY";
std::string EUTelFitHistograms::_beamShiftXYHistoName = "beamShiftXY";
std::string EUTelFitHistograms::_beamRotXHistoName    = "beamRotX";
std::string EUTelFitHistograms::_beamRotYHistoName    = "beamRotY";
std::string EUTelFitHistograms::_beamRotX2DHistoName  = "beamRotX2D";
std::string EUTelFitHistograms::_beamRotY2DHistoName  = "beamRotY2D";
std::string EUTelFitHistograms::_beamRot2XHistoName   = "beamRot2X";
std::string EUTelFitHistograms::_beamRot2YHistoName   = "beamRot2Y";

std::string EUTelFitHistograms::_relShiftXHistoName = "relShiftX";
std::string EUTelFitHistograms::_relShiftYHistoName = "relShiftY";
std::string EUTelFitHistograms::_relRotXHistoName   = "relRotX";
std::string EUTelFitHistograms::_relRotYHistoName   = "relRotY";
std::string EUTelFitHistograms::_relRotX2DHistoName   = "relRotX2D";
std::string EUTelFitHistograms::_relRotY2DHistoName   = "relRotY2D";

EUTelFitHistograms::EUTelFitHistograms() : Processor("EUTelFitHistograms") {

  // modify processor description
  _description = "Histogram track fit results" ;


  // register steering parameters:
  //       name, description, class-variable, default value

  // input collection first:

  registerInputCollection( LCIO::TRACK,
                           "InputCollectionName" ,
                           "Name of the input Track collection"  ,
                           _inputColName ,
                           std::string("testfittracks") ) ;

  // other processor parameters:


  registerProcessorParameter("HistoInfoFileName",
                             "Name of the histogram information file",
                             _histoInfoFileName, string( "histoinfo.xml" ) );

  registerProcessorParameter ("AlignCheckHistograms",
                              "Flag for producing additional histograms for alignment consistency check",
                              _alignCheckHistograms,  static_cast < bool > (false));

  registerProcessorParameter ("BeamReferenceID",
                              "ID of the layer used for beam based alignment check",
                              _BeamReferenceID,  static_cast < int > (0));

  std::vector<int > initLayerIDs;
  initLayerIDs.push_back(0);
  initLayerIDs.push_back(1);

  registerProcessorParameter ("TelescopeReferenceIDs",
                              "IDs of two layers used to check internal telescope alignment",
                              _TelescopeReferenceIDs, initLayerIDs);

  registerProcessorParameter ("DebugEventCount",
                              "Print out every DebugEnevtCount event",
                              _debugCount,  static_cast < int > (100));


}


void EUTelFitHistograms::init() {

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;


  // check if the GEAR manager pointer is not null!
  if ( Global::GEAR == 0x0 ) {
    streamlog_out( ERROR5 ) <<  "The GearMgr is not available, for an unknown reason." << endl;
    exit(-1);
  }

  // Read geometry information from GEAR
  streamlog_out ( MESSAGE5 )  << "Reading telescope geometry description from GEAR " << endl;

  _siPlanesParameters  = const_cast<gear::SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));


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

// Read position in Z (for sorting)

  _planeSort = new int[_nTelPlanes];
  _planePosition   = new double[_nTelPlanes];

  for(int ipl=0; ipl <  _siPlanesLayerLayout->getNLayers(); ipl++)
    {
      _planePosition[ipl]=_siPlanesLayerLayout->getLayerPositionZ(ipl);
      _planeSort[ipl]=ipl;
    }

  if(_iDUT>0)
    {
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
  _isActive        = new bool[_nTelPlanes];

  _beamID=-1;
  _referenceID0=-1;
  _referenceID1=-1;

// Fill remaining layer parameters

  for(int iz=0; iz < _nTelPlanes ; iz++)
    {
      int ipl=_planeSort[iz];

      double resolution;

      if(ipl != _iDUT )
        {
          _planeID[iz]=_siPlanesLayerLayout->getID(ipl);
          resolution = _siPlanesLayerLayout->getSensitiveResolution(ipl);
        }
      else
        {
          _planeID[iz]=_siPlanesLayerLayout->getDUTID();
          resolution = _siPlanesLayerLayout->getDUTSensitiveResolution();
        }

      _isActive[iz] = (resolution > 0);

      // Set local indexes of reference planes

      if(_planeID[iz]==_BeamReferenceID)_beamID=iz;
      if(_planeID[iz]==_TelescopeReferenceIDs.at(0))_referenceID0=iz;
      if(_planeID[iz]==_TelescopeReferenceIDs.at(1))_referenceID1=iz;

    }

  // Get new DUT position (after sorting)

  for(int iz=0;iz< _nTelPlanes ; iz++)
    if(_planeSort[iz]==_iDUT)
      {
        _iDUT=iz;
        break;
      }


  // Print out geometry information

  streamlog_out ( MESSAGE5 ) << "Telescope configuration with " << _nTelPlanes << " planes" << endl;


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
         << "  at Z [mm] = " << _planePosition[ipl];

      streamlog_out ( MESSAGE5 ) << ss.str() << endl;
    }


  // Allocate arrays for track fitting

  _isMeasured      = new bool[_nTelPlanes];
  _isFitted        = new bool[_nTelPlanes];

  _measuredX     = new double[_nTelPlanes];
  _measuredY     = new double[_nTelPlanes];
  _measuredQ     = new double[_nTelPlanes];
  _fittedX       = new double[_nTelPlanes];
  _fittedY       = new double[_nTelPlanes];


// Book histograms

  bookHistos();

}

void EUTelFitHistograms::processRunHeader( LCRunHeader* runHeader) {
  std::unique_ptr<EUTelRunHeaderImpl> eutelHeader = std::make_unique<EUTelRunHeaderImpl>(runHeader);
  eutelHeader->addProcessor(type());
  _nRun++ ;

  // Decode and print out Run Header information - just a check

  int runNr = runHeader->getRunNumber();

  streamlog_out( MESSAGE5 ) << "Processing run header " << _nRun
                           << ", run nr " << runNr << endl;

  const std::string detectorName = runHeader->getDetectorName();
  const std::string detectorDescription = runHeader->getDescription();
  const std::vector<std::string> * subDets = runHeader->getActiveSubdetectors();

  streamlog_out( MESSAGE5 ) << detectorName << " : " << detectorDescription << endl;

  int nDet = subDets->size();

  if(nDet) {
    streamlog_out( MESSAGE5 ) << nDet << " subdetectors defined :" << endl;
  }
  stringstream ss;
  for(int idet=0;idet<nDet;idet++)  {
    streamlog_out( MESSAGE5 )  << idet+1 << " : " << subDets->at(idet) << endl;
  }


}

void EUTelFitHistograms::processEvent( LCEvent * event ) {

  EUTelEventImpl * euEvent = static_cast<EUTelEventImpl*> ( event );
  if ( euEvent->getEventType() == kEORE ) {
    streamlog_out( DEBUG5 ) << "EORE found: nothing else to do." << endl;
    return;
  }

  bool debug = ( _debugCount>0 && _nEvt%_debugCount == 0);

  _nEvt ++ ;

  LCCollection* col;
  try {
    col = event->getCollection( _inputColName ) ;
  } catch (lcio::DataNotAvailableException& e) {
    streamlog_out( WARNING ) << "Not able to get collection "
                             << _inputColName
                             << " from event " << event->getEventNumber()
                             << " in run " << event->getRunNumber()  << endl;
    return;
  }

  // Loop over tracks in input collections

  int nTrack = col->getNumberOfElements()  ;

  if(debug) {
    streamlog_out ( TESTFITTERMESSAGE )  << "Total of " << nTrack << " tracks in input collection " << endl;
  }


  for(int itrack=0; itrack< nTrack ; itrack++)
    {
      Track * fittrack = dynamic_cast<Track*>( col->getElementAt(itrack) ) ;

      // Hit list assigned to track

      std::vector<EVENT::TrackerHit*>  trackhits = fittrack->getTrackerHits();

      // Copy hits assign to the track to local table
      // Assign hits to sensor planes


      int nHit =   trackhits.size();



      if(debug){
        streamlog_out ( TESTFITTERMESSAGE )  << "Track " << itrack << " with " << nHit << " hits " << endl;
      }


      // Clear plane tables

      for(int ipl=0;ipl<_nTelPlanes;ipl++)
        {
          _isMeasured[ipl]=false;
          _isFitted[ipl]=false;
        }

      // Loop over hits and fill hit tables

      // setup cellIdDecoder to decode the hit properties
      CellIDDecoder<TrackerHit>  hitCellDecoder(EUTELESCOPE::HITENCODING);

      for(int ihit=0; ihit< nHit ; ihit++)
        {
          TrackerHit * meshit = trackhits.at(ihit);

          // Hit position

          const double * pos = meshit->getPosition();

          // We find plane number of the hit
          // by looking at the Z position

          double distMin = 3.;
          int hitPlane = -1 ;

          for(int ipl=0;ipl<_nTelPlanes;ipl++)
            {
              double dist =  pos[2] - _planePosition[ipl] ;

              if(dist*dist < distMin*distMin)
                {
                  hitPlane=ipl;
                  distMin=dist;
                }
            }

          // Ignore hits not matched to any plane

          if(hitPlane<0)
            {
              streamlog_out ( ERROR5 )  << "Hit outside telescope plane at z [mm] = "  << pos[2] << endl;
              continue;
            }


          if( (hitCellDecoder(meshit)["properties"] & kFittedHit) == 0 )
            {
              // Measured hits

              _isMeasured[hitPlane]=true;

              _measuredX[hitPlane]=pos[0];
              _measuredY[hitPlane]=pos[1];

              // Get cluster charge
              EVENT::LCObjectVec rawdata =  meshit->getRawHits();

              if(rawdata.size()>0)
                {
                  EUTelVirtualCluster * cluster = new EUTelFFClusterImpl ( static_cast<TrackerDataImpl*> (rawdata.at(0))) ;

                  _measuredQ[hitPlane]=cluster->getTotalCharge();
                }
              else
                _measuredQ[hitPlane]=0.;



              if(debug) {
                streamlog_out ( TESTFITTERMESSAGE )  << "Measured hit in plane " << hitPlane << " at  X = "
                                                     << pos[0] << ", Y = " << pos[1] << ", Q = " << _measuredQ[hitPlane] << endl;
              }

            }
          else
            {
              // Fitted hits

              _isFitted[hitPlane]=true;

              _fittedX[hitPlane]=pos[0];
              _fittedY[hitPlane]=pos[1];

              if(debug) {
                streamlog_out ( TESTFITTERMESSAGE )   << "Fitted  hit  in plane " << hitPlane << " at  X = "
                                                      << pos[0] << ", Y = " << pos[1] << endl;
              }

            }

        }

      // Histograms of measured positions

      for(int ipl=0;ipl<_nTelPlanes; ipl++)
        {
          if(_isMeasured[ipl])
            {
              string tempHistoName;

              stringstream nam;
              nam << _MeasuredXHistoName << "_" << _planeID[ ipl ] ;
              tempHistoName=nam.str();
              (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredX[ipl]);

              stringstream nam2;
              nam2 << _MeasuredYHistoName << "_" << _planeID[ ipl ] ;
              tempHistoName=nam2.str();
              (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredY[ipl]);

              stringstream nam3;
              nam3 << _MeasuredXYHistoName << "_" << _planeID[ ipl ] ;
              tempHistoName=nam3.str();
              (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredX[ipl],_measuredY[ipl]);

              stringstream nam4;
              nam4 << _clusterSignalHistoName << "_" << _planeID[ ipl ] ;
              tempHistoName=nam4.str();
              (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredQ[ipl]);

              stringstream nam5;
              nam5 << _meanSignalXHistoName << "_" << _planeID[ ipl ] ;
              tempHistoName=nam5.str();
              (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredX[ipl],_measuredQ[ipl]);


              stringstream nam6;
              nam6 << _meanSignalYHistoName << "_" << _planeID[ ipl ] ;
              tempHistoName=nam6.str();
              (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredY[ipl],_measuredQ[ipl]);


              stringstream nam7;
              nam7 << _meanSignalXYHistoName << "_" << _planeID[ ipl ] ;
              tempHistoName=nam7.str();
              (dynamic_cast<AIDA::IProfile2D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredX[ipl],_measuredY[ipl],_measuredQ[ipl]);

              stringstream nam8;
              nam8 << _ShiftXvsYHistoName << "_" << _planeID[ ipl ] ;
              tempHistoName=nam8.str();
              (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredX[ipl], _measuredY[ipl] - _fittedY[ipl]);

              stringstream nam9;
              nam9 << _ShiftYvsXHistoName << "_" << _planeID[ ipl ] ;
              tempHistoName=nam9.str();
              (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredY[ipl], _measuredX[ipl] - _fittedX[ipl]);
            }
        }


      // Histograms of fitted positions

      for(int ipl=0;ipl<_nTelPlanes; ipl++)
        {
          if(_isFitted[ipl])
            {
              string tempHistoName;

              stringstream nam;
              nam << _FittedXHistoName << "_" << _planeID[ ipl ] ;
              tempHistoName=nam.str();
              AIDA::IHistogram1D* histo = (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]));
              if ( histo ) histo->fill(_fittedX[ipl]);
              else cout << tempHistoName << endl;

              stringstream nam2;
              nam2 << _FittedYHistoName << "_" << _planeID[ ipl ] ;
              tempHistoName=nam2.str();
              histo = (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]));
              if ( histo) histo->fill(_fittedY[ipl]);
              else cout << tempHistoName << endl;

              stringstream nam3;
              nam3 << _FittedXYHistoName << "_" << _planeID[ ipl ] ;
              tempHistoName=nam3.str();
              AIDA::IHistogram2D * histo2D = (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[tempHistoName]));
              if ( histo2D ) histo2D->fill(_fittedX[ipl],_fittedY[ipl]);
              else cout << tempHistoName << endl;

            }
        }


      // Histograms of incident particle angle

      for(int ipl=1;ipl<_nTelPlanes; ipl++)
        {
          if(_isFitted[ipl] && _isFitted[ipl-1])
            {
              string tempHistoName;

              double angleX=(_fittedX[ipl]-_fittedX[ipl-1])/
                (_planePosition[ipl]- _planePosition[ipl-1]);

              double angleY=(_fittedY[ipl]-_fittedY[ipl-1])/
                (_planePosition[ipl]- _planePosition[ipl-1]);

              stringstream nam;
              nam << _AngleXHistoName << "_" << _planeID[ ipl ] ;
              tempHistoName=nam.str();
              (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(angleX);

              stringstream nam2;
              nam2 << _AngleYHistoName << "_" << _planeID[ ipl ] ;
              tempHistoName=nam2.str();
              (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(angleY);

              stringstream nam3;
              nam3 << _AngleXYHistoName << "_" << _planeID[ ipl ] ;
              tempHistoName=nam3.str();
              (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[tempHistoName]))->fill(angleX,angleY);

            }
        }



      // Histograms of particle scattering angle

      for(int ipl=1;ipl<_nTelPlanes-1; ipl++)
        {
          if(_isFitted[ipl] && _isFitted[ipl+1] && _isFitted[ipl-1] )
            {
              string tempHistoName;

              double scatX=(_fittedX[ipl+1]-_fittedX[ipl])/
                (_planePosition[ipl+1]- _planePosition[ipl]);

              if(ipl>0)scatX-=(_fittedX[ipl]-_fittedX[ipl-1])/
                         (_planePosition[ipl]- _planePosition[ipl-1]);

              double scatY=(_fittedY[ipl+1]-_fittedY[ipl])/
                (_planePosition[ipl+1]- _planePosition[ipl]);

              if(ipl>0)scatY-=(_fittedY[ipl]-_fittedY[ipl-1])/
                         (_planePosition[ipl]- _planePosition[ipl-1]);

              stringstream nam;
              nam << _ScatXHistoName << "_" << _planeID[ ipl ] ;
              tempHistoName=nam.str();
              (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(scatX);

              stringstream nam2;
              nam2 << _ScatYHistoName << "_" << _planeID[ ipl ] ;
              tempHistoName=nam2.str();
              (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(scatY);

              stringstream nam3;
              nam3 << _ScatXYHistoName << "_" << _planeID[ ipl ] ;
              tempHistoName=nam3.str();
              (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[tempHistoName]))->fill(scatX,scatY);

            }
        }



      // Histograms of residuals

      for(int ipl=0;ipl<_nTelPlanes; ipl++)
        {
          if(_isMeasured[ipl] && _isFitted[ipl])
            {
              string tempHistoName;

              stringstream nam;
              nam << _ResidualXHistoName << "_" << _planeID[ ipl ] ;
              tempHistoName=nam.str();
              (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(_fittedX[ipl]-_measuredX[ipl]);

              stringstream nam2;
              nam2 << _ResidualYHistoName << "_" << _planeID[ ipl ] ;
              tempHistoName=nam2.str();
              (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(_fittedY[ipl]-_measuredY[ipl]);

              stringstream nam3;
              nam3 << _ResidualXYHistoName << "_" << _planeID[ ipl ] ;
              tempHistoName=nam3.str();
              (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[tempHistoName]))->fill(_fittedX[ipl]-_measuredX[ipl],_fittedY[ipl]-_measuredY[ipl]);

            }
        }



      // Alignment w.r.t the beam direction check histograms; only if requested

      if(_isMeasured[_beamID] && _alignCheckHistograms )
        {
          for(int ipl=0;ipl<_nTelPlanes; ipl++)
            {
              if(ipl!=_beamID && _isMeasured[ipl])
                {
                  string tempHistoName;

                  stringstream nam;
                  nam << _beamShiftXHistoName << "_" << _planeID[ ipl ] ;
                  tempHistoName=nam.str();
                  (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredX[ipl]-_measuredX[_beamID]);

                  stringstream nam1;
                  nam1 << _beamShiftYHistoName << "_" << _planeID[ ipl ] ;
                  tempHistoName=nam1.str();
                  (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredY[ipl]-_measuredY[_beamID]);

                  stringstream nam2;
                  nam2 << _beamShiftXYHistoName << "_" << _planeID[ ipl ] ;
                  tempHistoName=nam2.str();
                  (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredX[ipl]-_measuredX[_beamID],_measuredY[ipl]-_measuredY[_beamID]);

                  stringstream nam3;
                  nam3 << _beamRotXHistoName << "_" << _planeID[ ipl ] ;
                  tempHistoName=nam3.str();
                  (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredY[_beamID],_measuredX[ipl]-_measuredX[_beamID]);

                  stringstream nam4;
                  nam4 << _beamRotYHistoName << "_" << _planeID[ ipl ] ;
                  tempHistoName=nam4.str();
                  (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredX[_beamID],_measuredY[ipl]-_measuredY[_beamID]);

                  stringstream nam5;
                  nam5 << _beamRot2XHistoName << "_" << _planeID[ ipl ] ;
                  tempHistoName=nam5.str();
                  (dynamic_cast<AIDA::IProfile2D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredX[_beamID],_measuredY[_beamID],_measuredX[ipl]-_measuredX[_beamID]);

                  stringstream nam6;
                  nam6 << _beamRot2YHistoName << "_" << _planeID[ ipl ] ;
                  tempHistoName=nam6.str();
                  (dynamic_cast<AIDA::IProfile2D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredX[_beamID],_measuredY[_beamID],_measuredY[ipl]-_measuredY[_beamID]);

                  stringstream nam7;
                  nam7 << _beamRotX2DHistoName << "_" << _planeID[ ipl ] ;
                  tempHistoName=nam7.str();
                  (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredY[_beamID],_measuredX[ipl]-_measuredX[_beamID]);

                  stringstream nam8;
                  nam8 << _beamRotY2DHistoName << "_" << _planeID[ ipl ] ;
                  tempHistoName=nam8.str();
                  (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredX[_beamID],_measuredY[ipl]-_measuredY[_beamID]);

                }
            }
        }


      // Alignment check w.r.t two selected planes; only if requested

      if(_isMeasured[_referenceID0] && _isMeasured[_referenceID1]  && _alignCheckHistograms)
        {
          for(int ipl=0;ipl<_nTelPlanes; ipl++)
            {
              if(ipl!=_referenceID0 && ipl!=_referenceID1 && _isMeasured[ipl])
                {
                  string tempHistoName;

                  double lineX =
                    ( _measuredX[_referenceID0]*(_planePosition[_referenceID1]-_planePosition[ipl])
                      + _measuredX[_referenceID1]*(_planePosition[ipl]-_planePosition[_referenceID0]))/
                    (_planePosition[_referenceID1]- _planePosition[_referenceID0]);

                  double lineY =
                    ( _measuredY[_referenceID0]*(_planePosition[_referenceID1]-_planePosition[ipl])
                      + _measuredY[_referenceID1]*(_planePosition[ipl]-_planePosition[_referenceID0]))/
                    (_planePosition[_referenceID1]- _planePosition[_referenceID0]);

                  stringstream nam;
                  nam << _relShiftXHistoName << "_" << _planeID[ ipl ] ;
                  tempHistoName=nam.str();
                  (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredX[ipl]-lineX);

                  stringstream nam2;
                  nam2 << _relShiftYHistoName << "_" << _planeID[ ipl ] ;
                  tempHistoName=nam2.str();
                  (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredY[ipl]-lineY);

                  stringstream nam3;
                  nam3 << _relRotXHistoName << "_" << _planeID[ ipl ] ;
                  tempHistoName=nam3.str();
                  (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[tempHistoName]))->fill(lineY,_measuredX[ipl]-lineX);

                  stringstream nam4;
                  nam4 << _relRotYHistoName << "_" << _planeID[ ipl ] ;
                  tempHistoName=nam4.str();
                  (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[tempHistoName]))->fill(lineX,_measuredY[ipl]-lineY);
                  stringstream nam5;
                  nam5 << _relRotX2DHistoName << "_" << _planeID[ ipl ] ;
                  tempHistoName=nam5.str();
                  (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[tempHistoName]))->fill(lineY,_measuredX[ipl]-lineX);

                  stringstream nam6;
                  nam6 << _relRotY2DHistoName << "_" << _planeID[ ipl ] ;
                  tempHistoName=nam6.str();
                  (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[tempHistoName]))->fill(lineX,_measuredY[ipl]-lineY);
                }
            }
        }


      // End of loop over tracks
    }


  return;
}



void EUTelFitHistograms::check( LCEvent * /* evt  */ ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void EUTelFitHistograms::end(){

  //   std::cout << "EUTelFitHistograms::end()  " << name()
  //        << " processed " << _nEvt << " events in " << _nRun << " runs "
  //        << std::endl ;


  // Clean memory

  delete [] _planeSort ;
  delete [] _planePosition ;
  delete [] _planeID ;
  delete [] _isActive ;

  delete [] _isMeasured ;
  delete [] _isFitted ;
  delete [] _measuredX  ;
  delete [] _measuredY  ;
  delete [] _measuredQ  ;
  delete [] _fittedX ;
  delete [] _fittedY ;

}



void EUTelFitHistograms::bookHistos()
{

  streamlog_out ( MESSAGE5 ) << "Booking histograms " << endl;


  streamlog_out ( MESSAGE5 ) << "Histogram information searched in " << _histoInfoFileName << endl;

  std::unique_ptr<EUTelHistogramManager> histoMgr = std::make_unique<EUTelHistogramManager>(_histoInfoFileName);
  EUTelHistogramInfo    * histoInfo;
  bool                    isHistoManagerAvailable;

  try {
    isHistoManagerAvailable = histoMgr->init();
  } catch ( ios::failure& e) {
    streamlog_out ( ERROR5 ) << "I/O problem with " << _histoInfoFileName << "\n"
                            << "Continuing without histogram manager"    << endl;
    isHistoManagerAvailable = false;
  } catch ( ParseException& e ) {
    streamlog_out ( ERROR5 ) << e.what() << "\n"
                            << "Continuing without histogram manager" << endl;
    isHistoManagerAvailable = false;
  }



  //
  {
    int    shiftXNBin  = 200;
    double  shiftXMin  = -5.;
    double  shiftXMax  = 5.;

    string shiftXTitle = "Measured - fitted X position vs Y";

    if ( isHistoManagerAvailable )
      {
        histoInfo = histoMgr->getHistogramInfo(_ShiftXvsYHistoName);
        if ( histoInfo )
          {
            streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
            shiftXNBin = histoInfo->_xBin;
            shiftXMin  = histoInfo->_xMin;
            shiftXMax  = histoInfo->_xMax;
            if ( histoInfo->_title != "" ) shiftXTitle = histoInfo->_title;
          }
      }


    string tempHistoTitle;
    string tempHistoName;

    for(int ipl=0;ipl<_nTelPlanes; ipl++)    {
      if(_isActive[ipl]) {
        tempHistoName  = _ShiftXvsYHistoName +  "_"  + to_string( _planeID[ ipl ] );
        tempHistoTitle =  shiftXTitle + " for plane " + to_string( _planeID[ ipl ] );
        AIDA::IProfile1D * tempHisto = AIDAProcessor::histogramFactory(this)->createProfile1D( tempHistoName.c_str(), shiftXNBin,shiftXMin,shiftXMax);
      
        tempHisto->setTitle(tempHistoTitle.c_str());
        _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));
        
      }
    }
  }

  //
  {
    int    shiftYNBin  = 200;
    double  shiftYMin  = -5.;
    double  shiftYMax  = 5.;

    string shiftYTitle = "Measured - fitted Y position vs X";

    if ( isHistoManagerAvailable )
      {
        histoInfo = histoMgr->getHistogramInfo(_ShiftYvsXHistoName);
        if ( histoInfo )
          {
            streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
            shiftYNBin = histoInfo->_xBin;
            shiftYMin  = histoInfo->_xMin;
            shiftYMax  = histoInfo->_xMax;
            if ( histoInfo->_title != "" ) shiftYTitle = histoInfo->_title;
          }
      }


    string tempHistoTitle;
    string tempHistoName;

    for(int ipl=0;ipl<_nTelPlanes; ipl++)    {
      if(_isActive[ipl]) {
        tempHistoName  = _ShiftYvsXHistoName +  "_"  + to_string( _planeID[ ipl ] );
        tempHistoTitle =  shiftYTitle + " for plane " + to_string( _planeID[ ipl ] );
        AIDA::IProfile1D * tempHisto = AIDAProcessor::histogramFactory(this)->createProfile1D( tempHistoName.c_str(), shiftYNBin,shiftYMin,shiftYMax);
        
        tempHisto->setTitle(tempHistoTitle.c_str());
        _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));
        
      }
    }
  }





  // Measured position in X

  int    measXNBin  = 400;
  double measXMin   = -5.;
  double measXMax   = 5.;
  string measXTitle = "Measured particle position in X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_MeasuredXHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          measXNBin = histoInfo->_xBin;
          measXMin  = histoInfo->_xMin;
          measXMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) measXTitle = histoInfo->_title;
        }
    }


  string tempHistoTitle;
  string tempHistoName;

  for(int ipl=0;ipl<_nTelPlanes; ipl++)    {
    if(_isActive[ipl]) {
      tempHistoName  = _MeasuredXHistoName +  "_"  + to_string( _planeID[ ipl ] );
      tempHistoTitle =  measXTitle + " for plane " + to_string( _planeID[ ipl ] );
      AIDA::IHistogram1D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( tempHistoName.c_str(),measXNBin,measXMin,measXMax);
      tempHisto->setTitle(tempHistoTitle.c_str());
      _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));
    }

  }


  // Measured position in Y

  int    measYNBin  = 400;
  double measYMin   = -5.;
  double measYMax   = 5.;
  string measYTitle = "Measured particle position in Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_MeasuredYHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          measYNBin = histoInfo->_xBin;
          measYMin  = histoInfo->_xMin;
          measYMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) measYTitle = histoInfo->_title;
        }
    }


  for(int ipl=0;ipl<_nTelPlanes; ipl++) {
    if(_isActive[ipl]) {
      tempHistoName   =  _MeasuredYHistoName + "_" + to_string( _planeID[ ipl ] ) ;
      tempHistoTitle  =  measYTitle + " for plane " + to_string( _planeID[ ipl ] ) ;
      AIDA::IHistogram1D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( tempHistoName.c_str(),measYNBin,measYMin,measYMax);
      tempHisto->setTitle(tempHistoTitle.c_str());
      _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));
    }
  }

  // Measured position in X-Y
  measXNBin  = 100;
  measXMin   = -5.;
  measXMax   = 5.;
  measYNBin  = 100;
  measYMin   = -5.;
  measYMax   = 5.;
  string measXYTitle = "Measured particle position in XY";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_MeasuredXYHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          measXNBin = histoInfo->_xBin;
          measXMin  = histoInfo->_xMin;
          measXMax  = histoInfo->_xMax;
          measYNBin = histoInfo->_yBin;
          measYMin  = histoInfo->_yMin;
          measYMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) measXYTitle = histoInfo->_title;
        }
    }


  for(int ipl=0;ipl<_nTelPlanes; ipl++) {
    if(_isActive[ipl])   {
      tempHistoName   =  _MeasuredXYHistoName + "_" +  to_string( _planeID[ ipl ] ) ; ;
      tempHistoTitle  =   measXYTitle + " for plane " + to_string( _planeID[ ipl ] ) ;
      AIDA::IHistogram2D * tempHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(),measXNBin,measXMin,measXMax,measYNBin,measYMin,measYMax);
      tempHisto->setTitle(tempHistoTitle.c_str());
      _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));
    }
  }





  // Fitted position in X

  int    fitXNBin  = 400;
  double fitXMin   = -5.;
  double fitXMax   = 5.;
  string fitXTitle = "Fitted particle position in X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_FittedXHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          fitXNBin = histoInfo->_xBin;
          fitXMin  = histoInfo->_xMin;
          fitXMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) fitXTitle = histoInfo->_title;
        }
    }


  for(int ipl=0;ipl<_nTelPlanes; ipl++)    {
    if(_isActive[ipl])        {
      tempHistoName  = _FittedXHistoName + "_" +  to_string( _planeID[ ipl ] ) ;
      tempHistoTitle = fitXTitle + " for plane "  + to_string( _planeID[ ipl ] ) ;
      AIDA::IHistogram1D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( tempHistoName.c_str(),fitXNBin,fitXMin,fitXMax);
      tempHisto->setTitle(tempHistoTitle.c_str());
      _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));
    }
  }


  // Fitted position in Y

  int    fitYNBin  = 400;
  double fitYMin   = -5.;
  double fitYMax   = 5.;
  string fitYTitle = "Fitted particle position in Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_FittedYHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          fitYNBin = histoInfo->_xBin;
          fitYMin  = histoInfo->_xMin;
          fitYMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) fitYTitle = histoInfo->_title;
        }
    }


  for(int ipl=0;ipl<_nTelPlanes; ipl++){
    if(_isActive[ipl])        {
      tempHistoName  = _FittedYHistoName + "_" + to_string( _planeID[ ipl ] ) ;
      tempHistoTitle = fitYTitle + " for plane " + to_string( _planeID[ ipl ] ) ;
      AIDA::IHistogram1D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( tempHistoName.c_str(),fitYNBin,fitYMin,fitYMax);
      tempHisto->setTitle(tempHistoTitle.c_str());
      _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));
    }
  }



  // Fitted position in X-Y

  fitXNBin  = 100;
  fitXMin   = -5.;
  fitXMax   = 5.;
  fitYNBin  = 100;
  fitYMin   = -5.;
  fitYMax   = 5.;
  string fitXYTitle = "Fitted particle position in XY";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_FittedXYHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          fitXNBin = histoInfo->_xBin;
          fitXMin  = histoInfo->_xMin;
          fitXMax  = histoInfo->_xMax;
          fitYNBin = histoInfo->_yBin;
          fitYMin  = histoInfo->_yMin;
          fitYMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) fitXYTitle = histoInfo->_title;
        }
    }


  for(int ipl=0;ipl<_nTelPlanes; ipl++)    {
    if(_isActive[ipl])        {
      tempHistoName  = _FittedXYHistoName + "_" + to_string( _planeID[ ipl ] ) ;
      tempHistoTitle = fitXYTitle + " for plane " + to_string( _planeID[ ipl ] ) ;
      AIDA::IHistogram2D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(),fitXNBin,fitXMin,fitXMax,fitYNBin,fitYMin,fitYMax);
      tempHisto->setTitle(tempHistoTitle.c_str());
      _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));
    }
  }

  // Incident X angle

  int    angleXNBin  = 400;
  double angleXMin   = -0.04;
  double angleXMax   = 0.04;
  string angleXTitle = "Incident particle angle in X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_AngleXHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          angleXNBin = histoInfo->_xBin;
          angleXMin  = histoInfo->_xMin;
          angleXMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) angleXTitle = histoInfo->_title;
        }
    }


  for(int ipl=1;ipl<_nTelPlanes; ipl++)  {
    if(_isActive[ipl])  {
      tempHistoName  = _AngleXHistoName + "_" + to_string( _planeID[ ipl ] ) ;
      tempHistoTitle = angleXTitle + " for plane " + to_string( _planeID[ ipl ] ) ;
      AIDA::IHistogram1D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( tempHistoName.c_str(),angleXNBin,angleXMin,angleXMax);
      tempHisto->setTitle(tempHistoTitle.c_str());
      _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));
    }
  }


  // Incident angle in Y

  int    angleYNBin  = 400;
  double angleYMin   = -0.04;
  double angleYMax   = 0.04;
  string angleYTitle = "Incident particle angle in Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_AngleYHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          angleYNBin = histoInfo->_xBin;
          angleYMin  = histoInfo->_xMin;
          angleYMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) angleYTitle = histoInfo->_title;
        }
    }


  for(int ipl=1;ipl<_nTelPlanes; ipl++){
    if(_isActive[ipl])        {
      tempHistoName  = _AngleYHistoName + "_" + to_string( _planeID[ ipl ] ) ;
      tempHistoTitle = angleYTitle + " for plane " + to_string( _planeID[ ipl ] ) ;
      AIDA::IHistogram1D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( tempHistoName.c_str(),angleYNBin,angleYMin,angleYMax);
      tempHisto->setTitle(tempHistoTitle.c_str());
      _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));
    }
  }



  // Incident angle in X-Y

  angleXNBin  = 100;
  angleXMin   = -0.1;
  angleXMax   = 0.1;
  angleYNBin  = 100;
  angleYMin   = -0.1;
  angleYMax   = 0.1;
  string angleXYTitle = "Incident particle angle in XY";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_AngleXYHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          angleXNBin = histoInfo->_xBin;
          angleXMin  = histoInfo->_xMin;
          angleXMax  = histoInfo->_xMax;
          angleYNBin = histoInfo->_yBin;
          angleYMin  = histoInfo->_yMin;
          angleYMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) angleXYTitle = histoInfo->_title;
        }
    }


  for(int ipl=1;ipl<_nTelPlanes; ipl++)  {
    if(_isActive[ipl])   {
      tempHistoName  = _AngleXYHistoName + "_" + to_string( _planeID[ ipl ] ) ;
      tempHistoTitle = angleXYTitle + " for plane " + to_string( _planeID[ ipl ] ) ;
      AIDA::IHistogram2D * tempHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(),angleXNBin,angleXMin,angleXMax,angleYNBin,angleYMin,angleYMax);
      tempHisto->setTitle(tempHistoTitle.c_str());
      _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));
    }
  }


  // Scattering angle in X

  int    scatXNBin  = 400;
  double scatXMin   = -0.01;
  double scatXMax   = 0.01;
  string scatXTitle = "Particle scattering angle in X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ScatXHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          scatXNBin = histoInfo->_xBin;
          scatXMin  = histoInfo->_xMin;
          scatXMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) scatXTitle = histoInfo->_title;
        }
    }


  for(int ipl=1;ipl<_nTelPlanes-1; ipl++) {
    if(_isActive[ipl])        {
      tempHistoName  = _ScatXHistoName + "_" + to_string( _planeID[ ipl ] ) ;
      tempHistoTitle = scatXTitle + " for plane " + to_string( _planeID[ ipl ] ) ;
      AIDA::IHistogram1D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( tempHistoName.c_str(),scatXNBin,scatXMin,scatXMax);
      tempHisto->setTitle(tempHistoTitle.c_str());
      _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));
    }
  }


  // Scattering angle in Y

  int    scatYNBin  = 400;
  double scatYMin   = -0.01;
  double scatYMax   = 0.01;
  string scatYTitle = "Particle scattering angle in Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ScatYHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          scatYNBin = histoInfo->_xBin;
          scatYMin  = histoInfo->_xMin;
          scatYMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) scatYTitle = histoInfo->_title;
        }
    }


  for(int ipl=1;ipl<_nTelPlanes-1; ipl++) {
    if(_isActive[ipl]) {
      tempHistoName  = _ScatYHistoName + "_" + to_string( _planeID[ ipl ] ) ;
      tempHistoTitle = scatYTitle + " for plane " + to_string( _planeID[ ipl ] ) ;
      AIDA::IHistogram1D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( tempHistoName.c_str(),scatYNBin,scatYMin,scatYMax);
      tempHisto->setTitle(tempHistoTitle.c_str());
      _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));
    }
  }



  // Scattering angle in X-Y

  scatXNBin  = 100;
  scatXMin   = -0.01;
  scatXMax   = 0.01;
  scatYNBin  = 100;
  scatYMin   = -0.01;
  scatYMax   = 0.01;
  string scatXYTitle = "Particle scattering angle in XY";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ScatXYHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          scatXNBin = histoInfo->_xBin;
          scatXMin  = histoInfo->_xMin;
          scatXMax  = histoInfo->_xMax;
          scatYNBin = histoInfo->_yBin;
          scatYMin  = histoInfo->_yMin;
          scatYMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) scatXYTitle = histoInfo->_title;
        }
    }


  for(int ipl=1;ipl<_nTelPlanes-1; ipl++) {
      if(_isActive[ipl])         {
      tempHistoName  = _ScatXYHistoName + "_" + to_string( _planeID[ ipl ] ) ;
      tempHistoTitle = scatXYTitle + " for plane " + to_string( _planeID[ ipl ] ) ;
          AIDA::IHistogram2D * tempHisto =
            AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(),scatXNBin,scatXMin,scatXMax,scatYNBin,scatYMin,scatYMax);
          tempHisto->setTitle(tempHistoTitle.c_str());
          _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));
        }
    }



  // Residual position in X

  int    residXNBin  = 400;
  double residXMin   = -0.1;
  double residXMax   = 0.1;
  string residXTitle = "Residual in X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ResidualXHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          residXNBin = histoInfo->_xBin;
          residXMin  = histoInfo->_xMin;
          residXMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) residXTitle = histoInfo->_title;
        }
    }


  for(int ipl=0;ipl<_nTelPlanes; ipl++){
      if(_isActive[ipl]){
      tempHistoName  = _ResidualXHistoName + "_" + to_string( _planeID[ ipl ] ) ;
      tempHistoTitle = residXTitle + " for plane " + to_string( _planeID[ ipl ] ) ;
          AIDA::IHistogram1D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( tempHistoName.c_str(),residXNBin,residXMin,residXMax);
          tempHisto->setTitle(tempHistoTitle.c_str());
          _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));
        }
    }


  // Residual position in Y

  int    residYNBin  = 400;
  double residYMin   = -0.1;
  double residYMax   = 0.1;
  string residYTitle = "Residual in Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ResidualYHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          residYNBin = histoInfo->_xBin;
          residYMin  = histoInfo->_xMin;
          residYMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) residYTitle = histoInfo->_title;
        }
    }


  for(int ipl=0;ipl<_nTelPlanes; ipl++)  {
      if(_isActive[ipl])  {
      tempHistoName  = _ResidualYHistoName + "_" + to_string( _planeID[ ipl ] ) ;
      tempHistoTitle = residYTitle + " for plane " + to_string( _planeID[ ipl ] ) ;
          AIDA::IHistogram1D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( tempHistoName.c_str(),residYNBin,residYMin,residYMax);
          tempHisto->setTitle(tempHistoTitle.c_str());
          _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));
        }
    }



  // Residual position in X-Y

  residXNBin  = 100;
  residXMin   = -0.1;
  residXMax   = 0.1;
  residYNBin  = 100;
  residYMin   = -0.1;
  residYMax   = 0.1;
  string residXYTitle = "Residual in XY";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_ResidualXYHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          residXNBin = histoInfo->_xBin;
          residXMin  = histoInfo->_xMin;
          residXMax  = histoInfo->_xMax;
          residYNBin = histoInfo->_yBin;
          residYMin  = histoInfo->_yMin;
          residYMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) residXYTitle = histoInfo->_title;
        }
    }


  for(int ipl=0;ipl<_nTelPlanes; ipl++){
      if(_isActive[ipl]){

      tempHistoName  = _ResidualXYHistoName + "_" + to_string( _planeID[ ipl ] ) ;
      tempHistoTitle = residXYTitle + " for plane " + to_string( _planeID[ ipl ] ) ;
          AIDA::IHistogram2D * tempHisto = 
AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(),residXNBin,residXMin,residXMax,residYNBin,residYMin,residYMax);
          tempHisto->setTitle(tempHistoTitle.c_str());
          _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));
        }
    }





  //
  // Cluster signal distribution
  //

  int    clusterNBin  = 1000;
  double clusterMin   = 0.;
  double clusterMax   = 1000.;
  string clusterTitle = "Cluster spectrum with all pixels";

  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_clusterSignalHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          clusterNBin = histoInfo->_xBin;
          clusterMin  = histoInfo->_xMin;
          clusterMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) clusterTitle = histoInfo->_title;
        }
    }


  for(int ipl=0;ipl<_nTelPlanes; ipl++)
    {
      if(_isActive[ipl])
        {
          stringstream nam,tit;

          nam << _clusterSignalHistoName << "_" << _planeID[ ipl ] ;
          tempHistoName=nam.str();

          tit << clusterTitle << " for plane " << _planeID[ ipl ] ;
          tempHistoTitle=tit.str();

          AIDA::IHistogram1D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( tempHistoName.c_str(), clusterNBin,clusterMin,clusterMax);

          tempHisto->setTitle(tempHistoTitle.c_str());

          _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));

        }

    }

  //
  // Mean cluster signal vs X
  //

  int    meanXNBin  = 100;
  double meanXMin   = -5.;
  double meanXMax   = 5.;
  string meanXTitle = "Mean cluster signal vs X ";

  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_meanSignalXHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          meanXNBin = histoInfo->_xBin;
          meanXMin  = histoInfo->_xMin;
          meanXMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) meanXTitle = histoInfo->_title;
        }
    }


  for(int ipl=0;ipl<_nTelPlanes; ipl++)
    {
      if(_isActive[ipl])
        {
          stringstream nam,tit;

          nam << _meanSignalXHistoName << "_" << _planeID[ ipl ] ;
          tempHistoName=nam.str();

          tit << meanXTitle << " for plane " << _planeID[ ipl ] ;
          tempHistoTitle=tit.str();

          AIDA::IProfile1D * tempHisto = AIDAProcessor::histogramFactory(this)->createProfile1D( tempHistoName.c_str(), meanXNBin,meanXMin,meanXMax);

          tempHisto->setTitle(tempHistoTitle.c_str());

          _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));

        }

    }
  //
  // Mean cluster signal vs Y
  //

  int    meanYNBin  = 100;
  double meanYMin   = -5.;
  double meanYMax   = 5.;
  string meanYTitle = "Mean cluster signal vs Y ";

  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_meanSignalYHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          meanYNBin = histoInfo->_xBin;
          meanYMin  = histoInfo->_xMin;
          meanYMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) meanYTitle = histoInfo->_title;
        }
    }


  for(int ipl=0;ipl<_nTelPlanes; ipl++)
    {
      if(_isActive[ipl])
        {
          stringstream nam,tit;

          nam << _meanSignalYHistoName << "_" << _planeID[ ipl ] ;
          tempHistoName=nam.str();

          tit << meanYTitle << " for plane " << _planeID[ ipl ] ;
          tempHistoTitle=tit.str();

          AIDA::IProfile1D * tempHisto = AIDAProcessor::histogramFactory(this)->createProfile1D( tempHistoName.c_str(), meanYNBin,meanYMin,meanYMax);

          tempHisto->setTitle(tempHistoTitle.c_str());

          _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));

        }

    }

  //
  // Mean cluster signal vs XY
  //

  meanXNBin  = 100;
  meanXMin   = -5.;
  meanXMax   = 5.;
  meanYNBin  = 100;
  meanYMin   = -5.;
  meanYMax   = 5.;
  string meanXYTitle = "Mean cluster signal vs XY ";

  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_meanSignalXYHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          meanXNBin = histoInfo->_xBin;
          meanXMin  = histoInfo->_xMin;
          meanXMax  = histoInfo->_xMax;
          meanYNBin = histoInfo->_yBin;
          meanYMin  = histoInfo->_yMin;
          meanYMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) meanXYTitle = histoInfo->_title;
        }
    }


  for(int ipl=0;ipl<_nTelPlanes; ipl++)
    {
      if(_isActive[ipl])
        {
          stringstream nam,tit;

          nam << _meanSignalXYHistoName << "_" << _planeID[ ipl ] ;
          tempHistoName=nam.str();

          tit << meanXYTitle << " for plane " << _planeID[ ipl ] ;
          tempHistoTitle=tit.str();

          AIDA::IProfile2D * tempHisto = AIDAProcessor::histogramFactory(this)->createProfile2D( tempHistoName.c_str(), meanXNBin,meanXMin,meanXMax,  meanYNBin,meanYMin,meanYMax);

          tempHisto->setTitle(tempHistoTitle.c_str());

          _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));

        }

    }


  // 
  // Book alignment check histograms only if requested
  // =================================================

  if(_alignCheckHistograms){


  // Beam alignment histograms: shift in X

  int    shiftXNBin  = 400;
  double shiftXMin   = -2.;
  double shiftXMax   = 2.;
  string shiftXTitle = "Beam alignment shift in X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_beamShiftXHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          shiftXNBin = histoInfo->_xBin;
          shiftXMin  = histoInfo->_xMin;
          shiftXMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) shiftXTitle = histoInfo->_title;
        }
    }


  for(int ipl=0;ipl<_nTelPlanes; ipl++)
  {
      if(_isActive[ipl] && ipl!=_beamID)        {

        tempHistoName  = _beamShiftXHistoName + "_" + to_string( _planeID[ ipl ] ) ;
        tempHistoTitle = shiftXTitle + " for plane " + to_string( _planeID[ ipl ] ) + " w.r.t. " + to_string( _planeID[ _beamID ] ) ;
        AIDA::IHistogram1D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( tempHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);

        tempHisto->setTitle(tempHistoTitle.c_str());

          _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));

        }

    }


  // Beam alignment histograms: shift in Y

  int    shiftYNBin  = 400;
  double shiftYMin   = -2.;
  double shiftYMax   = 2.;
  string shiftYTitle = "Beam alignment shift in Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_beamShiftYHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          shiftYNBin = histoInfo->_xBin;
          shiftYMin  = histoInfo->_xMin;
          shiftYMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) shiftYTitle = histoInfo->_title;
        }
    }


  for(int ipl=0;ipl<_nTelPlanes; ipl++)
    {
      if(_isActive[ipl] && ipl!=_beamID)
        {
          stringstream nam,tit;

          nam << _beamShiftYHistoName << "_" << _planeID[ ipl ] ;
          tempHistoName=nam.str();

          tit << shiftYTitle <<  " for plane " << _planeID[ ipl ] << " w.r.t. plane " << _planeID[ _beamID ];
          tempHistoTitle=tit.str();

          AIDA::IHistogram1D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( tempHistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);

          tempHisto->setTitle(tempHistoTitle.c_str());

          _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));

        }

    }


  // Beam alignment: shift in X-Y

  shiftXNBin  = 100;
  shiftXMin   = -2.;
  shiftXMax   = 2.;
  shiftYNBin  = 100;
  shiftYMin   = -2.;
  shiftYMax   = 2.;
  string shiftXYTitle = "Beam alignment shift in XY";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_beamShiftXYHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          shiftXNBin = histoInfo->_xBin;
          shiftXMin  = histoInfo->_xMin;
          shiftXMax  = histoInfo->_xMax;
          shiftYNBin = histoInfo->_yBin;
          shiftYMin  = histoInfo->_yMin;
          shiftYMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) shiftXYTitle = histoInfo->_title;
        }
    }


  for(int ipl=0;ipl<_nTelPlanes; ipl++)
    {
      if(_isActive[ipl] && ipl!=_beamID)
        {
          stringstream nam,tit;

          nam << _beamShiftXYHistoName << "_" << _planeID[ ipl ] ;
          tempHistoName=nam.str();

          tit << shiftXYTitle <<  " for plane " << _planeID[ ipl ] << " w.r.t. plane " << _planeID[_beamID] ;
          tempHistoTitle=tit.str();

          AIDA::IHistogram2D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax,shiftYNBin,shiftYMin,shiftYMax);

          tempHisto->setTitle(tempHistoTitle.c_str());

          _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));

        }

    }



  // Beam alignment histograms: X vs Y to extract rotation in Z

  int    rotXNBin  = 100;
  double rotXMin   = -5.;
  double rotXMax   = 5.;
  double rotVMin   = -0.5;
  double rotVMax   = 0.5;
  string rotXTitle = "Beam alignment shift in X vs Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_beamRotXHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          rotXNBin = histoInfo->_xBin;
          rotXMin  = histoInfo->_xMin;
          rotXMax  = histoInfo->_xMax;
          rotVMin  = histoInfo->_yMin;
          rotVMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) rotXTitle = histoInfo->_title;
        }
    }

  for(int ipl=0;ipl<_nTelPlanes; ipl++)
    {
      if(_isActive[ipl] && ipl!=_beamID)
        {
          stringstream nam,tit;

          nam << _beamRotXHistoName << "_" << _planeID[ ipl ] ;
          tempHistoName=nam.str();

          tit << rotXTitle << " for plane " << _planeID[ ipl ] << " w.r.t. plane " << _planeID[ _beamID ] ;
          tempHistoTitle=tit.str();

          AIDA::IProfile1D * tempHisto = AIDAProcessor::histogramFactory(this)->createProfile1D( tempHistoName.c_str(),rotXNBin,rotXMin,rotXMax,rotVMin,rotVMax);

          tempHisto->setTitle(tempHistoTitle.c_str());

          _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));

        }

    }


  // Beam alignment histograms: Y vs X to extract rotation in Z

  int    rotYNBin  = 100;
  double rotYMin   = -5.;
  double rotYMax   = 5.;
  rotVMin   = -0.5;
  rotVMax   = 0.5;
  string rotYTitle = "Beam alignment shift in Y vs X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_beamRotYHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          rotYNBin = histoInfo->_xBin;
          rotYMin  = histoInfo->_xMin;
          rotYMax  = histoInfo->_xMax;
          rotVMin  = histoInfo->_yMin;
          rotVMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) rotYTitle = histoInfo->_title;
        }
    }

  for(int ipl=0;ipl<_nTelPlanes; ipl++)
    {
      if(_isActive[ipl] && ipl!=_beamID)
        {
          stringstream nam,tit;

          nam << _beamRotYHistoName << "_" << _planeID[ ipl ] ;
          tempHistoName=nam.str();

          tit << rotYTitle << " for plane " << _planeID[ ipl ] << " w.r.t. plane " << _planeID[ _beamID ];
          tempHistoTitle=tit.str();

          AIDA::IProfile1D * tempHisto = AIDAProcessor::histogramFactory(this)->createProfile1D( tempHistoName.c_str(),rotYNBin,rotYMin,rotYMax,rotVMin,rotVMax);

          tempHisto->setTitle(tempHistoTitle.c_str());

          _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));

        }

    }

  // Beam alignment histograms: X vs Y (2D plot)

  rotXNBin  = 100;
  rotXMin   = -5.;
  rotXMax   = 5.;
  int rotVNBin  = 100;
  rotVMin   = -0.5;
  rotVMax   = 0.5;
  rotXTitle = "Beam alignment shift in X vs Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_beamRotX2DHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          rotXNBin = histoInfo->_xBin;
          rotXMin  = histoInfo->_xMin;
          rotXMax  = histoInfo->_xMax;
          rotVNBin = histoInfo->_yBin;
          rotVMin  = histoInfo->_yMin;
          rotVMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) rotXTitle = histoInfo->_title;
        }
    }

  for(int ipl=0;ipl<_nTelPlanes; ipl++)
    {
      if(_isActive[ipl] && ipl!=_beamID)
        {
          stringstream nam,tit;

          nam << _beamRotX2DHistoName << "_" << _planeID[ ipl ] ;
          tempHistoName=nam.str();

          tit << rotXTitle << " for plane " << _planeID[ ipl ] << " w.r.t. plane " << _planeID[ _beamID ];
          tempHistoTitle=tit.str();

          AIDA::IHistogram2D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(),rotXNBin,rotXMin,rotXMax,rotVNBin,rotVMin,rotVMax);

          tempHisto->setTitle(tempHistoTitle.c_str());

          _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));

        }

    }


  // Beam alignment histograms: Y vs X (2D plot)

  rotYNBin  = 100;
  rotYMin   = -5.;
  rotYMax   = 5.;
  rotVNBin  = 100;
  rotVMin   = -0.5;
  rotVMax   = 0.5;
  rotYTitle = "Beam alignment shift in Y vs X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_beamRotY2DHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          rotYNBin = histoInfo->_xBin;
          rotYMin  = histoInfo->_xMin;
          rotYMax  = histoInfo->_xMax;
          rotVNBin = histoInfo->_yBin;
          rotVMin  = histoInfo->_yMin;
          rotVMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) rotYTitle = histoInfo->_title;
        }
    }

  for(int ipl=0;ipl<_nTelPlanes; ipl++)
    {
      if(_isActive[ipl] && ipl!=_beamID)
        {
          stringstream nam,tit;

          nam << _beamRotY2DHistoName << "_" << _planeID[ ipl ] ;
          tempHistoName=nam.str();

          tit << rotYTitle << " for plane " << _planeID[ ipl ] << " w.r.t. plane " << _planeID[ _beamID ];
          tempHistoTitle=tit.str();

          AIDA::IHistogram2D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(),rotYNBin,rotYMin,rotYMax,rotVNBin,rotVMin,rotVMax);

          tempHisto->setTitle(tempHistoTitle.c_str());

          _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));

        }

    }

  // Beam alignment histograms: X vs XY

  rotXNBin  = 50;
  rotXMin   = -5.;
  rotXMax   = 5.;
  rotYNBin  = 50;
  rotYMin   = -5.;
  rotYMax   = 5.;
  rotVMin   = -0.5;
  rotVMax   = 0.5;
  rotXTitle = "Beam alignment shift in X vs XY";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_beamRot2XHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          rotXNBin = histoInfo->_xBin;
          rotXMin  = histoInfo->_xMin;
          rotXMax  = histoInfo->_xMax;
          rotYNBin = histoInfo->_yBin;
          rotYMin  = histoInfo->_yMin;
          rotYMax  = histoInfo->_yMax;
          rotVMin  = histoInfo->_zMin;
          rotVMax  = histoInfo->_zMax;
          if ( histoInfo->_title != "" ) rotXTitle = histoInfo->_title;
        }
    }

  for(int ipl=0;ipl<_nTelPlanes; ipl++)
    {
      if(_isActive[ipl] && ipl!=_beamID)
        {
          stringstream nam,tit;

          nam << _beamRot2XHistoName << "_" << _planeID[ ipl ] ;
          tempHistoName=nam.str();

          tit << rotXTitle << " for plane " << _planeID[ ipl ] << " w.r.t. plane " << _planeID[ _beamID ];
          tempHistoTitle=tit.str();

          AIDA::IProfile2D * tempHisto = AIDAProcessor::histogramFactory(this)->createProfile2D( tempHistoName.c_str(),rotXNBin,rotXMin,rotXMax,rotYNBin,rotYMin,rotYMax,rotVMin,rotVMax);

          tempHisto->setTitle(tempHistoTitle.c_str());

          _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));

        }

    }




  // Beam alignment histograms: Y vs XY

  rotXNBin  = 50;
  rotXMin   = -5.;
  rotXMax   = 5.;
  rotYNBin  = 50;
  rotYMin   = -5.;
  rotYMax   = 5.;
  rotVMin   = -0.5;
  rotVMax   = 0.5;
  rotYTitle = "Beam alignment shift in Y vs XY";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_beamRot2YHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          rotXNBin = histoInfo->_xBin;
          rotXMin  = histoInfo->_xMin;
          rotXMax  = histoInfo->_xMax;
          rotYNBin = histoInfo->_yBin;
          rotYMin  = histoInfo->_yMin;
          rotYMax  = histoInfo->_yMax;
          rotVMin  = histoInfo->_zMin;
          rotVMax  = histoInfo->_zMax;
          if ( histoInfo->_title != "" ) rotYTitle = histoInfo->_title;
        }
    }

  for(int ipl=0;ipl<_nTelPlanes; ipl++)
    {
      if(_isActive[ipl] && ipl!=_beamID)
        {
          stringstream nam,tit;

          nam << _beamRot2YHistoName << "_" << _planeID[ ipl ] ;
          tempHistoName=nam.str();

          tit << rotYTitle << " for plane " << _planeID[ ipl ] << " w.r.t. plane " << _planeID[ _beamID ] ;
          tempHistoTitle=tit.str();

          AIDA::IProfile2D * tempHisto = AIDAProcessor::histogramFactory(this)->createProfile2D( tempHistoName.c_str(),rotXNBin,rotXMin,rotXMax,rotYNBin,rotYMin,rotYMax,rotVMin,rotVMax);

          tempHisto->setTitle(tempHistoTitle.c_str());

          _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));

        }

    }



  // Relative alignment histograms: shift in X

  shiftXNBin  = 400;
  shiftXMin   = -0.2;
  shiftXMax   = 0.2;
  shiftXTitle = "Relative alignment shift in X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_relShiftXHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          shiftXNBin = histoInfo->_xBin;
          shiftXMin  = histoInfo->_xMin;
          shiftXMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) shiftXTitle = histoInfo->_title;
        }
    }


  for(int ipl=0;ipl<_nTelPlanes; ipl++)
    {
      if(_isActive[ipl] && ipl!=_referenceID0 && ipl!=_referenceID1)
        {
          stringstream nam,tit;

          nam << _relShiftXHistoName << "_" << _planeID[ ipl ] ;
          tempHistoName=nam.str();

          tit << shiftXTitle << " for plane " << _planeID[ ipl ] << " w.r.t. planes "
              << _planeID[ _referenceID0 ] << " and " << _planeID[ _referenceID1 ] ;
          tempHistoTitle=tit.str();

          AIDA::IHistogram1D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( tempHistoName.c_str(),shiftXNBin,shiftXMin,shiftXMax);

          tempHisto->setTitle(tempHistoTitle.c_str());

          _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));

        }

    }


  // Relative alignment histograms: shift in Y

  shiftYNBin  = 400;
  shiftYMin   = -0.2;
  shiftYMax   = 0.2;
  shiftYTitle = "Relative alignment shift in Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_relShiftYHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          shiftYNBin = histoInfo->_xBin;
          shiftYMin  = histoInfo->_xMin;
          shiftYMax  = histoInfo->_xMax;
          if ( histoInfo->_title != "" ) shiftYTitle = histoInfo->_title;
        }
    }


  for(int ipl=0;ipl<_nTelPlanes; ipl++)
    {
      if(_isActive[ipl] && ipl!=_referenceID0 && ipl!=_referenceID1)
        {
          stringstream nam,tit;

          nam << _relShiftYHistoName << "_" << _planeID[ ipl ] ;
          tempHistoName=nam.str();

          tit << shiftYTitle <<  " for plane " << _planeID[ ipl ] << " w.r.t. planes "
              << _planeID[ _referenceID0 ] << " and " << _planeID[ _referenceID1 ];
          tempHistoTitle=tit.str();

          AIDA::IHistogram1D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( tempHistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);

          tempHisto->setTitle(tempHistoTitle.c_str());

          _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));

        }

    }


  // Relative alignment histograms: X vs Y to extract rotation in Z

  rotXNBin  = 200;
  rotXMin   = -5.;
  rotXMax   = 5.;
  rotVMin   = -0.1;
  rotVMax   = 0.1;
  rotXTitle = "Relative alignment shift in X vs Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_relRotXHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          rotXNBin = histoInfo->_xBin;
          rotXMin  = histoInfo->_xMin;
          rotXMax  = histoInfo->_xMax;
          rotVMin  = histoInfo->_yMin;
          rotVMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) rotXTitle = histoInfo->_title;
        }
    }

  for(int ipl=0;ipl<_nTelPlanes; ipl++)
    {
      if(_isActive[ipl] && ipl!=_referenceID0 && ipl!=_referenceID1)
        {
          stringstream nam,tit;

          nam << _relRotXHistoName << "_" << _planeID[ ipl ] ;
          tempHistoName=nam.str();

          tit << rotXTitle << " for plane " << _planeID[ ipl ] << " w.r.t. planes "
              << _planeID[ _referenceID0 ] << " and " << _planeID[ _referenceID1 ] ;
          tempHistoTitle=tit.str();

          AIDA::IProfile1D * tempHisto = AIDAProcessor::histogramFactory(this)->createProfile1D( tempHistoName.c_str(),rotXNBin,rotXMin,rotXMax,rotVMin,rotVMax);

          tempHisto->setTitle(tempHistoTitle.c_str());

          _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));

        }

    }


  // Relative alignment histograms: Y vs X to extract rotation in Z

  rotYNBin  = 200;
  rotYMin   = -5.;
  rotYMax   = 5.;
  rotVMin   = -0.1;
  rotVMax   = 0.1;
  rotYTitle = "Relative alignment shift in Y vs X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_relRotYHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          rotYNBin = histoInfo->_xBin;
          rotYMin  = histoInfo->_xMin;
          rotYMax  = histoInfo->_xMax;
          rotVMin  = histoInfo->_yMin;
          rotVMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) rotYTitle = histoInfo->_title;
        }
    }

  for(int ipl=0;ipl<_nTelPlanes; ipl++)
    {
      if(_isActive[ipl] && ipl!=_referenceID0 && ipl!=_referenceID1)
        {
          stringstream nam,tit;

          nam << _relRotYHistoName << "_" << _planeID[ ipl ] ;
          tempHistoName=nam.str();

          tit << rotYTitle << " for plane " << _planeID[ ipl ] << " w.r.t. planes "
              << _planeID[ _referenceID0 ] << " and " << _planeID[ _referenceID1 ] ;
          tempHistoTitle=tit.str();

          AIDA::IProfile1D * tempHisto = AIDAProcessor::histogramFactory(this)->createProfile1D( tempHistoName.c_str(),rotYNBin,rotYMin,rotYMax,rotVMin,rotVMax);

          tempHisto->setTitle(tempHistoTitle.c_str());

          _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));

        }

    }


  // Relative alignment histograms: X vs Y (2D plot)

  rotXNBin  = 100;
  rotXMin   = -5.;
  rotXMax   = 5.;
  rotVNBin  = 100;
  rotVMin   = -0.2;
  rotVMax   = 0.2;
  rotXTitle = "Relative alignment shift in X vs Y";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_relRotX2DHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          rotXNBin = histoInfo->_xBin;
          rotXMin  = histoInfo->_xMin;
          rotXMax  = histoInfo->_xMax;
          rotVNBin = histoInfo->_yBin;
          rotVMin  = histoInfo->_yMin;
          rotVMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) rotXTitle = histoInfo->_title;
        }
    }

  for(int ipl=0;ipl<_nTelPlanes; ipl++)
    {
      if(_isActive[ipl] && ipl!=_referenceID0 && ipl!=_referenceID1)
        {
          stringstream nam,tit;

          nam << _relRotX2DHistoName << "_" << _planeID[ ipl ] ;
          tempHistoName=nam.str();

          tit << rotXTitle << " for plane " << _planeID[ ipl ] << " w.r.t. planes "
              << _planeID[ _referenceID0 ] << " and " << _planeID[ _referenceID1 ];
          tempHistoTitle=tit.str();

          AIDA::IHistogram2D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(),rotXNBin,rotXMin,rotXMax,rotVNBin,rotVMin,rotVMax);

          tempHisto->setTitle(tempHistoTitle.c_str());

          _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));

        }

    }


  // Relative alignment histograms: Y vs X (2D plot)

  rotYNBin  = 100;
  rotYMin   = -5.;
  rotYMax   = 5.;
  rotVNBin  = 100;
  rotVMin   = -0.2;
  rotVMax   = 0.2;
  rotYTitle = "Relative alignment shift in Y vs X";


  if ( isHistoManagerAvailable )
    {
      histoInfo = histoMgr->getHistogramInfo(_relRotY2DHistoName);
      if ( histoInfo )
        {
          streamlog_out ( DEBUG5 ) << (* histoInfo ) << endl;
          rotYNBin = histoInfo->_xBin;
          rotYMin  = histoInfo->_xMin;
          rotYMax  = histoInfo->_xMax;
          rotVNBin = histoInfo->_yBin;
          rotVMin  = histoInfo->_yMin;
          rotVMax  = histoInfo->_yMax;
          if ( histoInfo->_title != "" ) rotYTitle = histoInfo->_title;
        }
    }

  for(int ipl=0;ipl<_nTelPlanes; ipl++)
    {
      if(_isActive[ipl] && ipl!=_referenceID0 && ipl!=_referenceID1)
        {
          stringstream nam,tit;

          nam << _relRotY2DHistoName << "_" << _planeID[ ipl ] ;
          tempHistoName=nam.str();

          tit << rotYTitle << " for plane " << _planeID[ ipl ] << " w.r.t. planes "
              << _planeID[ _referenceID0 ] << " and " << _planeID[ _referenceID1 ];
          tempHistoTitle=tit.str();

          AIDA::IHistogram2D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(),rotYNBin,rotYMin,rotYMax,rotVNBin,rotVMin,rotVMax);

          tempHisto->setTitle(tempHistoTitle.c_str());

          _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));

        }

    }

  } // end of alignment histogram booking if:  if(_alignCheckHistograms){


// List all booked histogram - check of histogram map filling

  streamlog_out ( MESSAGE5 ) <<  _aidaHistoMap.size() << " histograms booked" << endl;


  map<string, AIDA::IBaseHistogram *>::iterator mapIter;
  for(mapIter = _aidaHistoMap.begin(); mapIter != _aidaHistoMap.end() ; mapIter++ ) {
    streamlog_out ( DEBUG5 ) <<  mapIter->first << " : " <<  (mapIter->second)->title() << endl;
  }
  streamlog_out ( DEBUG5 ) << "Histogram booking completed \n\n" << endl;


  return;
}

#endif // GEAR && AIDA
