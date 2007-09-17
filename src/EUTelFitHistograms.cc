// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-

// Author: A.F.Zarnecki, University of Warsaw <mailto:zarnecki@fuw.edu.pl>
// Version: $Id: EUTelFitHistograms.cc,v 1.2 2007-09-17 22:24:23 zarnecki Exp $
// Date 2007.09.10

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope inlcudes
#include "EUTelFitHistograms.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelHistogramManager.h"
#include "EUTelExceptions.h"

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
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

std::string EUTelFitHistograms::_beamShiftXHistoName = "beamShiftX";
std::string EUTelFitHistograms::_beamShiftYHistoName = "beamShiftY";
std::string EUTelFitHistograms::_beamRotXHistoName   = "beamRotX";
std::string EUTelFitHistograms::_beamRotYHistoName   = "beamRotY";

std::string EUTelFitHistograms::_relShiftXHistoName = "relShiftX";
std::string EUTelFitHistograms::_relShiftYHistoName = "relShiftY";
std::string EUTelFitHistograms::_relRotXHistoName   = "relRotX";
std::string EUTelFitHistograms::_relRotYHistoName   = "relRotY";
#endif


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
         << "  at Z [mm] = " << _planePosition[ipl]; 
      
      message<MESSAGE> ( log() << ss.str() );
    }

  
    // Allocate arrays for track fitting

  _isMeasured      = new bool[_nTelPlanes];
  _isFitted        = new bool[_nTelPlanes];

  _measuredX     = new double[_nTelPlanes];
  _measuredY     = new double[_nTelPlanes];
  _fittedX       = new double[_nTelPlanes];
  _fittedY       = new double[_nTelPlanes];


// Book histograms

#ifdef MARLIN_USE_AIDA
    bookHistos();
#endif

}

void EUTelFitHistograms::processRunHeader( LCRunHeader* runHeader) { 

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

  if(nDet)message<MESSAGE> ( log() << nDet << " subdetectors defined :" );
  stringstream ss;
  for(int idet=0;idet<nDet;idet++)  message<MESSAGE> (log()  << idet+1 << " : " << subDets->at(idet) );


} 

void EUTelFitHistograms::processEvent( LCEvent * event ) { 
  
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
   throw SkipEventException(this);
  }
    

  // Loop over tracks in input collections

  int nTrack = col->getNumberOfElements()  ;

  if(debug)message<DEBUG> ( log() << "Total of " << nTrack << " tracks in input collection " );


  for(int itrack=0; itrack< nTrack ; itrack++)
    {
      Track * fittrack = dynamic_cast<Track*>( col->getElementAt(itrack) ) ;

      // Hit list asign to track

      std::vector<EVENT::TrackerHit*>  trackhits = fittrack->getTrackerHits();

  // Copy hits assign to the track to local table
  // Assign hits to sensor planes


   int nHit =   trackhits.size();

   if(debug)message<DEBUG> ( log() << "Track " << itrack << " with " << nHit << " hits " );


 // Clear plane tables

      for(int ipl=0;ipl<_nTelPlanes;ipl++)
	{
	  _isMeasured[ipl]=false;
	  _isFitted[ipl]=false;
        }


 // Loop over hits and fill hit tables

  for(int ihit=0; ihit< nHit ; ihit++)
    {
      TrackerHit * meshit = trackhits.at(ihit); 

      // Hit position

      const double * pos = meshit->getPosition();

      // We find plane number of the hit
      // by looking at the Z position

      double distMin = 1.;
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
	message<ERROR> ( log() << "Hit outside telescope plane at z [mm] = "  << pos[2] );
        continue;
        }


      if( meshit->getType() < 32 )
        {
	  // Measured hits

	  _isMeasured[hitPlane]=true;

          _measuredX[hitPlane]=pos[0];
          _measuredY[hitPlane]=pos[1];
	}
      else
	{
	  // Measured hits

	  _isFitted[hitPlane]=true;

          _fittedX[hitPlane]=pos[0];
          _fittedY[hitPlane]=pos[1];
	}

    }

  // Histograms of measured positions

   for(int ipl=0;ipl<_nTelPlanes; ipl++)
      {
      if(_isMeasured[ipl])
 	{
        string tempHistoName;

        stringstream nam;
        nam << _MeasuredXHistoName << "_" << ipl ;
        tempHistoName=nam.str();
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredX[ipl]);

        stringstream nam2;
        nam2 << _MeasuredYHistoName << "_" << ipl ;
        tempHistoName=nam2.str();
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredY[ipl]);

        stringstream nam3;
        nam3 << _MeasuredXYHistoName << "_" << ipl ;
        tempHistoName=nam3.str();
        (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredX[ipl],_measuredY[ipl]);

	}
      }



  // Histograms of fitted positions

   for(int ipl=0;ipl<_nTelPlanes; ipl++)
      {
      if(_isFitted[ipl])
 	{
        string tempHistoName;

        stringstream nam;
        nam << _FittedXHistoName << "_" << ipl ;
        tempHistoName=nam.str();
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(_fittedX[ipl]);

        stringstream nam2;
        nam2 << _FittedYHistoName << "_" << ipl ;
        tempHistoName=nam2.str();
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(_fittedY[ipl]);

        stringstream nam3;
        nam3 << _FittedXYHistoName << "_" << ipl ;
        tempHistoName=nam3.str();
        (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[tempHistoName]))->fill(_fittedX[ipl],_fittedY[ipl]);

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
        nam << _AngleXHistoName << "_" << ipl ;
        tempHistoName=nam.str();
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(angleX);

        stringstream nam2;
        nam2 << _AngleYHistoName << "_" << ipl ;
        tempHistoName=nam2.str();
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(angleY);

        stringstream nam3;
        nam3 << _AngleXYHistoName << "_" << ipl ;
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
        nam << _ScatXHistoName << "_" << ipl ;
        tempHistoName=nam.str();
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(scatX);

        stringstream nam2;
        nam2 << _ScatYHistoName << "_" << ipl ;
        tempHistoName=nam2.str();
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(scatY);

        stringstream nam3;
        nam3 << _ScatXYHistoName << "_" << ipl ;
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
        nam << _ResidualXHistoName << "_" << ipl ;
        tempHistoName=nam.str();
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(_fittedX[ipl]-_measuredX[ipl]);

        stringstream nam2;
        nam2 << _ResidualYHistoName << "_" << ipl ;
        tempHistoName=nam2.str();
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(_fittedY[ipl]-_measuredY[ipl]);

        stringstream nam3;
        nam3 << _ResidualXYHistoName << "_" << ipl ;
        tempHistoName=nam3.str();
        (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[tempHistoName]))->fill(_fittedX[ipl]-_measuredX[ipl],_fittedY[ipl]-_measuredY[ipl]);

	}
      }



  // Alignment w.r.t the beam direction

  if(_isMeasured[_beamID])
    {
    for(int ipl=0;ipl<_nTelPlanes; ipl++)
      {
      if(ipl!=_beamID && _isMeasured[ipl])
 	{
        string tempHistoName;

        stringstream nam;
        nam << _beamShiftXHistoName << "_" << ipl ;
        tempHistoName=nam.str();
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredX[ipl]-_measuredX[_beamID]);

        stringstream nam2;
        nam2 << _beamShiftYHistoName << "_" << ipl ;
        tempHistoName=nam2.str();
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredY[ipl]-_measuredY[_beamID]);

        stringstream nam3;
        nam3 << _beamRotXHistoName << "_" << ipl ;
        tempHistoName=nam3.str();
        (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredY[_beamID],_measuredX[ipl]-_measuredX[_beamID]);

        stringstream nam4;
        nam4 << _beamRotYHistoName << "_" << ipl ;
        tempHistoName=nam4.str();
        (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredX[_beamID],_measuredY[ipl]-_measuredY[_beamID]);
	}
      }
    }


  // Alignment w.r.t two selected planes

  if(_isMeasured[_referenceID0] && _isMeasured[_referenceID1] )
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
        nam << _relShiftXHistoName << "_" << ipl ;
        tempHistoName=nam.str();
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredX[ipl]-lineX);

        stringstream nam2;
        nam2 << _relShiftYHistoName << "_" << ipl ;
        tempHistoName=nam2.str();
        (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[tempHistoName]))->fill(_measuredY[ipl]-lineY);

        stringstream nam3;
        nam3 << _relRotXHistoName << "_" << ipl ;
        tempHistoName=nam3.str();
        (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[tempHistoName]))->fill(lineY,_measuredX[ipl]-lineX);

        stringstream nam4;
        nam4 << _relRotYHistoName << "_" << ipl ;
        tempHistoName=nam4.str();
        (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[tempHistoName]))->fill(lineX,_measuredY[ipl]-lineY);
	}
      }
    }


  // End of loop over tracks
    }


  return;
}



void EUTelFitHistograms::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void EUTelFitHistograms::end(){ 
  
  //   std::cout << "EUTelFitHistograms::end()  " << name() 
  // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
  // 	    << std::endl ;


  // Clean memory 

  delete [] _planeSort ;
  delete [] _planePosition ;
  delete [] _planeID ;
  delete [] _isActive ;

  delete [] _isMeasured ;
  delete [] _isFitted ;
  delete [] _measuredX  ;
  delete [] _measuredY  ;
  delete [] _fittedX ;
  delete [] _fittedY ;

}



void EUTelFitHistograms::bookHistos() 
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
	 message<DEBUG> ( log() << (* histoInfo ) );
	 measXNBin = histoInfo->_xBin;
	 measXMin  = histoInfo->_xMin;
	 measXMax  = histoInfo->_xMax;
	 if ( histoInfo->_title != "" ) measXTitle = histoInfo->_title;
         }
      }


    string tempHistoTitle;
    string tempHistoName;

    for(int ipl=0;ipl<_nTelPlanes; ipl++)
      {
      if(_isActive[ipl])
 	{
        stringstream nam,tit;

        nam << _MeasuredXHistoName << "_" << ipl ;
        tempHistoName=nam.str();

        tit << measXTitle << " for plane " << ipl ;
        tempHistoTitle=tit.str();

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
	 message<DEBUG> ( log() << (* histoInfo ) );
	 measYNBin = histoInfo->_xBin;
	 measYMin  = histoInfo->_xMin;
	 measYMax  = histoInfo->_xMax;
	 if ( histoInfo->_title != "" ) measYTitle = histoInfo->_title;
         }
      }


    for(int ipl=0;ipl<_nTelPlanes; ipl++)
      {
      if(_isActive[ipl])
 	{
        stringstream nam,tit;

        nam << _MeasuredYHistoName << "_" << ipl ;
        tempHistoName=nam.str();

        tit << measYTitle << " for plane " << ipl ;
        tempHistoTitle=tit.str();

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


    for(int ipl=0;ipl<_nTelPlanes; ipl++)
      {
      if(_isActive[ipl])
 	{
        stringstream nam,tit;

        nam << _MeasuredXYHistoName << "_" << ipl ;
        tempHistoName=nam.str();

        tit << measXYTitle << " for plane " << ipl ;
        tempHistoTitle=tit.str();

        AIDA::IHistogram2D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(),measXNBin,measXMin,measXMax,measYNBin,measYMin,measYMax);

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
	 message<DEBUG> ( log() << (* histoInfo ) );
	 fitXNBin = histoInfo->_xBin;
	 fitXMin  = histoInfo->_xMin;
	 fitXMax  = histoInfo->_xMax;
	 if ( histoInfo->_title != "" ) fitXTitle = histoInfo->_title;
         }
      }


    for(int ipl=0;ipl<_nTelPlanes; ipl++)
      {
      if(_isActive[ipl])
 	{
        stringstream nam,tit;

        nam << _FittedXHistoName << "_" << ipl ;
        tempHistoName=nam.str();

        tit << fitXTitle << " for plane " << ipl ;
        tempHistoTitle=tit.str();

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
	 message<DEBUG> ( log() << (* histoInfo ) );
	 fitYNBin = histoInfo->_xBin;
	 fitYMin  = histoInfo->_xMin;
	 fitYMax  = histoInfo->_xMax;
	 if ( histoInfo->_title != "" ) fitYTitle = histoInfo->_title;
         }
      }


    for(int ipl=0;ipl<_nTelPlanes; ipl++)
      {
      if(_isActive[ipl])
 	{
        stringstream nam,tit;

        nam << _FittedYHistoName << "_" << ipl ;
        tempHistoName=nam.str();

        tit << fitYTitle << " for plane " << ipl ;
        tempHistoTitle=tit.str();

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


    for(int ipl=0;ipl<_nTelPlanes; ipl++)
      {
      if(_isActive[ipl])
 	{
        stringstream nam,tit;

        nam << _FittedXYHistoName << "_" << ipl ;
        tempHistoName=nam.str();

        tit << fitXYTitle << " for plane " << ipl ;
        tempHistoTitle=tit.str();

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
	 message<DEBUG> ( log() << (* histoInfo ) );
	 angleXNBin = histoInfo->_xBin;
	 angleXMin  = histoInfo->_xMin;
	 angleXMax  = histoInfo->_xMax;
	 if ( histoInfo->_title != "" ) angleXTitle = histoInfo->_title;
         }
      }


    for(int ipl=1;ipl<_nTelPlanes; ipl++)
      {
      if(_isActive[ipl])
 	{
        stringstream nam,tit;

        nam << _AngleXHistoName << "_" << ipl ;
        tempHistoName=nam.str();

        tit << angleXTitle << " for plane " << ipl ;
        tempHistoTitle=tit.str();

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
	 message<DEBUG> ( log() << (* histoInfo ) );
	 angleYNBin = histoInfo->_xBin;
	 angleYMin  = histoInfo->_xMin;
	 angleYMax  = histoInfo->_xMax;
	 if ( histoInfo->_title != "" ) angleYTitle = histoInfo->_title;
         }
      }


    for(int ipl=1;ipl<_nTelPlanes; ipl++)
      {
      if(_isActive[ipl])
 	{
        stringstream nam,tit;

        nam << _AngleYHistoName << "_" << ipl ;
        tempHistoName=nam.str();

        tit << angleYTitle << " for plane " << ipl ;
        tempHistoTitle=tit.str();

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
	 message<DEBUG> ( log() << (* histoInfo ) );
	 angleXNBin = histoInfo->_xBin;
	 angleXMin  = histoInfo->_xMin;
	 angleXMax  = histoInfo->_xMax;
	 angleYNBin = histoInfo->_yBin;
	 angleYMin  = histoInfo->_yMin;
	 angleYMax  = histoInfo->_yMax;
	 if ( histoInfo->_title != "" ) angleXYTitle = histoInfo->_title;
         }
      }


    for(int ipl=1;ipl<_nTelPlanes; ipl++)
      {
      if(_isActive[ipl])
 	{
        stringstream nam,tit;

        nam << _AngleXYHistoName << "_" << ipl ;
        tempHistoName=nam.str();

        tit << angleXYTitle << " for plane " << ipl ;
        tempHistoTitle=tit.str();

        AIDA::IHistogram2D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(),angleXNBin,angleXMin,angleXMax,angleYNBin,angleYMin,angleYMax);

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
	 message<DEBUG> ( log() << (* histoInfo ) );
	 scatXNBin = histoInfo->_xBin;
	 scatXMin  = histoInfo->_xMin;
	 scatXMax  = histoInfo->_xMax;
	 if ( histoInfo->_title != "" ) scatXTitle = histoInfo->_title;
         }
      }


    for(int ipl=1;ipl<_nTelPlanes-1; ipl++)
      {
      if(_isActive[ipl])
 	{
        stringstream nam,tit;

        nam << _ScatXHistoName << "_" << ipl ;
        tempHistoName=nam.str();

        tit << scatXTitle << " for plane " << ipl ;
        tempHistoTitle=tit.str();

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
	 message<DEBUG> ( log() << (* histoInfo ) );
	 scatYNBin = histoInfo->_xBin;
	 scatYMin  = histoInfo->_xMin;
	 scatYMax  = histoInfo->_xMax;
	 if ( histoInfo->_title != "" ) scatYTitle = histoInfo->_title;
         }
      }


    for(int ipl=1;ipl<_nTelPlanes-1; ipl++)
      {
      if(_isActive[ipl])
 	{
        stringstream nam,tit;

        nam << _ScatYHistoName << "_" << ipl ;
        tempHistoName=nam.str();

        tit << scatYTitle << " for plane " << ipl ;
        tempHistoTitle=tit.str();

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
	 message<DEBUG> ( log() << (* histoInfo ) );
	 scatXNBin = histoInfo->_xBin;
	 scatXMin  = histoInfo->_xMin;
	 scatXMax  = histoInfo->_xMax;
	 scatYNBin = histoInfo->_yBin;
	 scatYMin  = histoInfo->_yMin;
	 scatYMax  = histoInfo->_yMax;
	 if ( histoInfo->_title != "" ) scatXYTitle = histoInfo->_title;
         }
      }


    for(int ipl=1;ipl<_nTelPlanes-1; ipl++)
      {
      if(_isActive[ipl])
 	{
        stringstream nam,tit;

        nam << _ScatXYHistoName << "_" << ipl ;
        tempHistoName=nam.str();

        tit << scatXYTitle << " for plane " << ipl ;
        tempHistoTitle=tit.str();

        AIDA::IHistogram2D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(),scatXNBin,scatXMin,scatXMax,scatYNBin,scatYMin,scatYMax);

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
	 message<DEBUG> ( log() << (* histoInfo ) );
	 residXNBin = histoInfo->_xBin;
	 residXMin  = histoInfo->_xMin;
	 residXMax  = histoInfo->_xMax;
	 if ( histoInfo->_title != "" ) residXTitle = histoInfo->_title;
         }
      }


    for(int ipl=0;ipl<_nTelPlanes; ipl++)
      {
      if(_isActive[ipl])
 	{
        stringstream nam,tit;

        nam << _ResidualXHistoName << "_" << ipl ;
        tempHistoName=nam.str();

        tit << residXTitle << " for plane " << ipl ;
        tempHistoTitle=tit.str();

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
	 message<DEBUG> ( log() << (* histoInfo ) );
	 residYNBin = histoInfo->_xBin;
	 residYMin  = histoInfo->_xMin;
	 residYMax  = histoInfo->_xMax;
	 if ( histoInfo->_title != "" ) residYTitle = histoInfo->_title;
         }
      }


    for(int ipl=0;ipl<_nTelPlanes; ipl++)
      {
      if(_isActive[ipl])
 	{
        stringstream nam,tit;

        nam << _ResidualYHistoName << "_" << ipl ;
        tempHistoName=nam.str();

        tit << residYTitle << " for plane " << ipl ;
        tempHistoTitle=tit.str();

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
	 message<DEBUG> ( log() << (* histoInfo ) );
	 residXNBin = histoInfo->_xBin;
	 residXMin  = histoInfo->_xMin;
	 residXMax  = histoInfo->_xMax;
	 residYNBin = histoInfo->_yBin;
	 residYMin  = histoInfo->_yMin;
	 residYMax  = histoInfo->_yMax;
	 if ( histoInfo->_title != "" ) residXYTitle = histoInfo->_title;
         }
      }


    for(int ipl=0;ipl<_nTelPlanes; ipl++)
      {
      if(_isActive[ipl])
 	{
        stringstream nam,tit;

        nam << _ResidualXYHistoName << "_" << ipl ;
        tempHistoName=nam.str();

        tit << residXYTitle << " for plane " << ipl ;
        tempHistoTitle=tit.str();

        AIDA::IHistogram2D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(),residXNBin,residXMin,residXMax,residYNBin,residYMin,residYMax);

         tempHisto->setTitle(tempHistoTitle.c_str());

        _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));

        }

      }




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
	 message<DEBUG> ( log() << (* histoInfo ) );
	 shiftXNBin = histoInfo->_xBin;
	 shiftXMin  = histoInfo->_xMin;
	 shiftXMax  = histoInfo->_xMax;
	 if ( histoInfo->_title != "" ) shiftXTitle = histoInfo->_title;
         }
      }


    for(int ipl=0;ipl<_nTelPlanes; ipl++)
      {
      if(ipl!=_beamID)
 	{
        stringstream nam,tit;

        nam << _beamShiftXHistoName << "_" << ipl ;
        tempHistoName=nam.str();

        tit << shiftXTitle << " for plane " << ipl << " w.r.t. plane " << _beamID ;
        tempHistoTitle=tit.str();

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
	 message<DEBUG> ( log() << (* histoInfo ) );
	 shiftYNBin = histoInfo->_xBin;
	 shiftYMin  = histoInfo->_xMin;
	 shiftYMax  = histoInfo->_xMax;
	 if ( histoInfo->_title != "" ) shiftYTitle = histoInfo->_title;
         }
      }


    for(int ipl=0;ipl<_nTelPlanes; ipl++)
      {
      if(ipl!=_beamID)
 	{
        stringstream nam,tit;

        nam << _beamShiftYHistoName << "_" << ipl ;
        tempHistoName=nam.str();

        tit << shiftYTitle <<  " for plane " << ipl << " w.r.t. plane " << _beamID ;
        tempHistoTitle=tit.str();

        AIDA::IHistogram1D * tempHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( tempHistoName.c_str(),shiftYNBin,shiftYMin,shiftYMax);

         tempHisto->setTitle(tempHistoTitle.c_str());

        _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));

        }

      }


 // Beam alignment histograms: X vs Y to extract rotation in Z

    int    rotXNBin  = 200;
    double rotXMin   = -5.;
    double rotXMax   = 5.;
    string rotXTitle = "Beam alignment shift in X vs Y";


    if ( isHistoManagerAvailable ) 
      {
      histoInfo = histoMgr->getHistogramInfo(_beamRotXHistoName);
      if ( histoInfo ) 
         {
	 message<DEBUG> ( log() << (* histoInfo ) );
	 rotXNBin = histoInfo->_xBin;
	 rotXMin  = histoInfo->_xMin;
	 rotXMax  = histoInfo->_xMax;
	 if ( histoInfo->_title != "" ) rotXTitle = histoInfo->_title;
         }
      }

    for(int ipl=0;ipl<_nTelPlanes; ipl++)
      {
      if(ipl!=_beamID)
 	{
        stringstream nam,tit;

        nam << _beamRotXHistoName << "_" << ipl ;
        tempHistoName=nam.str();

        tit << rotXTitle << " for plane " << ipl << " w.r.t. plane " << _beamID ;
        tempHistoTitle=tit.str();

        AIDA::IProfile1D * tempHisto = AIDAProcessor::histogramFactory(this)->createProfile1D( tempHistoName.c_str(),rotXNBin,rotXMin,rotXMax);

         tempHisto->setTitle(tempHistoTitle.c_str());

        _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));

        }

      }


 // Beam alignment histograms: Y vs X to extract rotation in Z

    int    rotYNBin  = 200;
    double rotYMin   = -5.;
    double rotYMax   = 5.;
    string rotYTitle = "Beam alignment shift in Y vs X";


    if ( isHistoManagerAvailable ) 
      {
      histoInfo = histoMgr->getHistogramInfo(_beamRotYHistoName);
      if ( histoInfo ) 
         {
	 message<DEBUG> ( log() << (* histoInfo ) );
	 rotYNBin = histoInfo->_xBin;
	 rotYMin  = histoInfo->_xMin;
	 rotYMax  = histoInfo->_xMax;
	 if ( histoInfo->_title != "" ) rotYTitle = histoInfo->_title;
         }
      }

    for(int ipl=0;ipl<_nTelPlanes; ipl++)
      {
      if(ipl!=_beamID)
 	{
        stringstream nam,tit;

        nam << _beamRotYHistoName << "_" << ipl ;
        tempHistoName=nam.str();

        tit << rotYTitle << " for plane " << ipl << " w.r.t. plane " << _beamID ;
        tempHistoTitle=tit.str();

        AIDA::IProfile1D * tempHisto = AIDAProcessor::histogramFactory(this)->createProfile1D( tempHistoName.c_str(),rotYNBin,rotYMin,rotYMax);

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
	 message<DEBUG> ( log() << (* histoInfo ) );
	 shiftXNBin = histoInfo->_xBin;
	 shiftXMin  = histoInfo->_xMin;
	 shiftXMax  = histoInfo->_xMax;
	 if ( histoInfo->_title != "" ) shiftXTitle = histoInfo->_title;
         }
      }


    for(int ipl=0;ipl<_nTelPlanes; ipl++)
      {
      if(ipl!=_referenceID0 && ipl!=_referenceID1)
 	{
        stringstream nam,tit;

        nam << _relShiftXHistoName << "_" << ipl ;
        tempHistoName=nam.str();

        tit << shiftXTitle << " for plane " << ipl << " w.r.t. planes "
	    << _referenceID0 << " and " << _referenceID1 ;
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
	 message<DEBUG> ( log() << (* histoInfo ) );
	 shiftYNBin = histoInfo->_xBin;
	 shiftYMin  = histoInfo->_xMin;
	 shiftYMax  = histoInfo->_xMax;
	 if ( histoInfo->_title != "" ) shiftYTitle = histoInfo->_title;
         }
      }


    for(int ipl=0;ipl<_nTelPlanes; ipl++)
      {
      if(ipl!=_referenceID0 && ipl!=_referenceID1)
      if(ipl!=_beamID)
 	{
        stringstream nam,tit;

        nam << _relShiftYHistoName << "_" << ipl ;
        tempHistoName=nam.str();

        tit << shiftYTitle <<  " for plane " << ipl << " w.r.t. planes "
	    << _referenceID0 << " and " << _referenceID1 ;
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
    rotXTitle = "Relative alignment shift in X vs Y";


    if ( isHistoManagerAvailable ) 
      {
      histoInfo = histoMgr->getHistogramInfo(_relRotXHistoName);
      if ( histoInfo ) 
         {
	 message<DEBUG> ( log() << (* histoInfo ) );
	 rotXNBin = histoInfo->_xBin;
	 rotXMin  = histoInfo->_xMin;
	 rotXMax  = histoInfo->_xMax;
	 if ( histoInfo->_title != "" ) rotXTitle = histoInfo->_title;
         }
      }

    for(int ipl=0;ipl<_nTelPlanes; ipl++)
      {
      if(ipl!=_beamID)
      if(ipl!=_referenceID0 && ipl!=_referenceID1)
 	{
        stringstream nam,tit;

        nam << _relRotXHistoName << "_" << ipl ;
        tempHistoName=nam.str();

        tit << rotXTitle << " for plane " << ipl << " w.r.t. planes "
	    << _referenceID0 << " and " << _referenceID1 ;
        tempHistoTitle=tit.str();

        AIDA::IProfile1D * tempHisto = AIDAProcessor::histogramFactory(this)->createProfile1D( tempHistoName.c_str(),rotXNBin,rotXMin,rotXMax);

         tempHisto->setTitle(tempHistoTitle.c_str());

        _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));

        }

      }


 // Relative alignment histograms: Y vs X to extract rotation in Z

    rotYNBin  = 200;
    rotYMin   = -5.;
    rotYMax   = 5.;
    rotYTitle = "Relative alignment shift in Y vs X";


    if ( isHistoManagerAvailable ) 
      {
      histoInfo = histoMgr->getHistogramInfo(_relRotYHistoName);
      if ( histoInfo ) 
         {
	 message<DEBUG> ( log() << (* histoInfo ) );
	 rotYNBin = histoInfo->_xBin;
	 rotYMin  = histoInfo->_xMin;
	 rotYMax  = histoInfo->_xMax;
	 if ( histoInfo->_title != "" ) rotYTitle = histoInfo->_title;
         }
      }

    for(int ipl=0;ipl<_nTelPlanes; ipl++)
      {
      if(ipl!=_referenceID0 && ipl!=_referenceID1)
 	{
        stringstream nam,tit;

        nam << _relRotYHistoName << "_" << ipl ;
        tempHistoName=nam.str();

        tit << rotYTitle << " for plane " << ipl << " w.r.t. planes "
	    << _referenceID0 << " and " << _referenceID1 ;
        tempHistoTitle=tit.str();

        AIDA::IProfile1D * tempHisto = AIDAProcessor::histogramFactory(this)->createProfile1D( tempHistoName.c_str(),rotYNBin,rotYMin,rotYMax);

         tempHisto->setTitle(tempHistoTitle.c_str());

        _aidaHistoMap.insert(make_pair(tempHistoName, tempHisto));

        }

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
