// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-

// Author: A.F.Zarnecki, University of Warsaw <mailto:zarnecki@fuw.edu.pl>
// Version: $Id: EUTelFitHistograms.cc,v 1.0 2007/09/10
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
std::string EUTelFitHistograms::_beamShiftXHistoName = "beamShiftX";
std::string EUTelFitHistograms::_beamShiftYHistoName = "beamShiftY";
std::string EUTelFitHistograms::_beamRotXHistoName   = "beamRotX";
std::string EUTelFitHistograms::_beamRotYHistoName   = "beamRotY";
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


  // beam alignment

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


    string tempHistoTitle;
    string tempHistoName;

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
