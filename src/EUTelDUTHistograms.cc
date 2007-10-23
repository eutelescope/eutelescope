// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-

// Author: A.F.Zarnecki, University of Warsaw <mailto:zarnecki@fuw.edu.pl>
// @version: $Id: EUTelDUTHistograms.cc,v 1.1 2007-10-23 21:29:02 zarnecki Exp $
// Date 2007.09.15

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope inlcudes
#include "EUTelDUTHistograms.h"
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


using namespace std;
using namespace lcio ;
using namespace marlin ;
using namespace eutelescope;

// definition of static members mainly used to name histograms
#ifdef MARLIN_USE_AIDA
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

std::string EUTelDUTHistograms::_NoiseXHistoName  = "DUTnoiseX";
std::string EUTelDUTHistograms::_NoiseYHistoName  = "DUTnoiseY";
std::string EUTelDUTHistograms::_NoiseXYHistoName = "DUTnoiseXY";

std::string EUTelDUTHistograms::_ShiftXHistoName       = "DUTshiftX";
std::string EUTelDUTHistograms::_ShiftYHistoName       = "DUTshiftY";
std::string EUTelDUTHistograms::_ShiftXYHistoName      = "DUTshiftXY";
std::string EUTelDUTHistograms::_ShiftXvsYHistoName    = "DUTshiftXvsY";
std::string EUTelDUTHistograms::_ShiftYvsXHistoName    = "DUTshiftYvsX";
std::string EUTelDUTHistograms::_ShiftXvsY2DHistoName    = "DUTshiftXvsY2D";
std::string EUTelDUTHistograms::_ShiftYvsX2DHistoName    = "DUTshiftYvsX2D";

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

  registerProcessorParameter ("DistMax",
			      "Maximum allowed distance between fit and matched DUT hit",
			      _distMax,  static_cast < double > (0.1));


  std::vector<float > initAlign;
  initAlign.push_back(0.);
  initAlign.push_back(0.);
  initAlign.push_back(0.);

  registerProcessorParameter ("DUTalignment",
     "Alignment corrections for DUT: shift in X, Y and rotation around Z",
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

#ifdef MARLIN_USE_AIDA
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
  
  EUTelEventImpl * euEvent = static_cast<EUTelEventImpl*> ( event );
  if ( euEvent->getEventType() == kEORE ) {
    message<DEBUG> ( "EORE found: nothing else to do." );
    return;
  }

  bool debug = ( _debugCount>0 && _nEvt%_debugCount == 0);

  _nEvt ++ ;
  int evtNr = event->getEventNumber();


  if(debug)message<DEBUG> ( log() << "Processing record " << _nEvt << " == event " << evtNr );

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
   throw SkipEventException(this);
  }
    

  LCCollection* hitcol;
  try {
    hitcol = event->getCollection( _inputHitColName ) ;
  } catch (lcio::DataNotAvailableException& e) {
    message<ERROR> ( log() << "Not able to get collection " 
		     << _inputHitColName 
		     << "\nfrom event " << event->getEventNumber()
		     << " in run " << event->getRunNumber()  );
   throw SkipEventException(this);
  }
    
  // Clear local storage

  _measuredX.clear();
  _measuredY.clear();

  _fittedX.clear();
  _fittedY.clear();


  // Loop over tracks in input track collection
  // Read fitted positions at DUT

  int nTrack = trackcol->getNumberOfElements()  ;


  if(debug)message<DEBUG> ( log() << "Total of " << nTrack << " tracks in input collection " );



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
            break;
	    }
	  }
	}

      // End of loop over fitted tracks
    }

  if(debug)message<DEBUG> ( log() << _fittedX.size()  << " fitted positions at DUT " );

  // Histograms of fitted positions

  for(int ifit=0;ifit<(int)_fittedX.size(); ifit++)
      {
      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_FittedXHistoName]))->fill(_fittedX[ifit]);


      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_FittedYHistoName]))->fill(_fittedY[ifit]);

      (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_FittedXYHistoName]))->fill(_fittedX[ifit],_fittedY[ifit]);


      if(debug)message<DEBUG> ( log() << "Fit " << ifit
		                << "   X = " << _fittedX[ifit] 
				<< "   Y = " << _fittedY[ifit]) ;
      }



  // Loop over hits in input collection
  // Read measured positions at DUT

  int nHit = hitcol->getNumberOfElements()  ;

  if(debug)message<DEBUG> ( log() << "Total of " << nHit << " tracker hits in input collection " );


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
	}
    }
 

  if(debug)message<DEBUG> ( log() << _measuredX.size() << " hits at DUT " );


  // Histograms of measured positions

  for(int ihit=0;ihit<(int)_measuredX.size(); ihit++)
      {
      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_MeasuredXHistoName]))->fill(_measuredX[ihit]);


      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_MeasuredYHistoName]))->fill(_measuredY[ihit]);

      (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_MeasuredXYHistoName]))->fill(_measuredX[ihit],_measuredY[ihit]);

      if(debug)message<DEBUG> ( log() << "Hit " << ihit
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

      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_MatchedXHistoName]))->fill(_measuredX[besthit]);


      (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_MatchedYHistoName]))->fill(_measuredY[besthit]);

      (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_MatchedXYHistoName]))->fill(_measuredX[besthit],_measuredY[besthit]);

      // Histograms of measured-fitted shifts

   (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_ShiftXHistoName]))->fill(_measuredX[besthit]-_fittedX[bestfit]);

   (dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[_ShiftYHistoName]))->fill(_measuredY[besthit]-_fittedY[bestfit]);

   (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_ShiftXYHistoName]))->fill(_measuredX[besthit]-_fittedX[bestfit],_measuredY[besthit]-_fittedY[bestfit]);

   (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_ShiftXvsYHistoName]))->fill(_fittedY[bestfit],_measuredX[besthit]-_fittedX[bestfit]);

   (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_ShiftYvsXHistoName]))->fill(_fittedX[bestfit],_measuredY[besthit]-_fittedY[bestfit]);

   (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_ShiftXvsY2DHistoName]))->fill(_fittedY[bestfit],_measuredX[besthit]-_fittedX[bestfit]);

   (dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[_ShiftYvsX2DHistoName]))->fill(_fittedX[bestfit],_measuredY[besthit]-_fittedY[bestfit]);

   // Efficiency plots

   (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_EfficiencyXHistoName]))->fill(_fittedX[bestfit],1.);

   (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_EfficiencyYHistoName]))->fill(_fittedY[bestfit],1.);

   (dynamic_cast<AIDA::IProfile2D*> ( _aidaHistoMap[_EfficiencyXYHistoName]))->fill(_fittedX[bestfit],_fittedY[bestfit],1.);


   // Noise plots

   (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_NoiseXHistoName]))->fill(_measuredX[bestfit],0.);

   (dynamic_cast<AIDA::IProfile1D*> ( _aidaHistoMap[_NoiseYHistoName]))->fill(_measuredY[bestfit],0.);

   (dynamic_cast<AIDA::IProfile2D*> ( _aidaHistoMap[_NoiseXYHistoName]))->fill(_measuredX[bestfit],_measuredY[bestfit],0.);


   // Remove matched entries from the list (so the next matching pair
   // can be looked for)

   _fittedX.erase(_fittedX.begin()+bestfit);
   _fittedY.erase(_fittedY.begin()+bestfit);

   _measuredX.erase(_measuredX.begin()+besthit);
   _measuredY.erase(_measuredY.begin()+besthit);

    }


  // End of loop of matching DUT hits to fitted positions 
   }
   while(distmin < _distMax*_distMax);


   if(debug)
     {
     message<DEBUG> ( log() << nMatch << " DUT hits matched to fitted tracks ");
     message<DEBUG> ( log() << _measuredX.size() << " DUT hits not matched to any track ");
     message<DEBUG> ( log() << _fittedX.size() << " fitted tracks not matched to any DUT hit ");
     }


   // Efficiency plots - unmatched tracks

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

  return;
}



void EUTelDUTHistograms::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void EUTelDUTHistograms::end(){ 
  
  // Nothing to do here
}



void EUTelDUTHistograms::bookHistos() 
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

    int    measYNBin  = 400;
    double measYMin   = -5.;
    double measYMax   = 5.;
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

    measXNBin  = 100;
    measXMin   = -5.;
    measXMax   = 5.;
    measYNBin  = 100;
    measYMin   = -5.;
    measYMax   = 5.;
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



    AIDA::IHistogram1D * fitXHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_FittedXHistoName.c_str(),fitXNBin,fitXMin,fitXMax);

    fitXHisto->setTitle(fitXTitle.c_str());

    _aidaHistoMap.insert(make_pair(_FittedXHistoName, fitXHisto));



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




    AIDA::IHistogram1D * fitYHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(_FittedYHistoName.c_str(),fitYNBin,fitYMin,fitYMax);

    fitYHisto->setTitle(fitYTitle.c_str());

    _aidaHistoMap.insert(make_pair(_FittedYHistoName, fitYHisto));




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



   AIDA::IHistogram2D * fitXYHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( _FittedXYHistoName.c_str(),fitXNBin,fitXMin,fitXMax,fitYNBin,fitYMin,fitYMax);

   fitXYHisto->setTitle(fitXYTitle.c_str());

   _aidaHistoMap.insert(make_pair(_FittedXYHistoName, fitXYHisto));




 // Measured - fitted position in X 

    int    shiftXNBin  = 400;
    double shiftXMin   = -0.1;
    double shiftXMax   = 0.1;
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



 // Measured - fitted position in Y 

    int    shiftYNBin  = 400;
    double shiftYMin   = -0.1;
    double shiftYMax   = 0.1;
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





 // Measured - fitted position in X  vs Y

    shiftXNBin  = 200;
    shiftXMin   = -5.;
    shiftXMax   = 5.;
    double shiftVMin   = -0.1;
    double shiftVMax   = 0.1;
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



 // Measured - fitted position in Y vs X 

    shiftYNBin  = 200;
    shiftYMin   = -5.;
    shiftYMax   = 5.;
    shiftVMin   = -0.1;
    shiftVMax   = 0.1;
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




 // Measured - fitted position in X  vs Y (2D plot)

    shiftXNBin  = 200;
    shiftXMin   = -5.;
    shiftXMax   = 5.;
    int shiftVNBin  = 200;
    shiftVMin   = -0.1;
    shiftVMax   = 0.1;
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



 // Measured - fitted position in Y vs X  (2D plot)

    shiftYNBin  = 200;
    shiftYMin   = -5.;
    shiftYMax   = 5.;
    shiftVNBin  = 200;
    shiftVMin   = -0.1;
    shiftVMax   = 0.1;
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




 // Measured - fitted position in X-Y 

    shiftXNBin  = 100;
    shiftXMin   = -0.1;
    shiftXMax   = 0.1;
    shiftYNBin  = 100;
    shiftYMin   = -0.1;
    shiftYMax   = 0.1;
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






 // Efficiency as a function of the fitted position in X 

    int    effiXNBin  = 100;
    double effiXMin   = -5.;
    double effiXMax   = 5.;
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



 // Efficiency as a function of the fitted position in Y 

    int    effiYNBin  = 100;
    double effiYMin   = -5.;
    double effiYMax   = 5.;
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




 // Efficiency as a function of the fitted position in X-Y 

    effiXNBin  = 50;
    effiXMin   = -5.;
    effiXMax   = 5.;
    effiYNBin  = 50;
    effiYMin   = -5.;
    effiYMax   = 5.;
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



 // Noise as a function of the measured position in X 

    int    noiseXNBin  = 100;
    double noiseXMin   = -5.;
    double noiseXMax   = 5.;
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

    int    noiseYNBin  = 100;
    double noiseYMin   = -5.;
    double noiseYMax   = 5.;
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

    noiseXNBin  = 50;
    noiseXMin   = -5.;
    noiseXMax   = 5.;
    noiseYNBin  = 50;
    noiseYMin   = -5.;
    noiseYMax   = 5.;
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
