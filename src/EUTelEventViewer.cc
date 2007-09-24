// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelEventViewer.cc,v 1.5 2007-09-24 01:20:06 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h" 
#include "EUTELESCOPE.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelEventViewer.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"

// MarlinUtil
#include "MarlinCED.h"
#include "ClusterShapes.h"
#include "ced_cli.h"

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>

// system includes <>
#include <vector>
#include <string>
#include <memory>

using namespace std ;
using namespace lcio ;
using namespace marlin ;
using namespace eutelescope;


EUTelEventViewer::EUTelEventViewer() : Processor("EUTelEventViewer") {

  _description = "Event display" ;
  
  vector<string > trackerHitCollectionNameVecExample;
  trackerHitCollectionNameVecExample.push_back("hit");
  trackerHitCollectionNameVecExample.push_back("testfithit");

  registerInputCollections( LCIO::TRACKERHIT, "TrackerHitCollections",
			    "Tracker hit collection names",  _trackerHitCollectionNameVec,
			    trackerHitCollectionNameVecExample );


  vector<string > trackCollectionNameVecExample;
  trackCollectionNameVecExample.push_back("testfittrack");
  
  registerInputCollections( LCIO::TRACK, "TrackCollections",
			    "Track collection names", _trackCollectionNameVec,
			    trackCollectionNameVecExample );

  registerProcessorParameter("LayerTrackerHit", 
			     "Layer for Tracker Hits",
			     _layerTrackerHit, 
			     (int)-1);

  registerProcessorParameter("LayerTracks", 
			     "Layer for Tracks",
			     _layerTrack, 
			     (int)-1);

  registerProcessorParameter("DetectorModel",
			     "Detector Model",
			     _detModel,
			     (int)99999);

  registerProcessorParameter("WaitForKeyboard",
			     "Wait for Keyboard before proceed",
			     _waitForKeyboard,
			     (int)1);

}

void EUTelEventViewer::init() {

  printParameters();

  MarlinCED::init(this) ;

  _iEvt = 0;;

}


void EUTelEventViewer::processRunHeader( LCRunHeader * rdr ) { 
  
  auto_ptr<EUTelRunHeaderImpl> runHeader ( new EUTelRunHeaderImpl( rdr ) );
  runHeader->addProcessor( type() );

} 

void EUTelEventViewer::processEvent( LCEvent * evt ) { 

  streamlog_out( MESSAGE4 ) << "Processing event " 
			    << setw(6) << setiosflags(ios::right) << evt->getEventNumber() << " in run "
			    << setw(6) << setiosflags(ios::right) << setfill('0')  << evt->getRunNumber() << setfill(' ')
			    << " (Total = " << setw(10) << _iEvt << ")" << resetiosflags(ios::left) << endl;

  EUTelEventImpl * event = static_cast<EUTelEventImpl *> ( evt );
  
  if ( event->getEventType() == kEORE ) {
    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return;
  } else if ( event->getEventType() == kUNKNOWN ) {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
			       << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }


  // Drawing Geometry

  MarlinCED::newEvent( this, _detModel ) ;

  if ( _layerTrackerHit > 0 ) {
    for ( unsigned int iCollection = 0; iCollection < _trackerHitCollectionNameVec.size(); iCollection++ ) {
      try {
	LCCollection * collection = evt->getCollection( _trackerHitCollectionNameVec[iCollection].c_str() );
	for ( int iHit = 0; iHit < collection->getNumberOfElements(); iHit++ ) {
	  TrackerHitImpl * hit = dynamic_cast<TrackerHitImpl *> ( collection->getElementAt(iHit) ) ;
	  float x = static_cast<float > ( hit->getPosition()[0] );
	  float y = static_cast<float > ( hit->getPosition()[1] );
	  float z = static_cast<float > ( hit->getPosition()[2] );
	  unsigned int color =  returnColor(iCollection);
	  ced_hit(x,y,z, _layerTrackerHit << CED_LAYER_SHIFT,2,color);
	}
      }  catch (DataNotAvailableException& e) {
	streamlog_out ( WARNING2 ) << "No input collection (" << _trackerHitCollectionNameVec[iCollection] << " found on " 
				   <<  event->getEventNumber() 
				   << " in run " << event->getRunNumber() << endl;
      }
    }
  }


  // Drawing  Tracks
  if (_layerTrack >= 0) {
    for ( unsigned int iCollection = 0; iCollection < _trackCollectionNameVec.size(); iCollection++ ) {
      try {
	LCCollection * collection =    evt->getCollection(_trackCollectionNameVec[iCollection].c_str());
	for ( int iTrack = 0; iTrack < collection->getNumberOfElements(); iTrack++ ) {
	  TrackImpl * track = dynamic_cast<TrackImpl*> ( collection->getElementAt(iTrack) );
	  TrackerHitVec hitvec = track->getTrackerHits();
	  int nHits = static_cast<int> (hitvec.size());
	  float * ah = new float[nHits];
	  float * xh = new float[nHits];
	  float * yh = new float[nHits];
	  float * zh = new float[nHits];
	  float zmin = 1.0E+10;
	  float zmax = -1.0E+10;
	  for (int iHit = 0; iHit < nHits; ++iHit) {
	    TrackerHit * hit = hitvec[iHit];
	    float x = (float)hit->getPosition()[0];
	    float y = (float)hit->getPosition()[1];
	    float z = (float)hit->getPosition()[2];
	    int kcol = returnColor(iTrack);
	    ced_hit(x,y,z, _layerTrack<<CED_LAYER_SHIFT,2,kcol);
	    ah[iHit] = 1.0;
	    xh[iHit] = x;
	    yh[iHit] = y;
	    zh[iHit] = z;
	    if (z < zmin)
	      zmin = z;
	    if (z > zmax)
	      zmax = z;
	  } 	
	  ClusterShapes * shapes = new ClusterShapes(nHits,ah,xh,yh,zh);
	  float dz = (zmax - zmin) / 500.;
	  float par[5];
	  float dpar[5];
	  float chi2;
	  float distmax;
	  shapes->FitHelix(500, 0, 1, par, dpar, chi2, distmax);
	  float x0 = par[0];
	  float y0 = par[1];
	  float r0 = par[2];
	  float bz = par[3];
	  float phi0 = par[4];
	  if (chi2 > 0. && chi2 < 10.) {
	    for (int iz(0); iz < 500; ++iz) {
	      float z1 = zmin + iz*dz;
	      float z2 = z1 + dz;
	      float x1 = x0 + r0*cos(bz*z1+phi0);
	      float y1 = y0 + r0*sin(bz*z1+phi0);
	      float x2 = x0 + r0*cos(bz*z2+phi0);
	      float y2 = y0 + r0*sin(bz*z2+phi0);
	      ced_line(x1,y1,z1,x2,y2,z2,_layerTrack<<CED_LAYER_SHIFT,2,0xFFFFFF);
	    }
	  }
	  //   ced_send_event();
	  delete shapes;
	  delete[] xh;
	  delete[] yh;
	  delete[] zh;
	  delete[] ah;
	}
      }
      catch(DataNotAvailableException &e){
	streamlog_out ( WARNING2 ) << "No input collection (" << _trackerHitCollectionNameVec[iCollection] << " found on " 
				   <<  event->getEventNumber() 
				   << " in run " << event->getRunNumber() << endl;
      }	
    }
  }
  
  MarlinCED::draw( this, _waitForKeyboard ) ;

}


void EUTelEventViewer::check( LCEvent * /* evt */ ) { }
  
void EUTelEventViewer::end(){ } 

int EUTelEventViewer::returnColor(int counter) {

  int icol =  counter % 16;
  int kcol =  0x00ff00;
  if (icol==1)  kcol = 0xAA00ff;
  if (icol==2)  kcol = 0xff0000;
  if (icol==3)  kcol = 0x00ffff;
  if (icol==4)  kcol = 0xffff00;
  if (icol==5)  kcol = 0xff00ff;
  if (icol==6)  kcol = 0xffffff;
  if (icol==7)  kcol = 0x0fffff;
  if (icol==8)  kcol = 0x000088;
  if (icol==9)  kcol = 0x008800;
  if (icol==10) kcol = 0x880000;
  if (icol==11) kcol = 0x008888;
  if (icol==12) kcol = 0x888800;
  if (icol==13) kcol = 0x880088;
  if (icol==14) kcol = 0x888888;
  if (icol==15) kcol = 0x00A888;

  return kcol;

}
