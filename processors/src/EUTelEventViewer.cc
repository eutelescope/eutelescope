// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// this processor is built only if GEAR and MARLINUTIL are available, otherwise
// it is simply skipped
#if defined(USE_GEAR) && defined(USE_MARLINUTIL)

// eutelescope includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelEventViewer.h"
#include "EUTelAlignmentConstant.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
# include "marlin/Global.h"

// MarlinUtil
#include "MarlinCED.h"
#include "ClusterShapes.h"
#include "ced_cli.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/LCEventImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCCollectionVec.h>

// system includes <>
#include <vector>
#include <string>
#include <memory>
#include <ctime>
#include <unistd.h>

using namespace std ;
using namespace lcio ;
using namespace marlin ;
using namespace eutelescope;


EUTelEventViewer::EUTelEventViewer()
 : Processor("EUTelEventViewer"),
   _trackerHitCollectionNameVec(),
   _trackCollectionNameVec(),
   _alignmentCollectionName(""),
   _layerTrackerHit(-1),
   _layerTrack(-1),
   _waitForKeyboard(true),
   _autoForwardDelay(4.0),
   _detModel(99999)
{

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
                             static_cast< int >(-1));

  registerProcessorParameter("LayerTracks",
                             "Layer for Tracks",
                             _layerTrack,
                             static_cast< int >(-1));

  registerProcessorParameter("DetectorModel",
                             "Detector Model:\n"
			     " 99999 to draw the ideal model taken from the GEAR description\n"
			     "    -1 to draw the model described by GEAR and corrected for alignment",
                             _detModel,
                             static_cast< int >(99999));

  registerProcessorParameter("WaitForKeyboard",
                             "Wait for Keyboard before proceed",
                             _waitForKeyboard,
                             static_cast< bool > (true) );

  registerProcessorParameter("AutoForwardDelay","This is the time in second between two following display\n"
                             "Enable only when WaitForKeybord is off", _autoForwardDelay, 4.0);

  registerOptionalParameter("AlignmentConstantName",
                            "Alignment constant from the condition file",
                            _alignmentCollectionName, string ("alignment"));



}

void EUTelEventViewer::init() {

  printParameters();

  MarlinCED::init(this) ;

}


void EUTelEventViewer::processRunHeader( LCRunHeader * rdr ) {
  std::unique_ptr<EUTelRunHeaderImpl> runHeader = std::make_unique<EUTelRunHeaderImpl>(rdr);
  runHeader->addProcessor(type());
}

void EUTelEventViewer::processEvent( LCEvent * evt ) {


  EUTelEventImpl * event = static_cast<EUTelEventImpl *> ( evt );

  if ( event->getEventType() == kEORE ) {
    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return;
  } else if ( event->getEventType() == kUNKNOWN ) {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }

  // for convenience
  int tempDetModel = _detModel;

  // try to load the alignment collections
  LCCollectionVec * alignmentCollectionVec = NULL;
  try {

    alignmentCollectionVec     = dynamic_cast < LCCollectionVec * > (evt->getCollection(_alignmentCollectionName));

  } catch ( DataNotAvailableException& e ) {

    if ( _detModel == -1 ) {
      // the user wants to draw the alignment corrected geometry, but
      // no alignment collection found
      streamlog_out ( WARNING2 ) << "No alignment condition found, using default ideal geometry model (99999)" << endl;
      tempDetModel = 99999;
    }
  }


  // Drawing Geometry
  MarlinCED::newEvent( this, tempDetModel ) ;

  // in the case the user wants to draw the alignment corrected
  // telescope planes ( tempDetModel == -1 ), then we have to draw
  // manually from here the geoboxes. For the time being only shifts
  // are applied. To apply rotations, one should rotated the box and
  // this is still not possible with CED.
  if ( tempDetModel == -1 ) {

    const gear::SiPlanesParameters&  siPlanesParameters  = Global::GEAR->getSiPlanesParameters();
    const gear::SiPlanesLayerLayout& siPlanesLayerLayout = siPlanesParameters.getSiPlanesLayerLayout();

    double * sizes  = new double[3];
    double * center = new double[3];
    unsigned int color = 0xFFFFFF;

    for ( int iLayer = 0 ; iLayer < siPlanesLayerLayout.getNLayers() ; iLayer++ ) {
      EUTelAlignmentConstant * alignment = static_cast< EUTelAlignmentConstant * >
        ( alignmentCollectionVec->getElementAt( iLayer ) );
      center[0] = siPlanesLayerLayout.getSensitivePositionX(iLayer) - alignment->getXOffset();
      center[1] = siPlanesLayerLayout.getSensitivePositionY(iLayer) - alignment->getYOffset();
      center[2] = siPlanesLayerLayout.getSensitivePositionZ(iLayer) - alignment->getZOffset();
      sizes[0]  = siPlanesLayerLayout.getSensitiveSizeX(iLayer);
      sizes[1]  = siPlanesLayerLayout.getSensitiveSizeY(iLayer);
      sizes[2]  = siPlanesLayerLayout.getSensitiveThickness(iLayer) ;
      ced_geobox( sizes, center, color );
    }
    delete [] center;
    delete [] sizes;

  }

  // Drawing hit collections
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
          unsigned int size  = 1;
          unsigned int type  = CED_HIT_STAR;
          ced_hit(x,y,z,type ,size , color);
        }
      }  catch (DataNotAvailableException& e) {
        streamlog_out ( WARNING2 ) << "No input collection (" << _trackerHitCollectionNameVec[iCollection] << " found on "
                                   <<  event->getEventNumber()
                                   << " in run " << event->getRunNumber() << endl;
      }
    }
  }

  // setup cellIdDecoder to decode the hit properties
  CellIDDecoder<TrackerHit>  hitCellDecoder(EUTELESCOPE::HITENCODING);

  // Drawing  Tracks
  if (_layerTrack >= 0) {
    for ( unsigned int iCollection = 0; iCollection < _trackCollectionNameVec.size(); iCollection++ ) {
      try {
        LCCollection * collection =    evt->getCollection(_trackCollectionNameVec[iCollection].c_str());
        for ( int iTrack = 0; iTrack < collection->getNumberOfElements(); iTrack++ ) {
          TrackImpl * track = dynamic_cast<TrackImpl*> ( collection->getElementAt(iTrack) );
          TrackerHitVec hitvec = track->getTrackerHits();

          unsigned int color =  returnColor(iCollection);

          // ok now I have everything I'm interested in. The idea
          // behing is that for each fitted hit I will draw a hit and
          // a line to the next plane.
          float xPrev = 0, yPrev = 0, zPrev = 0;
          float xNext, yNext, zNext;

          bool first = true;

          for ( size_t iHit = 0; iHit < hitvec.size()  ; ++iHit ) {

            TrackerHitImpl * hit;
            if ( (hit = dynamic_cast<TrackerHitImpl*> ( hitvec[ iHit ] )) != 0x0 ) {


              // show only hit resulting from fit
              if ( (hitCellDecoder(hit)["properties"] & kFittedHit) > 0 ) {
                if ( first ) {

                  first = false;
                  xPrev = static_cast< float > ( hit->getPosition()[0] );
                  yPrev = static_cast<float > ( hit->getPosition()[1] );
                  zPrev = static_cast<float > ( hit->getPosition()[2] );
                  ced_hit(xPrev, yPrev, zPrev, _layerTrackerHit << CED_LAYER_SHIFT,2,color);

                } else {

                  xNext = static_cast< float > ( hit->getPosition()[0] );
                  yNext = static_cast<float > ( hit->getPosition()[1] );
                  zNext = static_cast<float > ( hit->getPosition()[2] );
                  ced_hit(xNext, yNext, zNext, _layerTrackerHit << CED_LAYER_SHIFT,2,color);

                  ced_line( xPrev, yPrev, zPrev, xNext, yNext, zNext, _layerTrack << CED_LAYER_SHIFT, 2, color );

                  xPrev = xNext;
                  yPrev = yNext;
                  zPrev = zNext;


                }
              }

            }


          }
        }
      }  catch(DataNotAvailableException &e) {
        streamlog_out ( WARNING2 ) << "No input collection (" << _trackerHitCollectionNameVec[iCollection] << " found on "
                                   <<  event->getEventNumber()
                                   << " in run " << event->getRunNumber() << endl;
      }
    }
  }

  // draw the current event!
  MarlinCED::draw( this, _waitForKeyboard ) ;

  // in case we don't have to wait for the kepyboard, so wait
  // _autoForwardDelay second before continue
  if ( ! _waitForKeyboard ) usleep( static_cast< int > (_autoForwardDelay * 1000000) );

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


#endif // USE_GEAR && USE_CED
