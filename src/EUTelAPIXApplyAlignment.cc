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
 
#ifdef USE_GEAR
// eutelescope includes ".h"
#include "EUTelAPIXApplyAlignment.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelAlignmentConstant.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"

// lcio includes <.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes
#include <AIDA/IHistogram3D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ITree.h>
#include <AIDA/IHistogramFactory.h>
#include <marlin/AIDAProcessor.h>
#endif

// gear includes <.h>
#include "marlin/Global.h"
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// system includes <>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <memory>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;
using namespace gear;

EUTelAPIXApplyAlign::EUTelAPIXApplyAlign () :Processor("EUTelAPIXApplyAlignment") {
  _description = "Apply alignment constants to hit collection";

  registerInputCollection (LCIO::TRACKERHIT, "InputHitCollectionName",
                           "The name of the input hit collection",
                           _inputHitCollectionName, string ("hit"));
  registerInputCollection (LCIO::LCGENERICOBJECT, "AlignmentConstantName",
                           "Alignment constant from the condition file",
                           _alignmentCollectionName, string ("alignment"));
  registerOutputCollection (LCIO::TRACKERHIT, "OutputHitCollectionName",
                            "The name of the output hit collection",
                            _outputHitCollectionName, string("correctedHit"));
}


void EUTelAPIXApplyAlign::init () {
  // this method is called only once even when the rewind is active
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

  // check if the GEAR manager pointer is not null!
  if ( Global::GEAR == NULL ) {
    streamlog_out ( ERROR4 ) <<  "The GearMgr is not available, for an unknown reason." << endl;
    exit(-1);
  }

  _siPlanesParameters  = const_cast<SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

  _siPlaneZPosition = new double[ _siPlanesLayerLayout->getNLayers() ];
  for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++ ) {
    _siPlaneZPosition[ iPlane ] = _siPlanesLayerLayout->getLayerPositionZ(iPlane);
  }

  _isFirstEvent = true;
  fevent = true;
}

void EUTelAPIXApplyAlign::processRunHeader (LCRunHeader * rdr) {
  auto_ptr<EUTelRunHeaderImpl> runHeader ( new EUTelRunHeaderImpl( rdr ) ) ;
  runHeader->addProcessor( type() );
  ++_iRun;
}

void EUTelAPIXApplyAlign::processEvent (LCEvent * event) {
  ++_iEvt;
  if ( _iEvt % 10000 == 0 )
    streamlog_out ( MESSAGE4 ) << "Processing event "
                               << setw(6) << setiosflags(ios::right) << event->getEventNumber() << " in run "
                               << setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber()
                               << setfill(' ')
                               << " (Total = " << setw(10) << _iEvt << ")" << resetiosflags(ios::left) << endl;

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);

  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return;
  } else if ( evt->getEventType() == kUNKNOWN ) {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }

  try {
    LCCollectionVec * inputCollectionVec         = dynamic_cast < LCCollectionVec * > (evt->getCollection(_inputHitCollectionName));
    LCCollectionVec * alignmentCollectionVec     = dynamic_cast < LCCollectionVec * > (evt->getCollection(_alignmentCollectionName));
    
    if (fevent) {
      streamlog_out ( MESSAGE5 ) << "The alignment collection contains: " <<  alignmentCollectionVec->size()
                                << " planes " << endl;
      
      for ( size_t iPos = 0; iPos < alignmentCollectionVec->size(); ++iPos ) {
        EUTelAlignmentConstant * alignment = static_cast< EUTelAlignmentConstant * > ( alignmentCollectionVec->getElementAt( iPos ) );
        _lookUpTable[ alignment->getSensorID() ] = iPos;
      }

      map< int , int >::iterator mapIter = _lookUpTable.begin();
      for(;mapIter != _lookUpTable.end(); mapIter++){
        streamlog_out ( MESSAGE5 ) << "Sensor ID = " << mapIter->first
                                << " is in position " << mapIter->second << endl;
      }
      _isFirstEvent = false;
      fevent = false;
    }

    LCCollectionVec * outputCollectionVec = new LCCollectionVec(LCIO::TRACKERHIT);

    for (size_t iHit = 0; iHit < inputCollectionVec->size(); iHit++) {
      TrackerHitImpl   * inputHit   = dynamic_cast< TrackerHitImpl * >  ( inputCollectionVec->getElementAt( iHit ) ) ;
      // now we have to understand which layer this hit belongs to.
      int sensorID = guessSensorID( inputHit );

      // copy the input to the output, at least for the common part
      TrackerHitImpl   * outputHit  = new TrackerHitImpl;
      outputHit->setType( inputHit->getType() );
      outputHit->rawHits() = inputHit->getRawHits();

      // now that we know at which sensor the hit belongs to, we can
      // get the corresponding alignment constants
      map< int , int >::iterator  positionIter = _lookUpTable.find( sensorID );

      double * inputPosition      = const_cast< double * > ( inputHit->getPosition() ) ;
      double   outputPosition[3]  = { inputPosition[0],inputPosition[1], inputPosition[2] };

      if ( positionIter != _lookUpTable.end() ) {
        EUTelAlignmentConstant * alignment = static_cast< EUTelAlignmentConstant * >
          ( alignmentCollectionVec->getElementAt( positionIter->second ) );
	// Rotations
	outputPosition[0] = inputPosition[0] * (1.0 + alignment->getAlpha())  + alignment->getGamma() * inputPosition[1];
	outputPosition[1] = inputPosition[1] * (1.0 + alignment->getBeta())   - alignment->getGamma() * inputPosition[0];
	// Shifts
	outputPosition[0] -= alignment->getXOffset();
	outputPosition[1] -= alignment->getYOffset();
	outputPosition[2] -= alignment->getZOffset();

	streamlog_out( DEBUG5 ) << "position: " 
			       << inputPosition[0] << "," 
			       << inputPosition[1] << ","
			       << inputPosition[2] << " -> "
			       << outputPosition[0] << "," 
			       << outputPosition[1] << ","
			       << outputPosition[2] << endl;
      } else {
        // this hit belongs to a plane whose sensorID is not in the
        // alignment constants. So the idea is to eventually advice
        // the users if running in DEBUG and copy the not aligned hit
        // in the new collection.
        streamlog_out ( DEBUG5 ) << "Sensor ID " << sensorID << " not found. Skipping alignment for hit " << iHit << endl;
      }
      outputHit->setPosition( outputPosition );
      outputCollectionVec->push_back( outputHit );
    }
    evt->addCollection( outputCollectionVec, _outputHitCollectionName );
  } catch (DataNotAvailableException& e) {
    streamlog_out  ( WARNING2 ) <<  "No input collection " << _inputHitCollectionName << " found on event " << event->getEventNumber()
                                << " in run " << event->getRunNumber() << endl;
  }

}

void EUTelAPIXApplyAlign::end() {
  streamlog_out ( MESSAGE2 ) <<  "Successfully finished" << endl;
  delete [] _siPlaneZPosition;
}

int EUTelAPIXApplyAlign::guessSensorID( TrackerHitImpl * hit ) {
  int sensorID = -1;
  double minDistance =  numeric_limits< double >::max() ;
  double * hitPosition = const_cast<double * > (hit->getPosition());

  for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); ++iPlane ) {
    double distance = std::abs( hitPosition[2] - _siPlaneZPosition[ iPlane ] );
    if ( distance < minDistance ) {
      minDistance = distance;
      sensorID = _siPlanesLayerLayout->getID( iPlane );
    }
  }
  if  ( _siPlanesParameters->getSiPlanesType() == _siPlanesParameters->TelescopeWithDUT ) {
    double distance = std::abs( hitPosition[2] - _siPlanesLayerLayout->getDUTPositionZ() );
    if( distance < minDistance )
      {
        minDistance = distance;
        sensorID = _siPlanesLayerLayout->getDUTID();
      }
  }
  if ( minDistance > 5 /* mm */ ) {
    // advice the user that the guessing wasn't successful 
    streamlog_out( WARNING3 ) << "A hit was found " << minDistance << " mm far from the nearest plane\n"
      "Please check the consistency of the data with the GEAR file " << endl;
    throw SkipEventException(this);
  }
  return sensorID;
}


#endif
