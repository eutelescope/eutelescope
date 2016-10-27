// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelProcessorApplyAlignment.cc 2367 2013-02-12 15:41:08Z hperrey $
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
#include "EUTelProcessorApplyAlignment.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelAlignmentConstant.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

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

#include "TVector3.h"

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;
using namespace gear;

EUTelProcessorApplyAlign::EUTelProcessorApplyAlign () : Processor("EUTelProcessorApplyAlignment"),
_inputHitCollectionName("hit"),
_alignmentCollectionName("alignment"),
_outputHitCollectionName("correctedHit"),
_correctionMethod(0),
_iRun(0),
_iEvt(0),
_lookUpTable()
{
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


void EUTelProcessorApplyAlign::init () {
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

  _isFirstEvent = true;
}

void EUTelProcessorApplyAlign::processRunHeader (LCRunHeader * rdr) {
  std::unique_ptr<EUTelRunHeaderImpl> runHeader = std::make_unique<EUTelRunHeaderImpl>(rdr);
  runHeader->addProcessor( type() );
  ++_iRun;
}

void EUTelProcessorApplyAlign::processEvent (LCEvent * event) {
  ++_iEvt;
  
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
    
    if (_isFirstEvent) {
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
    }

    LCCollectionVec* outputCollectionVec = new LCCollectionVec(LCIO::TRACKERHIT);
    UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder( EUTELESCOPE::HITENCODING );      

    for (size_t iHit = 0; iHit < inputCollectionVec->size(); iHit++) {
      TrackerHitImpl   * inputHit   = dynamic_cast< TrackerHitImpl * >  ( inputCollectionVec->getElementAt( iHit ) ) ;
      // now we have to understand which layer this hit belongs to.

			UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder ( EUTELESCOPE::HITENCODING );
			const int sensorID = hitDecoder(inputHit)["sensorID"];

      // copy the input to the output, at least for the common part
      TrackerHitImpl   * outputHit  = new TrackerHitImpl;
      outputHit->setType( inputHit->getType() );
      outputHit->rawHits() = inputHit->getRawHits();
      FloatVec cov = inputHit->getCovMatrix();
      outputHit->setCovMatrix( cov );
      outputHit->setCellID0( inputHit->getCellID0() );
      outputHit->setCellID1( inputHit->getCellID1() );
      outputHit->setTime( inputHit->getTime() );

      // now that we know at which sensor the hit belongs to, we can
      // get the corresponding alignment constants
      map< int , int >::iterator  positionIter = _lookUpTable.find( sensorID );

      double * inputPosition      = const_cast< double * > ( inputHit->getPosition() ) ;
      double   outputPosition[3]  = { inputPosition[0],inputPosition[1], inputPosition[2] };

      if ( positionIter != _lookUpTable.end() ) {
        EUTelAlignmentConstant * alignment = static_cast< EUTelAlignmentConstant * >
          ( alignmentCollectionVec->getElementAt( positionIter->second ) );
	// Rotations
        double xPlaneCenter    = geo::gGeometry().siPlaneXPosition( sensorID );
        double yPlaneCenter    = geo::gGeometry().siPlaneYPosition( sensorID );
        double zPlaneThickness = geo::gGeometry().siPlaneZSize(sensorID) ;
        double zPlaneCenter    = geo::gGeometry().siPlaneZPosition( sensorID ) + zPlaneThickness / 2.;

        TVector3 inputVec( inputPosition[0] - xPlaneCenter, inputPosition[1] - yPlaneCenter, inputPosition[2] - zPlaneCenter );
        inputVec.RotateX( -alignment->getAlpha() );
        inputVec.RotateY( -alignment->getBeta() );
        inputVec.RotateZ( -alignment->getGamma() );
        

        outputPosition[0] = inputVec.X() + xPlaneCenter;
        outputPosition[1] = inputVec.Y() + yPlaneCenter;
        outputPosition[2] = inputVec.Z() + zPlaneCenter;
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

void EUTelProcessorApplyAlign::end() {
  streamlog_out ( MESSAGE2 ) <<  "Successfully finished" << endl;
}

#endif

