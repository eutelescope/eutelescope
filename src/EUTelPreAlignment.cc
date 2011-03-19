#ifdef USE_GEAR
// eutelescope includes ".h"
#include "EUTelPreAlignment.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelAlignmentConstant.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"

// lcio includes <.h>
#include <IO/LCWriter.h>
#include <UTIL/LCTime.h>
#include <EVENT/LCCollection.h>
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

EUTelPreAlign::EUTelPreAlign () :Processor("EUTelPreAlign") {
  _description = "Apply alignment constants to hit collection";

  registerInputCollection (LCIO::TRACKERHIT, "InputHitCollectionName",
                           "The name of the input hit collection",
                           _inputHitCollectionName, string ("hit"));
  registerOptionalParameter ("FixedPlane", "SensorID of fixed plane", _fixedID, 0);
  registerOptionalParameter("AlignmentConstantLCIOFile","Name of LCIO db file where alignment constantds will be stored", 
			    _alignmentConstantLCIOFile, std::string( "alignment.slcio" ) );
}


void EUTelPreAlign::init () {
  // this method is called only once even when the rewind is active
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;  _iEvt = 0;
  
  _siPlanesParameters  = const_cast<SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));
  for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++ ) {
    if(_siPlanesLayerLayout->getID(iPlane) == _fixedID){ 
      //Get Zpos of ref plane
      _fixedZ = _siPlanesLayerLayout->getSensitivePositionZ(iPlane); 
    } else {
      //Get 
      _preAligners.push_back( PreAligner(  _siPlanesLayerLayout->getSensitivePitchX(iPlane),
					   _siPlanesLayerLayout->getSensitivePitchY(iPlane),
					   _siPlanesLayerLayout->getSensitivePositionZ(iPlane),
					   _siPlanesLayerLayout->getID(iPlane)) );
    }
  }
}

void EUTelPreAlign::processRunHeader (LCRunHeader * rdr) {
  auto_ptr<EUTelRunHeaderImpl> runHeader ( new EUTelRunHeaderImpl( rdr ) ) ;
  runHeader->addProcessor( type() );
  ++_iRun;
}

void EUTelPreAlign::processEvent (LCEvent * event) {
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
    //Loop over hits in fixed plane
    for (size_t ref = 0; ref < inputCollectionVec->size(); ref++) {
      TrackerHitImpl   * refHit   = dynamic_cast< TrackerHitImpl * >  ( inputCollectionVec->getElementAt( ref ) ) ;
      const double * refPos = refHit->getPosition();
      if( std::fabs(refPos[2] - _fixedZ) > 2.5) { continue; }
      
      for (size_t iHit = 0; iHit < inputCollectionVec->size(); iHit++) {
	TrackerHitImpl   * hit   = dynamic_cast< TrackerHitImpl * >  ( inputCollectionVec->getElementAt( iHit ) ) ;
	const double * pos = hit->getPosition();
	if( std::fabs(pos[2] - _fixedZ) < 2.5) { continue; }
	bool gotIt(false);
	for(size_t ii = 0; ii < _preAligners.size(); ii++){
	  PreAligner& pa = _preAligners.at(ii);
	  if( std::fabs(pa.getZPos() - pos[2]) > 2.5  ) { continue; }
	  gotIt = true;
	  pa.addPoint(refPos[0] - pos[0], refPos[1] - pos[1]);
	  break;
	}
	if(not gotIt) {
	  streamlog_out ( ERROR ) << "Mismatched hit at " << pos[2] << endl;
	}
      }
    }
  }
  catch (DataNotAvailableException& e) { 
    streamlog_out  ( WARNING2 ) <<  "No input collection " << _inputHitCollectionName << " found on event " << event->getEventNumber()
                                << " in run " << event->getRunNumber() << endl;
  }
}
void EUTelPreAlign::end() {
  LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();
  try {
    lcWriter->open( _alignmentConstantLCIOFile, LCIO::WRITE_NEW );
  } catch ( IOException& e ) {
    streamlog_out ( ERROR4 ) << e.what() << endl;
    exit(-1);
  }

  std::cout << "Writing to " << _alignmentConstantLCIOFile << std::endl;

  LCRunHeaderImpl * lcHeader  = new LCRunHeaderImpl;
  lcHeader->setRunNumber( 0 );
  lcWriter->writeRunHeader(lcHeader);
  delete lcHeader;
  LCEventImpl * event = new LCEventImpl;
  event->setRunNumber( 0 );
  event->setEventNumber( 0 );
  LCTime * now = new LCTime;
  event->setTimeStamp( now->timeStamp() );
  delete now;

  LCCollectionVec * constantsCollection = new LCCollectionVec( LCIO::LCGENERICOBJECT );
  for(size_t ii = 0 ; ii < _preAligners.size(); ii++){
    EUTelAlignmentConstant* constant = new EUTelAlignmentConstant();
    constant->setXOffset( -1.0 * _preAligners.at(ii).getPeakX());
    constant->setYOffset( -1.0 * _preAligners.at(ii).getPeakY());
    constant->setSensorID( _preAligners.at(ii).getIden() );
    constantsCollection->push_back( constant );
    streamlog_out ( MESSAGE ) << (*constant) << endl;
    // std::cout << "Iden: " << _preAligners.at(ii).getIden()
    // 	      << " Aligned x: "
    // 	      << _preAligners.at(ii).getPeakX()
    // 	      << " Aligned y:" 
    // 	      << _preAligners.at(ii).getPeakY() << std::endl;
  }
  event->addCollection( constantsCollection, "alignment" );
  lcWriter->writeEvent( event );
  delete event;
  lcWriter->close();
}
#endif
