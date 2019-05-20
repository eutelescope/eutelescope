/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelHitCoordinateTransformer.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "CellIDReencoder.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"

// lcio includes <.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>

// system includes <>
#include <vector>

using namespace eutelescope;

EUTelHitCoordinateTransformer::EUTelHitCoordinateTransformer ():Processor ("EUTelHitCoordinateTransformer"),
_hitCollectionNameInput (), _hitCollectionNameOutput (),
_undoAlignment (false)
{

  _description =
    "EUTelHitCoordinateTransformer is responsible to change local "
    "coordinates to global using the EUTelGeometryClass.";

  registerInputCollection (LCIO::TRACKERHIT,
			   "hitCollectionNameInput",
			   "Input hit collection name",
			   _hitCollectionNameInput,
			   std::string ("local_hit"));

  registerOutputCollection (LCIO::TRACKERHIT,
			    "hitCollectionNameOutput",
			    "Output hit collection name",
			    _hitCollectionNameOutput,
			    std::string ("global_hit"));

  registerOptionalParameter ("UndoAlignment",
			     "Set to true to undo the alignment instead",
			     _undoAlignment, false);
}

void
EUTelHitCoordinateTransformer::init ()
{
  //initialize geometry 
  geo::gGeometry ().initializeTGeoDescription (EUTELESCOPE::GEOFILENAME,
					       EUTELESCOPE::DUMPGEOROOT);
}

void
EUTelHitCoordinateTransformer::processEvent (LCEvent * event)
{
  //check event type and if it is the last event
  EUTelEventImpl *evt = static_cast < EUTelEventImpl * >(event);
  if (evt->getEventType () == kEORE)
    {
      streamlog_out (MESSAGE5) << "EORE found: nothing else to do." << std::
	endl;
      return;
    }

  //opens collection for input
  LCCollection *inputCollection = nullptr;
  try
  {
    inputCollection = evt->getCollection (_hitCollectionNameInput);
  }
  catch (DataNotAvailableException & e)
  {
    streamlog_out (WARNING2) << _hitCollectionNameInput
      << " collection not available" << std::endl;
    return;
  }

  //opens collection for output
  LCCollectionVec *outputCollection = nullptr;
  try
  {
    outputCollection =
      static_cast <
      LCCollectionVec * >(event->getCollection (_hitCollectionNameOutput));
  }
  catch ( ...)
  {
    outputCollection = new LCCollectionVec (LCIO::TRACKERHIT);
  }

  //get encoding from input
  std::string encoding =
    inputCollection->getParameters ().getStringVal (LCIO::CellIDEncoding);
  if (encoding.empty ())
    {
      encoding = EUTELESCOPE::HITENCODING;
    }

  //get decoder
  lcio::CellIDDecoder < TrackerHitImpl > hitDecoder (encoding);
  lcio::UTIL::CellIDReencoder < TrackerHitImpl > cellReencoder (encoding,
								outputCollection);

  //[START] loop over hits
  for (int iHit = 0; iHit < inputCollection->getNumberOfElements (); ++iHit)
    {

      TrackerHitImpl *inputHit =
	static_cast <
	TrackerHitImpl * >(inputCollection->getElementAt (iHit));
      TrackerHitImpl *outputHit = new IMPL::TrackerHitImpl ();

      //get some basic information
      int properties = hitDecoder (inputHit)["properties"];
      int sensorID = hitDecoder (inputHit)["sensorID"];

      //use local2masterHit/master2localHit function in EUTelGeometryTelescopeDescription
      //to translate input/output position                
      const double *inputPos = inputHit->getPosition ();
      double outputPos[3];

      //transform local->global
      if (!(properties & kHitInGlobalCoord) && !_undoAlignment)
	{
	  streamlog_out (DEBUG0) << "Transforming hit from local to global!"
	    << std::endl;
	  geo::gGeometry ().local2Master (sensorID, inputPos, outputPos);
	}
      //transform global->local
      else if ((properties & kHitInGlobalCoord) && _undoAlignment)
	{
	  streamlog_out (DEBUG0) << "Transforming hit from global to local!"
	    << std::endl;
	  geo::gGeometry ().master2Local (sensorID, inputPos, outputPos);
	}
      //error case
      else
	{
	  std::cout << "Properties: " << properties << std::endl;
	  std::string errMsg;
	  if (!_undoAlignment)
	    {
	      errMsg =
		"Provided global hit, but trying to transform into global. Something is wrong!";
	    }
	  else
	    {
	      errMsg =
		"Provided local hit, but trying to transform into local. Something is wrong!";
	    }
	  throw InvalidGeometryException (errMsg);
	}

      //fill new outputHit with information
      outputHit->setPosition (outputPos);
      outputHit->setCovMatrix (inputHit->getCovMatrix ());
      outputHit->setType (inputHit->getType ());
      outputHit->setTime (inputHit->getTime ());
      outputHit->setCellID0 (inputHit->getCellID0 ());
      outputHit->setCellID1 (inputHit->getCellID1 ());
      outputHit->setQuality (inputHit->getQuality ());
      outputHit->rawHits () = inputHit->getRawHits ();

      //and reencode hit
      cellReencoder.readValues (outputHit);
      //^= is a bitwise XOR i.e. will switch the coordinate system
      cellReencoder["properties"] = properties ^= kHitInGlobalCoord;
      cellReencoder.setCellID (outputHit);
      //finally store it in collection
      outputCollection->push_back (outputHit);

    }				//[END] loop over hits

  //push the hit for this event onto the collection
  try
  {
    event->addCollection (outputCollection, _hitCollectionNameOutput);
  }
  catch ( ...)
  {
    streamlog_out (WARNING5) << "Problem with pushing collection onto event"
      << std::endl;
  }
}

void
EUTelHitCoordinateTransformer::end ()
{
  streamlog_out (MESSAGE4) << "Successfully finished" << std::endl;
}
