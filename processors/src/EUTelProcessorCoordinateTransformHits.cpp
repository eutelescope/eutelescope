//Writen by Alexander Morton
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelProcessorCoordinateTransformHits.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "CellIDReencoder.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"

//Standard C++ libraries 
#include <vector>

// lcio includes <.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>

using namespace eutelescope;

EUTelProcessorCoordinateTransformHits::EUTelProcessorCoordinateTransformHits():
Processor("EUTelProcessorCoordinateTransformHits"),
_hitCollectionNameInput(), 
_hitCollectionNameOutput(),
_undoAlignment(false)
{
		_description ="EUTelLocaltoGlobalHitMaker is responsible to change local coordinates to global. This is done using the EUTelGeometryClass";

		registerInputCollection(LCIO::TRACKERHIT, "hitCollectionNameInput", "Local input hit collection name", _hitCollectionNameInput, std::string("local_hit"));

		registerOutputCollection(LCIO::TRACKERHIT,"hitCollectionNameOutput", "Global output hit collection name", _hitCollectionNameOutput, std::string ("global_hit"));

		registerOptionalParameter("Undo Alignment (boolean)", "Set to true to undo the alignment instead", _undoAlignment, bool(false));
}

void EUTelProcessorCoordinateTransformHits::init() {
		geo::gGeometry().initializeTGeoDescription(EUTELESCOPE::GEOFILENAME, EUTELESCOPE::DUMPGEOROOT);
}

void EUTelProcessorCoordinateTransformHits::processRunHeader(LCRunHeader* rdr)
{
		std::unique_ptr<EUTelRunHeaderImpl> header = std::make_unique<EUTelRunHeaderImpl>(rdr);

		// this is the right place also to check the geometry ID. This is a
		// unique number identifying each different geometry used at the
		// beam test. The same number should be saved in the run header and
		// in the xml file. If the numbers are different, instead of barely
		// quitting ask the user what to do.
		if ( header->getGeoID() == 0 )
				streamlog_out ( WARNING0 )
						<<  "The geometry ID in the run header is set to zero." << std::endl
						<<  "This may mean that the GeoID parameter was not set" << std::endl;


		if ( (unsigned int)header->getGeoID() != geo::gGeometry().getSiPlanesLayoutID() ) {
				streamlog_out ( WARNING5 ) 
						<<  "Error during the geometry consistency check: " << std::endl
						<< "The run header says the GeoID is " << header->getGeoID() << std::endl
						<< "The GEAR description says is     " << geo::gGeometry().getSiPlanesLayoutID() << std::endl;
		}
}

void EUTelProcessorCoordinateTransformHits::processEvent(LCEvent* event)
{
		//Check the event type and if it is the last event.
		EUTelEventImpl* evt	= static_cast<EUTelEventImpl*>(event);				
		if( evt->getEventType() == kEORE )
		{
				streamlog_out( MESSAGE5 ) << "EORE found: nothing else to do." << std::endl;
				return;
		} 
		else if( evt->getEventType() == kUNKNOWN )
		{
				streamlog_out( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
						<< " is of unknown type. Continue considering it as a normal Data Event." << std::endl;
		}

		//Opens collection for input.
		LCCollection* inputCollection = nullptr;
		try
		{
				inputCollection = evt->getCollection(_hitCollectionNameInput);
		}
		catch (DataNotAvailableException& e)
		{
				streamlog_out( WARNING2 ) << _hitCollectionNameInput << " collection not available" << std::endl;
				return;
		}

		LCCollectionVec* outputCollection = nullptr;
		try
		{
				outputCollection  = static_cast<LCCollectionVec*> (event->getCollection( _hitCollectionNameOutput ));
		}
		catch(...)
		{
				outputCollection = new LCCollectionVec(LCIO::TRACKERHIT);
		}

		std::string encoding = inputCollection->getParameters().getStringVal( LCIO::CellIDEncoding );

		if(encoding.empty())
		{
			encoding = EUTELESCOPE::HITENCODING;
		}

		lcio::CellIDDecoder<TrackerHitImpl> hitDecoder ( encoding );
		lcio::UTIL::CellIDReencoder<TrackerHitImpl> cellReencoder( encoding, outputCollection );

		//Now get each individual hit LOOP OVER!
		for(int iHit = 0; iHit < inputCollection->getNumberOfElements(); ++iHit)
		{  
			TrackerHitImpl*	inputHit = static_cast<TrackerHitImpl*>(inputCollection->getElementAt(iHit));
			TrackerHitImpl* outputHit = new IMPL::TrackerHitImpl(); 

			//Call the local2masterHit/master2localHit function defined int EUTelGeometryTelescopeDescription
			int properties = hitDecoder(inputHit)["properties"];
			int sensorID = hitDecoder(inputHit)["sensorID"];
			
			const double* inputPos = inputHit->getPosition();
			double outputPos[3];

			if( !(properties & kHitInGlobalCoord) && !_undoAlignment )
			{
				streamlog_out(DEBUG5) << "Transforming hit from local to global!" << std::endl;
				geo::gGeometry().local2Master(sensorID, inputPos, outputPos);
			}
			else if( (properties & kHitInGlobalCoord) && _undoAlignment )
			{
				streamlog_out(DEBUG5) << "Transforming hit from global to local!" << std::endl;
				geo::gGeometry().master2Local(sensorID, inputPos, outputPos);
			}
			else
			{
				std::cout << "Properties: " << properties <<std::endl;
				std::string errMsg;
				if(!_undoAlignment) errMsg = "Provided global hit, but trying to transform into global. Something is wrong!";
				else errMsg = "Provided local hit, but trying to transform into local. Something is wrong!";
				throw InvalidGeometryException(errMsg);
			}
	
			//Fill the new outputHit with information
			outputHit->setPosition(outputPos);
			outputHit->setCovMatrix( inputHit->getCovMatrix());
			outputHit->setType( inputHit->getType() );
			outputHit->setTime( inputHit->getTime() );
			outputHit->setCellID0( inputHit->getCellID0() );
			outputHit->setCellID1( inputHit->getCellID1() );
			outputHit->setQuality( inputHit->getQuality() );
			outputHit->rawHits() = inputHit->getRawHits();

			cellReencoder.readValues(outputHit);
			//^= is a bitwise XOR i.e. we will switch the coordinate sytsem
			cellReencoder["properties"] = properties ^= kHitInGlobalCoord;
			cellReencoder.setCellID(outputHit);

			outputCollection->push_back(outputHit);
		}
	
		//Now push the hit for this event onto the collection
		try
		{	
				event->addCollection(outputCollection, _hitCollectionNameOutput );
		}
		catch(...)
		{
				streamlog_out ( WARNING5 )  << "Problem with pushing collection onto event" << std::endl;
		}
}

void EUTelProcessorCoordinateTransformHits::end()
{
	streamlog_out(MESSAGE4) << "Successfully finished" << std::endl;
}
