//Writen by Alexander Morton
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#include "EUTelProcessorCoordinateTransformHits.h"

// eutelescope includes ".h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"

//Standard C++ libraries 
#include <vector>

// lcio includes <.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>


using namespace std;
using namespace marlin;
using namespace eutelescope;

EUTelProcessorCoordinateTransformHits::EUTelProcessorCoordinateTransformHits():
Processor("EUTelProcessorCoordinateTransformHits"),
_hitCollectionNameInput(), 
_hitCollectionNameOutput()
{
		// modify processor description
		_description ="EUTelLocaltoGlobalHitMaker is responsible to change local coordinates to global. This is done using the EUTelGeometryClass";

		registerInputCollection(LCIO::TRACKERHIT, "hitCollectionNameInput", "Local input hit collection name", _hitCollectionNameInput, std::string("local_hit"));

		registerOutputCollection(LCIO::TRACKERHIT,"hitCollectionNameOutput", "Global output hit collection name", _hitCollectionNameOutput, std::string ("global_hit"));
}

void EUTelProcessorCoordinateTransformHits::init()
{
		std::string name("telescope_geometry.root");
		geo::gGeometry().initializeTGeoDescription(name,false);
}

void EUTelProcessorCoordinateTransformHits::processRunHeader(LCRunHeader* rdr)
{
		auto_ptr<EUTelRunHeaderImpl> header( new EUTelRunHeaderImpl(rdr) );

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
		LCCollection* inputCollection;
		try
		{
				inputCollection = evt->getCollection(_hitCollectionNameInput);
		}
		catch (DataNotAvailableException& e)
		{
				streamlog_out( WARNING2 ) << _hitCollectionNameInput << " collection not available" << std::endl;
				return;
		}

		LCCollectionVec * outputCollection = NULL;
		try{
				hitCollectionOutput  = static_cast<LCCollectionVec*> (event->getCollection( _hitCollectionNameOutput ));
		}
		catch(...){
				hitCollectionOutput = new LCCollectionVec(LCIO::TRACKERHIT);
				streamlog_out ( DEBUG5 ) << "Collection does not exist. Create new collection."<<endl;
		}
		UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder ( EUTELESCOPE::HITENCODING );   
		//Now get each individual hit LOOP OVER!
		for (int iHit = 0; iHit < collection->getNumberOfElements(); ++iHit) {  
				TrackerHitImpl*	hit_input = static_cast<TrackerHitImpl*>(collection->getElementAt(iHit)); //This will return a LCObject. Must cast to specify which type
				TrackerHitImpl* hit_output = new IMPL::TrackerHitImpl; 
				//Call the local2masterHit/master2localHit function defined int EUTelGeometryTelescopeDescription
				int properties = hitDecoder(static_cast< IMPL::TrackerHitImpl* >(hit_input))["properties"];
				if(properties == kHitInGlobalCoord){
						streamlog_out(DEBUG5) << " The properties cell ID is global. So will now change to local" << std::endl;
						geo::gGeometry().master2localHit(hit_input, hit_output, hitCollectionOutput);
				}
				else{
						streamlog_out(DEBUG5) << " The properties cell ID is not set so assume local. So will change to global now" << std::endl;
						geo::gGeometry().local2masterHit(hit_input, hit_output, hitCollectionOutput);
				}
				streamlog_out ( DEBUG5 )  << "New hit "<< iHit << " for event  "<< evt->getEventNumber() <<" created" << std::endl;

				try{
						streamlog_out ( DEBUG5 )  << "HIT OUTPUT CONTAINS!!!!!!!!!!!!!!!!!!!!" << std::endl;
						streamlog_out ( DEBUG5 )  << "Hit position: "<<*(hit_output->getPosition())<<" , " << *(hit_output->getPosition()+1)<<" , " << *(hit_output->getPosition()+2) << std::endl;
						streamlog_out ( DEBUG5 )  << "Type: "<< 	hit_output->getType() <<endl;
						UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder ( EUTELESCOPE::HITENCODING );
						streamlog_out ( DEBUG5 )  << "Sensor ID: " << 	hitDecoder(static_cast< IMPL::TrackerHitImpl* >(hit_output))["sensorID"] <<endl;
						streamlog_out ( DEBUG5 )  << "Properties: " << hitDecoder(static_cast< IMPL::TrackerHitImpl* >(hit_output))["properties"] <<endl;
				}
				catch(...){
						streamlog_out ( DEBUG5 )  << "There is a problem opening the new hit_output object" << std::endl;
				}
				hitCollectionOutput->push_back(hit_output);
		}
		streamlog_out ( DEBUG5 )  << "ALL HITS ON COLLECTIONVEC: Now for event "<< evt->getEventNumber() <<" push onto collection" << std::endl;
		//Now push the hit for this event onto the collection
		try{	
				event->addCollection(hitCollectionOutput, _hitCollectionNameOutput );
				streamlog_out ( DEBUG5 )  << "Pushed onto collection: " << _hitCollectionNameOutput <<endl;	
		}
		catch(...){
				streamlog_out ( MESSAGE5 )  << "Problem with pushing collection onto event"<<endl;
		}

}

void EUTelProcessorCoordinateTransformHits::end()
{
	streamlog_out(MESSAGE4) << "Successfully finished" << std::endl;
}
