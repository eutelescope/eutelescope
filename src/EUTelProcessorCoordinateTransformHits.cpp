//Writen by Alexander Morton <alexander.morton@desy.de>
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// built only if GEAR is available
#ifdef USE_GEAR

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

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

//Only create histograms if one of these are available 
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram2D.h>
#endif

//Standard C++ libraries 
#include <vector>

// lcio includes <.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>


using namespace std;
using namespace marlin;
using namespace gear;
using namespace eutelescope;

//Define constructor. This will take collection name as input 
EUTelProcessorCoordinateTransformHits::EUTelProcessorCoordinateTransformHits () : Processor("EUTelProcessorCoordinateTransformHits"),
_hitCollectionNameInput(), //Variables here initialized with the constructor of this class
_hitCollectionNameOutput() //This variable is the output of this processor
{
	// modify processor description
	_description ="EUTelLocaltoGlobalHitMaker is responsible to change local coordinates to global. This is done using the EUTelGeometryClass";

  registerInputCollection(LCIO::TRACKERHIT,"hitCollectionNameInput",
                           "Hit collection name",
                           _hitCollectionNameInput, string ( "" ));

  registerOutputCollection(LCIO::TRACKERHIT,"hitCollectionNameOutput",   //Must be specified to output so it can be read.
                           "Hit collection name",
                           _hitCollectionNameOutput, string ( "" ));

}

//This function is run at the start of the job only
void EUTelProcessorCoordinateTransformHits::init() {
		//geo::gGeometry automatically calls the constructor of the class which will fill all variables used by the gear file. Then initialise with create the volumes of the telescope using TGeo from ROOT.
    std::string name("test.root");
    geo::gGeometry().initializeTGeoDescription(name,false);
}


void EUTelProcessorCoordinateTransformHits::processRunHeader (LCRunHeader * rdr) {

  streamlog_out ( MESSAGE5 ) << "Initialize gear" << endl;

  auto_ptr<EUTelRunHeaderImpl> header ( new EUTelRunHeaderImpl (rdr) );

  // this is the right place also to check the geometry ID. This is a
  // unique number identifying each different geometry used at the
  // beam test. The same number should be saved in the run header and
  // in the xml file. If the numbers are different, instead of barely
  // quitting ask the user what to do.

  if ( header->getGeoID() == 0 )
    streamlog_out ( WARNING0 ) <<  "The geometry ID in the run header is set to zero." << endl
                               <<  "This may mean that the GeoID parameter was not set" << endl;


  if ( (unsigned int)header->getGeoID() != geo::gGeometry().getSiPlanesLayoutID() ) {
  	streamlog_out ( WARNING5 ) <<  "Error during the geometry consistency check: " << endl
                             << "The run header says the GeoID is " << header->getGeoID() << endl
                             << "The GEAR description says is     " << geo::gGeometry().getSiPlanesLayoutID() << endl;

//Here if the geometry number and gear file number don't match decide if you should continue
#ifdef EUTEL_INTERACTIVE
    string answer;
    while (true) {
      streamlog_out ( ERROR1 ) << "Type Q to quit now or C to continue using the actual GEAR description anyway [Q/C]" << endl;
      cin >> answer;
      // put the answer in lower case before making the comparison.
      transform( answer.begin(), answer.end(), answer.begin(), ::tolower );
      if ( answer == "q" ) {
        exit(-1);
      } else if ( answer == "c" ) {
        break;
      }
    }//End of while statement
#endif
  }


}//end of processRunHeader

void EUTelProcessorCoordinateTransformHits::check(LCEvent *event){
}



void EUTelProcessorCoordinateTransformHits::processEvent (LCEvent * event) {
	streamlog_out ( DEBUG5 ) << "Beginning event number " << event->getEventNumber() << " in run " << event->getRunNumber() <<std::endl;
	//Check the event type and if it is the last event.
	EUTelEventImpl * evt	= static_cast<EUTelEventImpl*> (event) ;				
	if ( evt->getEventType() == kEORE ) {
    streamlog_out ( MESSAGE5 ) << "EORE found: nothing else to do." << endl;
		return;
  } else if ( evt->getEventType() == kUNKNOWN ) {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }


	LCCollectionVec * hitCollectionOutput = NULL;
	try{
		hitCollectionOutput  = static_cast<LCCollectionVec*> (event->getCollection( _hitCollectionNameOutput ));
	}
	catch(...){
		hitCollectionOutput = new LCCollectionVec(LCIO::TRACKERHIT);
		streamlog_out ( DEBUG5 ) << "Collection does not exist. Create new collection."<<endl;
	}

	//Opens collection for input. If it can not find it then through exception
	streamlog_out(DEBUG1) << "Try to get input collection: " << _hitCollectionNameInput << std::endl;
	LCCollection* collection;
	try {
		collection = evt->getCollection(_hitCollectionNameInput);
	} catch (...) {
		streamlog_out(WARNING2) << _hitCollectionNameInput << " collection not available" << std::endl;
		return;
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
		streamlog_out ( DEBUG5 )  << "New hit "<< iHit << " for event  "<< evt->getEventNumber() <<" created" << endl;
	
		try{
			streamlog_out ( DEBUG5 )  << "HIT OUTPUT CONTAINS!!!!!!!!!!!!!!!!!!!!" << endl;
			streamlog_out ( DEBUG5 )  << "Hit position: "<<*(hit_output->getPosition())<<" , " << *(hit_output->getPosition()+1)<<" , " << *(hit_output->getPosition()+2) << endl;
 			streamlog_out ( DEBUG5 )  << "Type: "<< 	hit_output->getType() <<endl;
			UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder ( EUTELESCOPE::HITENCODING );
			streamlog_out ( DEBUG5 )  << "Sensor ID: " << 	hitDecoder(static_cast< IMPL::TrackerHitImpl* >(hit_output))["sensorID"] <<endl;
			streamlog_out ( DEBUG5 )  << "Properties: " << hitDecoder(static_cast< IMPL::TrackerHitImpl* >(hit_output))["properties"] <<endl;
		}
		catch(...){
			streamlog_out ( DEBUG5 )  << "There is a problem opening the new hit_output object" << endl;
		}
		hitCollectionOutput->push_back(hit_output);
	}
	streamlog_out ( DEBUG5 )  << "ALL HITS ON COLLECTIONVEC: Now for event "<< evt->getEventNumber() <<" push onto collection" << endl;
	//Now push the hit for this event onto the collection
	try{	
		event->addCollection(hitCollectionOutput, _hitCollectionNameOutput );
		streamlog_out ( DEBUG5 )  << "Pushed onto collection: " << _hitCollectionNameOutput <<endl;	
	}
	catch(...){
		streamlog_out ( MESSAGE5 )  << "Problem with pushing collection onto event"<<endl;
	}

}

void EUTelProcessorCoordinateTransformHits::end(){

	streamlog_out ( MESSAGE4 )  << "Successfully finished" << endl;

}

#endif
