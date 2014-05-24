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

#include "EUTelLocaltoGlobalHitMaker.h"

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


using namespace std;
using namespace marlin;
using namespace gear;
using namespace eutelescope;

//Define constructor. This will take collection name as input 
EUTelLocaltoGlobalHitMaker::EUTelLocaltoGlobalHitMaker () : Processor("EUTelLocaltoGlobalHitMaker"),
_hitCollectionNameInput(), //Variables here initialized with the constructor of this class
_hitCollectionNameOutput() //This variable is the output of this processor
{
	// modify processor description
	_description ="EUTelLocaltoGlobalHitMaker is responsible to change local coordinates to global. This is done using the GEAR geometry description";

  registerInputCollection(LCIO::TRACKERHIT,"hitCollectionNameInput",
                           "Hit collection name",
                           _hitCollectionNameInput, string ( "" ));

  registerOutputCollection(LCIO::TRACKERHIT,"hitCollectionNameOutput",   //Must be specified to Output so it can be read.
                           "Hit collection name",
                           _hitCollectionNameOutput, string ( "" ));

}

//This function is run at the start of the job
void EUTelLocaltoGlobalHitMaker::init() {
		//geo::gGeometry automatically calls the constructor of the class and then initializeTGeoDescription will create the TGeo object
    std::string name("test.C");
    geo::gGeometry().initializeTGeoDescription(name,false);
}


void EUTelLocaltoGlobalHitMaker::processRunHeader (LCRunHeader * rdr) {

  streamlog_out ( MESSAGE5 ) << "Initialize gear" << endl;

	//auto_ptr to ensure the object is deleted with the pointer after leaving block scope
  auto_ptr<EUTelRunHeaderImpl> header ( new EUTelRunHeaderImpl (rdr) );

  // this is the right place also to check the geometry ID. This is a
  // unique number identifying each different geometry used at the
  // beam test. The same number should be saved in the run header and
  // in the xml file. If the numbers are different, instead of barely
  // quitting ask the user what to do.

  if ( header->getGeoID() == 0 )
    streamlog_out ( WARNING0 ) <<  "The geometry ID in the run header is set to zero." << endl
                               <<  "This may mean that the GeoID parameter was not set" << endl;


  if ( header->getGeoID() != geo::gGeometry()._siPlanesParameters->getSiPlanesID() ) {
  	streamlog_out ( WARNING5 ) <<  "Error during the geometry consistency check: " << endl
                             << "The run header says the GeoID is " << header->getGeoID() << endl
                             << "The GEAR description says is     " << geo::gGeometry()._siPlanesParameters->getSiPlanesID() << endl;

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

void EUTelLocaltoGlobalHitMaker::check(LCEvent *event){
	cout << "Checking event!!!!!!!   "<<event->getEventNumber() << endl;
}



void EUTelLocaltoGlobalHitMaker::processEvent (LCEvent * event) {
	streamlog_out ( DEBUG5 ) << "Beginning event number " << event->getEventNumber() << " in run " << event->getRunNumber() <<std::endl;
	//Check if the event is the last or if the event is of a unknown type
	EUTelEventImpl * evt	= static_cast<EUTelEventImpl*> (event) ;				
	if ( evt->getEventType() == kEORE ) {
    streamlog_out ( MESSAGE5 ) << "EORE found: nothing else to do." << endl;
		return;
  } else if ( evt->getEventType() == kUNKNOWN ) {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }
	bool already_exists=false;
	LCCollectionVec * hitCollectionOutput = NULL;
	try{
		already_exists=true;
		hitCollectionOutput  = static_cast<LCCollectionVec*> (event->getCollection( _hitCollectionNameOutput ));
		streamlog_out ( MESSAGE5 ) << "Got collection: "<< _hitCollectionNameOutput<<endl;
	}
	catch(...){
		already_exists=false;
		hitCollectionOutput = new LCCollectionVec(LCIO::TRACKERHIT);
		streamlog_out ( MESSAGE5 ) << "Create Collection "<<endl;
	}

	//Opens collection for input. If it can not find it then through exception
	LCCollection* collection;
	try {
		collection = evt->getCollection(_hitCollectionNameInput);
	} catch (...) {
		streamlog_out(WARNING2) << _hitCollectionNameInput << " collection not available" << std::endl;
		return;
	}
	streamlog_out(DEBUG1) << "collection : " << _hitCollectionNameInput << " retrieved" << std::endl;
	//Create two pointers to the input and output hit
	TrackerHit* hit_input;
	TrackerHitImpl* hit_output_total = new IMPL::TrackerHitImpl[collection->getNumberOfElements()];    
	//Now get each individual hit LOOP OVER!
	for (int iHit = 0; iHit < collection->getNumberOfElements(); ++iHit) {  
		hit_input = static_cast<TrackerHit*>(collection->getElementAt(iHit)); //This wll return a LCObject. Must Cast to specify which object
		TrackerHitImpl* hit_output = hit_output_total+iHit;
		//Call the local2masterHit function defined int EUTelGeometryTelescopeDescription
		geo::gGeometry().local2masterHit(hit_input, hit_output, hitCollectionOutput);
		streamlog_out ( DEBUG5 )  << "New hit "<< iHit << " for event  "<< evt->getEventNumber() <<" created" << endl;
	
		try{
			streamlog_out ( DEBUG5 )  << "HIT OUTPUT CONTAINS!!!!!!!!!!!!!!!!!!!!" << endl;
			streamlog_out ( DEBUG5 )  << "Hit position: "<<*(hit_output->getPosition())<<" , " << *(hit_output->getPosition()+1)<<" , " << *(hit_output->getPosition()+2) << endl;
		//	streamlog_out ( DEBUG5 )  << "Matrix: "<<	*(hit_output->getCovMatrix()) <<endl;
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
	if(!already_exists){ 
		event->addCollection(hitCollectionOutput, _hitCollectionNameOutput );
		streamlog_out ( MESSAGE5 )  << "Pushed onto collection: " << _hitCollectionNameOutput <<endl;	
	}
	else 	streamlog_out ( MESSAGE5 )  << "The collection already exists. Did not add to event"<<endl;	
	}
	catch(...){
		streamlog_out ( MESSAGE5 )  << "Problem with pushing collection onto event"<<endl;
	}
	streamlog_out ( DEBUG5 )  << "Push Successful"<<endl;
	delete [] hit_output_total;
//	delete hitCollectionOutput;
	streamlog_out(DEBUG5) << "Delete successful" << endl;

}

void EUTelLocaltoGlobalHitMaker::end(){

	streamlog_out ( MESSAGE4 )  << "Successfully finished" << endl;

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////
////////////////////////////// HERE WE DEFINE THE FUNCTIONS THAT ARE USED IN THIS PROCESSOR.
/////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
