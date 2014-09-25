/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */



// personal includes ".h"
#include "ALIBAVA.h"
#include "AlibavaClusterCollectionMerger.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"

// eutelescope includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"

// marlin includes
#include "marlin/Global.h"
#include "marlin/Exceptions.h"
#include "marlin/Processor.h"
#include "marlin/DataSourceProcessor.h"
#include "marlin/ProcessorMgr.h"


// lcio includes
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/LCGenericObjectImpl.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/LCEventImpl.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>

// system includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <memory>
#include <stdlib.h>
#include <algorithm>

using namespace std;
using namespace marlin;
using namespace lcio;
using namespace alibava;
using namespace eutelescope;

AlibavaClusterCollectionMerger::AlibavaClusterCollectionMerger ():
DataSourceProcessor("AlibavaClusterCollectionMerger"),
//telescope
//_telescope_lcReader(),
_telescopeFileName(ALIBAVA::NOTSET),
_telescopePulseCollectionName(ALIBAVA::NOTSET),
_telescopeSparseCollectionName(ALIBAVA::NOTSET),
// alibava
//_alibava_lcReader(),
_alibavaFileName(ALIBAVA::NOTSET),
_alibavaPulseCollectionName(ALIBAVA::NOTSET),
_alibavaSparseCollectionName(ALIBAVA::NOTSET),
// output
_outputPulseCollectionName(ALIBAVA::NOTSET),
_outputSparseCollectionName(ALIBAVA::NOTSET)
{
	
	
	// initialize few variables
	
	_description = "Merges alibava and telescope cluster collections";
	
	///////////////
	// Telescope //
	///////////////
	
	// file name
	registerProcessorParameter("InputTelescopeFileName", "This is the input file name that telescope cluster collections stored",
										_telescopeFileName, string("runXXXXXX.slcio") );
	
	// pulse collection
	registerInputCollection (LCIO::TRACKERPULSE, "TelescopeClusterPulseCollectionName",
									 "Name of the cluster pulse collection of telescope data",
									 _telescopePulseCollectionName, string("telescope_cluster_pulse") );
	
	// sparse collection
	registerInputCollection (LCIO::TRACKERDATA, "TelescopeSparseClusterCollectionName",
									 "Name of the sparse cluster collection of telescope data",
									 _telescopeSparseCollectionName, string("telescope_sparse_cluster") );
	
	/////////////
	// Alibava //
	/////////////
	
	// file name
	registerProcessorParameter("InputAlibavaFileName", "This is the input file name that alibava cluster collections stored",
										_alibavaFileName, string("runXXXXXX.slcio") );
	
	// pulse collection
	registerInputCollection (LCIO::TRACKERPULSE, "AlibavaClusterPulseCollectionName",
									 "Name of the cluster pulse collection of alibava data",
									 _alibavaPulseCollectionName, string("alibava_cluster_pulse") );
	
	// sparse collection
	registerInputCollection (LCIO::TRACKERDATA, "AlibavaSparseClusterCollectionName",
									 "Name of the sparse cluster collection of alibava data",
									 _alibavaSparseCollectionName, string("alibava_sparse_cluster") );
	
	
	////////////
	// Output //
	////////////
	// pulse collection
	registerOutputCollection (LCIO::TRACKERPULSE, "OutputClusterPulseCollectionName",
									  "Name of the merged/output cluster pulse collection",
									  _outputPulseCollectionName, string("merged_cluster_pulse") );
	
	// sparse collection
	registerOutputCollection (LCIO::TRACKERDATA, "OutputSparseClusterCollectionName",
									  "Name of the merged/output sparse cluster collection. DO NOT Change this. This is hard coded in other  ",
									  _outputSparseCollectionName, string("original_zsdata") );
 
      registerProcessorParameter ("EventIDDifference",
                                                                                 "AlibavaEventNumber - TelescopeEventNumber",
                                                                                 _eventIDDiff , int(0));

 	
}

AlibavaClusterCollectionMerger * AlibavaClusterCollectionMerger::newProcessor () {
	return new AlibavaClusterCollectionMerger;
}



void AlibavaClusterCollectionMerger::init () {
	
	printParameters ();
}

void AlibavaClusterCollectionMerger::readDataSource(int /* numEvents */) {
	
	
	// open telescope file
	LCReader* telescope_lcReader = LCFactory::getInstance()->createLCReader();
	try {
		telescope_lcReader->open(_telescopeFileName );
	} catch( IOException& e ){
		streamlog_out ( ERROR1 ) << "Can't open the telescope file: " << e.what() << endl ;
	}
	
	// open alibava file
	LCReader* alibava_lcReader = LCFactory::getInstance()->createLCReader();
	try {
		alibava_lcReader->open( _alibavaFileName );
	} catch( IOException& e ){
		streamlog_out ( ERROR1 ) << "Can't open the alibava file: " << e.what() << endl ;
	}
	
	
	// we will copy alibava run header to as the header of output file.
	try {
		LCRunHeader* alibava_runHeader = alibava_lcReader->readNextRunHeader();
		ProcessorMgr::instance()->processRunHeader( alibava_runHeader ) ;
	} catch( IOException& e ){
		streamlog_out ( ERROR1 ) << "Can't access run header of the alibava file: " << e.what() << endl ;
	}
	
	int eventCounter=0;
	EUTelEventImpl*  telescopeEvent;
	AlibavaEventImpl* alibavaEvent;

	if (_eventIDDiff<0 ){
		for (int i=0; i<abs(_eventIDDiff); i++)
			telescopeEvent = static_cast<EUTelEventImpl*> ( telescope_lcReader->readNextEvent() );
	}
	else if (_eventIDDiff>0){
		for (int i=0; i<_eventIDDiff; i++)
			alibavaEvent = static_cast<AlibavaEventImpl*> ( alibava_lcReader->readNextEvent() );
	}

	while( ((telescopeEvent = static_cast<EUTelEventImpl*> ( telescope_lcReader->readNextEvent())) != 0 )
			&& ((alibavaEvent = static_cast<AlibavaEventImpl*> (alibava_lcReader->readNextEvent())) != 0 ) )
	{
		if (telescopeEvent->getEventType() == kEORE){ 
			streamlog_out ( MESSAGE5 ) << "Reached EORE of telescope data"<< endl;
			break;		
		}
		if ( eventCounter % 1000 == 0 )
			streamlog_out ( MESSAGE4 ) << "Looping events "<<alibavaEvent->getEventNumber() << endl;
		
		LCCollectionVec * alibavaPulseColVec = 0;
		LCCollectionVec * alibavaSparseColVec = 0;
		LCCollectionVec * telescopePulseColVec = 0;
		LCCollectionVec * telescopeSparseColVec = 0;
		try
		{
			// get alibava collections
			alibavaPulseColVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( _alibavaPulseCollectionName ) ) ;
			alibavaSparseColVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( _alibavaSparseCollectionName ) ) ;
			
			
		} catch ( IOException& e) {
			// do nothing again
			streamlog_out( ERROR5 ) << e.what() << endl;
		}
		
		try
		{
			// get telescope collections
			telescopePulseColVec = dynamic_cast< LCCollectionVec * > ( telescopeEvent->getCollection( _telescopePulseCollectionName ) ) ;
			telescopeSparseColVec = dynamic_cast< LCCollectionVec * > ( telescopeEvent->getCollection( _telescopeSparseCollectionName ) ) ;
			
		} catch ( IOException& e) {
			// do nothing again
			streamlog_out( ERROR5 ) << e.what() << endl;
		}
		
		// create output collections
		LCCollectionVec * outputPulseColVec = new LCCollectionVec(LCIO::TRACKERPULSE);
		LCCollectionVec * outputSparseColVec = new LCCollectionVec(LCIO::TRACKERDATA);
		
		// copy telescope clusters
		copyClustersInCollection(outputPulseColVec, outputSparseColVec, telescopePulseColVec, telescopeSparseColVec);
		// copy alibava cluster
		copyClustersInCollection(outputPulseColVec, outputSparseColVec, alibavaPulseColVec, alibavaSparseColVec);
		
		try
		{
			AlibavaEventImpl* outputEvent = new AlibavaEventImpl();
			
			outputEvent->setRunNumber( alibavaEvent->getRunNumber() );
			outputEvent->setEventNumber(eventCounter);
			outputEvent->setEventType( alibavaEvent->getEventType() );
			outputEvent->setEventSize( alibavaEvent->getEventSize() );
			outputEvent->setEventValue( alibavaEvent->getEventValue() );
			outputEvent->setEventTime( alibavaEvent->getEventTime() );
			outputEvent->setEventTemp( alibavaEvent->getEventTemp() );
			outputEvent->setCalCharge( alibavaEvent->getCalCharge() );
			outputEvent->setCalDelay( alibavaEvent->getCalDelay() );
			if (alibavaEvent->isEventMasked())
				outputEvent->maskEvent();
			else
				outputEvent->unmaskEvent();
						
			outputEvent->addCollection(outputPulseColVec, _outputPulseCollectionName);
			outputEvent->addCollection(outputSparseColVec, _outputSparseCollectionName);
			
			ProcessorMgr::instance()->processEvent( static_cast<LCEventImpl*> ( outputEvent ) ) ;
			// delete outputEvent;
			//	streamlog_out ( MESSAGE1 ) << "Successfully copied Alibava collections to output event" << endl ;
			
		} catch ( IOException& e) {
			// do nothing again
			streamlog_out( ERROR5 ) << e.what() << endl;
		}
		
      eventCounter++;
	}// end of loop over events
	
}


void AlibavaClusterCollectionMerger::copyClustersInCollection(LCCollectionVec * outputPulseColVec, LCCollectionVec * outputSparseColVec, LCCollectionVec * inputPulseColVec, LCCollectionVec * inputSparseColVec){
	
	// Here is the Cell ID Encodes for pulseFrame and sparseFrame
	// CellID Encodes are introduced in eutelescope::EUTELESCOPE
	
	// for sparseFrame (usually called cluster collection)
	CellIDEncoder<TrackerDataImpl> outputSparseColEncoder ( eutelescope::EUTELESCOPE::ZSCLUSTERDEFAULTENCODING, outputSparseColVec );
	// for pulseFrame
	CellIDEncoder<TrackerPulseImpl> outputPulseColEncoder ( eutelescope::EUTELESCOPE::PULSEDEFAULTENCODING, outputPulseColVec );
	
	// for input
	CellIDDecoder<TrackerPulseImpl> inputPulseColDecoder(inputPulseColVec);
	CellIDDecoder<TrackerDataImpl> inputSparseColDecoder(inputSparseColVec);
	
	unsigned int noOfClusters;
	// go through input clusters and copy them to output cluster collection
	
	noOfClusters = inputPulseColVec->getNumberOfElements();
	for ( size_t i = 0; i < noOfClusters; ++i ){
		TrackerPulseImpl * outputPulseFrame = new TrackerPulseImpl();
		TrackerDataImpl * outputSparseFrame = new TrackerDataImpl();
		
		TrackerPulseImpl* inputPulseFrame = dynamic_cast<TrackerPulseImpl*>(inputPulseColVec->getElementAt(i));
		TrackerDataImpl* inputSparseFrame = dynamic_cast<TrackerDataImpl*>(inputPulseFrame->getTrackerData());
		
		// set Cell ID for sparse collection
		outputSparseColEncoder["sensorID"] = static_cast<int>(inputSparseColDecoder(inputSparseFrame) ["sensorID"]);
		outputSparseColEncoder["sparsePixelType"] =static_cast<int>(inputSparseColDecoder(inputSparseFrame)["sparsePixelType"]);
		outputSparseColEncoder["quality"] = static_cast<int>(inputSparseColDecoder(inputSparseFrame)["quality"]);
		outputSparseColEncoder.setCellID( outputSparseFrame );
		
		// copy tracker data
		outputSparseFrame->setChargeValues(inputSparseFrame->getChargeValues());
		// add it to the cluster collection
		outputSparseColVec->push_back( outputSparseFrame );
		
		// prepare a pulse for this cluster
		outputPulseColEncoder["sensorID"] = static_cast<int> (inputPulseColDecoder(inputPulseFrame) ["sensorID"]);
		outputPulseColEncoder["type"] = static_cast<int>(inputPulseColDecoder(inputPulseFrame) ["type"]);
		outputPulseColEncoder.setCellID( outputPulseFrame );
		
		outputPulseFrame->setCharge( inputPulseFrame->getCharge() );
		outputPulseFrame->setTrackerData( outputSparseFrame);
		outputPulseColVec->push_back( outputPulseFrame );
		
	} // end of loop over input clusters
	
}


void AlibavaClusterCollectionMerger::end () {
	
	streamlog_out ( MESSAGE5 )  << "AlibavaClusterCollectionMerger Successfully finished" << endl;
}

