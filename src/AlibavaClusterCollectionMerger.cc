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
	streamlog_out ( MESSAGE1 ) << "Successfully opened Telescope file" << endl ;
	
	// open alibava file
	LCReader* alibava_lcReader = LCFactory::getInstance()->createLCReader();
	try {
		alibava_lcReader->open( _alibavaFileName );
	} catch( IOException& e ){
		streamlog_out ( ERROR1 ) << "Can't open the alibava file: " << e.what() << endl ;
	}
	streamlog_out ( MESSAGE1 ) << "Successfully opened Alibava file" << endl ;
	
	
	// we will copy alibava run header to as the header of output file.
	try {
		LCRunHeader* alibava_runHeader = alibava_lcReader->readNextRunHeader();
		ProcessorMgr::instance()->processRunHeader( alibava_runHeader ) ;
	} catch( IOException& e ){
		streamlog_out ( ERROR1 ) << "Can't access run header of the alibava file: " << e.what() << endl ;
	}
	streamlog_out ( MESSAGE1 ) << "Successfully copied Alibava run header to output run header" << endl ;
	
	
	LCEvent*  telescopeEvent;
	LCEvent* alibavaEvent;
	while( ((telescopeEvent = telescope_lcReader->readNextEvent()) != 0)
			&& ((alibavaEvent = alibava_lcReader->readNextEvent()) != 0 ) )
	{
		if ( alibavaEvent->getEventNumber() % 1000 == 0 )
			streamlog_out ( MESSAGE4 ) << "Looping events "<<alibavaEvent->getEventNumber() << endl;
		
		LCCollectionVec * alibavaPulseColVec;
		LCCollectionVec * alibavaSparseColVec;
		LCCollectionVec * telescopePulseColVec;
		LCCollectionVec * telescopeSparseColVec;
		try
		{
			// get alibava collections
			alibavaPulseColVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( _alibavaPulseCollectionName ) ) ;
			alibavaSparseColVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( _alibavaSparseCollectionName ) ) ;
			
			streamlog_out ( DEBUG1 ) << "Alibava collections successfully found" << endl ;
			
		} catch ( IOException& e) {
			// do nothing again
			streamlog_out( ERROR5 ) << e.what() << endl;
		}
		
		try
		{
			// get telescope collections
			telescopePulseColVec = dynamic_cast< LCCollectionVec * > ( telescopeEvent->getCollection( _telescopePulseCollectionName ) ) ;
			telescopeSparseColVec = dynamic_cast< LCCollectionVec * > ( telescopeEvent->getCollection( _telescopeSparseCollectionName ) ) ;
			
			streamlog_out ( DEBUG1 ) << "Telescope collections successfully found" << endl ;
			
		} catch ( IOException& e) {
			// do nothing again
			streamlog_out( ERROR5 ) << e.what() << endl;
		}
		
		// create output collections
		LCCollectionVec * outputPulseColVec = new LCCollectionVec(LCIO::TRACKERPULSE);
		LCCollectionVec * outputSparseColVec = new LCCollectionVec(LCIO::TRACKERDATA);
		
		// Here is the Cell ID Encodes for pulseFrame and sparseFrame
		// CellID Encodes are introduced in eutelescope::EUTELESCOPE
		
		// for sparseFrame (usually called cluster collection)
		CellIDEncoder<TrackerDataImpl> outputSparseColEncoder ( eutelescope::EUTELESCOPE::ZSCLUSTERDEFAULTENCODING, outputSparseColVec );
		// for pulseFrame
		CellIDEncoder<TrackerPulseImpl> outputPulseColEncoder ( eutelescope::EUTELESCOPE::PULSEDEFAULTENCODING, outputPulseColVec );
		
		// for telescope
		CellIDDecoder<TrackerPulseImpl> telescopePulseColDecoder(telescopePulseColVec);
		CellIDDecoder<TrackerDataImpl> telescopeSparseColDecoder(telescopeSparseColVec);

		// for alibava
		CellIDDecoder<TrackerPulseImpl> alibavaPulseColDecoder(alibavaPulseColVec);
		CellIDDecoder<TrackerDataImpl> alibavaSparseColDecoder(alibavaSparseColVec);
		
		unsigned int noOfClusters;
		try
		{
			// first go through telescope clusters and copy them to output cluster collection
			// for sure pulse and sparse collections have same number of clusters
			noOfClusters = telescopePulseColVec->getNumberOfElements();
			for ( size_t i = 0; i < noOfClusters; ++i ){
				TrackerPulseImpl * outputPulseFrame = new TrackerPulseImpl();
				TrackerDataImpl * outputSparseFrame = new TrackerDataImpl();
				
				TrackerPulseImpl* telescopePulseFrame = dynamic_cast<TrackerPulseImpl*>(telescopePulseColVec->getElementAt(i));
				TrackerDataImpl* telescopeSparseFrame = dynamic_cast<TrackerDataImpl*>(telescopePulseFrame->getTrackerData());
				
				// set Cell ID for sparse collection
				outputSparseColEncoder["sensorID"] = static_cast<int>(telescopeSparseColDecoder(telescopeSparseFrame) ["sensorID"]);
				outputSparseColEncoder["sparsePixelType"] =static_cast<int>(telescopeSparseColDecoder(telescopeSparseFrame)["sparsePixelType"]);
				outputSparseColEncoder["quality"] = static_cast<int>(telescopeSparseColDecoder(telescopeSparseFrame)["quality"]);
				outputSparseColEncoder.setCellID( outputSparseFrame );
				
				// copy tracker data
				outputSparseFrame->setChargeValues(telescopeSparseFrame->getChargeValues());
				// add it to the cluster collection
				outputSparseColVec->push_back( outputSparseFrame );
				
				// prepare a pulse for this cluster
				outputPulseColEncoder["sensorID"] = static_cast<int> (telescopePulseColDecoder(telescopePulseFrame) ["sensorID"]);
				outputPulseColEncoder["type"] = static_cast<int>(telescopePulseColDecoder(telescopePulseFrame) ["type"]);
				outputPulseColEncoder.setCellID( outputPulseFrame );
				
				outputPulseFrame->setCharge( telescopePulseFrame->getCharge() );
				outputPulseFrame->setTrackerData( outputSparseFrame);
				outputPulseColVec->push_back( outputPulseFrame );
				/*
				 // set the ID for this zsCluster
				 sparseColEncoder["sensorID"] = _sensorIDStartsFrom + chipnum;
				 sparseColEncoder["sparsePixelType"] = static_cast<int>( kEUTelGenericSparsePixel );
				 sparseColEncoder["quality"] = static_cast<int>(0);
				 sparseColEncoder.setCellID( sparseFrame );
				 
				 // add it to the cluster collection
				 sparseColVec->push_back( sparseFrame );
				 
				 // prepare a pulse for this cluster
				 pulseColEncoder["sensorID"] = _sensorIDStartsFrom + chipnum;
				 pulseColEncoder["type"] = static_cast<int>(kEUTelSparseClusterImpl);
				 pulseColEncoder.setCellID( pulseFrame );
				 
				 pulseFrame->setCharge( totalSignal ); // no need
				 pulseFrame->setTrackerData( sparseFrame );
				 pulseColVec->push_back( pulseFrame );
				 */
				
				//	outputSparseColVec->push_back( telescopeSparseColVec->getElementAt(i) );
				//	outputPulseColVec->push_back( telescopePulseColVec->getElementAt(i) );
				
				
				
			} // end of loop over telescope clusters
			streamlog_out ( MESSAGE1 ) << "Successfully copied Telescope collections to output event" << endl ;
			
			// now alibava clusters, copy them to output cluster collection
			// for sure pulse and sparse collections have same number of clusters
			noOfClusters = alibavaSparseColVec->getNumberOfElements();
		/*
			for ( size_t i = 0; i < noOfClusters; ++i ){
				outputSparseColVec->push_back( alibavaSparseColVec->getElementAt(i) );
				outputPulseColVec->push_back( alibavaPulseColVec->getElementAt(i) );
			} // end of loop over alibava clusters
		 */
			streamlog_out ( MESSAGE1 ) << "Successfully copied Alibava collections to output event" << endl ;
		} catch ( IOException& e) {
			// do nothing again
			streamlog_out( ERROR5 ) << e.what() << endl;
		}
		
		try
		{
			LCEventImpl* outputEvent = new LCEventImpl();
			outputEvent->addCollection(outputPulseColVec, _outputPulseCollectionName);
			outputEvent->addCollection(outputSparseColVec, _outputSparseCollectionName);
			
			ProcessorMgr::instance()->processEvent( static_cast<LCEventImpl*> ( outputEvent ) ) ;
			// delete outputEvent;
			streamlog_out ( MESSAGE1 ) << "Successfully copied Alibava collections to output event" << endl ;
			
		} catch ( IOException& e) {
			// do nothing again
			streamlog_out( ERROR5 ) << e.what() << endl;
		}
		
      
	}// end of loop over events
	
}


void AlibavaClusterCollectionMerger::end () {
	
	streamlog_out ( MESSAGE5 )  << "AlibavaClusterCollectionMerger Successfully finished" << endl;
}

