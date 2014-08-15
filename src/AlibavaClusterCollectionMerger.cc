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

// marlin includes
#include "marlin/Global.h"
#include "marlin/Exceptions.h"
#include "marlin/Processor.h"
#include "marlin/DataSourceProcessor.h"
#include "marlin/ProcessorMgr.h"


// lcio includes
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCGenericObjectImpl.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/LCEventImpl.h>
#include <UTIL/CellIDEncoder.h>

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
	
	
	LCEvent*  telescopeEvent;
	LCEvent* alibavaEvent;
	while( ((telescopeEvent = telescope_lcReader->readNextEvent()) != 0)
			&& ((alibavaEvent = alibava_lcReader->readNextEvent()) != 0 ) )
	{
		try
		{
			
			// get alibava collections
			LCCollectionVec * alibavaPulseColVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( _alibavaPulseCollectionName ) ) ;
			LCCollectionVec * alibavaSparseColVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( _alibavaSparseCollectionName ) ) ;
			
			// get telescope collections
			LCCollectionVec * telescopePulseColVec = dynamic_cast< LCCollectionVec * > ( telescopeEvent->getCollection( _telescopePulseCollectionName ) ) ;
			LCCollectionVec * telescopeSparseColVec = dynamic_cast< LCCollectionVec * > ( telescopeEvent->getCollection( _telescopeSparseCollectionName ) ) ;
			
			// create output collections
			LCCollectionVec * outputPulseColVec = new LCCollectionVec(LCIO::TRACKERPULSE);
			LCCollectionVec * outputSparseColVec = new LCCollectionVec(LCIO::TRACKERDATA);
			
			// Here is the Cell ID Encodes for pulseFrame and sparseFrame
			// CellID Encodes are introduced in eutelescope::EUTELESCOPE
			
			// for sparseFrame (usually called cluster collection)
			CellIDEncoder<TrackerDataImpl> sparseColEncoder ( eutelescope::EUTELESCOPE::ZSCLUSTERDEFAULTENCODING, outputSparseColVec );
			// for pulseFrame
			CellIDEncoder<TrackerPulseImpl> pulseColEncoder ( eutelescope::EUTELESCOPE::PULSEDEFAULTENCODING, outputPulseColVec );
			
			
			unsigned int noOfClusters;
			
			// first go through telescope clusters and copy them to output cluster collection
			// for sure pulse and sparse collections have same number of clusters
			noOfClusters = telescopeSparseColVec->getNumberOfElements();
			
			for ( size_t i = 0; i < noOfClusters; ++i ){
				outputSparseColVec->push_back( telescopeSparseColVec->getElementAt(i) );
				outputPulseColVec->push_back( telescopePulseColVec->getElementAt(i) );
			} // end of loop over telescope clusters
			
			// now alibava clusters, copy them to output cluster collection
			// for sure pulse and sparse collections have same number of clusters
			noOfClusters = alibavaSparseColVec->getNumberOfElements();
			
			for ( size_t i = 0; i < noOfClusters; ++i ){
				outputSparseColVec->push_back( alibavaSparseColVec->getElementAt(i) );
				outputPulseColVec->push_back( alibavaPulseColVec->getElementAt(i) );
			} // end of loop over alibava clusters
			
			
			LCEvent* outputEvent = new LCEvent();
			outputEvent->addCollection(outputPulseColVec, _outputPulseCollectionName);
			outputEvent->addCollection(outputSparseColVec, _outputSparseCollectionName);

			ProcessorMgr::instance()->processEvent( static_cast<LCEventImpl*> ( outputEvent ) ) ;
			delete outputEvent;
			
		} catch ( IOException& e) {
			// do nothing again
			streamlog_out( ERROR5 ) << e.what() << endl;
		}
		
      
	}
	
	
	/*
	 ///////////////////
	 // Process Event //
	 ///////////////////
	 
	 
	 // now write these to AlibavaEvent
	 AlibavaEventImpl* anEvent = new AlibavaEventImpl();
	 anEvent->setRunNumber(_runNumber);
	 anEvent->setEventNumber(eventCounter);
	 anEvent->setEventType(eventTypeCode);
	 anEvent->setEventSize(eventSize);
	 anEvent->setEventValue(value);
	 anEvent->setEventTime(tdc_time(tdcTime));
	 anEvent->setEventTemp(get_temperature(temp));
	 anEvent->setCalCharge(charge);
	 anEvent->setCalDelay(delay);
	 anEvent->unmaskEvent();
	 
	 
	 // creating LCCollection
	 LCCollectionVec* rawDataCollection = new LCCollectionVec(LCIO::TRACKERDATA);
	 CellIDEncoder<TrackerDataImpl> chipIDEncoder(ALIBAVA::ALIBAVADATA_ENCODE,rawDataCollection);
	 
	 // for this to work the _chipselection has to be sorted in ascending order!!!
	 for (unsigned int ichip=0; ichip<_chipSelection.size(); ichip++) {
	 FloatVec chipdata;
	 chipdata.clear();
	 
	 chipdata.insert(chipdata.end(), all_data.begin()+_chipSelection[ichip]*ALIBAVA::NOOFCHANNELS, all_data.begin()+(_chipSelection[ichip]+1)*ALIBAVA::NOOFCHANNELS);
	 TrackerDataImpl * arawdata = new TrackerDataImpl();
	 arawdata->setChargeValues(chipdata);
	 chipIDEncoder[ALIBAVA::ALIBAVADATA_ENCODE_CHIPNUM] = _chipSelection[ichip];
	 chipIDEncoder.setCellID(arawdata);
	 rawDataCollection->push_back(arawdata);
	 }
	 
	 anEvent->addCollection(rawDataCollection, _rawDataCollectionName);
	 
	 
	 if (_startEventNum!=-1 && eventCounter<_startEventNum) {
	 streamlog_out( MESSAGE5 )<<" Skipping event "<<eventCounter<<". StartEventNum is set to "<<_startEventNum<<endl;
	 eventCounter++;
	 continue;
	 }
	 
	 if (_stopEventNum!=-1 && eventCounter>_stopEventNum) {
	 streamlog_out( MESSAGE5 )<<" Reached StopEventNum: "<<_stopEventNum<<". Last saved event number is "<<eventCounter<<endl;
	 break;
	 }
	 
	 ProcessorMgr::instance()->processEvent( static_cast<LCEventImpl*> ( anEvent ) ) ;
	 eventCounter++;
	 
	 delete anEvent;
	 */
	
	
}


void AlibavaClusterCollectionMerger::end () {
	
	streamlog_out ( MESSAGE5 )  << "AlibavaClusterCollectionMerger Successfully finished" << endl;
}

