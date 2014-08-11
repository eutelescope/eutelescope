/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


// alibava includes ".h"
#include "AlibavaClusterConverter.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"
#include "AlibavaCluster.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <lcio.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/LCEventImpl.h>

// eutelescope includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelSparseClusterImpl.h"

// ROOT includes ".h"

// system includes <>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <memory>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;
using namespace eutelescope;


AlibavaClusterConverter::AlibavaClusterConverter () :
AlibavaBaseProcessor("AlibavaClusterConverter"),
_pulseCollectionName(ALIBAVA::NOTSET),
_sparseCollectionName(ALIBAVA::NOTSET),
_sensorIDStartsFrom(0),
_missingCorrdinateValue(0)
{
	
	// modify processor description
	_description =
	"AlibavaClusterConverter converts AlibavaClusters to EUTelSparseCluster and  :) ";
	
	
	// first of register the input collection
	registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
									 "Input alibava cluster collection name",
									 _inputCollectionName, string("alibava_clusters") );
	
	// if needed one can change these to optional parameters
	
	registerProcessorParameter ("OutputEUTelClusterPulseCollectionName",
										 "The collection name of cluster pulse.  This might be hardcoded in EUTelescope framework",
										 _pulseCollectionName , string("clustercollection") );
	
	registerProcessorParameter ("OutputEUTelSparseClusterCollectionName",
										 "The collection name of sparse cluster.  This might be hardcoded in EUTelescope framework",
										 _sparseCollectionName , string("original_zsdata") );
	
	
	registerProcessorParameter ("SensorIDStartsFrom",
										 "The sensor ID for the data. The actual sensorID will be stored as SensorIDStartsFrom + ChipNumber ",
										 _sensorIDStartsFrom, int(6) );
	
	
	// now the optional parameters
	registerOptionalParameter ("MissingCoordinateValue",
										"The value that should be stored in missing coordinate. This number has to be integer since it will be used as channel number of the missing coordinate",
										_missingCorrdinateValue, int(0) );
}


void AlibavaClusterConverter::init () {
	streamlog_out ( MESSAGE4 ) << "Running init" << endl;
	
	// this method is called only once even when the rewind is active
	// usually a good idea to
	printParameters ();
	
	/* To choose if processor should skip masked events
	 ex. Set the value to 0 for false, to 1 for true
	 */
	if (Global::parameters->isParameterSet(ALIBAVA::SKIPMASKEDEVENTS))
		_skipMaskedEvents = bool ( Global::parameters->getIntVal(ALIBAVA::SKIPMASKEDEVENTS) );
	else {
		streamlog_out ( MESSAGE4 ) << "The Global Parameter "<< ALIBAVA::SKIPMASKEDEVENTS <<" is not set! Masked events will be used!" << endl;
	}
	
	
}
void AlibavaClusterConverter::processRunHeader (LCRunHeader * rdr) {
	streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;
	
	// Add processor name to the runheader
	auto_ptr<AlibavaRunHeaderImpl> arunHeader ( new AlibavaRunHeaderImpl(rdr)) ;
	arunHeader->addProcessor(type());
	
	// get and set selected chips
	setChipSelection(arunHeader->getChipSelection());
	
	// if you want
	bookHistos();
	
	// set number of skipped events to zero (defined in AlibavaBaseProcessor)
	_numberOfSkippedEvents = 0;
	
}


void AlibavaClusterConverter::processEvent (LCEvent * anEvent) {
	
	if ( anEvent->getEventNumber() % 1000 == 0 )
		streamlog_out ( MESSAGE4 ) << "Looping events "<<anEvent->getEventNumber() << endl;
	
	AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
	
	if (_skipMaskedEvents && (alibavaEvent->isEventMasked()) ) {
		_numberOfSkippedEvents++;
		return;
	}
	
	LCCollectionVec * alibavaCluColVec; // alibava cluster collection vector
	
	// the pulse collection we need for EUTelescope
	LCCollectionVec * pulseColVec =new LCCollectionVec(LCIO::TRACKERPULSE);
	// the sparse collection needed for EUTelescope
	LCCollectionVec * sparseColVec  =  new LCCollectionVec(LCIO::TRACKERDATA);
	
	// Here is the Cell ID Encodes for pulseFrame and sparseFrame
	// CellID Encodes are introduced in eutelescope::EUTELESCOPE
	
	// for sparseFrame (usually called cluster collection)
	CellIDEncoder<TrackerDataImpl> sparseColEncoder ( EUTELESCOPE::ZSCLUSTERDEFAULTENCODING, sparseColVec );
	// for pulseFrame
	CellIDEncoder<TrackerPulseImpl> pulseColEncoder ( EUTELESCOPE::PULSEDEFAULTENCODING, pulseColVec );
	
	// HERE
	
	
	unsigned int noOfClusters;
	try
	{
		alibavaCluColVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( getInputCollectionName() ) ) ;
		noOfClusters = alibavaCluColVec->getNumberOfElements();
		
		for ( size_t i = 0; i < noOfClusters; ++i )
		{
			// get your data from the collection and do what ever you want
			TrackerDataImpl * alibavaClu = dynamic_cast< TrackerDataImpl * > ( alibavaCluColVec->getElementAt( i ) ) ;
			
			AlibavaCluster anAlibavaCluster(alibavaClu);
			int chipnum = anAlibavaCluster.getChipNum();
			
			// For each cluster we will have pulseFrame and sparseFrame
			lcio::TrackerPulseImpl * pulseFrame = new lcio::TrackerPulseImpl();
			lcio::TrackerDataImpl * sparseFrame = new lcio::TrackerDataImpl();
			
			// Form EUTelSparseClusterImpl
			EUTelSparseClusterImpl< EUTelGenericSparsePixel > *eutelPixelCluster = new EUTelSparseClusterImpl< eutelescope::EUTelGenericSparsePixel >(sparseFrame);
			
			// for each member (channel) in the alibava cluster create a EUTelSimpleSparsePixel
			// and add this EUTelGenericSparsePixel to the EUTelSparseClusterImpl
			int clusterSize = anAlibavaCluster.getClusterSize();
			for (int imember=0; imember < clusterSize; imember++) {
				EUTelGenericSparsePixel eutelPixel;
				
				int memberChanNum = anAlibavaCluster.getChanNum(imember);
				float memberSignal = anAlibavaCluster.getSignal(imember) * anAlibavaCluster.getSignalPolarity();
				
				if (anAlibavaCluster.getIsSensitiveAxisX() == true) {
					eutelPixel.setXCoord( memberChanNum );
					eutelPixel.setYCoord( _missingCorrdinateValue );
				}
				else {
					eutelPixel.setXCoord( _missingCorrdinateValue );
					eutelPixel.setYCoord( memberChanNum );
				}
				eutelPixel.setSignal( memberSignal ); // EUTelGenericSparsePixel::setSignal() gets float
				eutelPixel.setTime(0); // there is no time info for channels in Alibava
				eutelPixelCluster->addSparsePixel( new EUTelGenericSparsePixel(eutelPixel) );
			}
			
			// Now we have a EUTelSparseCluster

			// Fill pulse collection
			float totalSignal = anAlibavaCluster.getTotalSignal() *anAlibavaCluster.getSignalPolarity();
			
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
			
		} // end of loop over number of clusters
		
		alibavaEvent->addCollection(pulseColVec, _pulseCollectionName);
		alibavaEvent->addCollection(sparseColVec, _sparseCollectionName);
		
	} catch ( lcio::DataNotAvailableException ) {
		// do nothing again
		streamlog_out( ERROR5 ) << "Collection ("<<getInputCollectionName()<<") not found! " << endl;
	}
	
}

void AlibavaClusterConverter::check (LCEvent * /* evt */ ) {
	// nothing to check here - could be used to fill check plots in reconstruction processor
}


void AlibavaClusterConverter::end() {
	
	if (_numberOfSkippedEvents > 0)
		streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents<<" events skipped since they are masked" << endl;
	streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
	
}

void AlibavaClusterConverter::fillHistos(TrackerDataImpl * /* trkdata */){
	// nothing is done here
}


void AlibavaClusterConverter::bookHistos(){
	// nothing is done here
	//	streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << endl;
}


