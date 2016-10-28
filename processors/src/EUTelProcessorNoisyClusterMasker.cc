/*
 *   This processor masks noisy clustes, i.e. clusters containing 
 *   hot pixels as specified by the hot pixel collection.
 *
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelProcessorNoisyClusterMasker.h"
#include "EUTELESCOPE.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelTrackerDataInterfacerImpl.h"
#include "CellIDReencoder.h"
#include "EUTelUtility.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <LCIOTypes.h>
#include <IO/LCWriter.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <Exceptions.h>

// system includes
#include <memory>
#include <algorithm>

namespace eutelescope {

bool EUTelProcessorNoisyClusterMasker::_staticPrintedSummary = false;

EUTelProcessorNoisyClusterMasker::EUTelProcessorNoisyClusterMasker():
  Processor("EUTelProcessorNoisyClusterMasker"),
  _inputCollectionName(""),
  _iRun(0),
  _iEvt(0),
  _firstEvent(true),
  _dataFormatChecked(false),
  _wrongDataFormat(false)
{
  _description ="EUTelProcessorNoisyClusterMasker masks pulses which contain hot pixels. For this, the quality field of pulses is used to encode the kNoisyCluster enum provided by EUTelescope.";

  registerInputCollection (LCIO::TRACKERPULSE, "InputCollectionName", "Input of zero suppressed data, still containing hot pixels", _inputCollectionName, std::string("cluster") );

  registerOptionalParameter("HotPixelCollectionName", "Name of the hot pixel collection.",  _noisyPixelCollectionName, std::string("hotpixel"));

}

void EUTelProcessorNoisyClusterMasker::init () {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters();
  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;
}

void EUTelProcessorNoisyClusterMasker::processRunHeader(LCRunHeader* /*rdr*/) {
  // increment the run counter
  ++_iRun;
  // reset the event counter
  _iEvt = 0;
}

void EUTelProcessorNoisyClusterMasker::processEvent(LCEvent * event) {
	if(_firstEvent) {
		//The noisy pixel collection stores all thot pixels in event #1
		//Thus we have to read it in in that case
		_noisyPixelMap = Utility::readNoisyPixelList(event, _noisyPixelCollectionName);
		_firstEvent = false;
	}

 	// get the collection of interest from the event.
	LCCollectionVec* pulseInputCollectionVec = nullptr;

	try {
    		pulseInputCollectionVec  = dynamic_cast <LCCollectionVec*>( event->getCollection(_inputCollectionName) );
	} catch( lcio::DataNotAvailableException& e ) {
		return;
  	}

	// prepare decoder for input data
	CellIDDecoder<TrackerPulseImpl> cellDecoder( pulseInputCollectionVec );
	
	//read the encoding string from the input collection
	std::string encoding = pulseInputCollectionVec->getParameters().getStringVal( LCIO::CellIDEncoding );
	//and the encoder for the data
	lcio::UTIL::CellIDReencoder<TrackerPulseImpl> cellReencoder( encoding, pulseInputCollectionVec );
	
	//loop over all the pulses
	for ( size_t iPulse = 0 ; iPulse < pulseInputCollectionVec->size(); iPulse++ ) {
		//the vector contains tracker pulses
        	TrackerPulseImpl* pulseData = dynamic_cast<TrackerPulseImpl*> ( pulseInputCollectionVec->getElementAt( iPulse ) );
		int sensorID = cellDecoder(pulseData)["sensorID"];		
	
	        //get the noise vector for the given plane
		std::vector<int>* noiseVector = &(_noisyPixelMap[sensorID]);
		
		//each pulse has the tracker data attached to it
		TrackerDataImpl* trackerData = dynamic_cast<TrackerDataImpl*>( pulseData->getTrackerData() );
		//decoder for tracker data
		CellIDDecoder<TrackerDataImpl> trackerDecoder ( EUTELESCOPE::ZSCLUSTERDEFAULTENCODING );
		int pixelType = trackerDecoder(trackerData)["sparsePixelType"];

		//interface to sparsified data
        auto sparseData = Utility::getSparseData(trackerData, pixelType);
		bool noisy = false;
		
		//Loop over all hits!
		for(auto& pixelRef: *sparseData) {
			auto& pixel = pixelRef.get();
			if(std::binary_search( noiseVector->begin(), noiseVector->end(), Utility::cantorEncode(pixel.getXCoord(), pixel.getYCoord()) )) {
				noisy = true;
				break;
			}
		}

		if(noisy) {
			int quality = cellDecoder(pulseData)["quality"];
			quality = quality | kNoisyCluster;
			//next line actually copies the old values
			cellReencoder.readValues(pulseData);
			//then we overwrite the "quality" one
			cellReencoder["quality"] = quality;
			//and apply the changes
			cellReencoder.setCellID(pulseData);
			_maskedNoisyClusters[sensorID]++;
		}	
	}
}

void EUTelProcessorNoisyClusterMasker::end() 
{
	//Print out some stats for the user
	if(!_staticPrintedSummary) {
		streamlog_out ( MESSAGE4 ) << "Noisy cluster masker(s) successfully finished" << std::endl;
		streamlog_out ( MESSAGE4 ) << "Printing summary:" << std::endl;
		_staticPrintedSummary = true;
	}
	for(std::map<int,int>::iterator it = _maskedNoisyClusters.begin(); it != _maskedNoisyClusters.end(); ++it) {
  		streamlog_out ( MESSAGE4 ) << "Masked " << (*it).second << " noisy clusters on plane " << (*it).first << "." << std::endl;
	}
}

} //namespace eutelescope
