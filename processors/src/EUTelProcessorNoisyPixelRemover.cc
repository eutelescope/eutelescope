/*
 *   This processor removes noisy pixels from a TrackerData collection
 *   been masked as noisy and written out into a noisy pixel database
 *
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelProcessorNoisyPixelRemover.h"
#include "EUTelTrackerDataInterfacerImpl.h"
#include "EUTelUtility.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "EUTelRunHeaderImpl.h"

// lcio includes <.h>
#include <LCIOTypes.h>
#include <IO/LCWriter.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>

// system includes
#include <memory>
#include <algorithm>

namespace eutelescope {

EUTelProcessorNoisyPixelRemover::EUTelProcessorNoisyPixelRemover():
  Processor("EUTelProcessorNoisyPixelRemover"),
  _inputCollectionName(""),
  _outputCollectionName(""),
  _noisyPixelCollectionName("")
{
  _description ="EUTelProcessorNoisyPixelRemover removes noisy pixels (TrackerData) from a collection. This processor requires a noisy pixel collection.";

  registerInputCollection(LCIO::TRACKERDATA, "InputCollectionName", "Input collection containing noisy raw data", _inputCollectionName, std::string ("noisy_raw_data_collection"));
  registerOutputCollection(LCIO::TRACKERDATA, "OutputCollectionName", "Output collection where noisy pixels have been removed", _outputCollectionName, std::string("noisefree_raw_data_collection"));
  registerProcessorParameter("NoisyPixelCollectionName", "Name of the noisy pixel collection.",  _noisyPixelCollectionName, std::string("noisypixel"));
}

void EUTelProcessorNoisyPixelRemover::init() {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters();
}

void EUTelProcessorNoisyPixelRemover::processRunHeader(LCRunHeader* rdr){
	auto runHeader = std::make_unique<EUTelRunHeaderImpl>(rdr);
	runHeader->addProcessor(type());
}

void EUTelProcessorNoisyPixelRemover::processEvent(LCEvent* event) {

	if(_firstEvent) {
		//The noisy pixel collection stores all thot pixels in event #1
		//Thus we have to read it in in that case
		_noisyPixelMap = Utility::readNoisyPixelList(event, _noisyPixelCollectionName);
		_firstEvent = false;
	}

 	// get the collection of interest from the event.
	LCCollectionVec* inputCollection = nullptr;

	try {
    		inputCollection  = dynamic_cast<LCCollectionVec*>( event->getCollection(_inputCollectionName) );
	} catch( lcio::DataNotAvailableException& e ) {
		return;
  	}

	//now prepare output collection
	LCCollectionVec* outputCollection = nullptr;
	bool outputCollectionExists = false;
  	size_t initialOutputCollectionSize = 0;

  	try {
   		outputCollection = dynamic_cast<LCCollectionVec*>( event->getCollection( _outputCollectionName ) );
    		outputCollectionExists = true;
    		initialOutputCollectionSize = outputCollection->size();
  	} catch ( lcio::DataNotAvailableException& e ) {
    		outputCollection = new LCCollectionVec(LCIO::TRACKERDATA);
  	}
	
	// prepare decoder for input data
	CellIDDecoder<TrackerDataImpl> inputDataDecoder( EUTELESCOPE::ZSDATADEFAULTENCODING );

 	//read the encoding std::string from the input collection
	std::string encodingString = inputCollection->getParameters().getStringVal( LCIO::CellIDEncoding );	
	outputCollection->parameters().setValue(LCIO::CellIDEncoding, encodingString);
	
	auto trackerData = std::make_unique<lcio::TrackerDataImpl>();
	
	for ( size_t iEntry = 0; iEntry < inputCollection->size(); ++iEntry ){

		if(iEntry > 0) {
				outputCollection->push_back( trackerData.release() );
				auto newTrackerData = std::make_unique<lcio::TrackerDataImpl>();
				trackerData = std::move( newTrackerData );	
		}

        	TrackerDataImpl* inputData = dynamic_cast<TrackerDataImpl*>( inputCollection->getElementAt(iEntry) );
		
		int sensorID = inputDataDecoder(inputData)["sensorID"];		
		SparsePixelType pixelType = static_cast<SparsePixelType>(static_cast<int>(inputDataDecoder(inputData)["sparsePixelType"]));

		trackerData->setCellID0( inputData->getCellID0() );
		trackerData->setCellID1( inputData->getCellID1() );
		trackerData->setTime( inputData->getTime() );
				
		//get the noise vector for the given plane
		std::vector<int>* noiseVector = &(_noisyPixelMap[sensorID]);
		
		//interface to sparsified data
		auto sparseDataInterface = Utility::getSparseData(inputData, pixelType);
		auto sparseOutputData = Utility::getSparseData(trackerData.get(), pixelType);

		for(auto& pixelRef: *sparseDataInterface) {
			auto& pixel = pixelRef.get();
			if(!std::binary_search( noiseVector->begin(), noiseVector->end(), Utility::cantorEncode(pixel.getXCoord(), pixel.getYCoord()) )) {
					sparseOutputData->push_back(pixel);
			}
		}
	}	
	outputCollection->push_back( trackerData.release() );
	
	//add the collection if we created it and added elements
	if ( !outputCollectionExists ) {      
		if( outputCollection->size() != initialOutputCollectionSize ) {
			event->addCollection( outputCollection, _outputCollectionName );
		} else {
			delete outputCollection;
		}
	}
}

void EUTelProcessorNoisyPixelRemover::end() {
}

}
