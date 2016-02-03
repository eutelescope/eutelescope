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

  registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName", "Input of zero suppressed data, still containing hot pixels", _inputCollectionName, std::string("cluster") );

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
		readNoisyPixelList(event);
		_firstEvent = false;
	}

 	// get the collection of interest from the event.
	LCCollectionVec* pulseInputCollectionVec = NULL;

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

		//IMPORTANT: if we pass a nullptr to EUTelBaseSparsePixel* getSparsePixelAt( int index, EUTelBaseSparsePixel* pixel)
		//where pixel=NULL, this method will return a pointer to a new pixel of derived type. YOU have to delete it!
		EUTelBaseSparsePixel* pixel = nullptr;

		//Loop over all hits!
		for ( unsigned int iPixel = 0; iPixel < sparseData->size(); iPixel++ ) {
			//IF pixel=NULL this yields a new pixel
		        pixel = sparseData->getSparsePixelAt( iPixel, pixel );
			if(std::binary_search( noiseVector->begin(), noiseVector->end(), encode(pixel->getXCoord(), pixel->getYCoord()) )) {
				noisy=true;
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
		delete pixel;
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

int EUTelProcessorNoisyClusterMasker::encode(int X, int Y) {
	//Cantor pairing function
	return static_cast<int>( 0.5*(X+Y)*(X+Y+1)+Y );
} 

void EUTelProcessorNoisyClusterMasker::readNoisyPixelList(LCEvent* event) {
	//Preapare pointer to hot pixel collection
	LCCollectionVec* noisyPixelCollectionVec = nullptr;

	//Try to obtain the collection
	try  {
		noisyPixelCollectionVec = static_cast< LCCollectionVec*>( event->getCollection(_noisyPixelCollectionName) );
	} catch (...) {
		if (!_noisyPixelCollectionName.empty()) {
			streamlog_out ( WARNING1 ) << "_noisyPixelCollectionName " << _noisyPixelCollectionName.c_str() << " not found" << std::endl;
			streamlog_out ( WARNING1 ) << "READ CAREFULLY: This means that no noisy pixels will be removed, despite the processor successfully running!" << std::endl;
		}
		return;
    	}

	//Decoder to get sensor ID
	CellIDDecoder<TrackerDataImpl> cellDecoder( noisyPixelCollectionVec );
        EUTelBaseSparsePixel* pixel = nullptr;

	//Loop over all hot pixels
	for(int i=0; i<  noisyPixelCollectionVec->getNumberOfElements(); i++) {
		//Get the TrackerData for the sensor ID
		TrackerDataImpl* noisyPixelData = dynamic_cast< TrackerDataImpl *> ( noisyPixelCollectionVec->getElementAt( i ) );
		int sensorID = cellDecoder( noisyPixelData )["sensorID"];
		int pixelType = cellDecoder( noisyPixelData )["sparsePixelType"];
		
		//And get the corresponding noise vector for that plane
		std::vector<int>* noiseSensorVector = &(_noisyPixelMap[sensorID]);
		std::unique_ptr<EUTelTrackerDataInterfacer> noisyPixelDataInterface;

		if( pixelType == kEUTelGenericSparsePixel ) {
			noisyPixelDataInterface =  std::unique_ptr<EUTelTrackerDataInterfacer>( new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>(noisyPixelData) );
		} else {
			streamlog_out( ERROR5 ) << "The noisy pixel collection is corrupted, it does not contain the right pixel type. Something is wrong!" << std::endl;
		}
		
		//Store all the noisy pixels in the noise vector, use the provided encoding to map two int's to an unique int
		for ( unsigned int iPixel = 0; iPixel < noisyPixelDataInterface->size(); iPixel++ ) {
			pixel = noisyPixelDataInterface->getSparsePixelAt( iPixel, pixel);
			noiseSensorVector->push_back( encode(pixel->getXCoord(), pixel->getYCoord()) );

		}
	}
	delete pixel;

	for( std::map<int, std::vector<int> >::iterator it = _noisyPixelMap.begin(); it != _noisyPixelMap.end(); ++it) {
		//Sort the noisy pixel maps
		std::sort( (it->second).begin(), (it->second).end() );
		std::cout << "Read in " << (it->second).size() << " hot pixels on plane " << (it->first) << std::endl;
	}
}

} //namespace eutelescope
