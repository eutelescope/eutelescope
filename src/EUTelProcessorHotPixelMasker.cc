// Author Tobias Bisanz,  <tobias.bisanz@phys.uni-goettingen.de>
// Version $Id$
/*
 *   This processor removes hot pixels as specified by the hot pixel
 *   collection and removes them from any given other tracker collection.
 *   It is only campatible with the EUTelGenerisSparsePixel type, the
 *   hot pixel as well as the input collections must use it.
 *
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelProcessorHotPixelMasker.h"
#include "EUTELESCOPE.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelTrackerDataInterfacerImpl.h"

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

using namespace std;
using namespace marlin;
using namespace eutelescope;


EUTelProcessorHotPixelMasker::EUTelProcessorHotPixelMasker():
  Processor("EUTelProcessorHotPixelMasker"),
  _inputCollectionName(""),
  _iRun(0),
  _iEvt(0),
  _firstEvent(true),
  _dataFormatChecked(false),
  _wrongDataFormat(false)
{
  _description ="EUTelProcessorHotPixelMasker removes hot pixels from a tracker data collection, it reads in a hot pixel collection and an input collection. If any of the hits in the input collection is a hot pixel it gets removed.";

  registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName", "Input of zero suppressed data, still containing hot pixels", _inputCollectionName, string ("zsdata") );

  registerOutputCollection(LCIO::TRACKERDATA, "OutputCollectionName", "Hot Pixel free output collection name", _outputCollectionName, string("zsdatafree"));

  registerOptionalParameter("HotPixelCollectionName", "Name of the hot pixel collection.",  _hotPixelCollectionName, static_cast< string > ("hotpixel"));

}

void EUTelProcessorHotPixelMasker::init () 
{
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters();
  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;
}

void EUTelProcessorHotPixelMasker::processRunHeader(LCRunHeader* rdr){

  auto_ptr<EUTelRunHeaderImpl> runHeader ( new EUTelRunHeaderImpl(rdr) );
  runHeader->addProcessor(type()) ;
  // increment the run counter
  ++_iRun;
  // reset the event counter
  _iEvt = 0;
}

void EUTelProcessorHotPixelMasker::processEvent(LCEvent * event) 
{
	if(_firstEvent)
	{
		//The hot pixel collection stores all thot pixels in event #1
		//Thus we have to read it in in that case
		readHotPixelList(event);
		_firstEvent = false;
	}

 	// get the collection of interest from the event.
	LCCollectionVec* pulseInputCollectionVec = NULL;

	try
	{
    		pulseInputCollectionVec  = dynamic_cast <LCCollectionVec*>( event->getCollection(_inputCollectionName) );
	}
  	catch( lcio::DataNotAvailableException& e ) 
  	{
		return;
  	}

	// prepare decoder for input data
	CellIDDecoder<TrackerPulseImpl> cellDecoder( pulseInputCollectionVec );
	
	//read the encoding string from the input collection
	std::string encodingString = pulseInputCollectionVec->getParameters().getStringVal( LCIO::CellIDEncoding );
	//and the encoder for the output data
	CellIDEncoder<TrackerPulseImpl> cellEncoder( encodingString , pulseInputCollectionVec);

	//now prepare output collection
	LCCollectionVec* outputPulseCollection;
	bool outputPulseCollectionExists = false;
  	_initialoutputPulseCollectionSize = 0;

  	try 
  	{
   		outputPulseCollection = dynamic_cast< LCCollectionVec* > ( event->getCollection( _outputCollectionName ) );
    		outputPulseCollectionExists = true;
    		_initialoutputPulseCollectionSize = outputPulseCollection->size();
  	} 
  	catch ( lcio::DataNotAvailableException& e ) 
  	{
    		outputPulseCollection = new LCCollectionVec(LCIO::TRACKERPULSE);
  	}

	//and the decoder for output data
	CellIDEncoder<TrackerPulseImpl> outputEncoder(encodingString,  outputPulseCollection);
	
	//loop over all the pulses
	for ( size_t iPulse = 0 ; iPulse < pulseInputCollectionVec->size(); iPulse++ ) 
	{
        	TrackerPulseImpl* pulseData = dynamic_cast<TrackerPulseImpl*> ( pulseInputCollectionVec->getElementAt( iPulse ) );
		int sensorID = cellDecoder(pulseData)["sensorID"];		
	
	        //get the noise vector for the given plane
		std::vector<int>* noiseVector = &(_hotPixelMap[sensorID]);
		
		//get the tracker data from the pulse
		TrackerDataImpl* trackerData = dynamic_cast<TrackerDataImpl*>( pulseData->getTrackerData() );
		//decoder for tracker data
		CellIDDecoder<TrackerDataImpl> trackerDecoder ( EUTELESCOPE::ZSCLUSTERDEFAULTENCODING );
		int pixelType = trackerDecoder(trackerData)["sparsePixelType"];

		//Interface to sparsified data
                auto_ptr<EUTelTrackerDataInterfacer> sparseData = auto_ptr<EUTelTrackerDataInterfacer>();

		EUTelBaseSparsePixel* pixel = NULL;
		if( pixelType == kEUTelSimpleSparsePixel )
		{
			sparseData =  auto_ptr<EUTelTrackerDataInterfacer>( new EUTelTrackerDataInterfacerImpl<EUTelSimpleSparsePixel>(trackerData) );
		}
		else
                {
		}

		bool noisy = false;

		//Loop over all hits!
		for ( unsigned int iPixel = 0; iPixel < sparseData->size(); iPixel++ )
       		{
		        pixel = sparseData->getSparsePixelAt( iPixel, pixel );
			if(std::binary_search( noiseVector->begin(), noiseVector->end(), encode(pixel->getXCoord(), pixel->getYCoord()) ))
			{
				noisy=true;
				break;
			}
		}

		/*if(noisy)
		{
			int quality = cellDecoder(pulseData)["quality"];
			quality = quality || kNoisyCluster;
			cellEncoder["quality"] = quality;
			cellEncoder.setCellID(pulseData);
		}*/
                
		//TrackerPulseImpl for the output collection
		auto_ptr<TrackerPulseImpl> outputPulse ( new TrackerPulseImpl );
		//copy the information which is the same
		outputPulse->setCellID0( pulseData->getCellID0() );
		outputPulse->setCellID1( pulseData->getCellID1() );
		if(noisy)
		{
			//int quality = cellDecoder(outputPulse)["quality"];
			int quality = 8;
			outputEncoder["quality"] = quality;
			outputEncoder.setCellID(outputPulse.get());
		}
		outputPulse->setTime( pulseData->getTime() );
		outputPulse->setCharge( pulseData->getCharge() );
		outputPulse->setCovMatrix( pulseData->getCovMatrix() );
		
		outputPulse->setQuality( pulseData->getQuality() );
		//outputPulse->setTrackerData( pulseData->getTrackerData() );
		
		outputPulseCollection->push_back( outputPulse.release()  );
		
		delete pixel;
        }// loop over detectors

	//add the collection if necessary
	if ( !outputPulseCollectionExists && ( outputPulseCollection->size() != _initialoutputPulseCollectionSize )) 
	{
		event->addCollection( outputPulseCollection, _outputCollectionName );
	}

	if ( !outputPulseCollectionExists && ( outputPulseCollection->size() == _initialoutputPulseCollectionSize ) ) 
	{
		delete outputPulseCollection;
	}	
//rest of memory cleaned up by auto_ptrs
}

void EUTelProcessorHotPixelMasker::end() 
{
	//Print out some stats for the user
	streamlog_out ( MESSAGE4 ) << "Hot pixel remover successfully finished" << endl;
	streamlog_out ( MESSAGE4 ) << "Printing summary:" << endl;
	for(std::map<int,int>::iterator it = _removedHotPixels.begin(); it != _removedHotPixels.end(); ++it)
	{
  		streamlog_out ( MESSAGE4 ) << "Removed " << (*it).second << " hot pixels from plane " << (*it).first << "." << endl;
	}
}

int EUTelProcessorHotPixelMasker::encode(int X, int Y)
{
	//Cantor pairing function
	return static_cast<int>( 0.5*(X+Y)*(X+Y+1)+Y );
} 

void EUTelProcessorHotPixelMasker::readHotPixelList(LCEvent* event)
{
	//Preapare pointer to hot pixel collection
	LCCollectionVec* hotPixelCollectionVec = NULL;

	//Try to obtain the collection
	try 
	{
		hotPixelCollectionVec = static_cast< LCCollectionVec*>( event->getCollection(_hotPixelCollectionName) );
	}	
	catch (...)
    	{
		if (!_hotPixelCollectionName.empty())
		{
			streamlog_out ( WARNING1 ) << "_hotPixelCollectionName " << _hotPixelCollectionName.c_str() << " not found" << endl;
			streamlog_out ( WARNING1 ) << "READ CAREFULLY: This means that no hot pixels will be removed, despite the processor successfully running!" << endl;
		}
		return;
    	}

	//Decoder to get sensor ID
	CellIDDecoder<TrackerDataImpl> cellDecoder( hotPixelCollectionVec );

        EUTelBaseSparsePixel* pixel = NULL;

	//Loop over all hot pixels
	for(int i=0; i<  hotPixelCollectionVec->getNumberOfElements(); i++)
	{
		//Get the TrackerData for the sensor ID
		TrackerDataImpl* hotPixelData = dynamic_cast< TrackerDataImpl *> ( hotPixelCollectionVec->getElementAt( i ) );
		int sensorID = cellDecoder( hotPixelData )["sensorID"];
		int pixelType = cellDecoder( hotPixelData )["sparsePixelType"];
		//And get the corresponding noise vector for that plane
		std::vector<int>* noiseSensorVector = &(_hotPixelMap[sensorID]);

		auto_ptr<EUTelTrackerDataInterfacer> noisyPixelData = auto_ptr<EUTelTrackerDataInterfacer>();

		if( pixelType == kEUTelSimpleSparsePixel )
		{
			noisyPixelData =  auto_ptr<EUTelTrackerDataInterfacer>( new EUTelTrackerDataInterfacerImpl<EUTelSimpleSparsePixel>(hotPixelData) );
		}
		//else if
		else
		{
		}
		//Store all the noisy pixels in the noise vector, use the provided encoding to map two int's to an unique int
		for ( unsigned int iPixel = 0; iPixel < noisyPixelData->size(); iPixel++ ) 
		{
			pixel = noisyPixelData->getSparsePixelAt( iPixel, pixel);
			noiseSensorVector->push_back( encode(pixel->getXCoord(), pixel->getYCoord()) );

		}
	}

	delete pixel;

	for( std::map<int, std::vector<int> >::iterator it = _hotPixelMap.begin(); it != _hotPixelMap.end(); ++it)
	{
		//Sort the hot pixel maps
		std::sort( (it->second).begin(), (it->second).end() );
		std::cout << "Read in " << (it->second).size() << " hot pixels on plane " << (it->first) << std::endl;
	}
}

