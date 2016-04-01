/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelProcessorRawHistos.h"
#include "EUTELESCOPE.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelTrackerDataInterfacerImpl.h"

// eutelescope geometry
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelGenericPixGeoDescr.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include "marlin/AIDAProcessor.h"
#include <AIDA/ITree.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IHistogramFactory.h>
#endif

// lcio includes <.h>
#include <LCIOTypes.h>
#include <IO/LCWriter.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <Exceptions.h>


// system includes <>
#include <map>
#include <memory>
#include <cmath>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <algorithm>

using namespace std;
using namespace marlin;
using namespace eutelescope;

EUTelProcessorRawHistos::EUTelProcessorRawHistos(): 
  Processor("EUTelProcessorRawHistos"),
  _zsDataCollectionName(""),
  _iRun(0),
  _iEvt(0),
  _sensorIDVec()
{
  //processor description
  _description = "EUTelProcessorRawHistos computes the firing frequency of pixels and applies a cut on this value to mask (NOT remove) hot pixels.";

  registerInputCollection(LCIO::TRACKERDATA, "ZSDataCollectionName", "Input of Zero Suppressed data", _zsDataCollectionName, std::string("zsdata") );
  registerOptionalParameter("NoisyPixelCollectionName", "Name of the noisy pixel collection", _noisyPixCollectionName, std::string("noisypixel"));
}

void EUTelProcessorRawHistos::init() {
	// this method is called only once even when the rewind is active usually a good idea to
	printParameters ();

	// set to zero the run and event counters
	_iRun = 0;
	_iEvt = 0;
	_noOfEvents = 0;
	
	geo::gGeometry().initializeTGeoDescription(EUTELESCOPE::GEOFILENAME, EUTELESCOPE::DUMPGEOROOT);
	
	_sensorIDVec = geo::gGeometry().sensorIDsVec();
	bookHistos();	
	_treatNoise = true;
}

int EUTelProcessorRawHistos::cantorEncode(int X, int Y) {
	//Cantor pairing function
	return static_cast<int>( 0.5*(X+Y)*(X+Y+1)+Y );
}

void EUTelProcessorRawHistos::initialiseNoisyPixels( LCCollectionVec* const noisyPixCollectionVec) {

	//Decoder to get sensor ID
	CellIDDecoder<TrackerDataImpl> cellDecoder( noisyPixCollectionVec );
	
	//Loop over all hot pixels
	for(int i=0; i<  noisyPixCollectionVec->getNumberOfElements(); i++) {
		//Get the TrackerData for the sensor ID
		TrackerDataImpl* noisyTrackerData = dynamic_cast<TrackerDataImpl*>( noisyPixCollectionVec->getElementAt(i) );
		int sensorID = cellDecoder( noisyTrackerData )["sensorID"];
		int pixelType = cellDecoder( noisyTrackerData )["sparsePixelType"];
		//And get the corresponding noise vector for that plane
		std::vector<int>* noiseSensorVector = &(_noisyPixelVecMap[sensorID]);

		if( pixelType == kEUTelGenericSparsePixel ) {
			auto noisyPixelData =  std::make_unique<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>>(noisyTrackerData);
			auto& pixelVec = noisyPixelData->getPixels();
			//Store all the noisy pixels in the noise vector, use the provided encoding to map two int's to an unique int
			for( auto& pixel: pixelVec ) {	
				noiseSensorVector->push_back( cantorEncode(pixel.getXCoord(), pixel.getYCoord()) );
			}
		} else { /*PANIC*/ }
	}

	for(auto& i: _noisyPixelVecMap) {
		std::sort(i.second.begin(), i.second.end() );
	}
}
 
void EUTelProcessorRawHistos::processRunHeader(LCRunHeader* rdr) {
	unique_ptr<EUTelRunHeaderImpl> runHeader( new EUTelRunHeaderImpl(rdr) );
	runHeader->addProcessor(type());
	// increment the run counter
	++_iRun;
	// reset the event counter
	_iEvt = 0;
}

void EUTelProcessorRawHistos::processEvent (LCEvent* event) {
	if( event == nullptr ) {
		streamlog_out ( WARNING2 ) <<  "Event does not exist! Skipping!" <<  std::endl;       
		return;
	}

	if( _iEvt == 0 ) {
		try 
		{
			auto noisyPixelCollection = static_cast< LCCollectionVec*>( event->getCollection(_noisyPixCollectionName) );
			initialiseNoisyPixels( noisyPixelCollection );		
		} catch (...) {
			if (!_noisyPixCollectionName.empty())
			{
				streamlog_out( WARNING1 ) << "_noisyPixCollectionName " << _noisyPixCollectionName << " not found" << endl;
			}
			_treatNoise = false;
		}
	}

	EUTelEventImpl* evt = static_cast<EUTelEventImpl*> (event);
	if ( evt->getEventType() == kEORE ) {
		streamlog_out ( DEBUG4 ) <<  "EORE found: nothing else to do." <<  std::endl;
		return;
	} else if ( evt->getEventType() == kUNKNOWN ) {
		streamlog_out ( WARNING2 ) << "Event number " << event->getEventNumber() << " is of unknown type. Continue considering it as a normal Data Event." << std::endl;
	}

	try
	{
		LCCollectionVec* zsInputCollectionVec  = dynamic_cast < LCCollectionVec * > (evt->getCollection( _zsDataCollectionName ));
		CellIDDecoder<TrackerDataImpl> cellDecoder( zsInputCollectionVec );

		std::map<int, size_t> rawHitsPerPlane;
		std::map<int, size_t> rawHitsPerPlaneNoNoise;

		for( auto sensorID: _sensorIDVec) {
			rawHitsPerPlane[sensorID] = 0;
			rawHitsPerPlaneNoNoise[sensorID] = 0;
		
		}

		for ( size_t iDetector = 0 ; iDetector < zsInputCollectionVec->size(); iDetector++ ) {

			// get the TrackerData and guess which kind of sparsified data it contains.
			TrackerDataImpl* zsData = dynamic_cast< TrackerDataImpl* > ( zsInputCollectionVec->getElementAt( iDetector ) );
			int sensorID            = static_cast<int > ( cellDecoder( zsData )["sensorID"] );

			// now prepare the EUTelescope interface to sparsified data.  
			auto sparseData = std::make_unique<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>>( zsData );
			auto& pixelVec = sparseData->getPixels();

			for( auto& genericPixel: pixelVec ) {
				bool isNoisy = false;
				
				int xCo = genericPixel.getXCoord(); 
				int yCo = genericPixel.getYCoord();
				int encoded = cantorEncode(xCo, yCo);

				if (_treatNoise) isNoisy = std::binary_search(_noisyPixelVecMap.at(sensorID).begin(), _noisyPixelVecMap.at(sensorID).end(), encoded );


				rawHitsPerPlane[sensorID]++;
				_chargeHisto.at(sensorID)->fill(genericPixel.getSignal());
				_timeHisto.at(sensorID)->fill(genericPixel.getTime());		
		
				if(!isNoisy) {	
					rawHitsPerPlaneNoNoise[sensorID]++;
					_chargeHistoNoNoise.at(sensorID)->fill(genericPixel.getSignal());
					_timeHistoNoNoise.at(sensorID)->fill(genericPixel.getTime());		
				}
			}
		}

		for(auto& i: rawHitsPerPlane) {
			_countHisto.at(i.first)->fill(i.second);
		}	
		for(auto& i: rawHitsPerPlaneNoNoise) {
			_countHistoNoNoise.at(i.first)->fill(i.second);
		}
	} catch (lcio::DataNotAvailableException& e ) {
		streamlog_out ( WARNING2 )  << "Input collection not found in the current event. Skipping..." << e.what() << std::endl;
		return;    
	}
	//don't forget to increment the event counter
	_iEvt++;
}

void EUTelProcessorRawHistos::end() {
}

void EUTelProcessorRawHistos::bookHistos() {
	streamlog_out ( MESSAGE1 ) << "Booking histograms " << std::endl;

	for( auto sensorID: _sensorIDVec) {
		auto basePath = "detector_" + to_string(sensorID);			
		AIDAProcessor::tree(this)->mkdir(basePath.c_str());
		basePath.append("/");

		auto const countHistoName = "raw_hits_per_trigger_d"+to_string(sensorID);
		auto const chargeHistoName = "raw_hits_charge_d"+to_string(sensorID);
		auto const timeHistoName = "raw_hits_time_d"+to_string(sensorID);

		auto rawHitCountHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + countHistoName).c_str(), 20, -0.5, 19.5 );
		auto rawHitChargeHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + chargeHistoName).c_str(), 20, -0.5, 19.5 );
		auto rawHitTimeHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + timeHistoName).c_str(),20, -0.5, 19.5 );

		auto rawHitCountHistoNoNoise = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + countHistoName+"_noNoise").c_str(), 20, -0.5, 19.5 );
		auto rawHitChargeHistoNoNoise = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + chargeHistoName+"_noNoise").c_str(), 20, -0.5, 19.5 );
		auto rawHitTimeHistoNoNoise = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + timeHistoName+"_noNoise").c_str(),20, -0.5, 19.5 );

		_countHisto[sensorID] = rawHitCountHisto;
		_chargeHisto[sensorID] = rawHitChargeHisto;
		_timeHisto[sensorID] = rawHitTimeHisto;

		_countHistoNoNoise[sensorID] = rawHitCountHistoNoNoise;
		_chargeHistoNoNoise[sensorID] = rawHitChargeHistoNoNoise;
		_timeHistoNoNoise[sensorID] = rawHitTimeHistoNoNoise;
	}
}
