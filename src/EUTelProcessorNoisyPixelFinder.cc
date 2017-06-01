/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelProcessorNoisyPixelFinder.h"
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

namespace eutelescope {

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelProcessorNoisyPixelFinder::_firing2DHistoName = "Firing2D";
std::string EUTelProcessorNoisyPixelFinder::_firing1DHistoName = "Firing1D";
#endif

EUTelProcessorNoisyPixelFinder::EUTelProcessorNoisyPixelFinder(): 
  Processor("EUTelProcessorNoisyPixelFinder"),
  _zsDataCollectionName(""),
  _noisyPixelCollectionName(""),
  _excludedPlanes(),
  _noOfEvents(0),
  _maxAllowedFiringFreq(0.0),
  _iRun(0),
  _iEvt(0),
  _sensorIDVec(),
  _noisyPixelDBFile(""),
  _finished(false)
{
  //processor description
  _description = "EUTelProcessorNoisyPixelFinder computes the firing frequency of pixels and applies a cut on this value to mask (NOT remove) noisy pixels.";

  registerInputCollection (LCIO::TRACKERDATA, "ZSDataCollectionName", "Input of Zero Suppressed data",_zsDataCollectionName, std::string("zsdata") );

  registerProcessorParameter("NoOfEvents", "The number of events to be considered for each update cycle", _noOfEvents, static_cast<int>(100) );

  registerOptionalParameter("SensorIDVec", "The sensorID for the generated collection (one per detector)", _sensorIDVec, std::vector<int>() );

  registerProcessorParameter("MaxAllowedFiringFreq", "This float number [0,1] represents the maximum allowed firing frequency\n"
                             "within the selected number of event per cycle", _maxAllowedFiringFreq, static_cast<float>(0.2) );

  registerOptionalParameter("HotPixelDBFile","This is the name of the LCIO file name with the output noisyPixel db (add .slcio)",
                             _noisyPixelDBFile, std::string("noisyPixel.slcio"));

  registerOptionalParameter("ExcludedPlanes", "The list of sensor IDs that shall be excluded.", _excludedPlanes, std::vector<int> () );

  registerOptionalParameter("HotPixelCollectionName", "This is the name of the hot pixel collection to be saved into the output slcio file",
                             _noisyPixelCollectionName, std::string("noisyPixel"));
}

void EUTelProcessorNoisyPixelFinder::initializeHitMaps() {
	//it stored detectoID and itt stores the vector for the y-entries
	for(auto sensorID: _sensorIDVec) {
		try {
			//get the geoemtry description of the plane
			geo::EUTelGenericPixGeoDescr* geoDescr = geo::gGeometry().getPixGeoDescr(sensorID) ;

			//get the max/min pixel indices
			int minX, minY, maxX, maxY;
			minX = minY = maxX = maxY = 0;
			geoDescr->getPixelIndexRange( minX, maxX, minY, maxY );

			//a sensor structs holds all the helpful information
			sensor thisSensor;

			//store the information in them
			thisSensor.offX = minX;
			thisSensor.sizeX =  maxX - minX+1;
			thisSensor.offY = minY;
			thisSensor.sizeY = maxY - minY+1;

			//this is the 2-dimensional array used to store all the pixels, it is implemented as a vector of vectors
			//first the x-entries-vector
			std::vector<std::vector<int>> hitVecP = std::vector<std::vector<int>>( thisSensor.sizeX );

			//and the y-entries vector
			for(auto& i: hitVecP) {
				//resize it now
				i.resize(thisSensor.sizeY, 0);
			}

			//collection to later hold the hot pixels
			std::vector<EUTelGenericSparsePixel> noisyPixelMap;

			//store all the collections/pointers in the corresponding maps
		    	_sensorMap[sensorID] = thisSensor;
			_hitVecMap[sensorID] = hitVecP;
			_noisyPixelMap[sensorID] = noisyPixelMap;
		} catch(std::runtime_error& e) {
			streamlog_out ( ERROR0 ) << "Noisy pixel masker could not retrieve plane " << sensorID << std::endl;
			streamlog_out ( ERROR0 ) << e.what() << std::endl;
			throw marlin::StopProcessingException(this);
		}
	}
}

void EUTelProcessorNoisyPixelFinder::init() {
	// this method is called only once even when the rewind is active usually a good idea to
	printParameters ();

	// set to zero the run and event counters
	_iRun = 0;
	_iEvt = 0;

	geo::gGeometry().initializeTGeoDescription(EUTELESCOPE::GEOFILENAME, EUTELESCOPE::DUMPGEOROOT);

	//and use it to prepare the hit maps
	initializeHitMaps();
}

void EUTelProcessorNoisyPixelFinder::processRunHeader(LCRunHeader* /*rdr*/) {
	// increment the run counter
	++_iRun;
	// reset the event counter
	_iEvt = 0;
}

void EUTelProcessorNoisyPixelFinder::noisyPixelFinder(EUTelEventImpl* evt) {
	if (evt == nullptr) {
		return; //exit(-1);
	}
	try {
		// get the collections of interest from the event.
		LCCollectionVec* zsInputCollectionVec  = dynamic_cast < LCCollectionVec * > (evt->getCollection( _zsDataCollectionName ));
		// prepare some decoders
		CellIDDecoder<TrackerDataImpl> cellDecoder( zsInputCollectionVec );

		for( size_t iDetector = 0; iDetector < zsInputCollectionVec->size(); iDetector++ ) {
			// get the TrackerData and guess which kind of sparsified data it contains.
			TrackerDataImpl* zsData = dynamic_cast<TrackerDataImpl*>( zsInputCollectionVec->getElementAt(iDetector) );
			int sensorID            = static_cast<int>( cellDecoder(zsData)["sensorID"] );

			sensor* currentSensor = &_sensorMap[sensorID];
			std::vector<std::vector<int>>* hitArray = &_hitVecMap[sensorID];

			//if this is an excluded sensor go to the next element
			bool foundexcludedsensor = false;
			for(auto i : _excludedPlanes) {
				if(i == sensorID) {
					foundexcludedsensor = true;
				}
			}
			if(foundexcludedsensor) continue;

			// now prepare the EUTelescope interface to sparsified data.  
			int pixelType = cellDecoder(zsData)["sparsePixelType"];
			auto sparseData = Utility::getSparseData(zsData, pixelType);
			auto basePixelPtrVec = sparseData->getBasePixelPtrVec();

			// loop over all pixels in the sparseData object, these are the hit pixels!
			for(auto pixel: basePixelPtrVec) {
				//compute the address in the array-like-structure, any offset
				//has to be substracted (array index starts at 0)
				int indexX = pixel->getXCoord() - currentSensor->offX;
				int indexY = pixel->getYCoord() - currentSensor->offY;

				try {
					//increment the hit counter for this pixel
					(hitArray->at(indexX)).at(indexY)++;
				} catch(std::out_of_range& e) {
					streamlog_out ( ERROR5 )  << "Pixel: " << pixel->getXCoord() << "|" <<  pixel->getYCoord() << " on plane: " << sensorID << " fired." << std::endl 
						<< "This pixel is out of the range defined by the geometry. Either your data is corrupted or your pixel geometry not specified correctly!" << std::endl;
				}
			}
		}    
	} catch (lcio::DataNotAvailableException& e ) {
		streamlog_out ( WARNING2 )  << "Input collection not found in the current event. Skipping..." << e.what() << std::endl;
		return;    
	} 
}

void EUTelProcessorNoisyPixelFinder::processEvent (LCEvent * event) {
	//if we are over the number of events we need we just skip
	if(_noOfEvents < _iEvt) {
		++_iEvt;
		return;
	}

	if( event == nullptr ) {
		streamlog_out ( WARNING2 ) <<  "event does not exist!. skip " <<  std::endl;       
		return;
	}

	EUTelEventImpl* evt = static_cast<EUTelEventImpl*> (event);
	
	if( evt->getEventType() == kEORE ) 
	{
		streamlog_out ( DEBUG4 ) <<  "EORE found: nothing else to do." <<  std::endl;
		return;
	} else if ( evt->getEventType() == kUNKNOWN ) {
		streamlog_out ( WARNING2 ) << "Event number " << event->getEventNumber() << " is of unknown type. Continue considering it as a normal Data Event." << std::endl;
	}

	//if we don't skip the event, we pass it to NoisyPixelFinder() which does all the work
	try {
		noisyPixelFinder(evt);
	} catch (lcio::DataNotAvailableException& e ) {
		streamlog_out ( WARNING2 )  << "Input collection not found in the current event. Skipping..." << e.what() << std::endl;
		return;    
	} catch( marlin::ParseException& e ) {
		streamlog_out ( MESSAGE5 )  << "Input collection not found in the current event. Skipping..." << e.what() << std::endl;
		return;     
	}

	//don't forget to increment the event counter
	++_iEvt;
}

void EUTelProcessorNoisyPixelFinder::end() {
	if(_finished) {
		streamlog_out ( MESSAGE4 ) << "Noisy pixel finder has successfully finished!" << std::endl;
	} else {
		streamlog_out ( ERROR3 ) << "End of run reached before enough events were processed for noisy pixel finder, databases have NOT been written!" << std::endl;
		streamlog_out ( ERROR3 ) << "Increase the number or events to be processed for the run or decrease the number required for noisy pixel finding" << std::endl;
	}
}

void EUTelProcessorNoisyPixelFinder::check(LCEvent* /*event*/ ) {
	//only if the eventNo is the amount of events to be processed we analyse the data
	//since check() runs after the event and we increment the event counter before that
	//the comparison is valid. If we set _noOfEvents=1 we only want to have processed event
	//0 before calling this. Since we increment 0++ before calling this function, this 
	//criteria is fullfilled
	if( _iEvt == _noOfEvents) {
		streamlog_out ( MESSAGE4 ) << "Finished determining hot pixels, writing them out..." << std::endl;

		//iterate over all the sensors in our sensorMap
		for(auto& thisSensor: _sensorMap)
		{
			auto sensorID = thisSensor.first;
			streamlog_out ( MESSAGE3 ) << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
			streamlog_out ( MESSAGE3 ) << "Noisy pixels found on plane " << sensorID << " (max. fire freq set to: " << _maxAllowedFiringFreq << ")" << std::endl;
			streamlog_out ( MESSAGE3 ) << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

			//get the correpsonding hit array-like vector of vectors
			std::vector<std::vector<int>>* hitVector =  &_hitVecMap[sensorID];
			//and the sensor which stores offsets
			sensor* currentSensor = &_sensorMap[sensorID];

			//loop over all pixels
			for (std::vector<std::vector<int> >::iterator xIt = hitVector->begin(); xIt != hitVector->end(); ++xIt)
			{
				for(std::vector<int>::iterator yIt = xIt->begin(); yIt != xIt->end(); ++yIt)
				{
					//compute the firing frequency
					float fireFreq = (float)*yIt/(float)_iEvt;
					//if it is larger than the allowed one, we write this pixel into a collection
					if(fireFreq > _maxAllowedFiringFreq) {
						streamlog_out ( MESSAGE3 )	<< "Pixel: " << xIt-hitVector->begin() + currentSensor->offX << "|" 
										<< yIt - xIt->begin() + currentSensor->offY  << " fired " << fireFreq << std::endl;
						EUTelGenericSparsePixel pixel;
						pixel.setXCoord(xIt-hitVector->begin() + currentSensor->offX);
						pixel.setYCoord(yIt - xIt->begin() + currentSensor->offY );
						pixel.setSignal( fireFreq );
						//writing out is done here
						_noisyPixelMap[sensorID].push_back(pixel);
					}
				}
			}
		}

		//write out the databases and histograms
		noisyPixelDBWriter();
		bookAndFillHistos();
		//we reached enough events, wrote out noisy pixel db and are done now
		_finished = true;
	}
}

void EUTelProcessorNoisyPixelFinder::noisyPixelDBWriter() {    
	streamlog_out ( DEBUG5 ) << "Writing out hot pixel db into " << _noisyPixelDBFile.c_str() << std::endl;

	// reopen the LCIO file this time in append mode
	LCWriter* lcWriter = LCFactory::getInstance()->createLCWriter();
	std::unique_ptr<LCEventImpl> event (new LCEventImpl);
	// create new file if requested OR opening of existing file did not succeed
	try {
		lcWriter->open( _noisyPixelDBFile, LCIO::WRITE_NEW );
		LCRunHeaderImpl lcHeader;
		lcHeader.setRunNumber( 0 );
		lcWriter->writeRunHeader(&lcHeader);

		event->setRunNumber( 0 );
		event->setEventNumber( 0 );
		event->setDetectorName("EUTelNoisyPixel");
		
		streamlog_out ( DEBUG5 ) << "otPixelDB file: run header and event created ok"  << std::endl;       
		LCTime now;
		event->setTimeStamp( now.timeStamp() );        
	} catch ( IOException& e ) {
		streamlog_out ( ERROR4 ) << e.what() << std::endl << "Sorry, was not able to create new NoisyPixelDB file " << _noisyPixelDBFile <<", will quit now. " << std::endl;
		return;
		//exit(-1);
	}

	if(event == nullptr) {
		streamlog_out ( ERROR5 ) << "Problem opening NoisyPixelDB file, event is nullptr" << std::endl;
		return;  
	}

	// create main collection to be saved into the db file 
	LCCollectionVec* noisyPixelCollection;
	// create new or open existing noisyPixel collection
	try {
		noisyPixelCollection = static_cast< LCCollectionVec* > ( event->getCollection( _noisyPixelCollectionName  ) );
		std::cout << "noisyPixelCollection: " << _noisyPixelCollectionName << " found found with " << noisyPixelCollection->getNumberOfElements() << " elements " <<  std::endl; 
		noisyPixelCollection->clear();
		std::cout << "noisyPixelCollection: " << _noisyPixelCollectionName << " cleared: now " << noisyPixelCollection->getNumberOfElements() << " elements " <<  std::endl; 
	} catch ( lcio::DataNotAvailableException& e ) {
		noisyPixelCollection = new LCCollectionVec( lcio::LCIO::TRACKERDATA );
		event->addCollection( noisyPixelCollection, _noisyPixelCollectionName );
		std::cout << "noisyPixelCollection: " << _noisyPixelCollectionName << " created" <<  std::endl; 
	}

	streamlog_out( MESSAGE5 ) << "Noisy Pixel Finder summary:" << std::endl;
	//_noisyPixelMap holds the sensor if (first) and a vector of noisy pixels (second)
	for(auto& mapEntry: _noisyPixelMap) {
		CellIDEncoder< TrackerDataImpl > noisyPixelEncoder  ( EUTELESCOPE::ZSDATADEFAULTENCODING, noisyPixelCollection  );
		noisyPixelEncoder["sensorID"]        = mapEntry.first;
		noisyPixelEncoder["sparsePixelType"] = kEUTelGenericSparsePixel;

		// prepare a new TrackerData for the hot Pixel data
		std::unique_ptr<lcio::TrackerDataImpl> currentFrame( new lcio::TrackerDataImpl );
		noisyPixelEncoder.setCellID( currentFrame.get() );

		// this is the structure that will host the sparse pixel  
		std::unique_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>>
			sparseFrame( new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>(currentFrame.get()) );

		for( auto& pixel: mapEntry.second) {
			sparseFrame->addSparsePixel( pixel );                
		}
		noisyPixelCollection->push_back( currentFrame.release() );

		streamlog_out( MESSAGE5 ) << "Found " << mapEntry.second.size() << " noisy pixels on sensor: " << mapEntry.first << std::endl;
	}
	lcWriter->writeEvent( event.get() );
	lcWriter->close();
}

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

void EUTelProcessorNoisyPixelFinder::bookAndFillHistos() {

	streamlog_out ( MESSAGE1 ) << "Booking and filling histograms " << std::endl;
	std::string tempHistoName, basePath;

	for( auto det:_sensorIDVec) {   
		basePath = "detector_" + std::to_string(det) ;
		marlin::AIDAProcessor::tree(this)->mkdir(basePath.c_str());

		basePath.append("/");

		tempHistoName = _firing2DHistoName + "_d" + std::to_string(det);
		sensor* currentSensor = &_sensorMap[det];

		//determine range for 2D firing histo
		int     xBin = currentSensor->sizeX +1 ;
		double  xMin = static_cast<double >( currentSensor->offX ) - 0.5;
		double  xMax = static_cast<double >( currentSensor->offX + currentSensor->sizeX) + 0.5;

		int     yBin = currentSensor->sizeY +1 ;
		double  yMin = static_cast<double >( currentSensor->offY ) - 0.5;
		double  yMax = static_cast<double >( currentSensor->offY + currentSensor->sizeY) + 0.5;

		AIDA::IHistogram2D* firing2DHisto = marlin::AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(), xBin, xMin, xMax,yBin, yMin, yMax);

		if(!firing2DHisto) {
			streamlog_out ( ERROR5 )	<< "CreateHistogram2D for " <<  (basePath + tempHistoName).c_str()  << " failed " << std::endl
							<< "Execution stopped, check that your path (" << basePath.c_str() << ")exists  " << std::endl;
			return;
		}

		firing2DHisto->setTitle("Firing frequency map of hot pixels (in percent, rounded to the next integer); Pixel Index X; Pixel Index Y; Percent (%)");
		tempHistoName = _firing1DHistoName + "_d" + std::to_string(det);

		//range for 1D firing histo
		int nBin = 101;
		double min = -0.5;
		double max = 100.5;

		AIDA::IHistogram1D* firing1DHisto = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),nBin, min, max );
		firing1DHisto->setTitle("Firing frequency distribution of hot pixels (in percent, rounded to the next integer);Firing Frequency (%); Count (#)");

		//actually fill both histos for each detector
		for ( auto& pixel: _noisyPixelMap[det]) {
			firing2DHisto->fill(pixel.getXCoord(), pixel.getYCoord(), pixel.getSignal());
			firing1DHisto->fill(pixel.getSignal());
		}
	}//loop over dteectors
}
#endif


}//namespace eutelescope
