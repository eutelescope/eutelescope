/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelNoisyPixelFinder.h"
#include "EUTELESCOPE.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelTrackerDataInterfacerImpl.h"

// eutelescope geometry
#include "EUTelGenericPixGeoDescr.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

// marlin includes ".h"
#include "marlin/Exceptions.h"
#include "marlin/Processor.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include "marlin/AIDAProcessor.h"
#include <AIDA/IHistogram1D.h>
#include <AIDA/ICloud1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <IO/LCWriter.h>
#include <LCIOTypes.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTime.h>

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <Exceptions.h>

// system includes <>
#include <cmath>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <stdexcept>

namespace eutelescope
{

  EUTelNoisyPixelFinder::
    EUTelNoisyPixelFinder ():Processor ("EUTelNoisyPixelFinder"),
    _zsDataCollectionName (""), _noisyPixelCollectionName (""),
    _excludedPlanes (), _noOfEvents (0), _maxAllowedFiringFreq (0.0),
    _iRun (0), _iEvt (0), _sensorIDVec (), _noisyPixelDBFile (""),
    _finished (false)
  {

    _description = "EUTelNoisyPixelFinder computes the firing "
      "frequency of pixels and applies a cut on this value to "
      "mask (NOT remove) noisy pixels.";

    registerInputCollection (LCIO::TRACKERDATA,
			     "ZSDataCollectionName",
			     "Input of Zero Suppressed data",
			     _zsDataCollectionName, std::string ("zsdata"));

    registerProcessorParameter ("NoOfEvents",
				"The number of events to be considered for each update cycle",
				_noOfEvents, 100);

    registerProcessorParameter ("MaxAllowedFiringFreq",
				"This float number [0,1] represents the maximum "
				"allowed firing frequency\n"
				"within the selected number of event per cycle",
				_maxAllowedFiringFreq, 0.2f);

    registerOptionalParameter ("SensorIDVec",
			       "The sensorID for the generated collection (one per detector)",
			       _sensorIDVec, std::vector < int >());

      registerOptionalParameter ("HotPixelDBFile",
				 "This is the name of the LCIO file name with the "
				 "output noisyPixel db (add .slcio)",
				 _noisyPixelDBFile,
				 std::string ("noisyPixel.slcio"));

      registerOptionalParameter ("ExcludedPlanes",
				 "The list of sensor IDs that shall be excluded.",
				 _excludedPlanes, std::vector < int >());

      registerOptionalParameter ("HotPixelCollectionName",
				 "This is the name of the hot pixel collection to "
				 "be saved into the output slcio file",
				 _noisyPixelCollectionName,
				 std::string ("noisyPixel"));
  }

  void EUTelNoisyPixelFinder::initializeHitMaps ()
  {
    //store detectorID and vector for the y-entries
  for (auto sensorID:_sensorIDVec)
      {
	try
	{
	  //get the geometry description of the plane
	  geo::EUTelGenericPixGeoDescr * geoDescr =
	    geo::gGeometry ().getPixGeoDescr (sensorID);

	  //get the max/min pixel indices
	  int minX, minY, maxX, maxY;
	  minX = minY = maxX = maxY = 0;
	  geoDescr->getPixelIndexRange (minX, maxX, minY, maxY);

	  //a sensor structs holds all the helpful information
	  sensor thisSensor;
	  thisSensor.offX = minX;
	  thisSensor.sizeX = maxX - minX + 1;
	  thisSensor.offY = minY;
	  thisSensor.sizeY = maxY - minY + 1;

	  //this is the 2-dimensional array used to store all the pixels, it is
	  //implemented as a vector of vectors
	  //first the x-entries-vector
	  std::vector < std::vector < long int >>hitVecP =
	    std::vector < std::vector < long int >>(thisSensor.sizeX);

	  //and the y-entries vector
	for (auto & ihit:hitVecP)
	    {
	      //resize it now
	      ihit.resize (thisSensor.sizeY, 0);
	    }

	  //make vector for firing frequency
	  auto & firingFreqVec = _firingFreqForAllPixels[sensorID];
	  firingFreqVec.resize (thisSensor.sizeX * thisSensor.sizeY);

	  //collection to later hold the hot pixels
	  std::vector < EUTelGenericSparsePixel > noisyPixelMap;

	  //store all the collections/pointers in the corresponding maps
	  _sensorMap[sensorID] = thisSensor;
	  _hitVecMap[sensorID] = hitVecP;
	  _noisyPixelMap[sensorID] = noisyPixelMap;
	}
	catch (std::runtime_error & e)
	{
	  streamlog_out (ERROR0) <<
	    "Noisy pixel masker could not retrieve plane " << sensorID <<
	    std::endl;
	  streamlog_out (ERROR0) << e.what () << std::endl;
	  throw marlin::StopProcessingException (this);
	}
      }
  }

  void EUTelNoisyPixelFinder::init ()
  {

    printParameters ();

    //set to zero the run and event counters
    _iRun = 0;
    _iEvt = 0;

    geo::gGeometry ().initializeTGeoDescription (EUTELESCOPE::GEOFILENAME,
						 EUTELESCOPE::DUMPGEOROOT);

    //prepare the hit maps
    initializeHitMaps ();
  }

  void EUTelNoisyPixelFinder::processRunHeader (LCRunHeader * /*rdr */ )
  {
    //increment the run counter
    ++_iRun;
    //reset the event counter
    _iEvt = 0;
  }

  void EUTelNoisyPixelFinder::noisyPixelFinder (EUTelEventImpl * evt)
  {
    if (evt == nullptr)
      {
	return;
      }
    try
    {
      //get the collections of interest from the event
      LCCollectionVec *zsInputCollectionVec =
	dynamic_cast <
	LCCollectionVec * >(evt->getCollection (_zsDataCollectionName));
      //prepare some decoders
      CellIDDecoder < TrackerDataImpl > cellDecoder (zsInputCollectionVec);

      for (size_t iDetector = 0; iDetector < zsInputCollectionVec->size ();
	   iDetector++)
	{
	  //get TrackerData and guess which kind of sparsified data it contains
	  TrackerDataImpl *zsData =
	    dynamic_cast <
	    TrackerDataImpl *
	    >(zsInputCollectionVec->getElementAt (iDetector));
	  int sensorID =
	    static_cast < int >(cellDecoder (zsData)["sensorID"]);

	  sensor *currentSensor = &_sensorMap[sensorID];
	  std::vector < std::vector < long int >>*hitArray =
	    &_hitVecMap[sensorID];

	  //if this is an excluded sensor go to the next element
	  bool foundExcludedSensor = false;
	for (auto planeID:_excludedPlanes)
	    {
	      if (planeID == sensorID)
		foundExcludedSensor = true;
	    }
	  if (foundExcludedSensor)
	    continue;

	  //now prepare the EUTelescope interface to sparsified data
	  int pixelType = cellDecoder (zsData)["sparsePixelType"];
	  auto sparseData = Utility::getSparseData (zsData, pixelType);

	  //loop over all pixels in the sparseData object, these are the hit pixels
	for (auto & pixelRef:*sparseData)
	    {
	      auto & pixel = pixelRef.get ();

	      //compute the address in the array-like-structure, any offset
	      //has to be substracted (array index starts at 0)
	      int indexX = pixel.getXCoord () - currentSensor->offX;
	      int indexY = pixel.getYCoord () - currentSensor->offY;

	      try
	      {
		//increment the hit counter for this pixel
		(hitArray->at (indexX)).at (indexY)++;
	      } catch (std::out_of_range & e)
	      {
		streamlog_out (ERROR5)
		  << "Pixel: " << pixel.getXCoord () << "|" << pixel.
		  getYCoord () << " on plane: " << sensorID << " fired." <<
		  std::
		  endl <<
		  "This pixel is out of the range defined by the geometry. "
		  "Either your data is corrupted or your pixel geometry not "
		  "specified correctly!" << std::endl;
	      }
	    }
	}
    } catch (lcio::DataNotAvailableException & e)
    {
      streamlog_out (WARNING2)
	<< "Input collection not found in the current event. Skipping..."
	<< e.what () << std::endl;
      return;
    }
  }

  void EUTelNoisyPixelFinder::processEvent (LCEvent * event)
  {
    //if we are over the number of events we need, we just skip
    if (_noOfEvents < _iEvt)
      {
	++_iEvt;
	return;
      }
    //check if event exists
    if (event == nullptr)
      {
	streamlog_out (WARNING2) << "Event does not exist! Skip. " << std::
	  endl;
	return;
      }

    EUTelEventImpl *evt = static_cast < EUTelEventImpl * >(event);
    //check event type
    if (evt->getEventType () == kEORE)
      {
	streamlog_out (DEBUG4) << "EORE found: nothing else to do." << std::
	  endl;
	return;
      }

    //if we don't skip the event, we pass it to NoisyPixelFinder() which does
    //all the work
    try
    {
      noisyPixelFinder (evt);
    }
    catch (lcio::DataNotAvailableException & e)
    {
      streamlog_out (WARNING2)
	<< "Input collection not found in the current event. Skipping..."
	<< e.what () << std::endl;
      return;
    }
    catch (marlin::ParseException & e)
    {
      streamlog_out (MESSAGE5)
	<< "Input collection not found in the current event. Skipping..."
	<< e.what () << std::endl;
      return;
    }

    //increment the event counter
    ++_iEvt;
  }

  void EUTelNoisyPixelFinder::end ()
  {
    if (_finished)
      {
	streamlog_out (MESSAGE4) <<
	  "Noisy pixel finder has successfully finished!" << std::endl;
      }
    else
      {
	streamlog_out (ERROR3)
	  <<
	  "End of run reached before enough events were processed for noisy "
	  "pixel finder, databases have NOT been written!" << std::endl;
	streamlog_out (ERROR3) <<
	  "Increase the number of events to be processed for the run or "
	  "decrease the number required for noisy pixel finding." << std::
	  endl;
      }
  }

  void EUTelNoisyPixelFinder::check (LCEvent * /*event */ )
  {
    //only if the eventNo is the amount of events to be processed, we analyse
    //the data since check() runs after the event and we increment the event counter
    //before that the comparison is valid. If we set _noOfEvents=1 we only want to 
    //have processed event 0 before calling this. Since we increment 0++ before calling 
    //this function, this criteria is fullfilled
    if (_iEvt == _noOfEvents)
      {
	streamlog_out (MESSAGE4)
	  << "Finished determining hot pixels, writing them out..."
	  << std::endl;

	//[START] loop over all the sensors in sensorMap
      for (auto & thisSensor:_sensorMap)
	  {
	    auto sensorID = thisSensor.first;
	    streamlog_out (MESSAGE3) <<
	      "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	      "~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
	    streamlog_out (MESSAGE3) << "Noisy pixels found on plane " <<
	      sensorID << " (max. firing frequency set to: " <<
	      _maxAllowedFiringFreq << ")" << std::endl;
	    streamlog_out (MESSAGE3) <<
	      "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	      "~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

	    //get the corresponding hit array-like vector of vectors
	    std::vector < std::vector < long int >>*hitVector =
	      &_hitVecMap[sensorID];
	    //and the sensor which stores offsets
	    sensor *currentSensor = &_sensorMap[sensorID];
	    auto & firingFreqVec = _firingFreqForAllPixels[sensorID];

	    //[START] loop over all pixels
	    for (auto xIt = hitVector->begin (); xIt != hitVector->end ();
		 ++xIt)
	      {
		for (auto yIt = xIt->begin (); yIt != xIt->end (); ++yIt)
		  {
		    //compute the firing frequency
		    double fireFreq =
		      static_cast < double >(*yIt) / static_cast <
		      double >(_iEvt);
		    firingFreqVec.push_back (fireFreq);
		    //if larger than the allowed one, write pixel into a collection
		    if (fireFreq > _maxAllowedFiringFreq)
		      {
			streamlog_out (MESSAGE3)
			  << "Pixel: " << xIt - hitVector->begin () +
			  currentSensor->offX << "|" << yIt - xIt->begin () +
			  currentSensor->
			  offY << " fired " << fireFreq << std::endl;
			EUTelGenericSparsePixel pixel;
			pixel.setXCoord (xIt - hitVector->begin () +
					 currentSensor->offX);
			pixel.setYCoord (yIt - xIt->begin () +
					 currentSensor->offY);
			pixel.setSignal (fireFreq);
			//writing it out
			_noisyPixelMap[sensorID].push_back (pixel);
		      }
		  }
	      }			//[END] loop over pixel
	  }			//[END] loop over sensors

	//write out the databases and histograms
	noisyPixelDBWriter ();
	bookAndFillHistos ();
	_finished = true;
      }
  }

  void EUTelNoisyPixelFinder::noisyPixelDBWriter ()
  {

    streamlog_out (DEBUG5) << "Writing out hot pixel DB into "
      << _noisyPixelDBFile.c_str () << std::endl;

    //reopen LCIO file in append mode
    LCWriter *lcWriter = LCFactory::getInstance ()->createLCWriter ();
    std::unique_ptr < LCEventImpl > event (new LCEventImpl);
    //create new file if requested OR opening of existing file did not succeed
    try
    {
      lcWriter->open (_noisyPixelDBFile, LCIO::WRITE_NEW);
      LCRunHeaderImpl lcHeader;
      lcHeader.setRunNumber (0);
      lcWriter->writeRunHeader (&lcHeader);

      event->setRunNumber (0);
      event->setEventNumber (0);
      event->setDetectorName ("EUTelNoisyPixel");

      streamlog_out (DEBUG5) <<
	"hotPixelDB file: run header and event created" << std::endl;
      LCTime now;
      event->setTimeStamp (now.timeStamp ());
    } catch (IOException & e)
    {
      streamlog_out (ERROR4)
	<< e.what () << std::endl
	<< "Sorry, was not able to create new NoisyPixelDB file "
	<< _noisyPixelDBFile << ", will quit now. " << std::endl;
      return;
    }

    if (event == nullptr)
      {
	streamlog_out (ERROR5)
	  << "Problem opening NoisyPixelDB file, event is nullptr" << std::
	  endl;
	return;
      }

    //create main collection to be saved into the db file
    LCCollectionVec *noisyPixelCollection;
    //create new or open existing noisyPixel collection
    try
    {
      noisyPixelCollection =
	static_cast <
	LCCollectionVec * >(event->getCollection (_noisyPixelCollectionName));
      streamlog_out (MESSAGE4) << "noisyPixelCollection: " <<
	_noisyPixelCollectionName << " found with " << noisyPixelCollection->
	getNumberOfElements () << " elements " << std::endl;
      noisyPixelCollection->clear ();
      streamlog_out (MESSAGE4) << "noisyPixelCollection: " <<
	_noisyPixelCollectionName << " cleared: now " <<
	noisyPixelCollection->
	getNumberOfElements () << " elements " << std::endl;
    }
    catch (lcio::DataNotAvailableException & e)
    {
      noisyPixelCollection = new LCCollectionVec (lcio::LCIO::TRACKERDATA);
      event->addCollection (noisyPixelCollection, _noisyPixelCollectionName);
      streamlog_out (MESSAGE4) << "noisyPixelCollection: " <<
	_noisyPixelCollectionName << " created" << std::endl;
    }

    streamlog_out (MESSAGE5) << "Noisy Pixel Finder summary:" << std::endl;
    //_noisyPixelMap holds sensor id (first) and vector of noisy pixels (second)
  for (auto & mapEntry:_noisyPixelMap)
      {
	CellIDEncoder < TrackerDataImpl >
	  noisyPixelEncoder (EUTELESCOPE::ZSDATADEFAULTENCODING,
			     noisyPixelCollection);
	noisyPixelEncoder["sensorID"] = mapEntry.first;
	noisyPixelEncoder["sparsePixelType"] = kEUTelGenericSparsePixel;

	//prepare a new TrackerData for the hot pixel data
	std::unique_ptr < lcio::TrackerDataImpl >
	  currentFrame (new lcio::TrackerDataImpl);
	noisyPixelEncoder.setCellID (currentFrame.get ());

	//structure that will host the sparse pixel
	std::unique_ptr < EUTelTrackerDataInterfacerImpl <
	  EUTelGenericSparsePixel >> sparseFrame (new
						  EUTelTrackerDataInterfacerImpl
						  < EUTelGenericSparsePixel >
						  (currentFrame.get ()));

      for (auto & pixel:mapEntry.second)
	  {
	    sparseFrame->push_back (pixel);
	  }
	noisyPixelCollection->push_back (currentFrame.release ());

	streamlog_out (MESSAGE5) << "Found " << mapEntry.second.size ()
	  << " noisy pixels on sensor: " << mapEntry.first << std::endl;
      }
    lcWriter->writeEvent (event.get ());
    lcWriter->close ();
  }

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  void EUTelNoisyPixelFinder::bookAndFillHistos ()
  {

    streamlog_out (MESSAGE1) << "Booking histograms " << std::endl;

    //[START] loop of sensors
  for (auto det:_sensorIDVec)
      {

	sensor *currentSensor = &_sensorMap[det];

	//create folder for current detector
	std::string basePath = "detector_" + std::to_string (det);
	marlin::AIDAProcessor::tree (this)->mkdir (basePath.c_str ());

	//create scatter plot
	std::string dataPointName =
	  basePath + "/FiringFrequencyDependency_det" + std::to_string (det);
	AIDA::ICloud1D * dataPointSet =
	  marlin::AIDAProcessor::histogramFactory (this)->
	  createCloud1D (dataPointName);

	//create 1D histogram: noisy pixel vs noise cut           
	std::string histName_noisyPixelVsNoiseCut =
	  basePath + "/NoisyPixel_vs_NoiseCut_det" + std::to_string (det);
	int nbin = 100;
	double cutLow = 0;
	double cutHigh = 0.006;
	AIDA::IHistogram1D * hist1D_noisyPixelVsNoiseCut =
	  marlin::AIDAProcessor::
	  histogramFactory (this)->createHistogram1D
	  (histName_noisyPixelVsNoiseCut, nbin, cutLow, cutHigh);
	hist1D_noisyPixelVsNoiseCut->
	  setTitle
	  ("Number of noisy pixels for given noise cut; noise cut; #noisy pixels");

	//create 1D histogram: firing frequency
	std::string histName_firingFreq1D =
	  basePath + "/FiringFrequency1D_det" + std::to_string (det);
	AIDA::IHistogram1D * hist1D_firingFreq =
	  marlin::AIDAProcessor::
	  histogramFactory (this)->createHistogram1D (histName_firingFreq1D,
						      100, 0, 100);
	hist1D_firingFreq->
	  setTitle
	  ("Firing frequency distribution of hot pixels; Firing Frequency (%); Count (#)");

	//create 2D histogram: firing frequency
	std::string histName_firingFreq2D =
	  basePath + "/FiringFrequency2D_det";
	int xBin = currentSensor->sizeX + 1;
	double xMin = static_cast < double >(currentSensor->offX) - 0.5;
	double xMax =
	  static_cast <
	  double >(currentSensor->offX + currentSensor->sizeX) + 0.5;
	int yBin = currentSensor->sizeY + 1;
	double yMin = static_cast < double >(currentSensor->offY) - 0.5;
	double yMax =
	  static_cast <
	  double >(currentSensor->offY + currentSensor->sizeY) + 0.5;
	AIDA::IHistogram2D * hist2D_firingFreq =
	  marlin::AIDAProcessor::
	  histogramFactory (this)->createHistogram2D (histName_firingFreq2D,
						      xBin, xMin, xMax, yBin,
						      yMin, yMax);
	hist2D_firingFreq->
	  setTitle
	  ("Firing frequency map of hot pixels; Pixel Index X; Pixel Index Y; Percent (%)");

	//fill vectors for firing frequency and noise cut values
	auto & firingFreqVec = _firingFreqForAllPixels[det];
	std::sort (firingFreqVec.begin (), firingFreqVec.end (),
		   [](const double a, const double b)
		   {
		   return a > b;
		   }
	);
	std::vector < long double >cuts;
	long double cutsteps = cutHigh / static_cast < long double >(nbin);
	for (int ibin = 1; ibin <= nbin * 10; ++ibin)
	  {
	    cuts.emplace_back (ibin * cutsteps);
	  }

	//fill dataPointSet with noisy pixels in dependence of noise cut
	long counter = 0;
	for (auto it = firingFreqVec.begin (); it != firingFreqVec.end ();
	     ++it)
	  {
	    if (cuts.back () >= *it)
	      {
		while (cuts.back () >= *it)
		  {
		    dataPointSet->fill (cuts.back (), counter);
		    cuts.pop_back ();
		  }
	      }
	    counter++;
	  }
	dataPointSet->fillHistogram (*hist1D_noisyPixelVsNoiseCut);

	//fill firing frequency histograms      
      for (auto & pixel:_noisyPixelMap[det])
	  {
	    hist1D_firingFreq->fill (pixel.getSignal ());
	    hist2D_firingFreq->fill (pixel.getXCoord (), pixel.getYCoord (),
				     pixel.getSignal ());
	  }
      }				//[END] loop over detectors
  }
#endif

}				//namespace eutelescope
