/*
 *   Most of this code is based on the EUTelClusteringProcessor by
 *   Antonio Bulgheroni
 *
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

//eutelescope includes
#include "EUTelProcessorGeometricClustering.h"

#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelHistogramManager.h"

//eutel data specific
#include "EUTelTrackerDataInterfacerImpl.h"
#include "EUTelGenericSparseClusterImpl.h"

//eutel geometry
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelGenericPixGeoDescr.h"

//ROOT includes
#include "TGeoShape.h"
#include "TGeoBBox.h"

//marlin includes
#include "marlin/Processor.h"
#include "marlin/AIDAProcessor.h"
#include "marlin/Exceptions.h"

//lcio includes
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/LCCollectionVec.h>

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
//aida includes
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ITree.h>
#endif

//system includes
#include <string>
#include <vector>
#include <memory>
//#include <iostream>
#include <cmath>

using namespace lcio;
using namespace marlin;
using namespace eutelescope;

EUTelProcessorGeometricClustering::EUTelProcessorGeometricClustering(): 
  Processor("EUTelProcessorGeometricClustering"), 
  _zsDataCollectionName(""),
  _pulseCollectionName(""),
  _initialPulseCollectionSize(0),
  _iRun(0),
  _iEvt(0),
  _fillHistos(false),
  _histoInfoFileName(""),
  _cutT(0.0),
  _totClusterMap(),
  _noOfDetector(0),
  _ExcludedPlanes(),
  _clusterSignalHistos(),
  _clusterSizeXHistos(),
  _clusterSizeYHistos(),
  _seedSignalHistos(),
  _hitMapHistos(),
  _eventMultiplicityHistos(),
  _isGeometryReady(false),
  _sensorIDVec(),
  _zsInputDataCollectionVec(NULL),
  _pulseCollectionVec(NULL)
 {
  
  // modify processor description
  _description = "EUTelProcessorGeometricClustering is looking for clusters into a calibrated pixel matrix.";

  // first of all we need to register the input collection
  registerInputCollection (LCIO::TRACKERDATA, "ZSDataCollectionName", "Input of Zero Suppressed data",
                           _zsDataCollectionName, std::string ("zsdata") );

  registerOutputCollection(LCIO::TRACKERPULSE, "PulseCollectionName", "Cluster (output) collection name",
                           _pulseCollectionName, std::string("cluster"));

  // now the optional parameters
  registerProcessorParameter("TCut","Time cut in time units of your sensor",
                             _cutT, static_cast<float > ( std::numeric_limits<float>::max() ));

  registerProcessorParameter("HistoInfoFileName", "This is the name of the histogram information file",
                             _histoInfoFileName, std::string( "histoinfo.xml" ) );

  registerProcessorParameter("HistogramFilling","Switch on or off the histogram filling",
                             _fillHistos, static_cast< bool > ( true ) );

  registerOptionalParameter("ExcludedPlanes", "The list of sensor ids that have to be excluded from the clustering.",
                             _ExcludedPlanes, std::vector<int> () );

  		_isFirstEvent = true;
}

void EUTelProcessorGeometricClustering::init() {
	//this method is called only once even when the rewind is active, it is usually a good idea to
	printParameters ();

	//init new geometry
	std::string name("test.root");
	geo::gGeometry().initializeTGeoDescription(name,true);

	//set to zero the run and event counters
	_iRun = 0;
	_iEvt = 0;

	//the geometry is not yet initialized, so set the corresponding switch to false
	_isGeometryReady = false;
}

void EUTelProcessorGeometricClustering::processRunHeader (LCRunHeader * rdr) {

	std::auto_ptr<EUTelRunHeaderImpl> runHeader( new EUTelRunHeaderImpl( rdr ) );
	runHeader->addProcessor( type() );
  	//increment the run counter
	++_iRun;
}

void EUTelProcessorGeometricClustering::initializeGeometry( LCEvent * event ) throw ( marlin::SkipEventException ) {

	//set the total number of detector to zero. This number can be different from the one written in the gear description because
	//the input collection can contain only a fraction of all the sensors.

	_noOfDetector = 0;
	_sensorIDVec.clear();

	streamlog_out( DEBUG5 ) << "Initializing geometry" << std::endl;

  	try 
  	{
		_zsInputDataCollectionVec = dynamic_cast<LCCollectionVec*>( event->getCollection(_zsDataCollectionName) );
		_noOfDetector += _zsInputDataCollectionVec->getNumberOfElements();
		CellIDDecoder<TrackerDataImpl> cellDecoder(_zsInputDataCollectionVec);

		for ( size_t i = 0; i < _zsInputDataCollectionVec->size(); ++i ) 
		{
			TrackerDataImpl * data = dynamic_cast< TrackerDataImpl * > ( _zsInputDataCollectionVec->getElementAt( i ) ) ;
			_sensorIDVec.push_back( cellDecoder(data)["sensorID"] );
			_totClusterMap.insert( std::make_pair( cellDecoder(data)["sensorID"], 0) );
		}
	} 

	catch ( lcio::DataNotAvailableException ) 
	{
		streamlog_out( DEBUG5 ) << "Could not find the input collection: " << _zsDataCollectionName.c_str() << " !" << std::endl;
		return;
	}

    _isGeometryReady = true;
}

void EUTelProcessorGeometricClustering::modifyEvent( LCEvent * /* event */ )
{
	return;
}

void EUTelProcessorGeometricClustering::readCollections(LCEvent* event)
{
	try 
	{
		_zsInputDataCollectionVec = dynamic_cast< LCCollectionVec * > ( event->getCollection( _zsDataCollectionName ) ) ;
        	streamlog_out ( DEBUG4 ) << "_zsInputDataCollectionVec: " << _zsDataCollectionName.c_str() << " found " << std::endl;
	} 
	catch ( lcio::DataNotAvailableException )   // do nothing
	{
		streamlog_out ( DEBUG4 ) << "_zsInputDataCollectionVec: " << _zsDataCollectionName.c_str() << " not found " << std::endl;
	}

	try 
	{
		event->getCollection( _zsDataCollectionName ) ;
	} 
	catch (lcio::DataNotAvailableException& e ) 
	{
		streamlog_out(MESSAGE2) << "The current event doesn't contain nZS data collections: skip # " << event->getEventNumber() << std::endl;
		throw SkipEventException(this);
	}
}

void EUTelProcessorGeometricClustering::processEvent (LCEvent * event) 
{
	//increment event counter
	++_iEvt;

	//first of all we need to be sure that the geometry is properly initialized!
	if ( !_isGeometryReady ) 
	{
		initializeGeometry( event );
	}

	//read zsInputData collection
	readCollections(event);

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
	// book the histograms now
	if ( _fillHistos && isFirstEvent() ) 
	{
		bookHistos();
	}
#endif

	EUTelEventImpl* evt = static_cast<EUTelEventImpl*> (event);
	if ( evt->getEventType() == kEORE ) 
  	{
		streamlog_out ( DEBUG4 ) <<  "EORE found: nothing else to do." <<  std::endl;
		return;
	}
	else if ( evt->getEventType() == kUNKNOWN ) 
	{
		streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " is of unknown type. Continue considering it as a normal Data Event." << std::endl;
	}

	// prepare a pulse collection to add all clusters found this can be either a new collection or already existing in the event
	LCCollectionVec* pulseCollection;
	bool pulseCollectionExists = false;
	_initialPulseCollectionSize = 0;
	try 
	{
		pulseCollection = dynamic_cast< LCCollectionVec * > ( evt->getCollection( _pulseCollectionName ) );
		pulseCollectionExists = true;
		_initialPulseCollectionSize = pulseCollection->size();
	} 
	catch ( lcio::DataNotAvailableException& e ) 
	{
		pulseCollection = new LCCollectionVec(LCIO::TRACKERPULSE);
	}

	//HERE WE ACTUALLY CALL THE CLUSTERING ROUTINE:
	geometricClustering(evt, pulseCollection);

	// if the pulseCollection is not empty add it to the event
	if ( ! pulseCollectionExists && ( pulseCollection->size() != _initialPulseCollectionSize )) 
	{
		evt->addCollection( pulseCollection, _pulseCollectionName );
	}
  
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
	if ( (pulseCollection->size() != _initialPulseCollectionSize) && (_fillHistos) ) 
	{
		fillHistos(event);
	}
#endif

	if ( ! pulseCollectionExists && ( pulseCollection->size() == _initialPulseCollectionSize ) ) 
	{
		delete pulseCollection;
	}
	_isFirstEvent = false;
}

void EUTelProcessorGeometricClustering::geometricClustering(LCEvent * evt, LCCollectionVec * pulseCollection) 
{
	// prepare some decoders
	CellIDDecoder<TrackerDataImpl> cellDecoder( _zsInputDataCollectionVec );

	bool isDummyAlreadyExisting = false;
	LCCollectionVec* sparseClusterCollectionVec = NULL;

	try 
	{
		sparseClusterCollectionVec = dynamic_cast< LCCollectionVec* > ( evt->getCollection( "original_zsdata") );
		isDummyAlreadyExisting = true ;
	} 
	catch (lcio::DataNotAvailableException& e) 
	{
		sparseClusterCollectionVec =  new LCCollectionVec(LCIO::TRACKERDATA);
		isDummyAlreadyExisting = false;
	}

	CellIDEncoder<TrackerDataImpl> idZSClusterEncoder( EUTELESCOPE::ZSCLUSTERDEFAULTENCODING, sparseClusterCollectionVec  );

	// prepare an encoder also for the pulse collection
	CellIDEncoder<TrackerPulseImpl> idZSPulseEncoder(EUTELESCOPE::PULSEDEFAULTENCODING, pulseCollection);

	// in the _zsInputDataCollectionVec we should have one TrackerData for each 
	// detector working in ZS mode. We need to loop over all of them
	for ( unsigned int idetector = 0 ; idetector < _zsInputDataCollectionVec->size(); idetector++ ) 
	{
		// get the TrackerData and guess which kind of sparsified data it contains.
		TrackerDataImpl * zsData = dynamic_cast< TrackerDataImpl * > ( _zsInputDataCollectionVec->getElementAt( idetector ) );
		SparsePixelType   type   = static_cast<SparsePixelType> ( static_cast<int> (cellDecoder( zsData )["sparsePixelType"]) );
		int sensorID             = static_cast<int > ( cellDecoder( zsData )["sensorID"] );
    
		//get alle the plane relevant geo information, that is the plane name and the plane pix geometry
		std::string planePath = geo::gGeometry().getPlanePath( sensorID );
		geo::EUTelGenericPixGeoDescr* geoDescr =  ( geo::gGeometry().getPixGeoDescr( sensorID ) );

		//if this is an excluded sensor go to the next element
		bool foundexcludedsensor = false;
		for(size_t iexclude = 0; iexclude < _ExcludedPlanes.size(); ++iexclude)
		{
		    if(_ExcludedPlanes[iexclude] == sensorID)
		    {
		        foundexcludedsensor = true;
		    }
		}

		if(foundexcludedsensor)       
		{		
			continue;
		}

		//now that we know which is the sensorID, we can ask which are the minX, minY, maxX and maxY.
		int minX, minY, maxX, maxY;
		minX = minY = maxX = maxY = 0;
		geoDescr->getPixelIndexRange( minX, maxX, minY, maxY );

    		if ( type == kEUTelGenericSparsePixel ) 
		{

			// now prepare the EUTelescope interface to sparsified data.
			std::auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel > > sparseData( new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> ( zsData ) );

			streamlog_out ( DEBUG2 ) << "Processing sparse data on detector " << sensorID << " with " << sparseData->size() << " pixels " << std::endl;

			int hitPixelsInEvent = sparseData->size();
			std::vector<EUTelGeometricPixel> hitPixelVec;
			EUTelGenericSparsePixel* genericPixel = new EUTelGenericSparsePixel;

			//This for-loop loads all the hits of the given event and detector plane and stores them as GeometricPixels
			for(int i = 0; i < hitPixelsInEvent; ++i )
			{
				//Load the information of the hit pixel into genericPixel
				sparseData->getSparsePixelAt( i, genericPixel );
				EUTelGeometricPixel hitPixel( *genericPixel );

				//And get the path to the given pixel
				std::string pixelPath = geoDescr->getPixName(hitPixel.getXCoord(), hitPixel.getYCoord());

				//Then navigate to this pixel with the TGeo manager
				geo::gGeometry()._geoManager->cd( (planePath+pixelPath).c_str() );

				//get the imbedding box
				TGeoShape* currentShape =  geo::gGeometry()._geoManager->GetCurrentVolume()->GetShape();
				TGeoBBox* bbox = dynamic_cast<TGeoBBox*>( currentShape );
				//store the dimensions of this box in the GeometricPixel
				hitPixel.setBoundaryX( bbox->GetDX() );
				hitPixel.setBoundaryY( bbox->GetDY() );

				//Get how deep the node description goes (this is how often we have to transform to get coordinates in the local plane coordinate system)
				std::vector<std::string> split = Utility::stringSplit( planePath+pixelPath , "/", false);

				//Three recursions for the telescope/plane
				int recursionDepth = split.size() - 3;

				//The do the transformation
				Double_t origin_pt[3] = {0,0,0};
				Double_t transformed1_pt[3];
				Double_t transformed2_pt[3];
				gGeoManager->GetCurrentNode()->LocalToMaster(origin_pt, transformed1_pt);

				transformed2_pt[0] = transformed1_pt[0];
				transformed2_pt[1] = transformed1_pt[1];
				transformed2_pt[2] = transformed1_pt[2];

				//transform into local plane coordinate system
				for(int i = 1 ; i < recursionDepth; ++i)
				{
					gGeoManager->GetMother(i)->LocalToMaster(transformed1_pt, transformed2_pt);
					transformed1_pt[0] = transformed2_pt[0];
					transformed1_pt[1] = transformed2_pt[1];
					transformed1_pt[2] = transformed2_pt[2];
				}

				//store all the position information in the GeometricPixel
				hitPixel.setPosX( transformed2_pt[0] );
				hitPixel.setPosY( transformed2_pt[1] );
				//and push this pixel back
				hitPixelVec.push_back( hitPixel );
			}		

			std::vector<EUTelGeometricPixel> newlyAdded;
			//We now cluster those hits together
			while( !hitPixelVec.empty() )
			{
				// prepare a TrackerData to store the cluster candidate
				std::auto_ptr< TrackerDataImpl > zsCluster ( new TrackerDataImpl );
				// prepare a reimplementation of sparsified cluster
				std::auto_ptr<EUTelGenericSparseClusterImpl<EUTelGeometricPixel > > sparseCluster ( new EUTelGenericSparseClusterImpl<EUTelGeometricPixel >( zsCluster.get() ) );

				//First we need to take any pixel, so let's take the first one
				//Add it to the cluster as well as the newly added pixels
				newlyAdded.push_back( hitPixelVec.front() );
				sparseCluster->addSparsePixel( &hitPixelVec.front() );
				//And remove it from the original collection
				hitPixelVec.erase( hitPixelVec.begin() );
		
				//Now process all newly added pixels, initially this is the just previously added one
				//but in the process of neighbour finding we continue to add new pixels
				while( !newlyAdded.empty() )
				{
					bool newlyDone = true;
					float x1, x2, y1, y2, dX, dY, cx1, cy1, cx2, cy2, cutX, cutY, t1 , t2, dT;

					//check against all pixels in the hitPixelVec
					for( std::vector<EUTelGeometricPixel>::iterator hitVec = hitPixelVec.begin(); hitVec != hitPixelVec.end(); ++hitVec )
					{
						//get the relevant infos from the newly added pixel
						x1 = newlyAdded.front().getPosX();
						y1 = newlyAdded.front().getPosY();
						t1 =  newlyAdded.front().getTime();
						cx1 = newlyAdded.front().getBoundaryX();
						cy1 = newlyAdded.front().getBoundaryY();

						//and the pixel we test against
						x2 = hitVec->getPosX();
						y2 = hitVec->getPosY();
						t2 = hitVec->getTime();
						cx2 = hitVec->getBoundaryX();
						cy2 = hitVec->getBoundaryY();

						dX = x1 - x2;
						dY = y1 - y2;
						dT = t1 - t2;
						cutX = (cx1+cx2)*1.01; //this additional 1% is accounting for precision
						cutY = (cy1+cy2)*1.01; //uncertainty with the geo framework

						//if they pass the spatial and temporal cuts, we add them	
						if(	(dX*dX <= cutX*cutX) && (dY*dY <= cutY*cutY) && (dT*dT <= _cutT*_cutT) )
						{
							//add them to the cluster as well as to the newly added ones
							newlyAdded.push_back( *hitVec );
							sparseCluster->addSparsePixel( &(*hitVec) );
							//and remove it from the original collection
							hitPixelVec.erase( hitVec );
							//for the pixel we test there might be other neighbours, we still have to check
							newlyDone = false;
							break;
						}
					}
				
					//if no neighbours are found, we can delete the pixel from the newly added
					//we tested against _ALL_ non cluster pixels, there are no other pixels
					//which could be neighbours
					if(newlyDone) newlyAdded.erase( newlyAdded.begin() );
				}	

				//Now we need to process the found cluster
				if ( sparseCluster->size() > 0 ) 
				{
					// set the ID for this zsCluster
					idZSClusterEncoder["sensorID"]  = sensorID;
					idZSClusterEncoder["sparsePixelType"] = static_cast<int>( type );
					idZSClusterEncoder["quality"] = 0;
					idZSClusterEncoder.setCellID( zsCluster.get() );

					// add it to the cluster collection
					sparseClusterCollectionVec->push_back( zsCluster.get() );

					//int xSeed, ySeed, xSize, ySize;
					//sparseCluster->getClusterInfo(xSeed, ySeed, xSize, ySize);

					// prepare a pulse for this cluster
					std::auto_ptr<TrackerPulseImpl> zsPulse ( new TrackerPulseImpl );
					idZSPulseEncoder["sensorID"]  = sensorID;
					//idZSPulseEncoder["xSeed"]     = xSeed;
					//idZSPulseEncoder["ySeed"]     = ySeed;
					idZSPulseEncoder["type"]      = static_cast<int>(kEUTelGenericSparseClusterImpl);
					idZSPulseEncoder.setCellID( zsPulse.get() );

					zsPulse->setCharge( sparseCluster->getTotalCharge() );
					//zsPulse->setQuality( static_cast<int > (sparseCluster->getClusterQuality()) );
					zsPulse->setTrackerData( zsCluster.release() );
					pulseCollection->push_back( zsPulse.release() );

					// last but not least increment the totClusterMap
					_totClusterMap[ sensorID ] += 1;

				} //cluster processing if

				else 
				{
					//in the case the cluster candidate is not passing the threshold ...
					//forget about them, the memory should be automatically cleaned by std::auto_ptr's
				}
			} //loop over all found clusters
			
			delete genericPixel;
    		}	 
		else 
		{
    			throw UnknownDataTypeException("Unknown sparsified pixel");
    		}
	} // this is the end of the loop over all ZS detectors

	// if the sparseClusterCollectionVec isn't empty add it to the
	// current event. The pulse collection will be added afterwards
	if ( ! isDummyAlreadyExisting ) 
	{
		if ( sparseClusterCollectionVec->size() != 0 ) 
		{
			evt->addCollection( sparseClusterCollectionVec, "original_zsdata" );
		}
		else 
		{
			delete sparseClusterCollectionVec;
		}
	}
}

void EUTelProcessorGeometricClustering::check (LCEvent * /* evt */) {
  // nothing to check here - could be used to fill check plots in reconstruction processor
}

void EUTelProcessorGeometricClustering::end() {
	
	streamlog_out ( MESSAGE4 ) <<  "Successfully finished" << std::endl;
  
	std::map<int,int>::iterator iter = _totClusterMap.begin();
	while( iter != _totClusterMap.end() )
	{
		streamlog_out ( MESSAGE4 ) << "Found " << iter->second << " clusters on detector " << iter->first << std::endl;
		++iter;
	}
}

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
void EUTelProcessorGeometricClustering::fillHistos (LCEvent * evt) 
{
	EUTelEventImpl* eutelEvent = static_cast<EUTelEventImpl*> (evt);
	EventType type = eutelEvent->getEventType();

	if ( type == kEORE ) 
	{
	    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << std::endl;
	    return;
	}
	else if ( type == kUNKNOWN ) 
	{
		// if it is unknown we had already issued a warning to the user at
		// the beginning of the processEvent. If we get to here, it means
		// that the assumption that the event was a data event was
		// correct, so no harm to continue...
	}

	try 
	{
		LCCollectionVec* _pulseCollectionVec = dynamic_cast<LCCollectionVec*>  (evt->getCollection(_pulseCollectionName));
		CellIDDecoder<TrackerPulseImpl > cellDecoder(_pulseCollectionVec);

		std::map<int, int> eventCounterMap;

		for( int iPulse = _initialPulseCollectionSize; iPulse < _pulseCollectionVec->getNumberOfElements(); iPulse++ ) 
		{
			TrackerPulseImpl* pulse = dynamic_cast<TrackerPulseImpl*> ( _pulseCollectionVec->getElementAt(iPulse) );
			ClusterType type  = static_cast<ClusterType> ( static_cast<int> ( cellDecoder(pulse)["type"] ));
			int detectorID = static_cast<int>( cellDecoder(pulse)["sensorID"] ); 
			//TODO: do we need this check?
			//SparsePixelType pixelType = static_cast<SparsePixelType> (0);
			EUTelSimpleVirtualCluster* cluster;
	
			if( type == kEUTelGenericSparseClusterImpl ) 
			{
		    	cluster = new EUTelGenericSparseClusterImpl<EUTelGeometricPixel> ( static_cast<TrackerDataImpl*> ( pulse->getTrackerData() ) );
			}	 
			else 
			{
			    streamlog_out ( ERROR4 ) <<  "Unknown cluster type. Sorry for quitting" << std::endl;
			    throw UnknownDataTypeException("Cluster type unknown");
			}
	
			//if this key doesn't exist yet it will be value initialized, this is desired, for int this is 0!
			eventCounterMap[detectorID]++;

			//if this is an excluded sensor go to the next element
			bool foundexcludedsensor = false;
			for(size_t i = 0; i < _ExcludedPlanes.size(); ++i)
		    {
				if(_ExcludedPlanes[i] == detectorID)
		        {
					foundexcludedsensor = true;
					break;
		        }
		    }

			if(foundexcludedsensor) continue;

			// get the cluster size in X and Y separately and plot it:
			int xPos, yPos, xSize, ySize;
			cluster->getClusterInfo(xPos, yPos, xSize, ySize);
			float geoPosX, geoPosY, geoSizeX, geoSizeY;
			cluster->getClusterGeomInfo(geoPosX, geoPosY, geoSizeX, geoSizeY);

			//Do all the plots
			(dynamic_cast<AIDA::IHistogram1D*> (_clusterSizeXHistos[detectorID]))->fill(xSize);
			(dynamic_cast<AIDA::IHistogram1D*> (_clusterSizeYHistos[detectorID]))->fill(ySize);
			(dynamic_cast<AIDA::IHistogram2D*> (_hitMapHistos[detectorID]))->fill(static_cast<double >(xPos), static_cast<double >(yPos), 1.);
			(dynamic_cast<AIDA::IHistogram2D*> (_hitMapGeomHistos[detectorID]))->fill(geoPosX, geoPosY, 1.);
			(dynamic_cast<AIDA::IHistogram1D*> (_clusterSizeTotalHistos[detectorID]))->fill( static_cast<int>(cluster->size()) );
			(dynamic_cast<AIDA::IHistogram1D*> (_clusterSignalHistos[detectorID]))->fill(cluster->getTotalCharge());

			delete cluster;
		}

		//fill the event multiplicity here
		std::string tempHistoName;
		for ( int iDetector = 0; iDetector < _noOfDetector; iDetector++ ) 
		{
			AIDA::IHistogram1D * histo = dynamic_cast<AIDA::IHistogram1D*> ( _eventMultiplicityHistos[_sensorIDVec.at( iDetector)] );
			if ( histo ) 
			{
			    histo->fill( eventCounterMap[_sensorIDVec.at( iDetector)] );
			}
		}
	} 
	catch (lcio::DataNotAvailableException& e) 
	{
		return;
	}
}
#endif

#ifdef MARLIN_USE_AIDA
void EUTelProcessorGeometricClustering::bookHistos() {

  // histograms are grouped in loops and detectors
  streamlog_out ( DEBUG5 )  << "Booking histograms " << std::endl;
  std::auto_ptr<EUTelHistogramManager> histoMgr( new EUTelHistogramManager( _histoInfoFileName ) );
  EUTelHistogramInfo* histoInfo;
  bool isHistoManagerAvailable;

  try {
    isHistoManagerAvailable = histoMgr->init();
  } catch ( std::ios::failure& e) {
    streamlog_out ( WARNING2 ) << "I/O problem with " << _histoInfoFileName << "\n"
                               << "Continuing without histogram manager"  << std::endl;
    isHistoManagerAvailable = false;
  } catch ( ParseException& e ) {
    streamlog_out ( WARNING2 ) << e.what() << "\n"
                               << "Continuing without histogram manager" << std::endl;
    isHistoManagerAvailable = false;
  }

  // define our histogram names
  std::string _clusterSignalHistoName      = "clusterSignal";
  std::string _clusterSizeXHistoName       = "clusterSizeX";
  std::string _clusterSizeYHistoName       = "clusterSizeY";
  std::string _clusterSizeTotalHistoName   = "totalClusterSize";
  std::string _hitMapHistoName             = "hitMap";
  std::string _hitMapGeomHistoName         = "geometricHitMap";
  std::string _eventMultiplicityHistoName  = "eventMultiplicity";

  std::string tempHistoName;
  std::string basePath;

	for (size_t iDetector = 0; iDetector < _sensorIDVec.size(); iDetector++) 
	{
		int sensorID = _sensorIDVec.at( iDetector );
		geo::EUTelGenericPixGeoDescr* geoDescr =  ( geo::gGeometry().getPixGeoDescr( sensorID ) );

		//now that we know which is the sensorID, we can ask which are the minX, minY, maxX and maxY.
		int minX, minY, maxX, maxY;
		minX = minY = maxX = maxY = 0;
		float sizeX, sizeY, sizeZ;
		sizeX=sizeY=sizeZ=0;
		 
		geoDescr->getPixelIndexRange( minX, maxX, minY, maxY );
		geoDescr->getSensitiveSize(sizeX, sizeY, sizeZ);

		basePath = "detector_" + to_string( sensorID );
		AIDAProcessor::tree(this)->mkdir(basePath.c_str());
		basePath.append("/");

		int    clusterNBin  = 61;
		double clusterMin   = -0.5;
		double clusterMax   = 60.5;
		std::string clusterTitle = "Total Signal per Cluster (in Detector specific Charge Unit);Charge;Count [#]";
		if ( isHistoManagerAvailable ) 
		{
			histoInfo = histoMgr->getHistogramInfo( _clusterSizeTotalHistoName );
			if ( histoInfo ) 
			{
				streamlog_out ( DEBUG2 ) << (* histoInfo ) << std::endl;
				clusterNBin = histoInfo->_xBin;
				clusterMin  = histoInfo->_xMin;
				clusterMax  = histoInfo->_xMax;
				if ( histoInfo->_title != "" ) clusterTitle = histoInfo->_title;
			}
		}

		int    clusterTotBin  = 61;
		double clusterTotMin   = -0.5;
		double clusterTotMax   = 60.5;
		std::string clusterTotTitle = "Total size of Cluster (in Hit Pixels);No. of Pixel in Cluster [#];Count [#]";
		if ( isHistoManagerAvailable ) {
		  histoInfo = histoMgr->getHistogramInfo( _clusterSignalHistoName );
		  if ( histoInfo ) {
		    streamlog_out ( DEBUG2 ) << (* histoInfo ) << std::endl;
		    clusterTotBin = histoInfo->_xBin;
		    clusterTotMin  = histoInfo->_xMin;
		    clusterTotMax  = histoInfo->_xMax;
		    if ( histoInfo->_title != "" ) clusterTotTitle = histoInfo->_title;
		  }
		}

		int    clusterXNBin  = 21;
		double clusterXMin   = -0.5;
		double clusterXMax   = 20.5;
		std::string clusterXTitle = "Size in X-Direction of all Clusters in Pixels;Size [# of pixels];Count [#]";
		if ( isHistoManagerAvailable ) {
		  histoInfo = histoMgr->getHistogramInfo( _clusterSizeXHistoName );
		  if ( histoInfo ) {
		    streamlog_out ( DEBUG2 ) << (* histoInfo ) << std::endl;
		    clusterXNBin = histoInfo->_xBin;
		    clusterXMin  = histoInfo->_xMin;
		    clusterXMax  = histoInfo->_xMax;
		    if ( histoInfo->_title != "" ) clusterXTitle = histoInfo->_title;
		  }
		}


		int    clusterYNBin  = 21;
		double clusterYMin   = -0.5;
		double clusterYMax   = 20.5;
		std::string clusterYTitle = "Size in Y-Direction of all Clusters in Pixels;Size [# of pixels];Count [#]";
		if ( isHistoManagerAvailable ) {
		  histoInfo = histoMgr->getHistogramInfo( _clusterSizeYHistoName );
		  if ( histoInfo ) {
		    streamlog_out ( DEBUG2 ) << (* histoInfo ) << std::endl;
		    clusterYNBin = histoInfo->_xBin;
		    clusterYMin  = histoInfo->_xMin;
		    clusterYMax  = histoInfo->_xMax;
		    if ( histoInfo->_title != "" ) clusterYTitle = histoInfo->_title;
		  }
		}

		// cluster total size
		tempHistoName = _clusterSizeTotalHistoName + "_d" + to_string( sensorID );
		_clusterSizeTotalHistos.insert(std::make_pair(sensorID, 
						  AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
														clusterTotBin,clusterTotMin,clusterTotMax)));
		_clusterSizeTotalHistos[sensorID]->setTitle(clusterTotTitle.c_str());

		// cluster signal
		tempHistoName = _clusterSignalHistoName + "_d" + to_string( sensorID );
		_clusterSignalHistos.insert(std::make_pair(sensorID, 
						  AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
														clusterNBin,clusterMin,clusterMax)));
		_clusterSignalHistos[sensorID]->setTitle(clusterTitle.c_str());

		// cluster signal along X
		tempHistoName = _clusterSizeXHistoName + "_d" + to_string( sensorID );
		_clusterSizeXHistos.insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), clusterXNBin,clusterXMin,clusterXMax)));
		_clusterSizeXHistos[sensorID]->setTitle(clusterXTitle.c_str());

		// cluster signal along Y
		tempHistoName = _clusterSizeYHistoName + "_d" + to_string( sensorID );
		_clusterSizeYHistos.insert(std::make_pair(sensorID, 
						AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
													  clusterYNBin,clusterYMin,clusterYMax)));
		_clusterSizeYHistos[sensorID]->setTitle(clusterYTitle.c_str());



		tempHistoName = _hitMapHistoName + "_d" + to_string( sensorID );
		int     xBin = maxX - minX + 1;
		double  xMin = static_cast<double >( minX ) - 0.5 ;
		double  xMax = static_cast<double >( maxX ) + 0.5;
		int     yBin = maxY - minY + 1;
		double  yMin = static_cast<double >( minY ) - 0.5;
		double  yMax = static_cast<double >( maxY ) + 0.5;
		AIDA::IHistogram2D * hitMapHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(), xBin, xMin, xMax,yBin, yMin, yMax);
		_hitMapHistos.insert(std::make_pair(sensorID, hitMapHisto));
		hitMapHisto->setTitle("Pixel Index Hit Map;X Index [#];Y Index [#];Count [#]");


		float binX = (maxX-minX+0)*60/sizeX;
		float binY = (maxY-minY+0)*40/sizeY;

		tempHistoName = _hitMapGeomHistoName + "_d" + to_string( sensorID );
		int     xGeomBin = ceil(binX);
		double  xGeomMin = -30.25;
		double  xGeomMax =  29.75;
		int     yGeomBin = ceil(binY);
		double  yGeomMin = -20.25;
		double  yGeomMax = 19.75;
		AIDA::IHistogram2D * hitMapGeomHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(), xGeomBin, xGeomMin, xGeomMax, yGeomBin, yGeomMin, yGeomMax);
		_hitMapGeomHistos.insert(std::make_pair(sensorID, hitMapGeomHisto));
		hitMapGeomHisto->setTitle("Geometric Cluster Hit Map;X-Position [mm];Y-Position [mm];Count [#]");

		tempHistoName = _eventMultiplicityHistoName + "_d" + to_string( sensorID );
		int     eventMultiNBin  = 15;
		double  eventMultiMin   =  -0.5;
		double  eventMultiMax   = 14.5;
		std::string  eventMultiTitle = "Event multiplicity";
		if ( isHistoManagerAvailable ) {
		  histoInfo = histoMgr->getHistogramInfo(  _eventMultiplicityHistoName );
		  if ( histoInfo ) {
		    streamlog_out ( DEBUG2 ) << (* histoInfo ) << std::endl;
		    eventMultiNBin  = histoInfo->_xBin;
		    eventMultiMin   = histoInfo->_xMin;
		    eventMultiMax   = histoInfo->_xMax;
		    if ( histoInfo->_title != "" ) eventMultiTitle = histoInfo->_title;
		  }
		}
		AIDA::IHistogram1D * eventMultiHisto =
		  AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
		                                                            eventMultiNBin, eventMultiMin, eventMultiMax);
		_eventMultiplicityHistos.insert( std::make_pair(sensorID, eventMultiHisto) );
		eventMultiHisto->setTitle( eventMultiTitle.c_str() );

  }
  streamlog_out ( DEBUG5 )  << "end of Booking histograms " << std::endl; 
}
#endif
