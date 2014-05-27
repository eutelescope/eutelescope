// Version: $Id$
/*========================================================================*/
/*          CMSPixel ClusteringProcessor (clustering of zs data)          */
/*          Author: Simon Spannagel                                       */
/*                (simon.spannagel@student.kit.edu or s.spannagel@cern.ch)*/
/*          Created       14 mar 2012                                     */
/*          Last modified 24 apr 2012                                     */
/*========================================================================*/

// based on EUTelAPIXClusteringProcessor:
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $ $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
 
#ifdef USE_GEAR

// eutelescope includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "CMSPixelClusteringProcessor.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelExceptions.h"
#include "EUTelHistogramManager.h"
#include "EUTelMatrixDecoder.h"
#include "EUTelTrackerDataInterfacerImpl.h"
#include "EUTelSparseClusterImpl.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/AIDAProcessor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

// lcio includes <.h>
#include <UTIL/CellIDEncoder.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/LCCollectionVec.h>

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes <.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ITree.h>
#endif

// system includes <>
#ifdef MARLINDEBUG
#include <fstream>
#include <cassert>
#endif
#include <string>
#include <sstream>
#include <vector>
#include <memory>
#include <list>
#include <algorithm>
#include <set>
#include <stdio.h>
#include <iostream>


using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string CMSPixelClusteringProcessor::_clusterSignalHistoName      	= "clusterSignal";
std::string CMSPixelClusteringProcessor::_hitMapHistoName             	= "hitMap";
std::string CMSPixelClusteringProcessor::_chargeMapHistoName            = "chargeMap";
std::string CMSPixelClusteringProcessor::_pixelPerEventHistoName  	    = "pixelPerEvent";
std::string CMSPixelClusteringProcessor::_clusterPerEventHistoName  	= "clusterPerEvent";
std::string CMSPixelClusteringProcessor::_pixelSignalHistoName	    	= "pixelSignal";
std::string CMSPixelClusteringProcessor::_clusterSizeName	        	= "clustersize";
std::string CMSPixelClusteringProcessor::_clusterSizeVsChargeName      	= "clusterSizeVsCharge";
std::string CMSPixelClusteringProcessor::_clusterSizeXName	        	= "Xclusterwidth";
std::string CMSPixelClusteringProcessor::_clusterSizeYName	        	= "Yclusterwidth";
std::string CMSPixelClusteringProcessor::_cluster1pxHistoName           = "clusters1pixel";
std::string CMSPixelClusteringProcessor::_cluster2pxHistoName           = "clusters2pixel";
std::string CMSPixelClusteringProcessor::_cluster3pxHistoName           = "clusters3pixel";
std::string CMSPixelClusteringProcessor::_cluster4pxHistoName           = "clusters4pixel";
std::string CMSPixelClusteringProcessor::_cluster1_2pxHistoName         = "clusters1+2pixel";
std::string CMSPixelClusteringProcessor::_clusterMorepxHistoName        = "clustersMorePixel";
#endif

static const int NOCLUSTER=-1;

CMSPixelClusteringProcessor::CMSPixelClusteringProcessor () : Processor("CMSPixelClusteringProcessor"), _zsDataCollectionName(""), _clusterCollectionName(""), _iRun(0), _iEvt(0), _isFirstEvent(true), _iClusters(0), _iPlaneClusters(),  _initialClusterCollectionSize(0), _minNPixels(0), _minXDistance(0), _minYDistance(0), _minDiagDistance(0), _minCharge(0), _fillHistos(false), hotPixelCollectionVec(), _hitIndexMapVec(), _noOfDetector(0), _isGeometryReady(false), _sensorIDVec(), _siPlanesParameters(), _siPlanesLayerLayout(), _orderedSensorIDVec(), _histoInfoFileName(""), _hotPixelCollectionName(""), _clusterSpectraNVector(), _clusterSpectraNxNVector(), _aidaHistoMap() {
	 _description = "CMSPixelClusteringProcessor is searching clusters in zero suppressed data.";

	registerInputCollection (LCIO::TRACKERDATA, "ZSDataCollectionName", "LCIO converted data files", _zsDataCollectionName, string("zsdata_pixel"));
	registerOutputCollection(LCIO::TRACKERPULSE, "ClusterCollectionName", "Cluster (output) collection name", _clusterCollectionName, string("cluster_pixel"));
	
	registerProcessorParameter ("MinNumberOfPixels", "Minimum number of Pixels to build a cluster", _minNPixels, 1);
	registerProcessorParameter ("MinXDistance", "Minimum distance in X that pixels should have to build a cluster", _minXDistance, 1);
	registerProcessorParameter ("MinYDistance", "Minimum distance in Y that pixels should have to build a cluster", _minYDistance, 1);
	registerProcessorParameter ("MinDiagonalDistance", "Minimum diagonal distance that pixels should have to build a cluster", _minDiagDistance, 1);
	registerProcessorParameter ("MinCharge", "Minimum Charge (TOT) that clusters should have to build a cluster", _minCharge, 1);
	
	registerProcessorParameter("HistoInfoFileName", "This is the name of the histogram information file", _histoInfoFileName, string( "histoinfo.xml" ) );
	registerProcessorParameter("HistogramFilling","Switch on or off the histogram filling", _fillHistos, static_cast< bool > ( true ) );
	
    registerOptionalParameter("HotPixelCollectionName","This is the name of the hotpixel collection",
                             _hotPixelCollectionName, static_cast< string > ( "hotpixel" ) );
                             
	_isFirstEvent = true;
	
	if (_minNPixels < 1) _minNPixels = 1;
	if (_minXDistance < 1) _minXDistance = 1;
	if (_minYDistance < 1) _minYDistance = 1;
	if (_minDiagDistance < 1) _minDiagDistance = -1;
	if (_minCharge < 1) _minCharge = 0;

}


void CMSPixelClusteringProcessor::init () {
	streamlog_out( MESSAGE4 ) << "Initing" << endl;
	streamlog_out(MESSAGE4) << "Histofile is: " << _histoInfoFileName << endl;
  
    // this method is called only once even when the rewind is active
    // usually a good idea to
    printParameters ();

	_iRun = 0;
	_iEvt = 0;
	_iClusters = 0;
	_iPlaneClusters.clear();
	
    // reset hotpixel map vectors
    _hitIndexMapVec.clear();
    
    _isGeometryReady = false;
  
}

void CMSPixelClusteringProcessor::processRunHeader (LCRunHeader * rdr) {
	streamlog_out( MESSAGE4 ) << "Processing Run Header" << endl;
	auto_ptr<EUTelRunHeaderImpl> runHeader( new EUTelRunHeaderImpl( rdr ) );
	runHeader->addProcessor( type() );
	++_iRun;
}

void CMSPixelClusteringProcessor::initializeGeometry( LCEvent * event ) throw ( marlin::SkipEventException ) {
	streamlog_out( MESSAGE4 ) << "Initializing geometry" << endl;

	_noOfDetector = 0;
	_sensorIDVec.clear();

	try {
		LCCollectionVec * collection = dynamic_cast< LCCollectionVec * > ( event->getCollection( _zsDataCollectionName ) ) ;
		_noOfDetector += collection->getNumberOfElements();

		CellIDDecoder<TrackerDataImpl > cellDecoder( collection );
		for ( size_t i = 0; i < collection->size(); ++i ) {
			TrackerDataImpl * data = dynamic_cast< TrackerDataImpl * > ( collection->getElementAt( i ) ) ;
			_sensorIDVec.push_back( cellDecoder( data )[ "sensorID" ] );
    	}

	} catch ( lcio::DataNotAvailableException ) {
    	  cout << "Cannot find collection" << endl;
	}

	_siPlanesParameters  = const_cast< gear::SiPlanesParameters*  > ( &(Global::GEAR->getSiPlanesParameters()));
	_siPlanesLayerLayout = const_cast< gear::SiPlanesLayerLayout* > ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

	// now let's build a map relating the position in the layerindex with the sensorID.
	_layerIndexMap.clear();
	for ( int iLayer = 0; iLayer < _siPlanesLayerLayout->getNLayers(); ++iLayer ) {
		_layerIndexMap.insert( make_pair( _siPlanesLayerLayout->getID( iLayer ), iLayer ) );
		_iPlaneClusters.push_back(0);
	}

	if ( _noOfDetector == 0 ) {
		streamlog_out( WARNING2 ) << "Unable to initialize the geometry. Trying with the following event" << endl;
		_iPlaneClusters.clear();
		_isGeometryReady = false;
		throw SkipEventException( this );
	} else {
		_isGeometryReady = true;
	}
	streamlog_out(MESSAGE4) << "Active SensorPlanes: ";
		for (unsigned int i = 0; i < _sensorIDVec.size(); i++) {
			streamlog_out(MESSAGE4) << " " << _sensorIDVec.at(i);
		}
		streamlog_out(MESSAGE4) << endl;
	
}


void CMSPixelClusteringProcessor::initializeHotPixelMapVec(  )
{
    streamlog_out( MESSAGE5 ) << "CMSPixelClusteringProcessor::initializeHotPixelMapVec, hotPixelCollectionVec size = " << hotPixelCollectionVec->size() << endl;
    
    // prepare some decoders
    CellIDDecoder<TrackerDataImpl> cellDecoder( hotPixelCollectionVec );

    for ( unsigned int iDetector = 0 ; iDetector < hotPixelCollectionVec->size(); iDetector++ )         
    {
        if( _hitIndexMapVec.size() < iDetector+1 )            _hitIndexMapVec.resize(iDetector+1);

        TrackerDataImpl * hotData = dynamic_cast< TrackerDataImpl * > ( hotPixelCollectionVec->getElementAt( iDetector ) );
        unsigned int sensorID            = static_cast< unsigned int > ( cellDecoder( hotData )["sensorID"] );
 
 
        if ( _layerIndexMap.find( sensorID ) == _layerIndexMap.end()   )
        {
	  streamlog_out( WARNING5 ) << "sensor " << sensorID << " not found in the present Data, so skip the hotPixel info for this sensor " << endl;
	  continue;
        }
       
        // prepare the matrix decoder
        EUTelMatrixDecoder matrixDecoder( _siPlanesLayerLayout , iDetector );

        // now prepare the EUTelescope interface to sparsified data.  
        auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> > sparseData(new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>( hotData ));

        streamlog_out ( MESSAGE1 ) << "Processing hotpixel data on detector " << sensorID << " with "
                                 << sparseData->size() << " pixels " << endl;
        
        for ( unsigned int iPixel = 0; iPixel < sparseData->size(); iPixel++ ) 
        {
            // loop over all pixels in the sparseData object.      
            EUTelGenericSparsePixel *sparsePixel =  new EUTelGenericSparsePixel() ;

            sparseData->getSparsePixelAt( iPixel, sparsePixel );
            int decoded_XY_index = matrixDecoder.getIndexFromXY( sparsePixel->getXCoord(), sparsePixel->getYCoord() ); // unique pixel index !!

            streamlog_out ( DEBUG5 )   << " iPixel " << iPixel << " idet " << iDetector << " decoded_XY_index " << decoded_XY_index << endl;
            
            if( _hitIndexMapVec[iDetector].find( decoded_XY_index ) == _hitIndexMapVec[iDetector].end() )
            {
                _hitIndexMapVec[iDetector].insert ( make_pair ( decoded_XY_index, EUTELESCOPE::FIRINGPIXEL ) );               
                streamlog_out ( DEBUG5 ) 
                    << "adding hot pixel " << " Det." << iDetector << " [" << sparsePixel->getXCoord() << " "<< sparsePixel->getYCoord() << "]" << " status : " << EUTELESCOPE::FIRINGPIXEL << endl;
              }
              else
              {
                  streamlog_out ( ERROR5 ) << "hot pixel [index " << decoded_XY_index << "] reoccered ?!" << endl;
              }
           
        }
  }
    
}


void CMSPixelClusteringProcessor::modifyEvent( LCEvent * /* event */ ){
  return;
}

void CMSPixelClusteringProcessor::processEvent (LCEvent * event) {

	++_iEvt;

	
	try {
		event->getCollection( _zsDataCollectionName ) ;
	} catch (lcio::DataNotAvailableException& e ) {
		streamlog_out ( DEBUG4 ) << "No ZS data found in the event" << endl;
	}
	
	 if ( !_isGeometryReady ) {
		initializeGeometry( event ) ;
		
		hotPixelCollectionVec = 0;
        try 
        {
            hotPixelCollectionVec = static_cast< LCCollectionVec* >  (event->getCollection( _hotPixelCollectionName )) ;
            initializeHotPixelMapVec();
            streamlog_out ( MESSAGE5 ) << "hotPixelCollectionName: " << _hotPixelCollectionName.c_str() << " found " << endl;
        } 
        catch (lcio::DataNotAvailableException& e ) 
        {
            streamlog_out ( MESSAGE5 ) << "No hot pixel DB collection (" << _hotPixelCollectionName << ") found in the event" << endl;
        }
	}
	
	
	#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
	// Book the histograms:
	if ( _fillHistos && isFirstEvent()) bookHistos();
	#endif
  
  
	EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);
	if ( evt->getEventType() == kEORE ) {
		streamlog_out ( DEBUG4 ) <<  "EORE found: nothing else to do." <<  endl;
		return;
	} else if ( evt->getEventType() == kUNKNOWN ) {
		streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " is of unknown type. Continue considering it as a normal Data Event." << endl;
	}
	
	LCCollectionVec * clusterCollection;
	bool clusterCollectionExists = false;
	_initialClusterCollectionSize = 0;
	try {
		clusterCollection = dynamic_cast< LCCollectionVec * > ( evt->getCollection( _clusterCollectionName ) );
		clusterCollectionExists = true;
		_initialClusterCollectionSize = clusterCollection->size();
	} catch ( lcio::DataNotAvailableException& e ) {
		clusterCollection = new LCCollectionVec(LCIO::TRACKERPULSE);
	}

	
    // Start the clustering...	
	Clustering(evt, clusterCollection);
	
	
	// If we found some clusters (event not empty), we add the collection
	 if ( ! clusterCollectionExists && ( clusterCollection->size() != _initialClusterCollectionSize )) {
		evt->addCollection( clusterCollection, _clusterCollectionName );
	}
	
	
	#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
	if ( clusterCollection->size() != _initialClusterCollectionSize ) {
		if ( _fillHistos ) fillHistos(event);
	}
	#endif
  
  
	if ( ! clusterCollectionExists && ( clusterCollection->size() == _initialClusterCollectionSize ) ) {
		delete clusterCollection;
    	streamlog_out ( DEBUG5 ) << "delete clusterCollection;" <<  endl;
	}


	_isFirstEvent = false;
	streamlog_out ( DEBUG5 ) << "End of process event" <<  endl;
}


void CMSPixelClusteringProcessor::Clustering(LCEvent * evt, LCCollectionVec * clusterCollection) {
	
	//Data Input
	LCCollectionVec * zsInputCollectionVec  = dynamic_cast < LCCollectionVec * > (evt->getCollection( _zsDataCollectionName ));
	CellIDDecoder<TrackerDataImpl> cellDecoder( zsInputCollectionVec );
	
	//Data Output
    CellIDEncoder< TrackerPulseImpl > zsDataEncoder ( "sensorID:5,clusterID:12,xSeed:9,ySeed:10,xCluSize:9,yCluSize:9,type:5", clusterCollection );	

    LCCollectionVec * sparseClusterCollectionVec = NULL;
    sparseClusterCollectionVec =  new LCCollectionVec(LCIO::TRACKERDATA);
	  
    CellIDEncoder<TrackerDataImpl> idClusterEncoder( "sensorID:5,clusterID:12,sparsePixelType:5,type:6", sparseClusterCollectionVec  );
    
    
	
	// One loop for each sensor plane
	for ( unsigned int i = 0 ; i < zsInputCollectionVec->size(); i++ ) {
		TrackerDataImpl * zsData = dynamic_cast< TrackerDataImpl * > ( zsInputCollectionVec->getElementAt( i ) );
		SparsePixelType   type   = static_cast<SparsePixelType> ( static_cast<int> (cellDecoder( zsData )["sparsePixelType"]) );

		unsigned int sensorID              = static_cast<unsigned int > ( cellDecoder( zsData )["sensorID"] );
		streamlog_out ( DEBUG5 ) << "evt " << evt->getEventNumber() << " SensorID " << sensorID << endl;

        // prepare the matrix decoder
        EUTelMatrixDecoder matrixDecoder( _siPlanesLayerLayout , sensorID );
		
		if (type == kEUTelGenericSparsePixel  ) {
		    auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> > pixelData( new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>( zsData ));
			streamlog_out ( DEBUG5 ) << "Processing data on detector " << sensorID << ", " << pixelData->size() << " pixels " << endl;

			// Loop over all pixels in the sparseData object.
			std::vector<EUTelGenericSparsePixel*> PixelVec;
			EUTelGenericSparsePixel Pixel;

			 //Push all single Pixels of one plane in the PixelVec
			for ( unsigned int iPixel = 0; iPixel < pixelData->size(); iPixel++ ) {
				pixelData->getSparsePixelAt( iPixel, &Pixel);

                // HotPixel treatment: check if we read a hotpixel db:
				if(_hitIndexMapVec.size() > sensorID) {
                    // Now remove HotPixels
                    int index = matrixDecoder.getIndexFromXY( Pixel.getXCoord(), Pixel.getYCoord() );
                    if( _hitIndexMapVec[sensorID].find( index ) != _hitIndexMapVec[sensorID].end() )
                    {
                        streamlog_out ( DEBUG5 ) << "Detector " << sensorID << " Pixel " << Pixel.getXCoord() << " " << Pixel.getYCoord() << " -- HOTPIXEL, skipping... " << endl;
                        continue;
                    }
                }

				PixelVec.push_back(new EUTelGenericSparsePixel(Pixel));
			}
			
			streamlog_out ( DEBUG5 ) << "Hit Pixels: " << PixelVec.size() << endl;
			
			/* --- Here the real clustering happens --- */
			 int nClusters = 0;
			std::vector<int> clusterNumber;

			if (!PixelVec.empty()) clusterNumber.assign(PixelVec.size(), NOCLUSTER);
			
			
			streamlog_out( DEBUG5 ) << "Starting with clustering..." << endl;
			if (PixelVec.size() > 1) {
				//Now compare all pixel in one plane with each other
				for ( unsigned int aPixel=0;aPixel < PixelVec.size();++aPixel) {
					for ( unsigned int bPixel=aPixel+1; bPixel < PixelVec.size(); ++bPixel) {
						EUTelGenericSparsePixel *aPix = PixelVec.at(aPixel);
						EUTelGenericSparsePixel *bPix = PixelVec.at(bPixel);
						int xDist = abs( aPix->getXCoord() - bPix->getXCoord() );
						int yDist = abs( aPix->getYCoord() - bPix->getYCoord() );
						
						if ( ( xDist <= _minXDistance) && (yDist <= _minYDistance)  ) { // This means they are neighbours in x-direction && this means they are neighbours in y-direction
							bool skipPixel = false;
							if (_minDiagDistance != -1) {
								bool areDiagonalPartners =  (xDist != 0 && yDist != 0);  // If true then these pixels are diagonal partners
								int diagDist = max(xDist,yDist);
								if (areDiagonalPartners && diagDist > _minDiagDistance) skipPixel = true;
							}
						
							if (skipPixel == false) {
								
								if ( (clusterNumber.at(aPixel) == NOCLUSTER) && clusterNumber.at(bPixel) == NOCLUSTER) { // None of these pixels have been assigned to a cluster
									++nClusters;
									clusterNumber.at(aPixel) = nClusters;
									clusterNumber.at(bPixel) = nClusters;
									streamlog_out ( DEBUG5 ) << "assigning clusternumber: " << nClusters << endl;
								}else if ( (clusterNumber.at(aPixel) == NOCLUSTER) && clusterNumber.at(bPixel) != NOCLUSTER) { 
								    // b was assigned already, a not
									clusterNumber.at(aPixel) = clusterNumber.at(bPixel);
								}else if ( (clusterNumber.at(aPixel) != NOCLUSTER) && clusterNumber.at(bPixel) == NOCLUSTER) { 
								    // a was assigned already, b not
									clusterNumber.at(bPixel) = clusterNumber.at(aPixel);
								}else { 
								    // Both pixels have a clusternumber already
									int min = std::min(clusterNumber.at(aPixel), clusterNumber.at(bPixel));
									clusterNumber.at(aPixel) = min;
									clusterNumber.at(bPixel) = min;
								}
							} else { // These pixels are not neighboured due to diagonal cut
								if ( clusterNumber.at(aPixel) == NOCLUSTER ) {
									++nClusters;
									clusterNumber.at(aPixel) = nClusters;
								}
								if ( clusterNumber.at(bPixel) == NOCLUSTER ) {
									++nClusters;
									clusterNumber.at(bPixel) = nClusters;
								}
							}
						} else { // These pixels are not neighboured
							if ( clusterNumber.at(aPixel) == NOCLUSTER ) {
								++nClusters;
								clusterNumber.at(aPixel) = nClusters;
							}
							if ( clusterNumber.at(bPixel) == NOCLUSTER ) {
								++nClusters;
								clusterNumber.at(bPixel) = nClusters;
							}
						}
					}
				}
			} else { 
			    // You can't use the clustering algorithm with only one pixel
				if (PixelVec.size() == 1) clusterNumber.at(0) = 1;
			}
			
			// It happens, that nCluster is higher that the real number of clusters, as a clusternumber can be deleted in the else-block, so lets correct it
			std::set<int> clusterSet;
			for (unsigned int i=0; i<PixelVec.size();++i) {
				clusterSet.insert(clusterNumber.at(i));
				if ( clusterNumber.at(i) == NOCLUSTER) {streamlog_out(MESSAGE1) << "Cluster with ID -1 in " << i << endl;}
			}
			nClusters = clusterSet.size();
			
			
			if (nClusters != 0) streamlog_out( DEBUG5 ) << "Found " << nClusters << " clusters in sensor " << sensorID<< endl;
			
			/* --- Finished Clustering --- */
			
			/* --- Push back one Collection per Cluster --- */
			std::set<int>::iterator it;
			for (it = clusterSet.begin(); it != clusterSet.end(); ++it) {
				lcio::TrackerPulseImpl * pulseFrame = new lcio::TrackerPulseImpl();
				lcio::TrackerDataImpl * clusterFrame = new lcio::TrackerDataImpl();
				eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelGenericSparsePixel > *pixelCluster = new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelGenericSparsePixel >(clusterFrame);	
				for (unsigned int i=0;i< PixelVec.size();i++) {
					if(clusterNumber.at(i)== *it) { 
					    // Put only these pixels in that ClusterCollection that belong to that cluster
						pixelCluster->addSparsePixel(PixelVec.at(i));
						streamlog_out( DEBUG5 ) << "Adding Pixel " << i << " to cluster " << clusterNumber.at(i) << endl;
					}
				}
	            
	            streamlog_out( DEBUG5 ) << "size: " << pixelCluster->size() << ">=" << _minNPixels << " && charge: " << pixelCluster->getTotalCharge() << ">=" << _minCharge << endl;
                if ( (pixelCluster->size() >= static_cast< unsigned int >(_minNPixels)) && (pixelCluster->getTotalCharge() >= static_cast< unsigned int >(_minCharge)) ) {

					float x,y;
					int xsize,ysize;
					
                    pulseFrame->setCharge(pixelCluster->getTotalCharge());

					// This is using analogue hit information:
					pixelCluster->getCenterOfGravity(x,y);
					pixelCluster->getClusterSize(xsize,ysize);
					if (x >= 0 && x <= _siPlanesLayerLayout->getSensitiveNpixelX( _layerIndexMap[ sensorID ] ) && 
					    y >= 0 && y <= _siPlanesLayerLayout->getSensitiveNpixelY( _layerIndexMap[ sensorID ] )) {
						streamlog_out( DEBUG5 ) << "Clustervars: ROC" << sensorID << " Cl" << *it << " x" << x << " y" << y << " dx" <<xsize << " dy" << ysize << " " << type << endl;
						_iClusters++;
						_iPlaneClusters[sensorID]++;

						int clusterID = *it;
						zsDataEncoder["sensorID"]      = sensorID;
						zsDataEncoder["clusterID"]     = clusterID;
						zsDataEncoder["xSeed"]         = static_cast< long >(x);
						zsDataEncoder["ySeed"]         = static_cast< long >(y);
						zsDataEncoder["xCluSize"]      = xsize;
						zsDataEncoder["yCluSize"]      = ysize;
						zsDataEncoder["type"]          = static_cast<int>(kEUTelSparseClusterImpl);
						zsDataEncoder.setCellID(pulseFrame);
						pulseFrame->setTrackerData(clusterFrame);
						clusterCollection->push_back(pulseFrame);
						
						idClusterEncoder["sensorID"] 		= sensorID;
						idClusterEncoder["clusterID"] 		= clusterID;
						idClusterEncoder["sparsePixelType"] 	= static_cast<int>(type);
						idClusterEncoder["type"] 		= static_cast<int>(kEUTelSparseClusterImpl);
						idClusterEncoder.setCellID(clusterFrame);
						sparseClusterCollectionVec->push_back(clusterFrame);
					}
					else streamlog_out( DEBUG5 ) << "No cluster: ROC" << sensorID << " Cl" << *it << " x" << x << " y" << y << " dx" <<xsize << " dy" << ysize << " " << type << endl;

				}
            }
		 }
	}
    evt->addCollection( sparseClusterCollectionVec, "original_zsdata" );
}

void CMSPixelClusteringProcessor::check (LCEvent * /* evt */) {
    // Nothing to check here - could be used to fill check plots in reconstruction processor
}


void CMSPixelClusteringProcessor::end() {
	streamlog_out( MESSAGE4 ) << "Ending" << endl;
	for(unsigned int i = 0; i < _iPlaneClusters.size(); i++)
	    	streamlog_out( MESSAGE4 ) << "Found " << _iPlaneClusters[i] << " clusters on detector " << i << endl;
	streamlog_out( MESSAGE4 ) << "Clusters found in total: " << _iClusters << endl;
}


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
void CMSPixelClusteringProcessor::fillHistos (LCEvent * evt) {

	EUTelEventImpl * eutelEvent = static_cast<EUTelEventImpl*> (evt);
	EventType type              = eutelEvent->getEventType();

	if ( type == kEORE ) {
		streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
		return ;
	} else if ( type == kUNKNOWN ) {
	}

	  try {

		LCCollectionVec * clusterCollectionVec = dynamic_cast<LCCollectionVec*>  (evt->getCollection(_clusterCollectionName));
		CellIDDecoder<TrackerPulseImpl> cellDecoder(clusterCollectionVec);
		

		vector<unsigned short> eventCounterVec( _noOfDetector, 0 );
    	vector<unsigned short> noOfClusters( _noOfDetector, 0);
		string tempHistoName;

		for ( int iCluster = _initialClusterCollectionSize; iCluster < clusterCollectionVec->getNumberOfElements(); iCluster++ ) {
			TrackerPulseImpl * pulseFrame = dynamic_cast<TrackerPulseImpl*> ( clusterCollectionVec->getElementAt(iCluster) );
			TrackerDataImpl * clusterFrame = dynamic_cast<TrackerDataImpl*> ( pulseFrame->getTrackerData() );
			
			int sensorID = ( static_cast<int> ( cellDecoder(pulseFrame)["sensorID"] ));
			// FIXME: check if correctly working with "ClusterType":
			// SparsePixelType type = static_cast<SparsePixelType> ( static_cast<int> (cellDecoder( pulseFrame )["type"]) );
			ClusterType type = static_cast<ClusterType> ( static_cast<int> (cellDecoder( pulseFrame )["type"]) );
			eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelGenericSparsePixel > *pixelCluster = new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelGenericSparsePixel >(clusterFrame);	
			
			if (type == kEUTelSparseClusterImpl  ) {
			
				int size = pixelCluster->size();
				eventCounterVec[ sensorID ] += size;
				noOfClusters[ sensorID ]++;
				
				int xSize, ySize;
				pixelCluster->getClusterSize(xSize,ySize);


                if(size <= 2) {
				    tempHistoName = _cluster1_2pxHistoName + "_d" + to_string( sensorID ) ;
				    (dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(pixelCluster->getTotalCharge());
                }                
                else if(size > 2) {
				    tempHistoName = _clusterMorepxHistoName + "_d" + to_string( sensorID ) ;
				    (dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(pixelCluster->getTotalCharge());                
                }
                
                
                if(size == 1) {
				    tempHistoName = _cluster1pxHistoName + "_d" + to_string( sensorID ) ;
				    (dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(pixelCluster->getTotalCharge());
                }
                else if(size == 2) {
				    tempHistoName = _cluster2pxHistoName + "_d" + to_string( sensorID ) ;
				    (dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(pixelCluster->getTotalCharge());                
                }
                else if(size == 3) {
				    tempHistoName = _cluster3pxHistoName + "_d" + to_string( sensorID ) ;
				    (dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(pixelCluster->getTotalCharge());              
                }
                else if(size == 4) {
				    tempHistoName = _cluster4pxHistoName + "_d" + to_string( sensorID ) ;
				    (dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(pixelCluster->getTotalCharge());               
                }
                
				tempHistoName = _clusterSignalHistoName + "_d" + to_string( sensorID ) ;
				(dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(pixelCluster->getTotalCharge());
				
				tempHistoName = _clusterSizeName + "_d" + to_string( sensorID );
				(dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(size);
				
				tempHistoName = _clusterSizeXName + "_d" + to_string( sensorID );
				(dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(xSize);
				
				tempHistoName = _clusterSizeYName + "_d" + to_string( sensorID );
				(dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(ySize);
				
				tempHistoName = _clusterSizeVsChargeName + "_d" + to_string( sensorID );
				(dynamic_cast<AIDA::IHistogram2D*> (_aidaHistoMap[tempHistoName]))->fill(size, pixelCluster->getTotalCharge(), 1.);
				
				
				for (int iPixel=0; iPixel < size; iPixel++) {
					EUTelGenericSparsePixel Pixel;
					pixelCluster->getSparsePixelAt(iPixel, &Pixel);
					tempHistoName = _pixelSignalHistoName + "_d" + to_string( sensorID );
					(dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(Pixel.getSignal());
					
					tempHistoName = _chargeMapHistoName + "_d" + to_string( sensorID );
    				(dynamic_cast<AIDA::IHistogram2D*> (_aidaHistoMap[tempHistoName]))->fill(Pixel.getXCoord(), Pixel.getYCoord(), Pixel.getSignal());

				}
				
				tempHistoName = _hitMapHistoName + "_d" + to_string( sensorID );
				int xSeed, ySeed;
				pixelCluster->getCenterCoord(xSeed, ySeed);
				(dynamic_cast<AIDA::IHistogram2D*> (_aidaHistoMap[tempHistoName]))->fill(static_cast<double >(xSeed), static_cast<double >(ySeed), 1.);
								               
			}
			
			
			
		}
		
        for(unsigned int i = 0; i < _noOfDetector; i++) {
        
		    tempHistoName = _pixelPerEventHistoName + "_d" + to_string( i );
		    (dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(eventCounterVec[i]);        

            tempHistoName = _clusterPerEventHistoName + "_d" + to_string( i );
            (dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(noOfClusters[i]);

		}


	} catch (lcio::DataNotAvailableException& e) {
		return;
	}

}

#endif

#ifdef MARLIN_USE_AIDA
void CMSPixelClusteringProcessor::bookHistos() {
	
	streamlog_out ( MESSAGE0 )  << "Booking histograms for " << _sensorIDVec.size() << " detectors..." << endl;

	string tempHistoName;
	string basePath;
	for (size_t iDetector = 0; iDetector < _sensorIDVec.size(); iDetector++) {

		int sensorID = _sensorIDVec.at( iDetector );

		// The min and max information are taken from GEAR.
		int minX, minY, maxX, maxY;
		minX = 0;
		minY = 0;

		if ( _layerIndexMap.find( sensorID ) != _layerIndexMap.end() ){
			maxX = _siPlanesLayerLayout->getSensitiveNpixelX( _layerIndexMap[ sensorID ] ) - 1;
			maxY = _siPlanesLayerLayout->getSensitiveNpixelY( _layerIndexMap[ sensorID ] ) - 1;
		} else {
			throw  InvalidGeometryException ("Wrong geometry file? Unknown sensorID " + to_string( sensorID ));
		}

		basePath = "detector_" + to_string( sensorID );
		AIDAProcessor::tree(this)->mkdir(basePath.c_str());
		basePath.append("/");

		int    signalNBin  = 1999;
		double signalMin   = 0.;
		double signalMax   = 2000.;
		
		int clusterNBin = 10;
		double clusterMin = 0.;
		double clusterMax = 10.;
		
		string clusterSizeTitle = "Clustersize";
		tempHistoName = _clusterSizeName + "_d" + to_string( sensorID );
		AIDA::IHistogram1D * clusterSizeHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), clusterNBin,clusterMin,clusterMax);
		_aidaHistoMap.insert(make_pair(tempHistoName, clusterSizeHisto));
		clusterSizeHisto->setTitle(clusterSizeTitle.c_str());
		
		string clusterSizeXTitle = "Clusterwidth in X";
		tempHistoName = _clusterSizeXName + "_d" + to_string( sensorID );
		AIDA::IHistogram1D * clusterXSize = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), clusterNBin,clusterMin,clusterMax);
		_aidaHistoMap.insert(make_pair(tempHistoName, clusterXSize));
		clusterXSize->setTitle(clusterSizeXTitle.c_str());
		
		string clusterSizeYTitle = "Clusterwidth in Y";
		tempHistoName = _clusterSizeYName + "_d" + to_string( sensorID );
		AIDA::IHistogram1D * clusterYSize = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), clusterNBin,clusterMin,clusterMax);
		_aidaHistoMap.insert(make_pair(tempHistoName, clusterYSize));
		clusterYSize->setTitle(clusterSizeYTitle.c_str());

		string clusterEventTitle = "Clusters per event";
		tempHistoName = _clusterPerEventHistoName + "_d" + to_string( sensorID );
		AIDA::IHistogram1D * clusterEventHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), clusterNBin,clusterMin,clusterMax);
		_aidaHistoMap.insert(make_pair(tempHistoName, clusterEventHisto));
		clusterEventHisto->setTitle(clusterEventTitle.c_str());
		
		string pixelEventTitle = "Pixels per event";
		tempHistoName = _pixelPerEventHistoName + "_d" + to_string( sensorID );
		AIDA::IHistogram1D * pixelEventHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), clusterNBin,clusterMin,clusterMax);
		_aidaHistoMap.insert(make_pair(tempHistoName, pixelEventHisto));
		pixelEventHisto->setTitle(pixelEventTitle.c_str());

		string clusterTitle = "Cluster spectrum with all pixels";
		tempHistoName = _clusterSignalHistoName + "_d" + to_string( sensorID );
		AIDA::IHistogram1D * clusterSignalHisto =  AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), signalNBin,signalMin,signalMax);
		_aidaHistoMap.insert(make_pair(tempHistoName, clusterSignalHisto));
		clusterSignalHisto->setTitle(clusterTitle.c_str());


		clusterTitle = "1 pixel clusters";
		tempHistoName = _cluster1pxHistoName + "_d" + to_string( sensorID );
		AIDA::IHistogram1D * cluster1pxHisto =  AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), signalNBin,signalMin,signalMax);
		_aidaHistoMap.insert(make_pair(tempHistoName, cluster1pxHisto));
		cluster1pxHisto->setTitle(clusterTitle.c_str());

		clusterTitle = "2 pixel clusters";
		tempHistoName = _cluster2pxHistoName + "_d" + to_string( sensorID );
		AIDA::IHistogram1D * cluster2pxHisto =  AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), signalNBin,signalMin,signalMax);
		_aidaHistoMap.insert(make_pair(tempHistoName, cluster2pxHisto));
		cluster2pxHisto->setTitle(clusterTitle.c_str());

		clusterTitle = "3 pixel clusters";
		tempHistoName = _cluster3pxHistoName + "_d" + to_string( sensorID );
		AIDA::IHistogram1D * cluster3pxHisto =  AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), signalNBin,signalMin,signalMax);
		_aidaHistoMap.insert(make_pair(tempHistoName, cluster3pxHisto));
		cluster3pxHisto->setTitle(clusterTitle.c_str());

		clusterTitle = "4 pixel clusters";
		tempHistoName = _cluster4pxHistoName + "_d" + to_string( sensorID );
		AIDA::IHistogram1D * cluster4pxHisto =  AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), signalNBin,signalMin,signalMax);
		_aidaHistoMap.insert(make_pair(tempHistoName, cluster4pxHisto));
		cluster4pxHisto->setTitle(clusterTitle.c_str());

		clusterTitle = "1 and 2 pixel clusters";
		tempHistoName = _cluster1_2pxHistoName + "_d" + to_string( sensorID );
		AIDA::IHistogram1D * cluster1_2pxHisto =  AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), signalNBin,signalMin,signalMax);
		_aidaHistoMap.insert(make_pair(tempHistoName, cluster1_2pxHisto));
		cluster1_2pxHisto->setTitle(clusterTitle.c_str());

		clusterTitle = "Clusters with more than 2 pixels";
		tempHistoName = _clusterMorepxHistoName + "_d" + to_string( sensorID );
		AIDA::IHistogram1D * clusterMorepxHisto =  AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), signalNBin,signalMin,signalMax);
		_aidaHistoMap.insert(make_pair(tempHistoName, clusterMorepxHisto));
		clusterMorepxHisto->setTitle(clusterTitle.c_str());
	

		string singlePixelTitle = "Spectrum of single pixels";
		tempHistoName = _pixelSignalHistoName + "_d" + to_string( sensorID );
		AIDA::IHistogram1D * pixelSignalHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), signalNBin,signalMin,signalMax);
		_aidaHistoMap.insert(make_pair(tempHistoName, pixelSignalHisto));
		pixelSignalHisto->setTitle(singlePixelTitle.c_str());
	
	
    	string clusterSizeVsChargeTitle = "Cluster Size vs. Cluster Charge";
		tempHistoName = _clusterSizeVsChargeName + "_d" + to_string( sensorID );
		AIDA::IHistogram2D * clusterSizeVsChargeHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(), clusterNBin, clusterMin, clusterMax,signalNBin, signalMin, signalMax);
		_aidaHistoMap.insert(make_pair(tempHistoName, clusterSizeVsChargeHisto));
		clusterSizeVsChargeHisto->setTitle(clusterSizeVsChargeTitle.c_str());
	
	
		tempHistoName = _hitMapHistoName + "_d" + to_string( sensorID );
		int     xBin = maxX - minX + 1;
		double  xMin = static_cast<double >( minX ) - 0.5 ;
		double  xMax = static_cast<double >( maxX ) + 0.5;
		int     yBin = maxY - minY + 1;
		double  yMin = static_cast<double >( minY ) - 0.5;
		double  yMax = static_cast<double >( maxY ) + 0.5;
		AIDA::IHistogram2D * hitMapHisto =
		AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(), xBin, xMin, xMax,yBin, yMin, yMax);
		_aidaHistoMap.insert(make_pair(tempHistoName, hitMapHisto));
		hitMapHisto->setTitle("Hit map");
		
	
		tempHistoName = _chargeMapHistoName + "_d" + to_string( sensorID );
		AIDA::IHistogram2D * chargeMapHisto =
		AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(), xBin, xMin, xMax,yBin, yMin, yMax);
		_aidaHistoMap.insert(make_pair(tempHistoName, chargeMapHisto));
		hitMapHisto->setTitle("Charge map [in DAC]");
		
	}
	streamlog_out ( MESSAGE0 )  << "end of Booking histograms " << endl;
}
#endif
#endif // USE_GEAR
