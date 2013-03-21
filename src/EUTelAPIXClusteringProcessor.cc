// Version: $Id$
   
#ifdef USE_GEAR

// eutelescope includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelAPIXClusteringProcessor.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelExceptions.h"
#include "EUTelHistogramManager.h"
#include "EUTelMatrixDecoder.h"
#include "EUTelSparseDataImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelSparseData2Impl.h"
#include "EUTelSparseCluster2Impl.h"


// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/AIDAProcessor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h" // this is because we want to use GEAR.

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



using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;

#ifdef MARLINDEBUG
/// /* DEBUG */ ofstream logfile;
#endif

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelAPIXClusteringProcessor::_clusterSignalHistoName      	= "clusterSignal";
std::string EUTelAPIXClusteringProcessor::_seedSignalHistoName         	= "seedSignal";
std::string EUTelAPIXClusteringProcessor::_hitMapHistoName             	= "hitMap";
std::string EUTelAPIXClusteringProcessor::_pixelMapHistoName           	= "pixelMap";
std::string EUTelAPIXClusteringProcessor::_seedSNRHistoName            	= "seedSNR";
std::string EUTelAPIXClusteringProcessor::_clusterNoiseHistoName       	= "clusterNoise";
std::string EUTelAPIXClusteringProcessor::_clusterSNRHistoName         	= "clusterSNR";
std::string EUTelAPIXClusteringProcessor::_eventMultiplicityHistoName  	= "eventMultiplicity";
std::string EUTelAPIXClusteringProcessor::_pixelSignalHistoName		= "pixelSignal";
std::string EUTelAPIXClusteringProcessor::_clusterSizeName		= "clustersize";
std::string EUTelAPIXClusteringProcessor::_clusterSizeXName		= "Xclusterwidth";
std::string EUTelAPIXClusteringProcessor::_clusterSizeYName		= "Yclusterwidth";
std::string EUTelAPIXClusteringProcessor::_lvl1TriggerName		= "lvl1trigger";
std::string EUTelAPIXClusteringProcessor::_lvl1TriggerDiffName		= "lvl1triggerdiff";
#endif

static const int NOCLUSTER=-1;

EUTelAPIXClusteringProcessor::EUTelAPIXClusteringProcessor () : Processor("EUTelAPIXClusteringProcessor") {

        _description = "EUTelAPIXClusteringProcessor is looking for clusters.";

	registerInputCollection (LCIO::TRACKERDATA, "ZSDataCollectionName", "LCIO converted data files", _zsDataCollectionName, string("zsdata_apix"));
	registerOutputCollection(LCIO::TRACKERPULSE, "ClusterCollectionName", "Cluster (output) collection name", _clusterCollectionName, string("cluster_apix"));
	//registerOutputCollection(LCIO::TRACKERPULSE, "ClusteredZSDataCollectionName", "LCIO converted data file with cluster-IDs", _clusterCollectionName, string("clustered_zsdata_apix"));
	
	registerProcessorParameter ("MinNumberOfPixels", "Minimum number of Pixels to build a cluster", _minNPixels, 1);
	registerProcessorParameter ("MinXDistance", "Minimum distance in X that pixels should have to build a cluster", _minXDistance, 1);
	registerProcessorParameter ("MinYDistance", "Minimum distance in Y that pixels should have to build a cluster", _minYDistance, 1);
	registerProcessorParameter ("MinDiagonalDistance", "Minimum diagonal distance that pixels should have to build a cluster", _minDiagDistance, 1);
	registerProcessorParameter ("MinCharge", "Minimum Charge (TOT) that clusters should have to build a cluster", _minCharge, 1);
	registerProcessorParameter ("MinLVL1Difference", "Minimim Level 1 difference of pixels building one cluster", _minLVL1, -1);
        registerProcessorParameter ("MaxXSize","Maximim size along X direction", _maxXsize, 1000 );
        registerProcessorParameter ("MaxYSize","Maximim size along Y direction", _maxYsize, 1000 );
        registerProcessorParameter ("MinXSize","Minimim size along X direction", _minXsize,  0   );
        registerProcessorParameter ("MinYSize","Minimim size along Y direction", _minYsize,  0   );

        registerProcessorParameter ("ClusterLimits","Minimim/Maximum size along X and Y direction [ID xmin xmax ymin ymax]", _ClusterLimits,  std::vector<int> (static_cast<int> (5), -1) );
	
	registerProcessorParameter("HistoInfoFileName", "This is the name of the histogram information file", _histoInfoFileName, string( "histoinfo.xml" ) );
	registerProcessorParameter("HistogramFilling","Switch on or off the histogram filling", _fillHistos, static_cast< bool > ( true ) );
	_isFirstEvent = true;
	
	if (_minNPixels < 1) _minNPixels = 1;
	if (_minXDistance < 1) _minXDistance = 1;
	if (_minYDistance < 1) _minYDistance = 1;
	if (_minDiagDistance < 1) _minDiagDistance = -1;
	if (_minCharge < 1) _minCharge = 0;
	if (_minLVL1 > 15) _minLVL1 = -1;


}


void EUTelAPIXClusteringProcessor::init () {
	streamlog_out( MESSAGE4 ) << "Initing" << endl;
	streamlog_out(MESSAGE4) << "Histofile is: " << _histoInfoFileName << endl;
	_iRun = 0;
	_iEvt = 0;
  // this method is called only once even when the rewind is active
  // usually a good idea to

        _isGeometryReady = false;

    printParameters ();

    if( _ClusterLimits.size() % 5 != 0) 
    {
      streamlog_out ( ERROR5 ) << "Wrong number of input values in ClusterLimits; check your steering file !" <<endl;
    }
    else
    {
      streamlog_out ( MESSAGE5 ) << " ClusterLimits initialised with " << _ClusterLimits.size() << " elements" << endl;
    } 

    for(int i = 0; i < static_cast< int >(_ClusterLimits.size()); i=i+5)
    {
       int iSensor = _ClusterLimits.at(i);
       _ClusterLimitsMap[ iSensor ].push_back(_ClusterLimits.at(i+1));
       _ClusterLimitsMap[ iSensor ].push_back(_ClusterLimits.at(i+2));
       _ClusterLimitsMap[ iSensor ].push_back(_ClusterLimits.at(i+3));
       _ClusterLimitsMap[ iSensor ].push_back(_ClusterLimits.at(i+4));
    }  

}

void EUTelAPIXClusteringProcessor::processRunHeader (LCRunHeader * rdr) {
	streamlog_out( MESSAGE4 ) << "Processing Run Header" << endl;
	auto_ptr<EUTelRunHeaderImpl> runHeader( new EUTelRunHeaderImpl( rdr ) );
	runHeader->addProcessor( type() );
	++_iRun;
}

void EUTelAPIXClusteringProcessor::initializeGeometry( LCEvent * event ) throw ( marlin::SkipEventException ) {
	streamlog_out( MESSAGE4 ) << "Initializing geometry" << endl;

	_noOfDetector = 0;
	_sensorIDVec.clear();

	try {
		LCCollectionVec * collection = dynamic_cast< LCCollectionVec * > ( event->getCollection( _zsDataCollectionName ) ) ;
		_noOfDetector += collection->getNumberOfElements();
		//cout<<"collection->size()= "<<collection->size()<<endl;
		CellIDDecoder<TrackerDataImpl > cellDecoder( collection );
		for ( size_t i = 0; i < collection->size(); ++i ) {
			TrackerDataImpl * data = dynamic_cast< TrackerDataImpl * > ( collection->getElementAt( i ) ) ;
			_sensorIDVec.push_back( cellDecoder( data )[ "sensorID" ] );
			_totClusterMap.insert( make_pair( cellDecoder( data )[ "sensorID" ] , 0 ));
	}

	} catch ( lcio::DataNotAvailableException ) {
		// do nothing again
	  cout << "Cannot find collection" << endl;
	}

	_siPlanesParameters  = const_cast< gear::SiPlanesParameters*  > ( &(Global::GEAR->getSiPlanesParameters()));
	_siPlanesLayerLayout = const_cast< gear::SiPlanesLayerLayout* > ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

	// now let's build a map relating the position in the layerindex
	// with the sensorID.
	_layerIndexMap.clear();
	for ( int iLayer = 0; iLayer < _siPlanesLayerLayout->getNLayers(); ++iLayer ) {
		_layerIndexMap.insert( make_pair( _siPlanesLayerLayout->getID( iLayer ), iLayer ) );
	}

	// now another map relating the position in the ancillary
	_orderedSensorIDVec.clear();

	if ( _noOfDetector == 0 ) {
		streamlog_out( WARNING2 ) << "Unable to initialize the geometry. Trying with the following event" << endl;
		_isGeometryReady = false;
		throw SkipEventException( this );
	} else {
		_isGeometryReady = true;
	}
	streamlog_out(MESSAGE4) << "Active SensorPlanes: ";
		for (size_t i= 0; i< _sensorIDVec.size(); i++) {
			streamlog_out(MESSAGE4) << " " << _sensorIDVec.at(i);
		}
		streamlog_out(MESSAGE4) << endl;
	
}

void EUTelAPIXClusteringProcessor::modifyEvent( LCEvent * /* event */ ){
  return;
}

void EUTelAPIXClusteringProcessor::processEvent (LCEvent * event) {
	if (_iEvt % 50 == 0)
	streamlog_out( MESSAGE4 ) << "Processing event "
		<< setw(6) << setiosflags(ios::right) << event->getEventNumber() << " in run "
		<< setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber() << setfill(' ')
		<< " (Total = " << setw(10) << _iEvt << ")" << resetiosflags(ios::left) << endl;
	++_iEvt;
	
	/* -------------- */

        if ( !_isGeometryReady ) {
		initializeGeometry( event ) ;
	}
	#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
	// book the histograms now
	if ( _fillHistos && isFirstEvent() ) {
		bookHistos();
	}
	#endif
  
	/* -------------- */

	try {
		event->getCollection( _zsDataCollectionName ) ;
	} catch (lcio::DataNotAvailableException& e ) {
		streamlog_out ( DEBUG4 ) << "No ZS data found in the event " << _iEvt << endl;
          return;
	}
	
	/* -------------- */
	
	
	EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);
	if ( evt->getEventType() == kEORE ) {
		streamlog_out ( DEBUG4 ) <<  "EORE found: nothing else to do." <<  endl;
		return;
	} else if ( evt->getEventType() == kUNKNOWN ) {
		streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " is of unknown type. Continue considering it as a normal Data Event." << endl;
	}
	
	/* -------------- */
	
	LCCollectionVec * clusterCollection;
	bool clusterCollectionExists = false;
	_initialClusterCollectionSize = 0;
	try {
		clusterCollection = dynamic_cast< LCCollectionVec * > ( evt->getCollection( _clusterCollectionName ) );
		clusterCollectionExists = true;
		_initialClusterCollectionSize = clusterCollection->size();
	} catch ( lcio::DataNotAvailableException& e ) {
		clusterCollection = new LCCollectionVec(LCIO::TRACKERPULSE);
		//streamlog_out ( MESSAGE1 ) << "Creating new Collection " << _pulseCollectionName << endl;
	}
	
	/* -------------- */
	
	Clustering(evt, clusterCollection);
	
	/* -------------- */
	
        if ( ! clusterCollectionExists && ( clusterCollection->size() != _initialClusterCollectionSize ))  
        {
		evt->addCollection( clusterCollection, _clusterCollectionName );
		//streamlog_out ( MESSAGE1 ) << "Adding collection " << _clusterCollectionName <<  endl;
	}
	#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
	if ( clusterCollection->size() != _initialClusterCollectionSize ) {
		if ( _fillHistos ) fillHistos(event);
	}
	#endif
  
	if ( ! clusterCollectionExists && ( clusterCollection->size() == _initialClusterCollectionSize ) ) {
		delete clusterCollection;
	}


	_isFirstEvent = false;
	//streamlog_out ( MESSAGE1 ) << "End of process event" <<  endl;
}

void EUTelAPIXClusteringProcessor::Clustering(LCEvent * evt, LCCollectionVec * clusterCollection) {
	
	//Data Input
	LCCollectionVec * zsInputCollectionVec  = dynamic_cast < LCCollectionVec * > (evt->getCollection( _zsDataCollectionName ));
	CellIDDecoder<TrackerDataImpl> cellDecoder( zsInputCollectionVec );
	
	//Data Output
	//CellIDEncoder<TrackerPulseImpl> idClusterEncoder(EUTELESCOPE::ZSAPIXCLUSTERENCODING, clusterCollection);
        CellIDEncoder< TrackerPulseImpl > zsDataEncoder ( "sensorID:5,clusterID:12,xSeed:9,ySeed:10,xCluSize:9,yCluSize:9,type:5", clusterCollection );	

	//Dummy collection
	
	 bool isDummyAlreadyExisting = false;
	 LCCollectionVec * sparseClusterCollectionVec = NULL;
	 size_t dummyCollectionInitialSize = 0;
	 try 
         {
		 sparseClusterCollectionVec = dynamic_cast< LCCollectionVec* > ( evt->getCollection( "original_zsdata_apix") );
		 isDummyAlreadyExisting = true ;
	 }
         catch (lcio::DataNotAvailableException& e) 
         {
		//streamlog_out(MESSAGE4) << "Creating new dummy collection" << endl;
		sparseClusterCollectionVec =  new LCCollectionVec(LCIO::TRACKERDATA);
		isDummyAlreadyExisting = false;
	 }
	 dummyCollectionInitialSize = sparseClusterCollectionVec->size();
	  
	 CellIDEncoder<TrackerDataImpl> idClusterEncoder( "sensorID:5,clusterID:12,sparsePixelType:5,type:6", sparseClusterCollectionVec  );
	
	 //One loop for each sensorplane
	 for ( unsigned int i = 0 ; i < zsInputCollectionVec->size(); i++ ) 
         {
		TrackerDataImpl * zsData = dynamic_cast< TrackerDataImpl * > ( zsInputCollectionVec->getElementAt( i ) );
		SparsePixelType   type   = static_cast<SparsePixelType> ( static_cast<int> (cellDecoder( zsData )["sparsePixelType"]) );

		int sensorID             = static_cast<int > ( cellDecoder( zsData )["sensorID"] );
		//streamlog_out ( MESSAGE1 ) << "SensorID / Type is " << sensorID << " / " << type <<  endl;
		
	        if (type == kEUTelAPIXSparsePixel  ) 
                {
			auto_ptr<EUTelSparseDataImpl<EUTelAPIXSparsePixel > > apixData(new EUTelSparseDataImpl<EUTelAPIXSparsePixel> ( zsData ));
			//streamlog_out ( MESSAGE1 ) << "Processing APIX data on detector " << sensorID << " with " << apixData->size() << " pixels " << endl;
			// loop over all pixels in the sparseData object.
			std::vector<EUTelAPIXSparsePixel*> apixPixelVec;
			//auto_ptr<EUTelAPIXSparsePixel> apixPixel( new EUTelAPIXSparsePixel );
			EUTelAPIXSparsePixel apixPixel;
			 //Push all single Pixels of one plane in the apixPixelVec
			for ( unsigned int iPixel = 0; iPixel < apixData->size(); iPixel++ ) 
                        {
			  apixData->getSparsePixelAt( iPixel, &apixPixel);
                          EUTelAPIXSparsePixel *apix_temp = new EUTelAPIXSparsePixel(apixPixel);  
			  apixPixelVec.push_back( apix_temp );
//streamlog_out ( MESSAGE1 ) << "PixelInfo:  " << apixPixel.getXCoord() << " " << apixPixel.getYCoord() << " " << apixPixel.getSignal() << " " << apixPixel.getChip() << " " << apixPixel.getTime()<< endl;
			}
			
			//streamlog_out ( MESSAGE1 ) << "Hit Pixels: " << apixPixelVec.size() << endl;
			
			/* --- Here the real clustering happens --- */
			 int nClusters = 0;
			//int clusterNumArray[1000];
			std::vector<int> clusterNumber;
			if (!apixPixelVec.empty()) 
                        {
				clusterNumber.assign(apixPixelVec.size(), NOCLUSTER);
				//for (int i=0;i< apixPixelVec.size(); ++i) clusterNumArray[i]=NOCLUSTER;
			}
			
			
			//streamlog_out(MESSAGE1) << "Starting with clustering..." << endl;
			if (apixPixelVec.size() > 1) 
                        {
				//Now compare all pixel in one plane with each other
                		for ( unsigned int aPixel=0;aPixel < apixPixelVec.size();++aPixel) 
                                {
	 				EUTelAPIXSparsePixel *aPix = apixPixelVec.at(aPixel);
// printf("starting Pixel %5d at %5d %5d \n", aPixel, aPix->getXCoord(), aPix->getYCoord());
					for ( unsigned int bPixel=aPixel+1; bPixel < apixPixelVec.size(); ++bPixel) 
                                        {
						EUTelAPIXSparsePixel *bPix = apixPixelVec.at(bPixel);
// printf("-- neeting Pixel %5d at %5d %5d \n", bPixel, bPix->getXCoord(), bPix->getYCoord());
						//streamlog_out (MESSAGE1) << "looping..." << aPixel << " " << bPixel << endl;
						int xDist = abs( aPix->getXCoord() - bPix->getXCoord() );
						int yDist = abs( aPix->getYCoord() - bPix->getYCoord() );
						
						if ( ( xDist <= _minXDistance) && (yDist <= _minYDistance)  )
                                                { //this means they are neighbours in x-direction && / this means they are neighbours in y-direction
							bool skipPixel = false;
							if (_minDiagDistance != -1) 
                                                        {
								bool areDiagonalPartners =  (xDist != 0 && yDist != 0);  // if true then these pixels are diagonal partners
								int diagDist = max(xDist,yDist);
								if (areDiagonalPartners && diagDist > _minDiagDistance) skipPixel = true;
							}

							if (_minLVL1 != -1 )
                                                        {
								float LVL1Diff = abs(aPix->getTime() - bPix->getTime() );
								if (LVL1Diff > _minLVL1) skipPixel = true;
 							}
							
							if (skipPixel == false) 
                                                        {
								
								if ( (clusterNumber.at(aPixel) == NOCLUSTER) && clusterNumber.at(bPixel) == NOCLUSTER) 
                                                                { // none of these pixels have been assigned to a cluster
									++nClusters;
									clusterNumber.at(aPixel) = nClusters;
									clusterNumber.at(bPixel) = nClusters;
									//streamlog_out (MESSAGE1) << "assigning clusternumber: " << nClusters << endl;
								} else if ( (clusterNumber.at(aPixel) == NOCLUSTER) && clusterNumber.at(bPixel) != NOCLUSTER) { // b was assigned already, a not
									clusterNumber.at(aPixel) = clusterNumber.at(bPixel);
								} else if ( (clusterNumber.at(aPixel) != NOCLUSTER) && clusterNumber.at(bPixel) == NOCLUSTER) { // a was assigned already, b not
									clusterNumber.at(bPixel) = clusterNumber.at(aPixel);
								} else { //both pixels have a clusternumber already
									int min = std::min(clusterNumber.at(aPixel), clusterNumber.at(bPixel));
									clusterNumber.at(aPixel) = min;
									clusterNumber.at(bPixel) = min;
								}

							}
                                                        else
                                                        { // these pixels are not neighboured due to diagonal cut or too different LVL Triggers

								if ( clusterNumber.at(aPixel) == NOCLUSTER ) 
                                                                {
									++nClusters;
									clusterNumber.at(aPixel) = nClusters;
								}

								if ( clusterNumber.at(bPixel) == NOCLUSTER ) {
									++nClusters;
									clusterNumber.at(bPixel) = nClusters;
								}
							}

						}
                                                else
                                                { // these pixels are not neighboured

							if ( clusterNumber.at(aPixel) == NOCLUSTER ) 
                                                        {
								++nClusters;
								clusterNumber.at(aPixel) = nClusters;
							}

							if ( clusterNumber.at(bPixel) == NOCLUSTER ) 
                                                        {
								++nClusters;
								clusterNumber.at(bPixel) = nClusters;
							}

						}
					}
				}

			} 
                        else
                        { //You can't use the clustering algorithym with only one pixel
		
                          if (apixPixelVec.size() == 1) clusterNumber.at(0) = 1;
			}
			
			//it happens, that nCluster is higher that the real number of clusters, as a clusternumber can be deleted in the else-block, so lets correct it
			std::set<int> clusterSet;
			for (unsigned int i=0; i<apixPixelVec.size();++i) 
                        {
			  clusterSet.insert(clusterNumber.at(i));
			  if ( clusterNumber.at(i) == NOCLUSTER) 
                          {
                            streamlog_out(MESSAGE1) << "Cluster with ID -1 in " << i << endl;
                          }
			}
			nClusters = clusterSet.size();
			
			/* --- Finished Clustering --- */

			/* --- Push back one Collection per Cluster --- */
			std::set<int>::iterator it;
                        int icounter=0;
			for (it=clusterSet.begin();it!=clusterSet.end();++it) 
                        {
			  lcio::TrackerPulseImpl * pulseFrame   = new lcio::TrackerPulseImpl();
		          lcio::TrackerDataImpl  * clusterFrame = new lcio::TrackerDataImpl();
			  eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelAPIXSparsePixel >  
                                                 * apixCluster  = new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelAPIXSparsePixel >(clusterFrame);	

			  for (unsigned int i=0;i< apixPixelVec.size();i++) 
                          {
			    if(clusterNumber.at(i)== *it) 
                            { //Put only these pixels in that ClusterCollection that belong to that cluster
                              apixCluster->addSparsePixel(apixPixelVec.at(i));
			    }
  			  }

 			  int xsize0,ysize0;
			  apixCluster->getClusterSize(xsize0,ysize0);
 
			  if( xsize0 < _minXsize || xsize0 > _maxXsize || ysize0 > _maxYsize || ysize0 < _minYsize )
			  {
                            if( pulseFrame   != 0 ) delete pulseFrame;
                            if( clusterFrame != 0 ) delete clusterFrame;
                            if( apixCluster  != 0 ) delete apixCluster; 
                            continue;
		  	  }

                          bool bSensorID = true; // sensor is not defined in _ClusterLimitsMap
                          try
                          {
                            bSensorID = _ClusterLimits[0] > -1 ;
                          }
                          catch(...)
                          {
                            bSensorID = false;
                          }

                          if( bSensorID )
                          {
                            if(  
                                 xsize0 < _ClusterLimitsMap[sensorID][0] 
                                 ||
                                 xsize0 > _ClusterLimitsMap[sensorID][1] 
                                 ||
                                 ysize0 < _ClusterLimitsMap[sensorID][2] 
                                 ||
                                 ysize0 > _ClusterLimitsMap[sensorID][3] 
                               )
                                {
                                  if( pulseFrame   != 0 ) delete pulseFrame;
                                  if( clusterFrame != 0 ) delete clusterFrame;
                                  if( apixCluster  != 0 ) delete apixCluster; 
                    		  continue;
	                        }
                          }

                          if (   
                             ( apixCluster->size() >= static_cast< unsigned int >(_minNPixels) ) 
                             &&
                             ( apixCluster->getTotalCharge() >= static_cast< unsigned int >(_minCharge) )
                             )
                          {
			    float x = 0;
                            float y = 0;
		            int xsize = 0;
                            int ysize = 0;
			  
                            // This is only using digital hit information
			    // apixCluster->getCenterCoord(x,y); 
			    // This is using analogue hit information
			    apixCluster->getCenterOfGravity(x,y);
			    apixCluster->getClusterSize(xsize,ysize);
			
                            if (x > 0 && x < 10000 && y > 0 && y < 10000) 
                            {
                              // increment global cluster counter (per plane)
                              _totClusterMap[ sensorID ] += 1;
 
			      int clusterID = *it;
			      zsDataEncoder["sensorID"]      = sensorID;
			      zsDataEncoder["clusterID"]     = clusterID;
			      zsDataEncoder["xSeed"]         = static_cast< int >(x);
			      zsDataEncoder["ySeed"]         = static_cast< int >(y);
			      zsDataEncoder["xCluSize"]      = xsize;
			      zsDataEncoder["yCluSize"]      = ysize;
			      zsDataEncoder["type"]          = static_cast<int>(kEUTelAPIXClusterImpl);
			      zsDataEncoder.setCellID(pulseFrame);
			      pulseFrame->setTrackerData(clusterFrame);
			      clusterCollection->push_back(pulseFrame);
			      idClusterEncoder["sensorID"] 	     = sensorID;
			      idClusterEncoder["clusterID"] 	     = clusterID;
			      idClusterEncoder["sparsePixelType"]    = static_cast<int>(type);
			      idClusterEncoder["type"] 	             = static_cast<int>(kEUTelAPIXClusterImpl);
			      idClusterEncoder.setCellID(clusterFrame);
			      sparseClusterCollectionVec->push_back(clusterFrame);
			    }
                          }
                          if( apixCluster  != 0 ) delete apixCluster;
                          icounter++;   
                        }
		}
	}
	if ( !isDummyAlreadyExisting ) {
		if ( sparseClusterCollectionVec->size() != dummyCollectionInitialSize ) {
				evt->addCollection( sparseClusterCollectionVec, "original_zsdata_apix" );
		} else {
			delete sparseClusterCollectionVec;
		}
	}
}

void EUTelAPIXClusteringProcessor::check (LCEvent * /* evt */) {
  // nothing to check here - could be used to fill check plots in reconstruction processor
}


void EUTelAPIXClusteringProcessor::end() {
  streamlog_out ( MESSAGE2 ) <<  "Successfully finished" << endl;

  map< int, int >::iterator iter = _totClusterMap.begin();
  while ( iter != _totClusterMap.end() ) {

#ifdef MARLINDEBUG
    /// /* DEBUG */    message<DEBUG5> ( logfile << "Found " << iter->second << " clusters on detector " << iter->first );
#endif
    streamlog_out ( MESSAGE2 ) << "Found " << iter->second << " clusters on detector " << iter->first << endl;
    ++iter;

  }

}

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
void EUTelAPIXClusteringProcessor::fillHistos (LCEvent * evt) {

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

		for ( int iCluster = _initialClusterCollectionSize; iCluster < clusterCollectionVec->getNumberOfElements(); iCluster++ ) {
			TrackerPulseImpl * pulseFrame   = dynamic_cast<TrackerPulseImpl*> ( clusterCollectionVec->getElementAt(iCluster) );
			TrackerDataImpl  * clusterFrame = dynamic_cast<TrackerDataImpl*> ( pulseFrame->getTrackerData() );
			
			//ClusterType        type  = static_cast<ClusterType> ( static_cast<int> ( cellDecoder(clusterFrame)["type"] ));
			int		         sensorID  = ( static_cast<int> ( cellDecoder(pulseFrame)["sensorID"] ));
			SparsePixelType   type   = static_cast<SparsePixelType> ( static_cast<int> (cellDecoder( pulseFrame )["type"]) );
			eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelAPIXSparsePixel > *apixCluster = new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelAPIXSparsePixel >(clusterFrame);	
			
			if (static_cast< int >(type) == static_cast< int >(kEUTelAPIXClusterImpl)  ) {
			
				string tempHistoName;
				
				tempHistoName = _clusterSignalHistoName + "_d" + to_string( sensorID ) ;


				(dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(apixCluster->getTotalCharge());
 				
				int size = apixCluster->size();
				int xSize, ySize;
				apixCluster->getClusterSize(xSize,ySize);
				
				tempHistoName = _clusterSizeName + "_d" + to_string( sensorID );
				(dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(size);
				
				tempHistoName = _clusterSizeXName + "_d" + to_string( sensorID );
				(dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(xSize);
				
				tempHistoName = _clusterSizeYName + "_d" + to_string( sensorID );
				(dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(ySize);
				
				int lvl1min = -1;
				int lvl1max = -1;
				for (int iPixel=0; iPixel < size; iPixel++) {
					
					
					EUTelAPIXSparsePixel apixPixel;
					apixCluster->getSparsePixelAt(iPixel, &apixPixel);
					tempHistoName = _pixelSignalHistoName + "_d" + to_string( sensorID );
					(dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(apixPixel.getSignal());
					
					int lvl1 = static_cast< int >(apixPixel.getTime());
					
					tempHistoName = _lvl1TriggerName + "_d" + to_string( sensorID );
					(dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(lvl1);
					
					tempHistoName = _pixelMapHistoName + "_d" + to_string( sensorID );
				 	(dynamic_cast<AIDA::IHistogram2D*> (_aidaHistoMap[tempHistoName]))->fill(
                                                 static_cast<double >(apixPixel.getXCoord() ),
                                                 static_cast<double >(apixPixel.getYCoord() ),
                                                 static_cast<double >(apixPixel.getSignal() )  );
	
					if (lvl1min == -1 && lvl1max == -1) {
						lvl1min = lvl1;
						lvl1max =lvl1;
					}
					if (lvl1min > lvl1) lvl1min = lvl1;
					if (lvl1max < lvl1) lvl1max = lvl1;
				}
				
				tempHistoName = _lvl1TriggerDiffName + "_d" + to_string( sensorID );
				(dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(abs(lvl1max-lvl1min));
				
				
				tempHistoName = _hitMapHistoName + "_d" + to_string( sensorID );
				int xSeed, ySeed;
				apixCluster->getCenterCoord(xSeed, ySeed);
				(dynamic_cast<AIDA::IHistogram2D*> (_aidaHistoMap[tempHistoName]))->fill(static_cast<double >(xSeed), static_cast<double >(ySeed), 1.);
			}
			
			
			
		}

	} catch (lcio::DataNotAvailableException& e) {
		return;
	}

}

#endif

#ifdef MARLIN_USE_AIDA
void EUTelAPIXClusteringProcessor::bookHistos() {
	
	streamlog_out ( MESSAGE0 )  << "Booking histograms " << endl;
	auto_ptr<EUTelHistogramManager> histoMgr( new EUTelHistogramManager( _histoInfoFileName ));
	EUTelHistogramInfo    * histoInfo;
	bool isHistoManagerAvailable;

	try {
		isHistoManagerAvailable = histoMgr->init();
	} catch ( ios::failure& e) {
		streamlog_out ( WARNING2 ) << "I/O problem with " << _histoInfoFileName << "\n"<< "Continuing without histogram manager"  << endl;
	isHistoManagerAvailable = false;
	} catch ( ParseException& e ) {
		streamlog_out ( WARNING2 ) << e.what() << "\n"<< "Continuing without histogram manager" << endl;
		isHistoManagerAvailable = false;
	}

	string tempHistoName;
	string basePath;
	for (size_t iDetector = 0; iDetector < _sensorIDVec.size(); iDetector++) {

		int sensorID = _sensorIDVec.at( iDetector );

		// the min and max information are taken from GEAR.
		int minX, minY, maxX, maxY;
		minX = 0;
		minY = 0;

		if ( _layerIndexMap.find( sensorID ) != _layerIndexMap.end() ){
			maxX = _siPlanesLayerLayout->getSensitiveNpixelX( _layerIndexMap[ sensorID ] ) - 1;
			maxY = _siPlanesLayerLayout->getSensitiveNpixelY( _layerIndexMap[ sensorID ] ) - 1;
		} else {
			throw  InvalidGeometryException ("Unknown sensorID " + to_string( sensorID ));
		}

		basePath = "detector_" + to_string( sensorID );
		AIDAProcessor::tree(this)->mkdir(basePath.c_str());
		basePath.append("/");

		int    signalNBin  = 255;
		double signalMin   = 0.;
		double signalMax   = 256.;
		
		int clusterNBin = 10;
		double clusterMin = 0.;
		double clusterMax = 10.;
		
		int lvl1NBin = 16;
		int lvl1Min = 0;
		int lvl1Max = 15;
		
		
		string clusterTitle = "Cluster spectrum with all pixels";
		
		if ( isHistoManagerAvailable ) {
			histoInfo = histoMgr->getHistogramInfo( _clusterSignalHistoName );
			if ( histoInfo ) {
				streamlog_out ( DEBUG2 ) << (* histoInfo ) << endl;
				clusterNBin = histoInfo->_xBin;
				clusterMin  = histoInfo->_xMin;
				clusterMax  = histoInfo->_xMax;
				if ( histoInfo->_title != "" ) clusterTitle = histoInfo->_title;
			}
		}
		
	        tempHistoName = _clusterSignalHistoName + "_d" + to_string( sensorID );
		AIDA::IHistogram1D * clusterSignalHisto =  AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), signalNBin,signalMin,signalMax);
		_aidaHistoMap.insert(make_pair(tempHistoName, clusterSignalHisto));
		clusterSignalHisto->setTitle(clusterTitle.c_str());

		string singlePixelTitle = "Spectrum of single pixels";
		tempHistoName = _pixelSignalHistoName + "_d" + to_string( sensorID );
		AIDA::IHistogram1D * pixelSignalHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), signalNBin,signalMin,signalMax);
		_aidaHistoMap.insert(make_pair(tempHistoName, pixelSignalHisto));
		pixelSignalHisto->setTitle(singlePixelTitle.c_str());
		
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
		
		string lvl1triggerTitle = "LVL1 trigger distribution";
		tempHistoName = _lvl1TriggerName + "_d" + to_string( sensorID );
		AIDA::IHistogram1D * lvl1trigger = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), lvl1NBin,lvl1Min,lvl1Max);
		_aidaHistoMap.insert(make_pair(tempHistoName, lvl1trigger));
		lvl1trigger->setTitle(lvl1triggerTitle.c_str());
		
		string lvl1TriggerDiffTitle = "Max LVL1 trigger differece";
		tempHistoName = _lvl1TriggerDiffName + "_d" + to_string( sensorID );
		AIDA::IHistogram1D * lvl1triggerdiff = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), lvl1NBin,lvl1Min,lvl1Max);
		_aidaHistoMap.insert(make_pair(tempHistoName, lvl1triggerdiff));
		lvl1triggerdiff->setTitle(lvl1TriggerDiffTitle.c_str());

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
	
		tempHistoName = _pixelMapHistoName + "_d" + to_string( sensorID );
           	AIDA::IHistogram2D * pixelMapHisto =
		AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(), xBin, xMin, xMax,yBin, yMin, yMax);
		_aidaHistoMap.insert(make_pair(tempHistoName, pixelMapHisto));
		pixelMapHisto->setTitle("Pixel map");
	}
	streamlog_out ( MESSAGE0 )  << "end of Booking histograms " << endl;
}
#endif
#endif // USE_GEAR
