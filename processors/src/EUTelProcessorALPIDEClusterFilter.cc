// eutelescope includes ".h"
#include "EUTelProcessorALPIDEClusterFilter.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelSparseClusterImpl.h"

// lcio includes <.h>
#include <IMPL/TrackerPulseImpl.h>

using namespace std;
using namespace marlin;
using namespace eutelescope;
using namespace lcio;

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#endif

EUTelProcessorALPIDEClusterFilter::EUTelProcessorALPIDEClusterFilter () : Processor("EUTelProcessorALPIDEClusterFilter"),
_zsDataCollectionName(""),
_nDeep(2),
_Range(0.1),
_initialPulseCollectionSize(0),
_pulseCollectionName(""),
_sparseClusterCollectionName(""),
_totClusterMap()
{
    registerInputCollection (LCIO::TRACKERDATA, "ZSDataCollectionName",
                             "Input of Zero Suppressed data",
                             _zsDataCollectionName, string ("zsdata") );

    registerOutputCollection(LCIO::TRACKERPULSE, "PulseCollectionName",
                             "Cluster (output) collection name",
                             _pulseCollectionName, string("cluster"));

    registerOutputCollection(LCIO::TRACKERPULSE, "sparseClusterCollectionName",
                             "Cluster (output) collection name to _sparseClusterCollectionName",
                             _sparseClusterCollectionName, string("filtered_zsdata"));

}


void EUTelProcessorALPIDEClusterFilter::init(){
}


bool EUTelProcessorALPIDEClusterFilter::SameCluster(int iEvent, int iCluster, int jEvent,int jCluster)
{
	int nSame=0;
	for(unsigned int iPixel=0; iPixel<PixelsOfEvents[iEvent][iCluster].size(); iPixel++)
	{
		for(unsigned int jPixel=0; jPixel<PixelsOfEvents[jEvent][jCluster].size(); jPixel++)
		{
			if(PixelsOfEvents[iEvent][iCluster][iPixel][0]==PixelsOfEvents[jEvent][jCluster][jPixel][0] && 
			PixelsOfEvents[iEvent][iCluster][iPixel][1]==PixelsOfEvents[jEvent][jCluster][jPixel][1] &&
			PixelsOfEvents[iEvent][iCluster][iPixel][2]==PixelsOfEvents[jEvent][jCluster][jPixel][2] &&
			PixelsOfEvents[iEvent][iCluster][iPixel][3]==PixelsOfEvents[jEvent][jCluster][jPixel][3])
			{
				nSame++;
				break;
			}
		}
	}
	if(nSame>PixelsOfEvents[iEvent][iCluster].size() * _Range || nSame>PixelsOfEvents[jEvent][jCluster].size() * _Range) return true;
	return false;
}

void EUTelProcessorALPIDEClusterFilter::AddCluster(int iEvent, int iCluster, int jEvent,int jCluster)
{
	for(unsigned int jPixel=0; jPixel<PixelsOfEvents[jEvent][jCluster].size(); jPixel++)
	{
		bool samePixel=false;
		for(unsigned int iPixel=0; iPixel<PixelsOfEvents[iEvent][iCluster].size(); iPixel++)
		{
			if(PixelsOfEvents[iEvent][iCluster][iPixel][0]==PixelsOfEvents[jEvent][jCluster][jPixel][0] && 
			PixelsOfEvents[iEvent][iCluster][iPixel][1]==PixelsOfEvents[jEvent][jCluster][jPixel][1] &&
			PixelsOfEvents[iEvent][iCluster][iPixel][2]==PixelsOfEvents[jEvent][jCluster][jPixel][2] &&
			PixelsOfEvents[iEvent][iCluster][iPixel][3]==PixelsOfEvents[jEvent][jCluster][jPixel][3])
			{
				samePixel=true;
				break;
			}
		}
		if(!samePixel) PixelsOfEvents[iEvent][iCluster].push_back(PixelsOfEvents[jEvent][jCluster][jPixel]);
	}
}

void EUTelProcessorALPIDEClusterFilter::DeleteCluster(int jEvent, int jCluster){
	PixelsOfEvents[jEvent].erase(PixelsOfEvents[jEvent].begin()+jCluster);
}

void EUTelProcessorALPIDEClusterFilter::readCollections (LCCollectionVec * zsInputDataCollectionVec) {
        vector<vector<vector<int>>>EventPixels;
	for ( size_t actualCluster=0 ; actualCluster<zsInputDataCollectionVec->size(); actualCluster++) {
		CellIDDecoder<TrackerDataImpl> cellDecoder( zsInputDataCollectionVec );
		TrackerDataImpl * zsData = dynamic_cast< TrackerDataImpl * > ( zsInputDataCollectionVec->getElementAt(actualCluster) );
		SparsePixelType   type   = static_cast<SparsePixelType> ( static_cast<int> (cellDecoder( zsData )["sparsePixelType"]) );

		int clusterSize = zsData->getChargeValues().size()/4;
		vector<int> X(clusterSize);
		vector<int> Y(clusterSize);
		if ( type == kEUTelGenericSparsePixel )
		{
			vector<vector<int> > pixVector;
			auto sparseData = EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>(zsData);	
			for(size_t iPixel = 0; iPixel < sparseData.size(); iPixel++ )
			{
				auto& pixel = sparseData.at( iPixel );
				X[iPixel] = pixel.getXCoord();
				Y[iPixel] = pixel.getYCoord();
				vector<int> pix;
				pix.push_back(X[iPixel]);
				pix.push_back(Y[iPixel]);
				//pix[2] will be the sensor id.
				pix.push_back(static_cast<int>(cellDecoder(zsData)["sensorID"]));
				//pix[3] will be the pixel type.
				pix.push_back(static_cast<int>(cellDecoder( zsData )["sparsePixelType"]));
				//pix[4] will be the time (it is like an id) of the cluster
				pix.push_back(zsData->getTime());
				//pix[5] will bw Signal
				pix.push_back(pixel.getSignal());
				//pix[6] will be the time (from the pixel)
				pix.push_back(pixel.getTime());
				pixVector.push_back(pix);
			}
			EventPixels.push_back(pixVector);
		}
	}
	if(EventPixels.size()!=0) {
 	   	PixelsOfEvents.push_back(EventPixels);
        }

}

void EUTelProcessorALPIDEClusterFilter::writeCollection (LCCollectionVec * sparseClusterCollectionVec, LCCollectionVec * pulseCollection) {
	CellIDEncoder<TrackerDataImpl> idZSClusterEncoder( EUTELESCOPE::ZSCLUSTERDEFAULTENCODING, sparseClusterCollectionVec );
	CellIDEncoder<TrackerPulseImpl> idZSPulseEncoder(EUTELESCOPE::PULSEDEFAULTENCODING, pulseCollection);
	if(PixelsOfEvents.size()>_nDeep) {
		for(unsigned int iCluster=0; iCluster<PixelsOfEvents[0].size();iCluster++) {
			// prepare a TrackerData to store the cluster candidate
                        auto zsCluster = std::make_unique<TrackerDataImpl>();
                        // prepare a reimplementation of sparsified cluster
                        auto sparseCluster = std::make_unique<EUTelSparseClusterImpl<EUTelGenericSparsePixel>>(zsCluster.get());
			int sensorID = 0;
			int TYPE = 0;
			int TIME = 0;

			while(PixelsOfEvents[0][iCluster].size()>0)
			{
				EUTelGenericSparsePixel Pixel;
				Pixel.setXCoord(PixelsOfEvents[0][iCluster][0][0]);
				Pixel.setYCoord(PixelsOfEvents[0][iCluster][0][1]);
				Pixel.setTime(PixelsOfEvents[0][iCluster][0][6]);
				Pixel.setSignal(PixelsOfEvents[0][iCluster][0][5]);
				sensorID=PixelsOfEvents[0][iCluster][0][2];
				TYPE=PixelsOfEvents[0][iCluster][0][3];
				TIME=PixelsOfEvents[0][iCluster][0][4];
				PixelsOfEvents[0][iCluster].erase(PixelsOfEvents[0][iCluster].begin());
				sparseCluster->push_back( Pixel );
			}
			if ( sparseCluster->size() > 0)
			{
				// set the ID for this zsCluster
                                idZSClusterEncoder["sensorID"] = sensorID;
				idZSClusterEncoder["sparsePixelType"]= TYPE;
                	        idZSClusterEncoder["quality"] = 0;
    	                        idZSClusterEncoder.setCellID( zsCluster.get() );
                	        zsCluster->setTime(TIME);

				// add it to the cluster collection
    	                        sparseClusterCollectionVec->push_back( zsCluster.get() );

				// prepare a pulse for this cluster
                	        auto zsPulse = std::make_unique<TrackerPulseImpl>();
    	                        idZSPulseEncoder["sensorID"] = sensorID;
                	        idZSPulseEncoder["type"] = static_cast<int>(kEUTelSparseClusterImpl);
    	                        idZSPulseEncoder.setCellID( zsPulse.get() );
                	        zsPulse->setTime(TIME);
    	                        zsPulse->setTrackerData( zsCluster.release() );
                	        pulseCollection->push_back( zsPulse.release() );
				_totClusterMap[sensorID] +=1;
			}
		}
		PixelsOfEvents.erase(PixelsOfEvents.begin());
	}
}

void EUTelProcessorALPIDEClusterFilter::filter () {
	if(PixelsOfEvents.size()>_nDeep)
	{
		for(unsigned int iCluster=0; iCluster<PixelsOfEvents[0].size(); iCluster++)
		{
			for(unsigned int jEvent=1; jEvent<=_nDeep; jEvent++)
			{
				bool wasSameCluster=false;
				for(unsigned int jCluster=0; jCluster<PixelsOfEvents[jEvent].size(); jCluster++)
				{
					if(SameCluster(0,iCluster,jEvent,jCluster))
					{
						//AddCluster(0,iCluster,jEvent,jCluster);
						DeleteCluster(jEvent,jCluster);
						wasSameCluster=true;
					}
				}
				if(!wasSameCluster) break;
			}
		}
	}
}

void EUTelProcessorALPIDEClusterFilter::processEvent (LCEvent * evt) {
	EUTelEventImpl * event = static_cast<EUTelEventImpl*> (evt);
	if ( event->getEventType() == kEORE )
        {
                streamlog_out ( DEBUG4 ) <<  "EORE found: nothing else to do." <<  endl;
                return;
        }
        else if ( event->getEventType() == kUNKNOWN )
        {
                streamlog_out ( WARNING2 ) << "Event number " << event->getEventNumber()
                                   << " is of unknown type. Continue considering it as a normal Data Event." << endl;
        }
  	_clusterAvailable = true;
  	try {
	        zsInputDataCollectionVec = dynamic_cast< LCCollectionVec * > ( evt->getCollection( _zsDataCollectionName ) ) ;
	        streamlog_out ( DEBUG5 ) << "zsInputDataCollectionVec: " << _zsDataCollectionName.c_str() << " found " << endl;
		CellIDDecoder<TrackerDataImpl > cellDecoder( zsInputDataCollectionVec );
		for ( size_t i = 0; i < zsInputDataCollectionVec->size(); ++i )
                {
                       TrackerDataImpl * data = dynamic_cast< TrackerDataImpl * > ( zsInputDataCollectionVec->getElementAt( i ) ) ;
                       _totClusterMap.insert( make_pair( cellDecoder( data )[ "sensorID" ] , 0 ));
                }
  	} catch ( lcio::DataNotAvailableException ) {
	        _clusterAvailable = false;
  	}

	bool isDummyAlreadyExisting = false;
  	LCCollectionVec * sparseClusterCollectionVec = nullptr;
        try
        {
                sparseClusterCollectionVec = dynamic_cast< LCCollectionVec* > ( evt->getCollection( _sparseClusterCollectionName ) );
                isDummyAlreadyExisting = true ;
        }
        catch (lcio::DataNotAvailableException& e)
        {
                sparseClusterCollectionVec = new LCCollectionVec(LCIO::TRACKERDATA);
                isDummyAlreadyExisting = false;
        }
	LCCollectionVec * pulseCollection = nullptr;
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
	if(_clusterAvailable) {
		readCollections(zsInputDataCollectionVec);
		filter();
		writeCollection(sparseClusterCollectionVec, pulseCollection);
	}
	if ( ! isDummyAlreadyExisting )
        {
                if ( sparseClusterCollectionVec->size() != 0 )
                {
                        evt->addCollection( sparseClusterCollectionVec, _sparseClusterCollectionName );
                }
                else
                {
                        delete sparseClusterCollectionVec;
                }
        }
        if ( ! pulseCollectionExists && ( pulseCollection->size() != _initialPulseCollectionSize ))
        {
                evt->addCollection( pulseCollection, _pulseCollectionName );
        }

        if ( ! pulseCollectionExists && ( pulseCollection->size() == _initialPulseCollectionSize ) )
        {
                delete pulseCollection;
        }
}

void EUTelProcessorALPIDEClusterFilter::end() 
{
	map< int, int >::iterator iter = _totClusterMap.begin();
	while ( iter != _totClusterMap.end() ) {
    	        streamlog_out ( MESSAGE2 ) << "Found " << iter->second << " clusters on detector " << iter->first << endl;
    	        ++iter;
        }
}
