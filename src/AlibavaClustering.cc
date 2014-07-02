/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@cern.ch
 */


// eutelescope includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelSparseDataImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelSparseData2Impl.h"
#include "EUTelSparseCluster2Impl.h"


// alibava includes ".h"
#include "AlibavaClustering.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"
#include "AlibavaPedNoiCalIOManager.h"


// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <lcio.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/ClusterImpl.h>
#include <IMPL/LCEventImpl.h>


// ROOT includes ".h"
#include <TH1D.h>
#include <TH2D.h>
#include <TH1.h>
#include <TF1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>

// system includes <>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <memory>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;
using namespace eutelescope;


AlibavaClustering::AlibavaClustering () :
AlibavaBaseProcessor("AlibavaClustering"),
_clusterCollectionName(ALIBAVA::NOTSET),
_clustercharge(),
_clustercount(0),
_telescopecoordinate()
{

	// modify processor description
	_description =	"AlibavaClustering finds clusters ";

	// first of register the input collection
	registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName", "Input data collection name", _inputCollectionName, string("data") );

	// noise input
	registerProcessorParameter ("NoiseInputFile", "The filename where the final noise is stored", _pedestalFile , string("finalnoise.slcio"));

	// seed cut
	registerProcessorParameter ("SeedCut", "The seed snr cut", _seedcut , float(3.0));

	// cluster cut
	registerProcessorParameter ("ClusterCut", "The cluster snr cut", _clustercut , float(1.85));

	// polarity
	registerProcessorParameter ("Polarity", "The sensor polarity: -1 for negative cluster signals (p-type sensor), 1 for positive cluster signals (n-type sensor)", _polarity, int(-1));

	// now the optional parameters
	registerProcessorParameter ("NoiseCollectionName", "Noise collection name", _noiseCollectionName, string ("finalnoise"));

	// the name of the cluster collection
	registerProcessorParameter ("ClusterCollectionName", "Cluster collection name", _clusterCollectionName, string ("clustercollection"));

	// the unsensitve axis of the strip sensor
	registerProcessorParameter ("UnsensitiveAxis", "The unsensitive axis of our strip sensor", _nonsensitiveaxis, string ("x"));

	// the name of the sparse cluster collection
	registerProcessorParameter ("SparseClusterCollectionName", "Sparse cluster collection name, needs to be original_zsdata for hitmaker", _sparseclusterCollectionName, string ("original_zsdata"));
	
	// do we want to get the missing coordinate from the telescope here?
	registerProcessorParameter ("GetMissingCoordFromTelescope","If we want to get the missing coordinate from the telescope this should be set to the telescope plane id (usually 0-5). Set to -1 to disable this feature", _telescopePlane, int(2));

	registerOptionalParameter ("TelescopeCollectionName","Telescope collection name we want to get the unsensitive axis positions from", _telescopeCollectionName, string("cluster_m26"));

	registerOptionalParameter ("TelescopeFile", "The filename where the telescope data is stored", _telescopeFile , string("telescope.slcio"));

	registerOptionalParameter ("MaxCount", "The maximum number of events read from the telescope. This number should be below the actual number of events in the telescope. The Alibava event count can be higher.", _maxcount, int(500000));
}


void AlibavaClustering::init () {
	streamlog_out ( MESSAGE4 ) << "Running init" << endl;

	// this method is called only once even when the rewind is active
	// usually a good idea to
	printParameters ();

	
	/* To set of channels to be used
	 ex.The format should be like $ChipNumber:StartChannel-EndChannel$
	 ex. $0:5-20$ $0:30-100$ $1:50-70$
	 means from chip 0 channels between 5-20 and 30-100, from chip 1 channels between 50-70 will be used (all numbers included). the rest will be masked and not used
	 Note that the numbers should be in ascending order and there should be no space between two $ character
	 */
	if (Global::parameters->isParameterSet(ALIBAVA::CHANNELSTOBEUSED))
		Global::parameters->getStringVals(ALIBAVA::CHANNELSTOBEUSED,_channelsToBeUsed);
	else {
		streamlog_out ( MESSAGE4 ) << "The Global Parameter "<< ALIBAVA::CHANNELSTOBEUSED <<" is not set! All channels will be used!" << endl;
	}

	
	/* To choose if processor should skip masked events
	 ex. Set the value to 0 for false, to 1 for true
	 */
	if (Global::parameters->isParameterSet(ALIBAVA::SKIPMASKEDEVENTS))
		_skipMaskedEvents = bool ( Global::parameters->getIntVal(ALIBAVA::SKIPMASKEDEVENTS) );
	else {
		streamlog_out ( MESSAGE4 ) << "The Global Parameter "<< ALIBAVA::SKIPMASKEDEVENTS <<" is not set! Masked events will be used!" << endl;
	}
	
}
void AlibavaClustering::processRunHeader (LCRunHeader * rdr) {
	streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;

	// Add processor name to the runheader
	auto_ptr<AlibavaRunHeaderImpl> arunHeader ( new AlibavaRunHeaderImpl(rdr)) ;
	arunHeader->addProcessor(type());

	// get and set the number of chips
	setChipSelection( arunHeader->getChipSelection() );

	// set channels to be used (if it is defined)
	setChannelsToBeUsed();

	// set pedestal and noise values
	setPedestals();
	
	// if you want
	bookHistos();
	
	// if we are using the telescope for the missing coordinate
	if (_telescopePlane != -1)
	{
		_telescopeopen = false;
		streamlog_out (MESSAGE4) << "Looping telescope clusters on plane " << _telescopePlane << " up to " << _maxcount << " events!" << endl;

		// read in the telescope and get the coordinate we need
		for (int i=0;i<_maxcount;i++)
		{
			_telescopecoordinate.push_back(getTelescope());
			if ( i % 10000 == 0 )
				streamlog_out (MESSAGE4) << "Storing telescope coordinate in event " << i << " as " << _telescopecoordinate[i] << endl;
		}
		lcReader->close() ;
		delete lcReader ;

		streamlog_out (MESSAGE4) << "Finished looping telescope clusters on plane " << _telescopePlane << "!" << endl;
		streamlog_out (MESSAGE4) << "Found " << _telescopecoordinate.size() << " cluster elements in the telescope!" << endl;
	}

}


void AlibavaClustering::processEvent (LCEvent * anEvent) {

	if ( anEvent->getEventNumber() % 1000 == 0 )
		streamlog_out ( MESSAGE4 ) << "Looping events "<<anEvent->getEventNumber() << endl;

	AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
	
	if (_skipMaskedEvents && (alibavaEvent->isEventMasked()) ) {
		_numberOfSkippedEvents++;
		return;
	}

	// the collection we read
	LCCollectionVec * inputCollectionVec;

	// the collection we output
	LCCollectionVec * clusterCollection;

	// the sparse collection needed for the hitmaker step
	LCCollectionVec * sparseClusterCollectionVec = NULL;
	sparseClusterCollectionVec =  new LCCollectionVec(LCIO::TRACKERDATA);

	try {
		clusterCollection = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( _clusterCollectionName ) );
		streamlog_out ( DEBUG2 ) << "clusterCollection exists..." <<  endl;
	} catch ( lcio::DataNotAvailableException& e ) {
		clusterCollection = new LCCollectionVec(LCIO::TRACKERPULSE);
		streamlog_out ( DEBUG2 ) << "clusterCollection doesn't exist..." <<  endl;
	}

	// find the seed clusters on out data
	try
        {
		// give the collection vec it's data
		inputCollectionVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( getInputCollectionName() ) ) ;

		// loop over detectors
		int noOfDetector = inputCollectionVec->getNumberOfElements();
                for ( int i = 0; i < noOfDetector; ++i ) 
                {
			TrackerDataImpl * trkdata = dynamic_cast< TrackerDataImpl * > ( inputCollectionVec->getElementAt( i ) ) ;

			// The missing coordinate. This value will be added as other coordinate to all clusters found.
			float othercoordinate = 0.0;

			// if we are using the telescope for this value, the coordinate will be the same for all clusters in this event to keep sync.
			if (_telescopePlane != -1)
			{

				streamlog_out ( DEBUG0 )<< "Getting telescope in Event " << alibavaEvent->getEventNumber() << endl;
				if ((anEvent->getEventNumber()) < _telescopecoordinate.size())
				{
					othercoordinate = _telescopecoordinate[anEvent->getEventNumber()];
				}
				else
				{
					streamlog_out (ERROR1 )<< "Telescope doesn't have enought events! Setting the missing coordinate to zero from event " << alibavaEvent->getEventNumber() << " onwards!" << endl;
				}
			}

			// let the clustering begin!
			findSeedClusters(trkdata, clusterCollection, sparseClusterCollectionVec, alibavaEvent, othercoordinate);
		}

		// write an output event, even if there is no cluster in this event, this makes hitmaking and merging easier
		alibavaEvent->addCollection( sparseClusterCollectionVec, _sparseclusterCollectionName );
		alibavaEvent->addCollection( clusterCollection, getClusterCollectionName() );

        } catch ( lcio::DataNotAvailableException ) {
                // do nothing
                streamlog_out( ERROR5 ) << "Collection ("<<getInputCollectionName()<<") not found in event " << anEvent->getEventNumber() << " ! " << endl;
        }
}


void AlibavaClustering::findSeedClusters(TrackerDataImpl * trkdata, LCCollectionVec * clusterCollection, LCCollectionVec * sparseClusterCollectionVec, AlibavaEventImpl * alibavaEvent, float othercoordinate)
{

	// Magic is done here
	streamlog_out ( DEBUG0 )<< "Find Seed loop " << alibavaEvent->getEventNumber() << endl;

	// the chip number
	int chipnum = getChipNum(trkdata);

	// the data
	FloatVec datavec;
	datavec = trkdata->getChargeValues();
	
	// cluster count in this event
	int nClusters = 0;
	
	// a vector with the cluster info, 0 for no cluster and nClusters for the individual cluster, e.g. {0,0,0,0,0,0}, then  {0,1,0,0,0,0} after finding the first cluster, {0,1,0,2,2,0} after finding the second, etc.
	std::vector<int> clusterNumber;
	clusterNumber.assign(datavec.size(), 0);

	// go through all channels
	for (size_t ichan=0; ichan<datavec.size();ichan++)
	{
		streamlog_out ( DEBUG0 ) << "Channel loop " << ichan << endl;

		float noise = 0.0;
		float seedlevel = 0.0;
		float clusterlevel = 0.0;

		// this is the noise level we get from the previous steps
		if (isMasked(chipnum,ichan) == false)
		{
			// get the noise levels
			noise = getNoiseAtChannel(chipnum,ichan);
			streamlog_out ( DEBUG0 ) << "Channel " << ichan << " has a noise level of | " << noise << " | ADCs" << endl;
			streamlog_out ( DEBUG0 ) << "Channel " << ichan << " has a signal level of " << datavec[ichan] << " ADCs" << endl;

			// the level for a seed
			seedlevel = _seedcut * noise;
			streamlog_out ( DEBUG0 ) << "Seed level is " << seedlevel * _polarity << endl;

			// the level for a cluster
			clusterlevel = _clustercut * noise;
			streamlog_out ( DEBUG0 ) << "Cluster level is " << clusterlevel * _polarity << endl;
		}
		else {
			// all levels should be 0
			streamlog_out ( DEBUG0 ) << "Channel " << ichan << " is masked. Setting noise level to " << noise << endl;
		}
		
		// if a channel passes the seed cluster clut, use only positive ADCs to check! The noise levels from previous steps will always be positive!
		if ((datavec[ichan] * _polarity) >= seedlevel && ( isMasked(chipnum,ichan) == false) )
		{

			streamlog_out ( DEBUG0 ) << "Found seed on channel " << ichan << " with seed of " << datavec[ichan] * _polarity << endl;

			// now check if any neighbours pass the cluster cut
			// the cluster starts out with size 1, so nothing left or right
			int posclustersize = 0;
			int negclustersize = 0;

			// look right
			bool posdir = true;
			// does this channel exist?
			while ( (ichan+posclustersize+1) < datavec.size() && posdir == true )
			{
				if ( (datavec[ichan+posclustersize+1] * _polarity) >= (_clustercut * getNoiseAtChannel(chipnum,ichan+posclustersize+1)) && ( isMasked(chipnum,ichan+posclustersize+1) == false) )
				{

					posclustersize += 1;
					streamlog_out (DEBUG3) << "Found large cluster in event: " << alibavaEvent->getEventNumber() << endl;
					
					// if this happens, we have found a second seed in this cluster with a higher ADC count
					// since ichan goes to the right, this feature can only be found here.
					
					// don't compare ACDs, look for ADC/seedlevel
					if (((datavec[ichan+posclustersize+1] * _polarity)/(getNoiseAtChannel(chipnum,ichan+posclustersize+1) * _seedcut)  ) > ((datavec[ichan] * _polarity) / seedlevel) )
					{

						streamlog_out (DEBUG4) << "Found multiple seeds in a cluster in event " << alibavaEvent->getEventNumber() << " ... splitting!" << endl;
						streamlog_out (DEBUG4) << ((datavec[ichan+posclustersize+1] * _polarity)/(getNoiseAtChannel(chipnum,ichan+posclustersize+1) * _seedcut)  ) << " > " << ((datavec[ichan] * _polarity) / seedlevel) << endl;
						
						// two options: if they touch, we simply move the seed position to the right one and inc negativesize, if they don't touch, we split:
						if (posclustersize == 1)
						{

						  streamlog_out (DEBUG5) << "Moving the seed in event: " << alibavaEvent->getEventNumber() << endl;

						  // move ichan and negclustersize
						  ichan++;
						  negclustersize++;

						  // now posclustersize has to go down, since we shifted the cluster seed position 
						  posclustersize--;
						}
						
						// if there are clusters below the seed between two seeds, we give them to the highest and split
						if (posclustersize >= 2)
						{
							streamlog_out (DEBUG5) << "Splitting the cluster between seeds in event: " << alibavaEvent->getEventNumber() << endl;

							// setting this to zero, let the higher seed find the inbetween strip when looping
							posclustersize = 0;
							posdir = false;
						}
						
						if (posclustersize >= 3)
						{
							// just report 
							streamlog_out (WARNING5) << "Multiple clusters between seeds in event: " << alibavaEvent->getEventNumber() << " ! The noise levels could be false!" << endl;
						}
						
					}
				} else {
					posdir = false;
				}
			}

			// look left
			bool negdir = true;
			// does this channel exist?
			while ( negdir == true )
			{
				// the channel has to fullfill 3 conditions: a) be above cluster level, b) not be masked and c) not be above the seed cut. c) is so that we don't refind clusters after splitting.
				if ((datavec[ichan-negclustersize-1] * _polarity) >= (_clustercut * getNoiseAtChannel(chipnum,ichan-negclustersize-1)) && ( isMasked(chipnum,ichan-negclustersize-1) == false) && (datavec[ichan-negclustersize-1] * _polarity) < (_seedcut * getNoiseAtChannel(chipnum,ichan-negclustersize-1)))
				{
				  streamlog_out (DEBUG3) << "Found large cluster in event: " << alibavaEvent->getEventNumber() << endl;
					negclustersize += 1;
				} else {
					negdir = false;
				}
			}

			// increment cluster count
			nClusters++;

			// save the cluster information
			int lowerlimit = ichan - negclustersize;
			int upperlimit = ichan + posclustersize;
			int clusize = upperlimit - lowerlimit;

			// debug output
			if (clusize > 1)
			{
				streamlog_out ( DEBUG0 ) << "Clustersize: " << clusize << endl;
			}
			
			// fill the clusterNumber vector with the position of the individual clusters
			for (int icluster=lowerlimit; icluster<=upperlimit; icluster++)
			{
				clusterNumber.at(icluster) = nClusters;
				streamlog_out ( DEBUG0 ) << "Wrote cluster to clusterNumber at position(s): " <<  icluster << endl;
			}
			
			// record the charge to the left and right of the seed to spot asymetries:
			if (isMasked(chipnum,ichan-2) == false && isMasked(chipnum,ichan-1) == false && isMasked(chipnum,ichan+1) == false && isMasked(chipnum,ichan+2) == false)
			{
				_clustercharge[0] += fabs(datavec[ichan-2]);
				_clustercharge[1] += fabs(datavec[ichan-1]);
				_clustercharge[2] += fabs(datavec[ichan]);
				_clustercharge[3] += fabs(datavec[ichan+1]);
				_clustercharge[4] += fabs(datavec[ichan+2]);
				fillChargeDistHisto(datavec[ichan-3]*_polarity,datavec[ichan-2]*_polarity,datavec[ichan-1]*_polarity,datavec[ichan]*_polarity,datavec[ichan+1]*_polarity,datavec[ichan+2]*_polarity,datavec[ichan+3]*_polarity);
			}

			// let's calculate the eta distribution, only works for clustersize >1
			if (clusize > 1)
			{
				float etaleft = 0;
				float etaright = 0;
				float cogpos = 0;
				float totalsignal = 0;
				float weight = 0;
				for (int i = lowerlimit;i<=upperlimit;i++)
				{
					totalsignal += fabs(datavec[i]);
					weight += fabs(datavec[i]*i);
				}
				cogpos = weight / totalsignal;
				for (int j = lowerlimit;j<=upperlimit;j++)
				{
					if (float (j) <cogpos)
					{
					  etaleft += fabs(datavec[j]);
					}
					else if (float (j)>cogpos)
					{
					  etaright += fabs(datavec[j]);
					}
				}
				float etaratio = etaleft/(etaleft+etaright);
				fillEtaHisto(etaratio);
			}

			// do an alternative calculation of eta:
			// require good channels left and right so we don't bias
			if (isMasked(chipnum,ichan-1) == false && isMasked(chipnum,ichan+1) == false)
			{
				float etaleft_a = 0;
				float etaright_a = 0;
				float etaleft_b = 0;
				float etaright_b = 0;
				float etaratio = -1;

				// divison by noise can be added...
				etaleft_a = datavec[ichan-1] * _polarity;// / getNoiseAtChannel(chipnum,ichan-1);
				etaright_a = datavec[ichan] * _polarity;// / getNoiseAtChannel(chipnum,ichan);

				etaright_b = datavec[ichan+1] * _polarity;// / getNoiseAtChannel(chipnum,ichan+1);
				etaleft_b = datavec[ichan] * _polarity;// / getNoiseAtChannel(chipnum,ichan);

				if (etaleft_a > etaright_b)
				{
					etaratio = etaleft_a/(etaleft_a+etaright_a);
					fillEtaHisto2(etaratio);
					streamlog_out (DEBUG2) << "Eta2: left side ratio written: " << etaratio << endl; 
				} else if (etaright_b > etaleft_a)
				{
					etaratio = etaleft_b/(etaleft_b+etaright_b);
					fillEtaHisto2(etaratio);
					streamlog_out (DEBUG2) << "Eta2: right side ratio written:" << etaratio << endl;
				} else {
					streamlog_out (DEBUG2) << "Eta2: no eta here (" << alibavaEvent->getEventNumber()<< ") because fail! ela:" << etaleft_a << " era: " << etaright_a << " elb: "<< etaleft_b << " erb: " << etaright_b << endl;
				}
			}

			// fill the cluster signal and snr histos
			float clustersignal = 0.0;
			float clusternoise = 0.0;
			for (int i=(lowerlimit); i<=(upperlimit); i++)
			{
				clustersignal += _polarity*datavec[i];
				clusternoise += getNoiseAtChannel(chipnum,i);
			}
			fillSignalHisto(clustersignal);
			fillSNRHisto(clustersignal/clusternoise);
			
			// fill the seed charge histogram
			fillSeedChargeHisto(_polarity*datavec[ichan]);

			// the overall cluster count in this run
			_clustercount++;

			// fill the clustersize histo
			fillClusterHisto(clusize);

			// fill the hitmap histogram
			fillHitmapHisto(ichan, negclustersize, posclustersize);

			// fill the seed histogram
			fillSeedHisto(ichan);

			// final: move the ichan var so that we dont include a seed channel in multiple clusters
			if ( (ichan + posclustersize + 1) < datavec.size() )
			{
				// this channel is the first one next to a previous cluster and below
				// the clustercut so it will fail the  seed check in the next check
				ichan = ichan + posclustersize +1;
			}

		} // done going over all seeds
	} // done going over all channels in this event

	// more debug output
	if (nClusters >= 1)
	{
		streamlog_out(DEBUG0) << "Clusters in this event: " << nClusters << endl;
	}
	
	// now output the clusters we have found
	
	// data -> check fields FIXME
	
	// The implementation from cms pixel:
	//CellIDEncoder< TrackerPulseImpl > zsDataEncoder ( "sensorID:10,clusterID:12,xSeed:9,ySeed:10,xCluSize:9,yCluSize:9,type:5", clusterCollection );
	
	// The implementation from eutel cluster filter:
	//CellIDEncoder< TrackerPulseImpl > zsDataEncoder ( EUTELESCOPE::PULSEDEFAULTENCODING, clusterCollection );
	
	// The dump from a telescope cluster
	CellIDEncoder< TrackerPulseImpl > zsDataEncoder ( "sensorID:5,clusterID:12,xSeed:12,ySeed:12,xCluSize:5,yCluSize:5,type:5", clusterCollection );
	
	// The implementation from cms pixel:
	//CellIDEncoder<TrackerDataImpl> idClusterEncoder( "sensorID:10,clusterID:12,sparsePixelType:5,type:6", sparseClusterCollectionVec  );
	
	// The implementation from a telescope cluster:
	CellIDEncoder<TrackerDataImpl> idClusterEncoder( "sensorID:5,clusterID:12,sparsePixelType:5,quality:5", sparseClusterCollectionVec  );

	// What have we found?
	for (unsigned int icluster=0; icluster < clusterNumber.size(); icluster++)
	{
		streamlog_out (DEBUG0) << "Pos: " << icluster << " , Value: " << clusterNumber.at(icluster) << endl;
	}
	
	// decode datavec into pixelvec, each strip goes into a Pixel
	// move to positive ADCs only, since identification has already happend. This is for EUTelSimpleSparsePixel.
	std::vector<eutelescope::EUTelSimpleSparsePixel*> PixelVec;
	eutelescope::EUTelSimpleSparsePixel Pixel;
	auto_ptr<eutelescope::EUTelSparseDataImpl<eutelescope::EUTelSimpleSparsePixel > > pixelData(new eutelescope::EUTelSparseDataImpl<eutelescope::EUTelSimpleSparsePixel> ( trkdata ));
	for ( unsigned int iPixel = 0; iPixel < datavec.size(); iPixel++ )
	{
		pixelData->getSparsePixelAt( iPixel, &Pixel);

		// select the sensor orientation, we give each "pixel" the coordinate 1 on the unsensitive axis
		if (_nonsensitiveaxis == "x")
		{
			Pixel.setXCoord(othercoordinate);
			Pixel.setYCoord(iPixel);
		}
		if (_nonsensitiveaxis == "y")
		{
			Pixel.setXCoord(iPixel);
			Pixel.setYCoord(othercoordinate);
		}

		Pixel.setSignal(datavec[iPixel] * _polarity);
		PixelVec.push_back(new EUTelSimpleSparsePixel(Pixel));

		streamlog_out(DEBUG0) << "Decoded pixel with coordinates: ( " << Pixel.getXCoord() << " , " << Pixel.getYCoord() << " ) and charge: " << datavec[iPixel] << " !" << endl;
		streamlog_out(DEBUG0) << "Input was: " << PixelVec.at(iPixel) << endl;
	}

	// now let's move this out into trackerpulses and trackerdata
	// iterate over the clusters
	std::vector<int>::iterator it;
	for (it = clusterNumber.begin(); it != clusterNumber.end(); ++it)
	{
		lcio::TrackerPulseImpl * pulseFrame = new lcio::TrackerPulseImpl();
		lcio::TrackerDataImpl * clusterFrame = new lcio::TrackerDataImpl();
		eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelSimpleSparsePixel > *pixelCluster = new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelSimpleSparsePixel >(clusterFrame);
		for (unsigned int i=0;i< datavec.size();i++)
		{
			if((clusterNumber.at(i)== *it) && (clusterNumber.at(i) > 0))
			{
				// put only these pixels in that ClusterCollection that belong to that cluster
				pixelCluster->addSparsePixel(PixelVec.at(i));
				streamlog_out( DEBUG1 ) << "Adding strip " << i << " to cluster " << clusterNumber.at(i) << endl;
				++it;
			}
		
		}
		
		// now we have pixelClusters with the corresponding pixels (which have the xy and q info) in them
		// this if stops making the "zero" cluster with all strips not in a cluster...
		if ( (pixelCluster->size() >= 1) && (pixelCluster->getTotalCharge() >= 1) )
		{

			// make a frame of each cluster and give it position, charge, etc. Then push back into clusterCollection and sparseClusterCollectionVec.
			float x,y=0;
			int xsize,ysize=0;
			pulseFrame->setCharge(pixelCluster->getTotalCharge());
			pixelCluster->getCenterOfGravity(x,y);
			pixelCluster->getClusterSize(xsize,ysize);

			streamlog_out( DEBUG6 ) << "Cluster:: Cl: " << *it << ", Q: " << pixelCluster->getTotalCharge() << " , x: " << x << " , y: " << y << " , dx: " << xsize << " , dy: " << ysize << " in event: " << alibavaEvent->getEventNumber() << endl;

			// for telescope hitmaker step: set sensor ID to chipnum+6, telescope should be 0-5. This assumes sensorid 6+chipnum in the gear file...
			int clusterID = *it;
			zsDataEncoder["sensorID"] = chipnum+6;
			zsDataEncoder["clusterID"] = clusterID;
			zsDataEncoder["xSeed"] = static_cast< long >(x);
			zsDataEncoder["ySeed"] = static_cast< long >(y);
			zsDataEncoder["xCluSize"] = xsize;
			zsDataEncoder["yCluSize"] = ysize;
			zsDataEncoder["type"] = static_cast<int>(kEUTelSparseClusterImpl);
			zsDataEncoder.setCellID(pulseFrame);
			pulseFrame->setTrackerData(clusterFrame);
			clusterCollection->push_back(pulseFrame);

			idClusterEncoder["sensorID"] = chipnum+6;
			idClusterEncoder["clusterID"] = clusterID;
			
			// set this to 1 for hitmaker
			idClusterEncoder["sparsePixelType"] = static_cast<int>(1);
			//idClusterEncoder["type"] = static_cast<int>(kEUTelSparseClusterImpl);
			idClusterEncoder["quality"] = static_cast<int>(0);
			idClusterEncoder.setCellID(clusterFrame);
			sparseClusterCollectionVec->push_back(clusterFrame);

		}
		else // clean up
		{
			delete pulseFrame;
			delete clusterFrame;
		}
		
		delete pixelCluster;

	} // done clusterNumber iteration

	// we don't need the decoding any more, cleanup
	for (std::vector<eutelescope::EUTelSimpleSparsePixel*>::iterator pixit = PixelVec.begin(); pixit < PixelVec.end();pixit++)
	{
		delete *pixit;
	}

}

void AlibavaClustering::check (LCEvent * /* evt */ ) {
	// nothing to check here - could be used to fill check plots in reconstruction processor
}

void AlibavaClustering::fillHitmapHisto(int ichan, int negclustersize, int posclustersize)
{
	string tempHistoName = "Cluster Hitmap";
	if ( TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[tempHistoName]) )
	{
		for ( int i= ( ichan - negclustersize ) ; i <= ( ichan + posclustersize ) ; i++ )
		{
			histo->Fill(i);
		}
	}
}

void AlibavaClustering::fillClusterHisto(int clusize)
{
	string tempHistoName = "Cluster Size";
	if ( TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[tempHistoName]) )
	{
		histo->Fill(clusize);
	}
}

void AlibavaClustering::fillEtaHisto(float etaratio)
{
	string tempHistoName = "Cluster COG Eta";
	if ( TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[tempHistoName]) )
	{
		histo->Fill(etaratio);
	}
}

void AlibavaClustering::fillEtaHisto2(float etaratio)
{
	string tempHistoName = "Seed and Neighbour Eta";
	if ( TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[tempHistoName]) )
	{
		histo->Fill(etaratio);
	}
}

void AlibavaClustering::fillSeedHisto(int ichan)
{
	string tempHistoName = "Seed Hitmap";
	if ( TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[tempHistoName]) )
	{
		histo->Fill(ichan);
	}
}

void AlibavaClustering::fillChargeDistHisto(float a, float b, float c, float d, float e, float f, float g)
{
	string tempHistoName = "Neighbour Charge Distribution";
	if ( TH2D * histo = dynamic_cast<TH2D*> (_rootObjectMap[tempHistoName]) )
	{
		//histo->Fill(-3.0,a/d);
		histo->Fill(-2.0,b/d);
		histo->Fill(-1.0,c/d);
		histo->Fill(0.0,d/d);
		histo->Fill(1.0,e/d);
		histo->Fill(2.0,f/d);
		//histo->Fill(3.0,g/d);
	}
}

void AlibavaClustering::fillSignalHisto(float signal)
{
	string tempHistoName = "Signal from Clusters";
	if ( TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[tempHistoName]) )
	{
		histo->Fill(signal);
	}
}

void AlibavaClustering::fillSNRHisto(float signal)
{
	string tempHistoName = "SNR from Clusters";
	if ( TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[tempHistoName]) )
	{
		histo->Fill(signal);
	}
}

void AlibavaClustering::fillSeedChargeHisto(float signal)
{
	string tempHistoName = "Signal from Seeds";
	if ( TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[tempHistoName]) )
	{
		histo->Fill(signal);
	}
}

void AlibavaClustering::end()
{
	dolandaugausfit("Signal from Seeds");
	dolandaugausfit("Signal from Clusters");

	streamlog_out ( MESSAGE4 ) << "===============================" << endl;
	streamlog_out ( MESSAGE4 ) << "Left 2 : " << _clustercharge[0] / fabs(_clustercharge[2]) << endl;
	streamlog_out ( MESSAGE4 ) << "Left 1 : " << _clustercharge[1] / fabs(_clustercharge[2]) << endl;
	streamlog_out ( MESSAGE4 ) << "Seed   : " << _clustercharge[2] / _clustercharge[2] << endl;
	streamlog_out ( MESSAGE4 ) << "Right 1: " << _clustercharge[3] / fabs(_clustercharge[2]) << endl;
	streamlog_out ( MESSAGE4 ) << "Right 2: " << _clustercharge[4] / fabs(_clustercharge[2]) << endl;
	streamlog_out ( MESSAGE4 ) << "===============================" << endl;
	streamlog_out ( MESSAGE4 ) << "Clustercount:" << _clustercount << endl;
	streamlog_out ( MESSAGE4 ) << "===============================" << endl;
	streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
	
}

void AlibavaClustering::fillHistos()
{


}


void AlibavaClustering::bookHistos()
{
	// a histogram showing the clustersize
	string tempHistoName = "Cluster Size";
	stringstream tempHistoTitle;
	tempHistoTitle << tempHistoName << ";Clustersize;NumberofEntries";

	TH1D * clusterHisto = new TH1D (tempHistoName.c_str(),"",4,1,5);
	_rootObjectMap.insert(make_pair(tempHistoName, clusterHisto));
	string tmp_string = tempHistoTitle.str();
	clusterHisto->SetTitle(tmp_string.c_str());

	// a histogram showing the eta distribution
	string tempHistoName2 = "Cluster COG Eta";
	stringstream tempHistoTitle2;
	tempHistoTitle2 << tempHistoName2 << ";Eta;NumberofEntries";

	TH1D * etaHisto = new TH1D (tempHistoName2.c_str(),"",100,0,1);
	_rootObjectMap.insert(make_pair(tempHistoName2, etaHisto));
	string tmp_string2 = tempHistoTitle2.str();
	etaHisto->SetTitle(tmp_string2.c_str());

	// a hitmap histo
	string tempHistoName3 = "Cluster Hitmap";
	stringstream tempHistoTitle3;
	tempHistoTitle3 << tempHistoName3 << ";Channel;NumberofEntries";

	TH1D * hitmapHisto = new TH1D (tempHistoName3.c_str(),"",128,0,127);
	_rootObjectMap.insert(make_pair(tempHistoName3, hitmapHisto));
	string tmp_string3 = tempHistoTitle3.str();
	hitmapHisto->SetTitle(tmp_string3.c_str());

	// a seed histo
	string tempHistoName4 = "Seed Hitmap";
	stringstream tempHistoTitle4;
	tempHistoTitle4 << tempHistoName4 << ";Channel;NumberofEntries";

	TH1D * seedHisto = new TH1D (tempHistoName4.c_str(),"",128,0,127);
	_rootObjectMap.insert(make_pair(tempHistoName4, seedHisto));
	string tmp_string4 = tempHistoTitle4.str();
	seedHisto->SetTitle(tmp_string4.c_str());

	// a histogram showing the eta distribution - alternative calculation
	string tempHistoName5 = "Seed and Neighbour Eta";
	stringstream tempHistoTitle5;
	tempHistoTitle5 << tempHistoName5 << ";Eta;NumberofEntries";

	TH1D * etaHisto2 = new TH1D (tempHistoName5.c_str(),"",100,0,1);
	_rootObjectMap.insert(make_pair(tempHistoName5, etaHisto2));
	string tmp_string5 = tempHistoTitle5.str();
	etaHisto2->SetTitle(tmp_string5.c_str());

	// a histogram showing the charge distribution
	string tempHistoName6 = "Neighbour Charge Distribution";
	stringstream tempHistoTitle6;
	tempHistoTitle6 << tempHistoName6 << ";Distance to seed;Strip charge / Seed charge";

	TH2D * chargehisto = new TH2D (tempHistoName6.c_str(),"",7,-3.5,3.5,1000,-1,1);
	_rootObjectMap.insert(make_pair(tempHistoName6, chargehisto));
	string tmp_string6 = tempHistoTitle6.str();
	chargehisto->SetTitle(tmp_string6.c_str());

	// a histogram showing the cluster signals
	string tempHistoName7 = "Signal from Clusters";
	stringstream tempHistoTitle7;
	tempHistoTitle7 << tempHistoName7 << ";Cluster signal (ADCs) * (-1);Number of Entries";

	TH1D * clustersignalhisto = new TH1D (tempHistoName7.c_str(),"",1000,0,100);
	_rootObjectMap.insert(make_pair(tempHistoName7, clustersignalhisto));
	string tmp_string7 = tempHistoTitle7.str();
	clustersignalhisto->SetTitle(tmp_string7.c_str());

	// a histogram showing the cluster signal to noise ratio
	string tempHistoName8 = "SNR from Clusters";
	stringstream tempHistoTitle8;
	tempHistoTitle8 << tempHistoName8 << ";Cluster SNR;Number of Entries";

	TH1D * clustersnrhisto = new TH1D (tempHistoName8.c_str(),"",500,0,50);
	_rootObjectMap.insert(make_pair(tempHistoName8, clustersnrhisto));
	string tmp_string8 = tempHistoTitle8.str();
	clustersnrhisto->SetTitle(tmp_string8.c_str());

	// a histogram showing the seed charge
	string tempHistoName9 = "Signal from Seeds";
	stringstream tempHistoTitle9;
	tempHistoTitle9 << tempHistoName9 << ";Seed charge (ADCs) * (-1);Number of Entries";

	TH1D * seedchargehisto = new TH1D (tempHistoName9.c_str(),"",1000,0,100);
	_rootObjectMap.insert(make_pair(tempHistoName9, seedchargehisto));
	string tmp_string9 = tempHistoTitle9.str();
	seedchargehisto->SetTitle(tmp_string9.c_str());

	streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << endl;
}



// getter and setter for _clusterCollectionName
void AlibavaClustering::setClusterCollectionName(std::string clusterCollectionName){
	_clusterCollectionName = clusterCollectionName;
}
std::string AlibavaClustering::getClusterCollectionName(){
	return _clusterCollectionName;
}




float AlibavaClustering::getTelescope()
{

  	// the input collections
	LCCollectionVec * telescopeCollectionVec;

	// the telescope is read by the function
	LCEvent* evt= readTelescope();

	try{
	// and the secondary collections
	telescopeCollectionVec = dynamic_cast< LCCollectionVec * > ( evt->getCollection( _telescopeCollectionName ) ) ;
	
	} catch ( lcio::DataNotAvailableException& e ) {
	    
	    // FIXME: Sometimes an EOF is not caught, it doesn't affect the actual clustering.

                streamlog_out( DEBUG0 ) << "Collection fail" << endl;
		telescopeCollectionVec = new LCCollectionVec(LCIO::TRACKERPULSE);
		
		return(0.0);
        }
	
	int telescopesize = telescopeCollectionVec->getNumberOfElements();
	streamlog_out ( DEBUG0 ) << telescopesize << " Elements in this telescope event!" << endl;

	CellIDDecoder<TrackerPulseImpl> inputDecoder(telescopeCollectionVec);
	
	float telescope_x = 0.0;
	float telescope_y = 0.0;
	int planecount = 0;
	float position_x = 0.0;
	float position_y = 0.0;

	// loop over the telescope cluster event
	for (int j = 0; j<telescopesize; j++)
	{ 
		streamlog_out ( DEBUG0 ) << "Reading element " << j << " of " << telescopesize << " in this telescope event!" << endl;
		lcio::TrackerPulseImpl * input  = dynamic_cast< lcio::TrackerPulseImpl * > ( telescopeCollectionVec->getElementAt( j ) ) ;

		// decode our values into these ints ...
		int sensorID = inputDecoder(input)["sensorID"];
		int xSeed = inputDecoder(input)["xSeed"];
		int ySeed = inputDecoder(input)["ySeed"];

		// We only want one plane
		if (sensorID == _telescopePlane)
		{
			streamlog_out (DEBUG0) << "Telescope x read: " << xSeed << endl;
			streamlog_out (DEBUG0) << "Telescope y read: " << ySeed << endl;
			telescope_x += xSeed;
			telescope_y += ySeed;
			planecount++;
		}
	}

	if (planecount > 0)
	{
		streamlog_out (DEBUG0) << "The telescope plane " << _telescopePlane << " has " << planecount << " clusters on it!" << endl;
		position_x = telescope_x/planecount;
		position_y = telescope_y/planecount;
		streamlog_out (DEBUG0) << "Telescope x total: " << position_x << endl;
		streamlog_out (DEBUG0) << "Telescope y total: " << position_y << endl;
	}
	
	if ( _nonsensitiveaxis == "x" )
	{
		streamlog_out ( DEBUG4 ) << "Wrote telescope x position: " << position_x << endl;
		return(position_x) ;
	}
	else if ( _nonsensitiveaxis == "y" )
	{
		streamlog_out ( DEBUG4 ) << "Wrote telescope y position: " << position_y << endl;
		return(position_y) ;
	}
	else
	{
		streamlog_out ( ERROR1 ) << "_nonsensitiveaxis is set to an invalid value: " << _nonsensitiveaxis << " ! Returning 0 for this coordinate from the telescope!" << endl;
		return(0.0);
	}
}

// the telescope file is read here:
LCEvent *AlibavaClustering::readTelescope ()
{
	if (_telescopeopen == false)
	{
		lcReader = LCFactory::getInstance()->createLCReader( IO::LCReader::directAccess ) ;
		try
		{
			lcReader->open( _telescopeFile ) ;
			_telescopeopen = true;
		}
		catch( IOException& e )
		{
			streamlog_out ( ERROR1 ) << "Can't open the telescope file: " << e.what() << endl ;
		}
	}

	try
	{
		LCEvent *evt = lcReader->readNextEvent();
		if (evt == NULL)
		{
		  return(0);
		  streamlog_out ( ERROR1 ) << "FAIL: " << endl ;
		}
		return(evt);
	}
	catch ( IOException& e )
	{
		streamlog_out ( ERROR1 ) << "FAIL: " << e.what() << endl ;
		return(0);
	}
}


void AlibavaClustering::dolandaugausfit(string tempHistoName)
{
  	streamlog_out (DEBUG5) << "Fitting landau gaus on histo " << tempHistoName << endl;
	if ( TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[tempHistoName]) )
	{
		// Setting fit range and start values
		Double_t fr[2];
		Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
		fr[0]=0.7*histo->GetMean();
		fr[1]=5.0*histo->GetMean();
		pllo[0]=0.05; pllo[1]=0.50; pllo[2]=0.1; pllo[3]=0.04;
		plhi[0]=50.0; plhi[1]=500.0; plhi[2]=1000000.0; plhi[3]=50.0;
		sv[0]=1.8; sv[1]=20.0; sv[2]=50000.0; sv[3]=3.0;

		Double_t chisqr;
		Int_t    ndf;
		TF1 *fitsnr = langaufit(histo,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
		histo->Fit(fitsnr,"QR");

		Double_t SNRPeak, SNRFWHM;
		langaupro(fp,SNRPeak,SNRFWHM);

		streamlog_out ( MESSAGE4 ) << "===============================" << endl;
		streamlog_out ( DEBUG5 ) << "Done fitting landau gaus on histo" << tempHistoName << endl;
		streamlog_out ( MESSAGE4 ) << "Landau-Gaus Peak for : " << tempHistoName << " is: " << SNRPeak << endl;
		streamlog_out ( MESSAGE4 ) << "Landau-Gaus Sigma for: " << tempHistoName << " is: " << SNRFWHM << endl;
		streamlog_out ( MESSAGE4 ) << "===============================" << endl;
	}

}


Double_t langaufun(Double_t *x, Double_t *par)
{
	//Fit parameters:
	//par[0]=Width (scale) parameter of Landau density
	//par[1]=Most Probable (MP, location) parameter of Landau density
	//par[2]=Total area (integral -inf to inf, normalization constant)
	//par[3]=Width (sigma) of convoluted Gaussian function
	//
	//In the Landau distribution (represented by the CERNLIB approximation), 
	//the maximum is located at x=-0.22278298 with the location parameter=0.
	//This shift is corrected within this function, so that the actual
	//maximum is identical to the MP parameter.

	// Numeric constants
	Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
	Double_t mpshift  = -0.22278298;       // Landau maximum location

	// Control constants
	Double_t np = 100.0;      // number of convolution steps
	Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

	// Variables
	Double_t xx;
	Double_t mpc;
	Double_t fland;
	Double_t sum = 0.0;
	Double_t xlow,xupp;
	Double_t step;
	Double_t i;

	// MP shift correction
	mpc = par[1] - mpshift * par[0]; 

	// Range of convolution integral
	xlow = x[0] - sc * par[3];
	xupp = x[0] + sc * par[3];

	step = (xupp-xlow) / np;

	// Convolution integral of Landau and Gaussian by sum
	for(i=1.0; i<=np/2; i++) {
		xx = xlow + (i-.5) * step;
		fland = TMath::Landau(xx,mpc,par[0]) / par[0];
		sum += fland * TMath::Gaus(x[0],xx,par[3]);
		xx = xupp - (i-.5) * step;
		fland = TMath::Landau(xx,mpc,par[0]) / par[0];
		sum += fland * TMath::Gaus(x[0],xx,par[3]);
	}

	return (par[2] * step * sum * invsq2pi / par[3]);
}

TF1 *langaufit(TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
	// Variables for langaufit call:
	//   his             histogram to fit
	//   fitrange[2]     lo and hi boundaries of fit range
	//   startvalues[4]  reasonable start values for the fit
	//   parlimitslo[4]  lower parameter limits
	//   parlimitshi[4]  upper parameter limits
	//   fitparams[4]    returns the final fit parameters
	//   fiterrors[4]    returns the final fit errors
	//   ChiSqr          returns the chi square
	//   NDF             returns ndf

	Int_t i;
	Char_t FunName[100];

	sprintf(FunName,"Fitfcn_%s",his->GetName());

	TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
	if (ffitold) delete ffitold;

	TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
	ffit->SetParameters(startvalues);
	ffit->SetParNames("Width","MP","Area","GSigma");
	
	for (i=0; i<4; i++) {
	    ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
	}

	// fit within specified range, use ParLimits, do not plot
	his->Fit(FunName,"RB0Q");

	// obtain fit parameters
	ffit->GetParameters(fitparams);
	for (i=0; i<4; i++) {
		// obtain fit parameter errors
		fiterrors[i] = ffit->GetParError(i);
	}
	// obtain chi^2
	ChiSqr[0] = ffit->GetChisquare();
	// obtain ndf
	NDF[0] = ffit->GetNDF();
	// return fit function
	return (ffit);
}

Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM)
{
	// Seaches for the location (x value) at the maximum of the 
	// Landau-Gaussian convolute and its full width at half-maximum.
	// The search is probably not very efficient, but it's a first try.

	Double_t p,x,fy,fxr,fxl;
	Double_t step;
	Double_t l,lold;
	Int_t i = 0;
	Int_t MAXCALLS = 10000;

	// Search for maximum
	p = params[1] - 0.1 * params[0];
	step = 0.05 * params[0];
	lold = -2.0;
	l    = -1.0;

	while ( (l != lold) && (i < MAXCALLS) )
	{
		i++;
		lold = l;
		x = p + step;
		l = langaufun(&x,params);
      
		if (l < lold)
			step = -step/10;

		p += step;
	}

	if (i == MAXCALLS)
		return (-1);

	maxx = x;
	fy = l/2;

	// Search for right x location of fy
	p = maxx + params[0];
	step = params[0];
	lold = -2.0;
	l    = -1e300;
	i    = 0;

	while ( (l != lold) && (i < MAXCALLS) )
	{
		i++;
		lold = l;
		x = p + step;
		l = TMath::Abs(langaufun(&x,params) - fy);

		if (l > lold)
			step = -step/10;

		p += step;
	}

	if (i == MAXCALLS)
		return (-2);

	fxr = x;

	// Search for left x location of fy
	p = maxx - 0.5 * params[0];
	step = -params[0];
	lold = -2.0;
	l    = -1e300;
	i    = 0;

	while ( (l != lold) && (i < MAXCALLS) )
	{
		i++;

		lold = l;
		x = p + step;
		l = TMath::Abs(langaufun(&x,params) - fy);
      
		if (l > lold)
			step = -step/10;
      
		p += step;
	}

	if (i == MAXCALLS)
		return (-3);

	fxl = x;
	FWHM = fxr - fxl;
	return (0);
}