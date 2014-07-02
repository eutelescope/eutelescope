/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@cern.ch
 */


// alibava includes ".h"
#include "AlibavaMerger.h"
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
#include <IMPL/TrackerDataImpl.h>
#include <IO/LCWriter.h>
#include <IO/LCReader.h>
#include <Exceptions.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackerHitImpl.h>

// system includes <>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <memory>
#include <glob.h>
#include <vector>
#include <set>
#include <map>

// eutelescope includes ""
#include "anyoption.h"
#include "EUTELESCOPE.h"

// ROOT includes ".h"
#include "TH2D.h"
#include "TF1.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TSystem.h"


using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;
using namespace IMPL;


AlibavaMerger::AlibavaMerger () :
AlibavaBaseProcessor("AlibavaMerger")
{

	_description = "AlibavaMerger merges the Alibava data stream with the telescope data stream.";

	registerProcessorParameter ("AlibavaFile", "The filename where the alibava data is stored", _alibavaFile , string("alibava.slcio"));

	registerProcessorParameter ("TelescopeFile", "The filename where the telescope data is stored", _telescopeFile , string("telescope.slcio"));
	
	registerProcessorParameter ("AlibavaCollectionName", "The name of the alibava collection we want to merge", _alibavaCollectionName , string("original_zsdata"));
	
	registerProcessorParameter ("TelescopeCollectionName", "The name of the telescope collection we want to merge", _telescopeCollectionName , string("original_zsdata"));
	
	registerProcessorParameter ("AlibavaCollectionName2", "The name of the secondary alibava collection we want to merge", _alibavaCollectionName2 , string("clusters"));
	
	registerProcessorParameter ("TelescopeCollectionName2", "The name of the secondary telescope collection we want to merge", _telescopeCollectionName2 , string("cluster_m26"));
	
	registerProcessorParameter ("OutputCollectionName", "The name of the output collection we want to create", _outputCollectionName , string("original_zsdata"));

	registerProcessorParameter ("OutputCollectionName2", "The name of the secondary output collection we want to create", _outputCollectionName2 , string("combinedcluster"));
	
	registerProcessorParameter ("OutputCollectionName3", "The name of the tertiary output collection we want to create", _outputCollectionName3 , string("zsdata_m26"));
	
	registerProcessorParameter ("MergeType", "Do we want to merge hits (0) or clusters (1)? If (1), then the first collection type must be trackerdata, the second trackerpulse!", _mergetype , int(0));
	
	registerProcessorParameter ("OutputMode", "The verbosity of the merged hits: 0 to only write events where we have a hit in both systems, 1 for events where we at least have one alibava hit, 2 for events where we at least have one telescope hit and 3 for writing out all events, even if there are no hits in them", _outputmode, int(0));

	// the unsensitve axis of the strip sensor
	registerProcessorParameter ("UnsensitiveAxis", "The unsensitive axis of our strip sensor", _nonsensitiveaxis, string ("x"));

}


void AlibavaMerger::init () {
	streamlog_out ( MESSAGE4 ) << "Running init" << endl;

	// this method is called only once even when the rewind is active
	// usually a good idea to
	printParameters ();
}

void AlibavaMerger::processRunHeader (LCRunHeader * rdr) {
	streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;

	auto_ptr<AlibavaRunHeaderImpl> arunHeader ( new AlibavaRunHeaderImpl(rdr)) ;
	arunHeader->addProcessor(type());

	// so we only open the telescope file once...
	_telescopeopen = false;

	bookHistos();
}

// the telescope file is read here:
LCEvent *AlibavaMerger::readTelescope ()
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

void AlibavaMerger::addCorrelation(float ali_x, float ali_y, float ali_z, float tele_x, float tele_y, float tele_z, int event)
{
	if ( TH2D * signalHisto = dynamic_cast<TH2D*> (_rootObjectMap["Correlation_X"]) )
	{
		signalHisto->Fill(ali_x,tele_x);
	}
	if ( TH2D * signalHisto2 = dynamic_cast<TH2D*> (_rootObjectMap["Correlation_Y"]) )
	{
		signalHisto2->Fill(ali_y,tele_y);
	}
	if ( TH2D * signalHisto3 = dynamic_cast<TH2D*> (_rootObjectMap["Correlation_Z"]) )
	{
		signalHisto3->Fill(ali_z,tele_z);
	}
	if ( TH2D * signalHisto4 = dynamic_cast<TH2D*> (_rootObjectMap["Correlation_Event"]) )
	{
		if (_nonsensitiveaxis == "x")
		{
			signalHisto4->Fill(event,(tele_y/1152.0 - ali_y/256.0));
		}
		if (_nonsensitiveaxis == "y")
		{
			signalHisto4->Fill(event,(tele_x/1152.0 - ali_x/256.0));
		}
	}
}


void AlibavaMerger::processEvent (LCEvent * anEvent) {

	// the alibava event we are processing
	AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);

	// the input collections
	LCCollectionVec * alibavaCollectionVec;
	LCCollectionVec * telescopeCollectionVec;

	// let's merge hits
	if (_mergetype == 0)
	{
		// for hits we only need one output collection
		LCCollectionVec * outputCollectionVec = new LCCollectionVec(LCIO::TRACKERHIT);

		try
		{

			// the correlation output
			float alibava_x = 0.0;
			float alibava_y = 0.0;
			float alibava_z = 0.0;
			float telescope_x = 0.0;
			float telescope_y = 0.0;
			float telescope_z = 0.0;

			// the alibava data is passed from the event
			alibavaCollectionVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( _alibavaCollectionName ) ) ;
			int alibavasize = alibavaCollectionVec->getNumberOfElements();
			streamlog_out ( DEBUG0 ) << alibavasize << " Elements in Alibava event!" << endl;

			// the telescope is read by the function
			LCEvent* evt= readTelescope();
			telescopeCollectionVec = dynamic_cast< LCCollectionVec * > (evt->getCollection(_telescopeCollectionName)) ;
			int telescopesize = telescopeCollectionVec->getNumberOfElements();
			streamlog_out ( DEBUG0 ) << telescopesize << " Elements in Telescope event!" << endl;

			// loop over the alibava event
			for (int i = 0; i<alibavasize; i++)
			{
				streamlog_out ( DEBUG0 ) << "Reading alibava..." << endl;
				lcio::TrackerHitImpl * input  = dynamic_cast< lcio::TrackerHitImpl * > ( alibavaCollectionVec->getElementAt( i ) ) ;
				lcio::TrackerHitImpl * output = new lcio::TrackerHitImpl;
				output->setCellID0 ( input->getCellID0() ) ;
				output->setCellID1 ( input->getCellID1() ) ;
				output->setCovMatrix ( input->getCovMatrix() ) ;
				output->setEDep ( input->getEDep() ) ;
				output->setEDepError( input->getEDepError() ) ;
				output->setPosition ( input->getPosition() ) ;
				output->setQuality ( input->getQuality() ) ;
				output->setTime ( input->getTime() ) ;
				output->setType ( input->getType() ) ;
				outputCollectionVec->addElement(output);
				streamlog_out ( DEBUG0 ) << "Wrote alibava..." << endl;

				// get average position
				alibava_x += input->getPosition()[0];
				alibava_y += input->getPosition()[1];
				alibava_z += input->getPosition()[2];

			}

			// loop over the telescope event
			for (int j = 0; j<telescopesize; j++)
			{
				streamlog_out ( DEBUG0 ) << "Reading telescope..." << endl;
				lcio::TrackerHitImpl * input  = dynamic_cast< lcio::TrackerHitImpl * > ( telescopeCollectionVec->getElementAt( j ) ) ;
				lcio::TrackerHitImpl * output = new lcio::TrackerHitImpl;
				output->setCellID0 ( input->getCellID0() ) ;
				output->setCellID1 ( input->getCellID1() ) ;
				output->setCovMatrix ( input->getCovMatrix() ) ;
				output->setEDep ( input->getEDep() ) ;
				output->setEDepError( input->getEDepError() ) ;
				output->setPosition ( input->getPosition() ) ;
				output->setQuality ( input->getQuality() ) ;
				output->setTime ( input->getTime() ) ;
				output->setType ( input->getType() ) ;
				outputCollectionVec->addElement(output);
				streamlog_out ( DEBUG0 ) << "Wrote telescope..." << endl;

				// get average position
				telescope_x = input->getPosition()[0];
				telescope_y = input->getPosition()[1];
				telescope_z = input->getPosition()[2];

			}
			
			// correlation plot. Only makes sense if we have an event in both systems
			if (alibavasize > 0 && telescopesize > 0)
			{
				addCorrelation((alibava_x/alibavasize),(alibava_y/alibavasize),(alibava_x/alibavasize),(telescope_x/telescopesize),(telescope_y/telescopesize),(telescope_z/telescopesize),anEvent->getEventNumber());

				streamlog_out( DEBUG0 ) << "Filling histogram with: Ali_X: " << (alibava_x/alibavasize) << " Ali_Y: " << (alibava_y/alibavasize) << " Ali_Z: " << (alibava_z/alibavasize) << " Tele_X: " << (telescope_x/telescopesize) << " Tele_Y: " << (telescope_y/telescopesize) << " Tele_Z: " << (telescope_z/telescopesize) << endl;
			}

			// Now the output. Usually the telescope will have more events...
			if (_outputmode == 0)
			{
				if (alibavasize > 0 && telescopesize > 0)
				{
					streamlog_out ( DEBUG1 ) << "Mode 0. Writing out event " << anEvent->getEventNumber() << endl;
					alibavaEvent->addCollection(outputCollectionVec,_outputCollectionName);
				}
			} if (_outputmode == 1) {
				if (alibavasize > 0 )
				{
					streamlog_out ( DEBUG1 ) << "Mode 1. Writing out event " << anEvent->getEventNumber() << endl;
					alibavaEvent->addCollection(outputCollectionVec,_outputCollectionName);
				}
			} if (_outputmode == 2) {
				if (telescopesize > 0 )
				{
					streamlog_out ( DEBUG1 ) << "Mode 2. Writing out event " << anEvent->getEventNumber() << endl;
					alibavaEvent->addCollection(outputCollectionVec,_outputCollectionName);
				}
			} if (_outputmode == 3) {
				streamlog_out ( DEBUG1 ) << "Mode 3. Writing out event " << anEvent->getEventNumber() << endl;
				alibavaEvent->addCollection(outputCollectionVec,_outputCollectionName);
			}

		}
		catch ( lcio::DataNotAvailableException )
		{
			streamlog_out( DEBUG5 ) << "Collections ("<<_alibavaCollectionName<<") or ("<<_telescopeCollectionName<<") not found! " << endl;
		}
	
	}
	
	// let's merge clusters
	if (_mergetype == 1)
	{

		// the correlation output
		float alibava_x = 0.0;
		float alibava_y = 0.0;
		float alibava_z = 0.0;
		float telescope_x = 0.0;
		float telescope_y = 0.0;
		float telescope_z = 0.0;

		// here we have to merge two collections: the trackerdata and the trackerpulse. We also keep the zsdata for reference -> 3 collections in total

		// out output collections
		LCCollectionVec * outputCollectionVec = new LCCollectionVec(LCIO::TRACKERDATA);
		LCCollectionVec * outputCollectionVec2 = new LCCollectionVec(LCIO::TRACKERPULSE);
		LCCollectionVec * outputCollectionVec3 = new LCCollectionVec(LCIO::TRACKERDATA);

		// the input collections
		LCCollectionVec * alibavaCollectionVec2;
		LCCollectionVec * telescopeCollectionVec2;
		LCCollectionVec * telescopeCollectionVec3;
		
		// these encoders produce output - the encoding string should be the same as in the clustering and fit to the telescope
		CellIDEncoder< TrackerPulseImpl > zsDataEncoder ( "sensorID:5,clusterID:12,xSeed:12,ySeed:12,xCluSize:5,yCluSize:5,type:5", outputCollectionVec2 );	
		CellIDEncoder<TrackerDataImpl> idClusterEncoder( "sensorID:5,clusterID:12,sparsePixelType:5,quality:5", outputCollectionVec  );
		CellIDEncoder<TrackerDataImpl> idClusterEncoder3( "sensorID:5,clusterID:12,sparsePixelType:5,quality:5", outputCollectionVec3  );
		CellIDEncoder< TrackerPulseImpl > zsDataEncoder2 ( "sensorID:5,clusterID:12,xSeed:12,ySeed:12,xCluSize:5,yCluSize:5,type:5", outputCollectionVec2 );	
		CellIDEncoder<TrackerDataImpl> idClusterEncoder2( "sensorID:5,clusterID:12,sparsePixelType:5,quality:5", outputCollectionVec  );

		try
		{
			// the alibava data is passed from the event
			alibavaCollectionVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( _alibavaCollectionName ) ) ;
			int alibavasize = alibavaCollectionVec->getNumberOfElements();
			streamlog_out ( DEBUG0 ) << alibavasize << " Elements in Alibava event!" << endl;

			// so is the secondary collection of the same size (we hope)
			alibavaCollectionVec2 = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( _alibavaCollectionName2 ) ) ;

			// the telescope is read by the function
			LCEvent* evt= readTelescope();
			telescopeCollectionVec = dynamic_cast< LCCollectionVec * > ( evt->getCollection( _telescopeCollectionName ) ) ;
			int telescopesize = telescopeCollectionVec->getNumberOfElements();
			streamlog_out ( DEBUG0 ) << telescopesize << " Elements in Telescope event!" << endl;
			
			// and the secondary collections
			telescopeCollectionVec2 = dynamic_cast< LCCollectionVec * > ( evt->getCollection( _telescopeCollectionName2 ) ) ;
			int telescopesize2 = telescopeCollectionVec2->getNumberOfElements();
			streamlog_out ( DEBUG0 ) << telescopesize2 << " Elements in Telescope event - Collection 2!" << endl;
			
			// this guy can have less events and is not merged, but copied
			telescopeCollectionVec3 = dynamic_cast< LCCollectionVec * > ( evt->getCollection( _outputCollectionName3 ) ) ;
			int telescopesize3 = telescopeCollectionVec3->getNumberOfElements();
			streamlog_out ( DEBUG0 ) << telescopesize3 << " Elements in Telescope event - Collection 3!" << endl;

			// this decoder produces input
			CellIDDecoder<TrackerPulseImpl> inputDecoder(telescopeCollectionVec2);
			CellIDDecoder<TrackerPulseImpl> inputDecoder2(alibavaCollectionVec2);

			// first, loop over the alibava event
			for (int i = 0; i<alibavasize; i++)
			{
				streamlog_out ( DEBUG0 ) << "Reading alibava..." << endl;
				lcio::TrackerDataImpl * input  = dynamic_cast< lcio::TrackerDataImpl * > ( alibavaCollectionVec->getElementAt( i ) ) ;
				lcio::TrackerDataImpl * output = new lcio::TrackerDataImpl;
				output->setChargeValues ( input->getChargeValues() ) ;
				output->setCellID0 ( input->getCellID0() ) ;
				output->setCellID1 ( input->getCellID1() ) ;
				output->setTime ( input->getTime() ) ;
				outputCollectionVec->addElement(output);

				lcio::TrackerPulseImpl * input2  = dynamic_cast< lcio::TrackerPulseImpl * > ( alibavaCollectionVec2->getElementAt( i ) ) ;
				lcio::TrackerPulseImpl * output2 = new lcio::TrackerPulseImpl;
				output2->setCellID0 ( input2->getCellID0() ) ;
				output2->setCellID1 ( input2->getCellID1() ) ;
				output2->setTime ( input2->getTime() ) ;
				output2->setCharge ( input2->getCharge() ) ;
				output2->setQuality ( input2->getQuality() ) ;
				output2->setTrackerData ( input2->getTrackerData() ) ;
				outputCollectionVec2->addElement(output2);

				streamlog_out ( DEBUG0 ) << "Wrote alibava..." << endl;
				
				// output correlation into a histogram
				alibava_x += inputDecoder(input2)["xSeed"];
				alibava_y += inputDecoder(input2)["ySeed"];
				alibava_z += 0.0;
				
			}

			// loop over the telescope cluster event
			for (int j = 0; j<telescopesize; j++)
			{ 
				streamlog_out ( DEBUG0 ) << "Reading telescope..." << endl;
				lcio::TrackerDataImpl * input  = dynamic_cast< lcio::TrackerDataImpl * > ( telescopeCollectionVec->getElementAt( j ) ) ;
				lcio::TrackerPulseImpl * input2  = dynamic_cast< lcio::TrackerPulseImpl * > ( telescopeCollectionVec2->getElementAt( j ) ) ;
				lcio::TrackerPulseImpl * output2 = new lcio::TrackerPulseImpl;
				lcio::TrackerDataImpl * clusterFrame = new lcio::TrackerDataImpl();

				// decode our values into these ints ...
				int sensorID = inputDecoder(input2)["sensorID"];
				int clusterID = inputDecoder(input2)["clusterID"];
				int xSeed = inputDecoder(input2)["xSeed"];
				int ySeed = inputDecoder(input2)["ySeed"];
				int xCluSize = inputDecoder(input2)["xCluSize"];
				int yCluSize = inputDecoder(input2)["yCluSize"];
				int type = inputDecoder(input2)["type"];

				// ... and encode them out again
				zsDataEncoder["sensorID"] = sensorID;
				zsDataEncoder["clusterID"] = clusterID;
				zsDataEncoder["xSeed"] = xSeed;
				zsDataEncoder["ySeed"] = ySeed;
				zsDataEncoder["xCluSize"] = xCluSize;
				zsDataEncoder["yCluSize"] = yCluSize;
				zsDataEncoder["type"] = type;
				zsDataEncoder.setCellID(output2);
				output2->setTrackerData(clusterFrame);
				idClusterEncoder["sensorID"] = sensorID;
				idClusterEncoder["clusterID"] = clusterID;

				// set this to 1 for hitmaker
				idClusterEncoder["sparsePixelType"] = static_cast<int>(1);
				//idClusterEncoder["type"] = static_cast<int>(kEUTelSparseClusterImpl);
				idClusterEncoder["quality"] = static_cast<int>(0);
				clusterFrame->setChargeValues ( input->getChargeValues() ) ;
				idClusterEncoder.setCellID(clusterFrame);
				outputCollectionVec->push_back(clusterFrame);

				outputCollectionVec2->push_back(output2);

				streamlog_out ( DEBUG0 ) << "Wrote telescope..." << endl;
				
				// output correlation into a histogram
				telescope_x += xSeed;
				telescope_y += ySeed;
				telescope_z += sensorID;
			}

			// one more loop over the telescope, now the third collection
			for (int k = 0; k<telescopesize3; k++)
			{
				streamlog_out ( DEBUG0 ) << "Reading telescope again..." << endl;
				lcio::TrackerDataImpl * input  = dynamic_cast< lcio::TrackerDataImpl * > ( telescopeCollectionVec3->getElementAt( k ) ) ;
				lcio::TrackerDataImpl * output = new lcio::TrackerDataImpl;
				output->setChargeValues ( input->getChargeValues() ) ;
				output->setCellID0 ( input->getCellID0() ) ;
				output->setCellID1 ( input->getCellID1() ) ;
				output->setTime ( input->getTime() ) ;
				outputCollectionVec3->addElement(output);
				streamlog_out ( DEBUG0 ) << "Wrote telescope again..." << endl;
			}
			
			/*
			// loop over the alibava event and change the missing coordinate to the average telescope position
			for (int i = 0; i<alibavasize; i++)
			{
				streamlog_out ( DEBUG0 ) << "Reading alibava..." << endl;
				lcio::TrackerDataImpl * input  = dynamic_cast< lcio::TrackerDataImpl * > ( alibavaCollectionVec->getElementAt( i ) ) ;
				lcio::TrackerPulseImpl * input2  = dynamic_cast< lcio::TrackerPulseImpl * > ( alibavaCollectionVec2->getElementAt( i ) ) ;
				//lcio::TrackerDataImpl * output = new lcio::TrackerDataImpl;
				lcio::TrackerPulseImpl * output3 = new lcio::TrackerPulseImpl;
				lcio::TrackerDataImpl * clusterFrame2 = new lcio::TrackerDataImpl();



				// decode our values into these ints ...
				int sensorID = inputDecoder2(input2)["sensorID"];
				int clusterID = inputDecoder2(input2)["clusterID"];
				int xSeed = inputDecoder2(input2)["xSeed"];
				int ySeed = inputDecoder2(input2)["ySeed"];
				int xCluSize = inputDecoder2(input2)["xCluSize"];
				int yCluSize = inputDecoder2(input2)["yCluSize"];
				int type = inputDecoder2(input2)["type"];

				// ... and encode them out again
				zsDataEncoder2["sensorID"] = sensorID;
				zsDataEncoder2["clusterID"] = clusterID;
				
				if (_nonsensitiveaxis == "x" )
				{
					if (telescopesize > 0)
					{
						zsDataEncoder2["xSeed"] = telescope_x/telescopesize;
						streamlog_out (DEBUG2) << "Setting average alibava x to " << telescope_x/telescopesize << endl;
					}
					else
					{
						zsDataEncoder2["xSeed"] = xSeed;
						streamlog_out (DEBUG2) << "No telescope event!" << endl;
					}
					zsDataEncoder2["ySeed"] = ySeed;
				}
				if (_nonsensitiveaxis == "y" )
				{
					zsDataEncoder2["xSeed"] = xSeed;
					if (telescopesize > 0)
					{
						zsDataEncoder2["ySeed"] = telescope_y/telescopesize;
						streamlog_out (DEBUG2) << "Setting average alibava y to " << telescope_y/telescopesize << endl;
					}
					else
					{
						zsDataEncoder2["ySeed"] = ySeed;
						streamlog_out (DEBUG2) << "No telescope event!" << endl;
					}
				}
				zsDataEncoder2["xCluSize"] = xCluSize;
				zsDataEncoder2["yCluSize"] = yCluSize;
				zsDataEncoder2["type"] = type;
				zsDataEncoder2.setCellID(output3);
				output3->setTrackerData(clusterFrame2);
				idClusterEncoder2["sensorID"] = sensorID;
				idClusterEncoder2["clusterID"] = clusterID;

				// set this to 1 for hitmaker
				idClusterEncoder2["sparsePixelType"] = static_cast<int>(1);
				//idClusterEncoder["type"] = static_cast<int>(kEUTelSparseClusterImpl);
				idClusterEncoder2["quality"] = static_cast<int>(0);
				clusterFrame2->setChargeValues ( input->getChargeValues() ) ;
				idClusterEncoder2.setCellID(clusterFrame2);
				outputCollectionVec->push_back(clusterFrame2);

				outputCollectionVec2->push_back(output3);

				streamlog_out ( DEBUG0 ) << "Wrote alibava..." << endl;
				
				// output correlation into a histogram
				alibava_x += inputDecoder2(input2)["xSeed"];
				alibava_y += inputDecoder2(input2)["ySeed"];
				alibava_z += 0.0;
				
			}
			
			*/
			// correlation plot. Only makes sense if we have an event in both systems
			if (alibavasize > 0 && telescopesize > 0)
			{
				addCorrelation((alibava_x/alibavasize),(alibava_y/alibavasize),(alibava_z/alibavasize),(telescope_x/telescopesize),(telescope_y/telescopesize),(telescope_z/telescopesize),anEvent->getEventNumber());
			
				streamlog_out( DEBUG0 ) << "Filling histogram with: Ali_X: " << (alibava_x/alibavasize) << " Ali_Y: " << (alibava_y/alibavasize) << " Ali_Z: " << (alibava_z/alibavasize) << " Tele_X: " << (telescope_x/telescopesize) << " Tele_Y: " << (telescope_y/telescopesize) << " Tele_Z: " << (telescope_z/telescopesize) << endl;
			}

			// Now the output. Usually the telescope will have more events...
			if (_outputmode == 0)
			{
				if (alibavasize > 0 && telescopesize > 0)
				{
					streamlog_out ( DEBUG1 ) << "Mode 0. Writing out event " << anEvent->getEventNumber() << endl;
					alibavaEvent->addCollection(outputCollectionVec,_outputCollectionName);
					alibavaEvent->addCollection(outputCollectionVec2,_outputCollectionName2);
					alibavaEvent->addCollection(outputCollectionVec3,_outputCollectionName3);
				}
			} if (_outputmode == 1) {
				if (alibavasize > 0 )
				{
					streamlog_out ( DEBUG1 ) << "Mode 1. Writing out event " << anEvent->getEventNumber() << endl;
					alibavaEvent->addCollection(outputCollectionVec,_outputCollectionName);
					alibavaEvent->addCollection(outputCollectionVec2,_outputCollectionName2);
					alibavaEvent->addCollection(outputCollectionVec3,_outputCollectionName3);
				}
			} if (_outputmode == 2) {
				if (telescopesize > 0 )
				{
					streamlog_out ( DEBUG1 ) << "Mode 2. Writing out event " << anEvent->getEventNumber() << endl;
					alibavaEvent->addCollection(outputCollectionVec,_outputCollectionName);
					alibavaEvent->addCollection(outputCollectionVec2,_outputCollectionName2);
					alibavaEvent->addCollection(outputCollectionVec3,_outputCollectionName3);
				}
			} if (_outputmode == 3) {
				streamlog_out ( DEBUG1 ) << "Mode 3. Writing out event " << anEvent->getEventNumber() << endl;
				alibavaEvent->addCollection(outputCollectionVec,_outputCollectionName);
				alibavaEvent->addCollection(outputCollectionVec2,_outputCollectionName2);
				alibavaEvent->addCollection(outputCollectionVec3,_outputCollectionName3);
			}

		}
		catch ( lcio::DataNotAvailableException )
		{
			streamlog_out( DEBUG5 ) << "Collections ("<<_alibavaCollectionName<<") or ("<<_telescopeCollectionName<<") not found! " << endl;
		}
	
	}
}


void AlibavaMerger::check (LCEvent * /* evt */ )
{

}


void AlibavaMerger::end()
{
	// the telescope file is still open, we can now close it
	lcReader->close() ;
	delete lcReader ;
	streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
}


void AlibavaMerger::fillHistos()
{

}


void AlibavaMerger::bookHistos()
{
	stringstream tempHistoTitle;
	tempHistoTitle << "Correlation_X" << ";Alibava Average Channel;Telescope Average Pixel";

	TH2D * signalHisto = new TH2D ("Correlation_X","",256,0,255,1152,0,1152);
	_rootObjectMap.insert(make_pair("Correlation_X", signalHisto));
	string tmp_string = tempHistoTitle.str();
	signalHisto->SetTitle(tmp_string.c_str());


	stringstream tempHistoTitle2;
	tempHistoTitle2 << "Correlation_Y" << ";Alibava Average Channel;Telescope Average Pixel";

	TH2D * signalHisto2 = new TH2D ("Correlation_Y","",256,0,255,1152,0,1152);
	_rootObjectMap.insert(make_pair("Correlation_Y", signalHisto2));
	string tmp_string2 = tempHistoTitle2.str();
	signalHisto2->SetTitle(tmp_string2.c_str());


	stringstream tempHistoTitle3;
	tempHistoTitle3 << "Correlation_Z" << ";Alibava Average Channel;Telescope Average Pixel";

	TH2D * signalHisto3 = new TH2D ("Correlation_Z","",256,0,255,1152,0,1152);
	_rootObjectMap.insert(make_pair("Correlation_Z", signalHisto3));
	string tmp_string3 = tempHistoTitle3.str();
	signalHisto3->SetTitle(tmp_string3.c_str());
	
	
	stringstream tempHistoTitle4;
	tempHistoTitle4 << "Correlation_Event" << ";Event Nr.;Tele_Avg_XorY - Alib_Avg_XorY";

	TH2D * signalHisto4 = new TH2D ("Correlation_Event","",1000,0,500000,100,-1,1);
	_rootObjectMap.insert(make_pair("Correlation_Event", signalHisto4));
	string tmp_string4 = tempHistoTitle4.str();
	signalHisto4->SetTitle(tmp_string4.c_str());
	
	
	streamlog_out ( MESSAGE4 )  << "End of Booking histograms. " << endl;
}


