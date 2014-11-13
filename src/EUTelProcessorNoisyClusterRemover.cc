/*
 *   This processor removes noisy clusters, i.e. clusters which have
 *   been masked as noisy by the EUTelProcessorNoisyClusterMasker
 *
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelProcessorNoisyClusterRemover.h"
#include "EUTELESCOPE.h"
#include "EUTelTrackerDataInterfacerImpl.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "EUTelRunHeaderImpl.h"

// lcio includes <.h>
#include <LCIOTypes.h>
#include <IO/LCWriter.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerPulseImpl.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>

// system includes
#include <memory>
#include <algorithm>

using namespace std;
using namespace marlin;
using namespace eutelescope;


EUTelProcessorNoisyClusterRemover::EUTelProcessorNoisyClusterRemover():
  Processor("EUTelProcessorNoisyClusterRemover"),
  _inputCollectionName(""),
  _outputCollectionName(""),
  _iRun(0),
  _iEvt(0)
{
  _description ="EUTelProcessorNoisyClusterRemover removes masked noisy clusters (TrackerPulses) from a collection, please note that the clusters have to be masked previously. This processor does not read in a hot pixel collection, it simply removes previusly masked TrackerPulses.";

  registerInputCollection (LCIO::TRACKERPULSE, "InputCollectionName", "Input collection containing noise masked tracker pulse objects", _inputCollectionName, string ("noisy_cluster") );

  registerOutputCollection(LCIO::TRACKERPULSE, "OutputCollectionName", "Output collection where noisy clusters have been removed", _outputCollectionName, string("noisefree_clusters"));

}

void EUTelProcessorNoisyClusterRemover::init () 
{
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters();
  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;
}

void EUTelProcessorNoisyClusterRemover::processRunHeader(LCRunHeader* rdr){

  auto_ptr<EUTelRunHeaderImpl> runHeader ( new EUTelRunHeaderImpl(rdr) );
  runHeader->addProcessor(type()) ;
  // increment the run counter
  ++_iRun;
  // reset the event counter
  _iEvt = 0;
}

void EUTelProcessorNoisyClusterRemover::processEvent(LCEvent * event) 
{
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
	
	//now prepare output collection
	LCCollectionVec* outputCollection;
	bool outputCollectionExists = false;
  	_initialOutputCollectionSize = 0;

  	try 
  	{
   		outputCollection = dynamic_cast< LCCollectionVec* > ( event->getCollection( _outputCollectionName ) );
    		outputCollectionExists = true;
    		_initialOutputCollectionSize = outputCollection->size();
  	} 
  	catch ( lcio::DataNotAvailableException& e ) 
  	{
    		outputCollection = new LCCollectionVec(LCIO::TRACKERPULSE);
  	}

        //read the encoding string from the input collection
	std::string encodingString = pulseInputCollectionVec->getParameters().getStringVal( LCIO::CellIDEncoding );	
	//and the encoder for the output data
	CellIDEncoder<TrackerPulseImpl> idZSGenDataEncoder( encodingString , outputCollection);
	
	//loop over all pulses
	for( size_t iPulse = 0 ; iPulse < pulseInputCollectionVec->size(); iPulse++ ) 
	{
		//get the tracker pulse
        	TrackerPulseImpl* inputPulse = dynamic_cast<TrackerPulseImpl*> ( pulseInputCollectionVec->getElementAt(iPulse) );
		//and its quality
		int quality = cellDecoder(inputPulse)["quality"];
		int sensorID  = cellDecoder(inputPulse)["sensorID"];
		
		//if the kNoisyCluster flag is NOT set, we add the pulse to the output collection
		if(!(quality & kNoisyCluster))
		{
			//TrackerPulseImpl for the output collection
			auto_ptr<TrackerPulseImpl> outputPulse (new TrackerPulseImpl );

			//copy the information which is the same
			outputPulse->setCellID0( inputPulse->getCellID0() );
			outputPulse->setCellID1( inputPulse->getCellID1() );
			outputPulse->setTime( inputPulse->getTime() );
			outputPulse->setCharge( inputPulse->getCharge() );
			outputPulse->setCovMatrix( inputPulse->getCovMatrix() );
			outputPulse->setQuality( inputPulse->getQuality() );
			outputPulse->setTrackerData( inputPulse->getTrackerData() );
			
	                outputCollection->push_back( outputPulse.release() );
                        streamlog_out ( DEBUG4 ) << "elements in collection : " << outputCollection->size() << std::endl;
		}
		//if the cluster is noisy, thus removed, we count that for nice user output
		else
		{
			_removedNoisyPulses[sensorID]++;
		}
        }// loop over all pulses

	//add the collection if necessary
	if ( !outputCollectionExists && ( outputCollection->size() != _initialOutputCollectionSize )) 
	{
		event->addCollection( outputCollection, _outputCollectionName );
		streamlog_out ( DEBUG4 ) << "adding output collection: " << _outputCollectionName << " " <<  outputCollection->size() << std::endl;
	}

	if ( !outputCollectionExists && ( outputCollection->size() == _initialOutputCollectionSize ) ) 
	{
		delete outputCollection;
	}	
//rest of memory cleaned up by auto_ptrs
}

void EUTelProcessorNoisyClusterRemover::end() 
{
	//Print out some stats for the user
	streamlog_out ( MESSAGE4 ) << "Noisy cluster remover successfully finished" << endl;
	streamlog_out ( MESSAGE4 ) << "Printing summary:" << endl;
	for(std::map<int,int>::iterator it = _removedNoisyPulses.begin(); it != _removedNoisyPulses.end(); ++it)
	{
		streamlog_out ( MESSAGE4 ) << "Removed " << (*it).second << " noisy pulses from plane " << (*it).first << "." << endl;
	}
	
}

