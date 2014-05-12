// Author Tobias Bisanz,  <tobias.bisanz@phys.uni-goettingen.de>
// Version $Id$
/*
 *   This processor removes hot pixels as specified by the hot pixel
 *   collection and removes them from any given other tracker collection.
 *   It is only campatible with the EUTelGenerisSparsePixel type, the
 *   hot pixel as well as the input collections must use it.
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
#include "EUTelRunHeaderImpl.h"
#include "EUTelTrackerDataInterfacer.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <LCIOTypes.h>
#include <IO/LCWriter.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <Exceptions.h>

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
  _description ="EUTelProcessorNoisyClusterRemover removes hot pixels from a tracker data collection, it reads in a hot pixel collection and an input collection. If any of the hits in the input collection is a hot pixel it gets removed.";

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
        	TrackerPulseImpl* inputPulse = dynamic_cast<TrackerPulseImpl*> ( pulseInputCollectionVec->getElementAt( iPulse ) );
	
		//get the quality
		int quality = cellDecoder(inputPulse)["quality"];
		
		if(!(quality & kNoisyCluster))
		{
			
			//TrackerPulseImpl for the output collection
			auto_ptr<TrackerPulseImpl> outputPulse ( new TrackerPulseImpl );

			//copy the information which is the same
			outputPulse->setCellID0( inputPulse->getCellID0() );
			outputPulse->setCellID1( inputPulse->getCellID1() );
			outputPulse->setTime( inputPulse->getTime() );
			outputPulse->setCharge( inputPulse->getCharge() );
			outputPulse->setCovMatrix( inputPulse->getCovMatrix() );
			outputPulse->setQuality( inputPulse->getQuality() );
			outputPulse->setTrackerData( inputPulse->getTrackerData() );
			
			outputCollection->push_back( outputPulse.release() );
		}		
        }// loop over all pulses

	//add the collection if necessary
	if ( !outputCollectionExists && ( outputCollection->size() != _initialOutputCollectionSize )) 
	{
		event->addCollection( outputCollection, _outputCollectionName );
	}

	if ( !outputCollectionExists && ( outputCollection->size() == _initialOutputCollectionSize ) ) 
	{
		delete outputCollection;
	}	
//rest of memory cleaned up by auto_ptrs
}

void EUTelProcessorNoisyClusterRemover::end() 
{
	//maybe print some info for the user
}

