// eutelescope includes ".h"
#include "EUTelProcessorTransformFromGEAR.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"

// lcio includes <.h>
#include <EVENT/LCCollection.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

// system includes <>
#include <memory>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;
using namespace gear;

EUTelProcessorTransformFromGEAR::EUTelProcessorTransformFromGEAR () :Processor("EUTelProcessorTransformFromGEAR")
{
  _description = "Apply alignment constants from GEAR-File to hit collection";

  registerInputCollection(LCIO::TRACKERHIT, "InputHitCollectionName","The name of the input hit collection in the local frame", _inputHitCollectionName, std::string("inhit"));
  
  registerOutputCollection(LCIO::TRACKERHIT, "OutputHitCollectionName","The name of the output hit collection in the global frame", _outputHitCollectionName, std::string("outhit"));
}


void EUTelProcessorTransformFromGEAR::init()
{
	//for info, good idea to:
	printParameters();

	// set to zero the run and event counters
	_iRun = 0;  
	_iEvt = 0;
	std::vector<int> sensorVec = geo::gGeometry().sensorIDsVec();

	for( std::vector<int>::iterator it = sensorVec.begin(); it != sensorVec.end(); it++ ) 
	{
		int sensorID = *it;
		
		Eigen::Vector3d offVec = geo::gGeometry().getOffsetVector(sensorID);
		Eigen::Matrix3d rotMat = geo::gGeometry().rotationMatrixFromAngles(sensorID);
		Eigen::Matrix3d flipMat = geo::gGeometry().getFlipMatrix(sensorID);

		_rotMat[sensorID] = rotMat*flipMat;
		_offVec[sensorID] = offVec;

		std::cout << "Offset Vector for plane: " << sensorID << ": " << std::endl << _offVec[sensorID] << std::endl;
		std::cout << "Rot*FlipMat for plane: " << sensorID << ": " << std::endl << _rotMat[sensorID] << std::endl;
	}

}

void EUTelProcessorTransformFromGEAR::processRunHeader(LCRunHeader* rdr)
{
	std::auto_ptr<EUTelRunHeaderImpl>runHeader(new EUTelRunHeaderImpl( rdr ));
	runHeader->addProcessor(type());
	++_iRun;
}

void EUTelProcessorTransformFromGEAR::processEvent(LCEvent* event)
{
	_iEvt++;
 	LCCollectionVec* inputCollectionVec = NULL;
	
	try
	{
		inputCollectionVec = dynamic_cast<LCCollectionVec*>(event->getCollection(_inputHitCollectionName));

  	}

	catch(DataNotAvailableException& e)
	{
	       return;	
	}
	
	UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder( inputCollectionVec );
	
	//now prepare output collection
	LCCollectionVec* outputCollectionVec;
	bool outputCollectionVecExists = false;
  	_initialOutputCollectionSize = 0;

  	try 
  	{
   		outputCollectionVec = dynamic_cast< LCCollectionVec* > ( event->getCollection( _outputHitCollectionName ) );
    		outputCollectionVecExists = true;
    		_initialOutputCollectionSize = outputCollectionVec->size();
  	} 
  	catch ( lcio::DataNotAvailableException& e ) 
  	{
    		outputCollectionVec = new LCCollectionVec(LCIO::TRACKERHIT);
  	}

        //read the encoding string from the input collection
	std::string encodingString = inputCollectionVec->getParameters().getStringVal( LCIO::CellIDEncoding );	
	//and the encoder for the output data
	CellIDEncoder<TrackerHitImpl> idZSGenDataEncoder( encodingString , outputCollectionVec);


	//Loop over all hits and determine the hit on the fixed plane:
	for(size_t hitNo = 0; hitNo < inputCollectionVec->size(); hitNo++)
	{
		TrackerHitImpl* inputHit = dynamic_cast<TrackerHitImpl*>(inputCollectionVec->getElementAt(hitNo));
		
		int sensorID = hitDecoder(inputHit)["sensorID"];
		const double* hitPos = inputHit->getPosition();
		
		Eigen::Vector3d oldPos;
		oldPos << hitPos[0], hitPos[1], hitPos[2];

		Eigen::Vector3d newPos = _rotMat[sensorID]*oldPos+_offVec[sensorID];

		double newPosArray[3];

		newPosArray[0] = newPos(0);
		newPosArray[1] = newPos(1);
		newPosArray[2] = newPos(2); 

		//TrackerHitImpl for the output collection
		std::auto_ptr<TrackerHitImpl> outputHit ( new TrackerHitImpl );

		//copy the information which is the same
		outputHit->setCellID0( inputHit->getCellID0() );
		outputHit->setCellID1( inputHit->getCellID1() );
		outputHit->setTime( inputHit->getTime() );
		outputHit->setCovMatrix( inputHit->getCovMatrix() );
		outputHit->setQuality( inputHit->getQuality() );
		outputHit->rawHits() =  inputHit->getRawHits();
		outputHit->setPosition( &newPosArray[0] );

		outputCollectionVec->push_back( outputHit.release() );
	}


	//add the collection if necessary
	if ( !outputCollectionVecExists && ( outputCollectionVec->size() != _initialOutputCollectionSize )) 
	{
		event->addCollection( outputCollectionVec, _outputHitCollectionName );
	}

	if ( !outputCollectionVecExists && ( outputCollectionVec->size() == _initialOutputCollectionSize ) ) 
	{
		delete outputCollectionVec;
	}	
	//rest of memory cleaned up by auto_ptrs
}

void EUTelProcessorTransformFromGEAR::end()
{
	//NOP NOP NOP
}
