// eutelescope inlcudes
#include "EUTelAPIXTbTrackTuple.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelTrackerDataInterfacerImpl.h"
#include "EUTelExternalTrigger.h"
#include "EUTelTrackerDataTriggerInterfacer.h"

// eutelescope geometry
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelGenericPixGeoDescr.h"

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <UTIL/CellIDDecoder.h>

#include <algorithm>

using namespace eutelescope;

EUTelAPIXTbTrackTuple::EUTelAPIXTbTrackTuple()
: Processor("EUTelAPIXTbTrackTuple"),
  _inputTrackColName (""),
  _inputTrackerHitColName (""),
  _inputTelPulseCollectionName(""),
  _inputDutPulseCollectionName(""),
  _telZsColName(""),
  _dutZsColName(""),
  _path2file(""),
  _DUTIDs(std::vector<int>()),
  _nRun (0),
  _nEvt (0),
  _runNr(0),
  _evtNr(0),
  _isFirstEvent(false),
  _file(NULL),
  _eutracks(NULL),
  _nTrackParams(0),
  _xPos(NULL),
  _yPos(NULL),
  _dxdz(NULL),
  _dydz(NULL),
  _trackIden(NULL),
  _trackNum(NULL),
  _chi2(NULL),
  _ndof(NULL),
  _zstree(NULL),
  _nPixHits(0),
  p_col(NULL),
  p_row(NULL),
  p_tot(NULL),
  p_iden(NULL),
  p_lv1(NULL),
  p_hitTime(NULL),
  p_frameTime(NULL),
  p_TLU(NULL),
  _euhits(NULL),
  _nHits(0),
  _hitXPos(NULL),
  _hitYPos(NULL),
  _hitZPos(NULL),
  _hitSensorId(NULL),
  _triggers(NULL),
  _nTriggers(0),
  _nTLUTriggers(0),
  _nExtTriggers(0),
  _TLUTrigTime(NULL),
  _ExtTrigTime(NULL),
  _tots(NULL),
  _nToTs(0),
  _ToTTime(NULL),
  _ToTLength(NULL)
 {
  //processor description
  _description = "Prepare tbtrack style n-tuple with track fit results" ;


  registerInputCollection(LCIO::TRACK, "InputTrackCollectionName", "Name of the input Track collection",
		  		_inputTrackColName, std::string("fittracks"));
  	   
  registerInputCollection( LCIO::TRACKERHIT, "InputTrackerHitCollectionName", "Name of the plane-wide hit-data hit collection"  ,
		_inputTrackerHitColName, std::string("fitpoints") );
 
  registerProcessorParameter ("DutZsColName", "DUT zero supressed data colection name",
 		     _dutZsColName, std::string("zsdata_apix"));

  registerProcessorParameter ("TriggerColName", "External trigger data colection name",
 		     _triggerColName, std::string("eudet_triggers"));

  registerProcessorParameter ("ToTColName", "Tot data colection name",
 		     _totColName, std::string("eudet_tots"));
 
  registerProcessorParameter ("OutputPath", "Path/File where root-file should be stored",
			      _path2file, std::string("NTuple.root"));
  
  registerProcessorParameter ("DUTIDs", "Int std::vector containing the IDs of the DUTs",
		  		_DUTIDs, std::vector<int>());

}


void EUTelAPIXTbTrackTuple::init()
{
	// usually a good idea to
	printParameters();

	_isFirstEvent = true; 
 
	_nRun = 0;
	_nEvt = 0;

	// Prepare TTree/TFiles
	prepareTree();

	//init new geometry
	std::string name("test.root");
	geo::gGeometry().initializeTGeoDescription(name,true);


	for(std::vector<int>::iterator it = _DUTIDs.begin(); it != _DUTIDs.end(); it++)
	{
		//Later we need to shift the sensor since in EUTel centre of sensor is 0|0 while in TBmon(II) it is in the lower left corner
		geo::EUTelGenericPixGeoDescr* geoDescr = geo::gGeometry().getPixGeoDescr( *it ) ;
		float xSize,ySize;
		geoDescr->getSensitiveSize(xSize, ySize);

		_xSensSize[*it] = xSize;
		_ySensSize[*it] = ySize;
	}
}

void EUTelAPIXTbTrackTuple::processRunHeader( LCRunHeader* runHeader) 
{
	std::auto_ptr<EUTelRunHeaderImpl> eutelHeader( new EUTelRunHeaderImpl ( runHeader ) );
	eutelHeader->addProcessor( type() );
	_nRun++;
	
	// Decode and print out Run Header information - just a check
	_runNr = runHeader->getRunNumber();
}

void EUTelAPIXTbTrackTuple::processEvent( LCEvent * event )
{
	_nEvt ++;
	_evtNr = event->getEventNumber();
	EUTelEventImpl* euEvent = static_cast<EUTelEventImpl*> ( event );

	if( euEvent->getEventType() == kEORE )
       	{
		streamlog_out( DEBUG5) << "EORE found: nothing else to do." << std::endl;
    		return;
	}

	//Clear all event info containers
	clear();

	//try to read in hits (e.g. fitted hits in local frame)	
  	if( !readHits( _inputTrackerHitColName, event )) 
	{ 
		return; 
	}

	// read in triggers
	if( readTriggers(_triggerColName, event) )
	  {
	    _triggers->Fill();
	  }
  	
	//read in raw data	
	if(!readZsHits( _dutZsColName , event )) 
	{ 
		return; 
	}
  
  	//read in tracks
	if(!readTracks(event))
       	{
		return;
	}
 
        //fill the trees	
	_zstree->Fill();
	_eutracks->Fill();
	_euhits->Fill();

	
	

	// read in tots
	if( readToTs(_totColName, event) )
	  {
	    _tots->Fill();
	  }

	
	_isFirstEvent = false;
}

void EUTelAPIXTbTrackTuple::end()
{
	//write version number
	_versionNo->push_back(1.1);
	_versionTree->Fill();
	//Maybe some stats output?
	_file->Write();
}

//Read in TrackerHit(Impl) to later dump them
bool EUTelAPIXTbTrackTuple::readHits( std::string hitColName, LCEvent* event )
{
	LCCollection* hitCollection = NULL;
  
	try
       	{
		hitCollection = event->getCollection( hitColName ); 
  	} 
	catch(lcio::DataNotAvailableException& e) 
	{
		streamlog_out( DEBUG2 ) << "Hit collection " << hitColName << " not found in event " << event->getEventNumber()  << "!" << std::endl;
    		return false;
  	}
  	
	int nHit = hitCollection->getNumberOfElements();
	_nHits = nHit;
 
  	for(int ihit=0; ihit< nHit ; ihit++)
       	{
    		TrackerHitImpl* meshit = dynamic_cast<TrackerHitImpl*>( hitCollection->getElementAt(ihit) ) ;
    		const double* pos = meshit->getPosition();	
    		LCObjectVec clusterVec = (meshit->getRawHits());

		UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder ( EUTELESCOPE::HITENCODING );
    		int sensorID = hitDecoder(meshit)["sensorID"];

		//Only dump DUT hits
		if( std::find( _DUTIDs.begin(), _DUTIDs.end(), sensorID) == _DUTIDs.end() )
		{
			continue;
		}

    		double x = pos[0];
    		double y = pos[1];
    		double z = pos[2];

	       	//offset by half sensor/sensitive size
			_hitXPos->push_back(x + _xSensSize.at(sensorID)/2.0);
    		_hitYPos->push_back(y + _ySensSize.at(sensorID)/2.0);
    		_hitZPos->push_back(z);
    		_hitSensorId->push_back(sensorID);
	}

	return true;
}

//Read in TrackerHit to later dump
bool EUTelAPIXTbTrackTuple::readTracks(LCEvent* event)
{
	LCCollection* trackCol = NULL;

	try
	{
    		trackCol = event->getCollection( _inputTrackColName ) ;
  	}
	catch(lcio::DataNotAvailableException& e)
	{
		streamlog_out( DEBUG2 ) << "Track collection " << _inputTrackColName << " not found in event " << event->getEventNumber()  << "!" << std::endl;
		return false;
  	}

	// setup cellIdDecoder to decode the hit properties
	UTIL::CellIDDecoder<TrackerHitImpl>  hitCellDecoder( EUTELESCOPE::HITENCODING );

	int nTrackParams=0;
	for(int itrack=0; itrack< trackCol->getNumberOfElements(); itrack++)
       	{
		lcio::Track* fittrack = dynamic_cast<lcio::Track*>( trackCol->getElementAt(itrack) ) ;
		
		std::vector<EVENT::TrackerHit*> trackhits = fittrack->getTrackerHits();
		double chi2 = fittrack->getChi2();
		double ndof = fittrack->getNdf();
		double dxdz = fittrack->getOmega();
		double dydz = fittrack->getPhi();

 		/* Get the (fitted) hits belonging to this track, 
		   they are in global frame when coming from the straight track fitter */
		for(unsigned int ihit=0; ihit< trackhits.size(); ihit++)
	       	{
      			TrackerHitImpl* fittedHit = dynamic_cast<TrackerHitImpl*>( trackhits.at(ihit) );
      			const double* pos = fittedHit->getPosition();
      			if( (hitCellDecoder(fittedHit)["properties"] & kFittedHit) == 0)
		       	{
			       	continue;
		       	}

     			UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder ( EUTELESCOPE::HITENCODING );
      			int sensorID = hitDecoder(fittedHit)["sensorID"];

			//Dump the (fitted) hits for the DUTs
			if( std::find( _DUTIDs.begin(), _DUTIDs.end(), sensorID) == _DUTIDs.end() )
			{
				continue;
			}

      			nTrackParams++;
      			
			/* Transform to local coordinates */
			double pos_loc[3];
			geo::gGeometry().master2Local(sensorID, pos, pos_loc);
			
			double x = pos_loc[0];
			double y = pos_loc[1];
			
      			//double z = pos[2]; //not used!
			
				//eutrack tree
      			_xPos->push_back(x);
      			_yPos->push_back(y);
      			_dxdz->push_back(dxdz);
      			_dydz->push_back(dydz);
      			_trackIden->push_back(sensorID);
      			_trackNum->push_back(itrack);
      			_chi2->push_back(chi2);
      			_ndof->push_back(ndof);
    		}
  	}

	_nTrackParams = nTrackParams;
	return true;
}

//Read in raw (zs) TrackerData(Impl) to later dump
bool EUTelAPIXTbTrackTuple::readZsHits( std::string colName, LCEvent* event)
{
	LCCollectionVec* zsInputCollectionVec = NULL;
 
	try
	{
		zsInputCollectionVec = dynamic_cast<LCCollectionVec*>( event->getCollection(colName) );
  	}
       	catch(DataNotAvailableException& e)
	{
		streamlog_out( DEBUG2 ) << "Raw ZS data collection " << colName << " not found in event " << event->getEventNumber()  << "!" << std::endl;
    		return false;
  	}
	
	UTIL::CellIDDecoder<TrackerDataImpl> cellDecoder( zsInputCollectionVec );
	for( unsigned int plane = 0; plane < zsInputCollectionVec->size(); plane++ ) 
	{
		TrackerDataImpl* zsData = dynamic_cast< TrackerDataImpl * > ( zsInputCollectionVec->getElementAt( plane ) );
		SparsePixelType type = static_cast<SparsePixelType> ( static_cast<int> (cellDecoder( zsData )["sparsePixelType"]) );
    		int sensorID = cellDecoder( zsData )["sensorID"];
		
		std::auto_ptr<EUTelTrackerDataInterfacer> sparseData = std::auto_ptr<EUTelTrackerDataInterfacer>();

		if (type == kEUTelGenericSparsePixel  ) 
		  {
		   sparseData =  std::auto_ptr<EUTelTrackerDataInterfacer>( new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>(zsData) );
		   EUTelGenericSparsePixel apixPixel;
		   
		   for( unsigned int iHit = 0; iHit < sparseData->size(); iHit++ ) 
		     {
		       sparseData->getSparsePixelAt( iHit, &apixPixel);
		       _nPixHits++;
		       p_iden->push_back( sensorID );
		       p_row->push_back( apixPixel.getYCoord() );
		       p_col->push_back( apixPixel.getXCoord() );
		       p_tot->push_back( static_cast< int >(apixPixel.getSignal()) );
		       p_lv1->push_back( static_cast< int >(apixPixel.getTime()) );
		     }
		   
		  }
		else if( type == kEUTelMuPixel )
		  {
		    sparseData =  std::auto_ptr<EUTelTrackerDataInterfacer>( new EUTelTrackerDataInterfacerImpl<EUTelMuPixel>(zsData) );
		    EUTelMuPixel binaryPixel;
		    for( unsigned int iHit = 0; iHit < sparseData->size(); iHit++ ) 
		     {
		       sparseData->getSparsePixelAt( iHit, &binaryPixel);
		       _nPixHits++;
		       p_iden->push_back( sensorID );
		       p_row->push_back( binaryPixel.getYCoord() );
		       p_col->push_back( binaryPixel.getXCoord() );
		       p_hitTime->push_back( binaryPixel.getHitTime() );
		       p_frameTime->push_back( binaryPixel.getFrameTime() );
		       if ( _TLUTrigTime->size() != 0) 
			 p_TLU->push_back( _TLUTrigTime->at(0) );
		       else 
			 p_TLU->push_back( 0 );
		     }
		  }
		else
		  {
		    throw UnknownDataTypeException("Unknown sparsified pixel");
		  }
	}
  	return true;
} 


//Read in external triggers TrackerData(Impl) to later dump
bool EUTelAPIXTbTrackTuple::readTriggers( std::string colName, LCEvent* event)
{
	LCCollectionVec* triggerInputCollectionVec = NULL;
 
	try
	{
		triggerInputCollectionVec = dynamic_cast<LCCollectionVec*>( event->getCollection(colName) );
  	}
       	catch(DataNotAvailableException& e)
	{
		streamlog_out( DEBUG2 ) << "External trigger collection " << colName << " not found in event " << event->getEventNumber()  << "!" << std::endl;
    		return false;
  	}
	
	UTIL::CellIDDecoder<TrackerDataImpl> cellDecoder( triggerInputCollectionVec );
		
	TrackerDataImpl* triggerData = dynamic_cast< TrackerDataImpl * > ( triggerInputCollectionVec->getElementAt( 0 ) );
	std::auto_ptr<EUTelTrackerDataTriggerInterfacer> sparseData = std::auto_ptr<EUTelTrackerDataTriggerInterfacer>( );
	sparseData =  std::auto_ptr<EUTelTrackerDataTriggerInterfacer>( new EUTelTrackerDataTriggerInterfacer(triggerData) );

	EUTelExternalTrigger trigger;
	for( unsigned int iTrigger = 0; iTrigger < sparseData->size(); iTrigger++ ) 
	  {
	    sparseData->getExternalTriggerAt( iTrigger, &trigger);
	    _nTriggers++;
	    if ( trigger.getLabel() == 0x1 ) // TLU trigger
	      {
	      _TLUTrigTime->push_back( trigger.getTimestamp() );
	      _nTLUTriggers++;
	      }
	    else if ( trigger.getLabel() == 0xBA ) {  // coincidence trigger
	      if ( trigger.getTimestamp() != (unsigned int)trigger.getTimestamp() )
		std::cout << std::hex << trigger.getTimestamp() << " " <<   (unsigned int)trigger.getTimestamp() << std::endl;
	      _ExtTrigTime->push_back( (unsigned int)trigger.getTimestamp() );
	      _nExtTriggers++;
	    }
	    else
	      throw UnknownDataTypeException("Unknown trigger label");
	  }

	return true;
} 

// Read in ToTs saved as external trigger TrackerData 
bool EUTelAPIXTbTrackTuple::readToTs( std::string colName, LCEvent* event)
{
	LCCollectionVec* triggerInputCollectionVec = NULL;
 
	try
	{
		triggerInputCollectionVec = dynamic_cast<LCCollectionVec*>( event->getCollection(colName) );
  	}
       	catch(DataNotAvailableException& e)
	{
		streamlog_out( DEBUG2 ) << "Tot collection " << colName << " not found in event " << event->getEventNumber()  << "!" << std::endl;
    		return false;
  	}
	
	UTIL::CellIDDecoder<TrackerDataImpl> cellDecoder( triggerInputCollectionVec );
		
	TrackerDataImpl* triggerData = dynamic_cast< TrackerDataImpl * > ( triggerInputCollectionVec->getElementAt( 0 ) );
	std::auto_ptr<EUTelTrackerDataTriggerInterfacer> sparseData = std::auto_ptr<EUTelTrackerDataTriggerInterfacer>( );
	sparseData =  std::auto_ptr<EUTelTrackerDataTriggerInterfacer>( new EUTelTrackerDataTriggerInterfacer(triggerData) );

	EUTelExternalTrigger trigger;
	for( unsigned int iTrigger = 0; iTrigger < sparseData->size(); iTrigger++ ) 
	  {
	    sparseData->getExternalTriggerAt( iTrigger, &trigger);
	    _nToTs++;
	    if ( trigger.getLabel() == 0x2 ) // ToT information encoded as ExternalTrigger
	      {
		_ToTTime->push_back( (trigger.getTimestamp() >> 8) & 0xFFFFFFFFFFFF );
		_ToTLength->push_back( trigger. getTimestamp() & 0xFF );
	      }
	    else
	      throw UnknownDataTypeException("Unknown trigger / tot label");
	  }

	return true;
} 


void EUTelAPIXTbTrackTuple::clear()
{
	/* Clear zsdata */
	p_col->clear();
	p_row->clear();
	p_tot->clear();
	p_iden->clear();
	p_lv1->clear();
	p_hitTime->clear();
	p_frameTime->clear();
	p_TLU->clear();
	_nPixHits = 0;
	/* Clear hittrack */
	_xPos->clear();
	_yPos->clear();
	_dxdz->clear();
	_dydz->clear();
	_trackNum->clear();
	_trackIden->clear();
	_chi2->clear();
	_ndof->clear();
	//Clear hits
	_hitXPos->clear();
	_hitYPos->clear();
	_hitZPos->clear();
	_hitSensorId->clear();
	// Clear triggers
	_nTriggers = 0;
	_nTLUTriggers = 0;
	_nExtTriggers = 0;
	_TLUTrigTime->clear();
	_ExtTrigTime->clear();
	_ToTTime->clear();
	_ToTLength->clear();
	_nToTs = 0;
}

void EUTelAPIXTbTrackTuple::prepareTree()
{
	_file = new TFile(_path2file.c_str(),"RECREATE");

	_xPos = new std::vector<double>();	     
	_yPos = new std::vector<double>();	   
	_dxdz = new std::vector<double>();	   
	_dydz = new std::vector<double>();	   
	_trackIden  = new std::vector<int>();
	_trackNum = new std::vector<int>();
	_chi2 = new std::vector<double>();	   
	_ndof = new std::vector<double>();    

	p_col = new std::vector<int>();
	p_row = new std::vector<int>();
	p_tot = new std::vector<double>();
	p_iden = new std::vector<int>();
	p_lv1 = new std::vector<int>();
	p_hitTime = new std::vector<int>();
	p_frameTime = new std::vector<double>();
	p_TLU = new std::vector<unsigned int>();

	_hitXPos = new std::vector<double>();
	_hitYPos = new std::vector<double>();
	_hitZPos = new std::vector<double>();
	_hitSensorId  = new std::vector<int>();

	_TLUTrigTime = new std::vector<unsigned int>();
	_ExtTrigTime = new std::vector<unsigned int>();
	_ToTTime = new std::vector<double>();
	_ToTLength = new std::vector<int>();

	_versionNo = new std::vector<double>();
	_versionTree = new TTree("version","version");
	_versionTree->Branch("no", &_versionNo);

	_euhits = new TTree("fitpoints","fitpoints");
	_euhits->SetAutoSave(1000000000);
	
	_euhits->Branch("nHits", &_nHits);
	_euhits->Branch("xPos", &_hitXPos);
	_euhits->Branch("yPos", &_hitYPos);
	_euhits->Branch("zPos", &_hitZPos);
	_euhits->Branch("sensorId", &_hitSensorId);

	_zstree = new TTree("rawdata", "rawdata");
	_zstree->SetAutoSave(1000000000);
	_zstree->Branch("nPixHits", &_nPixHits);
	_zstree->Branch("euEvt",    &_nEvt);
	_zstree->Branch("col",      &p_col);
	_zstree->Branch("row",      &p_row);
	_zstree->Branch("tot",   "std::vector<double>",   &p_tot);
	_zstree->Branch("lv1",      &p_lv1);
	_zstree->Branch("iden",     &p_iden);
	_zstree->Branch("hitTime",  &p_hitTime);
	_zstree->Branch("frameTime",&p_frameTime);
	_zstree->Branch("TLU",&p_TLU);
	
	//Tree for storing all track param info
	_eutracks = new TTree("tracks", "tracks");
	_eutracks->SetAutoSave(1000000000);
	_eutracks->Branch("nTrackParams", &_nTrackParams);
	_eutracks->Branch("euEvt", &_nEvt);
	_eutracks->Branch("xPos", &_xPos);
	_eutracks->Branch("yPos", &_yPos);
	_eutracks->Branch("dxdz", &_dxdz);
	_eutracks->Branch("dydz", &_dydz);
	_eutracks->Branch("trackNum", &_trackNum);
	_eutracks->Branch("iden", &_trackIden);
	_eutracks->Branch("chi2", &_chi2);
	_eutracks->Branch("ndof", &_ndof);
	
	//Tree for storing external triggers
	_triggers = new TTree("triggers","triggers");
	_triggers->SetAutoSave(1000000000);
	_triggers->Branch("nTriggers", &_nTriggers);
	_triggers->Branch("nTLUTriggers", &_nTLUTriggers);
	_triggers->Branch("nExtTriggers", &_nExtTriggers);
	_triggers->Branch("TLUTrigTime", "std::vector<unsigned int>", &_TLUTrigTime);
	_triggers->Branch("ExtTrigTime", "std::vector<unsigned int>", &_ExtTrigTime);
	
	// Tree for storing ToTs
	_tots = new TTree("tots","tots");
	_tots->SetAutoSave(1000000000);
	_tots->Branch("nToTs", &_nToTs);
	_tots->Branch("ToTTime","std::vector<double>", &_ToTTime);
	_tots->Branch("ToTLength","std::vector<int>", &_ToTLength);

	_euhits->AddFriend(_zstree);
	_euhits->AddFriend(_eutracks);
	_euhits->AddFriend(_triggers);
	_euhits->AddFriend(_tots);
}
