// eutelescope inlcudes
#include "EUTelAPIXTbTrackTupleGBL.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelTrackerDataInterfacerImpl.h"

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

EUTelAPIXTbTrackTupleGBL::EUTelAPIXTbTrackTupleGBL()
: Processor("EUTelAPIXTbTrackTupleGBL"),
  _inputTrackColName (""),
  _inputTelPulseCollectionName(""),
  _inputDutPulseCollectionName(""),
  _telZsColName(""),
  _path2file(""),
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
  _euhits(NULL),
  _nHits(0),
  _hitXPos(NULL),
  _hitYPos(NULL),
  _hitZPos(NULL),
  _hitSensorId(NULL)
 {
  //processor description
  _description = "Prepare tbtrack style n-tuple with track fit results" ;


  registerInputCollection(LCIO::TRACK, "InputTrackCollectionName", "Name of the input Track collection",
		  		_inputTrackColName, std::string("fittracks"));
    registerProcessorParameter ("DutZsColName", "DUT zero surpressed data colection name",
                         _dutZsColName, std::string("zsdata_apix"));
  	   
  registerProcessorParameter ("OutputPath", "Path/File where root-file should be stored",
			      _path2file, std::string("NTuple.root"));

}


void EUTelAPIXTbTrackTupleGBL::init()
{
    EUTelExcludedPlanes();

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


	for(std::vector<int>::iterator it = EUTelExcludedPlanes::_senNoDeadMaterial.begin(); it != EUTelExcludedPlanes::_senNoDeadMaterial.end(); it++)
	{
		//Later we need to shift the sensor since in EUTel centre of sensor is 0|0 while in TBmon(II) it is in the lower left corner
        if(*it > 5){
            geo::EUTelGenericPixGeoDescr* geoDescr = geo::gGeometry().getPixGeoDescr( *it ) ;
            float xSize,ySize;
            geoDescr->getSensitiveSize(xSize, ySize);

            _xSensSize[*it] = xSize;
            _ySensSize[*it] = ySize;
        }
	}
}

void EUTelAPIXTbTrackTupleGBL::processRunHeader( LCRunHeader* runHeader) 
{
	std::auto_ptr<EUTelRunHeaderImpl> eutelHeader( new EUTelRunHeaderImpl ( runHeader ) );
	eutelHeader->addProcessor( type() );
	_nRun++;
	
	// Decode and print out Run Header information - just a check
	_runNr = runHeader->getRunNumber();
}

void EUTelAPIXTbTrackTupleGBL::processEvent( LCEvent * event )
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
    EUTelReaderGenericLCIO reader = EUTelReaderGenericLCIO();
    std::vector<EUTelTrack> tracks = reader.getTracks(event, _inputTrackColName );

	//try to read in hits (e.g. fitted hits in local frame)	
  	if( !readHits( tracks )) 
	{ 
		return; 
	}
  	
	//read in raw data	
    if(!readZsHits( _dutZsColName , event )) 
	{ 
		return; 
	}
  
  	//read in tracks
	if(!readTracks(tracks))
       	{
		return;
	}
 
        //fill the trees	
	_zstree->Fill();
	_eutracks->Fill();
	_euhits->Fill();

	_isFirstEvent = false;
}

void EUTelAPIXTbTrackTupleGBL::end()
{
	//write version number
	_versionNo->push_back(1.1);
	_versionTree->Fill();
	//Maybe some stats output?
	_file->Write();
}

//Read in TrackerHit(Impl) to later dump them
bool EUTelAPIXTbTrackTupleGBL::readHits( std::vector<EUTelTrack>& tracks )
{
 
  	for( std::vector<EUTelTrack>::iterator iTrack= tracks.begin() ; iTrack != tracks.end() ; iTrack++){
        for( std::vector<EUTelState>::iterator iState= iTrack->getStates().begin() ; iState != iTrack->getStates().end() ; iState++){
            int sensorID = iState->getLocation();
            //Only dump DUT hits
            if(  sensorID > 5 ){
                TVector3 pos;
                pos[0]=-999; pos[1]=-999; pos[2]=-999;
                if(iState->getStateHasHit()){
                    pos = iState->getHit().getPositionGlobal();	
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
        }
    }
	return true;
}

//Read in TrackerHit to later dump
bool EUTelAPIXTbTrackTupleGBL::readTracks(std::vector<EUTelTrack>& tracks)
{
    unsigned int trackCounter=0;
	unsigned int nTrackParams=0;

    for( std::vector<EUTelTrack>::iterator iTrack= tracks.begin() ; iTrack != tracks.end() ; iTrack++){
		double chi2 = iTrack->getChi2();
		double ndof = iTrack->getNdf();
        for( std::vector<EUTelState>::iterator iState= iTrack->getStates().begin() ; iState != iTrack->getStates().end() ; iState++){
            TVector3 pos = iState->getPositionGlobal();	

            double dxdz = iState->getSlopeXGlobal();
            double dydz = iState->getSlopeYGlobal();
            int sensorID = iState->getLocation();

            if( sensorID > 5 ){
                    nTrackParams++;
                    double x = pos[0];
                    double y = pos[1];
                    //eutrack tree
                    _xPos->push_back(x);
                    _yPos->push_back(y);
                    _dxdz->push_back(dxdz);
                    _dydz->push_back(dydz);
                    _trackIden->push_back(sensorID);
                    _trackNum->push_back(trackCounter);
                    _chi2->push_back(chi2);
                    _ndof->push_back(ndof);
            }
        }
        trackCounter++;
    }

	_nTrackParams = nTrackParams;
	return true;
}
bool EUTelAPIXTbTrackTupleGBL::readZsHits( std::string colName, LCEvent* event)
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

        if (type == kEUTelGenericSparsePixel  ) 
        {
            std::auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> > apixData( new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> ( zsData ));
            EUTelGenericSparsePixel apixPixel;

            for( unsigned int iHit = 0; iHit < apixData->size(); iHit++ ) 
            {
                apixData->getSparsePixelAt( iHit, &apixPixel);
                _nPixHits++;
                p_iden->push_back( sensorID );
                p_row->push_back( apixPixel.getYCoord() );
                p_col->push_back( apixPixel.getXCoord() );
                p_tot->push_back( static_cast< int >(apixPixel.getSignal()) );
                p_lv1->push_back( static_cast< int >(apixPixel.getTime()) );
            }
        }
        else
        {
        }
    }
        return true;
}


void EUTelAPIXTbTrackTupleGBL::clear()
{
	/* Clear zsdata */
	p_col->clear();
	p_row->clear();
	p_tot->clear();
	p_iden->clear();
	p_lv1->clear();
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
}

void EUTelAPIXTbTrackTupleGBL::prepareTree()
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
	p_tot = new std::vector<int>();
	p_iden = new std::vector<int>();
	p_lv1 = new std::vector<int>();

	_hitXPos = new std::vector<double>();
	_hitYPos = new std::vector<double>();
	_hitZPos = new std::vector<double>();
	_hitSensorId  = new std::vector<int>();

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
	_zstree->Branch("tot",      &p_tot);
	_zstree->Branch("lv1",      &p_lv1);
	_zstree->Branch("iden",     &p_iden);

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

	_euhits->AddFriend(_zstree);
	_euhits->AddFriend(_eutracks);
}
