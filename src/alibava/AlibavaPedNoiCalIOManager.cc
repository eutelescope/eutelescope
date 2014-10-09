/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */

// alibava includes ".h"
#include "ALIBAVA.h"
#include "AlibavaPedNoiCalIOManager.h"
#include "AlibavaRunHeaderImpl.h"

// lcio includes <.h>
#include <lcio.h>
#include "IO/LCReader.h"
#include "IO/LCWriter.h"
#include <UTIL/LCTOOLS.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <EVENT/LCParameters.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/TrackerDataImpl.h>

// system includes <>
#include <string>
#include <sys/stat.h>

using namespace std;
using namespace lcio;
using namespace alibava;

AlibavaPedNoiCalIOManager::AlibavaPedNoiCalIOManager(){
}

AlibavaPedNoiCalIOManager::~AlibavaPedNoiCalIOManager(){
}

EVENT::FloatVec AlibavaPedNoiCalIOManager::getDataFromEventForChip(LCEvent* evt, string collectionName, unsigned int chipnum){
	
	EVENT::FloatVec tmp_vec;
	tmp_vec.clear();
	
	LCCollectionVec* col = dynamic_cast< LCCollectionVec * > (evt->getCollection(collectionName)) ;
	
	int ielement = getElementNumberOfChip(col,chipnum);
	if (ielement!=-1) {
		TrackerDataImpl * trkdata = dynamic_cast< TrackerDataImpl * > ( col->getElementAt( ielement ) ) ;
		tmp_vec = trkdata->getChargeValues();
	}
	
	return tmp_vec;
}

int AlibavaPedNoiCalIOManager::getElementNumberOfChip(LCCollectionVec* col, int chipnum){
	int ielement = -1;
	CellIDDecoder<TrackerDataImpl> chipIDDecoder(col);
	unsigned int noOfChips = col->getNumberOfElements();
	for ( size_t i = 0; i < noOfChips; ++i )
	{
		
		TrackerDataImpl * trkdata = dynamic_cast< TrackerDataImpl * > ( col->getElementAt( i ) ) ;
		// check chip number
		const int ichip = static_cast<int> ( chipIDDecoder( trkdata )[ALIBAVA::ALIBAVADATA_ENCODE_CHIPNUM] );

		if (chipnum==ichip)
			ielement=i;
	}
	return ielement;
}


EVENT::FloatVec AlibavaPedNoiCalIOManager::getPedNoiCalForChip(string filename, string collectionName, unsigned int chipnum){
	
	// open pedestal file
	LCReader* lcReader = LCFactory::getInstance()->createLCReader() ;
	
	EVENT::FloatVec tmp_vec;
	tmp_vec.clear();
	
	try{
		lcReader->open( filename ) ;
		
		// check if there is only one run and only one event as it is supposed to
		if (lcReader->getNumberOfRuns() !=1 )
			streamlog_out( ERROR5 ) << " There are more than one run in AlibavaPedNoiCalFile: "<< filename<< endl ;
		if (lcReader->getNumberOfEvents() !=1 )
			streamlog_out( ERROR5 ) << " There are more than one event in AlibavaPedNoiCalFile: "<< filename<< endl ;
		
		LCEvent*  evt = lcReader->readNextEvent();
		
		tmp_vec = getDataFromEventForChip(evt,collectionName,chipnum);
		
		// if datavec is empty
		if (tmp_vec.size()==0)
			streamlog_out( ERROR5 ) <<"Trying to access"<<collectionName<<" for non existing chip ("<<chipnum<<")."<< endl;
		
		
		lcReader->close() ;
	}
	catch(IOException& e){
		streamlog_out( ERROR5 ) << " Unable to read the AlibavaPedNoiCal file - "<< filename<<e.what() << endl ;
	}
	
	//delete lcReader;
	return tmp_vec;
	
}

void AlibavaPedNoiCalIOManager::createFile(string filename, IMPL::LCRunHeaderImpl* runHeader){
	
	
	LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();
	
	try {
		lcWriter->open(filename, LCIO::WRITE_NEW);
	} catch (IOException& e) {
		cerr << e.what() << endl;
		return;
	}
	
	lcWriter->writeRunHeader(runHeader);
	
	lcWriter->close();
}




void AlibavaPedNoiCalIOManager::addToFile( string filename, string collectionName, int chipnum, EVENT::FloatVec datavec){

	// if file doesn't exist
	if (!doesFileExist(filename)) {
		streamlog_out( WARNING5 ) << " The AlibavaPedNoiCalFile: "<<filename<<" doesn't exist. "<< endl ;
		streamlog_out( WARNING5 ) << " Creating new AlibavaPedNoiCalFile "<<endl;
		createFile(filename);
	}

	
	LCRunHeaderImpl* runHeader = getRunHeader(filename);
	LCEventImpl*  evt = getEvent(filename);
	
	FloatVec tmp_vec;
	tmp_vec.clear();
		
	LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();
	// we will write a new lcio file with the copied run header and event
	try {
		lcWriter->open(filename, LCIO::WRITE_NEW);
		
		// first write runheader
		lcWriter->writeRunHeader(runHeader);
		
		// check if the collection exists
		LCCollectionVec* newCol = new LCCollectionVec(LCIO::TRACKERDATA);


		if (doesCollectionExist(evt,collectionName)){
			LCCollectionVec* col = dynamic_cast < LCCollectionVec * > (evt->getCollection(collectionName));
			*newCol = *col;
			evt->removeCollection(collectionName);
		}

		// check if the data exists for this chip in this event
		// if exists remove it
		int ielement=0;
		do {
			ielement= getElementNumberOfChip(newCol,chipnum);
			if (ielement!=-1)
				newCol->removeElementAt(ielement);
		} while (ielement!=-1);

				
		// now, add data vector to the collecton
		TrackerDataImpl * tmp_data = new TrackerDataImpl();
		tmp_data->setChargeValues(datavec);
		
		// set Cell ID encode
		CellIDEncoder<TrackerDataImpl> chipIDEncoder(ALIBAVA::ALIBAVADATA_ENCODE,newCol);
		
		chipIDEncoder[ALIBAVA::ALIBAVADATA_ENCODE_CHIPNUM] = chipnum;
		chipIDEncoder.setCellID(tmp_data);

		newCol->push_back(tmp_data);
		evt->addCollection(newCol, collectionName);
		
		lcWriter->writeEvent(evt);
		lcWriter->close();
	}
	catch (IOException& e) {
		cerr << e.what() << endl;
	}
	//delete lcWriter;
	
}

bool AlibavaPedNoiCalIOManager::doesCollectionExist(LCEvent* evt, string collectionName){

	bool it_exists = false;
//	vector<string> * colnames = (vector<string> *) evt->getCollectionNames();
	StringVec * colnames = const_cast<StringVec * > (evt->getCollectionNames());
	
	for (unsigned int i=0; i<colnames->size(); i++)
		if (colnames->at(i)==collectionName) it_exists = true;
	
	return it_exists;
}

bool AlibavaPedNoiCalIOManager::doesFileExist(string astring){
	struct stat buffer;
	return (stat (astring.c_str(), &buffer) == 0);
}

void AlibavaPedNoiCalIOManager::createFile(string filename){
	LCRunHeaderImpl* runHeader = new LCRunHeaderImpl();
	int runnumber = 0;
	runHeader->setRunNumber(runnumber);
	createFile(filename,runHeader);
	streamlog_out( WARNING5 ) << " New AlibavaPedNoiCalFile is created with empty header. File is: " << filename <<endl;
	
}

LCRunHeaderImpl* AlibavaPedNoiCalIOManager::getRunHeader(string filename){
	
	LCRunHeaderImpl* newRunHeader = new LCRunHeaderImpl();
	
	// first try to access file and get run header and the event
	LCReader* lcReader = LCFactory::getInstance()->createLCReader() ;
	
	try{
		// first try to access file
		lcReader->open( filename ) ;
		
		// read run header
		if ( (lcReader->getNumberOfRuns()) >0 ){
			if ( (lcReader->getNumberOfRuns()) >1 ){
				streamlog_out( ERROR5 ) << " There are more than one run in AlibavaPedNoiCalFile: "<< filename<<endl;
				streamlog_out( ERROR5 ) << " Using only one of them (probably first one)! Might cause problems!"<<endl ;
			}
			newRunHeader = dynamic_cast< LCRunHeaderImpl* > (lcReader->readNextRunHeader(LCIO::UPDATE));
		}
		else if ( (lcReader->getNumberOfRuns()) ==0 ){
			streamlog_out( ERROR5 ) << " There is no run header in AlibavaPedNoiCalFile: "<< filename<< endl;
		}
		
		lcReader->close() ;
	}
	catch(IOException& e){
		cerr << e.what() << endl;
	}
	//delete lcReader;
	return newRunHeader;
	
}

LCEventImpl* AlibavaPedNoiCalIOManager::getEvent(string filename){
	
	LCEventImpl* newEvent = new LCEventImpl();
	
	// first try to access file and get run header and the event
	LCReader* lcReader = LCFactory::getInstance()->createLCReader() ;
	
	try{
		// first try to access file
		lcReader->open( filename ) ;
		
		// read the event
		if (lcReader->getNumberOfEvents() >0 ){
			if (lcReader->getNumberOfEvents() >1 ){
				streamlog_out( ERROR5 ) << " There are more than one event in AlibavaPedNoiCalFile: "<< filename<<endl;
				streamlog_out( ERROR5 ) << " Using only one of them (probably first one)! Might cause problems!"<<endl;
			}
			newEvent = dynamic_cast< LCEventImpl* > ( lcReader->readNextEvent(LCIO::UPDATE) );
		}
		else if (lcReader->getNumberOfEvents() ==0 ){
			streamlog_out( MESSAGE5 ) << " New event created in AlibavaPedNoiCalFile: "<< filename<< endl ;
		}
		
		lcReader->close() ;
	}
	catch(IOException& e){
		cerr << e.what() << endl;
	}
	// clean up
	//delete lcReader;
	return newEvent;
	
}
