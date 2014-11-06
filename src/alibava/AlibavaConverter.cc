/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */



// personal includes ".h"
#include "ALIBAVA.h"
#include "AlibavaConverter.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"

// marlin includes
#include "marlin/Global.h"
#include "marlin/Exceptions.h"
#include "marlin/Processor.h"
#include "marlin/DataSourceProcessor.h"
#include "marlin/ProcessorMgr.h"


// lcio includes
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCGenericObjectImpl.h>
#include <IMPL/LCEventImpl.h>
#include <UTIL/CellIDEncoder.h>

// system includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <memory>
#include <stdlib.h>
#include <algorithm>

using namespace std;
using namespace marlin;
using namespace alibava;

AlibavaConverter::AlibavaConverter ():
DataSourceProcessor("AlibavaConverter"),
_fileName(ALIBAVA::NOTSET),
_geoID(0),
_formattedRunNumber("0"),
_runNumber(-1),
_rawDataCollectionName("rawdata"),
_chipSelection(),
_tiltAngle(0.0),
_sensorTemperature(111.111),
_startEventNum(-1),
_stopEventNum(-1),
_readChannelsReverse(false),
_storeHeaderPedestalNoise(false)
{
	
	
	// initialize few variables
	
	_description = "Reads data streams produced by Alibava and produces the corresponding LCIO output";
	
	registerProcessorParameter("InputFileName", "This is the input file name",
										_fileName, string("runXXXXXX.dat") );
	
	registerProcessorParameter("GeoID", "The geometry identification number", _geoID, static_cast<int> ( 0 ));
	
	registerProcessorParameter("RunNumber", "Run number of file (formatted)",
										_formattedRunNumber, string("0") );
	registerOutputCollection (LCIO::TRACKERDATA, "RawDataCollectionName",
									  "Name of the collection",
									  _rawDataCollectionName, string("rawdata") );
	
	
	// now optional parameters
	registerOptionalParameter("ReadChannelsReverse", "Alibava read channels from right to left if you want to revert this i.e. make it from left to right set this parameter to true. CAUTION: This will be applied first!",
									  _readChannelsReverse, bool(false) );
	
	registerOptionalParameter("ChipSelection", "Selection of chip that you want to store data from. Chip numbers start from 0. If not set, all data (i.e. chip 0 and 1) will be stored",
									  _chipSelection, EVENT::IntVec() );
	
	registerOptionalParameter("TiltAngle", "The tilt angle of the sensors",
									  _tiltAngle, float(0.0) );
	
	registerOptionalParameter("SensorTemperature", "The temperature of the sensors",
									  _sensorTemperature, float(111.111) );
	
	registerOptionalParameter("StartEventNum", "The event number that AlibavaConverter should start storing. Default value is -1, in this case it will store every event",
									  _startEventNum, int(-1) );
	
	registerOptionalParameter("StopEventNum", "The event number that AlibavaConverter should stop storing. Default value is -1, in this case it will store every event",
									  _stopEventNum, int(-1) );
	registerOptionalParameter("StoreHeaderPedestalNoise", "Alibava stores a pedestal and noise set in the run header. These values are not used in te rest of the analysis, so it is optional to store it. By default it will not be stored, but it you want you can set this variable to true to store it in the header of slcio file",
									  _storeHeaderPedestalNoise, bool(false) );
	
	
}

AlibavaConverter * AlibavaConverter::newProcessor () {
	return new AlibavaConverter;
}



void AlibavaConverter::init () {
	
	// print out warning if _readChannelsReverse
	if (_readChannelsReverse) {
		streamlog_out( WARNING5 )<< "You select ReadChannelsReverse="<<_readChannelsReverse<<endl;
		streamlog_out( WARNING5 )<<"This will be applied before chip selection and channel masking!"<<endl;
	}
	
	checkIfChipSelectionIsValid();
	
	
	if (_startEventNum!=-1 && _stopEventNum==-1){
		streamlog_out( WARNING5 )<< "First "<<_startEventNum <<" will be skipped!"<<endl;
	}
	else if (_startEventNum!=-1 && _stopEventNum!=-1){
		streamlog_out( WARNING5 )<< "Only events from "<<_startEventNum <<" to "<<_stopEventNum<<" will be stored! "<<endl;
	}
	else if (_startEventNum==-1 && _stopEventNum!=-1){
		streamlog_out( WARNING5 )<< "Events after event number: "<<_stopEventNum<<" will be skipped! "<<endl;
	}
	
	
	printParameters ();
}

void AlibavaConverter::readDataSource(int /* numEvents */) {
	
	// this event counter is used to stop the processing when it is
	// greater than numEvents.
	int eventCounter = 0;
	
	
	// this is to make the output messages nicer
	streamlog::logscope scope(streamlog::out);
	scope.setName(name());
	
	// print out a debug message
	streamlog_out( MESSAGE5 ) << "Reading " << _fileName << " with AlibavaConverter " << endl;
	_runNumber = atoi(_formattedRunNumber.c_str());
	
	/////////////////
	//  Open File  //
	/////////////////
	ifstream infile;
	infile.open(_fileName.c_str());
	if(!infile.is_open()) {
		streamlog_out( ERROR5 ) << "AlibavaConverter could not read the file "<<_fileName<<" correctly. Please check the path and file names that have been input" << endl;
		exit(-1);
	}
	else streamlog_out( MESSAGE4 )<<"Input file "<<_fileName<<" is opened!"<<endl;
	
	time_t date;
	int type;
	unsigned int lheader; // length of the header
	string header;
	int version; // Alibava firmware version
	
	
	/////////////////
	// Read Header //
	/////////////////
	infile.read(reinterpret_cast< char *> (&date), sizeof(time_t));
	infile.read(reinterpret_cast< char *> (&type), sizeof(int));
	
	header.clear();
	infile.read(reinterpret_cast< char *> (&lheader), sizeof(unsigned int)); //length of header
	for (unsigned int ic=0; ic<lheader; ic++){
		char tmp_c;
		infile.read(&tmp_c, sizeof(char));
		header.append(1, tmp_c);
	}
	
	header = trim_str(header);
	
	if (header[0]!='V' && header[0]!='v')
	{
		version = 0;
	}
	else
	{
		version = int(header[1]-'0');
		header = header.substr(5);
	}
	
	////////////////////////////////////
	// Read header pedestal and noise //
	////////////////////////////////////
	
	// Alibava stores a pedestal and noise set in the run header. These values are not used in te rest of the analysis, so it is optional to store it. By default it will not be stored, but it you want you can set _storeHeaderPedestalNoise variable to true.
	float tmp_float;
	FloatVec headerPedestal;
	FloatVec headerNoise;
	
	// first pedestal
	for (int ichan=0; ichan<ALIBAVA::NOOFCHIPS*ALIBAVA::NOOFCHANNELS; ichan++) {
		infile.read(reinterpret_cast< char *> (&tmp_float), sizeof(double));
		headerPedestal.push_back(tmp_float);
	}
	// now noise
	for (int ichan=0; ichan<ALIBAVA::NOOFCHIPS*ALIBAVA::NOOFCHANNELS; ichan++) {
		infile.read(reinterpret_cast< char *> (&tmp_float), sizeof(double));
		headerNoise.push_back(tmp_float);
	}
	
	////////////////////
	// Process Header //
	////////////////////
	
	LCRunHeaderImpl * arunHeader = new LCRunHeaderImpl();
	AlibavaRunHeaderImpl* runHeader = new AlibavaRunHeaderImpl(arunHeader);

	runHeader->setDetectorName(Global::GEAR->getDetectorName());
	runHeader->setHeader(header);
	runHeader->setHeaderVersion(version);
	runHeader->setDataType(type);
	runHeader->setDateTime(string(ctime(&date)));
	if (_storeHeaderPedestalNoise) {
		runHeader->setHeaderPedestal(headerPedestal);
		runHeader->setHeaderNoise(headerNoise);
	}
	runHeader->setRunNumber(_runNumber);
	runHeader->setTiltAngle(_tiltAngle);
	runHeader->setSensorTemperature(_sensorTemperature);
	runHeader->setChipSelection(_chipSelection);
	
	// get number of events from header
	string tmpstring = getSubStringUpToChar(header,";",0);
	int noofevents = atoi(tmpstring.c_str());
	runHeader->setNoOfEvents(noofevents);
	
	
	//runHeader->addProcessor(type());
	
	ProcessorMgr::instance()->processRunHeader( runHeader->lcRunHeader() ) ;
	
	delete arunHeader;
	delete runHeader;
	
	
	////////////////
	// Read Event //
	////////////////
	
	
	if (version<2) {
		// this code is not written for version<=1.
		streamlog_out( ERROR5 )<<" Unexpected data version found (version="<<version<<"<2). Data is not saved"<<endl;
		return;
	}
	
	do
	{
		
		if ( eventCounter % 1000 == 0 )
			streamlog_out ( MESSAGE4 ) << "Processing event "<< eventCounter << " in run " << _runNumber<<endl;
		
		
		unsigned int headerCode, eventSize, userEventTypeCode=0, eventTypeCode=0;
		do
		{
			infile.read(reinterpret_cast< char *> (&headerCode), sizeof(unsigned int));
			if (infile.bad() || infile.eof())
				return;
			
			eventTypeCode = (headerCode>>16) & 0xFFFF;
		} while ( eventTypeCode != 0xcafe );
		
		eventTypeCode = headerCode & 0x0fff;
		userEventTypeCode = headerCode & 0x1000;
		
		if (userEventTypeCode){
			streamlog_out( ERROR5 )<<" Unexpected data type found (type= User type). Data is not saved"<<endl;
			return;
		}
		
		infile.read(reinterpret_cast< char *> (&eventSize), sizeof(unsigned int));
		
		double value, charge, delay;
		infile.read(reinterpret_cast< char *> (&value), sizeof(double));
		
		//see AlibavaGUI.cc
		charge = int(value) & 0xff;
		delay = int(value) >> 16;
		charge = charge * 1024;
		
		unsigned int tdcTime;
		unsigned short temp;  // temperature measured on Daughter board
		//  ifile->read((char *)&_t0, sizeof(time_t));
		
		infile.read(reinterpret_cast< char *> (&tdcTime), sizeof(unsigned int));
		infile.read(reinterpret_cast< char *> (&temp), sizeof(unsigned short));
		
		
		unsigned short chipHeader[2][16];
		short tmp_short;
		FloatVec all_data;
		all_data.reserve(ALIBAVA::NOOFCHIPS*ALIBAVA::NOOFCHANNELS);
		FloatVec::iterator it;
		for (int ichip=0; ichip<ALIBAVA::NOOFCHIPS; ichip++)
		{
			infile.read(reinterpret_cast< char *> (chipHeader[ichip]), 16*sizeof(unsigned short));
			for (int ichan=0; ichan<ALIBAVA::NOOFCHANNELS; ichan++) {
				infile.read(reinterpret_cast< char *> (&tmp_short), sizeof(unsigned short));
				if (_readChannelsReverse) {
					it = all_data.begin();
					all_data.insert(it, float(tmp_short));
				}
				else
					all_data.push_back(float(tmp_short));
			}
		}
		
		///////////////////
		// Process Event //
		///////////////////
		
		
		// now write these to AlibavaEvent
		AlibavaEventImpl* anEvent = new AlibavaEventImpl();
		anEvent->setRunNumber(_runNumber);
		anEvent->setEventNumber(eventCounter);
		anEvent->setEventType(eventTypeCode);
		anEvent->setEventSize(eventSize);
		anEvent->setEventValue(value);
		anEvent->setEventTime(tdc_time(tdcTime));
		anEvent->setEventTemp(get_temperature(temp));
		anEvent->setCalCharge(charge);
		anEvent->setCalDelay(delay);
		anEvent->unmaskEvent();
		
		
		// creating LCCollection
		LCCollectionVec* rawDataCollection = new LCCollectionVec(LCIO::TRACKERDATA);
		CellIDEncoder<TrackerDataImpl> chipIDEncoder(ALIBAVA::ALIBAVADATA_ENCODE,rawDataCollection);
		
		// for this to work the _chipselection has to be sorted in ascending order!!!
		for (unsigned int ichip=0; ichip<_chipSelection.size(); ichip++) {
			FloatVec chipdata;
			chipdata.clear();
			
			chipdata.insert(chipdata.end(), all_data.begin()+_chipSelection[ichip]*ALIBAVA::NOOFCHANNELS, all_data.begin()+(_chipSelection[ichip]+1)*ALIBAVA::NOOFCHANNELS);
			TrackerDataImpl * arawdata = new TrackerDataImpl();
			arawdata->setChargeValues(chipdata);
			chipIDEncoder[ALIBAVA::ALIBAVADATA_ENCODE_CHIPNUM] = _chipSelection[ichip];
			chipIDEncoder.setCellID(arawdata);
			rawDataCollection->push_back(arawdata);
		}
		
		anEvent->addCollection(rawDataCollection, _rawDataCollectionName);
		
		
		if (_startEventNum!=-1 && eventCounter<_startEventNum) {
			streamlog_out( MESSAGE5 )<<" Skipping event "<<eventCounter<<". StartEventNum is set to "<<_startEventNum<<endl;
			eventCounter++;
			continue;
		}
		
		if (_stopEventNum!=-1 && eventCounter>_stopEventNum) {
			streamlog_out( MESSAGE5 )<<" Reached StopEventNum: "<<_stopEventNum<<". Last saved event number is "<<eventCounter<<endl;
			break;
		}
		
		ProcessorMgr::instance()->processEvent( static_cast<LCEventImpl*> ( anEvent ) ) ;
		eventCounter++;
		
		delete anEvent;
		
	} while ( !(infile.bad() || infile.eof()) );
	
	infile.close();
	
	if (_stopEventNum!=-1 && eventCounter<_stopEventNum)
		streamlog_out( MESSAGE5 )<<" Stooped before reaching StopEventNum: "<<_stopEventNum<<". The file has "<<eventCounter<<" events."<<endl;
	
	
	
	
}


void AlibavaConverter::end () {
	
	streamlog_out ( MESSAGE5 )  << "AlibavaConverter Successfully finished" << endl;
}

double AlibavaConverter::tdc_time(unsigned int tdcTime){
	unsigned short fpart = tdcTime & 0xffff;
	short ipart = (tdcTime & 0xffff0000)>>16;
	if (ipart<0)
		fpart *= -1;
	//double tt = 100.*(1. -(ipart + (fpart/65535.)));
	double tt = 100.0*(ipart + (fpart/65535.));
	return tt;
}

double AlibavaConverter::get_temperature(unsigned short temp){
	if (temp==0)
		return 9999.;
	else
		return 0.12*temp - 39.8;
}

void AlibavaConverter::checkIfChipSelectionIsValid(){
	
	bool resetChipSelection = false;
	
	// check if there is chip selection or if there is chip selection but not valid
	if (_chipSelection.size()==0){
		streamlog_out( WARNING5 )<< "You didn't select any chip"<<endl;
		resetChipSelection = true;
	}
	else if (int(_chipSelection.size())> ALIBAVA::NOOFCHIPS) {
		streamlog_out( WARNING5 )<< "The number of chips you selected ( "<<_chipSelection.size()<<" ) is more than an alibava daughter board can have ( "<<ALIBAVA::NOOFCHIPS<<" )"<<endl;
		resetChipSelection= true;
	}
	else{
		// first sort the chip numbers in ascending order
		sort(_chipSelection.begin(),_chipSelection.end());
		
		// check if the selected chips make sense
		for (int ichip=0; ichip<int(_chipSelection.size()); ichip++) {
			bool del_this_chip = false;
			
			if (_chipSelection[ichip]<0){
				streamlog_out( ERROR5 )<< "Selected chip cannot have negative value. "<< endl;
				del_this_chip =true;
			}
			else if (_chipSelection[ichip]>=ALIBAVA::NOOFCHIPS){
				streamlog_out( ERROR5 )<< "Chip numbering has to start from zero \"0\" and cannot be greater than " << ALIBAVA::NOOFCHIPS-1 << endl;
				del_this_chip = true;
			}
			
			if (del_this_chip) { // if this chip selection is not valid, delete it.
				streamlog_out( ERROR5 )<< "Chip "<<_chipSelection[ichip]<<" is deleted from the chip selection list"<< endl;
				_chipSelection.erase(_chipSelection.begin()+ichip);
				ichip = ichip-1;
			}
		}// end of ichip loop
		
		// check again if there is any selected chip left
		if (_chipSelection.size()==0) resetChipSelection= true;
		
		
	}
	
	if (resetChipSelection) {
		streamlog_out( WARNING5 )<< "I will save data from all chips"<< endl;
		
		_chipSelection.clear();
		
		for (int ichip=0; ichip<ALIBAVA::NOOFCHIPS; ichip++) {
			_chipSelection.push_back(ichip);
		}
	}
	
	
	// now there is valid chip selection!!!
	
	streamlog_out( WARNING5 )<< "Final applied chip selection: ";
	for (int ichip=0; ichip<int(_chipSelection.size()); ichip++)
		streamlog_out( WARNING5 )<< _chipSelection[ichip] << " ";
	streamlog_out( WARNING5 )<<endl;
	streamlog_out( WARNING5 )<< "Only data coming from these chips will be stored"<<endl;
	
}

