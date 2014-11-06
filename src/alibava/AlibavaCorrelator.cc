/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


// alibava includes ".h"
#include "AlibavaCorrelator.h"
#include "AlibavaBaseHistogramMaker.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"

// eutelescope includes ".h"
#include "EUTELESCOPE.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"
#include "marlin/tinyxml.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <lcio.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/TrackerHitImpl.h>


// ROOT includes ".h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

// system includes <>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <memory>
#include <cstdlib>
#include <ctime>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;


AlibavaCorrelator::AlibavaCorrelator () :
AlibavaBaseHistogramMaker("AlibavaCorrelator"),
// List of Histogram names, initialized here. As an example we put only 2
_hHitPosX("hHitPosX"),
_hHitPosY("hHitPosY"),
_hCorX("hCorX"),
_hCorY("hCorY"),
_hSyncX("hSyncX"),
_hSyncY("hSyncY"),
_detectorIDs()
{
	
	// modify processor description
	_description =
	"AlibavaCorrelator gets hit collection and plots correlation histograms ";
	
	
	// first of register the input collection
	registerInputCollection (LCIO::TRACKERHIT, "InputCollectionName",
									 "Input raw data collection name",
									 _inputCollectionName, string("rawdata") );
	// Details about HistoXMLFile
	registerProcessorParameter ("HistoXMLFile",
										 "The path of XML file where the histograms are defined",
										 _histoXMLFileName , string("AlibavaHistoList.xml"));
	
	registerProcessorParameter ("TopTagInXMLFile",
										 "The top tag in HistoXMLFile",
										 _topTag , string("AlibavaHistoList"));
	
	registerProcessorParameter ("TagToProcess",
										 "The tag in TopTagInXMLFile. This processor will only consider the histogram definitions inside this tag. This tag should be inside <TopTagInXMLFile> ... <TopTagInXMLFile/>",
										 _tagToProcess , string("myAlibavaCorrelator"));
	
	registerProcessorParameter ("DetectorIDs",
										 "The list of detector IDs",
										 _detectorIDs , IntVec());
	
	
}


void AlibavaCorrelator::init () {
	streamlog_out ( MESSAGE4 ) << "Running init" << endl;
	
	/* To choose if processor should skip masked events
	 ex. Set the value to 0 for false, to 1 for true
	 */
	if (Global::parameters->isParameterSet(ALIBAVA::SKIPMASKEDEVENTS))
		_skipMaskedEvents = bool ( Global::parameters->getIntVal(ALIBAVA::SKIPMASKEDEVENTS) );
	else {
		streamlog_out ( MESSAGE4 ) << "The Global Parameter "<< ALIBAVA::SKIPMASKEDEVENTS <<" is not set! Masked events will be used!" << endl;
	}
	
	// here sort _detectorIDs
	std::sort(_detectorIDs.begin(),_detectorIDs.end());
	
	// this method is called only once even when the rewind is active
	// usually a good idea to
	printParameters ();
	
	// if you want to change any variable defined in HistoXMLFile use this function
	// see example function below
	createRulesToChangeXMLValues();
}

void AlibavaCorrelator::createRulesToChangeXMLValues(){
	// here we define which variables in HistoXMLFile we want to change or modify
	// You can only change or modify these variables:
	// xBin, xMin, xMax, yBin, yMin, yMax, title, labelX and labelY
	// others cannot be changed!
	
	
	/*! If you want to replace one of these variables: xBin, xMin, xMax, yBin, yMin, yMax
	 *  We will use the function
	 *      		void changeXMLVariable(string histoName, string variableName, float newValue);
	 *          defined in AlibavaBaseHistogramMaker
	 *          Note that variableNames are case sensitive!
	 */
	
	//  For example usually it is good idea to replace maximum event number for example
	//  Lets say _someOtherHistoName histograms X axis should be total number of events
	//	changeXMLVariable(_someOtherHistoName, "xMax", float(_totalNumberOfEvents));
	// And _someOtherHistoName histograms Y axis should be 1000 since it is the max value it can get
	//	changeXMLVariable(_someHistoName, "yMax", 1000.0);
	// It is really bad idea, but for some reason if you want to hard code the binning of X axis, here how you can do. Note that float number will be changed to integer in AlibavaBaseHistogramMaker::updateXMLVariables() method but here we should give it as float
	//	changeXMLVariable(_someOtherHistoName, "xBin", float(100));
	
	/*! If you want to replace or modify one of these variables: title, labelX and labelY
	 *  We will use the function
	 *      		void addToXMLTitle(string histoName, string titleName, string whichSide, string stringToAdd);
	 *          defined in AlibavaBaseHistogramMaker
	 *          Note that titleName and whichSide are case sensitive!
	 */
	
	// Say that you have the signal in X axis of _someHistoName.
	// Assuming that "Signal" is written in HistoXMLFile as this histograms labelY, to get "(_multiplySignalby) Signal" you should add signalMultipliedby string to the left
	//		addToXMLTitle(_someHistoName, "labelY", "left", signalMultipliedby);
	// Again it is usually bad idea but if you want to replace the X label here how you can do
	//		addToXMLTitle(_someOtherHistoName, "labelX", "all", string("Event Number"));
	// And here how you add string to right of title
	//		addToXMLTitle(_someHistoName, "title", "right", string("Some thing"));
	
	
}


void AlibavaCorrelator::processRunHeader (LCRunHeader * rdr) {
	streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;
	
	// Add processor name to the runheader
	auto_ptr<AlibavaRunHeaderImpl> arunHeader ( new AlibavaRunHeaderImpl(rdr)) ;
	arunHeader->addProcessor(type());
	
	// set total number of events in this run
	_totalNumberOfEvents = arunHeader->getNoOfEvents();
	streamlog_out ( DEBUG1 ) << "N events "<<_totalNumberOfEvents << endl;
	
	// reads and creates the histograms defined in HistoXMLFile
	bookHistos();
	
	// set number of skipped events to zero (defined in AlibavaBaseProcessor)
	_numberOfSkippedEvents = 0;
	
	
}
std::string AlibavaCorrelator::getHistoNameForDetector(std::string name, int detID){
	stringstream s;
	s<< name <<"_d" << detID;
	return s.str();
	
}
std::string AlibavaCorrelator::getHistoNameForDetector(std::string name, int detID1, int detID2){
	stringstream s;
	s<< name <<"_d" << detID1 << "_d" << detID2;
	return s.str();
	
}

void AlibavaCorrelator::fillListOfHistos(){
	addToHistoCheckList(_hHitPosX);
	addToHistoCheckList(_hHitPosY);
	addToHistoCheckList(_hCorX);
	addToHistoCheckList(_hCorY);
	addToHistoCheckList(_hSyncX);
	addToHistoCheckList(_hSyncY);
	
	checkListOfHistosCreatedByXMLFile();
	
}

void AlibavaCorrelator::createClones_hHitPos(string histoName){
	
	AIDAProcessor::tree(this)->cd(this->name());
	AIDAProcessor::tree(this)->cd(getInputCollectionName().c_str());
	TH1D * h = dynamic_cast<TH1D*> (_rootObjectMap[histoName]);
	
	streamlog_out ( MESSAGE1 )  << "hist "<< histoName<<" "<<h << endl;
	AIDAProcessor::tree(this)->mkdir(histoName.c_str());
	AIDAProcessor::tree(this)->cd(histoName.c_str());
	
	for (unsigned int idet=0; idet<_detectorIDs.size(); idet++) {
		int detID = _detectorIDs[idet];
		string newHistoName = getHistoNameForDetector(histoName, detID);
		TH1D * hnew = (TH1D*)h->Clone( newHistoName.c_str() );
		
		string title = hnew->GetTitle();
		title = title + string(" (det ") + to_string(detID) + string(")");
		_rootObjectMap.insert(make_pair(newHistoName, hnew));
	}
	
}
void AlibavaCorrelator::createClones_hCor(string histoName){
	AIDAProcessor::tree(this)->cd(this->name());
	AIDAProcessor::tree(this)->cd(getInputCollectionName().c_str());
	TH2D * h = dynamic_cast<TH2D*> (_rootObjectMap[histoName]);
	
	AIDAProcessor::tree(this)->mkdir(histoName.c_str());
	AIDAProcessor::tree(this)->cd(histoName.c_str());
	
	// the _detectorIDs vector has to be sorted
	for (unsigned int idet=0; idet<_detectorIDs.size(); idet++) {
		int detID = _detectorIDs[idet];
		for (unsigned int iCorDet=idet+1; iCorDet<_detectorIDs.size(); iCorDet++) {
			int corDetID = _detectorIDs[iCorDet]; // correlated detector id
			
			string newHistoName = getHistoNameForDetector(histoName, detID, corDetID);
			TH2D * hnew = (TH2D*)h->Clone( newHistoName.c_str() );
			string title = hnew->GetXaxis()->GetTitle();
			title = title + string(" det") + to_string(detID);
			
			title = hnew->GetYaxis()->GetTitle();
			title = title + string(" det") + to_string(corDetID);
			
			_rootObjectMap.insert(make_pair(newHistoName, hnew));
		}
	}
}
void AlibavaCorrelator::createClones_hSync(string histoName){
	AIDAProcessor::tree(this)->cd(this->name());
	AIDAProcessor::tree(this)->cd(getInputCollectionName().c_str());
	TH2D * h = dynamic_cast<TH2D*> (_rootObjectMap[histoName]);
	
	AIDAProcessor::tree(this)->mkdir(histoName.c_str());
	AIDAProcessor::tree(this)->cd(histoName.c_str());
	
	// the _detectorIDs vector has to be sorted
	for (unsigned int idet=0; idet<_detectorIDs.size(); idet++) {
		int detID = _detectorIDs[idet];
		for (unsigned int iCorDet=idet+1; iCorDet<_detectorIDs.size(); iCorDet++) {
			int corDetID = _detectorIDs[iCorDet]; // correlated detector id
			string newHistoName = getHistoNameForDetector(histoName, detID, corDetID);
			TH2D * hnew = (TH2D*)h->Clone( newHistoName.c_str() );
			string title = hnew->GetYaxis()->GetTitle();
			title = title + string(" d") + to_string(detID) + string(" - d")+to_string(corDetID);
			
			_rootObjectMap.insert(make_pair(newHistoName, hnew));
		}
	}
}

void AlibavaCorrelator::bookHistos(){
	// create histograms defined in HistoXMLFile
	processHistoXMLFile();
	
	
	// hX plots
	createClones_hHitPos(_hHitPosX);
	// hY plots
	createClones_hHitPos(_hHitPosY);
	// hCorX
	createClones_hCor(_hCorX);
	// hCorY
	createClones_hCor(_hCorY);
	// hSyncX
	createClones_hSync(_hSyncX);
	// hSyncY
	createClones_hSync(_hSyncY);
	
	streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << endl;
}

bool AlibavaCorrelator::isInDetectorIDsList(int detID){
	for (unsigned int i=0; i<_detectorIDs.size(); i++) {
		if (detID == _detectorIDs[i]) return true;
	}
	return false;
}

void AlibavaCorrelator::processEvent (LCEvent * anEvent) {
	
	if ( anEvent->getEventNumber() % 1000 == 0 )
		streamlog_out ( DEBUG1 ) << "Looping events "<<anEvent->getEventNumber() << endl;
	
	AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
	int eventnum = alibavaEvent->getEventNumber();
	
	/////////////////////////////
	// Now loop ever detectors //
	
	LCCollectionVec * collectionVec;
	unsigned int noOfHits;
	try{
		collectionVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( getInputCollectionName() ) ) ;
		CellIDDecoder<TrackerHitImpl> hitDecoder ( eutelescope::EUTELESCOPE::HITENCODING );
		
		noOfHits = collectionVec->getNumberOfElements();
		
		for ( size_t ihit = 0; ihit < noOfHits; ++ihit ){
			TrackerHitImpl * ahit = dynamic_cast< TrackerHitImpl * > ( collectionVec->getElementAt( ihit ) ) ;
			int detID = hitDecoder( ahit )["sensorID"];
			if ( !isInDetectorIDsList(detID) ) continue;

			const double* pos = ahit->getPosition();
			
			string histoName;
			// fill hX
			histoName = getHistoNameForDetector(_hHitPosX, detID);
			TH1D * hX = dynamic_cast<TH1D*> (_rootObjectMap[ histoName ]);
			hX->Fill(pos[0]);
			// fill hY
			histoName = getHistoNameForDetector(_hHitPosY, detID);
			TH1D * hY = dynamic_cast<TH1D*> (_rootObjectMap[ histoName ]);
			hY->Fill(pos[1]);
			
			// correlation plots
			for (size_t i = 0; i < noOfHits; ++i ) {
				TrackerHitImpl * anotherHit = dynamic_cast< TrackerHitImpl * > ( collectionVec->getElementAt( i ) ) ;
				int anotherDetID = hitDecoder( anotherHit )["sensorID"];
				// only consider hits from other detectors
				if (detID >= anotherDetID) continue;
				
				const double* anotherPos = anotherHit->getPosition();
				
				histoName = getHistoNameForDetector(_hCorX, detID, anotherDetID);
				TH2D * hCorX = dynamic_cast<TH2D*> (_rootObjectMap[ histoName ]);
				hCorX->Fill(pos[0], anotherPos[0]);
				
				histoName = getHistoNameForDetector(_hCorY, detID, anotherDetID);
				TH2D * hCorY = dynamic_cast<TH2D*> (_rootObjectMap[ histoName ]);
				hCorY->Fill(pos[1], anotherPos[1]);
				
				histoName = getHistoNameForDetector(_hSyncX, detID, anotherDetID);
				TH2D * hSyncX = dynamic_cast<TH2D*> (_rootObjectMap[ histoName ]);
				hSyncX->Fill(eventnum, pos[0]-anotherPos[0]);
				
				histoName = getHistoNameForDetector(_hSyncY, detID, anotherDetID);
				TH2D * hSyncY = dynamic_cast<TH2D*> (_rootObjectMap[ histoName ]);
				hSyncY->Fill(eventnum, pos[1]-anotherPos[1]);
				
			}
			
		}
		
	} catch ( lcio::DataNotAvailableException ) {
		// do nothing again
		streamlog_out( ERROR5 ) << "Collection ("<<getInputCollectionName()<<") not found! " << endl;
	}
	
}


void AlibavaCorrelator::fillHistos(TrackerDataImpl * /* trkdata */ ){
	// nothing here
}

void AlibavaCorrelator::check (LCEvent * /* evt */ ) {
	// nothing to check here
}


void AlibavaCorrelator::end() {
	
	if (_numberOfSkippedEvents > 0)
		streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents<<" events skipped since they are masked" << endl;
	streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
	
}

// You have to define what bookEventHisto will do
void AlibavaCorrelator::bookEventHisto(int ){
	// does nothing
}

// You have to define what fillEventHisto will do
void AlibavaCorrelator::fillEventHisto(int , TrackerDataImpl * ){
	// does nothing
}
