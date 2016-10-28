/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


// alibava includes ".h"
#include "AlibavaBaseHistogramMaker.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"
#include "AlibavaPedNoiCalIOManager.h"


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
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCEventImpl.h>


// ROOT includes ".h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

// system includes <>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <memory>
#include <cstdlib>
#include <ctime>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;


AlibavaBaseHistogramMaker::AlibavaBaseHistogramMaker (std::string processorName) :
AlibavaBaseProcessor(processorName),
_histoXMLFileName("AlibavaHistoList.xml"),
_topTag("AlibavaHistoList"),
_tagToProcess( "myAlibavaBaseHistogramMaker"),
_plotPedestalAndNoise(false),
_plotEvents(),
_plotXPercentOfEvents(0),
_pedestalHistoName ("hPedestal"),
_noiseHistoName ("hNoise"),
_pedestalAndNoiseHistoDir("PedestalAndNoise"),
// Event plots
_eventHistoName("hEvent"),
_eventHistoDir("Events"),
_maskedEventsHistoName("hMaskedEvents"),
// Calibration
_calChargeHistoName("hCalibrationChargeValues"),
_delayHistoName("hDelayValues"),
_histoList(),
// Others
_multiplySignalby(1.0),
_totalNumberOfEvents(0),
_insertToXMLTitle(),
_replaceXMLValues()
{
	
	// modify processor description
	_description =
	"AlibavaBaseHistogramMaker is the base processor to produce histograms ";
	
	// Initialize random number generator
	srand(time(0));

}
void AlibavaBaseHistogramMaker::processHistoXMLFile(){
	
	streamlog_out ( MESSAGE1 )  << "Reading HistoXMLFile : "<< _histoXMLFileName << endl;
	
	string titleString;
	
	AIDAProcessor::tree(this)->cd(this->name());
	AIDAProcessor::tree(this)->mkdir(getInputCollectionName().c_str());
	AIDAProcessor::tree(this)->cd(getInputCollectionName().c_str());
		
	EVENT::IntVec chipSelection = getChipSelection();
	
	string signalMultipliedby;
	if (_multiplySignalby != 1.0){
		signalMultipliedby = to_string(_multiplySignalby);
		signalMultipliedby = string(" (")+signalMultipliedby+string(") ");
	}

	// first open XML document which contains list of histograms
	TiXmlDocument * histoXMLDoc = new TiXmlDocument;
	
	if ( histoXMLDoc->LoadFile( _histoXMLFileName ) ) {
		
		TiXmlHandle    handleHistoXMLDoc(histoXMLDoc);
		TiXmlHandle    handleTopTag(0);
		TiXmlElement * topTagElement;
		
		// go to topTag which is defined by _topTag
		topTagElement = handleHistoXMLDoc.FirstChildElement().Element();
		if ( !topTagElement ) {
			delete histoXMLDoc;
			throw ParseException( string( "AlibavaBaseHistogramMaker::bookHistos no <") + _topTag + string("> ... </")+ _topTag +string("> block found in ")  + _histoXMLFileName);
			
		} else {
			handleTopTag = TiXmlHandle(topTagElement);
		}
		
		// now go to the tag chosen to be processed
		TiXmlHandle handleTagToProcess = handleTopTag.FirstChild( _tagToProcess );
//		TiXmlNode * nHistosBlock = handleTagToProcess.ToNode();
		
		TiXmlNode * histosBlockNode = handleTagToProcess.ToNode();
		if ( !histosBlockNode ) {
			delete histoXMLDoc;
			throw ParseException( string( "AlibavaBaseHistogramMaker::bookHistos no <") + _tagToProcess +string("> ... </")+_tagToProcess+string("> block found in ")  + _histoXMLFileName);
		}
		
		int xBin, yBin;
		float xMin, xMax, yMin, yMax;
		string histoName, title, labelX, labelY, perEachChip ,histoType;
		string toSeparate(";");
		
		TH1D * histo1D;
		TH2D * histo2D;
		TProfile * profile;
		
		// now read histos
		TiXmlElement * histoNode = handleTagToProcess.FirstChild( "histo" ).Element();
		while (histoNode) {
			// read attributes
			histoName = histoNode->Attribute("name");
			histoType = histoNode->Attribute("type");
			perEachChip = histoNode->Attribute("perEachChip");
			histoNode->QueryIntAttribute("xBin", &xBin);
			histoNode->QueryFloatAttribute("xMin", &xMin);
			histoNode->QueryFloatAttribute("xMax", &xMax);
			histoNode->QueryIntAttribute("yBin", &yBin);
			histoNode->QueryFloatAttribute("yMin", &yMin);
			histoNode->QueryFloatAttribute("yMax", &yMax);
			
			title = histoNode->Attribute("title");
			labelX = histoNode->Attribute("labelX");
			labelY = histoNode->Attribute("labelY");
			
			
			// Here we will replace these values if any rule defined
			updateXMLVariables(histoName, &xBin, &xMin, &xMax, &yBin, &yMin, &yMax, &title, &labelX, &labelY);

			// if we should create this histogram for each chip
	
			if (perEachChip == string("true") || perEachChip == string("True")) {
				for (unsigned int i=0; i<chipSelection.size(); i++) {
					int ichip=chipSelection[i];
					string newHistoName = histoName;
					newHistoName =getHistoNameForChip(newHistoName,ichip);

					titleString = title + string(" (chip")+ to_string(ichip) + string(")");
					titleString = titleString + toSeparate + labelX + toSeparate + labelY;
					
					if (histoType == string("TH1D")){
						histo1D = new TH1D (newHistoName.c_str(),"",xBin, xMin, xMax);
						histo1D->SetTitle(titleString.c_str());
						_rootObjectMap.insert(make_pair(newHistoName, histo1D));
					}
					else if (histoType == string("TH2D")){
						histo2D = new TH2D (newHistoName.c_str(),"",xBin, xMin, xMax, yBin, yMin, yMax);
						histo2D->SetTitle(titleString.c_str());
						_rootObjectMap.insert(make_pair(newHistoName, histo2D));
					}
					else if (histoType == string("TProfile")){
						profile = new TProfile (newHistoName.c_str(),"",xBin, xMin, xMax, yMin, yMax);
						profile->SetTitle(titleString.c_str());
						_rootObjectMap.insert(make_pair(newHistoName, profile));
					}
					else
						streamlog_out (ERROR5) << "Histogram type of "<<newHistoName<<" is not valid!"<<endl;
				streamlog_out (DEBUG1) << "Histogram  "<< newHistoName <<" created!"<<endl;
				}
			}
			else { // just create one
				titleString = title + toSeparate + labelX + toSeparate + labelY;

				if (histoType == string("TH1D")){
					histo1D = new TH1D (histoName.c_str(),"",xBin, xMin, xMax);
					histo1D->SetTitle(titleString.c_str());
					_rootObjectMap.insert(make_pair(histoName, histo1D));
				}
				else if (histoType == string("TH2D")){
					histo2D = new TH2D (histoName.c_str(),"",xBin, xMin, xMax, yBin, yMin, yMax);
					histo2D->SetTitle(titleString.c_str());
					_rootObjectMap.insert(make_pair(histoName, histo2D));
				}
				else if (histoType == string("TProfile")){
					profile = new TProfile (histoName.c_str(),"",xBin, xMin, xMax, yMin, yMax);
					profile->SetTitle(titleString.c_str());
					_rootObjectMap.insert(make_pair(histoName, profile));
				}
				else
					streamlog_out (ERROR5) << "Histogram type of "<<histoName<<" is not valid!"<<endl;

				streamlog_out (DEBUG1) << "Histogram  "<< histoName <<" created!" <<endl;

			}
			
			// move to next element
			histoNode = histoNode->NextSiblingElement();
		}
		
	}
	if (checkListOfHistosCreatedByXMLFile()) {
		streamlog_out (MESSAGE1) << "All histograms needed are defined in "<<_histoXMLFileName<<endl;
	}
}

void AlibavaBaseHistogramMaker::addToHistoCheckList(string histoName){
	_histoList.push_back(histoName);
}
void AlibavaBaseHistogramMaker::addToHistoCheckList_PerChip(string histoName){
	IntVec chipSelection = getChipSelection();
	for (unsigned int i=0; i<chipSelection.size(); i++){
		int ichip=chipSelection[i];
		_histoList.push_back(getHistoNameForChip(histoName,ichip));
	}
}


bool AlibavaBaseHistogramMaker::checkListOfHistosCreatedByXMLFile(){
	// Checks if all the histograms needed by this processor is defined in _histoXMLFileName
	// Unfortunately histogram names are hard coded,  we are checking if all histo names exists in the _rootObjectMap
	
	// Now check if all histos exists
	for (unsigned int ihisto=0; ihisto<_histoList.size(); ihisto++) {
		string histoName = _histoList[ihisto];
		if (!doesRootObjectExists(histoName)) {
			streamlog_out (ERROR5) << "Histogram "<< histoName <<" doesn't exists"<<endl;
			streamlog_out (ERROR5) << "Define "<< _histoList[ihisto] <<" in "<< _histoXMLFileName <<endl;
			return false;
		} // end of if
	}

	// if we reach here, it means that all histograms exists, then return true
	return true;
}

///////////
// Event //
///////////

string AlibavaBaseHistogramMaker::getEventHistoName( int eventnum,  int ichip){
	stringstream s;
	s<< _eventHistoName<<"_"<< eventnum <<"_chip" << ichip;
	return s.str();
}

bool AlibavaBaseHistogramMaker::isEventToBePlotted(int eventnum){
		
	for (unsigned int iEvent=0; iEvent<_plotEvents.size(); iEvent++)
		if (_plotEvents[iEvent] == eventnum)
			return true;

	// generate a float number between 0-100 and
	float randnum = float(rand() % (100* 10000 ))/10000.0;
	if ( randnum <= _plotXPercentOfEvents )
		return true;
	// if event num is not in _plotEvents vector return false
	return false;
	
}

////////////////////////
// Pedestal and Noise //
////////////////////////

void AlibavaBaseHistogramMaker::plotPedestalAndNoise(){
	AIDAProcessor::tree(this)->cd(this->name());
	EVENT::IntVec chipSelection = getChipSelection();
	
	//the chipSelection should be in ascending order!
	//this is guaranteed with AlibavaConverter::checkIfChipSelectionIsValid()
	
	AIDAProcessor::tree(this)->mkdir(_pedestalAndNoiseHistoDir.c_str());
	AIDAProcessor::tree(this)->cd(_pedestalAndNoiseHistoDir.c_str());

	for (unsigned int i=0; i<chipSelection.size(); i++) {
		 int ichip=chipSelection[i];
		
		// pedestal
		if (_pedestalCollectionName != string(ALIBAVA::NOTSET) ) {
			TH1D * pedestalHisto = new TH1D (getPedestalHistoName(ichip).c_str(),"",ALIBAVA::NOOFCHANNELS, -0.5, ALIBAVA::NOOFCHANNELS-0.5);
			
			stringstream sp; //title string for pedestal histogram
			sp<< "Pedestal (chip "<<ichip<<");Channel Number;Pedestal (ADCs)";
			pedestalHisto->SetTitle((sp.str()).c_str());
			
			FloatVec pedVec = getPedestalOfChip(ichip);
			
			for (int ichan=0; ichan<ALIBAVA::NOOFCHANNELS; ichan++) {
				// if channel is masked, do not fill histo
				if (isMasked(ichip,ichan)) continue;
				else pedestalHisto->SetBinContent(ichan+1,pedVec[ichan]);
			}
			_rootObjectMap.insert(make_pair(getPedestalHistoName(ichip), pedestalHisto));
		}
		
		//noise
		if (_noiseCollectionName != string(ALIBAVA::NOTSET) ) {
			TH1D * noiseHisto = new TH1D (getNoiseHistoName(ichip).c_str(),"",ALIBAVA::NOOFCHANNELS, -0.5, ALIBAVA::NOOFCHANNELS-0.5);
			
			
			stringstream sn; //title string for noise histogram
			sn<< "Noise (chip "<<ichip<<");Channel Number;Pedestal (ADCs)";
			noiseHisto->SetTitle((sn.str()).c_str());
			
			FloatVec noiVec = getNoiseOfChip(ichip);
			for (int ichan=0; ichan<ALIBAVA::NOOFCHANNELS; ichan++) {
				// if channel is masked, do not fill histo
				if (isMasked(ichip,ichan)) continue;
				else noiseHisto->SetBinContent(ichan+1,noiVec[ichan]);
			}
			
			_rootObjectMap.insert(make_pair(getNoiseHistoName(ichip),noiseHisto));
			
		}
	} // end of loop over selected chips
}

string AlibavaBaseHistogramMaker::getPedestalHistoName( int ichip){
	return getHistoNameForChip(_pedestalHistoName, ichip);
}
string AlibavaBaseHistogramMaker::getNoiseHistoName( int ichip){
	return getHistoNameForChip(_noiseHistoName, ichip);
}

string AlibavaBaseHistogramMaker::getHistoNameForChip(string histoName, int ichip){
	stringstream s;
	s<< histoName <<"_chip" << ichip;
	return s.str();
}

/////////////////////////////
// To change XML variables //
/////////////////////////////

// changeXMLVariable method is used to replace xBin, xMin, xMax, yBin, yMin, yMax
// so variableName has to be one of these: xBin, xMin, xMax, yBin, yMin, yMax (case sensitive)
// for xBin and yBin newValue will changed to integer
void AlibavaBaseHistogramMaker::changeXMLVariable(string histoName, string variableName, float newValue){
	ReplaceXMLValues newValueToReplace;
	newValueToReplace._histoName = histoName;
	newValueToReplace._variableName = variableName;
	newValueToReplace._newValue = newValue;
	
	_replaceXMLValues.push_back(newValueToReplace);
}

// addToXMLTitle method is used to change title, labelX and labelY
// so titleName has to be one of these: title, labelX and labelY
// stringToAdd is the string that will added to the defined title
// whichSide defines which side of the title you want to add the "stringToAdd"
void AlibavaBaseHistogramMaker::addToXMLTitle(string histoName, string titleName, string whichSide, string stringToAdd){
	InsertToXMLTitle newTitleToInsert;
	newTitleToInsert._histoName = histoName;
	newTitleToInsert._titleName = titleName;
	newTitleToInsert._whichSide = whichSide;
	newTitleToInsert._stringToAdd = stringToAdd;
	
	_insertToXMLTitle.push_back(newTitleToInsert);
}
// Here we will update the variables read from XML file if needed
void AlibavaBaseHistogramMaker::updateXMLVariables(string histoName, int* xBin, float* xMin, float* xMax, int* yBin, float* yMin, float* yMax, string* title, string* labelX, string* labelY){
	
	// first loop over _insertToXMLTitle
	for (unsigned int i=0; i<_insertToXMLTitle.size(); i++) {
		InsertToXMLTitle newTitleToInsert = _insertToXMLTitle[i];
		// if there is a rule for this histogram name
		if (histoName == newTitleToInsert._histoName) {
			// find which title you want to change
			if (newTitleToInsert._titleName == string("title")) {
				if (newTitleToInsert._whichSide == string("left")) {
					*title = newTitleToInsert._stringToAdd + *title;
				}
				if (newTitleToInsert._whichSide == string("right")) {
					*title = *title + newTitleToInsert._stringToAdd;
				}
				if (newTitleToInsert._whichSide == string("all")) {
					*title = newTitleToInsert._stringToAdd;
				}
			}
			if (newTitleToInsert._titleName == string("labelX")) {
				if (newTitleToInsert._whichSide == string("left")) {
					*labelX = newTitleToInsert._stringToAdd + *labelX;
				}
				if (newTitleToInsert._whichSide == string("right")) {
					*labelX = *labelX + newTitleToInsert._stringToAdd;
				}
				if (newTitleToInsert._whichSide == string("all")) {
					*labelX = newTitleToInsert._stringToAdd;
				}
			}
			if (newTitleToInsert._titleName == string("labelY")) {
				if (newTitleToInsert._whichSide == string("left")) {
					*labelY = newTitleToInsert._stringToAdd + *labelY;
				}
				if (newTitleToInsert._whichSide == string("right")) {
					*labelY = *labelY + newTitleToInsert._stringToAdd;
				}
				if (newTitleToInsert._whichSide == string("all")) {
					*labelY = newTitleToInsert._stringToAdd;
				}
			}
		} // end of if histoName
	} // end of loop over _insertToXMLTitle
	
	
	// now loop over _replaceXMLValues
	for (unsigned int i=0; i<_replaceXMLValues.size(); i++) {
		ReplaceXMLValues newValueToReplace = _replaceXMLValues[i];

		// if there is a rule for this histogram name
		if (histoName == newValueToReplace._histoName) {
			// find which variable you want to change
			
			if (newValueToReplace._variableName == string("xBin")) {
				*xBin = int(newValueToReplace._newValue);
			}
			if (newValueToReplace._variableName == string("xMin")) {
				*xMin = newValueToReplace._newValue;
			}
			if (newValueToReplace._variableName == string("xMax")) {
				*xMax = newValueToReplace._newValue;
			}
			if (newValueToReplace._variableName == string("yBin")) {
				*yBin = int(newValueToReplace._newValue);
			}
			if (newValueToReplace._variableName == string("yMin")) {
				*yMin = newValueToReplace._newValue;
			}
			if (newValueToReplace._variableName == string("yMax")) {
				*yMax = newValueToReplace._newValue;
			}
			
		} // end of if histoName
	}// end of loop over _replaceXMLValues
}

// Helper Classes
ReplaceXMLValues::ReplaceXMLValues() :
_histoName(),
_variableName(),
_newValue(0)
{
}

InsertToXMLTitle::InsertToXMLTitle() :
_histoName(),
_titleName(),
_whichSide(),
_stringToAdd()
{
}






