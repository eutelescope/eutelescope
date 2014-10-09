/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


// alibava includes ".h"
#include "AlibavaPedestalNoiseProcessor.h"
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
#include <IMPL/TrackerDataImpl.h>

// ROOT includes ".h"
#include "TH1D.h"
#include "TF1.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TSystem.h"

// system includes <>
#include <string>
#include <iostream>
#include <sstream>
#include <memory>


using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;


AlibavaPedestalNoiseProcessor::AlibavaPedestalNoiseProcessor () :
AlibavaBaseProcessor("AlibavaPedestalNoiseProcessor"),
_pedestalHistoName ("hpedestal"),
_noiseHistoName ("hnoise"),
_temperatureHistoName("htemperature"),
_chanDataHistoName ("Data_chan"),
_chanDataFitName ("Fit_chan")
{
	
	// modify processor description
	_description =
	"AlibavaPedestalNoiseProcessor computes the pedestal and noise values of each channel";
	
	
	// first of register the input collection
	registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
									  "Input raw data collection name",
									  _inputCollectionName, string("rawdata") );

	registerProcessorParameter ("PedestalOutputFile",
										 "The filename to store the pedestal and noise values",
										 _pedestalFile , string("outputped.slcio"));

	// now the optional parameters
	registerOptionalParameter ("PedestalCollectionName",
										"Pedestal collection name, better not to change",
										_pedestalCollectionName, string ("pedestal"));
	
	registerOptionalParameter ("NoiseCollectionName",
										"Noise collection name, better not to change",
										_noiseCollectionName, string ("noise"));

}


void AlibavaPedestalNoiseProcessor::init () {
	streamlog_out ( MESSAGE4 ) << "Running init" << endl;


	/* To set of channels to be used
	 ex.The format should be like $ChipNumber:StartChannel-EndChannel$
	 ex. $0:5-20$ $0:30-100$ $1:50-70$
	 means from chip 0 channels between 5-20 and 30-100, from chip 1 channels between 50-70 will be used (all numbers included). the rest will be masked and not used
	 Note that the numbers should be in ascending order and there should be no space between two $ character
	 */
	if (Global::parameters->isParameterSet(ALIBAVA::CHANNELSTOBEUSED))
		Global::parameters->getStringVals(ALIBAVA::CHANNELSTOBEUSED,_channelsToBeUsed);
	else {
		streamlog_out ( MESSAGE4 ) << "The Global Parameter "<< ALIBAVA::CHANNELSTOBEUSED <<" is not set!" << endl;
	}
	
	/* To choose if processor should skip masked events
	 ex. Set the value to 0 for false, to 1 for true
	 */
	if (Global::parameters->isParameterSet(ALIBAVA::SKIPMASKEDEVENTS))
		_skipMaskedEvents = bool ( Global::parameters->getIntVal(ALIBAVA::SKIPMASKEDEVENTS) );
	else {
		streamlog_out ( MESSAGE4 ) << "The Global Parameter "<< ALIBAVA::SKIPMASKEDEVENTS <<" is not set! Masked events will be used!" << endl;
	}
	// this method is called only once even when the rewind is active
	// usually a good idea to
	printParameters ();

}
void AlibavaPedestalNoiseProcessor::processRunHeader (LCRunHeader * rdr) {
	streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;

	auto_ptr<AlibavaRunHeaderImpl> arunHeader ( new AlibavaRunHeaderImpl(rdr)) ;
	arunHeader->addProcessor(type());
	setChipSelection( arunHeader->getChipSelection() );
	setChannelsToBeUsed();
	
	AlibavaPedNoiCalIOManager man;
	man.createFile(_pedestalFile, arunHeader->lcRunHeader());
	
	bookHistos();

	// set number of skipped events to zero (defined in AlibavaBaseProcessor)
	_numberOfSkippedEvents = 0;

}


void AlibavaPedestalNoiseProcessor::processEvent (LCEvent * anEvent) {

	if ( anEvent->getEventNumber() % 1000 == 0 )
		streamlog_out ( MESSAGE4 ) << "Looping events "<<anEvent->getEventNumber() << endl;
	
	AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
	
	if (_skipMaskedEvents && (alibavaEvent->isEventMasked()) ) {
		_numberOfSkippedEvents++;
		return;
	}

	
	LCCollectionVec * collectionVec;
	unsigned int noOfDetector; // here it also stands for number or chips stored
	try
	{
		collectionVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( getInputCollectionName() ) ) ;
		noOfDetector = collectionVec->getNumberOfElements();
		
		// fill temperature histogram
		TH1D * temperatureHisto = dynamic_cast<TH1D*> (_rootObjectMap[_temperatureHistoName]);
		temperatureHisto->Fill(alibavaEvent->getEventTemp());
		
		
		for ( size_t i = 0; i < noOfDetector; ++i )
		{
			
			TrackerDataImpl * trkdata = dynamic_cast< TrackerDataImpl * > ( collectionVec->getElementAt( i ) ) ;
			fillHistos(trkdata);
			
		}
		
	} catch ( lcio::DataNotAvailableException ) {
		// do nothing again
		streamlog_out( ERROR5 ) << "Collection ("<<getInputCollectionName()<<") not found! " << endl;
	}
	
}

void AlibavaPedestalNoiseProcessor::check (LCEvent * /* evt */ ) {
	// nothing to check here - could be used to fill check plots in reconstruction processor
}


void AlibavaPedestalNoiseProcessor::end() {
	calculatePedestalNoise();

	if (_numberOfSkippedEvents > 0)
		streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents<<" events skipped since they are masked" << endl;
	streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
}

void AlibavaPedestalNoiseProcessor::calculatePedestalNoise(){
	string tempHistoName,tempFitName;
	TCanvas *cc = new TCanvas("cc","cc",800,600);
	
	EVENT::IntVec chipSelection = getChipSelection();
	for (unsigned int i=0; i<chipSelection.size(); i++) {
		unsigned int ichip=chipSelection[i];
		
		TH1D * hped = dynamic_cast<TH1D*> (_rootObjectMap[getPedestalHistoName(ichip)]);
		TH1D * hnoi = dynamic_cast<TH1D*> (_rootObjectMap[getNoiseHistoName(ichip)]);
		EVENT::FloatVec pedestalVec,noiseVec;
		
		
		for (int ichan=0; ichan<ALIBAVA::NOOFCHANNELS; ichan++) {
			double ped, noi;
			// if channel is masked, set pedestal and noise to 0
			if (isMasked(ichip,ichan)){
				ped=0; noi=0;
			}
			else {
				tempFitName = getChanDataFitName(ichip, ichan);
				tempHistoName = getChanDataHistoName(ichip, ichan);
				TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[tempHistoName]);
				TF1 * tempfit = dynamic_cast<TF1*> (_rootObjectMap[tempFitName]);
				histo->Fit(tempfit,"Q");
				ped = tempfit->GetParameter(1);
				noi = tempfit->GetParameter(2);
				hped->SetBinContent(ichan+1,ped);
				hnoi->SetBinContent(ichan+1,noi);
			}
			pedestalVec.push_back(ped);
			noiseVec.push_back(noi);
		}
		
		AlibavaPedNoiCalIOManager man;
		man.addToFile(_pedestalFile,_pedestalCollectionName, ichip, pedestalVec);
		man.addToFile(_pedestalFile,_noiseCollectionName, ichip, noiseVec);
	}
	delete cc;
}

string AlibavaPedestalNoiseProcessor::getChanDataHistoName(unsigned int ichip, unsigned int ichan){
	stringstream s;
	s<< _chanDataHistoName<<"_chip"<<ichip<<"_chan" << ichan;
	return s.str();
}
string AlibavaPedestalNoiseProcessor::getChanDataFitName(unsigned int ichip, unsigned int ichan){
	stringstream s;
	s<< _chanDataFitName<<"_chip"<<ichip<<"_chan"  << ichan;
	return s.str();
}
string AlibavaPedestalNoiseProcessor::getPedestalHistoName(unsigned int ichip){
	stringstream s;
	s<< _pedestalHistoName<<"_chip" << ichip;
	return s.str();
}
string AlibavaPedestalNoiseProcessor::getNoiseHistoName(unsigned int ichip){
	stringstream s;
	s<< _noiseHistoName<<"_chip" << ichip;
	return s.str();
}


void AlibavaPedestalNoiseProcessor::fillHistos(TrackerDataImpl * trkdata){
	
	FloatVec datavec;
	datavec = trkdata->getChargeValues();
	
	int chipnum = getChipNum(trkdata);
	
	for (size_t ichan=0; ichan<datavec.size();ichan++) {
		if(isMasked(chipnum, ichan)) continue;
		
		string tempHistoName = getChanDataHistoName(chipnum, ichan);
		if ( TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[tempHistoName]) )
			histo->Fill(datavec[ichan]);
	}

}


void AlibavaPedestalNoiseProcessor::bookHistos(){
	AIDAProcessor::tree(this)->cd(this->name());
	EVENT::IntVec chipSelection = getChipSelection();

	//the chipSelection should be in ascending order!
	//this is guaranteed with AlibavaConverter::checkIfChipSelectionIsValid()
	
	// temperature of event
	TH1D * temperatureHisto = new TH1D(_temperatureHistoName.c_str(),"Temperature",1000,-50,50);
	_rootObjectMap.insert(make_pair(_temperatureHistoName,temperatureHisto));
	
	for (unsigned int i=0; i<chipSelection.size(); i++) {
		unsigned int ichip=chipSelection[i];
		
		TH1D * pedestalHisto = new TH1D (getPedestalHistoName(ichip).c_str(),"",ALIBAVA::NOOFCHANNELS, -0.5, ALIBAVA::NOOFCHANNELS-0.5);
		_rootObjectMap.insert(make_pair(getPedestalHistoName(ichip), pedestalHisto));

		stringstream sp; //title string for pedestal histogram
		sp<< "Pedestal (chip "<<ichip<<");Channel Number;Pedestal (ADCs)";
		pedestalHisto->SetTitle((sp.str()).c_str());
		
		TH1D * noiseHisto = new TH1D (getNoiseHistoName(ichip).c_str(),"",ALIBAVA::NOOFCHANNELS, -0.5, ALIBAVA::NOOFCHANNELS-0.5);

		_rootObjectMap.insert(make_pair(getNoiseHistoName(ichip),noiseHisto));
		
		stringstream sn; //title string for noise histogram
		sn<< "Noise (chip "<<ichip<<");Channel Number;Pedestal (ADCs)";
		noiseHisto->SetTitle((sn.str()).c_str());
		
	}


	AIDAProcessor::tree(this)->mkdir(getInputCollectionName().c_str());
	AIDAProcessor::tree(this)->cd(getInputCollectionName().c_str());
		
	// here are the histograms used to calculate pedestal and noise for each channel
	string tempHistoName,tempFitName;
	for (unsigned int i=0; i<chipSelection.size(); i++) {
		unsigned int ichip=chipSelection[i];

		for ( int ichan=0; ichan<ALIBAVA::NOOFCHANNELS; ichan++) {
			if (isMasked(ichip,ichan)) continue;
			
			
			tempHistoName = getChanDataHistoName(ichip,ichan);
			tempFitName = getChanDataFitName(ichip,ichan);
			stringstream tempHistoTitle;
			tempHistoTitle<<tempHistoName<<";ADCs;NumberofEntries";
			
			TH1D * chanDataHisto =
			new TH1D (tempHistoName.c_str(),"",1000,0,1000);
			_rootObjectMap.insert(make_pair(tempHistoName, chanDataHisto));
			string tmp_string = tempHistoTitle.str();
			chanDataHisto->SetTitle(tmp_string.c_str());
			
			TF1 *chanDataFit = new TF1(tempFitName.c_str(),"gaus");
			_rootObjectMap.insert(make_pair(tempFitName, chanDataFit));
			
			
		}
	}

	streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << endl;
}


