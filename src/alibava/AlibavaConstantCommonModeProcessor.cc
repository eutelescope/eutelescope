/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@cern.ch
 *
 *  modified by: Eda Yildirim eda.yildirim@cern.ch
 */


// alibava includes ".h"
#include "AlibavaConstantCommonModeProcessor.h"
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
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>

// ROOT includes ".h"
#include "TH1D.h"
#include "TH2D.h"
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

// TODO::put masking option

AlibavaConstantCommonModeProcessor::AlibavaConstantCommonModeProcessor () :
AlibavaBaseProcessor("AlibavaConstantCommonModeProcessor"),
_commonmodeCollectionName(ALIBAVA::NOTSET),
_commonmodeerrorCollectionName(ALIBAVA::NOTSET),
_Niteration(3),
_NoiseDeviation(2.5),
_commonmodeHistoName ("hcommonmode"),
_commonmodeerrorHistoName ("hcommonmodeerror"),
_commonmode(),
_commonmodeerror()
{
	
	// modify processor description
	_description =
	"AlibavaConstantCommonModeProcessor computes the common mode values of each chip and their errors";
	
	
	// first of register the input collection
	registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
									  "Input data collection name (should be pedestal subtracted!)",
									  _inputCollectionName, string("recodata") );

	// now the optional parameters
	registerOptionalParameter ("CommonModeCollectionName",
										"Common mode collection name, better not to change",
										_commonmodeCollectionName, string ("commonmode"));
	registerOptionalParameter ("CommonModeErrorCollectionName",
										"Common mode error collection name, better not to change",
										_commonmodeerrorCollectionName, string ("commonmodeerror"));
	
	registerOptionalParameter ("CommonModeErrorCalculationIteration",
										"The number of iteration that should be used in common mode calculation",
										_Niteration, int(3) );

	registerOptionalParameter ("NoiseDeviation",
										"The limit to the deviation of noise. The data exceeds this deviation will be considered as signal and not be included in common mode error calculation",
										_NoiseDeviation, float(2.5) );

}


void AlibavaConstantCommonModeProcessor::init () {
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
		streamlog_out ( MESSAGE4 ) << "The Global Parameter "<< ALIBAVA::CHANNELSTOBEUSED <<" is not set! All channels will be used!" << endl;
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
void AlibavaConstantCommonModeProcessor::processRunHeader (LCRunHeader * rdr) {
	streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;

	auto_ptr<AlibavaRunHeaderImpl> arunHeader ( new AlibavaRunHeaderImpl(rdr)) ;
	arunHeader->addProcessor(type());

	setChipSelection( arunHeader->getChipSelection() );
	setChannelsToBeUsed();
		
	bookHistos();
	
	// set number of skipped events to zero (defined in AlibavaBaseProcessor)
	_numberOfSkippedEvents = 0;

}

void AlibavaConstantCommonModeProcessor::processEvent (LCEvent * anEvent) {

	if ( anEvent->getEventNumber() % 1000 == 0 )
		streamlog_out ( MESSAGE4 ) << "Looping events "<<anEvent->getEventNumber() << endl;
	
	AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
	
	if (_skipMaskedEvents && (alibavaEvent->isEventMasked()) ) {
		_numberOfSkippedEvents++;
		return;
	}

	LCCollectionVec * collectionVec;
	LCCollectionVec* commonCollection = new LCCollectionVec(LCIO::TRACKERDATA);
	LCCollectionVec* commerrCollection = new LCCollectionVec(LCIO::TRACKERDATA);

	unsigned int noOfDetector;
	int chipnum;
	try
	{
		collectionVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( getInputCollectionName() ) ) ;
		noOfDetector = collectionVec->getNumberOfElements();
		
		streamlog_out( DEBUG0 ) << "===============================================================================" << endl;
		streamlog_out( DEBUG0 ) << "Will now loop chips..."<< endl;
		streamlog_out( DEBUG0 ) << "===============================================================================" << endl;
		

		CellIDEncoder<TrackerDataImpl> commonCol_CellIDEncode (ALIBAVA::ALIBAVADATA_ENCODE,commonCollection);
		CellIDEncoder<TrackerDataImpl> commerrCol_CellIDEncode (ALIBAVA::ALIBAVADATA_ENCODE,commerrCollection);

		for ( size_t i = 0; i < noOfDetector; ++i )
		{

			// Get the data
			TrackerDataImpl * trkdata = dynamic_cast< TrackerDataImpl * > ( collectionVec->getElementAt( i ) ) ;

			chipnum = getChipNum(trkdata);
			
			// Find Commonmode
			calculateConstantCommonMode(trkdata);

			TrackerDataImpl * commonData = new TrackerDataImpl();
			commonData->setChargeValues(getCommonModeVec());
			commonCol_CellIDEncode[ALIBAVA::ALIBAVADATA_ENCODE_CHIPNUM] = chipnum;
			commonCol_CellIDEncode.setCellID(commonData);
			commonCollection->push_back(commonData);

			TrackerDataImpl * commerrData = new TrackerDataImpl();
			commerrData->setChargeValues(getCommonModeErrorVec());
			commerrCol_CellIDEncode[ALIBAVA::ALIBAVADATA_ENCODE_CHIPNUM] = chipnum;
			commerrCol_CellIDEncode.setCellID(commerrData);
			commerrCollection->push_back(commerrData);
			
			// Fill Histos
			fillHistos(commonData,anEvent->getEventNumber());
  
		}
		alibavaEvent->addCollection(commonCollection, getCommonModeCollectionName());
		alibavaEvent->addCollection(commerrCollection, getCommonModeErrorCollectionName());

	}
	catch ( lcio::DataNotAvailableException )
	{
		// do nothing again
		streamlog_out( ERROR5 ) << "Collection ("<<getInputCollectionName()<<") not found! " << endl;
	}
	
}

void AlibavaConstantCommonModeProcessor::check (LCEvent * /* evt */ ) {
	// nothing to check here - could be used to fill check plots in reconstruction processor
}


void AlibavaConstantCommonModeProcessor::end() {

	if (_numberOfSkippedEvents > 0)
		streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents<<" events skipped since they are masked" << endl;
	streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
	
}

void AlibavaConstantCommonModeProcessor::calculateConstantCommonMode(TrackerDataImpl *trkdata){
  
  	EVENT::FloatVec commonmodeVec;
	EVENT::FloatVec commonmodeerrorVec;

	FloatVec datavec;
	datavec = trkdata->getChargeValues();

	
	int chipnum = getChipNum(trkdata);
	
	streamlog_out( DEBUG0 ) << "Chip " << chipnum << " of " << getNumberOfChips() << ", now iterating..." << endl;
	
	
	double sig=0, tmpdouble=0;
	double mean_signal=0;
	double sigma_mean_signal=0;
	
	
	for (int i=0; i<_Niteration; i++) {
		int nchan=0;
		double total_signal=0;
		double total_signal_square=0;
		
		// find mean value of signals
		for ( int ichan =0; ichan < int (datavec.size()); ichan++ ){
			if ( isMasked(chipnum, ichan) )
				continue;

			sig = datavec[ichan];
			
			// First iteration: take everything
			if (i==0) {
				total_signal += sig;
				total_signal_square += sig*sig;
				nchan++;
			}
			else{ // exclude outliers
				tmpdouble = fabs((sig - mean_signal)/sigma_mean_signal);
				if (tmpdouble<_NoiseDeviation) {
					total_signal += sig;
					total_signal_square += sig*sig;
					nchan++;
				}
			}
		} // end of loop over channels
		// here find the deviation from mean value
		// standard deviation = SQRT( E[x^2] - E[x]^2 )
		// where E denotes average value
		if (nchan>0) {
			mean_signal=total_signal/nchan;
			sigma_mean_signal=sqrt(total_signal_square/nchan - mean_signal*mean_signal);
		}

	} // end of iterations
	
	streamlog_out( DEBUG0 ) << "===============================================================================" << endl;
	streamlog_out( DEBUG0 ) << "Chip " << chipnum << " : CommonModeCorrection = " << mean_signal << ", CommonModeCorrectionError = " << sigma_mean_signal << endl;

	streamlog_out( DEBUG0 ) << "===============================================================================" << endl;
	
	// Now reloop over all channels and create the output vector. This will be the same for all channels
	for ( int ichan = 0; ichan < ALIBAVA::NOOFCHANNELS ; ichan++ )
	{
		commonmodeVec.push_back(mean_signal);
		commonmodeerrorVec.push_back(sigma_mean_signal);
	}
	
	setCommonModeVec(commonmodeVec);
	setCommonModeErrorVec(commonmodeerrorVec);
	
}


string AlibavaConstantCommonModeProcessor::getCommonCorrectionName(){
	string s;
	s = "Common Mode Correction Values";
	return s;
}


void AlibavaConstantCommonModeProcessor::fillHistos(TrackerDataImpl * trkdata, int event){

	// Fill the histograms with the corrected data
	FloatVec datavec;
	datavec = trkdata->getChargeValues();
	
	int chipnum = getChipNum(trkdata);
	
	for ( size_t ichan = 0 ; ichan < datavec.size() ; ichan++ )
	{
		if ( isMasked(chipnum,ichan) ) continue;
		
		string tempHistoName = getCommonCorrectionName();
		if ( TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[tempHistoName]) )
			histo->Fill(datavec[ichan]);
		
		if (TH2D * histo = dynamic_cast<TH2D*> (_rootObjectMap["Common Mode Correction Values over Events"]) )
			histo->Fill(event,datavec[ichan]);
		
	}
}


// getter and setter for _commonmodeCollectionName
void AlibavaConstantCommonModeProcessor::setCommonModeCollectionName(std::string commonmodeCollectionName){
	_commonmodeCollectionName = commonmodeCollectionName;
}
std::string AlibavaConstantCommonModeProcessor::getCommonModeCollectionName(){
	return _commonmodeCollectionName;
}

// getter and setter for _commonmodeerrorCollectionName
void AlibavaConstantCommonModeProcessor::setCommonModeErrorCollectionName(std::string commonmodeerrorCollectionName){
	_commonmodeerrorCollectionName = commonmodeerrorCollectionName;
}
std::string AlibavaConstantCommonModeProcessor::getCommonModeErrorCollectionName(){
	return _commonmodeerrorCollectionName;
}


// getter and setter for _commonmode
void AlibavaConstantCommonModeProcessor::setCommonModeVec(EVENT::FloatVec common){
	if (common.size()!=unsigned(ALIBAVA::NOOFCHANNELS)) {
		streamlog_out( ERROR5 ) <<"Size of common mode vector set doesn't match with number of channels in the data stream."<< endl;
	}
	_commonmode = common;
}
EVENT::FloatVec AlibavaConstantCommonModeProcessor::getCommonModeVec(){
	return _commonmode;
}

// getter and setter for _commonmodeerror
void AlibavaConstantCommonModeProcessor::setCommonModeErrorVec(EVENT::FloatVec commonerror){
	if (commonerror.size()!=unsigned(ALIBAVA::NOOFCHANNELS)) {
		streamlog_out( ERROR5 ) <<"Size of common mode error vector set doesn't match with number of channels in the data stream."<< endl;
	}
	_commonmodeerror = commonerror;
}
EVENT::FloatVec AlibavaConstantCommonModeProcessor::getCommonModeErrorVec(){
	return _commonmodeerror;
}





void AlibavaConstantCommonModeProcessor::bookHistos(){

	string tempHistoName;

	// a histogram showing the corrected signals
	tempHistoName = getCommonCorrectionName();
	stringstream tempHistoTitle;
	tempHistoTitle << tempHistoName << ";ADCs;NumberofEntries";

	TH1D * signalHisto =
	new TH1D (tempHistoName.c_str(),"",1000,-500,500);
	_rootObjectMap.insert(make_pair(tempHistoName, signalHisto));
	string tmp_string = tempHistoTitle.str();
	signalHisto->SetTitle(tmp_string.c_str());
	
	
	// a histogram showing the corrected signals over events
	stringstream tempHistoTitle2;
	tempHistoTitle2 << "Common Mode Correction Values over Events" << ";ADCs;NumberofEntries";

	TH2D * signalHisto2 = new TH2D ("Common Mode Correction Values over Events","",5000,0,500000,1000,-500,500);
	_rootObjectMap.insert(make_pair("Common Mode Correction Values over Events", signalHisto2));
	string tmp_string2 = tempHistoTitle2.str();
	signalHisto2->SetTitle(tmp_string2.c_str());

	streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << endl; 
}
