/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@cern.ch
 */


// alibava includes ".h"
#include "AlibavaHeader.h"
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
#include <IMPL/LCEventImpl.h>


// ROOT includes ".h"
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include "TROOT.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TProfile.h"

// system includes <>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <memory>
#include <fstream>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;


AlibavaHeader::AlibavaHeader () :
AlibavaBaseProcessor("AlibavaHeader")
{

	_description = "AlibavaHeader does some analysis of the beetle header. Run once to get the header pedestals, then run again with the output headers in the steering file!";

	registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName", "Header name",  _inputCollectionName, string("Header") );

	registerProcessorParameter ("RawDataCollection", "The collection the raw data of channel 0 is stored in", _rawdatacollection , string("rawdata"));

	registerProcessorParameter ("Pedestal14", "The pedestal of header 14", _pedestal14, float (522.0));

	registerProcessorParameter ("Pedestal15", "The pedestal of header 15", _pedestal15, float (521.4));

	registerProcessorParameter ("PedestalChannel", "The pedestal of channel 0", _firstpedestal, float (516.7));

	registerOptionalParameter("OutputFile", "The filename to write coefficients to", _filterFileName, string("filtercoefficients.txt"));
}


void AlibavaHeader::init () {
	streamlog_out ( MESSAGE4 ) << "Running init" << endl;

	printParameters ();

	if (Global::parameters->isParameterSet(ALIBAVA::SKIPMASKEDEVENTS))
		_skipMaskedEvents = bool ( Global::parameters->getIntVal(ALIBAVA::SKIPMASKEDEVENTS) );
	else {
		streamlog_out ( MESSAGE4 ) << "The Global Parameter "<< ALIBAVA::SKIPMASKEDEVENTS <<" is not set! Masked events will be used!" << endl;
	}

}

void AlibavaHeader::processRunHeader (LCRunHeader * rdr)
{
	streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;

	// Add processor name to the runheader
	auto_ptr<AlibavaRunHeaderImpl> arunHeader ( new AlibavaRunHeaderImpl(rdr)) ;
	arunHeader->addProcessor(type());

	// get and set selected chips
	setChipSelection(arunHeader->getChipSelection());

	bookHistos();

	// set number of skipped events to zero (defined in AlibavaBaseProcessor)
	_numberOfSkippedEvents = 0;

}


void AlibavaHeader::processEvent (LCEvent * anEvent)
{

	if ( anEvent->getEventNumber() % 1000 == 0 )
		streamlog_out ( MESSAGE4 ) << "Looping events "<<anEvent->getEventNumber() << endl;

	AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);

	LCCollectionVec * collectionVec;
	LCCollectionVec * collectionVec2;
	try
	{
		// both chips have the same header -> element 0
		collectionVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( getInputCollectionName() ) ) ;
		collectionVec2 = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( _rawdatacollection ) ) ;
		TrackerDataImpl * trkdata = dynamic_cast< TrackerDataImpl * > ( collectionVec->getElementAt( 0 ) ) ;
		TrackerDataImpl * trkdata2 = dynamic_cast< TrackerDataImpl * > ( collectionVec2->getElementAt( 0 ) ) ;
		fillHistos(trkdata,trkdata2);
		correlateLastHeader(trkdata,trkdata2);

	} catch ( lcio::DataNotAvailableException ) {
		// do nothing again
		streamlog_out( ERROR5 ) << "Collection ("<<getInputCollectionName()<<") not found! " << endl;
	}
}


void AlibavaHeader::check (LCEvent * )
{

}


// end, calculate the pedestal of 16 headers and 1 channel, then finish
void AlibavaHeader::end()
{
	calculateHeaderNoise();
	streamlog_out ( MESSAGE4 ) << "You could run me again and set the pedestals to: header 14: " << _pedestal14 << " header 15: " << _pedestal15 << " channel 0: " << _firstpedestal << " !" << endl;
	streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;

	float x1=0;
	float x2=0;
	float x3=0;
	float x4=0;

	float y1=0;
	float y2=0;
	float y3=0;
	float y4=0;
	
	if ( TH2D * histo = dynamic_cast<TH2D*> (_rootObjectMap["LowLowHisto"]) )
	{
		x1 = histo->GetMean(1);
		y1 = histo->GetMean(2);
	}
	if ( TH2D * histo = dynamic_cast<TH2D*> (_rootObjectMap["LowHighHisto"]) )
	{
		x2 = histo->GetMean(1);
		y2 = histo->GetMean(2);
	}
	if ( TH2D * histo = dynamic_cast<TH2D*> (_rootObjectMap["HighLowHisto"]) )
	{
		x3 = histo->GetMean(1);
		y3 = histo->GetMean(2);
	}
	if ( TH2D * histo = dynamic_cast<TH2D*> (_rootObjectMap["HighHighHisto"]) )
	{
		x4 = histo->GetMean(1);
		y4 = histo->GetMean(2);
	}

	float c1=0;
	float c2=0;

	c1 = ( y1/x1 + y2/x2 + y3/x3 + y4/x4) / 4.0;
	// get a value for 14:
	float tempx = (fabs(x1)+fabs(x2)+fabs(x3)+fabs(x4))/4.0;
	// a delta y:
	float tempy = ( fabs(y2-y4) + fabs(y1-y3) )/2.0;
	c2 = tempy/tempx;

	streamlog_out ( MESSAGE4 ) << "Coefficients: c1: " << c1 << " c2: " << c2 << endl;

	ofstream filterFile;
	filterFile.open(_filterFileName.c_str());
	filterFile << c1 << endl;
	filterFile << c2 << endl;
	filterFile.close();

}


// fills histograms with the 16 headers and the first channel
void AlibavaHeader::fillHistos(TrackerDataImpl * trkdata, TrackerDataImpl * trkdata2)
{

	FloatVec datavec;
	datavec = trkdata->getChargeValues();
	FloatVec datavec2;
	datavec2 = trkdata2->getChargeValues();

	for (int i=0; i<16;i++)
	{
		string tempHistoName = getChanDataHistoName(i);
		if ( TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[tempHistoName]) )
		{
			histo->Fill(datavec[i]);
		}
	}
	if ( TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap["Channel1"]) )
	{
		histo->Fill(datavec2[0]);
	}
}


// returns histogram names for the header pedestal
string AlibavaHeader::getChanDataHistoName(unsigned int ichan)
{
	stringstream s;
	s<< "Header_" << ichan;
	return s.str();
}


// returns fit names for the header pedestal
string AlibavaHeader::getChanDataFitName(unsigned int ichan)
{
	stringstream s;
	s<< "Fit_" << ichan;
	return s.str();
}


// books histos
void AlibavaHeader::bookHistos()
{
	string tempHistoName,tempFitName;

	for ( int ichan=0; ichan<16; ichan++)
	{
		tempHistoName = getChanDataHistoName(ichan);
		tempFitName = getChanDataFitName(ichan);
		stringstream tempHistoTitle;
		tempHistoTitle<<tempHistoName<<";ADCs;NumberofEntries";

		TH1D * chanDataHisto = new TH1D (tempHistoName.c_str(),"",1000,0,1000);
		_rootObjectMap.insert(make_pair(tempHistoName, chanDataHisto));
		string tmp_string = tempHistoTitle.str();
		chanDataHisto->SetTitle(tmp_string.c_str());

		TF1 *chanDataFit = new TF1(tempFitName.c_str(),"gaus");
		_rootObjectMap.insert(make_pair(tempFitName, chanDataFit));

	}

	TH1D * chanDataHisto = new TH1D ("Channel1","",1000,0,1000);
	_rootObjectMap.insert(make_pair("Channel1", chanDataHisto));
	chanDataHisto->SetTitle("Channel1;ADCs;NumberofEntries");

	TF1 * chanDataFit = new TF1("Channel1Fit","gaus");
	_rootObjectMap.insert(make_pair("Channel1Fit", chanDataFit));

	TH2D * correlationHisto = new TH2D ("CorrelationToChan1","",1000,-500,500,1000,-500,500);
	_rootObjectMap.insert(make_pair("CorrelationToChan1", correlationHisto));
	correlationHisto->SetTitle("Correlation Last Header to Channel 1;Last Header [pedestal subtracted ADCs];First Channel [pedestal subtracted ADCs]");

	TProfile * lowprofile = new TProfile ("LowProfile","",100,-1000,1000,-1000,1000,"s");
	_rootObjectMap.insert(make_pair("LowProfile", lowprofile));
	lowprofile->SetTitle("Low Header Profile;Last Header [pedestal subtracted ADCs];First Channel [pedestal subtracted ADCs]");

	TProfile * highprofile = new TProfile ("HighProfile","",100,-1000,1000,-1000,1000,"s");
	_rootObjectMap.insert(make_pair("HighProfile", highprofile));
	highprofile->SetTitle("High Header Profile;Last Header [pedestal subtracted ADCs];First Channel [pedestal subtracted ADCs]");

	TH2D * correlation2Histo = new TH2D ("CorrelationLastHeaders","",1000,-500,500,1000,-500,500);
	_rootObjectMap.insert(make_pair("CorrelationLastHeaders", correlation2Histo));
	correlation2Histo->SetTitle("Correlation Last Header to Last-but-one Header;Last Header [pedestal subtracted ADCs];Last but one Header [pedestal subtracted ADCs]");

	TH2D * lowlowhisto = new TH2D ("LowLowHisto","",1000,-500,500,1000,-500,500);
	_rootObjectMap.insert(make_pair("LowLowHisto", lowlowhisto));
	lowlowhisto->SetTitle("Headers: 14Low 15Low;Last Header [pedestal subtracted ADCs];First Channel [pedestal subtracted ADCs]");

	TH2D * lowhighhisto = new TH2D ("LowHighHisto","",1000,-500,500,1000,-500,500);
	_rootObjectMap.insert(make_pair("LowHighHisto", lowhighhisto));
	lowhighhisto->SetTitle("Headers: 14Low 15High;Last Header [pedestal subtracted ADCs];First Channel [pedestal subtracted ADCs]");

	TH2D * highlowhisto = new TH2D ("HighLowHisto","",1000,-500,500,1000,-500,500);
	_rootObjectMap.insert(make_pair("HighLowHisto", highlowhisto));
	highlowhisto->SetTitle("Headers: 14High 15Low;Last Header [pedestal subtracted ADCs];First Channel [pedestal subtracted ADCs]");

	TH2D * highhighhisto = new TH2D ("HighHighHisto","",1000,-500,500,1000,-500,500);
	_rootObjectMap.insert(make_pair("HighHighHisto", highhighhisto));
	highhighhisto->SetTitle("Headers: 14High 15High;Last Header [pedestal subtracted ADCs];First Channel [pedestal subtracted ADCs]");

	streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << endl;
}


void AlibavaHeader::calculateHeaderNoise()
{
	string tempHistoName,tempFitName;
	TCanvas *cc = new TCanvas("cc","cc",800,600);

	EVENT::FloatVec pedestalVec,noiseVec;

	for (int ichan=0; ichan<16; ichan++)
	{
		double ped, noi;
		tempFitName = getChanDataFitName(ichan);
		tempHistoName = getChanDataHistoName(ichan);
		TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[tempHistoName]);
		TF1 * tempfit = dynamic_cast<TF1*> (_rootObjectMap[tempFitName]);
		histo->Fit(tempfit,"Q");
		ped = histo->GetMean(1);
		noi = histo->GetMean(11);
		pedestalVec.push_back(ped);
		noiseVec.push_back(noi);
		if (ichan == 14)
		{
			_pedestal14 = ped;
		}
		if (ichan == 15)
		{
			_pedestal15 = ped;
		}
	}
	TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap["Channel1"]);
	TF1 * tempfit = dynamic_cast<TF1*> (_rootObjectMap["Channel1Fit"]);
	histo->Fit(tempfit,"Q");
	_firstpedestal = histo->GetMean(1);

	delete cc;
}

// does some correlation plots
void AlibavaHeader::correlateLastHeader(TrackerDataImpl * trkdata, TrackerDataImpl * trkdata2)
{
	FloatVec datavec;
	datavec = trkdata->getChargeValues();
	FloatVec datavec2;
	datavec2 = trkdata2->getChargeValues();

	double header1 = datavec[14] - _pedestal14;
	double header = datavec[15] - _pedestal15;
	double channel = datavec2[0] - _firstpedestal;

	if ( TH2D * histo = dynamic_cast<TH2D*> (_rootObjectMap["CorrelationToChan1"]) )
	{
		histo->Fill(header,channel);
	}
	if ( TProfile* profile = dynamic_cast<TProfile*> (_rootObjectMap["LowProfile"]) )
	{
		if (header < 0)
		{
			profile->Fill(header,channel,1);
		}
	}
	if ( TProfile* profile = dynamic_cast<TProfile*> (_rootObjectMap["HighProfile"]) )
	{
		if (header > 0)
		{
			profile->Fill(header,channel,1);
		}
	}
	if ( TH2D * histo = dynamic_cast<TH2D*> (_rootObjectMap["CorrelationLastHeaders"]) )
	{
		histo->Fill(header,header1);
	}
	
	if (header1 < 0)
	{
		if (header < 0)
		{
			if ( TH2D * histo = dynamic_cast<TH2D*> (_rootObjectMap["LowLowHisto"]) )
			{
				histo->Fill(header,channel);
			}
		}
		if (header > 0)
		{
			if ( TH2D * histo = dynamic_cast<TH2D*> (_rootObjectMap["LowHighHisto"]) )
			{
				histo->Fill(header,channel);
			}
		}
	}
	if (header1 > 0)
	{
		if (header < 0)
		{
			if ( TH2D * histo = dynamic_cast<TH2D*> (_rootObjectMap["HighLowHisto"]) )
			{
				histo->Fill(header,channel);
			}
		}
		if (header > 0)
		{
			if ( TH2D * histo = dynamic_cast<TH2D*> (_rootObjectMap["HighHighHisto"]) )
			{
				histo->Fill(header,channel);
			}
		}
	}
}

