/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


// alibava includes ".h"
#include "AlibavaCrossTalkCalculator.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"
#include "AlibavaPedNoiCalIOManager.h"
#include "AlibavaCluster.h"


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
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"

// system includes <>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <memory>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;


AlibavaCrossTalkCalculator::AlibavaCrossTalkCalculator () :
AlibavaBaseProcessor("AlibavaCrossTalkCalculator"),
_recodataCollectionName("recodata"),
_crosstalkCollectionName("crosstalk"),
_maxNEvent(20000)
{
    
    // modify processor description
    _description =
    "AlibavaCrossTalkCalculator calculates the cross-talk correction values";
    
    
    // first of register the input collection
    registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
                             "Input collection name of alibava clusters",
                             _inputCollectionName, string("alibava_clusters") );
    
    registerProcessorParameter ("RecoDataCollectionName",
                                "The collection name of reconstructed data",
                                _recodataCollectionName, string ("recodata"));

    registerProcessorParameter ("DatabaseFile",
                                "The filename where the cross talk coefficients will be stored",
                                _databaseFile , string("database.slcio"));

    registerProcessorParameter ("CrossTalkCollectionName",
                                "The collection name of cross talk coefficients which will be stored in DatabaseFile",
                                _crosstalkCollectionName, string ("crosstalk"));
    
    registerProcessorParameter ("MaxNEvent",
                                "Max number of events to be used for cross talk calculation",
                                _maxNEvent, int (20000));
   
    
}


void AlibavaCrossTalkCalculator::init () {
    streamlog_out ( MESSAGE4 ) << "Running init" << endl;
    
    // this method is called only once even when the rewind is active
    // usually a good idea to
    printParameters ();
    
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
    
    
}
void AlibavaCrossTalkCalculator::processRunHeader (LCRunHeader * rdr) {
    streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;
    
    // Add processor name to the runheader
    auto_ptr<AlibavaRunHeaderImpl> arunHeader ( new AlibavaRunHeaderImpl(rdr)) ;
    arunHeader->addProcessor(type());
    
    // get and set selected chips
    setChipSelection(arunHeader->getChipSelection());
    
    // set channels to be used (if it is defined)
    setChannelsToBeUsed();
    
    // set pedestal and noise values
    setPedestals();
    
    // set number of skipped events to zero (defined in AlibavaBaseProcessor)
    _numberOfSkippedEvents = 0;
    
    
    bookHistos();
    fillHistos();
    calculateCrossTalk ();
    
    
}
void AlibavaCrossTalkCalculator::processEvent (LCEvent * /* anEvent */) {
 /*   AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
    
    if (_skipMaskedEvents && (alibavaEvent->isEventMasked()) ) {
        _numberOfSkippedEvents++;
        return;
    }
    LCCollectionVec * recodataColVec;
    
    LCCollectionVec* newDataColVec = new LCCollectionVec(LCIO::TRACKERDATA);
    CellIDEncoder<TrackerDataImpl> chipIDEncoder(ALIBAVA::ALIBAVADATA_ENCODE,newDataColVec);

    
    unsigned int noOfChips;
    try
    {
       recodataColVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( _recodataCollectionName ) ) ;
        
        noOfChips = recodataColVec->getNumberOfElements();
        
        for (unsigned int i=0; i< noOfChips; i++) {
            TrackerDataImpl * recodata = dynamic_cast< TrackerDataImpl * > ( recodataColVec->getElementAt( i ) ) ;
            int chipnum = getChipNum(recodata);
            FloatVec recodataVec = recodata->getChargeValues();
            
            FloatVec newRecoData;
            
            for (unsigned int ichan=0; ichan<recodataVec.size();ichan++) {
                newRecoData.push_back(recodataVec[ichan]);
            }
            
            float b1 = _crosstalkCoefficients[chipnum].at(0);
            float b2 = _crosstalkCoefficients[chipnum].at(1);

            for (unsigned int ichan=newRecoData.size()-1; ichan>1;ichan--) {
                // y[n] = x[n] + b1*x[n-1] + b2*x[n-2]
                // there will be no corrections for channel 0 and 1
                
                newRecoData[ichan]= newRecoData[ichan]+b1*newRecoData[ichan-1]+ b2*newRecoData[ichan-2];
                newRecoData[ichan-1]= newRecoData[ichan-1]-b1*newRecoData[ichan-1];
                newRecoData[ichan-2]= newRecoData[ichan-2]-b2*newRecoData[ichan-2];
            }

            
            TrackerDataImpl * newDataImpl = new TrackerDataImpl();
            newDataImpl->setChargeValues(newRecoData);
            chipIDEncoder[ALIBAVA::ALIBAVADATA_ENCODE_CHIPNUM] = chipnum;
            chipIDEncoder.setCellID(newDataImpl);
            newDataColVec->push_back(newDataImpl);
        }

        alibavaEvent->addCollection(newDataColVec, _outputDataCollectionName);
        
    } catch ( lcio::DataNotAvailableException ) {
        // do nothing again
        streamlog_out( DEBUG0 ) << "Collection ("<<getInputCollectionName()<<") not found in event number "<< alibavaEvent->getEventNumber()<< endl;
    }
  */
    
}


void AlibavaCrossTalkCalculator::fillHistos () {
    
    LCReader* lcReader = LCFactory::getInstance()->createLCReader();
    try {
        lcReader->open( Global::parameters->getStringVal("LCIOInputFiles") );
    } catch( IOException& e ){
        streamlog_out ( ERROR1 ) << "Can't open the file: " << e.what() << endl ;
    }
    
    AlibavaEventImpl* alibavaEvent;
    while( (alibavaEvent = static_cast<AlibavaEventImpl*> ( lcReader->readNextEvent() ) ) != 0 ){
        
        if (alibavaEvent->getEventNumber() > _maxNEvent && _maxNEvent!=0) {
            lcReader->close() ;
            return;
        }
        
        string histoName;
        LCCollectionVec * clusterColVec;
        LCCollectionVec * recodataColVec;
        unsigned int noOfClusters;
        unsigned int noOfChips;
        try
        {
            clusterColVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( getInputCollectionName() ) ) ;
            recodataColVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( _recodataCollectionName ) ) ;
            
            noOfClusters = clusterColVec->getNumberOfElements();
            CellIDDecoder<TrackerDataImpl> clusterIDDecoder(ALIBAVA::ALIBAVACLUSTER_ENCODE);
            
            noOfChips = recodataColVec->getNumberOfElements();
            vector<FloatVec> recodataVec;
            recodataVec.resize(noOfChips);
            
            for (unsigned int i=0; i< noOfChips; i++) {
                TrackerDataImpl * recodata = dynamic_cast< TrackerDataImpl * > ( recodataColVec->getElementAt( i ) ) ;
                int chipnum = getChipNum(recodata);
                recodataVec[chipnum] = recodata->getChargeValues();
            }
            
            for ( size_t i = 0; i < noOfClusters; ++i )
            {
                TrackerDataImpl * clusterdata = dynamic_cast< TrackerDataImpl * > ( clusterColVec->getElementAt( i ) ) ;
                AlibavaCluster anAlibavaCluster(clusterdata);
                int chipnum = anAlibavaCluster.getChipNum();
                int seedChan = anAlibavaCluster.getSeedChanNum();
                
                double seedSignal=0, leftSignal=0, rightSignal=0, leftleftSignal=0, rightrightSignal=0;
                
                seedSignal = recodataVec[chipnum].at(seedChan);
                
                if (seedChan-1>0 && !isMasked(chipnum,seedChan-1))
                    leftSignal = recodataVec[chipnum].at(seedChan-1);
                
                if (seedChan-2>0 && !isMasked(chipnum,seedChan-2))
                    leftleftSignal = recodataVec[chipnum].at(seedChan-2);
                
                if (seedChan+1<int(recodataVec[chipnum].size()) && !isMasked(chipnum,seedChan+1))
                    rightSignal = recodataVec[chipnum].at(seedChan+1);
                
                if (seedChan+2<int(recodataVec[chipnum].size()) && !isMasked(chipnum,seedChan+2))
                    rightrightSignal = recodataVec[chipnum].at(seedChan+2);
                
                histoName = string("b1_L_")+to_string(chipnum);
                if ( TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[histoName]) )
                    histo->Fill(leftSignal/seedSignal);
                
                histoName = string("b1_R_")+to_string(chipnum);
                if ( TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[histoName]) )
                    histo->Fill(rightSignal/seedSignal);
                
                histoName = string("b2_L_")+to_string(chipnum);
                if ( TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[histoName]) )
                    histo->Fill(leftleftSignal/seedSignal);
                
                histoName = string("b2_R_")+to_string(chipnum);
                if ( TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[histoName]) )
                    histo->Fill(rightrightSignal/seedSignal);
                
                
            }
            
        } catch ( lcio::DataNotAvailableException ) {
            // do nothing again
            streamlog_out( DEBUG0 ) << "Collection ("<<getInputCollectionName()<<") not found in event number "<< alibavaEvent->getEventNumber()<< endl;
        }
        
    }
    lcReader->close() ;
    
}

void AlibavaCrossTalkCalculator::check (LCEvent * /* evt */ ) {
    // nothing to check here - could be used to fill check plots in reconstruction processor
}


void AlibavaCrossTalkCalculator::end() {
}
void AlibavaCrossTalkCalculator::calculateCrossTalk() {
   	TCanvas *cc = new TCanvas("cc","cc",800,600);
 
    
    int nChips = getNumberOfChips();
    string histoName;
    
    for (int ichip=0; ichip<nChips; ichip++) {
        histoName = string("b1_L_")+to_string(ichip);
        TH1D * hb1_L = dynamic_cast<TH1D*> (_rootObjectMap[histoName]);
        histoName+=string("fit");
        TF1 *hb1_L_Fit=dynamic_cast<TF1*> (_rootObjectMap[histoName]);
	hb1_L->Fit(hb1_L_Fit,"QR");
        double b1_L = hb1_L_Fit->GetParameter(1);
        
        histoName = string("b1_R_")+to_string(ichip);
        TH1D * hb1_R = dynamic_cast<TH1D*> (_rootObjectMap[histoName]);
        histoName+=string("fit");
        TF1 *hb1_R_Fit=dynamic_cast<TF1*> (_rootObjectMap[histoName]);
        hb1_R->Fit(hb1_R_Fit,"QR");
        double b1_R = hb1_R_Fit->GetParameter(1);
        
        histoName = string("b2_L_")+to_string(ichip);
        TH1D * hb2_L = dynamic_cast<TH1D*> (_rootObjectMap[histoName]);
        histoName+=string("fit");
        TF1 *hb2_L_Fit=dynamic_cast<TF1*> (_rootObjectMap[histoName]);
        hb2_L->Fit(hb2_L_Fit,"QR");
        double b2_L = hb2_L_Fit->GetParameter(1);
        
        histoName = string("b2_R_")+to_string(ichip);
        TH1D * hb2_R = dynamic_cast<TH1D*> (_rootObjectMap[histoName]);
        histoName+=string("fit");
        TF1 *hb2_R_Fit=dynamic_cast<TF1*> (_rootObjectMap[histoName]);
        hb2_R->Fit(hb2_R_Fit,"QR");
        double b2_R = hb2_R_Fit->GetParameter(1);
        
        double b1 = b1_L - b1_R;
        double b2 = b2_L - b2_R;
        // we expect b1 and b2 to be negative
        // so we save them multiplied by -1
        b1 = -b1;
        b2 = -b2;
        
        streamlog_out ( MESSAGE5 ) << "Chip "<<ichip<<" b1 "<<b1<<" b2 "<<b2 << endl;
        _crosstalkCoefficients[ichip].push_back(b1);
        _crosstalkCoefficients[ichip].push_back(b2);
        
        AlibavaPedNoiCalIOManager man;
        man.addToFile(_databaseFile,_crosstalkCollectionName, ichip, _crosstalkCoefficients[ichip]);
        
    }
    if (_numberOfSkippedEvents > 0)
        streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents<<" events skipped since they are masked" << endl;
    streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
delete cc;    
}

void AlibavaCrossTalkCalculator::fillHistos(TrackerDataImpl * /* trkdata */){
    
}


void AlibavaCrossTalkCalculator::bookHistos(){
    AIDAProcessor::tree(this)->cd(this->name());
    AIDAProcessor::tree(this)->mkdir(getInputCollectionName().c_str());
    AIDAProcessor::tree(this)->cd(getInputCollectionName().c_str());
    
    int nChips = getNumberOfChips();
    string histoName;
    
    for (int ichip=0; ichip<nChips; ichip++) {
        histoName = string("b1_L_")+to_string(ichip);
        TH1D * b1_L = new TH1D (histoName.c_str(),"", 100, -2.0,2.0);
        _rootObjectMap.insert(make_pair(histoName, 	b1_L));
        b1_L->SetTitle("left signal / seed signal;b1_L;Number of Entries");
        histoName+=string("fit");
        TF1 *b1_L_Fit= new TF1(histoName.c_str(),"gaus");
        b1_L_Fit->SetRange(-0.2,0.2);
	_rootObjectMap.insert(make_pair(histoName, b1_L_Fit));
        
        histoName = string("b1_R_")+to_string(ichip);
        TH1D * b1_R = new TH1D (histoName.c_str(),"", 100, -2.0,2.0);
        _rootObjectMap.insert(make_pair(histoName, 	b1_R));
        b1_R->SetTitle("right signal / seed signal;b1_R;Number of Entries");
        histoName+=string("fit");
        TF1 *b1_R_Fit= new TF1(histoName.c_str(),"gaus");
        b1_R_Fit->SetRange(-0.2,0.2);
        _rootObjectMap.insert(make_pair(histoName, b1_R_Fit));
        
        histoName = string("b2_L_")+to_string(ichip);
        TH1D * b2_L = new TH1D (histoName.c_str(),"", 100, -2.0,2.0);
        _rootObjectMap.insert(make_pair(histoName, 	b2_L));
        b2_L->SetTitle("left left signal / seed signal;b2_L;Number of Entries");
        histoName+=string("fit");
        TF1 *b2_L_Fit= new TF1(histoName.c_str(),"gaus");
        b2_L_Fit->SetRange(-0.2,0.2);
        _rootObjectMap.insert(make_pair(histoName, b2_L_Fit));
        
        histoName = string("b2_R_")+to_string(ichip);
        TH1D * b2_R = new TH1D (histoName.c_str(),"", 100, -2.0,2.0);
        _rootObjectMap.insert(make_pair(histoName, 	b2_R));
        b2_R->SetTitle("left left signal / seed signal;b2_R;Number of Entries");
        histoName+=string("fit");
        TF1 *b2_R_Fit= new TF1(histoName.c_str(),"gaus");
        b2_R_Fit->SetRange(-0.2,0.2);
        _rootObjectMap.insert(make_pair(histoName, b2_R_Fit));
        
        histoName = string("eta_")+to_string(ichip);
        TH1D * heta = new TH1D (histoName.c_str(),"", 100, 0.0,1.0);
        _rootObjectMap.insert(make_pair(histoName, heta));
        heta->SetTitle("eta;eta;Number of Entries");
        
    }
    
    
    streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << endl;
}


