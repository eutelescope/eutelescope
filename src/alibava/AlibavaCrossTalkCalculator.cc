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
AlibavaBaseProcessor("AlibavaCrossTalkCalculator")
{
    
    // modify processor description
    _description =
    "AlibavaCrossTalkCalculator does whatever it does :) ";
    
    
    // first of register the input collection
    registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
                             "Input raw data collection name",
                             _inputCollectionName, string("rawdata") );
    
    // if needed one can change these to optional parameters
    
    registerProcessorParameter ("PedestalInputFile",
                                "The filename where the pedestal and noise values stored",
                                _pedestalFile , string("pedestal.slcio"));
    
    registerProcessorParameter ("CalibrationInputFile",
                                "The filename where the calibration values stored",
                                _calibrationFile , string("calibration.slcio"));
    
    
    // now the optional parameters
    registerProcessorParameter ("PedestalCollectionName",
                                "Pedestal collection name, better not to change",
                                _pedestalCollectionName, string ("pedestal"));
    
    registerProcessorParameter ("NoiseCollectionName",
                                "Noise collection name, better not to change",
                                _noiseCollectionName, string ("noise"));
    
    registerProcessorParameter ("ChargeCalibrationCollectionName",
                                "Charge calibration collection name, better not to change",
                                _chargeCalCollectionName, string ("chargeCal"));
    
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
    
    // if you want
    bookHistos();
    
    // set number of skipped events to zero (defined in AlibavaBaseProcessor)
    _numberOfSkippedEvents = 0;
    
}


void AlibavaCrossTalkCalculator::processEvent (LCEvent * anEvent) {
    
    AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
    
    if (_skipMaskedEvents && (alibavaEvent->isEventMasked()) ) {
        _numberOfSkippedEvents++;
        return;
    }
    string histoName;
    LCCollectionVec * collectionVec;
    unsigned int noOfClusters;
    try
    {
        collectionVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( getInputCollectionName() ) ) ;
        noOfClusters = collectionVec->getNumberOfElements();
        CellIDDecoder<TrackerDataImpl> clusterIDDecoder(ALIBAVA::ALIBAVACLUSTER_ENCODE);
        
        for ( size_t i = 0; i < noOfClusters; ++i )
        {
            TrackerDataImpl * trkdata = dynamic_cast< TrackerDataImpl * > ( collectionVec->getElementAt( i ) ) ;
            AlibavaCluster anAlibavaCluster(trkdata);
            int chipnum = anAlibavaCluster.getChipNum();
            int seedChan = anAlibavaCluster.getSeedChanNum();
            double signalPolarity = anAlibavaCluster.getSignalPolarity();
            
            int neighChan=0;
            double seedSignal=0.0, neighSignal=0.0;
            
            // using only clusters with cluster size = 2
            if (anAlibavaCluster.getClusterSize() == 2){
                for (int imember=0; imember<anAlibavaCluster.getClusterSize(); imember++) {
                    if (anAlibavaCluster.getChanNum(imember) == seedChan) {
                        seedSignal =anAlibavaCluster.getSignal(imember);
                    }else{
                        neighChan=anAlibavaCluster.getChanNum(imember);
                        neighSignal = anAlibavaCluster.getSignal(imember);
                    }
                }
            }
            double leftSignal =0 , rightSignal=0;
            if (neighChan<seedChan) {
                // neighbour is left
                leftSignal = neighSignal;
                rightSignal = seedSignal;
                histoName = string("b1_L_")+to_string(chipnum);
                if ( TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[histoName]) )
                    histo->Fill(leftSignal/seedSignal);
            }else{
                rightSignal = neighSignal;
                leftSignal = seedSignal;
                histoName = string("b1_R_")+to_string(chipnum);
                if ( TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[histoName]) )
                    histo->Fill(rightSignal/seedSignal);

            }
            
            double eta= leftSignal /(leftSignal+rightSignal);
            histoName = string("eta_")+to_string(chipnum);
            if ( TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[histoName]) )
                histo->Fill(eta);

        }
        
    } catch ( lcio::DataNotAvailableException ) {
        // do nothing again
//        streamlog_out( ERROR5 ) << "Collection ("<<getInputCollectionName()<<") not found! " << endl;
    }
    
}

void AlibavaCrossTalkCalculator::check (LCEvent * /* evt */ ) {
    // nothing to check here - could be used to fill check plots in reconstruction processor
}


void AlibavaCrossTalkCalculator::end() {
    
    if (_numberOfSkippedEvents > 0)
        streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents<<" events skipped since they are masked" << endl;
    streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
    
}

void AlibavaCrossTalkCalculator::fillHistos(TrackerDataImpl * /* trkdata */){
    // this is an example
    
    /*
     int ichip = getChipNum(trkdata);
     FloatVec datavec;
     datavec = trkdata->getChargeValues();
     for (size_t ichan=0; ichan<datavec.size();ichan++) {
     if(isMasked(ichip, ichan)) continue;
     
     string tempHistoName = getChanDataHistoName(ichip,ichan);
     if ( TH1D * histo = dynamic_cast<TH1D*> (_rootObjectMap[tempHistoName]) )
     histo->Fill(datavec[ichan]);
     }
     */
    
}


void AlibavaCrossTalkCalculator::bookHistos(){
    // this is an example
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
        
        histoName = string("b1_R_")+to_string(ichip);
        TH1D * b1_R = new TH1D (histoName.c_str(),"", 100, -2.0,2.0);
        _rootObjectMap.insert(make_pair(histoName, 	b1_R));
        b1_R->SetTitle("right signal / seed signal;b1_R;Number of Entries");
        
        histoName = string("eta_")+to_string(ichip);
        TH1D * heta = new TH1D (histoName.c_str(),"", 100, 0.0,1.0);
        _rootObjectMap.insert(make_pair(histoName, heta));
        heta->SetTitle("eta;eta;Number of Entries");
        
    }
    
    
    streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << endl;
}



















