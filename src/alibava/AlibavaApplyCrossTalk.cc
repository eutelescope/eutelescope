/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


// alibava includes ".h"
#include "AlibavaApplyCrossTalk.h"
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


AlibavaApplyCrossTalk::AlibavaApplyCrossTalk () :
AlibavaBaseProcessor("AlibavaApplyCrossTalk"),
//_inputCollectionName("recodata"),
//_outputCollectionName("recodata_xtalk"),
_crosstalkCollectionName(),
_databaseFile("database.slcio"),
_stopCorrectionValue(0)
{
    
    // modify processor description
    _description =
    "AlibavaApplyCrossTalk does whatever it does :) ";
    
    
    // first of register the input collection
    registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
                             "Input collection name of reconstructed data",
                             _inputCollectionName, string("recodata") );
    
    registerOutputCollection (LCIO::TRACKERDATA, "OutputCollectionName",
                             "Output collection name of cross talk corrected reconstructed data",
                             _outputCollectionName, string("recodata_xtalk") );
    
    registerProcessorParameter ("DatabaseFile",
                                "The filename where the cross talk coefficients are stored",
                                _databaseFile , string("database.slcio"));
    
    registerProcessorParameter ("CrossTalkCollectionName",
                                "The collection name of cross talk coefficients which will are stored in DatabaseFile",
                                _crosstalkCollectionName, EVENT::StringVec ());
    
    registerProcessorParameter ("StopCorrectionValue",
                                "The cross talk correction will stop if the next channel crosstalk (b1) is smaller than StopCorrectionValue",
                                _stopCorrectionValue, float (0.0));
    
}


void AlibavaApplyCrossTalk::init () {
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
    
    for (int i=0; i<ALIBAVA::NOOFCHIPS; i++) {
        _crosstalkCoefficients[i].clear();
    }
    
    
}
void AlibavaApplyCrossTalk::processRunHeader (LCRunHeader * rdr) {
    streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;
    
    // Add processor name to the runheader
    auto_ptr<AlibavaRunHeaderImpl> arunHeader ( new AlibavaRunHeaderImpl(rdr)) ;
    arunHeader->addProcessor(type());
    
    // get and set selected chips
    setChipSelection(arunHeader->getChipSelection());
    
    // set channels to be used (if it is defined)
    setChannelsToBeUsed();
    
    // set pedestal and noise values
    // setPedestals();
    
    // set number of skipped events to zero (defined in AlibavaBaseProcessor)
    _numberOfSkippedEvents = 0;
    
}
void AlibavaApplyCrossTalk::getCrossTalkCoefficients(){
    AlibavaPedNoiCalIOManager man;
    EVENT::IntVec used_chips = getChipSelection();

    for (unsigned int icor=0; icor<_crosstalkCollectionName.size(); icor++) {
        
        for (unsigned int ichip=0; ichip<used_chips.size(); ichip++) {
            int chipnum = used_chips[ichip];
            EVENT::FloatVec vec_float = man.getPedNoiCalForChip(_databaseFile,_crosstalkCollectionName[icor], chipnum);
            
            // check if vector size is 2
            if (vec_float.size() != 2) {
                streamlog_out(ERROR5)<<"Not enough values in the collection "<<_crosstalkCollectionName[icor]<< " in the database file "<<_databaseFile<<endl;
                streamlog_out(ERROR5)<<"Either b1 or b2 is missing!"<<endl;
            }
            else if(vec_float[0] > _stopCorrectionValue){
                _crosstalkCoefficients[chipnum].insert( _crosstalkCoefficients[chipnum].end(), vec_float.begin(), vec_float.end() );
            }
            else {
                streamlog_out(MESSAGE4)<<"The correction values: b1="<<vec_float[0]<<" b2="<<vec_float[1]<<" for chip"<<chipnum<<" will not be used since b1 is smaller than StopCorrectionValue. These correction values comes from "<<_crosstalkCollectionName[icor]<<" collection"<<endl;
            }
        }
            
    }
    
    for (unsigned int ichip=0; ichip<used_chips.size(); ichip++) {
        int chipnum = used_chips[ichip];
        
        streamlog_out(MESSAGE4)<< "Cross talk correction coefficients for chip"<< chipnum <<" :"<<endl;
        
        
        for (unsigned int icor=0; icor<_crosstalkCoefficients[chipnum].size() ; icor=icor+2) {
            streamlog_out(MESSAGE4)<< "     b1="<<_crosstalkCoefficients[chipnum].at(icor)<<"  b2="<< _crosstalkCoefficients[chipnum].at(icor+1)<<endl;
        }
    }

}

void AlibavaApplyCrossTalk::processEvent (LCEvent * anEvent) {
    AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
    
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
        recodataColVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( _inputCollectionName ) ) ;
        
        noOfChips = recodataColVec->getNumberOfElements();
        
        for (unsigned int i=0; i< noOfChips; i++) {
            TrackerDataImpl * recodata = dynamic_cast< TrackerDataImpl * > ( recodataColVec->getElementAt( i ) ) ;
            int chipnum = getChipNum(recodata);
            FloatVec recodataVec = recodata->getChargeValues();
            
            FloatVec newRecoData;
            
            for (unsigned int ichan=0; ichan<recodataVec.size();ichan++) {
                newRecoData.push_back(recodataVec[ichan]);
            }
            
            for (unsigned int icor=0; icor<_crosstalkCoefficients[chipnum].size(); icor=icor+2) {
                
                float b1 = _crosstalkCoefficients[chipnum].at(icor);
                float b2 = _crosstalkCoefficients[chipnum].at(icor+1);
                
                for (unsigned int ichan=newRecoData.size()-1; ichan>1;ichan--) {
                    // y[n] = x[n] - b1*x[n-1] - b2*x[n-2]
                    // there will be no corrections for channel 0 and 1
                    
                    newRecoData[ichan]= newRecoData[ichan]-b1*newRecoData[ichan-1]- b2*newRecoData[ichan-2];
                    newRecoData[ichan-1]= newRecoData[ichan-1]+b1*newRecoData[ichan-1];
                    newRecoData[ichan-2]= newRecoData[ichan-2]+b2*newRecoData[ichan-2];
                }
                
                
            }
            TrackerDataImpl * newDataImpl = new TrackerDataImpl();
            newDataImpl->setChargeValues(newRecoData);
            chipIDEncoder[ALIBAVA::ALIBAVADATA_ENCODE_CHIPNUM] = chipnum;
            chipIDEncoder.setCellID(newDataImpl);
            newDataColVec->push_back(newDataImpl);
        }
        
        alibavaEvent->addCollection(newDataColVec, _outputCollectionName);
        
    } catch ( lcio::DataNotAvailableException ) {
        // do nothing again
        streamlog_out( DEBUG0 ) << "Collection ("<<getInputCollectionName()<<") not found in event number "<< alibavaEvent->getEventNumber()<< endl;
    }
    
    
}

void AlibavaApplyCrossTalk::check (LCEvent * /* evt */ ) {
    // nothing to check here - could be used to fill check plots in reconstruction processor
}


void AlibavaApplyCrossTalk::end() {
}

void AlibavaApplyCrossTalk::fillHistos(TrackerDataImpl * /* trkdata */){
    
}


void AlibavaApplyCrossTalk::bookHistos(){
    
//    streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << endl;
}


