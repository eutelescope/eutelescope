/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


// alibava includes ".h"
#include "AlibavaDataTTreeMaker.h"
#include "AlibavaBaseHistogramMaker.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"


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
#include "TTree.h"
#include "TROOT.h"
#include "TBranch.h"

// system includes <>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <memory>
#include <cstdlib>
#include <ctime>
#include <vector>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;


AlibavaDataTTreeMaker::AlibavaDataTTreeMaker () :
AlibavaBaseProcessor("AlibavaDataTTreeMaker"),
// List of Histogram names, initialized here.
_treeName("dataTree"),
_runnumber(-1),
_eventnumber(-1),
_tdctime(0),
_temperature(0),
_chipNum(-1)
{
    
    // modify processor description
    _description =
    "AlibavaDataTTreeMaker creates a ROOT tree which stores data information ";
    
    
    // first of register the input collection
    registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
                             "Input Alibava data collection name",
                             _inputCollectionName, string("alibava_data") );
    
    
}


void AlibavaDataTTreeMaker::init () {
    streamlog_out ( MESSAGE4 ) << "Running init" << endl;
    
    // this method is called only once even when the rewind is active
    // usually a good idea to
    printParameters ();
    _treeName = _inputCollectionName;
    _treeName += string("_Tree");
    _tree = new TTree(_treeName.c_str(), "a ROOT tree which stores cluster information");
    
    _tree->Branch("runnumber",&_runnumber);
    _tree->Branch("eventnumber",&_eventnumber);
    _tree->Branch("tdctime",&_tdctime);
    _tree->Branch("temperature",&_temperature);
    _tree->Branch("chipNum",&_chipNum);
    _data = new std::vector<double>();
    _tree->Branch("data",&_data);
    
    _rootObjectMap.insert(make_pair(_treeName, _tree));
    
}


void AlibavaDataTTreeMaker::processRunHeader (LCRunHeader * rdr) {
    streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;
    
    // Add processor name to the runheader
    auto_ptr<AlibavaRunHeaderImpl> arunHeader ( new AlibavaRunHeaderImpl(rdr)) ;
    arunHeader->addProcessor(type());
    _runnumber = arunHeader->getRunNumber();
    
}


void AlibavaDataTTreeMaker::bookHistos(){
    // does nothing
    //	streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << endl;
}



void AlibavaDataTTreeMaker::processEvent (LCEvent * anEvent) { // HERE look for it
    
    
    AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
    _eventnumber = alibavaEvent->getEventNumber();
    _tdctime = alibavaEvent->getEventTime();
    _temperature = alibavaEvent->getEventTemp();
    /////////////////////////////
    // Now loop over clusters //
    LCCollectionVec * collectionVec;
    unsigned int noOfChips;
    try{
        
        collectionVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( getInputCollectionName() ) ) ;
        noOfChips = collectionVec->getNumberOfElements();
        
        for ( size_t i = 0; i < noOfChips; ++i ){
            TrackerDataImpl * trkdata = dynamic_cast< TrackerDataImpl * > ( collectionVec->getElementAt( i ) ) ;
            
	    _chipNum = getChipNum(trkdata);
            FloatVec datavec;
            datavec = trkdata->getChargeValues();
            _data->clear();
            // check if trkdata has size of ALIBAVA::NOOFCHANNELS
            if (datavec.size() != (unsigned int) ALIBAVA::NOOFCHANNELS) {
                streamlog_out( ERROR1 ) << "Collection ("<<getInputCollectionName()<<") doesn't have "<< ALIBAVA::NOOFCHANNELS<< "channels "<< endl;
                
            }
            else{
                for (unsigned int i=0; i<datavec.size(); i++){
                    _data->push_back(datavec[i]);
                    
                }
                
                _tree->Fill();
                
            }
            
        }    
        } catch ( lcio::DataNotAvailableException ) {
            // do nothing again
            streamlog_out( DEBUG1 ) << "Collection ("<<getInputCollectionName()<<") not found! " << endl;
        }
        
        
    }
    
    void AlibavaDataTTreeMaker::fillHistos(TrackerDataImpl * /* trkdata */){
        
    }
    
    void AlibavaDataTTreeMaker::check (LCEvent * /* evt */ ) {
        // nothing to check here
    }
    
    
    void AlibavaDataTTreeMaker::end() {
        
        if (_numberOfSkippedEvents > 0)
            streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents<<" events skipped since they are masked" << endl;
        streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
        
    }
