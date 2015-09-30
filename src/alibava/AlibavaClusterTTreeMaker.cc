/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


// alibava includes ".h"
#include "AlibavaClusterTTreeMaker.h"
#include "AlibavaBaseHistogramMaker.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"
#include "AlibavaCluster.h"


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


AlibavaClusterTTreeMaker::AlibavaClusterTTreeMaker () :
AlibavaBaseProcessor("AlibavaClusterTTreeMaker"),
// List of Histogram names, initialized here.
_treeName("clusterTree"),
_runnumber(-1),
_eventnumber(-1),
_chipNum(-1),
_clusterID(-1),
_clusterSize(-1),
_totalSignal(-1),
_totalSNR(-1),
_seedChanNum(-1),
_isSensitiveAxisX(true),
_signalPolarity(0)
{
	
	// modify processor description
	_description =
	"AlibavaClusterTTreeMaker creates a ROOT tree which stores cluster information ";
	
	
	// first of register the input collection
	registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
									 "Input ALIBAVA Cluster collection name",
									 _inputCollectionName, string("alibava_cluster") );
	
    // optional parameters, these parameters are defined in AlibavaBaseHistogramMaker
    registerOptionalParameter ("NoiseInputFile",
                               "The filename where the noise values stored",
                               _pedestalFile , string(ALIBAVA::NOTSET));
    
    registerOptionalParameter ("NoiseCollectionName",
                               "Noise collection name",
                               _noiseCollectionName, string (ALIBAVA::NOTSET));
    
	
}


void AlibavaClusterTTreeMaker::init () {
	streamlog_out ( MESSAGE4 ) << "Running init" << endl;
	
	// this method is called only once even when the rewind is active
	// usually a good idea to
	printParameters ();
    
    _tree = new TTree(_treeName.c_str(), "a ROOT tree which stores cluster information");
	
    _tree->Branch("runnumber",&_runnumber);
    _tree->Branch("eventnumber",&_eventnumber);
    _tree->Branch("chipNum",&_chipNum);
    _tree->Branch("clusterID",&_clusterID);
    _tree->Branch("clusterSize",&_clusterSize);
    _tree->Branch("totalSignal",&_totalSignal);
    _tree->Branch("totalSNR",&_totalSNR);
    _tree->Branch("seedChanNum",&_seedChanNum);
    _channums = new std::vector<int>();
    _tree->Branch("channums",&_channums);
    _signals = new std::vector<float>();
    _tree->Branch("signals",&_signals);
    _snrs = new std::vector<float>();
    _tree->Branch("snrs",&_snrs);
    _tree->Branch("isSensitiveAxisX",&_isSensitiveAxisX);
    _tree->Branch("signalPolarity",&_signalPolarity);
    
    _rootObjectMap.insert(make_pair(_treeName, _tree));

   gROOT->ProcessLine("#include <vector>"); 
}


void AlibavaClusterTTreeMaker::processRunHeader (LCRunHeader * rdr) {
	streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;
	
	// Add processor name to the runheader
	auto_ptr<AlibavaRunHeaderImpl> arunHeader ( new AlibavaRunHeaderImpl(rdr)) ;
	arunHeader->addProcessor(type());
    _runnumber = arunHeader->getRunNumber();
    
    // get and set selected chips
    setChipSelection(arunHeader->getChipSelection());

    // set pedestal and noise values (if it is defined)
    setPedestals();

}


void AlibavaClusterTTreeMaker::bookHistos(){
// does nothing
//	streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << endl;
}



void AlibavaClusterTTreeMaker::processEvent (LCEvent * anEvent) { // HERE look for it

	
	AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
	_eventnumber = alibavaEvent->getEventNumber();

	/////////////////////////////
	// Now loop over clusters //
	LCCollectionVec * collectionVec;
	unsigned int noOfCluster;
	try{
		collectionVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( getInputCollectionName() ) ) ;
		noOfCluster = collectionVec->getNumberOfElements();
		
		for ( size_t i = 0; i < noOfCluster; ++i ){
			TrackerDataImpl * trkdata = dynamic_cast< TrackerDataImpl * > ( collectionVec->getElementAt( i ) ) ;
			
            AlibavaCluster anAlibavaCluster(trkdata);
		_channums->clear();
		_signals->clear();
		_snrs->clear();
            _chipNum = anAlibavaCluster.getChipNum();
            _clusterID = anAlibavaCluster.getClusterID();
            _clusterSize = anAlibavaCluster.getClusterSize();
            _totalSignal = anAlibavaCluster.getTotalSignal();
            _totalSNR = anAlibavaCluster.getTotalSNR(getNoiseOfChip(_chipNum));
            _seedChanNum = anAlibavaCluster.getSeedChanNum();
std::vector<float> SNRS = anAlibavaCluster.getSNRs(getNoiseOfChip(_chipNum));
for (int i=0; i<_clusterSize; i++){
_channums->push_back(anAlibavaCluster.getChanNum(i));
_signals->push_back(anAlibavaCluster.getSignal(i));
_snrs->push_back(SNRS[i]);
}
         /*   
	_channums = anAlibavaCluster.getChanNums();
            _signals = anAlibavaCluster.getSignals();
            _snrs = anAlibavaCluster.getSNRs(getNoiseOfChip(_chipNum));
	*/
            _isSensitiveAxisX = anAlibavaCluster.getIsSensitiveAxisX();
            _signalPolarity = anAlibavaCluster.getSignalPolarity();
	    //TBranch *b_channums = (TBranch*)_tree->GetBranch("channums");
	    //b_channums->SetAddress(&_channums);
	    //TBranch *b_signals = (TBranch*)_tree->GetBranch("signals");
	    //b_signals->SetAddress(&_signals);
            //TBranch *b_snrs = (TBranch*)_tree->GetBranch("snrs");		
            //b_snrs->SetAddress(&_snrs);		
            _tree->Fill();

		}
		
	} catch ( lcio::DataNotAvailableException ) {
		// do nothing again
		streamlog_out( ERROR5 ) << "Collection ("<<getInputCollectionName()<<") not found! " << endl;
	}
	
}

void AlibavaClusterTTreeMaker::fillHistos(TrackerDataImpl * /* trkdata */){
	
}

void AlibavaClusterTTreeMaker::check (LCEvent * /* evt */ ) {
	// nothing to check here
}


void AlibavaClusterTTreeMaker::end() {
	
	if (_numberOfSkippedEvents > 0)
		streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents<<" events skipped since they are masked" << endl;
	streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
	
}
