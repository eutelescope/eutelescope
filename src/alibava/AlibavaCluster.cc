/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */

// alibava includes ".h"
#include "AlibavaCluster.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <lcio.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <IMPL/TrackerDataImpl.h>

// system includes <>
#include <iomanip>


using namespace std;
using namespace lcio;
using namespace alibava;

AlibavaCluster::AlibavaCluster() :
_channums(),
_signals(),
_chipNum(-1),
_seedChanNum(-1),
_clusterID(-1),
_isSensitiveAxisX(true),
_signalPolarity(-1.0)
{
    // does nothing
}

AlibavaCluster::AlibavaCluster(TrackerDataImpl* trkdata) :
_channums(),
_signals(),
_chipNum(-1),
_seedChanNum(-1),
_clusterID(-1),
_isSensitiveAxisX(true),
_signalPolarity(-1.0)
{
    CellIDDecoder<TrackerDataImpl> clusterIDDecoder(ALIBAVA::ALIBAVACLUSTER_ENCODE);
    _clusterID = static_cast<int> ( clusterIDDecoder( trkdata )[ALIBAVA::ALIBAVACLUSTER_ENCODE_CLUSTERID] );
    _seedChanNum = static_cast<int> ( clusterIDDecoder( trkdata )[ALIBAVA::ALIBAVACLUSTER_ENCODE_SEED] );
    _chipNum = static_cast<int> ( clusterIDDecoder( trkdata )[ALIBAVA::ALIBAVACLUSTER_ENCODE_CHIPNUM] );
    
    int sensitiveAxisX= static_cast<int> ( clusterIDDecoder( trkdata )[ALIBAVA::ALIBAVACLUSTER_ENCODE_ISSENSITIVEAXISX] );
    if (sensitiveAxisX == 0)
        _isSensitiveAxisX = false;
    else
        _isSensitiveAxisX = true;
    
    int negSignal = static_cast<int> ( clusterIDDecoder( trkdata )[ALIBAVA::ALIBAVACLUSTER_ENCODE_ISSIGNALNEGATIVE] );
    if (negSignal == 0)
        _signalPolarity = 1.0;
    else
        _signalPolarity = -1.0;
    
    FloatVec data = trkdata->getChargeValues();
    // first number is eta,
    // then it goes like channel number, signal // see createTrackerData
    // So data size has to be odd! (because of eta)
    if (data.size() % 2 != 0)
        streamlog_out (ERROR5) << "Size in TrackerData that stores cluster information is not even! Either channel number or signal information is missing!"<< endl;
    
    // numbers are channel number and corresponding signal
    for (unsigned int imember=0; imember<data.size(); imember++) {
        _channums.push_back( int(data[imember]) );
        imember++;
        _signals.push_back( data[imember] );
    }
}




AlibavaCluster::~AlibavaCluster(){
    
}

void AlibavaCluster::createTrackerData(lcio::TrackerDataImpl * alibavaCluster){
    // put channel information in TrackerData
    // first number we put is channel number then signal
    FloatVec dataToStore;
    
    // then channel number and signal for each member
    for (int imember=0; imember<getClusterSize(); imember++) {
        dataToStore.push_back( float(getChanNum(imember)) );
        dataToStore.push_back( getSignal(imember) );
    }
    alibavaCluster->setChargeValues(dataToStore);
    
    // Cell ID encoding will done in the main processor
    
}
float AlibavaCluster::getSignalOnChannel(int channelnum){
	for (int imember=0; imember<getClusterSize(); imember++ ){
		if(getChanNum(imember)==channelnum) return getSignal(imember);
	}
	
	// if reaches here, this means that the channel number is not in the cluster
	return _unrealisticSignal;
}

float AlibavaCluster::getEta(){
	 
	if (getClusterSize() == 1) return -1;	

	int seedChan = getSeedChanNum();
	int leftChan = seedChan - 1;
	int rightChan = seedChan+1;
	
	// we will multiply all signal values by _signalPolarity to work on positive signal always
	float seedSignal = _signalPolarity * getSignalOnChannel(seedChan);
	float leftSignal = getSignalOnChannel(leftChan);
	float rightSignal = getSignalOnChannel(rightChan);
	
	if (leftSignal != _unrealisticSignal) leftSignal = _signalPolarity * leftSignal;
	if (rightSignal!= _unrealisticSignal) rightSignal= _signalPolarity * rightSignal;

	// if both right anf left channel is masked. Simply return -1
	// this case should not be saved by clustering algorithm anyways
	if (rightSignal == _unrealisticSignal && leftSignal == _unrealisticSignal ) {
		streamlog_out (WARNING3) << "Both neighbours are masked!"<<endl;
		return -1;
	}
	
	float eta = -1;
	// compare left and right signals
	// here both signal has to be positive (if not noise)
	// if one of the channel is masked, since it will be set to unrealisticSignal other signal will always be higher then the masked one.

	// Eta calculation: chargeOnLeftChannel / (chargeOnLeftChannel + chargeOnRightChannel)
	if ( leftSignal > rightSignal) {
		// then seed channel is on the right
		eta = leftSignal / ( leftSignal + seedSignal );
	}
	else {
		// seed channel is on the left
		eta = seedSignal / (seedSignal + rightSignal);
	}
		
	return eta;
	
}



float AlibavaCluster::getCenterOfGravity(){
    if (getClusterSize()==1)
        return getChanNum(0); // channel number of the only member
    
    float totalsignal = 0;
    float weightedSum = 0;
    for (int imember=0; imember<getClusterSize(); imember++) {
        totalsignal += getSignal(imember);
        weightedSum += getSignal(imember)*getChanNum(imember);
    }
    
    return weightedSum/totalsignal;
}
std::vector<int> AlibavaCluster::getChanNums(){
    return _channums;
}
std::vector<float> AlibavaCluster::getSignals(){
    return _signals;
}

int AlibavaCluster::getChanNum(int imember) {
    if (imember >= getClusterSize()) {
        streamlog_out(ERROR5)<< "Not enough member in this AlibavaCluster"<<endl;
        return -1;
    }
    return _channums[imember];
}

float AlibavaCluster::getSignal(int imember) {
    if (imember >= getClusterSize()) {
        streamlog_out(ERROR5)<< "Not enough member in this AlibavaCluster"<<endl;
        return 0;
    }
    return _signals[imember];
}
float AlibavaCluster::getTotalSNR(FloatVec noiseVec) {
    float snr_total = 0;
    for (int imember=0; imember<getClusterSize(); imember++) {
        int channel = getChanNum(imember);
        if (noiseVec[channel]==0) {
            streamlog_out (ERROR5)<<"Channels noise is zero! Check this cluster! Returning zero!"<<endl;
            return 0;
        }
        snr_total += getSignal(imember)/noiseVec[channel];
    }
    return snr_total * getSignalPolarity();
}

std::vector<float> AlibavaCluster::getSNRs(FloatVec noiseVec) {
    std::vector<float> snrs;
    snrs.clear();
    for (int imember=0; imember<getClusterSize(); imember++) {
        int channel = getChanNum(imember);
        if (noiseVec[channel]==0) {
            streamlog_out (ERROR5)<<"Channels noise is zero! Check this cluster! Returning zero!"<<endl;
            snrs.clear();
            return snrs;
        }
        snrs.push_back( getSignalPolarity() * (getSignal(imember)/noiseVec[channel]) );
    }
    return snrs;
}

float AlibavaCluster::getTotalSignal() {
    float totalsignal=0;
    for(int imember=0; imember<getClusterSize(); imember++)
        totalsignal += _signals[imember];
    return totalsignal;
}

void AlibavaCluster::add(int achannum, float asignal){
    _channums.push_back(achannum);
    _signals.push_back(asignal);
}

int AlibavaCluster::getClusterSize(){
    if (_channums.size() != _signals.size()) {
        streamlog_out (ERROR5)<<"There is something wrong with this AlibavaCluster, number of channels is not consistent with number of signals, registered in this cluster!"<<endl;
    }
    return int (_channums.size());
}

bool AlibavaCluster::has_seed(){
    if (_seedChanNum != -1)
        return true;
    return false;
}

///////////////////////
// Setters - Getters //
///////////////////////

// setter / getter for _chipNum
int AlibavaCluster::getChipNum(){
    
    return _chipNum;
}
void AlibavaCluster::setChipNum(int chipnum){
    
    _chipNum = chipnum;
}

// setter / getter for _seedChanNum
int AlibavaCluster::getSeedChanNum(){
    if (!has_seed()) {
        streamlog_out(ERROR3)<<"Seed channel not defined!"<<endl;
    }
    return _seedChanNum;
}
void AlibavaCluster::setSeedChanNum(int seedChanNum){
    
    _seedChanNum = seedChanNum;
}

// setter / getter for _clusterID
int AlibavaCluster::getClusterID(){
    if (_clusterID<0) {
        streamlog_out(ERROR3)<<"Cluster ID is not defined!"<<endl;
    }
    return _clusterID;
}
void AlibavaCluster::setClusterID(int clusterID){
    
    _clusterID = clusterID;
}

// setter / getter for _isSensitiveAxisX
void AlibavaCluster::setIsSensitiveAxisX(bool isSensitiveAxisX){
    _isSensitiveAxisX = isSensitiveAxisX;
}
bool AlibavaCluster::getIsSensitiveAxisX(){
    return _isSensitiveAxisX;
}

// setter / getter for _signalPolarity
double AlibavaCluster::getSignalPolarity(){
    return _signalPolarity;
}
void AlibavaCluster::setSignalPolarity(double signalPolarity){
    _signalPolarity = signalPolarity;
}

void AlibavaCluster::print(){
    streamlog_out(MESSAGE5)<<"************ Cluster ************"<<endl;
    streamlog_out(MESSAGE5)<<"ClusterID: "<< getClusterID()<< " ChipNum: "<<getChipNum()<<endl;
    
    string axis;
    if (getIsSensitiveAxisX()) axis = string("X");
    else axis = string("Y");
    streamlog_out(MESSAGE5)<<"SensitiveAxis: "<< axis << " SignalPolarity: "<<getSignalPolarity()<<endl;
    
    streamlog_out(MESSAGE5)<<"SeedChannel: "<< getSeedChanNum() << " ClusterSize: "<<getClusterSize()<< " Eta: " <<getEta() <<endl;
    streamlog_out(MESSAGE5)<<"ClusterMembers: channelNum, signal";
    streamlog_out(MESSAGE5)<<"           ";
    for (int imember=0; imember<getClusterSize(); imember++)
        streamlog_out(MESSAGE5)<<getChanNum(imember)<<", "<<getSignal(imember)<< ", ";
    streamlog_out(MESSAGE5)<<endl;
    streamlog_out(MESSAGE5)<<"***********************************"<<endl;
    
    
    
}









































