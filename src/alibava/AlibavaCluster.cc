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
using namespace marlin;
using namespace alibava;

AlibavaCluster::AlibavaCluster() :
_channums(),
_signals(),
_eta(-1),
_chipNum(-1),
_seedChanNum(-1),
_clusterID(-1),
_isSensitiveAxisX(true),
_signalPolarity(-1)
{
// does nothing
}

AlibavaCluster::AlibavaCluster(TrackerDataImpl* trkdata) :
_channums(),
_signals(),
_eta(-1),
_chipNum(-1),
_seedChanNum(-1),
_clusterID(-1),
_isSensitiveAxisX(true),
_signalPolarity(-1)
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
		_signalPolarity = 1;
	else
		_signalPolarity = -1;
	
	FloatVec data = trkdata->getChargeValues();
	// first number is eta,
	// then it goes like channel number, signal // see createTrackerData
	// So data size has to be odd! (because of eta)
	if (data.size() % 2 != 1)
		streamlog_out (ERROR5) << "Size in TrackerData that stores cluster information is not even! Either channel number or signal information is missing!"<< endl;

	// first number is eta, others are channel number and corresponding signal
	setEta(data.at(0));
	for (unsigned int imember=1; imember<data.size(); imember++) {
		_channums.push_back( int(data[imember]) );
		imember++;
		_signals.push_back( data[imember] );
	}
}




AlibavaCluster::~AlibavaCluster(){
	
}

void AlibavaCluster::createTrackerData(TrackerDataImpl * alibavaCluster){
	// put channel information in TrackerData
	// first number we put is channel number then signal
	FloatVec dataToStore;
	
	// first store eta
	dataToStore.push_back(getEta());
	// then channel number and signal for each member
	for (int imember=0; imember<getClusterSize(); imember++) {
		dataToStore.push_back( float(getChanNum(imember)) );
		dataToStore.push_back( getSignal(imember) );
	}
	alibavaCluster->setChargeValues(dataToStore);

	// Cell ID encoding will done in the main processor
	
}
/*
float AlibavaCluster::getEta(){
	// if clustersize is 1, nothing to calculate, returns -1 
	if (getClusterSize()==1)
		return -1;
	
	float leftNeighs = 0;
	float rightNeighs = 0;
	float centerOfGravity = getCenterOfGravity();
	
	for (int imember=0; imember<getClusterSize(); imember++) {
		// since we will compare this number with float
		float memberChanNum = float ( getChanNum(imember) );
		if ( memberChanNum < centerOfGravity )
			leftNeighs += getSignal(imember);
		else if ( memberChanNum > centerOfGravity )
			rightNeighs += getSignal(imember);
	}
	
	float eta = leftNeighs/(leftNeighs+rightNeighs);
	return eta;
}
*/

float AlibavaCluster::getEta(){
	return _eta;
}

void AlibavaCluster::setEta(float eta){
	_eta = eta;
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
/*
// setter / getter for _sensorIDOffset
int AlibavaCluster::getSensorIDOffset(){
	return _sensorIDOffset;
}
void AlibavaCluster::setSensorIDOffset(int sensorIDOffset){
	_sensorIDOffset = sensorIDOffset;
}

// setter / getter for _missingCorrdinateValue
int AlibavaCluster::getMissingCorrdinateValue(){
	return _missingCorrdinateValue;
}
void AlibavaCluster::setMissingCorrdinateValue(int missingCorrdinateValue){
	_missingCorrdinateValue = missingCorrdinateValue;
}
*/
// setter / getter for _signalPolarity
int AlibavaCluster::getSignalPolarity(){
	return _signalPolarity;
}
void AlibavaCluster::setSignalPolarity(int signalPolarity){
	_signalPolarity = signalPolarity;
}

void AlibavaCluster::print(){
	streamlog_out(MESSAGE1)<<"************ Cluster ************"<<endl;
	streamlog_out(MESSAGE1)<<"ClusterID: "<< getClusterID()<< " ChipNum: "<<getChipNum()<<endl;

	string axis;
	if (getIsSensitiveAxisX()) axis = string("X");
	else axis = string("Y");
	streamlog_out(MESSAGE1)<<"SensitiveAxis: "<< axis << " SignalPolarity: "<<getSignalPolarity()<<endl;
	
	streamlog_out(MESSAGE1)<<"SeedChannel: "<< getSeedChanNum() << " ClusterSize: "<<getClusterSize()<< " Eta: " <<getEta() <<endl;
	streamlog_out(MESSAGE1)<<"ClusterMembers: channelNum, signal";
	streamlog_out(MESSAGE1)<<"           ";
	for (int imember=0; imember<getClusterSize(); imember++)
		streamlog_out(MESSAGE1)<<getChanNum(imember)<<", "<<getSignal(imember)<< ", ";
	streamlog_out(MESSAGE1)<<endl;
	streamlog_out(MESSAGE1)<<"***********************************"<<endl;
	
	
	
}









































