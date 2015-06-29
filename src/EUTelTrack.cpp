#include "EUTelTrack.h"
using namespace eutelescope;
EUTelTrack::EUTelTrack(){
} 
EUTelTrack::EUTelTrack(const EUTelTrack& track){
    _chi2=0;
    _nDF=0;
    _var=0;
	setChi2(track.getChi2());
	setNdf(track.getNdf());
    setStates(track.getStatesCopy());
    setTotalVariance(track.getTotalVariance());

}
EUTelTrack::EUTelTrack(const EUTelTrack& track, bool copyContents){
    _chi2=0;
    _nDF=0;
    _var=0;
	setChi2(track.getChi2());
	setNdf(track.getNdf());
    setTotalVariance(track.getTotalVariance());
}
//getters
float EUTelTrack::getChi2() const {
    return _chi2;
}

float EUTelTrack::getNdf() const {
    return _nDF;
}


std::vector<EUTelState>& EUTelTrack::getStates(){
	return _states;
}
std::vector<EUTelState> EUTelTrack::getStatesCopy() const {
	return _states;
}

unsigned int EUTelTrack::getNumberOfHitsOnTrack() const {
//	if(_states.size() == 0){
//		throw(lcio::Exception("The number of states is 0.")); 	
//	}
    unsigned int numHits = 0;
	for(unsigned int i = 0; i< _states.size();++i){
        if(_states.at(i).getStateHasHit()){
            numHits++;
        }
	}
	return numHits;
}

void EUTelTrack::print(){
	streamlog_out(DEBUG1) <<"TRACK==>"<< " Chi: "<<getChi2() <<" ndf: "<<getNdf() <<". Path total variance: " << _var << std::endl; 
    std::vector<EUTelState> states = getStates();
	streamlog_out(DEBUG1) <<"STATES:"<<std::endl;
	for(unsigned int i=0; i < states.size(); ++i){
        states.at(i).print();
	}	

}
//Setters
void EUTelTrack::setChi2(float chi2){
    _chi2 = chi2;
}
void EUTelTrack::setNdf(float nDF){
    _nDF=nDF;
}

void EUTelTrack::setTotalVariance(double rad){
_var = rad;

}

void EUTelTrack::setState(EUTelState state){
    _states.push_back(state);
}
void EUTelTrack::setStates(std::vector<EUTelState> states){
    _states = std::vector<EUTelState>(states);
}
std::vector<double> EUTelTrack::getLCIOOutput(){
    std::vector<double> output;
    output.push_back(static_cast<double>(getChi2()));
    output.push_back(static_cast<double>(getNdf()));
    output.push_back(static_cast<double>(getTotalVariance()));
    return output;


}
void EUTelTrack::setTrackFromLCIOVec(std::vector<double> input){
    setChi2(input.at(0));
    setNdf( input.at(1));
    setTotalVariance(input.at(2));

}

