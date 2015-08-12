#include "EUTelTrack.h"
using namespace eutelescope;
EUTelTrack::EUTelTrack(){
} 
EUTelTrack::EUTelTrack(const EUTelTrack& track):
_chi2(0),
_nDF(0),
_radPerTotal(0)
{
	setChi2(track.getChi2());
	setNdf(track.getNdf());
    setStates(track.getStatesCopy());
    setRadPerTotal(track.getRadPerTotal());
    setQOverP(track.getQOverP()); 


}
EUTelTrack::EUTelTrack(const EUTelTrack& track, bool copyContents):
_chi2(0),
_nDF(0),
_radPerTotal(0)
{
	setChi2(track.getChi2());
	setNdf(track.getNdf());
    setRadPerTotal(track.getRadPerTotal());
    setQOverP(track.getQOverP()); 

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
std::vector<EUTelHit> EUTelTrack::getHitsCopy() const {
    std::vector<EUTelHit> hits;
    for(std::vector<EUTelState>::const_iterator itSt = _states.begin(); itSt != _states.end(); ++itSt){
        if(itSt->getStateHasHit()){
            hits.push_back(itSt->getHitCopy());
        }
    }
    return hits;
}

std::vector<EUTelState> EUTelTrack::getStatesCopy() const {
	return _states;
}
std::vector<int>  EUTelTrack::getPlaIDs() const {
    //Using an iterator did not seem to work? 
    std::vector<int> planes;
    for(size_t i = 0; i < _states.size() ; ++i){
        EUTelState state = _states.at(i);
        planes.push_back(state.getLocation());
    }

    return planes;
}
std::vector<int>  EUTelTrack::getPlaIDDUTs() const {
    //Using an iterator did not seem to work? 
    std::vector<int> planes;
    for(size_t i = 0; i < _states.size() ; ++i){
        EUTelState state = _states.at(i);
        if(state.getLocation() > 5){
            planes.push_back(state.getLocation());
        }
    }

    return planes;
}



unsigned int EUTelTrack::getNumberOfHitsOnTrack() const {
    unsigned int numHits = 0;
	for(unsigned int i = 0; i< _states.size();++i){
        if(_states.at(i).getStateHasHit()){
            numHits++;
        }
	}
	return numHits;
}

void EUTelTrack::print(){
	streamlog_out(DEBUG1) <<"TRACK==>"<< " Chi: "<<getChi2() <<" ndf: "<<getNdf() <<". Path total variance: " << _radPerTotal <<" OoverP:" << getQOverP()  << std::endl; 
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
    output.push_back(static_cast<double>(getRadPerTotal()));
    output.push_back(static_cast<double>(getQOverP()));
    return output;


}
void EUTelTrack::setTrackUsingCorrection(TVectorD corrections){
//    std::cout<<"Q over P before: " << getQOverP() << std::endl;
//    std::cout<<"Corr: " << corrections[0] << std::endl;
   setQOverP(corrections[0] + getQOverP()); 
 //   std::cout<<"Q over P after: " << getQOverP() << std::endl;

}


void EUTelTrack::setTrackFromLCIOVec(std::vector<double> input){
    setChi2(input.at(0));
    setNdf( input.at(1));
    setRadPerTotal(input.at(2));
    setQOverP(input.at(3));

}

