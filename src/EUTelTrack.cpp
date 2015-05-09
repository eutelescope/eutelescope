#include "EUTelTrack.h"
using namespace eutelescope;
EUTelTrack::EUTelTrack(){
} 
EUTelTrack::EUTelTrack(const EUTelTrack& track): IMPL::TrackImpl(){
	setChi2(track.getChi2());
	setNdf(track.getNdf());
	for(size_t i=0; i<track.getTracks().size();++i){
		addTrack(track.getTracks().at(i));
	}
}
EUTelTrack::EUTelTrack(const EUTelTrack& track, bool copyContents){
	if( track.getChi2()== 0 or track.getNdf() == 0){
		streamlog_out(MESSAGE5)<<"Chi: "<<track.getChi2() <<" ndf: "<<track.getNdf() << std::endl;
		throw(lcio::Exception("You are trying to create a track that is empty. With another track that has no chi2 or degrees of freedom.")); 	
	}
	setChi2(track.getChi2());
	setNdf(track.getNdf());
}
//getters
std::vector<EUTelState> EUTelTrack::getStates(){
	std::vector<EUTelState> states;
	for(size_t i=0; i<getTracks().size();++i){
		EUTelState* state = static_cast<EUTelState*>(getTracks().at(i));
		states.push_back(*state);
	}
	return states;
}
std::vector<EUTelState*> EUTelTrack::getStatesPointers(){//This will return the pointers to the states. This is needed if we want to change the contents of the track and not just copyu like getStates()
	std::vector<EUTelState*> states;
	for(size_t i=0; i<getTracks().size();++i){
		EUTelState* state = static_cast<EUTelState*>(getTracks().at(i));
		states.push_back(state);
	}
	return states;
}
unsigned int EUTelTrack::getNumberOfHitsOnTrack() const {
	//streamlog_out( DEBUG1 ) << "EUTelTrack::EUTelTrack::getNumberOfHitsOnTrack()---------------------------BEGIN" << std::endl;
	unsigned int numberOfHitsOnTrack =0;
	const EVENT::TrackVec& states = getTracks();
	if(states.size() == 0){
		throw(lcio::Exception("The number of states is 0.")); 	
	}
	//streamlog_out(DEBUG0) <<"The number of states " << states.size()<<std::endl; 
	for(size_t i =0; i< states.size();++i){
		//streamlog_out(DEBUG0) <<"The states memory address for loop number "<<i<<" " << &states<<std::endl; 
		//const EVENT::TrackerHitVec& hit = states[i]->getTrackerHits();
		if(states[i]->getTrackerHits().size() == 0){
			continue;
		}
		if(states[i]->getTrackerHits().size()>1){
			streamlog_out( DEBUG1 ) << "The number of hits in one state is: "<< states[i]->getTrackerHits().size() << std::endl;
			throw(lcio::Exception("The number of hits for the state is greater than 1.")); 	
		}
		numberOfHitsOnTrack++;
	}
	//streamlog_out( DEBUG1 ) << "EUTelTrack::EUTelTrack::getNumberOfHitsOnTrack()---------------------------END" << std::endl;
	return numberOfHitsOnTrack;
}

void EUTelTrack::print(){
	streamlog_out(DEBUG1) <<"TRACK==>"<< " Chi: "<<getChi2() <<" ndf: "<<getNdf() << std::endl; 
    std::vector<EUTelState> states = getStates();
	streamlog_out(DEBUG1) <<"STATES:"<<std::endl;
	for(unsigned int i=0; i < states.size(); ++i){
        states.at(i).print();
	}	

}
//Setters
void EUTelTrack::setTotalVariance(double rad){
setPhi(rad);

}

