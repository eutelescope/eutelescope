#include "EUTelTrack.h"
using namespace eutelescope;
EUTelTrack::EUTelTrack(){
} 
EUTelTrack::EUTelTrack(const EUTelTrack& track){
	setChi2(track.getChi2());
	setNdf(track.getNdf());
	for(int i=0; i<track.getTracks().size();++i){
		addTrack(track.getTracks().at(i));
	}
}
EUTelTrack::EUTelTrack(const EUTelTrack& track, bool copyContents){
	if( track.getChi2()== 0 or track.getNdf() == 0){
		streamlog_out(MESSAGE5)<<"Chi: "<<track.getChi2() <<" ndf: "<<track.getNdf() <<endl;
		throw(lcio::Exception(Utility::outputColourString("You are trying to create a track that is empty. With another track that has not chi2 or degrees of freedom.", "RED"))); 	
	}
	setChi2(track.getChi2());
	setNdf(track.getNdf());
}
//getters
std::vector<EUTelState> EUTelTrack::getStates(){
	std::vector<EUTelState> states;
	for(int i=0; i<getTracks().size();++i){
		EUTelState* state = static_cast<EUTelState*>(getTracks().at(i));
		states.push_back(*state);
	}
	return states;
}
std::vector<EUTelState*> EUTelTrack::getStatesPointers(){//This will return the pointers to the states. This is needed if we want to change the contents of the track and not just copyu like getStates()
	std::vector<EUTelState*> states;
	for(int i=0; i<getTracks().size();++i){
		EUTelState* state = static_cast<EUTelState*>(getTracks().at(i));
		states.push_back(state);
	}
	return states;
}
int EUTelTrack::getNumberOfHitsOnTrack() const {
	//streamlog_out( DEBUG1 ) << "EUTelTrack::EUTelTrack::getNumberOfHitsOnTrack()---------------------------BEGIN" << std::endl;
	int numberOfHitsOnTrack =0;
	const EVENT::TrackVec& states = getTracks();
	if(states.size() == 0){
		throw(lcio::Exception(Utility::outputColourString("The number of states is 0.", "RED"))); 	
	}
	//streamlog_out(DEBUG0) <<"The number of states " << states.size()<<std::endl; 
	for(int i =0; i< states.size();++i){
		//streamlog_out(DEBUG0) <<"The states memory address for loop number "<<i<<" " << &states<<std::endl; 
		const EVENT::TrackerHitVec& hit = states[i]->getTrackerHits();
		if(states[i]->getTrackerHits().size() == 0){
			continue;
		}
		if(states[i]->getTrackerHits().size()>1){
			streamlog_out( DEBUG1 ) << "The number of hits in one state is: "<< states[i]->getTrackerHits().size() << std::endl;
			throw(lcio::Exception(Utility::outputColourString("The number of hits for the state is greater than 1.", "RED"))); 	
		}
		numberOfHitsOnTrack++;
	}
	//streamlog_out( DEBUG1 ) << "EUTelTrack::EUTelTrack::getNumberOfHitsOnTrack()---------------------------END" << std::endl;
	return numberOfHitsOnTrack;
}

void EUTelTrack::print(){
	streamlog_out(DEBUG1) <<"TRACK INFORMATION//////////////////////////////////////////////////////////////////////////START"<<std::endl;
	streamlog_out(DEBUG1) << "Track contains " << getTracks().size() << " states " << std::endl;
	for(int i=0; i < getTracks().size(); ++i){
		const  EVENT::Track* state = getTracks().at(i);
		//Z0 contains the location. Need to static cast to EUTelState to use getLocation
		streamlog_out(DEBUG1) <<"State memory location "<< state << " The sensor location of the state " <<state->getZ0()<<std::endl;
		if(!state->getTrackerHits().empty()){
			streamlog_out(DEBUG1) <<"The hit ID of the state is "<<state->getTrackerHits().at(0)->id()<<std::endl;
		}
	}	
	streamlog_out(DEBUG1) <<"TRACK INFORMATION///////////////////////////////////////////////////////////////////////////////END"<<std::endl;

}
