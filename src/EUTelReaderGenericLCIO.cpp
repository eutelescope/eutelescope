#include "EUTelReaderGenericLCIO.h"
using namespace eutelescope;

EUTelReaderGenericLCIO::EUTelReaderGenericLCIO(){
} 
void EUTelReaderGenericLCIO::getColVec(std::vector<EUTelTrack> tracks,LCEvent* evt ,std::string colName ){
    streamlog_out(DEBUG1)<<"CREATE GENERIC CONTAINER..." <<std::endl;

    LCCollectionVec* colTrackVec = new LCCollectionVec(LCIO::LCGENERICOBJECT);
    LCCollectionVec* colStateVec = new LCCollectionVec(LCIO::LCGENERICOBJECT);
    LCCollectionVec* colHitVec = new LCCollectionVec(LCIO::LCGENERICOBJECT);
    LCCollectionVec* relTrackStateVec = new LCCollectionVec(LCIO::LCRELATION);
    LCCollectionVec* relStateHitVec = new LCCollectionVec(LCIO::LCRELATION);


    for(size_t i=0 ; i < tracks.size(); i++){
        IMPL::LCGenericObjectImpl* conTrack = new  IMPL::LCGenericObjectImpl();  
        //Save everything as double and down cast later.
        for(size_t j=0; j < tracks.at(i).getLCIOOutput().size(); j++){
            streamlog_out(DEBUG1)<<"Fill number: " << j << " Value " << tracks.at(i).getLCIOOutput().at(j) <<std::endl;
            conTrack->setDoubleVal (j, tracks.at(i).getLCIOOutput().at(j));
        }
        streamlog_out(DEBUG1)<<"Fill all track information...        Double number: "<< conTrack->getNDouble()  <<std::endl;
        colTrackVec->push_back(static_cast<EVENT::LCGenericObject*>(conTrack));
        streamlog_out(DEBUG1)<<"Tracks filled" <<std::endl;
        std::vector<EUTelState> states = tracks.at(i).getStates();
        for(size_t j=0 ; j < states.size(); j++){
            streamlog_out(DEBUG1)<<"Fill all state information " << " state location " << states.at(j).getLocation() <<std::endl;
            IMPL::LCGenericObjectImpl* conState = new  IMPL::LCGenericObjectImpl();  
            streamlog_out(DEBUG1)<<"Empty state container created" <<std::endl;
            for(size_t k=0 ; k < states.at(j).getLCIOOutput().size(); k++){
                streamlog_out(DEBUG1)<<"Fill number: " << k << " Value " << tracks.at(i).getStates().at(j).getLCIOOutput().at(k) <<std::endl;
                conState->setDoubleVal (k,states.at(j).getLCIOOutput().at(k));
            }
            IMPL::LCRelationImpl *relTrackState = new IMPL::LCRelationImpl(conTrack,conState); 
            colStateVec->push_back(static_cast<EVENT::LCGenericObject*>(conState));
            relTrackStateVec->push_back(static_cast<EVENT::LCRelation*>(relTrackState));
            if(states.at(j).getStateHasHit()){
                IMPL::LCGenericObjectImpl* conHit = new  IMPL::LCGenericObjectImpl();  
                for(size_t k=0 ; k < states.at(j).getHit().getLCIOOutput().size(); k++){
                    conHit->setDoubleVal (k, states.at(j).getHit().getLCIOOutput().at(k));
                }
                IMPL::LCRelationImpl *relStateHit = new IMPL::LCRelationImpl(conState,conHit); 
                colHitVec->push_back(static_cast<EVENT::LCGenericObject*>(conHit));
                relStateHitVec->push_back(static_cast<EVENT::LCRelation*>(relStateHit));

            }

        }
    }
    streamlog_out(DEBUG1)<<"Add collection to event!" <<std::endl;
    evt->addCollection(colTrackVec,"TrackFOR" + colName   );
    evt->addCollection(colStateVec,"StatesFOR" + colName);
    evt->addCollection(colHitVec,"HitsFOR" + colName);
    evt->addCollection(relTrackStateVec,"TrackStateFOR" + colName);
    evt->addCollection(relStateHitVec,"StateHitFOR" + colName);

} 

std::vector<EUTelTrack> EUTelReaderGenericLCIO::getTracks( LCEvent* evt, std::string colName){
    std::vector<EUTelTrack> tracks; 
    streamlog_out(DEBUG1)<<"Open Collections... " <<std::endl;

    LCCollection* relTrackStates =  evt->getCollection("TrackStateFOR"+ colName);
    LCCollection* relStatesHits =  evt->getCollection("StateHitFOR"+ colName);
    streamlog_out(DEBUG1)<<"Open!" <<std::endl;

    std::vector<int> trackIDVec;
    //Loop a link between tracks->states. //Remember multiple tracks for each collection 
    for (int iCol = 0; iCol < relTrackStates->getNumberOfElements(); iCol++) {//Loop through each track->state link
        EVENT::LCRelation* relTrackState = static_cast<EVENT::LCRelation*>(relTrackStates->getElementAt(iCol));
        EVENT::LCGenericObject* track  =  static_cast<EVENT::LCGenericObject*>(relTrackState->getFrom());
        int trackID = track->id();
        if(std::find(trackIDVec.begin(), trackIDVec.end(), trackID) == trackIDVec.end()){//This is a list of tracks already created
            //If track is new enter here.
            trackIDVec.push_back(trackID);
            std::vector<double> trackInput;
            for(int i =0 ; i < track->getNDouble(); i++){
                trackInput.push_back(track->getDoubleVal(i)); 
            }
            EUTelTrack track;
            track.setTrackFromLCIOVec(trackInput);
            for (int jCol = 0; jCol < relTrackStates->getNumberOfElements(); jCol++) {//Loop through track->state look for state linked to that track
                EVENT::LCRelation* relTrackState = static_cast<EVENT::LCRelation*>(relTrackStates->getElementAt(jCol));
                EVENT::LCGenericObject* state  =  static_cast<EVENT::LCGenericObject*>(relTrackState->getTo());
                EVENT::LCGenericObject* trackCheck  =  static_cast<EVENT::LCGenericObject*>(relTrackState->getFrom());
                if(trackCheck->id() == trackID){//If this is true then the state must be part of this track.
                    int stateID = state->id();
                    std::vector<double> stateInput;
                    for(int i =0 ; i < state->getNDouble(); i++){
                        stateInput.push_back(state->getDoubleVal(i)); 
                    }
                    EUTelState state;
                    state.setTrackFromLCIOVec(stateInput);
                    for (int kCol = 0; kCol < relStatesHits->getNumberOfElements(); kCol++) {//Loop through all state->hit link
                        EVENT::LCRelation* relStateHit = static_cast<EVENT::LCRelation*>(relStatesHits->getElementAt(kCol));
                        EVENT::LCGenericObject* stateCheck  =  static_cast<EVENT::LCGenericObject*>(relStateHit->getFrom());
                        EVENT::LCGenericObject* hit  =  static_cast<EVENT::LCGenericObject*>(relStateHit->getTo());
                        if(stateCheck->id() == stateID){//If this is true then you have the correct hit.
                            streamlog_out(DEBUG1)<<"Found correct ID. Add hit now..." <<std::endl;
                            std::vector<double> hitInput;
                            for(int i =0 ; i < hit->getNDouble(); i++){
                                hitInput.push_back(hit->getDoubleVal(i)); 
                            }
                            EUTelHit hit;
                            hit.setTrackFromLCIOVec(hitInput);
                            state.setHit(hit);
                            streamlog_out(DEBUG1)<<"Hit added." <<std::endl;

                        }
                    }
                    track.setState(state);
                }
            }
            tracks.push_back(track);
        }
    }
    streamlog_out(DEBUG1)<<"Return "<< tracks.size() <<" tracks" <<std::endl;
    return tracks;
}



