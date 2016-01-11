/*
 * EUTelTrackSelection.cpp 
 * 
 * Created on: April 19th 2015 
 *     author:Alexander Morton 
 * 
 *          This class creates a track based selection object which will take a vector of tracks and output a subset of these based on an a combined series of cuts. 
 *          Cuts are made on the tracks which have particular hits and chi2 
 * 
 *
 */



#include "EUTelTrackSelection.h"

using namespace eutelescope;
EUTelTrackSelection::EUTelTrackSelection(){
} 
bool EUTelTrackSelection::removeTracksWithHitsOnPlane(EUTelTrack track,std::vector<int> sensors){
    if(sensors.empty()){
        return true;
    }
    std::vector<EUTelState> states =  track.getStates();
    for(std::vector<EUTelState>::iterator it = states.begin(); it != states.end(); ++it){
        //If there is a hit and the location of the sensor is not supposed to have a hit then we throw a skip event exception.
        //I am unsure if this is the best way of doing this. Advice?
        if(it->getStateHasHit() and std::find(sensors.begin(), sensors.end(), it->getLocation())!=sensors.end() ){
            return false;
        }
    }
    return true;
} 


