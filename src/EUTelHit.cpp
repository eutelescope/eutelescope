#include "EUTelHit.h"
using namespace eutelescope;

EUTelHit::EUTelHit(){
} 

EUTelHit::EUTelHit(EUTelHit* hit){
    _position[0] = hit->getPosition()[0];
    _position[1] = hit->getPosition()[1];
    _position[2] = hit->getPosition()[2];
    _id = hit->getID();
} 

const double* EUTelHit::getPosition() const {
    return &_position[0];
}
void EUTelHit::setPosition(const double * position){
    _position[0] = position[0];
    _position[1] = position[1];
    _position[2] = position[2];


}

int EUTelHit::getID() const {
    return _id;
}
void EUTelHit::setID(int id){
    _id = id;
}

std::vector<double> EUTelHit::getLCIOOutput(){
    std::vector<double> output;
    output.push_back(getID());
    output.push_back(getPosition()[0]);
    output.push_back(getPosition()[1]);
    output.push_back(getPosition()[2]);

    return output;
}
void EUTelHit::setTrackFromLCIOVec(std::vector<double> input){
    setID(input.at(0));
    const double pos[3] = {input.at(1),input.at(2),input.at(3)};
    setPosition(pos);
}
