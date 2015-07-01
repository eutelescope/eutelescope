#include "EUTelHit.h"
using namespace eutelescope;

EUTelHit::EUTelHit():
_locationKnown(false)
{
} 

EUTelHit::EUTelHit(EUTelHit* hit){
    _position[0] = hit->getPosition()[0];
    _position[1] = hit->getPosition()[1];
    _position[2] = hit->getPosition()[2];
    _id = hit->getID();
    std::cout<<"INside" << std::endl;
    if(hit->_locationKnown){
        _location = hit->getLocation();
        _locationKnown=true;
    }else{
        _locationKnown=false;
    }
    std::cout<<"INside2" << std::endl;


} 

EUTelHit::EUTelHit(EVENT::TrackerHit* hit){
    _position[0] = hit->getPosition()[0];
    _position[1] = hit->getPosition()[1];
    _position[2] = hit->getPosition()[2];
    _id = hit->id();
    const int hitLoc = Utility::getSensorIDfromHit( static_cast<IMPL::TrackerHitImpl*> (hit) );
    _location = hitLoc;
    _locationKnown=true;

} 

const double* EUTelHit::getPosition() const {
    return &_position[0];
}
int EUTelHit::getLocation() const {
    if(_locationKnown){
        return _location;
    }else{
        throw(std::string("The hit location is not known!!!")); 
    }

}


TVector3 EUTelHit::getPositionGlobal() const {
	const double* local =  getPosition();
	const double posLocal[3] = {local[0],local[1],local[2]};
  double posGlobal[3];
	geo::gGeometry().local2Master(getLocation() ,posLocal,posGlobal);
	TVector3 posGlobalVec(posGlobal[0],posGlobal[1],posGlobal[2]);
	return posGlobalVec;
}


void EUTelHit::setPosition(const double * position){
    _position[0] = position[0];
    _position[1] = position[1];
    _position[2] = position[2];


}
void EUTelHit::setLocation(int location){
    _locationKnown=true;
    _location = location;
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
    output.push_back(getLocation());

    return output;
}
void EUTelHit::setTrackFromLCIOVec(std::vector<double> input){
    setID(input.at(0));
    const double pos[3] = {input.at(1),input.at(2),input.at(3)};
    setPosition(pos);
    setLocation(input.at(4));
}
void EUTelHit::print(){
	streamlog_out(DEBUG1)<< std::scientific << "Hit Information:" << std::endl;
	streamlog_out(DEBUG1)<< std::scientific <<"Position Local:  " << _position[0] << " " << _position[1] <<" " << _position[2] <<" Global: " << getPositionGlobal()[0] <<" " <<getPositionGlobal()[1] << " " << getPositionGlobal()[2] <<std::endl;
	streamlog_out(DEBUG1)<< std::scientific <<"ID:  " << _id  <<std::endl;
	streamlog_out(DEBUG1)<< std::scientific <<"Location:  " << _location  <<std::endl;

}
bool EUTelHit::operator==(const EUTelHit compareHit ) const {
	if(getLocation() == compareHit.getLocation() and 	getPosition()[0] == compareHit.getPosition()[0] and	getPosition()[1] == compareHit.getPosition()[1] and 	getPosition()[2] == compareHit.getPosition()[2]){
		return true;
	}else{
		return false;
	}
}

