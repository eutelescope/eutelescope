#include "EUTelHit.h"
using namespace eutelescope;

EUTelHit::EUTelHit():
_locationKnown(false),
_cov(2,2)
{
    _cov.Zero();
} 

EUTelHit::EUTelHit(EUTelHit* hit):
_cov(2,2)
{
    _position[0] = hit->getPosition()[0];
    _position[1] = hit->getPosition()[1];
    _position[2] = hit->getPosition()[2];
    _id = hit->getID();
    _time = hit->getTime();
    if(hit->_locationKnown){
        _location = hit->getLocation();
        _locationKnown=true;
    }else{
        _locationKnown=false;
    }
    setCov(hit->getCov());
    ///All hits must have clustering information 
    _pulse = hit->getPulse();
} 

EUTelHit::EUTelHit(EVENT::TrackerHit* hit):
_cov(2,2)
{
    _position[0] = hit->getPosition()[0];
    _position[1] = hit->getPosition()[1];
    _position[2] = hit->getPosition()[2];
    _id = hit->id();
    _time = hit->getTime();
    const int hitLoc = Utility::getSensorIDfromHit( static_cast<IMPL::TrackerHitImpl*> (hit) );
    _location = hitLoc;
    _pulse =hit->getRawHits().at(0);    
    _locationKnown=true;
    setCov(hit->getCovMatrix());

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

int EUTelHit::getID() const {
    return _id;
}
EVENT::LCObject* EUTelHit::getPulse(){
    return _pulse;
}
void EUTelHit::getCov( double (&cov)[4] ) const {
	cov[0] = _cov[0][0];
	cov[1] = _cov[0][1];
	cov[2] = _cov[1][0];
	cov[3] = _cov[1][1];
}

void EUTelHit::setCov(const std::vector<double>& cov){
	_cov[0][0] = cov.at(0);
	_cov[0][1] = cov.at(1);
	_cov[1][0] = cov.at(2);
	_cov[1][1] = cov.at(3);
}
void EUTelHit::setCov(const std::vector<float>& cov){
	_cov[0][0] = static_cast<double>(cov.at(0));
	_cov[0][1] = static_cast<double>(cov.at(1));
	_cov[1][0] = static_cast<double>(cov.at(2));
	_cov[1][1] = static_cast<double>(cov.at(3));
}

void EUTelHit::setCov(TMatrixD cov){
    _cov = cov;
}
void EUTelHit::setPulse( EVENT::LCObject* pulse){
    _pulse = pulse;
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

void EUTelHit::setID(int id){
    _id = id;
}

 float EUTelHit::getTime() const {
     return _time;
 }
 void EUTelHit::setTime(float time){
     _time = time;
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
	streamlog_out(DEBUG1)<< std::scientific <<"Location:  " << _location  << " Pulse " << _pulse <<std::endl;

}
bool EUTelHit::operator==(const EUTelHit compareHit ) const {
	if(getLocation() == compareHit.getLocation() and 	getPosition()[0] == compareHit.getPosition()[0] and	getPosition()[1] == compareHit.getPosition()[1] and 	getPosition()[2] == compareHit.getPosition()[2]){
		return true;
	}else{
		return false;
	}
}

