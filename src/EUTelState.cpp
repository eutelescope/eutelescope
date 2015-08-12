#include "EUTelState.h"
#include "EUTelNav.h"

using namespace eutelescope;
EUTelState::EUTelState():
    _arcLength(0),
    _kinks(2),
    _kinksMedium1(2),
    _kinksMedium2(2),
    _cov(5,5)
{
    _kinks[0]=0;
    _kinks[1]=0;
    _kinksMedium1[0]=0;
    _kinksMedium1[1]=0;
    _kinksMedium2[0]=0;
    _kinksMedium2[1]=0;
    _cov.Zero();
    _stateHasHit = false;
} 

EUTelState::EUTelState(EUTelState *state):
    _arcLength(0),
    _kinks(2),
    _kinksMedium1(2),
    _kinksMedium2(2),
    _cov(5,5)
{
     _kinks = state->getKinks();
    _kinksMedium1 = state->getKinksMedium1();
    _kinksMedium2 = state->getKinksMedium2();
    _stateHasHit = false;
	setDimensionSize(state->getDimensionSize());
    setArcLengthToNextState(state->getArcLengthToNextState());
	setLocation(state->getLocation());  
	double position[] = {state->getPosition()[0],state->getPosition()[1],state->getPosition()[2]}; 
	setPositionLocal(position); 	
	setDirLocalX(state->getDirLocalX()); 
	setDirLocalY(state->getDirLocalY());    
	setDirLocalZ(state->getDirLocalZ());  
    setRadFrac(state->getRadFracSensor(), state->getRadFracAir());
    setCov(state->getCov());
	if(state->getStateHasHit()){
		setHit(state->getHit());
	}
    block = state->block;
}
//getters
double EUTelState::getRadFracAir() const{
    return _radFracAir;
}

double EUTelState::getRadFracSensor() const {
    return _radFracSensor;
}

EUTelHit& EUTelState::getHit(){
    if(_stateHasHit){
        return _hit;
    }else{
        throw(std::string("Trying to access a hit which is not there "));
    }
}
EUTelHit EUTelState::getHitCopy() const {
    if(_stateHasHit){
        return _hit;
    }else{
        throw(std::string("Trying to access a hit which is not there "));
    }
}

int EUTelState::getDimensionSize() const {
	return _dimension;
}
int EUTelState::getLocation() const {
	return _location;
}
const double* EUTelState::getPosition() const {
	return &_position[0];
}
TVector3 EUTelState::getDirLocal() const {
    TVector3 vecLoc;
    vecLoc[0] = _dirLocalX;
    vecLoc[1] = _dirLocalY;
    vecLoc[2] = _dirLocalZ;
	return vecLoc;
}
TVector3 EUTelState::getDirGlobal() const {
    const double dirLoc[3] = { _dirLocalX ,_dirLocalY, _dirLocalZ};
	double dirGlo[3];
	geo::gGeometry().local2MasterVec(getLocation(), dirLoc, dirGlo );
    TVector3 vecLoc;
    vecLoc[0] = dirGlo[0];
    vecLoc[1] = dirGlo[1];
    vecLoc[2] = dirGlo[2];
	return vecLoc;
}
Eigen::Vector3d EUTelState::getDirGlobalEig() const {
	streamlog_out( DEBUG1 ) << "EUTelState::getDirGlobalEig()------------------------BEGIN" << std::endl;

    Eigen::Vector3d dirEigen;
    TVector3 dir = getDirGlobal();
    dirEigen << dir[0] , dir[1] , dir[2];
	streamlog_out( DEBUG1 ) << "EUTelState::getDirGlobalEig()------------------------END" << std::endl;

    return dirEigen;
}


TVector3 EUTelState::getPositionGlobal() const {
	const double* local =  getPosition();
	const double posLocal[3] = {local[0],local[1],local[2]};
  double posGlobal[3];
	geo::gGeometry().local2Master(getLocation() ,posLocal,posGlobal);
	TVector3 posGlobalVec(posGlobal[0],posGlobal[1],posGlobal[2]);
	return posGlobalVec;
}

Eigen::Vector3d EUTelState::getPositionGlobalEig() const {
	streamlog_out( DEBUG1 ) << "EUTelState::getDirGlobalEig()------------------------BEGIN" << std::endl;

    Eigen::Vector3d posE;
    TVector3 pos = getPositionGlobal();
    posE << pos[0] , pos[1] , pos[2];
	streamlog_out( DEBUG1 ) << "EUTelState::getDirGlobalEig()------------------------END" << std::endl;

    return posE;
}


TVectorD EUTelState::getStateVec(){ 
	streamlog_out( DEBUG1 ) << "EUTelState::getTrackStateVec()------------------------BEGIN" << std::endl;
	TVectorD stateVec(5);
	stateVec[0] = 1.0;//This has been moved to track object.
	stateVec[1] = getSlopeX();
	stateVec[2] = getSlopeY(); 
	stateVec[3] = getPosition()[0]; 
	stateVec[4] = getPosition()[1];
			
	if(stateVec[0] == INFINITY or stateVec[0] == NAN ){
		throw(lcio::Exception("Passing a state vector where curvature is not defined")); 
	}

	streamlog_out( DEBUG1 ) << "EUTelState::getTrackStateVec()------------------------END" << std::endl;
 	return stateVec;
}
TMatrixDSym EUTelState::getScatteringVarianceInLocalFrame(double const& var ){
	streamlog_out( DEBUG1 ) << "EUTelState::getScatteringVarianceInLocalFrame(Sensor)----------------------------BEGIN" << std::endl;
	if(var == 0){
		throw(std::string("Variance passed is zero."));
	}
	float scatPrecision = 1.0/var;
    TMatrixD TRotMatrix = geo::gGeometry().getRotMatrix( getLocation() );
	TVector3 axisX;
	TVector3 axisY;
	axisX[0] = TRotMatrix[0][0];	axisY[0] = TRotMatrix[0][1];
	axisX[1] = TRotMatrix[1][0];	axisY[1] = TRotMatrix[1][1];
	axisX[2] = TRotMatrix[2][0];	axisY[2] = TRotMatrix[2][1];
	float c1 = getDirGlobal().Dot(axisX); float c2 = getDirGlobal().Dot(axisY);
	streamlog_out( DEBUG1 ) << "The component in the x/y direction: "<< c1 <<"  "<<c2 << std::endl;
	TMatrixDSym precisionMatrix(2);
	float factor = scatPrecision/pow((1-pow(c1,2)-pow(c2,2)),2);
	streamlog_out( DEBUG1 ) << "The factor: "<< factor << std::endl;
	precisionMatrix[0][0]=factor*(1-pow(c2,2));precisionMatrix[0][1]=factor*c1*c2;
    precisionMatrix[1][0]=factor*c1*c2;				precisionMatrix[1][1]=factor*(1-pow(c1,2));
	streamlog_out( DEBUG1 ) << "EUTelState::getScatteringVarianceInLocalFrame(Sensor)----------------------------END" << std::endl;
	return precisionMatrix;
}
bool EUTelState::getStateHasHit() const {
    return _stateHasHit;
}

///Global -> local(Measurement)
///Calculate local to global and then invert.
TMatrixD EUTelState::getProjectionMatrix() const {
	TMatrixD xyDir(2, 3);
	xyDir[0][0] = 1; xyDir[0][1]=0.0; xyDir[0][2]=-1.0*getSlopeXGlobal();  
	xyDir[1][0] = 0; xyDir[1][1]=1.0; xyDir[1][2]=-1.0*getSlopeYGlobal();  
    //Rotation needed from local->global
    TMatrixD TRotMatrix = geo::gGeometry().getRotMatrix( getLocation() );
	TMatrixD measDir(3,2);
	measDir[0][0] = TRotMatrix[0][0];	measDir[0][1] = TRotMatrix[0][1];
	measDir[1][0] = TRotMatrix[1][0];	measDir[1][1] = TRotMatrix[1][1];
	measDir[2][0] = TRotMatrix[2][0];	measDir[2][1] = TRotMatrix[2][1];
    streamlog_out( DEBUG0 ) << "CALCULATE LOCAL TO GLOBAL STATE TRANSFORMATION... " << std::endl;
    streamlog_out( DEBUG0 ) << "Inputs... " << std::endl;
    streamlog_out( DEBUG0 ) << "The (X,Y)-axis of the global frame relative to the local  " << std::endl;
    streamlog_message( DEBUG0, measDir.Print();, std::endl; );
    streamlog_out( DEBUG0 ) << "The propagator  (Dx,Dy)   " << std::endl;
    streamlog_message( DEBUG0, xyDir.Print();, std::endl; );
	TMatrixD proM2l(2,2);
	proM2l = (xyDir*measDir).Invert(); 
    return proM2l;
}
float EUTelState::getSlopeX() const {return _dirLocalX/_dirLocalZ;}
float EUTelState::getSlopeY() const {return _dirLocalY/_dirLocalZ;}

float EUTelState::getSlopeXGlobal() const {return getDirGlobal()[0]/getDirGlobal()[2];}
float EUTelState::getSlopeYGlobal() const {return getDirGlobal()[1]/getDirGlobal()[2];}


TVectorD EUTelState::getKinks() const {
	return _kinks;
}
TVectorD EUTelState::getKinksMedium1() const  {
	return _kinksMedium1;
}
TVectorD EUTelState::getKinksMedium2() const {
	return _kinksMedium2;
}

//setters
void EUTelState::setHit(EUTelHit hit){
    _stateHasHit=true;
    _hit = hit;
}
//Can set the hit using EUTelHit or LCIO hit.
void EUTelState::setHit(EVENT::TrackerHit* hit){
    _stateHasHit=true;
    _hit = EUTelHit(hit);
}

void EUTelState::setDimensionSize(int dimension){
    _dimension = dimension;
}
void EUTelState::setLocation(int location){
    _location = location;
}
void EUTelState::setDirFromGloSlope(std::vector<double> slopes){
    double incX = slopes.at(0);
    double incY = slopes.at(1);
    double norm = sqrt(pow(incX,2)+pow(incX,2) + 1);
    const double dirGlo[3] = {incX/norm  ,incY/norm, 1.0/norm};
	double dirLoc[3];
	geo::gGeometry().master2LocalVec(getLocation(), dirGlo, dirLoc );
	setDirLocalX(dirLoc[0]);
	setDirLocalY(dirLoc[1]);
	setDirLocalZ(dirLoc[2]);
}

void EUTelState::setDirFromLocSlope(std::vector<double> slopes){
    double incX = slopes.at(0);
    double incY = slopes.at(1);
    double norm = sqrt(pow(incX,2)+pow(incX,2) + 1);
    double dir[3] = {incX/norm,incY/norm, 1.0/norm};
	setDirLocalX(dir[0]);
	setDirLocalY(dir[1]);
	setDirLocalZ(dir[2]);
}

void EUTelState::setDirLocalX(double dirX){
    _dirLocalX = dirX;
}

//Note this is dy/dz in the LOCAL frame
void EUTelState::setDirLocalY(double dirY){
    _dirLocalY = dirY;
}	
//Note this is the dx/dz in the LOCAL frame
void EUTelState::setDirLocalZ(double dirZ){
    _dirLocalZ = dirZ;
}
TMatrixD EUTelState::getCov(){
    return _cov;
}

//This variable is the RESIDUAL (Measurements - Prediction) of the kink angle. 
//Our measurement is assumed 0 in all cases.
void EUTelState::setKinks(TVectorD kinks){
    _kinks = kinks;
}
void EUTelState::setKinksMedium1(TVectorD kinks){
    _kinksMedium1 = kinks;
}
void EUTelState::setKinksMedium2(TVectorD kinks){
    _kinksMedium2 = kinks;
}
void EUTelState::setPositionGlobal(double positionGlobal[]){
    float pos[] = { static_cast<float>(positionGlobal[0]), static_cast<float>(positionGlobal[1]), static_cast<float>(positionGlobal[2]) };
    setPositionGlobal(pos);
}


void EUTelState::setPositionGlobal(float positionGlobal[]){
	double localPosition [3];
	const double referencePoint[]	= {positionGlobal[0], positionGlobal[1],positionGlobal[2]};//Need this since geometry works with const doubles not floats 
	geo::gGeometry().master2Local( getLocation(), referencePoint, localPosition );
	float posLocal[] = { static_cast<float>(localPosition[0]), static_cast<float>(localPosition[1]), static_cast<float>(localPosition[2]) };
    _position[0] = posLocal[0];
    _position[1] = posLocal[1];
    _position[2] = posLocal[2];
    if(abs(_position[2])  >  1e-7){
//        streamlog_out(MESSAGE5) << "The local z position which is being set is not zero! " << std::endl;
 //       streamlog_out(MESSAGE5) << "Global position: "<< positionGlobal[0] <<","<< positionGlobal[1] <<"," <<  positionGlobal[2] << "  ID " << this->getLocation() << std::endl;
  //      streamlog_out(MESSAGE5) << "Local position: "<< _position[0] <<","<< _position[1] <<"," <<  _position[2] << std::endl;
        throw(std::string("The corrections increase the local z postition"));
    }
}
void EUTelState::setLocalDirGlobalDir(TVector3 dirIn){
	//Now calculate the dir in LOCAL coordinates.
	const double dir[]	= {dirIn[0], dirIn[1],dirIn[2]};//Need this since geometry works with const doubles not floats 
	double localDir[3];
	geo::gGeometry().master2LocalVec(getLocation(), dir, localDir );
//	streamlog_out(DEBUG5) << "The local dir (x,y,z) is: "<< localDir[0]<<","<< localDir[1] <<"," <<localDir[2] << std::endl;
	setDirLocalX(localDir[0]);
	setDirLocalY(localDir[1]);
	setDirLocalZ(localDir[2]);

}

void EUTelState::setCov(TMatrixD cov){
    _cov = cov;
}

void EUTelState::setStateUsingCorrection(TVectorD corrections){
    /// Get position in global frame.
    Eigen::Vector3d  gPosCorr;
    gPosCorr << corrections[3] , corrections[4] , 0; 
//    std::cout<< "Corrections for global XY " << gPosCorr[0]  << " " << gPosCorr[1] << " " <<gPosCorr[2] << " Sensor " << this->getLocation()  <<std::endl;
    Eigen::Vector3d normal  = geo::gGeometry().siPlaneNormalEig(this->getLocation());
    if(normal[2] < 0){
        normal = -1*normal;
    }
    Eigen::Vector3d dir    = this->getDirGlobalEig();
    double factor =  dir.transpose()*normal;
    Eigen::MatrixXd I(3,3);
    I.setIdentity();
    Eigen::MatrixXd drldm = I - (dir*normal.transpose())*(1.0/factor);  
    Eigen::Matrix3d rot  =  geo::gGeometry().getRotMatrixEig(this->getLocation());
    Eigen::Matrix3d rotInv = rot.transpose();
    ///operates on global change in position to get change on plane defined in global frame. 
    Eigen::Vector3d corrE = drldm*gPosCorr; 
//    std::cout<< "Correction on plane " << corrE[0]  << " " << corrE[1] << " " <<corrE[2] << " Sensor " << this->getLocation()  <<std::endl;
    Eigen::Vector3d  gPos = this->getPositionGlobalEig();
    gPos = gPos + gPosCorr;
    /// Add correction in global frame and transform back to local internally 
	double corr[3] = { gPos[0],gPos[1],gPos[2]};
    setPositionGlobal(corr);
    /// slopes in global frame corrected and transformed back
    std::vector<double> slopes;
    slopes.push_back(getSlopeXGlobal() + corrections[1]);
    slopes.push_back(getSlopeYGlobal() + corrections[2]);
    setDirFromGloSlope(slopes);
}
void EUTelState::setRadFrac(double plane, double air){
    _radFracAir = air;
    _radFracSensor = plane;
}
void EUTelState::setPositionLocal(float position[]){
    double pos[3] = {position[0],position[1],position[2]};
    setPositionLocal(static_cast<double*>(pos));
}

void EUTelState::setPositionLocal(double position[]){
    float pos[3] = {position[0],position[1],position[2]};
    _position[0] = pos[0];
    _position[1] = pos[1];
    _position[2] = pos[2];
}
//print
void EUTelState::print(){
	streamlog_out(DEBUG1)<< std::scientific << "STATE VECTOR:" << std::endl;
	streamlog_out(DEBUG1)<< std::scientific <<"State memory location "<< this << " The location  " <<getLocation() <<" Distance to next state: " <<getArcLengthToNextState() <<std::endl;
	streamlog_out(DEBUG1)<< std::scientific <<"(Radiation fraction of full system)*(track Variance)    Plane: "<< getRadFracSensor() <<" Air:  " <<getRadFracAir() <<std::endl;

    streamlog_out(DEBUG1)<< std::scientific <<"Position local (X,Y,Z): "<< getPosition()[0] << " " <<  getPosition()[1]<< " " <<  getPosition()[2]<<" Global: "<<getPositionGlobal()[0]<<"  "<<getPositionGlobal()[1]<<" "<<getPositionGlobal()[2]<<std::endl;
    streamlog_out(DEBUG1)<< std::scientific <<"Direntum local (X,Y,Z): "<< getDirLocal()[0] << " " << getDirLocal()[1] << " " << getDirLocal()[2]<<" Global: "<< getDirGlobal()[0]<<" "<<getDirGlobal()[1]<<" "<<getDirGlobal()[2] <<std::endl;
    streamlog_out(DEBUG0)<< std::scientific <<"Incidence local (dx/dz,dy/dz): "<< getSlopeX()<< " " <<  getSlopeY()<< " Global: " << getSlopeXGlobal() << "  " << getSlopeYGlobal()  <<std::endl;
	streamlog_out(DEBUG0)<< std::scientific <<"Kinks local (d(dx/dz),d(dy/dz)) "<< _kinks[0] <<" ,  " << _kinks[1]<<std::endl;
	streamlog_out(DEBUG0)<< std::scientific <<"Kinks local Medium1 (d(dx/dz),d(dy/dz)) "<< _kinksMedium1[0] <<" ,  " << _kinksMedium1[1]<<std::endl;
	streamlog_out(DEBUG0)<< std::scientific <<"Kinks local Medium2 (d(dx/dz),d(dy/dz)) "<< _kinksMedium2[0] <<" ,  " << _kinksMedium2[1]<<std::endl;

	if(getStateHasHit()){
		streamlog_out(DEBUG1) <<"This state has hit: " << std::endl;
        getHit().print();

	}else{
		streamlog_out(DEBUG1) <<"This state has no hit " << std::endl;
	}
}	
//clear all contain so can be reused. This is needed in pattern recognition. 
void EUTelState::clear(){
    _kinks[0]=0;
    _kinks[1]=0;
    _stateHasHit = false;
   //Remove hit information.
	setDimensionSize(0);
    setArcLengthToNextState(0);
	setLocation(0);  
	double position[] = {0,0,0}; 
	setPositionLocal(position); 	
	setDirLocalX(0); 
	setDirLocalY(0);    
	setDirLocalZ(0);  
    setRadFrac(0,0);
}

//Overload operators.
bool EUTelState::operator<(const EUTelState compareState ) const {
	return getPosition()[2]<compareState.getPosition()[2];
}

bool EUTelState::operator==(const EUTelState compareState ) const {
	if(getLocation() == compareState.getLocation() and 	getPosition()[0] == compareState.getPosition()[0] and	getPosition()[1] == compareState.getPosition()[1] and 	getPosition()[2] == compareState.getPosition()[2]){
		return true;
	}else{
		return false;
	}
}
bool EUTelState::operator!=(const EUTelState compareState ) const {
	if(getLocation() == compareState.getLocation() and 	getPosition()[0] == compareState.getPosition()[0] and	getPosition()[1] == compareState.getPosition()[1] and 	getPosition()[2] == compareState.getPosition()[2]){
		return false;
	}else{
		return true;
	}
}


std::vector<double> EUTelState::getLCIOOutput(){
    std::vector<double> output;
    output.push_back(getDimensionSize());
    output.push_back(getLocation());
    output.push_back(_dirLocalX);
    output.push_back(_dirLocalY);
    output.push_back(_dirLocalZ);
    output.push_back(getArcLengthToNextState());
    output.push_back(getPosition()[0]);
    output.push_back(getPosition()[1]);
    output.push_back(getPosition()[2]);
    if(getStateHasHit()){
        output.push_back(1);
    }else{
        output.push_back(0);
    }
    output.push_back(getKinks()[0]);
    output.push_back(getKinks()[1]);
    output.push_back(getRadFracAir());
    output.push_back(getRadFracSensor());
		//	TMatrixDSym getStateCov() const;
    output.push_back(getKinksMedium1()[0]);
    output.push_back(getKinksMedium1()[1]);
    output.push_back(getKinksMedium2()[0]);
    output.push_back(getKinksMedium2()[1]);

    return output;

}
void EUTelState::setTrackFromLCIOVec(std::vector<double> input){
    setDimensionSize(input.at(0));
    setLocation(input.at(1));
    setDirLocalX(input.at(2));
    setDirLocalY(input.at(3));
    setDirLocalZ(input.at(4));
    setArcLengthToNextState(input.at(5)); 
    double pos[3] = {input.at(6),input.at(7),input.at(8)};
    setPositionLocal(pos);
    TVectorD kinks(2);
    kinks[0] = input.at(10);
    kinks[1] = input.at(11);
    setKinks(kinks);
    setRadFrac(input.at(13), input.at(12));
    TVectorD kinksMedium1(2);
    kinksMedium1[0] = input.at(14);
    kinksMedium1[1] = input.at(15);
    setKinksMedium1(kinksMedium1);
    TVectorD kinksMedium2(2);
    kinksMedium2[0] = input.at(16);
    kinksMedium2[1] = input.at(17);
    setKinksMedium2(kinksMedium2);

}
