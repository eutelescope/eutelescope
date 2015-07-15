#include "EUTelState.h"
#include "EUTelNav.h"

using namespace eutelescope;
EUTelState::EUTelState()
{
     _kinks.ResizeTo(2);
    _kinks[0]=0;
    _kinks[1]=0;
    _kinksMedium1.ResizeTo(2);
    _kinksMedium1[0]=0;
    _kinksMedium1[1]=0;
    _kinksMedium2.ResizeTo(2);
    _kinksMedium2[0]=0;
    _kinksMedium2[1]=0;


    _stateHasHit = false;
} 

EUTelState::EUTelState(EUTelState *state){
     _kinks.ResizeTo(2);
     _kinks = state->getKinks();
    _kinksMedium1.ResizeTo(2);
    _kinksMedium1 = state->getKinksMedium1();
    _kinksMedium2.ResizeTo(2);
    _kinksMedium2 = state->getKinksMedium2();
    _stateHasHit = false;
	setDimensionSize(state->getDimensionSize());
    setArcLengthToNextState(state->getArcLengthToNextState());
	setLocation(state->getLocation());  
	double position[] = {state->getPosition()[0],state->getPosition()[1],state->getPosition()[2]}; 
	setPositionLocal(position); 	
	setMomLocalX(state->getMomLocalX()); 
	setMomLocalY(state->getMomLocalY());    
	setMomLocalZ(state->getMomLocalZ());  
    setRadFrac(state->getRadFracSensor(), state->getRadFracAir());
	if(state->getStateHasHit()){
		setHit(state->getHit());
	}
}
//getters
double EUTelState::getRadFracAir() const{
    return _radFracAir;
}

double EUTelState::getRadFracSensor() const {
    return _radFracSensor;
}

EUTelHit EUTelState::getHit(){
	return _hit;
}
int EUTelState::getDimensionSize() const {
	return _dimension;
}
int EUTelState::getLocation() const {
	return _location;
}
const float* EUTelState::getPosition() const {
	return &_position[0];
}
TVector3 EUTelState::getPositionGlobal() const {
	const float* local =  getPosition();
	const double posLocal[3] = {local[0],local[1],local[2]};
  double posGlobal[3];
	geo::gGeometry().local2Master(getLocation() ,posLocal,posGlobal);
	TVector3 posGlobalVec(posGlobal[0],posGlobal[1],posGlobal[2]);
	return posGlobalVec;
}
TVectorD EUTelState::getStateVec(){ 
	streamlog_out( DEBUG1 ) << "EUTelState::getTrackStateVec()------------------------BEGIN" << std::endl;
	TVectorD stateVec(5);
	stateVec[0] = -1.0/getMomLocal().Mag();
	stateVec[1] = getMomLocalX()/getMomLocalZ();
	stateVec[2] = getMomLocalY()/getMomLocalZ(); 
	stateVec[3] = getPosition()[0]; 
	stateVec[4] = getPosition()[1];
			
	if(stateVec[0] == INFINITY or stateVec[0] == NAN ){
		throw(lcio::Exception("Passing a state vector where curvature is not defined")); 
	}

	streamlog_out( DEBUG1 ) << "EUTelState::getTrackStateVec()------------------------END" << std::endl;
 	return stateVec;
}
TMatrixDSym EUTelState::getScatteringVarianceInLocalFrame(){
	streamlog_out( DEBUG1 ) << "EUTelState::getScatteringVarianceInLocalFrame(Sensor)----------------------------BEGIN" << std::endl;
	streamlog_out(DEBUG1) << "Variance (Sensor):  " << std::scientific << getRadFracSensor() << "  Plane: " << getLocation()  << std::endl;
	if(getRadFracSensor() == 0){
		throw(std::string("Radiation of sensor is zero. Something is wrong with radiation length calculation."));
	}
	float scatPrecision = 1.0/getRadFracSensor();
	//We need the track direction in the direction of x/y in the local frame. 
	//This will be the same as unitMomentum in the x/y direction
	TVector3 unitMomentumLocalFrame =	getMomLocal().Unit();
	//c1 and c2 come from Claus's paper GBL
	float c1 = 	unitMomentumLocalFrame[0]; float c2 =	unitMomentumLocalFrame[1];
	streamlog_out( DEBUG1 ) << "The component in the x/y direction: "<< c1 <<"  "<<c2 << std::endl;
	TMatrixDSym precisionMatrix(2);
	float factor = scatPrecision/pow((1-pow(c1,2)-pow(c2,2)),2);
	streamlog_out( DEBUG1 ) << "The factor: "<< factor << std::endl;
	precisionMatrix[0][0]=factor*(1-pow(c2,2));
  precisionMatrix[1][0]=factor*c1*c2;				precisionMatrix[1][1]=factor*(1-pow(c1,2));
	streamlog_out( DEBUG1 ) << "EUTelState::getScatteringVarianceInLocalFrame(Sensor)----------------------------END" << std::endl;
	return precisionMatrix;
}
TMatrixDSym EUTelState::getScatteringVarianceInLocalFrame(float  variance){
	streamlog_out( DEBUG1 ) << "EUTelState::getScatteringVarianceInLocalFrame(Scatter)----------------------------BEGIN" << std::endl;
	streamlog_out(DEBUG5)<<"Variance (AIR Fraction): " <<std::scientific  <<  variance <<std::endl; 
	float scatPrecision = 1.0 /variance;
	//We need the track direction in the direction of x/y in the local frame. 
	//This will be the same as unitMomentum in the x/y direction
	TVector3 unitMomentumLocalFrame =	getMomLocal().Unit();
	//c1 and c2 come from Claus's paper GBL
	float c1 = 	unitMomentumLocalFrame[0]; float c2 = 	unitMomentumLocalFrame[1];
	streamlog_out( DEBUG1 ) << "The component in the x/y direction: "<< c1 <<"  "<<c2 << std::endl;
	TMatrixDSym precisionMatrix(2);
	float factor = scatPrecision/pow((1-pow(c1,2)-pow(c2,2)),2);
	streamlog_out( DEBUG1 ) << "The factor: "<< factor << std::endl;
	precisionMatrix[0][0]=factor*(1-pow(c2,2));
  precisionMatrix[1][0]=factor*c1*c2;				precisionMatrix[1][1]=factor*(1-pow(c1,2));
	streamlog_out( DEBUG1 ) << "EUTelState::getScatteringVarianceInLocalFrame(Scatter)----------------------------END" << std::endl;

	return precisionMatrix;
}
TMatrixDSym EUTelState::getStateCov() const {

//	streamlog_out( DEBUG1 ) << "EUTelState::getTrackStateCov()----------------------------BEGIN" << std::endl;
	TMatrixDSym C(5);   
//	const EVENT::FloatVec& trkCov = getCovMatrix();        
	C.Zero();
            
//	C[0][0] = trkCov[0]; 
//	C[1][0] = trkCov[1];  C[1][1] = trkCov[2]; 
//	C[2][0] = trkCov[3];  C[2][1] = trkCov[4];  C[2][2] = trkCov[5]; 
//	C[3][0] = trkCov[6];  C[3][1] = trkCov[7];  C[3][2] = trkCov[8];  C[3][3] = trkCov[9]; 
//	C[4][0] = trkCov[10]; C[4][1] = trkCov[11]; C[4][2] = trkCov[12]; C[4][3] = trkCov[13]; C[4][4] = trkCov[14]; 
//        
//	if ( streamlog_level(DEBUG0) ){
//		streamlog_out( DEBUG0 ) << "Track state covariance matrix:" << std::endl;
//		C.Print();
//	}
//        
	return C;
//	streamlog_out( DEBUG1 ) << "EUTelState::getTrackStateCov()----------------------------END" << std::endl;
}
bool EUTelState::getStateHasHit() const {
    return _stateHasHit;
}

void EUTelState::getCombinedHitAndStateCovMatrixInLocalFrame( double (&cov)[4] ) const {
	cov[0] = _covCombinedMatrix[0];
	cov[1] = _covCombinedMatrix[1];
	cov[2] = _covCombinedMatrix[2];
	cov[3] = _covCombinedMatrix[3];
}
TMatrixD EUTelState::getProjectionMatrix() const {
	TMatrixD projection(5,5);
	projection.Zero();
	TMatrixD proM2l(2, 2);
	proM2l.UnitMatrix();
	projection.SetSub(3, 3, proM2l);
	return proM2l;
}
TVector3 EUTelState::getMomLocal(){
	TVector3 pVecUnitLocal;
	pVecUnitLocal[0] = getMomLocalX(); 	pVecUnitLocal[1] = getMomLocalY(); 	pVecUnitLocal[2] = getMomLocalZ(); 
	return pVecUnitLocal;
}
TVectorD EUTelState::getKinks() const {
	return _kinks;
}
TVectorD EUTelState::getKinksMedium1() const {
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
    const double * pos = hit->getPosition();
    _hit.setPosition(pos);
    _hit.setID(hit->id());
}

void EUTelState::setDimensionSize(int dimension){
    _dimension = dimension;
}
void EUTelState::setLocation(int location){
    _location = location;
}
void EUTelState::setMomLocalX(float momX){
    _momLocalX = momX;
}

//Note this is dy/dz in the LOCAL frame
void EUTelState::setMomLocalY(float momY){
    _momLocalY = momY;
}	
//Note this is the dx/dz in the LOCAL frame
void EUTelState::setMomLocalZ(float momZ){
    _momLocalZ = momZ;
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

void EUTelState::setPositionGlobal(float positionGlobal[]){
	double localPosition [3];
	const double referencePoint[]	= {positionGlobal[0], positionGlobal[1],positionGlobal[2]};//Need this since geometry works with const doubles not floats 
	geo::gGeometry().master2Local( getLocation(), referencePoint, localPosition );
	float posLocal[] = { static_cast<float>(localPosition[0]), static_cast<float>(localPosition[1]), static_cast<float>(localPosition[2]) };
    _position[0] = posLocal[0];
    _position[1] = posLocal[1];
    _position[2] = posLocal[2];


}
void EUTelState::setLocalMomentumGlobalMomentum(TVector3 momentumIn){
	//Now calculate the momentum in LOCAL coordinates.
	const double momentum[]	= {momentumIn[0], momentumIn[1],momentumIn[2]};//Need this since geometry works with const doubles not floats 
	double localMomentum [3];
	geo::gGeometry().master2LocalVec(getLocation(), momentum, localMomentum );
//	streamlog_out(DEBUG5) << "The local momentum (x,y,z) is: "<< localMomentum[0]<<","<< localMomentum[1] <<"," <<localMomentum[2] << std::endl;
	setMomLocalX(localMomentum[0]);
	setMomLocalY(localMomentum[1]);
	setMomLocalZ(localMomentum[2]);

}

void EUTelState::setCombinedHitAndStateCovMatrixInLocalFrame(double cov[4]){
	_covCombinedMatrix[0] = cov[0];
	_covCombinedMatrix[1] = cov[1];
	_covCombinedMatrix[2] = cov[2];
	_covCombinedMatrix[3] = cov[3];
}

void EUTelState::setStateUsingCorrection(TVectorD corrections){
	double referencePoint[] = { getPosition()[0]+corrections[3],getPosition()[1]+corrections[4],0};
	setPositionLocal(referencePoint);
    float charge = -1.0;
    float omegaNew =  charge/getMomLocal().Mag() + corrections[0];
    float newMom = charge/omegaNew; 
 //   std::cout << "Here " << getMomLocalX() << " " << getMomLocalY()<< " " << getMomLocalZ() <<"New mom " << newMom<<std::endl;
    setMomLocalX( newMom*(getMomLocalX()/getMomLocalZ() + corrections[1]));  
    setMomLocalY( newMom*(getMomLocalY()/getMomLocalZ() + corrections[2]));  
    setMomLocalZ(sqrt(pow(newMom,2) - pow(getMomLocalX(),2) - pow(getMomLocalY(),2)));

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

//find
bool EUTelState::findIntersectionWithCertainID(int nextSensorID, float intersectionPoint[], TVector3& momentumAtIntersection, float& arcLength, int& newNextPlaneID )
{
	TVector3 pVec = getMomGlobal();

	//TODO: better exception
	if(pVec.Mag() == 0)
	{
		throw(lcio::Exception( "The momentum is 0")); 
	}

	double posLocal[] =  {getPosition()[0],getPosition()[1],getPosition()[2] };
	double temp[] = {0.,0.,0.};

	//IMPORTANT:For strip sensors this will make the hit strip look like a pixel at (Xstriplocal,somevalue,somevalue).
	geo::gGeometry().local2Master(getLocation() , posLocal, temp);
    int charge = -1;
	float posGlobal[] = { static_cast<float>(temp[0]), static_cast<float>(temp[1]), static_cast<float>(temp[2]) };
	return  EUTelNav::findIntersectionWithCertainID(	posGlobal[0], posGlobal[1], posGlobal[2], 
								pVec[0], pVec[1], pVec[2], charge,
								nextSensorID, intersectionPoint, 
								momentumAtIntersection, arcLength, newNextPlaneID); 
}

//compute
TVector3 EUTelState::getMomGlobal() const {
    const double input[3]={getMomLocalX(),getMomLocalY(),getMomLocalZ()};
    double momentum[3];
    geo::gGeometry().local2MasterVec(getLocation(),input,momentum);
    return TVector3(momentum[0],momentum[1],momentum[2]);
}
//TMatrix EUTelState::computePropagationJacobianFromLocalStateToNextLocalState(TVector3 momentumEnd, float arcLength,float nextPlaneID) {
//	streamlog_out(DEBUG2) << "-------------------------------EUTelState::computePropagationJacobianFromStateToThisZLocation()-------------------------BEGIN" << std::endl;
//	if(arcLength == 0 or arcLength < 0 ){ 
//		throw(lcio::Exception( "The arc length is less than or equal to zero.")); 
//	}
//	TMatrixD curvilinearJacobian = EUTelNav::getPropagationJacobianCurvilinear(arcLength,getOmega(), computeCartesianMomentum().Unit(),momentumEnd.Unit());
//	streamlog_out(DEBUG0)<<"This is the curvilinear jacobian at sensor:" << std::endl; 
//	streamlog_message( DEBUG0, curvilinearJacobian.Print();, std::endl; );
//	streamlog_out(DEBUG0)<<"The state vector that create the curvilinear system is:" << std::endl; 
//	print();
//	TMatrixD localToCurvilinearJacobianStart =  EUTelNav::getLocalToCurvilinearTransformMatrix(computeCartesianMomentum(),getLocation() ,getBeamCharge() );
//	streamlog_out(DEBUG0)<<"This is the local to curvilinear jacobian at sensor : " << std::endl; 
//	streamlog_message( DEBUG0, localToCurvilinearJacobianStart.Print();, std::endl; );
//	TMatrixD localToCurvilinearJacobianEnd =  EUTelNav::getLocalToCurvilinearTransformMatrix(momentumEnd,nextPlaneID ,getBeamCharge() );
//	streamlog_out(DEBUG0)<<"This is the local to curvilinear jacobian at sensor at last next sensor : " << std::endl; 
//	streamlog_message( DEBUG0, localToCurvilinearJacobianEnd.Print();, std::endl; );
//	TMatrixD curvilinearToLocalJacobianEnd = localToCurvilinearJacobianEnd.Invert();
//	streamlog_out(DEBUG0)<<"This is the curvilinear to local jacobian at sensor : " << std::endl; 
//	streamlog_message( DEBUG0, curvilinearToLocalJacobianEnd.Print();, std::endl; );
//	TMatrixD localToNextLocalJacobian = curvilinearToLocalJacobianEnd*curvilinearJacobian*localToCurvilinearJacobianStart;
//	streamlog_out(DEBUG0)<<"This is the full jacobian : "<<  std::endl; 
//	streamlog_message( DEBUG0, localToNextLocalJacobian.Print();, std::endl; );
//
//	streamlog_out(DEBUG2) << "-------------------------------EUTelState::computePropagationJacobianFromStateToThisZLocation()-------------------------END" << std::endl;
//	return localToNextLocalJacobian;
//}
//THIS WILL RETURN THE TOTAL RADIATION LENGTH OF THE TELESCOPE SYSTEM STARTING AT THE STATE AND MOVING FORWARD. 
//WE ALSO GET THE FRACTION OF RADIATION LENGTH THAT EACH PLANE AND VOLUME OF AIR SHOULD GET. 
//THIS IS ASSOCIATED SO THE AIR INFRONT OF A SENSOR IS ASSOCIATED WITH IT.
//EXCLUDED PLANES ARE REDUCED TO MORE RADIATION LENGTH IN FRONT OF A NON EXCLUDED PLANE.
float EUTelState::computeRadLengthsToEnd( std::map<const int,double> & mapSensor, std::map<const int ,double> & mapAir){
	//Get the ID of the last sensor
	int lastPlaneID; // =  	geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at( geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size()-1 );
	float intersectionPoint[3];
	TVector3 momentumAtIntersection;
	float arcLength;
	int holder; //This is used to return the plane which is found.
	//DETERMINE INTERSECTION ON LAST PLANE USING STATE INFORMATION.
	findIntersectionWithCertainID(lastPlaneID, intersectionPoint, momentumAtIntersection,arcLength,holder );
	TVector3 gPos =  getPositionGlobal();
	//TO DO: At the moment we just use a straight through the sensor in all enviroments. The code is designed to extend this to any straight line but we see some addition of extra radiation length beyond what is expect. This will have to be looked into but not a huge issue at the moment.
	//NOTE THE Z VALUE FOR THESE ARE NOT USED IN calculateTotalRadiationLengthAndWeights
	const double start[] = {gPos[0],gPos[1],-0.025+gPos[2]};
	const double end[]   = {gPos[0],gPos[1],gPos[2]+0.025};//Must make sure we add all silicon.
	//NOW WE CALCULATE THE RADIATION LENGTH FOR THE FULL FLIGHT AND THEN SPLIT THESE INTO LINEAR  COMMPONENTS FOR SCATTERING ESTIMATION. 
	//We will return the radiation lengths associate with the planes and air. Note excluded planes volume should be added to the air in front of non excluded planes. 
	float rad =	geo::gGeometry().calculateTotalRadiationLengthAndWeights( start,end,  mapSensor, mapAir);
	return rad;
}


//print
void EUTelState::print(){
	streamlog_out(DEBUG1)<< std::scientific << "STATE VECTOR:" << std::endl;
	streamlog_out(DEBUG1)<< std::scientific <<"State memory location "<< this << " The location  " <<getLocation() <<" Distance to next state: " <<getArcLengthToNextState() <<std::endl;
	streamlog_out(DEBUG1)<< std::scientific <<"(Radiation fraction of full system)*(track Variance)    Plane: "<< getRadFracSensor() <<" Air:  " <<getRadFracAir() <<std::endl;

    streamlog_out(DEBUG1)<< std::scientific <<"Position local (X,Y,Z): "<< getPosition()[0] << " " <<  getPosition()[1]<< " " <<  getPosition()[2]<<" Global: "<<getPositionGlobal()[0]<<"  "<<getPositionGlobal()[1]<<" "<<getPositionGlobal()[2]<<std::endl;
    streamlog_out(DEBUG1)<< std::scientific <<"Momentum local (X,Y,Z): "<< getMomLocal()[0] << " " << getMomLocal()[1] << " " << getMomLocal()[2]<<" Global: "<< getMomGlobal()[0]<<" "<<getMomGlobal()[1]<<" "<<getMomGlobal()[2] <<std::endl;
    streamlog_out(DEBUG1)<< std::scientific <<"Incidence local (dx/dz,dy/dz): "<< getMomLocal()[0]/getMomLocal()[2]<< " " << getMomLocal()[1]/getMomLocal()[2] <<std::endl;
	streamlog_out(DEBUG1)<< std::scientific <<"Kinks local (d(dx/dz),d(dy/dz)) "<< _kinks[0] <<" ,  " << _kinks[1]<<std::endl;
	streamlog_out(DEBUG1)<< std::scientific <<"Kinks local Medium1 (d(dx/dz),d(dy/dz)) "<< _kinksMedium1[0] <<" ,  " << _kinksMedium1[1]<<std::endl;
	streamlog_out(DEBUG1)<< std::scientific <<"Kinks local Medium2 (d(dx/dz),d(dy/dz)) "<< _kinksMedium2[0] <<" ,  " << _kinksMedium2[1]<<std::endl;

	if(getStateHasHit()){
        streamlog_out(DEBUG1)<< std::scientific<<"The hit ID of the state is "<<getHit().getID()<<std::endl;
        streamlog_out(DEBUG1)<< std::scientific<<"Position hit local (X,Y,Z): "<<getHit().getPosition()[0]<<" "<<getHit().getPosition()[1]<<getHit().getPosition()[2]<<std::endl;

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
	setMomLocalX(0); 
	setMomLocalY(0);    
	setMomLocalZ(0);  
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
    output.push_back(_momLocalX);
    output.push_back(_momLocalY);
    output.push_back(_momLocalZ);
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
    setMomLocalX(input.at(2));
    setMomLocalY(input.at(3));
    setMomLocalZ(input.at(4));
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
