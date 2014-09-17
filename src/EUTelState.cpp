#include "EUTelState.h"

using namespace eutelescope;
EUTelState::EUTelState(){
setBeamCharge(-1.0);
setBeamEnergy(5.0);
} 

EUTelState::EUTelState(EUTelState *state){
	//cout<<"Entering the copy constructor"<<endl;
	setDimensionSize(state->getDimensionSize());
	setLocation(state->getLocation());//This stores the location as a float in Z0 since no location for track LCIO. This is preferable to problems with storing hits.  
	//cout<<"HEREoldstate: "<<state->getPosition()[0]<<","<<state->getPosition()[1]<<","<<state->getPosition()[2]<<","<<state->getLocation()<<endl;
	float position[] = {state->getPosition()[0],state->getPosition()[1],state->getPosition()[2]};//Z position is not a state parameter but should add here for use later. 
	setPositionLocal(position);//This will automatically take cartesian coordinate system and save it to reference point. //This is different from most LCIO applications since each track will have own reference point. 	
	//cout<<"HEREnewstate: "<<getPosition()[0]<<","<<getPosition()[1]<<","<<getPosition()[2]<<","<<getLocation()<<endl;
	setIntersectionLocalXZ(state->getIntersectionLocalXZ());    
	setIntersectionLocalYZ(state->getIntersectionLocalYZ());  
	setBeamCharge(state->getBeamCharge());//this is set for each state. to do: is there a more efficient way of doing this since we only need this stored once?
	setBeamEnergy(state->getBeamEnergy());//this is saved in gev. 
	initialiseCurvature(); //this will perform the calculation _beamq/_beame ad place in invp
	if(getIsThereAHit()){
		addHit(const_cast<EVENT::TrackerHit*>(state->getTrackerHits().at(0)));
	}
}
//getters
EVENT::TrackerHit* EUTelState::getHit(){
	return	const_cast<EVENT::TrackerHit*>(getTrackerHits().at(0));
}
int EUTelState::getDimensionSize() const {
	int dimension = static_cast<int>(getD0());
	return dimension;
}
int EUTelState::getLocation() const {
	int location = static_cast<int>(getZ0());
	return location;
}
float* EUTelState::getPosition() const {
	return const_cast<float*>(getReferencePoint());
}
TVector3 EUTelState::getPositionGlobal() const {
	float* local =  const_cast<float*>(getReferencePoint());
	const double posLocal[3] = {local[0],local[1],local[2]};
  double posGlobal[3];
	geo::gGeometry().local2Master(getLocation() ,posLocal,posGlobal);
	TVector3 posGlobalVec(posGlobal[0],posGlobal[1],posGlobal[2]);
	return posGlobalVec;
}
TVectorD EUTelState::getStateVec() const { 
	streamlog_out( DEBUG1 ) << "EUTelState::getTrackStateVec()------------------------BEGIN" << std::endl;
	TVector3 momentum =	computeCartesianMomentum();
	TVectorD stateVec(5);
	const float lambda = asin(momentum[2]/(momentum.Mag()));//This will be in radians.
	const float phi = asin(momentum[1]/(momentum.Mag())*cos(lambda));


	stateVec[0] = getOmega();
	stateVec[1] = getIntersectionLocalXZ();
	stateVec[2] = getIntersectionLocalYZ(); 
	stateVec[3] = getPosition()[0]; 
	stateVec[4] = getPosition()[1];
			
//	if ( streamlog_level(DEBUG0) ){
//		streamlog_out( DEBUG0 ) << "Track state:" << std::endl;
//		stateVec.Print();
//	}
	if(stateVec[0] == INFINITY or stateVec[0] == NAN ){
		throw(lcio::Exception( Utility::outputColourString("Passing a state vector where curvature is not defined","RED"))); 
	}

	streamlog_out( DEBUG1 ) << "EUTelState::getTrackStateVec()------------------------END" << std::endl;
 	return stateVec;
}
TMatrixDSym EUTelState::getScatteringVarianceInLocalFrame(){
	streamlog_out( DEBUG1 ) << "EUTelState::getScatteringVarianceInLocalFrame(Sensor)----------------------------BEGIN" << std::endl;
	//TO DO:Should make the reading for any planes not just si. 
	const double x0  = geo::gGeometry().siPlaneRadLength(getLocation());//This is in mm.
	const double thickness        = geo::gGeometry().siPlaneZSize(getLocation());//TO DO: Need to get the correct thickness
	streamlog_out(DEBUG5)<<"The radiation length for the sensor is: " << x0 <<std::endl; 
	if( x0== 0){
		throw(lcio::Exception(Utility::outputColourString("Radiation length is zero.", "RED"))); 	
	}
	streamlog_out(DEBUG5)<<"The thickness for the sensor is: " << thickness <<std::endl; 
	if(thickness == 0 ){ 
		throw(lcio::Exception(Utility::outputColourString("thickness is zero", "RED"))); 	
	}
	const double percentageOfRadiationLength  = thickness / x0; 
	streamlog_out(DEBUG5)<<"The  percentageOfRadiationLength for the sensor is: " <<  percentageOfRadiationLength<<std::endl; 
	if( percentageOfRadiationLength== 0){
		streamlog_out(MESSAGE0)<<"The thickness: " << thickness << "The radiation length is: "<< x0 <<std::endl; 
		throw(lcio::Exception(Utility::outputColourString("percentageOfRadiationLength is zero for the plane itself.", "RED"))); 	
	}
	if(getBeamEnergy() == 0 ){ 
		throw(lcio::Exception(Utility::outputColourString("Beam energy is zero", "RED"))); 	
	}
	const double scatteringVariance  = Utility::getThetaRMSHighland(getBeamEnergy(), percentageOfRadiationLength);
	streamlog_out(DEBUG5)<<"The scattering variance in radians: " <<  scatteringVariance <<std::endl; 

	float scatPrecisionSqrt = 1.0/scatteringVariance;
	//We need the track direction in the direction of x/y in the local frame. 
	//This will be the same as unitMomentum in the x/y direction
	TVector3 unitMomentumLocalFrame =	getIncidenceUnitMomentumVectorInLocalFrame();
	//c1 and c2 come from Claus's paper GBL
	float c1 = 	unitMomentumLocalFrame[0]; float c2 =	unitMomentumLocalFrame[1];
	streamlog_out( DEBUG1 ) << "The component in the x/y direction: "<< c1 <<"  "<<c2 << std::endl;
	TMatrixDSym precisionMatrix(2);
	float factor = pow(scatPrecisionSqrt,2)/pow((1-pow(c1,2)-pow(c2,2)),2);
	streamlog_out( DEBUG1 ) << "The factor: "<< factor << std::endl;
	precisionMatrix[0][0]=factor*(1-pow(c2,2));
  precisionMatrix[1][0]=factor*c1*c2;				precisionMatrix[1][1]=factor*(1-pow(c1,2));
	streamlog_out( DEBUG1 ) << "EUTelState::getScatteringVarianceInLocalFrame(Sensor)----------------------------END" << std::endl;
	return precisionMatrix;
}
TMatrixDSym EUTelState::getScatteringVarianceInLocalFrame(float  percentageOfRadiationLength){
	streamlog_out( DEBUG1 ) << "EUTelState::getScatteringVarianceInLocalFrame(Scatter)----------------------------BEGIN" << std::endl;
	const double scatteringVariance  = Utility::getThetaRMSHighland(getBeamEnergy(), percentageOfRadiationLength);
	streamlog_out(DEBUG5)<<"The scattering variance in radians: " <<  scatteringVariance <<std::endl; 
	float scatPrecisionSqrt = 1.0 /scatteringVariance;
	//We need the track direction in the direction of x/y in the local frame. 
	//This will be the same as unitMomentum in the x/y direction
	TVector3 unitMomentumLocalFrame =	getIncidenceUnitMomentumVectorInLocalFrame();
	//c1 and c2 come from Claus's paper GBL
	float c1 = 	unitMomentumLocalFrame[0]; float c2 = 	unitMomentumLocalFrame[1];
	streamlog_out( DEBUG1 ) << "The component in the x/y direction: "<< c1 <<"  "<<c2 << std::endl;
	TMatrixDSym precisionMatrix(2);
	float factor = pow(scatPrecisionSqrt,2)/pow((1-pow(c1,2)-pow(c2,2)),2);
	streamlog_out( DEBUG1 ) << "The factor: "<< factor << std::endl;
	precisionMatrix[0][0]=factor*(1-pow(c2,2));
  precisionMatrix[1][0]=factor*c1*c2;				precisionMatrix[1][1]=factor*(1-pow(c1,2));
	streamlog_out( DEBUG1 ) << "EUTelState::getScatteringVarianceInLocalFrame(Scatter)----------------------------END" << std::endl;

	return precisionMatrix;
}
TMatrixDSym EUTelState::getStateCov() const {

	streamlog_out( DEBUG1 ) << "EUTelState::getTrackStateCov()----------------------------BEGIN" << std::endl;
	TMatrixDSym C(5);   
	const EVENT::FloatVec& trkCov = getCovMatrix();        
	C.Zero();
            
	C[0][0] = trkCov[0]; 
	C[1][0] = trkCov[1];  C[1][1] = trkCov[2]; 
	C[2][0] = trkCov[3];  C[2][1] = trkCov[4];  C[2][2] = trkCov[5]; 
	C[3][0] = trkCov[6];  C[3][1] = trkCov[7];  C[3][2] = trkCov[8];  C[3][3] = trkCov[9]; 
	C[4][0] = trkCov[10]; C[4][1] = trkCov[11]; C[4][2] = trkCov[12]; C[4][3] = trkCov[13]; C[4][4] = trkCov[14]; 
        
	if ( streamlog_level(DEBUG0) ){
		streamlog_out( DEBUG0 ) << "Track state covariance matrix:" << std::endl;
		C.Print();
	}
        
	return C;
	streamlog_out( DEBUG1 ) << "EUTelState::getTrackStateCov()----------------------------END" << std::endl;
}
bool EUTelState::getIsThereAHit() const {
	if(!getTrackerHits().empty()){
		return true;
	}else{
		return false;
	}
}

void EUTelState::getCombinedHitAndStateCovMatrixInLocalFrame( double (&cov)[4] ) const {
	cov[0] = _covCombinedMatrix[0];
	cov[1] = _covCombinedMatrix[1];
	cov[2] = _covCombinedMatrix[2];
	cov[3] = _covCombinedMatrix[3];
}
//TO DO:This matrix will only work for no tilted sensors. Must determine the generic projection matrix
TMatrixD EUTelState::getProjectionMatrix() const {
	TMatrixD projection(5,5);
	projection.Zero();
	TMatrixD proM2l(2, 2);
	proM2l.UnitMatrix();
	projection.SetSub(3, 3, proM2l);
	return proM2l;
}
TVector3 EUTelState::getIncidenceUnitMomentumVectorInLocalFrame(){
	TVector3 pVec =	computeCartesianMomentum();
	TVector3 pVecUnit = pVec.Unit();//Make the vector unit.
    streamlog_out(DEBUG2) << "Momentum in global coordinates  Px,Py,Pz= " << pVec[0]<<","<<pVec[1]<<","<<pVec[2] << std::endl;
	double globalVec[] = { pVecUnit[0],pVecUnit[1],pVecUnit[2] };
	double localVec[3];
	geo::gGeometry().master2LocalVec( getLocation() ,globalVec, localVec );
	TVector3 pVecUnitLocal;
	pVecUnitLocal[0] = localVec[0]; 	pVecUnitLocal[1] = localVec[1]; 	pVecUnitLocal[2] = localVec[2]; 
  streamlog_out(DEBUG2) << "Momentum in local coordinates  Px,Py,Pz= " << pVecUnitLocal[0]<<","<<pVecUnitLocal[1]<<","<<pVecUnitLocal[2]<< std::endl;
	return pVecUnitLocal;
}
TVectorD EUTelState::getKinks(){
	EVENT::FloatVec kinksVec = getCovMatrix();
	TVectorD kinks(2);
	kinks(0) = kinksVec.at(0);
	kinks(1) = kinksVec.at(1);
	return kinks;
}
//setters
void EUTelState::setDimensionSize(int dimension){
	setD0(static_cast<float>(dimension));
}
void EUTelState::setLocation(int location){
	float locationFloat = static_cast<float>(location);
	setZ0(locationFloat);
}
void EUTelState::setBeamEnergy(float beamE){
	setdEdx(beamE);
}

void EUTelState::setBeamCharge(float beamQ){
	setdEdxError(beamQ);
}
//Note this is dy/dz in the LOCAL frame
void EUTelState::setIntersectionLocalYZ(float directionYZ){
	setTanLambda(directionYZ);
}	
//Note this is the dx/dz in the LOCAL frame
void EUTelState::setIntersectionLocalXZ(float directionXZ){
	setPhi(directionXZ);//This is only a container. So hold inverseTan(Phi)
}
void EUTelState::setPositionLocal(float position[]){
	setReferencePoint(position);
}
void EUTelState::setKinks(TVectorD kinks){
	EVENT::FloatVec kinksInput;
	kinksInput.push_back(kinks[0]);
	kinksInput.push_back(kinks[1]);
	setCovMatrix(kinksInput);
}
void EUTelState::setPositionGlobal(float positionGlobal[]){
	double localPosition [3];
	//TO DO:Fix master2Localtwo so that is is the only one. This is currently a hack
	const double referencePoint[]	= {positionGlobal[0], positionGlobal[1],positionGlobal[2]};//Need this since geometry works with const doubles not floats 
	geo::gGeometry().master2Localtwo( getLocation(), referencePoint, localPosition );
	float posLocal[] = { static_cast<float>(localPosition[0]), static_cast<float>(localPosition[1]), static_cast<float>(localPosition[2]) };
	setReferencePoint(posLocal);
}
//This sets the LOCAL frame intersection. Not the curvilinear frames intersection
void EUTelState::setLocalXZAndYZIntersectionAndCurvatureUsingGlobalMomentum(TVector3 momentumIn){
	streamlog_out(DEBUG5) << "EUTelState::setLocalXZAndYZIntersectionAndCurvatureUsingGlobalMomentum--------------------------BEGIN" << std::endl;

	//set the beam energy and curvature
	setBeamEnergy(momentumIn.Mag());
	initialiseCurvature();//You must set beam charge before you call this.
	//Now calculate the momentum in LOCAL coordinates.
	const double momentum[]	= {momentumIn[0], momentumIn[1],momentumIn[2]};//Need this since geometry works with const doubles not floats 
	double localMomentum [3];
	geo::gGeometry().master2LocalVec(getLocation(), momentum, localMomentum );
	//In the LOCAL coordinates this is just dx/dz and dy/dz in the LOCAL frame
	streamlog_out(DEBUG5) << "The local momentum (x,y,z) is: "<< localMomentum[0]<<","<< localMomentum[1] <<"," <<localMomentum[2] << std::endl;
	setIntersectionLocalXZ(localMomentum[0]/localMomentum[2]);
	setIntersectionLocalYZ(localMomentum[1]/localMomentum[2]);
	streamlog_out(DEBUG5) << "The XZ tilt is: "<< getIntersectionLocalXZ()<<" The YZ tilt is: "<<  getIntersectionLocalYZ()<< std::endl;
	streamlog_out(DEBUG5) << "EUTelState::setLocalXZAndYZIntersectionAndCurvatureUsingGlobalMomentum--------------------------END" << std::endl;
}

void EUTelState::setCombinedHitAndStateCovMatrixInLocalFrame(double cov[4]){
	_covCombinedMatrix[0] = cov[0];
	_covCombinedMatrix[1] = cov[1];
	_covCombinedMatrix[2] = cov[2];
	_covCombinedMatrix[3] = cov[3];
}

void EUTelState::setStateVec(TVectorD stateVec){
	float referencePoint[] = {stateVec[3],stateVec[4],0};
	setPositionLocal(referencePoint);
	setIntersectionLocalXZ(stateVec[1]);
	setIntersectionLocalYZ(stateVec[2]);
	setOmega(stateVec[0]);

}
///initialise
void EUTelState::initialiseCurvature(){
	if(getBeamCharge() == 0){
		throw(lcio::Exception( Utility::outputColourString("The beam charge is 0.Can not set curvature","RED"))); 
	}
	if(getBeamEnergy() == 0){
		throw(lcio::Exception( Utility::outputColourString("The beam Energy is zero. Can not set curvature","RED"))); 
	}
	
	setOmega(getBeamCharge()/getBeamEnergy());
}
//find
int EUTelState::findIntersectionWithCertainID(int nextSensorID, float intersectionPoint[], TVector3& momentumAtIntersection, float& arcLength ){
	streamlog_out(DEBUG5) << "-EUTelState::findIntersectionWithCertainID---------------------------BEGIN" << std::endl;
	TVector3 pVec = computeCartesianMomentum();
	streamlog_out(DEBUG5) << "Momentum (Global): " << pVec[0]<<","<<pVec[1]<<","<<pVec[2]<<","<<" Position (local): "<<getPosition()[0]<<","<<getPosition()[1]<<","<<getPosition()[2]<< std::endl;
	if(pVec.Mag() == 0){
		throw(lcio::Exception( Utility::outputColourString("The momentum is 0","RED"))); 
	}
	double posLocal[] =  {getPosition()[0],getPosition()[1],getPosition()[2] };
	double temp[] = {0.,0.,0.};
	geo::gGeometry().local2Master(getLocation() , posLocal, temp);//IMPORTANT:For strip sensors this will make the hit strip look like a pixel at (Xstriplocal,somevalue,somevalue).
	float posGlobal[] = { static_cast<float>(temp[0]), static_cast<float>(temp[1]), static_cast<float>(temp[2]) };
	int sensorID = geo::gGeometry().findIntersectionWithCertainID(posGlobal[0] ,  posGlobal[1] , posGlobal[2], pVec[0],pVec[1],pVec[2], getBeamCharge(), nextSensorID, intersectionPoint, momentumAtIntersection, arcLength ); 
	streamlog_out(DEBUG5) << "-EUTelState::findIntersectionWithCertainID--------------------------END" << std::endl;
	return sensorID;
}

//compute
TVector3 EUTelState::computeCartesianMomentum() const {
	streamlog_out(DEBUG2) << "EUTelState::computeCartesianMomentum()-------------------------BEGIN" << std::endl;
float tx = getIntersectionLocalXZ();float ty= getIntersectionLocalYZ(); float curvature = getOmega(); 
	streamlog_out(DEBUG2) << "Input parameters: tx,ty, beamq,invp "<<tx <<","<<ty<<","<<getBeamCharge()<<","<<curvature<<std::endl;
	if(getBeamCharge() == 0){
		throw(lcio::Exception( Utility::outputColourString("The beam charge is 0.","RED"))); 
	}
	if(curvature == INFINITY or curvature == NAN ){
		throw(lcio::Exception( Utility::outputColourString("The curvature is zero","RED"))); 
	}
	const double p  =  1. / (curvature *getBeamCharge() );//Must times but beam charge or we would be going in the wrong direction     
  const double pz = p/(sqrt(pow(tx,2)+pow(ty,2)+1));
  const double px = tx*pz;
  const double py = ty*pz;
	const double input[3]={px,py,pz};
	double momentum[3];
	geo::gGeometry().local2MasterVec(getLocation(),input,momentum);
	streamlog_out(DEBUG2) << "output global momentum: px,py, pz "<<momentum[0] <<","<<momentum[1]<<","<<momentum[2]<<","<<std::endl;
        
  streamlog_out(DEBUG2) << "-------------------------------EUTelState::computeCartesianMomentum()-------------------------END" << std::endl;
        
  return TVector3(momentum[0],momentum[1],momentum[2]);
}
TMatrix EUTelState::computePropagationJacobianFromLocalStateToNextLocalState(TVector3 positionEnd, TVector3 momentumEnd, float arcLength,float nextPlaneID) {
	streamlog_out(DEBUG2) << "-------------------------------EUTelState::computePropagationJacobianFromStateToThisZLocation()-------------------------BEGIN" << std::endl;
	if(arcLength == 0 or arcLength < 0 ){ 
		throw(lcio::Exception( Utility::outputColourString("The arc length is less than or equal to zero. ","RED"))); 
	}
	TMatrixD curvilinearJacobian = geo::gGeometry().getPropagationJacobianCurvilinear(arcLength,getOmega(), computeCartesianMomentum().Unit(),momentumEnd.Unit());
	streamlog_out(DEBUG0)<<"This is the curvilinear jacobian at sensor:" << std::endl; 
	streamlog_message( DEBUG0, curvilinearJacobian.Print();, std::endl; );
	streamlog_out(DEBUG0)<<"The state vector that create the curvilinear system is:" << std::endl; 
	print();
	TMatrixD localToCurvilinearJacobianStart =  geo::gGeometry().getLocalToCurvilinearTransformMatrix(computeCartesianMomentum(),getLocation() ,getBeamCharge() );
	streamlog_out(DEBUG0)<<"This is the local to curvilinear jacobian at sensor : " << std::endl; 
	streamlog_message( DEBUG0, localToCurvilinearJacobianStart.Print();, std::endl; );
	TMatrixD localToCurvilinearJacobianEnd =  geo::gGeometry().getLocalToCurvilinearTransformMatrix(momentumEnd,nextPlaneID ,getBeamCharge() );
	streamlog_out(DEBUG0)<<"This is the local to curvilinear jacobian at sensor at last next sensor : " << std::endl; 
	streamlog_message( DEBUG0, localToCurvilinearJacobianEnd.Print();, std::endl; );
	TMatrixD curvilinearToLocalJacobianEnd = localToCurvilinearJacobianEnd.Invert();
	streamlog_out(DEBUG0)<<"This is the curvilinear to local jacobian at sensor : " << std::endl; 
	streamlog_message( DEBUG0, curvilinearToLocalJacobianEnd.Print();, std::endl; );
	TMatrixD localToNextLocalJacobian = curvilinearToLocalJacobianEnd*curvilinearJacobian*localToCurvilinearJacobianStart;
	streamlog_out(DEBUG0)<<"This is the full jacobian : "<<  std::endl; 
	streamlog_message( DEBUG0, localToNextLocalJacobian.Print();, std::endl; );

	streamlog_out(DEBUG2) << "-------------------------------EUTelState::computePropagationJacobianFromStateToThisZLocation()-------------------------END" << std::endl;
	return localToNextLocalJacobian;
}
//print
void EUTelState::print(){
	streamlog_out(DEBUG2) << "The state vector//////////////////////////////////////////////////////" << endl;
	TVectorD stateVec = getStateVec();
	streamlog_message( DEBUG0, stateVec.Print();, std::endl; );
	TVectorD kinks = getKinks();
	streamlog_out(DEBUG1) <<"The kink of this state is: "<< kinks[0] <<" ,  " << kinks[1]<<std::endl;
	streamlog_out(DEBUG2) << "/////////////////////////////////////////////////////" << endl;
	streamlog_out(DEBUG1) <<"State memory location "<< this << " The sensor location of the state " <<getLocation()<<std::endl;
	if(getIsThereAHit()){
			streamlog_out(DEBUG1) <<"The hit ID of the state is "<<getTrackerHits().at(0)->id()<<std::endl;
	}else{
		streamlog_out(DEBUG1) <<"This state has no hit " <<endl;
	}
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



