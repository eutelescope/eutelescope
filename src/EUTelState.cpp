#include "EUTelState.h"

using namespace eutelescope;
EUTelState::EUTelState():
_beamQ(-1.0),
_beamE(5.0){
} 

EUTelState::EUTelState(const EUTelState& state){
	setDimensionSize(state.getDimensionSize());
	setLocation(state.getLocation());//This stores the location as a float in Z0 since no location for track LCIO. This is preferable to problems with storing hits.  
	setPosition(state.getPosition());//This will automatically take cartesian coordinate system and save it to reference point. //This is different from most LCIO applications since each track will have own reference point. 	
	setDirectionXY(state.getDirectionXY());    
	setDirectionYZ(state.getDirectionYZ());  
	setBeamCharge(state.getBeamCharge());//this is set for each state. to do: is there a more efficient way of doing this since we only need this stored once?
	setBeamEnergy(state.getBeamEnergy());//this is saved in gev. 
	initialiseCurvature(); //this will perform the calculation _beamq/_beame ad place in invp
	if(!state.getTrackerHits().empty()){
		addHit(const_cast<EVENT::TrackerHit*>(state.getTrackerHits().at(0)));
	}
}
//getters
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
TVectorD EUTelState::getTrackStateVec() const { 
	streamlog_out( DEBUG1 ) << "EUTelState::getTrackStateVec()------------------------BEGIN" << std::endl;
	TVectorD stateVec(5);
	stateVec[3] = getReferencePoint()[0];
	stateVec[4] = getReferencePoint()[1];
	stateVec[1] = getPhi();
	stateVec[2] = getTanLambda();
	stateVec[0] = getOmega();
			
	if ( streamlog_level(DEBUG0) ){
		streamlog_out( DEBUG0 ) << "Track state:" << std::endl;
		stateVec.Print();
	}
	if(stateVec[0] == INFINITY or stateVec[0] == NAN ){
		throw(lcio::Exception( Utility::outputColourString("Passing a state vector where curvature is not defined","RED"))); 
	}

	streamlog_out( DEBUG1 ) << "EUTelState::getTrackStateVec()------------------------END" << std::endl;
 	return stateVec;
}
TMatrixDSym EUTelState::getTrackStateCov() const {

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
	return projection;
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
	_beamE = beamE;
}

void EUTelState::setBeamCharge(float beamQ){
	_beamQ = beamQ;
}
void EUTelState::setDirectionYZ(float directionYZ){
	setTanLambda(directionYZ);
}	
void EUTelState::setDirectionXY(float directionXY){
	setPhi(directionXY);//This is only a container. So hold inverseTan(Phi)
}
void EUTelState::setPosition(float position[]){
	setReferencePoint(position);
}
void EUTelState::setCombinedHitAndStateCovMatrixInLocalFrame(double cov[4]){
	_covCombinedMatrix[0] = cov[0];
	_covCombinedMatrix[1] = cov[1];
	_covCombinedMatrix[2] = cov[2];
	_covCombinedMatrix[3] = cov[3];
}

void EUTelState::setTrackStateVecPlusZParameter(TVectorD stateVec,float zParameter){
	float referencePoint[] = {stateVec[3],stateVec[4],zParameter};
	setReferencePoint(referencePoint);
	setPhi(stateVec[1]);
	setTanLambda(stateVec[2]);
	setOmega(stateVec[0]);

}
///initialise
void EUTelState::initialiseCurvature(){
	if(_beamQ == 0){
		throw(lcio::Exception( Utility::outputColourString("The beam charge is 0.Can not set curvature","RED"))); 
	}
	if(_beamE == 0){
		throw(lcio::Exception( Utility::outputColourString("The beam Energy is zero. Can not set curvature","RED"))); 
	}
	
	setOmega(_beamQ/_beamE);
}
//find
int EUTelState::findIntersectionWithCertainID(int nextSensorID, float intersectionPoint[] ){
	streamlog_out(DEBUG5) << "-EUTelState::findIntersectionWithCertainID---------------------------BEGIN" << std::endl;
	TVector3 pVec = computeCartesianMomentum();
	streamlog_out(DEBUG5) << "Momentum: " << pVec[0]<<","<<pVec[1]<<","<<pVec[2]<<","<< std::endl;
	if(pVec.Mag() == 0){
		throw(lcio::Exception( Utility::outputColourString("The momentum is 0","RED"))); 
	}
	int sensorID = geo::gGeometry().findIntersectionWithCertainID( getReferencePoint()[0],  getReferencePoint()[1] , getReferencePoint()[2], pVec[0],pVec[1],pVec[2], _beamQ, nextSensorID, intersectionPoint ); 
	streamlog_out(DEBUG5) << "-EUTelState::findIntersectionWithCertainID--------------------------END" << std::endl;
	return sensorID;
}

//compute
TVector3 EUTelState::computeCartesianMomentum(){
	streamlog_out(DEBUG2) << "EUTelState::computeCartesianMomentum()-------------------------BEGIN" << std::endl;
float tx = getPhi();float ty= getTanLambda(); float curvature = getOmega(); 
	streamlog_out(DEBUG2) << "Input parameters: tx,ty, beamq,invp "<<tx <<","<<ty<<","<<_beamQ<<","<<curvature<<std::endl;
	if(_beamQ == 0){
		throw(lcio::Exception( Utility::outputColourString("The beam charge is 0.","RED"))); 
	}
	if(curvature == INFINITY or curvature == NAN ){
		throw(lcio::Exception( Utility::outputColourString("The curvature is zero","RED"))); 
	}
	const double p  =  1. / (curvature *_beamQ );     
  const double px = p*tx / sqrt( 1. + tx*tx + ty*ty );
  const double py = p*ty / sqrt( 1. + tx*tx + ty*ty );
  const double pz = p    / sqrt( 1. + tx*tx + ty*ty );

	streamlog_out(DEBUG2) << "Output parameters: px,py, pz "<<px <<","<<py<<","<<pz<<","<<std::endl;
        
  streamlog_out(DEBUG2) << "-------------------------------EUTelStateStateImpl::getPfromCartesianParameters()-------------------------END" << std::endl;
        
  return TVector3(px,py,pz);
}
TMatrix EUTelState::computePropagationJacobianFromStateToThisZLocation(float zPosition){
	float dz = zPosition - getReferencePoint()[3];
	TVector3 pVec = computeCartesianMomentum();
	TMatrix jacobian(5,5);
	jacobian.Zero();
	float x = getReferencePoint()[0]; float y = getReferencePoint()[1];	float z = getReferencePoint()[2];
	streamlog_out( DEBUG1 ) << "These are the parameters being used in the calculation of the Jacobian. "<< "X= "<< x << " Y= "<< y << "Z= "<< z << " Px= " << pVec[0] << " Py= " << pVec[1] << " Pz= " << pVec[2] << " Qbeam= " <<_beamQ  << std::endl;
	jacobian = geo::gGeometry().getPropagationJacobianF(  x, y, z, pVec[0],pVec[1],pVec[2], _beamQ, dz);

	return jacobian;

}
//print
void EUTelState::print(){
	streamlog_out(DEBUG2) << "The state vector//////////////////////////////////////////////////////" << endl;
	TVectorD stateVec = getTrackStateVec();
	stateVec.Print();
	streamlog_out(DEBUG2) << "/////////////////////////////////////////////////////" << endl;
	streamlog_out(DEBUG1) <<"State memory location "<< this << " The sensor location of the state " <<getLocation()<<std::endl;
	if(!getTrackerHits().empty()){
			streamlog_out(DEBUG1) <<"The hit ID of the state is "<<getTrackerHits().at(0)->id()<<std::endl;
	}
}	
//Overload operators.
bool EUTelState::operator<(const EUTelState compareState ) const {
	return getPosition()[2]<compareState.getPosition()[2];
}





