#include "EUTelTrack.h"
using namespace eutelescope;
EUTelTrack::EUTelTrack(){} 

EUTelTrack::EUTelTrack(const EUTelTrack& state){
streamlog_out( DEBUG0 ) << "Call copy constructor for Track" << std::endl;

setLocation(state.getLocation());//This stores the location as a float in Z0 since no location for track LCIO. This is preferable to problems with storing hits.  
setPosition(const_cast<float*>(state.getReferencePoint()));//This will automatically take cartesian coordinate system and save it to reference point. //This is different from most LCIO applications since each track will have own reference point. 	
setDirectionXY(	state.getPhi());  // seed with 0 at first hit. This one is stored as tan (x/y) rather than angle as in LCIO format TO DO:This really should be give as a beam direction.  
setDirectionYZ(state.getTanLambda());  // seed with 0 at first hit.
setBeamCharge(state.getBeamCharge());//this is set for each state. to do: is there a more efficient way of doing this since we only need this stored once?
setBeamEnergy(state.getBeamEnergy());//this is saved in gev. 
initialiseCurvature(); //this will perform the calculation _beamq/_beame ad place in invp
addHit(const_cast<EVENT::TrackerHit*>(state.getTrackerHits()[0]));
setDimensionSize(state.getDimensionSize());
}
//getters
int EUTelTrack::getDimensionSize() const {
	int dimension = static_cast<int>(getD0());
	return dimension;
}
int EUTelTrack::getLocation() const {
	int location = static_cast<int>(getZ0());
	return location;
}
TVectorD EUTelTrack::getTrackStateVec() const { 
	streamlog_out( DEBUG1 ) << "EUTelTrack::getTrackStateVec()------------------------BEGIN" << std::endl;
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

	streamlog_out( DEBUG1 ) << "EUTelTrack::getTrackStateVec()------------------------END" << std::endl;
 	return stateVec;
}
TMatrixDSym EUTelTrack::getTrackStateCov() const {

	streamlog_out( DEBUG1 ) << "EUTelTrack::getTrackStateCov()----------------------------BEGIN" << std::endl;
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
	streamlog_out( DEBUG1 ) << "EUTelTrack::getTrackStateCov()----------------------------END" << std::endl;
}
int EUTelTrack::getNumberOfHitsOnTrack() const {
	int numberOfHitsOnTrack =0;
	const EVENT::TrackVec& states = getTracks();
	if(states.size() == 0){
		throw(lcio::Exception(Utility::outputColourString("The number of states is 0.", "RED"))); 	
	}
	streamlog_out(DEBUG0) <<"The number of states " << states.size()<<std::endl; 
	for(int i =0; i< states.size();++i){
		streamlog_out(DEBUG0) <<"The states memory address for loop number "<<i<<" " << &states<<std::endl; 
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
	return numberOfHitsOnTrack;
}
void EUTelTrack::getCombinedHitAndStateCovMatrixInLocalFrame( double (&cov)[4] ) const {
	cov[0] = _covCombinedMatrix[0];
	cov[1] = _covCombinedMatrix[1];
	cov[2] = _covCombinedMatrix[2];
	cov[3] = _covCombinedMatrix[3];
}
//TO DO:This matrix will only work for no tilted sensors. Must determine the generic projection matrix
TMatrixD EUTelTrack::getProjectionMatrix() const {
	TMatrixD projection(5,5);
	projection.Zero();
	TMatrixD proM2l(2, 2);
	proM2l.UnitMatrix();
	projection.SetSub(3, 3, proM2l);
	return projection;
}
//setters
void EUTelTrack::setDimensionSize(int dimension){
	setD0(static_cast<float>(dimension));
}
void EUTelTrack::setLocation(int location){
	float locationFloat = static_cast<float>(location);
	setZ0(locationFloat);
}
void EUTelTrack::setBeamEnergy(float beamE){
	_beamE = beamE;
}

void EUTelTrack::setBeamCharge(float beamQ){
	_beamQ = beamQ;
}
void EUTelTrack::setDirectionYZ(float directionYZ){
	setTanLambda(directionYZ);
}	
void EUTelTrack::setDirectionXY(float directionXY){
	setPhi(directionXY);//This is only a container. So hold inverseTan(Phi)
}
void EUTelTrack::setPosition(float position[]){
	setReferencePoint(position);
}
void EUTelTrack::setCombinedHitAndStateCovMatrixInLocalFrame(double cov[4]){
	_covCombinedMatrix[0] = cov[0];
	_covCombinedMatrix[1] = cov[1];
	_covCombinedMatrix[2] = cov[2];
	_covCombinedMatrix[3] = cov[3];
}

void EUTelTrack::setTrackStateVecPlusZParameter(TVectorD stateVec,float zParameter){
	float referencePoint[] = {stateVec[3],stateVec[4],zParameter};
	setReferencePoint(referencePoint);
	setPhi(stateVec[1]);
	setTanLambda(stateVec[2]);
	setOmega(stateVec[0]);

}
///initialise
void EUTelTrack::initialiseCurvature(){
	if(_beamQ == 0){
		throw(lcio::Exception( Utility::outputColourString("The beam charge is 0.Can not set curvature","RED"))); 
	}
	if(_beamE == 0){
		throw(lcio::Exception( Utility::outputColourString("The beam Energy is zero. Can not set curvature","RED"))); 
	}
	
	setOmega(_beamQ/_beamE);
}
//find
int EUTelTrack::findIntersectionWithCertainID(int nextSensorID, float intersectionPoint[] ){
	streamlog_out(DEBUG5) << "-EUTelTrack::findIntersectionWithCertainID---------------------------BEGIN" << std::endl;
	TVector3 pVec = computeCartesianMomentum();
	streamlog_out(DEBUG5) << "Momentum: " << pVec[0]<<","<<pVec[1]<<","<<pVec[2]<<","<< std::endl;
	if(pVec.Mag() == 0){
		throw(lcio::Exception( Utility::outputColourString("The momentum is 0","RED"))); 
	}
	int sensorID = geo::gGeometry().findIntersectionWithCertainID( getReferencePoint()[0],  getReferencePoint()[1] , getReferencePoint()[2], pVec[0],pVec[1],pVec[2], _beamQ, nextSensorID, intersectionPoint ); 
	streamlog_out(DEBUG5) << "-EUTelTrack::findIntersectionWithCertainID--------------------------END" << std::endl;
	return sensorID;
}

//compute
TVector3 EUTelTrack::computeCartesianMomentum(){
	streamlog_out(DEBUG2) << "EUTelTrack::computeCartesianMomentum()-------------------------BEGIN" << std::endl;
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
        
  streamlog_out(DEBUG2) << "-------------------------------EUTelTrackStateImpl::getPfromCartesianParameters()-------------------------END" << std::endl;
        
  return TVector3(px,py,pz);
}
TMatrix EUTelTrack::computePropagationJacobianFromStateToThisZLocation(float zPosition){
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
void EUTelTrack::print(){
	TVectorD stateVec = getTrackStateVec();
	stateVec.Print();
}
//Overload operators.
bool EUTelTrack::operator<(const EUTelTrack compareState ) const {
	return getReferencePoint()[2]<compareState.getReferencePoint()[2];
}
