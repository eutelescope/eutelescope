#include "EUTelTrack.h"
using namespace eutelescope;
EUTelTrack::EUTelTrack(){} 

//getters
int EUTelTrack::getLocation(){
	int location = static_cast<int>(getZ0());
	return location;
}

//setters
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
	if(curvature == 0){
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


