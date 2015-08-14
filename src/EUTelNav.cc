#include "EUTelNav.h"	

namespace eutelescope 
{

TMatrixD EUTelNav::getMeasToGlobal(TVector3 t1w, TMatrixD TRotMatrix )
{
	TMatrixD transM2l(5,5);
	transM2l.UnitMatrix();
	std::vector<double> slope;
	slope.push_back(t1w[0]/t1w[2]); slope.push_back(t1w[1]/t1w[2]);
	double norm = std::sqrt(pow(slope.at(0),2) + pow(slope.at(1),2) + 1);//This works since we have in the curvinlinear frame (dx/dz)^2 +(dy/dz)^2 +1 so time through by dz^2
	TVector3 direction;
	direction[0] = (slope.at(0)/norm); direction[1] =(slope.at(1)/norm);	direction[2] = (1.0/norm);
	TMatrixD xyDir(2, 3);
	xyDir[0][0] = 1; xyDir[0][1]=0.0; xyDir[0][2]=-slope.at(0);  
	xyDir[1][0] = 0; xyDir[1][1]=1.0; xyDir[1][2]=-slope.at(1);  
	TVector3 normalVec;
	normalVec[0] = TRotMatrix[0][2];	normalVec[1] = TRotMatrix[1][2];	normalVec[2] = TRotMatrix[2][2];
	double cosInc = direction*normalVec;
//	std::cout<<"Here is cosInc " << cosInc <<std::endl;
	TMatrixD measDir(3,2);
	measDir[0][0] = TRotMatrix[0][0];	measDir[0][1] = TRotMatrix[0][1];
	measDir[1][0] = TRotMatrix[1][0];	measDir[1][1] = TRotMatrix[1][1];
	measDir[2][0] = TRotMatrix[2][0];	measDir[2][1] = TRotMatrix[2][1];
    streamlog_out( DEBUG0 ) << "CALCULATE LOCAL TO GLOBAL STATE TRANSFORMATION... " << std::endl;
    streamlog_out( DEBUG0 ) << "Inputs... " << std::endl;

    streamlog_out( DEBUG0 ) << "The (X,Y)-axis of the global frame relative to the local  " << std::endl;
    streamlog_message( DEBUG0, measDir.Print();, std::endl; );
    streamlog_out( DEBUG0 ) << "The propagator (unit axis) (Dx,Dy)   " << std::endl;
    streamlog_message( DEBUG0, xyDir.Print();, std::endl; );
    double scaleFactor = cosInc/direction[2];
    streamlog_out( DEBUG0 ) << "Scale factor (s) " << scaleFactor << std::endl;
	TMatrixD proM2l(2,2);
	proM2l = xyDir*measDir; 
    streamlog_out( DEBUG0 ) << "Propagators... " << std::endl;

    streamlog_out( DEBUG0 ) << "PROJECTION MATRIX (shifts) (Dx,Dy)x(X,Y)" << std::endl;
    streamlog_message( DEBUG0, proM2l.Print();, std::endl; );
    streamlog_out( DEBUG0 ) << "PROJECTION MATRIX (incidence) s(Dx,Dy)x(X,Y) " << std::endl;
    TMatrixD proM2lInc = scaleFactor*proM2l;
    streamlog_message( DEBUG0, proM2lInc.Print();, std::endl; );

	transM2l.SetSub(1,1,proM2lInc);
	transM2l.SetSub(3,3,proM2l);
    streamlog_out( DEBUG0 ) << "OUTPUT:(Local to Global): " << std::endl;
    streamlog_message( DEBUG0, transM2l.Print();, std::endl; );

	return transM2l;
}
///This function creates a jacobain which links one state to another in the EUTelGlobal frame. 
/**
 * \param [in] ds Arc length between two states. 
 * \param [in] t1w Momentum on the initial states 
 * \return Jacobain 5x5 which links the two states 
 */

TMatrixD EUTelNav::getPropagationJacobianGlobalToGlobal(float ds, TVector3 t1w)
{
	t1w.Unit();
	std::vector<double> slope;
	slope.push_back(t1w[0]/t1w[2]); slope.push_back(t1w[1]/t1w[2]);
	double norm = std::sqrt(pow(slope.at(0),2) + pow(slope.at(1),2) + 1);//not this works since we have in the curvinlinear frame (dx/dz)^2 +(dy/dz)^2 +1 so time through by dz^2
	TVector3 direction;
	direction[0] = (slope.at(0)/norm); direction[1] =(slope.at(1)/norm);	direction[2] = (1.0/norm);
//	std::cout <<"DIRECTION: "<< direction[0] <<"   " << direction[1]<< "  "<< direction[2] <<std::endl;
	double sinLambda = direction[2]; 
	const gear::BField& Bfield = geo::gGeometry().getMagneticField();
	gear::Vector3D vectorGlobal(0.1,0.1,0.1);


	const double Bx = (Bfield.at( vectorGlobal ).x());  
	const double By = (Bfield.at( vectorGlobal ).y());
	const double Bz = (Bfield.at( vectorGlobal ).z());

	TVector3 b(Bx, By, Bz);
	TVector3 BxT = b.Cross(direction);
//	std::cout << "BxT" << BxT[0] << "  ,  " <<  BxT[1] <<"   ,  " <<BxT[2] << std::endl;
	TMatrixD xyDir(2, 3);
//    streamlog_out( DEBUG0 ) << "CALCULATE GLOBAL TO GLOBAL STATE TRANSFORMATION... " << std::endl;

	xyDir[0][0] = 1.0; xyDir[0][1]=0.0; xyDir[0][2]=-slope.at(0);  
	xyDir[1][0] = 0; xyDir[1][1]=1.0; xyDir[1][2]=-slope.at(1);  
    streamlog_out( DEBUG0 ) << "The propagator (Dx,Dy)   " << std::endl;
    streamlog_message( DEBUG0, xyDir.Print();, std::endl; );


	TMatrixD bFac(2,1);
	TMatrixD BxTMatrix(3,1);
	BxTMatrix.Zero();
	BxTMatrix[0][0] =BxT[0];	BxTMatrix[1][0] =BxT[1];	BxTMatrix[2][0] =BxT[2]; 
 
//	bFac = -0.0002998 * (xyDir*BxTMatrix); 
	TVector3 bFacVec = getBFac();
//	std::cout <<std::scientific<< "bFac" << bFac[0][0] << "  ,  " <<  bFac[1][0] << std::endl;
	TMatrixD ajac(5, 5);
	ajac.UnitMatrix();
//	if(b.Mag() < 0.001 ){
//			ajac[3][2] = ds * std::sqrt(t1w[0] * t1w[0] + t1w[2] * t1w[2]);
//			ajac[4][1] = ds;
//	}else{
		ajac[1][0] = bFacVec[0]*ds/sinLambda;
//		std::cout<<"Jacobian entries;" << ajac[1][0] << std::endl;
//		std::cout<<"1,0 " << ajac[1][0] << std::endl;
		ajac[2][0] = bFacVec[1]*ds/sinLambda;
//		std::cout<< std::scientific <<"2,0 " << ajac[2][0] << std::endl;
		ajac[3][0] = 0.5*bFacVec[0]*ds*ds;
//		std::cout<<"3,0 " << ajac[3][0] << std::endl;
		ajac[4][0] = 0.5*bFacVec[1]*ds*ds;
//		std::cout<< std::scientific  <<"4,0 " << ajac[4][0] << std::endl;
		ajac[3][1] = ds*sinLambda; 
//		std::cout<<"3,1 " << ajac[3][1] << std::endl;
		ajac[4][2] = ds*sinLambda; 
//		std::cout<<"4,2 " << ajac[4][2] << std::endl;

//	}
    streamlog_out( DEBUG0 ) << "Global to Global jacobian: " << std::endl;
    streamlog_message( DEBUG0, ajac.Print();, std::endl; );
	return ajac;
}


void EUTelNav::getTrackAvePara(EUTelHit & hit1, EUTelHit & hit2, std::vector<double>& offset, std::vector<double>& trackSlope){
    //NOW CREATE TRACK CANDIDATE
    offset.push_back(hit1.getPositionGlobal()[0]); 
    offset.push_back(hit1.getPositionGlobal()[1]); 
    offset.push_back(hit1.getPositionGlobal()[2]); 
    offset.push_back(hit2.getPositionGlobal()[2]); 
    const double dz = offset.at(3) - offset.at(2);
    trackSlope.push_back((hit2.getPositionGlobal()[0] - offset.at(0))/dz);trackSlope.push_back((hit2.getPositionGlobal()[1] - offset.at(1))/dz);
	streamlog_out(DEBUG1) << "Track average parameters: " << std::endl;
	streamlog_out(DEBUG1) << "Offsets:  " << offset.at(0)<<" " << offset.at(1)<<" " << offset.at(2) <<" " <<offset.at(3) << std::endl;
	streamlog_out(DEBUG1) << "Slope  " << trackSlope.at(0)<<" " << trackSlope.at(1) << std::endl;

}
double EUTelNav::getCorr(EUTelHit & hit1, EUTelHit & hit2, EUTelHit & hit3, EUTelHit & hit4){
    std::vector<double> curvCorr;
    std::vector<double> slopesArmOne; 
    slopesArmOne.push_back((hit2.getPositionGlobal()[0]-hit1.getPositionGlobal()[0])/(hit2.getPositionGlobal()[2]-hit1.getPositionGlobal()[2]));
    slopesArmOne.push_back((hit2.getPositionGlobal()[1]-hit1.getPositionGlobal()[1])/(hit2.getPositionGlobal()[2]-hit1.getPositionGlobal()[2]));
    std::vector<double> slopesArmTwo; 
    slopesArmTwo.push_back((hit4.getPositionGlobal()[0]-hit3.getPositionGlobal()[0])/(hit4.getPositionGlobal()[2]-hit3.getPositionGlobal()[2]));
    slopesArmTwo.push_back((hit4.getPositionGlobal()[1]-hit3.getPositionGlobal()[1])/(hit4.getPositionGlobal()[2]-hit3.getPositionGlobal()[2]));

    double averageDistArmOne = (hit2.getPositionGlobal()[2] +  hit1.getPositionGlobal()[2])/2.0;
    double averageDistArmTwo = (hit4.getPositionGlobal()[2] +  hit3.getPositionGlobal()[2])/2.0;

    //Slope change with curvature constant.
    double dSlopeXDCurv = _bFac[0]*(averageDistArmOne-averageDistArmTwo);  
    double dSlopeYDCurv = _bFac[1]*(averageDistArmOne-averageDistArmTwo);  
    //correct curvature
    double corr = (dSlopeXDCurv*(slopesArmOne.at(0)- slopesArmTwo.at(0)) + dSlopeYDCurv*(slopesArmOne.at(1)- slopesArmTwo.at(1)))/(pow(dSlopeXDCurv,2)+pow(dSlopeYDCurv,2));
    streamlog_out(DEBUG0) << "Corrected q/p: " << corr  <<std::endl; 
    return corr;
}
void EUTelNav::getTrackPredictionFromParam(std::vector<double> const & offset, std::vector<double> const & trackSlope,double const & qOverP, double const & posZ, Eigen::Vector3d & posPred, std::vector<double>& slopePred){
    std::vector<double> curvCorr; curvCorr.push_back(qOverP*_bFac[0]);curvCorr.push_back(qOverP*_bFac[1]);
    double dz1 = posZ - offset.at(2);
    double dz2 = posZ - offset.at(3); 
    double dz = (dz1 + dz2)/2.0;
    double posX = offset.at(0) + dz1*trackSlope.at(0) + 0.5*dz1*dz2*curvCorr.at(0);
    double posY = offset.at(1) + dz1*trackSlope.at(1) + 0.5*dz1*dz2*curvCorr.at(1);
    posPred << posX, posY, posZ ;
    slopePred.push_back(trackSlope.at(0)+dz*curvCorr.at(0));
    slopePred.push_back(trackSlope.at(1)+dz*curvCorr.at(1));
}

///Curvature in x/y is passed to naviagation.
std::vector<double>  EUTelNav::getCurvXY(){
    //Defined the same as saved in track parameters.
    const float omega = -1.0/_intBeamE;
    TVector3 bFac = getBFac();
    streamlog_out(DEBUG0) << "BFac field unit: " << bFac[0] << "  " << bFac[1] <<"  "<< bFac[2] << "  Omega: " << omega << std::endl;
    //Note the cross product
    const double curvX = bFac[0]*omega; 
    const double curvY = bFac[1]*omega; 
    std::vector<double> curv;
    curv.push_back(curvX); curv.push_back(curvY);
    //streamlog_out(DEBUG0) << "The curvature calculated: " << curv.at(0) << "  " << curv.at(1) << std::endl;

    return curv;

}
///cZxB is passed to navigation.
TVector3  EUTelNav::getBFac(){
    const gear::BField& Bfield = geo::gGeometry().getMagneticField();
    gear::Vector3D vectorGlobal(0.1,0.1,0.1);
    const double Bx = (Bfield.at( vectorGlobal ).x());  
    const double By = (Bfield.at( vectorGlobal ).y());
    const double Bz = (Bfield.at( vectorGlobal ).z());
    TVector3 B(Bx, By, Bz );
    TVector3 H = (B.Unit());
    TVector3 bFac = 0.0003*(TVector3(0,0,1).Cross(H));
    return bFac;
}
void EUTelNav::init(double beamEnergy){
    _intBeamE = beamEnergy;
    _bFac = getBFac();
    _curv = getCurvXY();
}
        TVector3 EUTelNav::_bFac;
        std::vector<double> EUTelNav::_curv;
        double EUTelNav::_intBeamE;

} //namespace eutelescope
