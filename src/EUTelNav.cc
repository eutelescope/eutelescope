#include "EUTelNav.h"	

namespace eutelescope 
{

/*This function given position/momentum of a particle. Will give you the approximate jacobian at any point along the track. 
 * This effectively relates changes in the particle position/momentum at the original to some distant point. 
 * So if I change the initial position by x amount how much will all the other variables position/momentum at the new position change? 
 * This is what the Jacobian tells you.
 *
 * Calculate track parameters propagation jacobian for given track state
 *  and propagation distance. The expressions were derived in parabolic approximation
 *  valid for small values of propagation distance |dz| < 10cm. Can be iterated if necessary.
 * 
 * @param ts track state
 * @param dz propagation distance
 * @return 
 */
TMatrix EUTelNav::getPropagationJacobianF( float x0, float y0, float z0, float px, float py, float pz, float beamQ, float dz )
{
		// The formulas below are derived from equations of motion of the particle in
		// magnetic field under assumption |dz| small. Must be valid for |dz| < 10 cm

		//Here is the conversion from k=(GeV/c) KG mm^-1
		const double k = 0.299792458*pow(10,-4);

		const double pxMomentum = px;//change momentum to GeV/c 
		const double pyMomentum = py;//change momentum to GeV/c 
		const double pzMomentum = pz;//change momentum to GeV/c 

		TVector3 pVec(pxMomentum, pyMomentum, pzMomentum);	

		//These are the track parameters
		const double invP = beamQ/pVec.Mag();//This should be in 1/(GeV/c)
		const double tx0 = px/pVec.Mag();//Unitless  
		const double ty0 = py/pVec.Mag();

		//And the magnetic field vector (assume uniform field)
		gear::Vector3D vectorGlobal( x0, y0, z0 );
		const gear::BField&   B =  geo::gGeometry().getMagneticField();
		//times but 10 to convert from Tesla to KiloGauss. 1 T = 10^4 Gauss.
		const double Bx         = (B.at( vectorGlobal ).x())*10;
		const double By         = (B.at( vectorGlobal ).y())*10;
		const double Bz         = (B.at( vectorGlobal ).z())*10;

		const double sqrtFactor = sqrt( 1. + tx0*tx0 + ty0*ty0 );

		const double Ax = sqrtFactor * (  ty0 * ( tx0 * Bx + Bz ) - ( 1. + tx0*tx0 ) * By );
		const double Ay = sqrtFactor * ( -tx0 * ( ty0 * By + Bz ) + ( 1. + ty0*ty0 ) * Bx );

		// Partial derivatives
		const double dAxdty0 = ty0 * Ax / (sqrtFactor*sqrtFactor) + sqrtFactor*( tx0*Bx + Bz );
		const double dAydtx0 = tx0 * Ay / (sqrtFactor*sqrtFactor) + sqrtFactor*( -ty0*By - Bz );

		const double dxdtx0 = dz;
		const double dxdty0 = 0.5 * invP * k * dz*dz * dAxdty0;

		const double dydtx0 = 0.5 * invP * k * dz*dz * dAydtx0;
		const double dydty0 = dz;

		const double dtxdty0 = invP * k * dz * dAxdty0;
		const double dtydtx0 = invP * k * dz * dAydtx0;

		const double dxdinvP0 = 0.5 * k * dz*dz * Ax;
		const double dydinvP0 = 0.5 * k * dz*dz * Ay;

		const double dtxdinvP0 = k * dz * Ax;
		const double dtydinvP0 = k * dz * Ay;

		//Fill-in matrix elements
		TMatrix jacobianF(5,5);
		jacobianF.UnitMatrix();

		jacobianF[3][1] = dxdtx0;	jacobianF[3][2] = dxdty0;	jacobianF[3][0] = dxdinvP0;
		jacobianF[4][1] = dydtx0;	jacobianF[4][2] = dydty0;	jacobianF[4][0] = dydinvP0;
		jacobianF[1][2] = dtxdty0;	jacobianF[1][0] = dtxdinvP0;
		jacobianF[2][1] = dtydtx0;	jacobianF[2][0] = dtydinvP0;

		if( streamlog_level(DEBUG0) )
		{
				streamlog_out( DEBUG0 ) << "Propagation jacobian: " << std::endl;
				jacobianF.Print();
		}
		return jacobianF;

}   

/* Here we define the transformation between the curvilinear and the local frame. Note the local frame we have is defined as the local frame of the telescope. 
 * We do this since the local frame is arbitrary. However, we must describe the local frame in the curvilinear frame, since we connect the two frames via the same global frame. 
 * The transformations are the same as described below. This is unlike the curvilinear frame which can not have the particles moving in z due to construction. 
 * Therefore we follow the paper: Derivation of Jacobians for the propagation of covariance matrices of track parameters in homogeneous magnetic fields; 
 * which implies that anytime we are describing the curvilinear system using our global momentums then we need to transform the from our system to theirs.
 * Within the function this will be clearly labelled.
 * This is a simple transform our x becomes their(curvilinear y), our y becomes their z and z becomes x
 * However, this is ok since we never directly access the curvilinear system. It is only a bridge between two local systems.
 */ 
TMatrixD EUTelNav::getLocalToCurvilinearTransformMatrix(TVector3 globalMomentum, int  planeID, float charge)
{
		const gear::BField& Bfield = geo::gGeometry().getMagneticField();

		//Since field is homogeneous this seems silly but we need to specify a position to geometry to get B-field.
		gear::Vector3D vectorGlobal(0.1,0.1,0.1);
		
		//Magnetic field must be changed to curvilinear coordinate system, as this is used in the curvilinear jacobian
		//We times bu 0.3 due to units of other variables. It is related to the speed of light over 1 ns.  See paper. Must be Tesla
		//TO DO: CHECK THAT I CAN REMOVE THIS 0.3 TERM FROM B FIELDS NOW WITH CLAUS'S JACOBIAN LIMIT
//		const double Bx = (Bfield.at( vectorGlobal ).z())*0.3;
//		const double By = (Bfield.at( vectorGlobal ).y())*0.3;
//		const double Bz = (Bfield.at( vectorGlobal ).x())*0.3;
		const double Bx = (Bfield.at( vectorGlobal ).z());
		const double By = (Bfield.at( vectorGlobal ).x());
		const double Bz = (Bfield.at( vectorGlobal ).y());
	
		TVector3 B(Bx, By, Bz);
		TVector3 H = B.Unit();

		//We transform the momentum to curvilinear frame, as this is what is used to describe the curvilinear frame
		TVector3 curvilinearGlobalMomentum(globalMomentum[2], globalMomentum[0], globalMomentum[1]);
		//With no magnetic field this will point in z direction.	
		TVector3 T = curvilinearGlobalMomentum.Unit();
		//Note here we create the curvilinear frame using the implicit global frame which with zGlobalNormal. Note T,U, V is the actual curvilnear frame.	
		//TO DO: I did not change T[1] here. Do in need to? Not sure. I assume not since it should be equivalent to T[2]
		const float cosLambda = sqrt(T[0]*T[0] + T[1]*T[1]);
		//We use the global curvilinear frame z direction to create the other unit axis
		TVector3 zGlobalNormal(0, 0, 1); 
		TVector3 U = (zGlobalNormal.Cross(T)).Unit(); 
		TVector3 V = (T.Cross(U));
 	
		streamlog_out(DEBUG0)<<"The Z(V) axis of the curvilinear system in the global system (Note not the same as global telescope frame)."<< std::endl; 
		streamlog_message( DEBUG0, V.Print();, std::endl; );
		streamlog_out(DEBUG0)<<"The Y(U) axis of the curvilinear system in the global system (Note not the same as global telescope frame)"<< std::endl; 
		streamlog_message( DEBUG0, U.Print();, std::endl; );
		streamlog_out(DEBUG0)<<"The X(T) axis of the curvilinear system in the global system  (Note not the same as global telescope frame). Should be beam direction."<< std::endl; 
		streamlog_message( DEBUG0, T.Print();, std::endl; );

		TVector3 ITelescopeFrame;
		TVector3 KTelescopeFrame;
		TVector3 JTelescopeFrame;

		///This is the EUTelescope local z direction.
		//314 is the number we chose to specify a scattering plane.
		if(planeID != 314)
		{ 
				ITelescopeFrame = geo::gGeometry().siPlaneNormal(planeID);       
				KTelescopeFrame = geo::gGeometry().siPlaneXAxis(planeID);	
				JTelescopeFrame = geo::gGeometry().siPlaneYAxis(planeID);
		}
		else
		{
				ITelescopeFrame.SetXYZ(0,0,1);
				KTelescopeFrame.SetXYZ(1,0,0);
				JTelescopeFrame.SetXYZ(0,1,0);
		}

		TVector3 I(ITelescopeFrame[2], ITelescopeFrame[1], ITelescopeFrame[0]);
		TVector3 K(KTelescopeFrame[2], KTelescopeFrame[1], KTelescopeFrame[0]);
		TVector3 J(JTelescopeFrame[2], JTelescopeFrame[1], JTelescopeFrame[0]);

		streamlog_out(DEBUG0)<<"The Z(J) axis of local system in the global curvilinear system"<< std::endl; 
		streamlog_message( DEBUG0, J.Print();, std::endl; );
		streamlog_out(DEBUG0)<<"The Y(K) axis of local system in the global curvilinear system"<< std::endl; 
		streamlog_message( DEBUG0, K.Print();, std::endl; );
		streamlog_out(DEBUG0)<<"The X(I) axis of local system in the global curvilinear system"<< std::endl; 
		streamlog_message( DEBUG0, I.Print();, std::endl; );

		TVector3 N = (H.Cross(T)).Unit();
		
		const double alpha = (H.Cross(T)).Mag();
		const double Q = -(B.Mag())*(charge/(curvilinearGlobalMomentum.Mag()));//You could use end momentum since it must be constant
		
		const double TDotI = T.Dot(I);
		const double TDotJ = T.Dot(J);
		const double TDotK = T.Dot(K);
		const double VDotJ = V.Dot(J);
		const double VDotK = V.Dot(K);
		const double VDotN = V.Dot(N);
		const double UDotJ = U.Dot(J);
		const double UDotK = U.Dot(K);
		const double UDotN = U.Dot(N);
	
		TMatrix jacobian(5,5);
		jacobian.Zero();
	
		/*	Matrix has following (X) entries set:
 		 *	X 0 0 0 0
		 *	0 X X X X
		 *	0 X X X X
		 *	0 0 0 X X
		 *	0 0 0 X X
		 */ 	

		//First Row
		jacobian[0][0]=1;
		//Second Row 
		jacobian[1][1]=TDotI*VDotJ;
		jacobian[1][2]=TDotI*VDotK;
		jacobian[1][3]=-alpha*Q*TDotJ*VDotN;
		jacobian[1][4]=-alpha*Q*TDotK*VDotN;
		//Third Row
		jacobian[2][1]=(TDotI*UDotJ)/cosLambda;
		jacobian[2][2]=(TDotI*UDotK)/cosLambda;
		jacobian[2][3]=(-alpha*Q*TDotJ*UDotN)/cosLambda;
		jacobian[2][4]=(-alpha*Q*TDotK*UDotN)/cosLambda;
		//Forth Row
		jacobian[3][3]=UDotJ;
		jacobian[3][4]=UDotK;
		//Fifth Row		
		jacobian[4][3]=VDotJ;
		jacobian[4][4]=VDotK;
		
		return jacobian;
}
///This will relate infintesimal changes in the local state vector to the global one. 
/**
 * \param [in] t1w Momentum of the state
 * \return Jacobain 5x5 which links the local and global states.
 */


TMatrixD EUTelNav::getMeasToGlobal(TVector3 t1w, int  planeID)
{
//	std::cout<<"Plane ID " << planeID <<std::endl;
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
	TMatrixD TRotMatrix(3,3);
	TRotMatrix	=  geo::gGeometry().getRotMatrix( planeID );
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
//    streamlog_out( DEBUG0 ) << "The propagator (Dx,Dy)   " << std::endl;
 //   streamlog_message( DEBUG0, xyDir.Print();, std::endl; );


	TMatrixD bFac(2,1);
	TMatrixD BxTMatrix(3,1);
	BxTMatrix.Zero();
	BxTMatrix[0][0] =BxT[0];	BxTMatrix[1][0] =BxT[1];	BxTMatrix[2][0] =BxT[2]; 
 
	bFac = -0.0002998 * (xyDir*BxTMatrix); 
//	std::cout <<std::scientific<< "bFac" << bFac[0][0] << "  ,  " <<  bFac[1][0] << std::endl;
	TMatrixD ajac(5, 5);
	ajac.UnitMatrix();
	if(b.Mag() < 0.001 ){
			ajac[3][2] = ds * std::sqrt(t1w[0] * t1w[0] + t1w[2] * t1w[2]);
			ajac[4][1] = ds;
	}else{
		ajac[1][0] = bFac[0][0]*ds/sinLambda;
//		std::cout<<"Jacobian entries;" << ajac[1][0] << std::endl;
//		std::cout<<"1,0 " << ajac[1][0] << std::endl;
		ajac[2][0] = bFac[1][0]*ds/sinLambda;
//		std::cout<< std::scientific <<"2,0 " << ajac[2][0] << std::endl;
		ajac[3][0] = 0.5*bFac[0][0]*ds*ds;
//		std::cout<<"3,0 " << ajac[3][0] << std::endl;
		ajac[4][0] = 0.5*bFac[1][0]*ds*ds;
//		std::cout<< std::scientific  <<"4,0 " << ajac[4][0] << std::endl;
		ajac[3][1] = ds*sinLambda; 
//		std::cout<<"3,1 " << ajac[3][1] << std::endl;
		ajac[4][2] = ds*sinLambda; 
//		std::cout<<"4,2 " << ajac[4][2] << std::endl;

	}
    streamlog_out( DEBUG0 ) << "Global to Global jacobian: " << std::endl;
    streamlog_message( DEBUG0, ajac.Print();, std::endl; );
	return ajac;
}

//TO DO: This used Z Y X system while claus and other limit jacobian uses Z X Y. 
/* Note that the curvilinear frame that this jacobian has been derived in does not work for particles moving in z-direction.
 * This is due to a construction that assumes tha beam pipe is in the z-direction. Since it is used for collider experiements. 
 * We therefore have to change to the coordinate system used in paper before we apply this jacobian.
 * This is ok since we never access the curvilinear system directly, but always through the local system which is defined 
 * in the local frame of the telescope; i.e Telescope x becomes y, y becomes z and z becomes x.
 */
TMatrixD EUTelNav::getPropagationJacobianCurvilinear(float ds, float qbyp, TVector3 t1w, TVector3 t2w)
{
		//This is needed to change to claus's coordinate system
		TVector3 t1(t1w[2],t1w[1],t1w[0]);
		TVector3 t2(t2w[2],t2w[1],t2w[0]);
		t1.Unit();
		t2.Unit();

		//	if(t1.Mag() != 1.0 or t2.Mag() != 1.0){//TO DO: This statement does not work for some reason.
		//		streamlog_out(MESSAGE9) << "The magnitude of the two vectors is:  "<< t1.Mag()<< " , "<<t2.Mag() << std::endl;
		//		throw(lcio::Exception("The calculation of the jacobian must be performed by unit vectors!.")); 	
		//	}

        const gear::BField& Bfield = geo::gGeometry().getMagneticField();

		//Must transform the B field to be in the frame used in the curvilinear frame.
		///Since field is homogeneous this seems silly but we need to specify a position to geometry to get B-field.
		gear::Vector3D vectorGlobal(0.1,0.1,0.1);

		//Must also change the magnetic field to be in the correct coordinate system
		//We times but 0.3 due to units of other variables. See paper. Must be Tesla
//		const double Bx = (Bfield.at( vectorGlobal ).z())*0.3*pow(10,-3); 
//		const double By = (Bfield.at( vectorGlobal ).y())*0.3*pow(10,-3);
//		const double Bz = (Bfield.at( vectorGlobal ).x())*0.3*pow(10,-3);
		const double Bx = (Bfield.at( vectorGlobal ).z())*pow(10,1);     ////////HERE CHANGE TO EXPRESS THE FIELD IN KGAUSS. CHANGE! 
		const double By = (Bfield.at( vectorGlobal ).y())*pow(10,1);
		const double Bz = (Bfield.at( vectorGlobal ).x())*pow(10,1);
	
		TVector3 b(Bx, By, Bz);
		
		streamlog_out( DEBUG2 ) << "EUTelGeometryTelescopeGeoDescription::getPropagationJacobianCurvilinear()------BEGIN" << std::endl;
		streamlog_out( DEBUG2 ) <<"This is the input to the jacobian"<< std::endl;  
		streamlog_out( DEBUG2 ) <<"The arc length: " <<ds << std::endl;
		streamlog_out( DEBUG2 ) <<"The curvature: "<< qbyp << std::endl; 
		streamlog_out(DEBUG0)<<"The unit momentum start "<< std::endl; 
		streamlog_message( DEBUG0, t1.Print();, std::endl; );
		streamlog_out(DEBUG0)<<"The unit momentum end "<< std::endl; 
		streamlog_message( DEBUG0, t2.Print();, std::endl; );
		streamlog_out(DEBUG0)<<"The unit Magnetic field  "<< std::endl; 
		streamlog_message( DEBUG0, b.Print();, std::endl; );
		
		TMatrixD ajac(5, 5);
		//This is b*c. speed of light in 1 nanosecond
		TVector3 bc = b*0.3*pow(10,-3);     //CHANGE HERE.
		ajac.UnitMatrix(); 
		
		// -|B*c|
		const double qp = -bc.Mag();

		// Q
		const double q = qp * qbyp;
		
		//if q is zero -> line, otherwise a helix
		if (q == 0.)
		{
				ajac[3][2] = ds * sqrt(t1[0] * t1[0] + t1[1] * t1[1]);
				ajac[4][1] = ds;
		}
		else
		{
				// at start
				const double cosl1 = sqrt(t1[0] * t1[0] + t1[1] * t1[1]);
				// at end
				const double cosl2 = sqrt(t2[0] * t2[0] + t2[1] * t2[1]);
				const double cosl2Inv = 1. / cosl2;
				// magnetic field direction
				TVector3 hn(bc.Unit());
				// (signed) momentum
				const double pav = 1.0 / qbyp;
				
//				const double theta = q * ds;
				const double theta = q * ds*0.1;//CONVERT DS TO CENTIMETRES. CHANGE!!!!!!!!!!!

				const double sint = sin(theta);
				const double cost = cos(theta);
				// H*T
				const double gamma = hn.Dot(t2);
				// HxT0
				TVector3 an1 = hn.Cross(t1);
				// HxT
				TVector3 an2 = hn.Cross(t2);
				// U0, V0
				const double au1 = 1. / sqrt(t1[0]*t1[0]+t1[1]*t1[1]);
				TVector3 u1(-au1 * t1[1], au1 * t1[0], 0.);
				TVector3 v1(-t1[2] * u1[1], t1[2] * u1[0], t1[0] * u1[1] - t1[1] * u1[0]);
				// U, V
				const double au2 = 1. /sqrt(t2[0]*t2[0]+t2[1]*t2[1]);
				TVector3 u2(-au2 * t2[1], au2 * t2[0], 0.);
				TVector3 v2(-t2[2] * u2[1], t2[2] * u2[0], t2[0] * u2[1] - t2[1] * u2[0]);
				// N*V = -H*U
				const double anv = -hn.Dot(u2);
				// N*U = H*V
				const double anu = hn.Dot(v2);
				const double omcost = 1. - cost;
				const double tmsint = theta - sint;
				// M0-M
				TVector3 dx(	-(gamma * tmsint * hn[0] + sint * t1[0] + omcost * an1[0]) / q,
								-(gamma * tmsint * hn[1] + sint * t1[1] + omcost * an1[1]) / q,
								-(gamma * tmsint * hn[2] + sint * t1[2] + omcost * an1[2]) / q	);
				// HxU0
				TVector3 hu1 = hn.Cross(u1);
				// HxV0
				TVector3 hv1 = hn.Cross(v1);
				// some.Dot products
				const double u1u2 = u1.Dot(u2), u1v2 = u1.Dot(v2), v1u2 = v1.Dot(u2), v1v2 = v1.Dot(v2);
				const double hu1u2 = hu1.Dot(u2), hu1v2 = hu1.Dot(v2), hv1u2 = hv1.Dot(u2), hv1v2 = hv1.Dot(v2);
				const double hnu1 = hn.Dot(u1), hnv1 = hn.Dot(v1), hnu2 = hn.Dot(u2), hnv2 = hn.Dot(v2);
				const double t2u1 = t2.Dot(u1), t2v1 = t2.Dot(v1);
				const double t2dx = t2.Dot(dx), u2dx = u2.Dot(dx), v2dx = v2.Dot(dx);
				const double an2u1 = an2.Dot(u1), an2v1 = an2.Dot(v1);
				// jacobian
				// 1/P
				ajac[0][0] = 1.;
				// Lambda
				ajac[1][0] = -qp * anv * t2dx;
				ajac[1][1] = cost * v1v2 + sint * hv1v2 + omcost * hnv1 * hnv2 + anv * (-sint * t2v1 + omcost * an2v1 - gamma * tmsint * hnv1);
				ajac[1][2] = cosl1
						* (cost * u1v2 + sint * hu1v2 + omcost * hnu1 * hnv2 + anv * (-sint * t2u1 + omcost * an2u1 - gamma * tmsint * hnu1));
				ajac[1][3] = -q * anv * t2u1;
				ajac[1][4] = -q * anv * t2v1;
				// Phi
				ajac[2][0] = -qp * anu * t2dx * cosl2Inv;
				ajac[2][1] = cosl2Inv
						* (cost * v1u2 + sint * hv1u2 + omcost * hnv1 * hnu2 + anu * (-sint * t2v1 + omcost * an2v1 - gamma * tmsint * hnv1));
				ajac[2][2] = cosl2Inv * cosl1
						* (cost * u1u2 + sint * hu1u2 + omcost * hnu1 * hnu2 + anu * (-sint * t2u1 + omcost * an2u1 - gamma * tmsint * hnu1));
				ajac[2][3] = -q * anu * t2u1 * cosl2Inv;
				ajac[2][4] = -q * anu * t2v1 * cosl2Inv;
				// Xt
				ajac[3][0] = pav * u2dx;
				ajac[3][1] = (sint * v1u2 + omcost * hv1u2 + tmsint * hnu2 * hnv1) / q;
				ajac[3][2] = (sint * u1u2 + omcost * hu1u2 + tmsint * hnu2 * hnu1) * cosl1 / q;
				ajac[3][3] = u1u2;
				ajac[3][4] = v1u2;
				// Yt
				ajac[4][0] = pav * v2dx;
				ajac[4][1] = (sint * v1v2 + omcost * hv1v2 + tmsint * hnv2 * hnv1) / q;
				ajac[4][2] = (sint * u1v2 + omcost * hu1v2 + tmsint * hnv2 * hnu1) * cosl1 / q;
				ajac[4][3] = u1v2;
				ajac[4][4] = v1v2;
		}
		return ajac;
}

//This function determined the xyz position in global coordinates using the state and arc length of the track s.
TVector3 EUTelNav::getPositionfromArcLength(TVector3 pos, TVector3 pVec, float beamQ, double s)
{

		// Get magnetic field vector, assuming uniform magnetic field running along X direction
		gear::Vector3D vectorGlobal( pos[0], pos[1], pos[2] ); 
		const gear::BField&   B = geo::gGeometry().getMagneticField();
		const double bx         = B.at( vectorGlobal ).x();
		const double by         = B.at( vectorGlobal ).y();
		const double bz         = B.at( vectorGlobal ).z();

		TVector3 hVec(bx,by,bz);

		const double H = hVec.Mag();
		const double p = pVec.Mag();
		const double constant = -0.299792458; 
		const double mm = 1000;
		const double combineConstantsAndMagneticField = (constant*beamQ*H)/mm;
		const double k = combineConstantsAndMagneticField;
		const double rho = combineConstantsAndMagneticField/p; 

		if ( fabs( k ) > 0  )
		{
				// Non-zero magnetic field case
				TVector3 pCrossH = pVec.Cross(hVec.Unit());
				TVector3 pCrossHCrossH = pCrossH.Cross(hVec.Unit());
				const double pDotH = pVec.Dot(hVec.Unit());
				TVector3 temp1 = pCrossHCrossH;	temp1 *= ( (-1./k) * sin( rho * s ) );
				TVector3 temp2 = pCrossH;       temp2 *= ( (-1./k) * ( 1. - cos( rho * s ) ) );
				TVector3 temp3 = hVec;          temp3 *= ( (pDotH / p) * s );
				pos += temp1;
				pos += temp2;
				pos += temp3;
		}
		else
		{
				// Vanishing magnetic field case. Here you just determine the fraction of P in each 
				// direction and this must be the fraction of s that this direction gets

				// Calculate cos of the angle between Z(beam) and X(solenoid field axis) //NEED TO MAKE SURE THAT TX=PX/P
				const double cosA = pVec[0]/p;
				// Calculate cos of the angle between Z(beam) and Y
				const double cosB = pVec[1]/p; 
				pos.SetX( pos[0] + cosA * s );
				pos.SetY( pos[1] + cosB * s );
				pos.SetZ( pos[2] + 1./p * pVec.Z() * s );
		}
		return pos;
}
TVector3 EUTelNav::getMomentumfromArcLengthLocal(TVector3 pVec, TVector3 pos, float beamQ, float s, int planeID)
{
	TVector3	newMomentum = EUTelNav::getMomentumfromArcLength(pVec,beamQ,s);
	TVector3 pVecUnitLocal;
	//TO DO: This transform is used also in state. Can make generic transform like this for both.
	double globalVec[] = { newMomentum[0],newMomentum[1],newMomentum[2] };
	double localVec[3];
	if(planeID != 314){
	geo::gGeometry().master2LocalVec( planeID ,globalVec, localVec );
	pVecUnitLocal[0] = localVec[0]; 	pVecUnitLocal[1] = localVec[1]; 	pVecUnitLocal[2] = localVec[2]; 
	}else{
		pVecUnitLocal[0] = pVec[0]; 	pVecUnitLocal[1] = pVec[1]; 	pVecUnitLocal[2] = pVec[2]; 
	}
	return pVecUnitLocal;
}


//This will calculate the momentum at a arc length away given initial parameters.
TVector3 EUTelNav::getMomentumfromArcLength(TVector3 momentum, float charge, float arcLength)
{
		//This is one coordinate axis of curvilinear coordinate system.	
		TVector3 T = momentum.Unit();

		const gear::BField&   Bfield = geo::gGeometry().getMagneticField();
		//Since field is homogeneous this seems silly but we need to specify a position to geometry to get B-field.

		const double Bx = (Bfield.at( TVector3(1,1,1) ).x());//We times bu 0.3 due to units of other variables. See paper. Must be Tesla
		const double By = (Bfield.at( TVector3(1,1,1) ).y());
		const double Bz = (Bfield.at( TVector3(1,1,1) ).z());

		TVector3 B(Bx*0.3, By*0.3, Bz*0.3 );
		TVector3 H = (B.Unit());
		
		const float alpha = (H.Cross(T)).Mag();
		const float gamma = H.Dot(T);

		//You could use end momentum since it must be constant
		const float Q = -(B.Mag())*(charge/(momentum.Mag()));
		//divide by 1000 to convert to meters 
		float theta = (Q*arcLength)/1000;

		TVector3 N = (H.Cross(T)).Unit();
		const float cosTheta = cos(theta);
		const float sinTheta = sin(theta);
		const float oneMinusCosTheta = (1-cos(theta));
		TVector3 momentumEndUnit = gamma*oneMinusCosTheta*H+cosTheta*T+alpha*sinTheta*N;
		TVector3 momentumEnd = momentumEndUnit*(momentum.Mag());
		
		streamlog_out( DEBUG0 ) << "Momentum direction (Unit Vector): " << momentumEndUnit[0] << " , " << momentumEndUnit[1] << " , " << momentumEndUnit[2] << std::endl;
		streamlog_out( DEBUG0 ) << "Momentum: " <<  momentumEnd[0] << " , " << momentumEnd[1] << " , " << momentumEnd[2] << std::endl;

		return momentumEnd;
}

/**
* Find closest surface intersected by the track and return the position
*/
bool EUTelNav::findIntersectionWithCertainID(	float x0, float y0, float z0, 
										float px, float py, float pz, 
										float beamQ, int nextPlaneID, float outputPosition[],
										TVector3& outputMomentum, float& arcLength, int& newNextPlaneID)
{
	//positions are in mm
	TVector3 trkVec(x0,y0,z0);
	TVector3 pVec(px,py,pz);

	//Assuming uniform magnetic field running along X direction. 
	//Why do we need this assumption? Equations of motion do not seem to dictate this.
	gear::Vector3D vectorGlobal( x0, y0, z0 );      
	const gear::BField& B	= geo::gGeometry().getMagneticField();
	const double bx         = B.at( vectorGlobal ).x();
	const double by         = B.at( vectorGlobal ).y();
	const double bz         = B.at( vectorGlobal ).z();

	//B field is in units of Tesla
	TVector3 hVec(bx,by,bz);
	const double H = hVec.Mag();

	//Calculate track momentum from track parameters and fill some useful variables
	const double p = pVec.Mag();
	//This is a constant used in the derivation of this equation. This is the distance light travels in a nano second    
	const double constant =  -0.299792458; 
	const double mm = 1000;
	const double combineConstantsAndMagneticField = (constant*beamQ*H)/mm;
	const double rho = combineConstantsAndMagneticField/p; 
				      
	//Determine geometry of sensor to be used to determine the point of intersection.//////////////////////////////////////
    TVector3 norm = geo::gGeometry().siPlaneNormal( nextPlaneID  );       
	streamlog_out (DEBUG5) << "The normal of the plane is (x,y,z): "<<norm[0]<<","<<norm[1]<<","<<norm[2]<< std::endl;
	TVector3 sensorCenter( geo::gGeometry().siPlaneXPosition( nextPlaneID  ), geo::gGeometry().siPlaneYPosition( nextPlaneID  ), geo::gGeometry().siPlaneZPosition( nextPlaneID  ) );
//	streamlog_out (DEBUG5) << "The sensor centre (x,y,z): "<<sensorCenter[0]<<","<<sensorCenter[1]<<","<<sensorCenter[2] << std::endl;
	TVector3 delta = (trkVec - sensorCenter);//Must be in mm since momentum is.
	TVector3 pVecCrosH = pVec.Cross(hVec.Unit());

	{

		streamlog_out(DEBUG5) << "PROPAGATE FROM (" << trkVec[0] <<","<<trkVec[1]<<","<<trkVec[2]<<") USING: " << std::endl;
		streamlog_out(DEBUG5) << "INPUT: -------------------------------------------" << std::endl;
		streamlog_out(DEBUG5) << "STATE AND ENVIROMENT:---------------------------- " << std::endl;
		streamlog_out(DEBUG5) << "Magnetic field                      :(" << hVec[0] <<","<< hVec[1] <<","<< hVec[2] <<")"<< std::endl;
		streamlog_out(DEBUG5) << "Momentum                            :(" << pVec[0] <<","<< pVec[1] <<","<< pVec[2] <<")"<< std::endl;
		streamlog_out(DEBUG5) << "INFINITE PLANE WE ARE AIMING AT:---------------------------- " << std::endl;
		streamlog_out(DEBUG5) << "Sensor Centre:                      :(" << sensorCenter[0] <<","<< sensorCenter[1] <<","<< sensorCenter[2] <<")"<< std::endl;
		streamlog_out(DEBUG5) << "Normal vector                       :(" << norm[0] <<","<< norm[1] <<","<< norm[2] <<")"<< std::endl;
		streamlog_out(DEBUG5) << "CALCULATE.... " << std::endl;

		streamlog_out(DEBUG5) << "Rho [ (charge.magnetic)/momentum  ] :" << rho << std::endl;
      	streamlog_out(DEBUG5) << "Momentum x Magnetic Field           :(" << pVecCrosH[0] <<","<< pVecCrosH[1] <<","<< pVecCrosH[2] <<")"<< std::endl;
		streamlog_out(DEBUG5) << "Delta(start - sensor centre)        :(" << delta[0] <<","<< delta[1] <<","<< delta[2] <<")"<< std::endl;
 	}

	//Solution to the plane equation and the curved line intersection will be an 
	//quadratic with the coefficients. The solution is the arc length along the curve
	const double a = -0.5 * rho * ( norm.Dot( pVecCrosH ) ) / p;
	const double b = norm.Dot( pVec ) / p;
	const double c = norm.Dot( delta );

	//Solution must be in femto metres, vector of arc length
	std::vector<double> sol = Utility::solveQuadratic(a,b,c); 
	//solutions are sorted in ascending order, choose minimum arc length
	double solution = ( sol[0] > 0. ) ? sol[0] : ( ( sol[0] < 0. && sol[1] > 0. ) ? sol[1] : -1. );
    {
		streamlog_out(DEBUG5) << "FORM QUADRATIC AND SOLVE TO FIND DISTANCE TO INTERSECTION.... " << std::endl;
		streamlog_out(DEBUG5) << "Quadratic coefficients              :(" << a<<","<<b<<","<<c<<")" << std::endl;
      	streamlog_out(DEBUG5) << "Both Possible arc lengths           :(" <<  sol.at(0) <<","<< sol.at(1) <<")"<< std::endl;
		streamlog_out(DEBUG5) << "FINAL ARCLENGTH:                    :("<< solution << std::endl;
    }
	if(solution < 0)
	{
		return false;
	}
			
	//Determine the global position from arc length.             
	TVector3 newPos;
	TVector3 newMomentum;
	newPos = EUTelNav::getPositionfromArcLength(trkVec,pVec,beamQ,solution);
	newMomentum = EUTelNav::getMomentumfromArcLength(pVec, beamQ, solution);
	outputMomentum[0] = newMomentum[0];
	outputMomentum[1] = newMomentum[1];
	outputMomentum[2] = newMomentum[2];

	//Is the new point within the sensor. If not then we may have to propagate a little bit further to enter.
 	double pos[3] = {newPos[0],newPos[1],newPos[2]};
	int sensorIDCheck = geo::gGeometry().getSensorID(pos); 
	bool foundIntersection = false;
	if(sensorIDCheck == nextPlaneID){
		streamlog_out( DEBUG3 ) << "INTERSECTION FOUND! " << std::endl;
		foundIntersection = true;
		outputPosition[0] = newPos[0];
		outputPosition[1] = newPos[1];
		outputPosition[2] = newPos[2];
	}
	if(!foundIntersection)
	{
		foundIntersection = geo::gGeometry().findNextPlaneEntrance( newPos,  newMomentum, nextPlaneID, outputPosition);
	}
	if(!foundIntersection)
	{
		foundIntersection = geo::gGeometry().findNextPlaneEntrance( newPos,  -newMomentum, nextPlaneID, outputPosition);
	}

	arcLength = solution;		

	//We have still not found intersection
	if(!foundIntersection)
	{
		streamlog_out( DEBUG3 ) << "FINAL: NO INTERSECTION FOUND. " << std::endl;
		return false;
	}
	else
	{
		newNextPlaneID = nextPlaneID;
		return true;
	}
}

} //namespace eutelescope
