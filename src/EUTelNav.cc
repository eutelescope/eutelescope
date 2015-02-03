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
		//We times bu 0.3 due to units of other variables. See paper. Must be Tesla
		const double Bx = (Bfield.at( vectorGlobal ).z())*0.3;
		const double By = (Bfield.at( vectorGlobal ).y())*0.3;
		const double Bz = (Bfield.at( vectorGlobal ).x())*0.3;
		
		TVector3 B(Bx, By, Bz);
		TVector3 H = B.Unit();

		//We transform the momentum to curvilinear frame, as this is what is used to describe the curvilinear frame
		TVector3 curvilinearGlobalMomentum(globalMomentum[2], globalMomentum[1], globalMomentum[0]);
		//With no magnetic field this will point in z direction.	
		TVector3 T = curvilinearGlobalMomentum.Unit();
		
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
				KTelescopeFrame = geo::gGeometry().siPlaneYAxis(planeID);	
				JTelescopeFrame = geo::gGeometry().siPlaneXAxis(planeID);
		}
		else
		{
				ITelescopeFrame.SetXYZ(0,0,1);
				KTelescopeFrame.SetXYZ(0,1,0);
				JTelescopeFrame.SetXYZ(1,0,0);
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
		const double Bx = (Bfield.at( vectorGlobal ).z())*0.3*pow(10,-3); 
		const double By = (Bfield.at( vectorGlobal ).y())*0.3*pow(10,-3);
		const double Bz = (Bfield.at( vectorGlobal ).x())*0.3*pow(10,-3);
		
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
		TVector3 bc = b;
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
				
				const double theta = q * ds;
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
TVector3 EUTelNav::getXYZfromArcLength(TVector3 pos, TVector3 pVec, float beamQ, double s)
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

//This will calculate the momentum at a arc length away given initial parameters.
TVector3 EUTelNav::getXYZMomentumfromArcLength(TVector3 momentum, TVector3 globalPositionStart, float charge, float arcLength)
{
		//This is one coordinate axis of curvilinear coordinate system.	
		TVector3 T = momentum.Unit();

		const gear::BField&   Bfield = geo::gGeometry().getMagneticField();
		//Since field is homogeneous this seems silly but we need to specify a position to geometry to get B-field.
		gear::Vector3D vectorGlobal(globalPositionStart[0],globalPositionStart[1],globalPositionStart[2]);

		const double Bx = (Bfield.at( vectorGlobal ).x());//We times bu 0.3 due to units of other variables. See paper. Must be Tesla
		const double By = (Bfield.at( vectorGlobal ).y());
		const double Bz = (Bfield.at( vectorGlobal ).z());

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

} //namespace eutelescope
