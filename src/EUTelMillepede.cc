#include "EUTelMillepede.h"

namespace eutelescope {

	EUTelMillepede::EUTelMillepede(){}

	EUTelMillepede::EUTelMillepede(int alignmentMode){
	SetAlignmentMode(alignmentMode);
	}

	EUTelMillepede::~EUTelMillepede(){}

	void EUTelMillepede::SetAlignmentMode(int alignmentMode){

		//This is important since we need to know how millepede numbers theres axis. I think this is the reason for this step. Also set the matrix size //////BEGIN        
		if (alignmentMode==0) {
			streamlog_out(WARNING1) << "No alignment was chosen "<< std::endl;	
			_alignmentMode = Utility::noAlignment;
		} else if (alignmentMode==1) {
    	_alignmentMode = Utility::XYShift;
			_globalLabels->resize(2);
    	_jacobian->ResizeTo(2, 2);
    } else if (alignmentMode==2) {
    	_alignmentMode = Utility::XYShiftXYRot;
  		_globalLabels->resize(3);
    	_jacobian->ResizeTo(2, 3);
 
    } else if (alignmentMode==3) {
    	_alignmentMode = Utility::XYZShiftXYRot;
  	_globalLabels->resize(4);
    _jacobian->ResizeTo(2, 4);

   	} else if (alignmentMode==4) {
    	_alignmentMode = Utility::XYShiftYZRotXYRot;
		_globalLabels->resize(4);
		_jacobian->ResizeTo(2, 4);

   	} else if (alignmentMode==5) {
			_alignmentMode = Utility::XYShiftXZRotXYRot;
		_globalLabels->resize(4);
		_jacobian->ResizeTo(2, 4);

    } else if (alignmentMode==6) {
    	_alignmentMode = Utility::XYShiftXZRotYZRotXYRot;
		_globalLabels->resize(5);
		_jacobian->ResizeTo(2, 5);

    } else if (alignmentMode==7) {
    	_alignmentMode = Utility::XYZShiftXZRotYZRotXYRot;
		_globalLabels->resize(6);
 		_jacobian->ResizeTo(2, 6);

    }else {
    	streamlog_out(WARNING3) << "Alignment mode was not recognized:" << _alignmentMode << std::endl;
      streamlog_out(WARNING3) << "Alignment will not be performed" << std::endl;
      alignmentMode = Utility::noAlignment;
    }
  	_jacobian->Zero();
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////END	
}

	//This will take as input a EUTelState and output the correct jacobian. Need to set the aligment mode before.
	int EUTelMillepede::CreateAlignmentToMeasurementJacobian( EUTelTrackStateImpl* state, TMatrixD* Jacobian ){

		if(state->getHit() == NULL){
			streamlog_out( WARNING1 ) << "This state contains no hit. Will continue but this can not be used for alignment. So why do you want the matrix????" << std::endl;
		}

		if(_alignmentMode == Utility::noAlignment){
			streamlog_out( WARNING1 ) << "You have not set the alignment mode OR have set it to no alignment. Either way this function must end!" << std::endl;
			return -999;
		}

		TVector3 incidenceVecLocal = state->getIncidenceVectorInLocalFrame();

		float TxLocal =  incidenceVecLocal[0]/incidenceVecLocal[2];
		float TyLocal =  incidenceVecLocal[1]/incidenceVecLocal[2];
		const float* localpos = state->getReferencePoint();
		

		CreateAlignmentToMeasurementJacobian( *localpos,*(localpos+1), TxLocal, TyLocal,Jacobian );

	} 

// Jacobian from below looks like //   (-1   0  -y   -dx/dz   x*dx/dz   -y*dx/dz)( x local )
																  //   (0  -1 	x   -dy/dz   x*dy/dz   -y*dy/dz )( y local )
                                  //                                             ( anti clockwise rotation around z) ?? Strange but whatever
                                  //                                             (moving the plane in the z direction)
																	//                                             (this is clockwise rotation look in + y direction )
																	// 	                                           (this is clockwise rotations in the x direction  )
		
int EUTelMillepede::CreateAlignmentToMeasurementJacobian( float x,float y, float slopeXvsZ, float slopeYvsZ, TMatrixD* Jacobian ){			
		
	//////////////////////////////////////Moving the sensor in x and y. Obviously if the sensor move right the hit will appear to move left. Always create this!!!! BEGIN
	Jacobian[0][0] = -1.0; // dxh/dxs      dxh => change in hit position         dxs => Change in sensor position
	Jacobian[0][1] =  0.0; // dxh/dys
	Jacobian[1][0] =  0.0; // dyh/dxs
	Jacobian[1][1] = -1.0; // dyh/dys
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////END


	//////////////////////////////////////////////////////////////Rotations around the z axis of sensor. Only create this if you want rotations around z. BEGIN
	if (_alignmentMode == Utility::XYShiftXYRot
                || _alignmentMode == Utility::XYZShiftXYRot
                || _alignmentMode == Utility::XYShiftXZRotXYRot
                || _alignmentMode == Utility::XYShiftYZRotXYRot
                || _alignmentMode == Utility::XYShiftXZRotYZRotXYRot
                || _alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
		Jacobian[0][2] = -y; // dxh/rotzs   rotzs => rotation of sensor around z axis
		Jacobian[1][2] =  x; // dyh/rotzs
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////END

	///////////////////////////////////////////////////Moving the sensor in the z axis. Only create this if you want shifts in z BEGIN
	if (_alignmentMode == Utility::XYZShiftXYRot
                || _alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
  	Jacobian[0][3] =   -slopeXvsZ; // dxh/dzs
    Jacobian[1][3] =   -slopeYvsZ; // dyh/dzs
  }
	///////////////////////////////////////////////////////////////////////////////////////////END

		
	///////////////////////////////////////////////////////////////Rotations around x and y axis. Only do this if you want to move everything. WHY NOT ALSO PARTIAL SHIFTS?????? BEGIN.
	if (_alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
  	Jacobian[0][4] =   x*slopeXvsZ; // dxh/rotys
    Jacobian[1][4] =   x*slopeYvsZ; // dyh/rotys
    Jacobian[0][5] =  -y*slopeXvsZ; // dxh/rotxs          
    Jacobian[1][5] =  -y*slopeYvsZ; // dyh/rotxs         
  }
	///////////////////////////////////////////////////////////////////////////////////END


//This part is if there is only partial alignment. Therefore you need to overwrite some parts of the full size matrix we have just filled BEGIN
	/////////////////////////////rotation around y axis BEGIN
	if (_alignmentMode == Utility::XYShiftXZRotXYRot) {
		Jacobian[0][3] = x*slopeXvsZ; // dxh/rotys
   	Jacobian[1][3] = x*slopeYvsZ; // dyh/rotys
  }
	///////////////////////////////////////////////////////////////////////Rotation around x axis BEGIN
	if (_alignmentMode == Utility::XYShiftYZRotXYRot) {
  	Jacobian[0][3] = -y*slopeXvsZ; // dxh/rotxs  //Note if changed  the signs here since they were wrong I think. Should match smae calculation above
    Jacobian[1][3] = -y*slopeYvsZ; // dyh/rotxs  
  }
	///////////////////////////////////////////////////////////////////////////////////////////END

 	///////////////This does all rotations but not z shift////////////////////////BEGIN
	if (_alignmentMode == Utility::XYShiftXZRotYZRotXYRot) {
		Jacobian[0][3] =  x*slopeXvsZ; // dxh/rotys
		Jacobian[1][3] =  x*slopeYvsZ; // dyh/rotys
		Jacobian[0][4] = -y*slopeXvsZ; // dxh/rotxs
		Jacobian[1][4] = -y*slopeYvsZ; // dyh/rotxs
  }
	/////////////////////////////////////////////////////////////////////////////END
					//////////////////////////////////////////////////////////////////////////////////////////////////////////END OF PARTIAL MATRIX FILL		
}



} // namespace eutelescope

