#include "EUTelMillepede.h"

namespace eutelescope {

	EUTelMillepede::EUTelMillepede(){}

	EUTelMillepede::~EUTelMillepede(){}

	void EUTelMillepede::SetAlignmentMode(int alignmentMode){
        
		_alignmentMode = Utility::noAlignment;
		if (alignmentMode==0) {
			_alignmentMode = Utility::noAlignment;
		} else if (alignmentMode==1) {
    	_alignmentMode = Utility::XYShift;
    } else if (alignmentMode==2) {
    	_alignmentMode = Utility::XYShiftXYRot; 
    } else if (alignmentMode==3) {
    	_alignmentMode = Utility::XYZShiftXYRot;
   	} else if (alignmentMode==4) {
    	_alignmentMode = Utility::XYShiftYZRotXYRot;
   	} else if (alignmentMode==5) {
			_alignmentMode = Utility::XYShiftXZRotXYRot;
    } else if (alignmentMode==6) {
    	_alignmentMode = Utility::XYShiftXZRotYZRotXYRot;
    } else if (alignmentMode==7) {
    	_alignmentMode = Utility::XYZShiftXZRotYZRotXYRot;
    }else {
    	streamlog_out(WARNING3) << "Alignment mode was not recognized:" << _alignmentMode << std::endl;
      streamlog_out(WARNING3) << "Alignment will not be performed" << std::endl;
      alignmentMode = Utility::noAlignment;
    }
	}

	

	//This will take as input a EUTelState and output the correct jacobian. Need to set the aligment mode before.
	int EUTelMillepede::CreateAlignmentToMeasurementJacobian( EUTelTrackStateImpl* state, TMatrixD Jacobian ){
		
	/*		/////////////////////////////////////////////////////////////////////////////////////////BEGIN  Using the state information in local coordinates create the matrix
      alDer[0][0] = -1.0; // dx/dx
      alDer[0][1] =  0.0; // dx/dy
      alDer[1][0] =  0.0; // dy/dx
      alDer[1][1] = -1.0; // dy/dy
        globalLabels[0] = _parameterIdXShiftsMap[iPlane]; // dx
        globalLabels[1] = _parameterIdYShiftsMap[iPlane]; // dy


        if (_alignmentMode == Utility::XYShiftXYRot
                || _alignmentMode == Utility::XYZShiftXYRot
                || _alignmentMode == Utility::XYShiftXZRotXYRot
                || _alignmentMode == Utility::XYShiftYZRotXYRot
                || _alignmentMode == Utility::XYShiftXZRotYZRotXYRot
                || _alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
            alDer[0][2] = -predpos[1]; // dx/rot
            alDer[1][2] =  predpos[0]; // dy/rot
            globalLabels[2] = _parameterIdZRotationsMap[iPlane]; // rot z
        }


        if (_alignmentMode == Utility::XYZShiftXYRot
                || _alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
            alDer[0][3] =   -xSlope; // dx/dz
            alDer[1][3] =   -ySlope; // dy/dz
            globalLabels[3] = _parameterIdZShiftsMap[iPlane]; // dz
        }

        if (_alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
            alDer[0][4] =   predpos[0]*xSlope; // dx/rot y
            alDer[1][4] =   predpos[0]*ySlope; // dy/rot y
            globalLabels[4] = _parameterIdYRotationsMap[iPlane]; // drot y  - actually X?
            alDer[0][5] =  -predpos[1]*xSlope; // dx/rot x          
            alDer[1][5] =  -predpos[1]*ySlope; // dy/rot x         
            globalLabels[5] = _parameterIdXRotationsMap[iPlane]; // drot x  - actually Y?
        }


// partial alignment 
        if (_alignmentMode == Utility::XYShiftXZRotXYRot) {
            alDer[0][3] = predpos[0]*xSlope; // dx/rot y
            alDer[1][3] = predpos[0]*ySlope; // dy/rot y
            globalLabels[3] = _parameterIdYRotationsMap[iPlane]; // drot y
        }
        if (_alignmentMode == Utility::XYShiftYZRotXYRot) {
            alDer[0][3] = predpos[1]*xSlope; // dx/rot x
            alDer[1][3] = predpos[1]*ySlope; // dy/rot x
            globalLabels[3] = _parameterIdXRotationsMap[iPlane]; // drot x
        }
 
        if (_alignmentMode == Utility::XYShiftXZRotYZRotXYRot) {
            alDer[0][3] = predpos[0]*xSlope; // dx/rot y
            alDer[1][3] = predpos[0]*ySlope; // dy/rot y
            globalLabels[3] = _parameterIdYRotationsMap[iPlane]; // drot y
            alDer[0][4] = predpos[1]*xSlope; // dx/rot x
            alDer[1][4] = predpos[1]*ySlope; // dy/rot x
            globalLabels[4] = _parameterIdXRotationsMap[iPlane]; // drot x
        }
        if (_alignmentMode == Utility::XYShiftXZRotYZRotXYRot) {
            alDer[0][3] = predpos[0]*xSlope; // dx/rot y
            alDer[1][3] = predpos[0]*ySlope; // dy/rot y
            globalLabels[3] = _parameterIdYRotationsMap[iPlane]; // drot y
            alDer[0][4] = predpos[1]*xSlope; // dx/rot x
            alDer[1][4] = predpos[1]*ySlope; // dy/rot x
            globalLabels[4] = _parameterIdXRotationsMap[iPlane]; // drot x
        }
				//////////////////////////////////////////////////////////////////////////////////////////////////////////END OF CREATION OF MATRI

*/


}


} // namespace eutelescope

