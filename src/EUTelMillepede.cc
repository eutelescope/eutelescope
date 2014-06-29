#include "EUTelMillepede.h"

using namespace lcio;
using namespace std;
using namespace marlin;
using namespace eutelescope;



namespace eutelescope {

	//This constructor is useful for running the binary files execution
	EUTelMillepede::EUTelMillepede(){
	FillMilleParametersLabels();
	}

	//This constructor useful for mille binary output part
	EUTelMillepede::EUTelMillepede(int alignmentMode){
	SetAlignmentMode(alignmentMode);
	FillMilleParametersLabels();
	}

	EUTelMillepede::~EUTelMillepede(){}

	void EUTelMillepede::SetAlignmentMode(int alignmentMode){

		//This is important since we need to know how millepede numbers theres axis. I think this is the reason for this step. Also set the matrix size //////BEGIN        
		if (alignmentMode==0) {
			streamlog_out(WARNING1) << "No alignment was chosen "<< std::endl;	
			_alignmentMode = Utility::noAlignment;
		} else if (alignmentMode==1) {
    	_alignmentMode = Utility::XYShift;
			_globalLabels.resize(2);
    	_jacobian->ResizeTo(2, 2);
    } else if (alignmentMode==2) {
    	_alignmentMode = Utility::XYShiftXYRot;
  		_globalLabels.resize(3);
    	_jacobian->ResizeTo(2, 3);
 
    } else if (alignmentMode==3) {
    	_alignmentMode = Utility::XYZShiftXYRot;
  	_globalLabels.resize(4);
    _jacobian->ResizeTo(2, 4);

   	} else if (alignmentMode==4) {
    	_alignmentMode = Utility::XYShiftYZRotXYRot;
		_globalLabels.resize(4);
		_jacobian->ResizeTo(2, 4);

   	} else if (alignmentMode==5) {
			_alignmentMode = Utility::XYShiftXZRotXYRot;
		_globalLabels.resize(4);
		_jacobian->ResizeTo(2, 4);

    } else if (alignmentMode==6) {
    	_alignmentMode = Utility::XYShiftXZRotYZRotXYRot;
		_globalLabels.resize(5);
		_jacobian->ResizeTo(2, 5);

    } else if (alignmentMode==7) {
    	_alignmentMode = Utility::XYZShiftXZRotYZRotXYRot;
		_globalLabels.resize(6);
 		_jacobian->ResizeTo(2, 6);

    }else {
    	streamlog_out(WARNING3) << "Alignment mode was not recognized:" << _alignmentMode << std::endl;
      streamlog_out(WARNING3) << "Alignment will not be performed" << std::endl;
      alignmentMode = Utility::noAlignment;
    }
  	_jacobian->Zero();
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////END	
}

//Note here we label sensors and every alignment degree of freedom uniquely. Note that even if the sensor is to remain fixed. The fixing is done latter.
void EUTelMillepede::FillMilleParametersLabels() {

    int currentLabel = 0;
    const IntVec sensorIDsVec = geo::gGeometry().sensorIDsVec();
    IntVec::const_iterator itr;
    for( itr = sensorIDsVec.begin(); itr != sensorIDsVec.end(); ++itr ) {
        _xShiftsMap.insert( make_pair(*itr, ++currentLabel) );
    }
    for( itr = sensorIDsVec.begin(); itr != sensorIDsVec.end(); ++itr ) {
        _yShiftsMap.insert( make_pair(*itr, ++currentLabel) );
    }
    for( itr = sensorIDsVec.begin(); itr != sensorIDsVec.end(); ++itr ) {
        _zShiftsMap.insert( make_pair(*itr, ++currentLabel) );
    }
    for( itr = sensorIDsVec.begin(); itr != sensorIDsVec.end(); ++itr ) {
        _xRotationsMap.insert( make_pair(*itr, ++currentLabel) );
    }
    for( itr = sensorIDsVec.begin(); itr != sensorIDsVec.end(); ++itr ) {
        _yRotationsMap.insert( make_pair(*itr, ++currentLabel) );
    }
    for( itr = sensorIDsVec.begin(); itr != sensorIDsVec.end(); ++itr ) {
        _zRotationsMap.insert( make_pair(*itr, ++currentLabel) );
    }
}

	//This will take as input a EUTelState and output the correct jacobian. Need to set the aligment mode before.
	int EUTelMillepede::CreateAlignmentToMeasurementJacobian( EUTelTrackStateImpl* state){

		if(state->getHit() == NULL){
			streamlog_out( WARNING1 ) << "This state contains no hit. Will continue but this can not be used for alignment. So why do you want the matrix????" << std::endl;
		}

		if(_alignmentMode == Utility::noAlignment){
			streamlog_out( WARNING1 ) << "You have not set the alignment mode OR have set it to no alignment. Either way this function must end!" << std::endl;
			return -999;
		}

		TVector3 incidenceVecLocal = state->getIncidenceVectorInLocalFrame();

		streamlog_out( DEBUG0 ) << "Incidence vector in local frame"<<incidenceVecLocal[0] <<","<< incidenceVecLocal[1]<<","<<incidenceVecLocal[2] << std::endl;

		float TxLocal =  incidenceVecLocal[0]/incidenceVecLocal[2];
		float TyLocal =  incidenceVecLocal[1]/incidenceVecLocal[2];
		const float* localpos = state->getReferencePoint();

			streamlog_out( DEBUG0 ) << "This is px/pz, py/pz (local) "<< TxLocal <<","<< TyLocal << std::endl;
			streamlog_out( DEBUG0 ) << "Local frame position "<< *localpos<<","<<*(localpos+1)<<","<<*(localpos+2) << std::endl;
		

		CreateAlignmentToMeasurementJacobian( *localpos,*(localpos+1), TxLocal, TyLocal);

	} 

// _jacobian from below looks like //   (-1   0  -y   -dx/dz   x*dx/dz   -y*dx/dz)( x local )
																  //   (0  -1 	x   -dy/dz   x*dy/dz   -y*dy/dz )( y local )
                                  //                                             ( anti clockwise rotation around z) ?? Strange but whatever
                                  //                                             (moving the plane in the z direction)
																	//                                             (this is clockwise rotation look in + y direction )
																	// 	                                           (this is clockwise rotations in the x direction  )
		
int EUTelMillepede::CreateAlignmentToMeasurementJacobian( float x,float y, float slopeXvsZ, float slopeYvsZ){
	_jacobian->Zero();			
		
	//////////////////////////////////////Moving the sensor in x and y. Obviously if the sensor move right the hit will appear to move left. Always create this!!!! BEGIN
	_jacobian[0][0] = -1.0; // dxh/dxs      dxh => change in hit position         dxs => Change in sensor position
	_jacobian[0][1] =  0.0; // dxh/dys
	_jacobian[1][0] =  0.0; // dyh/dxs
	_jacobian[1][1] = -1.0; // dyh/dys
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////END


	//////////////////////////////////////////////////////////////Rotations around the z axis of sensor. Only create this if you want rotations around z. BEGIN
	if (_alignmentMode == Utility::XYShiftXYRot
                || _alignmentMode == Utility::XYZShiftXYRot
                || _alignmentMode == Utility::XYShiftXZRotXYRot
                || _alignmentMode == Utility::XYShiftYZRotXYRot
                || _alignmentMode == Utility::XYShiftXZRotYZRotXYRot
                || _alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
		_jacobian[0][2] = -y; // dxh/rotzs   rotzs => rotation of sensor around z axis
		_jacobian[1][2] =  x; // dyh/rotzs
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////END

	///////////////////////////////////////////////////Moving the sensor in the z axis. Only create this if you want shifts in z BEGIN
	if (_alignmentMode == Utility::XYZShiftXYRot
                || _alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
  	_jacobian[0][3] =   -slopeXvsZ; // dxh/dzs
    _jacobian[1][3] =   -slopeYvsZ; // dyh/dzs
  }
	///////////////////////////////////////////////////////////////////////////////////////////END

		
	///////////////////////////////////////////////////////////////Rotations around x and y axis. Only do this if you want to move everything. WHY NOT ALSO PARTIAL SHIFTS?????? BEGIN.
	if (_alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
  	_jacobian[0][4] =   x*slopeXvsZ; // dxh/rotys
    _jacobian[1][4] =   x*slopeYvsZ; // dyh/rotys
    _jacobian[0][5] =  -y*slopeXvsZ; // dxh/rotxs          
    _jacobian[1][5] =  -y*slopeYvsZ; // dyh/rotxs         
  }
	///////////////////////////////////////////////////////////////////////////////////END


//This part is if there is only partial alignment. Therefore you need to overwrite some parts of the full size matrix we have just filled BEGIN
	/////////////////////////////rotation around y axis BEGIN
	if (_alignmentMode == Utility::XYShiftXZRotXYRot) {
		_jacobian[0][3] = x*slopeXvsZ; // dxh/rotys
   	_jacobian[1][3] = x*slopeYvsZ; // dyh/rotys
  }
	///////////////////////////////////////////////////////////////////////Rotation around x axis BEGIN
	if (_alignmentMode == Utility::XYShiftYZRotXYRot) {
  	_jacobian[0][3] = -y*slopeXvsZ; // dxh/rotxs  //Note if changed  the signs here since they were wrong I think. Should match smae calculation above
    _jacobian[1][3] = -y*slopeYvsZ; // dyh/rotxs  
  }
	///////////////////////////////////////////////////////////////////////////////////////////END

 	///////////////This does all rotations but not z shift////////////////////////BEGIN
	if (_alignmentMode == Utility::XYShiftXZRotYZRotXYRot) {
		_jacobian[0][3] =  x*slopeXvsZ; // dxh/rotys
		_jacobian[1][3] =  x*slopeYvsZ; // dyh/rotys
		_jacobian[0][4] = -y*slopeXvsZ; // dxh/rotxs
		_jacobian[1][4] = -y*slopeYvsZ; // dyh/rotxs
  }
	/////////////////////////////////////////////////////////////////////////////END
					//////////////////////////////////////////////////////////////////////////////////////////////////////////END OF PARTIAL MATRIX FILL		
}

//This might not need to separate but just incase we want to fill this without the need for the jacobian. Unlikely but I like to keep different processes separated

void EUTelMillepede::CreateGlobalLabels(EUTelTrackStateImpl* state){
	CreateGlobalLabels( state->getLocation());
}

void EUTelMillepede::CreateGlobalLabels( int iPlane){
	_globalLabels.clear();

	_globalLabels[0] = _xShiftsMap[iPlane]; // dx
	_globalLabels[1] = _yShiftsMap[iPlane]; // dy


	if (_alignmentMode == Utility::XYShiftXYRot
                || _alignmentMode == Utility::XYZShiftXYRot
                || _alignmentMode == Utility::XYShiftXZRotXYRot
                || _alignmentMode == Utility::XYShiftYZRotXYRot
                || _alignmentMode == Utility::XYShiftXZRotYZRotXYRot
                || _alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
  	_globalLabels[2] = _zRotationsMap[iPlane]; // rot z
  }


	if (_alignmentMode == Utility::XYZShiftXYRot
                || _alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
		_globalLabels[3] = _zRotationsMap[iPlane]; // dz
  }

	if (_alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
		_globalLabels[4] = _yRotationsMap[iPlane]; // drot y  - actually X?       
		_globalLabels[5] = _xRotationsMap[iPlane]; // drot x  - actually Y?
  }


// partial alignment 
	if (_alignmentMode == Utility::XYShiftXZRotXYRot) {
		_globalLabels[3] = _yRotationsMap[iPlane]; // drot y
  }
	if (_alignmentMode == Utility::XYShiftYZRotXYRot) {
		_globalLabels[3] = _xRotationsMap[iPlane]; // drot x
	}
 
	if (_alignmentMode == Utility::XYShiftXZRotYZRotXYRot) {
		_globalLabels[3] = _yRotationsMap[iPlane]; // drot y
		_globalLabels[4] = _xRotationsMap[iPlane]; // drot x
  }

	if (_alignmentMode == Utility::XYShiftXZRotYZRotXYRot) {
		_globalLabels[3] = _yRotationsMap[iPlane]; // drot y
		_globalLabels[4] = _xRotationsMap[iPlane]; // drot x
  }


}
//////////////////////////////////////////////////////////////////////////////////////////////////////All these functions are used after binary file creation.
void EUTelMillepede::setXShiftFixed(lcio::IntVec xfixed){
	_fixedAlignmentXShfitPlaneIds = xfixed;
}

void EUTelMillepede::setYShiftFixed(lcio::IntVec yfixed){
	_fixedAlignmentYShfitPlaneIds = yfixed;
}

void EUTelMillepede::setZShiftFixed(lcio::IntVec zfixed){
	_fixedAlignmentZShfitPlaneIds = zfixed;
}

void EUTelMillepede::setXRotationsFixed(lcio::IntVec xRotfixed){
	_fixedAlignmentXRotationPlaneIds = xRotfixed;
}

void EUTelMillepede::setYRotationsFixed(lcio::IntVec yRotfixed){
	_fixedAlignmentYRotationPlaneIds = yRotfixed;
}

void EUTelMillepede::setZRotationsFixed(lcio::IntVec zRotfixed){
	_fixedAlignmentZRotationPlaneIds = zRotfixed;
}

void EUTelMillepede::setPlanesExclude(lcio::IntVec exclude){
	_alignmentPlaneIdsExclude = exclude;
}

void EUTelMillepede::setBinaryFileName(std::string binary){
	_milleBinaryFilename = binary;
}



void EUTelMillepede::setSteeringFileName(std::string name){

	_milleSteeringFilename = name;
}

void EUTelMillepede::setResultsFileName(std::string name){

	_milleResultFileName = name;

}




  
int EUTelMillepede::writeMilleSteeringFile(lcio::StringVec pedeSteerAddCmds){

streamlog_out(DEBUG2) << "EUTelMillepede::writeMilleSteeringFile------------------------------------BEGIN" << endl;

	// Prepare millepede steering files only if alignment was requested //////////////Check the alignment and that you can open the steering file name/////////BEGIN
	if (_alignmentMode == Utility::noAlignment) {
		streamlog_out(WARNING1) << "Alignment steering file will not be created" << endl;
		return -999;
	}
	
	ofstream steerFile;
  steerFile.open(_milleSteeringFilename.c_str());

	if (!steerFile.is_open()) {
		streamlog_out(ERROR2) << "Could not open steering file." << _milleSteeringFilename << endl;
		return -999;
  }
	//////////////////////////////////////////////////////////////////////////////////////////////////END 

	
	steerFile << "Cfiles" << endl;
  steerFile << _milleBinaryFilename << endl;
  steerFile << endl;
    //
  steerFile << "Parameter" << endl;

  int counter = 0;

////////////////////////////////////////////////////////////////////////// loop over all planes BEGIN
// @TODO assumes that planes have ids 0..._nplanes !generaly wrong            
	for (unsigned int help = 0; help < geo::gGeometry().nPlanes(); help++) {

///////////////////////////////////////////////////////////////////////////////////////////////////Determine if the plane should be excluded or if some of alignment parameters are fixed. BEGIN
	const int sensorId = geo::gGeometry().sensorZOrderToID(help);
  const bool isPlaneExcluded = std::find(_alignmentPlaneIdsExclude.begin(), _alignmentPlaneIdsExclude.end(), sensorId) == _alignmentPlaneIdsExclude.end();
        
  // check if plane has to be used as fixed
  const bool isFixedXShift = std::find(_fixedAlignmentXShfitPlaneIds.begin(), _fixedAlignmentXShfitPlaneIds.end(), sensorId) != _fixedAlignmentXShfitPlaneIds.end();
  const bool isFixedYShift = std::find(_fixedAlignmentYShfitPlaneIds.begin(), _fixedAlignmentYShfitPlaneIds.end(), sensorId) != _fixedAlignmentYShfitPlaneIds.end();
  const bool isFixedZShift = std::find(_fixedAlignmentZShfitPlaneIds.begin(), _fixedAlignmentZShfitPlaneIds.end(), sensorId) != _fixedAlignmentZShfitPlaneIds.end();
  const bool isFixedXRotation = std::find(_fixedAlignmentXRotationPlaneIds.begin(), _fixedAlignmentXRotationPlaneIds.end(), sensorId) != _fixedAlignmentXRotationPlaneIds.end();
  const bool isFixedYRotation = std::find(_fixedAlignmentYRotationPlaneIds.begin(), _fixedAlignmentYRotationPlaneIds.end(), sensorId) != _fixedAlignmentYRotationPlaneIds.end();
  const bool isFixedZRotation = std::find(_fixedAlignmentZRotationPlaneIds.begin(), _fixedAlignmentZRotationPlaneIds.end(), sensorId) != _fixedAlignmentZRotationPlaneIds.end();

////////////////////////////////////////////////////////////////////////////////////////////////////END

  // if plane not excluded
  if ( !isPlaneExcluded ) {
///////////////////////////////////////////////////////////////////////////////////////////Now fill string that will go into steering depending on if fixed or not BEGIN
		const string initUncertaintyXShift = (isFixedXShift) ? "-1." : "0.01";
    const string initUncertaintyYShift = (isFixedYShift) ? "-1." : "0.01";
    const string initUncertaintyZShift = (isFixedZShift) ? "-1." : "0.01";
    const string initUncertaintyXRotation = (isFixedXRotation) ? "-1." : "0.01";
    const string initUncertaintyYRotation = (isFixedYRotation) ? "-1." : "0.01";
    const string initUncertaintyZRotation = (isFixedZRotation) ? "-1." : "0.01";

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////END


/*  Not sure what this part does. Need to fix is later.        
    const double initXshift = (isFixedXShift) ? 0. : _seedAlignmentConstants._xResiduals[sensorId]/_seedAlignmentConstants._nxResiduals[sensorId];
    const double initYshift = (isFixedYShift) ? 0. : _seedAlignmentConstants._yResiduals[sensorId]/_seedAlignmentConstants._nyResiduals[sensorId];
  */          

const double initXshift =0; const double initYshift = 0;
            if( _alignmentMode==Utility::XYZShiftXYRot ) {
                steerFile << left << setw(25) << _xShiftsMap[sensorId] << setw(25) << -initXshift << setw(25) << initUncertaintyXShift
                           << setw(25) << " ! X shift " << setw(25) << sensorId << endl;
                steerFile << left << setw(25) << _yShiftsMap[sensorId] << setw(25) << -initYshift << setw(25) << initUncertaintyYShift
                           << setw(25) << " ! Y shift " << setw(25) << sensorId << endl;
                steerFile << left << setw(25) << _zShiftsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyZShift
                           << setw(25) << " ! Z shift " << setw(25) << sensorId << endl;
                steerFile << left << setw(25) << _zRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyZRotation
                          << setw(25) << " ! XY rotation " << sensorId << endl;
            } else if( _alignmentMode==Utility::XYShiftYZRotXYRot ) {
                steerFile << left << setw(25) << _xShiftsMap[sensorId] << setw(25) << -initXshift << setw(25) << initUncertaintyXShift
                          << setw(25) << " ! X shift " << sensorId << endl;
                steerFile << left << setw(25) << _yShiftsMap[sensorId] << setw(25)  << -initYshift << setw(25) << initUncertaintyYShift
                          << setw(25) << " ! Y shift " << sensorId << endl;
                steerFile << left << setw(25) << _xRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyXRotation
                          << setw(25) << " ! YZ rotation " << sensorId << endl;
                steerFile << left << setw(25) << _zRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyZRotation
                          << setw(25) << " ! XY rotation " << sensorId << endl;
            } else if( _alignmentMode==Utility::XYShiftXZRotXYRot ) {
                steerFile << left << setw(25) << _xShiftsMap[sensorId] << setw(25) << -initXshift << setw(25) << initUncertaintyXShift
                          << setw(25) << " ! X shift " << sensorId << endl;
                steerFile << left << setw(25) << _yShiftsMap[sensorId] << setw(25) << -initYshift << setw(25) << initUncertaintyYShift
                          << setw(25) << " ! Y shift " << sensorId << endl;
                steerFile << left << setw(25) << _yRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyYRotation
                          << setw(25) << " ! XZ rotation " << sensorId << endl;
                steerFile << left << setw(25) << _zRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyZRotation
                         << setw(25)  << " ! XY rotation " << sensorId << endl;
            } else if( _alignmentMode==Utility::XYShiftXZRotYZRotXYRot ) {
                steerFile << left << setw(25) << _xShiftsMap[sensorId] << setw(25) << -initXshift << setw(25) << initUncertaintyXShift
                          << setw(25) << " ! X shift " << sensorId << endl;
                steerFile << left << setw(25) << _yShiftsMap[sensorId] << setw(25) << -initYshift << setw(25) << initUncertaintyYShift
                          << setw(25) << " ! Y shift " << sensorId << endl;
                steerFile << left << setw(25) << _yRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyYRotation
                          << setw(25) << " ! XZ rotation " << sensorId << endl;
                steerFile << left << setw(25) << _xRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyXRotation
                          << setw(25) << " ! YZ rotation " << sensorId << endl;
                steerFile << left << setw(25) << _zRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyZRotation
                         << setw(25)  << " ! XY rotation " << sensorId << endl;
            } else if( _alignmentMode==Utility::XYZShiftXZRotYZRotXYRot ) {
                steerFile << left << setw(25) << _xShiftsMap[sensorId] << setw(25) << -initXshift << setw(25) << initUncertaintyXShift
                          << setw(25) << " ! X shift " << sensorId << endl;
                steerFile << left << setw(25) << _yShiftsMap[sensorId] << setw(25) << -initYshift << setw(25) << initUncertaintyYShift
                          << setw(25) << " ! Y shift " << sensorId << endl;
                steerFile << left << setw(25) << _zShiftsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyZShift
                          << setw(25) << " ! Z shift " << sensorId << endl;
                steerFile << left << setw(25) << _yRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyYRotation
                          << setw(25) << " ! XZ rotation " << sensorId << endl;
                steerFile << left << setw(25) << _xRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyXRotation
                          << setw(25) << " ! YZ rotation " << sensorId << endl;
                steerFile << left << setw(25) << _zRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyZRotation
                         << setw(25)  << " ! XY rotation " << sensorId << endl;
            } else if ( _alignmentMode==Utility::XYShiftXYRot ) {
                steerFile << left << setw(25) << _xShiftsMap[sensorId] << setw(25) << -initXshift << setw(25) << initUncertaintyXShift
                          << setw(25) << " ! X shift " << sensorId << endl;
                steerFile << left << setw(25) << _yShiftsMap[sensorId] << setw(25) << -initYshift << setw(25) << initUncertaintyYShift
                          << setw(25) << " ! Y shift " << sensorId << endl;
                steerFile << left << setw(25) << _zRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyZRotation
                          << setw(25) << " ! XY rotation " << sensorId << endl;
            } else if ( _alignmentMode==Utility::XYShift ) {
                steerFile << left << setw(25) << _xShiftsMap[sensorId] << setw(25) << -initXshift << setw(25) << initUncertaintyXShift
                          << setw(25) << " ! X shift " << sensorId << endl;
                steerFile << left << setw(25) << _yShiftsMap[sensorId] << setw(25) << -initYshift << setw(25) << initUncertaintyYShift
                          << setw(25) << " ! Y shift " << sensorId << endl;
                steerFile << left << setw(25) << _zRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << "-1.0"
                          << setw(25) << " ! XY rotation fixed" << sensorId << endl;
            }  

            counter++;

        } // end if plane not excluded

    } // end loop over all planes

//    steerFile << "method diagonalization 15 0.1" << endl;
//    steerFile << "hugecut 500." << endl;
//    steerFile << "!chiscut 50. 25." << endl;
//    steerFile << "outlierdownweighting 4" << endl;
//    steerFile << "dwfractioncut 0.2" << endl;

    steerFile << endl;
    for ( StringVec::iterator it = pedeSteerAddCmds.begin( ); it != pedeSteerAddCmds.end( ); ++it ) {
        // two backslashes will be interpreted as newline
        if ( *it == "\\\\" )
            steerFile << endl;
        else
            steerFile << *it << " ";
    }
    steerFile << endl;
    steerFile << "end" << endl;

    steerFile.close();

    if( _alignmentMode != Utility::noAlignment ) streamlog_out(MESSAGE5) << "File " << _milleSteeringFilename << " written." << endl;   

  
}

int EUTelMillepede::runPede(){
	// check if alignment was requested
	if ( _alignmentMode == Utility::noAlignment ) {
		streamlog_out( WARNING1 ) << "No alignment mode has been set. Must end" << endl;
        return -999;
  }

	std::string command = "pede " + _milleSteeringFilename;
	streamlog_out ( MESSAGE5 ) << "Starting pede...: " << command.c_str( ) << endl;

  // run pede and create a streambuf that reads its stdout and stderr
  redi::ipstream pede( command.c_str( ), redi::pstreams::pstdout | redi::pstreams::pstderr );

	if ( !pede.is_open( ) ) {
		streamlog_out( ERROR5 ) << "Pede cannot be executed: command not found in the path" << endl;
  } else {

		// Currently unused variable:
    //bool encounteredError = false;

    // output multiplexing: parse pede output in both stdout and stderr and echo messages accordingly
    char buf[1024];
		std::streamsize n;
    std::stringstream pedeoutput; // store stdout to parse later
    std::stringstream pedeerrors;
    bool finished[2] = { false, false };
    while ( !finished[0] || !finished[1] ) {
    	if ( !finished[0] ) {
    		while ( ( n = pede.err( ).readsome( buf, sizeof (buf ) ) ) > 0 ) {
    			streamlog_out( ERROR5 ).write( buf, n ).flush( );
          string error ( buf, n );
          pedeerrors << error;
                    //encounteredError = true;
        }
       if ( pede.eof( ) ) {
       	finished[0] = true;
       		if ( !finished[1] ) pede.clear( );
       }
    	}

      if ( !finished[1] ) {
                while ( ( n = pede.out( ).readsome( buf, sizeof (buf ) ) ) > 0 ) {
                    streamlog_out( MESSAGE4 ).write( buf, n ).flush( );
                    string output ( buf, n );
                    pedeoutput << output;
                }
                if ( pede.eof( ) ) {
                    finished[1] = true;
                    if ( !finished[0] )
                        pede.clear( );
                }
            }
        }
        // wait for the pede execution to finish
        pede.close( );
        
        // Parse and rename MILLEPEDE result file
      //  if ( parseMilleOutput( "millepede.res" ) ) //moveMilleResultFile( "millepede.res", _milleResultFileName );
    }//END OF IF STATEMENT

return 0;

}

bool EUTelMillepede::parseMilleOutput(std::string alignmentConstantLCIOFile, std::string gear_aligned_file){
   
    bool isOK = true;
    
    // Check if the file is avaliable
    ifstream file( _milleResultFileName.c_str() );
    if ( !file.good( ) ) {
        streamlog_out( WARNING2 ) << "Can't read/find " << _milleResultFileName << " in current directory." << endl;
        isOK = false;
        return isOK;
    }
 
   
    const string command = "parsemilleout.sh " + _milleSteeringFilename + " " + _milleResultFileName + " " + alignmentConstantLCIOFile + 
                           " " + Global::parameters->getStringVal("GearXMLFile" ) + " " + gear_aligned_file;
    streamlog_out ( MESSAGE5 ) << "Convering millepede results to LCIO collections... " << endl;
    streamlog_out ( MESSAGE5 ) << command << endl;

    // run pede and create a streambuf that reads its stdout and stderr
    redi::ipstream parsepede( command.c_str( ), redi::pstreams::pstdout | redi::pstreams::pstderr );

    if ( !parsepede.is_open( ) ) {
        streamlog_out( ERROR5 ) << "Pede cannot be executed: command not found in the path" << endl;
    } else {
        // Currently unused variable:
        // bool encounteredError = false;
        // output multiplexing: parse parsepede output in both stdout and stderr and echo messages accordingly
        char buf[1024];
        std::streamsize n;
        std::stringstream parsepedeoutput; // store stdout to parse later
        std::stringstream parsepedeerrors;
        bool finished[2] = { false, false };
        while ( !finished[0] || !finished[1] ) {
            if ( !finished[0] ) {
                while ( ( n = parsepede.err( ).readsome( buf, sizeof (buf ) ) ) > 0 ) {
                    streamlog_out( ERROR5 ).write( buf, n ).flush( );
                    string error ( buf, n );
                    parsepedeerrors << error;
                    //encounteredError = true;
                }
                if ( parsepede.eof( ) ) {
                    finished[0] = true;
                    if ( !finished[1] )
                        parsepede.clear( );
                }
            }

            if ( !finished[1] ) {
                while ( ( n = parsepede.out( ).readsome( buf, sizeof (buf ) ) ) > 0 ) {
                    streamlog_out( MESSAGE4 ).write( buf, n ).flush( );
                    string output ( buf, n );
                    parsepedeoutput << output;
                }
                if ( parsepede.eof( ) ) {
                    finished[1] = true;
                    if ( !finished[0] )
                        parsepede.clear( );
                }
            }
        }
        // wait for the parsepede execution to finish
        parsepede.close( );
    }
    
    return isOK;
    

}


} // namespace eutelescope




