#include "EUTelMillepede.h"

using namespace lcio;
using namespace std;
using namespace marlin;
using namespace eutelescope;


//TO DO: The last alignment mode does not work. I.e alignment mode 7. This is not a big deal since we will never align all degrees of freedom in one go.
namespace eutelescope {

	//This constructor is useful for running the binary files execution
	EUTelMillepede::EUTelMillepede() : 
	_milleGBL(NULL),
	_alignmentMode(Utility::noAlignment),
	_jacobian(5,5),
	_globalLabels(5)
 	{
	FillMilleParametersLabels();
	}

	//This constructor useful for mille binary output part
	EUTelMillepede::EUTelMillepede(int alignmentMode) :
	_milleGBL(NULL),
	_alignmentMode(Utility::noAlignment),
	_jacobian(5,5),
	_globalLabels(5)
	{
	SetAlignmentMode(alignmentMode);
	FillMilleParametersLabels();
	CreateBinary();
	}

	EUTelMillepede::~EUTelMillepede(){}

	//Depending on the number alignmentMode(This tell us what shift/rotation we can do), we create the alignment jacobain here with the correct dimensions. 
	//Furthermore we determine the size of the labels per plane.
	//Note this is used for all sensors but we can still fix some alignement parameters independent of this.
	void EUTelMillepede::SetAlignmentMode(int alignmentMode){

		//This is important since we need to know how millepede numbers theres axis. I think this is the reason for this step. Also set the matrix size //////BEGIN        
		if (alignmentMode==0) {
			streamlog_out(WARNING1) << "No alignment was chosen "<< std::endl;	
			_alignmentMode = Utility::noAlignment;
		} else if (alignmentMode==1) {
    	_alignmentMode = Utility::XYShift;
			_globalLabels.resize(2);
    	_jacobian.ResizeTo(2, 2);
    } else if (alignmentMode==2) {
    	_alignmentMode = Utility::XYShiftXYRot;
  		_globalLabels.resize(3);
    	_jacobian.ResizeTo(2, 3);
 
    } else if (alignmentMode==3) {
    	_alignmentMode = Utility::XYZShiftXYRot;
  	_globalLabels.resize(4);
    _jacobian.ResizeTo(2, 4);

   	} else if (alignmentMode==4) {
    	_alignmentMode = Utility::XYShiftYZRotXYRot;
		_globalLabels.resize(4);
		_jacobian.ResizeTo(2, 4);

   	} else if (alignmentMode==5) {
			_alignmentMode = Utility::XYShiftXZRotXYRot;
		_globalLabels.resize(4);
		_jacobian.ResizeTo(2, 4);

    } else if (alignmentMode==6) {
    	_alignmentMode = Utility::XYShiftXZRotYZRotXYRot;
		_globalLabels.resize(5);
		_jacobian.ResizeTo(2, 5);

    } else if (alignmentMode==7) {
    	_alignmentMode = Utility::XYZShiftXZRotYZRotXYRot;
		_globalLabels.resize(6);
 		_jacobian.ResizeTo(2, 6);

    }else {
			throw(lcio::Exception(Utility::outputColourString(" Alignment mode was not recognized. ", "RED")));
    }
  	_jacobian.Zero();
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////END	
}

//Note here we label sensors and every alignment degree of freedom uniquely. Note that even if the sensor is to remain fixed. The fixing is done latter.
void EUTelMillepede::FillMilleParametersLabels() {

    int currentLabel = 0;
    const IntVec sensorIDsVec = geo::gGeometry().sensorIDsVec();
    IntVec::const_iterator itr;
    for( itr = sensorIDsVec.begin(); itr != sensorIDsVec.end(); ++itr ) {//sensor 0 to 5 will have numbers 1 to 6 for this x shift
        _xShiftsMap.insert( make_pair(*itr, ++currentLabel) );
    }
    for( itr = sensorIDsVec.begin(); itr != sensorIDsVec.end(); ++itr ) {// sensor 0 to 5 will have numbers 7 to 12 for this y shift
        _yShiftsMap.insert( make_pair(*itr, ++currentLabel) );
    }
    for( itr = sensorIDsVec.begin(); itr != sensorIDsVec.end(); ++itr ) {// sensor 0 to 5 will have numbers 13 to 18 for this z shift
        _zShiftsMap.insert( make_pair(*itr, ++currentLabel) );
    }
    for( itr = sensorIDsVec.begin(); itr != sensorIDsVec.end(); ++itr ) {// sensor 0 to 5 will have numbers 19 to 24  for this x rotation 
        _xRotationsMap.insert( make_pair(*itr, ++currentLabel) );
    }
    for( itr = sensorIDsVec.begin(); itr != sensorIDsVec.end(); ++itr ) {// sensor 0 to 5 will have numbers 25 to 30  for this y rotation 
        _yRotationsMap.insert( make_pair(*itr, ++currentLabel) );
    }
    for( itr = sensorIDsVec.begin(); itr != sensorIDsVec.end(); ++itr ) {// sensor 0 to 5 will have numbers 31 to 36  for this z rotation 
        _zRotationsMap.insert( make_pair(*itr, ++currentLabel) );
    }
}
//This function calculates the alignment jacobain in the local frame of the telescope. Using the state parameters
int EUTelMillepede::computeAlignmentToMeasurementJacobian( EUTelState &state){
	if(_alignmentMode == Utility::noAlignment){
		throw(lcio::Exception(Utility::outputColourString("No alignment has been chosen.", "RED"))); 	
	}
	streamlog_out(DEBUG3) <<"State we arew about to add: "<< state.getLocation()<<endl; 
	state.print();
	float TxLocal =  state.getIntersectionLocalXZ();
	float TyLocal =  state.getIntersectionLocalYZ();
	float* localpos = state.getPosition();
	streamlog_out( DEBUG0 ) << "This is px/pz, py/pz (local) "<< TxLocal <<","<< TyLocal << std::endl;
	streamlog_out( DEBUG0 ) << "Local frame position "<< *localpos<<","<<*(localpos+1)<<","<<*(localpos+2) << std::endl;
	computeAlignmentToMeasurementJacobian( *localpos,*(localpos+1), TxLocal, TyLocal);
} 



//This function does the leg work of all the calculation for each state to determine the alignment jacobian
// _jacobian from below looks like //   (-1   0  -y   -dx/dz   x*dx/dz   -y*dx/dz)( x local )
																  //   (0  -1 	x   -dy/dz   x*dy/dz   -y*dy/dz )( y local )
                                  //                                             ( anti clockwise rotation around z) ?? Strange but whatever
                                  //                                             (moving the plane in the z direction)
																	//                                             (this is clockwise rotation look in + y direction )
																	// 	                                           (this is clockwise rotations in the x direction  )
int EUTelMillepede::computeAlignmentToMeasurementJacobian( float x,float y, float slopeXvsZ, float slopeYvsZ){
	_jacobian.Zero();
	streamlog_out(DEBUG0) << "This is the empty Alignment Jacobian" << std::endl;
	streamlog_message( DEBUG0, _jacobian.Print();, std::endl; );			
		
	//////////////////////////////////////Moving the sensor in x and y. Obviously if the sensor move right the hit will appear to move left. Always create this!!!! BEGIN
	_jacobian[0][0] = -1.0; // dxh/dxs      dxh => change in hit position         dxs => Change in sensor position
	_jacobian[0][1] = 0.0; // dxh/dys      //I have changed this from -1 to 1. It seems this is the sign convention that Millepede uses.
	_jacobian[1][0] = 0.0; // dyh/dxs
	_jacobian[1][1] = -1.0; // dyh/dys
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////END


	//////////////////////////////////////////////////////////////Rotations around the z axis of sensor. Only create this if you want rotations around z. BEGIN
	if (_alignmentMode == Utility::XYShiftXYRot
                || _alignmentMode == Utility::XYZShiftXYRot
                || _alignmentMode == Utility::XYShiftXZRotXYRot
                || _alignmentMode == Utility::XYShiftYZRotXYRot
                || _alignmentMode == Utility::XYShiftXZRotYZRotXYRot
                || _alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
		_jacobian[0][2] =  -y; // dxh/rotzs   rotzs => rotation of sensor around z axis
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

void EUTelMillepede::setGlobalLabels(EUTelState& state){
	setGlobalLabels( state.getLocation());
}
//Here depending on the palne the states is on we return a particular label for x shift, y shift.....
void EUTelMillepede::setGlobalLabels( int iPlane){
	//We alway add these to the first to places in the alignment matrix
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
		_globalLabels[3] = _zShiftsMap[iPlane]; // dz
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

			streamlog_out(DEBUG1) << "Output of global labels for plane "<<iPlane<<" The size of labels "<<_globalLabels.size() <<std::endl;
			for( std::vector<int>::const_iterator i = _globalLabels.begin(); i != _globalLabels.end(); ++i){
    		streamlog_out(DEBUG1) << *i << ' ';
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
	if(_alignmentMode == Utility::noAlignment){
		throw(lcio::Exception(Utility::outputColourString("No alignment has been chosen.", "RED"))); 	
	}	
	ofstream steerFile;
	steerFile.open(_milleSteeringFilename.c_str());//We open the text file se we can add text to it.
	if (!steerFile.is_open()) {
		throw(lcio::Exception(Utility::outputColourString("Could not open steering file.", "RED"))); 	
	}
	streamlog_out(DEBUG0) << "Millepede binary:" << _milleBinaryFilename << endl;
	steerFile << "Cfiles" << endl;
	steerFile << _milleBinaryFilename << endl;
	steerFile << endl;
	steerFile << "Parameter" << endl;
	//TO DO: There should be a test that all planes that are used have a state associated with them and that state has a hit
	for(int i =0 ; i < geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size(); ++i){
		int sensorId = geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(i); 
		///////////////////////////////////////////////////////////////////////////////////////////////////Determine if some of alignment parameters are fixed. BEGIN
		// check if plane has to be used as fixed
		// end points to one past the last element so this is ok.
		const bool isFixedXShift = std::find(_fixedAlignmentXShfitPlaneIds.begin(), _fixedAlignmentXShfitPlaneIds.end(), sensorId) != _fixedAlignmentXShfitPlaneIds.end();
		const bool isFixedYShift = std::find(_fixedAlignmentYShfitPlaneIds.begin(), _fixedAlignmentYShfitPlaneIds.end(), sensorId) != _fixedAlignmentYShfitPlaneIds.end();
		const bool isFixedZShift = std::find(_fixedAlignmentZShfitPlaneIds.begin(), _fixedAlignmentZShfitPlaneIds.end(), sensorId) != _fixedAlignmentZShfitPlaneIds.end();
		const bool isFixedXRotation = std::find(_fixedAlignmentXRotationPlaneIds.begin(), _fixedAlignmentXRotationPlaneIds.end(), sensorId) != _fixedAlignmentXRotationPlaneIds.end();
		const bool isFixedYRotation = std::find(_fixedAlignmentYRotationPlaneIds.begin(), _fixedAlignmentYRotationPlaneIds.end(), sensorId) != _fixedAlignmentYRotationPlaneIds.end();
		const bool isFixedZRotation = std::find(_fixedAlignmentZRotationPlaneIds.begin(), _fixedAlignmentZRotationPlaneIds.end(), sensorId) != _fixedAlignmentZRotationPlaneIds.end();
		////////////////////////////////////////////////////////////////////////////////////////////////////END
		//cout<<"(Z rotation)This is for sensor ID:  "<<sensorId<< " Found sensor ID using find   "<< *(std::find(_fixedAlignmentZRotationPlaneIds.begin(), _fixedAlignmentZRotationPlaneIds.end(), sensorId)) <<" Is it fixed? " << isFixedZRotation<<endl;
	//	streamlog_out(DEBUG0)<<"(Y shift) This is for sensor ID:  "<<sensorId<< " Found sensor ID using find   "<< *(std::find(_fixedAlignmentYShfitPlaneIds.begin(), _fixedAlignmentYShfitPlaneIds.end(), sensorId)) <<" Is it fixed? " << isFixedYShift<<endl;

		//TO DO: These uncertainties I believe come from the accuracy of the alignment jacobain. We currently just say this is 0.01. However it there a way to quantify this? 
		/////////////////////////////////////////////////////////////////////////////////////////////Now fill string that will go into steering depending on if fixed or not BEGIN
		const string initUncertaintyXShift = (isFixedXShift) ? "-1." : "1";//-1 means that this is fixed
		const string initUncertaintyYShift = (isFixedYShift) ? "-1." : "1";
		const string initUncertaintyZShift = (isFixedZShift) ? "-1." : "1";
		const string initUncertaintyXRotation = (isFixedXRotation) ? "-1." : "1";
		const string initUncertaintyYRotation = (isFixedYRotation) ? "-1." : "1";
		const string initUncertaintyZRotation = (isFixedZRotation) ? "-1." : "1";
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////END

		//TO DO: Determine is this initial shift is needed. 
		/*We can set the initial shift that millepede will work from. This would be the same as changing the gear file I think. I do not know why this is here.        
		const double initXshift = (isFixedXShift) ? 0. : _seedAlignmentConstants._xResiduals[sensorId]/_seedAlignmentConstants._nxResiduals[sensorId];
		const double initYshift = (isFixedYShift) ? 0. : _seedAlignmentConstants._yResiduals[sensorId]/_seedAlignmentConstants._nyResiduals[sensorId];
		*/          
		//Here we fill the steering file with:What planes are fixed,initial shifts,the uncertainties.
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

	} // end loop over all planes
	steerFile << endl;
	//Here we add some more paramter that millepede needs. This is involves: How is the solution found, How are outliers down weighted(These are hits that are very far from state hit) and chi2 cuts
	for ( StringVec::iterator it = pedeSteerAddCmds.begin( ); it != pedeSteerAddCmds.end( ); ++it ) {
		// two backslashes will be interpreted as newline
		if ( *it == "\\\\" ){
			steerFile << endl;
		}else{
			steerFile << *it << " ";
		}
	}
	steerFile << endl;
	steerFile << "end" << endl;
	steerFile.close();
}
//This function will use the mille binary file and steering and execute pede. Pede is the work horse of millepede. It does the actual minimisation procedure.
//It also write the results of this into a log file. This is very important since we need the information that this log file provides to determine what is the next step in out iterative alignment
//By this I mean if too many tracks were rejected by millepede then on the next iteration we need to increase increase the chi2 cut and increase the hit residual.
int EUTelMillepede::runPede(){
	std::string command = "pede " + _milleSteeringFilename;//This is just the same as running a command line command pede <steering file> the minimisation would still be done.
	streamlog_out ( MESSAGE5 ) << "Starting pede...: " << command.c_str( ) << endl;
  redi::ipstream pede( command.c_str( ), redi::pstreams::pstdout | redi::pstreams::pstderr );// run pede and create a streambuf that reads its stdout and stderr

	if ( !pede.is_open( ) ) {
		throw(lcio::Exception(Utility::outputColourString("The pede file could not be openned. ", "RED")));
  } else {
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
					streamlog_out( ERROR5 ) << error;
                    //encounteredError = true;
        }
				if ( pede.eof( ) ) {
					finished[0] = true;
					if ( !finished[1] ) pede.clear( );
				}
    	}

			if ( !finished[1] ) {
				while ( ( n = pede.out( ).readsome( buf, sizeof (buf ) ) ) > 0 ) {
					streamlog_out( MESSAGE0 ).write( buf, n ).flush( );
					string output ( buf, n );
					pedeoutput << output;
					streamlog_out( MESSAGE9 )  << output;
				}
				if ( pede.eof( ) ) {
					finished[1] = true;
					if ( !finished[0] )
					pede.clear( );
				}
			}
		}
		pede.close( );
		//TO DO: Surely we can just specify the directory that we want this placed in. Need to check     
		//  if ( parseMilleOutput( "millepede.res" ) ) //moveMilleResultFile( "millepede.res", _milleResultFileName );
	}//END OF IF STATEMENT
return 0;
}
//This part using the output of millepede will create a new gear file based on the alignment parameters that have just been determined
//It will also create LCIO file that will hold the alignment constants
bool EUTelMillepede::parseMilleOutput(std::string alignmentConstantLCIOFile, std::string gear_aligned_file){
	ifstream file( _milleResultFileName.c_str() );
	if ( !file.good( ) ) {
		throw(lcio::Exception(Utility::outputColourString("Can not open millepede results file. ", "RED")));
	}
	const string command = "parsemilleout.sh " + _milleSteeringFilename + " " + _milleResultFileName + " " + alignmentConstantLCIOFile + 
												 " " + Global::parameters->getStringVal("GearXMLFile" ) + " " + gear_aligned_file;
	streamlog_out ( MESSAGE5 ) << "Converting millepede results to LCIO collections... " << endl;
	streamlog_out ( MESSAGE5 ) << command << endl;
	// run pede and create a streambuf that reads its stdout and stderr
	redi::ipstream parsepede( command.c_str( ), redi::pstreams::pstdout | redi::pstreams::pstderr );
	if ( !parsepede.is_open( )){
		throw(lcio::Exception(Utility::outputColourString("Could not open the parsepede file. ", "RED")));
	}else{
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
				}
				if ( parsepede.eof( ) ) {
					finished[0] = true;
					if ( !finished[1] )	parsepede.clear( );
				}
			}

			if ( !finished[1] ) {
				while ( ( n = parsepede.out( ).readsome( buf, sizeof (buf ) ) ) > 0 ) {
					streamlog_out( MESSAGE9 ).write( buf, n ).flush( );
					string output ( buf, n );
					parsepedeoutput << output;
				}
				if ( parsepede.eof( ) ) {
					finished[1] = true;
					if ( !finished[0] ) parsepede.clear( );
				}
			}
		}
		parsepede.close( );
	}
	return true;
}

void EUTelMillepede::CreateBinary(){
        streamlog_out(DEBUG0) << "Initialising Mille..." << std::endl;
				streamlog_out(DEBUG0) << "Millepede binary:" << _milleBinaryFilename << endl;

        const unsigned int reserveSize = 80000;
				std::string string = "millepede.bin"; //TO DO:need to fix this. Not reading it correctly
        _milleGBL = new gbl::MilleBinary(string, reserveSize);

        if (_milleGBL == NULL) {
            streamlog_out(ERROR) << "Can't allocate an instance of mMilleBinary. Stopping ..." << std::endl;
            throw UnknownDataTypeException("MilleBinary was not created");
        }
}

void EUTelMillepede::testUserInput(){
	bool fixedGood=true;
	if(_fixedAlignmentXShfitPlaneIds.size()== 0){
		streamlog_out(MESSAGE9) <<	Utility::outputColourString("You have not fixed any X shifts. This is ill advised .", "YELLOW")<<std::endl;
		fixedGood=false;
	}
	if(_fixedAlignmentYShfitPlaneIds.size()== 0){
		streamlog_out(MESSAGE9) <<	Utility::outputColourString("You have not fixed any Y shifts. This is ill advised .", "YELLOW")<<std::endl;
		fixedGood=false;
	}
	if(_fixedAlignmentZShfitPlaneIds.size()== 0){
		streamlog_out(MESSAGE9) <<	Utility::outputColourString("You have not fixed any Z shifts. This is ill advised .", "YELLOW")<<std::endl;
		fixedGood=false;
	}
	if(_fixedAlignmentXRotationPlaneIds.size()== 0){
		streamlog_out(MESSAGE9) <<	Utility::outputColourString("You have not fixed any X rotations. This is ill advised .", "YELLOW")<<std::endl;
		fixedGood=false;
	}
	if(_fixedAlignmentYRotationPlaneIds.size()== 0){
		streamlog_out(MESSAGE9) <<	Utility::outputColourString("You have not fixed any Y rotations. This is ill advised .", "YELLOW")<<std::endl;
		fixedGood=false;
	}
	if(_fixedAlignmentZRotationPlaneIds.size()== 0){
		streamlog_out(MESSAGE9) <<	Utility::outputColourString("You have not fixed any Z rotations. This is ill advised .", "YELLOW")<<std::endl;
		fixedGood=false;
	}
	if(fixedGood){
		streamlog_out(MESSAGE9) <<	Utility::outputColourString("For all alignment parameters each has one fixed. GOOD! .", "GREEN")<<std::endl;
	}
}
void EUTelMillepede::printFixedPlanes(){
	streamlog_out(MESSAGE5)<<"These are the planes what have X shifts fixed"<<endl;
	for(int i=0;i<_fixedAlignmentXShfitPlaneIds.size();++i){
		streamlog_out(MESSAGE5)<<_fixedAlignmentXShfitPlaneIds.at(i)<<"  ";
	}
	streamlog_out(MESSAGE5)<<endl<<"These are the planes what have Y shifts fixed"<<endl;
	for(int i=0;i<_fixedAlignmentYShfitPlaneIds.size();++i){
		streamlog_out(MESSAGE5)<<_fixedAlignmentYShfitPlaneIds.at(i)<<"  ";
	}
	streamlog_out(MESSAGE5)<<endl<<"These are the planes what have Z shifts fixed"<<endl;
	for(int i=0;i<_fixedAlignmentZShfitPlaneIds.size();++i){
		streamlog_out(MESSAGE5)<<_fixedAlignmentZShfitPlaneIds.at(i)<<"  ";
	}
	streamlog_out(MESSAGE5)<<endl<<"These are the planes what have X Rotation fixed"<<endl;
	for(int i=0;i<_fixedAlignmentXRotationPlaneIds.size();++i){
		streamlog_out(MESSAGE5)<<_fixedAlignmentXRotationPlaneIds.at(i)<<"  ";
	}
	streamlog_out(MESSAGE5)<<endl<<"These are the planes what have Y Rotation fixed"<<endl;
	for(int i=0;i<_fixedAlignmentYRotationPlaneIds.size();++i){
		streamlog_out(MESSAGE5)<<_fixedAlignmentYRotationPlaneIds.at(i)<<"  ";
	}
	streamlog_out(MESSAGE5)<<endl<<"These are the planes what have Z Rotation fixed"<<endl;
	for(int i=0;i<_fixedAlignmentZRotationPlaneIds.size();++i){
		streamlog_out(MESSAGE5)<<_fixedAlignmentZRotationPlaneIds.at(i)<<"  ";
	}
	streamlog_out(MESSAGE5)<<endl<<"The planes we will align with are: "<<endl;
	for(int i =0 ; i < geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size(); ++i){
		streamlog_out(MESSAGE5)<<geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(i)<<"  ";
	}
	streamlog_out(MESSAGE5)<<endl;
}
} // namespace eutelescope






