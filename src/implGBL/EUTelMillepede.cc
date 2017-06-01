#include "EUTelMillepede.h"

using namespace lcio;
using namespace std;
using namespace marlin;
using namespace eutelescope;


namespace eutelescope {

//This constructor useful for mille binary output part
EUTelMillepede::EUTelMillepede(std::string& res, std::string& newGear) :
_milleGBL(NULL),
_jacobian(2,6),
_globalLabels(6),
_milleSteeringFilename("steer.txt"),
_milleSteerNameOldFormat("steer-iteration-0.txt"),
_iteration(1),
_newGear(newGear),
_milleResultFileName(res)   
    
{
	FillMilleParametersLabels();
	CreateBinary();
}

EUTelMillepede::~EUTelMillepede(){}

//Note here we label sensors and every alignment degree of freedom uniquely. Note that even if the sensor is to remain fixed. The fixing is done latter.
void EUTelMillepede::FillMilleParametersLabels() {

    const std::vector<int> sensorIDsVec = EUTelExcludedPlanes::_senNoDeadMaterial;
    size_t noOfSensors = sensorIDsVec.size(); 

    for( IntVec::const_iterator itr = sensorIDsVec.begin(); itr != sensorIDsVec.end(); ++itr ) {//sensor 0 to 5 will have numbers 1 to 6 for this x shift
        _xShiftsMap.insert( make_pair(*itr, (*itr)*10 + 1) );
        _yShiftsMap.insert( make_pair(*itr, (*itr)*10 + 2) );
        _zShiftsMap.insert( make_pair(*itr, (*itr)*10 + 3) );
        _xRotationsMap.insert( make_pair(*itr, (*itr)*10 + 4) );
        _yRotationsMap.insert( make_pair(*itr, (*itr)*10 + 5) );
        _zRotationsMap.insert( make_pair(*itr, (*itr)*10 + 6) );
    }
}
//This function calculates the alignment jacobain in the local frame of the telescope. Using the state parameters
void EUTelMillepede::computeAlignmentToMeasurementJacobian( EUTelState &state){
	streamlog_out(DEBUG3) <<"State we are about to add: "<< state.getLocation()<<endl; 
	state.print();
	float TxLocal =  state.getSlopeX();
	float TyLocal =  state.getSlopeY();
	const double* localpos = state.getPosition();
	streamlog_out( DEBUG0 ) << "This is px/pz, py/pz (local) "<< TxLocal <<","<< TyLocal << std::endl;
	streamlog_out( DEBUG0 ) << "This is px/pz, py/pz (Global) "<< state.getDirGlobal()[0]/state.getDirGlobal()[2]<<","<< state.getDirGlobal()[1]/state.getDirGlobal()[2]<< std::endl;
	streamlog_out( DEBUG0 ) << "Local frame position "<< *localpos<<","<<*(localpos+1)<<","<<*(localpos+2) << std::endl;
	computeAlignmentToMeasurementJacobian( *localpos,*(localpos+1), TxLocal, TyLocal);
} 
void EUTelMillepede::computeAlignmentGlobal( EUTelState &state){
	streamlog_out(DEBUG3) <<"State we are about to add: "<< state.getLocation()<<endl; 
	state.print();
    /// Everything direction related defined in the global frame.
    /// The position of the hits are before the offsets are taken into account. 
    /// This allows the correct magnitude of rotation to be determined.
    Eigen::Vector3d normal  = geo::gGeometry().siPlaneNormalEig(state.getLocation());
    Eigen::Vector3d offset =  geo::gGeometry().getOffsetVector(state.getLocation());
	streamlog_out( DEBUG0 ) << "Offset: "<< offset[0] <<  " " << offset[1] << " " <<offset[2] << std::endl;

    Eigen::Vector3d dir    = state.getDirGlobalEig();
    ///Deduct offset to remove set from gear.
    double relX = state.getPositionGlobal()[0] - offset[0];
    double relY = state.getPositionGlobal()[1] - offset[1];
    double relZ = state.getPositionGlobal()[2] - offset[2];
    Eigen::MatrixXd I(3,3);
    I.setIdentity();
    double factor =  dir.transpose()*normal;
    Eigen::MatrixXd drldm = I - (dir*normal.transpose())*(1.0/factor); //Residual change with change in global position 

	streamlog_out( DEBUG0 ) << "drldm: "<<endl <<  drldm << std::endl;

    //Change in global position with change in alignment parameters.  
    Eigen::MatrixXd dmdg(3,6);
    dmdg << 1 , 0,  0, 0,   relZ, -relY,
            0 ,  1,  0,-relZ,  0,    relX,
           0 ,   0,  1, relY, -relX, 0;

    Eigen::Matrix3d rot  =  geo::gGeometry().getRotMatrixEig(state.getLocation());
    Eigen::Matrix3d rotInv = rot.transpose();
    Eigen::MatrixXd gloMes = rotInv*drldm;
    Eigen::MatrixXd aliToLoc =  gloMes*dmdg;
	streamlog_out( DEBUG0 ) << "1)Align to global shifts: "<< endl<< dmdg << std::endl;
	streamlog_out( DEBUG0 ) << "Relates motion of the global position to position on plane: "<< endl << drldm << std::endl;
	streamlog_out( DEBUG0 ) << "2)Rotate plane to work in local coordinates: "<< endl << gloMes << std::endl;
	streamlog_out( DEBUG0 ) << "Times (1)/(2) together and we get the relation between alignment and local residual "<< endl << aliToLoc << std::endl;

	_jacobian.Zero();
    for(int i = 0 ; i < 6; ++i){
        _jacobian[0][i] = aliToLoc(0,i); 
        _jacobian[1][i] = aliToLoc(1,i); 
    }

// dxh/dxs
// dyh/dxs
// dxh/dxs      
// dyh/dxs
// dxh/dys     
// dyh/dys
// dxh/rotzs   
// dyh/rotzs
// dxh/dzs
// dyh/dzs
// dxh/rotyr
// dyh/rotyr
// dxh/rotxr          
// dyh/rotxr         


} 


//This function does the leg work of all the calculation for each state to determine the alignment jacobian
// _jacobian from below looks like //   (-1   0  -y   -dx/dz   x*dx/dz   -y*dx/dz)( x local )
																  //   (0  -1 	x   -dy/dz   x*dy/dz   -y*dy/dz )( y local )
                                  //                                             ( anti clockwise rotation around z) ?? Strange but whatever
                                  //                                             (moving the plane in the z direction)
																	//                                             (this is clockwise rotation look in + y direction )
																	// 	                                           (this is clockwise rotations in the x direction  )
void EUTelMillepede::computeAlignmentToMeasurementJacobian( float x,float y, float slopeXvsZ, float slopeYvsZ){
	_jacobian.Zero();
	streamlog_out(DEBUG0) << "This is the empty Alignment Jacobian" << std::endl;
	streamlog_message( DEBUG0, _jacobian.Print();, std::endl; );			
		
	_jacobian[0][0] = -1.0; // dxh/dxs      dxh => change in hit position         dxs => Change in sensor position
	_jacobian[1][0] = 0.0; // dyh/dxs
	_jacobian[0][1] = 0.0; // dxh/dys     
	_jacobian[1][1] = -1.0; // dyh/dys
	_jacobian[0][2] = y; // dxh/rotzs   
	_jacobian[1][2] = -x; // dyh/rotzs
	_jacobian[0][3] = slopeXvsZ; // dxh/dzs
	_jacobian[1][3] = slopeYvsZ; // dyh/dzs
	_jacobian[0][4] = -x*slopeXvsZ; // dxh/rotyr
	_jacobian[1][4] = -x*slopeYvsZ; // dyh/rotyr
	_jacobian[0][5] = -y*slopeXvsZ; // dxh/rotxr          
	_jacobian[1][5] = -y*slopeYvsZ; // dyh/rotxr         
}

void EUTelMillepede::setGlobalLabels(EUTelState& state){
	setGlobalLabels( state.getLocation());
}
//Here depending on the palne the states is on we return a particular label for x shift, y shift.....
void EUTelMillepede::setGlobalLabels( int iPlane){
	//We alway add these to the first to places in the alignment matrix
//	_globalLabels[0] = _xShiftsMap[iPlane]; // dx
//	_globalLabels[1] = _yShiftsMap[iPlane]; // dy
//  	_globalLabels[2] = _zRotationsMap[iPlane]; // rot z
//    _globalLabels[3] = _zShiftsMap[iPlane]; // dz
//    _globalLabels[4] = _yRotationsMap[iPlane]; // drot y         
//    _globalLabels[5] = _xRotationsMap[iPlane]; // drot x  
//
    _globalLabels[0] = _xShiftsMap[iPlane]; 
    _globalLabels[1] = _yShiftsMap[iPlane]; 
    _globalLabels[2] = _zShiftsMap[iPlane];    
    _globalLabels[3] = _xRotationsMap[iPlane];     
    _globalLabels[4] = _yRotationsMap[iPlane];   
    _globalLabels[5] = _zRotationsMap[iPlane]; 

    streamlog_out(DEBUG1) << "Output of global labels for plane "<<iPlane<<" The size of labels "<<_globalLabels.size() <<std::endl;
    for( std::vector<int>::const_iterator i = _globalLabels.begin(); i != _globalLabels.end(); ++i){
    streamlog_out(DEBUG1) << *i << ' ';
    }
    streamlog_out(DEBUG1) << std::endl;

}
  
void EUTelMillepede::writeMilleSteeringFile(lcio::StringVec pedeSteerAddCmds){
    streamlog_out(DEBUG2) << "EUTelMillepede::writeMilleSteeringFile------------------------------------BEGIN" << std::endl;
    ofstream steerFile;
    steerFile.open(_milleSteeringFilename.c_str());//We open the text file se we can add text to it.
    if (!steerFile.is_open()) {
        throw(lcio::Exception("Could not open steering file."));    
    }
    streamlog_out(DEBUG0) << "Millepede binary:" << _milleBinaryFilename << std::endl;
    steerFile << "Cfiles" << std::endl;
    steerFile << _milleBinaryFilename << std::endl;
    steerFile << std::endl;
    steerFile << "Parameter" << std::endl;
    for(size_t i =0 ; i < EUTelExcludedPlanes::_senNoDeadMaterial.size(); ++i){
        int sensorId = EUTelExcludedPlanes::_senNoDeadMaterial.at(i); 
        const bool isFixedXShift = std::find(_fixedAlignmentXShfitPlaneIds.begin(), _fixedAlignmentXShfitPlaneIds.end(), sensorId) != _fixedAlignmentXShfitPlaneIds.end();
        const bool isFixedYShift = std::find(_fixedAlignmentYShfitPlaneIds.begin(), _fixedAlignmentYShfitPlaneIds.end(), sensorId) != _fixedAlignmentYShfitPlaneIds.end();
        const bool isFixedZShift = std::find(_fixedAlignmentZShfitPlaneIds.begin(), _fixedAlignmentZShfitPlaneIds.end(), sensorId) != _fixedAlignmentZShfitPlaneIds.end();
        const bool isFixedXRotation = std::find(_fixedAlignmentXRotationPlaneIds.begin(), _fixedAlignmentXRotationPlaneIds.end(), sensorId) != _fixedAlignmentXRotationPlaneIds.end();
        const bool isFixedYRotation = std::find(_fixedAlignmentYRotationPlaneIds.begin(), _fixedAlignmentYRotationPlaneIds.end(), sensorId) != _fixedAlignmentYRotationPlaneIds.end();
        const bool isFixedZRotation = std::find(_fixedAlignmentZRotationPlaneIds.begin(), _fixedAlignmentZRotationPlaneIds.end(), sensorId) != _fixedAlignmentZRotationPlaneIds.end();
        const string initUncertaintyXShift = (isFixedXShift) ? "-1." : "1";//-1 means that this is fixed
        const string initUncertaintyYShift = (isFixedYShift) ? "-1." : "1";
        const string initUncertaintyZShift = (isFixedZShift) ? "-1." : "1";
        const string initUncertaintyXRotation = (isFixedXRotation) ? "-1." : "1";
        const string initUncertaintyYRotation = (isFixedYRotation) ? "-1." : "1";
        const string initUncertaintyZRotation = (isFixedZRotation) ? "-1." : "1";

        const double initXshift =0; const double initYshift = 0;
        steerFile << left << setw(25) << _xShiftsMap[sensorId] << setw(25) << -initXshift << setw(25) << initUncertaintyXShift
        << setw(25) << " ! X shift " << sensorId << std::endl;
        steerFile << left << setw(25) << _yShiftsMap[sensorId] << setw(25) << -initYshift << setw(25) << initUncertaintyYShift
        << setw(25) << " ! Y shift " << sensorId << std::endl;
        steerFile << left << setw(25) << _zShiftsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyZShift
        << setw(25) << " ! Z shift " << sensorId << std::endl;
        steerFile << left << setw(25) << _yRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyYRotation
        << setw(25) << " ! XZ rotation " << sensorId << std::endl;
        steerFile << left << setw(25) << _xRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyXRotation
        << setw(25) << " ! YZ rotation " << sensorId << std::endl;
        steerFile << left << setw(25) << _zRotationsMap[sensorId] << setw(25) << "0.0" << setw(25) << initUncertaintyZRotation
        << setw(25)  << " ! XY rotation " << sensorId << std::endl;

    } 
    steerFile << std::endl;
    //Here we add some more paramter that millepede needs. This involves: How is the solution found, How are outliers down weighted(These are hits that are very far from state hit) and chi2 cuts
 //   for ( StringVec::iterator it = pedeSteerAddCmds.begin( ); it != pedeSteerAddCmds.end( ); ++it ) {
        // two backslashes will be interpreted as newline
//        if ( *it == "\\\\" ){
//            steerFile << std::endl;
 //       }else{
  //          steerFile << *it << " ";
   //     }
 //   }
    steerFile <<  "method inversion 5 0.1" <<std::endl;
    steerFile <<  "outlierdownweighting 4" << std::endl; 
    steerFile << "chiscut 30 6" <<std::endl; 
    steerFile << "end" << std::endl;
    steerFile.close();
}
void EUTelMillepede::copyFile(std::string _milleSteeringFilename, std::string _milleSteerNameOldFormat){
	std::ifstream infile (_milleSteeringFilename.c_str(),std::ifstream::binary);
	std::ofstream outfile (_milleSteerNameOldFormat.c_str() ,std::ofstream::binary);

	// get size of file
	infile.seekg (0,infile.end);
	long size = infile.tellg();
	infile.seekg (0);

	// allocate memory for file content
	char* buffer = new char[size];

	// read content of infile
	infile.read (buffer,size);

	// write to outfile
	outfile.write (buffer,size);

	// release dynamically-allocated memory
	delete[] buffer;

	outfile.close();
	infile.close();
}

//This function will use the mille binary file and steering and execute pede. Pede is the work horse of millepede. It does the actual minimisation procedure.
//It also write the results of this into a log file. This is very important since we need the information that this log file provides to determine what is the next step in out iterative alignment
//By this I mean if too many tracks were rejected by millepede then on the next iteration we need to increase increase the chi2 cut and increase the hit residual.
bool EUTelMillepede::runPede(){
	std::string command = "pede " + _milleSteeringFilename;//This is just the same as running a command line command pede <steering file> the minimisation would still be done.
	streamlog_out ( MESSAGE5 ) << "Starting pede...: " << command.c_str( ) << std::endl;
  redi::ipstream pede( command.c_str( ), redi::pstreams::pstdout | redi::pstreams::pstderr );// run pede and create a streambuf that reads its stdout and stderr

	if ( !pede.is_open( ) ) {
		throw(lcio::Exception("The pede file could not be openned."));
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
		std::string output = pedeoutput.str();
		bool found = findTooManyRejects(output);
		//TO DO: Surely we can just specify the directory that we want this placed in. Need to check     
		//  if ( parseMilleOutput( "millepede.res" ) ) //moveMilleResultFile( "millepede.res", _milleResultFileName );
		return found;
	}//END OF IF STATEMENT
}
bool EUTelMillepede::findTooManyRejects(std::string output){
	int found = output.find("Too many rejects (>33.3%)");
	if (found == std::string::npos){
		streamlog_out(MESSAGE5)<<endl<<"Number of rejects low. Continue with alignment."<<endl;
		return false;

	}else{
		streamlog_out(MESSAGE5)<<endl<<"Number of rejects high. We can't use this binary for alignment"<<endl;
		return true;
	}
}
void EUTelMillepede::editSteerUsingRes(){
	ifstream file( _milleResultFileName.c_str() );
	if ( !file.good( ) ) {
		throw(lcio::Exception("Can not open millepede results file. In editSteerUsingRes()"));
	}
	const string command = "resIntoSteer.py " + _milleSteeringFilename + " " + _milleResultFileName;
	streamlog_out ( MESSAGE5 ) << "Results fill used to create new steering file: " << std::endl;
	streamlog_out ( MESSAGE5 ) << command << std::endl;
	// run pede and create a streambuf that reads its stdout and stderr
	redi::ipstream parsepede( command.c_str( ), redi::pstreams::pstdout | redi::pstreams::pstderr );


}

bool EUTelMillepede::converge(){
    bool converged;
    bool rejectsHigh;
    for(unsigned int i=0; i<2 ; i++){
		editSteerUsingRes();
		converged = checkConverged();//Will simply output the steering files used in each iteration. 
		rejectsHigh = runPede();	
	}
	streamlog_out ( MESSAGE5 ) << "The number if rejects pass/fail: " << rejectsHigh  << std::endl;
	//We do this to create a new steering file from all the iterations but do not run pede.
	editSteerUsingRes();
	converged = checkConverged();//Will simply output the steering files used in each iteration. 
	return converged;
}

bool EUTelMillepede::checkConverged(){
	outputSteeringFiles();
    return true;
}
void EUTelMillepede::outputSteeringFiles(){
	std::stringstream output;
	output << "steer-iteration-" << _iteration <<".txt";
	std:: string outputName = output.str();
	copyFile(_milleSteeringFilename, outputName);
	_iteration++;

}
void EUTelMillepede::getNewGear(){
    ///Must create string explicitly on right had side. Since c++ does not have + char operator.
    const string command = std::string("gearUpdate ") + " -og " + Global::parameters->getStringVal("GearXMLFile" )  + " -r " + _milleResultFileName + " -ng " +  _newGear;
	streamlog_out ( MESSAGE5 ) << command << std::endl;
	redi::ipstream updateGear( command.c_str( ), redi::pstreams::pstdout | redi::pstreams::pstderr );

	if ( !updateGear.is_open( )){
		throw(lcio::Exception("Could not open the updateGear file. "));
	}else{
	// output multiplexing: parse updateGear output in both stdout and stderr and echo messages accordingly
	char buf[1024];
		std::streamsize n;
		std::stringstream updateGearoutput; // store stdout to parse later
		std::stringstream updateGearerrors;
		bool finished[2] = { false, false };
		while ( !finished[0] || !finished[1] ) {
			if ( !finished[0] ) {
				while ( ( n = updateGear.err( ).readsome( buf, sizeof (buf ) ) ) > 0 ) {
					streamlog_out( ERROR5 ).write( buf, n ).flush( );
					string error ( buf, n );
					updateGearerrors << error;
				}
				if ( updateGear.eof( ) ) {
					finished[0] = true;
					if ( !finished[1] )	updateGear.clear( );
				}
			}

			if ( !finished[1] ) {
				while ( ( n = updateGear.out( ).readsome( buf, sizeof (buf ) ) ) > 0 ) {
					streamlog_out( MESSAGE9 ).write( buf, n ).flush( );
					string output ( buf, n );
					updateGearoutput << output;
				}
				if ( updateGear.eof( ) ) {
					finished[1] = true;
					if ( !finished[0] ) updateGear.clear( );
				}
			}
		}
		updateGear.close( );
	}


};

//This part using the output of millepede will create a new gear file based on the alignment parameters that have just been determined
//It will also create LCIO file that will hold the alignment constants
bool EUTelMillepede::parseMilleOutput(std::string alignmentConstantLCIOFile, std::string gear_aligned_file){
	ifstream file( _milleResultFileName.c_str() );
	if ( !file.good( ) ) {
		throw(lcio::Exception("Can not open millepede results file. parseMilleOutput()"));
	}
	ifstream file2( _milleSteerNameOldFormat.c_str() );
	if ( !file2.good( ) ) {
		throw(lcio::Exception("Can not open millepede old steer file. parseMilleOutput()"));
	}

	const string command = "parsemilleout.sh " + _milleSteerNameOldFormat + " " + _milleResultFileName + " " + alignmentConstantLCIOFile + 
												 " " + Global::parameters->getStringVal("GearXMLFile" ) + " " + gear_aligned_file;
	streamlog_out ( MESSAGE5 ) << "Converting millepede results to LCIO collections... " << std::endl;
	streamlog_out ( MESSAGE5 ) << command << std::endl;
	// run pede and create a streambuf that reads its stdout and stderr
	redi::ipstream parsepede( command.c_str( ), redi::pstreams::pstdout | redi::pstreams::pstderr );
	if ( !parsepede.is_open( )){
		throw(lcio::Exception("Could not open the parsepede file. "));
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
				streamlog_out(DEBUG0) << "Millepede binary:" << _milleBinaryFilename << std::endl;

        const unsigned int reserveSize = 0;//This is the number of elements the vector will have as start for alignment parameters and derivatives.
				//Can still push more onto the vector.
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
		streamlog_out(MESSAGE9) <<	"You have not fixed any X shifts. This is ill advised ." << std::endl;
		fixedGood=false;
	}
	if(_fixedAlignmentYShfitPlaneIds.size()== 0){
		streamlog_out(MESSAGE9) <<	"You have not fixed any Y shifts. This is ill advised." << std::endl;
		fixedGood=false;
	}
	if(_fixedAlignmentZShfitPlaneIds.size()== 0){
		streamlog_out(MESSAGE9) <<	"You have not fixed any Z shifts. This is ill advised ." << std::endl; 
		fixedGood=false;
	}
	if(_fixedAlignmentXRotationPlaneIds.size()== 0){
		streamlog_out(MESSAGE9) <<	"You have not fixed any X rotations. This is ill advised ." << std::endl;
		fixedGood=false;
	}
	if(_fixedAlignmentYRotationPlaneIds.size()== 0){
		streamlog_out(MESSAGE9) <<	"You have not fixed any Y rotations. This is ill advised ."<< std::endl;
		fixedGood=false;
	}
	if(_fixedAlignmentZRotationPlaneIds.size()== 0){
		streamlog_out(MESSAGE9) <<	"You have not fixed any Z rotations. This is ill advised ." << std::endl;
		fixedGood=false;
	}
	if(fixedGood){
		streamlog_out(MESSAGE9) <<	"For all alignment parameters each has one fixed. GOOD! ." << std::endl;
	}
}
void EUTelMillepede::printFixedPlanes(){
	streamlog_out(MESSAGE5)<<"These are the planes what have X shifts fixed"<<endl;
	for(size_t i=0;i<_fixedAlignmentXShfitPlaneIds.size();++i){
		streamlog_out(MESSAGE5)<<_fixedAlignmentXShfitPlaneIds.at(i)<<"  ";
	}
	streamlog_out(MESSAGE5)<<endl<<"These are the planes what have Y shifts fixed"<<endl;
	for(size_t i=0;i<_fixedAlignmentYShfitPlaneIds.size();++i){
		streamlog_out(MESSAGE5)<<_fixedAlignmentYShfitPlaneIds.at(i)<<"  ";
	}
	streamlog_out(MESSAGE5)<<endl<<"These are the planes what have Z shifts fixed"<<endl;
	for(size_t i=0;i<_fixedAlignmentZShfitPlaneIds.size();++i){
		streamlog_out(MESSAGE5)<<_fixedAlignmentZShfitPlaneIds.at(i)<<"  ";
	}
	streamlog_out(MESSAGE5)<<endl<<"These are the planes what have X Rotation fixed"<<endl;
	for(size_t i=0;i<_fixedAlignmentXRotationPlaneIds.size();++i){
		streamlog_out(MESSAGE5)<<_fixedAlignmentXRotationPlaneIds.at(i)<<"  ";
	}
	streamlog_out(MESSAGE5)<<endl<<"These are the planes what have Y Rotation fixed"<<endl;
	for(size_t i=0;i<_fixedAlignmentYRotationPlaneIds.size();++i){
		streamlog_out(MESSAGE5)<<_fixedAlignmentYRotationPlaneIds.at(i)<<"  ";
	}
	streamlog_out(MESSAGE5)<<endl<<"These are the planes what have Z Rotation fixed"<<endl;
	for(size_t i=0;i<_fixedAlignmentZRotationPlaneIds.size();++i){
		streamlog_out(MESSAGE5)<<_fixedAlignmentZRotationPlaneIds.at(i)<<"  ";
	}
//	streamlog_out(MESSAGE5)<<endl<<"The planes we will align with are: "<<endl;
//	for(size_t i =0 ; i < geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size(); ++i){
//		streamlog_out(MESSAGE5)<<geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(i)<<"  ";
//	}
//	streamlog_out(MESSAGE5)<<endl;
}
} // namespace eutelescope






