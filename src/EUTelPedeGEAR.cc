/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelPedeGEAR.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"
#include "EUTelPStream.h"
//#include "EUTelCDashMeasurement.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

// marlin includes ".h"
#include "marlin/Global.h"
#include "marlin/StringParameters.h"

// lcio includes <.h>
#include <Exceptions.h>

// system includes <>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace marlin;
using namespace eutelescope;

EUTelPedeGEAR::EUTelPedeGEAR() : Processor("EUTelPedeGEAR") {

  // modify processor description
  _description = "EUTelPedeGEAR calls PEDE to process a MILLE binary file and create an updated GEAR file with the updated MILLEPEDE II alignment constants.";

  registerOptionalParameter("ExcludePlanes","Exclude planes from fit according to their sensor ids.",_excludePlanes_sensorIDs ,IntVec());

  registerOptionalParameter("FixedPlanes","Fix sensor planes in the fit according to their sensor ids.",_FixedPlanes_sensorIDs ,IntVec());

  registerOptionalParameter("AlignMode","Number of alignment constants used. Available mode are: "
                            "\n1 - shifts in the X and Y directions and a rotation around the Z axis,"
                            "\n2 - only shifts in the X and Y directions"
                            "\n3 - (EXPERIMENTAL) shifts in the X,Y and Z directions and rotations around all three axis",
                            _alignMode, static_cast<int>(1));

  registerOptionalParameter("PedeSteerfileName","Name of the steering file for the pede program.",_pedeSteerfileName, std::string("steer_mille.txt"));
  
  registerOptionalParameter("NewGEARSuffix", "Suffix for the new GEAR file, set to empty string (this is not default!) to overwrite old GEAR file", _GEARFileSuffix, std::string("_aligned") );
}

void EUTelPedeGEAR::init() {
	//this method is called only once even when the rewind is active usually a good idea to
	printParameters ();

	// set to zero the run and event counters
	_iRun = 0;
	_iEvt = 0;

	//Getting access to geometry description
	geo::gGeometry().initializeTGeoDescription(EUTELESCOPE::GEOFILENAME, EUTELESCOPE::DUMPGEOROOT);

	//check if the GEAR manager pointer is not null!
	if( Global::GEAR == 0x0 ) {
		streamlog_out( ERROR2 ) << "The GearMgr is not available, for an unknown reason." << std::endl;
		throw InvalidGeometryException("GEAR manager is not initialised");
	}

	//the number of planes is got from the GEAR description and is
	//the sum of the telescope reference planes and the DUT (if any)
	_nPlanes = geo::gGeometry().nPlanes();
  
	//an associative std::map for getting also the sensorID ordered
	std::map<double, int> sensorIDMap;
	
	std::vector<int> sensorIDVec = geo::gGeometry().sensorIDsVec();
	//lets create an array with the z positions of each layer
	
	for( std::vector<int>::iterator it = sensorIDVec.begin(); it != sensorIDVec.end(); it++) {
		int sensorID = *it;
		sensorIDMap.insert( std::make_pair( geo::gGeometry().siPlaneZPosition(sensorID), sensorID ) );
	}
	
	//the user is giving sensor ids for the planes to be excluded. this
	//sensor ids have to be converted to a local index according to the
	//planes positions along the z axis.
	for(size_t i = 0; i < _FixedPlanes_sensorIDs.size(); i++) {
		std::map< double, int >::iterator iter = sensorIDMap.begin();
		int counter = 0;
		while( iter != sensorIDMap.end() ) {
			if( iter->second == _FixedPlanes_sensorIDs[i]) {
				_FixedPlanes.push_back(counter);
				break;
			}
			++iter;
			++counter;
		}
	}

	for(size_t i = 0; i < _excludePlanes_sensorIDs.size(); i++) {
		std::map< double, int >::iterator iter = sensorIDMap.begin();
		int counter = 0;
		while( iter != sensorIDMap.end() ) {
			if( iter->second == _excludePlanes_sensorIDs[i]) {
				_excludePlanes.push_back(counter);
				break;
			}
			++iter;
			++counter;
		}
	}

	//strip from the std::map the sensor id already sorted.
	std::map< double, int >::iterator iter = sensorIDMap.begin();
	unsigned int counter = 0;
	while( iter != sensorIDMap.end() ) {
		bool excluded = false;
		for(size_t i = 0; i < _excludePlanes.size(); i++) {
			if(_excludePlanes[i] == counter) {
				excluded = true;
				break;
			}
		}
		if(!excluded) _orderedSensorID_wo_excluded.push_back( iter->second );
		_orderedSensorID.push_back( iter->second );

		++iter;
		++counter;
	}

	//consistency
	if(sensorIDVec.size() != _nPlanes) {
		streamlog_out( ERROR2 ) << "the number of detected planes is " << _nPlanes << " but only " << _siPlaneZPosition.size() << " layer z positions were found!"  << std::endl;
		throw InvalidParameterException("number of layers and layer z positions mismatch");
	}
	// Initialize number of excluded planes
	_nExcludePlanes = _excludePlanes.size();

	streamlog_out( MESSAGE2 ) << "Number of planes excluded from the alignment fit: " << _nExcludePlanes << std::endl;
	streamlog_out( MESSAGE4 ) << "End of initialisation" << std::endl;
}

void EUTelPedeGEAR::processRunHeader(LCRunHeader* rdr) {
	auto header = std::make_unique<EUTelRunHeaderImpl>(rdr);
	header->addProcessor( type() ) ;

	// this is the right place also to check the geometry ID. This is a
	// unique number identifying each different geometry used at the
	// beam test. The same number should be saved in the run header and
	// in the xml file. If the numbers are different, instead of barely
	// quitting ask the user what to do.

	if( header->getGeoID() != (int)geo::gGeometry().getSiPlanesLayoutID() ) {
		streamlog_out( ERROR2 ) << "Error during the geometry consistency check: " << std::endl;
		streamlog_out( ERROR2 ) << "The run header says the GeoID is " << header->getGeoID() << std::endl;
		streamlog_out( ERROR2 ) << "The GEAR description says is     " << geo::gGeometry().getSiPlanesLayoutID() << std::endl;
	}

	// increment the run counter
	++_iRun;
}

void EUTelPedeGEAR::processEvent(LCEvent* /*event*/) {
	/*NOP NOP NOP*/
}

void EUTelPedeGEAR::end() {
	//TODO:  check if steering file exists
	std::ifstream pedeSteerFile(_pedeSteerfileName.c_str());
	if(pedeSteerFile.good()) {
		pedeSteerFile.close();
		streamlog_out( MESSAGE2 ) << "Found pede steer file, continuing ..." << std::endl;
	} else {
		pedeSteerFile.close();
		streamlog_out( ERROR5 ) << "Could not find pede steer file: " << _pedeSteerfileName << " EXITING!" << std::endl; 
		return;
	}

	std::string command = "pede " + _pedeSteerfileName;

	streamlog_out( MESSAGE5 ) << "Starting pede...: " << command.c_str() << std::endl;

	bool encounteredError = false;
      
	//run pede and create a streambuf that reads its stdout and stderr
	redi::ipstream pede( command.c_str(), redi::pstreams::pstdout|redi::pstreams::pstderr ); 
      
	if(!pede.is_open()) {
		streamlog_out( ERROR5 ) << "Pede cannot be executed: command not found in the path" << std::endl;
		encounteredError = true;	  
	} else {
		//output multiplexing: parse pede output in both stdout and stderr and echo messages accordingly
		char buf[1024];
		std::streamsize n;
		std::stringstream pedeoutput; // store stdout to parse later
		std::stringstream pedeerrors;
		bool finished[2] = { false, false };
		
		while(!finished[0] || !finished[1]) {
			if(!finished[0]) {
				while( (n = pede.err().readsome(buf, sizeof(buf))) > 0 ) {
					streamlog_out( ERROR5 ).write(buf, n).flush();
					std::string error(buf, n);
					pedeerrors << error;
					encounteredError = true;
				}
				if(pede.eof()) {
					finished[0] = true;
					if(!finished[1]) pede.clear();
				}
			}

			if(!finished[1]) {
				while((n = pede.out().readsome(buf, sizeof(buf))) > 0) {
					streamlog_out( MESSAGE4 ).write(buf, n).flush();
					std::string output (buf, n);
					pedeoutput << output;
				}
				if(pede.eof()) {
					finished[1] = true;
					if(!finished[0]) pede.clear();
				}
			}
		}

		// pede does not return exit codes on some errors (in V03-04-00)
		// check for some of those here by parsing the output
		const char* pch = strstr(pedeoutput.str().data(),"Too many rejects");
		if(pch) {
			streamlog_out( ERROR5 ) << "Pede stopped due to the large number of rejects. " << std::endl;
			encounteredError = true;
		}
	
		const char* pch0 = strstr(pedeoutput.str().data(),"Sum(Chi^2)/Sum(Ndf) = ");
		if(pch0 != 0) {
			streamlog_out( DEBUG5 ) << " Parsing pede output for final chi2/ndf result.. " << std::endl;
			//search for the equal sign after which the result for chi2/ndf is stated within the next 80 chars 
			//(with offset of 22 chars since pch points to beginning of "Sum(..." string just found)
			char* pch = (char*)((memchr (pch0+22, '=', 180)));

	    		if(pch != NULL) {
				char str[16];
				//now copy the numbers after the equal sign
				strncpy( str, pch+1, 15 );
				str[15] = '\0';   /* null character manually added */
				//TODO: monitor the chi2/ndf in CDash when running tests
				//CDashMeasurement meas_chi2ndf("chi2_ndf",atof(str));  cout << meas_chi2ndf; // output only if DO_TESTING is set
				streamlog_out( MESSAGE6 ) << "Final Sum(Chi^2)/Sum(Ndf) = " << str << std::endl;
			}
		}

		//wait for the pede execution to finish
		pede.close();

		//check the exit value of pede / react to previous errors
		if( pede.rdbuf()->status() == 0 && !encounteredError) {
			streamlog_out( MESSAGE7 ) << "Pede successfully finished" << std::endl;
		} else {
			streamlog_out( ERROR5 ) << "Problem during Pede execution, exit status: " << pede.rdbuf()->status() << ", error messages (repeated here): " << std::endl;
			streamlog_out( ERROR5 ) << pedeerrors.str() << std::endl;
			// TODO: decide what to do now; exit? and if, how?
			streamlog_out( ERROR5 ) << "Will exit now" << std::endl;
			 //exit(EXIT_FAILURE); // FIXME: can lead to (ROOT?) seg faults - points to corrupt memory? run valgrind...
			return; // does fine for now
		}

		//reading back the millepede.res file and getting the results.
		std::string millepedeResFileName = "millepede.res";

		streamlog_out( MESSAGE6 ) 	<< "Reading back the " << millepedeResFileName << std::endl;

		//open the millepede ASCII output file
		std::ifstream millepede( millepedeResFileName.c_str() );


		if( millepede.bad() || !millepede.is_open() ) {
			streamlog_out( ERROR4 )	<< "Error opening the " << millepedeResFileName << std::endl;
		} else {
			std::vector<double> tokens;
			std::stringstream tokenizer;
			std::string line;

			// get the first line and throw it away since it is a comment!
			std::getline( millepede, line );

			int counter = 0;
			int sensorID = _orderedSensorID.at(counter); 

			while( !millepede.eof() ) {
				bool goodLine = true;
				unsigned int numpars = 0;
				
				if(_alignMode != 3) {
					numpars = 3;
				} else {
					numpars = 6;
				}

				double xOff = 0;
				double yOff = 0;
				double zOff = 0;
			/*	double xOffErr = 0;
				double yOffErr = 0;
				double zOffErr = 0;  */
				double alpha = 0;
				double beta = 0;
				double gamma = 0;
			/*	double alphaErr = 0;
				double betaErr = 0;
				double gammaErr = 0; */
	
				for( unsigned int iParam = 0 ; iParam < numpars ; ++iParam ) {
					std::getline( millepede, line );

					if(line.empty()) {
						goodLine = false;
						continue;
					}

					tokens.clear();
					tokenizer.clear();
					tokenizer.str( line );

					double buffer;
					//check that all parts of the line are non zero
					while( tokenizer >> buffer ) {
						tokens.push_back( buffer ) ;
					}
					if( ( tokens.size() == 3 ) || ( tokens.size() == 6 ) || (tokens.size() == 5) ) {
						goodLine = true;
					} else {
						goodLine = false;
					}
				
				//Remove comments to read in uncertainty
				//	bool isFixed = (tokens.size() == 3);

					if(_alignMode != 3) {
						if( iParam == 0 ) {
							xOff			= tokens[1]/1000.;
				//			if(!isFixed) xOffErr	= tokens[4]/1000.;
						}
						if( iParam == 1 ) {
							yOff			= tokens[1]/1000.;
				//			if(!isFixed) yOffErr	= tokens[4]/1000.;
						}
						if( iParam == 2 ) {
							gamma			= tokens[1];
				//			if(!isFixed) gammaErr	= tokens[4];
						}
					} else {
						if( iParam == 0 ) {
							xOff			= tokens[1]/1000.;
				//			if(!isFixed) xOffErr	= tokens[4]/1000.;                    
						}
						if( iParam == 1 ) {
							yOff 			= tokens[1]/1000.;
				//			if(!isFixed) yOffErr	= tokens[4]/1000.;
						}
						if( iParam == 2 ) {
							zOff			= tokens[1]/1000.;
				//			if(!isFixed) zOffErr	= tokens[4]/1000.;
						}
						if( iParam == 3 ) {
							alpha			= tokens[1];
				//			if(!isFixed) alphaErr	= tokens[4];
						} 
						if( iParam == 4 ) {
							beta			= tokens[1];
				//			if(!isFixed) betaErr	= tokens[4];
						} 
						if( iParam == 5 ) {
							gamma			= tokens[1];
				//			if(!isFixed) gammaErr	= tokens[4];
						} 
					}
				}

				// right place to add the constant to the collection
				if( goodLine ) {
					sensorID = _orderedSensorID.at( counter );
					std::cout 	<< "Alignment on sensor " << sensorID << " determined to be: xOff: " << xOff << ", yOff: " << yOff << ", zOff: " << zOff << ", alpha: " 
							<< alpha << ", beta: " << beta << ", gamma: " << gamma << std::endl;

					//The old rotation matrix is well defined by GEAR file
					Eigen::Matrix3d rotOld = geo::gGeometry().rotationMatrixFromAngles( sensorID);
					//The new rotation matrix is obtained via the alpha, beta, gamma from MillepedeII
					Eigen::Matrix3d rotAlign = geo::gGeometry().rotationMatrixFromAngles( -alpha, -beta, -gamma);
					//The corrected rotation is given by: rotAlign*rotOld, from this rotation we can extract the
					//updated alpha', beta' and gamma'
					Eigen::Vector3d newCoeff = geo::gGeometry().getRotationAnglesFromMatrix(rotAlign*rotOld);

					//std::cout << "Old rotation matrix: " << rotOld << std::endl; 
					//std::cout << "Align rotation matrix: " << rotAlign << std::endl; 
					//std::cout << "Updated coefficients: " << newCoeff*57.29 << std::endl; 
					std::cout << "This results in the updated rotations (alpha', beta', gamma'): " << newCoeff[0] << ", " << newCoeff[1] << ", " << newCoeff[2] << std::endl;
					
					Eigen::Vector3d oldOffset;
					oldOffset << geo::gGeometry().siPlaneXPosition(sensorID), geo::gGeometry().siPlaneYPosition(sensorID), geo::gGeometry().siPlaneZPosition(sensorID);
					//Eigen::Vector3d newOffset = rotAlign*oldOffset;

					geo::gGeometry().setPlaneXPosition(sensorID, oldOffset[0]-xOff);
					geo::gGeometry().setPlaneYPosition(sensorID, oldOffset[1]-yOff);
					geo::gGeometry().setPlaneZPosition(sensorID, oldOffset[2]-zOff);
					
					geo::gGeometry().setPlaneXRotationRadians(sensorID,  newCoeff[0]);
					geo::gGeometry().setPlaneYRotationRadians(sensorID,  newCoeff[1]);
					geo::gGeometry().setPlaneZRotationRadians(sensorID,  newCoeff[2]);

					counter++;
				}
			}

		}
		millepede.close();
	}
	marlin::StringParameters* MarlinStringParams = marlin::Global::parameters;
	std::string outputFilename = (MarlinStringParams->getStringVal("GearXMLFile")).substr(0, (MarlinStringParams->getStringVal("GearXMLFile")).size()-4);
	std::cout << "GEAR Filename: " << outputFilename+_GEARFileSuffix+".xml" << std::endl;
	geo::gGeometry().writeGEARFile(outputFilename+_GEARFileSuffix+".xml");
	streamlog_out( MESSAGE2 ) << std::endl << "Successfully finished" << std::endl;
}
