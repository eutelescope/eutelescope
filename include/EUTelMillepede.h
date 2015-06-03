#ifndef EUTELMILLEPEDE_H
#define	EUTELMILLEPEDE_H

#include "EUTelUtility.h"
#include "EUTelTrack.h"
#include "EUTelState.h"


#include <fstream>
#include "EUTelPStream.h"

// MARLIN
#include "marlin/Exceptions.h"
#include "marlin/Global.h"
#include "marlin/Processor.h"
#include "marlin/VerbosityLevels.h"


// system includes <>
#include <map>
#include <fstream>      // std::ifstream, std::ofstream
#include <string>
#include <utility>
#include <vector>
#include <iomanip>
#include <iterator>
#include <algorithm>
#include <cstdio>

#include "include/MilleBinary.h"
#include "EUTelExceptions.h"

namespace eutelescope {

    class EUTelMillepede {
        
    private:
        DISALLOW_COPY_AND_ASSIGN(EUTelMillepede)        // prevent users from making (default) copies of processors

   public:
        EUTelMillepede();

        ~EUTelMillepede();

				//This take a state and outputs a its alignment jacobian given the alignment mode
				void computeAlignmentToMeasurementJacobian( EUTelState& state);

				void computeAlignmentToMeasurementJacobian( float x,float y, float slopeXvsZ, float slopeYvsZ);

				void setGlobalLabels(EUTelState& state);

				void setGlobalLabels(int iPlane);

				void FillMilleParametersLabels();

				void writeMilleSteeringFile(lcio::StringVec pedeSteerAddCmds);

				bool runPede();
	
				bool parseMilleOutput(std::string alignmentConstantLCIOFile, std::string gear_aligned_file);
				bool converge();
				bool checkConverged();
				void editSteerUsingRes();
				void copyFile(std::string input, std::string output);
				void outputSteeringFiles();


				void testUserInput();
				void printFixedPlanes();
				/////////////////////////set stuff!
				void setXShiftFixed(lcio::IntVec xfixed);
				void setYShiftFixed(lcio::IntVec yfixed);
				void setZShiftFixed(lcio::IntVec zfixed);
				void setXRotationsFixed(lcio::IntVec xRotfixed);
				void setYRotationsFixed(lcio::IntVec yRotfixed);
				void setZRotationsFixed(lcio::IntVec zRotfixed);
				void setPlanesExclude(lcio::IntVec exclude);
				void setSteeringFileName(std::string name);
				void setBinaryFileName(std::string binary);
				void setResultsFileName(std::string name);


				///////////////////////////////////////////get stuff
				TMatrixD& getAlignmentJacobian()  { return _jacobian; }
				std::vector<int> getGlobalParameters() { return _globalLabels; }
				///////find stuff
				bool findTooManyRejects(std::string output);


				gbl::MilleBinary * _milleGBL;


void CreateBinary();

		protected:
				TMatrixD _jacobian; //Remember you need to create the object before you point ot it
				std::vector<int> _globalLabels;
				std::map<int, int> _xShiftsMap;
				std::map<int, int> _yShiftsMap;
				std::map<int, int> _zShiftsMap;
      	std::map<int, int> _xRotationsMap;
      	std::map<int, int> _yRotationsMap;
      	std::map<int, int> _zRotationsMap;

        /** Mille steering filename */
				std::string _milleSteeringFilename;
				
				std::string _milleSteerNameOldFormat;
				int _iteration;

				std::string _milleBinaryFilename;
				//the results file
				std::string _milleResultFileName;

  		 /** Alignment X shift plane ids to be fixed */
			lcio::IntVec _fixedAlignmentXShfitPlaneIds;
        
        /** Alignment Y shift plane ids to be fixed */
				lcio::IntVec _fixedAlignmentYShfitPlaneIds;
        
        /** Alignment Z shift plane ids to be fixed */
				lcio::IntVec _fixedAlignmentZShfitPlaneIds;
        
        /** Alignment X rotation plane ids to be fixed */
				lcio::IntVec _fixedAlignmentXRotationPlaneIds;
        
        /** Alignment Y rotation plane ids to be fixed */
				lcio::IntVec _fixedAlignmentYRotationPlaneIds;
        
        /** Alignment Z rotation plane ids to be fixed */
				lcio::IntVec _fixedAlignmentZRotationPlaneIds;

        /** Alignment plane ids of planes that are excluded*/
				lcio::IntVec _alignmentPlaneIdsExclude;

    };

}



#endif	/* EUTelMillepede_H */
