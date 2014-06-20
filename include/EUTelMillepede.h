#ifndef EUTELMILLEPEDE_H
#define	EUTELMILLEPEDE_H

#include "EUTelUtility.h"
#include "EUTelTrackStateImpl.h"


// system includes <>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <iomanip>
#include <iterator>
#include <algorithm>

namespace eutelescope {

    class EUTelMillepede {
        
    private:
        DISALLOW_COPY_AND_ASSIGN(EUTelMillepede);        // prevent users from making (default) copies of processors

   public:
        EUTelMillepede();
				EUTelMillepede(int alignmentMode);

        ~EUTelMillepede();

				//This set the number given by the processor to a aligment mode string
				void SetAlignmentMode(int alignmentMode);
				//This take a state and outputs a its alignment jacobian given the alignment mode
				int CreateAlignmentToMeasurementJacobian( EUTelTrackStateImpl* state);

				int CreateAlignmentToMeasurementJacobian( float x,float y, float slopeXvsZ, float slopeYvsZ);

				void CreateGlobalLabels(EUTelTrackStateImpl* state);

				void CreateGlobalLabels(int iPlane);

				void FillMilleParametersLabels();

		protected:
			int alignmentMode;
			Utility::AlignmentMode _alignmentMode =  Utility::noAlignment;
			TMatrixD* _jacobian; 
			std::vector<int> _globalLabels;
			std::map<int, int> _xShiftsMap;
			std::map<int, int> _yShiftsMap;
			std::map<int, int> _zShiftsMap;
      std::map<int, int> _xRotationsMap;
      std::map<int, int> _yRotationsMap;
      std::map<int, int> _zRotationsMap;

    };

}



#endif	/* EUTelMillepede_H */
