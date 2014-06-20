#ifndef EUTELMILLEPEDE_H
#define	EUTELMILLEPEDE_H

#include "EUTelUtility.h"
#include "EUTelTrackStateImpl.h"


namespace eutelescope {

    class EUTelMillepede {
        
    private:
        DISALLOW_COPY_AND_ASSIGN(EUTelMillepede);        // prevent users from making (default) copies of processors

   public:
        EUTelMillepede();

        ~EUTelMillepede();

				//This set the number given by the processor to a aligment mode string
				void SetAlignmentMode(int alignmentMode);
				//This take a state and outputs a its alignment jacobian given the alignment mode
				int CreateAlignmentToMeasurementJacobian( EUTelTrackStateImpl* state, TMatrixD* Jacobian );

				int CreateAlignmentToMeasurementJacobian( float x,float y, float slopeXvsZ, float slopeYvsZ, TMatrixD* Jacobian );

		protected:
			int alignmentMode;
			Utility::AlignmentMode _alignmentMode =  Utility::noAlignment;

    };

}



#endif	/* EUTelMillepede_H */
