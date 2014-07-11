
#include "EUTelUtility.h"
// ROOT
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#endif
//lcio
#include "IMPL/TrackImpl.h"
#include "EUTelGeometryTelescopeGeoDescription.h"


using namespace lcio;
using namespace marlin;
using namespace std;

namespace eutelescope {

	class  EUTelTrack : public IMPL::TrackImpl{
		public: 
			EUTelTrack();
			//getters
			int	getLocation();

			//setters
			void setLocation(int location);
			void setBeamEnergy(float beamE);
			void setBeamCharge(float beamQ);
			void setDirectionYZ(float directionYZ);
			void setDirectionXY(float directionXY);
			void setPosition(float position[]);
			//initialise
			void initialiseCurvature();
			//find
			int findIntersectionWithCertainID(int nextsensorID, float intersectionPoint[] );
			//compute
			TVector3 computeCartesianMomentum();
			TMatrix computePropagationJacobianFromStateToThisZLocation(float zPosition);
  	private:
    	float _beamQ;
			float _beamE;
	};

}
