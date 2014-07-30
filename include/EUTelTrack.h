
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
			int getDimensionSize() const ;
			int	getLocation() const;
			TMatrixDSym getTrackStateCov() const;
			TVectorD getTrackStateVec() const ;
			int getNumberOfHitsOnTrack() const;
			void getCombinedHitAndStateCovMatrixInLocalFrame(double (&cov)[4]) const;
			TMatrixD getProjectionMatrix() const;
			//setters
			void setDimensionSize(int dimension);
			void setLocation(int location);
			void setBeamEnergy(float beamE);
			void setBeamCharge(float beamQ);
			void setDirectionYZ(float directionYZ);
			void setDirectionXY(float directionXY);
			void setPosition(float position[]);
			void setCombinedHitAndStateCovMatrixInLocalFrame(double cov[4]);
			void setTrackStateVecPlusZParameter(TVectorD stateVec, float zParameter);
			//initialise
			void initialiseCurvature();
			//find
			int findIntersectionWithCertainID(int nextsensorID, float intersectionPoint[] );
			//compute
			TVector3 computeCartesianMomentum();
			TMatrix computePropagationJacobianFromStateToThisZLocation(float zPosition);
			//print
			void print();
			bool operator<(const EUTelTrack compareState ) const;

  	private:
    	float _beamQ;
			float _beamE;
			float _covCombinedMatrix[4];
	};

}
