#ifndef EUTELSTATE_H
#define	EUTELSTATE_H

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

	class  EUTelState : public IMPL::TrackImpl{
		public: 
			EUTelState();
			EUTelState(EUTelState *state);
			//getters
			EVENT::TrackerHit* getHit();
			int getDimensionSize() const ;
			int	getLocation() const;
			TMatrixDSym getTrackStateCov() const;
			TVectorD getTrackStateVec() const ;
			inline float  getBeamCharge() const  { return getdEdxError();}
			inline float  getBeamEnergy() const {return getdEdx();}
			float getDirectionXY() const {return getPhi();}
			float getDirectionYZ() const {return getTanLambda();}
			float* getPosition() const; 
			void getCombinedHitAndStateCovMatrixInLocalFrame(double (&cov)[4]) const;
			bool getIsThereAHit() const;
			TMatrixD getProjectionMatrix() const;
			TVector3 getIncidenceUnitMomentumVectorInLocalFrame();
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
			bool operator<(const EUTelState compareState ) const;
			bool operator==(const EUTelState compareState ) const;

  	private:
			float _covCombinedMatrix[4];
	};

}
#endif
