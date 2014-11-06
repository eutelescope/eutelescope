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
			TMatrixDSym getStateCov() const;
			TVectorD getStateVec() const ;
			inline float  getBeamCharge() const  { return getdEdxError();}
			inline float  getBeamEnergy() const {return getdEdx();}
			float getIntersectionLocalXZ() const {return getPhi();}
			float getIntersectionLocalYZ() const {return getTanLambda();}
			float getArcLengthToNextState() const {return getChi2();} 
			float* getPosition() const; 
			TVector3 getPositionGlobal() const; 
			void getCombinedHitAndStateCovMatrixInLocalFrame(double (&cov)[4]) const;
			bool getIsThereAHit() const;
			TMatrixD getProjectionMatrix() const;
			TVector3 getIncidenceUnitMomentumVectorInLocalFrame();
			TMatrixDSym getScatteringVarianceInLocalFrame();
			TMatrixDSym getScatteringVarianceInLocalFrame(float percentageRadiationLength);
			TVectorD getKinks();
			//setters
			void setDimensionSize(int dimension);
			void setLocation(int location);
			void setBeamEnergy(float beamE);
			void setBeamCharge(float beamQ);
			void setIntersectionLocalYZ(float directionYZ);
			void setIntersectionLocalXZ(float directionXZ);
			void setLocalXZAndYZIntersectionAndCurvatureUsingGlobalMomentum(TVector3 momentumIn);
			void setPositionLocal(float position[]);
			void setPositionGlobal(float positionGlobal[]);
			void setCombinedHitAndStateCovMatrixInLocalFrame(double cov[4]);
			void setStateVec(TVectorD stateVec);
			void setArcLengthToNextState(float arcLength){setChi2(arcLength);} 
			void setKinks(TVectorD kinks);
			//initialise
			void initialiseCurvature();
			//find
			int findIntersectionWithCertainID(int nextsensorID, float intersectionPoint[],TVector3 & momentumAtIntersection, float & arcLength );
			//compute
			TVector3 computeCartesianMomentum() const ;
			TMatrix computePropagationJacobianFromLocalStateToNextLocalState(TVector3 positionEnd, TVector3 momentumEnd, float arcLength,float nextPlaneID);
			//print
			void print();
			bool operator<(const EUTelState compareState ) const;
			bool operator==(const EUTelState compareState ) const;

  	private:
			float _covCombinedMatrix[4];
	};

}
#endif
