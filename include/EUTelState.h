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
			TVectorD getStateVec();
            TVector3 getMomLocal();
			float getMomLocalX() const {return getdEdx();}
			float getMomLocalY() const {return getTanLambda();}
			float getMomLocalZ() const {return getPhi();}
			TVector3 getMomGlobal() const ;

			float getArcLengthToNextState() const {return getChi2();} 
			float* getPosition() const; 
			TVector3 getPositionGlobal() const; 
			void getCombinedHitAndStateCovMatrixInLocalFrame(double (&cov)[4]) const;
			bool getIsThereAHit() const;
			TMatrixD getProjectionMatrix() const;
			TVector3 getIncidenceUnitMomentumVectorInLocalFrame();
			TMatrixDSym getScatteringVarianceInLocalFrame();
			TMatrixDSym getScatteringVarianceInLocalFrame(float variance);
			TVectorD getKinks();
			TVectorD getRadFrac();
			//setters
			void setDimensionSize(int dimension);
			void setLocation(int location);
			void setMomLocalX(float momX);
			void setMomLocalY(float momY);
			void setMomLocalZ(float momZ);
			void setLocalMomentumGlobalMomentum(TVector3 momentumIn);
            //!Template input for setting local position of hit  
            /*!
             * @param position of hit on plane
             */
            template<class number>
			void setPositionLocal(number position[]);
			void setPositionGlobal(float positionGlobal[]);
			void setCombinedHitAndStateCovMatrixInLocalFrame(double cov[4]);
			void setStateUsingCorrection(TVectorD stateVec);
			void setArcLengthToNextState(float arcLength){setChi2(arcLength);} 
			void setKinks(TVectorD kinks);
			void setRadFrac(double plane, double air);

			//initialise
			void initialiseCurvature();
			//find
			bool findIntersectionWithCertainID(int nextsensorID, float intersectionPoint[], TVector3& momentumAtIntersection, float& arcLength, int& newNextPlaneID );
			//compute
//			TMatrix computePropagationJacobianFromLocalStateToNextLocalState(TVector3 momentumEnd, float arcLength,float nextPlaneID);
			float computeRadLengthsToEnd( std::map<const int,double> & mapSensor, std::map<const int ,double> & mapAir );

			//print
			void print();
			bool operator<(const EUTelState compareState ) const;
			bool operator==(const EUTelState compareState ) const;

  	private:
			float _covCombinedMatrix[4];
	};
    #include "EUTelState.tcc"
}
#endif
