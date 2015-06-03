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
#include "EUTelHit.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

namespace eutelescope {

	class  EUTelState{
		public: 
			EUTelState();
			EUTelState(EUTelState *state);
			//getters
			EUTelHit getHit();
			int getDimensionSize() const ;
			int	getLocation() const;
			TMatrixDSym getStateCov() const;
			TVectorD getStateVec();
            TVector3 getMomLocal();
			float getMomLocalX() const {return _momLocalX;}
			float getMomLocalY() const {return _momLocalY;}
			float getMomLocalZ() const {return _momLocalZ;}
			TVector3 getMomGlobal() const ;
            std::vector<double> getLCIOOutput();
			float getArcLengthToNextState() const {return _arcLength;} 
			const float* getPosition() const ; 
			TVector3 getPositionGlobal() const; 
			void getCombinedHitAndStateCovMatrixInLocalFrame(double (&cov)[4]) const;
			bool getStateHasHit() const;
			TMatrixD getProjectionMatrix() const;
			TVector3 getIncidenceUnitMomentumVectorInLocalFrame();
			TMatrixDSym getScatteringVarianceInLocalFrame();
			TMatrixDSym getScatteringVarianceInLocalFrame(float variance);
			TVectorD getKinks() const;
			TVectorD getKinksMedium1() const;
			TVectorD getKinksMedium2() const;
			double getRadFracAir() const ;
			double getRadFracSensor() const ;
			//setters
            void setHit(EUTelHit hit);
            void setHit(EVENT::TrackerHit* hit);
			void setDimensionSize(int dimension);
			void setLocation(int location);
			void setMomLocalX(float momX);
			void setMomLocalY(float momY);
			void setMomLocalZ(float momZ);
			void setLocalMomentumGlobalMomentum(TVector3 momentumIn);
            void setTrackFromLCIOVec(std::vector<double> input);
            //!Template input for setting local position of hit  
            /*!
             * @param position of hit on plane
             */
			void setPositionLocal(float position[]);
			void setPositionLocal(double position[]);
			void setPositionGlobal(float positionGlobal[]);
			void setCombinedHitAndStateCovMatrixInLocalFrame(double cov[4]);
			void setStateUsingCorrection(TVectorD stateVec);
			void setArcLengthToNextState(float arcLength){_arcLength = arcLength;} 
			void setKinks(TVectorD kinks);
			void setKinksMedium1(TVectorD kinks);
			void setKinksMedium2(TVectorD kinks);
			void setRadFrac(double plane, double air);

			//initialise
			void initialiseCurvature();
			//find
			bool findIntersectionWithCertainID(int nextsensorID, float intersectionPoint[], TVector3& momentumAtIntersection, float& arcLength, int& newNextPlaneID );
			//compute
//			TMatrix computePropagationJacobianFromLocalStateToNextLocalState(TVector3 momentumEnd, float arcLength,float nextPlaneID);
			float computeRadLengthsToEnd( std::map<const int,double> & mapSensor, std::map<const int ,double> & mapAir );
            //
            EUTelHit _hit;
            int _dimension;
            int _location; 
            float _position[3];
            bool _stateHasHit;
            TVectorD _kinks;
            TVectorD _kinksMedium1;
            TVectorD _kinksMedium2;
            float _momLocalX;
            float _momLocalY;
            float _momLocalZ; 
            double _radFracSensor;
            double _radFracAir;
            float _arcLength;

			//print
			void print();
            //clear
            void clear();

			bool operator<(const EUTelState compareState ) const;
			bool operator==(const EUTelState compareState ) const;
			bool operator!=(const EUTelState compareState ) const;

  	private:
			float _covCombinedMatrix[4];
	};
}
#endif
