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
			EUTelHit& getHit();
			int getDimensionSize() const ;
			TVectorD getStateVec();
			float getMomLocalX() const {return _momLocalX;}
			float getMomLocalY() const {return _momLocalY;}
			float getMomLocalZ() const {return _momLocalZ;}
			float getSlopeX() const; 
			float getSlopeY() const; 
			TVector3 getMomGlobal() const ;
            std::vector<double> getLCIOOutput();
			float getArcLengthToNextState() const {return _arcLength;} 
			TVector3 getPositionGlobal() const; 
			void getCombinedHitAndStateCovMatrixInLocalFrame(double (&cov)[4]) const;
			bool getStateHasHit() const;
			TMatrixD getProjectionMatrix() const;
			TVector3 getIncidenceUnitMomentumVectorInLocalFrame();
			TMatrixDSym getScatteringVarianceInLocalFrame();
			TMatrixDSym getScatteringVarianceInLocalFrame(float variance);
			double getRadFracAir() const ;
			double getRadFracSensor() const ;
            //STATE PARAMETERS:
			TVectorD getKinks() const;
			TVectorD getKinksMedium1() const;
			TVectorD getKinksMedium2() const;
			const float* getPosition() const ; 
            TVector3 getMomLocal() const ;
			int	getLocation() const;
			TMatrixDSym getStateCov() const;
            //END OF STATE PARAMETERS
			//setters
            void setHit(EUTelHit hit);
            void setHit(EVENT::TrackerHit* hit);
			void setDimensionSize(int dimension);
			void setLocation(int location);
			void setMomLocalX(float momX);
			void setMomLocalY(float momY);
			void setMomLocalZ(float momZ);
            void setMomGlobalIncEne(std::vector<float> slopes, float energy);
            void setMomGlobalIncEne(std::vector<float> slopes, double energy);
			void setLocalMomentumGlobalMomentum(TVector3 momentumIn);
            void setTrackFromLCIOVec(std::vector<double> input);
			void setPositionLocal(float position[]);
			void setPositionLocal(double position[]);
			void setPositionGlobal(float positionGlobal[]);
            void setPositionGlobal(double positionGlobal[]);
			void setCombinedHitAndStateCovMatrixInLocalFrame(double cov[4]);
			void setStateUsingCorrection(TVectorD corrections);
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
			float computeRadLengthsToEnd( std::map<const int,double> & mapSensor, std::map<const int ,double> & mapAir );
			//print
			void print();
            //clear
            void clear();

			bool operator<(const EUTelState compareState ) const;
			bool operator==(const EUTelState compareState ) const;
			bool operator!=(const EUTelState compareState ) const;

  	private:
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
			float _covCombinedMatrix[4];
	};
}
#endif
