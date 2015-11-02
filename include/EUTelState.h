#ifndef EUTELSTATE_H
#define	EUTELSTATE_H

#include "EUTelUtility.h"
#include "EUTelBlock.h"

// ROOT
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#endif
#include "EUTelHit.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
//Eigen
#include <Eigen/Core>

namespace eutelescope {

	class  EUTelState{
		public: 
            Block block;
            std::vector<unsigned int> GBLLabels;

			EUTelState();
			EUTelState(EUTelState *state);
			//getters
			EUTelHit& getHit();
            EUTelHit getHitCopy() const;

			int getDimensionSize() const ;
			TVectorD getStateVec();
			float getSlopeX() const; 
			float getSlopeY() const; 
            float getSlopeXGlobal() const; 
            float getSlopeYGlobal() const; 

            std::vector<double> getLCIOOutput();
			float getArcLengthToNextState() const {return _arcLength;} 
			TVector3 getPositionGlobal() const; 
            Eigen::Vector3d getPositionGlobalEig() const; 
			bool getStateHasHit() const;
            /// This will get the link between the local and global frames. 
            /// GBL Fitter expects the link to be give from global to local. 
            /// MUST INVERT FOR USE!!
            /**
             * \param [out] Projection Calculated using the incidence in the global frame and normal of the senso
             */

			TMatrixD getProjectionMatrix() const;
			TVector3 getIncidenceUnitMomentumVectorInLocalFrame();
			TMatrixDSym getScatteringVarianceInLocalFrame(double const & var );
			double getRadFracAir() const ;
			double getRadFracSensor() const ;
            TVector3 getDirLocal() const; 
            TVector3 getDirGlobal() const; 
            Eigen::Vector3d getDirGlobalEig() const; 
            //STATE PARAMETERS:
			TVectorD getKinks() const;
			TVectorD getKinksMedium1() const;
			TVectorD getKinksMedium2() const;
			const double* getPosition() const ; 
			int	getLocation() const;
			float getDirLocalX() const {return _dirLocalX;}
			float getDirLocalY() const {return _dirLocalY;}
			float getDirLocalZ() const {return _dirLocalZ;}
            TMatrixD getCov();
            //END OF STATE PARAMETERS
			//setters
            void setHit(EUTelHit hit);
            void setHit(EVENT::TrackerHit* hit);
			void setDimensionSize(int dimension);
			void setLocation(int location);
			void setDirLocalX(double dirX);
			void setDirLocalY(double dirY);
			void setDirLocalZ(double dirZ);
            void setDirFromGloSlope(std::vector<double> slopes);
            void setDirFromLocSlope(std::vector<double> slopes);
			void setLocalDirGlobalDir(TVector3 momentumIn);
            void setTrackFromLCIOVec(std::vector<double> input);
			void setPositionLocal(float position[]);
			void setPositionLocal(double position[]);
			void setPositionGlobal(float positionGlobal[]);
            void setPositionGlobal(double positionGlobal[]);
			void setStateUsingCorrection(TVectorD corrections);
			void setArcLengthToNextState(float arcLength){_arcLength = arcLength;} 
			void setKinks(TVectorD kinks);
			void setKinksMedium1(TVectorD kinks);
			void setKinksMedium2(TVectorD kinks);
			void setRadFrac(double plane, double air);
            void setCov(TMatrixD cov);

			//initialise
			void initialiseCurvature();
			//print
			void print();
            //clear
            void clear();

			bool operator<(const EUTelState compareState ) const;
			bool operator==(const EUTelState compareState ) const;
			bool operator!=(const EUTelState compareState ) const;

  	protected:
            EUTelHit _hit;
            int _dimension;
            int _location; 
            double _position[3];
            bool _stateHasHit;
            TVectorD _kinks;
            TVectorD _kinksMedium1;
            TVectorD _kinksMedium2;
            TMatrixD _cov;
            double _dirLocalX;
            double _dirLocalY;
            double _dirLocalZ; 
            double _radFracSensor;
            double _radFracAir;
            float _arcLength;
			float _covCombinedMatrix[4];
	};
}
#endif
