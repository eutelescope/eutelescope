/* 
 * File:   EUTelGBLFitter.h
 * Contact: alexander.morton975@gmail.com
 *
 */

#ifdef USE_GBL

#ifndef EUTELGBLFITTER_H
#define	EUTELGBLFITTER_H

// mother class:
#include "EUTelTrackFitter.h"

// eutelescope includes ".h"
#include "EUTelUtility.h"
#include "EUTelTrack.h"
#include "EUTelState.h"
#include "EUTelMillepede.h"

// EVENT includes
#include <IMPL/TrackerHitImpl.h>
#include <EVENT/LCCollection.h>

// LCIO includes
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include "lcio.h"
#include "LCIOTypes.h"

// ROOT
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TMatrixD.h"
#else
#error *** You need ROOT to compile this code.  *** 
#endif

// GBL
#include "include/GblTrajectory.h"
#include "include/MilleBinary.h"

// system includes <>
#include <map>
#include <string>
#include <algorithm>
#include <functional>
#include <iostream>

namespace eutelescope {

	class EUTelGBLFitter :  public EUTelTrackFitter {
        
    private:
			DISALLOW_COPY_AND_ASSIGN(EUTelGBLFitter)        // prevent users from making (default) copies of processors
      
    public:
			EUTelGBLFitter();
			~EUTelGBLFitter();
			//SET
			void setMomentsAndStartEndScattering(EUTelState& state);
			void setInformationForGBLPointList(EUTelTrack& track, std::vector< gbl::GblPoint >& pointList);
			void setMeasurementGBL(gbl::GblPoint& point, const double *hitPos, double statePos[3], double combinedCov[4], TMatrixD projection);
			void getKinkInformationToTrack(gbl::GblTrajectory* traj, std::vector< gbl::GblPoint >& pointList,EUTelTrack &track);
            TMatrixD getFullJacobian(TVector3 momStart, TVector3 momEnd, int locationStart, int locationEnd, double distance, double min );
			void setPointVec( std::vector< gbl::GblPoint >& pointList, gbl::GblPoint& point);
			void setPairAnyStateAndPointLabelVec(gbl::GblTrajectory*);
			void setPairMeasurementStateAndPointLabelVec(std::vector< gbl::GblPoint >& pointList);
			void setAlignmentToMeasurementJacobian(std::vector< gbl::GblPoint >& pointList);
			void setScattererGBL(gbl::GblPoint& point,EUTelState & state );
			void setScattererGBL(gbl::GblPoint& point,EUTelState & state,  float variance, TVectorD scat );
			void setLocalDerivativesToPoint(gbl::GblPoint& point, float distanceFromKinkTargetToNextPlane );
			void setPointListWithNewScatterers(std::vector< gbl::GblPoint >& pointList,EUTelState & state, std::vector<float> variance );
			void setMeasurementCov(EUTelState& state);
			inline void setAlignmentMode( int number) {
				this->_alignmentMode = number;
			}
			inline void setBeamCharge(double beamQ) {
				this->_beamQ = beamQ;
			}
      inline void setBeamEnergy(double beamE) {
				this->_eBeam = beamE;
			}
			void SetMilleBinary(gbl::MilleBinary* _mille) {
				this->_mille = _mille;
			}
			void setMillepede( EUTelMillepede* Mille ) { _MilleInterface =  Mille; }
			void SetJacobain(TMatrixD matrix ){
				_jacobianAlignment = matrix;
			}
			void SetAlignmentVariablesList(std::vector<int> globalLabels){
				_globalLabels = globalLabels;
			}
			void setParamterIdXResolutionVec( const std::vector<float>& );
			void setParamterIdYResolutionVec( const std::vector<float>& );
			void setExcludeFromFitPlanes(const std::vector<int>&);
			void setMEstimatorType( const std::string& _mEstimatorType );
			//GET
			float getPositionOfSecondScatter(float start, float end);
			gbl::GblPoint getLabelToPoint(std::vector<gbl::GblPoint> & pointList, unsigned  int label);
			void getResidualOfTrackandHits(gbl::GblTrajectory* traj, std::vector< gbl::GblPoint > pointList, EUTelTrack& track, std::map< int, std::map< float, float > > & SensorResidual, std::map< int, std::map< float, float > >& sensorResidualError);
			inline int getAlignmentMode() const {
				return _alignmentMode;
			}
			inline double getBeamEnergy() const {
				return _eBeam;
			}
			inline std::vector<double> getCorrectionsTotal() const {
				std::vector<double> correctionsTotal;
				correctionsTotal.push_back(_omegaCorrections);	
				correctionsTotal.push_back(_intersectionLocalXZCorrections);	
				correctionsTotal.push_back(	_intersectionLocalYZCorrections);	
				correctionsTotal.push_back(_localPosXCorrections);
				correctionsTotal.push_back(_localPosYCorrections);
				return correctionsTotal;
			}
			const std::vector<  float>& getParamterIdXResolutionVec() const;
			const std::vector<  float>& getParamterIdYResolutionVec() const;
			std::string getMEstimatorType() const;

			//TEST
			void testUserInput();
			void testTrack(EUTelTrack& track);
			void testDistanceBetweenPoints(double* position1,double* position2);
			//COMPUTE
			void computeTrajectoryAndFit(gbl::GblTrajectory* traj, double* chi2, int* ndf, int & ierr);
			std::vector<float> computeVarianceForEachScatterer(EUTelState & state);
			//OTHER FUNCTIONS
			void resetPerTrack();
			void findScattersZPositionBetweenTwoStates();
			TMatrixD findScattersJacobians(EUTelState state, EUTelState nextTrack);
			void updateTrackFromGBLTrajectory(gbl::GblTrajectory* traj, EUTelTrack& track, std::map<int,std::vector<double> >& mapSensorIDToCorrectionVec );
			void prepareLCIOTrack( gbl::GblTrajectory*, std::vector<const IMPL::TrackImpl*>::const_iterator&, double, int); 
			void prepareMilleOut( gbl::GblTrajectory* );
            TVector3 transVecGlobalToLocal(TVector3 input, int location);

			//VARIABLES
			int _alignmentMode;
			std::vector<EUTelState> _measurementStatesInOrder;
			std::vector<EUTelState> _statesInOrder;
			double _omegaCorrections;	
			double _intersectionLocalXZCorrections;	
			double _intersectionLocalYZCorrections;	
			double _localPosXCorrections;
			double _localPosYCorrections;
			double _beamQ;
			double _eBeam;
			float _normalMean;
			float _normalVariance;
			float _start; //This is the first scatterer location
			float _end; //This is arc length from the first plane to the second.
			/** Outlier downweighting option */
			std::string _mEstimatorType;
			/** Milipede binary file handle */
			gbl::MilleBinary* _mille;
			std::string _binaryname;
			TMatrixD _jacobianAlignment;
			std::vector<TMatrixD> _scattererJacobians;
			std::vector<float> _scattererPositions;
			std::vector<int> _globalLabels;
			/** Parameter resolutions */
			std::vector<int> _parameterIdPlaneVec;
			/** Parameter resolutions */
			std::map< int,  float> _parameterIdXResolutionVec;
			/** Parameter resolutions */
			std::map< int,  float> _parameterIdYResolutionVec;
			/** Parameter ids */
			std::map<int,int> _parameterIdXShiftsMap;
			/** Parameter ids */
			std::map<int,int> _parameterIdYShiftsMap;
			std::map<int,int> _parameterIdZShiftsMap;
			/** Parameter ids */
			std::map<int,int> _parameterIdXRotationsMap;
			/** Parameter ids */
			std::map<int,int> _parameterIdYRotationsMap;
			/** Parameter ids */
			std::map<int,int> _parameterIdZRotationsMap;
			std::vector< std::pair< EUTelState, unsigned int> > _vectorOfPairsMeasurementStatesAndLabels;//This is used within alignment since you need to associate MEASUREMENT states to  labels
			std::vector< std::pair< EUTelState, int> > _vectorOfPairsStatesAndLabels;//This is used in track fit since you want to associate ANY states to labels.
			unsigned int _counter_num_pointer;
			bool _kinkAngleEstimation; //This used to determine if the correction matrix from the GBL fit is 5 or 7 elements long. 
			EUTelMillepede* _MilleInterface;
			std::vector<int> _sensorIDVec;
        
    };
}
#endif	/* EUTELGBLFITTER_H */
#endif
