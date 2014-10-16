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
#include "EUTelUtilityRungeKutta.h"
#include "EUTelEquationsOfMotion.h"
#include "EUTelTrackStateImpl.h"
#include "EUTelTrack.h"
#include "EUTelState.h"
#include "EUTelTrackImpl.h"
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
        explicit EUTelGBLFitter(std::string name);
        virtual ~EUTelGBLFitter();

        // do some clean up of internal data structures
        // will be automatically run when calling EUTelGBLFitter::FitTracks()
        void resetPerTrack();

        /** Fit tracks */
        // public: 
				void findScattersZPositionBetweenTwoStates(EUTelState& state);
				TMatrixD findScattersJacobians(EUTelState state, EUTelState nextTrack);
				void setInformationForGBLPointList(EUTelTrack& track, std::vector< gbl::GblPoint >& pointList);
//				void setInformationForGBLPointListForAlignment(EUTelTrack& track, std::vector< gbl::GblPoint >& pointList);
				void testTrack(EUTelTrack& track);
				void changejacobainGBL(TMatrixD & input, TMatrixD & output);
				void FindHitIfThereIsOne(EUTelTrackImpl* trackimpl, EVENT::TrackerHit* hit, EUTelTrackStateImpl* state);
				void testDistanceBetweenPoints(double* position1,double* position2);
				void setMeasurementGBL(gbl::GblPoint& point, const double *hitPos, double statePos[3], double combinedCov[4], TMatrixD projection);
				void setKinkInformationToTrack(gbl::GblTrajectory* traj, std::vector< gbl::GblPoint >& pointList,EUTelTrack &track);



				void computeTrajectoryAndFit(std::vector< gbl::GblPoint >& pointList,  gbl::GblTrajectory* traj, double* chi2, int* ndf, int & ierr);

				void updateTrackFromGBLTrajectory(gbl::GblTrajectory* traj,std::vector< gbl::GblPoint >& pointList, EUTelTrack & track, map<int, vector<double> > &  mapSensorIDToCorrectionVec );
				void setPointVec( std::vector< gbl::GblPoint >& pointList, gbl::GblPoint& point);
				std::vector<EUTelState> _measurementStatesInOrder;
				std::vector<EUTelState> _statesInOrder;
				void setPairAnyStateAndPointLabelVec(std::vector< gbl::GblPoint >& pointList, gbl::GblTrajectory*);
				void setPairMeasurementStateAndPointLabelVec(std::vector< gbl::GblPoint >& pointList);
				gbl::GblPoint getLabelToPoint(std::vector<gbl::GblPoint> & pointList, int label);
				void setAlignmentToMeasurementJacobian(EUTelTrack& track, std::vector< gbl::GblPoint >& pointList);
				void setScattererGBL(gbl::GblPoint& point,EUTelState & state );
				void setScattererGBL(gbl::GblPoint& point,EUTelState & state,  float  percentageRadiationLength);
				void setPointListWithNewScatterers(std::vector< gbl::GblPoint >& pointList,EUTelState & state, float  percentageRadiationLength);
				void getResidualOfTrackandHits(gbl::GblTrajectory* traj, std::vector< gbl::GblPoint > pointList, EUTelTrack& track, map< int, map< float, float > > & SensorResidual, map< int, map< float, float > >& sensorResidualError);

        /*
         */  
        void TrackCandidatesToGBLTrajectory( vector<const IMPL::TrackImpl*>::const_iterator&  );

        /*
         */
        void PerformFitGBLTrajectory( gbl::GblTrajectory* ,  vector<const IMPL::TrackImpl*>::const_iterator& );
       

        void FitSingleTrackCandidate(EVENT::TrackVec::const_iterator& itTrkCand);

				void setMeasurementCov(EUTelState& state);

				inline void setAlignmentMode( int number) {
            this->_alignmentMode = number;
        }

        inline int getAlignmentMode() const {
            return _alignmentMode;
        }
				int _alignmentMode;
        inline void SetBeamCharge(double beamQ) {
            this->_beamQ = beamQ;
        }

        inline double getBeamCharge() const {
            return _beamQ;
        }

        inline void SetBeamEnergy(double beamE) {
            this->_eBeam = beamE;
        }

        inline double getBeamEnergy() const {
            return _eBeam;
        }





        void SetMilleBinary(gbl::MilleBinary* _mille) {
            this->_mille = _mille;
        }

				void setMillepede( EUTelMillepede* Mille ) { _MilleInterface =  Mille; }
			
        inline void setChi2Cut(double chi2cut) {
            this->_chi2cut = chi2cut;
        }

        inline double GetChi2Cut() const {
            return _chi2cut;
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


				void SetJacobain(TMatrixD matrix ){
					_jacobianAlignment = matrix;
				}
		
				void SetAlignmentVariablesList(std::vector<int> globalLabels){
					_globalLabels = globalLabels;
				}
 
        void setParamterIdPlaneVec( const std::vector<int>& );
 
        void setParamterIdXResolutionVec( const std::vector<float>& );

        void setParamterIdYResolutionVec( const std::vector<float>& );
       
        void setParamterIdXRotationsMap( const std::map<int, int>& );
        
        void setParamterIdYRotationsMap( const std::map<int, int>& );
        
        void setParamterIdZRotationsMap( const std::map<int, int>& );
        
        void setParamterIdZShiftsMap( const std::map<int, int>& );
        
        void setParamterIdYShiftsMap( const std::map<int, int>& );
        
        void setParamterIdXShiftsMap(const std::map<int, int>& );
				        
        const std::vector<  int>& getParamterIdPlaneVec() const;

        const std::vector<  float>& getParamterIdXResolutionVec() const;

        const std::vector<  float>& getParamterIdYResolutionVec() const;

        const std::map<int, int>& getParamterIdXRotationsMap() const;
        
        const std::map<int, int>& getParamterIdYRotationsMap() const;
        
        const std::map<int, int>& getParamterIdZRotationsMap() const;
        
        const std::map<int, int>& getParamterIdZShiftsMap() const;
        
        const std::map<int, int>& getParamterIdYShiftsMap() const;
        
        const std::map<int, int>& getParamterIdXShiftsMap() const;
        
        void setMEstimatorType( const std::string& _mEstimatorType );
        
        std::string getMEstimatorType() const;
        
        void setExcludeFromFitPlanes(const std::vector<int>&);
        
        std::vector<int> getExcludeFromFitPlanes() const;

    private:

        TMatrixD PropagatePar(  double, double, double, double, double, double, double );

	TVectorD getXYZfromDzNum( double, double, double, double, double, double, double ) const;
        
        TMatrixD propagatePar(double);

        void addMeasurementsGBL( gbl::GblPoint&, TVectorD&, TVectorD&, const double*, const double*, const EVENT::FloatVec&, TMatrixD& );
        
        void addGlobalParametersGBL( gbl::GblPoint&, TMatrixD&, std::vector<int>&, int, const double*, double, double );
        
        void pushBackPoint( std::vector< gbl::GblPoint >&, const gbl::GblPoint&, int );

        void pushBackPointMille( std::vector< gbl::GblPoint >&, const gbl::GblPoint&, int );
 
        void prepareLCIOTrack( gbl::GblTrajectory*, vector<const IMPL::TrackImpl*>::const_iterator&,
                                double, int); //, double, double, double, double, double );

        void prepareMilleOut( gbl::GblTrajectory* );


    private:
        vector<const IMPL::TrackImpl*> _trackCandidatesVec;


				double _omegaCorrections;	
				double _intersectionLocalXZCorrections;	
				double _intersectionLocalYZCorrections;	
				double _localPosXCorrections;
				double _localPosYCorrections;

    private:
        /** Beam charge in [e] */
        double _beamQ;

        /** Beam energy in [GeV] */
        double _eBeam;

    private:
        /** Alignment degrees of freedom */

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
        
        /** Parameter ids */
        std::map<int,int> _parameterIdZShiftsMap;
        
        /** Parameter ids */
        std::map<int,int> _parameterIdXRotationsMap;
        
        /** Parameter ids */
        std::map<int,int> _parameterIdYRotationsMap;
        
        /** Parameter ids */
        std::map<int,int> _parameterIdZRotationsMap;
				//this maps the original states from the pattern recognition with true and the other scatterers to take into account volumes inbetween them as false. This is needed since you may not know the number of scatterers inbetween in future versions. 

//Need to create a comparision function to order the points. Needs to be in class format. not sure why?
struct compare_points
{
  bool operator()(EUTelTrackStateImpl a, EUTelTrackStateImpl b) const
  {
 			streamlog_out(DEBUG0) << "Label 1 " << a.getLocation() <<" Label 2 " << b.getLocation() <<std::endl;
 			return a.getLocation() > b.getLocation();
  }
};
				std::map<EUTelTrackStateImpl, gbl::GblPoint*, compare_points > _PointToState;

void OutputMap(std::map<EUTelTrackStateImpl,gbl::GblPoint*, compare_points > _PointToState);

std::vector< pair< EUTelState, int> > _vectorOfPairsMeasurementStatesAndLabels;//This is used within alignment since you need to associate MEASUREMENT states to  labels
std::vector< pair< EUTelState, int> > _vectorOfPairsStatesAndLabels;//This is used in track fit since you want to associate ANY states to labels.
std::vector<gbl::GblPoint> _points;

unsigned int _counter_num_pointer = 1;
       
        /** Planes ids to be excluded from refit */
        std::vector< int > _excludeFromFit;
        
        // Track requirements for alignment step
        
        /** Maximum allowed track 2hi2 value*/
        double _chi2cut;

	/** ODE integrator for equations of motion */
        EUTelUtilityRungeKutta* _eomIntegrator;

	/** ODE for equations of motion */
        ODE* _eomODE;

				EUTelMillepede* _MilleInterface;
        
    };

}
#endif	/* EUTELGBLFITTER_H */


#endif
