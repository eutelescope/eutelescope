/**  EUTelGBLFitter
 * 
 *  This is the link between GBL and EUTelescope. 
 *  Everything you need to know about GBL can be found here. 
 *  An example of the use of this processor is in EUTelProcessorGBLFitter and EUTelProcessorGBLAlign.
 *  contact:alexander.morton975@gmail.com 
 */

#ifdef USE_GBL

#ifndef EUTELGBLFITTER_H
#define	EUTELGBLFITTER_H

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
#include "TMatrixD.h"

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

	class EUTelGBLFitter{
        
    private:
			DISALLOW_COPY_AND_ASSIGN(EUTelGBLFitter)        
      
    public:
			EUTelGBLFitter();
			~EUTelGBLFitter();
            /// This function will add the measurement to the GBL point. 
            /// Internally a projection matrix is used to link the global frame to the local one. 
            /**
             * \param [in] point GBL point which the measurement will be added to 
             * \param [in]  State  This is the EUTelState with measurement information
             */

			void setMeasurementGBL(gbl::GblPoint& point,EUTelState& );
            /// This function will add the scattering to the GBL point. 
            /// Internally a projection matrix is used to link the track frame to the local/global one 
            /**
             * \param [in] point GBL point which the measurement will be added to 
             * \param [in]  State  This is the EUTelState with measurement information
             */

			void setScattererGBL(gbl::GblPoint& point,EUTelState & state );
            /// Will fill track with information about the sensor and air radiation length 
            /**
             * \param [in] track 
             * \param [in] mapSensor Link between sensorID and radiation length of sensor. 
             * \param [in] mapAir Link between sensorID and radiation length of air after sensor.
             */
            void setRadLengths(EUTelTrack & track,std::map<const int,double>&  mapSensor, std::map<const int ,double>&  mapAir, double & rad );
            /// Will take a track and fill it with information about it's radiation length 
            /**
             * \param [in] track 
             */

            void setRad(EUTelTrack &);
            /// A simple way to set hits covariance matrix.
			void setMeasurementCov(EUTelState& state);
			inline void setBeamCharge(double beamQ) { this->_beamQ = beamQ; }
            inline void setBeamEnergy(double beamE) { this->_eBeam = beamE; }
            /// This links the binary created in millepede to GBLFitter.
			void SetMilleBinary(gbl::MilleBinary* _mille) { this->_mille = _mille; }
            ///Links the EUTelMillepede class to EUTelGBLFitter.
			void setMillepede( EUTelMillepede* Mille ) { _MilleInterface =  Mille; }
            /// Resolution of the points in the local frame X/Y
			void setParamterIdXResolutionVec( const std::vector<float>& );
			void setParamterIdYResolutionVec( const std::vector<float>& );
            /// This is the down weighting certain hits will get due to non Gaussian errors.
			void setMEstimatorType( const std::string& );
			/// GET
            ///This function will create the information needed to place the scatterers at the correct location and with the correct variance. 
            /**
             * \param [in] track 
             * \param [out]  scatPos Where we want to place the scatterers after the plane. 
             */

    		std::map< unsigned int, std::vector<double> > getScatPos(EUTelTrack&) const;
            /// This will create the point list which describes each points relation to one another.
            /// At this stage only the points with a local to global transform and the ones which do not must be saved. 
            /// This is done using the position of the scatterers which is calculated before. 
            /// There is no limit on the number of scatterers here since we can easily estimate the global direction at any position in the z axis.
            /**
             * \param [in] track This is the measurement and scattering planes. These must have Local->Global relation. 
             * \param [in] vectorScatPosAfterPlane Z-positions of scatterers in relation to the scatterer before.
             * \param [out] pointList This is the points which describe the EUTelTrack without any measurements (Scattering and residuals) 
             * \param [out] vecPairsMeasStatesAndLabels This is the link pointGBL.label->EUTelState.location (Denotes that GBL point has local->global transform) 
             */
            void getBasicList( EUTelTrack& ,std::vector< gbl::GblPoint >& pointList, std::map<  unsigned int, unsigned int>  & );   
            /// This will create the needed GBL points from EUTelTracks.
            /// This will create the GBL points from the states and the needed scatterering points.  
            /**
             * \param [in]   track EUTelTrack 
             * \param [out]  pointList vector of GBL points
             * \param [out]  linkGL Link between points with local->global transform and GBL point
             * \param [out]  linkMeas This is link between states with hits and GBL points. 
             */

			void getGBLPointsFromTrack(EUTelTrack&, std::vector< gbl::GblPoint >&, std::map< unsigned int, unsigned int > &, std::map< unsigned int,unsigned int > & );
            /// Add measurement to the trajectory. Measurement in local frame.
            /// Note local frame makes the description of the error easy to describe externally.  
            /**
             * \param [in]   track EUTelTrack 
             * \param [in]  pointList vector of GBL points
             * \param [in]  linkGL Link between points with local->global transform and GBL point 
             * \param [out]  linkMeas Link between states with a measurement for position and GBL point.
             */

            void getMeas(EUTelTrack&, std::vector< gbl::GblPoint >&,  std::map< unsigned int,unsigned int> &, std::map< unsigned int, unsigned int> & );
            /// Add Scattering information to GBL trajectories. Kinks are in local frame. 
            /// All planes with no global to local transforms will be assumed to have normal incidence on a plane. 
            /// Note here that kink angles are always measured in the local frame and therefore the precision matrix is also.
            /// The kink angle precision matrix is non diagonal in the local frame but in the frame in which the measurement of kink is on a plane perpendicular to the track
            /// it is diagonal.
            /// Internally this transform takes place. Note this transform is not between the local and global frame. It is between the frame perpendicular to the track 
            /// and plane of incidence. Transforms between frames will make no difference it is the relative incidence which counts.
            ///   
            /**
             * \param [in]   track EUTelTrack 
             * \param [in]  pointList vector of GBL points
             * \param [in]  linkGL Link between points with local->global transform and GBL point 
             */

            void getScat(EUTelTrack&, std::vector< gbl::GblPoint >&,  std::map< unsigned int,unsigned int> &);
            /// Get the corrections for the GBL trajectories for states. Also update track automatically. 
            ///   
            /**
             * \param [in]   GBLtraj GBL trajectory 
             * \param [in]   EUTelTrack track to update.
             * \param [in]   linkGL link to states we want to update. 
             * \param [out]  corrections This is a map from sensor ID to vector. 
             */
			void getCorr(gbl::GblTrajectory* , EUTelTrack& , std::map<unsigned int ,unsigned int >&  ,std::map<int,std::vector<double> >& );
            /// Get the corrections for the GBL trajectories for states. Also update track automatically. 
            ///   
            /**
             * \param [in]   GBLtraj GBL trajectory 
             * \param [in]   linkGL link to states we want to update.
             * \param [in]   pointList used to check that this list is correct. 
             * \param [out]  residuals Residuals error associated to the correct sensor ID
             * \param [out]  error this is the residual error associated to that sensorID. 
             */
			void getResLoc(gbl::GblTrajectory*, std::map< unsigned int, unsigned int> &, std::vector< gbl::GblPoint >, std::map< int, std::map< float, float > > & , std::map< int, std::map< float, float > >&);
             /// This is the function which links the GBL track to millepede.cc object. 
             /// Millepede has global labels internal. Jacobain is then calcualted and attached to point with labels . 
            /**
             * \param [in]   vector GBLpoints  
             * \param [in]   EUTelTrack track to update.
             * \param [in]   linkGl map from EUTelStates with hits to GBL label.
             */

			void getGloPar(std::vector< gbl::GblPoint >& , EUTelTrack&, std::map< unsigned int, unsigned int>  & );
			inline double getBeamEnergy() const { return _eBeam; }
			std::string getMEstimatorType() const;
			///TEST
			void testUserInput();
			void testTrack(EUTelTrack& track);
			///COMPUTE
            /// This will fit the trajectory and output the chi2 and ndf 
            ///   
            /**
             * \param [in]   GBLtraj GBL trajectory 
             * \param [out]  chi2 
             * \param [out]  ndf 
             * \param [out]  error code that tells you if there was a problem from GBL.
             */

			void computeTrajectoryAndFit(gbl::GblTrajectory* traj, double* chi2, int* ndf, int & ierr);
    protected:
            int _numberRadLoss;
			double _beamQ;
			double _eBeam;
			/** Outlier downweighting option */
			std::string _mEstimatorType;
			/** Milipede binary file handle */
			gbl::MilleBinary* _mille;
			std::string _binaryname;
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
			EUTelMillepede* _MilleInterface;
        
    };
}
#endif	/* EUTELGBLFITTER_H */
#endif
