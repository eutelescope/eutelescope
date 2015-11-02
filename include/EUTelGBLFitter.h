/**  EUTelGBLFitter
 * 
 *  This is the link between GBL and EUTelescope. 
 *  This is all the functions to create and analysis GBL trajectories directly.
 *  Everything you need to know about the GBL functions can be found here: http://www.desy.de/~kleinwrt/GBL/doc/cpp/html/index.html 
 *  The functions here describe how to create a GBL trajectory and then analysis the direct output. 
 *  Further analysis of the EUTelTrack objects, GBL tracks saved in a generic LCIO container.
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
#include "EUTelRadCal.h"
#include "EUTelTrackCreate.h"
#include "EUTelBlock.h"
#include "EUTelExcludedPlanes.h"
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
           void  initNav();
           ///\todo Internal parameterisation will only work for tracks with a hits on 0,2,3,5 
 
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

			void setScattererGBL(gbl::GblPoint& point,EUTelState & state,double & );
            void setScattererMedGBL(gbl::GblPoint& point,TVectorD& kinks, double& varMed );

            /// Will take a track and fill it with information about it's radiation length 
            /// This is split between radiation length of the included planes and radiation length in front of it.
            /// Excluded planes are included in the radiation length in front of the included sensors. 
            /**
             * \param [in] track 
             */

            double getRad(EUTelTrack &, std::map<const int,double>&, std::map<const int,double>&);
            std::vector<double> getZPosScat(EUTelState &);

            /// A simple way to set hits covariance matrix.
			void setMeasurementCov(EUTelState& state);
            inline void setBeamEnergy(double beamE) { this->_eBeam = beamE; }
            inline void setMode(int mode) { this->_mode = mode; }
            inline void setIncMed(int incMed) { this->_incMed = incMed; }

            /// This links the binary created in millepede to GBLFitter.
			void SetMilleBinary(gbl::MilleBinary* mille) { this->_mille = mille; }
            ///Links the EUTelMillepede class to EUTelGBLFitter.
			void setMillepede( EUTelMillepede* mille ) { _MilleInterface =  mille; }
            /// Resolution of the points in the local frame X/Y
			void setParamterIdXResolutionVec( const std::vector<float>& );
			void setParamterIdYResolutionVec( const std::vector<float>& );
            std::vector<double> getWeigMeanVar(double &, double &);


            /// This is the down weighting certain hits will get due to non Gaussian errors.
			void setMEstimatorType( const std::string& );
            ///This function will create the information needed to place the scatterers at the correct location and with the correct variance. 
            /**
             * \param [in] track 
             * \param [out]  scatPos Where we want to place the scatterers after the plane. 
             */

    		std::map< unsigned int, std::vector<double> > getScatPos(EUTelTrack&) const;
            /// This will create the point list which describes each points relation to one another.
            /// At this stage only the points with a local to global transform and the ones which do not must be saved. 
            /// This is done using the position of the scatterers which is calculated internally.
            /// The scatterers are all assumed  
            /**
             * \param [in] track This is the measurement and scattering planes. These must have Local->Global relation. 
             * \param [out] pointList This is the points which describe the EUTelTrack without any measurements (Scattering and residuals) 
             */
            void getBasicList( EUTelTrack & ,std::vector< gbl::GblPoint >& pointList);   
            /// This will create the needed GBL points from EUTelTracks.
            /// This will create the GBL points from the states and the needed scatterering points.  
            /**
             * \param [in]   track EUTelTrack 
             * \param [out]  pointList vector of GBL points
             */

			void getGBLPointsFromTrack(EUTelTrack&, std::vector< gbl::GblPoint >&);
            /// Add measurement to the trajectory. Measurement in local frame.
            /// Note local frame makes the description of the error easy to describe externally.  
            /**
             * \param [in]   track EUTelTrack 
             * \param [in]  pointList vector of GBL points
             * \param [in]  linkGL Link between points with local->global transform and GBL point 
             * \param [out]  linkMeas Link between states with a measurement for position and GBL point.
             */

            void getMeas(EUTelTrack&, std::vector< gbl::GblPoint >&);
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
             */

            void getScat(EUTelTrack&, std::vector< gbl::GblPoint >&);
            /// Get the corrections for the GBL trajectories for states. Also update track automatically. 
            ///   
            /**
             * \param [in]   GBLtraj GBL trajectory 
             * \param [in]   EUTelTrack track to update.
             * \param [out]  corrections This is a map from sensor ID to vector. 
             */
			void getCorr(gbl::GblTrajectory* , EUTelTrack& ,std::map<int,std::vector<double> >& );
            /// Get the corrections for the GBL trajectories for states. Also update track automatically. 
            ///   
            /**
             * \param [in]   GBLtraj GBL trajectory 
             * \param [in]   track EUTeltrack so we know which points are from the planes
             * \param [in]   pointList used to check that this list is correct. 
             * \param [out]  residuals Residuals error associated to the correct sensor ID
             * \param [out]  error this is the residual error associated to that sensorID. 
             */
			void getResLoc(gbl::GblTrajectory*, EUTelTrack &, std::vector< gbl::GblPoint >, std::map< int, std::map< float, float > > & , std::map< int, std::map< float, float > >&);
             /// This is the function which links the GBL track to millepede.cc object. 
             /// Millepede has global labels internal. Jacobain is then calcualted and attached to point with labels . 
            /**
             * \param [in]   vector GBLpoints  
             * \param [in]   EUTelTrack track to update.
             */

			void getGloPar(std::vector< gbl::GblPoint >& , EUTelTrack&);
            void getLocalKink(EUTelTrack& track, std::vector< gbl::GblPoint >& pointList);

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
            void setArcLengths(EUTelTrack & track);

    protected:
            int _numberRadLoss;
			double _eBeam;
            int _mode;
            int _incMed;
            double _dutDistCut;
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
