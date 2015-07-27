/** Basic navigation tools for general track fitting.  
 * TO DO: Need to change these to Eigen. Will change to TMatrix in GBL
 *
 *  The naviagation tools use an initial bfactor (ZxB) and energy (q/p). This is used to calculate the initial curvature.
 *  Hit information is added to update the q/p and output to parameterise the track. 
 *  Track parameterisation is also done on a state level also. 
 *  Jacobians exist for global to global linking of states (position and incidence on planes). This works for an homogeneous magnetic field. 
 *  There are also jacobains for linking local to local frames. This is not used in the GBL tracking where the linking from the global to local (Measurement) systems is
 *  done via a propagator dependent on the state incidence on the track. This propagator is part of the EUTelState.
 *  State defined as (q/p,slopeX,slopeY,x,y)
 *
 *  contact:alexander.morton975@gmail.com 
 */
///TO DO: Update to eigen and create TMatrix objects from utility.
#ifndef EUTELNAV_H
#define EUTELNAV_H

#include "EUTelGeometryTelescopeGeoDescription.h"
#include "TMatrix.h"
#include "TVector3.h"
#include "gear/BField.h"
#include "EUTelTrack.h"
#include "EUTelState.h"
#include "EUTelHit.h"


namespace eutelescope 
{

class EUTelNav
{
	public: 
        EUTelNav();
        ~EUTelNav();
        static void init( double beamEnergy);

        /// This function will will produce a 5x5 matrix to propagate a local state to global. Incidence information is propagated here too.
        /// To define this tranform the incidence in the global frame and rotation matrix local->global is passed. 
        /**
         * \param [in] tw1 TVector3  in the global frame. 
         * \param [in]  TRotMatrix Local->Global 
         * \return localToGlobal transforms the local state vector to the global state vector.  
         */

		static TMatrixD getMeasToGlobal(TVector3 t1w, TMatrixD TRotMatrix);
        /// This links states in the global frame. 
        /// The incidence and arc length to the next state to link to are needed. This works with/without magnetic field.
        /// Small bending limit is assumed in derivation q/p->0 
        /// 
        /**
         * \param [in] ds Arclength to next state. 
         * \param [in]  t1w direction in global frame.  
         * \return stateToState This jacobian which links two global states.  
         */

		static TMatrixD getPropagationJacobianGlobalToGlobal(float ds, TVector3 t1w);
        //! Parameterise track from two hits 
        /*! This will get the needed information to describe a track from two hits.
         *  Need 4 hits to determine qOverP change.  
         *
         *  @param[in] firstHit Plane upstream from endHit 
         *  @param[in] endHit last hit on trajectory. 
         *  @param[out] offset This is x,y,z of first hit and z of endHit. So (sx,sy,sz,ez) were s=>start and e => end.
         *  @param[out] trackSlope This is the track slope in X/Y.
         */

        static void getTrackAvePara(EUTelHit &, EUTelHit &, std::vector<double>& offset, std::vector<double>& trackSlope);
        /// This is the update to the track parameter which corresponds to curvature.
        /// Need 4 hits as a minimum to calculate this. 
        /// Hits are all in increasing z order!
        /// 
        /**
         * \param [in] hit1  
         * \param [in] hit2 
         * \param [in] hit3
         * \param [in] hit4
         * \return corr correction must be added to q/p. 
         */
        static void getTrackPredictionFromParam( std::vector<double> const & offset, std::vector<double> const & trackSlope,double const & qOverP, double const & posZ, Eigen::Vector3d & posPred, std::vector<double>& slopePred);

        static double getCorr(EUTelHit &, EUTelHit &, EUTelHit &, EUTelHit &);
        static std::vector<double>  getCurvXY();
        static TVector3 getBFac();
        static TVector3 _bFac;
        static std::vector<double> _curv;
        static double _intBeamE;
    ///\todo Need to initialised this in silly way. Problem is beam energy is passed to each processor. Should be the same for all

};


}
#endif
