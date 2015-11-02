#ifndef EUTELTRACKCREATE_H
#define EUTELTRACKCREATE_H

#include "EUTelGeometryTelescopeGeoDescription.h"
#include "TMatrix.h"
#include "TVector3.h"
#include "gear/BField.h"
#include "EUTelTrack.h"
#include "EUTelState.h"
#include "EUTelHit.h"
#include "EUTelNav.h"
#include "EUTelExcludedPlanes.h"


namespace eutelescope 
{

class EUTelTrackCreate
{
	public: 
    //! Create track with parameterisation.  
    /*! It will create a track from the input and add hits and states as required by excluded planes. 
     *  Offset and slope is defined from the hits at either end of the track. 
     *  Slope for curved(With magnetic field) tracks is defined from the central point between the hits. The slope change is taken into account from that point.
     *  
     *  qOverP is the curvature term. Note this is dependent on (cZxB) BFac to get the actual curvature term. However BFac is a constant for each track so qOverP
     *  is the variable of interest. This is calculated from the slope differences after the initial guessed curvature is deducted. 
     *  The slope differences left must be due to more/less energetic particles.
     *
     * \param [in] hits any number of hits.
     * \param [in] offset this is the (x,y,z) of first hit and z of last. 
     * \param [in] trackSlope slope at centre point of offset.
     * \param [in] qOverP This is the charge of the particle (-1 => electron) over the energy (GeV)
     */

    static EUTelTrack getTrack(std::vector<EUTelHit> hits, std::vector<double> offset, std::vector<double> trackSlope,double qOverP );

    ///  This function will take a vector with the planes 0,2,3,5 included as minimum. 
    /// These four hits are then used to determine correction to curvature and parameterisation of the track.
    /// HIT ORDER DOES NOT MATTER INTO THIS FUNCTION
    /**
     * \param [in] hits vector of 4 hits which form track a known track.
     * \return track EUTelTrack
     */

    static EUTelTrack getTrackFourHits(std::vector<EUTelHit> hits);

};

    ///\todo Must initialise EUTelNav to get the beam energy. Should initialise classes together in a more clever way.
    /// Exclude sensors class is also used here. Should really make sure this has been initialised. What is the best way to do this?

}
#endif
