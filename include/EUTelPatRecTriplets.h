#ifndef EUTelPatRecTriplets_H
#define	EUTelPatRecTriplets_H

// ROOT
#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"

// system includes <>
#include <iostream>
#include <functional>
#include <algorithm>
#include <string>
#include <vector>
#include <cmath>

// EUTELESCOPE
#include "EUTelUtility.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelTrack.h"
#include "EUTelState.h"
#include "EUTelNav.h"
//LCIO
#include "lcio.h"
#include "IMPL/TrackerHitImpl.h"
#include "IMPL/TrackImpl.h"
#include <UTIL/LCTOOLS.h>

//other
#include "streamlog/streamlog.h"
#include "gear/gearimpl/Vector3D.h"

namespace eutelescope {
  //! Triplet/doublet finder class. 
  /*! The class performed pattern recognition using triplet and doublets constructed from mimosa hits 
   *  Triplets are a collection of three hits and doublets two on different planes. 
   *  These hits are associated together using geometric distances between hits. 
   *  There are two modes:Alignment and high stats. 
   *
   *  The alignment mode will construct two triplets on each arm. Attach DUT hits then return the track.
   *  This will ensure that the tracks passed to alignment are of good quality but the number of tracks found will be smaller than expected.
   *
   *  High stats mode will construct doublets from two planes and then interpolate between these planes to search for new hits.
   *  This mode will produce track of a lower quality but will find more tracks.   
   *
   *  This form of pattern recognition is the limit of what can be done without use of measurement errors. 
   *  It has been test in high occupancy situation (SLAC testbeam) and performs well. 
   *  For more complex situations a Kalman filter of some form or other can be used. 
   *
   */ 

	class EUTelPatRecTriplets {
    public: 
    EUTelPatRecTriplets();
    ~EUTelPatRecTriplets();

    //! Struct object 
    /*! This contains the informtion needed to construct a track from two hits. 
     */

    struct doublets {
        std::vector<float> pos;
        std::vector<double> slope;
        std::vector<double> diff;
        std::vector<EUTelHit> hits;

    }; 
    //! Struct object 
    /*! This contains the informtion needed to construct a track from two hits. 
     *  It also contains information in associating triplet together
     */

    struct triplets {
        /// The central plane of the triplet (1 or 4) is saved for identification
        unsigned int cenPlane;
        ///How many times this triplet has been linked to another on the other arm
        unsigned int matches;
        ///Just an ID given to link two triplets together. 
        unsigned int fitID;

        std::vector<double> pos;
        std::vector<double> slope;
        std::vector<double> diff;
        std::vector<EUTelHit> hits;
    }; 
    //! Set the planes which should be excluded from the search.  
    /*! 
     *  Only tracks without a hit on this plane will be created. 
     *  Even if a hit is within the correct range of the cuts the hit is still not included.
     *  @param[in] planeIDs Planes to exclude. 
     */
    ///\todo Could infer this from resolution. 

    void setPlaneExclude(IntVec& planeIDs);  
    //! Set in the correct z order if we are dealing with a strip sensor or pixel sensor.  
    /*! 
     *  A strip sensor is 1 and pixel 2. See GBL example for use in EUTelescope.
     *  @param[in] planeDimenstions This is a vector of 1s and 2s. 1=>strip, 2=>pixel. 
     */

    void setPlaneDimensionsVec(EVENT::IntVec& planeDimensions);
    inline int getEventNumber()	const {
        return _eventNumber;
    }
    //! Create vector internally with ID and hits vector. 
    /*!  
     * @param[in] allHitsVec a simple vector of hits 
     */

    inline void setHitsVec(EVENT::TrackerHitVec& allHitsVec){
        _allHitsVec = allHitsVec;
    }
    void setHitsVecPerPlane();

    void setEventNumber(int eventNumber){
        _eventNumber = eventNumber;
    }
    inline void setTripletSlopeCuts(std::vector<float> cuts ){
        this->_tripletSlopeCuts = cuts;
    };

    inline void setDoubletDistCut(std::vector<float> cuts) {
        this->_doubletDistCut = cuts;
    }

    inline void setTripletConnectDistCut(std::vector<float> cuts) {
        this->_tripletConnectDistCut = cuts;
    }
    inline void setDoubletCenDistCut(std::vector<float> cuts) {
        this->_doubletCenDistCut = cuts;
    }

    inline void setBeamMomentum(double beam) {
        this->_beamE = beam;
    }

    //get public
    
    //! The function should be called in your processor to return EUTelTracks 
    /*! This function will do all the work for you and return the track reconstructed. 
     *  These are ready for use in GBL and can be saved to LCIO with EUTelReaderGenericLCIO. 
     *  See EUTelProcessorPatRecTriplets.cpp for details on how to construct and store tracks.
     *  @return tracks A vector of tracks. These tracks are ready to use for GBL fitting
     */

    std::vector<EUTelTrack> getTracks();
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

    EUTelTrack getTrack(std::vector<EUTelHit> hits, std::vector<double> offset, std::vector<double> trackSlope,double qOverP );

    private:
    //! This will place DUT hit in track. 
    /*! 
     *  @param[in] tripLeft Triplet from upstream planes. 
     *  @param[in] tripRight Triplet from downstream planes
     *  @return track EUTelTrack from triplets 
     */

    EUTelTrack getTrack(triplets tripLeft,triplets tripRight);
    //! Will create triplet from doublet and hit if passes cut.  
    /*! 
     *   Hit must come from plane in the middle of the doublet  
     *   Cut made on distance of central hit to doublet prediction
     *
     *  @param[in] doublet doublet to interpolate between hits 
     *  @param[in] hit possible hit to create triplet from.
     *  @return pass Has the cut been passed.
     */

    bool getTriplet(doublets&, EUTelHit &, triplets&  );
    //! Will create doublet from hits if passes distance cut.  
    /*! 
     * The doublet included all the information need to parameterise the track fro mthe two hits.
     *  @param[in] hit Hit on plane down stream 
     *  @param[in] hit Hit on plane upstream.
     *  @param[out] doublet doublet constructed from hits.
     *  @return pass Has the cut been passed.
     */

    bool getDoublet( EUTelHit&, EUTelHit&,doublets& );
    //! Predict global position using triplet at certain z position.  
    /*! 
     *  
     *  @param[in] trip Triplet to extrapolate from 
     *  @return posZ Position to extrapolate to.
     */

    std::vector<float>  getTripPosAtZ(triplets trip, float posZ );
    //! Predict global position using doublet at certain z position.  
    /*! 
     *  
     *  @param[in] doub Doublet to extrapolate from 
     *  @return posZ Position to extrapolate to.
     */
    std::vector<float>  getDoubPosAtZ(doublets doub, float posZ);

    //Get functions.

    //! Create vector of triplets from hits passed to fitter. 
    /*! Hits are linked by removing the difference in position relative to curvature for doublet and triplet creation.
     *  Cuts are performed on the distance between hits on planes. See GBL examples for more information. 
     */

    std::vector<EUTelPatRecTriplets::triplets> getTriplets();
    //! The triplets on each plane are passed. 
    /*! If they pass a extraplolated postion/slope comparison then form a track. 
     *  If there is more than 1 match for a triplet remove triplet.  
     *  @param[in] 
     *  @return map Map of triplet fit ID and vector of hits.
     */
    std::map<int,std::vector<EUTelHit> >  getTrackHitsFromTriplets(std::vector<EUTelPatRecTriplets::triplets>&);
    /// Will return hits in the correct z order. 
    /**
     * \param[in] hits Not correct in Z
     * \return -cZxB The correct order of hits returned. 
     */
    std::vector<EUTelHit> getCorrHitOrder(std::vector<EUTelHit> hits );
    //! This will place DUT hit in track. 
    /*! 
     *  @param[in] track track without DUT hits. 
     *  @param[in] hits DUT hits to add to hit.
     *  @return hits Hits in the correct order with DUT hits. 
     */

    std::vector<EUTelHit> getDUTHitsOrder(EUTelTrack track, std::vector<EUTelHit> dutHit );
    ///  This function will take a vector with the planes 0,2,3,5 included as minimum. 
    /// These four hits are then used to determine correction to curvature and parameterisation of the track.
    /// HIT ORDER DOES NOT MATTER INTO THIS FUNCTION
    /**
     * \param [in] hits vector of 4 hits which form track a known track.
     * \return track EUTelTrack
     */

    EUTelTrack getTrackFourHits(std::vector<EUTelHit> hits);
    //! Get all hits which the predicted doublet passes  
    /*! 
     *  A predicted track from the doublet is constructed. At each sensor the closest hit is taken if included in list in arguments.
     *  A cut is made on the number of hits in the final track. If less than this the bool returned in false.
     *  @param[in] doub Doublet to parameterise track from.
     *  @param[in] sen list of sensor to collect hits from.
     *  @param[out] newHits list of hits found. WILL NOT ATTACH ORIGINAL DOUBLET HITS
     *  @return pass If we have enough hits?
     */
    bool getDoubHitOnTraj(doublets& doub, std::vector<unsigned int> & sen,std::vector<EUTelHit>& newHits   );


    //! Get all hits which the predicted doublet passes  
    /*! 
     *  @param[in] itHit Iterator of hit object 
     *  @param[in] pos Position of prediction. 
     *  @return  dist This is the distance between the prediction and measurement in the local frame.
     */

    float getDistLocal(std::vector<EUTelHit>::iterator itHit, std::vector<float>& pos);
    inline double getBeamMomentum() const {
        return _beamE;
    }
    DISALLOW_COPY_AND_ASSIGN(EUTelPatRecTriplets) // prevent users from making (default) copies of processors
	public:
    //OTHER
    void printTrackQuality( std::vector<EUTelTrack>&);
    void printHits();
    void testHitsVecPerPlane();
    void testUserInput();
    protected:
    //VARIABLES
    int _eventNumber;
    int _totalNumberOfHits;
    std::map< int, int > _planeDimensions;
    std::map<int,int> _senZOrderToIDWithoutExcPla;
    unsigned int _numberTripletsLeft;
    unsigned int _numberTripletsRight;
    ///Member variables. public for now.
    std::vector<float> _doubletDistCut;
    std::vector<float> _doubletCenDistCut;
    std::vector<float> _tripletConnectDistCut;
    std::vector<float> _tripletSlopeCuts;
    std::map<int ,std::vector<EUTelHit> > _mapHitsVecPerPlane;
    ///This will store all the hits for a single event.
    EVENT::TrackerHitVec _allHitsVec; 
    /** Beam momentum [GeV/c] */
    double _beamE;
    public:
    unsigned int _numberOfTracksTotal;

		
	};

} // namespace eutelescope

#endif	/* EUTelPatRecTriplets_H */

