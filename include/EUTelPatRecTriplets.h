/**  EUTelPaRecTriplet
 *  This class performs pattern recognition in two ways. First triplets are formed from the full telescope system and then two hits on each track are used to propagate to the    *  DUT plane. Which two hit on the track can be selected in the steering file.  
 *  This way will have a very low fake rate which could hurt your statistics. 
 *
 *  Therefore a less stringent pattern recognition can be used which two planes' hits are searched. Tracks are parameterised and if there are enough hit then this is taken as a  *  final track. This is a simple straight line fit for non magnetic field situations. 
 *
 *  States are create on planes if included in fit. DUT hits are excluded if outside some cut range. 
 *  contact:alexander.morton975@gmail.com 
 */

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

	class EUTelPatRecTriplets {
    ///Doublet. This is a relation between the two outer planes of each arm of the telescope. 
    struct doublets {
        std::vector<float> pos;
        std::vector<double> slope;
        std::vector<double> diff;
        std::vector<EUTelHit> hits;

    }; 
    ///Triplet. This is the relationship between three hits on each arm. 
    struct triplets {
        unsigned int cenPlane;
        unsigned int matches;
        unsigned int fitID;

        std::vector<double> pos;
        std::vector<double> slope;
        std::vector<double> diff;
        std::vector<EUTelHit> hits
    }; 

	private:
		DISALLOW_COPY_AND_ASSIGN(EUTelPatRecTriplets) // prevent users from making (default) copies of processors
	public:
		EUTelPatRecTriplets();
		~EUTelPatRecTriplets();
        /// Create vector internally with ID and hits vector. 
        ///  
        /**
         * \param [in] allHitsVec a simple vector of hits 
         */

		void setHitsVec(EVENT::TrackerHitVec& allHitsVec){
			_allHitsVec = allHitsVec;
		}
        /// Create vector of triplets from hits passed to fitter. 
        /// Hits are linked by removing the difference in position relative to curvature for doublet and triplet creation.
        /// This is done assuming the particle begins to deflect out of the detector system. This initial start of the magnetic field is set with initial distance
        /// in geometry. This is done in EUTelProcessorPatRecTriplets as an example. 
        /// Cuts are performed on the distance between hits on planes. See GBL examples for more information. 
        /// 
        ///  
        /**
         * \return a vector of triplets which passed the cuts.
         */

        std::vector<EUTelPatRecTriplets::triplets> getTriplets();
        /// The triplets on each plane are passed. 
        /// If they pass a extraplolated postion/slope comparison then form a track. 
        /// If there are more than 1 match for a triplet remove triplet.  
        /// No curvature information is used at this here. This was all input in the creation of the triplets. 
        ///  
        /**
         * \return a vector of tracks. These tracks are ready to use for GBL fitting
         */
        std::map<int,std::vector<EUTelHit> >  getTrackHitsFromTriplets(std::vector<EUTelPatRecTriplets::triplets>&);

        std::vector<EUTelTrack> getTracks();
        ///Calculates the initial curvature fo the tracks. This is then passed to navigation.
        std::vector<double>  getCurvXY();
        ///cZxB=>This is used in the calculation of curvature.
        /// This is also needed for the determination of the update to the q/p parameter. All of which is dealt with in navigation (EUTelNav)
        TVector3  getBFac();
        /// 
        /// Will return hits in the correct z order. 
        ///  
        /**
         * \param [in] hits not the correct in Z
         * \return hits The correct order of hits returned. 
         */

        std::vector<EUTelHit> getCorrHitOrder(std::vector<EUTelHit> hits );

        std::vector<EUTelHit> getDUTHitsOrder(EUTelTrack track, std::vector<EUTelHit> dutHit );

        ///  This function will take a vector with the planes 0,2,3,5 included as minimum. 
        /// These four hits are then used to determine correction to curvature and parameterisation of the track.
        /**
         * \param [in] hits vector of 4 hits o track.
         * \return track EUTelTreack object
         */

        EUTelTrack getTrackFourHits(std::vector<EUTelHit> hits);
        ///  The function needs no minimum of hits. It will create a track from the input and add hits and states are required by excluded planes. 
        /**
         * \param [in] hits any number of hits.
         * \param [in] offset this is the (x,y,z) of first hit and z of last. 
         * \param [in] trackSlope slope at centre point of offset.
         * \param [in] curvCorr This is the co
         */

        EUTelTrack getTrack(std::vector<EUTelHit> hits, std::vector<double> offset, std::vector<double> trackSlope,double curvCorr );
        EUTelTrack getTrack(triplets tripLeft,triplets tripRight);
        void getTrackAvePara(EUTelHit& firstHit, EUTelHit& endHit, std::vector<double>& offset, std::vector<double>& trackSlope);
        bool getTriplet(doublets&, EUTelHit &, triplets&  );
        std::vector<float>  getTripPosAtZ(triplets trip, float posZ );
        std::vector<float>  getDoubPosAtZ(doublets doub, float posZ);
        float getDistLocal(int location, std::vector<float> pos, TVector3 hitPosGlo);
        bool getDoubHitOnTraj(const doublets doub, const std::vector<unsigned int> & sen,std::vector<EUTelHit>& newHits   );


        void setPlaneDimensionsVec(EVENT::IntVec& planeDimensions);
        void setPlaneExclude(IntVec& planeIDs);  

		inline int getEventNumber()	const {
			return _eventNumber;
		}

		inline double getBeamMomentum() const {
			return _beamE;
		}
		inline double getBeamCharge() const {
			return _beamQ;
		}

		inline int getNumberDoublets(){
			return _numberDoublets;
		}
		//SETTERS
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

		inline  void setPlanesToCreateSeedsFrom(EVENT::IntVec createSeedsFromPlanes){
			this-> _createSeedsFromPlanes = createSeedsFromPlanes;
		}

		inline void setBeamCharge(double q) {
			this->_beamQ = q;
		}
        void setArcLengths(EUTelTrack & track);

		//TEST
		void testUserInput();
//		void testTrackCandidates();
		
		//OTHER
		void testHitsVecPerPlane();
		void findTracksWithEnoughHits();
	    bool getDoublet( EUTelHit&, EUTelHit&,doublets& );
		void testTrackQuality( std::vector<EUTelTrack>&);
		//VARIABLES
		int _eventNumber;
		int _totalNumberOfHits;
		int _totalNumberOfSharedHits;
		std::map< int, int > _planeDimensions;
		bool _firstExecution;
		EVENT::IntVec _createSeedsFromPlanes;
		EVENT::FloatVec _excludePlanes;         
        std::vector<triplets> _tripletsVec;

		std::vector<EUTelTrack> _tracks;
		std::vector<EUTelTrack> _tracksWithDUTHit;
		std::vector<EUTelTrack> _tracksAfterEnoughHitsCut;
		std::vector<EUTelTrack>	_finalTracks;
        std::map<int,int> _senZOrderToIDWithoutExcPla;

		unsigned int _numberOfTracksTotal;
		unsigned int _numberOfTracksTotalWithDUT;
        unsigned int _tracksWithoutHit;
		unsigned int _numberTripletsLeft;
		unsigned int _numberTripletsRight;
		unsigned int _numberDoublets;
		void printHits();


        ///Member variables. public for now.
        std::vector<float> _doubletDistCut;
        std::vector<float> _doubletCenDistCut;
        std::vector<float> _tripletConnectDistCut;
        std::vector<float> _tripletSlopeCuts;

        std::vector<int> _senNotExc;
private:


		/** Update track state and it's cov matrix */
	
		/** Update track propagation matrix for a given step */
		
		/** Update Kalman gain matrix */

		/** Propagate track state */
		void propagateTrackState( EUTelState& );

		/** Sort hits according to particles propagation direction */
		bool sortHitsByMeasurementLayers( const EVENT::TrackerHitVec& );
		
		// Helper functions
private:
		
		/** Calculate track momentum from track parameters */
		
		/** Calculate position of the track in global 
		 * coordinate system for given arc length starting
		 * from track's ref. point*/

		void setNewState(float position[],float momentum[],  EUTelState& newState);
       /// This function will take tracks produced by findTrackFromTriplets and associate DUT tracks to the closest track to it.
       /// Each track will only add one DUT hit per track. Need to do track comparison afterwards. 
       /// Association is done in the local frame. So the track prediction is determined in the global frame and then transformed to the DUT local frame.
       /// After this the the distance between the hit and tracks is determined (1 or 2 directions dependent on strip/pixel) with the closest track being taken as correct.
        /**
         * \param [in] tracks This is the tracks to look through when the you associate a DUT hit to the closest track
         * \param [return] tracksDUT This is the full EUTelTrack with the hit attached.  
         */

        std::vector<EUTelTrack>	getDUTHit( std::vector<EUTelTrack> &);

		void setRadLengths(EUTelTrack & track,std::map<const int,double>&  mapSensor, std::map<const int ,double>&  mapAir, double & rad );


		
		/** Calculate position of the track in global 
		 * coordinate system for given arc length starting
		 * from track's ref. point */

		
		/** Find hit closest to the track */
		std::map<int ,std::vector<EUTelHit>> _mapHitsVecPerPlane;
	protected:
		EVENT::TrackerHitVec _allHitsVec;//This is all the hits for a single event. 
private:       
		
		/** Final set of tracks in cartesian parameterisation */
		std::map<int, std::vector<EUTelState> > _mapSensorIDToSeedStatesVec;

		// User supplied configuration of the fitter
private:
		/** Maximum number of sensitive planes to be considered for initial seed hits */
		int _planesForPR;

		/** Maximum number of missing on a track candidate */
		int _allowedMissingHits;

		/** Maximum number of track candidates to be stored */
		int _maxTrackCandidates;

		/** Maximal distance in XY plane */
		double _residualsRMax;
		
		/** Beam momentum [GeV/c] */
		double _beamE;
		
		/** Signed beam charge [e] */
		double _beamQ;

		/** Beam energy spread [%] */
		double _beamEnergyUncertainty;
	
		/** Beam angular spread (horizontal,vertical) [mr] */
		EVENT::FloatVec _beamAngularSpread;

		
private:
/** Track parameters propagation jacobian matrix */
TMatrixD _jacobianF;
		
		/** Track parameters covariance C(k,k-1) matrix */
		TMatrixD _trkParamCovCkkm1;
		
		/** Process noise matrix */
		TMatrixDSym _processNoiseQ;       
		
		/** Kalman residual covariance matrix */
		TMatrixD _residualCovR;
		
private:
	};

} // namespace eutelescope

#endif	/* EUTelPatRecTriplets_H */

