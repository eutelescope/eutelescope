//***************************************************
//                  EUTelPatRecTriplets            //
//This class create tracks using triplets formed   //
//on the telescope planes. DUT hits are then added //
//using the predicted track from the mimosas       //                                       
//***********************************************  //

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

    struct doublets {
        std::vector<float> pos;
        std::vector<double> slope;
        std::vector<double> diff;
    }; 
    struct triplets {
        unsigned int cenPlane;
        unsigned int matches;
        unsigned int fitID;

        std::vector<double> pos;
        std::vector<double> slope;
        std::vector<double> diff;
        std::vector<EUTelState> states;
    }; 

	private:
		DISALLOW_COPY_AND_ASSIGN(EUTelPatRecTriplets) // prevent users from making (default) copies of processors
	public:
		EUTelPatRecTriplets();
		~EUTelPatRecTriplets();

        //doublet distance cut
        std::vector<float> _doubletDistCut;
        std::vector<float> _doubletCenDistCut;
        std::vector<float> _tripletConnectDistCut;

        std::vector<EUTelPatRecTriplets::triplets> getTriplets();
        std::vector<EUTelTrack> getTracks( );
        std::vector<double>  getCurvXY();
        TVector3  getBFac();
        std::vector<EUTelHit> getDUTHitsOrder(EUTelTrack track, std::vector<EUTelHit> dutHit );


        EUTelTrack getTrackDUTHit(std::vector<EUTelTrack>::iterator itTrack, EUTelState stateDUT );
        EUTelTrack getTrack(std::vector<EUTelHit> hits);
        EUTelTrack getTrack(std::vector<EUTelHit> hits, std::vector<double> offset, std::vector<double> trackSlope,std::vector<double> curvCorr );
        EUTelTrack getTrack(triplets tripLeft,triplets tripRight);
        EUTelTrack getTrack(triplets tripLeft,triplets tripRight,EUTelState stateDUT );
        void getTrackAvePara(EUTelHit& firstHit, EUTelHit& endHit, std::vector<double>& offset, std::vector<double>& trackSlope);
        triplets getTriplet(EUTelState & left, EUTelState & cen,EUTelState & right, doublets& doublet );
        std::vector<double> getCorr(EUTelHit & hitArmOne1, EUTelHit & hitArmOne2, EUTelHit & hitArmTwo1, EUTelHit & hitArmTwo2);

        void setPlaneDimensionsVec(EVENT::IntVec& planeDimensions);



		inline int getEventNumber()	const {
			return _eventNumber;
		}

		inline int getAllowedMissingHits() const {
			return _allowedMissingHits;
		}

		inline double getWindowSize() const {
			return _residualsRMax;
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
        std::vector<float>  getTripPosAtZ(triplets trip, float posZ );

		//SETTERS
		void setHitsVecPerPlane();

		void setHitsVec(EVENT::TrackerHitVec& allHitsVec){
			_allHitsVec = allHitsVec;
		}

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


		//Here if the user does not set a create seeds from planes x. The we set it automatically to the first plane travelling as the beam travels. 
		//This has the best of both world. No reduction on functionality. User does not even know this is here. 	
		inline 	void setAutoPlanestoCreateSeedsFrom(){
			if(_createSeedsFromPlanes.size() == 0){
				_createSeedsFromPlanes.push_back(geo::gGeometry().sensorZOrdertoIDs().at(0));
			}
		}	
		//COMPUTE
		TVector3 computeInitialMomentumGlobal();
		//TEST
		void testUserInput();
//		void testTrackCandidates();
		
		//OTHER
		void testHitsVecPerPlane();
        std::vector<EUTelTrack>  findTrackFromTriplets(std::vector<EUTelPatRecTriplets::triplets>&);
		void findTracksWithEnoughHits();
	    doublets getDoublet( double hitLeftPos[3], double hitRightPos[3],double curvX,double curvY );
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

		unsigned int _numberOfTracksTotal;
		unsigned int _numberOfTracksTotalWithDUT;
        unsigned int _tracksWithoutHit;
		unsigned int _numberTripletsLeft;
		unsigned int _numberTripletsRight;
		unsigned int _numberDoublets;
		void printHits();

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
		std::map<int ,EVENT::TrackerHitVec> _mapHitsVecPerPlane;
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

        std::vector<float> _tripletSlopeCuts;

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

