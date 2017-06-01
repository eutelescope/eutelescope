// Version: $Id$
/**  EUTelProcessorPatRecTriplets
 * 
 *  If compiled with MARLIN_USE_AIDA 
 *  This will do pattern recognition using a collection of LCIO::TRACKERHIT objects.
 *  Triplets are formed from each arm of the telescope using some cuts DEPENDENT on geometry.
 *  An initial track prediction is calculated using the 6 hits deemed a track.
 *  Each DUT hit created is then associated to the closest track for each DUT. No cut applyed here.
 *  MUST EXCLUDE A SENSOR IF ANOTHER IS SIDE BY SIDE WITH IT. 
 * 
 */

#ifndef EUTelProcessorPatRecTriplets_h
#define EUTelProcessorPatRecTriplets_h 1

// AIDA
#include <AIDA/IBaseHistogram.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/IHistogramFactory.h>

// C++
#include <map>
#include <memory>
#include <string>
#include <vector>

// LCIO
#include <EVENT/LCCollection.h>
#include <EVENT/TrackerHit.h>
#include <IMPL/TrackImpl.h>
#include "IMPL/TrackerHitImpl.h"
#include "lcio.h"
#include "marlin/Processor.h"

// MARLIN
#include "marlin/Exceptions.h"
#include "marlin/Global.h"
#include "marlin/Processor.h"
#include "marlin/VerbosityLevels.h"

// EUTELESCOPE
#include "EUTelHistogramManager.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelUtility.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelPatRecTriplets.h"

// Cluster types
#include "EUTelSparseClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelReaderGenericLCIO.h"

namespace eutelescope {

	/** \class EUTelProcessorPatRecTriplets Pattern recognition processor
	 * 
     *  If compiled with MARLIN_USE_AIDA 
     *  This will do pattern recognition using a collection of LCIO::TRACKERHIT objects.
     *  Triplets are formed from each arm of the telescope using some cuts DEPENDENT on geometry.
     *  An initial track prediction is calculated using the 6 hits deemed a track.
     *  Each DUT hit created is then associated to the closest track for each DUT. No cut applyed here.
     *  MUST EXCLUDE A SENSOR IF ANOTHER IS SIDE BY SIDE WITH IT. 
	 */
	class EUTelProcessorPatRecTriplets : public marlin::Processor {
       
    public:

    virtual marlin::Processor* newProcessor() {
			return new EUTelProcessorPatRecTriplets;
		}
        /// Create object to begin triplet creation. This can be used from any processor.
		EUTelProcessorPatRecTriplets();

        /// Init objects for use in each event. Done once per jobsub.
		virtual void init();
        /// Initialise every run
		virtual void processRunHeader(lcio::LCRunHeader* run);
        /// Run each event.
		virtual void processEvent(lcio::LCEvent * evt);
        /// Performed at the end of each event.
		virtual void check(lcio::LCEvent * evt);
        ///called at the end of every run.
		virtual void end();

		/// Histogram booking 
		void bookHistograms();
        ///Place tracks in this function to plot the state parameters.
		void plotHistos( std::vector<EUTelTrack>&  trackCandidates );
        ///Use in each event to save these tracks to the lcio event.
		void outputLCIO(LCEvent* evt,std::vector< EUTelTrack >&);
		 /** AIDA histogram map
		 *  Instead of putting several pointers to AIDA histograms as
		 *  class members, histograms are booked in the init() method and
		 *  their pointers are inserted into this map keyed by their
		 *  names.
		 *  The histogram filling can proceed recalling an object through
		 *  its name
		 */
		std::map< std::string, AIDA::IHistogram1D* > _aidaHistoMap1D;

		struct _histName {
			static std::string _numberTracksCandidatesHistName;
			static std::string _HitOnTrackCandidateHistName;
			static std::string _HitOnTrackTimeHistName;
			static std::string _chi2CandidateHistName;
		};


    protected:
		/// Input TrackerHit collection name 
		std::string _hitInputCollectionName;

		/// Output Tracker collection name 
		std::string _trackCandidateHitsOutputCollectionName;

		/// Histogram info file name 
		std::string _histoInfoFileName;

		/// Track fitter
		EUTelPatRecTriplets* _trackFitter;

		/// Initial cut to create 2 hits associated.
        std::vector<float> _doubletDistCut;
		std::vector<float> _localDistDUT;
        /// Connect the central point to these two hits from the doublet.
        std::vector<float> _doubletCenDistCut;
        /// Connect triplets if under this slope cut comparision. 
        std::vector<float> _tripletSlopeCuts;
        ///Connect triplets if extrapolation under this cut in location prediction
        std::vector<float> _tripletConnectDistCut;
        ///Min hits needed.
        int _minHits;
        ///the dimention vector.
        IntVec _planeDimension;
		/// Number of events processed
		int _nProcessedRuns;
		/// Number of runs processed 
		int _nProcessedEvents;
        /// The number of skiped events due to data missing
        unsigned int _dataMissNumber;
        /// Initial displacement to estimate initial incidence at plane 0
		float _initialDisplacement;        
		/// Beam energy in [GeV]
		double _eBeam;

		/// Beam charge in [e] 
		double _qBeam;
	    /// Planes we have excluded.	
		EVENT::IntVec _excludePlanes;         
        int _mode;
        int _hitNum;
        double _dutDistCut;
		private:

    AIDA::IHistogram1D * _DoubletXseperationHistoRight;
    AIDA::IHistogram1D * _DoubletYseperationHistoRight;
    AIDA::IHistogram1D * _DoubletXseperationHistoLeft;
    AIDA::IHistogram1D * _DoubletYseperationHistoLeft;
    AIDA::IHistogram1D * _TripletXseperationHistoRight;
    AIDA::IHistogram1D * _TripletYseperationHistoRight;
    AIDA::IHistogram1D * _TripletXseperationHistoLeft;
    AIDA::IHistogram1D * _TripletYseperationHistoLeft;
    AIDA::IHistogram1D * _TripletDistCutXHisto;
    AIDA::IHistogram1D * _TripletDistCutYHisto;
    AIDA::IHistogram1D * _TripletSlopeHistoX ;
    AIDA::IHistogram1D * _TripletSlopeHistoY ;
    AIDA::IHistogram1D * _DUTWindowHisto;

        /// prevent users from making (default) copies of processors
		DISALLOW_COPY_AND_ASSIGN(EUTelProcessorPatRecTriplets)   
     
	};

	/** A global instance of the processor */
	EUTelProcessorPatRecTriplets gEUTelProcessorPatRecTriplets;

} // eutelescope

#endif 
