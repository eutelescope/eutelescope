// Version: $Id$
#ifndef EUTelProcessorPatRecTriplets_h
#define EUTelProcessorPatRecTriplets_h 1

// AIDA
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/IHistogramFactory.h>
#endif

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

	/** @class EUTelProcessorPatRecTriplets Pattern recognition processor
	 * 
	 *  This processor performs track pattern recognition step that is 
	 *  necessary for two-step (track-search/fit) tracking, for example GBL.
	 *  Pattern recognition uses simple equation of motion to estimate possible trajectory. 
	 */
	class EUTelProcessorPatRecTriplets : public marlin::Processor {
       
    public:

    virtual marlin::Processor* newProcessor() {
			return new EUTelProcessorPatRecTriplets;
		}

		EUTelProcessorPatRecTriplets();

		//Here we define the virtual function from marlin::Processor.
		/** Called at the begin of the job before anything is read.
		 * Use to initialize the processor, e.g. book histograms.
		 */
		virtual void init();

		/** Called for every run.
		 */
		virtual void processRunHeader(lcio::LCRunHeader* run);

		/** Called for every event - the working horse.
		 */
		virtual void processEvent(lcio::LCEvent * evt);

		virtual void check(lcio::LCEvent * evt);

		/** Called after data processing for clean up.
		 */
		virtual void end();

    //This is the function declarations.//////////////////////////////////////
		/** Histogram booking */
		void bookHistograms();

		void plotHistos( std::vector<EUTelTrack>&  trackCandidates );

		void outputLCIO(LCEvent* evt,std::vector< EUTelTrack >&);
		//Here we define all things to do with histogramming.
		#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
		 /** AIDA histogram map
		 *  Instead of putting several pointers to AIDA histograms as
		 *  class members, histograms are booked in the init() method and
		 *  their pointers are inserted into this map keyed by their
		 *  names.
		 *  The histogram filling can proceed recalling an object through
		 *  its name
		 */
		std::map< std::string, AIDA::IHistogram1D* > _aidaHistoMap1D;

		/** Names of histograms */
		struct _histName {
			static std::string _numberTracksCandidatesHistName;
			static std::string _numberOfHitOnTrackCandidateHistName;
			static std::string _HitOnTrackCandidateHistName;
			static std::string _chi2CandidateHistName;
		};
		#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)


    protected:
		//This is variables being declared here/////////////////////
		/** Input TrackerHit collection name */
		std::string _hitInputCollectionName;

		/** Output Tracker collection name */
		std::string _trackCandidateHitsOutputCollectionName;

		/* Histogram info file name */
		std::string _histoInfoFileName;

		/** Track fitter*/
		EUTelPatRecTriplets* _trackFitter;

        std::vector<float> _doubletDistCut;
        std::vector<float> _tripletSlopeCuts;
        std::vector<float> _tripletConnectDistCut;
        std::vector<float> _doubletCenDistCut;

		/** Maximal amount of missing hits per track candidate */
		/** Number of events processed */
		int _nProcessedRuns;
		/** Number of runs processed */
		int _nProcessedEvents;
        /** The number of skiped events due to data missing */
        unsigned int _dataMissNumber;
		/** Maximal amount of tracks per event */
		int _maxNTracks;

		float _initialDisplacement;        

		/** Maximal distance in XY plane */
		double _residualsRMax;
		
		/** Beam energy in [GeV] */
		double _eBeam;

		/** Beam charge in [e] */
		double _qBeam;
		
		EVENT::IntVec _createSeedsFromPlanes;
		EVENT::IntVec _excludePlanes;         
		EVENT::IntVec _planeDimension;

		private:
		DISALLOW_COPY_AND_ASSIGN(EUTelProcessorPatRecTriplets)   // prevent users from making (default) copies of processors
     
	};

	/** A global instance of the processor */
	EUTelProcessorPatRecTriplets gEUTelProcessorPatRecTriplets;

} // eutelescope

#endif 