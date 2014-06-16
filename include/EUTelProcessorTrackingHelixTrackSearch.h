// Version: $Id$
#ifndef EUTelTrackingHelixTrackSearch_h
#define EUTelTrackingHelixTrackSearch_h 1

// C++
#include <string>

// LCIO
#include "lcio.h"

#include "marlin/Processor.h"

#include "IMPL/TrackerHitImpl.h"

// AIDA
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#endif

// EUTELESCOPE
#include "EUTelTrackFitter.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelUtility.h"
#include "EUTelTrackImpl.h"

class EUTelTrackImpl;

namespace eutelescope {

    /** @class EUTelProcessorTrackingHelixTrackSearch Pattern recognition processor
     * 
     *  This processor performs track pattern recognition step that is 
     *  necessary for two-step (track-search/fit) tracking, for example GBL.
     *  Pattern recognition is based on brute-force enumeration of
     *  all possible hit combinations that fulfill some requirements (optional);
     * 
     *  @see EUTelTrackFinder
     */
  class EUTelProcessorTrackingHelixTrackSearch : public marlin::Processor {
    private:
        DISALLOW_COPY_AND_ASSIGN(EUTelProcessorTrackingHelixTrackSearch)   // prevent users from making (default) copies of processors
        
    public:

    virtual marlin::Processor* newProcessor() {
            return new EUTelProcessorTrackingHelixTrackSearch;
        }

        EUTelProcessorTrackingHelixTrackSearch();

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



    public:
        /** Histogram booking */
        void bookHistograms();

        void plotHistos( vector< EUTelTrackImpl* >&  trackCandidates );

				void outputLCIO(LCEvent* evt,std::vector< EUTelTrackImpl* >&);

				void cartesian2LCIOTrack( EUTelTrackImpl* track, IMPL::TrackImpl*) const ;


        // Processor parameters        

    public:
        /** Retuns amount of missing hits in a track candidate */
        int getAllowedMissingHits() const;

        /** Retuns maximal amount of track candidates per event */
        int getMaxTrackCandidates() const;


    public:
        /** Fills hits data structure for track finder */
        void FillHits(LCEvent*, LCCollection*, EVENT::TrackerHitVec&) const;
        
        /** Prepare LCIO data structure for dumping track
         * candidate hits into LCIO files
         */
        void addTrackCandidateToCollection(lcio::LCEvent*, const std::vector< EVENT::TrackerHitVec >&);


        /** Prepare LCIO data structure for dumping track hits from PR search
         *  into LCIO files
         */
        void addTrackCandidateHitFittedToCollection(LCEvent*,  EVENT::TrackerHitVec& );


        /** Prepare LCIO data structure for dumping track
         * candidate hits into LCIO files
         */
        void addTrackCandidateToCollection1(lcio::LCEvent* evt, std::vector< IMPL::TrackImpl* >&);

    protected:

        // Input/Output collections of the processor

        /** Input ZS Data collection name */
	std::string _zsDataCollectionName;

        /** Input TrackerHit collection name */
	std::string _hitInputCollectionName;

        /** Output TrackerHit collection name */
        string _hitFittedOutputCollectionName;

        /** Output TrackerHit collection name */
        string _trackCandidateHitsOutputCollectionName;



    protected:
        /** Hot pixel collection name 
         *  @TODO must be a part of a separate data structure
         */
//        string _hotPixelCollectionName;

        /** Map of hot pixels */
//        map<string, bool > _hotPixelMap;

    protected:

        /** Track fitter*/
        EUTelTrackFitter* _trackFitter;


    private:

        /** TGeo geometry file name */
				std::string _tgeoFileName;
        
        // Exhaustive finder state definition

        /** Maximal amount of missing hits per track candidate */
        int _maxMissingHitsPerTrackCand;

				//The number of allowed similar hits on track candidates of a single event
				int _AllowedSharedHitsOnTrackCandidate;

        /** Maximal amount of tracks per event */
        int _maxNTracks;
 
        /** Number of Planes to project starting with the first most upstream */
        int _planesProject;
        
        /** Maximal distance in XY plane */
        double _residualsRMax;
        
        /** Beam energy in [GeV] */
        double _eBeam;

        /** Beam charge in [e] */
        double _qBeam;
        
        /** Beam energy uncertainty [%] */
        double _eBeamUncertatinty;
        
        /** Beam spread in [mrad] */
        EVENT::FloatVec _beamSpread;

        /** Histogram info file name */
	std::string _histoInfoFileName;

    protected:

        // Statistics counters

        /** Number of events processed */
        int _nProcessedRuns;
        /** Number of runs processed */
        int _nProcessedEvents;

    public:
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
            static string _numberTracksCandidatesHistName;
            static string _numberOfHitOnTrackCandidateHistName;
            static string _HitOnTrackCandidateHistName;
						static string _chi2CandidateHistName;

        };

#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

    };

    /** A global instance of the processor */
    EUTelProcessorTrackingHelixTrackSearch gEUTelProcessorTrackingHelixTrackSearch;

} // eutelescope

#endif // EUTelTrackingExhaustiveTrackSearch_h
