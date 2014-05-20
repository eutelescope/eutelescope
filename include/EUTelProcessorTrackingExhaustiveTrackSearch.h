// Version: $Id$
#ifndef EUTelTrackingExhaustiveTrackSearch_h
#define EUTelTrackingExhaustiveTrackSearch_h 1

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
#include "EUTelTrackFinder.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelUtility.h"

namespace eutelescope {

    /** @class EUTelProcessorTrackingExhaustiveTrackSearch Pattern recognition processor
     * 
     *  This processor performs track pattern recognition step that is 
     *  necessary for two-step (track-search/fit) tracking, for example GBL.
     *  Pattern recognition is based on brute-force enumeration of
     *  all possible hit combinations that fulfill some requirements (optional);
     * 
     *  @see EUTelTrackFinder
     */
  class EUTelProcessorTrackingExhaustiveTrackSearch : public marlin::Processor {
    private:
        DISALLOW_COPY_AND_ASSIGN(EUTelProcessorTrackingExhaustiveTrackSearch)   // prevent users from making (default) copies of processors
        
    public:

    virtual marlin::Processor* newProcessor() {
            return new EUTelProcessorTrackingExhaustiveTrackSearch;
        }

        EUTelProcessorTrackingExhaustiveTrackSearch();

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


        // Processor parameters        

    public:
        /** Retuns amount of missing hits in a track candidate */
        int getAllowedMissingHits() const;

        /** Retuns maximal amount of track candidates per event */
        int getMaxTrackCandidates() const;

        /** Retuns maximal amount of track candidates per event */
        int getFinderMode() const;



    public:
        /** Fills hits data structure for track finder */
        int FillHits(lcio::LCEvent*, lcio::LCCollection*,
		     std::map< int, EVENT::TrackerHitVec >&, std::vector< EVENT::TrackerHitVec >&) const;

        /** Prepare LCIO data structure for dumping track
         * candidate hits into LCIO files
         */
        void addTrackCandidateToCollection(lcio::LCEvent*, const std::vector< EVENT::TrackerHitVec >&);


    protected:

        // Input/Output collections of the processor

        /** Input TrackerHit collection name */
	std::string _hitInputCollectionName;

        /** Output TrackerHit collection name */
	std::string _trackCandidateHitsOutputCollectionName;


    protected:
        /** Hot pixel collection name 
         *  @TODO must be a part of a separate data structure
         */
//        string _hotPixelCollectionName;

        /** Map of hot pixels */
//        map<string, bool > _hotPixelMap;

    protected:

        /** Track finder*/
        EUTelTrackFinder *_trackFinder;


    private:
        
        /** TGeo geometry file name */
	std::string _tgeoFileName;

        // Exhaustive finder state definition

        /** Maximal amount of missing hits per track candidate */
        int _maxMissingHitsPerTrackCand;

        /** Maximal amount of tracks per event */
        int _maxNTracks;

        /** Finder mode. Defines search algorithm */
        int _finderMode;

        lcio::FloatVec _residualsXMin; /** Maximal negative deviation in X direction */
        lcio::FloatVec _residualsYMin; /** Maximal negative deviation in Y direction */
        lcio::FloatVec _residualsXMax; /** Maximal positive deviation in X direction */
        lcio::FloatVec _residualsYMax; /** Maximal positive deviation in Y direction */
        lcio::FloatVec _residualsRMax; /** Maximal distance in XY plane */


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
	  static std::string _numberTracksCandidatesHistName;
	  static std::string _numberOfHitOnTrackCandidateHistName;
        };

#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

    };

    /** A global instance of the processor */
    EUTelProcessorTrackingExhaustiveTrackSearch gEUTelProcessorTrackingExhaustiveTrackSearch;

} // eutelescope

#endif // EUTelTrackingExhaustiveTrackSearch_h
