
#ifndef EUTELESCOPEPROCESSORTRACKINGGBLTRACKFIT_H
#define	EUTELESCOPEPROCESSORTRACKINGGBLTRACKFIT_H


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

#ifdef USE_GBL

// GBL
#include "include/MilleBinary.h"

// EUTELESCOPE
#include "EUTelTrackFitter.h"
#include "EUTelGBLFitter.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

using namespace lcio;
using namespace marlin;
using namespace std;

namespace eutelescope {

    /** @class EUTelProcessorTrackingGBLTrackFit 
     * 
     *  @see EUTelTrackFitter
     */
    class EUTelProcessorTrackingGBLTrackFit : public Processor {
    public:

        virtual Processor* newProcessor() {
            return new EUTelProcessorTrackingGBLTrackFit;
        }

        EUTelProcessorTrackingGBLTrackFit();

        /** Called at the begin of the job before anything is read.
         * Use to initialize the processor, e.g. book histograms.
         */
        virtual void init();

        /** Called for every run.
         */
        virtual void processRunHeader(LCRunHeader* run);

        /** Called for every event - the working horse.
         */
        virtual void processEvent(LCEvent * evt);

        virtual void check(LCEvent * evt);

        /** Called after data processing for clean up.
         */
        virtual void end();



    public:
        /** Histogram booking */
        void bookHistograms();


        // Processor parameters

    public:
        
        /** Mille binary filename */
        string _binaryFilename;

        /** TGeo geometry file name */
        string _tgeoFileName;
        
    protected:

        // Input/Output collections of the processor

        /** Input TrackerHit collection name */
        string _trackCandidateHitsInputCollectionName;

        /** Output Tracks collection name */
        string _tracksOutputCollectionName;

    protected:

        /** Track fitter */
        EUTelTrackFitter *_trackFitter;
        
        /** Mille */
        gbl::MilleBinary * _milleGBL;


    protected:

        // Geometry related information

        /** GEAR description */
        EUTelGeometryTelescopeGeoDescription* _geometry;


    private:


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

        map< string, AIDA::IHistogram1D* > _aidaHistoMap1D;
        map< string, AIDA::IHistogram2D* > _aidaHistoMap2D;

        /** Names of histograms */
        static string _chi2GblFitHistName;
        static string _probGblFitHistName;
        static string _residGblFitHistName;
        static string _residGblFitHistNameX;
        static string _residGblFitHistNameY;
        static string _resid2DGblFitHistNameXvsX;
        static string _resid2DGblFitHistNameXvsY;
        static string _resid2DGblFitHistNameYvsX;
        static string _resid2DGblFitHistNameYvsY;
        static string _kinkGblFitHistNameX;
        static string _kinkGblFitHistNameY;

#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

    };

    /** A global instance of the processor */
    EUTelProcessorTrackingGBLTrackFit gEUTelProcessorTrackingGBLTrackFit;

} // eutelescope

#endif // USE_GBL

#endif	/* EUTELESCOPEPROCESSORTRACKINGGBLTRACKFIT_H */

