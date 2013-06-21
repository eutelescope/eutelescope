
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
#include "EUTelUtility.h"

using namespace lcio;
using namespace marlin;
using namespace std;

namespace eutelescope {

    /** @class EUTelProcessorTrackingGBLTrackFit 
     * 
     *  @see EUTelTrackFitter
     */
    class EUTelProcessorTrackingGBLTrackFit : public Processor {

    private:
        DISALLOW_COPY_AND_ASSIGN(EUTelProcessorTrackingGBLTrackFit)     // prevent users from making (default) copies of processors
        
    public:

        virtual Processor* newProcessor() {
            return new EUTelProcessorTrackingGBLTrackFit;
        }

        EUTelProcessorTrackingGBLTrackFit();
        
    public:
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

    private:
        /** Prepare MILLEPEDE steering file */
        void writeMilleSteeringFile();

        /** Rename MILLEPEDE result file */
        void moveMilleResultFile( const string&, const string& );
        
        /** Parse MILLEPEDE result file
         * and create output collection
         */
        bool parseMilleOutput( const string& );
        
        /** Run pede executable to solve alignment problem */
        void runPede();

        /** Generate MILLEPEDE parameter labels */
        void fillMilleParametersLabels();
       
    public:
        /** Histogram booking */
        void bookHistograms();

        
        // Processor parameters

    public:

        // Necessary parameters

        /** Beam energy in [GeV] */
        double _eBeam;


        // Optional parameters

        /** Alignment mode */
        int _alignmentMode;
      
        /** GBL M-estimator option */
        string _mEstimatorType;
        
        /** Parameter ids */
        map<int,int> _xShiftsMap;
        
        /** Parameter ids */
        map<int,int> _yShiftsMap;
        
        /** Parameter ids */
        map<int,int> _zShiftsMap;
        
        /** Parameter ids */
        map<int,int> _xRotationsMap;
        
        /** Parameter ids */
        map<int,int> _yRotationsMap;
        
        /** Parameter ids */
        map<int,int> _zRotationsMap;
        
        /** Mille binary filename */
        string _milleBinaryFilename;

        /** Mille steering filename */
        string _milleSteeringFilename;

        /** Mille result filename */
        string _milleResultFileName;
        
        /** Allows user-added commands in the pede steering file */
        StringVec _pedeSteerAddCmds;
        
        /** Alignment plane ids*/
        IntVec _alignmentPlaneIds;
        
        /** Alignment X shift plane ids to be fixed */
        IntVec _fixedAlignmentXShfitPlaneIds;
        
        /** Alignment Y shift plane ids to be fixed */
        IntVec _fixedAlignmentYShfitPlaneIds;
        
        /** Alignment Z shift plane ids to be fixed */
        IntVec _fixedAlignmentZShfitPlaneIds;
        
        /** Alignment X rotation plane ids to be fixed */
        IntVec _fixedAlignmentXRotationPlaneIds;
        
        /** Alignment Y rotation plane ids to be fixed */
        IntVec _fixedAlignmentYRotationPlaneIds;
        
        /** Alignment Z rotation plane ids to be fixed */
        IntVec _fixedAlignmentZRotationPlaneIds;
        
        /** Automatic pede run flag*/
        bool _runPede;
        
        /** Alignment constants file name */
        string _alignmentConstantLCIOFile;
        
        /** Maximum value of track chi2 for millipede */
        double _maxChi2Cut;

        /** TGeo geometry file name */
        string _tgeoFileName;
        
        /** Histogram info file name */
        string _histoInfoFileName;

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

    private:

        struct SeedAlignmentConstants {
            SeedAlignmentConstants() : _xResiduals(), _nxResiduals(), _yResiduals(), _nyResiduals() {};
            std::map< int, double > _xResiduals;   //! sum all x residuals for given plane id
            std::map< int, int >    _nxResiduals;  //! number of residuals used to calculate mean for given plane id
            std::map< int, double > _yResiduals;   //! sum all y residuals for given plane id
            std::map< int, int >    _nyResiduals;  //! number of residuals used to calculate mean for given plane id
        };

        /** Initial alignment constants */
        SeedAlignmentConstants _seedAlignmentConstants;


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
        struct _histName {
            static string _chi2GblFitHistName;
            static string _probGblFitHistName;
            static string _residGblFitHistName;
            static string _normResidGblFitHistName;
            static string _residGblFitHistNameX;
            static string _residGblFitHistNameY;
            static string _normResidGblFitHistNameX;
            static string _normResidGblFitHistNameY;
            static string _resid2DGblFitHistNameXvsX;
            static string _resid2DGblFitHistNameXvsY;
            static string _resid2DGblFitHistNameYvsX;
            static string _resid2DGblFitHistNameYvsY;
            static string _normResid2DGblFitHistNameXvsX;
            static string _normResid2DGblFitHistNameXvsY;
            static string _normResid2DGblFitHistNameYvsX;
            static string _normResid2DGblFitHistNameYvsY;
            static string _kinkGblFitHistNameX;
            static string _kinkGblFitHistNameY;
        };

#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

    };

    /** A global instance of the processor */
    EUTelProcessorTrackingGBLTrackFit gEUTelProcessorTrackingGBLTrackFit;

} // eutelescope

#endif // USE_GBL

#endif	/* EUTELESCOPEPROCESSORTRACKINGGBLTRACKFIT_H */

